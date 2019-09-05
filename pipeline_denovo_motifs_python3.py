
"""=================
Denovo Motif calling
=================

:Author: Ian Sudbery and Jaime Alvarez-Benayas
:Release: $Id: pipeline_denovo_motifs.py v2 $
:Date: |today|
:Tags: Python

The de novo motifs pipeline takes a several set of :term:`bed` formatted
genomic intervals and annotates them attempts to identify enriched 
motifs within them. 

Motifs can be discoverd by one or several of:
* MEME - with a standard background model (meme) or a background model tailored to
  a background interval set (disc_meme)
* DREME, run either with a random background set (rand_dreme) or with a 
  negative interval set (disc_dreme)
* MEME-CHIP

Once motifs have been called they will be processed in several ways: 
* Motifs from each sample and approach combination will be clustered by
  comparing each motif found to all the others using TomTom and the seed
  motifs identified 
* Motifs identified in different samples will be compared to each other
  to identify motifs that are found in mulitple samples
* Motifs identified from different methods will be compared to find if
  the different methods result in the same motifs.
 
Implemented by not tested:
MAST - Motif alignment and search tool: searches sequences for predefined 
       sequences and sorts sequences by number of matches

Not Currently Implemented
* Bioprospector
* FIFO - Find sequences that contain matches to motifs
* GLAM2 - Find gapped motifs
* GLAM2Scan - Search for gapped motifs.

Usage 
=====

See :ref:`PipelineSettingUp` and :ref:`PipelineRunning` on general
information how to use CGAT pipelines.

Configuration
-------------

The pipeline requires a configured :file:`pipeline.ini` file. The
pipeline looks for a configuration file in several places:

   1. The default configuration in the :term:`code directory`.
   2. A shared configuration file :file:`../pipeline.ini`.
   3. A local configuration :file:`pipeline.ini`.

The order is as above. Thus, a local configuration setting will
override a shared configuration setting and a default configuration
setting.

Configuration files follow the ini format (see the python
`ConfigParser <http://docs.python.org/library/configparser.html>`
documentation).  The configuration file is organized by section and
the variables are documented within the file. In order to get a local
configuration file in the current directory, type::

    python <codedir>/pipeline_intervals.py config

The following sections and parameters probably should be changed from
the default values:

.. todo::
   describe important parameters

The sphinxreport report requires a :file:`conf.py` and
:file:`sphinxreport.ini` file (see :ref:`PipelineReporting`). To start
with, use the files supplied with the Example_ data.


Input
-----

Intervals
+++++++++

Input are :term:`bed`-formatted files of intervals. Intervals should
be at least bed4 formatted, i.e., each interval should be labelled
(uniquely).

Optional inputs
+++++++++++++++

optionally, bam files with the same names as the intervals can be 
provided. This will allow the quantification of the peaks to aid
in identifying the centre of the peaks and the most significant
peaks.

Reference motifs
++++++++++++++++

Reference motifs are described by fasta-formatted multiple alignments, see for
example Jaspar download. The motifs are build by running MEME on the file.

Reference motifs should end in the suffix ".motif.fasta", for example,
:file:`rxrvdr.motif.fasta`.

Requirements
------------

The pipeline requires the information from the following pipelines:

:doc:`pipeline_annotations`
   set the configuration variable :py:data:`annotations_database`
   and :py:data:`annotations_dir`.

Requirements:

* bedtools >= 2.21.0
* samtools >= 1.1


Pipline Output
==============

The results of the computation are all stored in an sqlite relational
database :file:`csvdb`.

Example
=======

Example data is available at
http://www.cgat.org/~andreas/sample_data/pipeline_chipseq.tgz.  To run
the example, simply unpack and untar::

   wget http://www.cgat.org/~andreas/sample_data/pipeline_chipseq.tgz
   tar -xvzf pipeline_chipseq.tgz
   cd pipeline_chipseq
   python <srcdir>/pipeline_chipseq.py make full

.. note::

   For the pipeline to run, install the :doc:`pipeline_annotations` as well.

Change Log
==========

26.2.2014 Andreas Heger
    added tag counting in intervals.

11.3.2014 Andreas Heger
    preprocess all intervals to make sure they are
    properly formatted:
    1. merge overlapping intervals
    2. sort intervals by coordinate

    As a consequence of this change, intervals and subsequent files
    will end up in the intervals.dir directory. Thus, existing
    pipelines will re-run.

Code
====

"""
import sys
import re
import glob
import os
import sqlite3
import pysam

from ruffus import *
from ruffus.combinatorics import *

from CGATCore import Experiment as E
from CGATCore import Pipeline as P
import CGATCore.IOTools as IOTools
import CGAT.Bed as Bed
import CGATPipelines.PipelinePeakcalling as PipelinePeakcalling
import PipelineDeNovoMotifs_python3 as PipelineMotifs
import CGATPipelines.PipelineTracks as PipelineTracks


###################################################
###################################################
###################################################
# Pipeline configuration
###################################################
P.get_parameters(
    ["%s/pipeline.yml" % os.path.splitext(__file__)[0],
     "../pipeline.yml",
     "pipeline.yml"],
    defaults={
        'paired_end': False})

PARAMS = P.PARAMS


PipelinePeakcalling.PARAMS = PARAMS
PipelineMotifs.PARAMS = PARAMS


###################################################################
###################################################################
###################################################################
# Helper functions mapping tracks to conditions, etc
###################################################################
# load all tracks - exclude input/control tracks
# determine the location of the input files (reads).
DATADIR = PARAMS.get('input', '.')
if not os.path.exists(DATADIR):
    raise OSError('data directory %s does not exists')

Sample = PipelineTracks.Sample

TRACKS = PipelineTracks.Tracks(Sample).loadFromDirectory(
    glob.glob(os.path.join(DATADIR, "*.bed.gz")),
    "(\S+).bed.gz")

BEDFILES = [os.path.join(
    DATADIR, "%s.bed.gz") % x for x in TRACKS]


# create an indicator target
@transform(BEDFILES, suffix(".gz"), ".gz")
def BedFiles(infile, outfile):
    pass


BAMFILES = glob.glob(os.path.join(DATADIR, "*.bam"))


def getAssociatedBAMFiles(track):
    '''return a list of BAM files associated with a track.

    By default, this method searches for ``track.bam`` file in the
    data directory and returns an offset of 0.

    Associations can be defined in the .ini file in the section
    [bams]. For example, the following snippet associates track
    track1 with the bamfiles :file:`track1.bam` and :file:`track2.bam`::

       [bams]
       track1=track1.bam,track2.bam

    Glob expressions are permitted.

    Offsets are used to shift tags in ChIP experiments. Offsets
    need to be defined in the [offsets] sections. If no offsets
    are defined, the method returns a list of 0 offsets.

    Offsets need to be defined in the same order as the bam files::

       [offsets]
       track1=120,200

    returns a list of BAM files and offsets.

    Default tracks and offsets can be specified using a placeholder ``%``. The
    following will associate all tracks with the same bam file::

        [bams]
        %=all.bam

    '''
    fn = os.path.basename(track.asFile())
    bamfiles = glob.glob("%s.bam" % fn)

    if bamfiles == []:
        if "bams_%s" % fn.lower() in PARAMS:
            for ff in P.as_list(PARAMS["bams_%s" % fn.lower()]):
                bamfiles.extend(glob.glob(ff))
#        else:
#            for pattern, value in P.CONFIG.items("bams"):
#                if "%" in pattern:
#                    p = re.sub("%", "\S+", pattern.lower())
#                    if re.search(p, fn.lower()):
#                        bamfiles.extend(glob.glob(value))

    offsets = []
    if "offsets_%s" % fn.lower() in PARAMS:
        offsets = map(int, P.as_list(PARAMS["offsets_%s" % fn.lower()]))
#    else:
#        for pattern, value in P.CONFIG.items("offsets"):
#            if "%" in pattern:
#                p = re.sub("%", "\S+", pattern)
#                if re.search(p, fn):
#                    offsets.extend(map(int, value.split(",")))

    if offsets == []:
        offsets = [0] * len(bamfiles)

    if len(bamfiles) != len(offsets):
        raise ValueError("number of BAM files %s is not the "
                         "same as number of offsets: %s" % (
                             str(bamfiles), str(offsets)))

    return bamfiles, offsets


def connect():
    '''connect to database.

    This method also attaches to helper databases.
    '''

    dbh = sqlite3.connect(re.sub("sqlite:///","", PARAMS["database_url"]))
    statement = '''ATTACH DATABASE '%s' as annotations''' % (
        PARAMS["annotations_database"])
    cc = dbh.cursor()
    cc.execute(statement)
    cc.close()

    return dbh


@transform('preprocess.dir/*.bed.gz',
           regex(".*/(.*).bed.gz"),
           r"%s/\1.bed.gz" % DATADIR)
def prepareIntervals(infile, outfile):
    '''optionally preprocess intervals.

    1. Sort intervals by position.
    2. Merge overlapping intervals.
    3. Rename intervals with a sequential number.
    '''

    statement = '''
    zcat %(infile)s
    | sort -k1,1 -k2,2n
    | bedtools merge
    | awk '{printf("%%s\\t%%s\\t%%s\\t%%i\n",$1,$2,$3,++a)}'
    | bgzip
    > %(outfile)s;
    tabix -p bed %(outfile)s
    '''
    P.run(statement)


@merge(prepareIntervals, "preprocess.tsv.gz")
def buildProcessingSummary(infiles, outfile):

    outf = IOTools.open_file(outfile, "w")

    for infile in infiles:
        before = os.path.join(
            'preprocess.dir',
            getNumLines(os.path.basename(infile)))
        after = getNumLines(infile)
        outf.write("%s\t%i\t%i\n" % (
            P.snip(os.path.basename(infile), ".bed.gz"),
            before,
            after))

    outf.close()


@follows(mkdir(os.path.join(PARAMS["exportdir"], "bed")))
@transform(BEDFILES,
           regex(r".*/(.*).bed.gz"),
           os.path.join(PARAMS["exportdir"], "bed", r"\1.bed.gz"))
def indexIntervals(infile, outfile):
    '''index intervals with tabix.
    '''
    statement = '''zcat %(infile)s
    | sort -k1,1 -k2,2n
    | bgzip > %(outfile)s;
    tabix -p bed %(outfile)s'''

    P.run(statement)


@transform(BEDFILES,
           suffix(".bed.gz"),
           "_intervals.load")
def loadIntervals(infile, outfile):
    '''load intervals from :term:`bed` formatted files into
    the database.

    If a :term:`bam` file is associated with a :term:`bed`
    file, re-evaluate the intervals by counting reads within
    the interval. In contrast to the initial pipeline, the
    genome is not binned.

       nprobes: number of reads in interval
       peakcenter: position with maximum number of reads in interval
       avgval: average coverage within interval
    '''

    tmpfile = P.get_temp_file(".")

    headers = ("avgval", "disttostart",
               "genelist", "length",
               "peakcenter", "peakval",
               "position", "interval_id",
               "npeaks", "nprobes",
               "contig", "start", "end", "score", "strand")

    tmpfile.write("\t".join(headers) + "\n")

    (avgval, contig, disttostart, end, genelist,
     length, peakcenter, peakval, position,
     start, interval_id, npeaks, nprobes) = \
        0, "", 0, 0, "", 0, 0, 0, 0, 0, 0, 0, 0

    track = Sample(filename=P.snip(infile, ".bed.gz"))

    bamfiles, offsets = getAssociatedBAMFiles(track)

    if bamfiles:
        E.info("%s: associated bamfiles = %s" % (track, bamfiles))
    else:
        E.info("%s: no bamfiles associated" % (track))

    # open all bamfiles
    samfiles = [pysam.Samfile(fn, "rb") for fn in bamfiles]

    c = E.Counter()

    # count tags
    for bed in Bed.iterator(IOTools.open_file(infile, "r")):

        c.input += 1

        if "name" not in bed:
            bed.name = c.input

        try:
            strand = bed["strand"]
        except IndexError:
            strand = "."
            
        # The fifth field of a bed file can be used to supply a
        # score. Our iterator returns the optional fields as a "fields
        # array". The first of these is the interval name, and the
        # second the score. The score may be more is better or less is
        # better.
        if len(bed.fields) > 1:
            value = bed.fields[1]
            if value != "":
                score = value
            else:
                score = 1
        else:
            score = 1

        if samfiles:
            npeaks, peakcenter, length, avgval, peakval, nprobes = \
                PipelinePeakcalling.countPeaks(
                    bed.contig,
                    bed.start,
                    bed.end,
                    samfiles,
                    offsets)
            if nprobes == 0:
                c.skipped_reads += 1

        else:
            # deal with bed12
            bed_intervals = bed.toIntervals()
            length = sum([e - s for s, e in bed_intervals])
            mid_point = length / 2
            for s, e in bed_intervals:
                peakcenter = s + mid_point
                if peakcenter >= e:
                    mid_point = peakcenter - e
                else:
                    break

            npeaks, avgval, peakval, nprobes = \
                (1,
                 1,
                 1,
                 1)

        c.output += 1
        tmpfile.write("\t".join(map(
            str,
            (avgval, disttostart, genelist, length,
             peakcenter, peakval, position, bed.name,
             npeaks, nprobes,
             bed.contig, bed.start, bed.end, score, strand))) + "\n")

    if c.output == 0:
        E.warn("%s - no aggregate intervals")

    tmpfile.close()

    P.load(tmpfile.name,
           outfile,
           tablename=os.path.basename("%s_intervals" % track.asTable()),
           options="--allow-empty-file "
           "--add-index=interval_id")

    os.unlink(tmpfile.name)

    E.info("%s\n" % str(c))


@follows(mkdir(os.path.join(PARAMS["exportdir"], "peaks")))
@transform(loadIntervals,
           regex(r"(.*)_intervals.load"),
           os.path.join(PARAMS["exportdir"], "peaks", r"\1.peak.bed.gz"))
def exportPeakLocations(infile, outfile):
    '''export peak locations
    '''

    dbh = connect()
    outf = IOTools.open_file(outfile, "w")
    cc = dbh.cursor()
    table = P.to_table(infile)
    for x in cc.execute(
            """SELECT contig, peakcenter, peakcenter+1, interval_id, peakval
            FROM %(table)s """ % locals()):
        outf.write("\t".join(map(str, x)) + "\n")
    outf.close()


@follows(loadIntervals)
def buildIntervals():
    pass


def exportIntervalSequences(infile, outfile, track, method):
    '''export sequences for motif discovery.

    This method requires the _interval tables.

    For motif discovery, only the sequences with the highest S/N ratio
    are supplied.

    1. The top *motifs_proportion* intervals sorted by peakval
    2. Only a region +/- *motifs_halfwidth* around the peak
    3. At least *motifs_min_sequences*. If there are not enough sequences
          to start with, all will be used.
    4. At most *motifs_max_size* sequences will be output.

    '''
    dbhandle = connect()

    try:
        halfwidth = int(PARAMS[method+"_halfwidth"])
        full = False
    except ValueError:
        full = True
        halfwidth = None

    try:
        maxsize = int(PARAMS[method+"_max_size"])
    except ValueError:
        maxsize = None

    nseq = PipelineMotifs.writeSequencesForIntervals(
        track,
        outfile,
        dbhandle,
        full=full,
        masker=P.as_list(PARAMS[method+'_masker']),
        halfwidth=halfwidth,
        maxsize=maxsize,
        num_sequences=PARAMS[method+"_num_sequences"],
        proportion=PARAMS[method+"_proportion"],
        min_sequences=PARAMS[method+"_min_sequences"],
        order=PARAMS[method+'_score'])

    if nseq == 0:
        E.warn("%s: no sequences - %s skipped" % (outfile, method))
        P.touch(outfile)


############################################################
############################################################
############################################################
@follows(mkdir("motifs"))
@transform(BEDFILES,
           regex(".*/(.*).bed.gz"),
           r"motifs/\1.foreground.fasta")
def exportMotifIntervalSequences(infile, outfile):
    '''export sequences for motif detection.

    This method requires the _interval tables.
    '''
    PipelineMotifs.exportSequencesFromBedFile(
        infile, outfile,
        masker=PARAMS['motifs_masker'])


@follows(mkdir("motifs"))
@transform(BEDFILES,
           regex(".*/(.*).bed.gz"),
           r"motifs/\1.control.fasta")
def exportMotifControlSequences(infile, outfile):
    '''for each interval, export the left and right
    sequence segment of the same size.
    '''
    PipelineMotifs.exportSequencesFromBedFile(
        infile, outfile,
        masker=PARAMS['motifs_masker'],
        mode="leftright")


############################################################
############################################################
############################################################
@active_if("meme" in P.as_list(PARAMS["methods"]) or
           "disc_meme" in P.as_list(PARAMS["methods"]))
@transform(loadIntervals,
           suffix("_intervals.load"),
           ".meme.fasta")
def exportMemeIntervalSequences(infile, outfile):
    
    track = os.path.basename(P.snip(infile, "_intervals.load"))

    exportIntervalSequences(infile, outfile, track, "meme")


############################################################
@follows(mkdir("meme.dir"))
@active_if("meme" in P.as_list(PARAMS["methods"]))
@transform(exportMemeIntervalSequences, regex("(.+).meme.fasta"),
           r"meme.dir/\1.meme")
def runMeme(infile, outfile):
    '''run MEME to find motifs.

    In order to increase the signal/noise ratio, MEME is not run on
    all intervals but only the top 10% of intervals (peakval) are
    used.  Also, only the segment of 200 bp around the peak is used
    and not the complete interval.

    * Softmasked sequence is converted to hardmasked
      sequence to avoid the detection of spurious motifs.

    * Sequence is run through dustmasker

    '''
    PipelineMotifs.runMEMEOnSequences(infile, outfile)


@transform(runMeme,
           regex("(.+)"),
           r"\1.self_matches")
def memeSelfMatches(infile, outfile):
    '''Compare output of MEME and DREME to itself so that similar 
    motifs within runs can be clustered '''

    PipelineMotifs.tomtom_comparison(infile, infile, outfile)


############################################################
@transform(memeSelfMatches,
           regex("(.+).self_matches"),
           add_inputs(r"\1"),
           r"\1.seeds")
def getMemeSeeds(infiles, outfile):
    ''' extract seed motifs for MEME and DREME output '''

    alignments, motifs = infiles
    PipelineMotifs.getSeedMotifs(motifs, alignments, outfile)

############################################################
@permutations(getMemeSeeds,
              formatter("meme.dir/(?P<TRACK>.+).meme.seeds"),
              2,
              r"meme.dir/{TRACK[0][0]}_to_{TRACK[1][0]}.tomtom")
def compareMemeTracks(infiles, outfile):
    '''Use tomtom to look if simlar motifs have been foundin more than one
    track '''

    track1, track2 = infiles

    PipelineMotifs.tomtom_comparison(track1, track2, outfile)


############################################################
@follows(compareMemeTracks)
def meme():
    pass

############################################################
def getdescMEMEfiles():

    try:
        design = IOTools.open_file("design.tsv")
    except OSError:
        raise "A design.tsv must be supplied for descriminative MEME analysis"

    design.readline()
    for line in design:

        pos, neg = line.strip().split("\t")
        posf = "%s.meme.fasta" % pos
        negf = "%s.meme.fasta" % neg
        out = "disc_meme.dir/%s_vs_%s.psp" %(pos,neg)
        yield ((posf, negf), out)


@active_if("disc_meme" in P.as_list(PARAMS["methods"]))
@follows(exportMemeIntervalSequences, mkdir("disc_meme.dir"))
@files(getdescMEMEfiles)
def getDiscMEMEPSPFile(infiles, outfile):
    '''Get the position specific prior file to allow
    discriminative motif finding with meme '''

    pos, neg = infiles
    PipelineMotifs.generatePSP(pos, neg, outfile)


############################################################
@transform(getDiscMEMEPSPFile,
           regex("disc_meme.dir/(.+)_vs_(.+).psp"),
           add_inputs(r"\1.meme.fasta"),
           r"disc_meme.dir/\1_vs_\2.meme")
def runDiscMEME(infiles, outfile):
    ''' Run MEME with PSP file, therefore making it 
    discrimenative'''

    psp, fasta = infiles
    PipelineMotifs.runMEMEOnSequences(fasta, outfile, psp=psp)


############################################################
@transform(runDiscMEME,
           regex("(.+)"),
           r"\1.self_matches")
def discMemeSelfMatches(infile, outfile):
    '''Compare output of MEME and DREME to itself so that similar 
    motifs within runs can be clustered '''

    PipelineMotifs.tomtom_comparison(infile, infile, outfile)


############################################################
@transform(discMemeSelfMatches,
           regex("(.+).self_matches"),
           add_inputs(r"\1"),
           r"\1.seeds")
def getDiscMemeSeeds(infiles, outfile):
    ''' extract seed motifs for MEME and DREME output '''

    alignments, motifs = infiles
    PipelineMotifs.getSeedMotifs(motifs, alignments, outfile)


############################################################
@permutations(getDiscMemeSeeds,
              formatter("disc_meme.dir/(?P<TRACK>.+).meme.seeds"),
              2,
              r"disc_meme.dir/{TRACK[0][0]}_to_{TRACK[1][0]}.tomtom")
def compareDiscMemeTracks(infiles, outfile):
    '''Use tomtom to look if simlar motifs have been foundin more than one
    track '''

    track1, track2 = infiles

    PipelineMotifs.tomtom_comparison(track1, track2, outfile)


############################################################
@follows(compareDiscMemeTracks)
def discMeme():
    pass


############################################################
# DREME
############################################################
@active_if("rand_dreme" in PARAMS["methods"] or
           "disc_dreme" in PARAMS["methods"])
@transform(loadIntervals,
           suffix("_intervals.load"),
           ".dreme.fasta")
def exportDremeIntervalSequences(infile, outfile):
    
    track = os.path.basename(P.snip(infile, "_intervals.load"))

    exportIntervalSequences(infile, outfile, track, "dreme")


############################################################
@active_if("rand_dreme" in PARAMS["methods"])
@follows(mkdir("rand_dreme.dir"))
@transform(exportDremeIntervalSequences,
           regex("(.+).dreme.fasta"),
           r"rand_dreme.dir/\1.dreme")
def runRandDreme(infile, outfile):
    '''Run DREME with a randomised negative set'''

    PipelineMotifs.runDREME(infile, outfile)


############################################################
@transform(runRandDreme,
           regex("(.+)"),
           r"\1.self_matches")
def randDremeSelfMatches(infile, outfile):
    '''Compare output of MEME and DREME to itself so that similar 
    motifs within runs can be clustered '''

    PipelineMotifs.tomtom_comparison(infile, infile, outfile)


############################################################
@transform(randDremeSelfMatches,
           regex("(.+).self_matches"),
           add_inputs(r"\1"),
           r"\1.seeds")
def getRandDremeSeeds(infiles, outfile):
    ''' extract seed motifs for MEME and DREME output '''

    alignments, motifs = infiles
    PipelineMotifs.getSeedMotifs(motifs, alignments, outfile)


############################################################
@permutations(getRandDremeSeeds,
              formatter("rand_dreme.dir/(?P<TRACK>.+).dreme.seeds"),
              2,
              r"rand_dreme.dir/{TRACK[0][0]}_to_{TRACK[1][0]}.tomtom")
def compareRandDremeTracks(infiles, outfile):
    '''Use tomtom to look if simlar motifs have been foundin more than one
    track '''

    track1, track2 = infiles

    PipelineMotifs.tomtom_comparison(track1, track2, outfile)


############################################################
@follows(compareRandDremeTracks)
def randDreme():
    pass


############################################################
def getDiscDREMEFiles():

    try:
        design = IOTools.open_file("design.tsv")
    except OSError:
        raise "A design.tsv must be supplied for descriminative DREME analysis"

    design.readline()
    for line in design:

        pos, neg = line.strip().split("\t")
        posf = "%s.dreme.fasta" % pos
        negf = "%s.dreme.fasta" % neg
        out = "disc_dreme.dir/%s_vs_%s.dreme" %(pos,neg)
        yield ((posf, negf), out)

@active_if("disc_dreme" in PARAMS["methods"])
@follows(mkdir("disc_dreme.dir"), exportDremeIntervalSequences)
@files(getDiscDREMEFiles)
def runDiscDREME(infiles, outfile):
    '''Run discriminative DREME using control file speicied
    in design.tsv'''

    infile, negatives = infiles
    PipelineMotifs.runDREME(infile, outfile, neg_file=negatives)


############################################################
@transform(runDiscDREME,
           regex("(.+)"),
           r"\1.self_matches")
def discDremeSelfMatches(infile, outfile):
    '''Compare output of MEME and DREME to itself so that similar 
    motifs within runs can be clustered '''

    PipelineMotifs.tomtom_comparison(infile, infile, outfile)


############################################################
@transform(discDremeSelfMatches,
           regex("(.+).self_matches"),
           add_inputs(r"\1"),
           r"\1.seeds")
def getDiscDremeSeeds(infiles, outfile):
    ''' extract seed motifs for MEME and DREME output '''

    alignments, motifs = infiles
    PipelineMotifs.getSeedMotifs(motifs, alignments, outfile)


############################################################
@permutations(getDiscDremeSeeds,
              formatter("disc_dreme.dir/(?P<TRACK>.+).dreme.seeds"),
              2,
              r"disc_dreme.dir/{TRACK[0][0]}_to_{TRACK[1][0]}.tomtom")
def compareDiscDremeTracks(infiles, outfile):
    '''Use tomtom to look if simlar motifs have been foundin more than one
    track '''

    track1, track2 = infiles

    PipelineMotifs.tomtom_comparison(track1, track2, outfile)


############################################################
@follows(compareDiscDremeTracks)
def discDreme():
    pass


############################################################
#MEME-ChIP
############################################################
@active_if("memechip" in PARAMS["methods"])
@transform(loadIntervals,
           suffix("_intervals.load"),
           ".memechip.fasta")
def exportMemeCHiPIntervalSequences(infile, outfile):
    
    track = os.path.basename(P.snip(infile, "_intervals.load"))

    exportIntervalSequences(infile, outfile, track, "memechip")


############################################################
@active_if("memechip" in PARAMS["methods"])
@follows(mkdir("memechip.dir"))
@transform(exportMemeCHiPIntervalSequences,
           regex("./(.+).memechip.fasta"),
           add_inputs("*.motif"),
           r"memechip.dir/\1.memechip")
def runMemeChIP(infiles, outfile):
    '''Run the MEME-ChIP pipeline, optionally with 
    reference motif sequences'''

    infile, motifs = infiles[0], infiles[1:]

    PipelineMotifs.runMemeCHIP(infile, outfile, motifs)


############################################################
@transform(runMemeChIP,
           regex("memechip.dir/(.+).memechip"),
           r"memechip.dir/\1.memechip.seeds")
def getMemeChipSeedMotifs(infile, outfile):
    ''' extract the seed motifs from the MEME-ChIP output'''

    motifs = infile
    track = os.path.basename(motifs)
    alignments = os.path.join(PARAMS["exportdir"],
                              "memechip.dir",
                              track,
                              "motif_alignment.txt")
                              
    PipelineMotifs.getSeedMotifs(motifs, alignments, outfile)


############################################################
@permutations(getMemeChipSeedMotifs,
              formatter("memechip.dir/(?P<TRACK>.+).memechip.seeds"),
              2,
              r"memechip.dir/{TRACK[0][0]}_to_{TRACK[1][0]}.tomtom")
def compareMemeChipTracks(infiles, outfile):
    '''Use tomtom to look if simlar motifs have been foundin more than one
    track '''

    track1, track2 = infiles

    PipelineMotifs.tomtom_comparison(track1, track2, outfile)


############################################################
@transform(getMemeChipSeedMotifs,
           regex("(.+).seeds"),
           r"\1.tsv")
def getMemeChipSeedTables(infile, outfile):

    code_dir = os.path.dirname(__file__)
    statement = ''' python %(code_dir)s/meme2table.py
                       -I %(infile)s
                       -S %(outfile)s
                       -L %(outfile)s.log'''
    P.run(statement)


############################################################
@subdivide(getMemeChipSeedTables,
           regex("memechip.dir/(.+).tsv"),
           add_inputs(r"memechip.dir/\1.seeds"),
           r"memechip.dir/\1*.png")
def getSeqlogosFromMemeChip(infiles, outfiles):
    '''get png seqlogos for each of the seed motifs found'''

    table, memefile = infiles
    track = re.match("memechip.dir/(.+).memechip.tsv",table).groups()[0]
    statement = []
    for nmotif, motif in enumerate(IOTools.open_file(table)):
        if motif.startswith("primary_id"):
            continue
        motif_id = motif.split("\t")[0]
        outfile = "memechip.dir/%(track)s_%(motif_id)s.png" % locals()
        statement.append('''ceqlogo -i%(nmotif)s %(memefile)s -h 7.5 -t "" -x "" -y ""
                                    | convert - %(outfile)s ''' % locals())

    statement = "; ".join(statement)
    P.run(statement)


############################################################
@follows(getSeqlogosFromMemeChip,
         compareMemeChipTracks)
def memechip():
    pass


############################################################
############################################################
@collate([compareMemeChipTracks,
          compareMemeTracks,
          compareDiscMemeTracks,
          compareRandDremeTracks,
          compareDiscDremeTracks],
         regex("(.+).dir/(.+).tomtom"),
         r"\1.dir/\1.track_comparisons.load")
def mergeAndLoadTrackComparisons(infiles, outfile):
    
    P.concatenate_and_load(infiles, outfile,
                         regex_filename=(".+/(.+)_to_(.+).tomtom"),
                         cat="track1,track2",
                         options="-i track1 -i track2")


############################################################
@transform([getMemeSeeds,
            getDiscMemeSeeds,
            getRandDremeSeeds,
            getDiscDremeSeeds],
           regex("(.+).seeds"),
           r"\1.tsv")
def getMemeDremeSeedTables(infile, outfile):

    code_dir = os.path.dirname(__file__)
    statement = '''python %(code_dir)s/meme2table.py
                       -I %(infile)s
                       -S %(outfile)s
                       -L %(outfile)s.log'''
    P.run(statement)


############################################################
@collate([getMemeChipSeedTables,
          getMemeDremeSeedTables],
         regex("(.+).dir/(.+).tsv"),
         r"\1.dir/\1_seeds.load")
def loadSeedTables(infiles, outfile):

    P.concatenate_and_load(infiles, outfile,
                         regex_filename=".+/(.+)\.(?:meme|dreme|memechip).tsv")


############################################################
@collate([runMeme,
          runDiscMEME,
          runRandDreme,
          runDiscDREME,
          runMemeChIP],
         regex("(?:.+_)?(.+).dir/(.+)"),
         r"\1_summary.load")
def loadMemeSummary(infiles, outfile):
    '''load information about motifs into database.'''

    outf = P.get_temp_file(".")

    outf.write("method\ttrack\n")

    for infile in infiles:
        if IOTools.is_empty(infile):
            continue
        method = re.match("(.+).dir/", infile).groups()[0]
        track = os.path.basename(".".join(infile.split(".")[:-1]))
        outf.write("%s\t%s\n" % (method,track))

    outf.close()

    P.load(outf.name, outfile)

    os.unlink(outf.name)


############################################################
@transform([exportMemeIntervalSequences,
            exportDremeIntervalSequences,
            exportMemeCHiPIntervalSequences],
           suffix(".fasta"),
           ".motifseq_stats.load")
def loadMotifSequenceComposition(infile, outfile):
    '''compute sequence composition of sequences used for ab-initio search.'''

    load_statement = P.build_load_statement(
        P.to_table(outfile))

    statement = '''
    cgat fasta2table
        --section=na
        --log=%(outfile)s
    < %(infile)s
    | %(load_statement)s
    > %(outfile)s'''

    P.run(statement)


############################################################
@transform([getMemeSeeds,
            getDiscMemeSeeds,
            getRandDremeSeeds,
            getDiscDremeSeeds],
           regex("(.+).seeds"),
           r"\1.tomtom")
def runTomTom(infile, outfile):
    '''compare ab-initio motifs against a databse of known motifs'''
    PipelineMotifs.runTomTom(infile, outfile)


############################################################
@transform(runTomTom, suffix(".tomtom"), "_tomtom.load")
def loadTomTom(infile, outfile):
    '''load tomtom results'''
    PipelineMotifs.loadTomTom(infile, outfile)


############################################################
@follows(runTomTom)
@transform([getMemeSeeds,
            getDiscMemeSeeds,
            getRandDremeSeeds,
            getDiscDremeSeeds],
           formatter("(?P<SAMPLE>.+)\.seeds"),
           add_inputs("{SAMPLE[0]}.tomtom"),
           "{SAMPLE[0]}_w_motif_enrichment.tomtom")
def outputTomTomWithMotifEnrichment(infiles, outfile):
    '''Original motif enrichment is not outputted in runTomTom output.
    This adds this data'''
    seed_result = infiles[0]

    tomtom_result = infiles[1]

    # Get the list of seed motifs
    PipelineMotifs.add_motif_enrichment_to_tomtom(seed_result,
                                                  tomtom_result,
                                                  outfile)




############################################################
@follows(meme,
         discMeme,
         randDreme,
         discDreme,
         memechip,
         loadTomTom,
         loadMotifSequenceComposition,
         loadMemeSummary,
         loadSeedTables,
         mergeAndLoadTrackComparisons)
def Motifs():
    pass

############################################################
############################################################
############################################################


@merge("*.motif", "motif_info.load")
def loadMotifInformation(infiles, outfile):
    '''load information about motifs into database.'''

    outf = P.get_temp_file(".")

    outf.write("motif\n")
 
    for infile in infiles:
        if IOTools.is_empty(infile):
            continue
        motif = P.snip(infile, ".motif")
        outf.write("%s\n" % motif)

    outf.close()


    P.load(outf.name, outfile, "--allow-empty-file")

    os.unlink(outf.name)

############################################################
############################################################
############################################################


@files_re((exportMotifIntervalSequences, exportMotifControlSequences),
          "(\S+).control.fasta",
          [r"\1.control.fasta", r"\1.foreground.fasta",
           glob.glob("*.motif")],
          r"\1.mast.gz")
def runMast(infiles, outfile):
    '''run mast on all intervals and motifs.

    Collect all results for an E-value up to 10000 so that all
    sequences are output and MAST curves can be computed.

    10000 is a heuristic.

    '''
    PipelineMotifs.runMAST(infiles, outfile)

############################################################
############################################################
############################################################


@transform(runMast,
           suffix(".mast.gz"),
           "_mast.load")
def loadMast(infile, outfile):
    '''parse mast file and load into database.

    Parse several motif runs and add them to the same
    table.

    Add columns for the control data as well.
    '''
    PipelineMotifs.loadMAST(infile, outfile)

############################################################
############################################################
############################################################


@follows(loadMotifInformation,
         mkdir(os.path.join(PARAMS["exportdir"],
                            "motifs")))
@merge(loadMast, "motifs.export")
def exportMotifLocations(infiles, outfile):
    '''export motif locations. There will be a bed-file per motif.

    Overlapping motif matches in different tracks will be merged.
    '''

    dbh = connect()
    cc = dbh.cursor()

    motifs = [x[0]
              for x in cc.execute("SELECT motif FROM motif_info").fetchall()]

    for motif in motifs:

        tmpf = P.get_temp_file(".")

        for infile in infiles:
            table = P.to_table(infile)
            track = P.snip(table, "_mast")
            for x in cc.execute(
                    """SELECT contig, start, end, '%(track)s', evalue
                    FROM %(table)s WHERE motif = '%(motif)s' AND
                    start IS NOT NULL""" % locals()):
                tmpf.write("\t".join(map(str, x)) + "\n")
        tmpf.close()

        outfile = os.path.join(
            PARAMS["exportdir"], "motifs", "%s.bed.gz" % motif)
        tmpfname = tmpf.name

        statement = '''mergeBed -i %(tmpfname)s -nms | gzip > %(outfile)s'''
        P.run(statement)

        os.unlink(tmpf.name)


@follows(Motifs)
def full():
    '''run the full pipeline.'''
    pass


@follows(mkdir("report"))
def build_report():
    '''build report from scratch.'''

    E.info("starting documentation build process from scratch")
    P.run_report(clean=True)


@follows(mkdir("report"))
def update_report():
    '''update report.'''

    E.info("updating documentation")
    P.run_report(clean=False)


@follows(update_report)
def publish():
    '''publish files.'''
    # publish web pages

    P.publish_report()


@merge(None, 'reset.log')
def reset(infile, outfile):
    '''reset pipeline to initial start.

    This will remove all results from running this
    pipeline!!!
    '''
    statement = '''
    rm -rf export motifs report;
    rm -rf _cache _static _templates;
    rm -f csvdb *.load;
    rm -f *.meme* *.tomtom* *.fasta*;
    rm -rf *.dir;
    '''
    P.run(statement)

def main(argv=None):
    if argv is None:
        argv = sys.argv
    P.main(argv)

if __name__ == "__main__":
    sys.exit(P.main(sys.argv))
