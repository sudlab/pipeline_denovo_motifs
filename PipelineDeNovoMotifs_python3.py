'''
pipeline_vitaminD_motifs.py - 
======================================================

:Author: Andreas Heger
:Release: $Id$
:Date: |today|
:Tags: Python

Purpose
-------

.. todo::
   
   describe purpose of the script.

Usage
-----

Example::

   python pipeline_vitaminD_motifs.py --help

Type::

   python pipeline_vitaminD_motifs.py --help

for command line help.

Documentation
-------------

Requirements:

* meme >= 4.9.1
* bioprospector >= 2004 (optional)

Code
----

'''
import re
import os
import tempfile
import collections
import shutil
import glob
import xml.etree.ElementTree
import random
import math

import pandas
import logging as L
from CGATCore import Experiment as E
from CGATCore import Pipeline as P
from CGATCore import IOTools as IOTools
import CGAT.IndexedFasta as IndexedFasta
import CGAT.Masker as Masker
import CGAT.Glam2Scan as Glam2Scan
import CGAT.MAST as MAST
import CGAT.Bed as Bed
import CGAT.FastaIterator as FastaIterator
from CGAT.MEME import MemeMotifFile, MotifCluster
from CGAT.FastaIterator import FastaRecord
import pandas as pd

# Set from importing module
PARAMS = {}


def maskSequences(sequences, masker=None):
    '''return a list of masked sequence.

    *masker* can be one of
        dust/dustmasker * run dustmasker on sequences
        softmask        * use softmask to hardmask sequences
    '''


    if masker in ("dust", "dustmasker"):
        masker_object = Masker.MaskerDustMasker()
    else:
        masker_object = None

    if masker == "softmask":
        # the genome sequence is repeat soft-masked
        masked_seq = sequences
    elif masker in ("dust", "dustmasker"):
        # run dust
        masked_seq = masker_object.maskSequences(
            [x.upper() for x in sequences])
    elif masker is None:
        masked_seq = [x.upper() for x in sequences]
    else:
        raise ValueError("unknown masker %s" % masker)

    # hard mask softmasked characters
    masked_seq = [re.sub("[a-z]", "N", x) for x in masked_seq]

    return masked_seq


def exportSequencesFromBedFile(infile, outfile, masker=None, mode="intervals"):
    '''export sequences for intervals in :term:`bed`-formatted *infile* 
    to :term:`fasta` formatted *outfile*
    '''

    track = P.snip(infile, ".bed.gz")

    fasta = IndexedFasta.IndexedFasta(
        os.path.join(PARAMS["genome_dir"], PARAMS["genome"]))
    outs = IOTools.open_file(outfile, "w")

    ids, seqs = [], []
    for bed in Bed.setName(Bed.iterator(IOTools.open_file(infile))):
        lcontig = fasta.getLength(bed.contig)

        if mode == "intervals":
            seqs.append(fasta.getSequence(bed.contig, "+", bed.start, bed.end))
            ids.append("%s_%s %s:%i..%i" %
                       (track, bed.name, bed.contig, bed.start, bed.end))

        elif mode == "leftright":
            l = bed.end - bed.start

            start, end = max(0, bed.start - l), bed.end - l
            ids.append("%s_%s_l %s:%i..%i" %
                       (track, bed.name, bed.contig, start, end))
            seqs.append(fasta.getSequence(bed.contig, "+", start, end))

            start, end = bed.start + l, min(lcontig, bed.end + l)
            ids.append("%s_%s_r %s:%i..%i" %
                       (track, bed.name, bed.contig, start, end))
            seqs.append(fasta.getSequence(bed.contig, "+", start, end))

    masked = maskSequences(seqs, masker)
    outs.write("\n".join([">%s\n%s" % (x, y) for x, y in zip(ids, masked)]))

    outs.close()


def bedsFromList(data):
    ''' takes a list of data and returns a bed object'''

    for interval in data:
        bed = Bed.Bed()
        try:
            bed.contig, bed.start, bed.end = \
                        interval[0], int(interval[1]), int(interval[2])
        except IndexError:
            raise ValueError("Insufficient fields to generate bed entry")
        except ValueError:
            raise ValueError("Fields 2 and 3 must be integer")
        bed.fields = interval[3:]

        yield bed


def shiftBeds(beds):
    '''Create two bed intervals on either side of the
    original interval'''

    for bed in beds:
        
        left = bed.copy()
        right = bed.copy()
        length = bed.end - bed.start

        left.start = bed.start - length
        left.end = bed.start
        left["name"] = str(left["name"]) + "_left"
        yield left

        right.start = bed.end
        right.end = bed.end + length
        right["name"] = str(bed["name"]) + "_right"
        yield right


def centreAndCrop(beds, halfwidth):
    '''Centre each bed around the peakcenter (given in the thickStart entry)
    and take "halfwidth" either side'''

    for bed in beds:
        try:
            bed.start = int(bed["thickStart"]) - halfwidth
            bed.end = int(bed["thickStart"]) + halfwidth
        except:
            print(bed.fields)
            raise

        yield bed


def truncateList(data, track, proportion=None, min_sequences=None, num_sequences=None,
                 shuffle=False):

    if shuffle:
        random.shuffle(data)

    if proportion:
        cutoff = int(math.ceil(len(data) * proportion))
        if min_sequences:
            cutoff = max(cutoff, min_sequences)
    elif num_sequences:
        cutoff = num_sequences
    else:
        cutoff = len(data)

    E.debug("Using %s sequences for track %s" % (cutoff, track))

    for dataum in data[:cutoff]:
        yield dataum


def getFASTAFromBed(beds, fasta, stranded, offset, maxsize):

    # get the sequences - cut at number of nucleotides
    current_size, nseq = 0, 0
    for bed in beds:
        lcontig = fasta.getLength(bed.contig)
        start, end = max(0, bed.start + offset), min(bed.end + offset, lcontig)

        if start >= end:
            E.info("writeSequencesForIntervals %s: sequence %s is empty: "
                   "start=%i, end=%i, offset=%i - ignored" %
                   (bed.track, bed.name, start, end, offset))
            continue

        if stranded:
            strand = bed.strand
        else:
            strand = "+"
        
        
        seq = ''.join(fasta.getSequence(bed.contig, 
                                        strand, interval[0], interval[1])
                      for interval in bed.toIntervals())

        current_size += len(seq)
        if maxsize and current_size >= maxsize:
            E.info("writeSequencesForIntervals %s: maximum size (%i) reached"
                   "- only %i sequences output" % (bed.track, maxsize, nseq))
            break

        nseq += 1

        yield FastaRecord("%s_%s %s:%i-%i" %
                          (bed.track, str(bed.name), bed.contig,
                           bed.start, bed.end),
                          seq)

        

def shuffleFasta(sequences):

    for fasta in sequences:
        sequence = list(fasta.sequence)
        random.shuffle(sequence)
        fasta.sequence = "".join(sequence)
        yield fasta

def writeSequencesForIntervals(track,
                               filename,
                               dbhandle,
                               full=False,
                               halfwidth=None,
                               maxsize=None,
                               proportion=None,
                               masker=[],
                               offset=0,
                               shuffled=False,
                               num_sequences=None,
                               min_sequences=None,
                               order="peakval",
                               shift=None,
                               stranded=False):
    '''build a sequence set for motif discovery. Intervals are taken from
    the table <track>_intervals in the database *dbhandle* and save to
    *filename* in :term:`fasta` format.

    If num_shuffles is set, shuffled copies are created as well with
    the shuffled number appended to the filename.

    The sequences are masked before shuffling (is this appropriate?)

    If *full* is set, the whole intervals will be output, otherwise
    only the region around the peak given by *halfwidth*

    If *maxsize* is set, the output is truncated at *maxsize* characters
    in order to create jobs that take too long.

    If proportion is set, only the top *proportion* intervals are output
    (sorted by peakval).

    If *num_sequences* is set, the first *num_sequences* will be used.

    *masker* can be a combination of
        * dust, dustmasker: apply dustmasker
        * softmask: mask softmasked genomic regions

    *order* is the order by which peaks should be sorted. Possible
    values are 'peakval' (peak value, descending order), 'score' (peak
    score, descending order)

    If *shift* is set, intervals will be shifted. ``leftright``
    creates two intervals on the left and right of the actual
    interval. The intervals will be centered around the mid-point and
    truncated the same way as the main intervals.

    '''
    cc = dbhandle.cursor()

    orderby = ""
    if order == "peakval":
        orderby = " ORDER BY peakval DESC"
    elif order == "max":
        orderby = " ORDER BY score DESC"
    elif order != "random":
        raise ValueError(
            "Unknown value passed as order parameter, check your ini file")

    tablename = "%s_intervals" % P.tablequote(track)
    statement = '''SELECT contig, start, end, interval_id, score, strand, peakcenter 
                       FROM %(tablename)s 
                       ''' % locals() + orderby

    cc.execute(statement)
    data = cc.fetchall()
    cc.close()

    E.debug("Got %s intervals for track %s" % ( len(data), track))
    if len(data) == 0:
        P.touch(filename)
        return

    data = truncateList(data, track,
                        proportion, min_sequences, num_sequences,
                        order == "random")

    beds = bedsFromList(data)

    L.info("writeSequencesForIntervals %s: masker=%s" % (track, str(masker)))

    fasta = IndexedFasta.IndexedFasta(
        os.path.join(PARAMS["genome_dir"], PARAMS["genome"]))
   
    # At the moment the pipeline retrieves from the DB the bed regions and they will
    # always be in the positive strand but if this were to change. The regions retrieved from 
    # the negative strand will be counted from the end of the chromosome and not the beginning without this.
    # This should be tested.
    fasta.setConverter(IndexedFasta.getConverter("zero-single-open"))

    # modify the ranges
    if shift == "leftright":
        beds = shitfBeds(beds)

    if halfwidth and not full:
        beds = centreAndCrop(beds, halfwidth)

    sequences = getFASTAFromBed(beds, fasta, stranded, offset, maxsize)

    if shuffled:
        sequences = shuffleFasta(sequences)

    c = E.Counter()
    outs = IOTools.open_file(filename, "w")
    for masker in masker:
        if masker not in ("unmasked", "none", None):
            ids, sequences = zip(*[(x.title, x.sequence) for x in sequences])
            sequences = maskSequences(sequences, masker)
            sequences = (FastaRecord(id, seq) for id, seq in zip(ids, sequences))


    with IOTools.open_file(filename, "w") as outs:

        for sequence in sequences:
            c.input += 1
            if len(sequence.sequence) == 0:
                c.empty += 1
                continue
            if len(sequence.sequence) < 0:
                c.too_short += 1
                continue

            outs.write(">%s\n%s\n" % (sequence.title, sequence.sequence))
            c.output += 1
        outs.close()

    E.info("%s" % c)

    return c.output


############################################################
############################################################
############################################################
def loadMAST(infile, outfile):
    '''parse mast file and load into database.

    Parse several motif runs and add them to the same
    table.

    Add columns for the control data as well.
    '''

    tablename = P.to_table(outfile)

    tmpfile = P.get_temp_file(".")

    tmpfile.write(MAST.Match().header +
                  "\tmotif\tcontig"
                  "\tl_evalue\tl_pvalue\tl_nmatches\tl_length\tl_start\tl_end"
                  "\tr_evalue\tr_pvalue\tr_nmatches\tr_length\tr_start\tr_end"
                  "\tmin_evalue\tmin_pvalue\tmax_nmatches" + "\n")

    lines = IOTools.open_file(infile).readlines()
    chunks = [x for x in range(len(lines)) if lines[x].startswith("::")]
    chunks.append(len(lines))

    def readChunk(lines, chunk):
        # use real file, as MAST parser can not deal with a
        # list of lines
        tmpfile2 = P.get_temp_file(".")
        try:
            motif, part = re.match(
                ":: motif = (\S+) - (\S+) ::", lines[chunks[chunk]]).groups()
        except AttributeError:
            raise ValueError(
                "parsing error in line '%s'" % lines[chunks[chunk]])

        E.info("reading %s - %s" % (motif, part))

        tmpfile2.write("".join(lines[chunks[chunk] + 1:chunks[chunk + 1]]))
        tmpfile2.close()

        mast = MAST.parse(IOTools.open_file(tmpfile2.name, "r"))

        os.unlink(tmpfile2.name)

        return motif, part, mast

    def splitId(s, mode):
        '''split background match id

        has three parts: track _ id _ pos

        track might contain '_'.
        '''
        d = match.id.split("_")
        if mode == "bg":
            return "_".join(d[:-2]), d[-2], d[-1]
        elif mode == "fg":
            return "_".join(d[:-1]), d[-1]

    for chunk in range(0, len(chunks) - 1, 2):

        motif_fg, part, mast_fg = readChunk(lines, chunk)
        assert part == "foreground"
        motif_bg, part, mast_bg = readChunk(lines, chunk + 1)
        assert part == "background"
        assert motif_fg == motif_bg

        # index control data
        controls = collections.defaultdict(dict)
        for match in mast_bg.matches:
            track, id, pos = splitId(match.id, "bg")
            controls[id][pos] = (
                match.evalue, match.pvalue, match.nmotifs, match.length, match.start, match.end)

        for match in mast_fg.matches:
            # remove track and pos
            track, match.id = splitId(match.id, "fg")
            # move to genomic coordinates
            contig, start, end = re.match(
                "(\S+):(\d+)..(\d+)", match.description).groups()
            if match.nmotifs > 0:
                start, end = int(start), int(end)
                match.start += start
                match.end += start
                match.positions = [x + start for x in match.positions]

            id = match.id
            if id not in controls:
                P.warn("no controls for %s - increase MAST evalue" % id)

            if "l" not in controls[id]:
                controls[id]["l"] = (
                    float(PARAMS["mast_evalue"]), 1, 0, 0, 0, 0)
            if "r" not in controls[id]:
                controls[id]["r"] = (
                    float(PARAMS["mast_evalue"]), 1, 0, 0, 0, 0)

            min_evalue = min(controls[id]["l"][0], controls[id]["r"][0])
            min_pvalue = min(controls[id]["l"][1], controls[id]["r"][1])
            max_nmatches = max(controls[id]["l"][2], controls[id]["r"][2])

            tmpfile.write(str(match) + "\t%s\t%s\t%s\t%s\t%s\t%s\t%s" %
                          (motif_fg, contig,
                           "\t".join(map(str, controls[id]["l"])),
                           "\t".join(map(str, controls[id]["r"])),
                           str(min_evalue),
                           str(min_pvalue),
                           str(max_nmatches),
                           ) + "\n")

    tmpfile.close()

    P.load(tmpfile.name,
           outfile,
           options="--add-index=id "
           "--add-index=motif "
           "--add-index=id,motif "
           "--allow-empty-file "
           "--map=base_qualities:text")

    os.unlink(tmpfile.name)



def runMAST(infiles, outfile):
    '''run mast on all intervals and motifs.

    Collect all results for an E-value up to 10000 so that all
    sequences are output and MAST curves can be computed.

    10000 is a heuristic.

    '''

    # job_options = "-l mem_free=8000M"

    controlfile, dbfile, motiffiles = infiles

    if IOTools.is_empty(dbfile):
        P.touch(outfile)
        return

    if not os.path.exists(controlfile):
        raise ValueError(
            "control file %s for %s does not exist" % (controlfile, dbfile))

    # remove previous results
    if os.path.exists(outfile):
        os.remove(outfile)

    tmpdir = P.get_temp_dir(".")
    tmpfile = P.get_temp_filename(".")

    for motiffile in motiffiles:
        if IOTools.is_empty(motiffile):
            L.info("skipping empty motif file %s" % motiffile)
            continue

        of = IOTools.open_file(tmpfile, "a")
        motif, x = os.path.splitext(motiffile)
        of.write(":: motif = %s - foreground ::\n" % motif)
        of.close()

        # mast bails if the number of nucleotides gets larger than
        # 2186800982?
        # To avoid this, run db and control file separately.
        statement = '''
        cat %(dbfile)s
        | mast %(motiffile)s - -nohtml -oc %(tmpdir)s -ev %(mast_evalue)f %(mast_options)s >> %(outfile)s.log 2>&1;
        cat %(tmpdir)s/mast.txt >> %(tmpfile)s 2>&1
        '''
        P.run(statement)

        of = IOTools.open_file(tmpfile, "a")
        motif, x = os.path.splitext(motiffile)
        of.write(":: motif = %s - background ::\n" % motif)
        of.close()

        statement = '''
        cat %(controlfile)s
        | mast %(motiffile)s - -nohtml -oc %(tmpdir)s -ev %(mast_evalue)f %(mast_options)s >> %(outfile)s.log 2>&1;
        cat %(tmpdir)s/mast.txt >> %(tmpfile)s 2>&1
        '''
        P.run(statement)

    statement = "gzip < %(tmpfile)s > %(outfile)s"
    P.run(statement)

    shutil.rmtree(tmpdir)
    os.unlink(tmpfile)


############################################################
############################################################
############################################################
def runGLAM2(infile, outfile, dbhandle):
    '''run glam2 on all intervals and motifs.

    In order to increase the signal/noise ratio,
    MEME is not run on all intervals but only the 
    top 10% of intervals (peakval) are used. 
    Also, only the segment of 200 bp around the peak
    is used and not the complete interval.

    * Softmasked sequence is converted to hardmasked
      sequence to avoid the detection of spurious motifs.

    * Sequence is run through dustmasker
    '''
    to_cluster = True

    target_path = os.path.join(
        os.path.abspath(PARAMS["exportdir"]), "glam2", outfile)
    track = infile[:-len(".fasta")]

    tmpdir = tempfile.mkdtemp()
    tmpfasta = os.path.join(tmpdir, "in.fa")

    nseq = PipelineMotifs.writeSequencesForIntervals(
        track, tmpfasta,
        dbhandle,
        full=False,
        halfwidth=int(
            PARAMS["meme_halfwidth"]),
        maxsize=int(
            PARAMS["meme_max_size"]),
        proportion=PARAMS["meme_proportion"])
    
    min_sequences = int(nseq / 10.0)
    statement = '''
    %(execglam2)s -2 -O %(tmpdir)s %(glam2_options)s -z %(min_sequences)i n %(tmpfasta)s > %(outfile)s.log
    '''
    P.run(statement)

    # copy over results
    try:
        os.makedirs(os.path.dirname(target_path))
    except OSError:
        # ignore "file exists" exception
        pass

    if os.path.exists(target_path):
        shutil.rmtree(target_path)
    shutil.move(tmpdir, target_path)

    shutil.copyfile(os.path.join(target_path, "glam2.txt"), outfile)


def collectMEMEResults(tmpdir, target_path, outfile,
                       method="meme"):
    '''collect output from a MEME run in tmpdir
    and copy all over to target_path

    convert images output by MEME (.eps files) to 
    .png files.'''

    # copy over results
    try:
        os.makedirs(os.path.dirname(target_path))
    except OSError:
        # ignore "file exists" exception
        pass

    if os.path.exists(target_path):
        shutil.rmtree(target_path)
    E.debug("tmpdir is" + tmpdir + "      target_path is " + target_path)
    shutil.move(tmpdir, target_path)

    if method == "dreme":
        shutil.copyfile(os.path.join(target_path, "dreme.txt"), outfile)
    elif method == "meme":
        shutil.copyfile(os.path.join(target_path, "meme.txt"), outfile)
    elif method == "memechip":
        try:
            shutil.copyfile(os.path.join(target_path, "combined.meme"), outfile)
        except IOError:
            E.warn ("%s: No motifs found")
            P.touch(outfile)

    # convert images to png
    epsfiles = glob.glob(os.path.join(target_path, "*.eps"))

    statement = []
    for epsfile in epsfiles:
        b, ext = os.path.splitext(epsfile)
        pngfile = b + ".png"
        statement.append('''convert %(epsfile)s %(pngfile)s ''')
    
    if len(statement) > 0:
        statement = " && ".join(statement)
        P.run(statement)



def runMEME(track, outfile, dbhandle):
    '''run MEME to find motifs.

    In order to increase the signal/noise ratio,
    MEME is not run on all intervals but only the
    top 10% of intervals (peakval) are used.
    Also, only the segment of 200 bp around the peak
    is used and not the complete interval.

    * Softmasked sequence is converted to hardmasked
      sequence to avoid the detection of spurious motifs.

    * Sequence is run through dustmasker

    This method is deprecated - use runMEMEOnSequences instead.
    '''
    # job_options = "-l mem_free=8000M"

    target_path = os.path.join(
        os.path.abspath(PARAMS["exportdir"]), outfile)

    fasta = IndexedFasta.IndexedFasta(
        os.path.join(PARAMS["genome_dir"], PARAMS["genome"]))

    tmpdir = P.get_temp_dir(".")
    tmpfasta = os.path.join(tmpdir, "in.fa")

    nseq = writeSequencesForIntervals(
        track, tmpfasta,
        dbhandle,
        full=False,
        masker=P.as_list(PARAMS['motifs_masker']),
        halfwidth=int(PARAMS["meme_halfwidth"]),
        maxsize=int(PARAMS["meme_max_size"]),
        proportion=PARAMS["meme_proportion"],
        min_sequences=PARAMS["meme_min_sequences"])

    if nseq == 0:
        E.warn("%s: no sequences - meme skipped" % outfile)
        P.touch(outfile)
    else:
        statement = '''
        meme %(tmpfasta)s -dna -revcomp -mod %(meme_model)s -nmotifs %(meme_nmotifs)s -oc %(tmpdir)s -maxsize %(meme_max_size)s %(meme_options)s > %(outfile)s.log
        '''
        P.run(statement)

        collectMEMEResults(tmpdir, target_path, outfile)


def generatePSP(positives, negatives, outfile):
    ''' generate a discrimitative PSP file from
    the positives and negatives that can be used
    to do descriminative MEME '''

    psp_options = PARAMS["psp_options"]
    
    nseqs_pos = int(FastaIterator.count(positives))
    nseqs_neg = int(FastaIterator.count(negatives))

    if nseqs_pos < 2 or nseqs_neg < 2:
        E.warn("%s: input files do not have sufficent sequences"
               "to run psp-gen, skipping" % outfile)
        P.touch(outfile)
        return

    # get appropriate options from meme options
    if PARAMS.get("meme_revcomp", True):
        psp_options += " -revcomp"

    statement = '''psp-gen -pos %(positives)s
                           -neg %(negatives)s
                           %(psp_options)s
                   > %(outfile)s '''

    P.run(statement)

        
def runMEMEOnSequences(infile, outfile, background=None,
                       psp=None):
    '''run MEME on fasta sequences to find motifs
   
    By defualt MEME calculates a zero-th order background
    model from the nucleotide frequencies in the input set.

    To use a different background set, a background
    file created by fasta-get-markov must be supplied.

    To perform descrimantive analysis a position specific
    prior (psp) file must be provided. This can be generated
    used generatePSP.

    '''
    # job_options = "-l mem_free=8000M"

    nseqs = int(FastaIterator.count(infile))
    if nseqs < 2:
        E.warn("%s: less than 2 sequences - meme skipped" % outfile)
        P.touch(outfile)
        return

    # Get the total length of the sequences to decide the memory
    total_seqs_length = 0

    with IOTools.open_file(infile, "r") as fasta_reader:

        iterator_fasta = FastaIterator.iterate(fasta_reader)

        for fasta_seq in iterator_fasta:
            total_seqs_length += len(fasta_seq.sequence)

    fasta_reader.close()

    # If the length of all sequences is higher than 160,000bp
    # Up the memory
    job_memory = "2G"

    if (total_seqs_length > 160000):
        job_memory = "4G"
    
    if PARAMS.get("meme_revcomp", True):
        revcomp = "-revcomp"
    else:
        revcomp = ""

    target_path = os.path.join(
        os.path.abspath(PARAMS["exportdir"]),  outfile)
    tmpdir = P.get_temp_dir(".")
    if background:
        background_model = "-bfile %s" % background
    else:
        background_model = ""

    if psp:
        E.info("Running MEME in descriminative mode")
        psp_file = "-psp %s" % psp
    else:
        psp_file = ""
    
    
    statement = '''
    meme %(infile)s -dna %(revcomp)s
    -p %(meme_threads)s
    -mod %(meme_model)s
    -nmotifs %(meme_nmotifs)s
    -oc %(tmpdir)s
    -maxsize %(meme_max_size)s
    %(background_model)s
    %(psp_file)s
    %(meme_options)s
       2> %(outfile)s.log
    '''

    # If running with more than one thread
    # http://git.net/ml/clustering.gridengine.users/2007-04/msg00058.html
    # specify "excl=false -w n -pe openmpi-ib num_threads" in cluster_options
    # through job_options
    if int(PARAMS["meme_threads"]) != 1:
        job_options = str(PARAMS["meme_job_options"])
        job_threads = int(PARAMS["meme_threads"])
        cluster_parallel_environment = str(PARAMS["meme_cluster_parallel_environment"])
    
     
    P.run(statement)

    collectMEMEResults(tmpdir, target_path, outfile)


def runMemeCHIP(infile, outfile, motifs=None):
    '''Run the MEME-CHiP pipeline on the input files.
    optional motifs files can be supplied as a list'''

    if motifs:
        motifs = " ".join("-db %s" % motif for motif in motifs)
    else:
        motifs = " "

    nseqs = int(FastaIterator.count(infile))
    if nseqs == 0:
        E.warn("%s: no sequences - meme-chip skipped")
        P.touch(outfile)
        return

    target_path = os.path.join(
        os.path.abspath(PARAMS["exportdir"]), outfile)
    tmpdir = P.get_temp_dir(".")

    statement = '''
    meme-chip %(infile)s
             -p %(meme_threads)s 
             -oc %(tmpdir)s
             -nmeme %(memechip_nmeme)s
             %(memechip_options)s     
             %(motifs)s > %(outfile)s.log '''
    
    # If running with more than one thread
    # http://git.net/ml/clustering.gridengine.users/2007-04/msg00058.html
    # specify "excl=false -w n -pe openmpi-ib num_threads" in cluster_options
    # through job_options
    if int(PARAMS["memechip_threads"]) != 1:
        job_options = str(PARAMS["memechip_job_options"])
        job_threads = int(PARAMS["memechip_threads"])
        cluster_parallel_environment = str(PARAMS["memechip_cluster_parallel_environment"])
     
    
    P.run(statement)
   

    collectMEMEResults(tmpdir, target_path, outfile, method="memechip")

    
def runTomTom(infile, outfile):
    '''compare ab-initio motifs against tomtom.'''

    tmpdir = P.get_temp_dir(".")
    databases = " ".join(P.as_list(PARAMS["tomtom_databases"]))


    target_path = os.path.join(
        os.path.abspath(PARAMS["exportdir"]), "tomtom", outfile)

    if IOTools.is_empty(infile):
        E.warn("input is empty - no computation performed")
        P.touch(outfile)
        return

    statement = '''
    tomtom %(tomtom_options)s -oc %(tmpdir)s %(infile)s %(databases)s > %(outfile)s.log
    '''

    P.run(statement)

    # copy over results
    try:
        os.makedirs(os.path.dirname(target_path))
    except OSError:
        # ignore "file exists" exception
        pass

    if os.path.exists(target_path):
        shutil.rmtree(target_path)
    shutil.move(tmpdir, target_path)

    shutil.copyfile(os.path.join(target_path, "tomtom.txt"), outfile)
def loadTomTom(infile, outfile):
    '''load tomtom results'''

    tablename = P.to_table(outfile)

    resultsdir = os.path.join(
        os.path.abspath(PARAMS["exportdir"]), "tomtom", infile)
    xml_file = os.path.join(resultsdir, "tomtom.xml")

    if not os.path.exists(xml_file):
        E.warn("no tomtom output - skipped loading ")
        P.touch(outfile)
        return

    # get the motif name from the xml file

    tree = xml.etree.ElementTree.ElementTree()
    tree.parse(xml_file)
    motifs = tree.find("targets")
    name2alt = {}
    for motif in motifs.getiterator("motif"):
        name = motif.get("id")
        alt = motif.get("alt")
        name2alt[name] = alt

    tmpfile = P.get_temp_file(".")

    # parse the text file
    for line in IOTools.open_file(infile):
        if line.startswith("#Query"):
            tmpfile.write('\t'.join(
                ("target_name", "query_id", "target_id",
                 "optimal_offset", "pvalue", "evalue",
                 "qvalue", "Overlap", "query_consensus",
                 "target_consensus", "orientation")) + "\n")
            continue
        data = line[:-1].split("\t")
        target_name = name2alt[data[1]]
        tmpfile.write("%s\t%s" % (target_name, line))
    tmpfile.close()

    P.load(tmpfile.name, outfile)

    os.unlink(tmpfile.name)


def runDREME(infile, outfile, neg_file = "", options = ""):
    ''' Run DREME on fasta file. If a neg_file is passed
    then DREME will use this as the negative set, otherwise
    the default is to shuffle the input '''

    nseqs_pos = int(FastaIterator.count(infile))
    if nseqs_pos < 2:
        E.warn("%s: less than 2 sequences - dreme skipped" % outfile)
        P.touch(outfile)
        return
    
    if neg_file:
        nseqs_neg = int(FastaIterator.count(neg_file))
        if nseqs_neg < 2:
            E.warn("%s: less than 2 sequences in negatives file - dreme skipped"
                   % outfile)
            P.touch(outfile)
            return
        else:
            neg_file = "-n %s" % neg_file

    logfile = outfile + ".log"
    target_path = os.path.join(
        os.path.abspath(PARAMS["exportdir"]), outfile)
    tmpdir = P.get_temp_dir(".")

    statement = '''
    dreme -p %(infile)s %(neg_file)s -png
        -oc %(tmpdir)s
            %(dreme_options)s
            %(options)s
       > %(logfile)s
    '''

    P.run(statement)

    collectMEMEResults(tmpdir, target_path, outfile, method="dreme")

def runFIMO(motifs, database, outfile, exportdir, options={}):
    '''run fimo to look for occurances of motifs supplied in sequence database.
    :param:`motifs` is the path to a MEME formated motif file.
    :param:`database` is a fasta file.
    :param:`outfile` is the text output from fimo
    :param:`exportdir` specifies the directory to put exported files (html,gff)
    :param:options is a dictionary: {'option':'value'} will be passed as
                    --option=value and will overwrite options specified in the
                     PARAMs'''


    # if the motifs file is empty, then fimo will return an error
    # this isn't very useful behavoir.

    inlines = IOTools.open_file(motifs).read()
    #print inlines
    if not re.search("MOTIF", inlines):
        E.warning("No motifs found in %s" % motifs)
        P.touch(outfile)
        return
    else:
        E.debug("%s: %i motifs found" % 
                (motifs, len(re.findall("MOTIF", inlines))))


    fimo_options = PARAMS.get("fimo_options", "")
    for option, value in options.iteritems():
        fimo_options = re.sub("%s=\S+" % option, "", fimo_options)
        if value is None:
            fimo_options += " --%s" % option
        else:
            fimo_options += " --%s=%s" % (option, value)

    tmpout = P.get_temp_filename()
    
    track = os.path.basename(outfile)
    exportdir = os.path.abspath(exportdir)

    xmlout = P.snip(outfile,".txt") + ".xml"
    logfile = P.snip(outfile,".txt") + ".log"
    gffout = os.path.join(exportdir, track + ".gff")
    htmlout = os.path.join(exportdir, track + ".html")
    
    statement = ''' fimo --oc %(tmpout)s
                         %(fimo_options)s
                         %(motifs)s
                         %(database)s &> %(logfile)s;
                     mv %(tmpout)s/fimo.txt %(outfile)s;
                     mv %(tmpout)s/fimo.xml %(xmlout)s;
                     mv %(tmpout)s/fimo.gff %(gffout)s;
                     mv %(tmpout)s/fimo.html %(htmlout)s;
                     rm -r %(tmpout)s '''

    P.run(statement)


def getSeedMotifs(motif_file, tomtom_file, outfile):

    ungrouped = MemeMotifFile(IOTools.open_file(motif_file))

    if len(ungrouped) == 0:
        with IOTools.open_file(outfile, "w") as outf:
            outf.write(str(ungrouped))
        return

    E.debug("%s: Loaded %i motifs" % (motif_file, len(ungrouped)))
    tomtom = pandas.read_csv(tomtom_file, sep="\t")
    tomtom["Query ID"] = tomtom["Query ID"].astype(str)
    tomtom["Target ID"] = tomtom["Target ID"].astype(str)
    
    all_clusters = ungrouped.keys()
    new_index = pandas.MultiIndex.from_product([all_clusters, all_clusters],
                                               names=["#Query ID", "Target ID"])
    tomtom = tomtom.set_index(["Query ID", "Target ID"])
    tomtom = tomtom.reindex(new_index)
    tomtom = tomtom.sort_index()
    ungrouped.sort("evalue")

    groups = []

    E.debug("%s: Clustering Motifs" % motif_file)
    while len(ungrouped) > 0:
        cur_cluster = MotifCluster(ungrouped.take(0))
        assert len(cur_cluster) == 1
        E.debug("%s: working on cluster %s, %i clusters remaining" %
                (motif_file, cur_cluster.seed.primary_id, len(ungrouped)))

        try:
            seed_distances = tomtom.loc[cur_cluster.seed.primary_id]
        except:
            print(tomtom)
            raise

        close_motifs = seed_distances[seed_distances["q-value"] < 0.05]

        E.debug("%s: Found %i similar motifs" %
                (motif_file, close_motifs.shape[0]))

        for motif in close_motifs.index.values:
            if motif == cur_cluster.seed.primary_id:
                continue
            
            try:
                cur_cluster.append(ungrouped.take(str(motif)))
            except KeyError:
                pass

        groups.append(cur_cluster)

    E.debug("%s: Got %i clusters, containing a total of %i motifs"
            % (motif_file, len(groups), sum(len(cluster) for cluster in groups)))

    groups.sort(key=lambda cluster: cluster.seed.evalue)

    merged_groups = []

    E.debug("%s: Merging groups with weak similarity" % motif_file)
    while len(groups) > 0:

        cur_cluster = groups.pop(0)
        
        to_merge = []
        
        distances = tomtom.loc[cur_cluster.seed.primary_id]
        for other_cluster in groups:

            qvals = distances.loc[other_cluster.keys()]["q-value"]

            if (qvals < 0.1).all():
                to_merge.append(other_cluster)

        for cluster in to_merge:
            cur_cluster.extend(cluster)
            groups.remove(cluster)

        merged_groups.append(cur_cluster)

    E.debug("%i final clusters found" % len(merged_groups))
    E.debug("%s bulding output" % motif_file)
    for group in merged_groups:
        group.seed.letter_probability_line += " nClustered= %i" % len(group)
        group.seed.letter_probability_line += " totalHits= %i" % sum(
            [motif.nsites for motif in group])

    output = MemeMotifFile(ungrouped)
    output.extend(cluster.seed for cluster in merged_groups)
    
    E.debug("%s outputting" % motif_file)
    with IOTools.open_file(outfile, "w") as outf:
        outf.write(str(output))


def tomtom_comparison(track1, track2, outfile):
    ''' Compare two meme files to each other using
    tomtom '''

    statement = ''' tomtom -verbosity 1 -text -thresh 0.05
                     %(track1)s
                     %(track2)s  2> %(outfile)s.log 
                 | sed 's/#//' > %(outfile)s '''
    P.run(statement)


def add_motif_enrichment_to_tomtom(seed_result, tomtom_result, outfile):

    seed_results = MemeMotifFile(IOTools.open_file(seed_result))

    tomtom_result_df = pd.read_table(tomtom_result,
                                       header=0,
                                       delimiter="\t")

    # Create a dataframe with the seeds and their corresponding motif enrichment E-value
    seed_results_df = pandas.DataFrame(columns=['#Query ID', 'Enrichment E-value'])

    for seed in seed_results.motifs:

        # Create a series of key and enrichment E-value
        seed_series = pd.Series([seed.primary_id, seed.evalue], index=['#Query ID', 'Enrichment E-value'])

        seed_results_df = seed_results_df.append(seed_series, ignore_index=True)

    # Merge the tables only when common keys are found
    merged_df = tomtom_result_df.merge(seed_results_df, on="#Query ID", how="inner")


    # Output the table
    merged_df.to_csv(outfile,
                     sep="\t",
                     header=True,
                     mode="w",
                     index_label=None,
                     index=False, na_rep="NA",
                     line_terminator="\n")


