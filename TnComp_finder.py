#!/usr/bin/env python3
#
# TnComp_finder: composite transposon finder
#
# Version 1.0.0 - November 11, 2019
#
# Copyright © 2019 Danillo Oliveira Alvarenga
#
# TnComp_finder is free software: you can redistribute it and/or modify it under
# the terms of the GNU Affero General Public License as published by the
# Free Software Foundation, either version 3 of the License, or (at your
# option) any later version.
#
# TnComp_finder is distributed in the hope that it will be useful, but WITHOUT
# ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
# FITNESS FOR A PARTICULAR PURPOSE. See the GNU Affero General Public
# License for more details.
#
# You should have received a copy of the GNU Affero GPL along with TnComp_finder.
# If not, see <http://www.gnu.org/licenses/agpl-3.0.html>.
#
import argparse
import os
import subprocess
import sys
from functools import partial
from multiprocessing import Pool
from time import strftime
from Bio import SeqFeature
from Bio import SeqIO
from Bio.Alphabet import IUPAC
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord


# Parse arguments from the command line.
parser = argparse.ArgumentParser(description="composite transposon finder",
                                 formatter_class=lambda prog:
                                 argparse.HelpFormatter(prog,
                                 max_help_position=10, width=100))

parser.add_argument("-v", "--version", action="version",
                    version="%(prog)s 1.0.0", help="show version and exit")
parser.add_argument("-f", "--files", metavar="sequences.fasta", nargs='+',
                    help="target sequences", required=True)
parser.add_argument("-o", "--out", metavar="directory",
                    help="output directory")

parser.add_argument("-p", "--processors", metavar="threads", type=int,
                    help="processor threads available for analyses")

parser.add_argument("-g", "--gbk", action="store_true",
                    help="write a genbank file with predictions")

parser.add_argument("-i", "--identity", metavar="%", type=float,
                    help="minimum identity with the database")
parser.add_argument("-c", "--coverage", metavar="%", type=float,
                    help="minimum coverage with the database")
parser.add_argument("-d", "--distance", metavar="bp", type=int,
                    help="maximum distance between transposons")

parser.add_argument("-e", "--extend", metavar="bp", type=int,
                    help="retrieve sequence with extended borders")
parser.add_argument("-k", "--flanking", action="store_true",
                    help="check for identical candidates with distinct flanking regions")

m_excl = parser.add_mutually_exclusive_group()

m_excl.add_argument("-s", "--same", action="store_true",
                    help="consider only copies on the same strand")
m_excl.add_argument("-t", "--different", action="store_true",
                    help="consider only copies on different strands")

parser.set_defaults(out=os.getcwd(), processors=1,
                    identity=90, coverage=95, distance=20000)

args = parser.parse_args()


# Align given file to the database using BLASTn.
def run_blastn(filename, outdir):

    # Stop analysis if previous results are found.
    outfile = '.'.join(os.path.basename(filename).split('.')[:-1]) + ".blastn"
    if outfile in os.listdir(outdir + "/blastn"):
        print(strftime("%c" + "\tPrevious analysis found for " + outfile + ". Passing."))
        pass

    # Run BLASTn if no .blastn file is found in the output directory.
    else:
        db = os.path.dirname(os.path.realpath(__file__)) + "/db/transposons.fna"
        print(strftime("%c") + "\tBlasting file " + filename + '.')
        subprocess.call(["blastn",
                        "-query", filename,
                        "-subject", db,
                        "-out", outdir + "/blastn/" + outfile,
                        "-outfmt", "6 qseqid qstart qend sseqid slen sstart send pident length",
                        "-dust", "no",
                        "-word_size", "7",
                        "-evalue", "1e-5"])


# Verify BLASTn results and check for candidates.
def parse_blastn(filename, outdir):

    infile = os.path.basename(filename)
    outfile = '.'.join(infile.split('.')[:-1]) + ".blastn"
    candidates = {}
    n = 1

    print(strftime("%c" + "\tParsing " + outfile + " results file."))

    with open(outdir + "/blastn/" + outfile, "rt") as out:
        header, transposon = '', ''
        start, end, strand = 0, 0, ''
        length, ident = 0, 0
        matches = []

        # Reverse sort BLAST results according to identity.
        for line in out:
            line = line.strip().split("\t")
            matches.append(line)
        matches = sorted(matches, key=lambda x: (x[0], int(x[1]), -float(x[7])))

        for item in matches:
            qseqid, qstart, qend = item[0], int(item[1]), int(item[2])
            sseqid, sstart, send = item[3], int(item[5]), int(item[6])

            stransp = sseqid.split('|')[1]
            qlen = abs(qstart - qend) + 1
            pident = float(item[7])
            pcvg = int(item[8])*100 / int(item[4])

            if sstart < send:
                qstrand = '+'
            else:
                qstrand = '-'

            # Ignore according to strands.
            if args.same and qstrand != strand:
                continue
            elif args.different and qstrand == strand:
                continue

            # Ignore if BLASTn coverage or identity is low.
            elif pcvg < args.coverage or pident < args.identity:
                continue

            # Ignore if current sequence or transposon is different
            # from the previous ones and get their information.
            elif qseqid != header or stransp != transposon:
                header, transposon = qseqid, stransp
                start, end, strand = qstart, qend, qstrand
                length, ident = qlen, pident
                avg_ident = (ident + pident) / 2
                continue

            # Ignore if the distance betwee the current and the previous
            # transposons is too long.
            elif (abs(max(start, end) - min(qstart, qend)) - 1) > args.distance:
                continue

            else:
                info = (header, transposon, start, end, strand, qstart, qend, qstrand, avg_ident)
                candidates.update({filename + '_' + str(n): info})
                header, transposon = '', ''
                n += 1

    candidates = merge_candidates(candidates)
    if candidates:
        write_report(candidates, outdir)


# Verify features in common between candidates and remove the ones presenting
# lower identities.
def merge_candidates(candidates):

    # Do nothing if there is only one or no candidates.
    if len(candidates) <= 1:
        return candidates

    merged_candidates = {}
    keys = sorted(candidates.items(), key=lambda d: d[1][8], reverse=True)

    for k in keys:
        if not merged_candidates:
            merged_candidates.update({k[0]: candidates[k[0]]})

        else:
            A_start = int(candidates[k[0]][2])
            A_end = int(candidates[k[0]][3])
            B_start = int(candidates[k[0]][5])
            B_end = int(candidates[k[0]][6])
            avg_ident = candidates[k[0]][8]
            for i in merged_candidates.items():
                merged_A_start = int(i[1][2])
                merged_A_end = int(i[1][3])
                merged_B_start = int(i[1][5])
                merged_B_end = int(i[1][6])
                merged_avg_ident = i[1][8]

                # Do not add if it starts and ends with coordinates close to
                # those of a previous candidate.
                overlapping = ((A_start-10 <= merged_A_start <= A_start+10) and \
                               (A_end-10 <= merged_A_end <= A_end+10) and \
                               (B_end-10 <= merged_B_end <=B_end+10) and \
                               (B_end-10 <= merged_B_end <=B_end+10) and \
                               (merged_avg_ident > avg_ident))
                if overlapping:
                    break

            else:
                merged_candidates.update({k[0]: candidates[k[0]]})

    return merged_candidates


# Retrieve candidate sequences.
def get_sequence(filename, header, start, end, strand):

    with open(filename, "rt") as fasta:
        seq = False
        complete_seq = ''
        for line in fasta:
            if seq:
                if '>' in line:
                    break
                else:
                    complete_seq += line.strip()
            elif header in line:
                seq = True

    sequence = complete_seq[start-1:end]
    if strand and strand == '-':
        sequence = reverse_complement(sequence)

    return sequence


# Produce the reverse complement of a DNA sequence.
def reverse_complement(sequence):

    sequence = sequence[::-1]
    reverse_complemented_sequence = ''

    for character in sequence:
        if character is 'T':
            character = 'A'
        elif character is 'A':
            character = 'T'
        elif character is 'C':
            character = 'G'
        elif character is 'G':
            character = 'C'
        elif character is 'N':
            character = 'N'
        reverse_complemented_sequence += character

    return reverse_complemented_sequence


# Use Prodigal for predicting open reading frames in the detected sequence.
def predict_orfs(sequence, start):

    sequence_orfs = []
    sequence = ">sequence\n" + sequence
    prodigal = subprocess.check_output(
                            ["prodigal", "-f", "sco", "-p", "meta"],
                             input=sequence, universal_newlines=True,
                             stderr=subprocess.DEVNULL)
    prodigal = [x.lstrip(">") for x in prodigal.split("\n") if ">" in x]

    for item in prodigal:
        item = item.split("_")
        if item[3] is '+':
            beginning = int(item[1])
            end = int(item[2])
        elif item[3] is '-':
            beginning = int(item[2])
            end = int(item[1])
        sequence_orfs.append((beginning, end))

    complete_orfs = [(x+start, y+start) for (x, y) in sequence_orfs]

    return sequence_orfs, complete_orfs


# Get candidate information and generate a report.
def write_report(candidates, outdir):

    basename = '.'.join(list(candidates.keys())[0].split('/')[-1].split('.')[:-1])
    seq_info = {}
    print(strftime("%c") + "\tWriting " + basename + " report.")

    outfile = outdir + '/' + basename + "_composite.txt"
    if os.path.exists(outfile):
        os.remove(outfile)

    n = 1
    for key in sorted(candidates.keys()):
        filename = '_'.join(key.split('_')[:-1])
        header = candidates[key][0]
        transposon = candidates[key][1]

        A_start = candidates[key][2]
        A_end = candidates[key][3]
        A_len = abs(A_end - A_start) + 1
        A_strand = candidates[key][4]
        A_seq = get_sequence(filename, header, A_start, A_end, A_strand)

        B_start = candidates[key][5]
        B_end = candidates[key][6]
        B_len = abs(B_end - B_start) + 1
        B_strand = candidates[key][7]
        B_seq = get_sequence(filename, header, B_start, B_end, B_strand)

        min_coord = min(A_start, A_end, B_start, B_end)
        max_coord = max(A_start, A_end, B_start, B_end)
        whole_seq = get_sequence(filename, header, min_coord, max_coord, False)
        if args.extend:
            ext_min_coord = min_coord - args.extend
            if ext_min_coord < 0:
                ext_min_coord = 1
            ext_max_coord = max_coord + args.extend
            extended_seq = get_sequence(filename, header, ext_min_coord, ext_max_coord, False)

        seq_info.update({key:[transposon, extended_seq]})

        total_length = max(A_start, A_end, B_start, B_end) - \
                       min(A_start, A_end, B_start, B_end) + 1
        middle_start = sorted([A_start, A_end, B_start, B_end])[1]
        middle_end = sorted([A_start, A_end, B_start, B_end])[2]
        total_distance = middle_end - middle_start - 1

        with open(filename, "rt") as genome:
            organism = next(genome)
            organism = organism.split(',')[0].split('|')[-1].strip(" >.\n")

        candidate = "Candidate " + str(n)

        with open(outfile, "at") as out:
            if os.path.getsize(outfile) == 0:
                out.write("QUERY: " + organism + "\n")
            out.write("\n*" + candidate + "*\n\n")
            out.write("Length: " + str(round(total_length / 1000, 1)) + " kbp\n")
            out.write("Distance: " + str(round(total_distance / 1000, 1)) + " kbp\n")
            out.write(("Feature\tPosition\tStrand\tLength(bp)\tType"
                      ).expandtabs(20) + "\n")

            out.write(("A transposon\t" + str(A_start) + ".." + str(A_end) +
                       "\t" + A_strand + "\t" + str(A_len) + "\t" +
                       transposon).expandtabs(20) + "\n")
            out.write(("B transposon\t" + str(B_start) + ".." + str(B_end) +
                       "\t" + B_strand + "\t" + str(B_len) + "\t" +
                       transposon).expandtabs(20) + "\n")

            out.write("\n>A transposon\n")
            out.write(A_seq + "\n")
            out.write("\n>B transposon\n")
            out.write(B_seq + "\n")
            out.write("\n>whole sequence\n")
            out.write(whole_seq + "\n")
            if args.extend:
                out.write("\n>extended sequence\n")
                out.write(extended_seq)

            whole_seq_orfs, complete_seq_orfs = predict_orfs(whole_seq, min_coord)
            out.write("\n\nPredicted ORFs in the whole sequence:\n")
            for orf in whole_seq_orfs:
                out.write(str(orf[0]) + ".." + str(orf[1]) + "\n")

            if args.extend:
                extended_seq_orfs, complete_extended_orfs = \
                            predict_orfs(extended_seq, ext_min_coord)
                out.write("\nPredicted ORFs in the extended sequence:\n")
                corrected_ext_orfs = [(x-1, y-1) for (x, y) in complete_extended_orfs]
                for orf in corrected_ext_orfs:
                    out.write(str(orf[0]) + ".." + str(orf[1]) + "\n")

            if args.gbk:
                candidate = candidate.replace(' ', '').lower() + ".gbk"
                if not args.extend:
                    name = args.out + '/' + basename + "_composite_" + candidate
                    write_gbk(whole_seq, whole_seq_orfs, name, organism,
                              A_start, A_end, A_strand,
                              B_start, B_end, B_strand,
                              transposon, min_coord)
                else:
                    name = args.out + '/' + basename + "_composite_extended_" + candidate
                    write_gbk(extended_seq, extended_seq_orfs, name, organism,
                              A_start, A_end, A_strand,
                              B_start, B_end, B_strand,
                              transposon, ext_min_coord)

        n += 1

    flanking = False
    if args.flanking:
        flanking = compare_flanking_regions(seq_info, args.extend, filename)
    if flanking:
        with open(outfile, "at") as out:
            for transposon in flanking:
                out.write("\n*NOTE*: composite " + transposon + " transposon " +
                          "duplicated in genome regions with different flanks.")


# Use a sequence and matched ORFs to generate a genbank file.
def write_gbk(sequence, matched_orfs, filename, organism,
              A_start, A_end, A_strand, B_start, B_end, B_strand,
              transposon, min_coord):

    date = strftime("%d-%b-%Y").upper()
    orfs = []
    features = []

    # Get the strands of predicted ORFs.
    for item in matched_orfs:
        if item[0] < item[1]:
            start = item[0]
            end = item[1]
            strand = 0
        else:
            start = item[1]
            end = item[0]
            strand = -1
        orfs.append((start, end, strand))

    # Get transposon strands.
    if A_strand == '+':
        A_strand = 0
    elif A_strand == '-':
        A_strand = -1
    if B_strand == '+':
        B_strand = 0
    elif B_strand == '-':
        B_strand = -1

    gbk_record = SeqRecord(Seq(sequence, IUPAC.unambiguous_dna),
                        description=organism+" predicted composite transposon",
                        annotations={"accession":'.', "version":'.',
                                     "organism":'.', "date":date,
                                     "data_file_division":"BCT"})

    # Create A transposon features.
    features.append(SeqFeature.SeqFeature(
            SeqFeature.FeatureLocation(
            A_start-min_coord, A_end-min_coord+1, strand=A_strand),
            type="misc_feature",
            qualifiers={"note":"insertion sequence " + transposon}))

    # Create B transposon features.
    features.append(SeqFeature.SeqFeature(
            SeqFeature.FeatureLocation(
            B_start-min_coord, B_end-min_coord+1, strand=B_strand),
            type="misc_feature",
            qualifiers={"note":"insertion sequence " + transposon}))

    # Create features for predicted ORFs.
    for item in orfs:
        features.append(SeqFeature.SeqFeature(
                SeqFeature.FeatureLocation(
                item[0]-1, item[1], strand=item[2]),
                type="CDS"))

    for item in features:
        gbk_record.features.append(item)

    SeqIO.write(gbk_record, filename, "gb")


# From sequences with extended borders, compare regions flanking composite
# transposons of the same type to check for different genomic contexts.
def compare_flanking_regions(seq_info, extension, filename):

    matches = {}
    for key in seq_info:
        transposon = seq_info[key][0]
        if transposon not in matches.keys():
            matches.update({transposon:[tuple([key] + seq_info[key][1:])]})
        else:
            matches.update({transposon:matches[transposon]+[tuple([key] + seq_info[key][1:])]})

    for key in list(matches.keys()):
        if len(matches[key]) < 2:
            matches.pop(key, None)

    for sequence in matches:
        seq_1 = matches[sequence][0]
        seq_2 = matches[sequence][1]

        # Create a temporary fasta file with the A transposon left flanking region.
        with open(filename+"left1.fasta.tmp", "wt") as left1_fasta:
            left1_flank = seq_1[1][0:extension]
            left1_fasta.write(">left1\n")
            left1_fasta.write(left1_flank)

        # Create a temporary fasta file with the B transposon left flanking region.
        with open(filename+"left2.fasta.tmp", "wt") as left2_fasta:
            left2_flank = seq_2[1][0:extension]
            left2_fasta.write(">left2\n")
            left2_fasta.write(left2_flank)

        # Create a temporary fasta file with the B transposon right flanking region.
        with open(filename+"right1.fasta.tmp", "wt") as right1_fasta:
            right1_flank = seq_1[1][len(seq_1)-extension:]
            right1_fasta.write(">right1\n")
            right1_fasta.write(right1_flank)

        # Create a temporary fasta file with the B transposon right flanking region.
        with open(filename+"right2.fasta.tmp", "wt") as right2_fasta:
            right2_flank = seq_2[1][len(seq_2)-extension:]
            right2_fasta.write(">right2\n")
            right2_fasta.write(right2_flank)

        # Align left flanking regions and check if they are similar.
        subprocess.call(["blastn",
                        "-query", filename+"left1.fasta.tmp",
                        "-subject", filename+"left2.fasta.tmp",
                        "-out", filename+"left1×2.blastn.tmp",
                        "-outfmt", "6 qseqid qstart qend sseqid slen sstart send pident length",
                        "-dust", "no",
                        "-word_size", "7",
                        "-evalue", "1e-5"])

        # Align right flanking regions and check if they are similar.
        subprocess.call(["blastn",
                        "-query", filename+"right1.fasta.tmp",
                        "-subject", filename+"right2.fasta.tmp",
                        "-out", filename+"right1×2.blastn.tmp",
                        "-outfmt", "6 qseqid qstart qend sseqid slen sstart send pident length",
                        "-dust", "no",
                        "-word_size", "7",
                        "-evalue", "1e-5"])

        with open(filename+"left1×2.blastn.tmp", "rt") as left_blastn, open(
                  filename+"right1×2.blastn.tmp", "rt") as right_blastn:
            left_different = True
            right_different = True

            for line in left_blastn:
                line = line.split("\t")
                pident = float(line[7])
                pcvg = int(line[8])*100 / int(line[4])
                length = int(line[8])

                if length >= (extension/2) and pcvg >= args.coverage and pident >= args.identity:
                    left_different = False

            for line in right_blastn:
                line = line.split("\t")
                pident = float(line[7])
                pcvg = int(line[8])*100 / int(line[4])
                length = int(line[8])

                if length >= (extension/2) and pcvg >= args.coverage and pident >= args.identity:
                    right_different = False

        for tmp in [filename+"left1.fasta.tmp",
                    filename+"left2.fasta.tmp",
                    filename+"left1×2.blastn.tmp",
                    filename+"right1.fasta.tmp",
                    filename+"right2.fasta.tmp",
                    filename+"right1×2.blastn.tmp"]:
            os.remove(tmp)

        if left_different or right_different:
            return matches
        else:
            return False


# Verify input and call functions accordingly.
def main():

    if args.flanking and not args.extend:
        name = os.path.basename(sys.argv[0])
        print(name + ": error: please provide a length to extend")
        sys.exit()

    try:
        print(strftime("%c" + "\tStarting analyses."))

        # Create necessary directories if they do not already exist.
        if not os.path.exists(args.out):
            os.mkdir(args.out)
        if not os.path.exists(args.out + "/blastn"):
            os.mkdir(args.out + "/blastn")

        # File the command line to better allow reproducibility.
        with open(args.out + "/info.txt", "wt") as info:
            info.write("time:\t" + strftime("%c") + "\n")
            info.write("working directory:\t" + os.getcwd().replace(' ', "\ ") + "\n")
            info.write("command line:\t" + subprocess.list2cmdline(sys.argv[0:]) + "\n")

        # Run multithreaded BLASTn analyses if no output files are found.
        Pool(args.processors).map(partial(run_blastn, outdir=args.out), args.files)

        # Parse BLASTn results and check for potential composite transposons.
        Pool(args.processors).map(partial(parse_blastn, outdir=args.out), args.files)

        print(strftime("%c" + "\tDone."))

    except KeyboardInterrupt:
        print(strftime("%c") + "\tProgram interrupted.")
    #except:
        #print(strftime("%c") + "\tError. Please check the input files.")


if __name__ == "__main__":
    main()
