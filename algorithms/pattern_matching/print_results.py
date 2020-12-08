import sys
import os
from pathlib import Path
from tqdm import tqdm

# 6WLJ


colors = ["#FF0000", "#00FF00", "#0000FF"]


def draw_structure(target, bps, path, varna):
    '''
    target : dot bracket seq
    bps : [[[1, 2], [5, 10]], ]
    path : outputpath
    vanra : path to varna
    '''
    if not Path(varna).exists():
        print("Cannot find {}".format(varna))
        return
    seq = " " * len(target)

    cmd = "java -cp {} fr.orsay.lri.varna.applications.VARNAcmd -sequenceDBN \"{}\" -structureDBN \"{}\" -o {}".format(
        varna, seq, target, path)
    cmd += " -bpStyle simple -algorithm radiate -resolution 4"

    auxhighlight = []

    for i in range(len(bps)):
        for j in range(len(bps[i])):
            auxhighlight.append(
                f"{bps[i][j][0]+1}-{bps[i][j][1]+1}:fill=#FFFFFF,outline={colors[i%3]},radius={16+j*2}")

    cmd += " -highlightRegion \"{}\"".format(";".join(auxhighlight))

    # print(cmd)
    os.popen(cmd).close()


def find_matching_par(start, seq):
    count = 0

    for i in range(start, len(seq)):
        if seq[i] == '(':
            count += 1
        elif seq[i] == ')':
            count -= 1

        if(count == 0):
            return i


def find_intervals(pattern, sequence, start):
    '''
    Takes a pattern like ((..).(*)), et db sequence, and the position of the pattenr in the seq
    Returns the intervals of the seq involved in the pattern
    '''
    start_sequence = start
    pos_sequence = start
    intervals = []

    for i in range(len(pattern)):
        if(pattern[i] == '*'):
            pos_sequence -= 1

            intervals.append([start_sequence, pos_sequence])
            pos_sequence = find_matching_par(pos_sequence, sequence)
            start_sequence = pos_sequence

        else:
            pos_sequence += 1

    intervals.append([start_sequence, pos_sequence-2])
    return intervals


def parse_seq(seq):
    def fix(s):
        pstack = []
        fixed = [c for c in s]

        for i, c in enumerate(s):
            if c == '(':
                pstack.append(i)
            elif c == ')':
                if len(pstack) == 0:
                    fixed[i] = '.'
                else:
                    pstack.pop()

        while len(pstack) > 0:
            fixed[pstack.pop()] = '.'

        return ''.join(fixed)
    chains = seq.strip().split('&')
    return [fix(c) for c in chains]


def contains_motif(line):
    return len(line.split('\t')) > 2


def treat_motif(dbn_path, out_path, motif_path, output_path):
    output = open(out_path)
    lignes_output = output.readlines()
    output.close()

    if(any([contains_motif(l) for l in lignes_output])):
        dbn = open(dbn_path)
        temp = dbn.readlines()
        structure = temp[2].strip()
        dbn.close()
        chaines = parse_seq(structure)

        motifs_file = open(motif_path)
        motifs = motifs_file.readlines()
        motifs_file.close()

        assert(len(chaines) == len(lignes_output))

        for i, ligne in enumerate(lignes_output):
            if contains_motif(ligne):
                name = ligne.strip().split('\t')[0]

                motifs_occurence = ligne.strip().split('\t')[1:]
                motifs_positions = []

                for occurence in motifs_occurence:
                    motif = motifs[int(occurence.strip().split(':')[0])]

                    positions = occurence.split(':')[-1].split(';')[:-1]

                    motifs_positions += [find_intervals(
                        motif, chaines[i], int(pos.split('-')[0])) for pos in positions]

                draw_structure(
                    chaines[i], motifs_positions, f'{output_path}/{name}.png', '/home/axel/Téléchargements/VARNAv3-93.jar')


def main():
    if len(sys.argv) < 5:
        print(
            f"Usage : python {sys.argv[0]} output_files dbn_files motifs_file output_path")
        return

    output_pathlist = [str(path)
                       for path in Path(sys.argv[1]).glob('**/*.out')]
    dbn_pathlist = [str(path) for path in Path(sys.argv[2]).glob('**/*.dbn')]

    for out_path in tqdm(output_pathlist):
        name = out_path.split('/')[-1].split('.')[0]

        dbn_paths = list(filter(lambda x: name in x, dbn_pathlist))
        if len(dbn_paths) > 0:
            dbn_path = dbn_paths[0]

            treat_motif(dbn_path, out_path, sys.argv[3], sys.argv[4])
        else:
            continue


if __name__ == '__main__':
    main()
