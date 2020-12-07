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
    cmd += " -bpStyle simple -algorithm line -resolution 4"

    auxhighlight = []

    for i in range(len(bps)):
        for j in range(len(bps[i])):
            auxhighlight.append(
                f"{bps[i][j][0]}-{bps[i][j][1]}:fill=#FFFFFF,outline={colors[i%3]},radius={16+j*2}")

    cmd += " -highlightRegion \"{}\"".format(";".join(auxhighlight))

    # print(cmd)
    os.popen(cmd).close()


def find_end_pattern(start, seq):
    count = 0

    for i in range(start, len(seq)):
        if seq[i] == '(':
            count += 1
        elif seq[i] == ')':
            count -= 1

        if(count == 0):
            return i


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


def treat_motif(dbn_path, out_path, output_path):
    output = open(out_path)
    lignes = output.readlines()
    output.close()

    if(any([contains_motif(l) for l in lignes])):
        dbn = open(dbn_path)
        temp = dbn.readlines()
        structure = temp[2].strip()
        dbn.close()

        motifs = parse_seq(structure)

        assert(len(motifs) == len(lignes))

        for i, ligne in enumerate(lignes):
            if contains_motif(ligne):
                name = ligne.strip().split('\t')[0]
                motifs_occurence = ligne.strip().split('\t')[1:]
                motifs_positions = []

                for occurence in motifs_occurence:
                    positions = occurence.strip('()').split(',')[1:]

                    motifs_positions.append(
                        [[int(pos), find_end_pattern(int(pos), motifs[i])] for pos in positions])

                # print(motifs_positions)

                draw_structure(
                    motifs[i], motifs_positions, f'{output_path}/{name}.png', '/home/axel/Téléchargements/VARNAv3-93.jar')


def main():
    if len(sys.argv) < 4:
        print(
            f"Usage : python {sys.argv[0]} output_files dbn_files output_path")
        return

    output_pathlist = [str(path)
                       for path in Path(sys.argv[1]).glob('**/*.out')]
    dbn_pathlist = [str(path) for path in Path(sys.argv[2]).glob('**/*.dbn')]

    for out_path in tqdm(output_pathlist):
        name = out_path.split('/')[-1].split('.')[0]

        dbn_paths = list(filter(lambda x: name in x, dbn_pathlist))
        if len(dbn_paths) > 0:
            dbn_path = dbn_paths[0]

            treat_motif(dbn_path, out_path, sys.argv[3])
        else:
            continue


if __name__ == '__main__':
    main()

    # because path is object not string

    # print(path_in_str)

# # pattern_file = open("all_sub_1.txt")
# # pattern = pattern_file.readlines()[44]
# # pattern_file.close()

# pattern_pos = 5

# draw_structure(structure, [[5, find_end_pattern(5, structure)], [10, 50]], './ouput_varna.png',
#                '/home/axel/Téléchargements/VARNAv3-93.jar')
