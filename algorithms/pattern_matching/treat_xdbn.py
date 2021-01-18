import sys
import os
from pathlib import Path
from tqdm import tqdm
import re

debug = True

def draw_structure(target, intervals, nonCanon, path, varna):
    '''
    target : dot bracket seq
    bps : [[[1, 2], [5, 10]], ]
    path : outputpath
    vanra : path to varna
    '''
    if not Path(varna).exists():
        print("Cannot find {}".format(varna))
        return

    seq = " "*len(target.strip())

    cmd = "java -cp {} fr.orsay.lri.varna.applications.VARNAcmd -sequenceDBN \"{}\" -structureDBN \"{}\" -o {}".format(
        varna, seq, target.strip(), path)
    cmd += " -bpStyle simple -algorithm radiate -resolution 6.0"

    auxhighlight = []
    colors = ["#FF0000", "#00FF00", "#0000FF"]

    for i in range(len(intervals)): #For each motif
        for j in range(len(intervals[i])): #For each occurence of the motif
            for k in range(len(intervals[i][j])): # For each intervl of the occ
                auxhighlight.append(f"{intervals[i][j][k][0]}-{intervals[i][j][k][1]}:fill=#FFFFFF,outline={colors[i%3]},radius={16+j*2}")

    cmd += " -highlightRegion \"{}\"".format(";".join(auxhighlight))

    auxbps = ["({},{}):color=#EE82EE".format(bp[0], bp[1]) for bp in nonCanon]
    cmd += " -auxBPs \"{}\"".format(";".join(auxbps))

    if(debug):
        print(cmd)
    os.popen(cmd).close()

def fix_parenthesis(s):
    """
    Fix the parenthesis notation if broken
    """
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


def find_matching_par(start, seq):
    count = 0

    for i in range(start-1, len(seq)):
        if seq[i] == '(':
            count += 1
        elif seq[i] == ')':
            count -= 1

        if(count == 0):
            return i+1


def find_nt_in_chain(chain, pos, nts):
    for i, nt in enumerate(nts):
        if(nt['chain'] == chain and nt['pos_chain'] == pos):
            return i


def break_chain(chain, chain_id, nts):
    subs = []  # Array to store all subchains in db notation
    sub = ""  # Store subchains
    subs_count = 1  # Counts the number of subchains
    pos_count = 1 # Count the overall position in the chain (used to find nucleotide)
    sub_pos = 1  # Position in the current subchain
    for c in chain:
        if c == '&':
            subs.append(fix_parenthesis(sub))
            subs_count += 1
            sub = ""
            sub_pos = 1
        else:
            nts[find_nt_in_chain(chain_id, pos_count, nts)]['subchain'] = subs_count
            nts[find_nt_in_chain(chain_id, pos_count, nts)]['pos_sub'] = sub_pos
            sub += c
            sub_pos += 1
            pos_count += 1

    subs.append(fix_parenthesis(sub))

    return subs


def parse_header(header):
    if header == "":
        return []
    positions = header.strip().split(';')[:-1]
    pairs = []

    for position in positions:
        nucs = position.split('-')
        chain1 = int(nucs[0].split('.')[0])
        pos1 = int(nucs[0].split('.')[1])
        chain2 = int(nucs[1].split('.')[0])
        pos2 = int(nucs[1].split('.')[1])
        pairs.append([[chain1, pos1], [chain2, pos2]])
    return pairs


def create_nts(dbns):
    """
    Create a list of nts represented as {chain subchain pos} given a list of chains in db notation
    """
    nts = []

    for chain_num, chain in enumerate(dbns):
        pos = 1
        for c in chain:
            if(c != '&'):
                nts.append({'chain': chain_num+1, 'subchain': -
                            1, 'pos_chain': pos, 'pos_sub': -1})
                pos += 1
    return nts


def parse_xdbn(xdbn_path):
    xdbn_file = open(xdbn_path, 'r')
    lines = xdbn_file.readlines()
    xdbn_file.close()
    header = lines[0].strip()
    nonCanonPairs_header = parse_header(header)

    dbns = [lines[2*i].strip() for i in range(1, len(lines)//2+1)]
    nts = create_nts(dbns)

    nonCanonPairs = []
    for pair in nonCanonPairs_header:
        if find_nt_in_chain(pair[0][0], pair[0][1], nts) is None or find_nt_in_chain(pair[1][0], pair[1][1], nts) is None:
            continue
        nonCanonPairs.append((find_nt_in_chain(pair[0][0], pair[0][1], nts), find_nt_in_chain(pair[1][0], pair[1][1], nts)))

    chains = []
    for i, ligne in enumerate(dbns):
        chains.append(break_chain(ligne, i+1, nts))

    return chains, nts, nonCanonPairs


def parse_output(output_path):
    output_file = open(output_path, 'r')
    lines = output_file.readlines()
    output_file.close()

    res = {}

    for line in lines:
        if(len(line.split('\t')) > 2):
            head = line.split('\t')[0]
            chain_nb = int(head.split('-')[1])
            sub_nb = int(head.split('-')[2])
            res[(chain_nb, sub_nb)] = []
            motifs_found = line.split('\t')[1:-1]
            for motif_found in motifs_found:
                motif = motif_found.split(':')[0]
                positions_string = motif_found.split(':')[-1]
                positions = [int(num)
                             for num in positions_string.split(';')[:-1]]
                res[(chain_nb, sub_nb)].append(
                    {'motif': motif, 'positions': positions})
    return res

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

    intervals.append([start_sequence, pos_sequence-1])
    return intervals

def check_nonCanon(nonCanon, intervals, log, name):
    for i in range(len(intervals)):
        for j in range(len(intervals[i])): # For each occurence of a motif
            res = False
            for k in range(len(intervals[i][j])): 
                inter = intervals[i][j][k]
                for nt1, nt2 in nonCanon:
                    if nt1 >= inter[0] and nt1 <= inter[1] or nt2 >= inter[0] and nt2 <= inter[1]:
                        res = True
            if not res:
                log.write(f'Motif without nonCanon found in RNA {name}\n')

                

def treat_ARN(name, out_path, xdbn_path, output_path):
    chains, nts, nonCanonPairs = parse_xdbn(xdbn_path)
    output_data = parse_output(out_path)
    log = open('motif_without_canon.txt', 'a')
    for chain_subchain, motifs_data in output_data.items():
        chain = chains[chain_subchain[0]-1][chain_subchain[1]-1] # The chain in which the motifs occur
        intervals = []
        for motif_data in motifs_data: # For each dict representing a motif and its positions
            motif = motif_data['motif']
            positions = motif_data['positions']
            inter_motif = []
            for pos in positions:
                inter_motif.append(find_intervals(motif, chain, pos))
            intervals.append(inter_motif)

        nonCanonPairs_subchain = [(nts[nt_pair[0]]['pos_sub'], nts[nt_pair[1]]['pos_sub']) for nt_pair in nonCanonPairs if nts[nt_pair[0]]['chain'] == chain_subchain[0] and nts[nt_pair[0]]['subchain'] == chain_subchain[1] and nts[nt_pair[1]]['chain'] == chain_subchain[0] and nts[nt_pair[1]]['subchain'] == chain_subchain[1]]
        
        draw_structure(chain, intervals, nonCanonPairs_subchain,  f'{output_path}/{name}-{chain_subchain[0]}-{chain_subchain[1]}.png', '/home/axel/Téléchargements/VARNAv3-93.jar')
        check_nonCanon(nonCanonPairs_subchain, intervals, log, name)
    log.close()

def test():
    out_path = '../../RNA_files/outputs/1AQO.out'
    xdbn_path = "../../RNA_files/xdbn/1AQO.xdbn"
    name = "1AQO"
    treat_ARN(name, out_path, xdbn_path, "bite")
    # chain = ".....(.((..(((...)..).)..))..."
    # motif = "((*)..)"
    # print(find_matching_par(14, chain))
    # print(find_intervals(motif, chain, 13))

def is_empty(path):
    f = open(path, 'r')
    return len(f.readlines()) == 0

def main():
    if len(sys.argv) < 4:
        print(f"Usage : python {sys.argv[0]} output_files xdbn_files output_path")
        return 
    
    already_treated = []
    if len(sys.argv) == 5:
        path_already = sys.argv[4]
        already_treated = [str(path).split('/')[-1].split('.')[0] for path in Path(path_already).glob('**/*.png')]
        
    output_pathlist = [str(path)
                       for path in Path(sys.argv[1]).glob('**/*.out')]
    xdbn_pathlist = [str(path) for path in Path(sys.argv[2]).glob('**/*.xdbn')]
    log_file = open('./log.txt', 'w')
    print(len(already_treated))

    for out_path in tqdm(output_pathlist):
        name = out_path.split('/')[-1].split('.')[0]
        if name in already_treated:
            continue
        if(debug):
            print(name)
        if(is_empty(out_path)):
            continue

        xdbn_paths = list(filter(lambda x: name in x, xdbn_pathlist))
        if len(xdbn_paths) > 0:
            xdbn_path = xdbn_paths[0]
            if(is_empty(xdbn_path)):
                continue
            treat_ARN(name, out_path, xdbn_path, sys.argv[3])
        else:
            continue
    log_file.close()

if __name__ == '__main__':
    # test()
    main()