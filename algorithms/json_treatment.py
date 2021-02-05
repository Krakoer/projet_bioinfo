import json
from tqdm import tqdm
import os
import sys
from pathlib import Path
import re

canonPairs = {'A': ['U'], 'U': ['A', 'G'], 'G': ['C', 'U'], 'C': ['G']}


def fix_par(s):
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


def isCanon(link):
    assert('bp' in link and 'LW' in link)
    bp1 = link['bp'][0].capitalize()
    bp2 = link['bp'][2].capitalize()
    try:
        answer = bp2 in canonPairs[bp1] and link["LW"] == 'cWW'
        return answer, True
    except:
        return False, False


def main():
    json_pathlist = [str(path) for path in Path("./json_files").glob('**/*.json')]
    log_file = open('./log-json-treatment.txt', 'w')

    if not os.path.exists('./xdbn_files'):
        os.popen("mkdir xdbn_files").close()

    print("Start Json treating : creating xdbn files...")
    for json_path in tqdm(json_pathlist):
        name = json_path.split('/')[-1].split('.')[0]
        with open(json_path) as f:
            try:
                data = json.load(f)
            except:
                log_file.write(json_path+"-JSON opening failure"'\n')
                continue

        f.close()
        out_file = open("./xdbn_files/"+name+".xdbn", "w")

        try:
            pairs = data["pairs"]
        except:
            log_file.write(json_path+"-pairs not in keys"'\n')
            continue

        try:
            chains = {chain_name[-1] : [chain['bseq'], chain['sstr']] for chain_name, chain in data['chains'].items()}
        except:
            log_file.write(json_path+"-Pb while parsing chains, check if chains key exixsts and if sstr and bseq keys are in chains objects"'\n')
            continue

        try:
            nts = {nt['nt_id'] : [nt['chain_name'], nt['index_chain']] for nt in data['nts']}
        except:
            log_file.write(name+"-Pb while parsing nts : check if nts is in keys and if index_chain, chains_name and nt_is are in nts objects"'\n')
            continue

        header = ""
        for i, pair in enumerate(pairs):
            canon, success = isCanon(pair)
            if(not success):
                log_file.write(f"{json_path}-Error at pair index {i+1} : must be a DNA strand"'\n')
                continue
            if not canon:
                try:
                    id_nt1 = pair['nt1']
                    id_nt2 = pair['nt2']
                    nt1 = nts[id_nt1]
                    nt2 = nts[id_nt2]
                    chain_nb1 = list(chains.keys()).index(nt1[0])+1
                    chain_nb2 = list(chains.keys()).index(nt2[0])+1
                    header += f"{chain_nb1}.{nt1[1]}-{chain_nb2}.{nt2[1]};"
                except ValueError:
                    continue

        out_file.write(header+'\n')
                
        for chain in chains.values():
            out_file.write(chain[0])
            out_file.write('\n')
            out_file.write(fix_par(chain[1]))
            out_file.write('\n')

        out_file.close()
    log_file.close()

    print("xdbn files treated succesfuly")
    


if __name__ == '__main__':
    main()