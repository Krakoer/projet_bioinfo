import sys
import os
from pathlib import Path
from tqdm import tqdm

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
    
def find_nt_in_chain(chain, pos, nts):
  for i, nt in enumerate(nts):
    if(nt['chain'] == chain and nt['pos_chain'] == pos):
      return i

def break_chain(chain, chain_id, nts):
  subs = [] # Array to store all subchains in db notation
  sub = "" # Store subchains
  subs_count = 1 # Counts the number of subchains
  pos_count = 1 # Count the overall position in the chain (used to find nucleotide)
  sub_pos = 1 # Position in the current subchain
  for c in chain:
    if c == '&':
      subs.append(fix_parenthesis(sub))
      subs_count +=1
      sub = ""
      sub_pos = 1
    else:
      nts[find_nt_in_chain(chain_id, pos_count, nts)]['subchain'] = subs_count
      nts[find_nt_in_chain(chain_id, pos_count, nts)]['pos_sub'] = sub_pos
      sub += c
      sub_pos+=1
      pos_count+=1

  subs.append(fix_parenthesis(sub))

  return subs


def parse_header(header):
    if header == "":
        return []
    positions = header.strip().split(';')[:-1]
    pairs = []

    for position in positions:
        nucs = position.split('-')
        chain1 = nucs[0].split('.')[0]
        pos1 = nucs[0].split('.')[1]
        chain2 = nucs[1].split('.')[0]
        pos2 = nucs[1].split('.')[1]
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
          nts.append({'chain': chain_num+1, 'subchain': -1, 'pos_chain': pos, 'pos_sub': -1})
          pos+=1
    return nts

def parse_xdbn(xdbn_path):
    xdbn_file = open(xdbn_path, 'r')
    lines = xdbn_file.readlines()
    header = lines[0]
    nonCanonPairs_header = parse_header(header)

    dbns = [lines[2*i] for i in range(1, len(lines)//2+1)]
    nts = create_nts(dbns)

    nonCanonPairs = []
    for pair in nonCanonPairs_header:
        nonCanonPairs.append((find_nt_in_chain(pair[0][0], pair[0][1], nts), find_nt_in_chain(pair[1][0], pair[1][1], nts)))

    chains = [] 
    for i, ligne in enumerate(dbns):
        chains.append(break_chain(ligne, i+1, nts))
    

    return chains, nts, nonCanonPairs

def parse_output(output_path):
    pass