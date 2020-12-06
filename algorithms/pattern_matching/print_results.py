import sys
import os
from pathlib import Path

# 6WLJ


def hex_to_RGB(hex):
    ''' "#FFFFFF" -> [255,255,255] '''
    # Pass 16 to the integer function for change of base
    return [int(hex[i:i+2], 16) for i in range(1, 6, 2)]


def RGB_to_hex(RGB):
    ''' [255,255,255] -> "#FFFFFF" '''
    # Components need to be integers for hex to make sense
    RGB = [int(x) for x in RGB]
    return "#"+"".join(["0{0:x}".format(v) if v < 16 else
                        "{0:x}".format(v) for v in RGB])


def color_dict(gradient):
    ''' Takes in a list of RGB sub-lists and returns dictionary of
      colors in RGB and hex form for use in a graphing function
      defined later on '''
    return {"hex": [RGB_to_hex(RGB) for RGB in gradient],
            "r": [RGB[0] for RGB in gradient],
            "g": [RGB[1] for RGB in gradient],
            "b": [RGB[2] for RGB in gradient]}


def linear_gradient(start_hex, finish_hex="#FFFFFF", n=10):
    ''' returns a gradient list of (n) colors between
      two hex colors. start_hex and finish_hex
      should be the full six-digit color string,
      inlcuding the number sign ("#FFFFFF") '''
    # Starting and ending colors in RGB form
    s = hex_to_RGB(start_hex)
    f = hex_to_RGB(finish_hex)
    # Initilize a list of the output colors with the starting color
    RGB_list = [s]
    # Calcuate a color at each evenly spaced value of t from 1 to n
    for t in range(1, n):
        # Interpolate RGB vector for color at the current value of t
        cv = [int(s[j] + (float(t)/(n-1))*(f[j]-s[j])) for j in range(3)]
        # Add it to our list of output colors
        RGB_list.append(cv)

    return color_dict(RGB_list)


def draw_structure(target, bps, path, varna):
    '''
    target : dot bracket seq
    bps : ??
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

    # colors = linear_gradient("#ffcccc", "#ff0000", len(bps)+1)["hex"][1:]
    # auxbps = ["({},{}):color={}".format(bp[0]+1, bp[1]+1, c)
    #           for bp, c in zip(bps, colors)]
    # cmd += " -auxBPs \"{}\"".format(";".join(auxbps))
    cmd += ' -highlightRegion "10-16:fill=#FF0000"'
    print(cmd)
    os.popen(cmd).close()


def find_end_pattern(start, pattern, seq):
    bp_nb = 0
    pos = start
    for car in pattern:
        if(car == '*'):
            while(seq[pos])


dbn = open("../dot_brackets/6WLJ.dbn")
lignes = dbn.readlines()
structure = lignes[2]
sequence = lignes[1]
dbn.close()

pattern_file = open("./all_sub_1.txt")
pattern = pattern_file.readlines()[44]
pattern_file.close()

pattern_pos = 5

draw_structure(structure, [], './ouput_varna.png',
               '/home/axel/Téléchargements/VARNAv3-93.jar')
