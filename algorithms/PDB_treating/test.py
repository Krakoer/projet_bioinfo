chains = {
    "chain_X": {
      "num_nts": 16,
      "bseq": "TGAGGTAGTAGGTTGT",
      "sstr": ".((((((((((((((.",
      "form": "...AAB.AA...B.B-",
      "helical_rise": 2.916,
      "helical_rise_std": 1.408,
      "helical_axis": [
        -0.993,
        0.075,
        -0.09
      ],
      "point1": [
        10.327,
        -25.34,
        -18.56
      ],
      "point2": [
        -33.835,
        -22.004,
        -22.563
      ],
      "num_chars": 46,
      "suite": "T!!G!!A1aG1aG1bT!!A!!G1aT1cA!!G!!G!!T!!T!!G!!T"
    },
    "chain_Y": {
      "num_nts": 14,
      "bseq": "CAACCUACUACCUC",
      "sstr": "))))))))))))))",
      "form": "AAAAAAAAAAAAA-",
      "helical_rise": 3.177,
      "helical_rise_std": 0.452,
      "helical_axis": [
        0.961,
        -0.011,
        -0.275
      ],
      "point1": [
        -28.677,
        -24.495,
        -15.235
      ],
      "point2": [
        11.26,
        -24.969,
        -26.68
      ],
      "num_chars": 40,
      "suite": "C__A1aA1fC1aC1aU1aA1cC1aU1aA1cC!!C1LU1aC"
    }
  }
new_chains = chains = {chain_name[-1] : [chain['bseq'], chain['sstr']] for chain_name, chain in chains.items()}
print(list(new_chains.keys()).index('X'))