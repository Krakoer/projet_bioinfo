import moderna

m = load_model('../pdb_files/1F1T.pdb', 'A')  # A is the chain
seq = m.get_sequence()
sec = m.get_secstruc()
