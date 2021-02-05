all: dotbpattern | out_files pdf_files
	python create_json.py ./PDB_files
	python json_treatment.py
	./dotbpattern ./xdbn_files ./motifs/all_p_0.01.csv ./out_files
	python print_results.py ./out_files ./xdbn_files ./images

dotbpattern: dotbpattern.cpp
	g++ -o dotbpattern dotbpattern.cpp -lboost_system -lboost_filesystem -L /usr/local/lib

out_files:
	mkdir out_files

pdf_files:
	mkdir pdf_files