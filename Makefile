all:

add:
	bash bash/add.sh
run:
	python	src/mainpb.py
prova:
	python src/thesis_discrepancy_compute.py
import-from-Freefem++:
	cp ~/Documents/FreeFem++cs/out/*.txt  ./runs/ff++/

