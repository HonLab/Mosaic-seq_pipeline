#!/bin/tcsh

foreach counter (1 2 3 4 5 6)
#foreach counter (3)

	mkdir sorted_dir_$counter

	cd sorted_dir_$counter
	
	mv ../log$counter.sh .
	cp ../*pl .
	cp ../*m .
	sbatch log$counter.sh

	cd ../

end
