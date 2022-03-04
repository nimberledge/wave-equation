make clean;
make;
if [ "$#" -ne 1 ]; then
	mpiexec -n 9 bin/parwave src/conf_file.txt && python3 tools/collate.py 9;
else
	mpiexec -n $1 bin/parwave src/conf_file.txt && python3 tools/collate.py $1;
fi
rm data/*txt;
echo "Cleared the data/ directory";
# rm images/img_dump/*png;
# echo "Cleared the images/img_dump/ directory";