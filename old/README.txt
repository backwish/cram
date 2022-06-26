./mkcram filename [b] [sbs] [bu]
makes an index for the file "filename".
Output filenames are cram1.dat and cram2.dat.
[b] is the size of large blocks (default: 1024)
[sbs] is the size of small blocks (default: 64)
[bu] is a parameter adjusting write speed and fitness to entropy changes (default: 4)
If [bu] is small, faster and slower adaptation.
This may not work for small files.


./cram_test filename cram1.dat cram2.dat
rewrites the index cram1.dat and cram2.dat by overwriting the file "filename".

./crambg_test filename1 filename2
constructs an index for filename1 using gzip, and overwrites filename2.
