Author: A. Jo Williams-Newkirk
Date of last changes: 07/10/2015
Email: IGY7@CDC.GOV

Applies to Prokka version 1.8.

Use case: I want to run Prokka on a single fasta file.
In the terminal: qsub -M <email> prokkasingle.sh <input fasta>
qsub submits to the cluster.
-M provides Sun Grid Engine (SGE; the cluser manager) an email to send a notification to when your job is done. Replace <email> with your email (no <>).
prokkasingle.sh is the name of the script you're running.
<input fasta> is the name of the fasta file you want to run Prokka on. Replace <input fasta> with the name (no <>).
A note on file locations: Linux looks in the directory that the terminal is currently "in" for files without paths specified. If they aren't there, it will say "file not found" and fail. If your files are inside a sub-directory of your current directory, you can tell Linux where to look from where you're at (eg. the terminal is in /scicomp/home/igy7/folder1 and I'm referring to a file in folder1a, which is inside of folder1, so I can use folder1a/myscript.sh instead of giving the absolute path starting with /scicomp). If your file is in another directory outside of the one you're in, then you must give the full, absolute path to it starting with /scicomp. Change directories in the terminal using cd </path/to/directory>. The tab key completes partial directory and file names. Cntrl + L from an explorer window opens up a window with the full path to the directory you're in.
If you place these scripts in the same directory as your fasta file, the run will fail with an error message in the error file about not being able to create a directory. That's because the folder is named from the input file, and Linux will balk at creating a folder and a file in the same directory with the exact same names. Make sure your script is somewhere different from your input file. 
Where are my output files? A folder was created in the directory where your script was located.

Use case: I have a whole folder of fasta files and I want to run Prokka on each of them. All at once. Really fast.
In the terminal: ./prokkabatch.sh <input directory> <email>
All of the same stipulations from the first use case still apply. Note that prokkabatch.sh works by submitting prokkasingle.sh to the cluster once for each file in your target directory. So you need both scripts and qstat should show you as many jobs as you have input files after running prokkabatch.sh.
