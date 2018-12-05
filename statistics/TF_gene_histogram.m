%----------------------------------------------------------
% Author: Mirko Palla.
% Date: August 8, 2006.
% For: Task 6 at the Church Lab - Genetics Department, 
% Harvard Medical School.
% 
% Purpose: This program takes a text files containing a list 
% of integers and displays the numeric distrobution via a 
% histogram with bin-size of 10 in MATLAB.
% 
% Input: a text file containing list of integers
% 
% Output: a histogram displaying TF-to-gene distance and 
% for all TF pairs in the yeast genome
%----------------------------------------------------------- 

% 0. Declare function and check IO validity.

function TF_gene_histogram(d_file)

if nargin < 1
  error('Error: not enough input args - usage: [distance_file]\n');
end

% 1. Open directory ".txt" files are located.

cd('/root/Church_lab/Polypromoter/task_6/project_14/statistics/');

%---------- Relative TF distance data ----------

% 2. Retrieve first argument and check for validity.

fid= fopen(d_file, 'r');
  if fid == -1
    disp('Error: file cannot be opened - d_file\n');
end

% 3. Read in TF-to-gene distance data into a 1-D array, hist_d.

hist_d=load(d_file);

% 4. Assign the bin size for 1000 bp region, i.e., promoter.

x=0:10:1000;

% 5. Create histogram of relative distance data.

hist(hist_d,x);

% 6. Set histogram axis parameters.

set(gca,'XTick', 0:10:1000);
set(gca,'XTicklabel', {''});
set(gca,'XTickLabel',0:1:9);

% 7. Set histogram title/label parameters.

xlabel({''; 'Transcription factor binding site distance to gene start'; ''; '[Bins with increments of 10 bp (up to 1000 bp)]'});
ylabel('TF-gene distance bin frequency');
grid;

title(['\bf\fontsize{12}Transcription factor binding site distances to gene start in all S. cerevisiae promoter regions']);

% 8. Change the bin color to green.
 
h = findobj(gca,'Type','patch');
set(h,'FaceColor','g');

plotedit on;
