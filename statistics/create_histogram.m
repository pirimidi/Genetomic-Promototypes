%----------------------------------------------------------
% Author: Mirko Palla.
% Date: August 4, 2006.
% For: Task 6 at the Church Lab - Genetics Department, 
% Harvard Medical School.
% 
% Purpose: This program takes two text files containing a list of integers) 
% and displays the numeric distrobution via a histogram with bin-size of 
% 10 and orientation types (++, --, -+, +-) respectively in MATLAB.
% 
% Input: two text files containing list of integers
% 
% Output: two histograms displaying relative distance and 
% orientation for all TF pairs in the yeast genome
%----------------------------------------------------------- 

% 0. Declare function and check IO validity.

function create_histogram(file1, file2)

if nargin < 2
  error('Error: not enough input args - usage: [distance_file], [orientation_file]\n');
end

% 1. Open directory ".txt" files are located.

cd('/root/Church_lab/Polypromoter/task_6/project_14/statistics/');

%---------- Relative TF distance data ----------

% 2. Retrieve first argument and check for validity.

fid= fopen(file1, 'r');
  if fid == -1
    disp('Error: file cannot be opened - file1\n');
end

% 3. Read in relative distance data into a 1-D array, hist_1.

hist_1=load(file1);

% 4. Assign the bin size for 1000 bp region, i.e., promoter.

x=0:10:1000;

% 5. Create histogram of relative distance data.

hist(hist_1,x);

% 6. Set histogram axis parameters.

set(gca,'XTick', 0:10:1000);
set(gca,'XTicklabel', {''});
set(gca,'XTickLabel',0:1:9);

% 7. Set histogram title/label parameters.

xlabel({''; 'Relative transcription factor binding site distance'; ''; '[Bins with increments of 10 bp (up to 1000 bp)]'});
ylabel('Relative-distance bin frequency');
grid;

title(['\bf\fontsize{12}Relative transcription factor binding site distances in all S. cerevisiae promoter regions']);

plotedit on;

%---------- Relative TF order data ----------

% 8. Create new display window for second histogram.

figure;

% 9. Retrieve second argument and check for validity.

fid= fopen(file2, 'r');
  if fid == -1
    disp('Error: file cannot be opened - file2\n');
end

% 10. Read in relative distance data into a 1-D array, hist_2.

hist_2=load(file2);

% 11. Assign the bin size for 1000 bp region, i.e., promoter.

y=1:1:4;

% 12. Create histogram of relative distance data.

hist(hist_2,y);

% 13. Set histogram axis parameters.

set(gca,'XTick', 1:1:4);
set(gca,'XTicklabel', {''});
set(gca,'XTickLabel',{'++';'+-';'-+';'--'});

% 14. Set histogram title/label parameters.

xlabel('Orientation types for transcription factor binding site pairs');
ylabel('Relative orientation type frequency');
grid;

title(['\bf\fontsize{12}Relative transcription factor binding site pair orientation in all S. cerevisiae promoter regions']);

% 15. Change the bin color to red.
 
h = findobj(gca,'Type','patch');
set(h,'FaceColor','r');

plotedit on;
