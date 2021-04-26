%% step1 : trim adaptor/adaptor dimer
clear all
Lib = 6;
Experiment = 21;
save(sprintf('E%d_L%d_ReadsLength.mat',Experiment, Lib),'Length');% a = dir(sprintf('G%d*.fastq',Lib));

for i = 1:length(a)
    filenames(i) = strread(a(i).name,'%s');
end

%write the forward reads and reverse reads inyo two files
fid_f = fopen(sprintf('E%d_TAdp_%d_R1.fastq',Experiment, Lib),'wt');
fid_r = fopen(sprintf('E%d_TAdp_%d_R2.fastq',Experiment, Lib),'wt');

%get the NGS data
fid1 = fopen(filenames{1});
fid2 = fopen(filenames{2});
currline_1 = fgets(fid1);  % get the first seq
currline_2 = fgets(fid2);
linecount = 1;

% Adaptor sequence
f_adp_seq = 'AGATCGGAAGAGCACACGTCTGAACTCC';
r_adp_seq = 'AGATCGGAAGAGCGTCGTGTAGGGAAAG';

tic
adp_dimer = 0;
seq = 0;
%loop all the seq in data
while ischar(currline_1)

    if mod(linecount,4) == 1 % line1: seq name
        find_adp = 0;
        f1 = currline_1;
        r1 = currline_2;

    elseif mod(linecount,4) == 2 % line2: reads
        f2 = currline_1;
        r2 = currline_2;
        f_adp_pos = strfind(currline_1,f_adp_seq);
        r_adp_pos = strfind(currline_2,r_adp_seq);

        if ~isempty(f_adp_pos) && ~isempty(r_adp_pos)
            fprintf(fid_f,'%s',f1);
            fprintf(fid_r,'%s',r1);
            fprintf(fid_f,'%s\n',currline_1(1:(f_adp_pos-1)));
            fprintf(fid_r,'%s\n',currline_2(1:(r_adp_pos-1)));
            find_adp = 1;
        end

    elseif mod(linecount,4) == 3 && find_adp == 1% line3: +
       fprintf(fid_f,'%s',currline_1);
       fprintf(fid_r,'%s',currline_2);
    elseif mod(linecount,4) == 0 && find_adp == 1 %line4: quality
       fprintf(fid_f,'%s\n',currline_1(1:(f_adp_pos-1)));
       fprintf(fid_r,'%s\n',currline_2(1:(r_adp_pos-1)));
    end

    currline_1 = fgets(fid1);
    currline_2 = fgets(fid2);
    linecount = linecount + 1;

    if mod(linecount,1000000) == 1
        linecount
        toc
    end
end
fprintf('total reads: %d\n',(linecount-1)/4);

%% step2: Quality control for pair end NGS raw data

% Input file
a = dir(sprintf('E%d_TAdp_%d*.fastq',Experiment, Lib));
for i = 1:length(a)
    filenames(i) = strread(a(i).name,'%s');
end

% write the forward reads and reverse reads into two files
fid_f = fopen(sprintf('E%d_Lib%d_QC_R1_Long.fastq',Experiment, Lib),'wt');
fid_r = fopen(sprintf('E%d_Lib%d_QC_R2_Long.fastq',Experiment, Lib),'wt');
fid_f2 = fopen(sprintf('E%d_Lib%d_QC_R1_short.fasta',Experiment, Lib),'wt');
fid_r2 = fopen(sprintf('E%d_Lib%d_QC_R2_short.fasta',Experiment, Lib),'wt');

% get the NGS data
fid1 = fopen(filenames{1});
fid2 = fopen(filenames{2});
currline_1 = fgets(fid1);  % get the first seq
currline_2 = fgets(fid2);
linecount = 1;
true_seq = 0;


Length = zeros(1,150);
tic
% %loop all the seq in data
while ischar(currline_1)
    if mod(linecount,4) == 1 % line1: seq name
        f1 = currline_1;
        r1 = currline_2;
    elseif mod(linecount,4) == 2 % line2: reads
        f2 = currline_1;
        r2 = currline_2;
    elseif mod(linecount,4) == 3 % line3: +
        f3 = currline_1;
        r3 = currline_2;
    else %line4: quality
        f4 = currline_1;
        r4 = currline_2;
        f_quality = length(find(double(f4)>63))/length(f4);
        r_quality = length(find(double(r4)>63))/length(r4);
        if f_quality > 0.7 && r_quality > 0.7 && ~contains(f2,'N') && ~contains(r2,'N') && length(f4) > 10 && length(r4)>10 
            idx1 = length(f4);
            idx2 = length(r4);
            Length(idx1) = Length(idx1) + 1;
            true_seq = true_seq + 1;
            
            % write into short reads file
            if idx1 >= 60 && idx2 >= 60
                fprintf(fid_f,'%s',f1);
                fprintf(fid_r,'%s',r1);
                fprintf(fid_f,'%s',f2);
                fprintf(fid_r,'%s',r2);
                fprintf(fid_f,'%s',f3);
                fprintf(fid_r,'%s',r3);
                fprintf(fid_f,'%s',f4);
                fprintf(fid_r,'%s',r4);
            % write into long reads file 
            elseif idx1 < 60 && idx2 < 60
                fprintf(fid_f2,'>%s',f1(2:end));
                fprintf(fid_r2,'>%s',r1(2:end));
                fprintf(fid_f2,'%s',f2);
                fprintf(fid_r2,'%s',r2);
            end
        end
    end
    
    currline_1 = fgets(fid1);
    currline_2 = fgets(fid2);
    linecount = linecount + 1;
    if mod(linecount,1000000) == 1
        linecount
        toc
    end
end
fprintf('reads after trim adaptor: %d\n',(linecount-1)/4);

