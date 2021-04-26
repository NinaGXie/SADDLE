clear all

index = 1:192; % Plex Number * 2
plex_num = length(index)/2;
Experiment = 21;
Library = 6;

%get the refernece primer seq
[num1,txt1] = xlsread('E9_96plex_5.xlsx');
fp = txt1(:,2);
rp = txt1(:,4);

ref1 = fp(1:plex_num);
ref2 = rp(1:plex_num);

for i = 1:plex_num
    p{i} = ref1{index(i)};
    p{i+plex_num} = ref2{index(i)};
end

%  Creat hash(dimension x2 in case of collision) for ref sequences
numref = length(p);
hashmod = 786433;
hash1 = zeros(3,hashmod);
errnum = 0;

% save all the primer's first 8 bases
for i = 1:numref
    temphash1 = PrimeNumberHash(lower(p{i}(1:8)), hashmod);
    isEmpty = length(find(hash1(:,temphash1)~=0));
    
    % if there's collision --> save to the next pos in the column
    if isEmpty ~= 0
        fprintf(['ReHash: ', num2str(temphash1), ' is already used by sequence ', num2str(hash1(1,temphash1)), '\n']);
        hash1(isEmpty+1,temphash1) = i;
       % i
    else
        hash1(1,temphash1) = i;
    end
    score_perfect(i) = swalign(p{i}, p{i},'Alphabet','nt');
end

%get the NGS data
fid1 = fopen(sprintf('E%d_Lib%d_DimerCand_R1.fasta',Experiment, Library),'r');
fid2 = fopen(sprintf('E%d_Lib%d_DimerCand_R2.fasta',Experiment, Library),'r');
currline1 = fgets(fid1);
currline2 = fgets(fid2);
linecount = 1;


% dimer & others matrix
dimer(1:length(index), 1:length(index)) = 0;
others(1:length(index), 1:length(index)) = 0;
Imperfect(1:length(index)/2) = 0;
dimer_count = 0;
others_count = 0;
Imperfect_count = 0;
Dimer_Length = zeros(1,150);
Others_Length = zeros(1,150);


tic
%loop all the seq in data
while ischar(currline1)
    if mod(linecount,2) == 0
        dimer_find = 0;
        others_find = 0;
        readF = currline1(1:end-1);
        readR = currline2(1:end-1);
        if length(readF) > 15
            curhash1 = PrimeNumberHash(lower(currline1(1:8)), hashmod);
            curhash2 = PrimeNumberHash(lower(currline2(1:8)), hashmod);
            idx1 = hash1(:,curhash1);
            idx2 = hash1(:,curhash2);
            
            % find primer no collision
            if length(find(idx1 ~= 0)) >= 1 && length(find(idx2 ~= 0)) >= 1
                
                for i = 1:length(find(idx1 ~= 0))
                    if dimer_find == 1 || others_find == 1
                        break                      % % break the outer loop
                    end
                    
                    for j = 1:length(find(idx2 ~= 0))
                        ID1  = idx1(i);
                        ID2  = idx2(j);
                        f_seq = p{ID1};
                        r_seq = p{ID2};
                        if  length(f_seq) < length(currline1) && length(r_seq) < length(currline2)
                            
                            f_score_temp = swalign(f_seq,currline1(1:length(f_seq)),'Alphabet','nt')/score_perfect(ID1);
                            r_score_temp = swalign(r_seq,currline2(1:length(r_seq)),'Alphabet','nt')/score_perfect(ID2);
                            
                            %  choose a threshold value--> however
                            if f_score_temp > 0.6 && r_score_temp > 0.6
                                if length(currline1)-1 < (length(f_seq) + length(r_seq) - 2)
                                    dimer_find = 1;
                                    dimer(ID1, ID2) = dimer(ID1, ID2) + 1;
                                    dimer_count = dimer_count + 1;
                                    Dimer_Length(length(currline1)) = Dimer_Length(length(currline1)) + 1;
                                    break
                                else
                                    others_find = 1;
                                    others_count = others_count + 1;
                                    others(ID1, ID2) = others(ID1, ID2) + 1;
                                    Others_Length(length(currline1)) = Others_Length(length(currline1)) + 1;
                                    break
                                end
                            end
                        end
                    end
                end
            end
        end
    end
    currline1 = fgets(fid1);
    currline2 = fgets(fid2);
    linecount = linecount + 1;
    if mod(linecount,100000) == 1
        linecount
        toc
    end
end

save(sprintf('E%d_L%d_NS.mat', Experiment, Library),'others');
save(sprintf('E%d_L%d_dimer.mat', Experiment, Library),'dimer');

sum(sum(dimer))
sum(sum(others))
save(sprintf('E%d_L%d_Dimer_ReadsLength.mat',Experiment, Library),'Dimer_Length');
save(sprintf('E%d_L%d_NS_ReadsLength.mat',Experiment, Library),'Others_Length');

