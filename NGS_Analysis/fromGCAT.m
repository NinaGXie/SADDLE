function [output]=fromGCAT(inputseq)
for i = 1:length(inputseq)
    if ((inputseq(i) == 'A')|(inputseq(i) == 'a'))
        output(i) = 1;
    elseif ((inputseq(i) == 'T')|(inputseq(i) == 't')|(inputseq(i) == 'U')|(inputseq(i) == 'u'))
        output(i) = 2;
    elseif ((inputseq(i) == 'C')|(inputseq(i) == 'c'))
        output(i) = 3;
    elseif ((inputseq(i) == 'G')|(inputseq(i) == 'g'))
        output(i) = 4;
	elseif (inputseq(i) == ' ')
		output(i) = 0;
    else
        inputseq(i)
        fprintf('Unexpected nucleotide!\n');
	end
	if (upper(inputseq(i)) == inputseq(i))
		output(i) = output(i) + 0.1;
	end
end