function [curhash] = PrimeNumberHash(text, hashmod)

rep = fromGCAT(text);

curhash = 0;
for i = 1:length(text)
    %rep(i)
	tempval = mod(rep(i) * 4^(i-1), hashmod);
	curhash = curhash + tempval;
	curhash = mod(curhash, hashmod);
end

curhash = curhash + 1;


