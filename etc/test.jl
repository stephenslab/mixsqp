using LowRankApprox

L = readdlm("sample5000x20.txt",' ');
out = mixSQP(L);
