function sim = evaluate(G,A)
%FUNCTION EVALUATE
%  Evaluate the clustering results.
%INPUT
%  G: the ground clusters.
%  A: the clustering results.
%OUTPUT
%  sim: the similarity between G and A

[G_num,~] = size(G);
[A_num,~] = size(A);

sum = 0;
SimGiAj = zeros(A_num,1);
for i = 1:G_num
    for j = 1:A_num
        %compute similarity
        GA = intersect(G{i},A{j});
        SimGiAj(j) = 2*length(GA) ./ (length(G{i})+length(A{j}));
    end
    sum = sum + max(SimGiAj);
end

sim = sum / G_num;

