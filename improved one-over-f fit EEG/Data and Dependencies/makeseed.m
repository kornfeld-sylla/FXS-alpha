function [seeds] = makeseed(B_L,sample_size,num_seeds)
%Function to make lists of seeds for non-parametric hierarchical 
%bootstrapping so results can be reproducible across multiple runs of code.
%Input is the number of values needed within a seed to match bootstrap
%repitions and the number of seeds needed based on how many hierarchical
%lists are bootstrapped

%calculate total number of values
B_L_total = B_L*sample_size;
%preallocate
seeds = cell(1,num_seeds);
internal = cell(1,B_L_total);

%generate seeds based on random integers
for mash = 1:num_seeds
    for mush = 1:B_L_total
        internal{mush} = rng(randi(B_L^3));
    end
    seeds{mash} = internal;
end

end

