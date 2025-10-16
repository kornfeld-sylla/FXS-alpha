function [re_ordered_corrected] = FDR_correct(uncorrected_values)
%Perform FDR correction on multiple p-values calculated from same figure
%   list of uncorrected p-values is the input argument
%   the function will rank order and correct scaled to the number of comparisons
%   then, it re-orders the values to match the input order

%sort uncorrected p-values in ascending order
[rank_order,indicies] = sort(uncorrected_values);

%preallocate
correction = ones(1,length(uncorrected_values));

%correction values are scaled by the number of comparisons, with the
%smallest uncorrected p-value multipled by the total number of comparisons
%and the largest uncorrected p-value multiplied by the total number of
%comparisons/total number of comparisons (i.e., 1).
for ie = 1:length(correction)
    correction(ie) = length(uncorrected_values)/ie;
end

%apply correction
rank_order_corrected = rank_order.*correction;

%preallocate
re_ordered_corrected = zeros(1,length(rank_order_corrected));

%sort back to input order
for ou = 1:length(indicies)
    re_ordered_corrected(indicies(ou)) = rank_order_corrected(ou);
end
end