%%% Biomass recruiting to size-class (g m-2 d-1)
function rec = sub_rec_larv(X,bio,RE)
  %X: repro rate
  %RE=rfrac: repro efficiency ~ sex ratio * egg survival
  
  rec = RE .* X .* bio;

end