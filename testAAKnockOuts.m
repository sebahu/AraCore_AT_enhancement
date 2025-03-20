function [ opt_obj, ko_prots ] = testAAKnockOuts( in_model, ko_prots )

    opt_obj = zeros(length(ko_prots),1);
    for i=1:length(ko_prots)
        opt_obj(i) = testKnockOut(in_model, ko_prots(i));
    end

end