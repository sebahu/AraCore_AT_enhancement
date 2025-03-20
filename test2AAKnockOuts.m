function [ opt_obj, ko_prots ] = test2AAKnockOuts( in_model, ko_prots1 )
    total = (length(ko_prots1)*(length(ko_prots1)-1))/2;
    opt_obj = zeros(total,1);
    ko_prots = strings(total,1);
    total_i=0;
    for i=1:(length(ko_prots1)-1)
        for j=((i+1):length(ko_prots1))
            total_i = total_i+1;
            ko_prots(total_i) = ko_prots1(i)+"___"+ko_prots1(j);
            opt_obj(total_i) = testKnockOut(in_model, [ko_prots1(i), ko_prots1(j)]);
        end
    end

end