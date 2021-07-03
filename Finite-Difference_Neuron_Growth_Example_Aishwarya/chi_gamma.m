function val = chi_gamma(gamma,p)
val = zeros(size(p));
for i = 1:length(p)
    if(p(i)~=0)
        val(i) = tanh(gamma*(p(i)))/(p(i));
    else
        val(i) = gamma;
    end
end
end