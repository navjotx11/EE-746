

function Iapp = IappCalculation(we,S,t,N)
    Io = 1e-12;
    tau = 15e-3;
    tau_s = tau/4;

    Iapp = zeros(1,length(t));
    for ii=1:N
        for i=1:length(t)
            for j=1:length(S{ii})
                if t(i)<S{ii}(j), break; end
                Iapp(i) = Iapp(i) + Io*we(ii)*(exp(-(t(i)-S{ii}(j))/tau)-exp(-(t(i)-S{ii}(j))/tau_s));
            end
        end
    end

end