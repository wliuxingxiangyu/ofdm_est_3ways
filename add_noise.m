%¼ÓÔëÂ•
function output=add_noise(sgma,input)
         [n,nl]=size(input);
         for k=1:nl
            for m=1:n
                noise=normrnd(0,sgma)+normrnd(0,sgma)*sqrt(-1);              
                output(m,k)=input(m,k)+noise;
            end
         end