function output=de_map1(input)
    [m,n]=size(input);
    output=zeros(m,n);
           for k=1:m
               for l=1:n
                if(real(input(k,l))>0)
                      output(k,l)=1;
                else
                     output(k,l)=-1;
               end
           end 
     end