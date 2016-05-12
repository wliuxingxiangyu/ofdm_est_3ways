%ȥѭ�h
function output=del_cp(input,cp_length)   
        [m,n]=size(input);
         output=zeros(m-cp_length,n);
         for j=1:n
             output(1:(m-cp_length),j)=input((cp_length+1:m),j);
         end