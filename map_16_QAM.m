%ßMÐÐ16QAMÕ{ÖÆÖ®º¯”µ
function output= map_16_QAM(matrix)
[N,NL]=size(matrix);
 N=N/4;
 output=zeros(N,NL);
  for j=1:NL
      for n=1:N
          for ic=1:4
              qam_matrix(ic)=matrix((n-1)*4+ic,j);  
          end
          output(n,j)=qam16(qam_matrix);           
      end
  end   
end