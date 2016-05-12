%ßMÐÐ4QAMÕ{ÖÆÖ®º¯”µ
function output= map_4_QAM(matrix)
[N,NL]=size(matrix);
 N=N/2;
 output=zeros(N,NL);
  for j=1:NL
      for n=1:N
          for ic=1:2
              qam_matrix(ic)=matrix((n-1)*2+ic,j);  
          end
          output(n,j)=qam4(qam_matrix);           
      end
  end   
end