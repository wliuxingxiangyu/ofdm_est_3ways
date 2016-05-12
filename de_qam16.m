%以下為把16qam調制的信號还原的子函數
function y=de_qam16(x)     
   y=real(x);
   y1=imag(x);
  if         (y>=0)&(y<=2)     y=1;
       elseif (y>2)             y=3;
       elseif (y<-2)            y=-3;
       else                     y=-1;
   end
   if         (y1>=0)&(y1<=2)   y1=1;
       elseif (y1>2)            y1=3;
       elseif (y1<-2)           y1=-3;
       else                     y1=-1;
    end
    x=complex(y,y1);
    if        x==1+j      y=[0 0 0 0];
      elseif  x==1-j      y=[0 0 1 0];
      elseif  x==-1+j     y=[1 0 0 0];
      elseif  x==-1-j     y=[1 0 1 0];
      elseif  x==3+j      y=[0 1 0 0];
      elseif  x==1+3*j    y=[0 0 0 1];
      elseif  x==3-j      y=[0 1 1 0];
      elseif  x==1-3*j    y=[0 0 1 1];
      elseif  x==-1+3*j   y=[1 0 0 1];
      elseif  x==-3+j     y=[1 1 0 0];
      elseif  x==-3-j     y=[1 1 1 0];
      elseif  x==-1-3*j   y=[1 0 1 1];
      elseif  x==3+3*j    y=[0 1 0 1];
      elseif  x==-3+3*j   y=[1 1 0 1];
      elseif  x==-3-3*j   y=[1 1 1 1];
      elseif  x==3-3*j    y=[0 1 1 1];
    end
