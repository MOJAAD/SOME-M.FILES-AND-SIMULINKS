%=======================================================================================%
%                                  IN THE NAME OF GOD                                   %
%                                 PROJECT OF  LOADFLOW                                  %
%                                BY: MOHAMMAD JAVAD ADEL                                %
%                                      9621010042                                       %
%                                    DATE: 98/11/7                                      %
%=======================================================================================%
clc        % Clear Command Window                                                       %
clear      % Remove items from workspace, freeing up system memory                      %
close all  % closes all figures                                                         %
%=======================================================================================%
%Putting Impedances of lines in Impedance Matrix:
z=[    Inf       0.077+0.308i      Inf            Inf            Inf       0.060+0.242i ;
   0.077+0.308i      Inf       0.083+0.353i       Inf            Inf       0.047+0.198i ;
       Inf       0.083+0.353i      Inf        0.055+0.253i       Inf           Inf      ;
       Inf           Inf       0.055+0.253i       Inf        0.103+0.411i  0.077+0.308i ;
       Inf           Inf           Inf        0.103+0.411i       Inf       0.041+0.176i ;
   0.060+0.242i  0.047+0.198i      Inf        0.077+0.308i   0.041+0.176i       Inf       
  ];
%=======================================================================================%
%Introdusing Admittance Matrix:
Y=zeros(6);
y=zeros(6);
for i=1:6
    for j=1:6
        y(i,j)=1/z(i,j);
        if i~=j
            Y(i,j)=-y(i,j);
        end
    end
end

for i=1:6
    for j=1:6
        if i~=j
        Y(i,i)=Y(i,i)+y(i,j);
        end
    end
end
%=======================================================================================%
                                   %Difinig Values
% The given parameters and initial conditions are
activepower  =[ -0.80   ;   -0.75   ;   -1.25   ;   1.70    ;   -0.85   ;   -1.10   ];
reactivepower=[ -0.30   ;   -0.25   ;   -0.50   ;   0.30    ;   -0.40   ;   -0.50   ];
mv           =[  1.05   ;     1     ;     1     ;   1.02    ;     1     ;     1     ];
v            =[  mv(1)  ;   mv(2)   ;   mv(3)   ;   mv(4)   ;   mv(5)   ;   mv(6)   ];
                                     
acc=0.001;%acceleration constant
stopcondition=1;
                           % The Gauss-Seidel iterations                              
while stopcondition>1e-5
                                    %For P-Q buses                                    
   for i=[2,3,5,6]
     temp1=(activepower(i)-j*reactivepower(i))/conj(v(i));
     temp2=0;
     for k=1:6
       if (i~=k)
           temp2=temp2+Y(i,k)*v(k);
       end
     end
     v_t=(temp1-temp2)/Y(i,i);
     v(i)=v(i)+acc*(v_t-v(i));
   end
                                        % P-V bus
     reactivepower(4)=-0.3;
     for i=1:6
       reactivepower(4)=reactivepower(4)+Y(4,i)*v(i);
     end
     reactivepower(4)=-imag(conj(v(4))*reactivepower(4));
     temp1=(activepower(4)-j*reactivepower(4))/conj(v(4));
     temp2=0;
     for k=[2,3,5,6]
       temp2=temp2+Y(4,k)*v(k);
     end
     v_t=(temp1-temp2)/Y(4,4);
     v(5)=abs(v(4))*v_t/abs(v_t);
                                   % Calculate P and Q
     for i=1:6
       sm=0;
       for k=1:6
           sm=sm+Y(i,k)*v(k);
       end
       s(i)=conj(v(i))*sm;
       activepower(i)=real(s(i));
       reactivepower(i)=imag(s(i));
     end
for i=2:6
     delp(i)=activepower(i)-activepower(i-1);
     delq(i)=reactivepower(i)-reactivepower(i-1);
end
     delpq=[delp(2:6);delq(2:6)];
     stopcondition=min(min(abs(delpq)));
end
%=======================================================================================%
                                   %Display Values
d2r=pi/180;w=100*pi;
disp('FINAL VOLTAGE MAGNITUDESARE:')
abs(v)
disp('FINAL ANGLES IN DEGREE ARE:')
angle(v)/d2r
disp('THE REAL POWERS IN EACH BUS IN MW ARE:')
(real(s)+[0 0 0 0.6 0 0])*100
disp('THE REACTIVE POWERS IN EACH BUS IN MVar ARE:')
(-imag(s)+[0 0 0 0.3 0 0])*100
%=======================================================================================%





