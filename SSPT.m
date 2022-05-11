%Phase transformation in Ti64
%@author: Xinyu Yang，March 2020
%Mechanical Engineering， NUI Galway
%PhD Project funded by I-Form, the SFI Research Centre for Advanced Manufacturing
%Grant number 16/RC/3872
close all
clear all
clc

%Material Props
%Load temperature history from excel file
t=xlsread('Temperature history.xls','a:a');
T=xlsread('Temperature history.xls','b:b');
K=T+273.15;

%Polynomial description of equilibrium alpha and belta phase fraction   
for i=1:numel(T)
if T(i)>=994
    xeqalpha(i)=0;
else
    if 650<=T(i)
xeqalpha(i)=(-31188.514)*(T(i)/1000).^8+170526.26*(T(i)/1000).^7+(-388991.69)*(T(i)/1000).^6+471927.45*(T(i)/1000).^5+(-315178.49)*(T(i)/1000).^4+99079.891*(T(i)/1000).^3+1667.1991*(T(i)/1000).^2+(-9726.8403)*(T(i)/1000).^1+1884.7280*(T(i)/1000).^0;
    else
        if T(i)<650
       xeqalpha(i)=0.93;    %need to adjust it to different AM process, i.e. fully alpha means 1
        end
    end
end
xeqbelta(i)=1-xeqalpha(i);
xeqalpha=xeqalpha';
xeqbelta=xeqbelta';
end

%Define kinetic parameters k, p from TTT curves
xi=0.01; %initial measured fraction
xf=0.5;  %final measured fraction
c=log(1-xf)/log(1-xi); %Calculate c

%Load grain boundary alpha TTT value;
%Calculate initial time
ti1_gb=xlsread('TTT value.xls','a:a');
Ti1_gb=xlsread('TTT value.xls','b:b');
ti_gb=interp1(Ti1_gb,ti1_gb,T,'PCHIP');
%Calculate final time
tf1_gb=xlsread('TTT value.xls','c:c');
Tf1_gb=xlsread('TTT value.xls','d:d');
tf_gb=interp1(Tf1_gb,tf1_gb,T,'PCHIP');

%Load widmanstatten alpha TTT value;
%Calculate initial time
ti1_w=xlsread('TTT value.xls','e:e');
Ti1_w=xlsread('TTT value.xls','f:f');
ti_w=interp1(Ti1_w,ti1_w,T,'PCHIP');
%Calculate final time
tf1_w=xlsread('TTT value.xls','g:g');
Tf1_w=xlsread('TTT value.xls','h:h');
tf_w=interp1(Tf1_w,tf1_w,T,'PCHIP');

%k_gb, p_gb of grain boundary alpha
p_gb=zeros(numel(T),1);
k_gb=zeros(numel(T),1);
for i=1:numel(T)
    if T(i)>986
        p_gb(i)=0;
        k_gb(i)=0;
    else
        if 378<=T(i)
        p_gb(i)=log(c)/log(tf_gb(i)/ti_gb(i));
        k_gb(i)=-log(1-xf)*tf_gb(i).^(-p_gb(i));
    else
        p_gb(i)=0;
        k_gb(i)=0;
        end
    end
end

%k_w, p_w of widmanstatten alpha
p_w=zeros(numel(T),1);
k_w=zeros(numel(T),1);
for i=1:numel(T)
    if T(i)>986
        p_w(i)=0;
        k_w(i)=0;
    else
        if 259<=T(i)
        p_w(i)=log(c)/log(tf_w(i)/ti_w(i));
        k_w(i)=-log(1-xf)*tf_w(i).^(-p_w(i));
    else
        p_w(i)=0;
        k_w(i)=0;
        end
    end
end

%k_m, p_m of martensite alpha
%load reference k1_m, p1_m
k1_m=xlsread('TTT value.xls','i:i');
p1_m=xlsread('TTT value.xls','j:j');
T1_m=xlsread('TTT value.xls','k:k');
k_m=interp1(T1_m,k1_m,T,'PCHIP');
p_m=interp1(T1_m,p1_m,T,'PCHIP');
%equilibrium martensite fraction
%load reference xeqalpha1_m
xeqalpha1_m=xlsread('TTT value.xls','l:l');
xeqalpha_m=interp1(T1_m,xeqalpha1_m,T,'PCHIP');        
        for i=1:numel(T)
            if T(i)>800
                p_m(i)=0;
                k_m(i)=0;
                xeqalpha_m(i)=0;
            else
                if T(i)<400
                    p_m(i)=0;
                    k_m(i)=0;
                    xeqalpha_m(i)=1;
                end
            end
        end

                
%Initialise varibles

%Define varibles
xbelta=zeros(numel(T),1);
xalpha_gb=zeros(numel(T),1);
xalpha_w=zeros(numel(T),1);
xalpha_m=zeros(numel(T),1);
xalpha_gbw=zeros(numel(T),1);
te_gb=zeros(numel(T),1);
te_w=zeros(numel(T),1);

xalpha_gbw(1)=0.95;  %Fraction of alpha_gb and alpha_w
xalpha_gb(1)=0.03;   %Fraction of alpha_gb
xalpha_w(1)=xalpha_gbw(1)-xalpha_gb(1);  %Fraction of alpha_w
xbelta(1)=0.05;  %Fraction of belta
xalpha_m(1)=0;  %Fraction of alpha_m


%Loop through values of temperature
for i=2:numel(T)
    %Temperature increment
    deltaT=abs(T(i)-T(i-1));
    %Time increment
    deltat=abs(t(i)-t(i-1));
    
    %Determine if heating or cooling
    if T(i)>=T(i-1)
        m=0; %Heating
    else
        m=1; %Cooling
    end
    
    if m==0 %Heating
        %Compare the T with martensite dissolution start temperature
        if T(i)<400
        xbelta(i)=xbelta(i-1);
        xalpha_gb(i)=xalpha_gb(i-1);
        xalpha_w(i)=xalpha_w(i-1);        
        xalpha_gbw(i)=xalpha_gb(i)+xalpha_w(i);
        xalpha_m(i)=xalpha_m(i-1);
        else 
            if T(i)<=700
                if  xalpha_m(i-1)==0                %confirm the martensite fraction
                    xbelta(i)=xbelta(i-1);
                    xalpha_gb(i)=xalpha_gb(i-1);
                    xalpha_w(i)=xalpha_w(i-1);        
                    xalpha_gbw(i)=xalpha_gb(i)+xalpha_w(i);
                    xalpha_m(i)=xalpha_m(i-1);
                else    
                %Dissolution of martensite alpha
                %Equivalent time of alpha_m         
                if xalpha_m(i-1)<=xeqalpha_m(i)   %modification1 compare with eqalpha_m
                te_m(i)=0;
                xbelta(i)=xbelta(i-1);
                xalpha_gb(i)=xalpha_gb(i-1);
                xalpha_w(i)=xalpha_w(i-1);        
                xalpha_gbw(i)=xalpha_gb(i)+xalpha_w(i);
                xalpha_m(i)=xalpha_m(i-1);
                else
                te_m(i)=((-log(1-((xalpha_m(i-1)-xeqalpha_m(i))/(xbelta(i-1)+xalpha_m(i-1)-xeqalpha_m(i)))))/(k_m(i))).^(1/p_m(i));
                %Dissolution of martensite alpha
                xalpha_m(i)=xeqalpha_m(i)-exp(-k_m(i)*(te_m(i)+deltat).^p_m(i))*(xalpha_m(i-1)+xbelta(i)-xeqalpha_m(i));
				if xalpha_m(i)<0                    %modification4 confirm retained xalpha_m
				xalpha_m(i)=xalpha_m(i-1);
				else
				xalpha_m(i)=xalpha_m(i);
				end
                %Fraction of alpha_gb at step n
                xalpha_gb(i)=xalpha_gb(i-1);
                %Fraction of alpha_w at step n
                xalpha_w(i)=xalpha_w(i-1)+(xalpha_m(i-1)-xalpha_m(i))*(1-xeqbelta(i));
                %Fraction of alpha_gb and alpha_w at step n
                xalpha_gbw(i)=xalpha_gb(i)+xalpha_w(i);
                %Fraction of belta at step n
                xbelta(i)=1-xalpha_gbw(i)-xalpha_m(i);
                end
                end
            else
                if T(i)<994
                %Determine martensite alpha
                xalpha_m(i)=xalpha_m(i-1);
                %Dissolution of alpha to belta
                %Dissolution function
                f_d(i)=2.2*10^(-31)*K(i).^9.89;             
                %Equivalent time of alpha dissolution
                te_d(i)=(xbelta(i-1)/(xeqbelta(i)*f_d(i))).^2;
                if (deltat+te_d(i))<f_d(i).^(-2)
                    if xalpha_w(i-1)>0
                        %Fraction of alpha_w at step n
                        if xbelta(i-1)>xeqbelta(i)
                            xbelta(i)=xbelta(i-1);
                            xalpha_gbw(i)=1-xbelta(i);
                            xalpha_gb(i)=xalpha_gb(i-1);
                            xalpha_w(i)=xalpha_gbw(i)-xalpha_gb(i);
                        else
                        xbelta(i)=xeqbelta(i)*f_d(i)*(deltat+te_d(i)).^(1/2);
                        xalpha_gbw(i)=1-xbelta(i);
                        xalpha_gb(i)=xalpha_gb(i-1);
                        xalpha_w(i)=xalpha_gbw(i)-xalpha_gb(i);
                        end
                    else
                        if xbelta(i-1)>xeqbelta(i)
                            xbelta(i)=xbelta(i-1);
                            xalpha_gbw(i)=1-xbelta(i);
                            xalpha_gb(i)=xalpha_gb(i-1);
                            xalpha_w(i)=xalpha_gbw(i)-xalpha_gb(i);
                        else
                        xbelta(i)=xeqbelta(i)*f_d(i)*(deltat+te_d(i)).^(1/2);
                        xalpha_gbw(i)=1-xbelta(i);
                        xalpha_w(i)=xalpha_w(i-1);
                        xalpha_gb(i)=xalpha_gbw(i)-xalpha_w(i);
                        end
                    end
                else
                    if xalpha_w(i-1)>0
                        %Fraction of alpha_w at step n
                        if xbelta(i-1)>xeqbelta(i)
                            xbelta(i)=xbelta(i-1);
                            xalpha_gbw(i)=1-xbelta(i);
                            xalpha_gb(i)=xalpha_gb(i-1);
                            xalpha_w(i)=xalpha_gbw(i)-xalpha_gb(i);
                        else
                        xbelta(i)=xeqbelta(i);
                        xalpha_gbw(i)=1-xbelta(i);
                        xalpha_gb(i)=xalpha_gb(i-1);
                        xalpha_w(i)=xalpha_gbw(i)-xalpha_gb(i);
                        end
                    else
                        if xbelta(i-1)>xeqbelta(i)
                            xbelta(i)=xbelta(i-1);
                            xalpha_gbw(i)=1-xbelta(i);
                            xalpha_gb(i)=xalpha_gb(i-1);
                            xalpha_w(i)=xalpha_gbw(i)-xalpha_gb(i);
                        else
                        xbelta(i)=xeqbelta(i);
                        xalpha_gbw(i)=1-xbelta(i);
                        xalpha_w(i)=xalpha_w(i-1);
                        xalpha_gb(i)=xalpha_gbw(i)-xalpha_w(i);
                        end
                    end
                end
                else
                    xalpha_gb(i)=0;
                    xalpha_w(i)=0;
                    xalpha_gbw(i)=xalpha_gb(i)+xalpha_w(i);
                    xbelta(i)=1-xalpha_gbw(i)-xalpha_m(i);
                end
            end
        end
        
    else %Cooling
    %Compare the T with T-Belta trans
    if T(i)>=986
        xbelta(i)=1;
        xalpha_gb(i)=0;
        xalpha_w(i)=0;        
        xalpha_gbw(i)=xalpha_gb(i)+xalpha_w(i);
        xalpha_m(i)=0;
    else
        %Compare the T with T-martensite start point
        if T(i)>=575
            %Formation of grain boundary and widmanstatten alpha   %modification2 retained alpha should compare with eqalpha in fear of incomplete transformation due to the short time
            if xalpha_gbw(i-1)>=xeqalpha(i)
                te_gb(i)=0;
                te_w(i)=0;
                xbelta(i)=xbelta(i-1);
                xalpha_gb(i)=xalpha_gb(i-1);
                xalpha_w(i)=xalpha_w(i-1);        
                xalpha_gbw(i)=xalpha_gb(i)+xalpha_w(i);
                xalpha_m(i)=xalpha_m(i-1);
            else
            %Equivalent time of alpha_gb
            te_gb(i)=((-log(1-((xalpha_gbw(i-1)/xeqalpha(i))/(xbelta(i-1)+xalpha_gbw(i-1)))))/(k_gb(i))).^(1/p_gb(i));
            %Equivalent time of alpha_w
            te_w(i)=((-log(1-((xalpha_gbw(i-1)/xeqalpha(i))/(xbelta(i-1)+xalpha_gbw(i-1)))))/(k_w(i))).^(1/p_w(i));
            %Fraction of alpha_gb at step n
            xalpha_gb(i)=(1-exp(-k_gb(i)*(te_gb(i)+deltat).^p_gb(i)))*(xbelta(i-1)+xalpha_gbw(i-1))*(xeqalpha(i))-xalpha_w(i-1);
            %Fraction of alpha_w at step n
            xalpha_w(i)=(1-exp(-k_w(i)*(te_w(i)+deltat).^p_w(i)))*(xbelta(i-1)+xalpha_gbw(i-1))*(xeqalpha(i))-xalpha_gb(i-1);
            %Fraction of alpha_gb and alpha_w at step n
            xalpha_gbw(i)=xalpha_gb(i)+xalpha_w(i);
            %Formation of martensite alpha
            xalpha_m(i)=xalpha_m(i-1);
            %Fraction of belta at step n
            xbelta(i)=1-xalpha_gbw(i)-xalpha_m(i);
            end
        else
            if T(i)>=378
            %Formation of grain boundary and widmanstatten alpha
            %Equivalent time of alpha_gb
            te_gb(i)=((-log(1-((xalpha_gbw(i-1)/xeqalpha(i))/(xbelta(i-1)+xalpha_gbw(i-1)))))/(k_gb(i))).^(1/p_gb(i));
            %Equivalent time of alpha_w
            te_w(i)=((-log(1-((xalpha_gbw(i-1)/xeqalpha(i))/(xbelta(i-1)+xalpha_gbw(i-1)))))/(k_w(i))).^(1/p_w(i));
            %Fraction of alpha_gb at step n
            xalpha_gb(i)=(1-exp(-k_gb(i)*(te_gb(i)+deltat).^p_gb(i)))*(xbelta(i-1)+xalpha_gbw(i-1))*(xeqalpha(i))-xalpha_w(i-1);
            %Fraction of alpha_w at step n
            xalpha_w(i)=(1-exp(-k_w(i)*(te_w(i)+deltat).^p_w(i)))*(xbelta(i-1)+xalpha_gbw(i-1))*(xeqalpha(i))-xalpha_gb(i-1);
            %Fraction of alpha_gb and alpha_w at step n
            xalpha_gbw(i)=xalpha_gb(i)+xalpha_w(i);
            %Formation of martensite alpha
            if deltaT/deltat>410 %Compare the cooling rate with 410 degree
                xalpha_m(i)=(1-exp(-0.005*(575-T(i))))*(xbelta(i-1)+xalpha_m(i-1));
                else
            if 20<deltaT/deltat<=410
                xalpha_m(i)=(1-exp(-0.005*(575-T(i))))*(xbelta(i-1)+xalpha_m(i-1)-xeqbelta(i));
            else
                xalpha_m(i)=xalpha_m(i-1);
            end
            end
            %Fraction of belta at step n
            xbelta(i)=1-xalpha_gbw(i)-xalpha_m(i);
            else
                if T(i)>=259
            %Formation of grain boundary and widmanstatten alpha
            %Fraction of alpha_gb at step n
            xalpha_gb(i)=xalpha_gb(i-1);
            %Equivalent time of alpha_w
            if ((xalpha_gbw(i-1)/xeqalpha(i))/(xbelta(i-1)+xalpha_gbw(i-1)))>1   %modification3 retained alpha_w should compare with eqalpha in fear of incomplete transformation due to the short time
                xalpha_w(i)=xalpha_w(i-1);
                xalpha_gbw(i)=xalpha_gb(i)+xalpha_w(i);
            %Formation of martensite alpha
            if deltaT/deltat>410 %Compare the cooling rate with 410 degree
                xalpha_m(i)=(1-exp(-0.005*(575-T(i))))*(xbelta(i-1)+xalpha_m(i-1));
                else
            if 20<deltaT/deltat<=410
                xalpha_m(i)=(1-exp(-0.005*(575-T(i))))*(xbelta(i-1)+xalpha_m(i-1)-xeqbelta(i));
            else
                xalpha_m(i)=xalpha_m(i-1);
            end
            end
            else
			%Formation of alpha_w
            te_w(i)=((-log(1-((xalpha_gbw(i-1)/xeqalpha(i))/(xbelta(i-1)+xalpha_gbw(i-1)))))/(k_w(i))).^(1/p_w(i));
            %Fraction of alpha_w at step n
            xalpha_w(i)=(1-exp(-k_w(i)*(te_w(i)+deltat).^p_w(i)))*(xbelta(i-1)+xalpha_gbw(i-1))*(xeqalpha(i))-xalpha_gb(i-1);
            %Fraction of alpha_gb and alpha_w at step n
            xalpha_gbw(i)=xalpha_gb(i)+xalpha_w(i);
            %Formation of martensite alpha
            if deltaT/deltat>410 %Compare the cooling rate with 410 degree
                xalpha_m(i)=(1-exp(-0.005*(575-T(i))))*(xbelta(i-1)+xalpha_m(i-1));
                else
            if 20<deltaT/deltat<=410
                xalpha_m(i)=(1-exp(-0.005*(575-T(i))))*(xbelta(i-1)+xalpha_m(i-1)-xeqbelta(i));
            else
                xalpha_m(i)=xalpha_m(i-1);
            end
            end
            end
            %Fraction of belta at step n
            xbelta(i)=1-xalpha_gbw(i)-xalpha_m(i);
                else
                    if T(i)<259
            %Formation of grain boundary and widmanstatten alpha
            %Fraction of alpha_gb at step n
            xalpha_gb(i)=xalpha_gb(i-1);
            %Fraction of alpha_w at step n
            xalpha_w(i)=xalpha_w(i-1);
            %Fraction of alpha_gb and alpha_w at step n
            xalpha_gbw(i)=xalpha_gb(i)+xalpha_w(i);
             %Formation of martensite alpha
            if deltaT/deltat>410 %Compare the cooling rate with 410 degree
                xalpha_m(i)=(1-exp(-0.005*(575-T(i))))*(xbelta(i-1)+xalpha_m(i-1));
                else
            if 20<deltaT/deltat<=410
                xalpha_m(i)=(1-exp(-0.005*(575-T(i))))*(xbelta(i-1)+xalpha_m(i-1)-xeqbelta(i));
            else
                xalpha_m(i)=xalpha_m(i-1);
            end
            end
            %Fraction of belta at step n
            xbelta(i)=1-xalpha_gbw(i)-xalpha_m(i);
                end
            end
        end
    end
    end
    end
end

te_m=te_m';
f_d=f_d';
te_d=te_d';
xeqalpha=xeqalpha';


%lath width prediction
xalpha=xalpha_gb+xalpha_w+xalpha_m;
w(1)=1;
%Loop through values of temperature
for i=2:numel(K)
weq(i)=1.42*exp(-294/K(i));
if xalpha(i)==0
    w(i)=0;
else
    w(i)=(w(i-1)*xalpha(i-1)+weq(i)*(xalpha(i)-xalpha(i-1)))/xalpha(i);
end
end
w=w';

plot(t,xalpha,'-black')


