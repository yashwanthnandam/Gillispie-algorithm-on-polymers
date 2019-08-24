clear all;
close all;


%j=different reaction channels
%tsteps=no. of reactions taking place in time

%fileID=fopen('test1_M2.dat','w');
fileID=fopen('test.dat','w');
fileID1=fopen('test1.dat','w');
fileID2=fopen('test2.dat','w');
fileID3=fopen('test3.dat','w');


nreal=10;
bias=0.9;
l_ini1=5;
l_ini=5; %Initially, all polymers are of length l_ini
conc=50; %Concentration of polymers of length l_ini at time=0
conc1=20;
lambda=1;
maxlength=l_ini*conc+l_ini1*conc1; %Maximum possible polymer length
mass=l_ini*conc+l_ini1*conc1; %%-------------------------
D=sqrt(lambda*(lambda+4*mass));
tsteps=100; %No. of time steps
nspecies=maxlength; %No. of species

y=zeros(maxlength,tsteps);
entropy=zeros(tsteps,1);
entropy_th=zeros(tsteps,1);
totconc=zeros(tsteps+1,1);
totconc_th=zeros(tsteps+1,1);


tobsmax=100; %Maximum observation time (give a safe value)
timebins=5000; %number of bins into which the observation time is divided

dt=tobsmax/timebins;

sumconc=zeros(timebins,1);  %Stores coarse-grained data
sumconc_th=zeros(timebins,1); %Stores coarse-grained data
jj=zeros(timebins+1,1);

for n=1:nreal
    
    % Total no. of forward reactions = 1+2+3+...+(maxlength-1)
    % Total no. of backward reactions is equal in number
    
    jmax1 = (maxlength-1)*maxlength/2; %No. of forward reactions
    jmax2=0;  %No. of backward reactions
    jmaxtot=jmax1+jmax2;
    
    x=zeros(nspecies,1);
    x(l_ini)=conc; %initial condition: 20 polymers of length l=5
    x(l_ini1)=conc1;
    %y=zeros(nspecies,1);
    %y=zeros(nspecies,1); %to simulate concentrations via rate eqn
    %y(l_ini)=conc;
    
    c(1:jmax1,1)=2*lambda;
    c(jmax1+1:jmaxtot,1)=2;
    
    %From j=1 to jmax1, all are forward reactions
    %Remaining equations are all backward reactions
    
    % x(time,species) --> matrix x
    % c(channel) --> vector c
    % tsteps=no. of random time steps (or no. of reactions)
    % a(j) = propensity function for jth reaction channel
    % a0 =  sum of a(j)
    % r1,r2 = two uniform random nos. in [0,1]
    % ajsum = sum of a(j) till the sum exceeds r2*a0
    % jmin = mininum value of j for which ajsum exceeds r2*a0
    % tau = random time step
    
    a=zeros(jmaxtot,1); t=zeros(tsteps,1);
    
    
    t(1)=0;
    
    for i=1:tsteps
        
        
        j=0;
        
        for l=2:maxlength
            for k=1:l-1
                j=j+1;
                if k==l-k
                    a(j)=0.5*c(j)*(x(k)-1)*x(l-k);
                else
                    a(j)=0.5*c(j)*x(k)*x(l-k);%Forward reactions (factor of 0.5 to remove double counting)
                end
              %  a(j+jmax1)=0.5*c(j+jmax1)*x(l); %For corresponding backward reactions
            end
        end
        
        
        a0=sum(a);
        
        
        % r1 and r2 are two random nos. in [0,1]
        r1=rand; r2=rand; ajsum=0;
        
        tau = (1/a0)*log(1/r1);
        
        ajsum=0;
        
        for j=1:jmaxtot
            if ajsum<=r2*a0
                ajsum=ajsum+a(j);
                jmin=j;
            end
        end
        
        
        
        t(i+1)=t(i)+tau;
        
        
        j=0;
        
        % Increments in no. of molecules of the species that appear in
        % jmin^{th} reaction
        
        
        
        for l=2:maxlength
            for k=1:l-1
                j=j+1;
                if j==jmin
                    if k>l-k && l-k>1 
                        x(k)=x(k)-1;
                        x(l-k)=x(l-k);
                        if randombias(bias)==1
                            x(k+1)=x(k+1)+1;
                            x(l-k-1)=x(l-k-1)+1;
                        elseif randombias(bias)==0
                            x(k-1)=x(k-1)+1;
                            x(l-k+1)=x(l-k+1)+1;
                        end
                    elseif k>=l-k && l-k>1 
                        x(k)=x(k)-1;
                        x(l-k)=x(l-k);
                        if randombias(bias)==0
                            x(k+1)=x(k+1)+1;
                            x(l-k-1)=x(l-k-1)+1;
                        elseif randombias(bias)==1
                            x(k-1)=x(k-1)+1;
                            x(l-k+1)=x(l-k+1)+1;
                        end   
                    end
 
                end
            end
        end
         
       
        totconc_th(i)=(D/2)*tanh(t(i)*D/2+atanh((2*conc+lambda)/D))-lambda/2;
        totconc(i)=sum(x);                 
            
        nbin=floor(t(i)/dt)+1;
        jj(nbin)=jj(nbin)+1;
        sumconc(nbin)=sumconc(nbin)+totconc(i);
        sumconc_th(nbin)=sumconc_th(nbin)+totconc_th(i);
        
%         for l=1:maxlength
%             if x(l)>0
%                 entropy(i)=entropy(i)-(x(l)/sum(x))*log(x(l)/sum(x));
%             end
%         end
        
        for l=1:maxlength
            y(l,i)=y(l,i)+x(l);
        end
                        
    end
        
end

for l=1:maxlength
   fprintf(fileID1,'%d %.3f\n',l,y(l,tsteps)/nreal);
end


for i=1:tsteps
    for l=1:maxlength
        totsum = sum(y(:,i));
        if y(l,i)>0
            entropy(i)=entropy(i)-(y(l,i)/totsum)*log(y(l,i)/totsum);
        end
    end
    fprintf(fileID2,'%.3f %.3f\n',t(i),entropy(i));
end


subplot(2,2,1) %Plot of concentration vs length

plot(y(:,tsteps)./nreal);
xlabel('length')
ylabel('concentration')
xlim([0 50])

i=1:tsteps;
nbin=1:timebins;

subplot(2,2,2) %Plot of numerical and theoretical values of total concentration

plot(i,totconc(i))
ylabel('totalconcentration')
xlabel('time')
hold on
%plot(i,totconc_th(i),'r')
hold off

subplot(2,2,3)  %Coarse-grained total concentrations (numerical & theoretical)

plot(nbin,sumconc./jj(1:nbin))
ylabel('coarse-grained concentrations')
xlabel('time')
hold on
plot(nbin,sumconc_th./jj(1:nbin),'r--')
hold off
xlim([0 50]);

i=1:tsteps;
subplot(2,2,4)  %Evolution of entropy

plot(t(1:tsteps),entropy);
ylabel('entropy')
xlabel('time')
hold on
%plot(i,entropy_th,'red');
hold off
%axis([0 300 0 10])
%xlim([0 300]);
%axis([0 30 -2 3]);