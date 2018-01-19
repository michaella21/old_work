%Basic Parameters
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
E_no=100;%number of total engines
rep=2;
fr=0.95;%required fill rate
sigma=0.0005;
theta=sigma*2.2;
F=0.07; %threshold
cycle=floor(1.1*F/theta);
decision=cycle+1;
L=ceil(cycle/3);%lead time
delta=cycle/7;
obs_prd=2*cycle;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
ED=zeros(rep,1);%mean
VD=zeros(rep,1);%variance
ex_signal=zeros(E_no,rep);%for example only
ex_visit=zeros(E_no,rep);%for example only
result_matrix=zeros(2,obs_prd,rep);%matrix of mean and variance of all replications
result_size={};
bs_level=zeros(rep,obs_prd);
bs_level2=zeros(rep,obs_prd);
phic={};upsilonc={};psic={};
kmatrix=zeros(rep,obs_prd);
k2matrix=zeros(rep,obs_prd-5);
for_table=zeros(5,rep);
for_table2={};
bsmax=0;
bsmin=1000;
frate=zeros(1,3,rep);
frate_mo=zeros(1,floor(obs_prd/4),rep);
frate_qr=zeros(1,floor(obs_prd/12),rep);
oh_inv=zeros(1,3,rep+1);
oh_inv_mo=zeros(1,floor(obs_prd/4),rep);
oh_inv_qr=zeros(1,floor(obs_prd/12),rep);
bs_mean=zeros(1,3,rep+1);
bs_mean_mo=zeros(1,floor(obs_prd/4),rep);
bs_mean_qr=zeros(1,floor(obs_prd/12),rep);
ini_oh=200;
inv=zeros(5,obs_prd,rep);
result_mo=zeros(3,floor(obs_prd/4));
result_qr=zeros(3,floor(obs_prd/12));
result_lt=zeros(3,3);

for i=1:1:rep
    inv(1,1,i)=ini_oh; %initial start on-hand inventory-set by trial and error
    inv(1,2,i)=ini_oh;
end
bo=zeros(2,obs_prd,rep);%backorders: first row how much occured, second row after filled how much left. 
order=zeros(1,obs_prd,rep);%how much on-order


%se=rng;

for ii=1:1:rep
    t=cputime;
    %rng(se);
    %rng
    engine_set=randi([1,cycle],E_no,1,ii);
    next_visit=zeros(E_no,1,ii);
    for i=1:1:E_no
        next_visit(i,1,ii)=engine_set(i,1,ii)+cycle;
    end
   
    %generate standard brownian error until the beggining of period k(=time
    %k+1)+L(lead time) to see the fluctuation of base stock level
    w=zeros(E_no,decision+obs_prd,ii);

    z=zeros(E_no,decision+obs_prd,ii);
    
    for i=1:1:E_no
        z(i,1,ii)=theta*1;
        for t=2:1:decision+obs_prd
        
            w(i,t,ii)=randn+w(i,t-1,ii);
            z(i,t,ii)=theta*t+sigma*w(i,t,ii);
        %while z(i,t)<0            
         %   w(i,t)=randn+w(i,t-1);
         %    z(i,t)=theta*t+sigma*w(i,t);            
        %end
            
        end
    end
    for i=1:1:E_no
        for j=1:1:decision+obs_prd
            if z(i,j,ii)<0
                z(i,j,ii)
            end
        end
    end
    %want to observe the value at time k+1(=cycle+1)
    ti=zeros(obs_prd,1,ii);
    ti(1,1,ii)=6;
    for i=1:1:obs_prd-6
        ti(i+1,1,ii)=ti(i,1,ii)+1;
    end
    
    base_signal=zeros(E_no,obs_prd,ii);
    for i=1:1:E_no
        for j=1:1:obs_prd
             base_signal(i,j,ii)=z(i,decision-engine_set(i,1)+j-1,ii);
             
        end
    end
    
    %categorize the engines into 3 groups: planned/both/unplanned only
    %while the time moment to estimate the future demand, engine groups
    %should be changed based on the same rule.
    ct=0;
    
    for j=1:1:obs_prd
        pr=0;
        upr=0;
        pr_wo_re=0;%planned maintenance without replacement
       %per period (each j)
        %update the next visit and reset the signal after the visit 
        %whether there was an either 'planned' or 'unplanned' maintenance
        %during the previous period, those engines' 'next visit' dates are
        %updated. 
       
          for i=1:1:E_no
              if base_signal(i,j,ii)>=F
                  if next_visit(i,1,ii)==decision+j-1 
                      pr=pr+1;
                      next_visit(i,1,ii)=decision+j-1+cycle;
                      base_signal(i,j,ii)=theta*1; %singal is reset from the current period.
                      w(i,j,ii)=0;
                      for s=j+1:1:obs_prd
                          base_signal(i,s,ii)=0;%first reset the rest of the base_signal for the future
                          w(i,s,ii)=randn+w(i,s-1,ii);
                          z(i,s,ii)=theta*(s-j)+sigma*w(i,s,ii);
                          base_signal(i,s,ii)=z(i,s,ii);%regeneratethe signal starting from j=theta (if signal>F in (j-1), took donw, replaced)
                      end
                  elseif next_visit(i,1,ii)>decision+j-1+floor(delta)
                      upr=upr+1;
                      next_visit(i,1,ii)=decision+j-1+cycle;
                      base_signal(i,j,ii)=theta*1;
                      w(i,j,ii)=0;
                      for s=j+1:1:obs_prd
                          base_signal(i,s,ii)=0;%first reset the rest of the base_signal for the future
                          w(i,s,ii)=randn+w(i,s-1,ii);
                          z(i,s,ii)=theta*(s-j)+sigma*w(i,s,ii);
                          if z(i,s,ii)>0
                              z(i,s,ii);
                          end
                          base_signal(i,s,ii)=z(i,s,ii);%regeneratethe signal starting from j=theta (if signal>F in (j-1), took donw, replaced)
                      end
                  end
              else
                  if next_visit(i,1,ii)==decision+j-1
                     pr_wo_re=pr_wo_re+1;
                     next_visit(i,1,ii)=decision+j-1+cycle;
                  end
              end
          end
          
       
              
            
               
      
          
  
        engine_group1={};
        engine_group2={};
        engine_group3={};
        c1=1;c2=1;c3=1;

        for i=1:1:E_no
            if next_visit(i,1,ii)>= decision+j-1 && next_visit(i,1,ii)<= floor(decision+j-1+delta)
                engine_group1{c1,1}=i;
                c1=c1+1;
            elseif next_visit(i,1,ii)>= ceil(decision+delta+j-1) && next_visit(i,1,ii)<= decision+j-1+L
                    engine_group2{c2,1}=i;
                    c2=c2+1;
            else
                 engine_group3{c3,1}=i;
                c3=c3+1;
            end
        end


    %compute the part probabilities for each engine in engine_groups
    remaining_life=@(z,tau,F,theta,sigma)(1-(normcdf(((F-z)-(theta*tau))/(sigma*sqrt(tau)),0,1))+(exp(2*theta*(F-z)/(sigma^2))*normcdf((-F+z-(theta*tau))/(sigma*sqrt(tau)),0,1)));
    %just to identify where NaN problems occur
    %remaining_life1=@(z,tau,F,theta,sigma)(2*theta*(F-z)/(sigma^2));
    %remaining_life2=@(z,tau,F,theta,sigma)((-F+z-(theta*tau))/(sigma*sqrt(tau)));
    
    egroup1=cell2mat(engine_group1);
    egroup2=cell2mat(engine_group2);
    egroup3=cell2mat(engine_group3);

    
    for i=1:1:length(egroup1)
         if base_signal(egroup1(i,1),j,ii)>= F
            phic{i,1}=1;
          else
            phic{i,1}=remaining_life(base_signal(egroup1(i,1),j,ii),L+1,F,theta,sigma)-remaining_life(base_signal(egroup1(i,1),j,ii),0,F,theta,sigma);
        end
    end

   
    for i=1:1:length(egroup2)
        if base_signal(egroup2(i,1),j,ii)>=F
            upsilonc{i,1}=1;
        else
            upsilonc{i,1}=remaining_life(base_signal(egroup2(i,1),j,ii),L+1,F,theta,sigma)-remaining_life(base_signal(egroup2(i,1),j,ii),0,F,theta,sigma);
        end
    end
  
   
    for i=1:1:length(egroup3)
        if base_signal(egroup3(i,1),j,ii)>=F
            psic{i,1}=1;
        else
           psic{i,1}=remaining_life(base_signal(egroup3(i,1),j,ii),min(L+1,next_visit(egroup3(i,1),1,ii)-ceil(delta)-decision-j+1),F,theta,sigma)-remaining_life(base_signal(egroup3(i,1),j,ii),0,F,theta,sigma);
           if psic{i,1}<0
              base_signal(egroup3(i,1),j,ii)
              min(L+1,next_visit(egroup3(i,1),1,ii)-ceil(delta)-decision-j+1)
              break;
           end
              
        end
    end
 %   
 %for b=1:1:length(egroup1)
  %   ex_signal(b,ii)=base_signal(egroup1(b,1),1);
   %  ex_visit(b,ii)=next_visit(egroup1(b,1),1);
 %end
 %for b=1:1:length(egroup2)
  %   ex_signal(b+length(egroup1),ii)=base_signal(egroup2(b,1),1);
  %  ex_visit(b+length(egroup1),ii)=next_visit(egroup2(b,1),1);
 %end
 %for b=1:1:length(egroup3)
  %   ex_signal(b+length(egroup1)+length(egroup2),ii)=base_signal(egroup3(b,1),1);
  %  ex_visit(b+length(egroup1)+length(egroup2),ii)=next_visit(egroup3(b,1),1);
% end
   
  phi=cell2mat(phic);
  upsilon=cell2mat(upsilonc);
  psi=cell2mat(psic);
  
  mu_L=(sum(phi)+sum(upsilon)+sum(psi));
  var=sum(phi.*(1-phi))+sum(upsilon.*(1-upsilon))+sum(psi.*(1-psi));
  std=sqrt(var);
  result_matrix(1,j,ii)=mu_L;
  result_matrix(2,j,ii)=std;
  if mu_L==0 || std ==0 
      j,ii
      break;
  end
  mu=mu_L/(L+1);
    for k=10:-1:-30
          if mu*(1-fr)-std*(normpdf(k,0,1)-k*(1-normcdf(k,0,1)))<0
              n=k;
              break;
          end
    end
    for kk=n:0.1:n+1
          if mu*(1-fr)-std*(normpdf(kk,0,1)-kk*(1-normcdf(kk,0,1)))>0
             nn=kk;
            break;
          end
    end
    for kkk=nn:-0.01:nn-0.1
          if mu*(1-fr)-std*(normpdf(kkk,0,1)-kkk*(1-normcdf(kkk,0,1)))<0
           
             bs_level(ii,j)=mu_L+kkk*std;
             kmatrix(ii,j)=kkk;
            break;
          end
    end
    %for k4=nnn:0.001:nnn+0.01
     %      if mu*(1-fr)-std*(normpdf(k4,0,1)-k4*(1-normcdf(k4,0,1)))>0
      %       bs_level(ii,j)=mu_L+k4*std;%bs_level newly calculated in period j
       %      kmatrix(ii,j)=k4;
        %    
           % break;
      
         %  end
    %end
      
  
    if  kkk==0;
        break;
    end
 
     
     %if j>5
        % bs_level2(ii,j)=bs_level(ii,j);
       %  k2matrix(ii,j-5)=mu/std;
     %end
    
    
  %   if j==6
   %      result_size{ii,1}=length(egroup1);
    %     result_size{ii,2}=length(egroup2);
     %    result_size{ii,3}=length(egroup3);
      %   for_table(1,ii)=mu_L;
      %   for_table(2,ii)=mu;
      %   for_table(3,ii)=std;
      %   for_table(4,ii)=mu_L+k4*std;
      %   for_table(5,ii)=k4;
         
          
     %end
       %periodical demand from plannde/unplanned maintenance
         inv(4,j,ii)=pr;
         inv(5,j,ii)=upr;
         onorder=0;
         onorder2=0;
         
         if j>1 && j<L+1 %before we get the first delivery
             inv(2,j,ii)=inv(1,j,ii); %starting OH inventory, there's no delivery during this period
             for q=1:1:j-1 %first fill the backordered items
                 onorder=onorder+order(1,q,ii);%every order upto L is on-order
             end
             inv(3,j,ii)=inv(2,j,ii)+onorder;%inventory position
             order(1,j,ii)=max(0,round(bs_level(ii,j-1))-inv(3,j,ii)); %since we order BFFORE filling the demand
             
             if inv(2,j,ii)>=inv(4,j,ii)+inv(5,j,ii) %have enough on-hand inventory to fill all the demand during j
                 bo(1,j,ii)=0;%no backorder
                inv(1,j+1,ii)=inv(2,j,ii)-inv(4,j,ii)-inv(5,j,ii); %starting on-hand at the beginning of j+1
             elseif inv(2,j,ii)<inv(4,j,ii)+inv(5,j,ii) && inv(2,j,ii)>0
              %if you don't have enough inventory to fill the entire demand
                inv(1,j+1,ii)=inv(2,j,ii)-inv(4,j,ii)-inv(5,j,ii);%next no on-hand inventory BEFORE delivery
                bo(1,j,ii)=inv(4,j,ii)+inv(5,j,ii)-inv(2,j,ii);
                bo(2,j,ii)=bo(1,j,ii);
             elseif inv(2,j,ii)<inv(4,j,ii)+inv(5,j,ii) && inv(2,j,ii)<=0
                 inv(1,j+1,ii)=inv(2,j,ii)-inv(4,j,ii)-inv(5,j,ii);
                 bo(1,j,ii)=inv(4,j,ii)+inv(5,j,ii);
                 bo(2,j,ii)=bo(1,j,ii);
             end
         elseif j==1
                 inv(2,j,ii)=inv(1,j,ii);
                 inv(3,j,ii)=inv(1,j,ii);
                 order(1,j,ii)=max(0,round(bs_level(ii,j))-inv(3,j,ii));
                 if inv(2,j,ii)>=inv(4,j,ii)+inv(5,j,ii) %have enough on-hand inventory to fill all the demand during j
                 bo(1,j,ii)=0;%no backorder
                inv(1,j+1,ii)=inv(2,j,ii)-inv(4,j,ii)-inv(5,j,ii); %starting on-hand at the beginning of j+1
                 elseif in(2,j,ii)<inv(4,j,ii)+inv(5,j,ii)
              %if you don't have enough inventory to fill the entire demand
                inv(1,j+1,ii)=0;%next no on-hand inventory BEFORE delivery
                bo(1,j,ii)=inv(4,j,ii)+inv(5,j,ii)-inv(2,j,ii);
                bo(2,j,ii)=bo(1,j,ii);
                end
                 
         elseif j>=L+1
            inv(2,j,ii)=inv(1,j,ii)+order(1,j-L,ii);
            for q=j-L+1:1:j-1
                onorder2=onorder2+order(1,q,ii);
            end
            sum_bo=0;
            for q1=1:1:j-1 %first fill the backordered items
                sum_bo=sum_bo+bo(2,q1,ii);
                if sum_bo<=inv(2,j,ii)
                   bo(2,q1,ii)=0;
                else
                  
                    break;
                end
             end
            inv(3,j,ii)=inv(2,j,ii)+onorder2;%inventory position 
             order(1,j,ii)=max(0,round(bs_level(ii,j-1))-inv(3,j,ii)); %since we order BFFORE filling the demand
             
             
              if inv(2,j,ii)>=inv(4,j,ii)+inv(5,j,ii) %have enough on-hand inventory to fill all the demand during j
                 bo(1,j,ii)=0;%no backorder
                inv(1,j+1,ii)=inv(2,j,ii)-inv(4,j,ii)-inv(5,j,ii); %starting on-hand at the beginning of j+1
             elseif inv(2,j,ii)<inv(4,j,ii)+inv(5,j,ii) && inv(2,j,ii)>0
              %if you don't have enough inventory to fill the entire demand
                inv(1,j+1,ii)=inv(2,j,ii)-inv(4,j,ii)-inv(5,j,ii);%next no on-hand inventory BEFORE delivery
                bo(1,j,ii)=inv(4,j,ii)+inv(5,j,ii)-inv(2,j,ii);
                bo(2,j,ii)=bo(1,j,ii);
             elseif inv(2,j,ii)<inv(4,j,ii)+inv(5,j,ii) && inv(2,j,ii)<=0
                 inv(1,j+1,ii)=inv(2,j,ii)-inv(4,j,ii)-inv(5,j,ii);
                 bo(1,j,ii)=inv(4,j,ii)+inv(5,j,ii);
                 bo(2,j,ii)=bo(1,j,ii);
             end
         end
     %frate(1,j,ii)=1-(bo(1,j,ii)/(inv(4,j,ii)+inv(5,j,ii)));
     
    end
    avgbs=sum(bs_level2(ii,:))/(obs_prd-5);
    
    if bsmax<max(bs_level2(ii,:))
        bsmax=max(bs_level2(ii,:));
    end
    if bsmin>min(bs_level2(ii,:))
        bsmin=min(bs_level2(ii,:));
    end
        
    bsgap=bsmax-bsmin;
  %cputime-t
  for_table2{1,ii}=pr;
  for_table2{2,ii}=upr;
  for_table2{3,ii}=pr_wo_re;
  for_table2{4,ii}=bsgap*100/avgbs;
  bos=0;dmd=0;
  ohi=0;bsm=0;
  for i=1:1:obs_prd
      bos=bos+bo(1,i,ii);
      dmd=dmd+inv(4,i,ii)+inv(5,i,ii);
      ohi=ohi+inv(1,i,ii);
      bsm=bsm+bs_level(ii,i);
  end
  oh_inv(1,1,ii)=ohi/obs_prd;
  bs_mean(1,1,ii)=bsm/obs_prd;
  ohi=0;bsm=0;
  bos_L=0;dmd_L=0;
  for i=L+1:1:obs_prd
      bos_L=bos_L+bo(1,i,ii);
      dmd_L=dmd_L+inv(4,i,ii)+inv(5,i,ii);
      ohi=ohi+inv(1,i,ii);
      bsm=bsm+bs_level(ii,i);
  end
  oh_inv(1,3,ii)=ohi/(obs_prd-L);
  bs_mean(1,3,ii)=bsm/(obs_prd-L);
  ohi=0;
  bos_L2=0;dmd_L2=0;
  for i=L/2+1:1:obs_prd
      bos_L2=bos_L2+bo(1,i,ii);
      dmd_L2=dmd_L2+inv(4,i,ii)+inv(5,i,ii);
      ohi=ohi+inv(1,i,ii);
      bsm=bsm+bs_level(ii,i);
  end
  oh_inv(1,2,ii)=ohi/(obs_prd-L/2);
  bs_mean(1,2,ii)=bsm/(obs_prd-L/2);
  ohi=0;bsm=0;
 frate(1,1,ii)=1-bos/dmd;
 frate(1,2,ii)=1-bos_L2/dmd_L2;
 frate(1,3,ii)=1-bos_L/dmd_L;
 
  bos_mo=0;dmd_mo=0;cot=1;
  for i=1:4:obs_prd-3    
      for j=0:1:3
           bos_mo=bos_mo+bo(1,i+j,ii);
           dmd_mo=dmd_mo+inv(4,i+j,ii)+inv(5,i+j,ii);
           ohi=ohi+inv(1,i+j,ii);
           bsm=bsm+bs_level(ii,i+j);
           if j==3
                frate_mo(1,cot,ii)=1-bos_mo/dmd_mo;
                oh_inv_mo(1,cot,ii)=ohi/4;
                bs_mean_mo(1,cot,ii)=bsm/4;
                bos_mo=0; dmd_mo=0;ohi=0;bsm=0;
                cot=cot+1;
           end
      end
  end
  bos_qr=0;dmd_qr=0;cot=1;ohi=0;bsm=0;
  for i=1:12:obs_prd-11    
      for j=0:1:11
           bos_qr=bos_qr+bo(1,i+j,ii);
           dmd_qr=dmd_qr+inv(4,i+j,ii)+inv(5,i+j,ii);
           ohi=ohi+inv(1,i+j,ii);
           bsm=bsm+bs_level(ii,i+j);
           if j==11
                frate_qr(1,cot,ii)=1-bos_qr/dmd_qr;
                oh_inv_qr(1,cot,ii)=ohi/12;
                bs_mean_qr(1,cot,ii)=bsm/12;
                bos_qr=0; dmd_qr=0;ohi=0;bsm=0;
                cot=cot+1;
           end
      end
  end
              
  
  
end
 %frate(1,rep+1)=sum(frate)/rep;
 
for i=1:1:length(frate_mo)
    result_mo(2,i)=sum(frate_mo(1,i,:))/rep;
    
    result_mo(3,i)=sum(oh_inv_mo(1,i,:))/rep;
    result_mo(1,i)=sum(bs_mean_mo(1,i,:))/rep;
 end
 for i=1:1:length(frate_qr)
     result_qr(2,i)=sum(frate_qr(1,i,:))/rep;
     result_qr(3,i)=sum(oh_inv_qr(1,i,:))/rep;
     result_qr(1,i)=sum(bs_mean_qr(1,i,:))/rep;
 end
 for i=1:1:3
     result_lt(2,i)=sum(frate(1,i,:))/rep;
     result_lt(3,i)=sum(oh_inv(1,i,:))/rep;
     result_lt(1,i)=sum(bs_mean(1,i,:))/rep;
 end
 
%bs_level_mean=zeros(obs_prd,1);
%k2matrix_50_4_7=zeros(obs_prd-5,1);
%for j=1:1:obs_prd
   % bs_level_mean(j,1)=mean(bs_level(:,j));
    %k2matrix_1500_3_10(j,1)=mean(k2matrix(:,j));
%end

table2=cell2mat(for_table2);
table3=zeros(4,1);
table3(1,1)=mean(table2(1,:));
table3(2,1)=mean(table2(2,:));
table3(3,1)=mean(table2(3,:));
table3(4,1)=mean(table2(4,:));



    
