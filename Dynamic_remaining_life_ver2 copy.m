%Basic Parameters
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
E_no=50;%number of total engines
rep=2;
fr=0.80;%required fill rate
csl=0.97;
sigma=0.0005;
theta=sigma*2.2;
F=0.07; %threshold
cycle=floor(1.1*F/theta);
decision=cycle+1;
L=ceil(cycle/3);%lead time
delta=ceil(cycle/7);
%delta=ceil(cycle/8);
%delta=ceil(cycle/10);
%delta=ceil(cycle/25);
%delta=0;
obs_prd=3*cycle;
ini_oh=10;
avg_dd=zeros(3,obs_prd);
avg_ms=zeros(2,obs_prd);
avg_dd_mo=zeros(6,floor(obs_prd/4));
avg_dd_qt=zeros(6,floor(obs_prd/12));
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
ED=zeros(rep,1);%mean
VD=zeros(rep,1);%variance

result_matrix=zeros(2,obs_prd,rep);%matrix of mean and variance of all replications
bs_level=zeros(rep,obs_prd);
kmatrix=zeros(rep,obs_prd);

frate=zeros(1,rep);
svl=zeros(1,rep);
frate_mo=zeros(1,floor(obs_prd/4),rep);
frate_qr=zeros(1,floor(obs_prd/12),rep);
oh_inv=zeros(1,3,rep+1);
oh_inv_mo=zeros(1,floor(obs_prd/4),rep);
oh_inv_qr=zeros(1,floor(obs_prd/12),rep);
bs_mean=zeros(1,3,rep+1);
bs_mean_mo=zeros(1,floor(obs_prd/4),rep);
bs_mean_qr=zeros(1,floor(obs_prd/12),rep);

inv=zeros(8,obs_prd,rep);
result_mo=zeros(3,floor(obs_prd/4));
result_qr=zeros(3,floor(obs_prd/12));
result_lt=zeros(3,3);

scheduled=zeros(E_no,20,rep);%initially scheduled maintenance
repaired=zeros(E_no,20,rep);%actual maintenance 

for i=1:1:rep
    inv(1,1,i)=ini_oh; %initial start on-hand inventory-set by trial and error
    inv(1,2,i)=ini_oh;
end
bo=zeros(2,obs_prd,rep);%backorders: first row how much occured, second row after filled how much left. 
order=zeros(1,obs_prd,rep);%how much on-order

engine_set=randi([1,cycle],E_no,1,rep);%first uniformly distribute over [1,70] which will be next visits
w=zeros(E_no,cycle+obs_prd,rep);
z=zeros(E_no,cycle+obs_prd,rep);
base_signal=zeros(E_no,obs_prd,rep);
next_visit=zeros(E_no,1,rep);
%se=rng;

for ii=1:1:rep
    t=cputime;
    %rng(se);
    %rng
    ii
    
    
    for i=1:1:E_no
        next_visit(i,1,ii)=engine_set(i,1,ii)+cycle;
        scheduled(i,1,ii)=engine_set(i,1,ii);
        scheduled(i,2,ii)=next_visit(i,1,ii);
    end
   
    %generate standard brownian error until the beggining of period k(=time
    %k+1)+L(lead time) to see the fluctuation of base stock level
   
    
    for i=1:1:E_no%First generate enough singlas for long time
        z(i,1,ii)=theta*1;
        for t=2:1:cycle+obs_prd
            w(i,t,ii)=randn+w(i,t-1,ii);
            z(i,t,ii)=theta*t+sigma*w(i,t,ii);
            if z(i,t,ii)<0
                z(i,t,ii)
            end
        end
    end

    %want to observe the value at time k+1(=cycle+1)  
    for i=1:1:E_no
        for j=1:1:obs_prd
             base_signal(i,j,ii)=z(i,decision-engine_set(i,1)+j,ii);%cycle-engine+j=how many days passed since last visit at the time of observing singal at j             
        end
    end
    
    scheduled_ct=ones(E_no,rep);
    repaired_ct=ones(E_no,rep);
    
    %categorize the engines into 3 groups: planned/both/unplanned only
    %while the time moment to estimate the future demand, engine groups
    %should be changed based on the same rule.
    
    for j=1:1:obs_prd%so we are observing 3 cycles from 71 to 280
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
                      repaired(i,repaired_ct(i,ii),ii)=decision+j-1;
                      scheduled(i,scheduled_ct(i,ii)+2,ii)=decision+j-1+cycle;
                      repaired_ct(i,ii)=repaired_ct(i,ii)+1;
                      scheduled_ct(i,ii)=scheduled_ct(i,ii)+1;
                      base_signal(i,j+1,ii)=0;
                      w(i,j+1,ii)=0;
                      base_signal(i,j+1,ii)=theta*1+w(i,j+1,ii); %singal is reset from the current period, so next period get just theta*1.
                      for s=j+2:1:obs_prd
                          base_signal(i,s,ii)=0;%first reset the rest of the base_signal for the future
                          w(i,s,ii)=randn+w(i,s-1,ii);
                          z(i,s,ii)=theta*(s-j)+sigma*w(i,s,ii);
                           if z(i,j,ii)<0
                              z(i,j,ii)
                           end
                          base_signal(i,s,ii)=z(i,s,ii);%regeneratethe signal starting from j=theta (if signal>F in (j-1), took donw, replaced)
                      end
                      
                  elseif next_visit(i,1,ii)>decision+j-1+floor(delta)%too far away from planned maintenance
                      upr=upr+1;
                      next_visit(i,1,ii)=decision+j-1+cycle;
                      repaired(i,repaired_ct(i,ii),ii)=decision+j-1;
                      scheduled(i,scheduled_ct(i,ii)+2,ii)=decision+j-1+cycle;
                      repaired_ct(i,ii)=repaired_ct(i,ii)+1;
                      scheduled_ct(i,ii)=scheduled_ct(i,ii)+1;
                      
                      w(i,j+1,ii)=0;
                      base_signal(i,j+1,ii)=theta*1+w(i,j+1,ii);
                      
                      for s=j+2:1:obs_prd
                          base_signal(i,s,ii)=0;%first reset the rest of the base_signal for the future
                          w(i,s,ii)=randn+w(i,s-1,ii);
                          z(i,s,ii)=theta*(s-j)+sigma*w(i,s,ii);
                           if z(i,j,ii)<0
                               z(i,j,ii)
                           end
                          base_signal(i,s,ii)=z(i,s,ii);%regeneratethe signal starting from j=theta (if signal>F in (j-1), took donw, replaced)
                      end
                  end
              else%if the base signal is not reaching the threshold level
                  if next_visit(i,1,ii)==decision+j-1
                     pr_wo_re=pr_wo_re+1;
                     next_visit(i,1,ii)=decision+j-1+cycle;
                     repaired(i,repaired_ct(i,ii),ii)=decision+j-1;
                     scheduled(i,scheduled_ct(i,ii)+2,ii)=decision+j-1+cycle;
                     repaired_ct(i,ii)=repaired_ct(i,ii)+1;
                     scheduled_ct(i,ii)=scheduled_ct(i,ii)+1;
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
    egroup1=[];egroup2=[];egroup3=[];
    egroup1=cell2mat(engine_group1);
    egroup2=cell2mat(engine_group2);
    egroup3=cell2mat(engine_group3);

   phic={};upsilonc={};psic={};
   phi=[];upsilon=[];psi=[];
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
             bs_level(ii,j)=mu_L+kk*std;
             kmatrix(ii,j)=kk;
            break;
          end
    end 
       %periodical demand from plannde/unplanned maintenance
         inv(4,j,ii)=pr;
         inv(5,j,ii)=upr;
         inv(6,j,ii)=mu_L;
         inv(7,j,ii)=std;
         inv(8,j,ii)=pr_wo_re;
         
         onorder=0;
         onorder2=0;
         
         if j>1 && j<L+1 %before we get the first delivery
             inv(2,j,ii)=inv(1,j,ii); %starting OH inventory, there's no delivery during this period
             for q=1:1:j-1
                 onorder=onorder+order(1,q,ii);%every order upto L is on-order
             end
             inv(3,j,ii)=inv(2,j,ii)+onorder;%inventory position
             order(1,j,ii)=max(0,ceil(bs_level(ii,j-1))-inv(3,j,ii)); %since we order BFFORE filling the demand
             
             if inv(2,j,ii)>=inv(4,j,ii)+inv(5,j,ii) %have enough on-hand inventory to fill all the demand during j
                %no backorder
                inv(1,j+1,ii)=inv(2,j,ii)-inv(4,j,ii)-inv(5,j,ii); %starting on-hand at the beginning of j+1
             elseif inv(2,j,ii)<inv(4,j,ii)+inv(5,j,ii) && inv(2,j,ii)>0
              %if you don't have enough inventory to fill the entire demand
                inv(1,j+1,ii)=inv(2,j,ii)-inv(4,j,ii)-inv(5,j,ii);%next no on-hand inventory BEFORE delivery
                bo(1,j,ii)=-inv(1,j+1,ii);
                bo(2,j,ii)=bo(1,j,ii);
             elseif inv(2,j,ii)<inv(4,j,ii)+inv(5,j,ii) && inv(2,j,ii)<=0
                 inv(1,j+1,ii)=inv(2,j,ii)-inv(4,j,ii)-inv(5,j,ii);
                 bo(1,j,ii)=inv(4,j,ii)+inv(5,j,ii);
                 bo(2,j,ii)=bo(1,j,ii);
             end
          elseif j==1
                 inv(2,j,ii)=inv(1,j,ii);
                 inv(3,j,ii)=inv(1,j,ii);
                 order(1,j,ii)=max(0,ceil(bs_level(ii,j))-inv(3,j,ii));
                 if inv(2,j,ii)>=inv(4,j,ii)+inv(5,j,ii) %have enough on-hand inventory to fill all the demand during j
                    %no backorder
                    inv(1,j+1,ii)=inv(2,j,ii)-inv(4,j,ii)-inv(5,j,ii); %starting on-hand at the beginning of j+1
                 elseif inv(2,j,ii)<inv(4,j,ii)+inv(5,j,ii)
                     %if you don't have enough inventory to fill the entire demand
                     inv(1,j+1,ii)=0;%next no on-hand inventory BEFORE delivery
                     bo(1,j,ii)=inv(4,j,ii)+inv(5,j,ii)-inv(2,j,ii);
                     bo(2,j,ii)=bo(1,j,ii);
                 end     
         else
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
            order(1,j,ii)=max(0,ceil(bs_level(ii,j-1))-inv(3,j,ii)); %since we order BFFORE filling the demand
             
             
              if inv(2,j,ii)>=inv(4,j,ii)+inv(5,j,ii) %have enough on-hand inventory to fill all the demand during j
                %no backorder
                inv(1,j+1,ii)=inv(2,j,ii)-inv(4,j,ii)-inv(5,j,ii); %starting on-hand at the beginning of j+1
             elseif inv(2,j,ii)<inv(4,j,ii)+inv(5,j,ii) && inv(2,j,ii)>0
              %if you don't have enough inventory to fill the entire demand
                inv(1,j+1,ii)=inv(2,j,ii)-inv(4,j,ii)-inv(5,j,ii);%next no on-hand inventory BEFORE delivery
                bo(1,j,ii)=-inv(1,j+1,ii);
                bo(2,j,ii)=bo(1,j,ii);
             elseif inv(2,j,ii)<inv(4,j,ii)+inv(5,j,ii) && inv(2,j,ii)<=0
                 inv(1,j+1,ii)=inv(2,j,ii)-inv(4,j,ii)-inv(5,j,ii);
                 bo(1,j,ii)=inv(4,j,ii)+inv(5,j,ii);
                 bo(2,j,ii)=bo(1,j,ii);
             end
         end
    end
  bos=0;dmd=0;
  ohi=0;bsm=0;
  
  bos_L=0;dmd_L=0;svl_ct=0;skip=0;
  for i=ceil(cycle)+5:1:obs_prd
      bos_L=bos_L+bo(1,i,ii);
      dmd_L=dmd_L+inv(4,i,ii)+inv(5,i,ii);
      ohi=ohi+inv(1,i,ii);
      bsm=bsm+bs_level(ii,i);
      skip=skip+inv(8,i,ii);
      if bo(1,i,ii)>0
         svl_ct=svl_ct+1;
      end
  end
  svl(1,ii)=1-(svl_ct/(obs_prd-(cycle+4)));
  
  oh_inv(1,3,ii)=ohi/(obs_prd-ceil(cycle/2));
  bs_mean(1,3,ii)=bsm/(obs_prd-ceil(cycle/2));

 
 frate(1,ii)=1-bos_L/dmd_L;
 
  bos_mo=0;dmd_mo=0;cot=1;bsm=0;ohi=0;
  for i=ceil(cycle)+5:4:obs_prd-3    
      for j=0:1:3
          if i+j>obs_prd
              break;
          else
            bos_mo=bos_mo+bo(1,i+j,ii);
            dmd_mo=dmd_mo+inv(4,i+j,ii)+inv(5,i+j,ii);
            ohi=ohi+inv(1,i+j,ii);
            bsm=bsm+bs_level(ii,i+j);
            
                if j==3
                    if dmd_mo==0;
                        frate_mo(1,cot,ii)=1;
                    else
                        frate_mo(1,cot,ii)=1-bos_mo/dmd_mo;
                        oh_inv_mo(1,cot,ii)=ohi/4;
                        bs_mean_mo(1,cot,ii)=bsm/4;
                        bos_mo=0; dmd_mo=0;ohi=0;bsm=0;
                        cot=cot+1;
                    end
                end
           end
      end
  end
  bos_qr=0;dmd_qr=0;cot=1;ohi=0;bsm=0;
  for i=ceil(cycle)+5:12:obs_prd-11    
      for j=0:1:11
          if i+j>obs_prd
              break;
          else
              bos_qr=bos_qr+bo(1,i+j,ii);
              dmd_qr=dmd_qr+inv(4,i+j,ii)+inv(5,i+j,ii);
              ohi=ohi+inv(1,i+j,ii);
              bsm=bsm+bs_level(ii,i+j);
             
             if j==11
                 if dmd_qr==0
                     frate_qr(1,cot,ii)=1;
                 else
                    frate_qr(1,cot,ii)=1-bos_qr/dmd_qr;
                    oh_inv_qr(1,cot,ii)=ohi/12;
                    bs_mean_qr(1,cot,ii)=bsm/12;
                    bos_qr=0; dmd_qr=0;ohi=0;bsm=0;
                    cot=cot+1;
                 end
             end
          end
      end
  end
  
end
 
    for i=1:1:obs_prd
     avg_dd(1,i)=mean(inv(4,i,:));%weekly average planned demand
     avg_dd(2,i)=mean(inv(5,i,:));%weekly average unplanned demand
     avg_dd(3,i)=mean(inv(4,i,:))+mean(inv(5,i,:));%weekly average total demand
     avg_ms(1,i)=mean(inv(6,i,:));%average leadtime demand ean over replications
     avg_ms(2,i)=mean(inv(7,i,:));%average leadtime demand std over replications
     avg_skip(1,i)=mean(inv(8,i,:));
     end 
     
 ct=1;
 for i=ceil(cycle)+5:4:obs_prd-3
     dd_mo=0;dd_mo2=0;ms_mo=0;ms_mo2=0;skip_mo=0;
     for j=0:1:3
         if i+j>obs_prd
              break;
         else
            dd_mo=dd_mo+avg_dd(1,i+j);
            dd_mo2=dd_mo2+avg_dd(2,i+j);
            ms_mo=ms_mo+avg_ms(1,i+j);
            ms_mo2=ms_mo2+avg_ms(2,i+j);
            skip_mo=skip_mo+avg_skip(1,i+j);
                if j==3
                    avg_dd_mo(1,ct)=dd_mo;
                    avg_dd_mo(2,ct)=dd_mo2;
                    avg_dd_mo(3,ct)=(dd_mo+dd_mo2);
                    avg_dd_mo(4,ct)=ms_mo/4;
                    avg_dd_mo(5,ct)=ms_mo2/4;
                    avg_dd_mo(6,ct)=skip_mo/4;
                    ct=ct+1;
                end
          end
      end
  end
  ct=1;
 for i=ceil(cycle)+5:12:obs_prd-11
     dd_qt=0;dd_qt2=0;ms_qt=0;ms_qt2=0;skip_qr=0;
     for j=0:1:11
         if i+j>obs_prd
              break;
         else
            dd_qt=dd_qt+avg_dd(1,i+j);
            dd_qt2=dd_qt2+avg_dd(2,i+j);
            ms_qt=ms_qt+avg_ms(1,i+j);
            ms_qt2=ms_qt2+avg_ms(2,i+j);
            skip_qr=skip_qr+avg_skip(1,i+j);
                if j==11
                    avg_dd_qt(1,ct)=dd_qt;
                    avg_dd_qt(2,ct)=dd_qt2;
                    avg_dd_qt(3,ct)=(dd_qt+dd_qt2);
                    avg_dd_qt(4,ct)=ms_qt/12;
                    avg_dd_qt(5,ct)=ms_qt2/12;
                    avg_dd_qt(6,ct)=skip_qr/12;
                    ct=ct+1;
                end
          end
      end
  end        
     
 %chi2gof(avg_dd_mo(1,:))
 %chi2gof(avg_dd_mo(3,:))
     
for i=1:1:length(result_mo)
    result_mo(2,i)=sum(frate_mo(1,i,:))/rep;
    
    result_mo(3,i)=sum(oh_inv_mo(1,i,:))/rep;
    result_mo(1,i)=sum(bs_mean_mo(1,i,:))/rep;
 end
 for i=1:1:length(result_qr)
     result_qr(2,i)=sum(frate_qr(1,i,:))/rep;
     result_qr(3,i)=sum(oh_inv_qr(1,i,:))/rep;
     result_qr(1,i)=sum(bs_mean_qr(1,i,:))/rep;
 end
 service_level=sum(svl(1,:))/rep
 for i=1:1:3
     %result_lt(2,i)=service_level;
     result_lt(2,i)=sum(frate(1,:))/rep;
     result_lt(3,i)=sum(oh_inv(1,i,:))/rep;
     result_lt(1,i)=sum(bs_mean(1,i,:))/rep;
 end
 
 
 




    
