load ~/M.mat

for icond = 1 : 4
  for isubj = 1 : 24
isubj
    for ivox = 1 : 90


    acf(:,ivox,isubj,icond) = autocorr(squeeze(M(isubj,icond,:,ivox)));
    
    end
  end
end

%%

t = 1:21;

for ivox = 1 : 90
  for isubj = 1 : 24
  
    x = acf(:,ivox,isubj,4); %if min(x)<0; x = x+abs(min(x)); end;
    f = @(lambda) sum((x'-exp(-lambda.*t)).^2);
    [l_hat(ivox,isubj,2),h]=fminsearch(f,[0.1]);


    x = acf(:,ivox,isubj,3); %if min(x)<0; x = x+abs(min(x)); end;
    f = @(lambda) sum((x'-exp(-lambda.*t)).^2);
    [l_hat(ivox,isubj,1),h]=fminsearch(f,[0.1]);
  
  end
end


% 
% plot(x); hold on
% plot(exp(-g.*t))

%% PERMUTATION

dat(:,:,:,1) = acf(:,:,:,3);
dat(:,:,:,2) = acf(:,:,:,4);

clear x 

nperm = 10000;

for iperm = 1 : nperm
  
  fprintf('Perm %d ...\n',iperm)
  
  idx1 = randi(2,[24,1]);
  idx2 = 3-idx1;
  
  for i = 1 : 24
    
    acf_perm(:,:,i,1) = dat(:,:,i,idx1(i));
    acf_perm(:,:,i,2) = dat(:,:,i,idx2(i));
    
  end
  
  t = 1:21;

  for ivox = 1 : 90
    for isubj = 1 : 24

      xy = acf_perm(:,ivox,isubj,2); %if min(x)<0; x = x+abs(min(x)); end;
      f = @(lambda) sum((xy'-exp(-lambda.*t)).^2);
      [l_hat_perm(ivox,isubj,2),h]=fminsearch(f,[0.1]);
      lag1(ivox,isubj,2) = xy(2);
      lag2(ivox,isubj,2) = xy(3);


      xy = acf_perm(:,ivox,isubj,1); %if min(x)<0; x = x+abs(min(x)); end;
      f = @(lambda) sum((xy'-exp(-lambda.*t)).^2);
      [l_hat_perm(ivox,isubj,1),h]=fminsearch(f,[0.1]);
      lag1(ivox,isubj,1) = xy(2);
      lag2(ivox,isubj,1) = xy(3);

    end
  end
  
  [~,~,~,s]=ttest(l_hat_perm(:,:,1),l_hat_perm(:,:,2),'dim',2);
  max_tstat(iperm) = max(s.tstat);
 
end

save('~/rudy_autocorr.mat','max_tstat','acf','lag1','lag2','l_hat');

%% POWER SPECTRA
fsample = 1 / 2.2;
freq = 0 : fsample / 211 : fsample/2;

for icond = 1 : 4
  for isubj = 1 : 24
    isubj
      
      tmp = abs(fft( squeeze( M(isubj,icond,:,:) ))).^2;
      pow(:,:,isubj,icond) = tmp(1:length(freq),:,:,:);
      
      for ivox = 1 : 90
        X = [ones(length(freq(6:end)),1) log10(freq(6:end))'];
        Y = log10(pow(6:end,ivox,isubj,icond));
        tmp_reg = X\Y;
        slp(ivox,isubj,icond) = tmp_reg(2);
      end
      
  end
end
      
 

      
      
      
      
      




