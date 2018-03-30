%% pupmod_src_powcorr

clear

% --------------------------------------------------------
% VERSION 1 - WEIGHTED AAL
% --------------------------------------------------------
v               = 1;
v_postproc      = 6;
fsample         = 400;
SUBJLIST        = [4 5 6 7 8 9 10 11 12 13 15 16 19 20 21 22 23 24 25 26 27 28 29 30 31 32 33 34];
allpara.filt    = 'jh_lcmv';
allpara.grid    = 'vtpm_6mm';
allpara.reg     = 0.05;
% --------------------------------------------------------

outdir   = '/home/tpfeffer/funchier/proc/';

if strcmp(allpara.grid,'xcoarse')
  v_grid = 2;
elseif strcmp(allpara.grid,'aal')
  v_grid = 4;
elseif strcmp(allpara.grid,'cortex')
  v_grid = 3;
elseif strcmp(allpara.grid,'medium')
  v_grid = 5;
elseif strcmp(allpara.grid,'aal_6mm')
  v_grid = 6;
elseif strcmp(allpara.grid,'aal_4mm')
  v_grid = 7;
elseif strcmp(allpara.grid,'m758_4mm')
  v_grid = 8;
elseif strcmp(allpara.grid,'cortex_lowres')
  v_grid = 9;
elseif strcmp(allpara.grid,'vtpm_4mm')
  v_grid = 10;
elseif strcmp(allpara.grid,'vtpm_6mm')
  v_grid = 11;
end

%% LOAD DATA COMPUTE SRC TIME COURSES

for isubj = SUBJLIST
  for m = 1:3

    if ~exist(sprintf([outdir 'funchier_granger_s%d_m%d_v%d_processing.txt'],isubj,m,v))
      system(['touch ' outdir sprintf('funchier_granger_s%d_m%d_v%d_processing.txt',isubj,m,v)]);
    else
      continue
    end
    %
    fprintf('Processing s%d m%d...\n', isubj,m)
    
    for iblock = 1:2
      
      fprintf('Loading MEG data ...\n');
      
      load(sprintf('/home/tpfeffer/pconn_cnt/proc/preproc/pconn_cnt_postpostproc_s%d_m%d_b%d_v%d.mat',isubj,m,iblock,1))
%       load(sprintf('/home/tpfeffer/pconn/proc/preproc/pconn_postpostproc_s%d_m%d_b%d_v%d.mat',isubj,m,iblock,1))
      
      [dat] = megdata2mydata(data); clear data
      
      pars      = [];
      pars.sa   = sprintf('~/pconn_cnt/proc/src/pconn_cnt_sa_s%d_m%d_b%d_v%d.mat',isubj,m,iblock,v_grid);
%       pars.sa   = sprintf('~/pconn/proc/src/pconn_sa_s%d_m%d_b%d_v%d.mat',isubj,m,iblock,v_grid);
      
      load(pars.sa);
      
      para = [];
      para.epleng   = size(dat,1);
      para.segleng  = 400;
      para.fsample  = 400;
      para.segshift = para.segleng / 2;
      para.grid     = allpara.grid;
           
      outp = funchier_computegranger(dat,sa,para);
      
      load
   
      save(sprintf([outdir 'funchier_granger_s%d_m%d_b%d_f%d_v%d.mat'],isubj,m,iblock,ifoi,v),'outp');
      
    end
  end
end



error('!')




