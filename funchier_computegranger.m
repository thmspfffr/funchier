function outp = funchier_computegranger(dat,sa,para)

% Computes Granger causality as described in Dhamala et al. (2008).
% The function first computes X using fieldtrips inbuilt function.
% Then it esitimates directed interactions (granger causality) from the
% cross spectra. The measure can be computed on band-limited or broadband
% signals. The latter is the default setting, as this is used in 
% Bastos et al. (2015) Neuron and Michalareas et al. (2016) Neuron.

% -----
% thms.pfffr@gmail.com, 03/2018 
% -----

f = 0: (para.fsample / para.segleng) : para.fsample/2;

para.segshift = para.segleng/2;
para.zeropad = 0;
para.mydetrend = 0;

para.maxfreqbin = find(f<=95,1,'last');
para.minfreqbin = find(f>=2,1,'first');
para.bsfreq     = [find(f>=45,1,'first') find(f<=55,1,'last')];

fprintf('Computing cross spectrum ...\n')

cs = data2cs_event(dat,para.segleng,para.segshift,para.epleng,para.maxfreqbin,para);
% cs = cs(:,:,1:end);

fprintf('Computing cross spectrum ... Done!\n')


% para      = [];
para.iscs = 1;
para.reg  = 0.05;

fprintf('Computing spatial filter ...\n')

if strcmp(para.grid,'vtpm_4mm')
  pos           = sa.grid_vtpm_4mm_indi;
  filt          = pconn_beamformer(nanmean(cs(:,:,para.minfreqbin:end),3),sa.L_vtpm_4mm,para);
  vtpm_4mm      = tp_create_grid('vtpm');
  sa.vtpm_label = vtpm_4mm.tissue_4mm(vtpm_4mm.tissue_4mm>0);
elseif strcmp(para.grid,'vtpm_6mm')
  pos           = sa.grid_vtpm_6mm_indi;
  filt          = pconn_beamformer(nanmean(cs(:,:,para.minfreqbin:end),3),sa.L_vtpm_6mm,para);
  vtpm          = tp_create_grid('vtpm');
  sa.vtpm_label = vtpm.tissue_6mm(vtpm.tissue_6mm>0);
end

fprintf('Computing spatial filter ... Done!\n')

for ifoi = 1 : size(cs,3)
  cs_src(:,:,ifoi) = filt'*cs(:,:,ifoi)*filt;
end

[H, Z] = funchier_sfactorization_wilson(cs_src,f(1:para.maxfreqbin));

foi_cnt = 0;

for ifoi = para.minfreqbin:para.maxfreqbin
  
  foi_cnt = foi_cnt+1;
  
  fprintf('Estimating granger causality for %d Hz...\n',f(ifoi))
  
  for iloc = 1 : size(cs_src,1)
    c = real(cs_src(iloc,iloc,ifoi));
    for jloc = 1 : size(cs_src,1)
   
      tmp_granger(iloc,jloc,foi_cnt) = log (  c / (c - (Z(jloc,jloc) - Z(iloc,jloc)^2 / Z(iloc,iloc) ) * abs( H(iloc,jloc,ifoi) )^2 ));
      
    end
  end
end

outp.granger = tmp_granger;
outp.granger(:,:,para.bsfreq(1):para.bsfreq(2)) = [];
outp.freq    = f(para.minfreqbin:para.maxfreqbin);
outp.freq(para.bsfreq(1):para.bsfreq(2)) = [];

 