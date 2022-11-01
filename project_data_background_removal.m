function IMDATA = project_data_background_removal(IMDATA, sig_readnoise)

nIMs = numel(IMDATA);

for iIM = 1 : nIMs
   IMDATA(iIM).Ib = IMDATA(iIM).FITSREADDATA - sig_readnoise^2;
end