% JBYRNE 
function showSpyrLevel(pyr, pind, level);

%% look at all orientation bands at one level (scale):
for b = 1:spyrNumBands(pind)
  band = spyrBand(pyr,pind,level,b);
  subplot(nrows,ceil(nfilts/nrows),b);
  showIm(band);
end

