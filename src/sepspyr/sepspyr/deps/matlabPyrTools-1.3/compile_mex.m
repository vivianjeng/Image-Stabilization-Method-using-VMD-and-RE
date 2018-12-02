% JEBYRNE
cd('MEX');

fprintf('[matlabPyrTools.%s]: compiling corrDn.c\n', mfilename)
mex corrDn.c wrap.c convolve.c edges.c

fprintf('[matlabPyrTools.%s]: compiling upConv.c\n', mfilename)
mex upConv.c wrap.c convolve.c edges.c

fprintf('[matlabPyrTools.%s]: compiling pointOp.c\n', mfilename)
mex pointOp.c

fprintf('[matlabPyrTools.%s]: compiling histo.c\n', mfilename)
mex histo.c

fprintf('[matlabPyrTools.%s]: compiling range2.c\n', mfilename)
mex range2.c

