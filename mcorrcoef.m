function [r, w] = mcorrcoef(X)
% MCORRCOEF Multiple correlation coefficients.
%   R=MCORRCOEF(X) calculates the r (MCC) and w (MUC) of correlation coefficients for an array X, in
%   which each row is an observation and each column is a variable.
%   References:
%      Cross Multivariate Correlation Coefficients as Screening Tool for Analysis of Concurrent EEG-fMRI Recordings
%      Wang, J., Zheng, N., 2014. Measures of linear correlation for multiple variables. arXiv preprint arXiv:1401.4827 .

[MF,MFp] = corrcoef(X);
w = sqrt(det(MF));
r = sqrt(1-w^2);