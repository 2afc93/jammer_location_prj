load('Ducati_jam2_LHCP.mat')

gain_aug = [];
for ii = 1:size(gainvals,1)
    for jj = 1:size(gainvals,2)
        gain_aug = [gain_aug gainvals/gainvals(ii,jj)];
    end
end
clear jj ii

cov = 1/numel(gain_aug)*sum(sum((gain_aug - mean(gain_aug)).*(gain_aug - mean(gain_aug))));
