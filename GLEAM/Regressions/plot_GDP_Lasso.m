function plot_GDP_Lasso(means, uppers, lowers, means_lasso,odds)

nonnull_idx=[100, 200,201, 300, 301,302, 400,401,402,403, 500,501,502,503,504];
null_idx=ones(size(means));
null_idx(nonnull_idx)=0;


subplot(2,1,1);
plot_betas(means(nonnull_idx), lowers(nonnull_idx), uppers(nonnull_idx) );

hold on
plot(means_lasso(nonnull_idx), 'g.')
plot(1:length(means_lasso(nonnull_idx)),repmat(log(odds),1,length(means_lasso(nonnull_idx))),'*');
hline=refline(0,0);
set(hline,'Color','k')
title('Non-null ')
hold off

subplot(2,1,2);

plot_betas(means(find(null_idx)),lowers(find(null_idx)), uppers(find(null_idx)));
hold on
plot(means_lasso(find(null_idx)), 'g.')
title('Null ')
hold off