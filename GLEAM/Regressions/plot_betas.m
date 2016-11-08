function plot_betas(means, lowers, uppers )
% plot_betas(means, lowers, uppers )
% plot posterior means and credible intervals

p = size(means,2);
y= zeros(2,p);
y(1,:) = uppers;
y(2,:) = lowers;
x = repmat(1:p,2,1);
h = plot(x,y);
xlim([0,p+2]);
set_linespec(h,'b-');
set(h,'linewidth',1)
hold on
h=plot(means,'r.');
set(h,'linewidth',2)
hold off
xlabel('Coefficients');
ylabel('Posterior summary')
%axis_pct
%set(gca,'ytick',1:p,'yticklabel',rownames,'ydir','reverse')
%set(gca,'tickdir','out')
