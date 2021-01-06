%% plotting nigative binomial distrebution pdf 
RR        = 1;  %avreage length of Type-L/H period
PnBin     = .2;
x = 0:1:25;
y = nbinpdf(x,RR,PnBin);
y1 = (PnBin).*((1-PnBin).^x);
figure;
hold on
for ind = 1:length(x)
plot([x(ind) x(ind)],[0 y(ind)],'-','Color',[.7 .7 .7],'LineWidth',1)
plot(x(ind),y(ind),'.','MarkerSize',20,'Color','k')
end
set(gca,'FontSize',20)
box on
axis([0 25 0 0.2])
xticks([0 5 10 15 20 25])
xticklabels([0 5 10 15 20])
yticks([0 .03 0.09 .15 0.2])
yticklabels([0 .03 0.09 .15 ])
ylabel('$f_{r,p}(x)$','Rotation',90)
xlabel('$x$')