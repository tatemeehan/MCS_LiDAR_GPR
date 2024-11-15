%% Determine Phase Wrapping due to LWC
clear
f = 1.3; % GHz
c = 0.3; % speedolight
lambda = c./f;
lwc = 0.0001:0.0001:0.1;
rho = 340;
hs = 1.65;
meaninc = 65;
meanphz = -0.75;
% AM Refraction
amLWC = 0;
amV = WetCrimVRMS(rho,amLWC);
amPerm = (c./amV).^2;
amindex = sqrt(amPerm);
theta2am = asind(sind(meaninc)./amindex);
deltaphz = zeros(numel(lwc),1);
for kk = 1:numel(lwc)
% PM refraction
pmV = WetCrimVRMS(rho,lwc(kk));
pmPerm = (c./pmV).^2;
pmindex = sqrt(pmPerm);
theta2pm = asind(sind(meaninc)./pmindex);
% deltaL = depth(ii).*(1./sind(theta2pm)-1./sind(theta2am(ii)));
l1 = hs./cosd(theta2am);
l2 = hs./cosd(theta2pm);
deltaL = l2-l1;
deltaT = ((l2./pmV)-(l1./amV));
delL(kk) = deltaL;
delT(kk) = deltaT;

% deltaL = hs.*(1./cosd(theta2pm)-1./cosd(theta2am));
% Phase Retrieval
% Wrapped
% deltaLphz = 2.*pi.*((sign(deltaL)).*(abs(deltaL./lambda)-abs(floor(deltaL./lambda))));
% % deltaLphz = wrapToPi(2.*pi.*(deltaL./lambda));
% deltaCphz = wrapToPi(-2.*pi.*f.*(deltaL.*(1./pmV-1./amV)));
% deltaphz(kk) = wrapToPi(deltaLphz+deltaCphz);
% UnWrapped
% deltaLphz = (2.*pi.*(deltaL./lambda));
% deltaCphz = (-2.*pi.*f.*(deltaL.*(1./pmV-1./amV)));
% Zach's Correciton
deltaLphz = 2.*pi.*(2.*deltaL./lambda);
deltaTphz = (2.*pi.*f.*(2.*deltaT));
% deltaCphz = (-2.*pi.*f.*((deltaL+(hs./cosd(theta2am))).*(1./pmV-1./amV)));
% deltaphz(kk) = deltaLphz+deltaCphz;
deltaphz(kk) = wrapToPi(deltaLphz+deltaTphz);
tmpT(kk) = wrapToPi(deltaTphz);
tmpL(kk) = wrapToPi(deltaLphz);
end
[~,lwcIx] = min(deltaphz);
ambgLWC = lwc(lwcIx);
% Make Figure of the Error Plot
tmperror= abs(unwrap(deltaphz(:)-meanphz));
% tmperror= abs((deltaphz(:)-meanphz));

[~,solnIx] = min(tmperror); soln = lwc(solnIx);
figure();plot(lwc.*100,tmperror,'k','LineWidth',2);grid on; grid minor;
xlabel('\Delta LWC (%)')
ylabel('Error |(Measured \Delta\phi - Modeled \Delta\phi)|')
set(gca,'fontname','serif','fontweight','bold','fontsize',12)
% Make Figure Describing the Ambiguity of the Retrieval
figure();
plot(lwc.*100,deltaphz,'k','LineWidth',2); grid on; grid minor;
xlabel('\Delta LWC (%)')
ylabel('Modeled \Delta\phi')
set(gca,'fontname','serif','fontweight','bold','fontsize',12)
ylim([-pi pi])

figure();
subplot(1,2,1)
plot(lwc.*100,tmperror,'k','LineWidth',2);grid on; grid minor;
xlabel('\Delta LWC (%)')
ylabel('Error |(Measured \Delta\phi - Modeled \Delta\phi)|')
set(gca,'fontname','serif','fontweight','bold','fontsize',12)
xlim([0 0.5])
axis square
% Make Figure Describing the Ambiguity of the Retrieval
subplot(1,2,2)
plot(lwc.*100,deltaphz,'k','LineWidth',2); grid on; grid minor;
xlabel('\Delta LWC (%)')
ylabel('Modeled \Delta\phi')
set(gca,'fontname','serif','fontweight','bold','fontsize',12)
ylim([-pi pi])
xlim([0 0.5])
axis square
set(gcf,"Position",[200 -54 1047 420])
annotation("textbox",[.13 1 0 0],'String','a)',"linestyle",'none','FitBoxToText','on','FontName','serif','FontSize',12,'FontWeight','bold')
annotation("textbox",[.5725 1 0 0],'String','b)',"linestyle",'none','FitBoxToText','on','FontName','serif','FontSize',12,'FontWeight','bold')

figure();plot(lwc,delL,'k');hold on; plot(lwc,delT,'r');
legend('\Delta L','\Delta T')
figure();plot(lwc,tmpL,'k');hold on;plot(lwc,tmpT,'r');plot(lwc,deltaphz,'b')
legend('\Delta \phi L','\Delta \phi T','\Delta \phi')