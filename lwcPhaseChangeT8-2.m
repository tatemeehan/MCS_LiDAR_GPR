%% Determine Phase Wrapping due to LWC
clear; close all;
f = 1.3; % GHz
c = 0.3; % speedolight
lambda = c./f;
lwc = 0.0001:0.0001:0.1;
rho = 340;
hs = 1.65;
meaninc = 45;
meanphz = -0.75;

% AM Refraction
amLWC = 0;
amV = WetCrimVRMS(rho,amLWC);
amPerm = (c./amV).^2;
amindex = sqrt(amPerm);
theta2am = asind(sind(meaninc)./amindex);

% PM refraction
pmV = WetCrimVRMS(rho,lwc);
pmPerm = (c./pmV).^2;
pmindex = sqrt(pmPerm);
theta2pm = asind(sind(meaninc)./pmindex);

% Pathlengths
l1 = 2.*hs./cosd(theta2am);
l2 = 2.*hs./cosd(theta2pm);

% Change in PathLength
deltaL = l2-l1;
deltaT = ((l2./pmV)-(l1./amV));

% Phase Retrieval
deltaLphz = 2.*pi.*(deltaL./lambda);
deltaTphz = (2.*pi.*f.*deltaT); % no-negative

% Phase Change is the sum of contributions
deltaphz = wrapToPi(deltaLphz+deltaTphz);
tmpT = wrapToPi(deltaTphz); % Individual Wrapped Phase L Contribution
tmpL = wrapToPi(deltaLphz); % Individual Wrapped Phase T Contributions

% Find Point of Phase Ambiguity
[~,lwcIx] = min(deltaphz);
ambgLWC = lwc(lwcIx);

% Calculate Solution Error
tmperror= abs(unwrap(deltaphz(:)-meanphz));
% Minimize the Error and find the Solution
[~,solnIx] = min(tmperror); soln = lwc(solnIx);

%% Make Figures
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
xlim([0 5])
axis square
% Make Figure Describing the Ambiguity of the Retrieval
subplot(1,2,2)
plot(lwc.*100,deltaphz,'k','LineWidth',2); grid on; grid minor;
xlabel('\Delta LWC (%)')
ylabel('Modeled \Delta\phi')
set(gca,'fontname','serif','fontweight','bold','fontsize',12)
ylim([-pi pi])
xlim([0 5])
axis square
set(gcf,"Position",[200 -54 1047 420])
annotation("textbox",[.13 1 0 0],'String','a)',"linestyle",'none','FitBoxToText','on','FontName','serif','FontSize',12,'FontWeight','bold')
annotation("textbox",[.5725 1 0 0],'String','b)',"linestyle",'none','FitBoxToText','on','FontName','serif','FontSize',12,'FontWeight','bold')

% Plot of Phase Contributions
figure();
subplot(2,1,1)
plot(100.*lwc,tmpL,'color',[0.625 0.625 0.625],'linewidth',2);hold on;plot(100.*lwc,tmpT,'k','linewidth',2);
% xlabel('\Delta LWC (%)')
title('$\theta_1 = 45^\circ$','interpreter','latex')
ylabel('\Delta\phi')
legend('$\Delta \phi_{\Delta \ell}$','$\Delta \phi_{\Delta t}$','interpreter','latex','Location','southeast')
set(gca,'fontname','serif','fontweight','bold','fontsize',12)
ylim([-pi pi])
yticks([-pi, -pi/2, 0, pi/2  pi])
% yticklabels({'-\pi','0','\pi'})
set(gca, 'TickLabelInterpreter', 'latex', 'YTickLabel', {'$-\pi$','$-\frac{\pi}{2}$','$0$','$\frac{\pi}{2}$','$\pi$'})
grid on
subplot(2,1,2)
plot(100.*lwc,deltaphz,'k','LineWidth',2)
xlabel('\Delta LWC (%)')
ylabel('\Delta\phi')
legend('$\Delta \phi$','interpreter','latex','Location','southeast')
set(gca,'fontname','serif','fontweight','bold','fontsize',12)
ylim([-pi pi])
yticks([-pi, -pi/2, 0, pi/2  pi])
% yticklabels({'-\pi','-\pi/2','0','\pi/2','\pi'})
grid on
set(gca, 'TickLabelInterpreter', 'latex', 'YTickLabel', {'$-\pi$','$-\frac{\pi}{2}$','$0$','$\frac{\pi}{2}$','$\pi$'})

annotation("textbox",[.13 1 0 0],'String','a)',"linestyle",'none','FitBoxToText','on','FontName','serif','FontSize',12,'FontWeight','bold')
annotation("textbox",[.13 0.525 0 0],'String','b)',"linestyle",'none','FitBoxToText','on','FontName','serif','FontSize',12,'FontWeight','bold')

