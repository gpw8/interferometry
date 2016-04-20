% ambientnoiseXC
% create animation to demonstrate direct wave correlation for GF
% This shows the relationship of direct arrivals from a distribution of
% sources to the cross correlation function

% Try a small number of sources first (Ns=20), so you can see what is going
% on. Then increase the number of sources (Ns=1000) for more effective
% interference. You can see how the source XC changes with increasing
% numbers of sources
clear, close all
Ns=400;

% First, define station locations:
xA=-100; % distance units can be anything here
xB=100;
yA=0;
yB=0;

%Next, set up an evenly distributed collection of sources ala Waapenar
% SEG 2005

% Here phi is the angle measured CCW from positive X.
phi=linspace(0,2*pi,Ns);
% Try removing a section of sources to see how it influences the final XC
% function
% phi=linspace(0.25,1.25*pi,Ns);

r=250;          % this is the radius of the circle of sources
xs=r*cos(phi);  % x position of sources
ys=r*sin(phi);  % y position of sources

figure
plot(xA,yA,'k^',xB,yB,'b^')
hold on
text(xA,yA+1,'A');
text(xB,yB+1,'B');

plot(xs,ys,'o','MarkerFaceColor','r','MarkerEdgeColor','r')
axis equal
xlabel('arbitrary distance')
ylabel('arbitrary distance')
legend('station A','station B','sources')
drawnow
%%
% define a wavelet for plotting - A ricker wavelet works well (wavetype 2)
nsamp=41;
delta=0.001;
s=stfunc(nsamp,delta,.02,.02,2);
t=(0:nsamp-1)*delta;
figure;
subplot(211)
plot(t,s)
xlabel('time')
ylabel('arbitrary amplitude')
title('source wavelet')

% plot the derivative of the source wavelet too
ds=diff(s);
ds(nsamp)=0;
subplot(212)
plot(t,ds)
xlabel('time')
ylabel('arbitrary amplitude')
title('d/dt of the source wavelet')

drawnow
%%
% now compute traveltimes for all points to A
dA=sqrt( (xA-xs).^2 + (yA-ys).^2); % dA is the distance between A and every source
c=1000; % wave speed equal 1000 (if everything is m and s, this is 1 km/s)
ttA=dA/c; % and the travel times

% and do the same for all points to station B
dB=sqrt( (xB-xs).^2 + (yB-ys).^2);
ttB=dB/c;

%% 
% and the seismograms for A and B
totsamp=ceil(max(ttB)/delta)+nsamp;
tottime=totsamp*delta;
seisA=zeros(Ns,totsamp);
seisB=seisA;

ttAsamp=round(ttA/delta);
ttBsamp=round(ttB/delta);
for n=1:Ns
    seisA(n,ttAsamp(n)+1:ttAsamp(n)+nsamp)=s;
    seisB(n,ttBsamp(n)+1:ttBsamp(n)+nsamp)=s;
end

% do not plot these if the number of sources is more than 500
if Ns<500
figure
offset=rad2deg(phi)'*ones(1,totsamp);
tvec=(0:totsamp-1)*delta;
normval=offset(2,1)-offset(1,1);
subplot(211)
plot(tvec,normval*seisA+offset,'k')
ylabel('center-to-source angle')
title('receiver A')

subplot(212)
plot(tvec,normval*seisB+offset,'k')
ylabel('center-to-source angle')
xlabel('time')
title('receiver B')
drawnow
end
%%
% now cross correlate
xc=zeros(Ns,(totsamp*2-1));
for n=1:Ns
    [xc(n,:),lag]=xcorr(seisA(n,:),seisB(n,:),'coeff');
end

if Ns<500

figure
offset2=rad2deg(phi)'*ones(1,size(xc,2));

subplot(4,1,1:3)
plot(lag*delta,normval*xc+offset2,'k')
ylabel('center-to-source angle')
xlabel('time lag')
axis tight

subplot(414)
plot(lag*delta,sum(xc)/Ns)
xlabel('time lag')
xlim([min(lag*delta),max(lag*delta)])
end

figure
plot(lag*delta,sum(xc)/Ns)
xlabel('time lag')
xlim([min(lag*delta),max(lag*delta)])
title('stack of all cross correlations')