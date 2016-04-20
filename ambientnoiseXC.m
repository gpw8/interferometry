% create animation to demonstrate direct wave correlation for GF
clear, close all
% First, define station locations:
xA=-100;
xB=100;
yA=0;
yB=0;

%Next, set up an evenly distributed collection of sources ala Waapenar
% SEG 2005

Ns=20;
phi=linspace(0,2*pi,Ns);
r=210;
xs=r*cos(phi);
ys=r*sin(phi);

figure
plot(xA,yA,'k^',xB,yB,'b^')
hold on
text(xA,yA+1,'A');
text(xB,yB+1,'B');

plot(xs,ys,'o','MarkerFaceColor','r','MarkerEdgeColor','r')
axis equal

%%
% define a wavelet for plotting
nsamp=21;
delta=0.001;
s=stfunc(nsamp,delta,.02,.01,3);
t=(0:nsamp-1)*delta;
figure;
plot(t,s)

%%
% now compute tt for all points to A
dA=sqrt( (xA-xs).^2 + (yA-ys).^2);
c=1000; % wave speed equal 1km/s
ttA=dA/c;

% and tt for all points to B
dB=sqrt( (xB-xs).^2 + (yB-ys).^2);
c=1000; % wave speed equal 1km/s
ttB=dB/c;

%% 
% and the seismograms for A and B
tottime=ceil(max(ttB)/delta);
totsamp=tottime/delta;
seisA=zeros(Ns,totsamp);
seisB=seisA;

ttAsamp=round(ttA/delta);
ttBsamp=round(ttB/delta);
for n=1:length(ttAsamp)
    seisA(n,ttAsamp+1:ttAsamp+nsamp)=s;
    seisB(n,ttBsamp+1:ttBsamp+nsamp)=s;
end

figure

%%


