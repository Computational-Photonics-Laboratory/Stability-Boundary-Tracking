function stable = fsolve_func(alpha,F)


clear all hidden; clc; close all hidden; drawnow

load hbd_Increas200_08_S2U.mat
%{
diffa = 4.4 - theta;
alpha = theta + diffa;
diffa = 2.3 - hbd0;
F = hbd0 + diffa;
%}
%fprintf('alpha = %15.10f\n',alpha);
%fprintf('F = %15.10f\n',F);

dL = L/N;
x = (-N/2:N/2-1)'*dL;

% Fourier spectral differentiation d_xx
%calculates second spatial derivative of del phi and del phi conj
[~, D2_fourdif] = fourdif(N, 2);
D2 = (2*pi/L)^2*D2_fourdif;
%----------------------------------------------------------------------  
options = optimset('Display','off','Jacobian','on','MaxIter',50);
uout = fsolve(@(u) LLPR_rhs(u,alpha,F,N,D2),uout,options);	% call fsolve
%----------------------------------------------------------------------
urout = uout(1:N);
uiout = uout(N+1:N*2);
%{
figure(15);
subplot(2,1,1);
plot(x,urout,'b');
ylabel('real part');

subplot(2,1,2);
plot(x,uiout,'b');
ylabel('imaginary part');
xlabel('x');
%}


%Figure A
%calculate phi(x)
uout = urout + 1j*uiout;
as_uout = abs(uout).^2;
%plot
%{
figure(1)
hold on
plot(x,as_uout,'Color','red')
xlim([-25 25])
ylim([0 12])
set(gca,'TickLabelInterpreter','LaTeX')
%x and y ticks
yticks([0:3:12])
yticklabels({'$0$','','$6$','','$12$'})
xticks([-20:10:20])
xticklabels({'$-20$','','$0$','','$20$'})
ylabel('$|\psi|^{2}$','interpreter','LaTeX',...
    'fontname','Times New Roman')
xlabel('$x$','interpreter','LaTeX',...
    'fontname','Times New Roman')
ax1 = gca;
ax1.FontSize = 32;
ax1.TickLength = [0.035 0.025];
%plotting continued
ax2 = axes('xlim',[-25 25],'ylim',[0 12],'color','none',...
    'YAxisLocation','right','XAxisLocation','top');
set(gca,'TickLabelInterpreter','LaTeX')
yticks([0:3:12])
yticklabels({'','','','',''})
%set(gca,'ytick',[])
xticks([-25 -20 -10 0 10 20 25])
xticklabels({'$-\pi$','','','$\theta$','','','$\pi$'})
%ax2 = gca;
ax2.FontSize = 32;
ax2.TickLength = [0.035 0.025];
saveas(gcf,'a.pdf','pdf')
%}

%Figure B
%calculate P_n
N_per = 8;
mode_range = 10;
x_dis = N_per * mode_range;
%take fft and abs sqrt
u_fft = ifftshift(ifft(uout));
u_absqrt = abs(u_fft).^2;
%normalize and convert to log form
u_norm = u_absqrt / max(u_absqrt);
u_log = 10*log10(u_norm);
%frequency domain
f = -N/2:N/2-1;
%plot
%{
figure(2)
stem(f,u_log,'Marker','None','Linewidth',2,'BaseValue',-250,'Color','red')
xlim([-x_dis,x_dis])
ylim([-60,5])
set(gca,'TickLabelInterpreter','LaTeX')
yticks([-60:10:0])
yticklabels({'$-60$','','','$-30$','','','$0$'})
xticks([-x_dis:x_dis/2:x_dis])
xticklabels({'$10$','','$0$','','$10$'})
ylabel('$P_{n}$ (dB)','interpreter','LaTeX',...
    'fontname','Times New Roman')
xlabel('Mode Number $n$','interpreter','LaTeX',...
    'fontname','Times New Roman')
ax1 = gca;
ax1.FontSize = 32;
ax1.TickLength = [0.035 0.025];
saveas(gcf,'b.pdf','pdf')
%}



%Figure C

%calculate eigenvalues
phi_0 = uout .* eye(N);
%create L matrix
L11 = i*(D2+2*abs(phi_0)^2-alpha*eye(N))-eye(N);
L12 = i*phi_0^2;
L21 = conj(L12);
L22 = conj(L11);
L = [L11 L12; L21 L22];
[V,D] = eig(L);
D = diag(D);
%take real e-values of L
e_r = real(D);
e_i = imag(D);


%plot
%{
y_line = [-20:0.5:20];
y_zeros = zeros(1,length(y_line));
x_line = [-2.25:0.05:0.25];
x_zeros = zeros(1,length(x_line));
figure(3)
hold on
scatter(e_r,e_i,'filled')
plot(y_zeros,y_line,'--k')
plot(x_line,x_zeros,'--k')
hold off
xlim([-2.25 0.25])
ylim([-20 20])
set(gca,'TickLabelInterpreter','LaTeX')
yticks([-20:10:20])
yticklabels({'$-20$','','$0$','','$20$'})
xticks([-2:0.5:0])
xticklabels({'$-2$','','$-1$','','$0$'})
ylabel('Im($\lambda$)','interpreter','LaTeX',...
    'fontname','Times New Roman')
xlabel('Re($\lambda$)','interpreter','LaTeX',...
    'fontname','Times New Roman')
ax1 = gca;
ax1.FontSize = 32;
ax1.TickLength = [0.035 0.025];
box on
saveas(gcf,'c.pdf','pdf')
%}


stable = 1
end
