clear all; close all; clc;
%Heisenberg-Weil sequence, 
%Delay-Doppler Channel Estimation in Almost Linear Complexity, TIT, 2013,
%Section IV
N=37; %sequence length

%Heisenberg sequence a
hsc=0; %[0,N-1] integer, slope of the curtain
hsb=1; %[0,N-1] integer
HSa= 1/sqrt(N)*exp(pi*1i*(hsc*(0:N-1).^2+2*hsb*(0:N-1))/(N)).';

%Heisenberg sequence b
hsc=1; %[0,N-1] integer, slope of the curtain
hsb=0; %[0,N-1] integer
HSb= 1/sqrt(N)*exp(pi*1i*(hsc*(0:N-1).^2+2*hsb*(0:N-1))/(N)).';

%Weil sequence a, belongs to WS type 1, IV-B.3
%fix a generator (primitive root) r
r=2;%37 has 12 primitive roots, they are 2, 5, 13, 15, 17, 18, 19, 20, 22, 24, 32, and 35.
k=18; %[0,N-2] integer
xi=exp(2*pi*1i*k/(N-1));
WSa=[0;(1/sqrt(N-1)*(xi.^(log(1:N-1)/log(r)))).'];

%WS type 2, IV-B.4
% WSo=[0 (1/sqrt(N-1)*(xi.^(log(1:N-1)/log(r))))];
% wsc=17;%[0,N-1] integer
% WS=(exp(pi*1i*(wsc*(0:N-1).^2)/N).*WSo).';

% Weil sequence b, belongs to WS type 3, IV-B.6
WSo=[0 (1/sqrt(N-1)*(exp(2*pi*1i*9/(N-1)).^(log(1:N-1)/log(2))))];
WSoo=exp(pi*1i*((0:N-1).^2)/N).*fft(exp(pi*1i*((0:N-1).^2)/N).*WSo)/sqrt(N);
WSb=(-WSoo).';

%%%%numerical values
% HeisenbergSeqa=[0.164398987305357 + 0.00000000000000i;0.162034257533211 + 0.0277835637146866i;0.155008097186865 + 0.0547678448867089i;0.143522636106984 + 0.0801765548769046i;0.127908290097311 + 0.103278731360378i;0.108614255463576 + 0.123409766781726i;0.0861955864536317 + 0.139990527904360i;0.0612972273573252 + 0.152544016419299i;0.0346354586335294 + 0.160709091318046i;0.00697729082546770 + 0.164250858261879i;-0.0208816009361213 + 0.163067427064302i;-0.0481397670578131 + 0.157192842886203i;-0.0740130397291050 + 0.146796106818556i;-0.0977570920025927 + 0.132176314028738i;-0.118688850725055 + 0.113754049337123i;-0.136206147308404 + 0.0920592877575547i;-0.149805041023599 + 0.0677161480811250i;-0.159094316457896 + 0.0414249381148817i;-0.163806738069993 + 0.0139420081012684i;-0.163806738069993 - 0.0139420081012683i;-0.159094316457896 - 0.0414249381148817i;-0.149805041023599 - 0.0677161480811250i;-0.136206147308404 - 0.0920592877575546i;-0.118688850725055 - 0.113754049337123i;-0.0977570920025927 - 0.132176314028738i;-0.0740130397291050 - 0.146796106818556i;-0.0481397670578131 - 0.157192842886203i;-0.0208816009361214 - 0.163067427064302i;0.00697729082546766 - 0.164250858261879i;0.0346354586335294 - 0.160709091318046i;0.0612972273573252 - 0.152544016419299i;0.0861955864536316 - 0.139990527904361i;0.108614255463576 - 0.123409766781726i;0.127908290097311 - 0.103278731360378i;0.143522636106984 - 0.0801765548769046i;0.155008097186864 - 0.0547678448867089i;0.162034257533211 - 0.0277835637146867i];
% HeisenbergSeqb=[0.164398987305357 + 0.00000000000000i;0.163806738069993 + 0.0139420081012683i;0.155008097186865 + 0.0547678448867089i;0.118688850725055 + 0.113754049337123i;0.0346354586335294 + 0.160709091318046i;-0.0861955864536317 + 0.139990527904361i;-0.163806738069993 + 0.0139420081012684i;-0.0861955864536318 - 0.139990527904360i;0.108614255463576 - 0.123409766781726i;0.136206147308405 + 0.0920592877575546i;-0.0977570920025926 + 0.132176314028738i;-0.108614255463576 - 0.123409766781726i;0.155008097186864 - 0.0547678448867091i;-0.0346354586335293 + 0.160709091318046i;-0.0977570920025930 - 0.132176314028738i;0.159094316457896 + 0.0414249381148817i;-0.159094316457896 + 0.0414249381148820i;0.136206147308404 - 0.0920592877575550i;-0.118688850725055 + 0.113754049337123i;0.118688850725055 - 0.113754049337123i;-0.136206147308404 + 0.0920592877575551i;0.159094316457896 - 0.0414249381148826i;-0.159094316457896 - 0.0414249381148816i;0.0977570920025931 + 0.132176314028738i;0.0346354586335286 - 0.160709091318046i;-0.155008097186865 + 0.0547678448867087i;0.108614255463576 + 0.123409766781726i;0.0977570920025916 - 0.132176314028739i;-0.136206147308405 - 0.0920592877575536i;-0.108614255463576 + 0.123409766781726i;0.0861955864536318 + 0.139990527904360i;0.163806738069993 - 0.0139420081012695i;0.0861955864536309 - 0.139990527904361i;-0.0346354586335289 - 0.160709091318046i;-0.118688850725056 - 0.113754049337122i;-0.155008097186864 - 0.0547678448867089i;-0.163806738069993 - 0.0139420081012665i];
% WeilSeqa=[0.00000000000000 + 0.00000000000000i;0.166666666666667 + 0.00000000000000i;-0.166666666666667 + 2.04107799857892e-17i;0.0439599026793467 - 0.160764749662979i;0.166666666666667 - 4.08215599715784e-17i;0.0884504430517670 + 0.141259678966519i;-0.0439599026793467 + 0.160764749662979i;-0.137063549309063 + 0.0948227885509584i;-0.166666666666667 + 6.12323399573677e-17i;-0.143476990143735 - 0.0848064329934490i;-0.0884504430517669 - 0.141259678966519i;-0.0211840958931224 - 0.165314886985319i;0.0439599026793469 - 0.160764749662978i;0.0981503799230996 - 0.134700707862761i;0.137063549309063 - 0.0948227885509586i;0.159587098770155 - 0.0480597095694775i;0.166666666666667 - 8.16431199431569e-17i;0.160414486202856 + 0.0452213488747468i;0.143476990143735 + 0.0848064329934488i;0.118615944746902 + 0.117081319729399i;0.0884504430517669 + 0.141259678966519i;0.0553131694513685 + 0.157220326494452i;0.0211840958931226 + 0.165314886985319i;-0.0123257481326882 + 0.166210269570648i;-0.0439599026793471 + 0.160764749662978i;-0.0727848961540201 + 0.149933774279267i;-0.0981503799230998 + 0.134700707862761i;-0.119646716960677 + 0.116027759171284i;-0.137063549309063 + 0.0948227885509587i;-0.150351366238695 + 0.0719183178886675i;-0.159587098770155 + 0.0480597095694775i;-0.164944124296634 + 0.0239000761043643i;-0.166666666666667 + 1.02053899928946e-16i;-0.165048343293420 - 0.0231694228213646i;-0.160414486202856 - 0.0452213488747468i;-0.153107809992317 - 0.0658466118880414i;-0.143476990143735 - 0.0848064329934488i];
% WeilSeqb=[-0.0806087187076362 + 0.0274664415081806i;-0.0687736189375451 - 0.00388600385253095i;-0.0457130473346778 - 0.0330428848720121i;-0.00741679564283217 - 0.0524953718492844i;0.0446561799899346 - 0.0492378061060088i;0.0965784623869390 - 0.00880987026196951i;0.118895297954699 + 0.0716700980539499i;0.0792372756381589 + 0.166762114002962i;-0.0270459174946510 + 0.222925834561872i;-0.153486753448656 + 0.192894681602810i;-0.221183373173322 + 0.0858225549672516i;-0.187748006565121 - 0.0116572102244877i;-0.133266696910033 - 0.0106190832233501i;-0.177211329997751 + 0.00325510777379904i;-0.192009572585984 - 0.0934285683434404i;-0.105578671899379 - 0.100345901160834i;-0.155450817799634 - 0.0842642599877488i;-0.0996024536875197 - 0.168413004471898i;-0.0919525402136634 - 0.105265357827785i;-0.0749810348565582 - 0.191715211925602i;-0.0542190956936035 - 0.124291869486906i;-0.0233041622066714 - 0.200057251466911i;-0.0405372474094591 - 0.143581360176980i;0.0345838757695672 - 0.173872798076888i;-0.00722086883802803 - 0.188115736220010i;0.0300620883688700 - 0.134201625337415i;0.0774614129207872 - 0.168979321658486i;0.0497000213438803 - 0.182255176649127i;0.0500936651189732 - 0.138113409199274i;0.0942965852383237 - 0.117482473465816i;0.124917524820508 - 0.133580147259422i;0.125260607679630 - 0.150478488253253i;0.112399576481292 - 0.150019277489997i;0.100514855650296 - 0.135144619133643i;0.0932304660849020 - 0.112676461822293i;0.0892573075391755 - 0.0864652773207572i;0.0860889724909810 - 0.0578977687946020i];

HWSa=(HSa+WSa)/sqrt(2);
HWSb=(HSb+WSb)/sqrt(2);

omegamax =N-1;taumax = N-1;osr = 1;
doppler0 = -omegamax:1/osr:omegamax;
delay0 = -taumax:taumax;
[dopplerx1,delayx1]=meshgrid(doppler0,delay0);

DDHSa = getaf(HSa.'/sqrt(2),HSa.'/sqrt(2),omegamax,osr,taumax,1,N).';
DDHSb = getaf(HSb.'/sqrt(2),HSb.'/sqrt(2),omegamax,osr,taumax,1,N).';
DDWSa = getaf(WSa.'/sqrt(2),WSa.'/sqrt(2),omegamax,osr,taumax,1,N).';
DDWSb = getaf(WSb.'/sqrt(2),WSb.'/sqrt(2),omegamax,osr,taumax,1,N).';
DDHWSa = getaf(HWSa.',HWSa.',omegamax,osr,taumax,1,N).';
DDHWSb = getaf(HWSb.',HWSb.',omegamax,osr,taumax,1,N).';


figure(1)
surf(delayx1, dopplerx1,DDHSa, 'FaceColor', 'interp', 'EdgeColor', 'none');
axis tight; grid on; colorbar;
xlabel('Delay', 'fontsize', 22); 
ylabel('Doppler', 'fontsize', 22); 
title('AAF of HSa', 'fontsize', 22);
set(gca, 'FontSize', 20, 'FontName', 'Arial'); 
colormap(jet); 
caxis([0 0.5]);
zlim([0 1]);

figure(2)
surf(delayx1, dopplerx1,DDHSb, 'FaceColor', 'interp', 'EdgeColor', 'none');
axis tight; grid on; colorbar;
xlabel('Delay', 'fontsize', 22); 
ylabel('Doppler', 'fontsize', 22); 
title('AAF of HSb', 'fontsize', 22);
set(gca, 'FontSize', 20, 'FontName', 'Arial'); 
colormap(jet); 
caxis([0 0.5]);
zlim([0 1]);

figure(3)
surf(delayx1, dopplerx1,DDWSa, 'FaceColor', 'interp', 'EdgeColor', 'none');
axis tight; grid on; colorbar;
xlabel('Delay', 'fontsize', 22); 
ylabel('Doppler', 'fontsize', 22); 
title('AAF of WSa', 'fontsize', 22);
set(gca, 'FontSize', 20, 'FontName', 'Arial'); 
colormap(jet); 
caxis([0 0.5]);
zlim([0 1]);

figure(4)
surf(delayx1, dopplerx1,DDWSb, 'FaceColor', 'interp', 'EdgeColor', 'none');
axis tight; grid on; colorbar;
xlabel('Delay', 'fontsize', 22); 
ylabel('Doppler', 'fontsize', 22); 
title('AAF of WSb', 'fontsize', 22);
set(gca, 'FontSize', 20, 'FontName', 'Arial'); 
colormap(jet); 
caxis([0 0.5]);
zlim([0 1]);

figure(5)
surf(delayx1, dopplerx1,DDHWSa, 'FaceColor', 'interp', 'EdgeColor', 'none');
axis tight; grid on; colorbar;
xlabel('Delay', 'fontsize', 22); 
ylabel('Doppler', 'fontsize', 22);
title('AAF of HWSa', 'fontsize', 22);
set(gca, 'FontSize', 20, 'FontName', 'Arial'); 
colormap(jet); 
caxis([0 1]);
zlim([0 1]);

figure(6)
surf(delayx1, dopplerx1,DDHWSb, 'FaceColor', 'interp', 'EdgeColor', 'none');
axis tight; grid on; colorbar;
xlabel('Delay', 'fontsize', 22); 
ylabel('Doppler', 'fontsize', 22); 
title('AAF of HWSb', 'fontsize', 22);
set(gca, 'FontSize', 20, 'FontName', 'Arial'); 
colormap(jet); 
caxis([0 1]);
zlim([0 1]);


function af = getaf(seq1,seq2,omegamax,omosr,taumax,tauosr,OLength) 
    u_basic1=seq1;
    u_ref=seq2;
    m_basic=length(u_basic1);
    F=omegamax;
    K=F*omosr;
    df=F/K/OLength;
    T=taumax/m_basic;
    N=taumax*tauosr;
    sr=1;
    r=ceil(sr*(N)/T/m_basic);

    if r==1
       dt=1;
       m=m_basic;
       u=u_basic1;
       u2=u_ref;
    else                               % i.e., several samples within a bit
       dt=1/r;	                       % interval between samples
       ud=diag(u_basic1);
       ud2=diag(u_ref);
       ao=ones(r,m_basic);
       ao2=[ones(1,m_basic);zeros(r-1,m_basic)];
       m=m_basic*r;
       u_basic1=reshape(ao*ud,1,m);    % u_basic with each element repeated r times
       u_ref=reshape(ao2*ud2,1,m);
       u=u_basic1;
       u2=sparse(u_ref);
    end

    t=[0:r*m_basic-1]/r;
    dtau=ceil(T*m)*dt/(N);
    tau=round([-N:1:N]*dtau/dt)*dt;
    f=[0:1:K]*df; 
    f=[-fliplr(f) f];

    % calculate ambiguity function using sparse matrix manipulations (no loops)
    mat1=spdiags(u2',0,m,m);

    % define a convolution sparse matrix based on the signal samples u1 u2 u3 ... um
    u_padded=[u(m-ceil(T*m)+1:m),u,u(1:ceil(T*m))];%periodic
    %u_padded=[zeros(1,ceil(T*m)),u,zeros(1,ceil(T*m))];%aperiodic

    % define column indexing and row indexing vectors
    cidx=[1:m];
    ridx=round(tau/dt)';

    % define indexing matrix with Nused+1 rows and m+ceil(T*m) columns 
    % where each element is the index of the correct place in the padded version of u
    index = cidx(ones(2*N+1,1),:) + ridx(:,ones(1,m))+max(ridx);

    % calculate matrix
    mat2 = sparse(u_padded(index)); 

    % calculate the ambiguity matrix for positive delays given by  
    uu_pos=mat2*mat1;
    clear mat2 mat1

    % calculate exponent matrix for full calculation of ambiguity function. the exponent
    e=exp(-1i*2*pi*f'*t);

    % calculate ambiguity function for positive delays by calculating the integral for each
    a_pos=abs(e*uu_pos');
    %a_pos=a_pos/max(max(a_pos));
    
%     % normalize ambiguity function to have a maximal value of 1
    a=flipud(a_pos);
    % exclude the zero Delay that was taken twice
    a=a([1:K (K+2):2*(K+1)],:);
    af = a;
end
