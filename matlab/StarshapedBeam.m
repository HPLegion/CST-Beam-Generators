%% Natural Constants

mp = 1.6726219e-27; %kg Proton mass
qe = 1.60217662e-19; %C elementary charge
c = 299792458; %m/s speed of light

%% Setup
Ekin = 6 * 30000; %eV kinetic energy 
radius =20 /1000; %m radius of generated ring
num = 16; % number of particles on said ring
offsetangle = 0; %deg first particle will be at (1,0) if zero otherwise rotated
m = 12*mp; %particle mass
q = 6*qe; %particle charge
i_c = 1;%A %particle "current" for space charge runs

Ekin = Ekin*qe; %J transform to SI units
pnorm = sqrt(2*Ekin/m)/c; %beta*gamma normalised momentum CST Stlyle

%% Generate Distribution
Lm = ones(num,1)*m;
Li = ones(num,1)*i_c;
Lq = ones(num,1)*q;
Lx = zeros(num,1);
Ly = zeros(num,1);
Lz = zeros(num,1);
Lpx = zeros(num,1);
Lpy = zeros(num,1);
Lpz = ones(num,1)*pnorm;



for k=1:num
    angle = offsetangle+(k-1)*360/num;
    Lx(k) = radius*cosd(angle);
    Ly(k) = radius*sind(angle);
end

% %%Transform distribution (optional, needed this at some point )
% theta = 160;
% R = [cosd(theta), 0, sind(theta);...
%     0, 1, 0;...
%     -sind(theta), 0, cosd(theta)];
% tr = [Lx,Ly,Lz]*R';
% tr = tr + repmat([-86.26/1000, 0, 501/1000],num,1);
% trp = [Lpx,Lpy,Lpz]*R';
% if num >1
%     Lx=tr(1,:);
%     Ly=tr(2,:);
%     Lz=tr(3,:);
%     
%     Lpx=trp(1,:);
%     Lpy=trp(2,:);
%     Lpz=trp(3,:);
% else if num == 1
%         Lx=tr(1);
%         Ly=tr(2);
%         Lz=tr(3);
%         
%         Lpx=trp(1);
%         Lpy=trp(2);
%         Lpz=trp(3);
%     end
% end
%     %%TransformEnd
    

 %% save to L.txt in current path 
 % will be improved in the future   
    L=table(Lx,Ly,Lz,Lpx,Lpy,Lpz,Lm,Lq,Li);
    writetable(L,'L','Delimiter',' ');
