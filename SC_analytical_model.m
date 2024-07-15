clear;clc;close

%Initial parameter from data sheet of the CoolSiC FF2MR12W3M1H B11.
g = 44;
Vth = 5.0;
Cgs = 28e-9;
Vdc = 700;
Ls = 0.5e-9;
Rg = 2.0;
Rg_soft_off = 20.0;
VCC = 15;
VEE = -8;
L_sigma = 38e-9;
t_step = 1e-10;


Tsw = 50e-6; %switching period of converter
TW = 3e-6; % Clamping time

vgs_on = VEE:1e-2:(VCC-1e-3);
NUM_ON = length(vgs_on);
dids_dt_on = zeros(1,NUM_ON);
ids_on = zeros(1,NUM_ON);
t_on = zeros(1,NUM_ON);
for k = 1:1:(NUM_ON)
    if k==1
        dids_dt_on(1) = 0;
        ids_on(1) = 0;
        t_on(1) = 0;
    elseif k<NUM_ON
        if vgs_on(k)>= Vth
            ids_on(k) = g*(vgs_on(k)-Vth)^2;
            t_on(k) = -2*g*Ls*(vgs_on(k)-Vth) + (Rg*Cgs - 2*g*Ls*Vth + 2*g*Ls*VCC)*log((VCC-Vth)/(VCC-vgs_on(k))) + Rg*Cgs*log((VCC-VEE)/(VCC-Vth));
            dids_dt_on(k) = (ids_on(k+1) - ids_on(k))/(t_on(k+1)-t_on(k));
        else
            ids_on(k) = 0;
            dids_dt_on(k) = 0;        
            t_on(k) = (Rg*Cgs)*log((VCC-VEE)/(VCC-vgs_on(k)));
        end
    elseif k==NUM_ON
        ids_on(k) = g*(vgs_on(k)-Vth)^2;
        t_on(k) = -2*g*Ls*(vgs_on(k)-Vth) + (Rg*Cgs - 2*g*Ls*Vth + 2*g*Ls*VCC)*log((VCC-Vth)/(VCC-vgs_on(k))) + Rg*Cgs*log((VCC-VEE)/(VCC-Vth));
        dids_dt_on(k) = dids_dt_on(k-1);
    end

end

vgs_off = VCC:-5e-5:VEE;
NUM_OFF = length(vgs_off);
dids_dt_off = zeros(1,NUM_OFF);
ids_off = zeros(1,NUM_OFF);
t_off = zeros(1,NUM_OFF);
for k = 1:1:(NUM_OFF)
    if k==1
        ids_off(k) = g*(VCC-Vth)^2;
        t_off(k) = 0;
    elseif k<NUM_OFF
        if vgs_off(k)>= Vth
            ids_off(k) = g*(vgs_off(k)-Vth)^2;
            t_off(k) = 2*g*Ls*(VCC - vgs_off(k)) + (Rg_soft_off*Cgs - 2*g*Ls*Vth + 2*g*Ls*VEE)*log((VCC-VEE)/(vgs_off(k)-VEE));
        else
            ids_off(k) = 0;
            t_off(k) = (Rg_soft_off*Cgs)*log((Vth-VEE)/(vgs_off(k)-VEE)) + 2*g*Ls*(VCC-Vth) + (2*g*Ls*VEE + Rg_soft_off*Cgs - 2*g*Ls*Vth)*log((VCC-VEE)/(Vth-VEE));
        end
    elseif k==NUM_OFF
        ids_off(k) = 0;
        t_off(k) = t_off(k-1)+5e-9;
    end
end

for k = 1:1:(NUM_OFF)
    if k==1
        dids_dt_off(k) = 0;
    elseif k<NUM_OFF
        if vgs_off(k)>= Vth
            dids_dt_off(k) = (ids_off(k+1) - ids_off(k))/(t_off(k+1)-t_off(k));
        else
            dids_dt_off(k) = 0;        
        end
    elseif k==NUM_OFF
            dids_dt_off(k) = dids_dt_off(k-1);

    end

end

%Length before SC happens
t_0 = -2e-6:t_step:0; % before 2us.
vgs_0 = VEE*ones(1,length(t_0));
ids_0 = 0*ones(1,length(t_0));
vds_0 = Vdc*ones(1,length(t_0));


%define a new 1-D data with linear fit 
t_on_1 = min(t_on):t_step:max(t_on); 
vgs_on_1 = interp1(t_on, vgs_on,t_on_1);
ids_on_1 = interp1(t_on, ids_on,t_on_1);
dids_dt_on_1 = interp1(t_on, dids_dt_on, t_on_1);


t_off_1 = min(t_off):t_step:max(t_off); 
vgs_off_1 = interp1(t_off, vgs_off,t_off_1);
ids_off_1 = interp1(t_off, ids_off,t_off_1);
dids_dt_off_1 = interp1(t_off, dids_dt_off, t_off_1);


%insert middle gap data
t_middle = max(t_on_1):t_step:(max(t_on_1)+TW);
vgs_middle_1 = ones(1,length(t_middle))*VCC;
ids_middle_1 = ones(1,length(t_middle))*g*(VCC-Vth)^2;
dids_dt_middle_1 = zeros(1,length(t_middle));


t_1 = [t_0,t_on_1,t_middle,t_off_1 + t_on_1(NUM_ON-1)+max(t_middle)];
vgs_1 = [vgs_0,vgs_on_1,vgs_middle_1, vgs_off_1];
ids_1 = [ids_0,ids_on_1,ids_middle_1, ids_off_1];
dids_dt_1 = [dids_dt_on_1,dids_dt_middle_1,dids_dt_off_1];
vds_1 = [vds_0, Vdc - L_sigma.*dids_dt_1];

Esc = sum(ids_1.*vds_1)*t_step;
Vpk = max(vds_1);

subplot(3,1,1);
plot(t_1,vgs_1);
grid on;
subplot(3,1,2);
plot(t_1,ids_1);
grid on;
subplot(3,1,3);
plot(t_1,vds_1);
grid on;

 
 
 
 
