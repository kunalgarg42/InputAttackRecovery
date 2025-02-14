clc
clearvars -except u2 u4 attack_start0 attack_length no_attack detect_delay
close all

global k_s1 k_s2 k_s3 dt kr kt m Ix Iy Iz options g b l d ptm u_av u_M

options = optimoptions('quadprog','Display','off');

param_generator;

num = 50000;

dt = 0.001;

j11 = 0.177;
j22 = 0.177;
j33 = 0.334;

g = 9.8;
m = 4.493;
kt = 1;
kr = 1.5;

Ix = j11;
Iy = j22;
Iz = j33;

% g = 9.8;

% m = 1;

b = 1;
l = 0.1;
d = 0.0024;

x = zeros(12,num);
% x(:,1) = zeros(12,1);
x(1:2:5,1) = [0.01; 0.01; -0.1];

u = zeros(4,num);
xd = zeros(12,1);
xd(5) = -5;

k_s1 = 0;
k_s2 = 0;
k_s3 = 0;


% kt = 1;
% kr = 1;
% u4 = 0;

phid = 0;
thetad = 0;
ptm = 0.5;


u_av = 11.08;

u_M = 2.5*u_av;


% delta = 0.1;
BDm = -100; %*ones(20,10);
Bm = -100; % *ones(20,10);

no_Detect = 0; % zeros(20,10);

delta = 0.1;
delta2 = 0.01;

gamma = 5;

delta_2 = 3.437/4*ptm^2;


for k = 1:1

    %     delta_2 = ptm^2*(0.8+k/50);

    for j = 1:1
%         j
%         k
        attack = 0;
        last_attack = 0;
        attack_start = attack_start0;

        attack_signal = zeros(num,1);
        detect_signal = zeros(num,1);
        delta = 0.02*j;

        for i = 1:num
            phi = x(7,i);
            theta = x(9,i);
            psi = x(11,i);

            cp = cos(phi);
            ct = cos(theta);
            cs = cos(psi);
            sp = sin(phi);
            st = sin(theta);
            ss = sin(psi);

            F = [x(2,i);
                0-kt*x(2,i);
                x(4,i);
                0-kt*x(4,i);
                x(6,i);
                g-kt*x(6,i);
                x(8,i)+x(10,i)*sp*tan(theta)+x(12,i)*cp*tan(theta);
                (Iy-Iz)/Ix*x(10,i)*x(12,i)-kr*x(8,i);
                x(10,i)*cp-x(12,i)*sp;
                (Iz-Ix)/Iy*x(8,i)*x(12,i)-kr*x(10,i);
                1/ct*(x(12,i)*cp+x(10,i)*sp);
                (Ix-Iy)/Iz*x(8,i)*x(10,i)-kr*x(12,i);
                ];


            G1 = [0 0 0 0;
                -1/m*(sp*ss+cp*cs*st) 0 0 0;
                0 0 0 0;
                -1/m*(cp*ss*st-cs*sp) 0 0 0;
                0 0 0 0;
                -1/m*(cp*ct) 0 0 0;
                0 0 0 0;
                0 1/Ix 0 0;
                0 0 0 0;
                0 0 1/Iy 0;
                0 0 0 0;
                0 0 0 1/Iz];

            %     G2 = [b  b  b  b;
            %         -b*l b*l b*l -b*l;
            %         -b*l -b*l b*l b*l;
            %         -d   d  -d  d];


            G2 = [b  b  b  b;
                0  -b*l 0 b*l;
                -b*l 0 b*l 0;
                d   -d  d  -d];

            G = G1*G2;
            %     [u_n1, u_n2, u_n3] = input_assign(x(:,i),xd);

            %     [u_n1, u_n2, u_n3,phid,thetad] = input_assign_PID(x(:,i),xd,phid,thetad);

            if (i>last_attack+no_attack && attack == 0)
                attack = 1;

                %         if attack_start<last_attack
                if i == last_attack+no_attack+1
                    attack_start = i;
                    last_attack = attack_start+attack_length;
                end
            end

            if attack>0
                if i>=attack_start+attack_length/2
                    attack = 2;
                end
            end

            if i>attack_start+attack_length
                attack = 0;
            end

            if i>1
                B1_old = B1;
                B2_old = B2;
                B3_old = B3;
            end

            B1 = x(5,i);
            B2 = x(7,i)^2-ptm^2;
            B3 = x(9,i)^2-ptm^2;

            if max([B1, B2, B3])>Bm(j,k)
                Bm(j,k) = max([B1, B2, B3]);
            end
            if i>1
                B1d = (B1-B1_old)/dt;
                B2d = (B2-B2_old)/dt;
                B3d = (B3-B3_old)/dt;
            else
                B1d = 0;
                B2d = 0;
                B3d = 0;
            end

            if max([B1d, B2d, B3d])>BDm(j,k)
                BDm(j,k) = max([B1d, B2d, B3d]);
            end

            if attack > 0
                if  detect_signal(i-1) == 1 || ((B1d>-delta*B1-10*dt && B1>-delta_2) || (B2d>-delta*B2-10*dt && B2>-delta_2) || (B3d>-delta*B3-10*dt && B3>-delta_2)) %i<attack_start+detect_delay
                    if i == 1
                        u1 = zeros(7,1);
                    end

                    detect_signal(i) = 1;

                    %             u1 = input_assign_QP_two_motors(x(:,i),xd,u1);
                    %             u(:,i) = [u1(1);u2(i);u1(2);u4(i)];

                    u1 = input_assign_QP(x(:,i),xd,u1);
                    u(:,i) = [u1(1);u1(2);u1(3);u4(i)];
                else
                    if i == 1
                        u1 = zeros(6,1);
                    end
                    detect_signal(i) = 0;
                    u1 = input_assign_QP_no_attack(x(:,i),xd,u1);
                    u(:,i) = [u1(1:3,1);u4(i)];

                end
                if attack>1
                    if  detect_signal(i-1) == 2 || ((B1d>=-delta2*B1-10*dt-gamma && B1>=-delta_2) || (B2d>=-delta2*B2-10*dt-gamma && B2>=-delta_2) || (B3d>=-delta2*B3-10*dt-gamma && B3>=-delta_2))
                        if i == 1
                            u1 = zeros(7,1);
                        end
                        detect_signal(i) = 2;

                        u1 = input_assign_QP_two_motors(x(:,i),xd,u1);
                        u(:,i) = [u1(1);u2(i);u1(2);u4(i)];


                    else
                        if i == 1
                            u1 = zeros(6,1);
                        end
                        detect_signal(i) = 1;
                        u1 = input_assign_QP_smart_attack(x(:,i),xd,u1);
                        u(:,i) = [u1(1);u2(i);u1(3);u1(4)];
                    end
                end
            else
                if i == 1
                    u1 = zeros(7,1);
                end
                detect_signal(i) = 0;

                u1 = input_assign_QP_no_attack(x(:,i),xd,u1);
                u(:,i) = u1(1:4);
            end


            %     u1 = [u_n1;u_n2;u_n3];
            %     if sum(u1(1:3))>9
            %         i
            %     end


            attack_signal(i) = attack;
            x(:,i+1) = x(:,i) + (F+G*u(:,i))*dt;
            if x(5,i+1)>0
                break
            end
        end
        no_Detect(j,k) = sum(attack_signal-detect_signal);

        %     plot(detect_signal)
        %     hold on

    end
end
% plot(attack_signal)
%     hold on
%%
% for k = 1:10
%     figure(1)
%     plot(BDm(:,k))
%     hold on
% 
%     figure(2)
%     plot(Bm(:,k))
%     hold on
% 
%     figure(3)
%     plot(no_Detect(:,k))
%     hold on
% end

% plot(x(5,1:i),'linewidth',3)
% set(gca, 'fontsize',30)
% legend('$$z$$','interpreter','latex')
%
% figure
%
%
% plot(x(8,1:i),'linewidth',2)
% hold on
% plot(x(10,1:i),'linewidth',2)
% legend('$$\dot \phi$$','$$\dot \theta$$','interpreter','latex')
%
% figure
%
% subplot(3,1,1)
% plot(x(7,1:i),'linewidth',2)
% hold on
% plot(0,0,'linewidth',2,'color','r')
% plot(0,0,'linewidth',2,'color','k')
% set(gca, 'fontsize',30)
% legend('$$\phi$$','$$\theta$$','$$\psi$$','interpreter','latex')
% subplot(3,1,2)
% plot(x(9,1:i),'linewidth',2,'color','r')
% set(gca, 'fontsize',30)
% subplot(3,1,3)
% plot(x(11,1:i),'linewidth',2,'color','k')
% set(gca, 'fontsize',30)
%
%
% figure
% subplot(4,1,1)
% plot(u(1,1:i),'linewidth',2)
% axis([ 0 i 0 7])
% hold on
% plot(0,0,'linewidth',2,'color','r')
% plot(0,0,'linewidth',2,'color','k')
% plot(0,0,'linewidth',2,'color','g')
% set(gca, 'fontsize',30)
% legend('$$u_1$$','$$u_2$$','$$u_3$$','$$u_4$$','interpreter','latex')
% subplot(4,1,2)
% plot(u(2,1:i),'linewidth',2,'color','r')
% axis([ 0 i 0 7])
% set(gca, 'fontsize',30)
% % legend('$$u_2$$','interpreter','latex')
% subplot(4,1,3)
% plot(u(3,1:i),'linewidth',2,'color','k')
% axis([ 0 i 0 7])
% set(gca, 'fontsize',30)
% % legend('$$u_3$$','interpreter','latex')
% subplot(4,1,4)
% plot(u(4,1:i),'linewidth',2,'color','g')
% set(gca, 'fontsize',30)
% axis([ 0 i 0 7])
%
% % figure
% % plot(attack_signal(1:i),'linewidth',2)
% % hold on
% % plot(detect_signal(1:i),'linewidth',2)
% % legend('$$u_4$$','interpreter','latex')
% x_data = x;
%%

%
% x = x_data(:,1:50:i);
% lambda = 0:0.01:2*pi;
% write_vid = 1;
%
% y = x(3,1:end-1);
% z = -x(5,1:end-1);
% phi = x(7,1:end-1);
% theta = x(9,1:end-2);
% psi = x(11,1:end-1);
% x = x(1,1:end-1);
%
% animatesimulation
%
% x = x_data(:,1:i);