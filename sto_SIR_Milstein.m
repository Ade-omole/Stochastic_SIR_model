randn('state',100)

lambda = 0.2; beta = 0.4; epsilon = 0.1; mu = 0.2; 
gamma = 0.2; sigma_1 = 0.01; sigma_2 = 0.02; sigma_3 = 0.01; Szero = 0.7;
Izero = 0.2; Rzero = 0.1; 
T = 100; N = 2^8; dt = T/N;

dW = sqrt(dt)*randn(1,N);       
R = 4; Dt = R*dt; L = N/R;        
S = zeros(1,L); I = zeros(1,L); R = zeros(1,L); 

S(1) = Szero; I(1) = Izero; R(1) = Rzero;

for j = 1:L-1

    Winc = sum(dW(R*(j-1)+1:R*j));

    S(j+1) = S(j) + (lambda-beta*S(j).*I(j)-mu*S(j))*Dt +...
            sigma_1*S(j).*Winc + 0.5*sigma_1^2*S(j)*(Winc.^2 - Dt);

    I(j+1) = I(j) + (beta*I(j)-(mu+epsilon+gamma)*I(j))*Dt...
            + sigma_2*I(j).*Winc+ 0.5*sigma_2^2*I(j).*(Winc.^2 - Dt);

    R(j+1) = R(j) + (gamma*I(j)-mu*R(j))*Dt + sigma_3*R(j).*Winc...
            + 0.5*sigma_3^2*R(j).*(Winc.^2 - Dt);
end

figure(1)
plot([0:Dt:T], [Szero,S],'b-*');
xlabel('$t$','Interpreter','latex','FontSize',14 );
ylabel('$S(t)$','Interpreter','latex','FontSize',14 );
title('Susceptible','Interpreter','latex','FontSize',16);

figure(2)
plot([0:Dt:T], [Izero,I],'r--*');
xlabel('$t$','Interpreter','latex','FontSize',14 );
ylabel('$I(t)$','Interpreter','latex','FontSize',14 );
title('Infected','Interpreter','latex','FontSize',16);

figure(3)
plot([0:Dt:T], [Rzero,R],'g--*');
xlabel('$t$','Interpreter','latex','FontSize',14 );
ylabel('$R(t)$','Interpreter','latex','FontSize',14 );
title('Removed','Interpreter','latex','FontSize',16);
