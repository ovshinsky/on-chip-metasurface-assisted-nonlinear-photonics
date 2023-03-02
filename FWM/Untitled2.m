N = 500;
a = zeros(1,N);
b = zeros(1,N);
c = zeros(1,N);
x_axis = zeros(1,N);
beta_max = 1e4;

for k = 1:N
    k
[a(k),b(k),c(k)] = FWM_coupled_equations_v3(beta_max/N*(k-N/2));
x_axis(k) = beta_max/N*(k-N/2);
end
figure;
hold on;
plot(x_axis,b);
legend("P {idler}")
xlabel("delta beta (m^{-1})")
ylabel("Power {idler}(mW)")
