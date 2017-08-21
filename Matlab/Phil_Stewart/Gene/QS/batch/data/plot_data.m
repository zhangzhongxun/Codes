%plot batch QS data
figure(2)

hold off;

subplot(2,3,1);

plot(t_int,S,'-o');
xlabel('Time (hours)','fontsize',14);
ylabel('Oxygen','fontsize',14)

subplot(2,3,2);

semilogy(t_int,X,'-o');
xlabel('Time (hours)','fontsize',14);
ylabel('Cell density (log)','fontsize',14)

subplot(2,3,3);

plot(t_int,QS,'-o');
xlabel('Time (hours)','fontsize',14);
ylabel('Quorum Sensing Molecule','fontsize',14)

subplot(2,3,4);

plot(t_int,WL,'-o');
xlabel('Time (hours)','fontsize',14);
ylabel('lasI','fontsize',14)

subplot(2,3,5);

plot(t_int,WR,'-o');
xlabel('Time (hours)','fontsize',14);
ylabel('rsaL','fontsize',14)

fid = fopen('QS_batch.data','w');

fprintf(fid,'%12.10e ', S);
fprintf(fid,'\n');

fprintf(fid,'%12.10e ', X);
fprintf(fid,'\n');

fprintf(fid,'%12.10e ', QS);
fprintf(fid,'\n');

fprintf(fid,'%12.10e ', WL);
fprintf(fid,'\n');

fprintf(fid,'%12.10e ', WR);
fprintf(fid,'\n');

fclose(fid);