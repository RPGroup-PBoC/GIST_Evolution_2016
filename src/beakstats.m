%Want to organize beak data into four separate matrices

%Bird 1 is in odd columns. Bird 2 is in even columns. 
%Rows 1-10 are length. Rows 11-20 are width. 

bird1_length = beaks(1:10, 1:2:length(beaks));
bird2_length = beaks(1:10, 2:2:length(beaks));

bird1_width = beaks(11:20, 1:2:length(beaks));
bird2_width = beaks(11:20, 2:2:length(beaks));

bird1_length_var = var(bird1_length)
bird1_length_total_var = var(bird1_length(:))
%%

histogram(bird1_length(:),50);
hold on
histogram(bird1_length(:,1),5);
hold off
legend('Single Group', 'All Groups')
xlabel('Beak length [mm]')
ylabel('Frequency')

%%

boxplot(bird1_length)
xlabel('Team')
ylabel('Beak length [mm]')

%%
hold on
for i=1:length(bird1_length)
    ecdf(bird1_length(:, i))
end
ylim([0, 1.1])

%%

fano_bird1 = var(bird1_length) ./ mean(bird1_length);
fano_bird2 = var(bird2_length) ./ mean(bird2_length);

plot(fano_bird1, fano_bird2, 'o')
ylim([0, 0.1])
