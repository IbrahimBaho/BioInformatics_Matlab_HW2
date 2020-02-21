


[ ma60, mis60] = calc_matrecies (0.6);
[ ma70, mis70] = calc_matrecies (0.7);
[ ma80, mis80] = calc_matrecies (0.8);
[ ma90, mis90] = calc_matrecies (0.9);

array60 = [ma60,mis60,mis60,mis60 ;mis60 ,ma60, mis60, mis60;...
    mis60,mis60, ma60, mis60; mis60,mis60 ,mis60, ma60];

array70 = [ma70,mis70,mis70,mis70 ;mis70 ,ma70, mis70, mis70;...
    mis70,mis70, ma70, mis70; mis70,mis70 ,mis70, ma70];

array80 = [ma80,mis80,mis80,mis80 ;mis80 ,ma80, mis80, mis80;...
    mis80,mis80, ma80, mis80; mis80,mis80 ,mis80, ma80];

array90 = [ma90,mis90,mis90,mis90 ;mis90 ,ma90, mis90, mis90;...
    mis90, mis90, ma90, mis90; mis90,mis90 ,mis90, ma90];

rowNames = {'A', 'C','G','T'};
colNames = {'A', 'C','G','T'};
sTable60 = array2table(array60,'RowNames',...
    rowNames,'VariableNames',colNames);
disp(sTable60)
sTable70 = array2table(array70,'RowNames',...
    rowNames,'VariableNames',colNames);
disp(sTable70)
sTable80 = array2table(array80,'RowNames',...
    rowNames,'VariableNames',colNames);
disp(sTable80)
sTable90 = array2table(array90,'RowNames',...
    rowNames,'VariableNames',colNames);
disp(sTable90)

% (B) calculate E-value for score of 20
E_value = calc_E_value(20 , 100, 100);
fprintf("for score = 20 the E-value: \n");
disp( E_value);

% (C) calculate E-value for score of 5
E_value = calc_E_value(5 , 100, 100);
fprintf("for score = 5 the E-value: \n");
disp( E_value);

% (D) calculate E-value against a database of 1mil sequences of length 100bp
DB_E_value = multi_testing_E_value( 20, 1000000, 100);
fprintf("for score = 20 against a DB of size 1mil, E-value: \n");
disp(DB_E_value);

Scores  = zeros(4,4);
k = 1;
% (E) Smith-Waterman algorithm to align seqences
seqs = fastaread("hw2q2.fa");
for i = 1:2:length(seqs)
       seq1 = seqs(i).Sequence;
       seq2 = seqs(i+1).Sequence;
       % disp([seqs(i).Header, ' and  ', seqs(i+1).Header]);
       % perform alignment between seq1 and seq2 and output aligment results
       
       [Score6, Alignment] = swalign(seq1 ,seq2, 'alphabet','nt', 'ScoringMatrix', array60, 'GapOpen', 5);
       %fprintf("The score for %s and %s is:\n", seqs(i).Header,seqs(i+1).Header)
       %disp(Score6)
       Scores( k,1) = Score6;
       
       [Score7, Alignment] = swalign(seq1 ,seq2, 'alphabet','nt', 'ScoringMatrix', array70, 'GapOpen', 5);
       %fprintf("The score for %s and %s is:\n", seqs(i).Header,seqs(i+1).Header)
       %disp(Score7)
       Scores( k,2) = Score7;
       
       [Score8, Alignment] = swalign(seq1 ,seq2, 'alphabet','nt', 'ScoringMatrix', array80, 'GapOpen', 5);
       %fprintf("The score for %s and %s is:\n", seqs(i).Header,seqs(i+1).Header)
       %disp(Score8)
       Scores( k,3) = Score8;
       
       [Score9, Alignment] = swalign(seq1 ,seq2, 'alphabet','nt', 'ScoringMatrix', array90, 'GapOpen', 5);
       %fprintf("The score for %s and %s is:\n", seqs(i).Header,seqs(i+1).Header)
       %disp(Score9)
       Scores( k,4) = Score9;
       
       k = k + 1;
end
columnsNames = {'60%', '70%','80%','90%'};
scoreTable = array2table(Scores,'VariableNames',columnsNames);
disp("The Score Matrix: ")
disp(scoreTable)

P_vals = arrayfun( @(x) P_val_func(x), Scores);
E_vals = arrayfun( @(x) calc_E_value(x, 500, 500), Scores);
disp(P_vals);

figure(1)
x = [60 , 70, 80, 90];
plot(x, Scores', '-o')
%grid
legend({'seq %ident60','seq %ident70','seq %ident80','seq %ident90'},...
    'Location','northwest');
xlabel('Scoring Matrix %ident')
ylabel('Allignment Score')

figure(2)
semilogy(x, E_vals', '-o')
xlabel('Scoring Matrix %ident')
ylabel('E-value')

function P_value = P_val_func( Score)
    temp = calc_E_value(Score, 500,500);
    P_value = 1 - exp ( -1 * temp);
end

function [match, mismatch] = calc_matrecies (percentage )
    qLL = percentage * 0.25;
    rest = (1 - percentage)/12;
    match = log (qLL / (0.25*0.25));
    mismatch = log (rest / (0.25*0.25));
end

function E_value = calc_E_value (S, M, N)
    % E(S) = K*M*N exp{- lambda * S}
    E_value = 0.1 * N * M * exp(-1 * S);
end

function DB_E_value = multi_testing_E_value( S, DB_size, seq_len)
    % E(S) = K * seq_len * (DB_size* seq_length) * exp(-lambda * S)
    DB_E_value = 0.1 * seq_len * seq_len * DB_size * exp(-1* S);
end

