clear;

bc_file = 'CELL_BARCODE_COUNT_FILE';
size_file = 'NODUP_SIZE_FILE';
name = 'NAME';

% - Read full file.
fid  = fopen( bc_file, 'r' ) ;
data = textscan( fid, '%s%s%f%s%s%s%s%s%s%s%s%s%s', 'Delimiter', '\t' ) ;
fclose( fid ) ;

num_reads = load(size_file);

x = cumsum(data{3}) / num_reads;
plot(x, 'r-');
xlim([0 EST_NUM_CELL_BCS]);
ylim([0 1]);
title(name);
y = diff(x);               
i = find (y <= 1/EST_NUM_CELL_BCS / 2);    
hold on;
plot(i(1), x(i(1)), 'ko');
text(i(1), x(i(1)) - 0.1, sprintf('%d', i(1)));
hold off;
axis square;

xlabel('STAMP');
ylabel('cumulative fraction of reads');

print -dpdf OUT_FILE
