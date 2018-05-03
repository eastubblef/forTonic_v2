function [sink] = TNC_ExtTrigWinsEAS(source,indices,window)

numEvents   = round(numel(indices));
range       = -window(1,1):window(1,2);
numRange    = round(numel(range));

if indices(numEvents)+window(1,2) > numel(source)
    numEvents = numEvents-1;
end
if numel(indices) > numEvents      %EAS added these lines for 1.18.17 dataset
    sink.wins = NaN(numEvents,numRange);
    sink.range = range;
    return
end

sink.wins = zeros(numEvents,numRange);
sink.range = range;

for i=1:numEvents
        sink.wins(i,:) = source(indices(i)-window(1,1):indices(i)+window(1,2));
end

sink.avg = mean(sink.wins);
sink.err = std(sink.wins) ./ sqrt(numEvents-1);