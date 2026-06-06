v = ver;
for i = 1:length(v)
    fprintf('%s | %s\n', v(i).Name, v(i).Version);
end
