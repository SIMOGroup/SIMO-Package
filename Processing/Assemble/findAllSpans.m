function idcs = findAllSpans(n, p, KntVect, NEl)

idcs = zeros(1, NEl); % span index
el = 1;
for i = p + 1 : n
    if (abs(KntVect(i) - KntVect(i + 1)) > sqrt(eps))
        idcs(el) = i;
        el = el + 1;
    end
end
end