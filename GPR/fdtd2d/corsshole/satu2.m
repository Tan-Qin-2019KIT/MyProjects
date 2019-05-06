function dsa = satu2(d,sa)
sa=sa*max(max(d));
d=((d-sa)-abs(d-sa))/2+sa;
d=(abs(d+sa)+(d+sa))/2-sa;
dsa=d/sa;
return

