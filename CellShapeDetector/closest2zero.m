function C=closest2zero(A,B)
A(isnan(A))=Inf;
B(isnan(B))=Inf;
Csign=sign(A+B).*sign(A.*B);
C=Csign.*min(abs(A),abs(B));