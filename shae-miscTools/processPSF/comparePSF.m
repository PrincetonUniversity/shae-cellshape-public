function statsOut = comparePSF(psf1,psf2);

% measures that expect peak of unity
psf1 = psf1./max(psf1(:));
psf2 = psf2./max(psf2(:));

% difference in total intensity
statsOut.intUnitPeak1 = sum(psf1(:));
statsOut.intUnitPeak2 = sum(psf2(:));

% difference in volume above 0.5
statsOut.volPeak1_50 = nnz(psf1(:)>0.5);
statsOut.volPeak2_50 = nnz(psf2(:)>0.5);

% difference in volume above 0.05
statsOut.volPeak1_5 = nnz(psf1(:)>0.05);
statsOut.volPeak2_5 = nnz(psf2(:)>0.05);

% root mean squared difference
statsOut.rmsd  = sqrt(mean((psf1(:)-psf2(:)).^2));

% compare just the region that is above 5% in either psf
centralVox = or(psf1(:)>0.05,psf2(:)>0.05);
statsOut.rmsdCentral = sqrt(mean((psf1(centralVox)-psf2(centralVox)).^2));

% get ready for KL divergence
P = psf1(:)./sum(psf1(:));
Q = psf2(:)./sum(psf2(:));
statsOut.rmsdScale = sqrt(mean((P(:)-Q(:)).^2));

% truncate less than zero
P2 = P;
P2(P2(:)<0) = 0;
P2 = P2./(sum(P2(:)));
Q2 = Q;
Q2(Q2(:)<0) = 0;
Q2 = Q2./(sum(Q2(:)));

% invert less than zero
P3 = abs(P);
P3 = P3./sum(P3(:));
Q3 = abs(Q);
Q3 = Q3./sum(Q3(:));

% KL divergence of truncaated
KLTempPQ2 = P2./Q2;
KLTempPQ2(not(isfinite(KLTempPQ2)))=1;
KLTempPQ2(KLTempPQ2==0) = 1;

KLTempQP2 = Q2./P2;
KLTempQP2(not(isfinite(KLTempQP2)))=1;
KLTempQP2(KLTempQP2==0) = 1;


statsOut.KL_PQ_2 = sum(log(KLTempPQ2(:))./log(2).*P2(:));
statsOut.KL_QP_2 = sum(log(KLTempQP2(:))./log(2).*Q2(:));

% KL divergence of abs
KLTempPQ3 = P3./Q3;
KLTempPQ3(not(isfinite(KLTempPQ3)))=1;
KLTempPQ3(KLTempPQ3==0) = 1;


KLTempQP3 = Q3./P3;
KLTempQP3(not(isfinite(KLTempQP3)))=1;
KLTempPQ3(KLTempQP3==0) = 1;

statsOut.KL_PQ_3 = sum(log(KLTempPQ3(:))./log(2).*P3(:));
statsOut.KL_QP_3 = sum(log(KLTempQP3(:))./log(2).*Q3(:));

% noise floor? something about resampling as lower bit depth image and when
% that changes

for iBitDepth = 1:24
   tempMax = 2^iBitDepth-1;
   PDepth = round(tempMax.*P2)/tempMax;
   PDepth = PDepth./sum(PDepth(:));
   KLTempP = P3./PDepth;
   KLTempP(not(isfinite(KLTempP)))=1;
   statsOut.KL_PDepth(iBitDepth) = sum(log(KLTempP(:))./log(2).*P3(:));
 
   QDepth = round(tempMax.*Q2)/tempMax;
   QDepth = QDepth./sum(QDepth(:));
   KLTempQ = Q3./QDepth;
   KLTempQ(not(isfinite(KLTempQ)))=1;
   statsOut.KL_QDepth(iBitDepth) = sum(log(KLTempQ(:))./log(2).*Q3(:));
end
P2(P2==0) = 1;
Q2(Q2==0) = 1;
statsOut.PSelfInfo = -1*sum(log(P2(:))./log(2).*P2(:));
statsOut.QSelfInfo = -1*sum(log(Q2(:))./log(2).*Q2(:));






