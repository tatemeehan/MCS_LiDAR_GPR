mask = isnan(C);
tmpdiff = gradient(mask);
binix = find(abs(tmpdiff)>0);
tmp = zeros(size(C));
tmp(binix) = 1;
tmp2 = imtranslate(tmp,[1,1]);
tmp3 = imtranslate(tmp,[-1,-1]);
tmp4 = imtranslate(tmp,[1,-1]);
tmp5 = imtranslate(tmp,[-1,1]);
tmp6 = imtranslate(tmp,[2,2]);
tmp7 = imtranslate(tmp,[-2,-2]);
tmp8 = imtranslate(tmp,[2,-2]);
tmp9 = imtranslate(tmp,[-2,2]);
tmp10 = imtranslate(tmp,[3,3]);
tmp11 = imtranslate(tmp,[-3,-3]);
tmp12 = imtranslate(tmp,[3,-3]);
tmp13 = imtranslate(tmp,[-3,3]);
tmp14 = imtranslate(tmp,[1,0]);
tmp15 = imtranslate(tmp,[-1,0]);
tmp16 = imtranslate(tmp,[0,1]);
tmp17 = imtranslate(tmp,[0,-1]);

wtf = tmp+tmp2+tmp3+tmp4+tmp5+tmp6+tmp7+tmp8+tmp9+tmp10+tmp11+tmp12+tmp13+tmp14+tmp15+tmp16+tmp17;
wtf2 = imtranslate(wtf,[0 -1]);
wtf3 = imtranslate(wtf,[1,0]);
wtf4 = imtranslate(wtf,[-1,0]);
wtf5 = imtranslate(wtf,[0,1]);


wwtf = wtf+wtf2+wtf3+wtf4+wtf5;