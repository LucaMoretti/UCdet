CHPin=sum(costvals{1}./Fuelsall{2,2})+sum(costvals{2}./Fuelsall{1,2})
CHPout=sum(Pmat{1,3}(:))
CHPetaelmean=CHPout/CHPin

COMPRchillin=sum(Pmat{1,13})

CHPthout=sum(sum(Pmat{2,3}(1:2,:)))
CHPetathmean=CHPthout/CHPin

LTGBout=sum(Pmat{2,3}(3,:))

ABSchillin=sum(Pmat{2,13}(Dall{2}(3,:)>0))

COMPRchillout=sum(Pmat{3,3}(2,:))
ABSchillout=Pmat{3,3}(1,:);
meanchillout=mean(ABSchillout(and(ABSchillout>1,Dall{2}(3,:)>0)))
ABSchillout=sum(ABSchillout(Dall{2}(3,:)>0))


COPcompr=COMPRchillout/COMPRchillin
COPabs=ABSchillout/ABSchillin

ABSprof=Pmat{3,3}(1,:);
 mean(ABSprof(and(ABSprof>1,Dall{2}(3,:)>0)))