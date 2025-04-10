library(ggplot2)
library(data.table)
library(UpSetR)
library(pheatmap)
library(Seurat)

cicero = read.csv("../../Oligo_NN_cicero_0.1.csv",header =F)

nrow(cicero)

cicero_gr = GRanges(
  seqnames = Rle(sapply(strsplit(as.character(cicero$V2), "_"), `[`, 1)),
  ranges = IRanges(start = as.numeric(sapply(strsplit(as.character(cicero$V2), "_"), `[`, 2)), 
                   end = as.numeric(sapply(strsplit(as.character(cicero$V2), "_"), `[`, 3))), fnames = cicero$V2,
    tnames = cicero$V3, score = cicero$V4
)

cicero_gr

DARs = read.csv("../combined_diff/diff_csvs/diff_peaks_Oligo_NN_2vs18_iter.csv")

diffs = DARs[which(DARs$adjusted.p.value<0.05),]
diffs$chr = sapply(strsplit(as.character(diffs$feature.name), ":"), `[`, 1)
diffs$loc = sapply(strsplit(as.character(diffs$feature.name), ":"), `[`, 2)
head(diffs)

dar_gr = GRanges(
  seqnames = Rle(diffs$chr),
  ranges = IRanges(start = as.numeric(sapply(strsplit(as.character(diffs$loc), "-"), `[`, 1)), 
                   end = as.numeric(sapply(strsplit(as.character(diffs$loc), "-"), `[`, 2))), 
    names = diffs$feature.name, fc = -diffs$log2.fold_change.)

dar_gr

gtf = read.table("../../mm10.gencode.vM16.gene.bed")

rownames(gtf) = make.unique(gtf$V4)

gtf$V2 = gtf$V2-2000
gtf$V3 = gtf$V3+2000


genes2 = c('Ptgds','Plp1','Mal','Cd81','Cnp','Sbf1','Mag','Gsn','Cldn11','Tuba1a','Pcyt2','Olig1','Tmsb4x','Kcna1','Srsf2','Cryab','Scd2','Sirt2','Pllp','Fth1','Rps8','Vmp1','Nisch','Reep5','Slc38a2','Hspa8','AY036118','Tubb4a','Cyld','Tsc22d4','Rpl17','Igsf8','Ubb','Josd2','Shisa4','Actb','Eid1','Tspan3','Rps16','Hsp90ab1','Ftl1','Lgi3','Iffo1','Ubc','Cd9','Cers2','Glt8d1','Rpl41','Taldo1','Mast3','Gabbr1','Gas5','Msmo1','Itm2b','Tjap1','Adssl1','Gal3st1','Hsp90aa1','Fabp5','Rps21','Puf60','Atp6v0c','Cox4i1','Ptma','Eef1a1','Tecr','Rplp1','Slc22a17','Klf9','Rpl35a','Tpt1','Map7d1','Bcar1','Cox8a','Rps2','H2afj','Gnas','Wscd1','Tmbim1','Tmem59','Ss18l2','Npc2','Usp16','Cpt1c','Vcp','Rpl13','Cmtm5','0610012G03Rik','Ssbp1','Rtn1','Snrpn','Litaf','Rps5','Hcn2','Oaz1','Psap','Degs1','Gde1','S100a16','Arpc1a','Ppib','Lgmn','Rpl9','Plpbp','Rps23','Cab39l','Pabpc1','Rnasek','Rpl11','Glul','Rpl24','Aplp1','Trf','Apod','Lamp1','Ckb','Stmn4','Mog','Tmbim6','Gapdh','Rhog','Opalin','Evi2a','Rpl18','Spcs1','Gamt','Abca2','Itih3','Eif1','Ctsl','Ctsb','Gjb1','Tmem50a','Dbndd2','Emc7','Fau','Bsg','Dbi','Emc10','Rpsa','Rpl28','Rabac1','Dad1','Prnp','Pcsk1n','Ppia','Aldoa','Rps29','Saraf','Tkt','Selenos','Stmn1','Serf2','Gstm5','Cox6a1','Fdps','Atp6v0b','Rps3','Ctsa','Slc25a4','Atraid','Cfl1','Selenok','Arl6ip1','Elovl1','Klk6','Cd63','Actg1','Rpl8','Rpl6','Atp6ap1','Gpr62','Ndufs7','Serpinb1a','Bcap31','Rpl27a','Rpl29','Sept4','Rpl15','Tm2d2','Rnf130','Ndufa13','Gstp1','Prdx1','Calr','Rpl7a','Sec11c','Gnb2','Pdlim2','Rpl38','Serinc1','Desi1','Rps15','Rps26','Rpl37','Elob','Fkbp1a','Cst3','Ndufb11','1700047M11Rik','Slc48a1','Cyc1','Rnf7','Tmem33','Sec62','Rpl18a','Rhob','Pebp1','Omg','Gabarap','Tmem9b','Ppp1r14a','Atp1a1','Cox5b','Gpx4','Atp5b','Cisd1','Sys1','Calm1','Rpl34','Npdc1','Ndufa3','Vkorc1','Ndufc2','Pdgfa','Ddrgk1','Smdt1','Tmbim4','Rpl19','Phgdh','Pabpn1','Edf1','Rps27','Tmem258','Faim2','Rpl14','Ap2a1','Rplp0','1810037I17Rik','Tmem125','Rps10','Myl6','Pcbp1','Rpl7','Tle5','mt-Co2','Ilvbl','Btf3','Slc35a2','Prpf4b','Itm2c','Rps28','Qdpr','Hist3h2ba','Sccpdh','Rps27a','Rps17','Rplp2','Dctn1','Gjc2','Mrfap1','BC031181','Psenen','Tmed10','Selenof','Rpl21','Sh3bgrl3','Arf5','Cdh2','Ndufb8','Unc50','Hsbp1','Cfl2','Polr3e','Rps9','M6pr','Eef2','Rpl36','Efnb3','Ap2s1','Rpl4','Timp2','Gnai2','Slc25a3','Rpl23','Mlf2','Rpl30','H3f3b','Apoe','Rps11','Mif','Nkx6-2','Fis1','Clcn2','Gpr37','Rps25','Tomm20','Cox7c','Mcam','Rpl3','Tmem88b','H2afz','Psmb6','Rps3a1','Ndufa4','Ypel3','Rpl37a','mt-Atp6','Pmp22','Pim3','Cox5a','Rps4x','Jund','Selenop','Rps6','Cdk5','Abcb8','Rac1','Gabarapl2','Erp29','Rps7','Snrnp70','Rpl23a','Pigyl','Rpl10','Ndrg1','Commd3','Imp3','Dalrd3','Rack1','Hint1','Arl2','Olig2','Rpl10a','Golga7','Atp5j','Son','Naxe','Hmgb1','Hps4','Rpl5','Chchd2','Canx','Mt1','Cox6c','Tmem151a','Cox7b','Ndufb3','Rps12','2810468N07Rik','Gm48512','mt-Co3','Rpl26','Tsn','Tppp3','Gm11659','Pgp','Sar1b','Car14','Rps14','Mapk8ip1','Gm9895','Eif4a1','Sox10','Ndufb7','Pcbp4','Atp5f1','Mtch1','Cox6b1','Tmem189','Psmb3','Atp5h','Clcn4','Atp5a1','Sc5d','Rpl27','Jkamp','Cope','Cnot3','Gltp','Csrp1','Sdf4','H3f3a','Spcs2','Cops5','Tnni1','Uqcrq','Dvl1','Rpl31','Psma7','Llgl1','Ddx5','Rps13','Cyhr1','Map1lc3b','Ndufs5','Tma16','Psmc1','Vbp1','Tmed9','Tspan15','Rab4b','Erg28','Psmb1','Slc4a2','Tcp1','Atg9a','Hsd17b12','Rps15a','Ssr3','Clptm1','Gdi1','Eif3c','Pnkd','Ktn1','Tmed7','Sdhb','Selenom','Rpl13a','Rpl32','Selenow','Rps24','3110021N24Rik','Tmem19','Tmem223','Sft2d1','AC149090.1','Carhsp1','Glod4','Rps20','Clpp','Scrn1','Stard3','Sypl','Cct6a','Tmem229a','Trir','Mdh1','Rrbp1','Hsp90b1','Psmb7','Mtch2','Kmt5a','Clta','Psmb2','Ndufa7','Ndufaf4','Alkbh1','Uqcrh','Rhoa','Cdr2l','mt-Co1','Tob2','Fnta','Pip4p1','Ybx1','Adamts4','Rps19','S100b','Npc1','Tmed5','Atp5g2','Set','Pcp4','Atp6v0e','Anapc11','Acaa1a','Nmral1','Idh3b','Atp6v1e1','Gng11','Stambp','Agpat3','4930581F22Rik','Atp1b3','Mrpl41','mt-Nd1','Ermn','Pon2','Agpat4','Hist1h2bc','B230219D22Rik','Rpl12','Fscn1','Bcat1','Tubb3','Prrt1','Psmd7','Naxd','Ndufb9','Hapln2','Hnrnph2','Ddx39b','Tmed2','Plod1','Ankra2','Mpc2','Fuca1','Brk1','Pkm','4931406P16Rik','Mbp','9330159F19Rik','Zfp445','Atp11a','Smad7','Hepacam','Metap1d','Acot7','Dusp7','Acox1','Nudcd3','Kctd3','Gjc3','Snx32','Ppp1r10','Tmem165','Dusp11','Cd82','Paip2','Map1a','Dpm1','Aifm1','Mrpl48','Snhg20','Las1l','Dusp15','Mtif2','Zfp281','Ino80dos','Ppt1','Plekhb1','Dynlrb1','Srp14','Ndufb5','Atp5g1','Atp5d','Ndufa11','Atp5o','Rnpc3','Cox7a2','Atp6v1f','Uchl1','Atp5j2','Plaat3','Pea15a','Tprn','Ndufb10','Ndufa10','Ndufa12','Pfdn5','Naca','Timp3','Swi5','Atp5l','Ift43','Dynll1','A430106G13Rik','Nsdhl','2900097C17Rik','Preb','Peg3','Fn3k','Rnf25','Smco3','Snx24','Arhgef1','Snw1','Tmem106b','Zfp706','Ywhab','Tmem41b','Pkig')

genes1 = c('Man1a','Synpr','Pcdh7','Nckap5','Fam13c','Anks1b','Kndc1','Kank1','Plekha1','Psd3','Auts2','Gm4221','Csmd2','Scn8a','Tjp1','Tmod2','Diaph3','Add1','Phyhipl','Luzp2','Cdh19','Large1','Negr1','Cntn1','9630013A20Rik','Erbb4','Klf13','Tspan2','Asxl2','Tnik','Phldb1','Ppfia2','Ust','A230057D06Rik','Aatk','Rims2','Tmem163','Syt11','Chsy3','Serinc5','Syne1','Pde7a','Mobp','Pcdh17','Slk','Mical3','Mgat5','Pcdh11x','Lcorl','Dennd5a','Zfp608','Usp6nl','Cnksr3','Jam3','Ankrd44','Ptn','Srpk2','Dscam','Wasf1','Sh3gl3','Gm47283','Sema5a','Znrf1','Mtmr2','Cmss1','Ablim2','Ubl3','Fam214a','Far1','Tsc22d1','9530059O14Rik','Nudt4','A330076H08Rik','Mitf','Acsl3','Clmn','Cdk14','Smurf1','Cldn14','Fry','Sytl2','Gramd3','Gng12','R3hdm1','Gas7','Dcc','Pik3c2a','Gnb4','Gpr137c','Nalcn','Dusp10','Opcml','Ezh2','Ank3','Rapgef2','Gkap1','Astn2','Elf1','Tmod1','Brinp3','Plcb1','Strn','Rdx','Sema4d','Fyn','Ptpdc1','Map2','Cux2','Fchsd2','Pde4d','Etv6','Vps37a','Trim24','Pik3r3','Zfp827','Pdgfc','Itpr2','Rapgef1','Zdhhc9','Cask','Map7d2','Sptlc2','Fhit','Pdk3','Tns3','Gm42418','Sytl5','Gm5087','Ptpn13','Cyp51','Jph4','Hmgcs1','Lima1','Tmem30a','Kcna2','Nrg1','Snap91','Sh3rf1','Dock4','Plpp3','Chn2','Plxnb3','Ppp1r16b','Gria2','Reep3','Bmp1','Atat1','Ccnjl','Srpk3','Cadm4','Tulp4','Frmd4a','Tusc3','Lpar1','Brd9','Ctr9','Arhgap24','Hmgxb3','Csmd1','Tsc22d3','Nrd1','Pcdh10','Phactr3','Sv2a','Abca8b','Slc6a9','Acss2','Maml2','Ipo13','Fdft1','Prickle1','Supt5','Bcas1','Itpk1','Klhl5','2610316D01Rik','Bace1','Sorl1','Pcdh9','Foxn3','Tmeff2','Daam1','Zfp536','Dock10','Zbtb16','Nfasc','Ralgps1','Pak1','Myo18a','Exoc6b','Ube2e2','Acap2','Gm37494','Tafa5','Slc44a1','Rere','Rabgap1l','Pak7','Lars2','Bcl2l1','Rbpj','Srcin1','Sec14l5','Bin1','Nrxn1','Tcf4','Plekhh1','Ptpra','Rictor','Mink1','Jmjd1c','Fam102a','5031439G07Rik','Txndc16','Zmynd8','Snx29','Kazn','Prkcz','Rffl','E130308A19Rik','Stk39','Fnip2','Pum1','Abl1','Ldlrad4','Wdr20','Zzz3','Pnpla7','Arhgap35','Sept7','Ints6','Nfix','Evi5','Med13l','Soga1','Tbc1d16','Cpm','Sema3d','Fgf12','Spsb1','March8','Tecpr2','Arpp21','Sppl3','Csnk1g3','Jazf1','Camk1d','Mmp16','Tmem108','Pak3','Creb3l2','E2f3','Ntm','Gpbp1l1','Ptpro','Nedd4l','Lrrn1','Rasgef1b','Zmiz1','Gpc6','Zfhx2os','Cacna2d4','Ubac2','Mpped2','Rab6a','Marcks','Rassf8','Ptprs','Cdyl','Map1b','Gm20754','Zfp62','Mgat3','Ttc13','Dclk1','Prkg1','Rora','Nrxn2','Slc12a2','AC163347.1','Ttll7','Nkain1','Ncor2','Nectin3','Fam228b','Hook3','Tpd52','Eif4g2','Galnt2','Mad1l1','Arid1a','Cept1','Unc5b','Hs3st1','Enpp6','Trip12','Ppp1r9a','Lrp12','Stx6','Adcy5','Zer1','Adcy9','Stam','Prex1','Nectin1','Rab7','Tfeb','Ski','Tnr','Lhfpl3','Ppfibp1','Sox6','Tcf7l2','Abtb2','Cdh13','Mgll','Trio','Gna12','Srek1','Hp1bp3','Yeats2','Usp25','Tet3','Ccser2','Gm19951','Pdxk','Arhgef10','Adamts2','Gm15564','Tspan9','Hebp1','Gpc5','Slc15a2','Kit','Tenm4','Hdgfl2','Sgip1','Cyth1','Kcnj10','Hdac5','Aff3','Mb21d2','Nxph1','2010315B03Rik','Zfyve1','Spata13','Stimate','Gfra1','Fgd3','Paqr8','Gpr155','Tmlhe','Padi2','Rab6b','Nceh1','Prepl','Ptpre','Capzb')

genes  = c('Robo1','Neat1','Dpyd','Spock3','9330199G10Rik','9330111N05Rik','March1','Csmd3','Mdga2','Nav3','Lrrc4c','Kcnq3','Snca','Sorcs1','7630403G23Rik','Eya4','Grm3','Edil3','A330015K06Rik','Zbtb20','4930420G21Rik','A230001M10Rik','Lrp1b','Cacna2d3','Sox2ot','Pam','Dgki','Ctnna2','Il33','Kcnab1','Lrmda','Nbas','Adamtsl1','Sema6d','Slc4a10','Rtn4','Plxdc2','B3galt1','App','Schip1','Nf1','Ccser1','Akap6','Pkp4','Gnao1','Spock1','Chd2','Myo1d','Acaca','Mir100hg','Ano4','Stag1','Phactr1','Foxp1','Bmp2k','Tenm2','Rasgrp3','Dscaml1','Cog5','Ncoa2','Crppa','Neo1','Enox1','Snd1','Gnaq','Dennd1a','Cdh10','Cdh20','Shtn1','Ptprm','Immp2l','Otud7a','Slain1','Slc8a3','Adgrb3','Apbb2','Grid2','Zfp532','Ppp3ca','Fmn1','Thsd7a','Wwox','Kif13a','Adgrl3','Dock1','Slc8a1','Rcan2','Acyp2','Pbx1','Herc4','Iqsec1','Sik2','Igf1r','Fggy','Usp31','Gm36855','Cadm1','Gnai1','Setbp1','Daam2','Tanc2','Rapgef6','Cobll1','A330049N07Rik','Numb','Creb1','Fnip1','Ralyl','Meis1','Fut8','Gmds','Akt3','Aff2','Smyd3','Fam13b','Eml6','Ppp1r21','Aig1','Tmem132d','Add3','Mapre2','Lncpint','Alcam','Prkca','Mllt3','Slc20a2','Pip4k2a','Ephb1','4732471J01Rik','Supt3','Slc16a7','Nptn','Trps1','Gm16759','Sgce','Glcci1','Magi1','Prkcb','Hibadh','Fmnl2','Pik3r1','Phf3','Sorcs3','Zhx3','Thsd4','Grm7','Zmat1','Srgap1','Ttc7b','Magi2','Nkain2','Slc24a2','Plcl1','A530053G22Rik','Kctd16','St6galnac3','Tmem178b','Phlpp1','Gm49375','Grik4','St18','Lypd6b','Rnf220','Cadm2','Farp1','Cdk19','Gm816','Epdr1','Tbc1d5','Zeb2','Nfia','Hecw2','Fut9','Ctnna3','Ncam1','Nr3c2','Thrb','Dlg2','Frmd5','St8sia4','Ank2','Dmd','Zdhhc14','Dlg1','Elmo1','Tcf12','5033421B08Rik','Pvt1','Sh3kbp1','Cadps2','Gm29455','Agap1','Efemp1','Mbnl1','Gphn','Npas3','Jph1','Abca1','Naaladl2','Ubash3b','Nbea','Maml3','Skap2','Tanc1','A330023F24Rik','Cntnap2','Pde8a','Megf10','Ankrd28','Sybu','2610035D17Rik','Nlgn1','Tmtc2','Pard3','Frmd4b','Tspan5','Pakap','D7Ertd443e','Tmem117','Gramd1b','Galnt7','Gse1','Atxn1','Cpeb1','Gpr158','Arnt2','Apbb1ip','Sik3','Plekha5','Agmo','Plekhh2','Osbp2','Lpgat1','AW554918','A330072L02Rik','Gria4','Gm48371','Prrg1','Peli2','Exoc4','Ppfibp2','Nsmce2','Slco3a1','Mettl15','Trim2','Epb41l4a','Nos1ap','Gab2','Tmem132b','Gab1','Rftn1','Il1rap','Zcchc7','Cdk5rap2','Xkr4','Dgkb','9330117O12Rik','Lama2','Slc25a13','Macrod2','Dip2c','Prox1','St7','Hipk2','Ctnnd2','Plekhg1','Prkn','Grik2','Tbl1xr1','Aven','Hivep2','Sgms1','Retreg1','Pcca','Aopep','Heg1','4930488L21Rik','Pgbd5','Bach2','Macf1','Zfp462','Cep85l','Samd4','Dock9','Kcnk13','Ccdc88a','Phactr2','Ldlrad3','Gm50323','Tmem196','Tln2','Ttc28','Adam10','E130307A14Rik','Cdkl5','Plekhm3','B3galt5','Fam172a','Myef2','Thada','Tmem135','Pakap.1','Stxbp6','Ncald','Cacna1a','Elp4','Phf21a','Zfp438','Nova1','Pola1','Fign','Epn2','Mtss1','Rps6ka5','Folh1','Ncoa1','Gtdc1','Lrba','Asxl3','1700025G04Rik','Gm19500','Phc3','Amph','Nfe2l3','Inpp5a','9330182L06Rik','Cpeb3','Vwa8','Cpeb4','Diaph2','Lgr4','Megf9','Asap1','Wdpcp','St6gal1','Dpy19l1','Zcchc24','Klhl4','Col4a3bp','BC051408','Arid4b','Pls1','Gm10701','Asah2','Gca','Kif16b','Fbxo11','Cadps','St8sia3os','Zfand3','Efr3b','Atg10','Ercc6l2','Dyrk1a','Nr6a1','Sh3gl2','Ankub1','Prickle2','Unc5c','Il1rapl1','Ttll5','2210408I21Rik','Ssh2','Rybp','Kif13b','Picalm','Gpm6b','Ptprk','Ppp2r3a','Ehbp1','Ugt8a','Baz2b','Ncam2','Tmeff1','Dnm3','Tmcc1','Bicd1','Pacs2','Rnf13','Lrrc8d','Fbxl17','Gm10863','Zswim6','Wdr33','Fam168a','Pex5l','Sash1','Elovl7','Ank','Msi2','Osbpl1a','Zfp148','Eml1','Cox16','Rsf1','Prr5l','Arhgef12','Eri3','Pja2','Fer','Orc4','Tfdp2','Esyt2','Igsf11','Vgll4','Atf7','Cmip','Fars2','Chm','Jakmip2','Xrcc4','Itgb8','8030462N17Rik','Cntn3','Mtmr10','Csrnp3','Ppp3cb','Qk','C4b','Erbin','Mgat4c','Vti1a','Osbpl7','Mob3b','Ppp2r2a','Dnajc6','Nfkb1','Hs6st3','Enah','Myo1e','Lrch3','Pcdh15','Clcn5','Ppm1b','Pcsk6','Bche','Aim2','Wipf1','Sertad2','Fam171a1','Adamtsl4','Pcolce2','Adam19','Dennd5b','Vps13b','Cacnb4','Tnrc6a','Pde4b','Mast4','Fnbp1','Adipor2','Limch1','Mbd5','Dip2b','Otud7b','Rab10','Arid1b','Strn3','Tmcc3','Tbl1x','Grb14','Meis2','Ube2e3','Camkmt','Mecp2','Dock3','Rev3l','Cdkal1','Ankrd11','Mllt10','Zfp407','Dtna','Bcas3','Mtmr3','Rbfox2','Znrf3','Dusp16','Dtnb','Stox2','Gigyf2','Dnajc13','Mir99ahg','Creb5','Fto','Arih1','Gm16233','Gm6145','Cdh18','Slc35f1','Itgb5','Ogdhl','Adgrg2','Itpr1','Htr2c','Amy1','Sema3e','Zc3h6','Dennd1b','Gm2164','Rasal2','9930021J03Rik','Rfx7','Sbf2','Cyp2j12','Pbx3','Gm42495','Gm33948','Chpt1','Gk','Gm32647','Tank','Arhgap15','Tsc22d2','Atp11c','Dlc1','Tbc1d4','Depdc7','Gas2','Gm28322','Mindy3','Ptpn4','Abca6','Taf3','Agps','Ccdc187','Cep128','Crem','Snx13','Rlf','Psmd14','Spef2','Pou2f1','Phka1','Pter','Tgfa','Sgk3','2210408F21Rik','Cd46','BC052040','1110059E24Rik','Tcf20','Antxr1','Umad1','Hdx','Malt1','Zc3h12c','Slf1','Zfp110','Kpna3','Wars2','Fem1c','Btbd10','Sel1l3','A230009B12Rik','Gpd2','Tmem144','Xaf1','Itgav','Shroom3','Nxpe4','Zbtb38','Akap7','Prox1os','Foxo1','Rsu1','Tpk1','Xrcc5','Cachd1','Pag1','Mthfs','B230206H07Rik','Tsga10','Spg20','D3Ertd751e','Hdac4','Dennd2a','Tmem161b','Garnl3')

genes  = c("Mog","Mal","Olig2")

degs = read.csv("../../female_RNA/DEG_results_.01_.01/Oligo_NN.csv")

head(degs)

degs$avg_log2FC =  -degs$avg_log2FC

genes = degs$X[which(degs$p_val_adj<0.05 & degs$avg_log2FC<(-.1))]


genes = degs$X

length(genes)

#genes = genes2

agl = gtf[genes,]
agl = agl[-which(is.na(agl$V1)),]

gr <- GRanges(
  seqnames = Rle(agl[,1]),
  ranges = IRanges(start = agl[,2], end = agl[,3]),gene = agl$V4
)

gr

fd = findOverlaps(gr, dar_gr)

out_dar = dar_gr[fd@to]

out_dar

fo = findOverlaps(gr, cicero_gr)

out_cicero = cicero_gr[fo@to]
out_cicero

diffs$formatted_name = paste(diffs$chr, diffs$loc, sep = "_")
diffs$formatted_name = gsub("-", "_", diffs$formatted_name)

which(out_cicero$tnames %in% diffs$formatted_name)

nrow(diffs[which( diffs$formatted_name %in% c(out_cicero$fnames, out_cicero$tnames)),])

nrow(diffs[which( diffs$feature.name %in% c(out_dar$names)),])

diffs$start = as.numeric(sapply(strsplit(as.character(diffs$loc), "-"), `[`, 1))
diffs$end = as.numeric(sapply(strsplit(as.character(diffs$loc), "-"), `[`, 2))

diff_ov = diffs[c(which(diffs$formatted_name %in% c(out_cicero$fnames, out_cicero$tnames)), which( diffs$feature.name %in% c(out_dar$names))),]

diff_ov=diff_ov[-which(duplicated(diff_ov$feature.name)),]

nrow(diff_ov)

#diffs$log2.fold_change. = -diffs$log2.fold_change.

length(which(diff_ov$log2.fold_change.>0))/nrow(diff_ov)

bed = diff_ov[,c("chr", "start", "end")]
write.table(bed, "../Oligo_overlap_RNA_down_.1_DAR_cicero.bed", sep = "\t", row.names = F, col.names = F, quote=F)

# findMotifsGenome.pl Oligo_overlap_RNA_cluster2_DAR_cicero.bed mm10 Oligo_overlap_RNA_cluster2_DAR_cicero_motifs -size 250  -p 15 
# findMotifsGenome.pl Oligo_overlap_RNA_cluster1_DAR_cicero.bed mm10 Oligo_overlap_RNA_cluster1_DAR_cicero_motifs -size 250  -p 15 
# findMotifsGenome.pl Oligo_overlap_RNA_cluster3_DAR_cicero.bed mm10 Oligo_overlap_RNA_cluster3_DAR_cicero_motifs -size 250  -p 15 
# findMotifsGenome.pl Oligo_overlap_RNA_up_DAR_cicero.bed mm10 Oligo_overlap_RNA_up_DAR_cicero_motifs -size 250  -p 15  -bg Oligo_overlap_RNA_all_DAR_cicero.bed
# findMotifsGenome.pl Oligo_overlap_RNA_down_DAR_cicero.bed mm10 Oligo_overlap_RNA_down_DAR_cicero_motifs -size 250  -p 15 -bg Oligo_overlap_RNA_all_DAR_cicero.bed
findMotifsGenome.pl Oligo_overlap_RNA_cluster2_.1_DAR_cicero.bed mm10  Oligo_overlap_RNA_cluster2_.1_DAR_cicero_motifs -size 250  -p 15 
findMotifsGenome.pl Oligo_overlap_RNA_down_.1_DAR_cicero.bed mm10  Oligo_overlap_RNA_down_.1_DAR_cicero_motifs -size 250  -p 15 
findMotifsGenome.pl Oligo_overlap_RNA_up_.1_DAR_cicero.bed mm10  Oligo_overlap_RNA_up_.1_DAR_cicero_motifs -size 250  -p 15 -bg Oligo_overlap_RNA_all_.1_DAR_cicero.bed
findMotifsGenome.pl Oligo_overlap_RNA_down_.1_DAR_cicero.bed mm10  Oligo_overlap_RNA_down_.1_DAR_cicero_motifs_bg -size 250  -p 15 -bg Oligo_overlap_RNA_all_.1_DAR_cicero.bed

findMotifsGenome.pl Oligo_overlap_RNA_cluster2_.1_DAR_cicero.bed mm10 Oligo_overlap_RNA_cluster2_.1_DAR_cicero_bg -size 250  -p 15 -bg Oligo_overlap_RNA_all_.1_DAR_cicero.bed

which(out_cicero$fnames %in% diffs$formatted_name)

cicero_gr$fnames

out_cicero$tnames


