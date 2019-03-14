% Calculate matrix of centrifugal and coriolis load on the joints for
% S6RPPRRR8
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d4,d5,d6,theta3]';
% m_mdh [7x1]
%   mass of all robot links (including the base)
% rSges [7x3]
%   center of mass of all robot links (in body frames)
%   rows: links of the robot (starting with base)
%   columns: x-, y-, z-coordinates
% Icges [7x6]
%   inertia of all robot links about their respective center of mass, in body frames
%   rows: links of the robot (starting with base)
%   columns: xx, yy, zz, xy, xz, yz (see inertiavector2matrix.m)
% 
% Output:
% Cq [6x6]
%   matrix of coriolis and centrifugal joint torques.
%   Gives coriolis joint torques when multiplied with joint velocities

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 02:36
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Cq = S6RPPRRR8_coriolismatJ_fixb_slag_vp1(qJ, qJD, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(10,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPPRRR8_coriolismatJ_fixb_slag_vp1: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RPPRRR8_coriolismatJ_fixb_slag_vp1: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RPPRRR8_coriolismatJ_fixb_slag_vp1: pkin has to be [10x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RPPRRR8_coriolismatJ_fixb_slag_vp1: m has to be [7x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [7,3]), ...
  'S6RPPRRR8_coriolismatJ_fixb_slag_vp1: rSges has to be [7x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [7 6]), ...
  'S6RPPRRR8_coriolismatJ_fixb_slag_vp1: Icges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From coriolismat_joint_fixb_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 02:35:31
% EndTime: 2019-03-09 02:36:18
% DurationCPUTime: 41.99s
% Computational Cost: add. (144543->1175), mult. (143521->1666), div. (0->0), fcn. (156581->10), ass. (0->731)
t1077 = m(6) * qJD(1);
t1239 = -t1077 / 0.4e1;
t756 = pkin(10) + qJ(4);
t745 = cos(t756);
t763 = sin(qJ(1));
t1000 = t745 * t763;
t744 = sin(t756);
t759 = qJ(5) + qJ(6);
t746 = sin(t759);
t983 = t763 * t746;
t747 = cos(t759);
t765 = cos(qJ(1));
t991 = t747 * t765;
t666 = -t744 * t983 + t991;
t992 = t747 * t763;
t995 = t746 * t765;
t667 = t744 * t992 + t995;
t516 = Icges(7,5) * t667 + Icges(7,6) * t666 - Icges(7,3) * t1000;
t1059 = Icges(7,4) * t667;
t519 = Icges(7,2) * t666 - Icges(7,6) * t1000 + t1059;
t656 = Icges(7,4) * t666;
t522 = Icges(7,1) * t667 - Icges(7,5) * t1000 + t656;
t668 = t744 * t995 + t992;
t669 = t744 * t991 - t983;
t999 = t745 * t765;
t316 = t516 * t999 + t668 * t519 - t669 * t522;
t1247 = t316 * t763;
t1246 = t316 * t765;
t764 = cos(qJ(5));
t978 = t764 * t765;
t762 = sin(qJ(5));
t982 = t763 * t762;
t694 = -t744 * t982 + t978;
t981 = t763 * t764;
t988 = t762 * t765;
t695 = t744 * t981 + t988;
t539 = Icges(6,5) * t695 + Icges(6,6) * t694 - Icges(6,3) * t1000;
t1062 = Icges(6,4) * t695;
t542 = Icges(6,2) * t694 - Icges(6,6) * t1000 + t1062;
t679 = Icges(6,4) * t694;
t545 = Icges(6,1) * t695 - Icges(6,5) * t1000 + t679;
t696 = t744 * t988 + t981;
t697 = t744 * t978 - t982;
t353 = t539 * t999 + t696 * t542 - t697 * t545;
t1245 = t353 * t763;
t1244 = t353 * t765;
t541 = -Icges(6,5) * t697 + Icges(6,6) * t696 + Icges(6,3) * t999;
t1022 = t541 * t744;
t681 = Icges(6,4) * t697;
t544 = Icges(6,2) * t696 + Icges(6,6) * t999 - t681;
t680 = Icges(6,4) * t696;
t546 = Icges(6,1) * t697 - Icges(6,5) * t999 - t680;
t1234 = t544 * t762 + t546 * t764;
t381 = t1234 * t745 - t1022;
t518 = -Icges(7,5) * t669 + Icges(7,6) * t668 + Icges(7,3) * t999;
t1025 = t518 * t744;
t658 = Icges(7,4) * t669;
t521 = Icges(7,2) * t668 + Icges(7,6) * t999 - t658;
t657 = Icges(7,4) * t668;
t523 = Icges(7,1) * t669 - Icges(7,5) * t999 - t657;
t1235 = t521 * t746 + t523 * t747;
t343 = t1235 * t745 - t1025;
t1243 = -2 * m(7);
t1074 = m(7) * qJD(1);
t1228 = -t1074 / 0.4e1;
t819 = -t521 * t668 - t523 * t669;
t947 = t666 * t519 + t667 * t522;
t1242 = t947 + (-t516 * t763 - t518 * t765) * t745 + t819;
t816 = -t544 * t696 - t546 * t697;
t946 = t694 * t542 + t695 * t545;
t1241 = t946 + (-t539 * t763 - t541 * t765) * t745 + t816;
t1183 = t697 * rSges(6,1) - t696 * rSges(6,2);
t552 = -rSges(6,3) * t999 + t1183;
t842 = rSges(6,1) * t764 - rSges(6,2) * t762;
t645 = rSges(6,3) * t744 + t745 * t842;
t1240 = t552 * t744 + t645 * t999;
t1231 = -m(7) / 0.4e1;
t1005 = t744 * t765;
t879 = -t1005 / 0.4e1;
t1198 = 2 * m(7);
t1227 = t1198 / 0.4e1;
t1233 = 2 * m(6);
t1230 = t1233 / 0.4e1;
t528 = -t669 * rSges(7,1) + t668 * rSges(7,2) + rSges(7,3) * t999;
t841 = rSges(7,1) * t747 - rSges(7,2) * t746;
t634 = rSges(7,3) * t744 + t745 * t841;
t598 = t634 * t999;
t1237 = -t528 * t744 + t598;
t565 = rSges(7,1) * t668 + rSges(7,2) * t669;
t683 = t696 * pkin(5);
t1213 = -t565 - t683;
t564 = t666 * rSges(7,1) - t667 * rSges(7,2);
t682 = t694 * pkin(5);
t932 = -t564 - t682;
t1236 = t1213 * t763 + t765 * t932;
t352 = -t1000 * t541 + t694 * t544 - t546 * t695;
t315 = -t1000 * t518 + t666 * t521 - t523 * t667;
t1232 = -m(6) / 0.4e1;
t1144 = m(6) / 0.4e1;
t1143 = m(7) / 0.4e1;
t1132 = t744 / 0.2e1;
t1131 = t745 / 0.2e1;
t1129 = t763 / 0.2e1;
t1128 = t763 / 0.4e1;
t1127 = t765 / 0.2e1;
t1126 = t765 / 0.4e1;
t1199 = 2 * m(5);
t1226 = t1199 / 0.4e1;
t830 = Icges(7,5) * t747 - Icges(7,6) * t746;
t628 = Icges(7,3) * t744 + t745 * t830;
t1057 = Icges(7,4) * t747;
t834 = -Icges(7,2) * t746 + t1057;
t630 = Icges(7,6) * t744 + t745 * t834;
t1058 = Icges(7,4) * t746;
t837 = Icges(7,1) * t747 - t1058;
t632 = Icges(7,5) * t744 + t745 * t837;
t399 = t628 * t999 + t630 * t668 - t632 * t669;
t1223 = t399 * t744;
t317 = t518 * t999 - t819;
t827 = t317 * t765 - t1247;
t163 = t745 * t827 + t1223;
t314 = -t1000 * t516 + t947;
t967 = -t314 + t1242;
t57 = -t1223 + (t765 * t967 + t1247) * t745;
t1225 = t163 + t57;
t192 = t317 * t763 + t1246;
t70 = t763 * t967 - t1246;
t1224 = t192 + t70;
t831 = Icges(6,5) * t764 - Icges(6,6) * t762;
t639 = Icges(6,3) * t744 + t745 * t831;
t1060 = Icges(6,4) * t764;
t835 = -Icges(6,2) * t762 + t1060;
t641 = Icges(6,6) * t744 + t745 * t835;
t1061 = Icges(6,4) * t762;
t838 = Icges(6,1) * t764 - t1061;
t643 = Icges(6,5) * t744 + t745 * t838;
t416 = t639 * t999 + t641 * t696 - t643 * t697;
t1222 = t744 * t416;
t1006 = t744 * t763;
t613 = rSges(6,3) * t1006 + (rSges(6,1) * t981 - rSges(6,2) * t982) * t745;
t690 = pkin(4) * t1000 + pkin(8) * t1006;
t928 = -t613 - t690;
t1221 = t763 * t928;
t1220 = t763 * t932;
t592 = t630 * t763;
t594 = t632 * t763;
t820 = -t519 * t746 + t522 * t747;
t794 = t628 * t763 - t820;
t254 = (-t592 * t746 + t594 * t747 + t516) * t745 + t794 * t744;
t1017 = t628 * t744;
t627 = Icges(7,3) * t745 - t744 * t830;
t993 = t747 * t632;
t996 = t746 * t630;
t809 = t993 - t996;
t790 = t627 - t809;
t1173 = t745 * t790 - t1017;
t629 = Icges(7,6) * t745 - t744 * t834;
t631 = Icges(7,5) * t745 - t744 * t837;
t288 = -t1173 * t763 + t629 * t666 + t631 * t667;
t1218 = t254 + t288;
t593 = t630 * t765;
t595 = t632 * t765;
t793 = -t628 * t765 + t1235;
t255 = (t593 * t746 - t595 * t747 + t518) * t745 + t793 * t744;
t289 = t1173 * t765 + t629 * t668 - t631 * t669;
t1217 = t255 + t289;
t1027 = t516 * t744;
t342 = t745 * t820 + t1027;
t397 = -t1000 * t628 + t630 * t666 + t632 * t667;
t1216 = t342 + t397;
t1215 = -t343 + t399;
t1064 = Icges(5,4) * t744;
t836 = Icges(5,2) * t745 + t1064;
t648 = Icges(5,6) * t765 + t763 * t836;
t724 = Icges(5,4) * t1000;
t650 = Icges(5,1) * t1006 + Icges(5,5) * t765 + t724;
t807 = -t745 * t648 - t744 * t650;
t1214 = t765 * t807;
t924 = t632 + (-Icges(7,2) * t747 - t1058) * t745;
t1014 = t645 * t763;
t711 = t745 * pkin(4) + t744 * pkin(8);
t699 = t763 * t711;
t561 = t699 + t1014;
t1021 = t561 * t763;
t766 = -pkin(9) - pkin(8);
t1078 = -pkin(8) - t766;
t742 = pkin(5) * t764 + pkin(4);
t1080 = pkin(4) - t742;
t618 = t1078 * t744 - t1080 * t745;
t926 = t618 + t634;
t856 = t926 * t763;
t479 = t699 + t856;
t1030 = t479 * t763;
t706 = rSges(5,1) * t745 - rSges(5,2) * t744;
t757 = t763 ^ 2;
t758 = t765 ^ 2;
t913 = t757 + t758;
t1188 = t706 * t913;
t481 = (t711 + t926) * t765;
t1028 = t481 * t765;
t467 = -0.2e1 * t1028;
t563 = (t645 + t711) * t765;
t1019 = t563 * t765;
t535 = -0.2e1 * t1019;
t889 = (t467 - 0.2e1 * t1030) * t1143 + (t535 - 0.2e1 * t1021) * t1144 - t1188 * t1199 / 0.4e1;
t955 = (t467 + 0.2e1 * t1028) * t1143 + (t535 + 0.2e1 * t1019) * t1144;
t1212 = t889 - t955;
t1020 = t563 * t763;
t1029 = t481 * t763;
t466 = -0.2e1 * t1029;
t534 = -0.2e1 * t1020;
t956 = (t466 + 0.2e1 * t1029) * t1143 + (t534 + 0.2e1 * t1020) * t1144;
t1147 = 0.2e1 * t765;
t957 = (t1147 * t479 + t466) * t1143 + (t1147 * t561 + t534) * t1144;
t1211 = t956 - t957;
t526 = t667 * rSges(7,1) + t666 * rSges(7,2) - rSges(7,3) * t1000;
t1181 = t1080 * t744;
t730 = pkin(8) * t1000;
t901 = pkin(5) * t988;
t998 = t745 * t766;
t536 = t901 + t730 + (t998 - t1181) * t763;
t940 = t526 + t536;
t339 = t744 * t940 + t745 * t856;
t1045 = t339 * t763;
t578 = t618 * t999;
t733 = pkin(8) * t999;
t854 = pkin(4) * t1005 - t733;
t902 = pkin(5) * t982;
t916 = t742 * t1005 + t765 * t998;
t537 = t854 + t902 - t916;
t938 = t528 + t537;
t340 = -t744 * t938 + t578 + t598;
t213 = t340 * t1147 + 0.2e1 * t1045;
t551 = rSges(6,1) * t695 + rSges(6,2) * t694 - rSges(6,3) * t1000;
t468 = t1000 * t645 + t551 * t744;
t1033 = t468 * t763;
t358 = t1147 * t1240 + 0.2e1 * t1033;
t965 = t1143 * t213 + t1144 * t358;
t1148 = -0.2e1 * t765;
t341 = -t537 * t744 + t1237 + t578;
t973 = (t1148 * t341 - 0.2e1 * t1045 + t213) * t1143 + (t1148 * t1240 - 0.2e1 * t1033 + t358) * t1144;
t1210 = t965 - t973;
t1044 = t339 * t765;
t1177 = t340 * t763 - t1044;
t212 = 0.2e1 * t1177;
t1032 = t468 * t765;
t1179 = t1240 * t763 - t1032;
t357 = 0.2e1 * t1179;
t966 = t1231 * t212 + t1232 * t357;
t1149 = 0.2e1 * t763;
t972 = (t1149 * t341 - 0.2e1 * t1044 - t212) * t1143 + (t1149 * t1240 - 0.2e1 * t1032 - t357) * t1144;
t1209 = t966 - t972;
t1069 = rSges(6,3) * t745;
t1082 = t744 * pkin(4);
t1194 = pkin(1) * t763;
t1079 = -pkin(7) - qJ(3);
t760 = sin(pkin(10));
t1084 = pkin(3) * t760;
t748 = t765 * qJ(2);
t885 = t763 * t1079 + t765 * t1084 + t748;
t802 = t885 - t1194;
t472 = (-t1069 + t1082) * t765 - t733 + t802 + t1183;
t1083 = pkin(5) * t762;
t884 = pkin(1) + t1083;
t431 = -t763 * t884 - t528 + t885 + t916;
t740 = t765 * t1079;
t860 = qJ(2) + t1084;
t432 = -t740 + t884 * t765 + (t744 * t742 + t860 + t998) * t763 + t526;
t585 = rSges(6,1) * t694 - rSges(6,2) * t695;
t586 = rSges(6,1) * t696 + rSges(6,2) * t697;
t659 = (-rSges(7,1) * t746 - rSges(7,2) * t747) * t745;
t853 = t1083 * t745 - t659;
t587 = t853 * t763;
t588 = t853 * t765;
t685 = (-rSges(6,1) * t762 - rSges(6,2) * t764) * t745;
t1081 = t765 * pkin(1);
t473 = t1081 - t730 - t740 + (t860 + t1082) * t763 + t551;
t821 = t472 * t763 - t473 * t765;
t1206 = (t1213 * t479 - t431 * t587 + t432 * t588 + t481 * t932) * t1227 + (-t561 * t586 - t563 * t585 + t685 * t821) * t1230;
t1004 = t744 * t766;
t707 = t742 * t1000;
t571 = -t1004 * t763 - t690 + t707;
t617 = t1078 * t745 + t1181;
t1002 = t745 * t746;
t899 = rSges(7,2) * t1002;
t915 = t745 * rSges(7,1) * t992 + rSges(7,3) * t1006;
t599 = -t763 * t899 + t915;
t633 = rSges(7,3) * t745 - t744 * t841;
t886 = t633 * t1000 + t745 * t526 + t744 * t599;
t238 = (t617 * t763 + t536) * t745 + (t571 - t856) * t744 + t886;
t597 = t633 * t999;
t1189 = t634 * t765;
t914 = -pkin(4) * t999 - pkin(8) * t1005;
t929 = (-t745 * t742 + t1004) * t765 - t914 - t1189;
t239 = t597 + (t617 * t765 - t938) * t745 + (-t765 * t926 - t929) * t744;
t644 = -t744 * t842 + t1069;
t1015 = t644 * t763;
t384 = (t551 + t1015) * t745 + (t613 - t1014) * t744;
t385 = (t644 * t765 + t552) * t745;
t511 = t707 + (-t899 - t1004) * t763 + t915;
t512 = ((rSges(7,3) - t766) * t744 + (t742 + t841) * t745) * t765;
t1013 = t645 * t765;
t533 = -t914 + t1013;
t1205 = -(t1240 * t533 + t384 * t473 + t385 * t472 - t468 * t928) * t1233 / 0.4e1 - (t238 * t432 + t239 * t431 + t339 * t511 + t340 * t512) * t1198 / 0.4e1;
t960 = t340 - t341;
t1204 = 0.2e1 * t481 * t960 * t1231;
t649 = -Icges(5,6) * t763 + t765 * t836;
t1063 = Icges(5,4) * t745;
t839 = Icges(5,1) * t744 + t1063;
t651 = -Icges(5,5) * t763 + t765 * t839;
t674 = -Icges(5,2) * t1006 + t724;
t702 = -Icges(5,2) * t744 + t1063;
t675 = t702 * t765;
t704 = Icges(5,1) * t745 - t1064;
t677 = t704 * t763;
t678 = t704 * t765;
t1203 = -(t763 * (t649 - t678) - t765 * (t648 - t677)) * t744 + (t763 * (t651 + t675) - t765 * (t650 + t674)) * t745;
t653 = (-Icges(7,5) * t746 - Icges(7,6) * t747) * t745;
t1008 = t744 * t653;
t1202 = t1008 / 0.2e1 - t924 * t1002 / 0.2e1;
t1196 = -m(5) / 0.2e1;
t1195 = t757 / 0.2e1;
t1130 = -t763 / 0.2e1;
t1193 = t315 * t763;
t1192 = t315 * t765;
t1191 = t352 * t763;
t1190 = t352 * t765;
t1187 = (t745 * t649 + t744 * t651) * t765;
t1186 = t238 + t588;
t1185 = t239 - t587;
t538 = t765 * t565;
t446 = -t564 * t763 + t538;
t452 = t1000 * t634 + t526 * t744;
t1035 = t452 * t765;
t1178 = t1237 * t763 - t1035;
t554 = Icges(7,5) * t666 - Icges(7,6) * t667;
t942 = -Icges(7,2) * t667 + t522 + t656;
t944 = -Icges(7,1) * t666 + t1059 + t519;
t244 = -t1000 * t554 + t666 * t942 - t667 * t944;
t555 = Icges(7,5) * t668 + Icges(7,6) * t669;
t941 = Icges(7,2) * t669 - t523 + t657;
t943 = -Icges(7,1) * t668 + t521 - t658;
t245 = -t1000 * t555 + t666 * t941 - t667 * t943;
t140 = t244 * t765 + t245 * t763;
t246 = t554 * t999 + t668 * t942 + t669 * t944;
t247 = t555 * t999 + t668 * t941 + t669 * t943;
t141 = t246 * t765 + t247 * t763;
t971 = t1127 * t140 + t1129 * t141;
t994 = t747 * t631;
t997 = t746 * t629;
t303 = (t628 + t994 - t997) * t745 + t790 * t744;
t422 = t745 * t809 + t1017;
t826 = -t342 * t763 - t343 * t765;
t1176 = (t422 * t744 + t745 * t826) * t1131 + ((-t254 * t763 + t255 * t765 + t422) * t745 + (t303 - t826) * t744) * t1132;
t1175 = t422 * t1131 + t1132 * t303;
t810 = t585 * t763 - t586 * t765;
t952 = (t1213 * t765 - t1220) * t1227 + t810 * t1230;
t811 = -t585 * t765 - t586 * t763;
t953 = t1227 * t1236 + t811 * t1230;
t295 = (-t763 * t938 - t765 * t940) * t745;
t813 = -t564 * t765 - t565 * t763;
t429 = t813 * t745;
t485 = t659 * t1000 + t744 * t564;
t624 = t659 * t999;
t486 = -t565 * t744 + t624;
t847 = t295 * t429 + t339 * t485 + t340 * t486;
t1016 = t639 * t744;
t638 = Icges(6,3) * t745 - t744 * t831;
t979 = t764 * t643;
t989 = t762 * t641;
t808 = t979 - t989;
t789 = t638 - t808;
t1174 = t745 * t789 - t1016;
t791 = -t639 * t765 + t1234;
t1172 = t745 * t791 - t1022;
t1024 = t539 * t744;
t817 = -t542 * t762 + t545 * t764;
t792 = t639 * t763 - t817;
t1171 = t745 * t792 - t1024;
t1170 = t745 * t793 - t1025;
t1169 = t745 * t794 - t1027;
t844 = rSges(5,1) * t744 + rSges(5,2) * t745;
t783 = -t763 * rSges(5,3) + t765 * t844;
t568 = t783 + t802;
t569 = -t740 + (rSges(5,3) + pkin(1)) * t765 + (t844 + t860) * t763;
t684 = t706 * t763;
t686 = t706 * t765;
t1168 = -m(5) * (t568 * t686 + t569 * t684) - (-t993 / 0.2e1 + t996 / 0.2e1 + t627 / 0.2e1 - t979 / 0.2e1 + t989 / 0.2e1 + t638 / 0.2e1 - t704 / 0.2e1 + t836 / 0.2e1) * t744;
t804 = t684 * t763 + t686 * t765;
t887 = (t511 * t763 + t512 * t765) * t1227 + (t533 * t765 - t1221) * t1230 + t804 * t1226;
t805 = -t684 * t765 + t686 * t763;
t888 = (-t511 * t765 + t512 * t763) * t1227 + (t533 * t763 + t765 * t928) * t1230 + t805 * t1226;
t673 = (-Icges(6,2) * t764 - t1061) * t745;
t676 = (-Icges(6,1) * t762 - t1060) * t745;
t1167 = -t762 * (t643 / 0.2e1 + t673 / 0.2e1) + t764 * (t676 / 0.2e1 - t641 / 0.2e1);
t579 = Icges(6,5) * t694 - Icges(6,6) * t695;
t935 = -Icges(6,2) * t695 + t545 + t679;
t937 = -Icges(6,1) * t694 + t1062 + t542;
t293 = t579 * t744 + (-t762 * t935 - t764 * t937) * t745;
t580 = Icges(6,5) * t696 + Icges(6,6) * t697;
t934 = Icges(6,2) * t697 - t546 + t680;
t936 = -Icges(6,1) * t696 + t544 - t681;
t294 = t580 * t744 + (-t762 * t934 - t764 * t936) * t745;
t670 = (-Icges(6,5) * t762 - Icges(6,6) * t764) * t745;
t922 = t643 + t673;
t923 = t641 - t676;
t347 = -t1000 * t670 + t694 * t922 - t695 * t923;
t348 = t670 * t999 + t696 * t922 + t697 * t923;
t1164 = t1206 + (t294 + t348) * t1128 + (t293 + t347) * t1126;
t272 = t555 * t744 + (-t746 * t941 - t747 * t943) * t745;
t1049 = t272 * t763;
t271 = t554 * t744 + (-t746 * t942 - t747 * t944) * t745;
t1050 = t271 * t765;
t655 = (-Icges(7,1) * t746 - t1057) * t745;
t925 = t630 - t655;
t310 = -t1000 * t653 + t666 * t924 - t667 * t925;
t311 = t653 * t999 + t668 * t924 + t669 * t925;
t855 = t1050 / 0.4e1 + t1049 / 0.4e1 + t311 * t1128 + t310 * t1126;
t393 = t397 * t744;
t970 = t317 + t1242;
t56 = t393 + (-t763 * t970 + t1192) * t745;
t976 = t765 * t163;
t828 = -t314 * t763 + t1192;
t162 = t745 * t828 + t393;
t986 = t763 * t162;
t1163 = t57 * t1126 + t56 * t1128 + t976 / 0.4e1 - t986 / 0.4e1;
t351 = -t1000 * t539 + t946;
t219 = t351 * t765 + t1191;
t354 = t541 * t999 - t816;
t220 = t354 * t763 + t1244;
t414 = -t1000 * t639 + t641 * t694 + t643 * t695;
t411 = t414 * t744;
t964 = t354 + t1241;
t73 = t411 + (-t763 * t964 + t1190) * t745;
t961 = -t351 + t1241;
t74 = -t1222 + (t765 * t961 + t1245) * t745;
t86 = t765 * t964 + t1191;
t866 = t999 / 0.4e1;
t868 = -t999 / 0.4e1;
t87 = t763 * t961 - t1244;
t872 = -t1000 / 0.4e1;
t824 = t354 * t765 - t1245;
t181 = t745 * t824 + t1222;
t975 = t765 * t181;
t825 = -t351 * t763 + t1190;
t180 = t745 * t825 + t411;
t985 = t763 * t180;
t1162 = -t985 / 0.4e1 + t975 / 0.4e1 + t219 * t868 + t73 * t1128 + t74 * t1126 + t86 * t866 - t1204 + (t220 + t87) * t872;
t608 = t641 * t763;
t610 = t643 * t763;
t277 = (-t608 * t762 + t610 * t764 + t539) * t745 + t792 * t744;
t609 = t641 * t765;
t611 = t643 * t765;
t278 = (t609 * t762 - t611 * t764 + t541) * t745 + t791 * t744;
t640 = Icges(6,6) * t745 - t744 * t835;
t642 = Icges(6,5) * t745 - t744 * t838;
t304 = -t1174 * t763 + t640 * t694 + t642 * t695;
t305 = t1174 * t765 + t640 * t696 - t642 * t697;
t980 = t764 * t642;
t990 = t762 * t640;
t331 = (t639 + t980 - t990) * t745 + t789 * t744;
t380 = t745 * t817 + t1024;
t444 = t745 * t808 + t1016;
t881 = t1006 / 0.4e1;
t1161 = t444 * t1131 + t331 * t1132 - t1205 + (t380 + t414) * t881 + (-t381 + t416) * t879 + (t277 + t304) * t872 + (t278 + t305) * t866;
t1154 = 0.2e1 * t452;
t1153 = 0.2e1 * t1237;
t1152 = -0.2e1 * t481;
t1151 = 0.2e1 * t486;
t1150 = -0.2e1 * t565;
t1146 = -m(4) / 0.4e1;
t1145 = -m(5) / 0.4e1;
t930 = -t571 - t599;
t945 = t526 * t1005 + t528 * t1006;
t197 = (t536 * t765 + t537 * t763) * t744 + (-t763 * t929 + t765 * t930) * t745 + t945;
t318 = (t1189 * t763 - t599 * t765) * t745 + t945;
t366 = -t1006 * t634 + t886;
t367 = -t528 * t745 + t597;
t412 = (-t526 * t765 - t528 * t763) * t745;
t1141 = (t1237 * t239 + t197 * t412 + t238 * t452 + t295 * t318 + t339 * t366 + t340 * t367) * t1198;
t1139 = t452 * t960 * t1243;
t664 = t765 * t854;
t689 = pkin(4) * t1006 - t730;
t287 = -t664 + t938 * t765 + (-t689 - t940) * t763;
t225 = 0.2e1 * t429 * t287;
t378 = t485 * t1152;
t379 = t479 * t1151;
t890 = t225 + t378 + t379;
t904 = 0.2e1 * t659;
t1137 = m(7) * (t1177 * t904 + 0.2e1 * t295 * t446 + t890);
t396 = t683 * t765 + t1220 + t538;
t1135 = m(7) * (-t1153 * t587 + t1154 * t588 + 0.2e1 * t396 * t412 + t890);
t1134 = -t180 / 0.2e1;
t1133 = t294 / 0.2e1;
t736 = 0.2e1 * t913;
t803 = t736 / 0.4e1;
t1107 = (t1237 * t512 + t366 * t432 + t367 * t431 + t452 * t511) * t1198;
t345 = 0.2e1 * t485 * t432;
t346 = t431 * t1151;
t958 = t345 + t346;
t1106 = m(7) * (t1150 * t340 + 0.2e1 * t339 * t564 + t958);
t1103 = m(7) * (t1153 * t1213 - t1154 * t932 + t958);
t822 = t763 * t431 - t432 * t765;
t1101 = m(7) * (t1150 * t479 + t1152 * t564 + t822 * t904);
t1036 = t452 * t763;
t324 = t1147 * t1237 + 0.2e1 * t1036;
t1100 = m(7) * (t1148 * t1237 - 0.2e1 * t1036 + t324);
t323 = 0.2e1 * t1178;
t1099 = m(7) * (t1149 * t1237 - 0.2e1 * t1035 - t323);
t1096 = m(7) * t323;
t1095 = m(7) * t324;
t1086 = t813 * t1198;
t1085 = t446 * t1243;
t1076 = m(6) * qJD(4);
t1075 = m(6) * qJD(5);
t1073 = m(7) * qJD(4);
t1072 = m(7) * qJD(5);
t1071 = m(7) * qJD(6);
t1041 = t366 * t765;
t1040 = t367 * t763;
t1039 = t384 * t765;
t1038 = t385 * t763;
t1007 = t744 * t670;
t1001 = t745 * t747;
t927 = t617 + t633;
t912 = qJD(1) * t745;
t911 = qJD(5) * t745;
t861 = t1195 + t758 / 0.2e1;
t199 = (t1041 / 0.2e1 - t1040 / 0.2e1 + t861 * t659) * m(7);
t910 = t199 * qJD(2);
t235 = 0.2e1 * t366 * t763 + 0.2e1 * t367 * t765;
t909 = t235 * qJD(3);
t494 = (-m(7) / 0.2e1 - m(6) / 0.2e1 + t1196 - m(4) / 0.2e1) * t736;
t908 = t494 * qJD(1);
t903 = rSges(4,3) + pkin(1) + qJ(3);
t365 = (t1008 + (-t746 * t924 - t747 * t925) * t745) * t744;
t867 = t999 / 0.2e1;
t873 = -t1000 / 0.2e1;
t94 = t310 * t744 + (-t244 * t763 + t245 * t765) * t745;
t95 = t311 * t744 + (-t246 * t763 + t247 * t765) * t745;
t900 = t95 * t867 + t94 * t873 + (t365 + (-t271 * t763 + t272 * t765) * t745) * t1132;
t896 = t235 * t1143;
t894 = t1072 / 0.2e1;
t893 = t1071 / 0.2e1;
t892 = t1134 + t73 / 0.2e1;
t891 = -t74 / 0.2e1 - t181 / 0.2e1;
t832 = Icges(5,5) * t744 + Icges(5,6) * t745;
t458 = t648 * t1000 + t650 * t1006 + t765 * (Icges(5,3) * t765 + t763 * t832);
t647 = -Icges(5,3) * t763 + t765 * t832;
t459 = -t649 * t1000 - t651 * t1006 - t765 * t647;
t882 = t1006 / 0.2e1;
t880 = -t1005 / 0.2e1;
t875 = -t1001 / 0.2e1;
t874 = t1001 / 0.2e1;
t871 = t1000 / 0.2e1;
t870 = t1000 / 0.4e1;
t869 = -t999 / 0.2e1;
t865 = -t587 / 0.2e1 - t239 / 0.2e1;
t864 = t588 / 0.2e1 - t238 / 0.2e1;
t858 = 0.2e1 * t1143;
t743 = t745 ^ 2;
t419 = t682 * t744 - t743 * t902 + t485;
t497 = t1000 * t685 + t585 * t744;
t851 = m(6) * t497 + m(7) * t419;
t420 = t1213 * t744 - t743 * t901 + t624;
t498 = -t586 * t744 + t685 * t999;
t850 = m(6) * t498 + m(7) * t420;
t710 = t745 * pkin(8) - t1082;
t698 = t763 * t710;
t478 = t763 * t927 + t698;
t560 = t698 + t1015;
t849 = m(6) * t560 + m(7) * t478;
t480 = (-t710 - t927) * t765;
t562 = (-t644 - t710) * t765;
t848 = m(6) * t562 + m(7) * t480;
t846 = t1225 * t873 + t162 * t869 + t56 * t867;
t845 = rSges(4,1) * t760 + rSges(4,2) * cos(pkin(10));
t840 = t630 * t875 + t655 * t874 + t1202;
t833 = Icges(5,5) * t745 - Icges(5,6) * t744;
t81 = (t763 * t865 - t765 * t864) * m(7) + (t1039 / 0.2e1 - t1038 / 0.2e1 + t861 * t685) * m(6);
t784 = (t384 * t763 + t385 * t765) * t1233;
t85 = -t784 / 0.4e1 + (t763 * t864 + t765 * t865) * m(7);
t829 = -t81 * qJD(2) - t85 * qJD(3);
t823 = -t380 * t763 - t381 * t765;
t814 = t551 * t765 - t552 * t763;
t812 = t568 * t763 - t569 * t765;
t231 = -t1169 * t763 + t592 * t666 + t594 * t667;
t232 = -t1170 * t763 - t593 * t666 - t595 * t667;
t42 = (-t231 * t763 + t232 * t765 + t397) * t745 + (t288 - t828) * t744;
t233 = t1169 * t765 + t592 * t668 - t594 * t669;
t234 = t1170 * t765 - t593 * t668 + t595 * t669;
t43 = (-t233 * t763 + t234 * t765 + t399) * t745 + (t289 - t827) * t744;
t801 = t162 * t882 + t163 * t880 + t42 * t873 + t43 * t867 + t1176;
t800 = t1139 / 0.4e1 + t846;
t615 = -t763 * t903 + t765 * t845 + t748;
t616 = t903 * t765 + (qJ(2) + t845) * t763;
t799 = -0.4e1 * t812 * t1145 + 0.4e1 * (-t615 * t763 + t616 * t765) * t1146;
t273 = -t1000 * t579 + t694 * t935 - t695 * t937;
t274 = -t1000 * t580 + t694 * t934 - t695 * t936;
t166 = t273 * t765 + t274 * t763;
t275 = t579 * t999 + t696 * t935 + t697 * t937;
t276 = t580 * t999 + t696 * t934 + t697 * t936;
t167 = t275 * t765 + t276 * t763;
t788 = t1127 * t166 + t1129 * t167;
t786 = t1141 / 0.4e1 + t801;
t785 = t1237 * t486 + t412 * t429 + t452 * t485;
t191 = t314 * t765 + t1193;
t69 = t765 * t970 + t1193;
t782 = t1224 * t872 + t191 * t868 + t69 * t866 + t1163;
t781 = 0.4e1 * (t568 * t765 + t569 * t763) * t1145 + 0.4e1 * (t615 * t765 + t616 * t763) * t1146 - m(3) * ((rSges(3,3) * t765 - t1194 + t748) * t765 + (t1081 + (rSges(3,3) + qJ(2)) * t763) * t763);
t780 = t630 * t874 + t655 * t875 - t1202;
t461 = -t763 * t647 + t1187;
t779 = t192 / 0.2e1 + t70 / 0.2e1 + t220 / 0.2e1 + t87 / 0.2e1 + t647 * t1195 + t461 * t1129 + (t459 + (t647 - t807) * t765 + t1214) * t1127;
t778 = t69 / 0.2e1 - t191 / 0.2e1 + t86 / 0.2e1 - t219 / 0.2e1 - t458 * t765 / 0.2e1 + t459 * t1130 + (t459 - t1214) * t1129 + (t461 - t1187 + (t647 + t807) * t763 + t458) * t1127;
t122 = t231 * t765 + t232 * t763;
t123 = t233 * t765 + t234 * t763;
t777 = t42 * t1127 + t43 * t1129 + t122 * t873 + t123 * t867 + (t254 * t765 + t255 * t763) * t1132 + t191 * t882 + t192 * t880 + (t342 * t765 - t343 * t763) * t1131 - t971;
t776 = t56 * t869 + t365 + (t271 + t310) * t873 + t1225 * t871 + (t162 + t272 + t311) * t867;
t775 = t1215 * t879 + t1216 * t881 + t1217 * t866 + t1218 * t872 + t1175;
t774 = t42 * t871 + t43 * t869 + t140 * t873 + t141 * t867 - t744 * t986 / 0.2e1 + t95 * t1129 + t94 * t1127 - t1176 + (t976 + t1049 + t1050) * t1132;
t771 = t994 / 0.2e1 - t997 / 0.2e1 + t628 / 0.2e1 + t980 / 0.2e1 - t990 / 0.2e1 + t639 / 0.2e1 - t839 / 0.2e1 - t702 / 0.2e1;
t769 = t1224 * t870 + t191 * t866 + t69 * t868 - t1163 + t775 + t855;
t768 = -t1175 + t782 + t855 - t1216 * t1006 / 0.4e1 + t1215 * t1005 / 0.4e1 + t1218 * t870 + t1217 * t868;
t767 = t775 + t782 - t855;
t672 = t833 * t765;
t671 = t763 * t833;
t665 = t765 * t914;
t493 = (t1231 + t1232 + t1145 + t1146) * t736 + (m(4) + m(5) + m(6) + m(7)) * t803;
t457 = t811 * t745;
t439 = t1085 / 0.4e1;
t438 = t1086 / 0.4e1;
t425 = t814 * t745;
t421 = -t1013 * t765 + t1221 + t665;
t395 = -t552 * t765 - t664 + (-t551 - t689) * t763;
t386 = (t1007 + (-t762 * t922 - t764 * t923) * t745) * t744;
t383 = t1236 * t745;
t377 = 0.2e1 * t485 * t763 + 0.2e1 * t486 * t765;
t376 = -0.2e1 * t485 * t765 + 0.2e1 * t486 * t763;
t373 = t377 * t893;
t372 = t376 * t893;
t360 = 0.4e1 * t765 * t472 + 0.4e1 * t473 * t763;
t359 = 0.4e1 * t821;
t336 = (t1013 * t763 - t613 * t765) * t745 + t814 * t744;
t321 = t1095 / 0.4e1;
t320 = -t1096 / 0.4e1;
t319 = t665 + t929 * t765 + (-t690 + t930) * t763;
t313 = 0.4e1 * t765 * t431 + 0.4e1 * t432 * t763;
t312 = 0.4e1 * t822;
t292 = -0.4e1 * t472 * t586 + 0.4e1 * t473 * t585;
t281 = 0.4e1 * t472 * t533 - 0.4e1 * t473 * t928;
t262 = -0.4e1 * t431 * t565 + 0.4e1 * t432 * t564;
t261 = t1172 * t765 - t609 * t696 + t611 * t697;
t260 = t1171 * t765 + t608 * t696 - t610 * t697;
t259 = -t1172 * t763 - t609 * t694 - t611 * t695;
t258 = -t1171 * t763 + t608 * t694 + t610 * t695;
t253 = 0.4e1 * t431 * t512 + 0.4e1 * t432 * t511;
t243 = 0.4e1 * t1213 * t431 - 0.4e1 * t432 * t932;
t227 = -0.4e1 * t810 * t395 + 0.4e1 * (t1019 + t1021) * t685;
t226 = qJD(4) * t896;
t201 = t1099 / 0.4e1;
t200 = t1100 / 0.4e1;
t198 = (t1040 - t1041) * t858 + m(7) * t659 * t803;
t184 = t444 * t744 + t745 * t823;
t182 = t1143 * t262 + t840;
t173 = t1101 / 0.4e1;
t169 = -t1143 * t312 - t1144 * t359 - t799;
t168 = 0.4e1 * t287 * t446 + 0.4e1 * (t1028 + t1030) * t659;
t161 = 0.4e1 * t785;
t147 = t1143 * t313 + t1144 * t360 - t781;
t146 = t260 * t765 + t261 * t763;
t145 = t258 * t765 + t259 * t763;
t143 = 0.4e1 * t287 * t396 - 0.4e1 * t479 * t587 - 0.4e1 * t481 * t588;
t132 = t1103 / 0.4e1;
t131 = t321 + t200 - t1086 / 0.4e1;
t130 = t320 + t201 - t1085 / 0.4e1;
t129 = t439 + t320 - t1099 / 0.4e1;
t128 = t439 + t201 + t1096 / 0.4e1;
t127 = t438 + t321 - t1100 / 0.4e1;
t126 = t438 + t200 - t1095 / 0.4e1;
t112 = t348 * t744 + (-t275 * t763 + t276 * t765) * t745;
t111 = t347 * t744 + (-t273 * t763 + t274 * t765) * t745;
t110 = 0.4e1 * t1240 * t385 - 0.4e1 * t336 * t425 + 0.4e1 * t384 * t468;
t109 = t956 + t957 - t888;
t108 = t888 + t1211;
t107 = t888 - t1211;
t105 = t1106 / 0.4e1;
t103 = t889 + t955 - t887;
t102 = t887 - t1212;
t101 = t887 + t1212;
t99 = t1107 / 0.4e1;
t98 = 0.4e1 * t1237 * t367 + 0.4e1 * t318 * t412 + 0.4e1 * t366 * t452;
t97 = 0.4e1 * t847;
t96 = t243 * t1143 + t1007 / 0.2e1 + t292 * t1144 + t1167 * t745 + t840;
t84 = t784 / 0.4e1 + (t1185 * t765 + t1186 * t763) * t858;
t80 = (t1038 - t1039) * t1230 + m(6) * t685 * t803 + (t1185 * t763 - t1186 * t765) * t858;
t79 = (-t277 * t763 + t278 * t765 + t444) * t745 + (t331 - t823) * t744;
t77 = 0.4e1 * t339 * t960;
t76 = t1135 / 0.4e1;
t75 = t253 * t1143 + t281 * t1144 + t771 * t745 - t1168;
t64 = (-t260 * t763 + t261 * t765 + t416) * t745 + (t305 - t824) * t744;
t63 = (-t258 * t763 + t259 * t765 + t414) * t745 + (t304 - t825) * t744;
t50 = t1137 / 0.4e1;
t35 = 0.4e1 * t197 * t295 + 0.4e1 * t238 * t339 + 0.4e1 * t239 * t340;
t32 = t965 + t973 - t953;
t31 = t966 + t972 - t952;
t30 = t952 - t1209;
t29 = t952 + t1209;
t28 = t953 + t1210;
t27 = t953 - t1210;
t26 = t1143 * t168 + t971;
t21 = t1143 * t161 + t900;
t20 = t21 * qJD(6);
t19 = t1143 * t97 + t900;
t18 = t1143 * t143 + t1144 * t227 + t788 + t971;
t16 = t1143 * t98 + t801;
t15 = t132 - t1106 / 0.4e1 + t800;
t14 = t105 - t1103 / 0.4e1 + t800;
t13 = t105 + t132 - t1139 / 0.4e1 + t776;
t12 = t763 * t778 + t765 * t779;
t11 = -t77 * t1143 + (t763 * t891 + t765 * t892) * t745 + t846;
t10 = t76 - t1137 / 0.4e1 + t786;
t9 = t50 - t1135 / 0.4e1 + t786;
t8 = t173 + t99 + t769;
t7 = t173 + t768 - t1107 / 0.4e1;
t6 = t767 + t99 - t1101 / 0.4e1;
t5 = t35 * t1143 + t110 * t1144 + (t184 / 0.2e1 + t63 * t1130 + t64 * t1127) * t745 + (t79 / 0.2e1 + t985 / 0.2e1 - t975 / 0.2e1) * t744 + t801;
t4 = t76 + t774 + t50 - t1141 / 0.4e1;
t3 = t1161 + t769 + (-t73 / 0.4e1 + t180 / 0.4e1) * t763 + (-t181 / 0.4e1 - t74 / 0.4e1) * t765 + ((t219 / 0.4e1 - t86 / 0.4e1) * t765 + (t220 / 0.4e1 + t87 / 0.4e1) * t763) * t745 + t1164 + t1204;
t2 = t1161 + t1162 + t767 + (-t348 / 0.4e1 - t294 / 0.4e1) * t763 + (-t347 / 0.4e1 - t293 / 0.4e1) * t765 - t1206;
t1 = t1162 + t768 + (-t444 / 0.2e1 + (-t305 / 0.4e1 - t278 / 0.4e1) * t765 + (t304 / 0.4e1 + t277 / 0.4e1) * t763) * t745 + (-t331 / 0.2e1 + (t416 / 0.4e1 - t381 / 0.4e1) * t765 + (-t414 / 0.4e1 - t380 / 0.4e1) * t763) * t744 + t1164 + t1205;
t17 = [t147 * qJD(2) + t169 * qJD(3) + t75 * qJD(4) + t96 * qJD(5) + t182 * qJD(6), qJD(1) * t147 + qJD(3) * t493 + qJD(4) * t107 + qJD(5) * t28 + qJD(6) * t127, qJD(1) * t169 + qJD(2) * t493 + qJD(4) * t101 + qJD(5) * t29 + qJD(6) * t129, t75 * qJD(1) + t107 * qJD(2) + t101 * qJD(3) + t3 * qJD(5) + t8 * qJD(6) + (t431 * t478 + t432 * t480 + t479 * t512 - t481 * t511) * t1073 + (t472 * t560 + t473 * t562 + t533 * t561 + t928 * t563) * t1076 + ((t277 / 0.2e1 + t304 / 0.2e1 + t288 / 0.2e1 + t254 / 0.2e1 - t832 * t1127 + (-t648 / 0.2e1 + t677 / 0.2e1) * t745 + (-t650 / 0.2e1 - t674 / 0.2e1) * t744 - t779) * t765 + (t305 / 0.2e1 + t289 / 0.2e1 + t278 / 0.2e1 + t255 / 0.2e1 - t832 * t1129 + (t649 / 0.2e1 - t678 / 0.2e1) * t745 + (t651 / 0.2e1 + t675 / 0.2e1) * t744 - t778) * t763 + (t706 * t805 - t812 * t844) * m(5)) * qJD(4), t96 * qJD(1) + t28 * qJD(2) + t29 * qJD(3) + t3 * qJD(4) + (t386 + t776) * qJD(5) + t13 * qJD(6) + (t77 / 0.4e1 - t339 * t932 + t340 * t1213 + t419 * t432 + t420 * t431) * t1072 + (-t1240 * t586 + t468 * t585 + t472 * t498 + t473 * t497) * t1075 + ((t1133 + t348 / 0.2e1 - t892) * t765 + (-t293 / 0.2e1 - t347 / 0.2e1 - t891) * t763) * t911, t182 * qJD(1) + t127 * qJD(2) + t129 * qJD(3) + t8 * qJD(4) + t13 * qJD(5) + t776 * qJD(6) + (t452 * t564 - t1237 * t565 + t345 / 0.2e1 + t346 / 0.2e1) * t1071; t781 * qJD(1) + t494 * qJD(3) + t108 * qJD(4) + t27 * qJD(5) + t126 * qJD(6) + t313 * t1228 + t1239 * t360, 0, t908, t108 * qJD(1) + (t1196 * t736 * t844 + t849 * t763 - t848 * t765) * qJD(4) + t80 * qJD(5) + t198 * qJD(6), t27 * qJD(1) + t80 * qJD(4) + (t763 * t850 - t765 * t851) * qJD(5) + t372, t126 * qJD(1) + t198 * qJD(4) + t376 * t894 + t372; t799 * qJD(1) - t494 * qJD(2) + t102 * qJD(4) + t30 * qJD(5) + t128 * qJD(6) + t312 * t1074 / 0.4e1 + t359 * t1077 / 0.4e1, -t908, 0, t102 * qJD(1) + (t763 * t848 + t765 * t849) * qJD(4) + t84 * qJD(5) + qJD(6) * t896, t30 * qJD(1) + t84 * qJD(4) + (t763 * t851 + t765 * t850) * qJD(5) + t373, t128 * qJD(1) + t377 * t894 + t226 + t373; t1168 * qJD(1) + t109 * qJD(2) + t103 * qJD(3) + t12 * qJD(4) + t1 * qJD(5) + t7 * qJD(6) + t253 * t1228 + t281 * t1239 - t771 * t912, qJD(1) * t109 + qJD(5) * t81 + qJD(6) * t199, t103 * qJD(1) + t85 * qJD(5) - t235 * t1071 / 0.4e1, t12 * qJD(1) + t18 * qJD(5) + t26 * qJD(6) + ((t287 * t319 + t478 * t479 - t480 * t481) * m(7) + (t395 * t421 + t560 * t561 - t562 * t563) * m(6) + m(5) * (-t844 * t1188 - (-t765 * t783 + (-t765 * rSges(5,3) - t763 * t844) * t763) * t804) + (t123 + t146 - t757 * t672 + (t763 * t671 + t1203) * t765) * t1129 + (t122 + t145 + t758 * t671 + (-t765 * t672 - t1203) * t763) * t1127) * qJD(4), t1 * qJD(1) + t18 * qJD(4) + t4 * qJD(6) + (-t184 / 0.2e1 + (-t64 / 0.2e1 + t167 / 0.2e1) * t765 + (t63 / 0.2e1 - t166 / 0.2e1) * t763) * t911 + (t287 * t383 + t295 * t396 + t339 * t588 - t340 * t587 - t419 * t481 + t420 * t479 - t35 / 0.4e1) * t1072 + (t457 * t395 + t425 * t810 - t497 * t563 + t498 * t561 - t110 / 0.4e1 + t1179 * t685) * t1075 - t829 + (t111 * t1127 + t112 * t1129 + t774 + (-t79 / 0.2e1 + (t293 / 0.2e1 + t181 / 0.2e1) * t765 + (t1133 + t1134) * t763) * t744) * qJD(5), t7 * qJD(1) + t910 + t26 * qJD(4) + t4 * qJD(5) + t774 * qJD(6) + (-t909 / 0.4e1 + (t412 * t446 + t225 / 0.2e1 + t378 / 0.2e1 + t379 / 0.2e1 - t98 / 0.4e1 + t1178 * t659) * qJD(6)) * m(7); (-t1007 / 0.2e1 + t780) * qJD(1) + t32 * qJD(2) + t31 * qJD(3) + t2 * qJD(4) + t11 * qJD(5) + t14 * qJD(6) + (-t243 / 0.4e1 - t960 * t432) * t1074 + t292 * t1239 - t1167 * t912, qJD(1) * t32 - qJD(4) * t81, qJD(1) * t31 - qJD(4) * t85, t2 * qJD(1) + (t777 + (t380 * t765 - t381 * t763) * t1131 + (t277 * t765 + t278 * t763) * t1132 + t63 * t1127 + t64 * t1129 + t145 * t873 + t219 * t882 + t146 * t867 + t220 * t880 - t788) * qJD(4) + t5 * qJD(5) + t9 * qJD(6) + (-t143 / 0.4e1 + t197 * t287 - t238 * t481 + t239 * t479 + t295 * t319 + t339 * t480 + t340 * t478) * t1073 + (t336 * t395 - t384 * t563 + t385 * t561 - t421 * t425 + t468 * t562 + t1240 * t560 - t227 / 0.4e1) * t1076 + t829, t11 * qJD(1) + t5 * qJD(4) + ((t295 * t383 + t339 * t419 + t340 * t420) * m(7) + (t1240 * t498 - t425 * t457 + t468 * t497) * m(6) + t386 * t1132 + (t111 * t1130 + t112 * t1127 + (-t293 * t763 + t294 * t765) * t1132) * t745 + t900) * qJD(5) + t19 * qJD(6), t14 * qJD(1) + t9 * qJD(4) + t19 * qJD(5) + t900 * qJD(6) + (-t161 / 0.4e1 + t785 + t847) * t1071; t780 * qJD(1) + t131 * qJD(2) + t130 * qJD(3) + t6 * qJD(4) + t15 * qJD(5) + qJD(6) * t846 + t1228 * t262, qJD(1) * t131 - qJD(4) * t199, qJD(1) * t130 + t226, t6 * qJD(1) - t910 + t777 * qJD(4) + t10 * qJD(5) + t16 * qJD(6) + (t909 / 0.4e1 + (t287 * t318 + t319 * t412 - t366 * t481 + t367 * t479 + t452 * t480 + t1237 * t478 - t168 / 0.4e1) * qJD(4)) * m(7), t15 * qJD(1) + t10 * qJD(4) + t900 * qJD(5) + t20 + (t383 * t412 + t419 * t452 + t420 * t1237 - t97 / 0.4e1 + t847) * t1072, qJD(1) * t846 + qJD(4) * t16 + qJD(5) * t21 + t20;];
Cq  = t17;