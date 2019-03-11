% Calculate matrix of centrifugal and coriolis load on the joints for
% S6RPRPRP10
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d3,d5]';
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
% Datum: 2019-03-09 03:32
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Cq = S6RPRPRP10_coriolismatJ_fixb_slag_vp1(qJ, qJD, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(8,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRPRP10_coriolismatJ_fixb_slag_vp1: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RPRPRP10_coriolismatJ_fixb_slag_vp1: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S6RPRPRP10_coriolismatJ_fixb_slag_vp1: pkin has to be [8x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RPRPRP10_coriolismatJ_fixb_slag_vp1: m has to be [7x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [7,3]), ...
  'S6RPRPRP10_coriolismatJ_fixb_slag_vp1: rSges has to be [7x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [7 6]), ...
  'S6RPRPRP10_coriolismatJ_fixb_slag_vp1: Icges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From coriolismat_joint_fixb_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 03:31:00
% EndTime: 2019-03-09 03:31:55
% DurationCPUTime: 48.94s
% Computational Cost: add. (42412->1133), mult. (105732->1529), div. (0->0), fcn. (114491->6), ass. (0->655)
t732 = sin(qJ(3));
t736 = cos(qJ(1));
t965 = t732 * t736;
t718 = pkin(8) * t965;
t733 = sin(qJ(1));
t721 = t736 * qJ(2);
t735 = cos(qJ(3));
t960 = t735 * t736;
t901 = pkin(3) * t965 - qJ(4) * t960;
t874 = t721 + t901;
t1121 = pkin(1) + pkin(7);
t891 = pkin(4) + t1121;
t753 = -t733 * t891 + t718 + t874;
t1088 = rSges(7,1) + pkin(5);
t1221 = rSges(7,3) + qJ(6);
t731 = sin(qJ(5));
t734 = cos(qJ(5));
t637 = t731 * t733 - t734 * t960;
t964 = t733 * t734;
t638 = t731 * t960 + t964;
t715 = rSges(7,2) * t965;
t927 = -t1088 * t638 - t1221 * t637 + t715;
t1260 = t753 + t927;
t1087 = rSges(7,2) + pkin(8);
t962 = t734 * t735;
t639 = t731 * t736 + t733 * t962;
t961 = t734 * t736;
t963 = t733 * t735;
t640 = -t731 * t963 + t961;
t757 = t1088 * t640 + t1221 * t639;
t703 = qJ(4) * t963;
t781 = t736 * t891 - t703;
t342 = (qJ(2) + (pkin(3) + t1087) * t732) * t733 + t757 + t781;
t1261 = t1260 * t637 - t342 * t639;
t1265 = m(7) * qJD(1) * t1261;
t532 = -rSges(6,1) * t637 - rSges(6,2) * t638;
t537 = -rSges(6,1) * t639 - rSges(6,2) * t640;
t793 = t733 * t532 + t537 * t736;
t1048 = m(6) * t793;
t1093 = t735 / 0.2e1;
t1181 = m(7) * t735;
t1223 = -t1181 / 0.2e1;
t415 = -t1088 * t639 + t1221 * t640;
t926 = -t1088 * t637 + t1221 * t638;
t1259 = t415 * t736 + t733 * t926;
t946 = t1048 * t1093 - t1223 * t1259;
t1183 = m(6) * t735;
t829 = rSges(7,1) * t731 - rSges(7,3) * t734;
t1142 = -rSges(7,2) * t735 - t732 * t829;
t827 = pkin(5) * t731 - qJ(6) * t734;
t910 = t827 * t732 - t1142;
t847 = t910 * t736;
t1243 = -t732 * t847 + t735 * t927;
t967 = t732 * t733;
t1233 = -rSges(7,2) * t967 - t757;
t846 = t910 * t733;
t320 = t1233 * t735 + t732 * t846;
t192 = t1243 * t736 - t320 * t733;
t714 = rSges(6,3) * t965;
t504 = t638 * rSges(6,1) - t637 * rSges(6,2) - t714;
t830 = rSges(6,1) * t731 + rSges(6,2) * t734;
t758 = t830 * t732;
t610 = rSges(6,3) * t735 + t758;
t1228 = -t735 * t504 - t610 * t965;
t831 = t640 * rSges(6,1) - t639 * rSges(6,2);
t752 = rSges(6,3) * t967 + t831;
t408 = t610 * t967 - t735 * t752;
t297 = t1228 * t736 - t408 * t733;
t949 = t192 * t1223 - t297 * t1183 / 0.2e1;
t19 = t949 - t946;
t1264 = qJD(1) * t19;
t1123 = m(7) / 0.2e1;
t945 = -t1259 * t1123 - t1048 / 0.2e1;
t1126 = m(6) / 0.2e1;
t948 = t1123 * t192 + t1126 * t297;
t27 = t948 - t945;
t1263 = qJD(1) * t27;
t1128 = m(5) / 0.2e1;
t695 = pkin(3) * t735 + qJ(4) * t732;
t671 = t733 * t695;
t828 = rSges(5,2) * t735 - rSges(5,3) * t732;
t564 = -t733 * t828 + t671;
t566 = (t695 - t828) * t736;
t1163 = -t564 * t736 + t566 * t733;
t716 = pkin(8) * t963;
t905 = t671 + t716;
t975 = t610 * t733;
t478 = t905 + t975;
t850 = pkin(8) * t735 + t695;
t480 = (t610 + t850) * t736;
t1164 = -t478 * t736 + t480 * t733;
t429 = t846 + t905;
t431 = (t850 + t910) * t736;
t1168 = -t429 * t736 + t431 * t733;
t692 = rSges(5,2) * t732 + rSges(5,3) * t735;
t890 = rSges(5,1) + t1121;
t463 = -t692 * t736 - t733 * t890 + t874;
t1017 = pkin(3) * t732;
t902 = -rSges(5,2) * t967 - rSges(5,3) * t963;
t464 = -t703 + (qJ(2) + t1017) * t733 + t890 * t736 + t902;
t835 = t463 * t967 - t464 * t965;
t1155 = t753 - t504;
t1086 = rSges(6,3) + pkin(8);
t383 = (qJ(2) + (pkin(3) + t1086) * t732) * t733 + t781 + t831;
t837 = t1155 * t967 - t383 * t965;
t334 = t1260 * t967;
t838 = -t342 * t965 + t334;
t879 = (t1168 * t735 + t838) * t1123 + (t1164 * t735 + t837) * t1126 + (t1163 * t735 + t835) * t1128;
t665 = pkin(3) * t963 + qJ(4) * t967;
t543 = -rSges(5,2) * t963 + rSges(5,3) * t967 + t665;
t900 = -pkin(3) * t960 - qJ(4) * t965;
t544 = -t736 * t828 - t900;
t792 = t543 * t736 - t733 * t544;
t883 = t732 * t964;
t884 = t731 * t967;
t560 = rSges(6,1) * t884 + rSges(6,2) * t883 + rSges(6,3) * t963;
t875 = t716 + t665;
t461 = t875 + t560;
t462 = (t1086 * t735 + t758) * t736 - t900;
t800 = t461 * t736 - t733 * t462;
t1160 = rSges(7,2) * t963 + t1088 * t884;
t849 = t1221 * t734;
t417 = -t849 * t967 + t1160 + t875;
t418 = (t1087 * t735 + (t1088 * t731 - t849) * t732) * t736 - t900;
t803 = t417 * t736 - t733 * t418;
t880 = (t735 * t803 + t838) * t1123 + (t735 * t800 + t837) * t1126 + (t735 * t792 + t835) * t1128;
t1196 = t879 - t880;
t1262 = t1196 * qJD(1);
t897 = qJD(1) * t735;
t1170 = -t1260 * t736 - t342 * t733;
t520 = -Icges(6,5) * t637 - Icges(6,6) * t638;
t522 = -Icges(7,4) * t637 + Icges(7,6) * t638;
t1258 = t520 + t522;
t521 = -Icges(6,5) * t639 - Icges(6,6) * t640;
t523 = -Icges(7,4) * t639 + Icges(7,6) * t640;
t1257 = t521 + t523;
t628 = Icges(6,4) * t639;
t502 = Icges(6,1) * t640 + Icges(6,5) * t967 - t628;
t929 = -Icges(6,2) * t640 + t502 - t628;
t997 = Icges(7,5) * t639;
t499 = Icges(7,1) * t640 + Icges(7,4) * t967 + t997;
t931 = -Icges(7,3) * t640 + t499 + t997;
t1256 = -t929 - t931;
t626 = Icges(6,4) * t637;
t501 = -Icges(6,1) * t638 + Icges(6,5) * t965 + t626;
t930 = -Icges(6,2) * t638 - t501 - t626;
t623 = Icges(7,5) * t637;
t497 = Icges(7,1) * t638 - Icges(7,4) * t965 + t623;
t932 = -Icges(7,3) * t638 + t497 + t623;
t1255 = -t930 - t932;
t1001 = Icges(6,4) * t640;
t496 = -Icges(6,2) * t639 + Icges(6,6) * t967 + t1001;
t933 = -Icges(6,1) * t639 - t1001 - t496;
t625 = Icges(7,5) * t640;
t487 = Icges(7,6) * t967 + Icges(7,3) * t639 + t625;
t935 = -Icges(7,1) * t639 + t487 + t625;
t1254 = t933 + t935;
t627 = Icges(6,4) * t638;
t494 = -Icges(6,2) * t637 - Icges(6,6) * t965 + t627;
t934 = -Icges(6,1) * t637 - t494 - t627;
t624 = Icges(7,5) * t638;
t486 = Icges(7,6) * t965 - Icges(7,3) * t637 - t624;
t936 = -Icges(7,1) * t637 - t486 + t624;
t1253 = t934 + t936;
t817 = Icges(6,5) * t731 + Icges(6,6) * t734;
t587 = Icges(6,3) * t735 + t732 * t817;
t1000 = Icges(6,4) * t731;
t821 = Icges(6,2) * t734 + t1000;
t593 = Icges(6,6) * t735 + t732 * t821;
t966 = t732 * t734;
t707 = Icges(6,4) * t966;
t968 = t731 * t732;
t998 = Icges(6,5) * t735;
t599 = Icges(6,1) * t968 + t707 + t998;
t364 = t587 * t967 - t593 * t639 + t599 * t640;
t1226 = t639 * t494 + t640 * t501;
t489 = -Icges(6,5) * t638 + Icges(6,6) * t637 + Icges(6,3) * t965;
t259 = -t489 * t967 - t1226;
t490 = Icges(6,5) * t640 - Icges(6,6) * t639 + Icges(6,3) * t967;
t260 = t490 * t967 - t639 * t496 + t640 * t502;
t808 = -t259 * t736 + t260 * t733;
t1244 = t364 * t735 + t732 * t808;
t706 = Icges(7,5) * t968;
t991 = Icges(7,6) * t735;
t585 = -Icges(7,3) * t966 + t706 + t991;
t820 = Icges(7,4) * t731 - Icges(7,6) * t734;
t591 = Icges(7,2) * t735 + t732 * t820;
t996 = Icges(7,5) * t734;
t824 = Icges(7,1) * t731 - t996;
t597 = Icges(7,4) * t735 + t732 * t824;
t363 = t585 * t639 + t591 * t967 + t597 * t640;
t1227 = t639 * t486 - t640 * t497;
t492 = -Icges(7,4) * t638 + Icges(7,2) * t965 - Icges(7,6) * t637;
t257 = -t492 * t967 - t1227;
t493 = Icges(7,4) * t640 + Icges(7,2) * t967 + Icges(7,6) * t639;
t258 = t639 * t487 + t493 * t967 + t640 * t499;
t809 = -t257 * t736 + t258 * t733;
t1245 = t363 * t735 + t732 * t809;
t865 = t1244 / 0.2e1 + t1245 / 0.2e1;
t1124 = -m(7) / 0.2e1;
t1127 = -m(6) / 0.2e1;
t1224 = -m(5) / 0.2e1;
t697 = rSges(4,1) * t735 - rSges(4,2) * t732;
t666 = t697 * t733;
t669 = t697 * t736;
t842 = t803 * t1124 + t800 * t1127 + t792 * t1224 + m(4) * (-t666 * t736 + t733 * t669) / 0.2e1;
t877 = t1124 * t1168 + t1127 * t1164 + t1163 * t1224;
t66 = t877 - t842;
t1249 = qJD(1) * t66;
t798 = t486 * t734 + t497 * t731;
t980 = t492 * t735;
t302 = t732 * t798 - t980;
t796 = t494 * t734 - t501 * t731;
t982 = t489 * t735;
t305 = t732 * t796 - t982;
t1248 = -t302 - t305;
t797 = -t487 * t734 + t499 * t731;
t979 = t493 * t735;
t304 = t732 * t797 + t979;
t795 = t496 * t734 + t502 * t731;
t981 = t490 * t735;
t307 = t732 * t795 + t981;
t1247 = t304 + t307;
t1246 = (-Icges(5,4) + Icges(4,5)) * t735 + (Icges(5,5) - Icges(4,6)) * t732;
t1242 = t259 * t733 + t260 * t736;
t1241 = t257 * t733 + t258 * t736;
t1159 = -t637 * t733 - t639 * t736;
t1222 = m(7) * t1159;
t483 = -t1222 / 0.2e1;
t481 = t1222 / 0.2e1;
t729 = t733 ^ 2;
t730 = t736 ^ 2;
t899 = t729 + t730;
t1240 = t732 * t899;
t1239 = t1253 * t638 + t1255 * t637 - t1258 * t965;
t1238 = t1254 * t638 + t1256 * t637 - t1257 * t965;
t1237 = t1253 * t640 + t1255 * t639 + t1258 * t967;
t1236 = t1254 * t640 + t1256 * t639 + t1257 * t967;
t649 = (Icges(7,4) * t734 + Icges(7,6) * t731) * t732;
t645 = (Icges(7,3) * t731 + t996) * t732;
t918 = t597 - t645;
t655 = Icges(7,1) * t966 + t706;
t922 = t585 + t655;
t277 = -t637 * t918 + t638 * t922 - t649 * t965;
t646 = (Icges(6,5) * t734 - Icges(6,6) * t731) * t732;
t650 = -Icges(6,2) * t968 + t707;
t917 = t599 + t650;
t656 = (Icges(6,1) * t734 - t1000) * t732;
t921 = -t593 + t656;
t278 = -t637 * t917 + t638 * t921 - t646 * t965;
t1235 = t277 + t278;
t279 = -t639 * t918 + t640 * t922 + t649 * t967;
t280 = -t639 * t917 + t640 * t921 + t646 * t967;
t1234 = t280 + t279;
t1003 = Icges(4,4) * t732;
t823 = Icges(4,2) * t735 + t1003;
t596 = -Icges(4,6) * t733 + t736 * t823;
t1002 = Icges(4,4) * t735;
t826 = Icges(4,1) * t732 + t1002;
t602 = -Icges(4,5) * t733 + t736 * t826;
t723 = Icges(5,6) * t732;
t990 = Icges(5,3) * t735;
t813 = t723 + t990;
t603 = Icges(5,5) * t733 + t736 * t813;
t993 = Icges(5,6) * t735;
t680 = Icges(5,2) * t732 + t993;
t605 = Icges(5,4) * t733 + t680 * t736;
t785 = t603 * t735 + t605 * t732;
t1229 = (-t596 * t735 - t602 * t732 - t785) * t736;
t595 = Icges(4,6) * t736 + t733 * t823;
t708 = Icges(4,4) * t963;
t601 = Icges(4,1) * t967 + Icges(4,5) * t736 + t708;
t604 = Icges(5,5) * t736 - t733 * t813;
t606 = Icges(5,4) * t736 - t680 * t733;
t814 = -Icges(5,3) * t732 + t993;
t641 = t814 * t733;
t642 = t814 * t736;
t995 = Icges(5,2) * t735;
t815 = -t723 + t995;
t643 = t815 * t733;
t644 = t815 * t736;
t653 = -Icges(4,2) * t967 + t708;
t685 = -Icges(4,2) * t732 + t1002;
t654 = t685 * t736;
t687 = Icges(4,1) * t735 - t1003;
t657 = t687 * t733;
t658 = t687 * t736;
t1225 = -((-t601 - t653 + t606 - t641) * t736 + (t602 + t654 + t605 + t642) * t733) * t735 + ((-t595 + t657 + t604 + t643) * t736 + (t596 - t658 + t603 - t644) * t733) * t732;
t253 = -t637 * t486 + t492 * t965 + t638 * t497;
t255 = t489 * t965 - t637 * t494 - t638 * t501;
t1098 = -t732 / 0.2e1;
t1156 = t1159 * t1181;
t458 = -t1156 / 0.2e1;
t460 = t1156 / 0.2e1;
t1074 = m(5) * (t463 * t736 + t464 * t733);
t1220 = t253 * t733;
t1219 = t253 * t736;
t1218 = t255 * t733;
t1217 = t255 * t736;
t361 = -t637 * t585 + t591 * t965 - t638 * t597;
t1216 = t361 * t735;
t362 = t587 * t965 + t637 * t593 - t638 * t599;
t1215 = t362 * t735;
t1212 = t610 * t736;
t816 = Icges(7,5) * t731 - Icges(7,3) * t734;
t1144 = t732 * t816 + t991;
t548 = t1144 * t736;
t556 = t597 * t736;
t768 = t591 * t736 - t798;
t213 = -t768 * t735 + (t548 * t734 - t556 * t731 + t492) * t732;
t554 = t593 * t736;
t825 = Icges(6,1) * t731 + Icges(6,4) * t734;
t1143 = t732 * t825 + t998;
t558 = t1143 * t736;
t766 = -t587 * t736 + t796;
t215 = t766 * t735 + (-t554 * t734 - t558 * t731 + t489) * t732;
t1211 = -t213 - t215;
t547 = t1144 * t733;
t555 = t597 * t733;
t767 = -t591 * t733 - t797;
t214 = -t767 * t735 + (-t547 * t734 + t555 * t731 - t493) * t732;
t553 = t593 * t733;
t557 = t1143 * t733;
t765 = t587 * t733 + t795;
t216 = t765 * t735 + (t553 * t734 + t557 * t731 - t490) * t732;
t1210 = t214 + t216;
t228 = t522 * t735 + (t731 * t936 + t734 * t932) * t732;
t230 = t520 * t735 + (t731 * t934 + t734 * t930) * t732;
t1209 = -t228 - t230;
t229 = t523 * t735 + (t731 * t935 + t734 * t931) * t732;
t231 = t521 * t735 + (t731 * t933 + t734 * t929) * t732;
t1208 = t229 + t231;
t592 = -Icges(7,2) * t732 + t735 * t820;
t789 = -t585 * t734 + t597 * t731;
t764 = -t592 - t789;
t976 = t591 * t735;
t1149 = t732 * t764 - t976;
t586 = -Icges(7,6) * t732 + t735 * t816;
t598 = -Icges(7,4) * t732 + t735 * t824;
t239 = -t1149 * t733 + t586 * t639 + t598 * t640;
t588 = -Icges(6,3) * t732 + t735 * t817;
t788 = t593 * t734 + t599 * t731;
t763 = t588 + t788;
t977 = t587 * t735;
t1150 = t732 * t763 + t977;
t594 = -Icges(6,6) * t732 + t735 * t821;
t600 = -Icges(6,5) * t732 + t735 * t825;
t240 = t1150 * t733 - t594 * t639 + t600 * t640;
t1207 = t239 + t240;
t241 = t1149 * t736 + t637 * t586 + t638 * t598;
t242 = -t1150 * t736 - t637 * t594 + t638 * t600;
t1206 = t241 + t242;
t295 = -t764 * t735 + (-t734 * t586 + t731 * t598 - t591) * t732;
t296 = t763 * t735 + (t734 * t594 + t731 * t600 - t587) * t732;
t1205 = t295 + t296;
t1204 = -t361 - t362;
t1203 = t363 + t364;
t392 = t732 * t789 + t976;
t393 = t732 * t788 + t977;
t1202 = t392 + t393;
t683 = Icges(5,4) * t732 + Icges(5,5) * t735;
t607 = Icges(5,1) * t733 + t683 * t736;
t579 = t736 * t607;
t818 = Icges(4,5) * t732 + Icges(4,6) * t735;
t590 = -Icges(4,3) * t733 + t736 * t818;
t1201 = -t736 * t590 - t596 * t963 - t602 * t967 - t733 * t785 + t579;
t787 = -t595 * t735 - t601 * t732;
t1200 = t736 * t787;
t1199 = t1246 * t733;
t1198 = t1246 * t736;
t1197 = t1247 * t733 + t1248 * t736;
t739 = t757 * t736 + (t715 - t927) * t733;
t923 = t1142 * t736 - t827 * t965;
t924 = -t1221 * t883 + t1160;
t171 = (t733 * t923 + t736 * t924) * t732 + t739 * t735;
t244 = t739 * t732;
t262 = (-t1212 * t733 + t560 * t736) * t732 + (t733 * t504 + t736 * t752) * t735;
t313 = t1243 * t967;
t365 = (t736 * t831 + (t504 + t714) * t733) * t732;
t389 = t1228 * t967;
t1008 = rSges(6,3) * t732;
t612 = t735 * t830 - t1008;
t329 = (t560 - t975) * t735 + ((-t612 - t1008) * t733 - t831) * t732;
t330 = (-t612 * t736 + t504) * t732;
t805 = -t329 * t736 + t330 * t733;
t1010 = rSges(7,2) * t732;
t909 = t1010 + (-t827 - t829) * t735;
t222 = (t924 - t846) * t735 + ((t909 - t1010) * t733 - t757) * t732;
t223 = (-t847 - t923) * t735 + (t736 * t909 - t927) * t732;
t812 = -t222 * t736 + t223 * t733;
t987 = t408 * t736;
t1015 = (t313 + (t320 * t736 + t171) * t732 + (t244 - t812) * t735) * t1123 + (t389 + (t262 + t987) * t732 + (t365 - t805) * t735) * t1126;
t1171 = t899 * t735;
t286 = -t415 * t733 + t736 * t926;
t388 = t532 * t736 - t733 * t537;
t663 = (rSges(6,1) * t734 - rSges(6,2) * t731) * t732;
t907 = (t1088 * t734 + t1221 * t731) * t732;
t514 = t907 * t733;
t516 = t907 * t736;
t794 = t514 * t733 + t516 * t736;
t947 = (t732 * t286 - t735 * t794) * t1123 + (-t1171 * t663 + t388 * t732) * t1126;
t1193 = t947 - t1015;
t580 = t736 * (Icges(5,1) * t736 - t733 * t683);
t1190 = -t1229 + t580 + (-t590 + t607) * t733;
t1169 = -t1155 * t736 - t383 * t733;
t1096 = t733 / 0.2e1;
t1187 = -t735 / 0.2e1;
t1186 = -t736 / 0.2e1;
t1090 = t736 / 0.2e1;
t1180 = t1235 * t735 + (t1238 * t733 - t1239 * t736) * t732;
t1179 = t1234 * t735 + (t1236 * t733 - t1237 * t736) * t732;
t1122 = m(7) / 0.4e1;
t1125 = m(6) / 0.4e1;
t1172 = (m(5) / 0.4e1 + t1125 + t1122) * (0.1e1 - t899) * t735 * t732;
t1167 = -t637 * t487 - t638 * t499;
t1166 = t637 * t496 - t638 * t502;
t1158 = t216 / 0.2e1 + t214 / 0.2e1;
t1157 = t215 / 0.2e1 + t213 / 0.2e1;
t958 = (t1171 * t244 + t320 * t965 + t313) * t1123 + (t1171 * t365 + t408 * t965 + t389) * t1126;
t1153 = (-t478 * t532 - t480 * t537 + (t1155 * t733 - t383 * t736) * t663) * t1127 + (t1260 * t514 - t342 * t516 - t415 * t431 - t429 * t926) * t1124;
t1152 = (t1243 * t418 + t1260 * t223 + t222 * t342 - t320 * t417) * t1124 + (t1155 * t330 + t1228 * t462 + t329 * t383 - t408 * t461) * t1127;
t1148 = t732 * t765 + t981;
t1147 = t732 * t766 - t982;
t1146 = t732 * t767 - t979;
t1145 = t732 * t768 + t980;
t856 = t599 / 0.2e1 + t597 / 0.2e1;
t858 = t593 / 0.2e1 - t585 / 0.2e1;
t1141 = t731 * (t656 / 0.2e1 + t655 / 0.2e1 - t858) + t734 * (t650 / 0.2e1 - t645 / 0.2e1 + t856);
t1136 = t731 * t856 + t734 * t858 + t588 / 0.2e1 + t592 / 0.2e1 - t814 / 0.2e1 - t680 / 0.2e1 - t685 / 0.2e1 - t826 / 0.2e1;
t1135 = t731 * (t600 / 0.2e1 + t598 / 0.2e1) + (t594 / 0.2e1 - t586 / 0.2e1) * t734 - t587 / 0.2e1 - t591 / 0.2e1 + t723 + t990 / 0.2e1 - t995 / 0.2e1 + t823 / 0.2e1 - t687 / 0.2e1;
t1133 = 0.4e1 * qJD(1);
t1132 = 2 * qJD(3);
t1131 = 4 * qJD(3);
t1130 = 2 * qJD(5);
t1129 = 4 * qJD(5);
t1118 = m(6) * (t1228 * t330 + t262 * t365 - t329 * t408);
t881 = t732 * t961;
t312 = t320 * t881;
t1115 = m(7) * (t222 * t637 + t223 * t639 - t312 + (-t244 * t735 + (-t1243 * t733 - t171) * t732) * t734);
t1113 = m(7) * (t1243 * t223 + t171 * t244 - t222 * t320);
t664 = pkin(3) * t967 - t703;
t906 = -pkin(4) * t736 - pkin(8) * t967 - t664;
t619 = t736 * t901;
t908 = t736 * (pkin(4) * t733 - t718) - t619;
t220 = -t927 * t736 + (t906 + t1233) * t733 + t908;
t473 = t1159 * t732;
t972 = t639 * t733;
t973 = t637 * t736;
t513 = -t972 + t973;
t728 = t732 ^ 2;
t545 = t639 * t735 + t728 * t964;
t546 = -t735 * t637 + t728 * t961;
t1111 = m(7) * (-t1243 * t883 - t220 * t473 + t244 * t513 + t429 * t546 - t431 * t545 - t312);
t1106 = m(7) * (t429 * t640 - t431 * t638 + t514 * t639 - t516 * t637 + (t220 * t731 - t286 * t734) * t732);
t1105 = m(7) * (t220 * t286 + t429 * t514 + t431 * t516);
t1104 = m(7) * (-t1243 * t637 + t1260 * t546 - t639 * t320 + t342 * t545);
t1095 = t733 / 0.4e1;
t1089 = t736 / 0.4e1;
t1085 = m(3) * ((rSges(3,3) * t736 + t721) * t736 + (rSges(3,3) + qJ(2)) * t729);
t832 = rSges(4,1) * t732 + rSges(4,2) * t735;
t754 = -t733 * rSges(4,3) + t736 * t832;
t540 = -t1121 * t733 + t721 + t754;
t541 = (rSges(4,3) + t1121) * t736 + (qJ(2) + t832) * t733;
t1084 = m(4) * (t540 * t669 + t541 * t666);
t1083 = m(4) * (t540 * t736 + t541 * t733);
t379 = t692 * t730 - t619 + (-t664 - t902) * t733;
t833 = t564 * t967 + t566 * t965;
t1078 = m(5) * (t1171 * t379 + t833);
t1077 = m(5) * (t463 * t544 + t464 * t543);
t1075 = t735 * t1074;
t271 = t736 * t504 + (-t752 + t906) * t733 + t908;
t1067 = m(6) * (t271 * t388 + (t478 * t733 + t480 * t736) * t663);
t834 = t478 * t967 + t480 * t965;
t1062 = m(6) * (t1171 * t271 + t834);
t1059 = m(6) * (t1155 * t462 + t383 * t461);
t1058 = m(6) * (-t1155 * t532 + t383 * t537);
t1057 = t1169 * t1183;
t1056 = m(6) * t1169;
t1044 = m(7) * (t1260 * t640 + t342 * t638 + t415 * t637 - t639 * t926);
t779 = -t334 * t734 + t342 * t881;
t1043 = m(7) * (t417 * t637 + t418 * t639 + t779);
t1042 = m(7) * (-t637 * t429 - t639 * t431 + t779);
t836 = t429 * t967 + t431 * t965;
t1040 = m(7) * (t1171 * t220 + t836);
t1037 = m(7) * (-t1260 * t926 + t342 * t415);
t1036 = m(7) * (t1260 * t418 + t342 * t417);
t1032 = t1170 * t1181;
t1031 = m(7) * t1170;
t582 = t639 * t967;
t904 = t962 * t1240;
t1025 = m(7) * (t582 + (-0.2e1 * t962 - t973) * t732 + t904);
t1024 = m(7) * (t513 * t732 + t904);
t1023 = m(7) * (t582 + (-t1171 * t734 - t973) * t732);
t1016 = pkin(8) * t732;
t1013 = m(7) * qJD(3);
t1012 = m(7) * qJD(5);
t1011 = m(7) * qJD(6);
t740 = t1123 * t812 + t1126 * t805;
t741 = t1126 * t663 * t899 + t1123 * t794;
t70 = -t740 + t741;
t970 = t70 * qJD(2);
t969 = t728 * t731;
t956 = (t492 * t733 + t493 * t736) * t732 + t257 + t1167 + t1227;
t954 = (t489 * t733 + t490 * t736) * t732 + t259 + t1166 + t1226;
t898 = qJD(1) * t732;
t895 = qJD(3) * t732;
t894 = qJD(3) * t735;
t893 = qJD(5) * t732;
t892 = qJD(5) * t735;
t181 = -t1145 * t733 - t548 * t639 - t556 * t640;
t182 = -t1146 * t733 + t547 * t639 + t555 * t640;
t183 = t1147 * t733 + t554 * t639 - t558 * t640;
t184 = t1148 * t733 - t553 * t639 + t557 * t640;
t888 = ((-t181 - t183) * t736 + (t182 + t184) * t733 - t1203) * t732 / 0.2e1 + (t808 + t809 + t1207) * t1093;
t185 = t1145 * t736 - t637 * t548 - t638 * t556;
t186 = t1146 * t736 + t637 * t547 + t638 * t555;
t187 = -t1147 * t736 + t637 * t554 - t638 * t558;
t188 = -t1148 * t736 - t637 * t553 + t638 * t557;
t256 = -t490 * t965 - t1166;
t810 = t733 * t256 - t1217;
t254 = -t493 * t965 - t1167;
t811 = t733 * t254 - t1219;
t887 = (t811 + t810 + t1206) * t1187 + ((-t185 - t187) * t736 + (t186 + t188) * t733 - t1204) * t1098;
t886 = (t1197 + t1205) * t1187 + (t1210 * t733 + t1211 * t736 - t1202) * t1098;
t394 = t736 * (Icges(4,3) * t736 + t733 * t818) + t595 * t963 + t601 * t967;
t873 = t967 / 0.4e1;
t872 = -t965 / 0.4e1;
t868 = t1090 * t1238 + t1096 * t1239;
t867 = t1090 * t1236 + t1096 * t1237;
t122 = t732 * t811 - t1216;
t123 = t732 * t810 - t1215;
t866 = t122 / 0.2e1 + t123 / 0.2e1;
t864 = t1098 * t1197 + t1187 * t1202;
t153 = t254 * t736 + t1220;
t154 = t256 * t736 + t1218;
t863 = -t153 / 0.2e1 - t154 / 0.2e1;
t861 = t228 / 0.2e1 + t230 / 0.2e1;
t860 = t229 / 0.2e1 + t231 / 0.2e1;
t854 = t646 / 0.2e1 + t649 / 0.2e1;
t853 = t683 / 0.2e1 - t818 / 0.2e1;
t844 = t899 * t832;
t335 = t415 * t735 - t907 * t967;
t336 = -t516 * t732 - t735 * t926;
t804 = t335 * t736 - t336 * t733;
t691 = qJ(4) * t735 - t1017;
t670 = t733 * t691;
t428 = t670 + (-t909 - t1016) * t733;
t430 = t718 + (-t691 + t909) * t736;
t802 = t428 * t733 - t430 * t736;
t434 = t537 * t735 - t663 * t967;
t435 = -t735 * t532 - t663 * t965;
t801 = t434 * t736 - t435 * t733;
t477 = t670 + (t612 - t1016) * t733;
t479 = t718 + (-t612 - t691) * t736;
t799 = t477 * t733 - t479 * t736;
t791 = t545 * t736 - t546 * t733;
t563 = t692 * t733 + t670;
t565 = (-t691 - t692) * t736;
t790 = t563 * t733 - t565 * t736;
t784 = -t604 * t735 - t606 * t732;
t782 = t638 * t736 - t640 * t733;
t745 = (-t473 * t732 + t735 * t791) * t1123;
t746 = m(7) * (t735 * t782 + t969);
t249 = t745 - t746 / 0.2e1;
t755 = t791 * t1124;
t759 = m(7) * t782;
t346 = t755 + t759 / 0.2e1;
t780 = t346 * qJD(2) + t249 * qJD(4);
t620 = t736 * t900;
t778 = -pkin(8) * t1171 + t620;
t776 = -t868 - t887;
t775 = t867 - t888;
t45 = t1216 + (t733 * t956 + t1219) * t732;
t46 = t1215 + (t733 * t954 + t1217) * t732;
t774 = t46 / 0.2e1 + t45 / 0.2e1 + t866;
t762 = -t607 + t784;
t747 = -t1153 + (-t1209 + t1235) * t1095 + (t1208 + t1234) * t1089;
t401 = t733 * t784 + t580;
t57 = t736 * t956 - t1220;
t58 = t736 * t954 - t1218;
t744 = t58 / 0.2e1 + t57 / 0.2e1 + t729 * t590 / 0.2e1 - t863 + (t733 * t762 + t1190 - t401) * t1096 + ((t590 - t787) * t736 - t579 + t1200 + t1201) * t1090;
t743 = -t1201 * t733 / 0.2e1 + (t394 + t401) * t1186 + (t736 * t762 - t1200 + t1201) * t1096 + ((t590 + t787) * t733 + t394 + t1190 + t1229) * t1090;
t738 = (t153 + t154 + t57 + t58) * t873 + (t122 + t123 + t45 + t46) * t1089 + (t872 + t965 / 0.4e1) * (t1241 + t1242) + (-t733 / 0.4e1 + t1095) * (t1245 + t1244);
t737 = -t1152 + t1202 * t1098 + t1205 * t1093 + (t1207 + t1210) * t873 + (t1206 - t1211) * t872 + (t1203 + t1247) * t963 / 0.4e1 - (t1204 - t1248) * t960 / 0.4e1;
t621 = t899 * t966;
t574 = t637 * t881;
t451 = 0.2e1 * (t1128 + t1126 + t1123) * t1240;
t448 = t574 + (t962 - t972) * t966;
t420 = t1023 / 0.2e1;
t413 = t1024 / 0.2e1;
t412 = t431 * t881;
t406 = t637 * t638 + t639 * t640 - t734 * t969;
t405 = -t543 * t733 + t828 * t730 + t620;
t386 = t1025 / 0.2e1;
t376 = t793 * t732;
t373 = 0.4e1 * t1172;
t345 = t755 - t759 / 0.2e1;
t344 = -t1212 * t736 + (-t560 - t665) * t733 + t778;
t340 = (t735 * t646 + (t731 * t921 + t734 * t917) * t732) * t735;
t339 = (t735 * t649 + (t731 * t922 + t734 * t918) * t732) * t735;
t328 = 0.2e1 * t483;
t327 = t481 + t483;
t326 = 0.2e1 * t481;
t310 = 0.2e1 * t460;
t309 = t458 + t460;
t308 = 0.2e1 * t458;
t269 = t923 * t736 + (-t665 - t924) * t733 + t778;
t261 = t1259 * t732;
t248 = t745 + t746 / 0.2e1;
t234 = t420 + t386 - t1024 / 0.2e1;
t233 = t413 + t420 - t1025 / 0.2e1;
t232 = t413 + t386 - t1023 / 0.2e1;
t133 = t1042 / 0.2e1;
t131 = t220 * t513 - t429 * t883 - t412;
t128 = t1043 / 0.2e1;
t126 = t1044 / 0.2e1;
t104 = t1032 + t1057 - t1075;
t101 = t1243 * t546 - t244 * t473 - t320 * t545;
t99 = t187 * t733 + t188 * t736;
t98 = t185 * t733 + t186 * t736;
t97 = t183 * t733 + t184 * t736;
t96 = t181 * t733 + t182 * t736;
t91 = t1104 / 0.2e1;
t89 = -t1031 - t1056 + t1074 + t1083 + t1085;
t84 = t1106 / 0.2e1;
t72 = t1040 + t1062 + t1078;
t71 = t1141 * t732 + t854 * t735 + t1037 + t1058;
t69 = t740 + t741;
t64 = t842 + t877;
t63 = t1111 / 0.2e1;
t32 = t1115 / 0.2e1;
t31 = t1135 * t732 + t1136 * t735 + t1036 + t1059 + t1077 + t1084;
t30 = t128 - t1042 / 0.2e1;
t29 = t133 + t128;
t28 = t133 - t1043 / 0.2e1;
t25 = t945 + t948;
t22 = t126 - t1104 / 0.2e1;
t21 = t91 + t126;
t20 = t91 - t1044 / 0.2e1;
t17 = t949 + t946;
t15 = t879 + t880;
t13 = t84 + t32 - t1111 / 0.2e1;
t12 = t63 + t84 - t1115 / 0.2e1;
t11 = t63 + t32 - t1106 / 0.2e1;
t10 = t958 - t1193;
t9 = t947 + t1015 - t958;
t8 = t958 + t1193;
t7 = t733 * t868 + t736 * t867 + t1067 + t1105;
t6 = t774 * t967;
t5 = t733 * t743 + t736 * t744;
t4 = t1118 + t1113 + (t733 * t865 - t736 * t866 - t886) * t735 + (t733 * t888 + t736 * t887 + t864) * t732;
t3 = t747 + (-t58 / 0.4e1 - t57 / 0.4e1 - t154 / 0.4e1 - t153 / 0.4e1) * t967 + t737 + (-t123 / 0.4e1 - t122 / 0.4e1 - t46 / 0.4e1 - t45 / 0.4e1) * t736;
t2 = t738 + t737 + (-t278 / 0.4e1 - t277 / 0.4e1 - t230 / 0.4e1 - t228 / 0.4e1) * t733 + (-t280 / 0.4e1 - t279 / 0.4e1 - t231 / 0.4e1 - t229 / 0.4e1) * t736 + t1153;
t1 = t747 + t738 + (-t296 / 0.2e1 - t295 / 0.2e1 + (-t362 / 0.4e1 - t361 / 0.4e1 + t305 / 0.4e1 + t302 / 0.4e1) * t736 + (-t364 / 0.4e1 - t363 / 0.4e1 - t307 / 0.4e1 - t304 / 0.4e1) * t733) * t735 + (t393 / 0.2e1 + t392 / 0.2e1 + (t242 / 0.4e1 + t241 / 0.4e1 + t215 / 0.4e1 + t213 / 0.4e1) * t736 + (-t216 / 0.4e1 - t214 / 0.4e1 - t240 / 0.4e1 - t239 / 0.4e1) * t733) * t732 + t1152;
t14 = [t89 * qJD(2) + t31 * qJD(3) + t104 * qJD(4) + t71 * qJD(5) - t1011 * t1261, qJD(1) * t89 + qJD(3) * t64 + qJD(5) * t25 + qJD(6) * t327, t31 * qJD(1) + t64 * qJD(2) + t15 * qJD(4) + t3 * qJD(5) + t29 * qJD(6) + ((t463 * t563 + t464 * t565 - t543 * t566 + t544 * t564) * t1128 + (t1155 * t477 + t383 * t479 - t461 * t480 + t462 * t478) * t1126 + (t1260 * t428 + t342 * t430 - t417 * t431 + t418 * t429) * t1123) * t1132 + ((m(4) * (t541 * t832 - t666 * t697) + t239 / 0.2e1 + t240 / 0.2e1 + t853 * t736 - t744 + t1158) * qJD(3) + (t604 / 0.2e1 + t643 / 0.2e1 - t595 / 0.2e1 + t657 / 0.2e1) * t894 + (t606 / 0.2e1 - t641 / 0.2e1 - t601 / 0.2e1 - t653 / 0.2e1) * t895) * t736 + ((t241 / 0.2e1 + t242 / 0.2e1 + m(4) * (-t540 * t832 + t669 * t697) + t853 * t733 - t743 + t1157) * qJD(3) + (t603 / 0.2e1 - t644 / 0.2e1 + t596 / 0.2e1 - t658 / 0.2e1) * t894 + (t605 / 0.2e1 + t642 / 0.2e1 + t602 / 0.2e1 + t654 / 0.2e1) * t895) * t733, qJD(1) * t104 + qJD(3) * t15 + qJD(5) * t17 + qJD(6) * t309, t71 * qJD(1) + t25 * qJD(2) + t3 * qJD(3) + t17 * qJD(4) + (t339 + t340) * qJD(5) + t21 * qJD(6) + ((-t1243 * t926 + t1260 * t336 - t320 * t415 + t335 * t342) * t1123 + (t1155 * t435 - t1228 * t532 + t383 * t434 - t408 * t537) * t1126) * t1130 + ((-t278 / 0.2e1 - t277 / 0.2e1 - t861) * t736 + (t280 / 0.2e1 + t279 / 0.2e1 - t774 + t860) * t733) * t893, t327 * qJD(2) + t29 * qJD(3) + t309 * qJD(4) + t21 * qJD(5) - t1265; -t66 * qJD(3) - t27 * qJD(5) + t326 * qJD(6) + (t1031 / 0.4e1 + t1056 / 0.4e1 - t1074 / 0.4e1 - t1083 / 0.4e1 - t1085 / 0.4e1) * t1133, 0, -t1249 + (t790 * t1128 + t799 * t1126 + t802 * t1123 - m(4) * t844 / 0.2e1) * t1132 + t451 * qJD(4) + t69 * qJD(5) - t621 * t1011, t451 * qJD(3), -t1263 + t69 * qJD(3) + (t1124 * t804 + t1127 * t801) * t1130 + t345 * qJD(6), t326 * qJD(1) + t345 * qJD(5) - t1013 * t621; t66 * qJD(2) + t5 * qJD(3) + t1196 * qJD(4) + t1 * qJD(5) + t28 * qJD(6) + (-t1084 / 0.4e1 - t1077 / 0.4e1 - t1059 / 0.4e1 - t1036 / 0.4e1) * t1133 - t1136 * t897 - t1135 * t898, qJD(5) * t70 + t1249, t5 * qJD(1) + t72 * qJD(4) + t7 * qJD(5) + t131 * t1011 + (m(7) * (t220 * t269 + t428 * t429 - t430 * t431) + m(6) * (t271 * t344 + t477 * t478 - t479 * t480) + m(5) * (t379 * t405 + t563 * t564 - t565 * t566) + m(4) * ((-t736 * t754 + (-t736 * rSges(4,3) - t733 * t832) * t733) * (-t733 * t666 - t669 * t736) - t697 * t844) + (t99 + t98 - t1198 * t729 + (t1199 * t733 - t1225) * t736) * t1096 + (t97 + t96 + t1199 * t730 + (-t1198 * t736 + t1225) * t733) * t1090) * qJD(3), t72 * qJD(3) - 0.4e1 * qJD(4) * t1172 + t8 * qJD(5) + t233 * qJD(6) + t1262, t1 * qJD(1) + t970 + t7 * qJD(3) + t8 * qJD(4) + (t1090 * t1179 + t1096 * t1180) * qJD(5) + t12 * qJD(6) + (-t1118 / 0.4e1 - t1113 / 0.4e1) * t1129 + ((t376 * t271 + t365 * t388 - t434 * t480 + t435 * t478 + (t1228 * t733 + t987) * t663) * t1126 + (t1243 * t514 + t220 * t261 + t244 * t286 + t320 * t516 - t335 * t431 + t336 * t429) * t1123) * t1130 + ((t860 + t866) * t736 + (t861 - t865) * t733 + t886) * t892 + (t733 * t775 + t736 * t776 - t864) * t893, t28 * qJD(1) + t131 * t1013 + t233 * qJD(4) + t12 * qJD(5) + (t574 + (-t513 - t972) * t966 - t448) * t1011; -t1196 * qJD(3) - t19 * qJD(5) + t308 * qJD(6) + (-t1032 / 0.4e1 - t1057 / 0.4e1 + t1075 / 0.4e1) * t1133, 0, -t1262 + t373 * qJD(4) + t9 * qJD(5) + t232 * qJD(6) + (-t1040 / 0.4e1 - t1062 / 0.4e1 - t1078 / 0.4e1) * t1131 + ((t732 * t269 + t836) * t1123 + (t732 * t344 + t834) * t1126 + (t732 * t405 + t833) * t1128 + ((t220 - t802) * t1123 + (t271 - t799) * t1126 + (t379 - t790) * t1128) * t735) * t1132, t373 * qJD(3), -t1264 + t9 * qJD(3) + ((t376 * t732 + t735 * t801) * t1126 + (t261 * t732 + t735 * t804) * t1123) * t1130 + t248 * qJD(6), qJD(1) * t308 + qJD(3) * t232 + qJD(5) * t248; t27 * qJD(2) + t2 * qJD(3) + t19 * qJD(4) + t6 * qJD(5) + t20 * qJD(6) - t854 * t897 + (-t1037 / 0.4e1 - t1058 / 0.4e1) * t1133 - t1141 * t898, -qJD(3) * t70 + qJD(6) * t346 + t1263, t2 * qJD(1) - t970 + t10 * qJD(4) + t4 * qJD(5) + t11 * qJD(6) + (-t1105 / 0.4e1 - t1067 / 0.4e1) * t1131 + ((t1243 * t428 + t171 * t220 - t222 * t431 + t223 * t429 + t244 * t269 - t320 * t430) * t1123 + (t1228 * t477 + t262 * t271 - t329 * t480 + t330 * t478 + t344 * t365 - t408 * t479) * t1126) * t1132 + ((-t304 / 0.2e1 - t307 / 0.2e1 - t99 / 0.2e1 - t98 / 0.2e1) * t736 + (t97 / 0.2e1 + t96 / 0.2e1 - t302 / 0.2e1 - t305 / 0.2e1) * t733) * t895 + (((t863 + t1158) * t735 - t775) * t736 + ((t1241 / 0.2e1 + t1242 / 0.2e1 + t1157) * t735 + t776) * t733) * qJD(3), qJD(3) * t10 + qJD(6) * t249 + t1264, t6 * qJD(1) + t4 * qJD(3) + (t339 / 0.2e1 + t340 / 0.2e1) * t892 + ((t1243 * t336 + t244 * t261 - t320 * t335) * t1122 + (t1228 * t435 + t365 * t376 - t408 * t434) * t1125) * t1129 + t101 * t1011 + (t1179 * t1096 + (t1208 * t733 + t1209 * t736) * t1093 + t1180 * t1186) * t893, t20 * qJD(1) + t11 * qJD(3) + t101 * t1012 + (t473 * t966 + t545 * t637 + t546 * t639 - t406) * t1011 + t780; t328 * qJD(2) + t30 * qJD(3) + t310 * qJD(4) + t22 * qJD(5) + t1265, qJD(1) * t328 - qJD(5) * t346, t30 * qJD(1) + (t428 * t639 + t430 * t637 - t412 + (-t735 * t220 + (-t429 * t733 - t269) * t732) * t734 - t131) * t1013 + t234 * qJD(4) + t13 * qJD(5) + t448 * t1011, qJD(1) * t310 + qJD(3) * t234 - qJD(5) * t249, t22 * qJD(1) + t13 * qJD(3) + (-t320 * t638 + t1243 * t640 + t335 * t637 + t336 * t639 + (t244 * t731 - t261 * t734) * t732 - t101) * t1012 + t406 * t1011 - t780, 0.4e1 * (t448 * qJD(3) / 0.4e1 + t406 * qJD(5) / 0.4e1) * m(7);];
Cq  = t14;
