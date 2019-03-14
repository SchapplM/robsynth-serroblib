% Calculate matrix of centrifugal and coriolis load on the joints for
% S6RPPRRR3
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d4,d5,d6,theta2]';
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
% Datum: 2019-03-09 02:24
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Cq = S6RPPRRR3_coriolismatJ_fixb_slag_vp1(qJ, qJD, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(10,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPPRRR3_coriolismatJ_fixb_slag_vp1: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RPPRRR3_coriolismatJ_fixb_slag_vp1: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RPPRRR3_coriolismatJ_fixb_slag_vp1: pkin has to be [10x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RPPRRR3_coriolismatJ_fixb_slag_vp1: m has to be [7x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [7,3]), ...
  'S6RPPRRR3_coriolismatJ_fixb_slag_vp1: rSges has to be [7x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [7 6]), ...
  'S6RPPRRR3_coriolismatJ_fixb_slag_vp1: Icges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From coriolismat_joint_fixb_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 02:23:02
% EndTime: 2019-03-09 02:23:44
% DurationCPUTime: 39.10s
% Computational Cost: add. (146200->1112), mult. (138800->1594), div. (0->0), fcn. (152143->10), ass. (0->665)
t1147 = -m(6) * qJD(1) / 0.4e1;
t1017 = m(7) * qJD(1);
t1139 = -t1017 / 0.4e1;
t1014 = m(7) * qJD(6);
t708 = qJ(1) + pkin(10);
t704 = sin(t708);
t705 = cos(t708);
t710 = qJ(5) + qJ(6);
t707 = cos(t710);
t706 = sin(t710);
t712 = sin(qJ(4));
t947 = t706 * t712;
t607 = -t704 * t947 + t705 * t707;
t943 = t707 * t712;
t608 = t704 * t943 + t705 * t706;
t715 = cos(qJ(4));
t959 = t704 * t715;
t480 = Icges(7,5) * t608 + Icges(7,6) * t607 - Icges(7,3) * t959;
t1001 = Icges(7,4) * t608;
t483 = Icges(7,2) * t607 - Icges(7,6) * t959 + t1001;
t594 = Icges(7,4) * t607;
t486 = Icges(7,1) * t608 - Icges(7,5) * t959 + t594;
t609 = t704 * t707 + t705 * t947;
t610 = -t704 * t706 + t705 * t943;
t950 = t705 * t715;
t285 = t480 * t950 + t609 * t483 - t610 * t486;
t1155 = t285 * t704;
t1154 = t285 * t705;
t714 = cos(qJ(5));
t711 = sin(qJ(5));
t939 = t711 * t712;
t628 = -t704 * t939 + t705 * t714;
t928 = t712 * t714;
t953 = t705 * t711;
t629 = t704 * t928 + t953;
t496 = Icges(6,5) * t629 + Icges(6,6) * t628 - Icges(6,3) * t959;
t1004 = Icges(6,4) * t629;
t499 = Icges(6,2) * t628 - Icges(6,6) * t959 + t1004;
t621 = Icges(6,4) * t628;
t502 = Icges(6,1) * t629 - Icges(6,5) * t959 + t621;
t630 = t704 * t714 + t705 * t939;
t961 = t704 * t711;
t631 = t705 * t928 - t961;
t308 = t496 * t950 + t630 * t499 - t631 * t502;
t1153 = t308 * t704;
t1152 = t308 * t705;
t1013 = rSges(6,1) * t714;
t795 = -rSges(6,2) * t711 + t1013;
t647 = rSges(6,3) * t712 + t715 * t795;
t1102 = t631 * rSges(6,1) - t630 * rSges(6,2);
t507 = -rSges(6,3) * t950 + t1102;
t932 = t507 * t712;
t1151 = t647 * t950 + t932;
t623 = Icges(6,4) * t631;
t501 = Icges(6,2) * t630 + Icges(6,6) * t950 - t623;
t622 = Icges(6,4) * t630;
t503 = Icges(6,1) * t631 - Icges(6,5) * t950 - t622;
t1142 = t501 * t711 + t503 * t714;
t498 = -Icges(6,5) * t631 + Icges(6,6) * t630 + Icges(6,3) * t950;
t972 = t498 * t712;
t343 = t1142 * t715 - t972;
t596 = Icges(7,4) * t610;
t485 = Icges(7,2) * t609 + Icges(7,6) * t950 - t596;
t595 = Icges(7,4) * t609;
t487 = Icges(7,1) * t610 - Icges(7,5) * t950 - t595;
t1143 = t485 * t706 + t487 * t707;
t482 = -Icges(7,5) * t610 + Icges(7,6) * t609 + Icges(7,3) * t950;
t975 = t482 * t712;
t319 = t1143 * t715 - t975;
t717 = -pkin(9) - pkin(8);
t1021 = -pkin(8) - t717;
t701 = pkin(5) * t714 + pkin(4);
t1022 = pkin(4) - t701;
t606 = t1021 * t712 - t1022 * t715;
t571 = t606 * t950;
t794 = rSges(7,1) * t707 - rSges(7,2) * t706;
t620 = rSges(7,3) * t712 + t715 * t794;
t579 = t620 * t950;
t491 = -t610 * rSges(7,1) + t609 * rSges(7,2) + rSges(7,3) * t950;
t692 = pkin(8) * t950;
t952 = t705 * t712;
t802 = pkin(4) * t952 - t692;
t850 = pkin(5) * t961;
t920 = t715 * t717;
t1122 = t712 * t701 + t920;
t865 = t1122 * t705;
t510 = t802 + t850 - t865;
t888 = t491 + t510;
t333 = -t712 * t888 + t571 + t579;
t1145 = -t491 * t712 + t579;
t971 = t510 * t712;
t334 = t1145 + t571 - t971;
t906 = t333 - t334;
t1150 = m(7) * t906;
t774 = -t485 * t609 - t487 * t610;
t898 = t607 * t483 + t608 * t486;
t1149 = t898 + (-t480 * t704 - t482 * t705) * t715 + t774;
t771 = -t501 * t630 - t503 * t631;
t896 = t628 * t499 + t629 * t502;
t1148 = t896 + (-t496 * t704 - t498 * t705) * t715 + t771;
t855 = 0.2e1 * m(7);
t810 = t855 / 0.4e1;
t827 = -t952 / 0.4e1;
t1141 = 0.2e1 * m(6);
t1138 = t1141 / 0.4e1;
t999 = Icges(7,4) * t707;
t787 = -Icges(7,2) * t706 + t999;
t616 = Icges(7,6) * t712 + t715 * t787;
t1000 = Icges(7,4) * t706;
t790 = Icges(7,1) * t707 - t1000;
t618 = Icges(7,5) * t712 + t715 * t790;
t783 = Icges(7,5) * t707 - Icges(7,6) * t706;
t614 = Icges(7,3) * t712 + t715 * t783;
t923 = t715 * t614;
t385 = t609 * t616 - t610 * t618 + t705 * t923;
t1135 = t385 * t712;
t286 = t482 * t950 - t774;
t780 = t286 * t705 - t1155;
t154 = t715 * t780 + t1135;
t283 = -t480 * t959 + t898;
t914 = -t283 + t1149;
t1137 = t154 - t1135 + (t705 * t914 + t1155) * t715;
t522 = rSges(7,1) * t609 + rSges(7,2) * t610;
t625 = t630 * pkin(5);
t1125 = -t522 - t625;
t521 = t607 * rSges(7,1) - t608 * rSges(7,2);
t624 = t628 * pkin(5);
t882 = -t521 - t624;
t1144 = t1125 * t704 + t705 * t882;
t307 = -t498 * t959 + t628 * t501 - t503 * t629;
t284 = -t482 * t959 + t607 * t485 - t487 * t608;
t1069 = m(6) / 0.4e1;
t1067 = m(7) / 0.4e1;
t1053 = t704 / 0.2e1;
t1052 = t704 / 0.4e1;
t1051 = t705 / 0.2e1;
t1050 = t705 / 0.4e1;
t1049 = t712 / 0.2e1;
t1048 = t715 / 0.2e1;
t179 = t286 * t704 + t1154;
t60 = t704 * t914 - t1154;
t1136 = t179 + t60;
t1002 = Icges(6,4) * t714;
t788 = -Icges(6,2) * t711 + t1002;
t635 = Icges(6,6) * t712 + t715 * t788;
t1003 = Icges(6,4) * t711;
t791 = Icges(6,1) * t714 - t1003;
t637 = Icges(6,5) * t712 + t715 * t791;
t784 = Icges(6,5) * t714 - Icges(6,6) * t711;
t633 = Icges(6,3) * t712 + t715 * t784;
t922 = t715 * t633;
t410 = t630 * t635 - t631 * t637 + t705 * t922;
t1131 = t712 * t410;
t547 = t616 * t704;
t549 = t618 * t704;
t775 = -t483 * t706 + t486 * t707;
t752 = t614 * t704 - t775;
t243 = (-t547 * t706 + t549 * t707 + t480) * t715 + t752 * t712;
t613 = Icges(7,3) * t715 - t712 * t783;
t944 = t707 * t618;
t948 = t706 * t616;
t765 = t944 - t948;
t748 = t613 - t765;
t967 = t614 * t712;
t1095 = t715 * t748 - t967;
t615 = Icges(7,6) * t715 - t712 * t787;
t617 = Icges(7,5) * t715 - t712 * t790;
t289 = -t1095 * t704 + t607 * t615 + t608 * t617;
t1130 = t243 + t289;
t548 = t616 * t705;
t550 = t618 * t705;
t751 = -t614 * t705 + t1143;
t244 = (t548 * t706 - t550 * t707 + t482) * t715 + t751 * t712;
t290 = t1095 * t705 + t609 * t615 - t610 * t617;
t1129 = t244 + t290;
t976 = t480 * t712;
t318 = t715 * t775 + t976;
t383 = t607 * t616 + t608 * t618 - t704 * t923;
t1128 = t318 + t383;
t1127 = -t319 + t385;
t1006 = Icges(5,4) * t712;
t789 = Icges(5,2) * t715 + t1006;
t589 = Icges(5,6) * t705 + t704 * t789;
t679 = Icges(5,4) * t959;
t960 = t704 * t712;
t591 = Icges(5,1) * t960 + Icges(5,5) * t705 + t679;
t767 = -t589 * t715 - t591 * t712;
t1126 = t705 * t767;
t869 = t618 + (-Icges(7,2) * t707 - t1000) * t715;
t686 = t715 * pkin(4) + t712 * pkin(8);
t871 = t606 + t620;
t462 = (t686 + t871) * t705;
t978 = t462 * t704;
t449 = -0.2e1 * t978;
t544 = (t647 + t686) * t705;
t970 = t544 * t704;
t528 = -0.2e1 * t970;
t903 = (t449 + 0.2e1 * t978) * t1067 + (t528 + 0.2e1 * t970) * t1069;
t1070 = 0.2e1 * t705;
t663 = t704 * t686;
t807 = t871 * t704;
t460 = t663 + t807;
t542 = t647 * t704 + t663;
t904 = (t1070 * t460 + t449) * t1067 + (t1070 * t542 + t528) * t1069;
t1124 = t903 - t904;
t489 = t608 * rSges(7,1) + t607 * rSges(7,2) - rSges(7,3) * t959;
t1101 = t1022 * t712;
t688 = pkin(8) * t959;
t849 = pkin(5) * t953;
t509 = t849 + t688 + (t920 - t1101) * t704;
t890 = t489 + t509;
t332 = t712 * t890 + t715 * t807;
t988 = t332 * t704;
t224 = t333 * t1070 + 0.2e1 * t988;
t506 = rSges(6,1) * t629 + rSges(6,2) * t628 - rSges(6,3) * t959;
t933 = t712 * t506;
t444 = t647 * t959 + t933;
t979 = t444 * t704;
t339 = t1070 * t1151 + 0.2e1 * t979;
t909 = t1067 * t224 + t1069 * t339;
t1071 = -0.2e1 * t705;
t918 = (t1071 * t334 + t224 - 0.2e1 * t988) * t1067 + (t1071 * t1151 + t339 - 0.2e1 * t979) * t1069;
t1123 = t909 - t918;
t1011 = rSges(6,3) * t715;
t1024 = t712 * pkin(4);
t1066 = pkin(2) + pkin(7);
t814 = -sin(qJ(1)) * pkin(1) + t705 * qJ(3);
t739 = -t1066 * t704 + t814;
t426 = (-t1011 + t1024) * t705 - t692 + t739 + t1102;
t806 = pkin(5) * t711 + t1066;
t1099 = -t704 * t806 - t491 + t814 + t865;
t1023 = cos(qJ(1)) * pkin(1);
t397 = t1023 + t806 * t705 + (qJ(3) + t1122) * t704 + t489;
t427 = t1023 - t688 + t1066 * t705 + (qJ(3) + t1024) * t704 + t506;
t536 = rSges(6,1) * t628 - rSges(6,2) * t629;
t537 = rSges(6,1) * t630 + rSges(6,2) * t631;
t655 = (-rSges(7,1) * t706 - rSges(7,2) * t707) * t715;
t938 = t711 * t715;
t801 = pkin(5) * t938 - t655;
t553 = t801 * t704;
t554 = t801 * t705;
t661 = (-rSges(6,1) * t711 - rSges(6,2) * t714) * t715;
t1121 = (-t1099 * t553 + t1125 * t460 + t397 * t554 + t462 * t882) * t810 + (-t536 * t544 - t537 * t542 + (t426 * t704 - t427 * t705) * t661) * t1138;
t649 = pkin(4) * t959 + pkin(8) * t960;
t921 = t715 * t701;
t664 = t704 * t921;
t927 = t712 * t717;
t538 = -t704 * t927 - t649 + t664;
t605 = t1021 * t715 + t1101;
t946 = t706 * t715;
t847 = rSges(7,2) * t946;
t942 = t707 * t715;
t864 = t704 * rSges(7,1) * t942 + rSges(7,3) * t960;
t551 = -t704 * t847 + t864;
t619 = rSges(7,3) * t715 - t712 * t794;
t836 = t715 * t489 + t712 * t551 + t619 * t959;
t245 = t509 * t715 + t538 * t712 + (t605 * t715 - t712 * t871) * t704 + t836;
t578 = t619 * t950;
t1107 = t620 * t705;
t863 = -pkin(4) * t950 - pkin(8) * t952;
t539 = (-t921 + t927) * t705 - t863;
t879 = t539 - t1107;
t246 = t578 + (t605 * t705 - t888) * t715 + (-t705 * t871 - t879) * t712;
t565 = -rSges(6,2) * t704 * t938 + rSges(6,3) * t960 + t959 * t1013;
t646 = -t712 * t795 + t1011;
t762 = t646 * t715 - t647 * t712;
t365 = t506 * t715 + t565 * t712 + t704 * t762;
t1106 = t647 * t705;
t366 = t1106 * t712 + t507 * t715 + t705 * t762;
t477 = t664 + (-t847 - t927) * t704 + t864;
t478 = ((rSges(7,3) - t717) * t712 + (t701 + t794) * t715) * t705;
t513 = -t863 + t1106;
t878 = -t565 - t649;
t1120 = -(t1151 * t513 + t365 * t427 + t366 * t426 - t444 * t878) * t1141 / 0.4e1 - (t1099 * t246 + t245 * t397 + t332 * t477 + t333 * t478) * t855 / 0.4e1;
t1119 = -t462 * t1150 / 0.2e1;
t652 = (-Icges(7,5) * t706 - Icges(7,6) * t707) * t715;
t931 = t712 * t652;
t1118 = t931 / 0.2e1 - t869 * t946 / 0.2e1;
t813 = 0.2e1 * t1069;
t702 = t704 ^ 2;
t1112 = t702 / 0.2e1;
t1054 = -t704 / 0.2e1;
t1111 = t284 * t704;
t1110 = t284 * t705;
t1109 = t307 * t704;
t1108 = t307 * t705;
t590 = -Icges(5,6) * t704 + t705 * t789;
t1005 = Icges(5,4) * t715;
t792 = Icges(5,1) * t712 + t1005;
t592 = -Icges(5,5) * t704 + t705 * t792;
t1104 = (t590 * t715 + t592 * t712) * t705;
t968 = t1106 * t715;
t969 = t565 * t715;
t313 = (t933 - t969) * t705 + (-t932 + t968) * t704;
t494 = t705 * t522;
t350 = t705 * t625 + t882 * t704 + t494;
t276 = (-t704 * t888 - t705 * t890) * t715;
t769 = -t521 * t705 - t522 * t704;
t390 = t769 * t715;
t447 = t712 * t521 + t655 * t959;
t586 = t655 * t950;
t448 = -t522 * t712 + t586;
t799 = t276 * t390 + t332 * t447 + t333 * t448;
t515 = Icges(7,5) * t607 - Icges(7,6) * t608;
t892 = -Icges(7,2) * t608 + t486 + t594;
t894 = -Icges(7,1) * t607 + t1001 + t483;
t217 = -t515 * t959 + t607 * t892 - t608 * t894;
t516 = Icges(7,5) * t609 + Icges(7,6) * t610;
t891 = Icges(7,2) * t610 - t487 + t595;
t893 = -Icges(7,1) * t609 + t485 - t596;
t218 = -t516 * t959 + t607 * t891 - t608 * t893;
t118 = t217 * t705 + t218 * t704;
t219 = t515 * t950 + t609 * t892 + t610 * t894;
t220 = t516 * t950 + t609 * t891 + t610 * t893;
t119 = t219 * t705 + t220 * t704;
t919 = t1051 * t118 + t1053 * t119;
t945 = t707 * t617;
t949 = t706 * t615;
t340 = (t614 + t945 - t949) * t715 + t748 * t712;
t436 = t715 * t765 + t967;
t1098 = t436 * t1048 + t1049 * t340;
t777 = -t318 * t704 - t319 * t705;
t1097 = ((-t243 * t704 + t244 * t705 + t436) * t715 + (t340 - t777) * t712) * t1049 + (t436 * t712 + t715 * t777) * t1048;
t768 = -t536 * t705 - t537 * t704;
t902 = t768 * t1138 + t1144 * t810;
t632 = Icges(6,3) * t715 - t712 * t784;
t925 = t714 * t637;
t940 = t711 * t635;
t764 = t925 - t940;
t747 = t632 - t764;
t966 = t633 * t712;
t1096 = t715 * t747 - t966;
t749 = -t633 * t705 + t1142;
t1094 = t715 * t749 - t972;
t772 = -t499 * t711 + t502 * t714;
t750 = t633 * t704 - t772;
t973 = t496 * t712;
t1093 = t715 * t750 - t973;
t1092 = t715 * t751 - t975;
t1091 = t715 * t752 - t976;
t797 = rSges(5,1) * t712 + rSges(5,2) * t715;
t737 = -t704 * rSges(5,3) + t705 * t797;
t524 = t737 + t739;
t525 = t1023 + (rSges(5,3) + t1066) * t705 + (qJ(3) + t797) * t704;
t677 = rSges(5,1) * t715 - rSges(5,2) * t712;
t644 = t677 * t704;
t645 = rSges(5,1) * t950 - rSges(5,2) * t952;
t674 = Icges(5,1) * t715 - t1006;
t1090 = -m(5) * (t524 * t645 + t525 * t644) - (-t944 / 0.2e1 + t948 / 0.2e1 + t613 / 0.2e1 - t925 / 0.2e1 + t940 / 0.2e1 + t632 / 0.2e1 - t674 / 0.2e1 + t789 / 0.2e1) * t712;
t763 = -t644 * t705 + t645 * t704;
t837 = (-t477 * t705 + t478 * t704) * t810 + (t513 * t704 + t705 * t878) * t1138 + m(5) * t763 / 0.2e1;
t659 = (-Icges(6,2) * t714 - t1003) * t715;
t660 = (-Icges(6,1) * t711 - t1002) * t715;
t1089 = -t711 * (t637 / 0.2e1 + t659 / 0.2e1) + t714 * (t660 / 0.2e1 - t635 / 0.2e1);
t672 = -Icges(5,2) * t712 + t1005;
t641 = t672 * t705;
t643 = t674 * t705;
t1088 = t712 * (t590 - t643) - t715 * (t592 + t641);
t640 = -Icges(5,2) * t960 + t679;
t642 = t674 * t704;
t1087 = t712 * (t589 - t642) - (t591 + t640) * t715;
t530 = Icges(6,5) * t628 - Icges(6,6) * t629;
t885 = -Icges(6,2) * t629 + t502 + t621;
t887 = -Icges(6,1) * t628 + t1004 + t499;
t273 = t530 * t712 + (-t711 * t885 - t714 * t887) * t715;
t531 = Icges(6,5) * t630 + Icges(6,6) * t631;
t884 = Icges(6,2) * t631 - t503 + t622;
t886 = -Icges(6,1) * t630 + t501 - t623;
t274 = t531 * t712 + (-t711 * t884 - t714 * t886) * t715;
t658 = (-Icges(6,5) * t711 - Icges(6,6) * t714) * t715;
t866 = t637 + t659;
t867 = t635 - t660;
t346 = t628 * t866 - t629 * t867 - t658 * t959;
t347 = t630 * t866 + t631 * t867 + t658 * t950;
t1086 = t1121 + (t274 + t347) * t1052 + (t273 + t346) * t1050;
t654 = (-Icges(7,1) * t706 - t999) * t715;
t870 = t616 - t654;
t311 = t607 * t869 - t608 * t870 - t652 * t959;
t312 = t609 * t869 + t610 * t870 + t652 * t950;
t258 = t516 * t712 + (-t706 * t891 - t707 * t893) * t715;
t991 = t258 * t704;
t257 = t515 * t712 + (-t706 * t892 - t707 * t894) * t715;
t992 = t257 * t705;
t805 = t992 / 0.4e1 + t991 / 0.4e1 + t312 * t1052 + t311 * t1050;
t381 = t383 * t712;
t781 = -t283 * t704 + t1110;
t153 = t715 * t781 + t381;
t917 = t286 + t1149;
t49 = t381 + (-t704 * t917 + t1110) * t715;
t1085 = t49 * t1052 - t704 * t153 / 0.4e1 + t1137 * t1050;
t306 = -t496 * t959 + t896;
t190 = t306 * t705 + t1109;
t309 = t498 * t950 - t771;
t191 = t309 * t704 + t1152;
t408 = t628 * t635 + t629 * t637 - t704 * t922;
t406 = t408 * t712;
t913 = t309 + t1148;
t65 = t406 + (-t704 * t913 + t1108) * t715;
t910 = -t306 + t1148;
t66 = -t1131 + (t705 * t910 + t1153) * t715;
t76 = t705 * t913 + t1109;
t77 = t704 * t910 - t1152;
t822 = t950 / 0.4e1;
t824 = -t950 / 0.4e1;
t831 = -t959 / 0.4e1;
t778 = t309 * t705 - t1153;
t163 = t715 * t778 + t1131;
t956 = t705 * t163;
t779 = -t306 * t704 + t1108;
t162 = t715 * t779 + t406;
t963 = t704 * t162;
t1084 = -t963 / 0.4e1 + t956 / 0.4e1 + t190 * t824 + t65 * t1052 + t66 * t1050 + t76 * t822 - t1119 + (t191 + t77) * t831;
t561 = t635 * t704;
t563 = t637 * t704;
t266 = (-t561 * t711 + t563 * t714 + t496) * t715 + t750 * t712;
t562 = t635 * t705;
t564 = t637 * t705;
t267 = (t562 * t711 - t564 * t714 + t498) * t715 + t749 * t712;
t634 = Icges(6,6) * t715 - t712 * t788;
t636 = Icges(6,5) * t715 - t712 * t791;
t324 = -t1096 * t704 + t628 * t634 + t629 * t636;
t325 = t1096 * t705 + t630 * t634 - t631 * t636;
t342 = t715 * t772 + t973;
t926 = t714 * t636;
t941 = t711 * t634;
t364 = (t633 + t926 - t941) * t715 + t747 * t712;
t452 = t715 * t764 + t966;
t833 = t960 / 0.4e1;
t1083 = t452 * t1048 + t364 * t1049 - t1120 + (t342 + t408) * t833 + (t266 + t324) * t831 + (-t343 + t410) * t827 + (t267 + t325) * t822;
t703 = t705 ^ 2;
t422 = t489 * t712 + t620 * t959;
t1076 = 0.2e1 * t422;
t1075 = 0.2e1 * t1145;
t1074 = 0.2e1 * t448;
t1073 = -0.2e1 * t462;
t1072 = -0.2e1 * t522;
t1068 = m(7) / 0.2e1;
t880 = -t538 - t551;
t470 = t491 * t960;
t471 = t489 * t952;
t895 = t470 + t471;
t186 = (t509 * t705 + t510 * t704) * t712 + (-t704 * t879 + t705 * t880) * t715 + t895;
t291 = (t1107 * t704 - t551 * t705) * t715 + t895;
t351 = -t620 * t960 + t836;
t352 = -t491 * t715 + t578;
t371 = (-t489 * t705 - t491 * t704) * t715;
t1064 = (t1145 * t246 + t186 * t371 + t245 * t422 + t276 * t291 + t332 * t351 + t333 * t352) * t855;
t1062 = -0.2e1 * t422 * t1150;
t405 = -t704 * t521 + t494;
t603 = t705 * t802;
t648 = pkin(4) * t960 - t688;
t264 = -t603 + t888 * t705 + (-t648 - t890) * t704;
t202 = 0.2e1 * t390 * t264;
t360 = t447 * t1073;
t361 = t460 * t1074;
t838 = t202 + t360 + t361;
t853 = 0.2e1 * t655;
t1060 = m(7) * (0.2e1 * t405 * t276 + (-t332 * t705 + t333 * t704) * t853 + t838);
t1058 = m(7) * (-t1075 * t553 + t1076 * t554 + 0.2e1 * t350 * t371 + t838);
t1057 = (t1099 * t352 + t1145 * t478 + t351 * t397 + t422 * t477) * t855;
t1056 = -t162 / 0.2e1;
t1055 = t274 / 0.2e1;
t322 = 0.2e1 * t447 * t397;
t323 = t1099 * t1074;
t907 = t322 + t323;
t1037 = m(7) * (t1072 * t333 + 0.2e1 * t332 * t521 + t907);
t1035 = m(7) * (t1075 * t1125 - t1076 * t882 + t907);
t1033 = m(7) * (t460 * t1072 + t521 * t1073 + (t1099 * t704 - t397 * t705) * t853);
t981 = t422 * t704;
t317 = t1070 * t1145 + 0.2e1 * t981;
t1032 = m(7) * (t1071 * t1145 + t317 - 0.2e1 * t981);
t1030 = m(7) * t317;
t1025 = t769 * t855;
t1019 = m(6) * qJD(4);
t1018 = m(6) * qJD(5);
t1016 = m(7) * qJD(4);
t1015 = m(7) * qJD(5);
t985 = t351 * t705;
t984 = t352 * t704;
t983 = t365 * t705;
t982 = t366 * t704;
t930 = t712 * t658;
t851 = -0.2e1 * t950;
t852 = -0.2e1 * t959;
t391 = t521 * t851 + t522 * t852;
t872 = t605 + t619;
t862 = t702 + t703;
t861 = qJD(1) * t715;
t860 = qJD(5) * t715;
t368 = (t931 + (-t706 * t869 - t707 * t870) * t715) * t712;
t823 = t950 / 0.2e1;
t832 = -t959 / 0.2e1;
t89 = t311 * t712 + (-t217 * t704 + t218 * t705) * t715;
t90 = t312 * t712 + (-t219 * t704 + t220 * t705) * t715;
t848 = t90 * t823 + t89 * t832 + (t368 + (-t257 * t704 + t258 * t705) * t715) * t1049;
t842 = t1015 / 0.2e1;
t841 = t1014 / 0.2e1;
t840 = t1056 + t65 / 0.2e1;
t839 = -t66 / 0.2e1 - t163 / 0.2e1;
t785 = Icges(5,5) * t712 + Icges(5,6) * t715;
t416 = t705 * (Icges(5,3) * t705 + t704 * t785) + t589 * t959 + t591 * t960;
t588 = -Icges(5,3) * t704 + t705 * t785;
t417 = -t705 * t588 - t590 * t959 - t592 * t960;
t834 = t960 / 0.2e1;
t830 = t959 / 0.2e1;
t829 = t959 / 0.4e1;
t828 = -t952 / 0.2e1;
t825 = -t950 / 0.2e1;
t819 = -t942 / 0.2e1;
t818 = t942 / 0.2e1;
t815 = t1112 + t703 / 0.2e1;
t811 = t855 / 0.2e1;
t804 = -t1107 * t852 + t551 * t851 + 0.2e1 * t470 + 0.2e1 * t471;
t803 = 0.2e1 * t862;
t798 = t1137 * t832 + t153 * t825 + t49 * t823;
t793 = t616 * t819 + t654 * t818 + t1118;
t786 = Icges(5,5) * t715 - Icges(5,6) * t712;
t82 = ((-t554 / 0.2e1 + t245 / 0.2e1) * t705 + (-t553 / 0.2e1 - t246 / 0.2e1) * t704) * m(7) + (t983 / 0.2e1 - t982 / 0.2e1 + t815 * t661) * m(6);
t718 = (t804 + (t509 * t712 - t538 * t715) * t1070 + 0.2e1 * (-t539 * t715 + t971) * t704) * t1067;
t727 = t350 * t855;
t91 = t718 - t727 / 0.4e1 + ((t933 / 0.2e1 - t969 / 0.2e1 - t537 / 0.2e1) * t705 + (-t932 / 0.2e1 + t968 / 0.2e1 + t536 / 0.2e1) * t704) * m(6);
t782 = t91 * qJD(2) - t82 * qJD(3);
t776 = -t342 * t704 - t343 * t705;
t421 = -t536 * t704 + t537 * t705;
t761 = t803 / 0.4e1;
t208 = (t985 / 0.2e1 - t984 / 0.2e1 + t815 * t655) * m(7);
t743 = t804 * t1067;
t744 = t405 * t855;
t232 = t743 - t744 / 0.4e1;
t760 = -t232 * qJD(2) + t208 * qJD(3);
t211 = -t1091 * t704 + t547 * t607 + t549 * t608;
t212 = -t1092 * t704 - t548 * t607 - t550 * t608;
t39 = (-t211 * t704 + t212 * t705 + t383) * t715 + (t289 - t781) * t712;
t213 = t1091 * t705 + t547 * t609 - t549 * t610;
t214 = t1092 * t705 - t548 * t609 + t550 * t610;
t40 = (-t213 * t704 + t214 * t705 + t385) * t715 + (t290 - t780) * t712;
t759 = t153 * t834 + t154 * t828 + t39 * t832 + t40 * t823 + t1097;
t758 = t1062 / 0.4e1 + t798;
t757 = -m(5) * (t524 * t705 + t525 * t704) - m(4) * ((rSges(4,3) * t705 + t814) * t705 + (t1023 + (rSges(4,3) + qJ(3)) * t704) * t704);
t247 = -t530 * t959 + t628 * t885 - t629 * t887;
t248 = -t531 * t959 + t628 * t884 - t629 * t886;
t134 = t247 * t705 + t248 * t704;
t249 = t530 * t950 + t630 * t885 + t631 * t887;
t250 = t531 * t950 + t630 * t884 + t631 * t886;
t135 = t249 * t705 + t250 * t704;
t746 = t1051 * t134 + t1053 * t135;
t742 = t1064 / 0.4e1 + t759;
t738 = t1145 * t448 + t371 * t390 + t422 * t447;
t178 = t283 * t705 + t1111;
t59 = t705 * t917 + t1111;
t736 = t1136 * t831 + t178 * t824 + t59 * t822 + t1085;
t735 = t616 * t818 + t654 * t819 - t1118;
t419 = -t704 * t588 + t1104;
t734 = t588 * t1112 + t419 * t1053 + t179 / 0.2e1 + t60 / 0.2e1 + t77 / 0.2e1 + t191 / 0.2e1 + (t417 + (t588 - t767) * t705 + t1126) * t1051;
t733 = -t416 * t705 / 0.2e1 + t417 * t1054 + (t417 - t1126) * t1053 + (t419 - t1104 + (t588 + t767) * t704 + t416) * t1051 + t59 / 0.2e1 - t178 / 0.2e1 - t190 / 0.2e1 + t76 / 0.2e1;
t731 = t49 * t825 + t368 + (t257 + t311) * t832 + t1137 * t830 + (t153 + t258 + t312) * t823;
t110 = t211 * t705 + t212 * t704;
t111 = t213 * t705 + t214 * t704;
t730 = t39 * t1051 + t40 * t1053 + t110 * t832 + t111 * t823 + (t243 * t705 + t244 * t704) * t1049 + t178 * t834 + t179 * t828 + (t318 * t705 - t319 * t704) * t1048 - t919;
t729 = t1127 * t827 + t1128 * t833 + t1129 * t822 + t1130 * t831 + t1098;
t728 = t39 * t830 + t40 * t825 + t118 * t832 + t119 * t823 + (t991 + t992) * t1049 - t153 * t960 / 0.2e1 + t154 * t952 / 0.2e1 + t90 * t1053 + t89 * t1051 - t1097;
t722 = -t926 / 0.2e1 + t941 / 0.2e1 - t945 / 0.2e1 + t949 / 0.2e1 + t792 / 0.2e1 + t672 / 0.2e1 - t633 / 0.2e1 - t614 / 0.2e1;
t721 = t1136 * t829 + t178 * t822 + t59 * t824 - t1085 + t729 + t805;
t720 = -t1098 + t736 + t805 - t1128 * t960 / 0.4e1 + t1130 * t829 + t1127 * t952 / 0.4e1 + t1129 * t824;
t719 = t729 + t736 - t805;
t709 = t715 ^ 2;
t685 = t715 * pkin(8) - t1024;
t662 = t704 * t685;
t639 = t786 * t705;
t638 = t704 * t786;
t604 = t705 * t863;
t543 = (-t646 - t685) * t705;
t541 = t646 * t704 + t662;
t466 = -t537 * t712 + t661 * t950;
t465 = t536 * t712 + t661 * t959;
t461 = (-t685 - t872) * t705;
t459 = t704 * t872 + t662;
t411 = t768 * t715;
t403 = (t930 + (-t711 * t866 - t714 * t867) * t715) * t712;
t400 = t1025 / 0.4e1;
t393 = t1125 * t712 - t709 * t849 + t586;
t392 = t624 * t712 - t709 * t850 + t447;
t388 = (-t506 * t705 + t507 * t704) * t715;
t387 = t391 * t841;
t382 = -t1106 * t705 + t704 * t878 + t604;
t354 = -t507 * t705 - t603 + (-t506 - t648) * t704;
t345 = -0.2e1 * t447 * t705 + 0.2e1 * t448 * t704;
t341 = t1144 * t715;
t338 = t345 * t841;
t321 = 0.4e1 * t426 * t705 + 0.4e1 * t427 * t704;
t316 = t1030 / 0.4e1;
t287 = t604 + t879 * t705 + (-t649 + t880) * t704;
t278 = 0.4e1 * t1099 * t705 + 0.4e1 * t397 * t704;
t268 = -0.4e1 * t426 * t537 + 0.4e1 * t427 * t536;
t265 = 0.4e1 * t426 * t513 - 0.4e1 * t427 * t878;
t242 = -0.4e1 * t1099 * t522 + 0.4e1 * t397 * t521;
t241 = t1094 * t705 - t562 * t630 + t564 * t631;
t240 = t1093 * t705 + t561 * t630 - t563 * t631;
t239 = -t1094 * t704 - t562 * t628 - t564 * t629;
t238 = -t1093 * t704 + t561 * t628 + t563 * t629;
t233 = t744 / 0.4e1 + t743;
t228 = 0.4e1 * t1099 * t478 + 0.4e1 * t397 * t477;
t221 = 0.4e1 * t1099 * t1125 - 0.4e1 * t397 * t882;
t216 = 0.4e1 * t421 * t354 + 0.4e1 * (t542 * t704 + t544 * t705) * t661;
t207 = (t984 - t985) * t810 + m(7) * t655 * t761;
t204 = t1032 / 0.4e1;
t187 = t452 * t712 + t715 * t776;
t180 = t1067 * t242 + t793;
t164 = t1033 / 0.4e1;
t156 = 0.4e1 * t405 * t264 + 0.4e1 * (t460 * t704 + t462 * t705) * t655;
t155 = 0.4e1 * t738;
t137 = t1067 * t278 + t1069 * t321 - t757;
t131 = 0.4e1 * t264 * t350 - 0.4e1 * t460 * t553 - 0.4e1 * t462 * t554;
t130 = t240 * t705 + t241 * t704;
t129 = t238 * t705 + t239 * t704;
t128 = t316 + t204 - t1025 / 0.4e1;
t127 = t400 + t316 - t1032 / 0.4e1;
t126 = t400 + t204 - t1030 / 0.4e1;
t124 = t1035 / 0.4e1;
t108 = 0.4e1 * t1151 * t366 + 0.4e1 * t313 * t388 + 0.4e1 * t365 * t444;
t105 = t903 + t904 - t837;
t104 = t837 + t1124;
t103 = t837 - t1124;
t102 = t347 * t712 + (-t249 * t704 + t250 * t705) * t715;
t101 = t346 * t712 + (-t247 * t704 + t248 * t705) * t715;
t99 = t1037 / 0.4e1;
t96 = t1057 / 0.4e1;
t95 = 0.4e1 * t1145 * t352 + 0.4e1 * t291 * t371 + 0.4e1 * t351 * t422;
t94 = 0.4e1 * t799;
t93 = t221 * t1067 + t930 / 0.2e1 + t268 * t1069 + t1089 * t715 + t793;
t92 = t727 / 0.4e1 + t718 + (t421 + t313) * t813;
t81 = (t982 - t983) * t813 + m(6) * t661 * t761 + ((-t245 - t554) * t705 + (t246 - t553) * t704) * t810;
t80 = t228 * t1067 + t265 * t1069 - t722 * t715 - t1090;
t79 = (-t266 * t704 + t267 * t705 + t452) * t715 + (t364 - t776) * t712;
t78 = 0.4e1 * t332 * t906;
t72 = t1058 / 0.4e1;
t62 = (-t240 * t704 + t241 * t705 + t410) * t715 + (t325 - t778) * t712;
t61 = (-t238 * t704 + t239 * t705 + t408) * t715 + (t324 - t779) * t712;
t53 = t1060 / 0.4e1;
t32 = 0.4e1 * t186 * t276 + 0.4e1 * t245 * t332 + 0.4e1 * t246 * t333;
t29 = t909 + t918 - t902;
t28 = t902 - t1123;
t27 = t902 + t1123;
t24 = t1067 * t156 + t919;
t21 = t1067 * t155 + t848;
t20 = t21 * qJD(6);
t19 = t1067 * t94 + t848;
t18 = t1067 * t131 + t1069 * t216 + t746 + t919;
t16 = t1067 * t95 + t759;
t15 = t124 - t1037 / 0.4e1 + t758;
t14 = t99 - t1035 / 0.4e1 + t758;
t13 = t99 + t124 - t1062 / 0.4e1 + t731;
t12 = t704 * t733 + t705 * t734;
t11 = -t78 * t1067 + (t704 * t839 + t705 * t840) * t715 + t798;
t10 = t72 - t1060 / 0.4e1 + t742;
t9 = t53 - t1058 / 0.4e1 + t742;
t8 = t32 * t1067 + t108 * t1069 + (t187 / 0.2e1 + t61 * t1054 + t62 * t1051) * t715 + (t79 / 0.2e1 + t963 / 0.2e1 - t956 / 0.2e1) * t712 + t759;
t7 = t96 + t721 + t164;
t6 = t96 + t719 - t1033 / 0.4e1;
t5 = t720 - t1057 / 0.4e1 + t164;
t4 = t728 + t72 + t53 - t1064 / 0.4e1;
t3 = t1083 + (t162 / 0.4e1 - t65 / 0.4e1 + (t191 / 0.4e1 + t77 / 0.4e1) * t715) * t704 + (-t163 / 0.4e1 - t66 / 0.4e1 + (t190 / 0.4e1 - t76 / 0.4e1) * t715) * t705 + t721 + t1086 + t1119;
t2 = t1084 + (-t364 / 0.2e1 + (t410 / 0.4e1 - t343 / 0.4e1) * t705 + (-t408 / 0.4e1 - t342 / 0.4e1) * t704) * t712 + (-t452 / 0.2e1 + (-t267 / 0.4e1 - t325 / 0.4e1) * t705 + (t266 / 0.4e1 + t324 / 0.4e1) * t704) * t715 + t720 + t1086 + t1120;
t1 = t1084 + t1083 + (-t346 / 0.4e1 - t273 / 0.4e1) * t705 + (-t347 / 0.4e1 - t274 / 0.4e1) * t704 + t719 - t1121;
t17 = [t137 * qJD(3) + t80 * qJD(4) + t93 * qJD(5) + t180 * qJD(6), 0, qJD(1) * t137 + qJD(4) * t103 + qJD(5) * t27 + qJD(6) * t127, t80 * qJD(1) + t103 * qJD(3) + t3 * qJD(5) + t7 * qJD(6) + (t1099 * t459 + t397 * t461 + t460 * t478 - t462 * t477) * t1016 + (t426 * t541 + t427 * t543 + t513 * t542 + t878 * t544) * t1019 + ((t324 / 0.2e1 + t289 / 0.2e1 + t243 / 0.2e1 + t266 / 0.2e1 - t785 * t1051 + (-t589 / 0.2e1 + t642 / 0.2e1) * t715 + (-t591 / 0.2e1 - t640 / 0.2e1) * t712 - t734) * t705 + (t267 / 0.2e1 + t325 / 0.2e1 + t290 / 0.2e1 + t244 / 0.2e1 - t785 * t1053 + (t590 / 0.2e1 - t643 / 0.2e1) * t715 + (t592 / 0.2e1 + t641 / 0.2e1) * t712 - t733) * t704 + (t763 * t677 - (t524 * t704 - t525 * t705) * t797) * m(5)) * qJD(4), t93 * qJD(1) + t27 * qJD(3) + t3 * qJD(4) + (t403 + t731) * qJD(5) + t13 * qJD(6) + (t78 / 0.4e1 - t332 * t882 + t333 * t1125 + t392 * t397 + t393 * t1099) * t1015 + (-t1151 * t537 + t426 * t466 + t427 * t465 + t444 * t536) * t1018 + ((t347 / 0.2e1 + t1055 - t840) * t705 + (-t346 / 0.2e1 - t273 / 0.2e1 - t839) * t704) * t860, t180 * qJD(1) + t127 * qJD(3) + t7 * qJD(4) + t13 * qJD(5) + t731 * qJD(6) + (t422 * t521 - t1145 * t522 + t322 / 0.2e1 + t323 / 0.2e1) * t1014; 0, 0, 0 ((-m(5) * t645 - m(6) * t1106 + t811 * t879) * t705 + (-m(5) * t644 - m(6) * t565 + t811 * t880) * t704 + 0.2e1 * (m(6) / 0.2e1 + t1068) * (-t704 * t649 + t604)) * qJD(4) + t92 * qJD(5) + t233 * qJD(6), t92 * qJD(4) + (t391 * t1068 + (-m(6) * t536 - m(7) * t624) * t950 + (-m(6) * t537 - m(7) * t625) * t959) * qJD(5) + t387, t233 * qJD(4) + t391 * t842 + t387; t757 * qJD(1) + t104 * qJD(4) + t28 * qJD(5) + t126 * qJD(6) + t278 * t1139 + t1147 * t321, 0, 0, t104 * qJD(1) + ((-m(6) * t543 - m(7) * t461) * t705 - m(5) * t797 * t803 / 0.2e1 + (m(6) * t541 + m(7) * t459) * t704) * qJD(4) + t81 * qJD(5) + t207 * qJD(6), t28 * qJD(1) + t81 * qJD(4) + ((-m(6) * t465 - m(7) * t392) * t705 + (m(6) * t466 + m(7) * t393) * t704) * qJD(5) + t338, t126 * qJD(1) + t207 * qJD(4) + t345 * t842 + t338; qJD(1) * t1090 + t105 * qJD(3) + t12 * qJD(4) + t2 * qJD(5) + t5 * qJD(6) + t228 * t1139 + t265 * t1147 + t722 * t861, -qJD(5) * t91 - qJD(6) * t232, qJD(1) * t105 + qJD(5) * t82 + qJD(6) * t208, t12 * qJD(1) + t18 * qJD(5) + t24 * qJD(6) + ((t264 * t287 + t459 * t460 - t461 * t462) * m(7) + (t354 * t382 + t541 * t542 - t543 * t544) * m(6) + m(5) * ((-t705 * t737 + (-t705 * rSges(5,3) - t704 * t797) * t704) * (-t644 * t704 - t645 * t705) - t862 * t677 * t797) + (t111 + t130 - t702 * t639 + (t1087 * t705 + (-t1088 + t638) * t704) * t705) * t1053 + (t110 + t129 + t703 * t638 + (t1088 * t704 + (-t1087 - t639) * t705) * t704) * t1051) * qJD(4), t2 * qJD(1) + t18 * qJD(4) + t4 * qJD(6) + (-t187 / 0.2e1 + (-t62 / 0.2e1 + t135 / 0.2e1) * t705 + (t61 / 0.2e1 - t134 / 0.2e1) * t704) * t860 + (t264 * t341 + t276 * t350 + t332 * t554 - t333 * t553 - t392 * t462 + t393 * t460 - t32 / 0.4e1) * t1015 + (t354 * t411 + t388 * t421 - t465 * t544 + t466 * t542 - t108 / 0.4e1 + (t1151 * t704 - t444 * t705) * t661) * t1018 - t782 + (t101 * t1051 + t102 * t1053 + t728 + (-t79 / 0.2e1 + (t273 / 0.2e1 + t163 / 0.2e1) * t705 + (t1055 + t1056) * t704) * t712) * qJD(5), t5 * qJD(1) + t24 * qJD(4) + t4 * qJD(5) + t728 * qJD(6) + (t371 * t405 + t202 / 0.2e1 + t360 / 0.2e1 + t361 / 0.2e1 - t95 / 0.4e1 + (t1145 * t704 - t422 * t705) * t655) * t1014 + t760; (-t930 / 0.2e1 + t735) * qJD(1) + t29 * qJD(3) + t1 * qJD(4) + t11 * qJD(5) + t14 * qJD(6) + (-t221 / 0.4e1 - t906 * t397) * t1017 + t268 * t1147 - t1089 * t861, qJD(4) * t91, qJD(1) * t29 - qJD(4) * t82, t1 * qJD(1) + (t730 + t130 * t823 + t129 * t832 + (t342 * t705 - t343 * t704) * t1048 + (t266 * t705 + t267 * t704) * t1049 + t61 * t1051 + t62 * t1053 + t191 * t828 + t190 * t834 - t746) * qJD(4) + t8 * qJD(5) + t9 * qJD(6) + (t186 * t264 - t245 * t462 + t246 * t460 + t276 * t287 + t332 * t461 + t333 * t459 - t131 / 0.4e1) * t1016 + (t313 * t354 - t365 * t544 + t366 * t542 + t388 * t382 + t444 * t543 + t1151 * t541 - t216 / 0.4e1) * t1019 + t782, t11 * qJD(1) + t8 * qJD(4) + ((t276 * t341 + t332 * t392 + t333 * t393) * m(7) + t403 * t1049 + (t1151 * t466 + t388 * t411 + t444 * t465) * m(6) + (t101 * t1054 + t102 * t1051 + (-t273 * t704 + t274 * t705) * t1049) * t715 + t848) * qJD(5) + t19 * qJD(6), t14 * qJD(1) + t9 * qJD(4) + t19 * qJD(5) + t848 * qJD(6) + (-t155 / 0.4e1 + t738 + t799) * t1014; t735 * qJD(1) + t128 * qJD(3) + t6 * qJD(4) + t15 * qJD(5) + qJD(6) * t798 + t1139 * t242, qJD(4) * t232, qJD(1) * t128 - qJD(4) * t208, t6 * qJD(1) + t730 * qJD(4) + t10 * qJD(5) + t16 * qJD(6) + (t264 * t291 + t287 * t371 - t351 * t462 + t352 * t460 + t422 * t461 + t1145 * t459 - t156 / 0.4e1) * t1016 - t760, t15 * qJD(1) + t10 * qJD(4) + t848 * qJD(5) + t20 + (t371 * t341 + t422 * t392 + t1145 * t393 - t94 / 0.4e1 + t799) * t1015, qJD(1) * t798 + qJD(4) * t16 + qJD(5) * t21 + t20;];
Cq  = t17;