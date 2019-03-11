% Calculate vector of centrifugal and Coriolis load on the joints for
% S6RPRRRR2
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [11x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d3,d4,d5,d6,theta2]';
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
% tauc [6x1]
%   joint torques required to compensate Coriolis and centrifugal load

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 06:59
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauc = S6RPRRRR2_coriolisvecJ_fixb_slag_vp1(qJ, qJD, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(11,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRRRR2_coriolisvecJ_fixb_slag_vp1: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RPRRRR2_coriolisvecJ_fixb_slag_vp1: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6RPRRRR2_coriolisvecJ_fixb_slag_vp1: pkin has to be [11x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RPRRRR2_coriolisvecJ_fixb_slag_vp1: m has to be [7x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [7,3]), ...
  'S6RPRRRR2_coriolisvecJ_fixb_slag_vp1: rSges has to be [7x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [7 6]), ...
  'S6RPRRRR2_coriolisvecJ_fixb_slag_vp1: Icges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From coriolisvec_joint_fixb_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 06:57:01
% EndTime: 2019-03-09 06:58:20
% DurationCPUTime: 66.16s
% Computational Cost: add. (76522->1478), mult. (60697->1982), div. (0->0), fcn. (56688->12), ass. (0->727)
t638 = sin(qJ(1));
t1063 = pkin(1) * t638;
t633 = qJ(1) + pkin(11);
t624 = cos(t633);
t643 = -pkin(8) - pkin(7);
t600 = t624 * t643;
t640 = cos(qJ(3));
t620 = pkin(3) * t640 + pkin(2);
t623 = sin(t633);
t924 = -t623 * t620 - t600;
t804 = t924 - t1063;
t635 = qJ(3) + qJ(4);
t628 = cos(t635);
t982 = t623 * t628;
t891 = rSges(5,1) * t982;
t1143 = -t891 + t804;
t639 = cos(qJ(5));
t962 = t628 * t639;
t636 = sin(qJ(5));
t973 = t624 * t636;
t494 = t623 * t962 - t973;
t469 = Icges(6,4) * t494;
t963 = t628 * t636;
t971 = t624 * t639;
t493 = t623 * t963 + t971;
t626 = sin(t635);
t984 = t623 * t626;
t294 = -Icges(6,2) * t493 + Icges(6,6) * t984 + t469;
t468 = Icges(6,4) * t493;
t298 = -Icges(6,1) * t494 - Icges(6,5) * t984 + t468;
t1141 = t294 * t636 + t298 * t639;
t291 = Icges(6,5) * t494 - Icges(6,6) * t493 + Icges(6,3) * t984;
t133 = -t1141 * t626 - t291 * t628;
t634 = qJ(5) + qJ(6);
t627 = cos(t634);
t965 = t627 * t628;
t625 = sin(t634);
t977 = t624 * t625;
t452 = t623 * t965 - t977;
t437 = Icges(7,4) * t452;
t969 = t625 * t628;
t975 = t624 * t627;
t451 = t623 * t969 + t975;
t249 = -Icges(7,2) * t451 + Icges(7,6) * t984 + t437;
t436 = Icges(7,4) * t451;
t253 = -Icges(7,1) * t452 - Icges(7,5) * t984 + t436;
t1142 = t249 * t625 + t253 * t627;
t246 = Icges(7,5) * t452 - Icges(7,6) * t451 + Icges(7,3) * t984;
t127 = -t1142 * t626 - t246 * t628;
t631 = qJD(5) + qJD(6);
t816 = t628 * t631;
t1109 = qJD(1) - t816;
t115 = t246 * t984 - t249 * t451 - t253 * t452;
t983 = t623 * t627;
t453 = -t624 * t969 + t983;
t985 = t623 * t625;
t454 = t624 * t965 + t985;
t976 = t624 * t626;
t248 = Icges(7,5) * t454 + Icges(7,6) * t453 + Icges(7,3) * t976;
t1032 = Icges(7,4) * t454;
t251 = Icges(7,2) * t453 + Icges(7,6) * t976 + t1032;
t438 = Icges(7,4) * t453;
t254 = Icges(7,1) * t454 + Icges(7,5) * t976 + t438;
t116 = t248 * t984 - t451 * t251 + t452 * t254;
t772 = Icges(7,5) * t627 - Icges(7,6) * t625;
t422 = -Icges(7,3) * t628 + t626 * t772;
t1030 = Icges(7,4) * t627;
t776 = -Icges(7,2) * t625 + t1030;
t424 = -Icges(7,6) * t628 + t626 * t776;
t1031 = Icges(7,4) * t625;
t782 = Icges(7,1) * t627 - t1031;
t426 = -Icges(7,5) * t628 + t626 * t782;
t168 = t422 * t984 - t424 * t451 + t426 * t452;
t603 = qJD(3) * t623;
t526 = qJD(4) * t623 + t603;
t903 = qJD(5) * t626;
t447 = t624 * t903 + t526;
t900 = qJD(6) * t626;
t391 = t624 * t900 + t447;
t632 = qJD(3) + qJD(4);
t527 = t624 * t632;
t448 = -t623 * t903 + t527;
t392 = -t623 * t900 + t448;
t60 = t1109 * t168 - t115 * t392 + t116 * t391;
t121 = t291 * t984 - t294 * t493 - t298 * t494;
t979 = t623 * t639;
t495 = -t624 * t963 + t979;
t981 = t623 * t636;
t496 = t624 * t962 + t981;
t293 = Icges(6,5) * t496 + Icges(6,6) * t495 + Icges(6,3) * t976;
t1035 = Icges(6,4) * t496;
t296 = Icges(6,2) * t495 + Icges(6,6) * t976 + t1035;
t470 = Icges(6,4) * t495;
t299 = Icges(6,1) * t496 + Icges(6,5) * t976 + t470;
t122 = t293 * t984 - t493 * t296 + t494 * t299;
t774 = Icges(6,5) * t639 - Icges(6,6) * t636;
t455 = -Icges(6,3) * t628 + t626 * t774;
t1033 = Icges(6,4) * t639;
t778 = -Icges(6,2) * t636 + t1033;
t457 = -Icges(6,6) * t628 + t626 * t778;
t1034 = Icges(6,4) * t636;
t784 = Icges(6,1) * t639 - t1034;
t459 = -Icges(6,5) * t628 + t626 * t784;
t172 = t455 * t984 - t457 * t493 + t459 * t494;
t902 = qJD(5) * t628;
t592 = qJD(1) - t902;
t68 = -t121 * t448 + t122 * t447 + t172 * t592;
t123 = t291 * t976 + t495 * t294 - t298 * t496;
t124 = t293 * t976 + t495 * t296 + t496 * t299;
t173 = t455 * t976 + t457 * t495 + t459 * t496;
t69 = -t123 * t448 + t124 * t447 + t173 * t592;
t117 = t246 * t976 + t453 * t249 - t253 * t454;
t118 = t248 * t976 + t453 * t251 + t454 * t254;
t169 = t422 * t976 + t424 * t453 + t426 * t454;
t61 = t1109 * t169 - t117 * t392 + t118 * t391;
t645 = qJD(1) ^ 2;
t909 = qJD(1) * t626;
t869 = t623 * t909;
t968 = t626 * t632;
t687 = t1109 * t627 + t625 * t968;
t908 = qJD(1) * t628;
t817 = -t631 + t908;
t208 = t624 * t687 + t817 * t985;
t686 = t1109 * t625 - t627 * t968;
t209 = t624 * t686 - t817 * t983;
t964 = t628 * t632;
t882 = t624 * t964;
t881 = t209 * rSges(7,1) + t208 * rSges(7,2) + rSges(7,3) * t882;
t145 = -rSges(7,3) * t869 + t881;
t1060 = pkin(9) * t626;
t619 = pkin(5) * t639 + pkin(4);
t1056 = pkin(4) - t619;
t1111 = t1056 * t628;
t531 = pkin(9) * t882;
t851 = t1056 * t626;
t642 = -pkin(10) - pkin(9);
t1047 = pkin(5) * qJD(5);
t887 = t639 * t1047;
t895 = pkin(5) * t973;
t871 = qJD(1) * t895 + t623 * t887 + t642 * t869;
t888 = t636 * t1047;
t911 = qJD(1) * t623;
t961 = t628 * t642;
t164 = -t531 + (t1060 + t1111) * t911 + (-t628 * t888 + (t851 - t961) * t632) * t624 + t871;
t1134 = t145 + t164;
t789 = rSges(7,1) * t452 - rSges(7,2) * t451;
t258 = rSges(7,3) * t984 + t789;
t1061 = pkin(4) * t628;
t1055 = pkin(9) + t642;
t1113 = t1055 * t626;
t520 = t619 * t982;
t284 = t895 - t520 + (t1061 + t1113) * t623;
t1133 = t258 - t284;
t260 = t454 * rSges(7,1) + t453 * rSges(7,2) + rSges(7,3) * t976;
t974 = t624 * t628;
t568 = pkin(4) * t974;
t486 = pkin(9) * t976 + t568;
t591 = pkin(5) * t981;
t927 = t619 * t974 + t591;
t966 = t626 * t642;
t285 = -t624 * t966 - t486 + t927;
t1132 = t260 + t285;
t788 = rSges(7,1) * t627 - rSges(7,2) * t625;
t428 = -rSges(7,3) * t628 + t626 * t788;
t700 = -t1055 * t628 + t851;
t1131 = -t700 + t428;
t1022 = Icges(5,6) * t624;
t415 = Icges(5,4) * t982 - Icges(5,2) * t984 - t1022;
t618 = Icges(5,4) * t628;
t541 = Icges(5,1) * t626 + t618;
t1130 = -t541 * t623 - t415;
t779 = -Icges(5,2) * t626 + t618;
t416 = Icges(5,6) * t623 + t624 * t779;
t1129 = -t541 * t624 - t416;
t1036 = Icges(5,4) * t626;
t542 = Icges(5,1) * t628 - t1036;
t418 = Icges(5,5) * t623 + t542 * t624;
t539 = Icges(5,2) * t628 + t1036;
t1128 = -t539 * t624 + t418;
t1127 = t541 + t779;
t1028 = Icges(5,5) * t624;
t560 = Icges(5,4) * t984;
t417 = Icges(5,1) * t982 - t1028 - t560;
t996 = t415 * t626;
t750 = -t417 * t628 + t996;
t712 = t750 * t623;
t538 = Icges(5,5) * t628 - Icges(5,6) * t626;
t414 = Icges(5,3) * t623 + t538 * t624;
t945 = t623 * t414 + t418 * t974;
t1126 = -t712 - t945;
t910 = qJD(1) * t624;
t848 = t910 / 0.2e1;
t854 = t964 / 0.2e1;
t1125 = t623 * t854 + t626 * t848;
t1124 = -t869 / 0.2e1 + t624 * t854;
t867 = t624 * t909;
t884 = t623 * t964;
t703 = t867 + t884;
t793 = rSges(6,1) * t494 - rSges(6,2) * t493;
t300 = rSges(6,3) * t984 + t793;
t792 = rSges(6,1) * t639 - rSges(6,2) * t636;
t462 = -rSges(6,3) * t628 + t626 * t792;
t1123 = -t300 * t592 - t448 * t462;
t1122 = -t1109 * t258 + t284 * t592 - t392 * t428 + t448 * t700;
t1121 = 0.2e1 * qJD(3);
t210 = t623 * t687 - t817 * t977;
t211 = t623 * t686 + t817 * t975;
t790 = rSges(7,1) * t211 + rSges(7,2) * t210;
t146 = rSges(7,3) * t703 + t790;
t411 = -t1111 - t1113;
t885 = t623 * t968;
t530 = pkin(4) * t885;
t886 = t619 * t968;
t165 = -t624 * t887 + t530 + (-t886 + (-t1055 * t632 - t888) * t628) * t623 + (t411 * t624 + t591) * qJD(1);
t897 = qJD(1) * qJD(3);
t595 = t623 * t897;
t896 = qJD(1) * qJD(4);
t517 = t623 * t896 + t595;
t850 = qJD(1) * t903;
t901 = qJD(5) * t632;
t858 = t628 * t901;
t351 = t623 * t858 + t624 * t850 + t517;
t219 = qJD(6) * t703 + t351;
t723 = t788 * t628;
t787 = -rSges(7,1) * t625 - rSges(7,2) * t627;
t240 = t632 * t723 + (rSges(7,3) * t632 + t631 * t787) * t626;
t311 = t411 * t632 - t626 * t888;
t533 = t631 * t626;
t497 = t632 * t533;
t323 = pkin(9) * t703 + qJD(1) * t568 - t530;
t546 = t1060 + t1061;
t502 = t546 * t632;
t545 = pkin(4) * t626 - pkin(9) * t628;
t637 = sin(qJ(3));
t905 = qJD(3) * t637;
t890 = pkin(3) * t905;
t571 = t623 * t890;
t641 = cos(qJ(1));
t630 = t641 * pkin(1);
t892 = t645 * t630;
t960 = t640 * qJD(3) ^ 2;
t693 = -pkin(3) * t624 * t960 + qJD(1) * t571 - t892;
t1057 = pkin(2) - t620;
t615 = t623 * pkin(7);
t920 = t643 * t911 + t571;
t341 = (-t1057 * t624 - t615) * qJD(1) - t920;
t536 = t624 * pkin(2) + t615;
t519 = t536 * qJD(1);
t947 = -t341 - t519;
t654 = t517 * t545 + (-t323 + t947) * qJD(1) - t527 * t502 + t693;
t859 = t626 * t901;
t35 = -t1109 * t146 - t592 * t165 + t219 * t428 - t392 * t240 - t497 * t258 + t284 * t859 - t448 * t311 - t351 * t700 + t654;
t562 = t624 * t620;
t821 = -t623 * t643 + t562;
t406 = t821 - t536;
t841 = t536 + t630;
t805 = t406 + t841;
t683 = (t486 + t805) * qJD(1) - t526 * t545 - t571;
t98 = t1109 * t260 + t285 * t592 - t391 * t428 + t447 * t700 + t683;
t1120 = qJD(1) * t98 + t35;
t616 = t624 * pkin(7);
t535 = pkin(2) * t623 - t616;
t405 = t535 + t924;
t522 = qJD(1) * t535;
t1119 = qJD(1) * t405 - t522;
t862 = t624 * t905;
t815 = pkin(3) * t862;
t1118 = -t527 * t545 - t815;
t611 = t623 * rSges(4,3);
t970 = t624 * t640;
t972 = t624 * t637;
t450 = rSges(4,1) * t970 - rSges(4,2) * t972 + t611;
t1117 = t450 + t841;
t1116 = -rSges(5,2) * t984 - t624 * rSges(5,3);
t837 = t624 * rSges(3,1) - rSges(3,2) * t623;
t1115 = t630 + t837;
t381 = t462 * t623;
t1051 = rSges(6,3) * t626;
t724 = t792 * t628;
t463 = t724 + t1051;
t503 = t545 * t911;
t483 = t545 * t623;
t819 = qJD(1) * t483 - t527 * t546;
t861 = t623 * t902;
t1114 = t300 * t903 - t381 * t592 + t448 * t463 + t503 - t819 + (-t861 + t911) * t462;
t629 = Icges(4,4) * t640;
t780 = -Icges(4,2) * t637 + t629;
t579 = Icges(4,1) * t637 + t629;
t419 = t891 + t1116;
t610 = t623 * rSges(5,3);
t922 = rSges(5,1) * t974 + t610;
t420 = -rSges(5,2) * t976 + t922;
t906 = qJD(3) * t624;
t853 = -t405 * t603 + t406 * t906 + qJD(2);
t150 = t419 * t526 + t420 * t527 + t853;
t543 = rSges(5,1) * t626 + rSges(5,2) * t628;
t705 = -t527 * t543 - t815;
t842 = -t535 - t1063;
t806 = t405 + t842;
t174 = (-t419 + t806) * qJD(1) + t705;
t175 = -t526 * t543 - t571 + (t420 + t805) * qJD(1);
t481 = t543 * t623;
t482 = t543 * t624;
t1052 = rSges(5,2) * t626;
t544 = rSges(5,1) * t628 - t1052;
t1103 = -t174 * (qJD(1) * t481 - t527 * t544) - t150 * (-t526 * t481 - t482 * t527) - t175 * (-qJD(1) * t482 - t526 * t544);
t1019 = Icges(4,3) * t624;
t978 = t623 * t640;
t980 = t623 * t637;
t441 = Icges(4,5) * t978 - Icges(4,6) * t980 - t1019;
t1029 = Icges(4,5) * t624;
t588 = Icges(4,4) * t980;
t445 = Icges(4,1) * t978 - t1029 - t588;
t1023 = Icges(4,6) * t624;
t443 = Icges(4,4) * t978 - Icges(4,2) * t980 - t1023;
t993 = t443 * t637;
t747 = -t445 * t640 + t993;
t184 = -t441 * t624 - t623 * t747;
t1018 = Icges(5,3) * t624;
t537 = Icges(5,5) * t626 + Icges(5,6) * t628;
t989 = t537 * t624;
t995 = t416 * t626;
t1102 = -t632 * t989 + (-t418 * t628 - t538 * t623 + t1018 + t995) * qJD(1);
t914 = qJD(1) * t414;
t990 = t537 * t623;
t1101 = qJD(1) * t750 - t632 * t990 + t914;
t576 = Icges(4,5) * t640 - Icges(4,6) * t637;
t575 = Icges(4,5) * t637 + Icges(4,6) * t640;
t708 = qJD(3) * t575;
t1037 = Icges(4,4) * t637;
t580 = Icges(4,1) * t640 - t1037;
t446 = Icges(4,5) * t623 + t580 * t624;
t444 = Icges(4,6) * t623 + t624 * t780;
t992 = t444 * t637;
t746 = -t446 * t640 + t992;
t1100 = -t624 * t708 + (-t576 * t623 + t1019 + t746) * qJD(1);
t442 = Icges(4,3) * t623 + t576 * t624;
t913 = qJD(1) * t442;
t1099 = qJD(1) * t747 - t623 * t708 + t913;
t743 = t539 * t626 - t541 * t628;
t1098 = qJD(1) * t743 + t538 * t632;
t577 = Icges(4,2) * t640 + t1037;
t742 = t577 * t637 - t579 * t640;
t1097 = qJD(1) * t742 + t576 * qJD(3);
t931 = -Icges(4,2) * t978 + t445 - t588;
t933 = t579 * t623 + t443;
t1096 = -t637 * t931 - t640 * t933;
t302 = t496 * rSges(6,1) + t495 * rSges(6,2) + rSges(6,3) * t976;
t120 = t302 * t592 - t447 * t462 + t683;
t967 = t626 * t636;
t685 = t592 * t639 + t632 * t967;
t814 = -qJD(5) + t908;
t268 = t624 * t685 + t814 * t981;
t684 = t592 * t636 - t639 * t968;
t269 = t624 * t684 - t814 * t979;
t877 = t269 * rSges(6,1) + t268 * rSges(6,2) + rSges(6,3) * t882;
t162 = -rSges(6,3) * t869 + t877;
t791 = -rSges(6,1) * t636 - rSges(6,2) * t639;
t310 = t632 * t724 + (rSges(6,3) * t632 + qJD(5) * t791) * t626;
t596 = t624 * t897;
t518 = t624 * t896 + t596;
t872 = t624 * t858 + t518;
t352 = -t623 * t850 + t872;
t883 = t624 * t968;
t702 = -t623 * t908 - t883;
t322 = pkin(4) * t702 - pkin(9) * t869 + t531;
t601 = pkin(7) * t910;
t340 = -t815 - t601 + (t1057 * t623 - t600) * qJD(1);
t893 = t645 * t1063;
t807 = qJD(1) * (-pkin(2) * t911 + t601) - t893;
t655 = qJD(1) * t340 + (-t596 * t637 - t623 * t960) * pkin(3) + t807;
t652 = qJD(1) * t322 - t526 * t502 - t518 * t545 + t655;
t74 = t162 * t592 + t302 * t859 - t310 * t447 - t352 * t462 + t652;
t270 = t623 * t685 - t814 * t973;
t271 = t623 * t684 + t814 * t971;
t794 = rSges(6,1) * t271 + rSges(6,2) * t270;
t163 = rSges(6,3) * t703 + t794;
t75 = -t163 * t592 - t300 * t859 - t310 * t448 + t351 * t462 + t654;
t1095 = (qJD(1) * t120 + t75) * t624 + t74 * t623;
t713 = t772 * t628;
t748 = -t424 * t625 + t426 * t627;
t755 = -t251 * t625 + t254 * t627;
t1094 = t391 * (-t422 * t624 - t755) - t392 * (-t422 * t623 + t1142) + t1109 * (Icges(7,3) * t626 + t713 - t748);
t775 = -Icges(7,2) * t627 - t1031;
t1093 = t391 * (-Icges(7,2) * t454 + t254 + t438) - t392 * (-Icges(7,2) * t452 - t253 - t436) + t1109 * (t775 * t626 + t426);
t714 = t774 * t628;
t744 = -t457 * t636 + t459 * t639;
t752 = -t296 * t636 + t299 * t639;
t1092 = t447 * (-t455 * t624 - t752) - t448 * (-t455 * t623 + t1141) + t592 * (Icges(6,3) * t626 + t714 - t744);
t777 = -Icges(6,2) * t639 - t1034;
t1091 = t447 * (-Icges(6,2) * t496 + t299 + t470) - t448 * (-Icges(6,2) * t494 - t298 - t468) + t592 * (t777 * t626 + t459);
t1090 = qJD(1) * t1127 + t526 * t1128 - t527 * (-Icges(5,2) * t982 + t417 - t560);
t1089 = t219 / 0.2e1;
t220 = qJD(6) * t882 - t631 * t869 + t872;
t1088 = t220 / 0.2e1;
t1087 = t351 / 0.2e1;
t1086 = t352 / 0.2e1;
t1085 = -t391 / 0.2e1;
t1084 = t391 / 0.2e1;
t1083 = -t392 / 0.2e1;
t1082 = t392 / 0.2e1;
t1081 = -t447 / 0.2e1;
t1080 = t447 / 0.2e1;
t1079 = -t448 / 0.2e1;
t1078 = t448 / 0.2e1;
t1077 = t517 / 0.2e1;
t1076 = t518 / 0.2e1;
t1075 = -t1109 / 0.2e1;
t1074 = t1109 / 0.2e1;
t1073 = -t526 / 0.2e1;
t1072 = t526 / 0.2e1;
t1071 = -t527 / 0.2e1;
t1070 = t527 / 0.2e1;
t1069 = -t592 / 0.2e1;
t1068 = t592 / 0.2e1;
t1067 = t623 / 0.2e1;
t1066 = -t624 / 0.2e1;
t1065 = -t628 / 0.2e1;
t1064 = -rSges(6,3) - pkin(9);
t1062 = pkin(3) * t637;
t1059 = -qJD(1) / 0.2e1;
t1058 = qJD(1) / 0.2e1;
t1054 = rSges(4,1) * t640;
t1049 = rSges(7,3) * t626;
t140 = Icges(7,5) * t211 + Icges(7,6) * t210 + Icges(7,3) * t703;
t142 = Icges(7,4) * t211 + Icges(7,2) * t210 + Icges(7,6) * t703;
t144 = Icges(7,1) * t211 + Icges(7,4) * t210 + Icges(7,5) * t703;
t41 = (-t1142 * t632 - t140) * t628 + (t246 * t632 + (-t249 * t631 + t144) * t627 + (t253 * t631 - t142) * t625) * t626;
t1046 = t41 * t392;
t701 = -t869 + t882;
t139 = Icges(7,5) * t209 + Icges(7,6) * t208 + Icges(7,3) * t701;
t141 = Icges(7,4) * t209 + Icges(7,2) * t208 + Icges(7,6) * t701;
t143 = Icges(7,1) * t209 + Icges(7,4) * t208 + Icges(7,5) * t701;
t42 = (t632 * t755 - t139) * t628 + (t248 * t632 + (-t251 * t631 + t143) * t627 + (-t254 * t631 - t141) * t625) * t626;
t1045 = t42 * t391;
t157 = Icges(6,5) * t271 + Icges(6,6) * t270 + Icges(6,3) * t703;
t159 = Icges(6,4) * t271 + Icges(6,2) * t270 + Icges(6,6) * t703;
t161 = Icges(6,1) * t271 + Icges(6,4) * t270 + Icges(6,5) * t703;
t48 = (-t1141 * t632 - t157) * t628 + (-t159 * t636 + t161 * t639 + t291 * t632 + (-t294 * t639 + t298 * t636) * qJD(5)) * t626;
t1044 = t48 * t448;
t156 = Icges(6,5) * t269 + Icges(6,6) * t268 + Icges(6,3) * t701;
t158 = Icges(6,4) * t269 + Icges(6,2) * t268 + Icges(6,6) * t701;
t160 = Icges(6,1) * t269 + Icges(6,4) * t268 + Icges(6,5) * t701;
t49 = (t632 * t752 - t156) * t628 + (-t158 * t636 + t160 * t639 + t293 * t632 + (-t296 * t639 - t299 * t636) * qJD(5)) * t626;
t1043 = t49 * t447;
t1042 = -rSges(7,3) + t642;
t182 = -t422 * t628 + t626 * t748;
t771 = -Icges(7,5) * t625 - Icges(7,6) * t627;
t237 = t632 * t713 + (Icges(7,3) * t632 + t631 * t771) * t626;
t715 = t776 * t628;
t238 = t632 * t715 + (Icges(7,6) * t632 + t631 * t775) * t626;
t717 = t782 * t628;
t781 = -Icges(7,1) * t625 - t1030;
t239 = t632 * t717 + (Icges(7,5) * t632 + t631 * t781) * t626;
t91 = (t632 * t748 - t237) * t628 + (t422 * t632 + (-t424 * t631 + t239) * t627 + (-t426 * t631 - t238) * t625) * t626;
t1041 = t1109 * t91 + t182 * t497;
t773 = -Icges(6,5) * t636 - Icges(6,6) * t639;
t304 = t632 * t714 + (Icges(6,3) * t632 + qJD(5) * t773) * t626;
t716 = t778 * t628;
t305 = t632 * t716 + (Icges(6,6) * t632 + qJD(5) * t777) * t626;
t718 = t784 * t628;
t783 = -Icges(6,1) * t636 - t1033;
t306 = t632 * t718 + (Icges(6,5) * t632 + qJD(5) * t783) * t626;
t101 = (t632 * t744 - t304) * t628 + (-t305 * t636 + t306 * t639 + t455 * t632 + (-t457 * t639 - t459 * t636) * qJD(5)) * t626;
t188 = -t455 * t628 + t626 * t744;
t1040 = t101 * t592 + t188 * t859;
t484 = t546 * t623;
t711 = t526 * t484 + t486 * t527 + t853;
t81 = t258 * t391 + t260 * t392 - t284 * t447 + t285 * t448 + t711;
t1015 = qJD(1) * t81;
t1011 = t120 * t623;
t1008 = t127 * t219;
t128 = -t248 * t628 + t626 * t755;
t1007 = t128 * t220;
t1006 = t133 * t351;
t134 = -t293 * t628 + t626 * t752;
t1005 = t134 * t352;
t917 = rSges(4,2) * t980 + t624 * rSges(4,3);
t449 = rSges(4,1) * t978 - t917;
t581 = rSges(4,1) * t637 + rSges(4,2) * t640;
t863 = t581 * t906;
t256 = -t863 + (-t449 + t842) * qJD(1);
t1002 = t256 * t623;
t1001 = t256 * t624;
t865 = t581 * t603;
t257 = qJD(1) * t1117 - t865;
t511 = t581 * t624;
t1000 = t257 * t511;
t988 = t539 * t632;
t987 = t575 * t623;
t986 = t575 * t624;
t959 = t146 * t976 + t258 * t882;
t957 = t260 * t968 + t428 * t869;
t956 = -t240 - t311;
t949 = -t302 - t486;
t948 = -t310 - t502;
t413 = Icges(5,5) * t982 - Icges(5,6) * t984 - t1018;
t946 = -t623 * t413 - t417 * t974;
t943 = -t623 * t405 + t624 * t406;
t942 = t623 * t419 + t624 * t420;
t941 = -t623 * t441 - t445 * t970;
t940 = t623 * t442 + t446 * t970;
t934 = t623 * t484 + t624 * t486;
t932 = -t579 * t624 - t444;
t930 = -t577 * t624 + t446;
t925 = rSges(5,2) * t869 + rSges(5,3) * t910;
t923 = t562 + t630;
t907 = qJD(1) * t637;
t921 = rSges(4,2) * t623 * t907 + rSges(4,3) * t910;
t919 = -t577 + t580;
t918 = t579 + t780;
t110 = t300 * t447 + t302 * t448 + t711;
t916 = qJD(1) * t110;
t912 = qJD(1) * t576;
t904 = qJD(3) * t640;
t242 = -t623 * t743 - t989;
t899 = t242 * qJD(1);
t332 = -t623 * t742 - t986;
t898 = t332 * qJD(1);
t894 = pkin(5) * t967;
t889 = pkin(3) * t904;
t280 = rSges(5,1) * t702 - rSges(5,2) * t882 + t925;
t281 = -t632 * t481 + (t544 * t624 + t610) * qJD(1);
t880 = t624 * t280 + t623 * t281 + t419 * t910;
t879 = -t502 + t956;
t878 = -t486 - t1132;
t876 = t624 * t322 + t623 * t323 + t484 * t910;
t875 = t624 * t340 + t623 * t341 - t405 * t910;
t874 = t1131 * t911 + t503;
t870 = t1064 * t626;
t866 = t624 * t907;
t860 = t624 * t902;
t857 = t984 / 0.2e1;
t856 = t976 / 0.2e1;
t855 = t968 / 0.2e1;
t852 = -pkin(2) - t1054;
t849 = t911 / 0.2e1;
t847 = -t603 / 0.2e1;
t844 = t906 / 0.2e1;
t840 = -t543 - t1062;
t839 = -t545 - t1062;
t836 = (-t623 ^ 2 - t624 ^ 2) * t637;
t318 = -rSges(7,1) * t451 - rSges(7,2) * t452;
t319 = rSges(7,1) * t453 - rSges(7,2) * t454;
t835 = t391 * t318 + t319 * t392;
t490 = t787 * t626;
t834 = t1109 * t319 - t391 * t490;
t833 = (-t623 * t779 + t1022) * qJD(1) + t1128 * t632;
t832 = qJD(1) * t416 + t417 * t632 - t623 * t988;
t831 = (-t542 * t623 + t1028) * qJD(1) + t1129 * t632;
t830 = qJD(1) * t418 + t1130 * t632;
t829 = -t1109 * t318 - t392 * t490;
t485 = t545 * t624;
t828 = -t526 * t483 - t485 * t527;
t402 = t446 * t978;
t827 = t442 * t624 - t402;
t826 = -t413 + t995;
t825 = -t441 + t992;
t824 = -qJD(1) * t485 - t526 * t546;
t823 = t1127 * t632;
t822 = t542 * t632 - t988;
t813 = t628 * t146 + t240 * t984 + t428 * t703;
t812 = t623 * t300 + t624 * t302 + t934;
t501 = t544 * t632;
t800 = -t501 - t889;
t799 = -t502 - t889;
t365 = t418 * t982;
t798 = t416 * t984 - t365;
t532 = rSges(3,1) * t623 + rSges(3,2) * t624;
t795 = -rSges(4,2) * t637 + t1054;
t770 = t115 * t624 - t116 * t623;
t769 = t115 * t623 + t116 * t624;
t768 = t117 * t624 - t118 * t623;
t767 = t117 * t623 + t118 * t624;
t667 = (-t484 + t806) * qJD(1) + t1118;
t119 = t1123 + t667;
t766 = t119 * t624 + t1011;
t765 = t121 * t624 - t122 * t623;
t764 = t121 * t623 + t122 * t624;
t763 = t123 * t624 - t124 * t623;
t762 = t123 * t623 + t124 * t624;
t761 = t127 * t624 - t128 * t623;
t760 = t127 * t623 + t128 * t624;
t759 = t133 * t624 - t134 * t623;
t758 = t133 * t623 + t134 * t624;
t757 = -t174 * t624 - t175 * t623;
t754 = -t257 * t623 - t1001;
t751 = t300 * t624 - t302 * t623;
t225 = t415 * t628 + t417 * t626;
t286 = t443 * t640 + t445 * t637;
t287 = t444 * t640 + t446 * t637;
t745 = t449 * t623 + t450 * t624;
t741 = t839 - t1131;
t740 = -t310 + t799;
t737 = t624 * t162 + t623 * t163 + t300 * t910 + t876;
t736 = t1132 * t624 + t1133 * t623 + t934;
t735 = -t619 * t628 - t1049 - t620;
t510 = t581 * t623;
t695 = t340 * t906 + t341 * t603 - t405 * t596 - t406 * t595;
t659 = t322 * t527 + t526 * t323 + t518 * t484 - t486 * t517 + t695;
t25 = t145 * t392 + t146 * t391 + t164 * t448 + t165 * t447 - t219 * t260 + t220 * t258 - t284 * t352 - t285 * t351 + t659;
t722 = t25 * t258 * t976 + t35 * (t628 * t258 + t428 * t984);
t97 = t1122 + t667;
t721 = -t1132 * t81 - t700 * t97;
t720 = t799 + t956;
t710 = qJD(3) * t579;
t709 = qJD(3) * t577;
t185 = -t444 * t980 - t827;
t707 = (-t184 * t624 + t185 * t623) * qJD(3);
t186 = -t443 * t972 - t941;
t187 = -t444 * t972 + t940;
t706 = (-t186 * t624 + t187 * t623) * qJD(3);
t699 = -t291 * t448 + t293 * t447 + t455 * t592;
t698 = (-Icges(7,5) * t451 - Icges(7,6) * t452) * t392 - (Icges(7,5) * t453 - Icges(7,6) * t454) * t391 - t771 * t626 * t1109;
t697 = (-Icges(6,5) * t493 - Icges(6,6) * t494) * t448 - (Icges(6,5) * t495 - Icges(6,6) * t496) * t447 - t773 * t626 * t592;
t696 = qJD(1) * t538 - t526 * t989 + t527 * t990;
t694 = -t637 * t930 + t640 * t932;
t692 = t876 + t1133 * t910 + t1134 * t624 + (t146 + t165) * t623;
t691 = t626 * t698;
t690 = t626 * t697;
t676 = (-t637 * t918 + t640 * t919) * qJD(1);
t673 = (Icges(7,1) * t453 - t1032 - t251) * t391 - (-Icges(7,1) * t451 - t249 - t437) * t392 + (t781 * t626 - t424) * t1109;
t671 = (Icges(6,1) * t495 - t1035 - t296) * t447 - (-Icges(6,1) * t493 - t294 - t469) * t448 + (t783 * t626 - t457) * t592;
t335 = -rSges(4,2) * t624 * t904 + (-t640 * t911 - t862) * rSges(4,1) + t921;
t336 = -qJD(3) * t510 + (t624 * t795 + t611) * qJD(1);
t669 = t335 * t624 + t336 * t623 + (t449 * t624 - t450 * t623) * qJD(1);
t668 = t1129 * t526 - t1130 * t527 + (-t539 + t542) * qJD(1);
t666 = t1118 + t1119 + (-t1063 - t484) * qJD(1);
t36 = t140 * t976 + t142 * t453 + t144 * t454 + t208 * t249 - t209 * t253 + t246 * t701;
t37 = t139 * t976 + t141 * t453 + t143 * t454 + t208 * t251 + t209 * t254 + t248 * t701;
t38 = t140 * t984 - t142 * t451 + t144 * t452 + t210 * t249 - t211 * t253 + t246 * t703;
t39 = t139 * t984 - t141 * t451 + t143 * t452 + t210 * t251 + t211 * t254 + t248 * t703;
t77 = t208 * t424 + t209 * t426 + t237 * t976 + t238 * t453 + t239 * t454 + t422 * t701;
t7 = t1109 * t77 + t117 * t219 + t118 * t220 + t169 * t497 - t36 * t392 + t37 * t391;
t70 = t1109 * t182 - t127 * t392 + t128 * t391;
t78 = t210 * t424 + t211 * t426 + t237 * t984 - t238 * t451 + t239 * t452 + t422 * t703;
t8 = t1109 * t78 + t115 * t219 + t116 * t220 + t168 * t497 - t38 * t392 + t39 * t391;
t665 = ((t632 * t767 - t77) * t628 + (qJD(1) * t768 + t169 * t632 + t36 * t623 + t37 * t624) * t626) * t1084 + (t1007 + t1008 + t1041 + t1045 - t1046) * t1065 + ((t632 * t769 - t78) * t628 + (qJD(1) * t770 + t168 * t632 + t38 * t623 + t39 * t624) * t626) * t1083 + t8 * t857 + (t1093 * t453 + t673 * t454 - t624 * t691) * t1085 + (-t1093 * t451 + t452 * t673 - t623 * t691) * t1082 + t7 * t856 + (t698 * t628 + (-t1093 * t625 + t627 * t673) * t626) * t1075 + (-t168 * t628 + t626 * t769) * t1089 + (-t169 * t628 + t626 * t767) * t1088 + t70 * t855 + t497 * (-t182 * t628 + t626 * t760) / 0.2e1 + ((t632 * t760 - t91) * t628 + (qJD(1) * t761 + t182 * t632 + t41 * t623 + t42 * t624) * t626) * t1074 + t1124 * t61 + t1125 * t60;
t664 = qJD(1) * t413 - t626 * t832 + t628 * t830;
t663 = -t626 * t833 + t628 * t831 + t914;
t662 = qJD(1) * t537 - t626 * t823 + t628 * t822;
t329 = qJD(1) * t444 - t623 * t709;
t331 = qJD(1) * t446 - t623 * t710;
t658 = qJD(1) * t441 - qJD(3) * t286 - t329 * t637 + t331 * t640;
t328 = -t624 * t709 + (-t623 * t780 + t1023) * qJD(1);
t330 = -t624 * t710 + (-t580 * t623 + t1029) * qJD(1);
t657 = -qJD(3) * t287 - t328 * t637 + t330 * t640 + t913;
t549 = t780 * qJD(3);
t550 = t580 * qJD(3);
t656 = qJD(1) * t575 - t549 * t637 + t550 * t640 + (-t577 * t640 - t579 * t637) * qJD(3);
t337 = t700 * t623;
t359 = t428 * t623;
t429 = t723 + t1049;
t466 = t623 * t816;
t653 = t1109 * t359 - t258 * t533 + t284 * t903 - t337 * t592 - t392 * t429 - t448 * t411 + t466 * t428 - t700 * t861 + t819;
t651 = t1092 * t626;
t650 = -t1090 * t626 + t668 * t628;
t382 = t462 * t624;
t649 = t110 * (t300 * t860 - t302 * t861 - t447 * t381 - t382 * t448 + t828) + t120 * (t302 * t903 - t592 * t382 - t447 * t463 - t462 * t860 + t824);
t648 = (t1109 * t422 - t246 * t392 + t248 * t391) * t628 + t1094 * t626;
t176 = -t413 * t624 - t712;
t177 = -t414 * t624 - t798;
t108 = -t176 * t527 + t177 * t526 + t899;
t178 = -t415 * t976 - t946;
t179 = -t416 * t976 + t945;
t243 = -t624 * t743 + t990;
t231 = t243 * qJD(1);
t109 = -t178 * t527 + t179 * t526 + t231;
t129 = t626 * t830 + t628 * t832;
t44 = t157 * t976 + t159 * t495 + t161 * t496 + t268 * t294 - t269 * t298 + t291 * t701;
t45 = t156 * t976 + t158 * t495 + t160 * t496 + t268 * t296 + t269 * t299 + t293 * t701;
t87 = t268 * t457 + t269 * t459 + t304 * t976 + t305 * t495 + t306 * t496 + t455 * t701;
t13 = t123 * t351 + t124 * t352 + t173 * t859 - t44 * t448 + t447 * t45 + t592 * t87;
t130 = t626 * t831 + t628 * t833;
t131 = t1098 * t623 + t662 * t624;
t132 = -t1098 * t624 + t662 * t623;
t46 = t157 * t984 - t159 * t493 + t161 * t494 + t270 * t294 - t271 * t298 + t291 * t703;
t47 = t156 * t984 - t158 * t493 + t160 * t494 + t270 * t296 + t271 * t299 + t293 * t703;
t88 = t270 * t457 + t271 * t459 + t304 * t984 - t305 * t493 + t306 * t494 + t455 * t703;
t14 = t121 * t351 + t122 * t352 + t172 * t859 + t447 * t47 - t448 * t46 + t592 * t88;
t226 = t416 * t628 + t418 * t626;
t355 = t424 * t623;
t356 = t424 * t624;
t357 = t426 * t623;
t358 = t426 * t624;
t372 = t457 * t623;
t373 = t457 * t624;
t374 = t459 * t623;
t375 = t459 * t624;
t425 = Icges(7,6) * t626 + t715;
t427 = Icges(7,5) * t626 + t717;
t458 = Icges(6,6) * t626 + t716;
t460 = Icges(6,5) * t626 + t718;
t467 = t624 * t816;
t80 = -t133 * t448 + t134 * t447 + t188 * t592;
t92 = t1101 * t623 + t664 * t624;
t93 = t1102 * t623 + t663 * t624;
t94 = -t1101 * t624 + t664 * t623;
t95 = -t1102 * t624 + t663 * t623;
t647 = (t108 + t68 + t60) * t849 - (t623 * t68 + t624 * t69) * t902 / 0.2e1 + (t109 + t69 + t61) * t848 - t80 * t903 / 0.2e1 + (qJD(1) * t132 + t176 * t517 + t177 * t518 + t526 * t95 - t527 * t94 + t14 + t8) * t1066 - t759 * t859 / 0.2e1 - t497 * t761 / 0.2e1 - t352 * t763 / 0.2e1 - t351 * t765 / 0.2e1 + (t127 * t466 + t128 * t467 + t182 * t533 - t1094 * t628 + ((t356 * t625 - t358 * t627 + t248) * t391 - (t355 * t625 - t357 * t627 + t246) * t392 + (-t425 * t625 + t427 * t627 + t422) * t1109) * t626) * t1075 + ((t356 * t451 - t358 * t452) * t391 + t116 * t467 - (t355 * t451 - t357 * t452) * t392 + t115 * t466 + (-t425 * t451 + t427 * t452) * t1109 + t168 * t533 + t648 * t623) * t1082 + ((-t356 * t453 - t358 * t454) * t391 + t118 * t467 - (-t355 * t453 - t357 * t454) * t392 + t117 * t466 + (t425 * t453 + t427 * t454) * t1109 + t169 * t533 + t648 * t624) * t1085 + (qJD(1) * t131 + t178 * t517 + t179 * t518 + t526 * t93 - t527 * t92 + t13 + t7) * t1067 + (t1090 * t628 + t668 * t626) * t1059 - t220 * t768 / 0.2e1 - t219 * t770 / 0.2e1 - t533 * t70 / 0.2e1 + (((t373 * t636 - t375 * t639 + t293) * t447 - (t372 * t636 - t374 * t639 + t291) * t448 + (-t458 * t636 + t460 * t639 + t455) * t592 + t188 * qJD(5)) * t626 + (qJD(5) * t758 - t1092) * t628) * t1069 + ((t373 * t493 - t375 * t494) * t447 - (t372 * t493 - t374 * t494) * t448 + (-t458 * t493 + t460 * t494) * t592 + (t122 * t974 + t172 * t626) * qJD(5) + ((qJD(5) * t121 + t699) * t628 + t651) * t623) * t1078 + ((-t373 * t495 - t375 * t496) * t447 - (-t372 * t495 - t374 * t496) * t448 + (t458 * t495 + t460 * t496) * t592 + (t123 * t982 + t173 * t626) * qJD(5) + ((qJD(5) * t124 + t699) * t628 + t651) * t624) * t1081 + (-t129 * t624 + t130 * t623 + (t225 * t623 + t226 * t624) * qJD(1)) * t1058 + (qJD(1) * t758 - t48 * t624 + t49 * t623) * t1068 + (t623 * t650 - t624 * t696) * t1070 + (t623 * t95 - t624 * t94 + (t176 * t623 + t177 * t624) * qJD(1)) * t1071 + (t623 * t93 - t624 * t92 + (t178 * t623 + t179 * t624) * qJD(1)) * t1072 + (t623 * t696 + t624 * t650) * t1073 + (qJD(1) * t760 - t41 * t624 + t42 * t623) * t1074 + (-t178 * t624 + t179 * t623) * t1076 + (-t176 * t624 + t177 * t623) * t1077 + (qJD(1) * t764 - t46 * t624 + t47 * t623) * t1079 + (qJD(1) * t762 - t44 * t624 + t45 * t623) * t1080 + (qJD(1) * t769 - t38 * t624 + t39 * t623) * t1083 + (qJD(1) * t767 - t36 * t624 + t37 * t623) * t1084 - t466 * t60 / 0.2e1 - t467 * t61 / 0.2e1;
t338 = t700 * t624;
t360 = t428 * t624;
t646 = t81 * (t467 * t258 - t260 * t466 - t284 * t860 - t285 * t861 + t447 * t337 + t338 * t448 - t391 * t359 - t360 * t392 + t828) + t98 * (-t1109 * t360 + t533 * t260 + t285 * t903 + t592 * t338 - t391 * t429 - t411 * t447 - t428 * t467 + t700 * t860 + t824);
t552 = t795 * qJD(3);
t515 = t791 * t626;
t472 = t495 * pkin(5);
t471 = t493 * pkin(5);
t350 = rSges(6,1) * t495 - rSges(6,2) * t496;
t349 = -rSges(6,1) * t493 - rSges(6,2) * t494;
t333 = -t624 * t742 + t987;
t324 = t333 * qJD(1);
t241 = qJD(3) * t745 + qJD(2);
t181 = -t892 - t552 * t906 + (-t336 - t519 + t865) * qJD(1);
t180 = -t552 * t603 + (t335 - t863) * qJD(1) + t807;
t167 = -t1097 * t624 + t656 * t623;
t166 = t1097 * t623 + t656 * t624;
t149 = -qJD(3) * t746 + t328 * t640 + t330 * t637;
t148 = -qJD(3) * t747 + t329 * t640 + t331 * t637;
t147 = t669 * qJD(3);
t126 = -t501 * t527 + t517 * t543 + (-t281 + t947) * qJD(1) + t693;
t125 = qJD(1) * t280 - t501 * t526 - t518 * t543 + t655;
t112 = t324 + t706;
t111 = t707 + t898;
t90 = t280 * t527 + t281 * t526 + t419 * t518 - t420 * t517 + t695;
t40 = t162 * t448 + t163 * t447 + t300 * t352 - t302 * t351 + t659;
t34 = t1109 * t145 + t592 * t164 - t220 * t428 - t391 * t240 + t497 * t260 + t285 * t859 - t447 * t311 + t352 * t700 + t652;
t1 = [(-(t1123 - t119 + t666) * t120 + t75 * (-t793 + t804) + t119 * (t530 - t794 + t920) + t74 * (t923 - t949) + t120 * (-pkin(4) * t883 + t531 - t815 + t877) + (-t74 * t643 + t75 * t870 + (t1064 * t119 * t632 - t75 * pkin(4)) * t628) * t623 + ((-t119 * t641 - t120 * t638) * pkin(1) + (t119 * (-t546 - t620 - t1051) - t120 * t643) * t624 + (-t620 + t870 - t1061) * t1011) * qJD(1)) * m(6) + (t231 + (t177 + (t414 + t996) * t624 + t798 + t946) * t527 + (-t624 * t826 - t1126 + t176) * t526) * t1070 + (t243 + t226) * t1076 + (t242 + t225) * t1077 + (-(t1122 + t666 - t97) * t98 + t35 * (-t520 - t789 + t804) + t97 * (-t790 + t920) + t34 * (t260 + t923 + t927) + t98 * (t871 + t881) + (-t34 * t966 + t98 * (-t632 * t961 - t886 - t890) + (t35 * t636 + (t639 * t97 - t963 * t98) * qJD(5)) * pkin(5)) * t624 + (-t34 * t643 + t97 * (t1042 * t632 + t888) * t628 + (t97 * t619 * t632 + t1042 * t35) * t626) * t623 + ((-t638 * t98 - t641 * t97) * pkin(1) + (-t97 * pkin(5) * t636 + t735 * t98) * t623 + (t97 * (t735 + t966) - t98 * t643) * t624) * qJD(1)) * m(7) + ((t286 + t332) * t623 + (t287 + t333) * t624) * t897 / 0.2e1 + (-qJD(3) * t742 + t549 * t640 + t550 * t637 + t626 * t822 + t628 * t823) * qJD(1) + (t181 * (t623 * t852 - t1063 + t616 + t917) + t180 * t1117 + t257 * (t601 + t921) + (t1002 * t581 - t1000) * qJD(3) + ((-t256 * t641 - t257 * t638) * pkin(1) + (-pkin(2) - t795) * t1001 + (t256 * (-rSges(4,3) - pkin(7)) + t257 * t852) * t623) * qJD(1) - (-t863 - t256 - t522 + (-t449 - t1063) * qJD(1)) * t257) * m(4) + t1041 + t1007 / 0.2e1 + (t88 + t69) * t1079 + (t324 + ((t185 - t402 + (t442 + t993) * t624 + t941) * t624 + t940 * t623) * qJD(3)) * t844 + (t131 + t130) * t1072 + m(3) * ((-t532 * t645 - t893) * t1115 + (-t892 + (-0.2e1 * t837 - t630 + t1115) * t645) * (-t532 - t1063)) + t1005 / 0.2e1 + t1006 / 0.2e1 + (t78 + t61) * t1083 + (t126 * (-t1116 + t1143) + (-t1052 * t624 + t630 + t821 + t922) * t125 + (rSges(5,1) * t885 + rSges(5,2) * t884 + t920 + (-t610 - t630 + (-t544 - t620) * t624) * qJD(1)) * t174 + (t925 + (-rSges(5,1) * t968 - rSges(5,2) * t964 - t890) * t624 + t174 - t705 - t1119 + (t1063 + t419 + t1143) * qJD(1)) * t175) * m(5) + t61 * t1082 + t69 * t1078 - t1046 / 0.2e1 + (-t899 + (t179 + t1126) * t527 + (t623 * t826 + t178 - t365) * t526 + ((t414 + t750) * t526 + t826 * t527) * t624 + t108) * t1073 + t1008 / 0.2e1 + t1040 + t1045 / 0.2e1 - (t148 + t167 + t112) * t906 / 0.2e1 + (t149 + t166) * t603 / 0.2e1 + (t111 - t898 + ((t624 * t825 + t187 - t940) * t624 + (t623 * t825 + t186 + t827) * t623) * qJD(3)) * t847 + t87 * t1080 + t77 * t1084 + t173 * t1086 + t172 * t1087 + t169 * t1088 + t168 * t1089 + t1043 / 0.2e1 - t1044 / 0.2e1 + (t132 + t109 + t129) * t1071; m(4) * t147 + m(5) * t90 + m(6) * t40 + m(7) * t25; ((-t603 * t986 + t912) * t623 + (t676 + (-t1096 * t624 + (t987 + t694) * t623) * qJD(3)) * t624) * t847 + (-t148 * t624 + t149 * t623 + (t286 * t623 + t287 * t624) * qJD(1)) * t1058 + ((-t906 * t987 - t912) * t624 + (t676 + (t694 * t623 + (-t1096 + t986) * t624) * qJD(3)) * t623) * t844 + t647 + ((t637 * t919 + t640 * t918) * qJD(1) + ((t623 * t930 - t624 * t931) * t640 + (t623 * t932 + t624 * t933) * t637) * qJD(3)) * t1059 + (qJD(1) * t166 + (-(t1099 * t623 + t658 * t624) * t624 + (t1100 * t623 + t657 * t624) * t623 + (t186 * t623 + t187 * t624) * qJD(1)) * t1121) * t1067 + (qJD(1) * t167 + (-(-t1099 * t624 + t658 * t623) * t624 + (-t1100 * t624 + t657 * t623) * t623 + (t184 * t623 + t185 * t624) * qJD(1)) * t1121) * t1066 + (t111 + t707) * t849 + (t112 + t706) * t848 + (t97 * t874 + t25 * (t736 + t943) + t81 * (t692 + t875) + (t1120 * t741 + t720 * t97) * t624 + (t34 * t741 + t98 * t720 + (-t406 + t878) * t1015) * t623 - t653 * t97 - (-t98 * t866 + ((-t623 * t98 - t624 * t97) * t640 + t81 * t836) * qJD(3)) * pkin(3) - t646) * m(7) + (-(-t120 * t866 + (t110 * t836 - t640 * t766) * qJD(3)) * pkin(3) - t649 + t40 * (t812 + t943) + t110 * (t737 + t875) + (t120 * t740 + (-t406 + t949) * t916) * t623 + t1095 * (-t462 + t839) + (t624 * t740 + t1114) * t119) * m(6) + (-(-t175 * t866 + (t150 * t836 + t640 * t757) * qJD(3)) * pkin(3) + t90 * (t942 + t943) + t150 * (t875 + t880) + (t174 * t800 + (qJD(1) * t175 + t126) * t840) * t624 + (t125 * t840 + t175 * t800 + (t174 * t543 + t150 * (-t406 - t420)) * qJD(1)) * t623 + t1103) * m(5) + (-(t256 * t510 - t1000) * qJD(1) - (t241 * (-t510 * t623 - t511 * t624) + t754 * t795) * qJD(3) + t147 * t745 + t241 * t669 + t754 * t552 + (-t180 * t623 - t181 * t624 + (-t257 * t624 + t1002) * qJD(1)) * t581) * m(4); t647 + (t25 * t736 + t81 * t692 + (t1015 * t878 + t879 * t98) * t623 - t646 + (t1120 * t624 + t34 * t623) * (-t545 - t1131) + (t624 * t879 - t653 + t874) * t97) * m(7) + (-t649 + t40 * t812 + t110 * t737 + (t120 * t948 + t916 * t949) * t623 + t1095 * (-t462 - t545) + (t624 * t948 + t1114) * t119) * m(6) + (t90 * t942 + t150 * (-t420 * t911 + t880) + t757 * t501 + (-t125 * t623 - t126 * t624 + (t174 * t623 - t175 * t624) * qJD(1)) * t543 + t1103) * m(5); (t1005 + t1006 + t1040 + t1043 - t1044) * t1065 + t665 + (t697 * t628 + (-t1091 * t636 + t639 * t671) * t626) * t1069 + t13 * t856 + ((t632 * t758 - t101) * t628 + (qJD(1) * t759 + t188 * t632 + t48 * t623 + t49 * t624) * t626) * t1068 + ((t632 * t762 - t87) * t628 + (qJD(1) * t763 + t173 * t632 + t44 * t623 + t45 * t624) * t626) * t1080 + (-t1091 * t493 + t494 * t671 - t623 * t690) * t1078 + t14 * t857 + (t1091 * t495 + t671 * t496 - t624 * t690) * t1081 + ((t632 * t764 - t88) * t628 + (qJD(1) * t765 + t172 * t632 + t46 * t623 + t47 * t624) * t626) * t1079 + (-t172 * t628 + t626 * t764) * t1087 + (-t173 * t628 + t626 * t762) * t1086 + (qJD(5) * (-t188 * t628 + t626 * t758) + t80) * t855 + t1124 * t69 + t1125 * t68 + (-t97 * (t448 * t894 + t471 * t592 + t829) - t98 * (t447 * t894 + t472 * t592 + t834) - t81 * (-t447 * t471 + t448 * t472 + t835) + t97 * t813 + t98 * t957 + t81 * t959 + (-t35 * t284 + t97 * t165 - t34 * t1132 - t98 * t1134 + ((-t1131 * t98 - t81 * t284) * t624 + t721 * t623) * t632) * t628 + ((-t1133 * t97 + t98 * t285) * t632 + (qJD(1) * t721 - t1131 * t34 + t81 * t165 - t25 * t284 + t956 * t98) * t624 + (-t35 * t700 + t97 * t311 - t25 * t1132 - t81 * t1134 + (-t1133 * t81 - t700 * t98) * qJD(1)) * t623) * t626 + t722) * m(7) + (-t119 * (-t349 * t592 - t448 * t515) - t120 * (t350 * t592 - t447 * t515) - t110 * (t349 * t447 + t350 * t448) + (t119 * t163 - t120 * t162 + t75 * t300 - t74 * t302 + (t110 * t751 + (t119 * t623 - t120 * t624) * t462) * t632) * t628 + (t119 * (-t300 * t632 + t310 * t623) + t120 * (t302 * t632 - t310 * t624) + t40 * t751 + t110 * (-t162 * t623 + t163 * t624 - t300 * t911 - t302 * t910) + (qJD(1) * t766 + t75 * t623 - t74 * t624) * t462) * t626) * m(6); t665 + (t34 * (-t260 * t628 - t428 * t976) - t25 * t260 * t984 + t722 + (-t145 * t628 + (-t240 * t626 - t428 * t964) * t624 + t957 - t834) * t98 + (-t258 * t968 + t813 - t829) * t97 + (-t260 * t867 + t959 + (-t260 * t964 + (-qJD(1) * t258 - t145) * t626) * t623 - t835) * t81) * m(7);];
tauc  = t1(:);
