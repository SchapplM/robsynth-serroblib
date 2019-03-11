% Calculate vector of inverse dynamics joint torques for
% S6RPRPPR5
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% qJDD [6x1]
%   Generalized joint accelerations
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d3,d6,theta2,theta5]';
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
% tau [6x1]
%   joint torques of inverse dynamics (contains inertial, gravitational coriolis and centrifugal forces)

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 02:52
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tau = S6RPRPPR5_invdynJ_fixb_slag_vp1(qJ, qJD, qJDD, g, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(3,1),zeros(10,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRPPR5_invdynJ_fixb_slag_vp1: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RPRPPR5_invdynJ_fixb_slag_vp1: qJD has to be [6x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [6 1]), ...
  'S6RPRPPR5_invdynJ_fixb_slag_vp1: qJDD has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RPRPPR5_invdynJ_fixb_slag_vp1: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RPRPPR5_invdynJ_fixb_slag_vp1: pkin has to be [10x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RPRPPR5_invdynJ_fixb_slag_vp1: m has to be [7x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [7,3]), ...
  'S6RPRPPR5_invdynJ_fixb_slag_vp1: rSges has to be [7x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [7 6]), ...
  'S6RPRPPR5_invdynJ_fixb_slag_vp1: Icges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 02:50:05
% EndTime: 2019-03-09 02:51:36
% DurationCPUTime: 83.58s
% Computational Cost: add. (31346->1332), mult. (37305->1637), div. (0->0), fcn. (33228->10), ass. (0->627)
t1094 = Icges(5,4) - Icges(4,5);
t1093 = Icges(5,5) - Icges(4,6);
t1092 = Icges(5,1) + Icges(4,3);
t571 = pkin(9) + qJ(3);
t551 = sin(t571);
t553 = cos(t571);
t1061 = t1093 * t551 - t1094 * t553;
t930 = Icges(4,4) * t551;
t425 = Icges(4,2) * t553 + t930;
t918 = Icges(5,6) * t551;
t681 = Icges(5,3) * t553 + t918;
t1091 = t425 + t681;
t533 = Icges(4,4) * t553;
t427 = Icges(4,1) * t551 + t533;
t917 = Icges(5,6) * t553;
t683 = Icges(5,2) * t551 + t917;
t1090 = t427 + t683;
t581 = cos(qJ(1));
t1089 = t1092 * t581;
t580 = sin(qJ(1));
t892 = t553 * t580;
t896 = t551 * t580;
t1064 = -t1093 * t896 + t1094 * t892 + t1089;
t1075 = t1061 * t581 + t1092 * t580;
t919 = Icges(4,6) * t581;
t311 = Icges(4,4) * t892 - Icges(4,2) * t896 - t919;
t499 = Icges(5,6) * t896;
t928 = Icges(5,4) * t581;
t318 = Icges(5,2) * t892 - t499 + t928;
t1088 = t311 * t551 - t318 * t553;
t428 = Icges(4,1) * t553 - t930;
t314 = Icges(4,5) * t580 + t428 * t581;
t682 = -Icges(5,3) * t551 + t917;
t315 = Icges(5,5) * t580 - t581 * t682;
t1087 = -t314 * t892 - t315 * t896;
t923 = Icges(5,5) * t581;
t316 = Icges(5,6) * t892 - Icges(5,3) * t896 + t923;
t1086 = t311 + t316;
t691 = -Icges(4,2) * t551 + t533;
t312 = Icges(4,6) * t580 + t581 * t691;
t1085 = t312 - t315;
t504 = Icges(4,4) * t896;
t924 = Icges(4,5) * t581;
t313 = Icges(4,1) * t892 - t504 - t924;
t1084 = t313 + t318;
t895 = t551 * t581;
t500 = Icges(5,6) * t895;
t891 = t553 * t581;
t929 = Icges(5,4) * t580;
t317 = -Icges(5,2) * t891 + t500 + t929;
t1083 = -t314 + t317;
t421 = Icges(4,5) * t551 + Icges(4,6) * t553;
t689 = Icges(5,4) * t551 + Icges(5,5) * t553;
t1082 = t421 - t689;
t684 = Icges(5,2) * t553 - t918;
t1079 = t428 + t684;
t1078 = t1091 * qJD(3);
t1077 = t1090 * qJD(3);
t1057 = -t313 * t553 + t316 * t551 + t1088;
t1076 = t682 + t691;
t1021 = -t1075 * t580 - t314 * t891 - t315 * t895;
t1030 = -t312 * t895 - t317 * t891 - t1021;
t576 = cos(pkin(10));
t887 = t576 * t581;
t574 = sin(pkin(10));
t890 = t574 * t580;
t395 = t551 * t887 - t890;
t889 = t574 * t581;
t792 = t551 * t889;
t888 = t576 * t580;
t396 = t792 + t888;
t186 = Icges(6,5) * t396 + Icges(6,6) * t395 + Icges(6,3) * t891;
t189 = Icges(6,4) * t396 + Icges(6,2) * t395 + Icges(6,6) * t891;
t192 = Icges(6,1) * t396 + Icges(6,4) * t395 + Icges(6,5) * t891;
t65 = t186 * t891 + t395 * t189 + t396 * t192;
t1009 = t65 + t1030;
t657 = t425 * t551 - t427 * t553;
t1059 = -t551 * t681 + t553 * t683 - t657;
t1074 = -t1075 * t581 - t1087;
t1031 = -t312 * t896 - t317 * t892 + t1074;
t397 = t551 * t888 + t889;
t789 = t551 * t890;
t398 = -t789 + t887;
t67 = t186 * t892 + t397 * t189 - t398 * t192;
t1011 = t67 + t1031;
t1073 = t1064 * t580 - t313 * t891 + t316 * t895;
t1072 = t312 * t551 + t317 * t553;
t1071 = t311 * t895 - t318 * t891 + t1073;
t1070 = t1084 * t551 + t1086 * t553;
t1069 = t1083 * t551 - t1085 * t553;
t1068 = t1078 * t581 + (t1076 * t580 - t919 + t923) * qJD(1);
t1067 = qJD(1) * t1085 - t1078 * t580;
t1066 = -t1077 * t581 + (-t1079 * t580 + t924 - t928) * qJD(1);
t1065 = t1077 * t580 + (-t581 * t684 - t314 + t929) * qJD(1);
t1063 = t1076 * qJD(3);
t1062 = t1079 * qJD(3);
t1060 = t1082 * qJD(3);
t1058 = -t1090 * t551 - t1091 * t553;
t1056 = -t314 * t553 - t315 * t551 + t1072;
t1019 = t1082 * t580;
t686 = Icges(6,5) * t574 + Icges(6,6) * t576;
t293 = Icges(6,3) * t551 - t553 * t686;
t688 = Icges(6,4) * t574 + Icges(6,2) * t576;
t610 = -Icges(6,6) * t551 + t553 * t688;
t693 = Icges(6,1) * t574 + Icges(6,4) * t576;
t612 = -Icges(6,5) * t551 + t553 * t693;
t898 = t421 * t581;
t1036 = t293 * t892 - t397 * t610 + t398 * t612 - t580 * t657 - t898;
t373 = t689 * t581;
t152 = t681 * t896 - t683 * t892 - t373;
t1055 = t152 - t1036;
t188 = -Icges(6,5) * t398 + Icges(6,6) * t397 + Icges(6,3) * t892;
t191 = -Icges(6,4) * t398 + Icges(6,2) * t397 + Icges(6,6) * t892;
t194 = -Icges(6,1) * t398 + Icges(6,4) * t397 + Icges(6,5) * t892;
t671 = t191 * t397 - t194 * t398;
t1012 = -t1057 * t580 + t1064 * t581 + t188 * t892 + t671;
t1054 = t1075 * qJD(1);
t1008 = t1059 * t581 + t293 * t891 - t395 * t610 - t396 * t612 + t1019;
t570 = pkin(10) + qJ(6);
t550 = sin(t570);
t552 = cos(t570);
t894 = t552 * t580;
t349 = t550 * t581 + t551 * t894;
t893 = t552 * t581;
t350 = -t550 * t896 + t893;
t861 = t350 * rSges(7,1) - t349 * rSges(7,2);
t179 = rSges(7,3) * t892 - t861;
t715 = rSges(7,1) * t550 + rSges(7,2) * t552;
t1023 = t553 * t715;
t283 = rSges(7,3) * t551 - t1023;
t813 = qJD(6) * t553;
t820 = qJD(3) * t581;
t414 = -t580 * t813 + t820;
t814 = qJD(6) * t551;
t520 = qJD(1) + t814;
t960 = pkin(5) * t574;
t801 = t553 * t960;
t578 = -pkin(8) - qJ(5);
t884 = -qJ(5) - t578;
t305 = t551 * t884 - t801;
t433 = pkin(3) * t551 - qJ(4) * t553;
t911 = qJ(5) * t551;
t754 = -t433 - t911;
t724 = -t305 + t754;
t1053 = t179 * t520 + t283 * t414 - t724 * t820;
t1052 = qJD(1) * t1082 + t1058 * qJD(3) + t1062 * t553 - t1063 * t551;
t1051 = qJD(1) * t1064 + qJD(3) * t1070 + t1065 * t553 + t1067 * t551;
t1050 = qJD(3) * t1069 + t1066 * t553 + t1068 * t551 + t1054;
t1049 = -t1085 * t580 + t1086 * t581;
t1048 = -t1076 - t1090;
t1047 = t1079 - t1091;
t1046 = (t499 + t504 + (Icges(4,2) + Icges(5,3)) * t892 - t1084) * t581 + (-Icges(5,3) * t891 - t425 * t581 - t1083 - t500) * t580;
t1045 = qJD(1) * t1059 - qJD(3) * t1061;
t1044 = -t1060 * t581 + (-t1061 * t580 + t1056 + t1089) * qJD(1);
t1043 = qJD(1) * t1057 - t1060 * t580 + t1054;
t1042 = t1011 * t580 - t1012 * t581;
t999 = -t395 * t191 - t396 * t194;
t66 = t188 * t891 - t999;
t937 = t581 * t66;
t1041 = t1009 * t580 + t1071 * t581 - t937;
t416 = t581 * pkin(4) - qJ(5) * t892;
t577 = cos(pkin(9));
t539 = pkin(2) * t577 + pkin(1);
t579 = -pkin(7) - qJ(2);
t544 = t581 * t579;
t837 = -t580 * t539 - t544;
t857 = t398 * rSges(6,1) - t397 * rSges(6,2);
t1040 = t416 + t837 + t857;
t538 = pkin(5) * t576 + pkin(4);
t886 = t578 * t580;
t841 = t581 * t538 + t553 * t886;
t1039 = t837 + t841 + t861;
t209 = rSges(6,3) * t892 - t857;
t535 = t553 * rSges(7,3);
t282 = t551 * t715 + t535;
t1038 = t1008 * qJD(1);
t1037 = t1055 * qJD(1);
t821 = qJD(3) * t580;
t413 = t581 * t813 + t821;
t347 = -t550 * t580 + t551 * t893;
t348 = t550 * t895 + t894;
t164 = Icges(7,5) * t348 + Icges(7,6) * t347 + Icges(7,3) * t891;
t927 = Icges(7,4) * t348;
t167 = Icges(7,2) * t347 + Icges(7,6) * t891 + t927;
t331 = Icges(7,4) * t347;
t170 = Icges(7,1) * t348 + Icges(7,5) * t891 + t331;
t52 = t164 * t891 + t347 * t167 + t348 * t170;
t166 = -Icges(7,5) * t350 + Icges(7,6) * t349 + Icges(7,3) * t892;
t333 = Icges(7,4) * t350;
t169 = Icges(7,2) * t349 + Icges(7,6) * t892 - t333;
t332 = Icges(7,4) * t349;
t171 = Icges(7,1) * t350 - Icges(7,5) * t892 - t332;
t53 = t166 * t891 + t347 * t169 - t171 * t348;
t685 = Icges(7,5) * t550 + Icges(7,6) * t552;
t607 = -Icges(7,3) * t551 + t553 * t685;
t926 = Icges(7,4) * t550;
t687 = Icges(7,2) * t552 + t926;
t609 = -Icges(7,6) * t551 + t553 * t687;
t925 = Icges(7,4) * t552;
t692 = Icges(7,1) * t550 + t925;
t611 = -Icges(7,5) * t551 + t553 * t692;
t88 = -t347 * t609 - t348 * t611 - t607 * t891;
t13 = t413 * t52 - t414 * t53 + t88 * t520;
t54 = t164 * t892 + t349 * t167 - t350 * t170;
t55 = t166 * t892 + t169 * t349 + t171 * t350;
t89 = -t349 * t609 + t350 * t611 - t607 * t892;
t14 = t413 * t54 - t414 * t55 + t520 * t89;
t674 = t169 * t552 - t171 * t550;
t64 = t166 * t551 - t553 * t674;
t532 = t551 * qJ(4);
t962 = -rSges(6,3) - pkin(3);
t1029 = t553 * t962 - t532;
t568 = t580 * pkin(4);
t415 = qJ(5) * t891 + t568;
t773 = t551 * t821;
t468 = qJ(5) * t773;
t816 = qJD(5) * t580;
t245 = qJD(1) * t415 + t553 * t816 - t468;
t815 = qJD(5) * t581;
t484 = t553 * t815;
t824 = qJD(1) * t581;
t549 = pkin(4) * t824;
t772 = t551 * t820;
t825 = qJD(1) * t580;
t775 = t553 * t825;
t631 = t772 + t775;
t990 = -qJ(5) * t631 + t549;
t246 = t484 + t990;
t530 = qJD(5) * t553;
t494 = qJ(4) * t892;
t378 = -pkin(3) * t896 + t494;
t497 = qJ(4) * t891;
t382 = -pkin(3) * t895 + t497;
t531 = qJD(4) * t551;
t780 = t378 * t821 + t382 * t820 + t531;
t767 = t553 * t820;
t469 = qJ(4) * t767;
t817 = qJD(4) * t581;
t486 = t551 * t817;
t776 = t551 * t825;
t184 = -pkin(3) * t631 - qJ(4) * t776 + t469 + t486;
t768 = t553 * t821;
t357 = t551 * t824 + t768;
t476 = pkin(3) * t773;
t818 = qJD(4) * t580;
t769 = t551 * t818;
t774 = t553 * t824;
t185 = pkin(3) * t774 + qJ(4) * t357 - t476 + t769;
t537 = t553 * pkin(3);
t994 = t537 + t532;
t381 = t994 * t580;
t787 = t581 * t184 + t580 * t185 + t381 * t824;
t1028 = t580 * t245 + t581 * t246 - t416 * t824 - t530 - t780 + t787;
t559 = t581 * qJ(2);
t488 = pkin(1) * t580 - t559;
t307 = t488 + t837;
t459 = qJD(1) * t488;
t995 = qJD(1) * t307 - t459;
t991 = -qJD(1) * t381 + t995;
t1027 = -qJD(1) * t416 + t469 - t991;
t806 = qJD(3) * qJD(5);
t1026 = qJDD(5) * t553 - 0.2e1 * t551 * t806;
t1025 = t1064 + t1072;
t807 = qJD(3) * qJD(4);
t1024 = qJDD(4) * t551 + t553 * t807;
t1022 = t1061 * qJD(1);
t1020 = t898 - t373;
t1018 = qJD(3) * t1041 + t1038;
t1017 = qJD(3) * t1042 - t1037;
t259 = -qJD(1) * t397 + t576 * t767;
t729 = t574 * t767;
t260 = qJD(1) * t398 + t729;
t103 = Icges(6,5) * t260 + Icges(6,6) * t259 - Icges(6,3) * t631;
t105 = Icges(6,4) * t260 + Icges(6,2) * t259 - Icges(6,6) * t631;
t107 = Icges(6,1) * t260 + Icges(6,4) * t259 - Icges(6,5) * t631;
t672 = t189 * t576 + t192 * t574;
t1016 = (-t105 * t576 - t107 * t574 - t1068) * t553 + (t103 + t1066) * t551 + (t186 * t553 + t551 * t672 - t1056) * qJD(3);
t257 = qJD(1) * t395 + t576 * t768;
t258 = qJD(1) * t396 + t574 * t768;
t356 = -t773 + t774;
t102 = Icges(6,5) * t258 + Icges(6,6) * t257 + Icges(6,3) * t356;
t104 = Icges(6,4) * t258 + Icges(6,2) * t257 + Icges(6,6) * t356;
t106 = Icges(6,1) * t258 + Icges(6,4) * t257 + Icges(6,5) * t356;
t670 = t191 * t576 + t194 * t574;
t1015 = (t104 * t576 + t106 * t574 - t1067) * t553 + (-t102 + t1065) * t551 + (-t553 * t188 - t551 * t670 + t1057) * qJD(3);
t292 = Icges(6,3) * t553 + t551 * t686;
t271 = t292 * qJD(3);
t294 = Icges(6,6) * t553 + t551 * t688;
t272 = t294 * qJD(3);
t296 = Icges(6,5) * t553 + t551 * t693;
t273 = t296 * qJD(3);
t1014 = t1045 * t581 + t1052 * t580 - t257 * t610 - t258 * t612 + t271 * t892 + t272 * t397 - t273 * t398 + t293 * t356;
t1013 = -t1045 * t580 + t1052 * t581 - t259 * t610 - t260 * t612 + t271 * t891 + t272 * t395 + t273 * t396 - t293 * t631;
t1010 = t66 - t1071;
t1007 = t188 * t551 - t553 * t670 + t1070;
t1006 = t186 * t551 - t553 * t672 - t1069;
t667 = -t574 * t612 - t576 * t610;
t586 = (t292 + t667) * qJD(1) + (t580 * t672 - t581 * t670) * qJD(3);
t826 = qJD(1) * t551;
t1005 = -t293 * t826 + t553 * t586 + (t1047 * t553 + t1048 * t551) * qJD(1);
t1004 = t1049 * t553 + (-t186 * t580 + t188 * t581 - t1046) * t551;
t1003 = (-t102 * t891 - t104 * t395 + t1051 * t581 - t106 * t396 + t188 * t631 - t191 * t259 - t194 * t260) * t581 + (t103 * t891 + t105 * t395 + t107 * t396 - t186 * t631 + t189 * t259 + t192 * t260 + t1044 * t580 + (-t1043 + t1050) * t581) * t580;
t1002 = (-t102 * t892 - t104 * t397 + t1043 * t581 + t106 * t398 - t188 * t356 - t191 * t257 - t194 * t258) * t581 + (t103 * t892 + t105 * t397 - t107 * t398 + t186 * t356 + t189 * t257 + t192 * t258 + t1050 * t580 + (-t1044 + t1051) * t581) * t580;
t335 = (-rSges(7,1) * t552 + rSges(7,2) * t550) * t553;
t156 = qJD(3) * t282 + qJD(6) * t335;
t897 = t551 * t574;
t492 = pkin(5) * t897;
t632 = t884 * t553 + t492;
t284 = t632 * qJD(3);
t819 = qJD(4) * t553;
t334 = qJD(3) * t994 - t819;
t529 = qJD(5) * t551;
t822 = qJD(3) * t553;
t630 = -qJ(5) * t822 - t334 - t529;
t910 = qJ(5) * t553;
t753 = -t994 - t910;
t1001 = -qJD(3) * (t753 - t632) - t156 - t284 + t630;
t393 = t433 * t825;
t1000 = t551 * t815 + t393;
t727 = g(1) * t581 + g(2) * t580;
t997 = t551 * t727;
t534 = t551 * rSges(5,3);
t945 = rSges(5,2) * t553;
t714 = t534 - t945;
t558 = t580 * qJ(2);
t490 = t581 * pkin(1) + t558;
t575 = sin(pkin(9));
t947 = rSges(3,2) * t575;
t950 = rSges(3,1) * t577;
t365 = t580 * rSges(3,3) + (-t947 + t950) * t581;
t322 = (Icges(7,2) * t550 - t925) * t553;
t600 = t413 * (-Icges(7,2) * t348 + t170 + t331) - t414 * (Icges(7,2) * t350 - t171 + t332) + t520 * (-t611 + t322);
t323 = (-Icges(7,1) * t552 + t926) * t553;
t601 = t413 * (-Icges(7,1) * t347 + t167 + t927) - t414 * (-Icges(7,1) * t349 + t169 - t333) + t520 * (-t609 - t323);
t981 = 0.2e1 * t551;
t980 = 0.2e1 * t553;
t979 = m(5) / 0.2e1;
t978 = m(6) / 0.2e1;
t977 = m(7) / 0.2e1;
t976 = -m(6) - m(7);
t808 = qJD(1) * qJD(3);
t452 = qJDD(3) * t580 + t581 * t808;
t803 = qJDD(6) * t553;
t224 = -qJD(6) * t631 + t581 * t803 + t452;
t975 = t224 / 0.2e1;
t453 = -qJDD(3) * t581 + t580 * t808;
t225 = qJD(6) * t356 + t580 * t803 + t453;
t974 = t225 / 0.2e1;
t399 = qJD(3) * t813 + qJDD(6) * t551 + qJDD(1);
t973 = t399 / 0.2e1;
t972 = -t413 / 0.2e1;
t971 = t413 / 0.2e1;
t970 = -t414 / 0.2e1;
t969 = t414 / 0.2e1;
t968 = t452 / 0.2e1;
t967 = t453 / 0.2e1;
t966 = -t520 / 0.2e1;
t965 = t520 / 0.2e1;
t964 = t580 / 0.2e1;
t963 = -t581 / 0.2e1;
t961 = -rSges(7,3) - pkin(3);
t959 = g(1) * t580;
t120 = t468 + (t551 * t578 + t801) * t821 + ((-pkin(4) + t538) * t580 + t632 * t581) * qJD(1);
t653 = t530 + t531;
t557 = qJD(2) * t581;
t408 = qJD(1) * t490 - t557;
t525 = t579 * t825;
t952 = pkin(1) - t539;
t870 = t525 - (-t581 * t952 - t558) * qJD(1) - t408;
t786 = -t185 + t870;
t603 = -t580 * t653 - t245 + t786;
t809 = qJD(1) * qJD(2);
t833 = qJDD(2) * t580 + t581 * t809;
t697 = t1024 * t581 + t453 * t433 + t833;
t622 = t1026 * t581 + t453 * t911 + t697;
t788 = qJD(3) ^ 2 * t910;
t625 = -t788 + (-t284 - t334) * qJD(3);
t745 = pkin(5) * t789;
t223 = t416 + t745 - t841;
t869 = t307 - t488;
t782 = -t381 + t869;
t740 = t416 + t782;
t696 = -t223 + t740;
t742 = qJD(6) + t826;
t614 = t581 * t742 + t768;
t655 = t550 * t520;
t145 = t552 * t614 - t580 * t655;
t654 = t552 * t520;
t146 = t550 * t614 + t580 * t654;
t717 = rSges(7,1) * t146 + rSges(7,2) * t145;
t86 = rSges(7,3) * t356 + t717;
t4 = t625 * t581 + t696 * qJDD(1) + t453 * t305 + (-t120 + t603) * qJD(1) - t520 * t86 + t622 - t414 * t156 - t399 * t179 + t225 * t283;
t956 = t4 * t580;
t735 = pkin(5) * t729 + t538 * t824 + t578 * t631;
t121 = -qJD(1) * t745 + t735 - t990;
t177 = t348 * rSges(7,1) + t347 * rSges(7,2) + rSges(7,3) * t891;
t885 = t578 * t581;
t652 = pkin(5) * t792 + t580 * t538 - t553 * t885;
t222 = t652 - t415;
t495 = qJ(4) * t895;
t385 = pkin(3) * t891 + t495;
t508 = t581 * t539;
t748 = -t579 * t580 + t508;
t308 = t748 - t490;
t545 = qJ(2) * t824;
t556 = qJD(2) * t580;
t832 = t545 + t556;
t648 = -qJDD(2) * t581 + qJD(1) * (-pkin(1) * t825 + t832) + qJDD(1) * t490 + t580 * t809;
t623 = qJD(1) * (-t545 + (t580 * t952 - t544) * qJD(1)) + qJDD(1) * t308 + t648;
t599 = qJDD(1) * t385 + t623 + t1024 * t580 + (t184 + t486) * qJD(1);
t591 = qJDD(1) * t415 + t599 + t1026 * t580 + (t246 + t484) * qJD(1);
t613 = -t580 * t742 + t767;
t147 = t552 * t613 - t581 * t655;
t148 = t550 * t613 + t581 * t654;
t883 = t148 * rSges(7,1) + t147 * rSges(7,2);
t87 = -rSges(7,3) * t631 + t883;
t5 = qJD(1) * t121 + qJDD(1) * t222 - t413 * t156 + t399 * t177 - t224 * t283 + t452 * t724 + t520 * t87 + t580 * t625 + t591;
t955 = t5 * t581;
t274 = Icges(7,3) * t553 + t551 * t685;
t321 = (-Icges(7,5) * t552 + Icges(7,6) * t550) * t553;
t153 = qJD(3) * t274 + qJD(6) * t321;
t276 = Icges(7,6) * t553 + t551 * t687;
t154 = qJD(3) * t276 + qJD(6) * t322;
t278 = Icges(7,5) * t553 + t551 * t692;
t155 = qJD(3) * t278 + qJD(6) * t323;
t668 = -t550 * t611 - t552 * t609;
t29 = (qJD(3) * t668 + t153) * t551 + (-qJD(3) * t607 - t154 * t552 - t155 * t550 + (-t550 * t609 + t552 * t611) * qJD(6)) * t553;
t99 = -t551 * t607 - t553 * t668;
t951 = t29 * t520 + t99 * t399;
t949 = rSges(4,1) * t553;
t946 = rSges(5,2) * t551;
t944 = rSges(6,2) * t576;
t675 = t167 * t552 + t170 * t550;
t81 = Icges(7,5) * t148 + Icges(7,6) * t147 - Icges(7,3) * t631;
t83 = Icges(7,4) * t148 + Icges(7,2) * t147 - Icges(7,6) * t631;
t85 = Icges(7,1) * t148 + Icges(7,4) * t147 - Icges(7,5) * t631;
t10 = (qJD(3) * t675 + t81) * t551 + (qJD(3) * t164 - t550 * t85 - t552 * t83 + (t167 * t550 - t170 * t552) * qJD(6)) * t553;
t942 = t10 * t413;
t80 = Icges(7,5) * t146 + Icges(7,6) * t145 + Icges(7,3) * t356;
t82 = Icges(7,4) * t146 + Icges(7,2) * t145 + Icges(7,6) * t356;
t84 = Icges(7,1) * t146 + Icges(7,4) * t145 + Icges(7,5) * t356;
t11 = (qJD(3) * t674 + t80) * t551 + (qJD(3) * t166 - t550 * t84 - t552 * t82 + (t169 * t550 + t171 * t552) * qJD(6)) * t553;
t941 = t11 * t414;
t536 = t553 * rSges(6,3);
t565 = t580 * rSges(5,1);
t563 = t580 * rSges(4,3);
t63 = t164 * t551 - t553 * t675;
t936 = t63 * t224;
t935 = t64 * t225;
t933 = rSges(7,3) - t578;
t435 = rSges(4,1) * t551 + rSges(4,2) * t553;
t722 = -t435 * t820 + t556;
t340 = rSges(4,1) * t892 - rSges(4,2) * t896 - t581 * rSges(4,3);
t783 = -t340 + t869;
t124 = qJD(1) * t783 + t722;
t908 = t124 * t580;
t341 = rSges(4,1) * t891 - rSges(4,2) * t895 + t563;
t868 = t308 + t490;
t125 = -t435 * t821 - t557 + (t341 + t868) * qJD(1);
t384 = t435 * t581;
t907 = t125 * t384;
t878 = t179 + t223;
t877 = t260 * rSges(6,1) + t259 * rSges(6,2);
t862 = -t714 * qJD(3) - t334;
t342 = -rSges(5,2) * t891 + rSges(5,3) * t895 + t565;
t860 = -t342 - t385;
t859 = t580 * t381 + t581 * t385;
t858 = qJD(1) * t382 + t553 * t818;
t286 = t365 + t490;
t856 = -t385 - t415;
t855 = -t433 * t821 - t557;
t849 = t1023 * t580;
t848 = t1023 * t581;
t713 = rSges(5,3) * t553 + t946;
t847 = -t433 + t713;
t846 = -t994 - t714;
t791 = t553 * t890;
t797 = t553 * t944;
t845 = rSges(6,1) * t791 + t580 * t797;
t790 = t553 * t889;
t844 = rSges(6,1) * t790 + t581 * t797;
t843 = pkin(5) * t791 + t551 * t886;
t842 = pkin(5) * t790 + t551 * t885;
t840 = rSges(4,2) * t776 + rSges(4,3) * t824;
t839 = t486 + t556;
t800 = t580 * t950;
t526 = t580 * t947;
t834 = t581 * rSges(3,3) + t526;
t364 = t800 - t834;
t838 = -t488 - t364;
t379 = rSges(5,2) * t896 + rSges(5,3) * t892;
t383 = rSges(5,2) * t895 + rSges(5,3) * t891;
t836 = rSges(3,3) * t824 + qJD(1) * t526;
t835 = t525 + t557;
t831 = t580 ^ 2 + t581 ^ 2;
t823 = qJD(3) * t551;
t810 = -m(5) + t976;
t802 = -qJ(5) + t962;
t793 = qJ(5) * t895;
t207 = t396 * rSges(6,1) + t395 * rSges(6,2) + rSges(6,3) * t891;
t785 = -t207 + t856;
t784 = -t222 + t856;
t781 = t385 + t868;
t779 = -t468 + t855;
t778 = t476 + t835;
t777 = t484 + t839;
t766 = -pkin(1) - t950;
t763 = t824 / 0.2e1;
t762 = -t821 / 0.2e1;
t761 = t821 / 0.2e1;
t760 = -t820 / 0.2e1;
t759 = t820 / 0.2e1;
t757 = -qJ(4) - t960;
t756 = rSges(5,1) * t581 - rSges(5,3) * t896;
t747 = qJD(3) * t862;
t746 = qJD(3) * t831;
t487 = t553 * t817;
t744 = -qJD(1) * t378 + t487;
t343 = rSges(5,2) * t892 + t756;
t741 = t343 + t782;
t739 = t415 + t781;
t736 = t581 * t415 - t580 * t416 + t859;
t732 = rSges(5,1) * t824 + rSges(5,2) * t631 + rSges(5,3) * t767;
t731 = t580 * t802;
t730 = t802 * t581;
t728 = t831 * t911;
t718 = rSges(6,1) * t574 + t944;
t301 = rSges(6,3) * t551 - t553 * t718;
t725 = -t301 + t754;
t721 = t381 * t821 + t385 * t820 - t819;
t491 = rSges(2,1) * t581 - rSges(2,2) * t580;
t489 = rSges(2,1) * t580 + rSges(2,2) * t581;
t438 = -rSges(4,2) * t551 + t949;
t720 = rSges(6,1) * t258 + rSges(6,2) * t257;
t705 = t52 * t581 + t53 * t580;
t704 = t52 * t580 - t53 * t581;
t703 = t54 * t581 + t55 * t580;
t702 = t54 * t580 - t55 * t581;
t701 = t580 * t64 + t581 * t63;
t700 = t580 * t63 - t581 * t64;
t695 = -t209 + t740;
t676 = -t124 * t581 - t125 * t580;
t673 = t177 * t580 - t179 * t581;
t210 = -rSges(4,1) * t631 - rSges(4,2) * t767 + t840;
t380 = t435 * t580;
t211 = -qJD(3) * t380 + (t438 * t581 + t563) * qJD(1);
t669 = t210 * t581 + t211 * t580;
t661 = t340 * t580 + t341 * t581;
t651 = -t283 + t724;
t650 = t748 + t385;
t649 = -t492 - t532;
t61 = (-qJD(3) * t301 + t653) * t580 + (t207 + t739) * qJD(1) + t779;
t647 = t61 * t725;
t636 = -qJDD(4) * t553 + t184 * t820 + t185 * t821 + t452 * t381 + t551 * t807;
t633 = t820 * t847 + t839;
t629 = -t164 * t413 + t166 * t414 + t520 * t607;
t628 = (Icges(7,5) * t347 - Icges(7,6) * t348) * t413 - (Icges(7,5) * t349 + Icges(7,6) * t350) * t414 + t321 * t520;
t627 = t415 * t820 - t416 * t821 + t529 + t721;
t280 = (t551 * t718 + t536) * qJD(3);
t626 = -t788 + (-t280 - t334) * qJD(3);
t602 = qJDD(5) * t551 + t245 * t821 + t246 * t820 - t452 * t416 + t553 * t806 + t636;
t3 = t602 + t784 * t453 + (t120 * t580 + t121 * t581) * qJD(3) - t177 * t225 + t179 * t224 + t223 * t452 + t413 * t86 + t414 * t87;
t46 = qJD(1) * t696 - t1053 + t777;
t47 = t177 * t520 - t283 * t413 + (-qJD(3) * t305 + t653) * t580 + (t222 + t739) * qJD(1) + t779;
t624 = t46 * t820 + t47 * t821 - t3;
t118 = rSges(6,3) * t356 + t720;
t119 = -rSges(6,3) * t631 + t877;
t12 = t209 * t452 + (t118 * t580 + t119 * t581) * qJD(3) + t785 * t453 + t602;
t60 = qJD(1) * t695 + t725 * t820 + t777;
t621 = t60 * t820 + t61 * t821 - t12;
t618 = -t539 - t994 - t534;
t617 = -t280 + t630;
t616 = t553 * t628;
t604 = -t535 + t649 - t537;
t33 = t177 * t414 + t179 * t413 + (t222 * t581 + t223 * t580) * qJD(3) + t627;
t598 = t33 * t673 + (-t46 * t580 + t47 * t581) * t283;
t588 = (t607 * t581 + t675) * t413 - (t607 * t580 + t674) * t414 + (t274 + t668) * t520;
t587 = t588 * t553;
t19 = qJD(1) * t119 + qJDD(1) * t207 + t452 * t725 + t580 * t626 + t591;
t20 = t453 * t301 + t626 * t581 + t695 * qJDD(1) + (-t118 + t603) * qJD(1) + t622;
t62 = (t207 * t581 + t209 * t580) * qJD(3) + t627;
t584 = (qJD(3) * t62 + t19 * t580 + t20 * t581 - t60 * t825 + t61 * t824) * t978 + (qJD(3) * t33 + t4 * t581 - t46 * t825 + t47 * t824 + t5 * t580) * t977;
t410 = t438 * qJD(3);
t359 = t767 - t776;
t355 = t553 * t746;
t354 = t551 * t746;
t262 = t793 + t842;
t261 = qJ(5) * t896 + t843;
t256 = -rSges(6,3) * t895 + t844;
t255 = -rSges(6,3) * t896 + t845;
t254 = t612 * t581;
t253 = t612 * t580;
t252 = t610 * t581;
t251 = t610 * t580;
t244 = -rSges(7,3) * t895 + t848;
t243 = -rSges(7,3) * t896 + t849;
t242 = qJD(1) * t286 - t557;
t241 = qJD(1) * t838 + t556;
t240 = t611 * t581;
t239 = t611 * t580;
t238 = t609 * t581;
t237 = t609 * t580;
t221 = rSges(7,1) * t349 + rSges(7,2) * t350;
t220 = rSges(7,1) * t347 - rSges(7,2) * t348;
t213 = -rSges(5,3) * t776 + t732;
t212 = t713 * t821 + (t581 * t714 + t565) * qJD(1);
t180 = t661 * qJD(3);
t123 = qJDD(1) * t365 + qJD(1) * (-qJD(1) * t800 + t836) + t648;
t122 = t838 * qJDD(1) + (-qJD(1) * t365 - t408) * qJD(1) + t833;
t100 = (t342 * t581 - t343 * t580) * qJD(3) + t721;
t94 = (qJD(3) * t713 + t531) * t580 + (t342 + t781) * qJD(1) + t855;
t93 = qJD(1) * t741 + t633;
t57 = qJD(1) * t210 + qJDD(1) * t341 - t410 * t821 - t435 * t452 + t623;
t56 = -t410 * t820 + t435 * t453 + t783 * qJDD(1) + (-t211 + t870) * qJD(1) + t833;
t32 = -t343 * t452 + t860 * t453 + (t212 * t580 + t213 * t581) * qJD(3) + t636;
t31 = qJD(1) * t213 + qJDD(1) * t342 + t452 * t847 + t580 * t747 + t599;
t30 = -t713 * t453 + t581 * t747 + t741 * qJDD(1) + (-t212 - t769 + t786) * qJD(1) + t697;
t25 = -t147 * t609 - t148 * t611 + t153 * t891 + t154 * t347 + t155 * t348 + t607 * t631;
t24 = -t145 * t609 - t146 * t611 + t153 * t892 + t154 * t349 - t155 * t350 - t356 * t607;
t21 = t413 * t63 - t414 * t64 + t520 * t99;
t9 = t147 * t169 - t148 * t171 - t166 * t631 + t347 * t82 + t348 * t84 + t80 * t891;
t8 = t147 * t167 + t148 * t170 - t164 * t631 + t347 * t83 + t348 * t85 + t81 * t891;
t7 = t145 * t169 - t146 * t171 + t166 * t356 + t349 * t82 - t350 * t84 + t80 * t892;
t6 = t145 * t167 + t146 * t170 + t164 * t356 + t349 * t83 - t350 * t85 + t81 * t892;
t2 = t224 * t52 + t225 * t53 + t25 * t520 + t399 * t88 + t413 * t8 - t414 * t9;
t1 = t224 * t54 + t225 * t55 + t24 * t520 + t399 * t89 + t413 * t6 - t414 * t7;
t15 = [(-(qJD(1) * t343 + t633 - t93 + t991) * t94 + t93 * t778 + t94 * (-pkin(3) * t772 + t469 + t732 + t839) + t93 * (-t531 + (-t946 + (-rSges(5,3) - qJ(4)) * t553) * qJD(3)) * t580 + ((-t93 * rSges(5,1) + t618 * t94) * t580 + (t93 * (t618 + t945) - t94 * t579) * t581) * qJD(1) + (t31 - g(2)) * (t342 + t650) + (t30 - g(1)) * ((-t532 + (rSges(5,2) - pkin(3)) * t553) * t580 + t756 + t837)) * m(5) + (-(-qJD(1) * t340 - t124 + t722 + t995) * t125 + t124 * t835 + t125 * (t556 + t840) + (t435 * t908 - t907) * qJD(3) + ((-t124 * rSges(4,3) + t125 * (-t539 - t949)) * t580 + (t124 * (-t438 - t539) - t125 * t579) * t581) * qJD(1) + (t57 - g(2)) * (t341 + t748) + (t56 - g(1)) * (-t340 + t837)) * m(4) + (t241 * t557 + t242 * (t832 + t836) + (t241 * (t766 + t947) * t581 + (t241 * (-rSges(3,3) - qJ(2)) + t242 * t766) * t580) * qJD(1) - (-qJD(1) * t364 - t241 - t459 + t556) * t242 + (t123 - g(2)) * t286 + (t122 - g(1)) * (t580 * t766 + t559 + t834)) * m(3) + (-t647 * t820 + (t1029 * t580 + t1040) * t20 + (t468 - t720 + t778 + ((rSges(6,3) * qJD(3) - qJD(4)) * t551 + (-qJ(4) * qJD(3) - qJD(5)) * t553) * t580 + (t553 * t730 - t495 - t508 - t568) * qJD(1)) * t60 - g(1) * t1040 - t1029 * t959 + (t730 * t823 + t549 + t60 + t877 + (-qJ(4) * t896 + t553 * t731 + t209 + t837) * qJD(1) + t1027) * t61 + (t19 - g(2)) * (t650 + t207 + t415)) * m(6) + (m(2) * (t489 ^ 2 + t491 ^ 2) + Icges(3,2) * t577 ^ 2 + (Icges(3,1) * t575 + 0.2e1 * Icges(3,4) * t577) * t575 + t293 * t551 - t553 * t667 + Icges(2,3) - t1058) * qJDD(1) + t13 * t969 + (t1013 + t1016) * t761 + (t1014 - t1015 + t1018) * t760 + (t1007 + t1036) * t967 - t453 * t152 / 0.2e1 + (-g(1) * t1039 - (t551 * t757 + t553 * t961) * t959 + (-g(2) + t5) * (t650 + t652 + t177) + (t580 * t604 + t1039) * t4 + (-t717 + t778 + (-t539 + t649 + (-pkin(3) - t933) * t553) * t824 + (-t653 + (t551 * t933 + t553 * t757) * qJD(3) - t538 * qJD(1)) * t580) * t46 + (t961 * t823 * t581 + t46 + t735 + t883 + (-t544 + (-t539 + t604) * t580 + t223) * qJD(1) + t1027 + t1053) * t47) * m(7) + t89 * t974 + t88 * t975 + t25 * t971 + (t1006 + t1008) * t968 + ((-t272 * t576 - t273 * t574 + t1063) * t553 + (t271 + t1062) * t551 + (t293 * t553 + t551 * t667 + t1059) * qJD(3)) * qJD(1) - m(2) * (-g(1) * t489 + g(2) * t491) + t935 / 0.2e1 + t936 / 0.2e1 - t941 / 0.2e1 + t942 / 0.2e1 + t951 + (t24 + t13) * t970 + (((t1025 * t581 + t1021 + t1030 + t671) * t581 + (t1025 * t580 + t1010 - t1074 - t67 + t999) * t580) * qJD(3) + t1017 + t1037) * t762 + ((-t937 + ((t1075 + t1088) * t581 + t1031 + t1073 + t1087) * t581 + (t65 - t1021) * t580) * qJD(3) + t1038) * t759; (-m(3) - m(4) + t810) * (-g(2) * t581 + t959) + 0.2e1 * (t956 / 0.2e1 - t955 / 0.2e1) * m(7) + 0.2e1 * (t19 * t963 + t20 * t964) * m(6) + 0.2e1 * (t30 * t964 + t31 * t963) * m(5) + 0.2e1 * (t56 * t964 + t57 * t963) * m(4) + 0.2e1 * (t122 * t964 + t123 * t963) * m(3); (t13 * t581 + t14 * t580) * t814 / 0.2e1 + (-t598 * t814 - g(1) * (t497 + t842 + t848) - g(2) * (t494 + t843 + t849) - g(3) * (-t553 * t578 + t282 + t492 + t994) - t961 * t997 + t3 * t736 + (t4 * t651 + t3 * (t177 + t222)) * t581 + (t3 * t878 + t5 * t651) * t580 + (-t177 * t813 - t244 * t520 + t282 * t413 - t858 - (-qJ(5) * t824 - t816) * t551 + t1001 * t580 + (t581 * t651 - t262) * qJD(1)) * t47 + (t179 * t813 + t243 * t520 + t282 * t414 - t744 + t1001 * t581 + (t261 + (t283 + t305) * t580) * qJD(1) + t1000) * t46 + (-t243 * t413 - t244 * t414 - (t261 * t580 + t262 * t581 - t728) * qJD(3) + (qJD(1) * t878 + t121 + t87) * t581 + (t120 + t86 + (-t177 + t784) * qJD(1)) * t580 + t1028) * t33) * m(7) + (t12 * t736 + (qJD(1) * t647 + t12 * t207 + t20 * t725) * t581 + (t12 * t209 + t19 * t725) * t580 - g(1) * (t497 + t844) - g(2) * (t494 + t845) - (g(1) * t730 + g(2) * t731) * t551 + (t551 * t816 - t858 - (t256 - t793) * qJD(1) + t617 * t580) * t61 + (t617 * t581 - t487 + (t301 * t580 + t255 + t378) * qJD(1) + t1000) * t60 + (-(t580 * t61 + t581 * t60) * qJD(3) + g(3)) * (-rSges(6,1) * t897 - t551 * t944 - t536 + t753) + (-(t255 * t580 + t256 * t581 - t728) * qJD(3) + (qJD(1) * t209 + t119) * t581 + (qJD(1) * t785 + t118) * t580 + t1028) * t62) * m(6) + ((t1009 * t581 + t1010 * t580) * qJD(1) + t1003) * t761 + ((t1011 * t581 + t1012 * t580) * qJD(1) + t1002) * t760 + (qJD(1) * t1013 + qJD(3) * t1003 + qJDD(1) * t1008 + t1009 * t452 + t1010 * t453 + t2) * t964 + (t1015 * t581 + t1016 * t580 + (t1006 * t581 + t1007 * t580) * qJD(1)) * qJD(1) / 0.2e1 + (t14 + t1017) * t825 / 0.2e1 + (t13 + t1018) * t763 + ((t252 * t395 + t254 * t396) * t821 + (t294 * t395 + t296 * t396) * qJD(1) + (-t1020 * t821 + t1022) * t580 + ((t1019 * t580 - t251 * t395 - t253 * t396 + t1004) * qJD(3) + t1005) * t581) * t762 + (-(t251 * t397 - t253 * t398) * t820 + (t294 * t397 - t296 * t398) * qJD(1) + (-t1019 * t820 - t1022) * t581 + ((t1020 * t581 + t252 * t397 - t254 * t398 + t1004) * qJD(3) + t1005) * t580) * t759 + (t1014 * qJD(1) + t1002 * qJD(3) - qJDD(1) * t1055 + t1011 * t452 + t1012 * t453 + t1) * t963 - (t586 * t551 + (t1049 * t551 + ((-t252 * t576 - t254 * t574 + t186) * t580 - (-t251 * t576 - t253 * t574 + t188) * t581 + t1046) * t553) * qJD(3) + (t1047 * t551 + (-t294 * t576 - t296 * t574 - t1048 + t293) * t553) * qJD(1)) * qJD(1) / 0.2e1 + t700 * t973 + t702 * t974 + t704 * t975 + ((t238 * t349 - t240 * t350) * t413 - (t237 * t349 - t239 * t350) * t414 + (t276 * t349 - t278 * t350) * t520 + (-t54 * t895 + t553 * t89) * qJD(6) + ((-qJD(6) * t55 + t629) * t551 + t587) * t580) * t969 + (qJD(1) * t703 + t580 * t6 - t581 * t7) * t970 + (qJD(1) * t705 + t580 * t8 - t581 * t9) * t971 + ((t238 * t347 + t240 * t348) * t413 - (t237 * t347 + t239 * t348) * t414 + (t276 * t347 + t278 * t348) * t520 + (-t53 * t896 + t553 * t88) * qJD(6) + ((-qJD(6) * t52 + t629) * t551 + t587) * t581) * t972 + (qJD(1) * t701 + t10 * t580 - t11 * t581) * t965 + (-(t124 * t380 - t907) * qJD(1) - (t180 * (-t380 * t580 - t384 * t581) + t676 * t438) * qJD(3) + g(1) * t384 + g(2) * t380 - g(3) * t438 + (qJD(3) * t669 + t340 * t452 - t341 * t453) * t661 + t180 * ((t340 * t581 - t341 * t580) * qJD(1) + t669) + t676 * t410 + (-t56 * t581 - t57 * t580 + (-t125 * t581 + t908) * qJD(1)) * t435) * m(4) + (t1006 * t580 - t1007 * t581) * qJDD(1) / 0.2e1 + t1041 * t968 + t1042 * t967 - t21 * t813 / 0.2e1 + (((-t238 * t552 - t240 * t550 + t164) * t413 - (-t237 * t552 - t239 * t550 + t166) * t414 + (-t276 * t552 - t278 * t550 - t607) * t520 + t99 * qJD(6)) * t553 + (-qJD(6) * t701 + t588) * t551) * t966 + (-g(1) * (t382 + t383) - g(2) * (t378 + t379) + g(3) * t846 + t93 * t393 + t32 * t859 + t100 * t787 + (t30 * t847 + t93 * t862 + t32 * t342 + t100 * t213 + (-t100 * t343 + t847 * t94) * qJD(1)) * t581 + (t31 * t847 + t94 * t862 - t32 * t343 + t100 * t212 + (t100 * t860 - t713 * t93) * qJD(1)) * t580 - t93 * (-qJD(1) * t379 + t744) - t94 * (qJD(1) * t383 + t858) - t100 * t780 - ((t100 * t383 + t846 * t93) * t581 + (t100 * t379 + t846 * t94) * t580) * qJD(3)) * m(5); t810 * (-g(3) * t553 + t997) - m(5) * (t100 * t354 + t357 * t94 + t359 * t93) - m(6) * (t354 * t62 + t357 * t61 + t359 * t60) - m(7) * (t33 * t354 + t357 * t47 + t359 * t46) + ((t820 * t93 + t821 * t94 - t32) * t979 + t621 * t978 + t624 * t977) * t980 + ((qJD(3) * t100 + t30 * t581 + t31 * t580 + t824 * t94 - t825 * t93) * t979 + t584) * t981; t976 * (g(3) * t551 + t553 * t727) - m(6) * (t355 * t62 + t356 * t61 - t60 * t631) - m(7) * (t33 * t355 + t356 * t47 - t46 * t631) + (-m(6) * t621 / 0.2e1 - m(7) * t624 / 0.2e1) * t981 + t584 * t980; t2 * t891 / 0.2e1 + (t551 * t88 + t553 * t705) * t975 + ((-qJD(3) * t705 + t25) * t551 + (-qJD(1) * t704 + qJD(3) * t88 + t580 * t9 + t581 * t8) * t553) * t971 + t1 * t892 / 0.2e1 + (t551 * t89 + t553 * t703) * t974 + ((-qJD(3) * t703 + t24) * t551 + (-qJD(1) * t702 + qJD(3) * t89 + t580 * t7 + t581 * t6) * t553) * t970 + t21 * t822 / 0.2e1 + t551 * (t935 + t936 - t941 + t942 + t951) / 0.2e1 + (t551 * t99 + t553 * t701) * t973 + ((-qJD(3) * t701 + t29) * t551 + (-qJD(1) * t700 + qJD(3) * t99 + t10 * t581 + t11 * t580) * t553) * t965 + (t600 * t347 - t348 * t601 + t581 * t616) * t972 + (t349 * t600 + t350 * t601 + t580 * t616) * t969 + (t628 * t551 + (t601 * t550 - t552 * t600) * t553) * t966 + (t551 * t762 + t553 * t763) * t14 + (-t775 / 0.2e1 + t551 * t760) * t13 + ((qJD(3) * t598 + t5 * t177 - t4 * t179 - t46 * t86 + t47 * t87) * t551 + (t46 * (-qJD(3) * t179 + t156 * t580) + t47 * (qJD(3) * t177 - t156 * t581) - t3 * t673 + t33 * (-t177 * t824 - t179 * t825 - t580 * t87 + t581 * t86) + (t956 - t955 + (t46 * t581 + t47 * t580) * qJD(1)) * t283) * t553 - t46 * (-t221 * t520 - t335 * t414) - t47 * (t220 * t520 - t335 * t413) - t33 * (t220 * t414 + t221 * t413) - g(1) * t220 - g(2) * t221 - g(3) * t335) * m(7);];
tau  = t15;
