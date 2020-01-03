% Calculate vector of inverse dynamics joint torques for
% S5RRPRR15
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% qJDD [5x1]
%   Generalized joint accelerations
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d2,d4,d5]';
% m_mdh [6x1]
%   mass of all robot links (including the base)
% rSges [6x3]
%   center of mass of all robot links (in body frames)
%   rows: links of the robot (starting with base)
%   columns: x-, y-, z-coordinates
% Icges [6x6]
%   inertia of all robot links about their respective center of mass, in body frames
%   rows: links of the robot (starting with base)
%   columns: xx, yy, zz, xy, xz, yz (see inertiavector2matrix.m)
% 
% Output:
% tau [5x1]
%   joint torques of inverse dynamics (contains inertial, gravitational coriolis and centrifugal forces)

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 20:43
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tau = S5RRPRR15_invdynJ_fixb_slag_vp1(qJ, qJD, qJDD, g, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(5,1),zeros(3,1),zeros(8,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRPRR15_invdynJ_fixb_slag_vp1: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RRPRR15_invdynJ_fixb_slag_vp1: qJD has to be [5x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [5 1]), ...
  'S5RRPRR15_invdynJ_fixb_slag_vp1: qJDD has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RRPRR15_invdynJ_fixb_slag_vp1: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RRPRR15_invdynJ_fixb_slag_vp1: pkin has to be [8x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RRPRR15_invdynJ_fixb_slag_vp1: m has to be [6x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [6,3]), ...
  'S5RRPRR15_invdynJ_fixb_slag_vp1: rSges has to be [6x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [6 6]), ...
  'S5RRPRR15_invdynJ_fixb_slag_vp1: Icges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 20:41:12
% EndTime: 2019-12-31 20:43:07
% DurationCPUTime: 102.57s
% Computational Cost: add. (26397->1459), mult. (48679->1885), div. (0->0), fcn. (45354->8), ass. (0->668)
t1146 = Icges(4,4) - Icges(3,5);
t1145 = Icges(4,5) - Icges(3,6);
t1144 = Icges(4,1) + Icges(3,3);
t600 = sin(qJ(2));
t603 = cos(qJ(2));
t1117 = t1145 * t600 - t1146 * t603;
t963 = Icges(3,4) * t600;
t479 = Icges(3,2) * t603 + t963;
t948 = Icges(4,6) * t600;
t706 = Icges(4,3) * t603 + t948;
t1134 = -t479 - t706;
t583 = Icges(3,4) * t603;
t481 = Icges(3,1) * t600 + t583;
t947 = Icges(4,6) * t603;
t708 = Icges(4,2) * t600 + t947;
t1143 = t481 + t708;
t604 = cos(qJ(1));
t1142 = t1144 * t604;
t601 = sin(qJ(1));
t918 = t601 * t603;
t921 = t600 * t601;
t1120 = -t1145 * t921 + t1146 * t918 + t1142;
t1128 = t1117 * t604 + t1144 * t601;
t949 = Icges(3,6) * t604;
t343 = Icges(3,4) * t918 - Icges(3,2) * t921 - t949;
t539 = Icges(4,6) * t921;
t961 = Icges(4,4) * t604;
t352 = Icges(4,2) * t918 - t539 + t961;
t1141 = t343 * t600 - t352 * t603;
t482 = Icges(3,1) * t603 - t963;
t348 = Icges(3,5) * t601 + t482 * t604;
t707 = -Icges(4,3) * t600 + t947;
t349 = Icges(4,5) * t601 - t604 * t707;
t1140 = -t348 * t918 - t349 * t921;
t953 = Icges(4,5) * t604;
t350 = Icges(4,6) * t918 - Icges(4,3) * t921 + t953;
t1139 = t343 + t350;
t716 = -Icges(3,2) * t600 + t583;
t344 = Icges(3,6) * t601 + t604 * t716;
t1138 = t344 - t349;
t544 = Icges(3,4) * t921;
t954 = Icges(3,5) * t604;
t347 = Icges(3,1) * t918 - t544 - t954;
t1137 = t347 + t352;
t920 = t600 * t604;
t540 = Icges(4,6) * t920;
t915 = t603 * t604;
t962 = Icges(4,4) * t601;
t351 = -Icges(4,2) * t915 + t540 + t962;
t1136 = t348 - t351;
t475 = Icges(3,5) * t600 + Icges(3,6) * t603;
t714 = Icges(4,4) * t600 + Icges(4,5) * t603;
t1135 = t475 - t714;
t709 = Icges(4,2) * t603 - t948;
t1132 = t482 + t709;
t1131 = t1134 * qJD(2);
t1130 = t1143 * qJD(2);
t1113 = -t347 * t603 + t350 * t600 + t1141;
t1129 = t707 + t716;
t683 = t479 * t600 - t481 * t603;
t1115 = -t600 * t706 + t603 * t708 - t683;
t1127 = t1128 * t604 + t1140;
t1058 = t1128 * t601 + t348 * t915 + t349 * t920;
t1126 = t1120 * t601 - t347 * t915 + t350 * t920;
t1125 = t344 * t600 + t351 * t603;
t1076 = -t1113 * t601 + t1120 * t604;
t1075 = -t344 * t921 - t351 * t918 - t1127;
t1074 = -t343 * t920 + t352 * t915 - t1126;
t1073 = -t344 * t920 - t351 * t915 + t1058;
t1072 = t1137 * t600 + t1139 * t603;
t1071 = t1136 * t600 + t1138 * t603;
t1124 = t1131 * t604 + (-t1129 * t601 + t949 - t953) * qJD(1);
t1123 = t1138 * qJD(1) + t1131 * t601;
t1122 = -t1130 * t604 + (-t1132 * t601 + t954 - t961) * qJD(1);
t1121 = t1130 * t601 + (-t604 * t709 - t348 + t962) * qJD(1);
t1119 = t1129 * qJD(2);
t1118 = t1132 * qJD(2);
t1116 = t1135 * qJD(2);
t1114 = t1134 * t603 - t1143 * t600;
t1112 = t348 * t603 + t349 * t600 - t1125;
t1056 = t1135 * t601;
t928 = t475 * t604;
t206 = -t601 * t683 - t928;
t419 = t714 * t604;
t209 = t706 * t921 - t708 * t918 - t419;
t1111 = t206 - t209;
t1070 = t1115 * t604 + t1056;
t1110 = t1128 * qJD(1);
t574 = qJD(4) * t600;
t1042 = qJD(1) + t574;
t598 = qJ(4) + qJ(5);
t577 = sin(t598);
t578 = cos(t598);
t370 = t577 * t604 + t578 * t921;
t371 = -t577 * t921 + t578 * t604;
t890 = t371 * rSges(6,1) - t370 * rSges(6,2);
t185 = rSges(6,3) * t918 - t890;
t559 = pkin(7) * t918;
t456 = pkin(3) * t604 - t559;
t599 = sin(qJ(4));
t819 = t599 * t921;
t602 = cos(qJ(4));
t1001 = pkin(4) * t602;
t567 = pkin(3) + t1001;
t605 = -pkin(8) - pkin(7);
t914 = t603 * t605;
t860 = t604 * t567 + t601 * t914;
t259 = pkin(4) * t819 + t456 - t860;
t741 = rSges(6,1) * t577 + rSges(6,2) * t578;
t316 = rSges(6,3) * t600 - t603 * t741;
t838 = t604 * qJD(2);
t842 = qJD(4) * t603;
t443 = -t601 * t842 + t838;
t841 = qJD(5) * t603;
t336 = -t601 * t841 + t443;
t1002 = pkin(4) * t599;
t993 = -pkin(7) - t605;
t379 = -t1002 * t603 + t600 * t993;
t459 = qJD(5) * t600 + t1042;
t843 = qJD(3) * t604;
t533 = t600 * t843;
t1000 = pkin(7) * t600;
t489 = pkin(2) * t600 - qJ(3) * t603;
t781 = -t489 - t1000;
t651 = t781 * t838 + t533;
t1109 = t1042 * t259 + t185 * t459 + t336 * t316 + t443 * t379 - t651;
t919 = t601 * t602;
t408 = t599 * t604 + t600 * t919;
t916 = t602 * t604;
t409 = -t819 + t916;
t878 = t409 * rSges(5,1) - t408 * rSges(5,2);
t234 = rSges(5,3) * t918 - t878;
t987 = rSges(5,2) * t602;
t744 = rSges(5,1) * t599 + t987;
t364 = rSges(5,3) * t600 - t603 * t744;
t1108 = -t1042 * t234 - t443 * t364 + t651;
t1107 = t604 ^ 2;
t1106 = t1135 * qJD(1) + t1114 * qJD(2) + t1118 * t603 - t1119 * t600;
t1105 = -t1071 * qJD(2) + t1122 * t603 - t1124 * t600 + t1110;
t1104 = t1120 * qJD(1) + t1072 * qJD(2) + t1121 * t603 + t1123 * t600;
t1103 = -t1138 * t601 + t1139 * t604;
t1102 = t1073 * t601 - t1074 * t604;
t1101 = t1075 * t601 - t1076 * t604;
t1100 = t1129 + t1143;
t1099 = t1132 + t1134;
t1098 = (t539 + t544 + (Icges(3,2) + Icges(4,3)) * t918 - t1137) * t604 + (-Icges(4,3) * t915 - t479 * t604 + t1136 - t540) * t601;
t1097 = t1115 * qJD(1) - t1117 * qJD(2);
t1096 = t1113 * qJD(1) - t1116 * t601 + t1110;
t1095 = -t1116 * t604 + (-t1117 * t601 - t1112 + t1142) * qJD(1);
t1094 = t1070 * qJD(1);
t593 = t604 * pkin(6);
t1093 = t593 + t878;
t1092 = t593 + t860 + t890;
t1091 = t1111 * qJD(1);
t587 = t603 * rSges(6,3);
t315 = t600 * t741 + t587;
t1003 = -rSges(6,3) - pkin(2);
t782 = -qJ(3) - t1002;
t671 = t782 * t600 - pkin(1);
t1090 = t1003 * t603 + t671;
t710 = Icges(6,5) * t577 + Icges(6,6) * t578;
t635 = -Icges(6,3) * t600 + t603 * t710;
t956 = Icges(6,4) * t577;
t712 = Icges(6,2) * t578 + t956;
t637 = -Icges(6,6) * t600 + t603 * t712;
t955 = Icges(6,4) * t578;
t717 = Icges(6,1) * t577 + t955;
t639 = -Icges(6,5) * t600 + t603 * t717;
t109 = -t370 * t637 + t371 * t639 - t635 * t918;
t531 = t604 * t842;
t576 = qJD(2) * t601;
t442 = t576 + t531;
t335 = t604 * t841 + t442;
t368 = -t577 * t601 + t578 * t920;
t369 = t577 * t920 + t578 * t601;
t174 = Icges(6,5) * t369 + Icges(6,6) * t368 + Icges(6,3) * t915;
t957 = Icges(6,4) * t369;
t177 = Icges(6,2) * t368 + Icges(6,6) * t915 + t957;
t324 = Icges(6,4) * t368;
t180 = Icges(6,1) * t369 + Icges(6,5) * t915 + t324;
t74 = t174 * t918 + t370 * t177 - t371 * t180;
t176 = -Icges(6,5) * t371 + Icges(6,6) * t370 + Icges(6,3) * t918;
t326 = Icges(6,4) * t371;
t179 = Icges(6,2) * t370 + Icges(6,6) * t918 - t326;
t325 = Icges(6,4) * t370;
t181 = Icges(6,1) * t371 - Icges(6,5) * t918 - t325;
t75 = t176 * t918 + t179 * t370 + t181 * t371;
t37 = t109 * t459 + t335 * t74 - t336 * t75;
t711 = Icges(5,5) * t599 + Icges(5,6) * t602;
t636 = -Icges(5,3) * t600 + t603 * t711;
t959 = Icges(5,4) * t599;
t713 = Icges(5,2) * t602 + t959;
t638 = -Icges(5,6) * t600 + t603 * t713;
t958 = Icges(5,4) * t602;
t718 = Icges(5,1) * t599 + t958;
t640 = -Icges(5,5) * t600 + t603 * t718;
t122 = -t408 * t638 + t409 * t640 - t636 * t918;
t820 = t600 * t916;
t923 = t599 * t601;
t406 = t820 - t923;
t824 = t599 * t920;
t407 = t824 + t919;
t221 = Icges(5,5) * t407 + Icges(5,6) * t406 + Icges(5,3) * t915;
t960 = Icges(5,4) * t407;
t224 = Icges(5,2) * t406 + Icges(5,6) * t915 + t960;
t381 = Icges(5,4) * t406;
t227 = Icges(5,1) * t407 + Icges(5,5) * t915 + t381;
t81 = t221 * t918 + t408 * t224 - t409 * t227;
t223 = -Icges(5,5) * t409 + Icges(5,6) * t408 + Icges(5,3) * t918;
t383 = Icges(5,4) * t409;
t226 = Icges(5,2) * t408 + Icges(5,6) * t918 - t383;
t382 = Icges(5,4) * t408;
t228 = Icges(5,1) * t409 - Icges(5,5) * t918 - t382;
t82 = t223 * t918 + t226 * t408 + t228 * t409;
t44 = t1042 * t122 + t442 * t81 - t443 * t82;
t121 = -t406 * t638 - t407 * t640 - t636 * t915;
t79 = t221 * t915 + t406 * t224 + t407 * t227;
t80 = t223 * t915 + t406 * t226 - t228 * t407;
t43 = t1042 * t121 + t442 * t79 - t443 * t80;
t108 = -t368 * t637 - t369 * t639 - t635 * t915;
t72 = t174 * t915 + t368 * t177 + t369 * t180;
t73 = t176 * t915 + t368 * t179 - t181 * t369;
t36 = t108 * t459 + t335 * t72 - t336 * t73;
t697 = t226 * t602 - t228 * t599;
t88 = t223 * t600 - t603 * t697;
t700 = t179 * t578 - t181 * t577;
t78 = t176 * t600 - t603 * t700;
t1089 = qJD(2) * t1101 + t1091;
t1088 = qJD(2) * t1102 + t1094;
t1087 = t1112 * qJD(2) + t1122 * t600 + t1124 * t603;
t1086 = t1113 * qJD(2) + t1121 * t600 - t1123 * t603;
t1078 = -t1097 * t601 + t1106 * t604;
t1077 = t1097 * t604 + t1106 * t601;
t1069 = -t1098 * t600 + t1103 * t603;
t1068 = (t1099 * t603 - t1100 * t600) * qJD(1);
t1067 = t1096 * t1107 + (t1105 * t601 + (-t1095 + t1104) * t604) * t601;
t1066 = t1104 * t1107 + (t1095 * t601 + (-t1096 + t1105) * t604) * t601;
t785 = -t838 / 0.2e1;
t849 = qJD(1) * t603;
t807 = t601 * t849;
t1065 = t600 * t785 - t807 / 0.2e1;
t789 = -t576 / 0.2e1;
t848 = qJD(1) * t604;
t790 = t848 / 0.2e1;
t1064 = t600 * t789 + t603 * t790;
t1063 = t1120 + t1125;
t836 = qJD(2) * qJD(3);
t1062 = qJDD(3) * t600 + t603 * t836;
t803 = t600 * t838;
t1061 = t803 + t807;
t926 = t578 * t603;
t927 = t577 * t603;
t1060 = rSges(6,1) * t927 + rSges(6,2) * t926;
t1059 = t1117 * qJD(1);
t1057 = t928 - t419;
t183 = t369 * rSges(6,1) + t368 * rSges(6,2) + rSges(6,3) * t915;
t455 = t601 * pkin(3) + pkin(7) * t915;
t912 = t604 * t605;
t674 = pkin(4) * t824 + t601 * t567 - t603 * t912;
t258 = t674 - t455;
t579 = t600 * qJ(3);
t1046 = t603 * pkin(2) + t579;
t430 = t1046 * t601;
t435 = pkin(2) * t915 + qJ(3) * t920;
t844 = qJD(3) * t603;
t749 = t430 * t576 + t435 * t838 - t844;
t663 = t455 * t838 - t456 * t576 + t749;
t52 = t183 * t336 + t185 * t335 + t258 * t443 + t259 * t442 + t663;
t903 = t185 + t259;
t1054 = t52 * t903;
t502 = pkin(1) * t601 - t593;
t876 = -t430 - t502;
t812 = t456 + t876;
t751 = t812 * qJD(1);
t66 = -t1109 + t751;
t891 = t316 + t379;
t1053 = t66 * t891;
t586 = t601 * rSges(4,1);
t366 = -rSges(4,2) * t915 + rSges(4,3) * t920 + t586;
t503 = t604 * pkin(1) + t601 * pkin(6);
t759 = t435 + t503;
t236 = t366 + t759;
t465 = qJD(1) * t502;
t1050 = -qJD(1) * t430 - t465;
t1049 = t455 + t759;
t780 = -pkin(7) * t603 - t1046;
t834 = qJDD(4) * t603;
t1048 = qJD(1) * t531 + t601 * t834;
t924 = t599 * t600;
t558 = pkin(4) * t924;
t1047 = t993 * t603 + t558;
t584 = t600 * rSges(4,3);
t740 = -rSges(4,2) * t603 + t584;
t783 = -pkin(1) - t579;
t1004 = -rSges(5,3) - pkin(2);
t833 = -pkin(7) + t1004;
t1045 = t833 * t603 + t783;
t850 = qJD(1) * t601;
t1041 = qJD(1) * t456 + t1050;
t1040 = (g(1) * t604 + g(2) * t601) * t600;
t804 = t600 * t576;
t522 = pkin(7) * t804;
t846 = qJD(2) * t603;
t801 = t601 * t846;
t845 = qJD(2) * t605;
t802 = t600 * t845;
t123 = t601 * t802 + t522 + (qJD(4) * t408 + t599 * t801) * pkin(4) + ((-pkin(3) + t567) * t601 + t1047 * t604) * qJD(1);
t571 = pkin(3) * t848;
t500 = pkin(4) * t820;
t822 = t599 * t915;
t501 = pkin(4) * t822;
t720 = qJD(2) * t501 + qJD(4) * t500 + t567 * t848 + t604 * t802 + t605 * t807;
t851 = qJD(1) * t600;
t763 = qJD(4) + t851;
t124 = pkin(7) * t803 - t571 + (pkin(7) * t849 - t1002 * t763) * t601 + t720;
t595 = qJD(4) + qJD(5);
t925 = t595 * t600;
t721 = qJD(2) * t925;
t837 = qJD(1) * qJD(2);
t457 = qJDD(2) * t601 + t604 * t837;
t809 = t604 * t834 + t457;
t147 = -t604 * t721 + (qJDD(5) * t604 - t595 * t850) * t603 + t809;
t565 = t601 * t837;
t148 = t565 + (qJD(1) * t841 - qJDD(2)) * t604 + (qJDD(5) * t603 - t721) * t601 + t1048;
t260 = -qJD(4) * t1061 + t809;
t458 = -qJDD(2) * t604 + t565;
t261 = -qJD(4) * t804 + t1048 + t458;
t321 = qJD(1) * t455 - t522;
t322 = -pkin(7) * t1061 + t571;
t808 = t600 * t850;
t800 = t603 * t838;
t859 = qJ(3) * t800 + t533;
t230 = -pkin(2) * t1061 - qJ(3) * t808 + t859;
t404 = t600 * t848 + t801;
t523 = pkin(2) * t804;
t575 = qJD(3) * t600;
t797 = t601 * t575;
t806 = t603 * t848;
t231 = pkin(2) * t806 + qJ(3) * t404 - t523 + t797;
t656 = -qJDD(3) * t603 + t230 * t838 + t231 * t576 + t457 * t430 + t600 * t836;
t875 = t435 + t455;
t624 = t321 * t576 + t322 * t838 - t457 * t456 - t458 * t875 + t656;
t653 = -t804 + t806;
t767 = t595 + t851;
t643 = t604 * t767 + t801;
t768 = qJD(1) + t925;
t682 = t577 * t768;
t151 = t578 * t643 - t601 * t682;
t681 = t578 * t768;
t152 = t577 * t643 + t601 * t681;
t743 = rSges(6,1) * t152 + rSges(6,2) * t151;
t95 = rSges(6,3) * t653 + t743;
t642 = -t601 * t767 + t800;
t153 = t578 * t642 - t604 * t682;
t154 = t577 * t642 + t604 * t681;
t911 = t154 * rSges(6,1) + t153 * rSges(6,2);
t96 = -rSges(6,3) * t1061 + t911;
t11 = t123 * t442 + t124 * t443 + t147 * t185 - t148 * t183 - t258 * t261 + t259 * t260 + t335 * t95 + t336 * t96 + t624;
t904 = t183 + t258;
t1033 = t11 * t904 + t52 * (t124 + t96);
t373 = (Icges(6,2) * t577 - t955) * t603;
t628 = t335 * (-Icges(6,2) * t369 + t180 + t324) - t336 * (Icges(6,2) * t371 - t181 + t325) + t459 * (-t639 + t373);
t374 = (-Icges(6,1) * t578 + t956) * t603;
t629 = t335 * (-Icges(6,1) * t368 + t177 + t957) - t336 * (-Icges(6,1) * t370 + t179 - t326) + t459 * (-t637 - t374);
t417 = (Icges(5,2) * t599 - t958) * t603;
t626 = t442 * (-Icges(5,2) * t407 + t227 + t381) - t443 * (Icges(5,2) * t409 - t228 + t382) + t1042 * (-t640 + t417);
t422 = (-Icges(5,1) * t602 + t959) * t603;
t627 = t442 * (-Icges(5,1) * t406 + t224 + t960) - t443 * (-Icges(5,1) * t408 + t226 - t383) + t1042 * (-t638 - t422);
t1030 = m(4) / 0.2e1;
t1029 = m(5) / 0.2e1;
t1028 = m(6) / 0.2e1;
t1027 = t147 / 0.2e1;
t1026 = t148 / 0.2e1;
t1025 = t260 / 0.2e1;
t1024 = t261 / 0.2e1;
t441 = qJD(2) * t842 + qJDD(4) * t600 + qJDD(1);
t307 = qJD(2) * t841 + qJDD(5) * t600 + t441;
t1023 = t307 / 0.2e1;
t1022 = -t335 / 0.2e1;
t1021 = t335 / 0.2e1;
t1020 = -t336 / 0.2e1;
t1019 = t336 / 0.2e1;
t1018 = t441 / 0.2e1;
t1017 = -t442 / 0.2e1;
t1016 = t442 / 0.2e1;
t1015 = -t443 / 0.2e1;
t1014 = t443 / 0.2e1;
t1013 = t457 / 0.2e1;
t1012 = t458 / 0.2e1;
t1011 = -t459 / 0.2e1;
t1010 = t459 / 0.2e1;
t1009 = -t1042 / 0.2e1;
t1008 = t1042 / 0.2e1;
t1007 = t600 / 0.2e1;
t998 = g(1) * t601;
t992 = rSges(3,1) * t603;
t990 = rSges(6,1) * t578;
t989 = rSges(4,2) * t600;
t985 = pkin(4) * qJD(4);
t701 = t177 * t578 + t180 * t577;
t90 = Icges(6,5) * t154 + Icges(6,6) * t153 - Icges(6,3) * t1061;
t92 = Icges(6,4) * t154 + Icges(6,2) * t153 - Icges(6,6) * t1061;
t94 = Icges(6,1) * t154 + Icges(6,4) * t153 - Icges(6,5) * t1061;
t20 = (qJD(2) * t701 + t90) * t600 + (qJD(2) * t174 + (-t180 * t595 - t92) * t578 + (t177 * t595 - t94) * t577) * t603;
t984 = t20 * t335;
t89 = Icges(6,5) * t152 + Icges(6,6) * t151 + Icges(6,3) * t653;
t91 = Icges(6,4) * t152 + Icges(6,2) * t151 + Icges(6,6) * t653;
t93 = Icges(6,1) * t152 + Icges(6,4) * t151 + Icges(6,5) * t653;
t21 = (qJD(2) * t700 + t89) * t600 + (qJD(2) * t176 + (t181 * t595 - t91) * t578 + (t179 * t595 - t93) * t577) * t603;
t983 = t21 * t336;
t641 = -t601 * t763 + t800;
t678 = t1042 * t599;
t204 = t602 * t641 - t604 * t678;
t679 = t602 * t1042;
t205 = t599 * t641 + t604 * t679;
t112 = Icges(5,5) * t205 + Icges(5,6) * t204 - Icges(5,3) * t1061;
t114 = Icges(5,4) * t205 + Icges(5,2) * t204 - Icges(5,6) * t1061;
t116 = Icges(5,1) * t205 + Icges(5,4) * t204 - Icges(5,5) * t1061;
t698 = t224 * t602 + t227 * t599;
t28 = (qJD(2) * t698 + t112) * t600 + (qJD(2) * t221 - t114 * t602 - t116 * t599 + (t224 * t599 - t227 * t602) * qJD(4)) * t603;
t982 = t28 * t442;
t677 = t763 * t604;
t202 = t602 * t677 + (t602 * t846 - t678) * t601;
t203 = t601 * t679 + (t677 + t801) * t599;
t111 = Icges(5,5) * t203 + Icges(5,6) * t202 + Icges(5,3) * t653;
t113 = Icges(5,4) * t203 + Icges(5,2) * t202 + Icges(5,6) * t653;
t115 = Icges(5,1) * t203 + Icges(5,4) * t202 + Icges(5,5) * t653;
t29 = (qJD(2) * t697 + t111) * t600 + (qJD(2) * t223 - t113 * t602 - t115 * t599 + (t226 * t599 + t228 * t602) * qJD(4)) * t603;
t981 = t29 * t443;
t585 = t601 * rSges(3,3);
t588 = t603 * rSges(5,3);
t77 = t174 * t600 - t603 * t701;
t976 = t77 * t147;
t975 = t78 * t148;
t87 = t221 * t600 - t603 * t698;
t974 = t87 * t260;
t973 = t88 * t261;
t971 = -rSges(4,3) - qJ(3);
t694 = -t577 * t639 - t578 * t637;
t126 = -t600 * t635 - t603 * t694;
t309 = Icges(6,3) * t603 + t600 * t710;
t372 = (-Icges(6,5) * t578 + Icges(6,6) * t577) * t603;
t169 = qJD(2) * t309 + t372 * t595;
t311 = Icges(6,6) * t603 + t600 * t712;
t170 = qJD(2) * t311 + t373 * t595;
t313 = Icges(6,5) * t603 + t600 * t717;
t171 = qJD(2) * t313 + t374 * t595;
t51 = (qJD(2) * t694 + t169) * t600 + (-qJD(2) * t635 + (t595 * t639 - t170) * t578 + (-t595 * t637 - t171) * t577) * t603;
t970 = t126 * t307 + t51 * t459;
t969 = t123 + t95;
t693 = -t599 * t640 - t602 * t638;
t133 = -t600 * t636 - t603 * t693;
t337 = Icges(5,3) * t603 + t600 * t711;
t414 = (-Icges(5,5) * t602 + Icges(5,6) * t599) * t603;
t237 = qJD(2) * t337 + qJD(4) * t414;
t341 = Icges(5,6) * t603 + t600 * t713;
t240 = qJD(2) * t341 + qJD(4) * t417;
t345 = Icges(5,5) * t603 + t600 * t718;
t243 = qJD(2) * t345 + qJD(4) * t422;
t55 = (qJD(2) * t693 + t237) * t600 + (-qJD(2) * t636 - t240 * t602 - t243 * t599 + (-t599 * t638 + t602 * t640) * qJD(4)) * t603;
t967 = t1042 * t55 + t133 * t441;
t966 = t183 * t804 + t95 * t915;
t941 = qJD(2) * t66;
t491 = rSges(3,1) * t600 + rSges(3,2) * t603;
t805 = t491 * t838;
t857 = rSges(3,2) * t921 + t604 * rSges(3,3);
t363 = rSges(3,1) * t918 - t857;
t882 = -t363 - t502;
t192 = qJD(1) * t882 - t805;
t938 = t192 * t601;
t937 = t192 * t604;
t365 = rSges(3,1) * t915 - rSges(3,2) * t920 + t585;
t284 = t365 + t503;
t193 = qJD(1) * t284 - t491 * t576;
t434 = t491 * t604;
t936 = t193 * t434;
t470 = t595 * t603;
t917 = t602 * t603;
t913 = t603 * qJD(2) ^ 2;
t173 = (rSges(6,2) * t577 - t990) * t470 + t315 * qJD(2);
t910 = t173 * t918 + t316 * t806;
t825 = t602 * t985;
t281 = qJD(2) * t1047 - t603 * t825;
t909 = -t173 - t281;
t902 = t205 * rSges(5,1) + t204 * rSges(5,2);
t219 = t368 * rSges(6,1) - t369 * rSges(6,2);
t220 = t370 * rSges(6,1) + t371 * rSges(6,2);
t881 = -t366 - t435;
t880 = t601 * t430 + t604 * t435;
t384 = qJD(2) * t1046 - t844;
t879 = -t740 * qJD(2) - t384;
t537 = qJ(3) * t915;
t432 = -pkin(2) * t920 + t537;
t877 = qJD(1) * t432 + t601 * t844;
t437 = t489 * t850;
t524 = pkin(7) * t808;
t874 = t437 + t524;
t570 = pkin(6) * t848;
t873 = qJD(1) * (-pkin(1) * t850 + t570) + qJDD(1) * t503;
t872 = t1060 * t601;
t871 = t1060 * t604;
t823 = t599 * t918;
t830 = rSges(5,2) * t917;
t866 = rSges(5,1) * t823 + t601 * t830;
t865 = rSges(5,1) * t822 + t604 * t830;
t739 = rSges(4,3) * t603 + t989;
t864 = -t489 + t739;
t863 = -t1046 - t740;
t862 = pkin(4) * t823 + t605 * t921;
t861 = t600 * t912 + t501;
t858 = rSges(3,2) * t808 + rSges(3,3) * t848;
t428 = rSges(4,2) * t921 + rSges(4,3) * t918;
t433 = rSges(4,2) * t920 + rSges(4,3) * t915;
t386 = t408 * pkin(4);
t856 = t601 ^ 2 + t1107;
t847 = qJD(2) * t600;
t832 = pkin(4) * t917;
t826 = t599 * t985;
t436 = t489 * t576;
t630 = qJD(1) * t1049 - t436 - t522 + t797;
t67 = t1042 * t258 + t183 * t459 - t316 * t335 - t379 * t442 + t630;
t818 = t67 * t848;
t232 = t407 * rSges(5,1) + t406 * rSges(5,2) + rSges(5,3) * t915;
t84 = t1042 * t232 - t364 * t442 + t630;
t817 = t84 * t848;
t816 = t604 * t230 + t601 * t231 + t430 * t848;
t815 = t1062 * t604 + t458 * t489;
t535 = qJ(3) * t918;
t427 = -pkin(2) * t921 + t535;
t814 = t427 * t576 + t432 * t838 + t575;
t779 = rSges(4,1) * t604 - rSges(4,3) * t921;
t367 = rSges(4,2) * t918 + t779;
t813 = t367 + t876;
t810 = t570 + t859;
t362 = rSges(5,1) * t924 + t600 * t987 + t588;
t796 = t918 / 0.2e1;
t795 = t915 / 0.2e1;
t794 = -pkin(1) - t992;
t788 = t576 / 0.2e1;
t787 = t846 / 0.2e1;
t784 = t838 / 0.2e1;
t778 = t219 * t336 + t335 * t220;
t508 = rSges(6,2) * t927;
t376 = -rSges(6,1) * t926 + t508;
t777 = t459 * t219 - t335 * t376;
t776 = -t220 * t459 - t336 * t376;
t771 = qJD(2) * t879;
t770 = t856 * qJD(2);
t769 = -qJD(1) * t427 + t603 * t843;
t761 = t604 * t455 - t601 * t456 + t880;
t760 = rSges(4,1) * t848 + rSges(4,2) * t1061 + rSges(4,3) * t800;
t385 = -pkin(4) * t923 + t500;
t754 = -t364 + t781;
t753 = (rSges(4,2) - pkin(2)) * t603 - pkin(1);
t752 = -pkin(7) * t846 - t384;
t750 = t524 + t769;
t748 = t601 * t780;
t747 = t780 * t604;
t496 = rSges(2,1) * t604 - rSges(2,2) * t601;
t492 = rSges(2,1) * t601 + rSges(2,2) * t604;
t495 = -rSges(3,2) * t600 + t992;
t746 = rSges(5,1) * t203 + rSges(5,2) * t202;
t733 = t601 * t73 + t604 * t72;
t732 = t601 * t72 - t604 * t73;
t731 = t601 * t75 + t604 * t74;
t730 = t601 * t74 - t604 * t75;
t729 = t601 * t78 + t604 * t77;
t728 = t601 * t77 - t604 * t78;
t727 = t601 * t80 + t604 * t79;
t726 = t601 * t79 - t604 * t80;
t725 = t601 * t82 + t604 * t81;
t724 = t601 * t81 - t604 * t82;
t723 = t601 * t88 + t604 * t87;
t722 = t601 * t87 - t604 * t88;
t699 = -t193 * t601 - t937;
t696 = t232 * t601 - t234 * t604;
t253 = -rSges(3,1) * t1061 - rSges(3,2) * t800 + t858;
t429 = t491 * t601;
t254 = -qJD(2) * t429 + (t495 * t604 + t585) * qJD(1);
t695 = t253 * t604 + t254 * t601;
t687 = t363 * t601 + t365 * t604;
t680 = t781 - t891;
t431 = (-rSges(5,1) * t602 + rSges(5,2) * t599) * t603;
t252 = qJD(4) * t431 + (t600 * t744 + t588) * qJD(2);
t676 = -t252 + t752;
t675 = -pkin(1) - t1046;
t454 = t503 * qJD(1);
t673 = -t231 - t454 - t797;
t672 = t601 * t321 + t604 * t322 - t456 * t848 + t816;
t666 = t838 * t864 + t533;
t665 = t67 * t379 - t1054;
t664 = t752 + t909;
t655 = qJDD(1) * t435 + t873 + t1062 * t601 + (t230 + t533) * qJD(1);
t650 = (Icges(6,5) * t368 - Icges(6,6) * t369) * t335 - (Icges(6,5) * t370 + Icges(6,6) * t371) * t336 + t372 * t459;
t649 = t1042 * t636 - t221 * t442 + t223 * t443;
t648 = (Icges(5,5) * t406 - Icges(5,6) * t407) * t442 - (Icges(5,5) * t408 + Icges(5,6) * t409) * t443 + t414 * t1042;
t645 = t603 * t650;
t644 = t603 * t648;
t613 = qJD(1) * t322 + qJDD(1) * t455 + (-t457 * t600 - t601 * t913) * pkin(7) - t457 * t489 - t384 * t576 + t655;
t15 = t1042 * t124 - t147 * t316 - t335 * t173 + t307 * t183 + t441 * t258 - t260 * t379 - t442 * t281 + t459 * t96 + t613;
t631 = t11 * t185 * t915 + t15 * t600 * t183 + t67 * (t1061 * t316 + t183 * t846 + t600 * t96);
t16 = t151 * t177 + t152 * t180 + t174 * t653 + t370 * t92 - t371 * t94 + t90 * t918;
t17 = t151 * t179 - t152 * t181 + t176 * t653 + t370 * t91 - t371 * t93 + t89 * t918;
t18 = -t1061 * t174 + t153 * t177 + t154 * t180 + t368 * t92 + t369 * t94 + t90 * t915;
t19 = -t1061 * t176 + t153 * t179 - t154 * t181 + t368 * t91 + t369 * t93 + t89 * t915;
t41 = t126 * t459 + t335 * t77 - t336 * t78;
t45 = -t151 * t637 - t152 * t639 + t169 * t918 + t170 * t370 - t171 * t371 - t635 * t653;
t46 = t1061 * t635 - t153 * t637 - t154 * t639 + t169 * t915 + t170 * t368 + t171 * t369;
t5 = t109 * t307 + t147 * t74 + t148 * t75 + t16 * t335 - t17 * t336 + t45 * t459;
t6 = t108 * t307 + t147 * t72 + t148 * t73 + t18 * t335 - t19 * t336 + t459 * t46;
t625 = ((-qJD(2) * t733 + t46) * t600 + (-qJD(1) * t732 + qJD(2) * t108 + t18 * t604 + t19 * t601) * t603) * t1021 + (t628 * t368 - t369 * t629 + t604 * t645) * t1022 + (t370 * t628 + t371 * t629 + t601 * t645) * t1019 + ((-qJD(2) * t731 + t45) * t600 + (-qJD(1) * t730 + qJD(2) * t109 + t16 * t604 + t17 * t601) * t603) * t1020 + (t650 * t600 + (t629 * t577 - t578 * t628) * t603) * t1011 + t5 * t796 + (t108 * t600 + t603 * t733) * t1027 + (t109 * t600 + t603 * t731) * t1026 + t6 * t795 + t41 * t787 + (t126 * t600 + t603 * t729) * t1023 + ((-qJD(2) * t729 + t51) * t600 + (-qJD(1) * t728 + qJD(2) * t126 + t20 * t604 + t21 * t601) * t603) * t1010 + (t970 + t975 + t976 - t983 + t984) * t1007 + t1064 * t37 + t1065 * t36;
t76 = t232 * t443 + t234 * t442 + t663;
t83 = t1108 + t751;
t623 = t76 * t696 + (-t601 * t83 + t604 * t84) * t364;
t612 = (t635 * t604 + t701) * t335 - (t635 * t601 + t700) * t336 + (t309 + t694) * t459;
t611 = (t636 * t604 + t698) * t442 - (t636 * t601 + t697) * t443 + (t337 + t693) * t1042;
t610 = t611 * t603;
t609 = t458 * t1000 + t812 * qJDD(1) + (-pkin(7) * t913 - qJD(2) * t384) * t604 + (-t321 + t673) * qJD(1) + t815;
t608 = (-t174 * t335 + t176 * t336 + t459 * t635) * t600 + t612 * t603;
t451 = t495 * qJD(2);
t405 = t800 - t808;
t403 = t856 * t847;
t402 = t604 * t925;
t401 = t601 * t925;
t305 = pkin(7) * t920 + t861;
t304 = pkin(7) * t921 + t862;
t295 = -rSges(5,3) * t920 + t865;
t294 = -rSges(5,3) * t921 + t866;
t290 = t640 * t604;
t289 = t640 * t601;
t288 = t638 * t604;
t287 = t638 * t601;
t282 = t316 * t918;
t277 = -rSges(6,3) * t920 + t871;
t276 = -rSges(6,3) * t921 + t872;
t275 = t639 * t604;
t274 = t639 * t601;
t273 = t637 * t604;
t272 = t637 * t601;
t269 = rSges(5,1) * t408 + rSges(5,2) * t409;
t268 = rSges(5,1) * t406 - rSges(5,2) * t407;
t256 = -rSges(4,3) * t808 + t760;
t255 = t739 * t576 + (t604 * t740 + t586) * qJD(1);
t186 = t687 * qJD(2);
t129 = -t436 + (qJD(2) * t739 + t575) * t601 + t236 * qJD(1);
t128 = qJD(1) * t813 + t666;
t125 = (t366 * t604 - t367 * t601) * qJD(2) + t749;
t118 = -rSges(5,3) * t1061 + t902;
t117 = rSges(5,3) * t653 + t746;
t107 = qJD(1) * t253 + qJDD(1) * t365 - t451 * t576 - t457 * t491 + t873;
t106 = -t451 * t838 + t458 * t491 + t882 * qJDD(1) + (-t254 - t454) * qJD(1);
t65 = qJD(1) * t256 + qJDD(1) * t366 + t457 * t864 + t601 * t771 + t655;
t64 = -t458 * t739 + t604 * t771 + t813 * qJDD(1) + (-t255 + t673) * qJD(1) + t815;
t53 = -t367 * t457 + t881 * t458 + (t255 * t601 + t256 * t604) * qJD(2) + t656;
t50 = t1061 * t636 - t204 * t638 - t205 * t640 + t237 * t915 + t240 * t406 + t243 * t407;
t49 = -t202 * t638 - t203 * t640 + t237 * t918 + t240 * t408 - t243 * t409 - t636 * t653;
t47 = t1042 * t133 + t442 * t87 - t443 * t88;
t39 = t1042 * t118 + t232 * t441 - t252 * t442 - t260 * t364 + t613;
t38 = -t1042 * t117 - t234 * t441 - t252 * t443 + t261 * t364 + t609;
t27 = -t1061 * t223 + t111 * t915 + t113 * t406 + t115 * t407 + t204 * t226 - t205 * t228;
t26 = -t1061 * t221 + t112 * t915 + t114 * t406 + t116 * t407 + t204 * t224 + t205 * t227;
t25 = t111 * t918 + t113 * t408 - t115 * t409 + t202 * t226 - t203 * t228 + t223 * t653;
t24 = t112 * t918 + t114 * t408 - t116 * t409 + t202 * t224 + t203 * t227 + t221 * t653;
t23 = t117 * t442 + t118 * t443 - t232 * t261 + t234 * t260 + t624;
t14 = -t1042 * t123 + t148 * t316 - t336 * t173 - t307 * t185 - t441 * t259 + t261 * t379 - t443 * t281 - t459 * t95 + t609;
t10 = t1042 * t50 + t121 * t441 + t26 * t442 + t260 * t79 + t261 * t80 - t27 * t443;
t9 = t1042 * t49 + t122 * t441 + t24 * t442 - t25 * t443 + t260 * t81 + t261 * t82;
t1 = [(t1115 * qJD(2) + t1118 * t600 + t1119 * t603) * qJD(1) - t981 / 0.2e1 + (((t1063 * t604 - t1058 + t1073) * t604 + (t1063 * t601 + t1074 + t1127) * t601) * qJD(2) + t1089 - t1091) * t789 + (t1070 + t1071) * t1013 + (t206 + t1072) * t1012 + t973 / 0.2e1 + (t49 + t43) * t1015 + (-g(1) * (t456 + t1093) - (t1004 * t603 + t783) * t998 - (t1041 + t1108 - t83) * t84 + t38 * (-t559 + t1093) + t83 * (t522 + t523 - t746) + t84 * (t571 + t810 + t902) + (qJD(1) * t1045 * t83 + t84 * t833 * t847 + t38 * pkin(3)) * t604 + (t38 * (t675 - t588) + t83 * (rSges(5,3) * t847 - qJ(3) * t846 - t575) + (t83 * (-pkin(3) - pkin(6)) + t1045 * t84) * qJD(1)) * t601 + (-g(2) + t39) * (t232 + t1049)) * m(5) + (-g(1) * t1092 - t1090 * t998 + t782 * t941 * t918 + (t1090 * t601 + t1092) * t14 + (t523 - t743 + (-t826 + ((t605 + t1003) * t603 + t671) * qJD(1)) * t604 + ((rSges(6,3) * qJD(2) - qJD(3) - t825 - t845) * t600 + (-pkin(6) - t567) * qJD(1)) * t601) * t66 + (t720 + t810 + t911 + t1003 * t847 * t604 + (-t826 + (-t558 + t675 - t587) * qJD(1)) * t601 - t1041 + t66 + t1109) * t67 + (-g(2) + t15) * (t674 + t759 + t183)) * m(6) + t36 * t1019 + t43 * t1014 + (m(2) * (t492 ^ 2 + t496 ^ 2) + Icges(2,3) - t1114) * qJDD(1) + ((t523 + (-t575 + (t603 * t971 - t989) * qJD(2)) * t601 + ((t600 * t971 + t753) * t604 + (-rSges(4,1) - pkin(6)) * t601) * qJD(1)) * t128 + (-g(2) + t65) * t236 + (-g(1) + t64) * (t593 + (t753 - t579) * t601 + t779) + (-qJD(1) * t367 - t1050 + t128 - t666 - pkin(2) * t803 + t760 + t810 + (t675 - t584) * t850) * t129) * m(4) + (-(-qJD(1) * t363 - t192 - t465 - t805) * t193 + t193 * (t570 + t858) + (t491 * t938 - t936) * qJD(2) + ((-pkin(1) - t495) * t937 + (t192 * (-rSges(3,3) - pkin(6)) + t193 * t794) * t601) * qJD(1) + (t107 - g(2)) * t284 + (t106 - g(1)) * (t794 * t601 + t593 + t857)) * m(3) + (t45 + t36) * t1020 + t982 / 0.2e1 + t970 + t967 - t983 / 0.2e1 + t984 / 0.2e1 + (t1078 + t1087) * t788 + (t1077 - t1086 + t1088) * t785 + t974 / 0.2e1 + t975 / 0.2e1 + t976 / 0.2e1 + t50 * t1016 + t46 * t1021 + t122 * t1024 + t121 * t1025 + t109 * t1026 + t108 * t1027 - t458 * t209 / 0.2e1 + ((t1058 * t601 + ((t1128 + t1141) * t604 + t1075 + t1126 + t1140) * t604) * qJD(2) + t1094) * t784 - m(2) * (-g(1) * t492 + g(2) * t496); ((t273 * t370 - t275 * t371) * t335 - t74 * t402 - (t272 * t370 - t274 * t371) * t336 - t75 * t401 + (t311 * t370 - t313 * t371) * t459 + t109 * t470 + t608 * t601) * t1019 + (-t66 * (-t1042 * t304 - t1047 * t443 - t185 * t470 - t259 * t842 - t276 * t459 - t315 * t336 - t316 * t401 + t750) - t67 * (t1042 * t305 - t1047 * t442 + t183 * t470 + t258 * t842 + t277 * t459 - t315 * t335 + t316 * t402 + t877) - t52 * (t183 * t401 - t185 * t402 + t276 * t335 + t277 * t336 + t304 * t442 + t305 * t443 + t814) - (t66 * t747 + t67 * t748) * qJD(2) - ((-t52 * t770 - t818) * pkin(7) + (t52 * (t258 * t601 - t259 * t604) + (-t601 * t66 + t604 * t67) * t379) * qJD(4)) * t600 + t66 * t874 + t11 * t761 + t52 * t672 + (t14 * t680 + t66 * t664 + (t67 * t680 + t1054) * qJD(1) + t1033) * t604 + (t15 * t680 + t67 * t664 + t11 * t903 + t52 * t969 + (t1053 + t52 * (-t875 - t904)) * qJD(1)) * t601 - g(1) * (t537 + t861 + t871) - g(2) * (t535 + t862 + t872) - g(3) * (t315 + t1046 + t558 - t914) - t1003 * t1040) * m(6) + ((t273 * t368 + t275 * t369) * t335 - t72 * t402 - (t272 * t368 + t274 * t369) * t336 - t73 * t401 + (t311 * t368 + t313 * t369) * t459 + t108 * t470 + t608 * t604) * t1022 + (t1071 * t601 - t1072 * t604) * qJDD(1) / 0.2e1 + ((t1073 * t604 + t1074 * t601) * qJD(1) + t1066) * t788 + (((-t288 * t602 - t290 * t599 + t221) * t442 - (-t287 * t602 - t289 * t599 + t223) * t443 + (-t341 * t602 - t345 * t599 - t636) * t1042 + t133 * qJD(4)) * t603 + (-qJD(4) * t723 + t611) * t600) * t1009 + (t126 * t470 - t78 * t401 - t77 * t402 + ((-t273 * t578 - t275 * t577 + t174) * t335 - (-t272 * t578 - t274 * t577 + t176) * t336 + (-t311 * t578 - t313 * t577 - t635) * t459) * t603 + t612 * t600) * t1011 + (t43 * t604 + t44 * t601) * t574 / 0.2e1 - t47 * t842 / 0.2e1 - (t1077 * qJD(1) + t1067 * qJD(2) + qJDD(1) * t1111 + t1075 * t457 + t1076 * t458 + t5 + t9) * t604 / 0.2e1 + (t1086 * t604 + t1087 * t601 + (t1071 * t604 + t1072 * t601) * qJD(1)) * qJD(1) / 0.2e1 + (t43 + t36 + t1088) * t790 + (t37 + t44 + t1089) * t850 / 0.2e1 + ((-t1056 * t838 - t1059) * t604 + ((t1057 * t604 + t1069) * qJD(2) + t1068) * t601) * t784 + ((-t1057 * t576 + t1059) * t601 + ((t1056 * t601 + t1069) * qJD(2) + t1068) * t604) * t789 + ((t1075 * t604 + t1076 * t601) * qJD(1) + t1067) * t785 + (qJD(1) * t1078 + qJD(2) * t1066 + qJDD(1) * t1070 + t1073 * t457 + t1074 * t458 + t10 + t6) * t601 / 0.2e1 + (g(1) * t434 + g(2) * t429 - g(3) * t495 + (qJD(2) * t695 + t363 * t457 - t365 * t458) * t687 + t186 * ((t363 * t604 - t365 * t601) * qJD(1) + t695) + t699 * t451 + (-t106 * t604 - t107 * t601 + (-t193 * t604 + t938) * qJD(1)) * t491 - (t192 * t429 - t936) * qJD(1) - (t186 * (-t429 * t601 - t434 * t604) + t699 * t495) * qJD(2)) * m(3) + t1101 * t1012 + t1102 * t1013 - ((t1098 * t603 + t1103 * t600) * qJD(2) + (t1099 * t600 + t1100 * t603) * qJD(1)) * qJD(1) / 0.2e1 + (qJD(1) * t729 + t20 * t601 - t21 * t604) * t1010 + (qJD(1) * t725 + t24 * t601 - t25 * t604) * t1015 + (qJD(1) * t727 + t26 * t601 - t27 * t604) * t1016 + t722 * t1018 + (t128 * t437 + t53 * t880 + t125 * t816 + (t64 * t864 + t128 * t879 + t53 * t366 + t125 * t256 + (-t125 * t367 + t129 * t864) * qJD(1)) * t604 + (t65 * t864 + t129 * t879 - t53 * t367 + t125 * t255 + (t125 * t881 - t128 * t739) * qJD(1)) * t601 - g(1) * (t432 + t433) - g(2) * (t427 + t428) + g(3) * t863 - t128 * (-qJD(1) * t428 + t769) - t129 * (qJD(1) * t433 + t877) - t125 * t814 - ((t125 * t433 + t128 * t863) * t604 + (t125 * t428 + t129 * t863) * t601) * qJD(2)) * m(4) + (-t83 * (-t1042 * t294 - t234 * t842 - t362 * t443 + t750) - t84 * (t1042 * t295 + t232 * t842 - t362 * t442 + t877) - t76 * (t294 * t442 + t295 * t443 + t814) - (t747 * t83 + t748 * t84) * qJD(2) - ((-t76 * t770 - t817) * pkin(7) + t623 * qJD(4)) * t600 + t83 * t874 + t23 * t761 + t76 * t672 + (t38 * t754 + t83 * t676 + t23 * t232 + t76 * t118 + (t76 * t234 + t754 * t84) * qJD(1)) * t604 + (t39 * t754 + t84 * t676 + t23 * t234 + t76 * t117 + (t83 * t364 + t76 * (-t232 - t875)) * qJD(1)) * t601 - g(1) * (t537 + t865) - g(2) * (t535 + t866) - g(3) * (t362 - t780) - t833 * t1040) * m(5) + ((t288 * t408 - t290 * t409) * t442 - (t287 * t408 - t289 * t409) * t443 + (t341 * t408 - t345 * t409) * t1042 + (t122 * t603 - t81 * t920) * qJD(4) + ((-qJD(4) * t82 + t649) * t600 + t610) * t601) * t1014 + ((t288 * t406 + t290 * t407) * t442 - (t287 * t406 + t289 * t407) * t443 + (t341 * t406 + t345 * t407) * t1042 + (t121 * t603 - t80 * t921) * qJD(4) + ((-qJD(4) * t79 + t649) * t600 + t610) * t604) * t1017 + (qJD(1) * t723 + t28 * t601 - t29 * t604) * t1008 + (qJD(1) * t731 + t16 * t601 - t17 * t604) * t1020 + (qJD(1) * t733 + t18 * t601 - t19 * t604) * t1021 + t728 * t1023 + t724 * t1024 + t726 * t1025 + t730 * t1026 + t732 * t1027 + t401 * t37 / 0.2e1 + t402 * t36 / 0.2e1 - t470 * t41 / 0.2e1; (-m(4) - m(5) - m(6)) * (-g(3) * t603 + t1040) - m(4) * (t125 * t403 + t128 * t405 + t129 * t404) - m(5) * (t403 * t76 + t404 * t84 + t405 * t83) - m(6) * (t403 * t52 + t404 * t67 + t405 * t66) + 0.2e1 * ((t128 * t838 + t129 * t576 - t53) * t1030 + (t576 * t84 + t83 * t838 - t23) * t1029 + (t576 * t67 + t66 * t838 - t11) * t1028) * t603 + 0.2e1 * ((qJD(2) * t125 - t128 * t850 + t129 * t848 + t601 * t65 + t604 * t64) * t1030 + (qJD(2) * t76 + t38 * t604 + t39 * t601 - t83 * t850 + t817) * t1029 + (qJD(2) * t52 + t14 * t604 + t15 * t601 - t66 * t850 + t818) * t1028) * t600; ((-qJD(2) * t723 + t55) * t600 + (-qJD(1) * t722 + qJD(2) * t133 + t28 * t604 + t29 * t601) * t603) * t1008 + ((-qJD(2) * t725 + t49) * t600 + (-qJD(1) * t724 + qJD(2) * t122 + t24 * t604 + t25 * t601) * t603) * t1015 + ((-qJD(2) * t727 + t50) * t600 + (-qJD(1) * t726 + qJD(2) * t121 + t26 * t604 + t27 * t601) * t603) * t1016 + t625 + (t121 * t600 + t603 * t727) * t1025 + (t122 * t600 + t603 * t725) * t1024 + (t133 * t600 + t603 * t723) * t1018 + t10 * t795 + (t648 * t600 + (t627 * t599 - t602 * t626) * t603) * t1009 + t9 * t796 + t47 * t787 + (t408 * t626 + t409 * t627 + t601 * t644) * t1014 + (t967 + t973 + t974 - t981 + t982) * t1007 + (t626 * t406 - t407 * t627 + t604 * t644) * t1017 + t1064 * t44 + t1065 * t43 + (-g(1) * (t385 + t219) - g(2) * (t386 + t220) - g(3) * (t508 + (-t990 - t1001) * t603) - t66 * (-t1042 * t386 + t443 * t832 + t776) - t67 * (t1042 * t385 + t442 * t832 + t777) - t52 * (t385 * t443 + t386 * t442 + t778) + t14 * t282 + t66 * t910 + t52 * t966 + (-t14 * t903 - t66 * t969 + t15 * t258 + t67 * t124 + (t665 * t604 + (t52 * t258 - t1053) * t601) * qJD(2)) * t600 + ((t67 * t258 - t66 * t903) * qJD(2) + (-t15 * t891 + t67 * t909 + t11 * t259 + t52 * t123 + (t66 * t379 - t52 * t904) * qJD(1)) * t604 + (qJD(1) * t665 + t14 * t379 + t66 * t281 - t1033) * t601) * t603 + t631) * m(6) + (-t83 * (-t1042 * t269 - t431 * t443) - t84 * (t1042 * t268 - t431 * t442) - t76 * (t268 * t443 + t269 * t442) + (qJD(2) * t623 - t83 * t117 + t84 * t118 + t39 * t232 - t38 * t234) * t600 + (t83 * (-qJD(2) * t234 + t252 * t601) + t84 * (qJD(2) * t232 - t252 * t604) - t23 * t696 + t76 * (t117 * t604 - t118 * t601 - t232 * t848 - t234 * t850) + (t38 * t601 - t39 * t604 + (t601 * t84 + t604 * t83) * qJD(1)) * t364) * t603 - g(1) * t268 - g(2) * t269 - g(3) * t431) * m(5); t625 + (t14 * (-t185 * t600 + t282) + t52 * (-t185 * t803 + t966) + (-t185 * t941 + (-qJD(1) * t183 * t52 - t15 * t316 - t173 * t67) * t604 + (-t11 * t183 + t52 * (-qJD(1) * t185 - t96)) * t601) * t603 + t631 - t52 * t778 - t67 * t777 - g(1) * t219 - g(2) * t220 - g(3) * t376 + (-t316 * t804 - t600 * t95 - t776 + t910) * t66) * m(6);];
tau = t1;