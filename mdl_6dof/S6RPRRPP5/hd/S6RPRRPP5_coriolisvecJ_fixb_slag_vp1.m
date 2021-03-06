% Calculate vector of centrifugal and Coriolis load on the joints for
% S6RPRRPP5
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d3,d4,theta2]';
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
% Datum: 2019-03-09 04:45
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauc = S6RPRRPP5_coriolisvecJ_fixb_slag_vp1(qJ, qJD, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(9,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRRPP5_coriolisvecJ_fixb_slag_vp1: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RPRRPP5_coriolisvecJ_fixb_slag_vp1: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S6RPRRPP5_coriolisvecJ_fixb_slag_vp1: pkin has to be [9x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RPRRPP5_coriolisvecJ_fixb_slag_vp1: m has to be [7x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [7,3]), ...
  'S6RPRRPP5_coriolisvecJ_fixb_slag_vp1: rSges has to be [7x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [7 6]), ...
  'S6RPRRPP5_coriolisvecJ_fixb_slag_vp1: Icges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From coriolisvec_joint_fixb_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 04:42:43
% EndTime: 2019-03-09 04:44:33
% DurationCPUTime: 101.42s
% Computational Cost: add. (39790->1279), mult. (57867->1583), div. (0->0), fcn. (55102->8), ass. (0->619)
t1148 = Icges(7,4) + Icges(6,5);
t1090 = -Icges(6,1) - Icges(7,1);
t1083 = Icges(5,1) - t1090;
t1145 = Icges(7,2) + Icges(6,3);
t1109 = Icges(6,4) + Icges(5,5) - Icges(7,5);
t1108 = Icges(5,6) - Icges(6,6) + Icges(7,6);
t582 = cos(qJ(4));
t1147 = t1148 * t582;
t580 = sin(qJ(4));
t1146 = (Icges(5,4) - t1148) * t580;
t1144 = -t1145 * t580 - t1147;
t1082 = Icges(5,2) + t1145;
t1143 = Icges(6,2) + Icges(5,3) + Icges(7,3);
t1142 = t1108 * t580 - t1109 * t582;
t1141 = t1083 * t582 - t1146;
t576 = pkin(9) + qJ(3);
t563 = cos(t576);
t583 = cos(qJ(1));
t899 = t582 * t583;
t581 = sin(qJ(1));
t901 = t580 * t581;
t481 = t563 * t901 + t899;
t898 = t583 * t580;
t900 = t581 * t582;
t482 = t563 * t900 - t898;
t562 = sin(t576);
t906 = t562 * t581;
t240 = Icges(5,5) * t482 - Icges(5,6) * t481 + Icges(5,3) * t906;
t246 = Icges(6,4) * t482 + Icges(6,2) * t906 + Icges(6,6) * t481;
t1140 = t240 + t246;
t905 = t562 * t582;
t1139 = t1148 * t905;
t940 = Icges(5,4) * t582;
t702 = -Icges(5,2) * t580 + t940;
t1138 = -t1108 * t563 + t562 * t702;
t450 = Icges(6,5) * t482;
t238 = -Icges(6,6) * t906 - Icges(6,3) * t481 - t450;
t456 = Icges(5,4) * t482;
t249 = -Icges(5,2) * t481 + Icges(5,6) * t906 + t456;
t1137 = t238 + t249;
t449 = Icges(6,5) * t481;
t255 = Icges(6,1) * t482 + Icges(6,4) * t906 + t449;
t455 = Icges(5,4) * t481;
t259 = -Icges(5,1) * t482 - Icges(5,5) * t906 + t455;
t1136 = t255 - t259;
t1111 = t1142 * t562 + t1143 * t563;
t1135 = -t1142 * t563 + t1143 * t562;
t907 = t562 * t580;
t1107 = t1145 * t907 - t1138 + t1139;
t1047 = (-t702 - t1144) * t563 - t1108 * t562;
t1098 = -t1109 * t563 + t1141 * t562;
t1046 = t1109 * t562 + t1141 * t563;
t1134 = (t1082 * t582 + t1146) * t562;
t1133 = (t1108 * t582 + t1109 * t580) * t562;
t484 = t563 * t899 + t901;
t451 = Icges(6,5) * t484;
t483 = t563 * t898 - t900;
t904 = t562 * t583;
t239 = Icges(6,6) * t904 + Icges(6,3) * t483 + t451;
t942 = Icges(5,4) * t484;
t251 = -Icges(5,2) * t483 + Icges(5,6) * t904 + t942;
t1132 = t239 - t251;
t242 = Icges(5,5) * t484 - Icges(5,6) * t483 + Icges(5,3) * t904;
t248 = Icges(6,4) * t484 + Icges(6,2) * t904 + Icges(6,6) * t483;
t1131 = t242 + t248;
t933 = Icges(6,5) * t483;
t257 = Icges(6,1) * t484 + Icges(6,4) * t904 + t933;
t457 = Icges(5,4) * t483;
t260 = Icges(5,1) * t484 + Icges(5,5) * t904 - t457;
t1130 = t257 + t260;
t1129 = t1140 * t906;
t1128 = t1140 * t904;
t1127 = t1136 * t482 - t1137 * t481 + t1129;
t1126 = -t1136 * t484 + t1137 * t483 - t1128;
t1094 = t1130 * t482 + t1131 * t906 + t1132 * t481;
t236 = Icges(7,5) * t484 + Icges(7,6) * t483 - Icges(7,3) * t904;
t454 = Icges(7,4) * t484;
t245 = Icges(7,2) * t483 - Icges(7,6) * t904 + t454;
t938 = Icges(7,4) * t483;
t254 = Icges(7,1) * t484 - Icges(7,5) * t904 + t938;
t893 = -t481 * t245 - t482 * t254;
t72 = -t236 * t906 - t893;
t1106 = t72 + t1094;
t1093 = t1130 * t484 + t1131 * t904 + t1132 * t483;
t891 = t483 * t245 + t484 * t254;
t78 = -t236 * t904 + t891;
t1105 = t78 + t1093;
t829 = qJD(4) * t583;
t788 = t563 * t829;
t834 = qJD(3) * t583;
t794 = t562 * t834;
t831 = qJD(4) * t581;
t215 = qJD(1) * t481 - t582 * t788 + (t794 - t831) * t580;
t832 = qJD(4) * t563;
t541 = qJD(1) - t832;
t681 = t541 * t580;
t759 = qJD(1) * t563 - qJD(4);
t216 = t583 * t681 + (-t581 * t759 - t794) * t582;
t791 = t563 * t834;
t840 = qJD(1) * t581;
t799 = t562 * t840;
t641 = -t791 + t799;
t109 = Icges(7,5) * t216 - Icges(7,6) * t215 + Icges(7,3) * t641;
t113 = Icges(5,5) * t216 + Icges(5,6) * t215 - Icges(5,3) * t641;
t117 = Icges(6,4) * t216 - Icges(6,2) * t641 - Icges(6,6) * t215;
t1125 = -t109 + t113 + t117;
t836 = qJD(3) * t581;
t795 = t562 * t836;
t830 = qJD(4) * t582;
t839 = qJD(1) * t583;
t903 = t563 * t581;
t217 = t830 * t903 - t582 * t840 + (t563 * t839 - t795 - t829) * t580;
t835 = qJD(3) * t582;
t218 = t759 * t899 + (-t562 * t835 + t681) * t581;
t642 = t562 * t839 + t563 * t836;
t110 = Icges(7,5) * t218 + Icges(7,6) * t217 - Icges(7,3) * t642;
t114 = Icges(5,5) * t218 - Icges(5,6) * t217 + Icges(5,3) * t642;
t118 = Icges(6,4) * t218 + Icges(6,2) * t642 + Icges(6,6) * t217;
t1124 = -t110 + t114 + t118;
t111 = Icges(6,5) * t216 - Icges(6,6) * t641 - Icges(6,3) * t215;
t115 = Icges(7,4) * t216 - Icges(7,2) * t215 + Icges(7,6) * t641;
t119 = Icges(5,4) * t216 + Icges(5,2) * t215 - Icges(5,6) * t641;
t1123 = t111 + t115 - t119;
t112 = Icges(6,5) * t218 + Icges(6,6) * t642 + Icges(6,3) * t217;
t116 = Icges(7,4) * t218 + Icges(7,2) * t217 - Icges(7,6) * t642;
t120 = Icges(5,4) * t218 - Icges(5,2) * t217 + Icges(5,6) * t642;
t1122 = t112 + t116 - t120;
t121 = Icges(7,1) * t216 - Icges(7,4) * t215 + Icges(7,5) * t641;
t123 = Icges(6,1) * t216 - Icges(6,4) * t641 - Icges(6,5) * t215;
t125 = Icges(5,1) * t216 + Icges(5,4) * t215 - Icges(5,5) * t641;
t1121 = t121 + t123 + t125;
t122 = Icges(7,1) * t218 + Icges(7,4) * t217 - Icges(7,5) * t642;
t124 = Icges(6,1) * t218 + Icges(6,4) * t642 + Icges(6,5) * t217;
t126 = Icges(5,1) * t218 - Icges(5,4) * t217 + Icges(5,5) * t642;
t1120 = t122 + t124 + t126;
t1016 = t1098 * t484 + t1107 * t483 - t1111 * t904;
t1119 = qJD(3) * t1135 - qJD(4) * t1133;
t1118 = qJD(3) * t1047 + qJD(4) * t1134;
t440 = (-Icges(5,1) * t580 - t940) * t562;
t833 = qJD(4) * t562;
t1117 = qJD(4) * t440 + (t1090 * t580 + t1147) * t833 + t1046 * qJD(3);
t235 = -Icges(7,5) * t482 - Icges(7,6) * t481 + Icges(7,3) * t906;
t1116 = t235 + t1140;
t1115 = t236 - t1131;
t453 = Icges(7,4) * t482;
t244 = -Icges(7,2) * t481 + Icges(7,6) * t906 - t453;
t1081 = t244 + t1137;
t1114 = t245 + t1132;
t452 = Icges(7,4) * t481;
t252 = Icges(7,1) * t482 - Icges(7,5) * t906 + t452;
t1113 = t252 + t1136;
t1112 = t254 + t1130;
t1017 = -t1098 * t482 - t1107 * t481 + t1111 * t906;
t1110 = t1098 * t582 + t1107 * t580;
t1092 = rSges(7,3) + qJ(6);
t1060 = t1081 * t215 + t1113 * t216 - t1116 * t641 + t1120 * t484 + t1122 * t483 + t1124 * t904;
t1059 = t1112 * t216 - t1114 * t215 + t1115 * t641 + t1121 * t484 + t1123 * t483 + t1125 * t904;
t1058 = -t1081 * t217 + t1113 * t218 + t1116 * t642 + t1120 * t482 + t1122 * t481 + t1124 * t906;
t1057 = t1112 * t218 + t1114 * t217 - t1115 * t642 + t1121 * t482 + t1123 * t481 + t1125 * t906;
t1088 = t1098 * t216 - t1107 * t215 + t1111 * t641 + t1117 * t484 + t1118 * t483 + t1119 * t904;
t1087 = t1098 * t218 + t1107 * t217 - t1111 * t642 + t1117 * t482 + t1118 * t481 + t1119 * t906;
t693 = -t244 * t481 + t252 * t482;
t71 = t235 * t906 + t693;
t1085 = t71 + t1127;
t892 = t244 * t483 - t484 * t252;
t77 = t235 * t904 - t892;
t1084 = t77 - t1126;
t1069 = t563 * t235;
t692 = -t244 * t580 + t252 * t582;
t92 = t562 * t692 - t1069;
t1072 = t246 * t563;
t695 = -t238 * t580 + t255 * t582;
t94 = t562 * t695 - t1072;
t1075 = t240 * t563;
t690 = -t249 * t580 - t259 * t582;
t96 = t562 * t690 - t1075;
t1104 = t92 + t94 + t96;
t691 = t245 * t580 + t254 * t582;
t93 = t236 * t563 + t562 * t691;
t694 = t239 * t580 + t257 * t582;
t95 = -t248 * t563 + t562 * t694;
t689 = -t251 * t580 + t260 * t582;
t97 = -t242 * t563 + t562 * t689;
t1103 = t93 + t95 + t97;
t1015 = t1110 * t562 + t1111 * t563;
t1099 = t1144 * t562 + t1138;
t493 = t562 * t829 + t836;
t494 = -t562 * t831 + t834;
t1097 = (-t1110 + t1135) * t541 + (-t1111 * t581 + t690 + t692 + t695) * t494 + (t1111 * t583 - t689 - t691 - t694) * t493;
t1096 = t1016 * t541 + t1126 * t494;
t1095 = t1017 * t541 + t1127 * t494;
t1076 = rSges(7,1) + pkin(5);
t1091 = -pkin(4) - t1076;
t983 = pkin(3) * t563;
t504 = pkin(8) * t562 + t983;
t476 = t504 * t581;
t579 = -pkin(7) - qJ(2);
t556 = t583 * t579;
t731 = rSges(7,1) * t582 + rSges(7,2) * t580;
t1067 = -pkin(5) * t905 - t1092 * t563 - t562 * t731;
t446 = qJD(5) * t483;
t133 = t216 * pkin(4) - qJ(5) * t215 + t446;
t448 = t481 * qJ(5);
t316 = pkin(4) * t482 + t448;
t1021 = qJD(1) * t476;
t567 = t583 * qJ(2);
t524 = pkin(1) * t581 - t567;
t578 = cos(pkin(9));
t553 = pkin(2) * t578 + pkin(1);
t850 = -t581 * t553 - t556;
t398 = t524 + t850;
t372 = qJD(1) * t398;
t503 = pkin(3) * t562 - pkin(8) * t563;
t796 = t503 * t834;
t508 = qJD(1) * t524;
t564 = qJD(2) * t581;
t845 = t564 - t508;
t632 = t372 - t796 + t845 - t1021;
t518 = pkin(8) * t791;
t851 = t518 + t564;
t1089 = t316 * t541 + t133 - t446 - t632 + t851;
t565 = qJD(2) * t583;
t566 = t581 * qJ(2);
t844 = t583 * pkin(1) + t566;
t1038 = qJD(1) * t844 - t565;
t977 = rSges(3,2) * sin(pkin(9));
t979 = rSges(3,1) * t578;
t680 = t581 * rSges(3,3) + (-t977 + t979) * t583;
t324 = qJD(1) * t680 + t1038;
t1054 = t1106 * t493 - t494 * t71 - t1095;
t1053 = t1105 * t493 - t494 * t77 + t1096;
t1086 = (qJD(3) * t1110 - t1119) * t563 + (t1117 * t582 + t1118 * t580 + (-t1098 * t580 + t1107 * t582) * qJD(4) - t1111 * qJD(3)) * t562;
t1079 = t1133 * t541 + (-t1108 * t482 - t1109 * t481) * t494 + (t1108 * t484 + t1109 * t483) * t493;
t1078 = -(t235 * t581 + t236 * t583) * t562 + t891;
t730 = pkin(4) * t582 + qJ(5) * t580;
t1035 = t730 * t562;
t322 = t484 * pkin(4) + qJ(5) * t483;
t902 = t563 * t583;
t543 = pkin(3) * t902;
t478 = pkin(8) * t904 + t543;
t535 = t583 * t553;
t762 = -t579 * t581 + t535;
t805 = qJD(1) * (t762 - t844) + t1038;
t709 = qJD(1) * t478 - t503 * t836 + t805;
t828 = qJD(5) * t481;
t1065 = t1035 * t493 - t541 * t322 - t709 - t828;
t826 = qJD(6) * t562;
t519 = t581 * t826;
t1027 = t483 * rSges(7,2) + t1076 * t484;
t880 = -t1092 * t904 + t1027;
t53 = t1067 * t493 + t541 * t880 - t1065 - t519;
t1077 = 0.2e1 * qJD(3);
t1066 = -rSges(7,2) * t481 + t1092 * t906;
t409 = t1035 * t840;
t464 = t730 * t563;
t480 = t503 * t840;
t1064 = -t1035 * t563 * t831 + t494 * t464 + t409 + t480;
t822 = qJD(3) * qJD(4);
t780 = t563 * t822;
t363 = qJD(1) * t493 + t581 * t780;
t364 = qJD(1) * t494 + t583 * t780;
t781 = t562 * t822;
t1063 = t1016 * t781 + t1059 * t493 - t1060 * t494 + t1084 * t363 + t1088 * t541 + t1105 * t364;
t1062 = -t1017 * t781 + t1057 * t493 - t1058 * t494 + t1085 * t363 + t1087 * t541 + t1106 * t364;
t25 = (qJD(3) * t692 + t110) * t563 + (qJD(3) * t235 + t116 * t580 + t122 * t582 + (-t244 * t582 - t252 * t580) * qJD(4)) * t562;
t27 = (qJD(3) * t695 - t118) * t563 + (qJD(3) * t246 + t112 * t580 + t124 * t582 + (-t238 * t582 - t255 * t580) * qJD(4)) * t562;
t29 = (qJD(3) * t690 - t114) * t563 + (qJD(3) * t240 - t120 * t580 + t126 * t582 + (-t249 * t582 + t259 * t580) * qJD(4)) * t562;
t1056 = t25 + t27 + t29;
t26 = (qJD(3) * t691 + t109) * t563 + (-qJD(3) * t236 + t115 * t580 + t121 * t582 + (t245 * t582 - t254 * t580) * qJD(4)) * t562;
t28 = (qJD(3) * t694 - t117) * t563 + (qJD(3) * t248 + t111 * t580 + t123 * t582 + (t239 * t582 - t257 * t580) * qJD(4)) * t562;
t30 = (qJD(3) * t689 - t113) * t563 + (qJD(3) * t242 - t119 * t580 + t125 * t582 + (-t251 * t582 - t260 * t580) * qJD(4)) * t562;
t1055 = t26 + t28 + t30;
t1052 = t1015 * t541 + t1103 * t493 - t1104 * t494;
t1051 = t1099 * t581;
t1050 = t1099 * t583;
t1049 = t1098 * t581;
t1048 = t1098 * t583;
t1045 = t1097 * t562;
t1044 = -t1111 * t541 - t1115 * t493 - t1116 * t494;
t1010 = t1103 * t583 + t1104 * t581;
t1043 = t1103 * t581 - t1104 * t583;
t1009 = t1084 * t581 + t1105 * t583;
t1042 = -t1084 * t583 + t1105 * t581;
t1008 = t1085 * t581 + t1106 * t583;
t1041 = -t1085 * t583 + t1106 * t581;
t736 = rSges(5,1) * t482 - rSges(5,2) * t481;
t269 = rSges(5,3) * t906 + t736;
t735 = rSges(5,1) * t582 - rSges(5,2) * t580;
t394 = -rSges(5,3) * t563 + t562 * t735;
t1039 = -t269 * t541 - t394 * t494;
t1037 = t893 + (-t235 * t583 + t236 * t581) * t562;
t672 = t481 * t541 + t494 * t907;
t1032 = -t215 + t672;
t671 = -t483 * t541 + t493 * t907;
t1031 = t217 + t671;
t883 = t1076 * t482 - t1066;
t402 = t583 * t1035;
t1030 = t322 * t833 - t541 * t402;
t1028 = -t1035 * t494 + t446;
t492 = qJD(3) * t504;
t551 = qJD(6) * t563;
t761 = -t492 - t551;
t1026 = -rSges(7,2) * t217 + t519;
t552 = Icges(4,4) * t563;
t703 = -Icges(4,2) * t562 + t552;
t499 = Icges(4,1) * t562 + t552;
t1024 = -t215 * rSges(7,2) + t1076 * t216 + t1092 * t799;
t274 = t484 * rSges(6,1) + rSges(6,2) * t904 + t483 * rSges(6,3);
t732 = rSges(6,1) * t582 + rSges(6,3) * t580;
t393 = -rSges(6,2) * t563 + t562 * t732;
t68 = t274 * t541 - t393 * t493 - t1065;
t827 = qJD(5) * t580;
t520 = t562 * t827;
t837 = qJD(3) * t563;
t639 = t562 * t830 + t580 * t837;
t640 = t563 * t835 - t580 * t833;
t205 = pkin(4) * t640 + qJ(5) * t639 + t520;
t1014 = -qJD(5) * t215 + t1035 * t363 - t494 * t205;
t1013 = (-t1098 + t1134) * t541 + (-t1082 * t482 + t1113 + t449 + t452 - t455) * t494 + (t1082 * t484 - t1112 + t457 - t933 - t938) * t493;
t1012 = (t1090 * t907 + t1107 + t1139 + t440) * t541 + (t1083 * t481 + t1081 - t450 - t453 + t456) * t494 + (-t1083 * t483 + t1114 + t451 + t454 - t942) * t493;
t1011 = t1079 * t562;
t1007 = t1015 * t781 + t1086 * t541;
t923 = Icges(4,3) * t583;
t403 = Icges(4,5) * t903 - Icges(4,6) * t906 - t923;
t533 = Icges(4,4) * t906;
t935 = Icges(4,5) * t583;
t407 = Icges(4,1) * t903 - t533 - t935;
t927 = Icges(4,6) * t583;
t405 = Icges(4,4) * t903 - Icges(4,2) * t906 - t927;
t912 = t405 * t562;
t684 = -t407 * t563 + t912;
t149 = -t403 * t583 - t581 * t684;
t496 = Icges(4,5) * t563 - Icges(4,6) * t562;
t495 = Icges(4,5) * t562 + Icges(4,6) * t563;
t651 = qJD(3) * t495;
t943 = Icges(4,4) * t562;
t500 = Icges(4,1) * t563 - t943;
t408 = Icges(4,5) * t581 + t500 * t583;
t406 = Icges(4,6) * t581 + t583 * t703;
t911 = t406 * t562;
t683 = -t408 * t563 + t911;
t1006 = -t583 * t651 + (-t496 * t581 + t683 + t923) * qJD(1);
t404 = Icges(4,3) * t581 + t496 * t583;
t843 = qJD(1) * t404;
t1005 = qJD(1) * t684 - t581 * t651 + t843;
t497 = Icges(4,2) * t563 + t943;
t682 = t497 * t562 - t499 * t563;
t1004 = qJD(1) * t682 + t496 * qJD(3);
t889 = t217 * qJ(5) + t828;
t982 = pkin(4) * t218;
t134 = t889 + t982;
t521 = qJD(5) * t905;
t522 = t563 * t827;
t1003 = qJD(3) * t522 + qJD(4) * t521 + t493 * t134 + t364 * t316;
t1002 = t581 * (-t497 * t583 + t408) - t583 * (-Icges(4,2) * t903 + t407 - t533);
t998 = -m(7) / 0.2e1;
t997 = m(7) / 0.2e1;
t996 = t363 / 0.2e1;
t995 = t364 / 0.2e1;
t994 = -t493 / 0.2e1;
t993 = t493 / 0.2e1;
t992 = -t494 / 0.2e1;
t991 = t494 / 0.2e1;
t990 = -t541 / 0.2e1;
t989 = t541 / 0.2e1;
t987 = t581 / 0.2e1;
t986 = -t583 / 0.2e1;
t985 = -rSges(6,2) - pkin(8);
t984 = -rSges(5,3) - pkin(8);
t980 = pkin(1) - t553;
t978 = rSges(4,1) * t563;
t976 = rSges(6,2) * t562;
t973 = rSges(5,3) * t562;
t159 = t481 * t493 + t483 * t494;
t862 = t476 * t836 + t478 * t834;
t756 = t493 * t316 + t520 + t862;
t810 = t322 + t880;
t49 = t493 * t883 + t494 * t810 + t551 + t756;
t970 = t159 * t49;
t969 = t25 * t494;
t968 = t26 * t493;
t967 = t27 * t494;
t966 = t28 * t493;
t965 = t29 * t494;
t964 = t30 * t493;
t816 = t216 * rSges(5,1) + t215 * rSges(5,2) + rSges(5,3) * t791;
t129 = -rSges(5,3) * t799 + t816;
t397 = t563 * t735 + t973;
t461 = (-rSges(5,1) * t580 - rSges(5,2) * t582) * t562;
t223 = qJD(3) * t397 + qJD(4) * t461;
t275 = t484 * rSges(5,1) - t483 * rSges(5,2) + rSges(5,3) * t904;
t800 = t503 * t839;
t660 = -t492 * t581 - t800;
t643 = -t563 * t840 - t794;
t289 = pkin(3) * t643 - pkin(8) * t799 + t518;
t557 = qJ(2) * t839;
t824 = qJD(1) * qJD(2);
t846 = t557 + t564;
t856 = qJD(1) * (-pkin(1) * t840 + t846) + t581 * t824;
t806 = qJD(1) * (-t557 + (t581 * t980 - t556) * qJD(1)) + t856;
t754 = qJD(1) * t289 + t806;
t41 = t129 * t541 - t223 * t493 - t364 * t394 + (t275 * t833 + t660) * qJD(3) + t754;
t963 = t41 * t583;
t737 = rSges(5,1) * t218 - rSges(5,2) * t217;
t132 = rSges(5,3) * t642 + t737;
t517 = pkin(3) * t795;
t290 = pkin(8) * t642 + qJD(1) * t543 - t517;
t545 = t579 * t840;
t872 = t545 - (-t583 * t980 - t566) * qJD(1) - t1038;
t807 = -t290 + t872;
t555 = t583 * t824;
t823 = qJD(1) * qJD(3);
t782 = t581 * t823;
t859 = t503 * t782 + t555;
t644 = qJD(1) * t807 + t859;
t910 = t492 * t583;
t42 = -t132 * t541 - t223 * t494 + t363 * t394 + (-t269 * t833 - t910) * qJD(3) + t644;
t962 = t42 * t581;
t571 = t581 * rSges(4,3);
t88 = t275 * t541 - t394 * t493 + t709;
t957 = t581 * t88;
t733 = rSges(6,1) * t482 + rSges(6,3) * t481;
t268 = rSges(6,2) * t906 + t733;
t879 = t274 + t322;
t64 = t268 * t493 + t494 * t879 + t756;
t956 = t64 * t268;
t867 = t398 - t524;
t661 = t564 + (-t476 + t867) * qJD(1);
t605 = t661 - t796;
t882 = -t268 - t316;
t67 = -t393 * t494 + t541 * t882 + t1028 + t605;
t955 = t67 * t393;
t954 = t92 * t363;
t953 = t93 * t364;
t952 = t94 * t363;
t951 = t95 * t364;
t950 = t96 * t363;
t949 = t97 * t364;
t410 = rSges(4,1) * t903 - rSges(4,2) * t906 - t583 * rSges(4,3);
t501 = rSges(4,1) * t562 + rSges(4,2) * t563;
t789 = t501 * t834;
t740 = t564 - t789;
t153 = (-t410 + t867) * qJD(1) + t740;
t918 = t153 * t581;
t411 = rSges(4,1) * t902 - rSges(4,2) * t904 + t571;
t797 = t501 * t836;
t154 = qJD(1) * t411 - t797 + t805;
t465 = t501 * t583;
t917 = t154 * t465;
t909 = t495 * t581;
t908 = t495 * t583;
t897 = t134 * t904 + t316 * t791;
t896 = -rSges(7,3) * t791 + (-qJ(6) * t837 - t826) * t583 + t1024;
t814 = t216 * rSges(6,1) + rSges(6,2) * t791 - t215 * rSges(6,3);
t128 = -rSges(6,2) * t799 + t814;
t895 = t128 + t133;
t894 = t1076 * t218 - t1092 * t642 - t1026;
t312 = -pkin(4) * t481 + qJ(5) * t482;
t890 = t493 * t312 + t521;
t396 = t563 * t732 + t976;
t460 = (-rSges(6,1) * t580 + rSges(6,3) * t582) * t562;
t222 = qJD(3) * t396 + qJD(4) * t460;
t887 = -t205 - t222;
t395 = -rSges(7,3) * t562 + t563 * t731;
t459 = (-rSges(7,1) * t580 + rSges(7,2) * t582) * t562;
t838 = qJD(3) * t562;
t886 = pkin(5) * t640 - qJ(6) * t838 + qJD(3) * t395 + qJD(4) * t459 + t551;
t885 = -t223 - t492;
t318 = -pkin(4) * t483 + qJ(5) * t484;
t884 = qJD(5) * t482 + t541 * t318;
t878 = t322 * t838 + t562 * t409;
t877 = t1035 * t906 + t563 * t316;
t876 = t1067 * t581;
t875 = t1067 * t583;
t874 = -t581 * t403 - t407 * t902;
t873 = t581 * t404 + t408 * t902;
t870 = -t393 - t1035;
t869 = -t394 - t503;
t868 = -pkin(5) * t563 * t582 + qJ(6) * t562 - t395;
t475 = t503 * t581;
t477 = t503 * t583;
t863 = -t475 * t836 - t477 * t834;
t860 = t581 * t476 + t583 * t478;
t458 = (-pkin(4) * t580 + qJ(5) * t582) * t562;
t858 = qJD(5) * t484 - t494 * t458;
t855 = -t497 + t500;
t854 = t499 + t703;
t852 = rSges(4,2) * t799 + rSges(4,3) * t839;
t546 = t581 * t977;
t849 = rSges(3,3) * t839 + qJD(1) * t546;
t848 = t545 + t565;
t847 = t583 * rSges(3,3) + t546;
t842 = qJD(1) * t475;
t841 = qJD(1) * t496;
t179 = -t581 * t682 - t908;
t825 = t179 * qJD(1);
t820 = t581 * t979;
t818 = t133 + t896;
t817 = t494 * t316;
t813 = -t205 - t886;
t812 = -t492 + t887;
t811 = -t316 - t883;
t809 = t290 * t836 + (t1021 + t289) * t834;
t808 = t583 * t289 + t581 * t290 + t476 * t839;
t804 = -t1035 + t1067;
t803 = -t503 + t870;
t802 = t517 + t848;
t801 = t535 + t478;
t787 = t583 * t826;
t783 = -pkin(1) - t979;
t778 = t839 / 0.2e1;
t777 = t838 / 0.2e1;
t776 = -t836 / 0.2e1;
t775 = t836 / 0.2e1;
t773 = t834 / 0.2e1;
t770 = -t553 - t983;
t769 = t49 * t883;
t52 = (-qJD(3) * t503 - t826) * t583 + t1067 * t494 + t811 * t541 + t661 + t1028;
t768 = t52 * t1067;
t767 = t68 * t870;
t356 = t408 * t903;
t764 = t404 * t583 - t356;
t763 = -t403 + t911;
t760 = pkin(3) * t794;
t758 = t52 * (-pkin(8) + t1092);
t757 = t1035 * t642 + t563 * t134 + t205 * t906;
t755 = -t492 + t813;
t753 = t581 * t316 + t583 * t322 + t860;
t752 = -t503 + t804;
t745 = t53 * t804;
t742 = qJD(4) * t777;
t424 = qJD(1) * t477;
t741 = -t504 * t836 - t424;
t738 = -rSges(4,2) * t562 + t978;
t734 = rSges(6,1) * t218 + rSges(6,3) * t217;
t131 = rSges(6,2) * t642 + t734;
t630 = -t478 * t782 + t809;
t10 = t131 * t493 + t268 * t364 - t363 * t879 + t494 * t895 + t1003 + t630;
t729 = t10 * t268 + t64 * t131;
t401 = t581 * t1035;
t710 = t316 * t788 - t493 * t401 + t522 + t863;
t708 = t802 - t889;
t696 = -t153 * t583 - t154 * t581;
t688 = t269 * t583 - t275 * t581;
t182 = t405 * t563 + t407 * t562;
t183 = t406 * t563 + t408 * t562;
t679 = t583 * t133 + t581 * t134 + t316 * t839 + t808;
t463 = t501 * t581;
t658 = -t64 * t879 + t955;
t657 = t68 * t274 + t67 * t882;
t656 = -t504 * t834 + t842;
t654 = t322 + t801;
t653 = qJD(3) * t499;
t652 = qJD(3) * t497;
t150 = -t406 * t906 - t764;
t650 = (-t149 * t583 + t150 * t581) * qJD(3);
t151 = -t405 * t904 - t874;
t152 = -t406 * t904 + t873;
t649 = (-t151 * t583 + t152 * t581) * qJD(3);
t184 = (t410 * t581 + t411 * t583) * qJD(3);
t648 = qJD(5) * t217 + t541 * t133 + t322 * t781 + t754;
t646 = -t504 - t976;
t645 = -t504 - t973;
t7 = t894 * t493 + t883 * t364 + (-t478 * t840 - t826) * qJD(3) + t818 * t494 - t810 * t363 + t809 + t1003;
t631 = t49 * t894 + t7 * t883;
t628 = t405 * t583 - t406 * t581;
t625 = -t49 * t810 - t768;
t624 = t52 * t811 + t53 * t880;
t623 = -t520 + t761;
t606 = (-t562 * t854 + t563 * t855) * qJD(1);
t87 = t1039 + t605;
t89 = t269 * t493 + t275 * t494 + t862;
t594 = t89 * t688 + (t581 * t87 - t583 * t88) * t394;
t229 = qJD(1) * t406 - t581 * t652;
t231 = qJD(1) * t408 - t581 * t653;
t593 = qJD(1) * t403 - qJD(3) * t182 - t229 * t562 + t231 * t563;
t228 = -t583 * t652 + (-t581 * t703 + t927) * qJD(1);
t230 = -t583 * t653 + (-t500 * t581 + t935) * qJD(1);
t592 = -qJD(3) * t183 - t228 * t562 + t230 * t563 + t843;
t486 = t703 * qJD(3);
t487 = t500 * qJD(3);
t591 = qJD(1) * t495 - t486 * t562 + t487 * t563 + (-t497 * t563 - t499 * t562) * qJD(3);
t590 = (t767 + t956) * t583 + t658 * t581;
t589 = -t1002 * t562 + t628 * t563;
t588 = (t745 + t769) * t583 + t625 * t581;
t489 = t738 * qJD(3);
t426 = t820 - t847;
t354 = t394 * t583;
t353 = t393 * t583;
t351 = t394 * t581;
t350 = t393 * t581;
t323 = t564 + (-t426 - t524) * qJD(1);
t321 = -rSges(5,1) * t483 - rSges(5,2) * t484;
t320 = -rSges(6,1) * t483 + rSges(6,3) * t484;
t319 = -rSges(7,1) * t483 + rSges(7,2) * t484;
t315 = -rSges(5,1) * t481 - rSges(5,2) * t482;
t314 = -rSges(6,1) * t481 + rSges(6,3) * t482;
t313 = -rSges(7,1) * t481 + rSges(7,2) * t482;
t277 = t316 * t904;
t233 = -qJD(3) * t463 + (t583 * t738 + t571) * qJD(1);
t232 = rSges(4,1) * t643 - rSges(4,2) * t791 + t852;
t220 = -qJD(1) * t324 + t555;
t219 = qJD(1) * (-qJD(1) * t820 + t849) + t856;
t180 = -t583 * t682 + t909;
t178 = t180 * qJD(1);
t99 = -t489 * t834 + t555 + (-t233 + t797 + t872) * qJD(1);
t98 = -t489 * t836 + (t232 - t789) * qJD(1) + t806;
t91 = -qJD(3) * t683 + t228 * t563 + t230 * t562;
t90 = -qJD(3) * t684 + t229 * t563 + t231 * t562;
t86 = -t1004 * t583 + t591 * t581;
t85 = t1004 * t581 + t591 * t583;
t70 = t178 + t649;
t69 = t650 + t825;
t31 = t129 * t494 + t132 * t493 + t269 * t364 - t275 * t363 + t630;
t12 = -t222 * t494 + t363 * t393 + (-t131 - t134) * t541 + (t833 * t882 - t910) * qJD(3) + t644 + t1014;
t11 = t128 * t541 + t887 * t493 + t870 * t364 + (t274 * t833 + t660) * qJD(3) + t648;
t9 = -t886 * t494 - t1067 * t363 + (-t134 - t894) * t541 + (t519 + t807) * qJD(1) + (t583 * t761 + t811 * t833) * qJD(3) + t859 + t1014;
t8 = -qJD(1) * t787 + t896 * t541 + t813 * t493 + t804 * t364 + (t581 * t761 + t833 * t880 - t800) * qJD(3) + t648;
t1 = [-(t90 + t86 + t70) * t834 / 0.2e1 + (t220 * (t581 * t783 + t567 + t847) + t323 * t565 + t219 * (t680 + t844) + t324 * (t846 + t849) + (t323 * (t783 + t977) * t583 + (t323 * (-rSges(3,3) - qJ(2)) + t324 * t783) * t581) * qJD(1) - (-qJD(1) * t426 - t323 + t845) * t324) * m(3) + (t85 + t91) * t775 + (-t49 * t817 - (-t316 * t49 + t745) * t494 + t8 * (t654 + t1027) + (-t8 * t579 + t758 * t837) * t581 + (t1091 * t482 + t1066 - t448 - t476 + t850) * t9 + (t758 * qJD(1) - t1092 * t8) * t904 + (t1091 * t218 + t770 * t839 + t1026 + t53 + t708) * t52 + (t787 + t883 * t541 + t1024 + (-t504 - t553) * t840 + (-t1092 * t837 + (-pkin(3) * qJD(3) - qJD(6)) * t562 - t579 * qJD(1)) * t583 + t1089) * t53) * m(7) + (t42 * (-t736 + t850) + t87 * (-t737 + t802) + t41 * (t801 + t275) + t88 * (-t760 + t816 + t851) + (t837 * t87 * t984 - t41 * t579 + t42 * t645) * t581 + ((t562 * t984 + t770) * t957 + (t87 * (-t553 + t645) - t88 * t579) * t583) * qJD(1) - (t1039 + t632 - t87) * t88) * m(5) + ((t72 + t892 + t1037) * t494 + (t1081 * t481 - t1113 * t482 + t1078 + t1085 + t1093 - t1129) * t493 + t1096) * t991 + (t1072 + t1075 + t1069 + (t1081 * t580 - t1113 * t582) * t562 + t1104) * t493 * t990 + t1007 + ((-t1078 + t693 + t78) * t494 + (t1081 * t483 - t1113 * t484 + t1037 + t1084 - t1094 - t1128) * t493 + t1054 + t1095) * t994 + (t178 + ((t150 - t356 + (t404 + t912) * t583 + t874) * t583 + t873 * t581) * qJD(3)) * t773 + (-(-qJD(1) * t410 - t153 + t372 - t508 + t740) * t154 + t99 * (-t410 + t850) + t153 * t848 + t98 * (t411 + t762) + t154 * (t564 + t852) + (t501 * t918 - t917) * qJD(3) + ((-t153 * rSges(4,3) + t154 * (-t553 - t978)) * t581 + (t153 * (-t553 - t738) - t154 * t579) * t583) * qJD(1)) * m(4) + ((t182 + t179) * t581 + (t183 + t180) * t583) * t823 / 0.2e1 - t1017 * t996 + t1016 * t995 + t1088 * t993 + (t1087 + t1053) * t992 - t969 / 0.2e1 + (-qJD(3) * t682 + t486 * t563 + t487 * t562) * qJD(1) + t953 / 0.2e1 + t954 / 0.2e1 + t951 / 0.2e1 + t952 / 0.2e1 + t949 / 0.2e1 + t950 / 0.2e1 - t967 / 0.2e1 + t968 / 0.2e1 - t965 / 0.2e1 + t966 / 0.2e1 + t964 / 0.2e1 + (t12 * (-t316 - t733 + t850) + t11 * (t654 + t274) + (-t11 * t579 + t12 * t646) * t581 - t64 * t817 - (-t316 * t64 + t767) * t494 + (t708 - t734 - t982 + t837 * t985 * t581 + (-t553 + t646) * t839) * t67 + (-t760 + t814 + ((t562 * t985 + t770) * t581 - t556) * qJD(1) + t268 * t541 + t67 + t1089) * t68) * m(6) + (-t825 + ((t583 * t763 + t152 - t873) * t583 + (t581 * t763 + t151 + t764) * t581) * qJD(3) + t69) * t776; 0.2e1 * (t8 * t986 + t9 * t987) * m(7) + 0.2e1 * (t11 * t986 + t12 * t987) * m(6) + 0.2e1 * (-t963 / 0.2e1 + t962 / 0.2e1) * m(5) + 0.2e1 * (t98 * t986 + t987 * t99) * m(4) + 0.2e1 * (t219 * t986 + t220 * t987) * m(3); -qJD(1) * ((t562 * t855 + t563 * t854) * qJD(1) + (t1002 * t563 + t628 * t562) * qJD(3)) / 0.2e1 + ((-t834 * t909 - t841) * t583 + (t606 + (t583 * t908 + t589) * qJD(3)) * t581) * t773 + ((-t836 * t908 + t841) * t581 + (t606 + (t581 * t909 + t589) * qJD(3)) * t583) * t776 + qJD(1) * (t581 * t91 - t583 * t90 + (t182 * t581 + t183 * t583) * qJD(1)) / 0.2e1 + (-(t562 * t624 + t563 * t588) * qJD(4) - (t710 - t826 + (-t402 + t875) * t494 + t876 * t493) * t49 - (-t424 + t623 * t581 + t875 * t541 + (-t464 + t868) * t493 + t1030) * t53 + t7 * t753 + t49 * t679 + (t9 * t752 + t7 * t880 + t49 * t896 + (t53 * t752 + t769) * qJD(1)) * t583 + (t8 * t752 + t53 * t755 + (-t768 + t49 * (-t478 - t810)) * qJD(1) + t631) * t581 + (-t842 - (t401 - t876) * t541 - t868 * t494 + (-t623 + t755) * t583 + t1064) * t52) * m(7) + (-t68 * (-t353 * t541 - t520 * t581 + t741 + (-t396 - t464) * t493 + t1030) - t64 * (-t350 * t493 + t710 + (-t353 - t402) * t494) - (t562 * t657 + t563 * t590) * qJD(4) + t10 * t753 + t64 * t679 + (t12 * t803 + t10 * t274 + t64 * t128 + (t68 * t803 + t956) * qJD(1)) * t583 + (t11 * t803 + t68 * t812 + (t955 + t64 * (-t478 - t879)) * qJD(1) + t729) * t581 + (t396 * t494 - t656 - (t350 + t401) * t541 + (t520 + t812) * t583 + t1064) * t67) * m(6) + (-t87 * (t351 * t541 - t397 * t494 + t656) - t88 * (-t354 * t541 - t397 * t493 + t741) - ((-t269 * t87 + t275 * t88) * t562 + t594 * t563) * qJD(4) + t87 * t480 + t31 * t860 + (t31 * t275 + t87 * t885 + (qJD(1) * t88 + t42) * t869) * t583 + (qJD(1) * t394 * t87 + t31 * t269 + t41 * t869 + t88 * t885) * t581 + (t351 * t493 + t354 * t494 - t863 + t808 + (qJD(1) * t269 + t129) * t583 + (t132 + (-t275 - t478) * qJD(1)) * t581) * t89) * m(5) + (0.2e1 * t184 * (t232 * t583 + t233 * t581 + (t410 * t583 - t411 * t581) * qJD(1)) + t696 * t489 + (-t98 * t581 - t99 * t583 + (-t154 * t583 + t918) * qJD(1)) * t501 - (t153 * t463 - t917) * qJD(1) - (t184 * (-t463 * t581 - t465 * t583) + t696 * t738) * qJD(3)) * m(4) + t1041 * t996 + t1042 * t995 + ((t1016 * t562 + t1084 * t903) * qJD(4) + ((qJD(4) * t1105 + t1044) * t563 + t1045) * t583 + (t1046 * t484 + t1047 * t483) * t541 + (t1049 * t484 - t1051 * t483) * t494 + (-t1048 * t484 + t1050 * t483) * t493) * t994 + (qJD(1) * t1009 + t1059 * t581 - t1060 * t583) * t993 + (qJD(1) * t1008 + t1057 * t581 - t1058 * t583) * t992 + ((-t1017 * t562 + t1106 * t902) * qJD(4) + ((qJD(4) * t1085 + t1044) * t563 + t1045) * t581 + (t1046 * t482 + t1047 * t481) * t541 + (t1049 * t482 - t1051 * t481) * t494 + (-t1048 * t482 + t1050 * t481) * t493) * t991 + ((qJD(4) * t1010 - t1097) * t563 + ((t1046 * t582 + t1047 * t580 - t1111) * t541 + (t1049 * t582 - t1051 * t580 - t1116) * t494 + (-t1048 * t582 + t1050 * t580 - t1115) * t493 + t1015 * qJD(4)) * t562) * t990 + (qJD(1) * t1010 + t1055 * t581 - t1056 * t583) * t989 - t1052 * t833 / 0.2e1 + t1043 * t742 + (qJD(1) * t85 + (t581 * (t1006 * t581 + t592 * t583) - t583 * (t1005 * t581 + t593 * t583) + (t151 * t581 + t152 * t583) * qJD(1)) * t1077 + t1063) * t987 + (qJD(1) * t86 + (t581 * (-t1006 * t583 + t592 * t581) - t583 * (-t1005 * t583 + t593 * t581) + (t149 * t581 + t150 * t583) * qJD(1)) * t1077 + t1062) * t986 + (t650 + t69 + t1054) * t840 / 0.2e1 + (t649 + t70 + t1053) * t778 - (t1053 * t583 + t1054 * t581) * t832 / 0.2e1; (t9 * t877 + t52 * t757 + t53 * t878 + t7 * t277 + t49 * t897 + (qJD(3) * t588 + t52 * t894 - t53 * t818 - t8 * t810 + t883 * t9) * t563 + (t624 * qJD(3) + (qJD(1) * t625 + t53 * t813 + t8 * t804 + t631) * t583 + (-t9 * t1067 + t52 * t886 - t7 * t810 - t49 * t818 + (-t1067 * t53 + t49 * t811) * qJD(1)) * t581) * t562 - t52 * (-t459 * t494 + t858 + (-t312 - t313) * t541) - t53 * (t319 * t541 + t884 + (-t458 - t459) * t493) - t49 * (t313 * t493 + t890 + (t318 + t319) * t494) - (t52 * t672 + t53 * t671 - t970) * pkin(5)) * m(7) + (-t67 * (-t460 * t494 + (-t312 - t314) * t541 + t858) - t68 * (t320 * t541 + (-t458 - t460) * t493 + t884) - t64 * (t314 * t493 + (t318 + t320) * t494 + t890) + t12 * t877 + t67 * t757 + t68 * t878 + t10 * t277 + t64 * t897 + (qJD(3) * t590 - t11 * t879 + t12 * t268 + t67 * t131 - t68 * t895) * t563 + (t657 * qJD(3) + (qJD(1) * t658 + t11 * t870 + t68 * t887 + t729) * t583 + (t12 * t393 + t67 * t222 - t10 * t879 - t64 * t895 + (t68 * t393 + t64 * t882) * qJD(1)) * t581) * t562) * m(6) + (-t87 * (-t315 * t541 - t461 * t494) - t88 * (t321 * t541 - t461 * t493) - t89 * (t315 * t493 + t321 * t494) + (qJD(3) * t594 - t88 * t129 + t87 * t132 + t42 * t269 - t41 * t275) * t563 + (t87 * (-qJD(3) * t269 + t223 * t581) + t88 * (qJD(3) * t275 - t223 * t583) + t31 * t688 + t89 * (-t129 * t581 + t132 * t583 - t269 * t840 - t275 * t839) + (-t963 + t962 + (t583 * t87 + t957) * qJD(1)) * t394) * t562) * m(5) + (t1008 * t562 + t1017 * t563) * t996 + (t1009 * t562 - t1016 * t563) * t995 + (-t1011 * t583 + t1012 * t484 + t1013 * t483) * t994 + ((qJD(3) * t1009 - t1088) * t563 + (-qJD(1) * t1042 + t1016 * qJD(3) + t1059 * t583 + t1060 * t581) * t562) * t993 + ((qJD(3) * t1008 - t1087) * t563 + (-qJD(1) * t1041 - t1017 * qJD(3) + t1057 * t583 + t1058 * t581) * t562) * t992 + (-t1011 * t581 + t1012 * t482 + t1013 * t481) * t991 + (t1079 * t563 + (t1012 * t582 + t1013 * t580) * t562) * t990 + ((qJD(3) * t1010 - t1086) * t563 + (-qJD(1) * t1043 + t1015 * qJD(3) + t1055 * t583 + t1056 * t581) * t562) * t989 - (t953 + t954 + t968 - t969 + t951 + t952 + t966 - t967 + t949 + t950 + t964 - t965 + t1007) * t563 / 0.2e1 + t1062 * t906 / 0.2e1 + t1063 * t904 / 0.2e1 + t1052 * t777 + (t1010 * t562 - t1015 * t563) * t742 + (t562 * t778 + t563 * t775) * t1054 + (-t799 / 0.2e1 + t563 * t773) * t1053; (t1031 * t53 + t1032 * t52 + t481 * t8 + t483 * t9 + t49 * t639 + t7 * t907 - t970) * m(7) + (t10 * t907 + t11 * t481 + t12 * t483 + t1031 * t68 + t1032 * t67 + (-t159 + t639) * t64) * m(6); 0.2e1 * ((-t52 * t834 - t53 * t836 + t7) * t997 + (-t493 * t53 - t494 * t52) * t998) * t563 + 0.2e1 * ((-qJD(3) * t49 + t52 * t840 - t53 * t839 - t8 * t581 - t583 * t9) * t997 + (t49 * (-t493 * t581 - t494 * t583) + (t52 * t581 - t53 * t583) * t541) * t998) * t562;];
tauc  = t1(:);
