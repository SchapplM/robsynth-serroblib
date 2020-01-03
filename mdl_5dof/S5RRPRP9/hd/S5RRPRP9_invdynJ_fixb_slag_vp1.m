% Calculate vector of inverse dynamics joint torques for
% S5RRPRP9
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
%   pkin=[a2,a3,a4,a5,d1,d2,d4,theta3]';
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
% Datum: 2019-12-31 20:08
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tau = S5RRPRP9_invdynJ_fixb_slag_vp1(qJ, qJD, qJDD, g, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(5,1),zeros(3,1),zeros(8,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRPRP9_invdynJ_fixb_slag_vp1: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RRPRP9_invdynJ_fixb_slag_vp1: qJD has to be [5x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [5 1]), ...
  'S5RRPRP9_invdynJ_fixb_slag_vp1: qJDD has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RRPRP9_invdynJ_fixb_slag_vp1: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RRPRP9_invdynJ_fixb_slag_vp1: pkin has to be [8x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RRPRP9_invdynJ_fixb_slag_vp1: m has to be [6x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [6,3]), ...
  'S5RRPRP9_invdynJ_fixb_slag_vp1: rSges has to be [6x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [6 6]), ...
  'S5RRPRP9_invdynJ_fixb_slag_vp1: Icges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 20:05:42
% EndTime: 2019-12-31 20:07:27
% DurationCPUTime: 92.19s
% Computational Cost: add. (29001->1302), mult. (47620->1632), div. (0->0), fcn. (45109->8), ass. (0->597)
t1071 = Icges(5,1) + Icges(6,1);
t1088 = Icges(6,4) + Icges(5,5);
t1087 = Icges(5,6) - Icges(6,6);
t566 = pkin(8) + qJ(4);
t548 = sin(t566);
t1109 = (Icges(5,4) - Icges(6,5)) * t548;
t1070 = Icges(5,2) + Icges(6,3);
t1108 = Icges(6,2) + Icges(5,3);
t549 = cos(t566);
t1107 = -t1087 * t548 + t1088 * t549;
t1106 = t1071 * t549 - t1109;
t570 = sin(qJ(2));
t572 = cos(qJ(2));
t903 = Icges(5,4) * t549;
t677 = -Icges(5,2) * t548 + t903;
t1105 = -t1087 * t572 + t570 * t677;
t872 = t549 * t570;
t492 = Icges(6,5) * t872;
t874 = t548 * t570;
t1086 = Icges(6,3) * t874 - t1105 + t492;
t896 = Icges(6,5) * t549;
t673 = Icges(6,3) * t548 + t896;
t1020 = (t673 - t677) * t572 - t1087 * t570;
t1104 = t1107 * t572 + t1108 * t570;
t1078 = -t1088 * t572 + t1106 * t570;
t1019 = t1088 * t570 + t1106 * t572;
t1103 = (t1070 * t549 + t1109) * t570;
t1102 = (-t1087 * t549 - t1088 * t548) * t570;
t571 = sin(qJ(1));
t573 = cos(qJ(1));
t864 = t573 * t548;
t387 = -t571 * t549 + t572 * t864;
t865 = t572 * t573;
t785 = t549 * t865;
t388 = t548 * t571 + t785;
t867 = t570 * t573;
t182 = Icges(5,5) * t388 - Icges(5,6) * t387 + Icges(5,3) * t867;
t185 = Icges(6,4) * t388 + Icges(6,2) * t867 + Icges(6,6) * t387;
t1081 = t182 + t185;
t898 = Icges(6,5) * t387;
t191 = Icges(6,1) * t388 + Icges(6,4) * t867 + t898;
t360 = Icges(5,4) * t387;
t194 = Icges(5,1) * t388 + Icges(5,5) * t867 - t360;
t1093 = t191 + t194;
t357 = Icges(6,5) * t388;
t179 = Icges(6,6) * t867 + Icges(6,3) * t387 + t357;
t905 = Icges(5,4) * t388;
t188 = -Icges(5,2) * t387 + Icges(5,6) * t867 + t905;
t1095 = t179 - t188;
t1038 = t1081 * t867 + t1093 * t388 + t1095 * t387;
t866 = t571 * t572;
t385 = t548 * t866 + t549 * t573;
t386 = t549 * t866 - t864;
t868 = t570 * t571;
t180 = Icges(5,5) * t386 - Icges(5,6) * t385 + Icges(5,3) * t868;
t183 = Icges(6,4) * t386 + Icges(6,2) * t868 + Icges(6,6) * t385;
t1026 = t180 + t183;
t356 = Icges(6,5) * t386;
t178 = -Icges(6,6) * t868 - Icges(6,3) * t385 - t356;
t359 = Icges(5,4) * t386;
t186 = -Icges(5,2) * t385 + Icges(5,6) * t868 + t359;
t1027 = t178 + t186;
t355 = Icges(6,5) * t385;
t189 = Icges(6,1) * t386 + Icges(6,4) * t868 + t355;
t358 = Icges(5,4) * t385;
t193 = -Icges(5,1) * t386 - Icges(5,5) * t868 + t358;
t1094 = t189 - t193;
t1039 = t1026 * t867 - t1027 * t387 + t1094 * t388;
t800 = qJD(4) * t573;
t808 = qJD(2) * t571;
t450 = t570 * t800 + t808;
t802 = qJD(4) * t571;
t806 = qJD(2) * t573;
t451 = -t570 * t802 + t806;
t801 = qJD(4) * t572;
t539 = qJD(1) - t801;
t1080 = t1107 * t570 - t1108 * t572;
t986 = t1078 * t388 + t1080 * t867 + t1086 * t387;
t1051 = t1038 * t450 - t1039 * t451 + t539 * t986;
t1040 = t1081 * t868 + t1093 * t386 + t1095 * t385;
t1041 = t1026 * t868 - t1027 * t385 + t1094 * t386;
t987 = t1078 * t386 + t1080 * t868 + t1086 * t385;
t1052 = t1040 * t450 - t1041 * t451 + t987 * t539;
t1059 = rSges(6,1) + pkin(4);
t762 = t570 * t806;
t171 = qJD(1) * t385 - qJD(4) * t785 + (t762 - t802) * t548;
t651 = t548 * t539;
t811 = qJD(1) * t572;
t727 = -qJD(4) + t811;
t172 = t573 * t651 + (-t571 * t727 - t762) * t549;
t760 = t572 * t806;
t812 = qJD(1) * t571;
t766 = t570 * t812;
t432 = t760 - t766;
t87 = Icges(6,5) * t172 + Icges(6,6) * t432 - Icges(6,3) * t171;
t93 = Icges(5,4) * t172 + Icges(5,2) * t171 + Icges(5,6) * t432;
t1101 = t87 - t93;
t763 = t570 * t808;
t810 = qJD(1) * t573;
t173 = (t572 * t802 - t812) * t549 + (t572 * t810 - t763 - t800) * t548;
t174 = t571 * t651 + (t573 * t727 - t763) * t549;
t807 = qJD(2) * t572;
t761 = t571 * t807;
t431 = t570 * t810 + t761;
t88 = Icges(6,5) * t174 + Icges(6,6) * t431 + Icges(6,3) * t173;
t94 = Icges(5,4) * t174 - Icges(5,2) * t173 + Icges(5,6) * t431;
t1100 = t88 - t94;
t89 = Icges(5,5) * t172 + Icges(5,6) * t171 + Icges(5,3) * t432;
t91 = Icges(6,4) * t172 + Icges(6,2) * t432 - Icges(6,6) * t171;
t1099 = t89 + t91;
t90 = Icges(5,5) * t174 - Icges(5,6) * t173 + Icges(5,3) * t431;
t92 = Icges(6,4) * t174 + Icges(6,2) * t431 + Icges(6,6) * t173;
t1098 = t90 + t92;
t95 = Icges(6,1) * t172 + Icges(6,4) * t432 - Icges(6,5) * t171;
t97 = Icges(5,1) * t172 + Icges(5,4) * t171 + Icges(5,5) * t432;
t1097 = t95 + t97;
t96 = Icges(6,1) * t174 + Icges(6,4) * t431 + Icges(6,5) * t173;
t98 = Icges(5,1) * t174 - Icges(5,4) * t173 + Icges(5,5) * t431;
t1096 = t96 + t98;
t1037 = rSges(6,3) + qJ(5);
t567 = sin(pkin(8));
t568 = cos(pkin(8));
t426 = t567 * t866 + t568 * t573;
t863 = t573 * t567;
t427 = t568 * t866 - t863;
t255 = Icges(4,5) * t427 - Icges(4,6) * t426 + Icges(4,3) * t868;
t889 = Icges(3,3) * t573;
t389 = Icges(3,5) * t866 - Icges(3,6) * t868 - t889;
t527 = Icges(3,4) * t868;
t901 = Icges(3,5) * t573;
t393 = Icges(3,1) * t866 - t527 - t901;
t893 = Icges(3,6) * t573;
t391 = Icges(3,4) * t866 - Icges(3,2) * t868 - t893;
t878 = t391 * t570;
t657 = -t393 * t572 + t878;
t258 = Icges(4,4) * t427 - Icges(4,2) * t426 + Icges(4,6) * t868;
t261 = Icges(4,1) * t427 - Icges(4,4) * t426 + Icges(4,5) * t868;
t664 = -t258 * t426 + t261 * t427;
t1036 = t255 * t868 - t389 * t573 - t571 * t657 + t664;
t556 = Icges(3,4) * t572;
t679 = -Icges(3,2) * t570 + t556;
t392 = Icges(3,6) * t571 + t573 * t679;
t906 = Icges(3,4) * t570;
t479 = Icges(3,1) * t572 - t906;
t394 = Icges(3,5) * t571 + t479 * t573;
t323 = t394 * t866;
t475 = Icges(3,5) * t572 - Icges(3,6) * t570;
t390 = Icges(3,3) * t571 + t475 * t573;
t732 = t390 * t573 - t323;
t138 = -t392 * t868 - t732;
t428 = t568 * t571 - t572 * t863;
t870 = t567 * t571;
t429 = t568 * t865 + t870;
t257 = Icges(4,5) * t429 + Icges(4,6) * t428 + Icges(4,3) * t867;
t260 = Icges(4,4) * t429 + Icges(4,2) * t428 + Icges(4,6) * t867;
t263 = Icges(4,1) * t429 + Icges(4,4) * t428 + Icges(4,5) * t867;
t75 = t257 * t868 - t426 * t260 + t427 * t263;
t1035 = t138 + t75;
t849 = t571 * t390 + t394 * t865;
t140 = -t392 * t867 + t849;
t77 = t257 * t867 + t428 * t260 + t429 * t263;
t1033 = t140 + t77;
t1092 = t1020 * qJD(2) + t1103 * qJD(4);
t1091 = t1104 * qJD(2) + t1102 * qJD(4);
t400 = (-Icges(5,1) * t548 - t903) * t570;
t803 = qJD(4) * t570;
t1090 = (-Icges(6,1) * t548 + t896) * t803 + qJD(4) * t400 + t1019 * qJD(2);
t1089 = t1078 * t549 + t1086 * t548;
t675 = Icges(4,5) * t568 - Icges(4,6) * t567;
t377 = -Icges(4,3) * t572 + t570 * t675;
t678 = Icges(4,4) * t568 - Icges(4,2) * t567;
t379 = -Icges(4,6) * t572 + t570 * t678;
t682 = Icges(4,1) * t568 - Icges(4,4) * t567;
t381 = -Icges(4,5) * t572 + t570 * t682;
t476 = Icges(3,2) * t572 + t906;
t478 = Icges(3,1) * t570 + t556;
t653 = t476 * t570 - t478 * t572;
t474 = Icges(3,5) * t570 + Icges(3,6) * t572;
t875 = t474 * t573;
t1029 = t377 * t868 - t379 * t426 + t381 * t427 - t571 * t653 - t875;
t876 = t474 * t571;
t1028 = t377 * t867 + t379 * t428 + t381 * t429 - t573 * t653 + t876;
t709 = t386 * rSges(5,1) - t385 * rSges(5,2);
t198 = -rSges(5,3) * t868 - t709;
t928 = rSges(5,1) * t549;
t708 = -rSges(5,2) * t548 + t928;
t347 = -rSges(5,3) * t572 + t570 * t708;
t1085 = t198 * t539 - t347 * t451;
t354 = qJD(5) * t387;
t1073 = -t1037 * t385 - t1059 * t386;
t860 = rSges(6,2) * t868 - t1073;
t1084 = -t539 * t860 + t354;
t1058 = t1026 * t432 + t1027 * t171 + t1094 * t172 + t1096 * t388 + t1098 * t867 + t1100 * t387;
t1057 = t1081 * t432 + t1093 * t172 - t1095 * t171 + t1097 * t388 + t1099 * t867 + t1101 * t387;
t1056 = t1026 * t431 - t1027 * t173 + t1094 * t174 + t1096 * t386 + t1098 * t868 + t1100 * t385;
t1055 = t1081 * t431 + t1093 * t174 + t1095 * t173 + t1097 * t386 + t1099 * t868 + t1101 * t385;
t1047 = t1078 * t172 + t1080 * t432 - t1086 * t171 + t1090 * t388 + t1091 * t867 + t1092 * t387;
t1046 = t1078 * t174 + t1080 * t431 + t1086 * t173 + t1090 * t386 + t1091 * t868 + t1092 * t385;
t670 = -t178 * t548 + t189 * t549;
t70 = -t183 * t572 + t570 * t670;
t668 = -t186 * t548 - t193 * t549;
t72 = -t180 * t572 + t570 * t668;
t1083 = t70 + t72;
t669 = t179 * t548 + t191 * t549;
t71 = -t185 * t572 + t570 * t669;
t667 = -t188 * t548 + t194 * t549;
t73 = -t182 * t572 + t570 * t667;
t1082 = t71 + t73;
t985 = -t1080 * t572 + t1089 * t570;
t1079 = -t570 * t673 + t1105;
t850 = -t571 * t389 - t393 * t865;
t139 = -t391 * t867 - t850;
t998 = t428 * t258 + t429 * t261;
t76 = t255 * t867 + t998;
t917 = t573 * t76;
t1077 = t1033 * t571 - t139 * t573 - t917;
t1076 = t1035 * t571 - t1036 * t573;
t1075 = (-t1089 + t1104) * t539 + (t1080 * t571 + t668 + t670) * t451 + (-t1080 * t573 - t667 - t669) * t450;
t927 = rSges(4,2) * t567;
t929 = rSges(4,1) * t568;
t1001 = t927 - t929;
t1074 = t1001 * t572;
t836 = -t1037 * t872 + t1059 * t874;
t1072 = (t1089 * qJD(2) - t1091) * t572 + (t1090 * t549 + t1092 * t548 + (-t1078 * t548 + t1086 * t549) * qJD(4) + t1080 * qJD(2)) * t570;
t1069 = t1028 * qJD(1);
t1068 = -t1102 * t539 + (-t1087 * t386 - t1088 * t385) * t451 + (t1087 * t388 + t1088 * t387) * t450;
t564 = t573 * pkin(6);
t542 = pkin(3) * t568 + pkin(2);
t493 = t572 * t542;
t793 = pkin(3) * t863;
t715 = -t571 * t493 + t793;
t652 = t564 + t715;
t1066 = t652 + t1073;
t1065 = t1029 * qJD(1);
t797 = qJD(1) * qJD(2);
t465 = qJDD(2) * t571 + t573 * t797;
t794 = qJDD(4) * t570;
t283 = qJD(4) * t432 + t573 * t794 + t465;
t449 = qJD(2) * t803 - qJDD(4) * t572 + qJDD(1);
t552 = t570 * qJ(3);
t932 = pkin(2) - t542;
t752 = t932 * t572;
t1002 = -t752 - t552;
t504 = qJ(3) * t760;
t569 = -pkin(7) - qJ(3);
t992 = t570 * t932;
t629 = -t569 * t572 + t992;
t610 = t629 * t573;
t825 = qJD(1) * t793 + t569 * t766;
t126 = qJD(2) * t610 - t1002 * t812 - t504 + t825;
t540 = pkin(2) * t865;
t444 = qJ(3) * t867 + t540;
t521 = pkin(3) * t870;
t650 = t542 * t865 - t569 * t867 + t521;
t268 = t650 - t444;
t796 = qJD(2) * qJD(3);
t1003 = qJDD(3) * t570 + t572 * t796;
t804 = qJD(3) * t573;
t519 = t570 * t804;
t619 = -t571 * t811 - t762;
t269 = pkin(2) * t619 - qJ(3) * t766 + t504 + t519;
t499 = t573 * pkin(1) + t571 * pkin(6);
t547 = pkin(6) * t810;
t828 = qJD(1) * (-pkin(1) * t812 + t547) + qJDD(1) * t499;
t620 = qJDD(1) * t444 + t828 + t1003 * t571 + (t269 + t519) * qJD(1);
t563 = t572 * pkin(2);
t704 = t563 + t552;
t805 = qJD(3) * t572;
t422 = qJD(2) * t704 - t805;
t862 = qJ(3) + t569;
t993 = t570 * t862;
t852 = -(-t752 - t993) * qJD(2) - t422;
t730 = qJD(2) * t852;
t326 = t572 * t862 - t992;
t484 = pkin(2) * t570 - qJ(3) * t572;
t848 = -t326 - t484;
t581 = qJD(1) * t126 + qJDD(1) * t268 + t465 * t848 + t571 * t730 + t620;
t705 = rSges(6,1) * t549 + rSges(6,3) * t548;
t1000 = (-pkin(4) * t549 - qJ(5) * t548 - t705) * t570;
t845 = -rSges(6,2) * t572 - t1000;
t858 = rSges(6,2) * t867 + t1037 * t387 + t1059 * t388;
t799 = qJD(5) * t548;
t482 = t570 * t799;
t559 = t570 * rSges(6,2);
t618 = t548 * t807 + t549 * t803;
t861 = t482 + t618 * qJ(5) + (-t548 * t803 + t549 * t807) * pkin(4) + (-rSges(6,1) * t548 + rSges(6,3) * t549) * t803 + (t572 * t705 + t559) * qJD(2);
t995 = rSges(6,2) * t760 - t1037 * t171 + t1059 * t172 + t354;
t931 = -rSges(6,2) * t766 + t995;
t7 = qJD(5) * t173 + qJDD(5) * t385 - t283 * t845 + t449 * t858 - t450 * t861 + t539 * t931 + t581;
t1062 = -g(2) + t7;
t466 = -qJDD(2) * t573 + t571 * t797;
t284 = qJD(4) * t431 + t571 * t794 + t466;
t1061 = t1038 * t283 + t1039 * t284 + t1047 * t539 + t1057 * t450 - t1058 * t451 + t986 * t449;
t1060 = t1040 * t283 + t1041 * t284 + t1046 * t539 + t1055 * t450 - t1056 * t451 + t449 * t987;
t17 = (qJD(2) * t670 - t92) * t572 + (qJD(2) * t183 + t548 * t88 + t549 * t96 + (-t178 * t549 - t189 * t548) * qJD(4)) * t570;
t19 = (qJD(2) * t668 - t90) * t572 + (qJD(2) * t180 - t548 * t94 + t549 * t98 + (-t186 * t549 + t193 * t548) * qJD(4)) * t570;
t1054 = t17 + t19;
t18 = (qJD(2) * t669 - t91) * t572 + (qJD(2) * t185 + t548 * t87 + t549 * t95 + (t179 * t549 - t191 * t548) * qJD(4)) * t570;
t20 = (qJD(2) * t667 - t89) * t572 + (qJD(2) * t182 - t548 * t93 + t549 * t97 + (-t188 * t549 - t194 * t548) * qJD(4)) * t570;
t1053 = t18 + t20;
t1050 = t1082 * t450 - t1083 * t451 + t539 * t985;
t1049 = qJD(2) * t1076 + t1065;
t1048 = qJD(2) * t1077 + t1069;
t310 = qJD(1) * t428 + t567 * t763;
t311 = qJD(1) * t429 - t568 * t763;
t130 = Icges(4,5) * t311 + Icges(4,6) * t310 + Icges(4,3) * t431;
t132 = Icges(4,4) * t311 + Icges(4,2) * t310 + Icges(4,6) * t431;
t134 = Icges(4,1) * t311 + Icges(4,4) * t310 + Icges(4,5) * t431;
t624 = qJD(2) * t476;
t274 = qJD(1) * t392 - t571 * t624;
t625 = qJD(2) * t478;
t276 = qJD(1) * t394 - t571 * t625;
t997 = t258 * t567 - t261 * t568;
t1045 = (t130 - t274) * t572 + (t132 * t567 - t134 * t568 - t276) * t570 + (-t255 * t570 + t572 * t997 + t657) * qJD(2);
t308 = qJD(1) * t426 + t567 * t762;
t309 = -qJD(1) * t427 - t568 * t762;
t129 = Icges(4,5) * t309 + Icges(4,6) * t308 + Icges(4,3) * t432;
t131 = Icges(4,4) * t309 + Icges(4,2) * t308 + Icges(4,6) * t432;
t133 = Icges(4,1) * t309 + Icges(4,4) * t308 + Icges(4,5) * t432;
t273 = -t573 * t624 + (-t571 * t679 + t893) * qJD(1);
t275 = -t573 * t625 + (-t479 * t571 + t901) * qJD(1);
t877 = t392 * t570;
t656 = -t394 * t572 + t877;
t662 = -t260 * t567 + t263 * t568;
t1044 = (-t129 + t273) * t572 + (-t131 * t567 + t133 * t568 + t275) * t570 + (t257 * t570 + t572 * t662 - t656) * qJD(2);
t378 = Icges(4,3) * t570 + t572 * t675;
t328 = t378 * qJD(2);
t380 = Icges(4,6) * t570 + t572 * t678;
t329 = t380 * qJD(2);
t382 = Icges(4,5) * t570 + t572 * t682;
t330 = t382 * qJD(2);
t453 = t679 * qJD(2);
t454 = t479 * qJD(2);
t654 = t476 * t572 + t478 * t570;
t582 = qJD(1) * t474 - qJD(2) * t654 - t453 * t570 + t454 * t572;
t962 = t653 * qJD(1) + t475 * qJD(2);
t1043 = t308 * t379 + t309 * t381 + t328 * t867 + t329 * t428 + t330 * t429 + t377 * t432 + t571 * t962 + t582 * t573;
t1042 = t310 * t379 + t311 * t381 + t328 * t868 - t329 * t426 + t330 * t427 + t377 * t431 + t582 * t571 - t573 * t962;
t1034 = t139 + t76;
t221 = t391 * t572 + t393 * t570;
t1032 = -t255 * t572 - t570 * t997 + t221;
t222 = t392 * t572 + t394 * t570;
t1031 = -t257 * t572 + t570 * t662 + t222;
t1024 = t1079 * t571;
t1023 = t1079 * t573;
t1022 = t1078 * t571;
t1021 = t1078 * t573;
t1017 = t482 + t852 - t861;
t984 = -t427 * rSges(4,1) + t426 * rSges(4,2);
t1016 = t564 + t984;
t821 = t478 + t679;
t822 = -t476 + t479;
t658 = -t379 * t567 + t381 * t568;
t999 = qJD(2) * (-t571 * t662 - t573 * t997) + (t378 - t658) * qJD(1);
t1015 = t377 * t811 + t570 * t999 + (-t570 * t821 + t572 * t822) * qJD(1);
t1014 = t1075 * t570;
t607 = t391 * t573 - t392 * t571;
t959 = t571 * (-t476 * t573 + t394) - t573 * (-Icges(3,2) * t866 + t393 - t527);
t1013 = -t570 * t959 + (-t255 * t573 + t257 * t571 + t607) * t572;
t1012 = -t1026 * t451 + t1080 * t539 + t1081 * t450;
t510 = pkin(2) * t763;
t771 = -t431 * t569 - t542 * t763;
t127 = -qJ(3) * t761 + t510 + (t1002 * t573 + t521) * qJD(1) + t771;
t267 = (t563 + t993) * t571 + t715;
t522 = qJ(3) * t866;
t285 = t571 * t629 - t522;
t524 = qJ(3) * t865;
t286 = -t524 + t610;
t439 = -pkin(2) * t868 + t522;
t442 = -pkin(2) * t867 + t524;
t551 = qJD(3) * t570;
t772 = t439 * t808 + t442 * t806 + t551;
t517 = t571 * t551;
t818 = t517 - t510;
t270 = qJ(3) * t431 + qJD(1) * t540 + t818;
t441 = t704 * t571;
t779 = t573 * t269 + t571 * t270 + t441 * t810;
t1011 = t573 * t126 + t571 * t127 - t267 * t810 - t285 * t808 - t286 * t806 - t772 + t779;
t977 = t1082 * t573 + t1083 * t571;
t1010 = t1082 * t571 - t1083 * t573;
t976 = t1038 * t573 + t1039 * t571;
t1009 = t1038 * t571 - t1039 * t573;
t975 = t1040 * t573 + t1041 * t571;
t1008 = t1040 * t571 - t1041 * t573;
t814 = qJD(1) * t390;
t583 = -qJD(2) * t222 - t273 * t570 + t275 * t572 + t814;
t584 = qJD(1) * t389 - qJD(2) * t221 - t274 * t570 + t276 * t572;
t623 = qJD(2) * t474;
t963 = qJD(1) * t657 - t571 * t623 + t814;
t964 = -t573 * t623 + (-t475 * t571 + t656 + t889) * qJD(1);
t1007 = (-t130 * t868 + t132 * t426 - t134 * t427 - t255 * t431 - t258 * t310 - t261 * t311 + t573 * t963) * t573 + (t129 * t868 - t131 * t426 + t133 * t427 + t257 * t431 + t260 * t310 + t263 * t311 + t583 * t571 + (-t584 - t964) * t573) * t571;
t1006 = (-t130 * t867 - t132 * t428 - t134 * t429 - t255 * t432 - t258 * t308 - t261 * t309 - t584 * t573) * t573 + (t129 * t867 + t131 * t428 + t133 * t429 + t257 * t432 + t260 * t308 + t263 * t309 + t571 * t964 + (t583 - t963) * t573) * t571;
t558 = t570 * rSges(4,3);
t1005 = -t704 - t558;
t352 = qJD(5) * t385;
t819 = -t484 * t808 + t517;
t829 = t444 + t499;
t606 = (t268 + t829) * qJD(1) - t326 * t808 + t819;
t50 = -t450 * t845 + t539 * t858 + t352 + t606;
t737 = t50 * t845;
t714 = t441 * t808 + t444 * t806 - t805;
t626 = -t267 * t808 + t268 * t806 + t714;
t42 = t450 * t860 + t451 * t858 + t482 + t626;
t739 = t42 * t860;
t1004 = t737 - t739;
t264 = rSges(4,3) * t868 - t984;
t996 = -t1037 * t173 - t1059 * t174 - t352;
t990 = t806 * t848;
t266 = t429 * rSges(4,1) + t428 * rSges(4,2) + rSges(4,3) * t867;
t153 = t266 + t829;
t498 = pkin(1) * t571 - t564;
t468 = qJD(1) * t498;
t983 = -qJD(1) * t441 - t468;
t981 = (-t1078 + t1103) * t539 + (-t1070 * t386 + t1094 + t355 - t358) * t451 + (t1070 * t388 - t1093 + t360 - t898) * t450;
t980 = (-Icges(6,1) * t874 + t1086 + t400 + t492) * t539 + (t1071 * t385 + t1027 - t356 + t359) * t451 + (-t1071 * t387 + t1095 + t357 - t905) * t450;
t979 = t1068 * t570;
t970 = t1072 * t539 + t449 * t985;
t969 = g(1) * t573 + g(2) * t571;
t966 = t969 * t570;
t622 = -qJDD(3) * t572 + t269 * t806 + t270 * t808 + t465 * t441 + t570 * t796;
t830 = -t444 - t268;
t586 = t126 * t806 + t127 * t808 - t465 * t267 + t466 * t830 + t622;
t910 = rSges(6,2) * t431 - t996;
t5 = qJD(5) * t618 + qJDD(5) * t874 + t283 * t860 - t284 * t858 + t450 * t910 + t451 * t931 + t586;
t961 = t42 * t931 + t5 * t858;
t956 = m(4) / 0.2e1;
t955 = m(5) / 0.2e1;
t954 = m(6) / 0.2e1;
t953 = t283 / 0.2e1;
t952 = t284 / 0.2e1;
t951 = t449 / 0.2e1;
t950 = -t450 / 0.2e1;
t949 = t450 / 0.2e1;
t948 = -t451 / 0.2e1;
t947 = t451 / 0.2e1;
t946 = t465 / 0.2e1;
t945 = t466 / 0.2e1;
t944 = -t539 / 0.2e1;
t943 = t539 / 0.2e1;
t938 = g(1) * t571;
t935 = g(3) * t569;
t930 = rSges(3,1) * t572;
t926 = t17 * t451;
t925 = t18 * t450;
t924 = t19 * t451;
t923 = t20 * t450;
t557 = t570 * rSges(5,3);
t560 = t571 * rSges(3,3);
t632 = t519 + t990;
t831 = -t441 - t498;
t770 = t267 + t831;
t592 = qJD(1) * t770 + t632;
t58 = t1085 + t592;
t918 = t573 * t58;
t916 = t70 * t284;
t915 = t71 * t283;
t914 = t72 * t284;
t913 = t73 * t283;
t911 = -rSges(4,3) - qJ(3);
t485 = rSges(3,1) * t570 + rSges(3,2) * t572;
t764 = t485 * t806;
t816 = rSges(3,2) * t868 + t573 * rSges(3,3);
t413 = rSges(3,1) * t866 - t816;
t835 = -t413 - t498;
t219 = qJD(1) * t835 - t764;
t883 = t219 * t571;
t882 = t219 * t573;
t414 = rSges(3,1) * t865 - rSges(3,2) * t867 + t560;
t321 = t414 + t499;
t220 = qJD(1) * t321 - t485 * t808;
t443 = t485 * t573;
t881 = t220 * t443;
t873 = t548 * t572;
t871 = t549 * t572;
t869 = t569 * t570;
t857 = t1037 * t386 - t1059 * t385;
t856 = t1037 * t388 - t1059 * t387;
t855 = -t266 - t444;
t531 = rSges(6,2) * t866;
t854 = t1000 * t571 + t531;
t537 = rSges(6,2) * t865;
t853 = t1000 * t573 + t537;
t446 = t484 * t812;
t851 = t326 * t812 + t446;
t847 = -t493 + t869;
t846 = -(t558 - t1074) * qJD(2) - t422;
t844 = -t1037 * t873 - t1059 * t871 - t559;
t383 = -rSges(4,3) * t572 - t1001 * t570;
t841 = -t383 - t484;
t840 = t1005 + t1074;
t834 = t571 * t441 + t573 * t444;
t832 = qJD(1) * t442 + t571 * t805;
t789 = rSges(5,2) * t874;
t827 = rSges(5,3) * t866 + t571 * t789;
t826 = rSges(5,3) * t865 + t573 * t789;
t790 = t570 * t927;
t824 = rSges(4,3) * t866 + t571 * t790;
t823 = rSges(4,3) * t865 + t573 * t790;
t820 = rSges(3,2) * t766 + rSges(3,3) * t810;
t817 = t519 + t547;
t813 = qJD(1) * t475;
t809 = qJD(2) * t570;
t792 = t570 * t929;
t791 = rSges(5,1) * t872;
t784 = t569 * t866;
t783 = t569 * t865;
t782 = t172 * rSges(5,1) + t171 * rSges(5,2) + rSges(5,3) * t760;
t410 = (-rSges(5,1) * t548 - rSges(5,2) * t549) * t570;
t211 = qJD(4) * t410 + (t572 * t708 + t557) * qJD(2);
t780 = -t211 + t852;
t778 = -t264 + t831;
t776 = qJD(1) * t286 + t832;
t775 = t309 * rSges(4,1) + t308 * rSges(4,2) + rSges(4,3) * t760;
t774 = -t347 + t848;
t773 = t1003 * t573 + t466 * t484;
t200 = t388 * rSges(5,1) - t387 * rSges(5,2) + rSges(5,3) * t867;
t768 = t519 + t983;
t767 = -pkin(3) * t567 - pkin(6);
t754 = -pkin(1) - t930;
t749 = t810 / 0.2e1;
t747 = -t808 / 0.2e1;
t746 = t808 / 0.2e1;
t745 = -t806 / 0.2e1;
t744 = t806 / 0.2e1;
t741 = -pkin(1) - t493;
t740 = -pkin(1) + t869;
t49 = -t451 * t845 + t1084 + t592;
t738 = t49 * t845;
t113 = qJD(1) * t153 - t383 * t808 + t819;
t735 = t113 * t841;
t734 = t571 * t847;
t733 = t847 * t573;
t731 = -t389 + t877;
t729 = qJD(2) * t846;
t728 = -qJD(1) * t439 + t572 * t804;
t724 = -t571 * t267 + t573 * t268 + t834;
t723 = -t845 + t848;
t722 = -t517 - t771;
t489 = rSges(2,1) * t573 - rSges(2,2) * t571;
t486 = rSges(2,1) * t571 + rSges(2,2) * t573;
t488 = -rSges(3,2) * t570 + t930;
t713 = rSges(4,1) * t311 + rSges(4,2) * t310;
t710 = rSges(5,1) * t174 - rSges(5,2) * t173;
t666 = -t198 * t573 - t200 * t571;
t665 = -t220 * t571 - t882;
t279 = rSges(3,1) * t619 - rSges(3,2) * t760 + t820;
t440 = t485 * t571;
t280 = -qJD(2) * t440 + (t488 * t573 + t560) * qJD(1);
t661 = t279 * t573 + t280 * t571;
t655 = t413 * t571 + t414 * t573;
t349 = rSges(5,1) * t871 - rSges(5,2) * t873 + t557;
t649 = t741 - t559;
t648 = t741 - t557;
t462 = t499 * qJD(1);
t647 = -t270 - t462 - t517;
t628 = -qJD(1) * t285 + t728;
t621 = t570 * t911 - pkin(1) - t563;
t617 = t650 + t499;
t616 = t42 * t910 + t5 * t860;
t609 = -t49 * t860 + t50 * t858;
t608 = -t42 * t858 + t738;
t587 = -t542 * t762 - t569 * t760 + t817 + t825;
t57 = -t198 * t450 + t200 * t451 + t626;
t59 = t200 * t539 - t347 * t450 + t606;
t585 = t57 * t666 + (t571 * t58 - t573 * t59) * t347;
t580 = -t1004 * t573 + t608 * t571;
t578 = t466 * t326 + t770 * qJDD(1) + (-t127 + t647) * qJD(1) + t573 * t730 + t773;
t455 = t488 * qJD(2);
t430 = (t571 ^ 2 + t573 ^ 2) * t809;
t319 = -t573 * t792 + t823;
t318 = -t571 * t792 + t824;
t317 = t381 * t573;
t316 = t381 * t571;
t315 = t379 * t573;
t314 = t379 * t571;
t302 = -t573 * t791 + t826;
t300 = -t571 * t791 + t827;
t252 = -rSges(5,1) * t387 - rSges(5,2) * t388;
t247 = -rSges(5,1) * t385 - rSges(5,2) * t386;
t232 = qJD(1) * t267;
t212 = t655 * qJD(2);
t136 = rSges(4,3) * t431 + t713;
t135 = -rSges(4,3) * t766 + t775;
t112 = qJD(1) * t778 + t806 * t841 + t519;
t105 = qJD(1) * t279 + qJDD(1) * t414 - t455 * t808 - t465 * t485 + t828;
t104 = -t455 * t806 + t466 * t485 + t835 * qJDD(1) + (-t280 - t462) * qJD(1);
t102 = rSges(5,3) * t431 + t710;
t100 = -rSges(5,3) * t766 + t782;
t80 = (t264 * t571 + t266 * t573) * qJD(2) + t714;
t44 = qJD(1) * t135 + qJDD(1) * t266 + t465 * t841 + t571 * t729 + t620;
t43 = t383 * t466 + t573 * t729 + t778 * qJDD(1) + (-t136 + t647) * qJD(1) + t773;
t35 = t264 * t465 + t855 * t466 + (t135 * t573 + t136 * t571) * qJD(2) + t622;
t22 = t100 * t539 + t200 * t449 - t211 * t450 - t283 * t347 + t581;
t21 = -t102 * t539 + t198 * t449 - t211 * t451 + t284 * t347 + t578;
t8 = t100 * t451 + t102 * t450 - t198 * t283 - t200 * t284 + t586;
t6 = -qJD(5) * t171 + qJDD(5) * t387 + t284 * t845 - t449 * t860 - t451 * t861 - t539 * t910 + t578;
t1 = [(t1046 + t1051) * t948 + (t1043 + t1044) * t746 + t1047 * t949 + (t1042 - t1045 + t1048) * t745 + (t1028 + t1031) * t946 + (t1029 + t1032) * t945 + ((-t917 + (t138 - t323 + (t390 + t878) * t573 + t850) * t573 + (t77 + t849) * t571) * qJD(2) + t1069) * t744 + (m(2) * (t486 ^ 2 + t489 ^ 2) - t377 * t572 + t570 * t658 + t654 + Icges(2,3)) * qJDD(1) + (-(-qJD(1) * t413 - t219 - t468 - t764) * t220 + t220 * (t547 + t820) + (t485 * t883 - t881) * qJD(2) + ((-pkin(1) - t488) * t882 + (t219 * (-rSges(3,3) - pkin(6)) + t220 * t754) * t571) * qJD(1) + (t105 - g(2)) * t321 + (t104 - g(1)) * (t754 * t571 + t564 + t816)) * m(3) + t970 + t986 * t953 + t987 * t952 - t924 / 0.2e1 + t925 / 0.2e1 + t923 / 0.2e1 + (-(t1085 + t232 - t58 + t632 + t983) * t59 - (-pkin(1) + (-rSges(5,3) + t569) * t570) * t938 + t58 * (-t710 + t722) + t59 * (t587 + t782) + (t21 * (t740 - t557) - t58 * rSges(5,3) * t807) * t571 + (t648 * t918 + (t58 * t767 + t59 * t648) * t571) * qJD(1) + (-g(2) + t22) * (t617 + t200) + (-g(1) + t21) * (t652 - t709)) * m(5) + (t737 * t451 - g(1) * t1066 - (-pkin(1) + (-rSges(6,2) + t569) * t570) * t938 + ((t740 - t559) * t571 + t1066) * t6 + (t649 * t812 - t1084 - t232 + t587 - t768 - t990 + t995) * t50 + (t722 + t996 - rSges(6,2) * t761 + (t767 * t571 + t649 * t573) * qJD(1) + t50) * t49 + t1062 * (t617 + t858)) * m(6) - m(2) * (-g(1) * t486 + g(2) * t489) + t916 / 0.2e1 + t913 / 0.2e1 + t914 / 0.2e1 + t915 / 0.2e1 + ((-t328 + t453) * t572 + (-t329 * t567 + t330 * t568 + t454) * t570 + (t377 * t570 + t572 * t658 - t653) * qJD(2)) * qJD(1) + t1051 * t947 + (((t573 * t731 + t140 + t664 - t849) * t573 + (t571 * t731 + t1034 + t732 - t75 - t998) * t571) * qJD(2) + t1049 - t1065) * t747 + (t1026 * t572 + (t1027 * t548 - t1094 * t549) * t570 + t1083) * t450 * t944 + (-t735 * t806 - g(1) * t1016 - t621 * t938 + (-t713 - t818 + t621 * t810 + (-pkin(6) * qJD(1) + t911 * t807) * t571) * t112 + (t621 * t571 + t1016) * t43 + (-g(2) + t44) * t153 + (-pkin(2) * t762 + t112 + t504 - t768 + t775 + t817 + (t264 + (-pkin(1) + t1005) * t571) * qJD(1)) * t113) * m(4) - t926 / 0.2e1; (qJD(1) * t976 + t1057 * t571 - t1058 * t573) * t949 - (qJD(1) * t1042 + qJD(2) * t1007 + qJDD(1) * t1029 + t1035 * t465 + t1036 * t466 + t1060) * t573 / 0.2e1 + (qJD(1) * t1043 + qJD(2) * t1006 + qJDD(1) * t1028 + t1033 * t465 + t1034 * t466 + t1061) * t571 / 0.2e1 - t1050 * t803 / 0.2e1 + (t1048 + t1051) * t749 + (t1049 + t1052) * t812 / 0.2e1 + (qJD(1) * t977 + t1053 * t571 - t1054 * t573) * t943 + (qJD(1) * t975 + t1055 * t571 - t1056 * t573) * t948 + (t1045 * t573 + t1044 * t571 + (t1031 * t573 + t1032 * t571) * qJD(1)) * qJD(1) / 0.2e1 + (t1031 * t571 - t1032 * t573) * qJDD(1) / 0.2e1 + ((t1033 * t573 + t1034 * t571) * qJD(1) + t1006) * t746 + ((t1035 * t573 + t1036 * t571) * qJD(1) + t1007) * t745 - (t1051 * t573 + t1052 * t571) * t801 / 0.2e1 + (-t58 * (-t300 * t539 - t349 * t451 + t628) - t59 * (t302 * t539 - t349 * t450 + t776) - (t58 * t733 + t59 * t734) * qJD(2) - ((t198 * t58 + t200 * t59) * t570 + t585 * t572) * qJD(4) + t58 * t851 + t8 * t724 + (t8 * t200 + t58 * t780 + (qJD(1) * t59 + t21) * t774) * t573 + (qJD(1) * t347 * t58 - t198 * t8 + t22 * t774 + t59 * t780) * t571 - g(1) * (-t783 + t826) - g(2) * (-t784 + t827) - g(3) * (t349 + t493) - (-t935 + t969 * (-t542 - t928)) * t570 + (-t300 * t450 - t302 * t451 + (-qJD(1) * t198 + t100) * t573 + (t102 + (-t200 + t830) * qJD(1)) * t571 + t1011) * t57) * m(5) + ((t1039 * t866 + t570 * t986) * qJD(4) + ((qJD(4) * t1038 + t1012) * t572 + t1014) * t573 + (t1019 * t388 + t1020 * t387) * t539 + (t1022 * t388 - t1024 * t387) * t451 + (-t1021 * t388 + t1023 * t387) * t450) * t950 + ((t1040 * t865 + t570 * t987) * qJD(4) + ((qJD(4) * t1041 + t1012) * t572 + t1014) * t571 + (t1019 * t386 + t1020 * t385) * t539 + (t1022 * t386 - t1024 * t385) * t451 + (-t1021 * t386 + t1023 * t385) * t450) * t947 + ((qJD(4) * t977 - t1075) * t572 + ((t1019 * t549 + t1020 * t548 + t1080) * t539 + (t1022 * t549 - t1024 * t548 - t1026) * t451 + (-t1021 * t549 + t1023 * t548 + t1081) * t450 + t985 * qJD(4)) * t570) * t944 + t1008 * t952 + t1009 * t953 + t1010 * t951 + (g(1) * t443 + g(2) * t440 - g(3) * t488 - (t219 * t440 - t881) * qJD(1) - (t212 * (-t440 * t571 - t443 * t573) + t665 * t488) * qJD(2) + (qJD(2) * t661 + t413 * t465 - t414 * t466) * t655 + t212 * ((t413 * t573 - t414 * t571) * qJD(1) + t661) + t665 * t455 + (-t104 * t573 - t105 * t571 + (-t220 * t573 + t883) * qJD(1)) * t485) * m(3) + (-(t314 * t426 - t316 * t427) * t806 + (-t380 * t426 + t382 * t427) * qJD(1) + (-t806 * t876 - t813) * t573 + ((t315 * t426 - t317 * t427 + t573 * t875 + t1013) * qJD(2) + t1015) * t571) * t744 + ((-t315 * t428 - t317 * t429) * t808 + (t380 * t428 + t382 * t429) * qJD(1) + (-t808 * t875 + t813) * t571 + ((t314 * t428 + t316 * t429 + t571 * t876 + t1013) * qJD(2) + t1015) * t573) * t747 - ((t570 * t822 + t572 * t821) * qJD(1) + (t607 * t570 + t572 * t959) * qJD(2) - t999 * t572 + ((-t380 * t567 + t382 * t568 + t377) * qJD(1) + ((t315 * t567 - t317 * t568 + t257) * t571 - (t314 * t567 - t316 * t568 + t255) * t573) * qJD(2)) * t570) * qJD(1) / 0.2e1 + (t112 * t446 + t35 * t834 + t80 * t779 + (t43 * t841 + t112 * t846 + t35 * t266 + t80 * t135 + (t80 * t264 + t735) * qJD(1)) * t573 + (t44 * t841 + t113 * t846 + t35 * t264 + t80 * t136 + (t112 * t383 + t80 * t855) * qJD(1)) * t571 - g(1) * (t524 + t823) - g(2) * (t522 + t824) + g(3) * t840 - (-pkin(2) - t929) * t966 - t112 * (-qJD(1) * t318 + t728) - t113 * (qJD(1) * t319 + t832) - t80 * t772 - ((t112 * t840 + t80 * t319) * t573 + (t113 * t840 + t80 * t318) * t571) * qJD(2)) * m(4) + (-(t570 * t609 + t572 * t580) * qJD(4) - g(1) * (t537 - t783) - g(2) * (t531 - t784) - g(3) * (t493 - t844) - (-t935 + t969 * (-t1037 * t548 - t1059 * t549 - t542)) * t570 + t5 * t724 + (qJD(1) * t739 + t6 * t723 + t961) * t573 + (qJD(1) * t738 + t7 * t723 + t616) * t571 + (-t734 * qJD(2) + t1017 * t571 - t844 * t450 - t853 * t539 + t723 * t810 - t776) * t50 + (-t572 * t799 - t853 * t451 - t854 * t450 + (t830 - t858) * t812 + t1011) * t42 + (-t733 * qJD(2) + t1017 * t573 - t844 * t451 + t854 * t539 - t628 + t851) * t49) * m(6) + t1076 * t945 + t1077 * t946; (-m(4) - m(5) - m(6)) * (-g(3) * t572 + t966) - m(4) * (t112 * t432 + t113 * t431 + t430 * t80) - m(5) * (t430 * t57 + t431 * t59 + t432 * t58) - m(6) * (t42 * t430 + t431 * t50 + t432 * t49) + 0.2e1 * ((t112 * t806 + t113 * t808 - t35) * t956 + (t58 * t806 + t59 * t808 - t8) * t955 + (t49 * t806 + t50 * t808 - t5) * t954) * t572 + 0.2e1 * ((qJD(2) * t80 - t112 * t812 + t113 * t810 + t43 * t573 + t44 * t571) * t956 + (qJD(2) * t57 + t21 * t573 + t22 * t571 - t58 * t812 + t59 * t810) * t955 + (qJD(2) * t42 - t49 * t812 + t50 * t810 + t571 * t7 + t573 * t6) * t954) * t570; (t570 * t976 - t572 * t986) * t953 + (t570 * t975 - t572 * t987) * t952 + (t570 * t977 - t572 * t985) * t951 + (t387 * t981 + t388 * t980 - t573 * t979) * t950 + ((qJD(2) * t976 - t1047) * t572 + (-qJD(1) * t1009 + t986 * qJD(2) + t1057 * t573 + t1058 * t571) * t570) * t949 + ((qJD(2) * t975 - t1046) * t572 + (-qJD(1) * t1008 + t987 * qJD(2) + t1055 * t573 + t1056 * t571) * t570) * t948 + (t385 * t981 + t386 * t980 - t571 * t979) * t947 + (t1068 * t572 + (t548 * t981 + t549 * t980) * t570) * t944 + ((qJD(2) * t977 - t1072) * t572 + (-qJD(1) * t1010 + t985 * qJD(2) + t1053 * t573 + t1054 * t571) * t570) * t943 - (t915 + t916 + t925 - t926 + t913 + t914 + t923 - t924 + t970) * t572 / 0.2e1 + t1060 * t868 / 0.2e1 + t1061 * t867 / 0.2e1 + t1050 * t809 / 0.2e1 + (-g(1) * t856 - g(2) * t857 + g(3) * t836 - (t386 * t50 + t388 * t49 + t42 * t872) * qJD(5) - (-t49 * t857 + t50 * t856) * t539 - (t42 * t856 + t49 * t836) * t451 - (t42 * t857 + t50 * t836) * t450 + (qJD(2) * t580 + t49 * t910 - t50 * t931 + t6 * t860 - t7 * t858) * t572 + (t609 * qJD(2) + (qJD(1) * t608 - t50 * t861 - t7 * t845 + t616) * t573 + (qJD(1) * t1004 + t49 * t861 + t6 * t845 - t961) * t571) * t570) * m(6) + ((qJD(2) * t585 - t59 * t100 + t58 * t102 - t198 * t21 - t22 * t200) * t572 + (t58 * (qJD(2) * t198 + t211 * t571) + t59 * (qJD(2) * t200 - t211 * t573) + t8 * t666 + t57 * (-t100 * t571 + t102 * t573 + t198 * t812 - t200 * t810) + (t21 * t571 - t22 * t573 + (t571 * t59 + t918) * qJD(1)) * t347) * t570 - t58 * (-t247 * t539 - t410 * t451) - t59 * (t252 * t539 - t410 * t450) - t57 * (t247 * t450 + t252 * t451) - g(1) * t252 - g(2) * t247 - g(3) * t410) * m(5) + t1052 * (t570 * t749 + t572 * t746) + t1051 * (t572 * t744 - t766 / 0.2e1); (-t171 * t49 + t173 * t50 + t618 * t42 + (t450 * t50 + t451 * t49 - g(3) + t5) * t874 + (-t42 * t451 - t50 * t539 - g(1) + t6) * t387 + (-t42 * t450 + t49 * t539 + t1062) * t385) * m(6);];
tau = t1;
