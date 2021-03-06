% Calculate vector of inverse dynamics joint torques for
% S5RRRPP6
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
%   pkin=[a2,a3,a4,a5,d1,d2,d3,theta4]';
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
% Datum: 2019-12-31 21:03
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tau = S5RRRPP6_invdynJ_fixb_slag_vp1(qJ, qJD, qJDD, g, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(5,1),zeros(3,1),zeros(8,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRRPP6_invdynJ_fixb_slag_vp1: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RRRPP6_invdynJ_fixb_slag_vp1: qJD has to be [5x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [5 1]), ...
  'S5RRRPP6_invdynJ_fixb_slag_vp1: qJDD has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RRRPP6_invdynJ_fixb_slag_vp1: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RRRPP6_invdynJ_fixb_slag_vp1: pkin has to be [8x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RRRPP6_invdynJ_fixb_slag_vp1: m has to be [6x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [6,3]), ...
  'S5RRRPP6_invdynJ_fixb_slag_vp1: rSges has to be [6x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [6 6]), ...
  'S5RRRPP6_invdynJ_fixb_slag_vp1: Icges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 20:59:59
% EndTime: 2019-12-31 21:01:58
% DurationCPUTime: 103.63s
% Computational Cost: add. (32797->1394), mult. (57086->1738), div. (0->0), fcn. (54303->8), ass. (0->633)
t1136 = Icges(5,1) + Icges(6,1);
t1147 = Icges(6,4) + Icges(5,5);
t1146 = Icges(5,6) - Icges(6,6);
t614 = qJ(3) + pkin(8);
t596 = sin(t614);
t1170 = (Icges(5,4) - Icges(6,5)) * t596;
t1135 = Icges(5,2) + Icges(6,3);
t597 = cos(t614);
t1169 = t1136 * t597 - t1170;
t1168 = Icges(6,2) + Icges(4,3) + Icges(5,3);
t616 = sin(qJ(3));
t619 = cos(qJ(3));
t1167 = Icges(4,5) * t619 - Icges(4,6) * t616 - t1146 * t596 + t1147 * t597;
t617 = sin(qJ(2));
t620 = cos(qJ(2));
t980 = Icges(5,4) * t597;
t743 = -Icges(5,2) * t596 + t980;
t1166 = -t1146 * t620 + t617 * t743;
t948 = t597 * t617;
t535 = Icges(6,5) * t948;
t950 = t596 * t617;
t1145 = Icges(6,3) * t950 - t1166 + t535;
t973 = Icges(6,5) * t597;
t739 = Icges(6,3) * t596 + t973;
t1118 = (t739 - t743) * t620 - t1146 * t617;
t1143 = -t1147 * t620 + t1169 * t617;
t1117 = t1147 * t617 + t1169 * t620;
t1165 = (t1135 * t597 + t1170) * t617;
t1164 = t1167 * t620 + t1168 * t617;
t1163 = (-Icges(4,5) * t616 - Icges(4,6) * t619 - t1146 * t597 - t1147 * t596) * t617;
t1139 = t1167 * t617 - t1168 * t620;
t618 = sin(qJ(1));
t621 = cos(qJ(1));
t937 = t621 * t596;
t416 = -t618 * t597 + t620 * t937;
t938 = t620 * t621;
t417 = t596 * t618 + t597 * t938;
t983 = Icges(4,4) * t619;
t744 = -Icges(4,2) * t616 + t983;
t422 = -Icges(4,6) * t620 + t617 * t744;
t984 = Icges(4,4) * t616;
t748 = Icges(4,1) * t619 - t984;
t426 = -Icges(4,5) * t620 + t617 * t748;
t936 = t621 * t616;
t470 = t618 * t619 - t620 * t936;
t944 = t616 * t618;
t471 = t619 * t938 + t944;
t942 = t617 * t621;
t1060 = t1139 * t942 + t1143 * t417 + t1145 * t416 + t422 * t470 + t426 * t471;
t195 = Icges(5,5) * t417 - Icges(5,6) * t416 + Icges(5,3) * t942;
t198 = Icges(6,4) * t417 + Icges(6,2) * t942 + Icges(6,6) * t416;
t276 = Icges(4,5) * t471 + Icges(4,6) * t470 + Icges(4,3) * t942;
t1140 = t195 + t198 + t276;
t975 = Icges(6,5) * t416;
t204 = Icges(6,1) * t417 + Icges(6,4) * t942 + t975;
t392 = Icges(5,4) * t416;
t207 = Icges(5,1) * t417 + Icges(5,5) * t942 - t392;
t1154 = t204 + t207;
t389 = Icges(6,5) * t417;
t192 = Icges(6,6) * t942 + Icges(6,3) * t416 + t389;
t982 = Icges(5,4) * t417;
t201 = -Icges(5,2) * t416 + Icges(5,6) * t942 + t982;
t1156 = t192 - t201;
t985 = Icges(4,4) * t471;
t279 = Icges(4,2) * t470 + Icges(4,6) * t942 + t985;
t458 = Icges(4,4) * t470;
t282 = Icges(4,1) * t471 + Icges(4,5) * t942 + t458;
t1094 = t1140 * t942 + t1154 * t417 + t1156 * t416 + t470 * t279 + t471 * t282;
t940 = t618 * t620;
t414 = t596 * t940 + t597 * t621;
t415 = t597 * t940 - t937;
t943 = t617 * t618;
t193 = Icges(5,5) * t415 - Icges(5,6) * t414 + Icges(5,3) * t943;
t196 = Icges(6,4) * t415 + Icges(6,2) * t943 + Icges(6,6) * t414;
t468 = t616 * t940 + t619 * t621;
t939 = t619 * t620;
t469 = t618 * t939 - t936;
t274 = Icges(4,5) * t469 - Icges(4,6) * t468 + Icges(4,3) * t943;
t1093 = t193 + t196 + t274;
t388 = Icges(6,5) * t415;
t191 = -Icges(6,6) * t943 - Icges(6,3) * t414 - t388;
t391 = Icges(5,4) * t415;
t199 = -Icges(5,2) * t414 + Icges(5,6) * t943 + t391;
t1125 = t191 + t199;
t387 = Icges(6,5) * t414;
t202 = Icges(6,1) * t415 + Icges(6,4) * t943 + t387;
t390 = Icges(5,4) * t414;
t206 = -Icges(5,1) * t415 - Icges(5,5) * t943 + t390;
t1155 = t202 - t206;
t457 = Icges(4,4) * t469;
t277 = -Icges(4,2) * t468 + Icges(4,6) * t943 + t457;
t456 = Icges(4,4) * t468;
t281 = -Icges(4,1) * t469 - Icges(4,5) * t943 + t456;
t1095 = t1093 * t942 - t1125 * t416 + t1155 * t417 + t470 * t277 - t281 * t471;
t870 = qJD(3) * t621;
t876 = qJD(2) * t618;
t493 = t617 * t870 + t876;
t872 = qJD(3) * t618;
t874 = qJD(2) * t621;
t494 = -t617 * t872 + t874;
t871 = qJD(3) * t620;
t582 = qJD(1) - t871;
t1101 = t1060 * t582 + t1094 * t493 - t1095 * t494;
t1061 = t1139 * t943 + t1143 * t415 + t1145 * t414 - t422 * t468 + t426 * t469;
t1096 = t1140 * t943 + t1154 * t415 + t1156 * t414 - t468 * t279 + t469 * t282;
t1097 = t1093 * t943 - t1125 * t414 + t1155 * t415 - t277 * t468 - t281 * t469;
t1102 = t1061 * t582 + t1096 * t493 - t1097 * t494;
t825 = t620 * t870;
t829 = t617 * t874;
t181 = qJD(1) * t414 - t597 * t825 + (t829 - t872) * t596;
t879 = qJD(1) * t620;
t799 = -qJD(3) + t879;
t1053 = t618 * t799 + t829;
t716 = t621 * t582;
t182 = -t1053 * t597 + t596 * t716;
t827 = t620 * t874;
t880 = qJD(1) * t618;
t834 = t617 * t880;
t677 = t827 - t834;
t92 = Icges(6,5) * t182 + Icges(6,6) * t677 - Icges(6,3) * t181;
t98 = Icges(5,4) * t182 + Icges(5,2) * t181 + Icges(5,6) * t677;
t1160 = t92 - t98;
t830 = t617 * t876;
t878 = qJD(1) * t621;
t183 = (t620 * t872 - t880) * t597 + (t620 * t878 - t830 - t870) * t596;
t714 = t799 * t621;
t1116 = t714 - t830;
t717 = t618 * t582;
t184 = t1116 * t597 + t596 * t717;
t875 = qJD(2) * t620;
t828 = t618 * t875;
t678 = t617 * t878 + t828;
t93 = Icges(6,5) * t184 + Icges(6,6) * t678 + Icges(6,3) * t183;
t99 = Icges(5,4) * t184 - Icges(5,2) * t183 + Icges(5,6) * t678;
t1159 = t93 - t99;
t100 = Icges(6,1) * t182 + Icges(6,4) * t677 - Icges(6,5) * t181;
t102 = Icges(5,1) * t182 + Icges(5,4) * t181 + Icges(5,5) * t677;
t1158 = t100 + t102;
t101 = Icges(6,1) * t184 + Icges(6,4) * t678 + Icges(6,5) * t183;
t103 = Icges(5,1) * t184 - Icges(5,4) * t183 + Icges(5,5) * t678;
t1157 = t101 + t103;
t1153 = t1118 * qJD(2) + t1165 * qJD(3);
t435 = (-Icges(5,1) * t596 - t980) * t617;
t873 = qJD(3) * t617;
t1152 = (-Icges(6,1) * t596 + t973) * t873 + qJD(3) * t435 + t1117 * qJD(2);
t242 = t1053 * t616 + t619 * t716;
t715 = t582 * t616;
t243 = -t1053 * t619 + t621 * t715;
t127 = Icges(4,5) * t243 + Icges(4,6) * t242 + Icges(4,3) * t677;
t94 = Icges(5,5) * t182 + Icges(5,6) * t181 + Icges(5,3) * t677;
t96 = Icges(6,4) * t182 + Icges(6,2) * t677 - Icges(6,6) * t181;
t1151 = t127 + t94 + t96;
t244 = -t1116 * t616 + t619 * t717;
t877 = qJD(2) * t617;
t245 = t619 * t714 + (-t619 * t877 + t715) * t618;
t128 = Icges(4,5) * t245 + Icges(4,6) * t244 + Icges(4,3) * t678;
t95 = Icges(5,5) * t184 - Icges(5,6) * t183 + Icges(5,3) * t678;
t97 = Icges(6,4) * t184 + Icges(6,2) * t678 + Icges(6,6) * t183;
t1150 = t128 + t95 + t97;
t1149 = t1164 * qJD(2) + t1163 * qJD(3);
t1148 = t1143 * t597 + t1145 * t596 - t422 * t616 + t426 * t619;
t1127 = rSges(6,1) + pkin(4);
t1126 = rSges(6,3) + qJ(5);
t1144 = -t617 * t739 + t1166;
t130 = Icges(4,4) * t245 + Icges(4,2) * t244 + Icges(4,6) * t678;
t132 = Icges(4,1) * t245 + Icges(4,4) * t244 + Icges(4,5) * t678;
t1108 = t1093 * t677 + t1125 * t181 + t1150 * t942 + t1155 * t182 + t1157 * t417 + t1159 * t416 + t130 * t470 + t132 * t471 + t242 * t277 - t243 * t281;
t129 = Icges(4,4) * t243 + Icges(4,2) * t242 + Icges(4,6) * t677;
t131 = Icges(4,1) * t243 + Icges(4,4) * t242 + Icges(4,5) * t677;
t1107 = t1140 * t677 + t1151 * t942 + t1154 * t182 - t1156 * t181 + t1158 * t417 + t1160 * t416 + t129 * t470 + t131 * t471 + t242 * t279 + t243 * t282;
t1106 = t1093 * t678 - t1125 * t183 + t1150 * t943 + t1155 * t184 + t1157 * t415 + t1159 * t414 - t130 * t468 + t132 * t469 + t244 * t277 - t245 * t281;
t1105 = t1140 * t678 + t1151 * t943 + t1154 * t184 + t1156 * t183 + t1158 * t415 + t1160 * t414 - t129 * t468 + t131 * t469 + t244 * t279 + t245 * t282;
t423 = Icges(4,6) * t617 + t620 * t744;
t475 = (-Icges(4,2) * t619 - t984) * t617;
t291 = qJD(2) * t423 + qJD(3) * t475;
t427 = Icges(4,5) * t617 + t620 * t748;
t478 = (-Icges(4,1) * t616 - t983) * t617;
t294 = qJD(2) * t427 + qJD(3) * t478;
t1100 = t1139 * t678 + t1143 * t184 + t1145 * t183 + t1149 * t943 + t1152 * t415 + t1153 * t414 + t244 * t422 + t245 * t426 - t291 * t468 + t294 * t469;
t1098 = t1139 * t677 + t1143 * t182 - t1145 * t181 + t1149 * t942 + t1152 * t417 + t1153 * t416 + t242 * t422 + t243 * t426 + t291 * t470 + t294 * t471;
t736 = -t191 * t596 + t202 * t597;
t75 = -t196 * t620 + t617 * t736;
t734 = -t199 * t596 - t206 * t597;
t77 = -t193 * t620 + t617 * t734;
t731 = -t277 * t616 - t281 * t619;
t84 = -t274 * t620 + t617 * t731;
t1142 = t75 + t77 + t84;
t735 = t192 * t596 + t204 * t597;
t76 = -t198 * t620 + t617 * t735;
t733 = -t201 * t596 + t207 * t597;
t78 = -t195 * t620 + t617 * t733;
t730 = -t279 * t616 + t282 * t619;
t85 = -t276 * t620 + t617 * t730;
t1141 = t76 + t78 + t85;
t1059 = -t1139 * t620 + t1148 * t617;
t1138 = (-t1148 + t1164) * t582 + (t1139 * t618 + t731 + t734 + t736) * t494 + (-t1139 * t621 - t730 - t733 - t735) * t493;
t1078 = -t469 * rSges(4,1) + t468 * rSges(4,2);
t286 = -rSges(4,3) * t943 + t1078;
t1014 = rSges(4,2) * t616;
t1016 = rSges(4,1) * t619;
t781 = -t1014 + t1016;
t442 = -rSges(4,3) * t620 + t617 * t781;
t1137 = -t286 * t582 + t442 * t494;
t779 = t415 * rSges(5,1) - t414 * rSges(5,2);
t211 = -rSges(5,3) * t943 - t779;
t1123 = -t1126 * t414 - t1127 * t415;
t933 = rSges(6,2) * t943 - t1123;
t1134 = (t1148 * qJD(2) - t1149) * t620 + (-t291 * t616 + t294 * t619 + t1152 * t597 + t1153 * t596 + (-t1143 * t596 + t1145 * t597 - t422 * t619 - t426 * t616) * qJD(3) + t1139 * qJD(2)) * t617;
t1133 = -t1163 * t582 + (-Icges(4,5) * t468 - Icges(4,6) * t469 - t1146 * t415 - t1147 * t414) * t494 + (-Icges(4,5) * t470 + Icges(4,6) * t471 + t1146 * t417 + t1147 * t416) * t493;
t865 = qJD(1) * qJD(2);
t507 = qJDD(2) * t618 + t621 * t865;
t863 = qJDD(3) * t617;
t308 = qJD(3) * t677 + t621 * t863 + t507;
t491 = qJD(2) * t873 - qJDD(3) * t620 + qJDD(1);
t864 = qJD(2) * qJD(4);
t1114 = qJDD(4) * t617 + t620 * t864;
t593 = pkin(3) * t619 + pkin(2);
t1019 = pkin(2) - t593;
t556 = pkin(7) * t827;
t615 = -qJ(4) - pkin(7);
t868 = qJD(4) * t621;
t564 = t617 * t868;
t1013 = pkin(3) * qJD(3);
t851 = t619 * t1013;
t860 = pkin(3) * t936;
t794 = qJD(1) * t860 + t615 * t834 + t618 * t851 + t564;
t852 = t616 * t1013;
t125 = -t556 + (pkin(7) * t880 + t1019 * t874) * t617 + ((-qJD(2) * t615 - t852) * t621 + t1019 * t880) * t620 + t794;
t588 = pkin(2) * t938;
t487 = pkin(7) * t942 + t588;
t583 = pkin(3) * t944;
t713 = t593 * t938 - t615 * t942 + t583;
t284 = t713 - t487;
t679 = -t618 * t879 - t829;
t321 = pkin(2) * t679 - pkin(7) * t834 + t556;
t609 = t617 * pkin(7);
t611 = t620 * pkin(2);
t1075 = t611 + t609;
t503 = t1075 * qJD(2);
t540 = pkin(2) * t617 - pkin(7) * t620;
t543 = t621 * pkin(1) + t618 * pkin(6);
t595 = pkin(6) * t878;
t893 = qJD(1) * (-pkin(1) * t880 + t595) + qJDD(1) * t543;
t643 = qJD(1) * t321 + qJDD(1) * t487 - t503 * t876 - t507 * t540 + t893;
t630 = qJD(1) * t564 + t1114 * t618 + t582 * t125 + t491 * t284 + t643;
t1018 = pkin(7) + t615;
t1071 = t1019 * t617;
t365 = t1018 * t620 - t1071;
t774 = rSges(6,1) * t597 + rSges(6,3) * t596;
t1112 = (-pkin(4) * t597 - qJ(5) * t596 - t774) * t617;
t909 = -rSges(6,2) * t620 - t1112;
t837 = -t365 - t909;
t1072 = t1018 * t617;
t819 = t1019 * t620;
t869 = qJD(4) * t620;
t272 = -t617 * t852 - t869 + (-t819 - t1072) * qJD(2);
t867 = qJD(5) * t617;
t527 = t596 * t867;
t606 = t617 * rSges(6,2);
t676 = t596 * t875 + t597 * t873;
t928 = -(-rSges(6,1) * t596 + rSges(6,3) * t597) * t873 - (t620 * t774 + t606) * qJD(2) - t527 - t676 * qJ(5) - (-t596 * t873 + t597 * t875) * pkin(4);
t843 = -t272 + t928;
t930 = rSges(6,2) * t942 + t1126 * t416 + t1127 * t417;
t386 = qJD(5) * t416;
t1085 = rSges(6,2) * t827 - t1126 * t181 + t1127 * t182 + t386;
t992 = -rSges(6,2) * t834 + t1085;
t9 = qJD(5) * t183 + qJDD(5) * t414 + t308 * t837 + t491 * t930 + t493 * t843 + t582 * t992 + t630;
t1128 = t9 - g(2);
t1122 = t1144 * t618;
t1121 = t1144 * t621;
t1120 = t1143 * t618;
t1119 = t1143 * t621;
t1115 = -t503 + t843 + t527;
t347 = t365 * t880;
t545 = t620 * t593;
t946 = t615 * t617;
t366 = -t1075 + t545 - t946;
t490 = t540 * t880;
t585 = pkin(7) * t940;
t484 = -pkin(2) * t943 + t585;
t685 = -qJD(1) * t484 - t1075 * t874;
t1113 = -t618 * t365 * t871 + t494 * t366 - t620 * t868 + t347 + t490 - t685;
t612 = t621 * pkin(6);
t1111 = t612 + t1078;
t508 = -qJDD(2) * t621 + t618 * t865;
t309 = qJD(3) * t678 + t618 * t863 + t508;
t1110 = t1060 * t491 + t1094 * t308 + t1095 * t309 + t1098 * t582 + t1107 * t493 - t1108 * t494;
t1109 = t1061 * t491 + t1096 * t308 + t1097 * t309 + t1100 * t582 + t1105 * t493 - t1106 * t494;
t21 = (qJD(2) * t736 - t97) * t620 + (qJD(2) * t196 + t101 * t597 + t596 * t93 + (-t191 * t597 - t202 * t596) * qJD(3)) * t617;
t23 = (qJD(2) * t734 - t95) * t620 + (qJD(2) * t193 + t103 * t597 - t596 * t99 + (-t199 * t597 + t206 * t596) * qJD(3)) * t617;
t29 = (qJD(2) * t731 - t128) * t620 + (qJD(2) * t274 - t130 * t616 + t132 * t619 + (-t277 * t619 + t281 * t616) * qJD(3)) * t617;
t1104 = t21 + t23 + t29;
t22 = (qJD(2) * t735 - t96) * t620 + (qJD(2) * t198 + t100 * t597 + t596 * t92 + (t192 * t597 - t204 * t596) * qJD(3)) * t617;
t24 = (qJD(2) * t733 - t94) * t620 + (qJD(2) * t195 + t102 * t597 - t596 * t98 + (-t201 * t597 - t207 * t596) * qJD(3)) * t617;
t30 = (qJD(2) * t730 - t127) * t620 + (qJD(2) * t276 - t129 * t616 + t131 * t619 + (-t279 * t619 - t282 * t616) * qJD(3)) * t617;
t1103 = t22 + t24 + t30;
t1099 = t1059 * t582 + t1141 * t493 - t1142 * t494;
t1092 = t1138 * t617;
t1091 = -t1093 * t494 + t1139 * t582 + t1140 * t493;
t1057 = t1141 * t621 + t1142 * t618;
t1056 = t1094 * t621 + t1095 * t618;
t1055 = t1096 * t621 + t1097 * t618;
t1090 = t1141 * t618 - t1142 * t621;
t1089 = t1094 * t618 - t1095 * t621;
t1088 = t1096 * t618 - t1097 * t621;
t385 = qJD(5) * t414;
t841 = t284 + t930;
t1087 = -t493 * t909 + t841 * t582 + t385;
t1086 = -t1126 * t183 - t1127 * t184 - t385;
t459 = t468 * pkin(3);
t924 = t1126 * t415 - t1127 * t414;
t1080 = t459 - t924;
t460 = t470 * pkin(3);
t922 = t1126 * t417 - t1127 * t416;
t1079 = t460 + t922;
t894 = t487 + t543;
t1077 = t1126 * t948;
t599 = qJD(4) * t617;
t562 = t618 * t599;
t1076 = -t493 * t365 + t562;
t603 = Icges(3,4) * t620;
t745 = -Icges(3,2) * t617 + t603;
t521 = Icges(3,1) * t617 + t603;
t1074 = (-t1143 + t1165) * t582 + (-t1135 * t415 + t1155 + t387 - t390) * t494 + (t1135 * t417 - t1154 + t392 - t975) * t493;
t1073 = (-Icges(6,1) * t950 + t1145 + t435 + t535) * t582 + (t1136 * t414 + t1125 - t388 + t391) * t494 + (-t1136 * t416 + t1156 + t389 - t982) * t493;
t1064 = g(1) * t621 + g(2) * t618;
t1058 = t1133 * t617;
t1054 = t1059 * t491 + t1134 * t582;
t966 = Icges(3,3) * t621;
t420 = Icges(3,5) * t940 - Icges(3,6) * t943 - t966;
t569 = Icges(3,4) * t943;
t978 = Icges(3,5) * t621;
t428 = Icges(3,1) * t940 - t569 - t978;
t970 = Icges(3,6) * t621;
t424 = Icges(3,4) * t940 - Icges(3,2) * t943 - t970;
t955 = t424 * t617;
t724 = -t428 * t620 + t955;
t146 = -t420 * t621 - t618 * t724;
t518 = Icges(3,5) * t620 - Icges(3,6) * t617;
t517 = Icges(3,5) * t617 + Icges(3,6) * t620;
t681 = qJD(2) * t517;
t986 = Icges(3,4) * t617;
t522 = Icges(3,1) * t620 - t986;
t429 = Icges(3,5) * t618 + t522 * t621;
t425 = Icges(3,6) * t618 + t621 * t745;
t954 = t425 * t617;
t723 = -t429 * t620 + t954;
t1052 = -t621 * t681 + (-t518 * t618 + t723 + t966) * qJD(1);
t421 = Icges(3,3) * t618 + t518 * t621;
t882 = qJD(1) * t421;
t1051 = qJD(1) * t724 - t618 * t681 + t882;
t519 = Icges(3,2) * t620 + t986;
t720 = t519 * t617 - t521 * t620;
t1050 = qJD(1) * t720 + t518 * qJD(2);
t1049 = t1064 * t617;
t1048 = t618 * (-t519 * t621 + t429) - t621 * (-Icges(3,2) * t940 + t428 - t569);
t1044 = t493 * (-Icges(4,2) * t471 + t282 + t458) - t494 * (-Icges(4,2) * t469 - t281 - t456) + t582 * (t426 + t475);
t1043 = m(5) / 0.2e1;
t1042 = m(6) / 0.2e1;
t1041 = t308 / 0.2e1;
t1040 = t309 / 0.2e1;
t1039 = t491 / 0.2e1;
t1038 = -t493 / 0.2e1;
t1037 = t493 / 0.2e1;
t1036 = -t494 / 0.2e1;
t1035 = t494 / 0.2e1;
t1034 = t507 / 0.2e1;
t1033 = t508 / 0.2e1;
t1032 = -t582 / 0.2e1;
t1031 = t582 / 0.2e1;
t1026 = -rSges(4,3) - pkin(7);
t1025 = pkin(3) * t616;
t1024 = g(1) * t618;
t1021 = g(3) * t615;
t1020 = g(3) * t617;
t1017 = rSges(3,1) * t620;
t1015 = rSges(5,1) * t597;
t1012 = t21 * t494;
t1011 = t22 * t493;
t1010 = t23 * t494;
t1009 = t24 * t493;
t1008 = t29 * t494;
t1007 = t30 * t493;
t509 = t618 * t545;
t787 = -t509 + t860;
t283 = (t611 + t1072) * t618 + t787;
t485 = t1075 * t618;
t898 = t485 * t876 + t487 * t874;
t708 = -t493 * t283 - t869 + t898;
t213 = t417 * rSges(5,1) - t416 * rSges(5,2) + rSges(5,3) * t942;
t929 = t213 + t284;
t56 = -t211 * t493 + t494 * t929 + t708;
t1000 = t56 * t211;
t605 = t617 * rSges(4,3);
t604 = t617 * rSges(5,3);
t607 = t618 * rSges(3,3);
t778 = -rSges(5,2) * t596 + t1015;
t380 = -rSges(5,3) * t620 + t617 * t778;
t831 = t540 * t874;
t541 = pkin(1) * t618 - t612;
t886 = -t541 - t485;
t665 = qJD(1) * t886 - t831;
t644 = -t494 * t365 + t564 + t665;
t932 = t211 + t283;
t63 = -t380 * t494 + t582 * t932 + t644;
t999 = t63 * t380;
t998 = t75 * t309;
t997 = t76 * t308;
t996 = t77 * t309;
t995 = t78 * t308;
t994 = t84 * t309;
t993 = t85 * t308;
t991 = rSges(6,2) * t678 - t1086;
t529 = rSges(3,1) * t617 + rSges(3,2) * t620;
t832 = t529 * t874;
t884 = rSges(3,2) * t943 + t621 * rSges(3,3);
t443 = rSges(3,1) * t940 - t884;
t900 = -t443 - t541;
t233 = qJD(1) * t900 - t832;
t959 = t233 * t618;
t958 = t233 * t621;
t445 = rSges(3,1) * t938 - rSges(3,2) * t942 + t607;
t349 = t445 + t543;
t234 = qJD(1) * t349 - t529 * t876;
t483 = t529 * t621;
t957 = t234 * t483;
t952 = t517 * t618;
t951 = t517 * t621;
t949 = t596 * t620;
t947 = t597 * t620;
t945 = t616 * t617;
t846 = t182 * rSges(5,1) + t181 * rSges(5,2) + rSges(5,3) * t827;
t105 = -rSges(5,3) * t834 + t846;
t935 = t105 + t125;
t554 = pkin(2) * t830;
t800 = t620 * t852;
t706 = -t593 * t830 - t615 * t678 - t618 * t800 - t621 * t851 + t562;
t126 = -pkin(7) * t828 + t554 + (t583 + (-t819 - t609) * t621) * qJD(1) + t706;
t934 = t126 * t942 - t283 * t827;
t777 = -rSges(5,1) * t596 - rSges(5,2) * t597;
t439 = t777 * t617;
t221 = qJD(3) * t439 + (t620 * t778 + t604) * qJD(2);
t927 = -t221 - t272;
t926 = t284 * t877 + t617 * t347;
t925 = -t620 * t283 + t365 * t943;
t265 = -t414 * rSges(5,1) - t415 * rSges(5,2);
t923 = -t265 + t459;
t270 = -t416 * rSges(5,1) - t417 * rSges(5,2);
t921 = t270 + t460;
t481 = (-rSges(4,1) * t616 - rSges(4,2) * t619) * t617;
t297 = qJD(3) * t481 + (t620 * t781 + t605) * qJD(2);
t918 = -t297 - t503;
t322 = pkin(7) * t678 + qJD(1) * t588 - t554;
t504 = qJD(1) * t543;
t917 = -t322 - t504;
t573 = rSges(6,2) * t940;
t916 = t1112 * t618 + t573;
t580 = rSges(6,2) * t938;
t915 = t1112 * t621 + t580;
t913 = -t618 * t420 - t428 * t938;
t912 = t618 * t421 + t429 * t938;
t861 = pkin(3) * t945;
t911 = t582 * t460 + t493 * t861;
t910 = -t365 - t380;
t908 = -t1126 * t949 - t1127 * t947 - t606;
t902 = t1127 * t950 - t1077;
t901 = -t442 - t540;
t587 = pkin(7) * t938;
t486 = -pkin(2) * t942 + t587;
t899 = t484 * t876 + t486 * t874;
t896 = t618 * t485 + t621 * t487;
t856 = rSges(5,2) * t950;
t892 = rSges(5,3) * t940 + t618 * t856;
t891 = rSges(5,3) * t938 + t621 * t856;
t890 = -t519 + t522;
t889 = t521 + t745;
t857 = rSges(4,2) * t945;
t888 = rSges(4,3) * t940 + t618 * t857;
t887 = rSges(4,3) * t938 + t621 * t857;
t885 = rSges(3,2) * t834 + rSges(3,3) * t878;
t883 = t612 - t509;
t881 = qJD(1) * t518;
t246 = -t618 * t720 - t951;
t866 = t246 * qJD(1);
t859 = t617 * t1016;
t858 = rSges(5,1) * t948;
t850 = t125 + t992;
t849 = t615 * t940;
t848 = t615 * t938;
t844 = t494 * t283;
t842 = t283 - t933;
t840 = -t503 + t927;
t839 = t243 * rSges(4,1) + t242 * rSges(4,2) + rSges(4,3) * t827;
t838 = t621 * t321 + t618 * t322 + t485 * t878;
t836 = -t540 + t910;
t287 = t471 * rSges(4,1) + t470 * rSges(4,2) + rSges(4,3) * t942;
t835 = -pkin(6) - t1025;
t821 = -pkin(1) - t1017;
t816 = t878 / 0.2e1;
t814 = -t876 / 0.2e1;
t813 = t876 / 0.2e1;
t812 = -t874 / 0.2e1;
t811 = t874 / 0.2e1;
t808 = -pkin(1) - t545;
t807 = -pkin(1) + t946;
t44 = t493 * t933 + t494 * t841 + t527 + t708;
t806 = t44 * t933;
t52 = -t494 * t909 + t582 * t842 + t386 + t644;
t805 = t52 * t909;
t489 = t540 * t876;
t719 = qJD(1) * t894 - t489;
t680 = t719 + t1076;
t64 = -t380 * t493 + t582 * t929 + t680;
t804 = t64 * t910;
t359 = t429 * t940;
t802 = t421 * t621 - t359;
t801 = -t420 + t954;
t798 = t620 * t126 + t272 * t943 + t365 * t678;
t796 = -t618 * t283 + t621 * t284 + t896;
t795 = -t540 + t837;
t53 = t1087 + t680;
t786 = t53 * t837;
t784 = qJD(1) * t486 - t1075 * t876;
t532 = rSges(2,1) * t621 - rSges(2,2) * t618;
t530 = rSges(2,1) * t618 + rSges(2,2) * t621;
t531 = -rSges(3,2) * t617 + t1017;
t783 = rSges(4,1) * t245 + rSges(4,2) * t244;
t780 = rSges(5,1) * t184 - rSges(5,2) * t183;
t107 = rSges(5,3) * t678 + t780;
t712 = t321 * t874 + t322 * t876 + t507 * t485 - t508 * t487;
t642 = -qJDD(4) * t620 + t493 * t126 - t308 * t283 + t617 * t864 + t712;
t10 = t107 * t493 - t211 * t308 - t309 * t929 + t494 * t935 + t642;
t772 = -t10 * t211 + t56 * t107;
t238 = t425 * t620 + t429 * t617;
t682 = qJD(2) * t519;
t292 = -t621 * t682 + (-t618 * t745 + t970) * qJD(1);
t683 = qJD(2) * t521;
t295 = -t621 * t683 + (-t522 * t618 + t978) * qJD(1);
t632 = -qJD(2) * t238 - t292 * t617 + t295 * t620 + t882;
t237 = t424 * t620 + t428 * t617;
t293 = qJD(1) * t425 - t618 * t682;
t296 = qJD(1) * t429 - t618 * t683;
t633 = qJD(1) * t420 - qJD(2) * t237 - t293 * t617 + t296 * t620;
t771 = -(t1051 * t618 + t633 * t621) * t621 + (t1052 * t618 + t632 * t621) * t618;
t770 = -(-t1051 * t621 + t633 * t618) * t621 + t618 * (-t1052 * t621 + t632 * t618);
t688 = -t615 * t620 + t1071;
t317 = t618 * t688 - t585;
t751 = -t283 * t825 + t493 * t317 + t599 + t899;
t750 = t595 + t794;
t147 = -t425 * t943 - t802;
t738 = -t146 * t621 + t147 * t618;
t148 = -t424 * t942 - t913;
t149 = -t425 * t942 + t912;
t737 = -t148 * t621 + t149 * t618;
t732 = -t234 * t618 - t958;
t729 = -t286 * t621 - t287 * t618;
t298 = rSges(3,1) * t679 - rSges(3,2) * t827 + t885;
t482 = t529 * t618;
t299 = -qJD(2) * t482 + (t531 * t621 + t607) * qJD(1);
t728 = t298 * t621 + t299 * t618;
t722 = t443 * t618 + t445 * t621;
t721 = t519 * t620 + t521 * t617;
t718 = t612 + t787;
t444 = rSges(4,1) * t939 - t1014 * t620 + t605;
t382 = rSges(5,1) * t947 - rSges(5,2) * t949 + t604;
t711 = t808 - t606;
t710 = t808 - t604;
t512 = qJD(1) * t541;
t709 = -qJD(1) * t485 - t512 - t831;
t707 = t621 * t125 + t618 * t126 - t283 * t878 + t838;
t687 = -t56 * t929 + t999;
t686 = t64 * t213 + t63 * t932;
t684 = t1026 * t617 - pkin(1) - t611;
t675 = t713 + t543;
t7 = qJD(5) * t676 + qJDD(5) * t950 + t308 * t933 - t309 * t841 + t493 * t991 + t494 * t850 + t642;
t674 = t44 * t991 + t7 * t933;
t667 = t582 * t283 + t564 + t709;
t318 = t621 * t688 - t587;
t666 = t284 * t873 + t582 * t318 + t618 * t869 + t784;
t664 = t424 * t621 - t425 * t618;
t662 = -t44 * t841 + t805;
t661 = t52 * t842 + t53 * t930;
t660 = qJDD(1) * t886 - t503 * t874 + t508 * t540;
t645 = (-t617 * t889 + t620 * t890) * qJD(1);
t641 = -t593 * t877 - t615 * t875 - t800;
t636 = (Icges(4,1) * t470 - t279 - t985) * t493 - (-Icges(4,1) * t468 - t277 - t457) * t494 + (-t422 + t478) * t582;
t114 = -t1137 + t665;
t115 = t287 * t582 - t442 * t493 + t719;
t83 = -t286 * t493 + t287 * t494 + t898;
t634 = t83 * t729 + (t114 * t618 - t115 * t621) * t442;
t496 = t745 * qJD(2);
t497 = t522 * qJD(2);
t631 = qJD(1) * t517 - qJD(2) * t721 - t496 * t617 + t497 * t620;
t629 = (t804 - t1000) * t621 + t687 * t618;
t628 = -t1048 * t617 + t664 * t620;
t627 = t309 * t365 + (-t562 + t917) * qJD(1) - t494 * t272 + t660 + t1114 * t621;
t626 = (t786 + t806) * t621 + t662 * t618;
t498 = t531 * qJD(2);
t357 = -t621 * t859 + t887;
t356 = -t618 * t859 + t888;
t355 = t426 * t621;
t354 = t426 * t618;
t353 = t422 * t621;
t352 = t422 * t618;
t346 = t494 * t620 - t582 * t943;
t345 = t493 * t620 + t582 * t942;
t343 = t493 * t459;
t342 = -t621 * t858 + t891;
t340 = -t618 * t858 + t892;
t320 = rSges(4,1) * t470 - rSges(4,2) * t471;
t319 = -rSges(4,1) * t468 - rSges(4,2) * t469;
t302 = (t493 * t618 + t494 * t621) * t617;
t247 = -t621 * t720 + t952;
t232 = t247 * qJD(1);
t231 = t283 * t942;
t226 = t722 * qJD(2);
t134 = rSges(4,3) * t678 + t783;
t133 = -rSges(4,3) * t834 + t839;
t110 = qJD(1) * t298 + qJDD(1) * t445 - t498 * t876 - t507 * t529 + t893;
t109 = -t498 * t874 + t508 * t529 + t900 * qJDD(1) + (-t299 - t504) * qJD(1);
t91 = -t1050 * t621 + t631 * t618;
t90 = t1050 * t618 + t631 * t621;
t89 = -qJD(2) * t723 + t292 * t620 + t295 * t617;
t88 = -qJD(2) * t724 + t293 * t620 + t296 * t617;
t66 = qJD(2) * t737 + t232;
t65 = qJD(2) * t738 + t866;
t46 = t133 * t582 + t287 * t491 - t297 * t493 - t308 * t442 + t643;
t45 = qJD(1) * t917 - t134 * t582 + t286 * t491 - t297 * t494 + t309 * t442 + t660;
t35 = t133 * t494 + t134 * t493 - t286 * t308 - t287 * t309 + t712;
t20 = t105 * t582 + t213 * t491 + t308 * t910 + t493 * t927 + t630;
t19 = -t221 * t494 + t309 * t380 + (-t107 - t126) * t582 + t932 * t491 + t627;
t8 = -qJD(5) * t181 + qJDD(5) * t416 + t928 * t494 + t909 * t309 + (-t126 - t991) * t582 + t842 * t491 + t627;
t1 = [(-t866 + ((t621 * t801 + t149 - t912) * t621 + (t618 * t801 + t148 + t802) * t618) * qJD(2) + t65) * t814 + (t88 + t91 + t66) * t812 + (m(2) * (t530 ^ 2 + t532 ^ 2) + Icges(2,3) + t721) * qJDD(1) + (t246 + t237) * t1033 + (t1093 * t620 + (t1125 * t596 - t1155 * t597 - t731) * t617 + t1142) * t493 * t1032 + (-(-qJD(1) * t443 - t233 - t512 - t832) * t234 + t234 * (t595 + t885) + (t529 * t959 - t957) * qJD(2) + ((-pkin(1) - t531) * t958 + (t233 * (-rSges(3,3) - pkin(6)) + t234 * t821) * t618) * qJD(1) + (t110 - g(2)) * t349 + (t109 - g(1)) * (t618 * t821 + t612 + t884)) * m(3) + t1054 + ((t554 - t783 + t684 * t878 + (-pkin(6) * qJD(1) + t1026 * t875) * t618) * t114 - g(1) * t1111 - t684 * t1024 + (t618 * t684 + t1111) * t45 + (t46 - g(2)) * (t287 + t894) + (t114 - t709 - pkin(2) * t829 + t556 + t595 + t839 + (-pkin(1) - t1075 - t605) * t880 + t1137) * t115) * m(4) + (-g(1) * (t718 - t779) - (-pkin(1) + (-rSges(5,3) + t615) * t617) * t1024 + t19 * (-t779 + t883) + t63 * (-t706 - t780) + t64 * (t750 + t846) + (t1025 * t19 + t64 * t641) * t621 + (t19 * (t807 - t604) - t63 * rSges(5,3) * t875) * t618 + (t63 * t710 * t621 + (t63 * t835 + t64 * t710) * t618) * qJD(1) + t56 * t844 - (t211 * t582 - t63 + t667) * t64 - (t56 * t283 + t804) * t494 + (-g(2) + t20) * (t675 + t213)) * m(5) + (-g(1) * (t718 + t1123) - (-pkin(1) + (-rSges(6,2) + t615) * t617) * t1024 + t44 * t844 - (t44 * t283 + t786) * t494 + (t883 + t1123 + t1025 * t621 + (t807 - t606) * t618) * t8 + (t582 * t933 + t621 * t641 + t711 * t880 + t1085 - t386 - t667 + t750) * t53 + (-rSges(6,2) * t828 + t1076 + t1086 - t489 - t706 + (t618 * t835 + t621 * t711 + t894) * qJD(1) + t1087) * t52 + t1128 * (t675 + t930)) * m(6) + t1098 * t1037 + (t1100 + t1101) * t1036 - m(2) * (-g(1) * t530 + g(2) * t532) + t994 / 0.2e1 + t995 / 0.2e1 + t993 / 0.2e1 + t1101 * t1035 + (t89 + t90) * t813 + t1060 * t1041 + t1061 * t1040 + t1011 / 0.2e1 - t1012 / 0.2e1 + (t232 + ((t147 - t359 + (t421 + t955) * t621 + t913) * t621 + t912 * t618) * qJD(2)) * t811 + t998 / 0.2e1 + t1009 / 0.2e1 - t1010 / 0.2e1 + (t247 + t238) * t1034 + t996 / 0.2e1 + t997 / 0.2e1 + t1007 / 0.2e1 - t1008 / 0.2e1 + (-qJD(2) * t720 + t496 * t620 + t497 * t617) * qJD(1); ((-t874 * t952 - t881) * t621 + (t645 + (t621 * t951 + t628) * qJD(2)) * t618) * t811 + ((-t876 * t951 + t881) * t618 + (t645 + (t618 * t952 + t628) * qJD(2)) * t621) * t814 - (t1101 * t621 + t1102 * t618) * t871 / 0.2e1 + (qJD(1) * t1056 + t1107 * t618 - t1108 * t621) * t1037 - (qJD(1) * t91 + qJD(2) * t770 + qJDD(1) * t246 + t146 * t508 + t147 * t507 + t1109) * t621 / 0.2e1 + (qJD(1) * t90 + qJD(2) * t771 + qJDD(1) * t247 + t148 * t508 + t149 * t507 + t1110) * t618 / 0.2e1 + (t114 * t490 + t35 * t896 + (t114 * t918 + t35 * t287 + (qJD(1) * t115 + t45) * t901) * t621 + (t114 * t442 * qJD(1) + t115 * t918 - t286 * t35 + t46 * t901) * t618 - t114 * (-t356 * t582 - t444 * t494 + t685) - t115 * (t357 * t582 - t444 * t493 + t784) - ((t114 * t286 + t115 * t287) * t617 + t634 * t620) * qJD(3) - g(1) * (t587 + t887) - g(2) * (t585 + t888) - g(3) * (t444 + t1075) - (-pkin(2) - t1016) * t1049 + (t838 + (-qJD(1) * t286 + t133) * t621 + (t134 + (-t287 - t487) * qJD(1)) * t618 - t356 * t493 - t357 * t494 - t899) * t83) * m(4) + (t7 * t796 + (t806 * qJD(1) + t7 * t930 + t8 * t795) * t621 + (t805 * qJD(1) + t9 * t795 + t674) * t618 - g(1) * (t580 - t848) - g(2) * (t573 - t849) - g(3) * (t545 - t908) - (-t1021 + t1064 * (-t1126 * t596 - t1127 * t597 - t593)) * t617 - (t617 * t661 + t620 * t626) * qJD(3) + (t795 * t878 - t666 - t915 * t582 - (-t366 + t908) * t493 + t1115 * t618) * t53 + (t707 + t992 * t621 + (-t487 - t841) * t880 - qJD(5) * t949 - t751 - (t318 + t915) * t494 - t916 * t493) * t44 + (-(-t317 - t916) * t582 - t908 * t494 + t1115 * t621 + t1113) * t52) * m(6) + (-t64 * (t342 * t582 + t666 + (-t366 - t382) * t493) - t56 * (t340 * t493 + t751 + (t318 + t342) * t494) - (t617 * t686 + t620 * t629) * qJD(3) + t10 * t796 + t56 * t707 + (t19 * t836 + t10 * t213 + t56 * t105 + (t64 * t836 - t1000) * qJD(1)) * t621 + (t20 * t836 + t64 * t840 + (t999 + t56 * (-t487 - t929)) * qJD(1) + t772) * t618 - g(1) * (-t848 + t891) - g(2) * (-t849 + t892) - g(3) * (t382 + t545) - (-t1021 + t1064 * (-t593 - t1015)) * t617 + (t382 * t494 - (-t317 - t340) * t582 + t840 * t621 + t1113) * t63) * m(5) + ((qJD(2) * t728 + t443 * t507 - t445 * t508) * t722 + t226 * ((t443 * t621 - t445 * t618) * qJD(1) + t728) + t732 * t498 + (-t109 * t621 - t110 * t618 + (-t234 * t621 + t959) * qJD(1)) * t529 + g(1) * t483 + g(2) * t482 - g(3) * t531 - (t233 * t482 - t957) * qJD(1) - (t226 * (-t482 * t618 - t483 * t621) + t732 * t531) * qJD(2)) * m(3) + t1088 * t1040 + t1089 * t1041 + t1090 * t1039 + ((t1060 * t617 + t1095 * t940) * qJD(3) + ((qJD(3) * t1094 + t1091) * t620 + t1092) * t621 + (t1117 * t417 + t1118 * t416 + t423 * t470 + t427 * t471) * t582 + (t1120 * t417 - t1122 * t416 + t352 * t470 + t354 * t471) * t494 + (-t1119 * t417 + t1121 * t416 - t353 * t470 - t355 * t471) * t493) * t1038 + ((t1061 * t617 + t1096 * t938) * qJD(3) + ((qJD(3) * t1097 + t1091) * t620 + t1092) * t618 + (t1117 * t415 + t1118 * t414 - t423 * t468 + t427 * t469) * t582 + (t1120 * t415 - t1122 * t414 - t352 * t468 + t354 * t469) * t494 + (-t1119 * t415 + t1121 * t414 + t353 * t468 - t355 * t469) * t493) * t1035 + ((qJD(3) * t1057 - t1138) * t620 + ((t1117 * t597 + t1118 * t596 - t423 * t616 + t427 * t619 + t1139) * t582 + (t1120 * t597 - t1122 * t596 - t352 * t616 + t354 * t619 - t1093) * t494 + (-t1119 * t597 + t1121 * t596 + t353 * t616 - t355 * t619 + t1140) * t493 + t1059 * qJD(3)) * t617) * t1032 - t1099 * t873 / 0.2e1 + (t66 + t1101) * t816 + (t65 + t1102) * t880 / 0.2e1 + (qJD(1) * t1057 + t1103 * t618 - t1104 * t621) * t1031 + (qJD(1) * t1055 + t1105 * t618 - t1106 * t621) * t1036 + qJDD(1) * (-t237 * t621 + t238 * t618) / 0.2e1 - qJD(1) * ((t617 * t890 + t620 * t889) * qJD(1) + (t1048 * t620 + t664 * t617) * qJD(2)) / 0.2e1 + qJD(1) * (t618 * t89 - t621 * t88 + (t237 * t618 + t238 * t621) * qJD(1)) / 0.2e1 + ((t146 * t618 + t147 * t621) * qJD(1) + t770) * t812 + ((t148 * t618 + t149 * t621) * qJD(1) + t771) * t813 + t738 * t1033 + t737 * t1034; t1109 * t943 / 0.2e1 + t1110 * t942 / 0.2e1 + t1099 * t877 / 0.2e1 + (t1056 * t617 - t1060 * t620) * t1041 + (t1055 * t617 - t1061 * t620) * t1040 + (t1057 * t617 - t1059 * t620) * t1039 + (t1044 * t470 - t1058 * t621 + t1073 * t417 + t1074 * t416 + t636 * t471) * t1038 + ((qJD(2) * t1056 - t1098) * t620 + (-qJD(1) * t1089 + t1060 * qJD(2) + t1107 * t621 + t1108 * t618) * t617) * t1037 + ((qJD(2) * t1055 - t1100) * t620 + (-qJD(1) * t1088 + t1061 * qJD(2) + t1105 * t621 + t1106 * t618) * t617) * t1036 + (-t1044 * t468 - t1058 * t618 + t1073 * t415 + t1074 * t414 + t469 * t636) * t1035 + (t1133 * t620 + (-t1044 * t616 + t1073 * t597 + t1074 * t596 + t619 * t636) * t617) * t1032 + ((qJD(2) * t1057 - t1134) * t620 + (-qJD(1) * t1090 + t1059 * qJD(2) + t1103 * t621 + t1104 * t618) * t617) * t1031 - (t1011 + t997 + t998 - t1012 + t1009 + t995 + t996 - t1010 + t1007 + t993 + t994 - t1008 + t1054) * t620 / 0.2e1 + (-t52 * qJD(5) * t417 - t53 * (qJD(5) * t415 + t911) - t44 * (t597 * t867 - t343) - (t44 * t924 + t53 * t902) * t493 - (t1080 * t52 + t53 * t922) * t582 - (t52 * (t861 + t902) + t44 * t1079) * t494 + t8 * t925 + t52 * t798 + t53 * t926 - t7 * t231 + t44 * t934 + (qJD(2) * t626 + t52 * t991 - t53 * t850 + t8 * t933 - t841 * t9) * t620 + (t661 * qJD(2) + (qJD(1) * t662 + t53 * t843 + t837 * t9 + t674) * t621 + (t8 * t909 - t52 * t928 - t7 * t841 - t44 * t850 + (t44 * t842 + t53 * t909) * qJD(1)) * t618) * t617 - g(1) * t1079 + g(2) * t1080 - g(3) * t1077 - (-t1127 * t596 - t1025) * t1020) * m(6) + (t19 * t925 + t63 * t798 + t64 * t926 - t10 * t231 + t56 * t934 + (qJD(2) * t629 + t63 * t107 - t19 * t211 - t20 * t929 - t64 * t935) * t620 + (t686 * qJD(2) + (qJD(1) * t687 + t20 * t910 + t64 * t927 + t772) * t621 + (t19 * t380 + t63 * t221 - t10 * t929 - t56 * t935 + (t64 * t380 + t56 * t932) * qJD(1)) * t618) * t617 - t64 * (-t439 * t493 + t911) - t56 * (t265 * t493 - t343) - (t64 * t270 + t63 * t923) * t582 - (t63 * (-t439 + t861) + t56 * t921) * t494 - g(1) * t921 + g(2) * t923 - (t777 - t1025) * t1020) * m(5) + ((qJD(2) * t634 + t114 * t134 - t115 * t133 - t286 * t45 - t46 * t287) * t620 + (t114 * (qJD(2) * t286 + t297 * t618) + t115 * (qJD(2) * t287 - t297 * t621) + t35 * t729 + t83 * (-t133 * t618 + t134 * t621 + t286 * t880 - t287 * t878) + (t45 * t618 - t46 * t621 + (t114 * t621 + t115 * t618) * qJD(1)) * t442) * t617 - g(1) * t320 - g(2) * t319 - g(3) * t481 - t114 * (-t319 * t582 - t481 * t494) - t115 * (t320 * t582 - t481 * t493) - t83 * (t319 * t493 + t320 * t494)) * m(4) + (t617 * t816 + t620 * t813) * t1102 + (-t834 / 0.2e1 + t620 * t811) * t1101; (-m(5) - m(6)) * (-g(3) * t620 + t1049) - m(5) * (t302 * t56 + t345 * t64 + t346 * t63) - m(6) * (t302 * t44 + t345 * t53 + t346 * t52) + 0.2e1 * ((t63 * t874 + t64 * t876 - t10) * t1043 + (t52 * t874 + t53 * t876 - t7) * t1042) * t620 + 0.2e1 * ((qJD(2) * t56 + t19 * t621 + t20 * t618 - t63 * t880 + t64 * t878) * t1043 + (qJD(2) * t44 - t52 * t880 + t53 * t878 + t618 * t9 + t621 * t8) * t1042) * t617; (-t181 * t52 + t183 * t53 + t676 * t44 + (t53 * t493 + t52 * t494 - g(3) + t7) * t950 + (-t44 * t494 - t53 * t582 - g(1) + t8) * t416 + (-t44 * t493 + t52 * t582 + t1128) * t414) * m(6);];
tau = t1;
