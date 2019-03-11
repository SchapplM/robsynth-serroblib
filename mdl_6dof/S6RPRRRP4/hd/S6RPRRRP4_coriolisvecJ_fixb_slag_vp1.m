% Calculate vector of centrifugal and Coriolis load on the joints for
% S6RPRRRP4
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d3,d4,d5,theta2]';
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
% Datum: 2019-03-09 06:09
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauc = S6RPRRRP4_coriolisvecJ_fixb_slag_vp1(qJ, qJD, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(10,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRRRP4_coriolisvecJ_fixb_slag_vp1: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RPRRRP4_coriolisvecJ_fixb_slag_vp1: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RPRRRP4_coriolisvecJ_fixb_slag_vp1: pkin has to be [10x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RPRRRP4_coriolisvecJ_fixb_slag_vp1: m has to be [7x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [7,3]), ...
  'S6RPRRRP4_coriolisvecJ_fixb_slag_vp1: rSges has to be [7x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [7 6]), ...
  'S6RPRRRP4_coriolisvecJ_fixb_slag_vp1: Icges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From coriolisvec_joint_fixb_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 06:07:04
% EndTime: 2019-03-09 06:08:32
% DurationCPUTime: 76.23s
% Computational Cost: add. (58809->1266), mult. (56774->1650), div. (0->0), fcn. (52403->10), ass. (0->642)
t1188 = -Icges(6,4) - Icges(7,4);
t1148 = Icges(6,1) + Icges(7,1);
t1160 = Icges(6,5) + Icges(7,5);
t1159 = -Icges(6,2) - Icges(7,2);
t1158 = Icges(6,6) + Icges(7,6);
t617 = cos(qJ(5));
t1187 = t1188 * t617;
t615 = sin(qJ(5));
t1186 = t1188 * t615;
t1185 = -t1158 * t615 + t1160 * t617;
t1184 = t1159 * t615 - t1187;
t1183 = t1148 * t617 + t1186;
t1182 = Icges(6,3) + Icges(7,3);
t607 = pkin(10) + qJ(3);
t585 = qJ(4) + t607;
t568 = sin(t585);
t569 = cos(t585);
t1162 = -t1158 * t569 + t1184 * t568;
t1151 = -t1160 * t569 + t1183 * t568;
t1181 = t1185 * t569;
t1180 = t1184 * t569;
t1179 = t1183 * t569;
t1178 = -t1158 * t617 - t1160 * t615;
t1177 = t1159 * t617 + t1186;
t1176 = -t1148 * t615 + t1187;
t1153 = -t1182 * t569 + t1185 * t568;
t618 = cos(qJ(1));
t932 = t618 * t615;
t616 = sin(qJ(1));
t934 = t616 * t617;
t486 = -t569 * t932 + t934;
t933 = t617 * t618;
t935 = t615 * t616;
t487 = t569 * t933 + t935;
t947 = t568 * t618;
t1084 = t1151 * t487 + t1153 * t947 + t1162 * t486;
t270 = Icges(7,5) * t487 + Icges(7,6) * t486 + Icges(7,3) * t947;
t273 = Icges(6,5) * t487 + Icges(6,6) * t486 + Icges(6,3) * t947;
t1154 = t270 + t273;
t464 = Icges(7,4) * t486;
t282 = Icges(7,1) * t487 + Icges(7,5) * t947 + t464;
t467 = Icges(6,4) * t486;
t285 = Icges(6,1) * t487 + Icges(6,5) * t947 + t467;
t1163 = t282 + t285;
t995 = Icges(7,4) * t487;
t276 = Icges(7,2) * t486 + Icges(7,6) * t947 + t995;
t998 = Icges(6,4) * t487;
t279 = Icges(6,2) * t486 + Icges(6,6) * t947 + t998;
t1164 = t276 + t279;
t1120 = t1154 * t947 + t1163 * t487 + t1164 * t486;
t484 = t569 * t935 + t933;
t462 = Icges(7,4) * t484;
t485 = t569 * t934 - t932;
t948 = t568 * t616;
t281 = -Icges(7,1) * t485 - Icges(7,5) * t948 + t462;
t465 = Icges(6,4) * t484;
t284 = -Icges(6,1) * t485 - Icges(6,5) * t948 + t465;
t1080 = t281 + t284;
t463 = Icges(7,4) * t485;
t274 = -Icges(7,2) * t484 + Icges(7,6) * t948 + t463;
t466 = Icges(6,4) * t485;
t277 = -Icges(6,2) * t484 + Icges(6,6) * t948 + t466;
t1165 = t274 + t277;
t268 = Icges(7,5) * t485 - Icges(7,6) * t484 + Icges(7,3) * t948;
t271 = Icges(6,5) * t485 - Icges(6,6) * t484 + Icges(6,3) * t948;
t1166 = t268 + t271;
t1121 = -t1080 * t487 + t1165 * t486 + t1166 * t947;
t587 = qJD(3) * t616;
t544 = qJD(4) * t616 + t587;
t871 = qJD(5) * t618;
t459 = t568 * t871 + t544;
t608 = qJD(3) + qJD(4);
t545 = t608 * t618;
t872 = qJD(5) * t616;
t460 = -t568 * t872 + t545;
t874 = qJD(5) * t569;
t541 = qJD(1) - t874;
t1126 = t1084 * t541 + t1120 * t459 - t1121 * t460;
t1085 = t1151 * t485 + t1153 * t948 - t1162 * t484;
t1122 = t1154 * t948 + t1163 * t485 - t1164 * t484;
t1123 = -t1080 * t485 - t1165 * t484 + t1166 * t948;
t1127 = t1085 * t541 + t1122 * t459 - t1123 * t460;
t787 = qJD(1) * t569 - qJD(5);
t854 = t568 * t545;
t1062 = t616 * t787 + t854;
t714 = t617 * t541;
t249 = t1062 * t615 + t618 * t714;
t713 = t541 * t615;
t250 = -t1062 * t617 + t618 * t713;
t878 = qJD(1) * t616;
t837 = t568 * t878;
t853 = t569 * t545;
t674 = -t837 + t853;
t143 = Icges(7,5) * t250 + Icges(7,6) * t249 + Icges(7,3) * t674;
t145 = Icges(6,5) * t250 + Icges(6,6) * t249 + Icges(6,3) * t674;
t1175 = t143 + t145;
t937 = t608 * t616;
t855 = t568 * t937;
t251 = t616 * t714 + (-t618 * t787 + t855) * t615;
t950 = t568 * t608;
t252 = t787 * t933 + (-t617 * t950 + t713) * t616;
t877 = qJD(1) * t618;
t675 = t568 * t877 + t569 * t937;
t144 = Icges(7,5) * t252 + Icges(7,6) * t251 + Icges(7,3) * t675;
t146 = Icges(6,5) * t252 + Icges(6,6) * t251 + Icges(6,3) * t675;
t1174 = t144 + t146;
t147 = Icges(7,4) * t250 + Icges(7,2) * t249 + Icges(7,6) * t674;
t149 = Icges(6,4) * t250 + Icges(6,2) * t249 + Icges(6,6) * t674;
t1173 = t147 + t149;
t148 = Icges(7,4) * t252 + Icges(7,2) * t251 + Icges(7,6) * t675;
t150 = Icges(6,4) * t252 + Icges(6,2) * t251 + Icges(6,6) * t675;
t1172 = t148 + t150;
t151 = Icges(7,1) * t250 + Icges(7,4) * t249 + Icges(7,5) * t674;
t153 = Icges(6,1) * t250 + Icges(6,4) * t249 + Icges(6,5) * t674;
t1171 = t151 + t153;
t152 = Icges(7,1) * t252 + Icges(7,4) * t251 + Icges(7,5) * t675;
t154 = Icges(6,1) * t252 + Icges(6,4) * t251 + Icges(6,5) * t675;
t1170 = t152 + t154;
t1169 = t1181 * t608 + (qJD(5) * t1178 + t1182 * t608) * t568;
t1168 = t1180 * t608 + (qJD(5) * t1177 + t1158 * t608) * t568;
t1167 = t1179 * t608 + (t1176 * qJD(5) + t1160 * t608) * t568;
t1161 = t1151 * t617 - t1162 * t615;
t1133 = -t1080 * t250 + t1165 * t249 + t1166 * t674 + t1170 * t487 + t1172 * t486 + t1174 * t947;
t1132 = t1154 * t674 + t1163 * t250 + t1164 * t249 + t1171 * t487 + t1173 * t486 + t1175 * t947;
t1131 = -t1080 * t252 + t1165 * t251 + t1166 * t675 + t1170 * t485 - t1172 * t484 + t1174 * t948;
t1130 = t1154 * t675 + t1163 * t252 + t1164 * t251 + t1171 * t485 - t1173 * t484 + t1175 * t948;
t1095 = t1151 * t250 + t1153 * t674 + t1162 * t249 + t1167 * t487 + t1168 * t486 + t1169 * t947;
t1157 = t1151 * t252 + t1153 * t675 + t1162 * t251 + t1167 * t485 - t1168 * t484 + t1169 * t948;
t728 = -t274 * t615 - t281 * t617;
t120 = -t268 * t569 + t568 * t728;
t726 = -t277 * t615 - t284 * t617;
t122 = -t271 * t569 + t568 * t726;
t1156 = t120 + t122;
t727 = -t276 * t615 + t282 * t617;
t121 = -t270 * t569 + t568 * t727;
t725 = -t279 * t615 + t285 * t617;
t123 = -t273 * t569 + t568 * t725;
t1155 = t121 + t123;
t1083 = -t1153 * t569 + t1161 * t568;
t613 = -qJ(6) - pkin(9);
t1023 = pkin(9) + t613;
t579 = pkin(5) * t617 + pkin(4);
t1024 = pkin(4) - t579;
t761 = rSges(7,1) * t617 - rSges(7,2) * t615;
t915 = (-rSges(7,3) + t1023) * t569 + (-t1024 + t761) * t568;
t1150 = (t1182 * t568 - t1161 + t1181) * t541 + (t1153 * t616 + t726 + t728) * t460 + (-t1153 * t618 - t725 - t727) * t459;
t1076 = (t1159 * t485 - t1080 - t462 - t465) * t460 - (t1159 * t487 + t1163 + t464 + t467) * t459 - (t1177 * t568 + t1151) * t541;
t1149 = (t1161 * t608 - t1169) * t569 + (t1167 * t617 - t1168 * t615 + t1153 * t608 + (-t1151 * t615 - t1162 * t617) * qJD(5)) * t568;
t945 = t569 * t618;
t1147 = t487 * rSges(7,1) + t486 * rSges(7,2) + rSges(7,3) * t947 + pkin(5) * t935 + t579 * t945;
t1065 = t1023 * t568;
t348 = -t1024 * t569 - t1065;
t1146 = t761 * t569 + t348;
t946 = t569 * t616;
t1145 = -rSges(7,1) * t485 + rSges(7,2) * t484 - t579 * t946;
t1144 = (-t1158 * t485 - t1160 * t484) * t460 + (t1158 * t487 - t1160 * t486) * t459 - t1178 * t541 * t568;
t868 = qJD(6) * t618;
t528 = t568 * t868;
t1012 = pkin(5) * qJD(5);
t856 = t617 * t1012;
t1028 = pkin(5) * t615;
t858 = qJD(1) * t1028;
t1143 = t250 * rSges(7,1) + t249 * rSges(7,2) + rSges(7,3) * t853 + t613 * t837 + t616 * t856 + t618 * t858 + t528;
t1029 = pkin(4) * t569;
t927 = pkin(5) * t932 + (t1029 + t1065) * t616 - rSges(7,3) * t948 + t1145;
t1139 = t541 * t927 + t528;
t864 = qJD(1) * qJD(3);
t573 = t616 * t864;
t863 = qJD(1) * qJD(4);
t520 = t616 * t863 + t573;
t341 = qJD(5) * t675 + t520;
t574 = t618 * t864;
t521 = t618 * t863 + t574;
t342 = qJD(5) * t674 + t521;
t873 = qJD(5) * t608;
t830 = t568 * t873;
t1136 = t1084 * t830 + t1095 * t541 + t1120 * t342 + t1121 * t341 + t1132 * t459 - t1133 * t460;
t1135 = t1085 * t830 + t1122 * t342 + t1123 * t341 + t1130 * t459 - t1131 * t460 + t1157 * t541;
t1134 = rSges(7,1) + pkin(5);
t37 = (t608 * t728 - t144) * t569 + (-t148 * t615 + t152 * t617 + t268 * t608 + (-t274 * t617 + t281 * t615) * qJD(5)) * t568;
t39 = (t608 * t726 - t146) * t569 + (-t150 * t615 + t154 * t617 + t271 * t608 + (-t277 * t617 + t284 * t615) * qJD(5)) * t568;
t1129 = t37 + t39;
t38 = (t608 * t727 - t143) * t569 + (-t147 * t615 + t151 * t617 + t270 * t608 + (-t276 * t617 - t282 * t615) * qJD(5)) * t568;
t40 = (t608 * t725 - t145) * t569 + (-t149 * t615 + t153 * t617 + t273 * t608 + (-t279 * t617 - t285 * t615) * qJD(5)) * t568;
t1128 = t38 + t40;
t1125 = t1083 * t541 + t1155 * t459 - t1156 * t460;
t517 = pkin(4) * t855;
t763 = rSges(7,1) * t252 + rSges(7,2) * t251;
t944 = t579 * t608;
t789 = -qJD(6) + t944;
t857 = t615 * t1012;
t931 = t517 + (qJD(1) * t348 - t856) * t618 + (t858 - t789 * t568 + (-t1023 * t608 - t857) * t569) * t616 + rSges(7,3) * t675 + t763;
t518 = pkin(9) * t853;
t938 = t608 * t613;
t930 = -rSges(7,3) * t837 - t518 + (pkin(9) * t878 + t1024 * t545) * t568 + ((-t857 - t938) * t618 + t1024 * t878) * t569 + t1143;
t543 = pkin(4) * t945;
t458 = pkin(9) * t947 + t543;
t926 = -t613 * t947 + t1147 - t458;
t1119 = t915 * t616;
t1118 = t915 * t618;
t1015 = rSges(7,3) * t568;
t1117 = -t1015 - t1146;
t1116 = t1162 * t616;
t1115 = t1162 * t618;
t1114 = t1151 * t616;
t1113 = t1151 * t618;
t1112 = t1158 * t568 + t1180;
t1111 = t1160 * t568 + t1179;
t985 = Icges(5,6) * t618;
t412 = Icges(5,4) * t946 - Icges(5,2) * t948 - t985;
t561 = Icges(5,4) * t569;
t497 = Icges(5,1) * t568 + t561;
t1110 = -t497 * t616 - t412;
t752 = -Icges(5,2) * t568 + t561;
t413 = Icges(5,6) * t616 + t618 * t752;
t1109 = -t497 * t618 - t413;
t999 = Icges(5,4) * t568;
t498 = Icges(5,1) * t569 - t999;
t415 = Icges(5,5) * t616 + t498 * t618;
t495 = Icges(5,2) * t569 + t999;
t1108 = -t495 * t618 + t415;
t1107 = -t495 + t498;
t1106 = t497 + t752;
t1105 = t1150 * t568;
t1104 = t1153 * t541 + t1154 * t459 - t1166 * t460;
t1073 = t1155 * t618 + t1156 * t616;
t1103 = -t1155 * t616 + t1156 * t618;
t1072 = t1120 * t618 + t1121 * t616;
t1102 = -t1120 * t616 + t1121 * t618;
t1071 = t1122 * t618 + t1123 * t616;
t1101 = -t1122 * t616 + t1123 * t618;
t502 = pkin(9) * t568 + t1029;
t456 = t502 * t616;
t612 = cos(pkin(10));
t570 = t612 * pkin(2) + pkin(1);
t584 = cos(t607);
t522 = pkin(3) * t584 + t570;
t614 = -pkin(7) - qJ(2);
t606 = -pkin(8) + t614;
t895 = -t616 * t522 - t618 * t606;
t936 = t614 * t618;
t343 = t570 * t616 + t895 + t936;
t504 = t618 * t522;
t553 = t618 * t570;
t883 = -t606 + t614;
t344 = t616 * t883 + t504 - t553;
t876 = qJD(3) * t618;
t919 = -t343 * t587 + t344 * t876;
t709 = t544 * t456 + t458 * t545 + t919;
t870 = qJD(6) * t569;
t66 = -t459 * t927 + t460 * t926 + t709 - t870;
t501 = pkin(4) * t568 - pkin(9) * t569;
t792 = -t614 * t616 + t553;
t842 = t344 + t792;
t1013 = pkin(3) * qJD(3);
t583 = sin(t607);
t860 = t583 * t1013;
t535 = t616 * t860;
t589 = qJD(2) * t618;
t888 = t535 + t589;
t659 = (t458 + t842) * qJD(1) - t544 * t501 - t888;
t869 = qJD(6) * t616;
t827 = t568 * t869;
t79 = -t459 * t915 + t541 * t926 + t659 + t827;
t809 = t79 * t915;
t1100 = t66 * t927 + t809;
t766 = rSges(6,1) * t485 - rSges(6,2) * t484;
t287 = rSges(6,3) * t948 + t766;
t765 = rSges(6,1) * t617 - rSges(6,2) * t615;
t398 = -rSges(6,3) * t569 + t568 * t765;
t1099 = -t287 * t541 - t398 * t460;
t1098 = 0.2e1 * qJD(3);
t1093 = rSges(4,2) * t583;
t527 = t569 * t869;
t676 = -t569 * t878 - t854;
t294 = pkin(4) * t676 - pkin(9) * t837 + t518;
t438 = t502 * t608;
t833 = t583 * t876;
t788 = pkin(3) * t833;
t891 = t522 - t570;
t263 = -t788 + (-t616 * t891 + t618 * t883) * qJD(1);
t577 = qJ(2) * t877;
t1025 = pkin(1) - t570;
t691 = t1025 * t616 - t936;
t865 = qJD(1) * qJD(2);
t588 = qJD(2) * t616;
t884 = t577 + t588;
t897 = qJD(1) * (-pkin(1) * t878 + t884) + t616 * t865;
t840 = qJD(1) * (qJD(1) * t691 - t577) + t897;
t939 = t584 * qJD(3) ^ 2;
t633 = qJD(1) * t263 + (-t574 * t583 - t616 * t939) * pkin(3) + t840;
t628 = qJD(1) * t294 - t544 * t438 - t521 * t501 + t633;
t760 = -rSges(7,1) * t615 - rSges(7,2) * t617;
t929 = -t870 + t1146 * t608 + (rSges(7,3) * t608 + qJD(5) * t760 - t857) * t568;
t24 = t608 * t527 + t930 * t541 - t929 * t459 - t915 * t342 + (qJD(1) * t868 + t873 * t926) * t568 + t628;
t1092 = t24 * t616;
t529 = t569 * t868;
t576 = t618 * t865;
t694 = -pkin(3) * t618 * t939 + qJD(1) * t535 + t576;
t664 = -t545 * t438 + t520 * t501 + t694;
t295 = pkin(9) * t675 + qJD(1) * t543 - t517;
t557 = t606 * t878;
t562 = t614 * t878;
t264 = t877 * t891 - t535 - t557 + t562;
t590 = t616 * qJ(2);
t548 = t618 * pkin(1) + t590;
t491 = qJD(1) * t548 - t589;
t909 = t562 - (-t1025 * t618 - t590) * qJD(1) - t491;
t846 = -t264 + t909;
t786 = -t295 + t846;
t806 = t927 * t568;
t25 = -t931 * t541 - t929 * t460 + t915 * t341 + (qJD(5) * t806 + t529) * t608 + (t786 - t827) * qJD(1) + t664;
t1090 = qJD(1) * t79 + t25;
t592 = t618 * qJ(2);
t419 = -t592 + t691;
t546 = pkin(1) * t616 - t592;
t526 = qJD(1) * t546;
t1079 = qJD(1) * t419 - t526;
t600 = t616 * rSges(4,3);
t940 = t584 * t618;
t942 = t583 * t618;
t434 = rSges(4,1) * t940 - rSges(4,2) * t942 + t600;
t1078 = t434 + t792;
t1019 = rSges(3,2) * sin(pkin(10));
t1022 = rSges(3,1) * t612;
t711 = t616 * rSges(3,3) + (-t1019 + t1022) * t618;
t1077 = t548 + t711;
t567 = Icges(4,4) * t584;
t753 = -Icges(4,2) * t583 + t567;
t1075 = (t1176 * t568 - t1162) * t541 + (t1148 * t484 + t1165 + t463 + t466) * t460 + (t1148 * t486 - t1164 - t995 - t998) * t459;
t1074 = t1144 * t568;
t362 = t398 * t616;
t1017 = rSges(6,3) * t568;
t693 = t765 * t569;
t400 = t693 + t1017;
t461 = t501 * t878;
t455 = t501 * t616;
t790 = qJD(1) * t455 - t545 * t502;
t829 = t569 * t872;
t875 = qJD(5) * t568;
t1068 = t287 * t875 - t362 * t541 + t460 * t400 + t461 - t790 + (-t829 + t878) * t398;
t511 = Icges(4,1) * t583 + t567;
t1064 = t1083 * t830 + t1149 * t541;
t416 = rSges(5,1) * t946 - rSges(5,2) * t948 - t618 * rSges(5,3);
t599 = t616 * rSges(5,3);
t417 = rSges(5,1) * t945 - rSges(5,2) * t947 + t599;
t126 = t416 * t544 + t417 * t545 + t919;
t712 = t588 - t788;
t902 = t419 - t546;
t843 = t343 + t902;
t499 = rSges(5,1) * t568 + rSges(5,2) * t569;
t953 = t499 * t545;
t135 = -t953 + (-t416 + t843) * qJD(1) + t712;
t136 = -t499 * t544 + (t417 + t842) * qJD(1) - t888;
t453 = t499 * t616;
t454 = t499 * t618;
t1020 = rSges(5,1) * t569;
t500 = -rSges(5,2) * t568 + t1020;
t1063 = -t135 * (qJD(1) * t453 - t545 * t500) - t126 * (-t544 * t453 - t454 * t545) - t136 * (-qJD(1) * t454 - t500 * t544);
t941 = t584 * t616;
t943 = t583 * t616;
t982 = Icges(4,3) * t618;
t423 = Icges(4,5) * t941 - Icges(4,6) * t943 - t982;
t552 = Icges(4,4) * t943;
t992 = Icges(4,5) * t618;
t427 = Icges(4,1) * t941 - t552 - t992;
t986 = Icges(4,6) * t618;
t425 = Icges(4,4) * t941 - Icges(4,2) * t943 - t986;
t957 = t425 * t583;
t719 = -t427 * t584 + t957;
t173 = -t423 * t618 - t616 * t719;
t494 = Icges(5,5) * t569 - Icges(5,6) * t568;
t493 = Icges(5,5) * t568 + Icges(5,6) * t569;
t954 = t493 * t618;
t959 = t413 * t568;
t981 = Icges(5,3) * t618;
t1061 = -t608 * t954 + (-t415 * t569 - t494 * t616 + t959 + t981) * qJD(1);
t534 = Icges(5,4) * t948;
t991 = Icges(5,5) * t618;
t414 = Icges(5,1) * t946 - t534 - t991;
t721 = t412 * t568 - t414 * t569;
t411 = Icges(5,3) * t616 + t494 * t618;
t881 = qJD(1) * t411;
t955 = t493 * t616;
t1060 = qJD(1) * t721 - t608 * t955 + t881;
t508 = Icges(4,5) * t584 - Icges(4,6) * t583;
t507 = Icges(4,5) * t583 + Icges(4,6) * t584;
t679 = qJD(3) * t507;
t1000 = Icges(4,4) * t583;
t512 = Icges(4,1) * t584 - t1000;
t428 = Icges(4,5) * t616 + t512 * t618;
t426 = Icges(4,6) * t616 + t618 * t753;
t956 = t426 * t583;
t718 = -t428 * t584 + t956;
t1059 = -t618 * t679 + (-t508 * t616 + t718 + t982) * qJD(1);
t424 = Icges(4,3) * t616 + t508 * t618;
t880 = qJD(1) * t424;
t1058 = qJD(1) * t719 - t616 * t679 + t880;
t717 = t495 * t568 - t497 * t569;
t1057 = qJD(1) * t717 + t494 * t608;
t509 = Icges(4,2) * t584 + t1000;
t716 = t509 * t583 - t511 * t584;
t1056 = qJD(1) * t716 + t508 * qJD(3);
t1055 = t616 * (-t509 * t618 + t428) - t618 * (-Icges(4,2) * t941 + t427 - t552);
t291 = t487 * rSges(6,1) + t486 * rSges(6,2) + rSges(6,3) * t947;
t103 = t291 * t541 - t398 * t459 + t659;
t849 = t250 * rSges(6,1) + t249 * rSges(6,2) + rSges(6,3) * t853;
t156 = -rSges(6,3) * t837 + t849;
t764 = -rSges(6,1) * t615 - rSges(6,2) * t617;
t207 = t608 * t693 + (rSges(6,3) * t608 + qJD(5) * t764) * t568;
t41 = t156 * t541 - t207 * t459 + t291 * t830 - t342 * t398 + t628;
t767 = rSges(6,1) * t252 + rSges(6,2) * t251;
t158 = rSges(6,3) * t675 + t767;
t42 = qJD(1) * t786 - t158 * t541 - t207 * t460 - t287 * t830 + t341 * t398 + t664;
t1054 = (qJD(1) * t103 + t42) * t618 + t41 * t616;
t1051 = qJD(1) * t1106 + t544 * t1108 - t545 * (-Icges(5,2) * t946 + t414 - t534);
t1050 = -m(7) / 0.2e1;
t1049 = m(7) / 0.2e1;
t1048 = t341 / 0.2e1;
t1047 = t342 / 0.2e1;
t1046 = -t459 / 0.2e1;
t1045 = t459 / 0.2e1;
t1044 = -t460 / 0.2e1;
t1043 = t460 / 0.2e1;
t1042 = t520 / 0.2e1;
t1041 = t521 / 0.2e1;
t1040 = -t541 / 0.2e1;
t1039 = t541 / 0.2e1;
t1038 = -t544 / 0.2e1;
t1037 = t544 / 0.2e1;
t1036 = -t545 / 0.2e1;
t1035 = t545 / 0.2e1;
t1033 = t616 / 0.2e1;
t1032 = -t618 / 0.2e1;
t1031 = -rSges(6,3) - pkin(9);
t1030 = pkin(3) * t583;
t1027 = -qJD(1) / 0.2e1;
t1026 = qJD(1) / 0.2e1;
t1021 = rSges(4,1) * t584;
t1018 = rSges(4,2) * t584;
t1011 = t37 * t460;
t1010 = t38 * t459;
t1009 = t39 * t460;
t1008 = t40 * t459;
t1007 = t41 * t618;
t1006 = t42 * t616;
t1005 = -rSges(7,3) + t613;
t978 = qJD(1) * t66;
t95 = t287 * t459 + t291 * t460 + t709;
t976 = qJD(1) * t95;
t975 = t103 * t616;
t970 = t120 * t341;
t969 = t121 * t342;
t968 = t122 * t341;
t967 = t123 * t342;
t966 = t135 * t616;
t887 = rSges(4,2) * t943 + t618 * rSges(4,3);
t433 = rSges(4,1) * t941 - t887;
t513 = rSges(4,1) * t583 + t1018;
t831 = t513 * t876;
t770 = t588 - t831;
t177 = (-t433 + t902) * qJD(1) + t770;
t963 = t177 * t616;
t834 = t513 * t587;
t178 = qJD(1) * t1078 - t589 - t834;
t482 = t513 * t618;
t962 = t178 * t482;
t410 = Icges(5,5) * t946 - Icges(5,6) * t948 - t981;
t960 = t410 * t618;
t952 = t507 * t616;
t951 = t507 * t618;
t408 = t545 * t501;
t949 = t568 * t613;
t928 = -t207 - t438;
t921 = -t291 - t458;
t918 = -t616 * t343 + t618 * t344;
t917 = -rSges(7,2) * t485 - t1134 * t484;
t916 = -rSges(7,2) * t487 + t1134 * t486;
t914 = -t616 * t410 - t414 * t945;
t913 = t616 * t411 + t415 * t945;
t911 = -t616 * t423 - t427 * t940;
t910 = t616 * t424 + t428 * t940;
t905 = t616 * t416 + t618 * t417;
t900 = t616 * t456 + t618 * t458;
t893 = -t509 + t512;
t892 = t511 + t753;
t890 = rSges(5,2) * t837 + rSges(5,3) * t877;
t889 = rSges(4,3) * t877 + t1093 * t878;
t563 = t616 * t1019;
t886 = rSges(3,3) * t877 + qJD(1) * t563;
t885 = t618 * rSges(3,3) + t563;
t879 = qJD(1) * t508;
t560 = qJD(6) * t568;
t180 = -t616 * t717 - t954;
t867 = t180 * qJD(1);
t209 = -t616 * t716 - t951;
t866 = t209 * qJD(1);
t862 = t616 * t1022;
t859 = t584 * t1013;
t852 = -t438 - t929;
t247 = rSges(5,1) * t676 - rSges(5,2) * t853 + t890;
t248 = -t608 * t453 + (t500 * t618 + t599) * qJD(1);
t851 = t618 * t247 + t616 * t248 + t416 * t877;
t848 = t618 * t263 + t616 * t264 - t343 * t877;
t847 = -t458 - t926;
t845 = t618 * t294 + t616 * t295 + t456 * t877;
t844 = t878 * t915 + t461;
t839 = t557 + t888;
t838 = t1031 * t568;
t835 = t583 * t877;
t828 = t569 * t871;
t822 = t569 * t608 / 0.2e1;
t821 = -pkin(1) - t1022;
t820 = t878 / 0.2e1;
t819 = t877 / 0.2e1;
t818 = -t587 / 0.2e1;
t815 = t876 / 0.2e1;
t812 = -t499 - t1030;
t811 = -t501 - t1030;
t810 = t66 * t926;
t807 = -t570 - t1021;
t805 = (-t616 ^ 2 - t618 ^ 2) * t583;
t804 = (-t616 * t752 + t985) * qJD(1) + t1108 * t608;
t803 = qJD(1) * t413 + t414 * t608 - t495 * t937;
t802 = (-t498 * t616 + t991) * qJD(1) + t1109 * t608;
t801 = qJD(1) * t415 + t1110 * t608;
t365 = t415 * t946;
t800 = t411 * t618 - t365;
t374 = t428 * t941;
t799 = t424 * t618 - t374;
t457 = t501 * t618;
t798 = -t544 * t455 - t457 * t545;
t797 = -t410 + t959;
t796 = -t423 + t956;
t795 = t1106 * t608;
t794 = t1107 * t608;
t793 = -qJD(1) * t457 - t502 * t544;
t785 = t616 * t287 + t618 * t291 + t900;
t776 = (t1028 - t760) * t568;
t435 = t500 * t608;
t772 = -t435 - t859;
t771 = -t438 - t859;
t768 = t1021 - t1093;
t640 = (-t456 + t843) * qJD(1) - t408 + t712;
t102 = t1099 + t640;
t743 = t102 * t618 + t975;
t730 = -t135 * t618 - t136 * t616;
t729 = -t177 * t618 - t178 * t616;
t724 = t287 * t618 - t291 * t616;
t196 = t412 * t569 + t414 * t568;
t225 = t425 * t584 + t427 * t583;
t226 = t426 * t584 + t428 * t583;
t715 = t811 - t915;
t710 = -t207 + t771;
t707 = t618 * t156 + t616 * t158 + t287 * t877 + t845;
t706 = -t616 * t927 + t618 * t926 + t900;
t705 = -t569 * t579 - t1015 - t522;
t481 = t513 * t616;
t690 = t771 - t929;
t682 = t721 * t616;
t681 = qJD(3) * t511;
t680 = qJD(3) * t509;
t174 = -t426 * t943 - t799;
t678 = (-t173 * t618 + t174 * t616) * qJD(3);
t175 = -t425 * t942 - t911;
t176 = -t426 * t942 + t910;
t677 = (-t175 * t618 + t176 * t616) * qJD(3);
t253 = (t433 * t616 + t434 * t618) * qJD(3);
t668 = qJD(1) * t494 - t544 * t954 + t545 * t955;
t667 = t263 * t876 + t264 * t587 - t343 * t574 - t344 * t573;
t78 = -t460 * t915 + t1139 + t640;
t666 = t78 * t915 - t810;
t665 = t425 * t618 - t426 * t616;
t663 = qJD(1) * t343 + t1079 + t712;
t662 = t616 * t931 + t618 * t930 - t877 * t927 + t845;
t650 = (-t583 * t892 + t584 * t893) * qJD(1);
t647 = -qJD(1) * t456 - t408 + t663;
t641 = qJD(1) * t1107 + t1109 * t544 - t1110 * t545;
t639 = t1117 * t460 + t1119 * t541 + t829 * t915 + t529 + t790;
t636 = qJD(1) * t410 - t568 * t803 + t569 * t801;
t635 = -t568 * t804 + t569 * t802 + t881;
t634 = qJD(1) * t493 - t568 * t795 + t569 * t794;
t632 = t294 * t545 + t544 * t295 + t521 * t456 - t458 * t520 + t667;
t301 = qJD(1) * t426 - t616 * t680;
t303 = qJD(1) * t428 - t616 * t681;
t631 = qJD(1) * t423 - qJD(3) * t225 - t301 * t583 + t303 * t584;
t300 = -t618 * t680 + (-t616 * t753 + t986) * qJD(1);
t302 = -t618 * t681 + (-t512 * t616 + t992) * qJD(1);
t630 = -qJD(3) * t226 - t300 * t583 + t302 * t584 + t880;
t489 = t753 * qJD(3);
t490 = t512 * qJD(3);
t629 = qJD(1) * t507 - t489 * t583 + t490 * t584 + (-t509 * t584 - t511 * t583) * qJD(3);
t627 = -t1055 * t583 + t665 * t584;
t624 = -t1051 * t568 + t641 * t569;
t364 = t398 * t618;
t623 = t103 * (t291 * t875 - t541 * t364 - t398 * t828 - t400 * t459 + t793) + t95 * (t287 * t828 - t291 * t829 - t459 * t362 - t364 * t460 + t798);
t104 = t1057 * t616 + t634 * t618;
t105 = -t1057 * t618 + t634 * t616;
t116 = t568 * t801 + t569 * t803;
t117 = t568 * t802 + t569 * t804;
t167 = -t682 - t960;
t168 = -t413 * t948 - t800;
t169 = -t412 * t947 - t914;
t170 = -t413 * t947 + t913;
t197 = t413 * t569 + t415 * t568;
t82 = t1060 * t616 + t636 * t618;
t83 = t1061 * t616 + t635 * t618;
t84 = -t1060 * t618 + t636 * t616;
t85 = -t1061 * t618 + t635 * t616;
t93 = -t167 * t545 + t168 * t544 + t867;
t181 = -t618 * t717 + t955;
t179 = t181 * qJD(1);
t94 = -t169 * t545 + t170 * t544 + t179;
t622 = (t616 * t624 - t618 * t668) * t1035 + (t616 * t85 - t618 * t84 + (t167 * t616 + t168 * t618) * qJD(1)) * t1036 + (t616 * t83 - t618 * t82 + (t169 * t616 + t170 * t618) * qJD(1)) * t1037 + (t616 * t668 + t618 * t624) * t1038 + (-t169 * t618 + t170 * t616) * t1041 + (-t167 * t618 + t168 * t616) * t1042 + (-t116 * t618 + t117 * t616 + (t196 * t616 + t197 * t618) * qJD(1)) * t1026 + (t1051 * t569 + t641 * t568) * t1027 - t1101 * t341 / 0.2e1 - t1102 * t342 / 0.2e1 + ((t1084 * t568 + t1121 * t946) * qJD(5) + ((qJD(5) * t1120 + t1104) * t569 + t1105) * t618 + (t1111 * t487 + t1112 * t486) * t541 + (t1114 * t487 + t1116 * t486) * t460 + (-t1113 * t487 - t1115 * t486) * t459) * t1046 + (qJD(1) * t1072 + t1132 * t616 - t1133 * t618) * t1045 + (qJD(1) * t1071 + t1130 * t616 - t1131 * t618) * t1044 + ((t1085 * t568 + t1122 * t945) * qJD(5) + ((qJD(5) * t1123 + t1104) * t569 + t1105) * t616 + (t1111 * t485 - t1112 * t484) * t541 + (t1114 * t485 - t1116 * t484) * t460 + (-t1113 * t485 + t1115 * t484) * t459) * t1043 + ((qJD(5) * t1073 - t1150) * t569 + ((t1111 * t617 - t1112 * t615 + t1153) * t541 + (t1114 * t617 - t1116 * t615 - t1166) * t460 + (-t1113 * t617 + t1115 * t615 + t1154) * t459 + t1083 * qJD(5)) * t568) * t1040 + (qJD(1) * t1073 + t1128 * t616 - t1129 * t618) * t1039 + (t93 + t1127) * t820 + (t94 + t1126) * t819 + (qJD(1) * t104 + t169 * t520 + t170 * t521 + t544 * t83 - t545 * t82 + t1136) * t1033 + (qJD(1) * t105 + t167 * t520 + t168 * t521 + t544 * t85 - t545 * t84 + t1135) * t1032 - (t1103 * t608 + t1125) * t875 / 0.2e1 - (t1126 * t618 + t1127 * t616) * t874 / 0.2e1;
t621 = t66 * (-t1118 * t460 - t1119 * t459 - t828 * t927 + t560 + t798) + t79 * (t1117 * t459 - t1118 * t541 + t875 * t926 + t527 + t793) + (t78 * t806 + (-t616 * t810 - t618 * t809) * t569) * qJD(5);
t492 = t768 * qJD(3);
t474 = t862 - t885;
t452 = t764 * t568;
t346 = qJD(1) * t1077 - t589;
t345 = t588 + (-t474 - t546) * qJD(1);
t337 = rSges(6,1) * t486 - rSges(6,2) * t487;
t335 = -rSges(6,1) * t484 - rSges(6,2) * t485;
t305 = -qJD(3) * t481 + (t618 * t768 + t600) * qJD(1);
t304 = -t876 * t1018 + (-t584 * t878 - t833) * rSges(4,1) + t889;
t293 = t576 + (-qJD(1) * t711 - t491) * qJD(1);
t292 = qJD(1) * (-qJD(1) * t862 + t886) + t897;
t210 = -t618 * t716 + t952;
t199 = t210 * qJD(1);
t134 = -t492 * t876 + t576 + (-t305 + t834 + t909) * qJD(1);
t133 = -t492 * t587 + (t304 - t831) * qJD(1) + t840;
t125 = -qJD(3) * t718 + t300 * t584 + t302 * t583;
t124 = -qJD(3) * t719 + t301 * t584 + t303 * t583;
t119 = -t1056 * t618 + t629 * t616;
t118 = t1056 * t616 + t629 * t618;
t101 = t199 + t677;
t100 = t678 + t866;
t99 = -t435 * t545 + t499 * t520 + (-t248 + t846) * qJD(1) + t694;
t98 = qJD(1) * t247 - t435 * t544 - t499 * t521 + t633;
t75 = t247 * t545 + t248 * t544 + t416 * t521 - t417 * t520 + t667;
t28 = t156 * t460 + t158 * t459 + t287 * t342 - t291 * t341 + t632;
t19 = -t341 * t926 - t342 * t927 + t459 * t931 + t460 * t930 + t560 * t608 + t632;
t1 = [-(t124 + t119 + t101) * t876 / 0.2e1 + (t125 + t118) * t587 / 0.2e1 + t1095 * t1045 + (t1166 * t569 + (t1080 * t617 + t1165 * t615) * t568 + t1156) * t459 * t1040 + (t180 + t196) * t1042 + (t181 + t197) * t1041 + (-(-qJD(1) * t416 - t135 + t663 - t953) * t136 + t99 * (-t416 + t895) + t135 * t839 + t98 * (-t606 * t616 + t417 + t504) + t136 * (t712 + t890) + (-t136 * t454 + t499 * t966) * t608 + ((-t135 * rSges(5,3) + t136 * (-t522 - t1020)) * t616 + (t135 * (-t500 - t522) - t136 * t606) * t618) * qJD(1)) * m(5) + (t179 + (t168 + (t412 * t618 + t413 * t616) * t568 + t800 + t914) * t545 + (-t414 * t946 + t960 + t167 + (t412 * t616 - t413 * t618) * t568 + t913) * t544) * t1035 - t1009 / 0.2e1 + t1010 / 0.2e1 - t1011 / 0.2e1 + (-qJD(3) * t716 + t489 * t584 + t490 * t583 + t568 * t794 + t569 * t795) * qJD(1) + t1126 * t1043 + t1084 * t1047 + t1085 * t1048 + (t100 - t866 + ((t618 * t796 + t176 - t910) * t618 + (t616 * t796 + t175 + t799) * t616) * qJD(3)) * t818 + (t42 * (-t766 + t895) + t102 * (t517 - t767 + t839) + t41 * (t504 - t921) + t103 * (-pkin(4) * t854 + t518 + t712 + t849) + (-t41 * t606 + t42 * t838 + (t102 * t1031 * t608 - t42 * pkin(4)) * t569) * t616 + ((-t522 + t838 - t1029) * t975 + (t102 * (-t502 - t522 - t1017) - t103 * t606) * t618) * qJD(1) - (-t102 + t1099 + t647) * t103) * m(6) + (t199 + ((t174 - t374 + (t424 + t957) * t618 + t911) * t618 + t910 * t616) * qJD(3)) * t815 + ((t210 + t226) * t618 + (t225 + t209) * t616) * t864 / 0.2e1 + t969 / 0.2e1 + t970 / 0.2e1 + (t105 + t116 + t94) * t1036 + (-t867 + (t170 - t682 - t913) * t545 + (t616 * t797 + t169 - t365) * t544 + ((t411 + t721) * t544 + t797 * t545) * t618 + t93) * t1038 + t1008 / 0.2e1 + (-(-qJD(1) * t433 + t1079 - t177 + t770) * t178 + t134 * (t616 * t807 + t887 - t936) + t177 * (t562 + t589) + t133 * t1078 + t178 * (t588 + t889) + (t513 * t963 - t962) * qJD(3) + ((-t177 * rSges(4,3) + t178 * t807) * t616 + (t177 * (-t570 - t768) - t178 * t614) * t618) * qJD(1)) * m(4) + (-(-qJD(1) * t474 - t345 - t526 + t588) * t346 + t293 * (t616 * t821 + t592 + t885) + t345 * t589 + t292 * t1077 + t346 * (t884 + t886) + (t345 * (t821 + t1019) * t618 + (t345 * (-rSges(3,3) - qJ(2)) + t346 * t821) * t616) * qJD(1)) * m(3) + t1064 + (t25 * (t895 + t1145) + t78 * (-t763 + t839) + t24 * (t504 + t1147) + t79 * (t588 + t1143) + (-t24 * t949 + t79 * (-t568 * t944 - t569 * t938 - t860) + (t25 * t615 + (-t569 * t615 * t79 + t617 * t78) * qJD(5)) * pkin(5) + (t78 * (t705 + t949) - t79 * t606) * qJD(1)) * t618 + (-t24 * t606 + t78 * (t1005 * t608 + t857) * t569 + (t1005 * t25 + t78 * t789) * t568 + (-t1028 * t78 + t705 * t79) * qJD(1)) * t616 - (t647 - t78 + t1139) * t79 + t809 * t460) * m(7) + t967 / 0.2e1 + t968 / 0.2e1 + (t1157 + t1126) * t1044 + (t104 + t117) * t1037; 0.2e1 * (t1032 * t24 + t1033 * t25) * m(7) + 0.2e1 * (-t1007 / 0.2e1 + t1006 / 0.2e1) * m(6) + 0.2e1 * (t1032 * t98 + t1033 * t99) * m(5) + 0.2e1 * (t1032 * t133 + t1033 * t134) * m(4) + 0.2e1 * (t1032 * t292 + t1033 * t293) * m(3); ((-t876 * t952 - t879) * t618 + (t650 + (t618 * t951 + t627) * qJD(3)) * t616) * t815 + ((t583 * t893 + t584 * t892) * qJD(1) + (t1055 * t584 + t665 * t583) * qJD(3)) * t1027 + ((-t587 * t951 + t879) * t616 + (t650 + (t616 * t952 + t627) * qJD(3)) * t618) * t818 + (-t124 * t618 + t125 * t616 + (t225 * t616 + t226 * t618) * qJD(1)) * t1026 + t622 + (qJD(1) * t118 + (t616 * (t1059 * t616 + t630 * t618) - t618 * (t1058 * t616 + t631 * t618) + (t175 * t616 + t176 * t618) * qJD(1)) * t1098) * t1033 + (qJD(1) * t119 + (t616 * (-t1059 * t618 + t630 * t616) - t618 * (-t1058 * t618 + t631 * t616) + (t173 * t616 + t174 * t618) * qJD(1)) * t1098) * t1032 + (t678 + t100) * t820 + (t677 + t101) * t819 + (t78 * t844 + t19 * (t706 + t918) + t66 * (t662 + t848) + (t1090 * t715 + t690 * t78) * t618 + (t24 * t715 + t79 * t690 + (-t344 + t847) * t978) * t616 - t639 * t78 - (-t79 * t835 + ((-t616 * t79 - t618 * t78) * t584 + t66 * t805) * qJD(3)) * pkin(3) - t621) * m(7) + (t28 * (t785 + t918) + t95 * (t707 + t848) + (t103 * t710 + (-t344 + t921) * t976) * t616 - (-t103 * t835 + (-t584 * t743 + t805 * t95) * qJD(3)) * pkin(3) - t623 + t1054 * (-t398 + t811) + (t618 * t710 + t1068) * t102) * m(6) + (t75 * (t905 + t918) + t126 * (t848 + t851) + (t135 * t772 + (qJD(1) * t136 + t99) * t812) * t618 + (t98 * t812 + t136 * t772 + (t135 * t499 + t126 * (-t344 - t417)) * qJD(1)) * t616 - (-t136 * t835 + (t126 * t805 + t584 * t730) * qJD(3)) * pkin(3) + t1063) * m(5) + (-(t177 * t481 - t962) * qJD(1) - (t253 * (-t481 * t616 - t482 * t618) + t729 * t768) * qJD(3) + 0.2e1 * t253 * (t304 * t618 + t305 * t616 + (t433 * t618 - t434 * t616) * qJD(1)) + t729 * t492 + (-t133 * t616 - t134 * t618 + (-t178 * t618 + t963) * qJD(1)) * t513) * m(4); t622 + (t19 * t706 + t66 * t662 + (t79 * t852 + t847 * t978) * t616 - t621 + (t1090 * t618 + t1092) * (-t501 - t915) + (t618 * t852 - t639 + t844) * t78) * m(7) + (-t623 + t28 * t785 + t95 * t707 + (t103 * t928 + t921 * t976) * t616 + t1054 * (-t398 - t501) + (t618 * t928 + t1068) * t102) * m(6) + (t75 * t905 + t126 * (-t417 * t878 + t851) + t730 * t435 + (-t98 * t616 - t99 * t618 + (-t136 * t618 + t966) * qJD(1)) * t499 + t1063) * m(5); t1135 * t948 / 0.2e1 + t1136 * t947 / 0.2e1 + (t1071 * t568 - t1085 * t569) * t1048 + (t1072 * t568 - t1084 * t569) * t1047 + (-t1074 * t618 + t1075 * t487 - t1076 * t486) * t1046 + ((t1072 * t608 - t1095) * t569 + (qJD(1) * t1102 + t1084 * t608 + t1132 * t618 + t1133 * t616) * t568) * t1045 + ((t1071 * t608 - t1157) * t569 + (qJD(1) * t1101 + t1085 * t608 + t1130 * t618 + t1131 * t616) * t568) * t1044 + (-t1074 * t616 + t1075 * t485 + t1076 * t484) * t1043 + (t1144 * t569 + (t1075 * t617 + t1076 * t615) * t568) * t1040 + ((t1073 * t608 - t1149) * t569 + (qJD(1) * t1103 + t1083 * t608 + t1128 * t618 + t1129 * t616) * t568) * t1039 - (t1010 + t969 + t970 - t1011 + t1008 + t967 + t968 - t1009 + t1064) * t569 / 0.2e1 + ((-t25 * t927 + t78 * t931 - t24 * t926 - t79 * t930 + (-t1100 * t618 + t666 * t616) * t608) * t569 + ((t78 * t927 + t79 * t926) * t608 + (qJD(1) * t666 - t19 * t927 - t24 * t915 + t66 * t931 - t79 * t929) * t618 + (qJD(1) * t1100 - t19 * t926 + t25 * t915 - t66 * t930 + t78 * t929) * t616) * t568 - (-t78 * t917 + t79 * t916) * t541 - (t66 * t916 + t776 * t78) * t460 - (t66 * t917 + t776 * t79) * t459) * m(7) + (-t102 * (-t335 * t541 - t452 * t460) - t103 * (t337 * t541 - t452 * t459) - t95 * (t335 * t459 + t337 * t460) + (t102 * t158 - t103 * t156 + t42 * t287 - t41 * t291 + (t95 * t724 + (t102 * t616 - t103 * t618) * t398) * t608) * t569 + (t102 * (t207 * t616 - t287 * t608) + t103 * (-t207 * t618 + t291 * t608) + t28 * t724 + t95 * (-t156 * t616 + t158 * t618 - t287 * t878 - t291 * t877) + (qJD(1) * t743 + t1006 - t1007) * t398) * t568) * m(6) + ((t1073 * t568 - t1083 * t569) * qJD(5) + t1125) * t950 / 0.2e1 + t1127 * (t568 * t819 + t616 * t822) + t1126 * (-t837 / 0.2e1 + t618 * t822); 0.2e1 * ((t545 * t78 + t79 * t937 - t19) * t1049 + (t459 * t79 + t460 * t78) * t1050) * t569 + 0.2e1 * ((t25 * t618 + t608 * t66 - t78 * t878 + t79 * t877 + t1092) * t1049 + (t66 * (t459 * t616 + t460 * t618) + (-t616 * t78 + t618 * t79) * t541) * t1050) * t568;];
tauc  = t1(:);
