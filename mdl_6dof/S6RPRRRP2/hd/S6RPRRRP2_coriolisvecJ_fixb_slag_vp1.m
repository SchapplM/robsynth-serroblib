% Calculate vector of centrifugal and Coriolis load on the joints for
% S6RPRRRP2
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
% Datum: 2019-03-09 06:01
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauc = S6RPRRRP2_coriolisvecJ_fixb_slag_vp1(qJ, qJD, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(10,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRRRP2_coriolisvecJ_fixb_slag_vp1: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RPRRRP2_coriolisvecJ_fixb_slag_vp1: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RPRRRP2_coriolisvecJ_fixb_slag_vp1: pkin has to be [10x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RPRRRP2_coriolisvecJ_fixb_slag_vp1: m has to be [7x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [7,3]), ...
  'S6RPRRRP2_coriolisvecJ_fixb_slag_vp1: rSges has to be [7x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [7 6]), ...
  'S6RPRRRP2_coriolisvecJ_fixb_slag_vp1: Icges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From coriolisvec_joint_fixb_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 05:58:53
% EndTime: 2019-03-09 06:00:54
% DurationCPUTime: 107.98s
% Computational Cost: add. (69382->1490), mult. (71207->1918), div. (0->0), fcn. (67224->10), ass. (0->711)
t1252 = -Icges(6,4) - Icges(7,4);
t1217 = Icges(6,1) + Icges(7,1);
t1249 = Icges(6,5) + Icges(7,5);
t1222 = Icges(6,2) + Icges(7,2);
t1239 = Icges(6,6) + Icges(7,6);
t645 = qJ(4) + qJ(5);
t638 = cos(t645);
t1251 = t1252 * t638;
t637 = sin(t645);
t1250 = t1252 * t637;
t1248 = Icges(6,3) + Icges(7,3);
t1247 = -t1239 * t637 + t1249 * t638;
t1246 = -t1222 * t637 - t1251;
t1245 = t1217 * t638 + t1250;
t647 = sin(qJ(3));
t650 = cos(qJ(3));
t1204 = t1247 * t647 - t1248 * t650;
t1244 = t1247 * t650 + t1248 * t647;
t1224 = -t1239 * t650 + t1246 * t647;
t1143 = t1239 * t647 + t1246 * t650;
t1202 = t1245 * t647 - t1249 * t650;
t1142 = t1245 * t650 + t1249 * t647;
t1243 = (-t1239 * t638 - t1249 * t637) * t647;
t1242 = (-t1222 * t638 + t1250) * t647;
t1241 = (-t1217 * t637 + t1251) * t647;
t644 = qJ(1) + pkin(10);
t635 = sin(t644);
t1002 = t635 * t647;
t993 = t637 * t650;
t636 = cos(t644);
t998 = t636 * t638;
t470 = t635 * t993 + t998;
t992 = t638 * t650;
t999 = t636 * t637;
t471 = t635 * t992 - t999;
t252 = Icges(7,5) * t471 - Icges(7,6) * t470 + Icges(7,3) * t1002;
t453 = Icges(7,4) * t471;
t258 = -Icges(7,2) * t470 + Icges(7,6) * t1002 + t453;
t452 = Icges(7,4) * t470;
t265 = -Icges(7,1) * t471 - Icges(7,5) * t1002 + t452;
t1004 = t635 * t638;
t472 = -t636 * t993 + t1004;
t1005 = t635 * t637;
t473 = t636 * t992 + t1005;
t996 = t636 * t647;
t105 = t252 * t996 + t472 * t258 - t265 * t473;
t255 = Icges(6,5) * t471 - Icges(6,6) * t470 + Icges(6,3) * t1002;
t456 = Icges(6,4) * t471;
t261 = -Icges(6,2) * t470 + Icges(6,6) * t1002 + t456;
t455 = Icges(6,4) * t470;
t268 = -Icges(6,1) * t471 - Icges(6,5) * t1002 + t455;
t107 = t255 * t996 + t472 * t261 - t268 * t473;
t1214 = t107 + t105;
t254 = Icges(7,5) * t473 + Icges(7,6) * t472 + Icges(7,3) * t996;
t1049 = Icges(7,4) * t473;
t260 = Icges(7,2) * t472 + Icges(7,6) * t996 + t1049;
t454 = Icges(7,4) * t472;
t266 = Icges(7,1) * t473 + Icges(7,5) * t996 + t454;
t106 = t254 * t996 + t472 * t260 + t473 * t266;
t257 = Icges(6,5) * t473 + Icges(6,6) * t472 + Icges(6,3) * t996;
t1052 = Icges(6,4) * t473;
t263 = Icges(6,2) * t472 + Icges(6,6) * t996 + t1052;
t457 = Icges(6,4) * t472;
t269 = Icges(6,1) * t473 + Icges(6,5) * t996 + t457;
t108 = t257 * t996 + t472 * t263 + t473 * t269;
t1213 = t108 + t106;
t643 = qJD(4) + qJD(5);
t831 = t643 * t650;
t1129 = qJD(1) - t831;
t927 = qJD(3) * t647;
t694 = t1129 * t638 + t637 * t927;
t929 = qJD(1) * t650;
t833 = -t643 + t929;
t223 = t1005 * t833 + t636 * t694;
t693 = t1129 * t637 - t638 * t927;
t224 = -t1004 * t833 + t636 * t693;
t926 = qJD(3) * t650;
t882 = t636 * t926;
t930 = qJD(1) * t647;
t893 = t635 * t930;
t712 = t882 - t893;
t131 = Icges(7,5) * t224 + Icges(7,6) * t223 + Icges(7,3) * t712;
t133 = Icges(6,5) * t224 + Icges(6,6) * t223 + Icges(6,3) * t712;
t1238 = t131 + t133;
t225 = t635 * t694 - t833 * t999;
t226 = t635 * t693 + t833 * t998;
t887 = t635 * t926;
t891 = t636 * t930;
t714 = t887 + t891;
t132 = Icges(7,5) * t226 + Icges(7,6) * t225 + Icges(7,3) * t714;
t134 = Icges(6,5) * t226 + Icges(6,6) * t225 + Icges(6,3) * t714;
t1237 = t132 + t134;
t135 = Icges(7,4) * t224 + Icges(7,2) * t223 + Icges(7,6) * t712;
t137 = Icges(6,4) * t224 + Icges(6,2) * t223 + Icges(6,6) * t712;
t1236 = t135 + t137;
t136 = Icges(7,4) * t226 + Icges(7,2) * t225 + Icges(7,6) * t714;
t138 = Icges(6,4) * t226 + Icges(6,2) * t225 + Icges(6,6) * t714;
t1235 = t136 + t138;
t139 = Icges(7,1) * t224 + Icges(7,4) * t223 + Icges(7,5) * t712;
t141 = Icges(6,1) * t224 + Icges(6,4) * t223 + Icges(6,5) * t712;
t1234 = t139 + t141;
t140 = Icges(7,1) * t226 + Icges(7,4) * t225 + Icges(7,5) * t714;
t142 = Icges(6,1) * t226 + Icges(6,4) * t225 + Icges(6,5) * t714;
t1233 = t140 + t142;
t1176 = t1202 * t473 + t1204 * t996 + t1224 * t472;
t1232 = t252 + t255;
t1206 = t254 + t257;
t1231 = t258 + t261;
t1230 = t260 + t263;
t1229 = -t265 - t268;
t1228 = t266 + t269;
t1227 = qJD(3) * t1244 + t1243 * t643;
t1226 = qJD(3) * t1143 + t1242 * t643;
t1225 = qJD(3) * t1142 + t1241 * t643;
t1223 = t1202 * t638 - t1224 * t637;
t101 = t1002 * t252 - t258 * t470 - t265 * t471;
t103 = t1002 * t255 - t261 * t470 - t268 * t471;
t1216 = t101 + t103;
t102 = t254 * t1002 - t470 * t260 + t471 * t266;
t104 = t257 * t1002 - t470 * t263 + t471 * t269;
t1215 = t102 + t104;
t1177 = t1002 * t1204 + t1202 * t471 - t1224 * t470;
t1195 = t261 * t637 + t268 * t638;
t1196 = t258 * t637 + t265 * t638;
t1221 = t1195 + t1196;
t649 = cos(qJ(4));
t988 = t649 * t650;
t646 = sin(qJ(4));
t997 = t636 * t646;
t508 = t635 * t988 - t997;
t494 = Icges(5,4) * t508;
t991 = t646 * t650;
t995 = t636 * t649;
t507 = t635 * t991 + t995;
t305 = -Icges(5,2) * t507 + Icges(5,6) * t1002 + t494;
t493 = Icges(5,4) * t507;
t309 = -Icges(5,1) * t508 - Icges(5,5) * t1002 + t493;
t1194 = t305 * t646 + t309 * t649;
t302 = Icges(5,5) * t508 - Icges(5,6) * t507 + Icges(5,3) * t1002;
t123 = -t1194 * t647 - t302 * t650;
t1191 = t1229 * t224 + t1231 * t223 + t1232 * t712 + t1233 * t473 + t1235 * t472 + t1237 * t996;
t1190 = t1206 * t712 + t1228 * t224 + t1230 * t223 + t1234 * t473 + t1236 * t472 + t1238 * t996;
t1189 = t1002 * t1237 + t1229 * t226 + t1231 * t225 + t1232 * t714 + t1233 * t471 - t1235 * t470;
t1188 = t1002 * t1238 + t1206 * t714 + t1228 * t226 + t1230 * t225 + t1234 * t471 - t1236 * t470;
t622 = qJD(3) * t635;
t925 = qJD(4) * t647;
t542 = t636 * t925 + t622;
t923 = qJD(5) * t647;
t435 = t636 * t923 + t542;
t928 = qJD(3) * t636;
t543 = -t635 * t925 + t928;
t436 = -t635 * t923 + t543;
t1151 = t1129 * t1176 + t1213 * t435 - t1214 * t436;
t1220 = t1202 * t224 + t1204 * t712 + t1224 * t223 + t1225 * t473 + t1226 * t472 + t1227 * t996;
t1219 = t1002 * t1227 + t1202 * t226 + t1204 * t714 + t1224 * t225 + t1225 * t471 - t1226 * t470;
t1218 = (qJD(3) * t1223 - t1227) * t650 + ((-t1224 * t643 + t1225) * t638 + (-t1202 * t643 - t1226) * t637 + t1204 * qJD(3)) * t647;
t1184 = t252 * t650;
t117 = -t1196 * t647 - t1184;
t1182 = t255 * t650;
t119 = -t1195 * t647 - t1182;
t1212 = t119 + t117;
t762 = -t260 * t637 + t266 * t638;
t118 = -t254 * t650 + t647 * t762;
t760 = -t263 * t637 + t269 * t638;
t120 = -t257 * t650 + t647 * t760;
t1211 = t120 + t118;
t1175 = -t1204 * t650 + t1223 * t647;
t1073 = rSges(7,3) * t647;
t640 = t649 * pkin(4);
t634 = t640 + pkin(3);
t575 = pkin(5) * t638 + t640;
t912 = pkin(3) + t575;
t836 = t634 - t912;
t749 = t836 * t650;
t796 = rSges(7,1) * t638 - rSges(7,2) * t637;
t652 = -pkin(9) - pkin(8);
t642 = -qJ(6) + t652;
t935 = t652 - t642;
t844 = t935 * t647;
t1210 = t650 * t796 + t1073 - t749 + t844;
t1209 = (-t1239 * t471 - t1249 * t470) * t436 + (t1239 * t473 - t1249 * t472) * t435 - t1243 * t1129;
t1208 = t1129 * (-t1223 + t1244) - (-t1204 * t635 + t1221) * t436 + (-t1204 * t636 - t760 - t762) * t435;
t1152 = t1129 * t1177 + t1215 * t435 - t1216 * t436;
t1084 = pkin(5) * t637;
t1085 = pkin(4) * t646;
t574 = t1084 + t1085;
t819 = t650 * t912;
t1205 = t473 * rSges(7,1) + t472 * rSges(7,2) + rSges(7,3) * t996 + t635 * t574 + t636 * t819;
t1201 = -rSges(7,1) * t471 + rSges(7,2) * t470 - t635 * t819;
t680 = t647 * t836 + t650 * t935;
t1171 = t650 * rSges(7,3) - t647 * t796 + t680;
t1173 = (-t1222 * t471 + t1229 - t452 - t455) * t436 + (t1222 * t473 - t1228 - t454 - t457) * t435 + (-t1202 - t1242) * t1129;
t113 = t1002 * t302 - t305 * t507 - t309 * t508;
t1003 = t635 * t646;
t510 = t636 * t988 + t1003;
t1001 = t635 * t649;
t736 = t636 * t991 - t1001;
t304 = Icges(5,5) * t510 - Icges(5,6) * t736 + Icges(5,3) * t996;
t1055 = Icges(5,4) * t510;
t307 = -Icges(5,2) * t736 + Icges(5,6) * t996 + t1055;
t495 = Icges(5,4) * t736;
t310 = Icges(5,1) * t510 + Icges(5,5) * t996 - t495;
t114 = t304 * t1002 - t507 * t307 + t508 * t310;
t785 = Icges(5,5) * t649 - Icges(5,6) * t646;
t511 = -Icges(5,3) * t650 + t647 * t785;
t1053 = Icges(5,4) * t649;
t788 = -Icges(5,2) * t646 + t1053;
t513 = -Icges(5,6) * t650 + t647 * t788;
t1054 = Icges(5,4) * t646;
t792 = Icges(5,1) * t649 - t1054;
t515 = -Icges(5,5) * t650 + t647 * t792;
t176 = t1002 * t511 - t507 * t513 + t508 * t515;
t924 = qJD(4) * t650;
t623 = qJD(1) - t924;
t67 = -t113 * t543 + t114 * t542 + t176 * t623;
t115 = t302 * t996 - t305 * t736 - t309 * t510;
t116 = t304 * t996 - t307 * t736 + t510 * t310;
t177 = t510 * t515 + t511 * t996 - t513 * t736;
t68 = -t115 * t543 + t116 * t542 + t177 * t623;
t922 = qJD(6) * t647;
t596 = t636 * t922;
t986 = t650 * t634;
t557 = t635 * t986;
t853 = -t574 + t1085;
t980 = -t636 * t853 + t1201 + t557 + (-rSges(7,3) - t935) * t1002;
t1197 = t1129 * t980 + t596;
t653 = qJD(1) ^ 2;
t919 = qJD(1) * qJD(3);
t617 = t635 * t919;
t918 = qJD(3) * qJD(4);
t866 = t650 * t918;
t868 = qJD(1) * t925;
t428 = t635 * t866 + t636 * t868 + t617;
t311 = qJD(5) * t714 + t428;
t618 = t636 * t919;
t940 = t636 * t866 + t618;
t312 = qJD(5) * t882 - t643 * t893 + t940;
t588 = t643 * t647;
t563 = qJD(3) * t588;
t1193 = t1220 * t1129 + t1176 * t563 + t1190 * t435 - t1191 * t436 + t1213 * t312 + t1214 * t311;
t1192 = t1129 * t1219 + t1177 * t563 + t1188 * t435 - t1189 * t436 + t1215 * t312 + t1216 * t311;
t35 = (-qJD(3) * t1196 - t132) * t650 + (qJD(3) * t252 + (-t258 * t643 + t140) * t638 + (t265 * t643 - t136) * t637) * t647;
t37 = (-qJD(3) * t1195 - t134) * t650 + (qJD(3) * t255 + (-t261 * t643 + t142) * t638 + (t268 * t643 - t138) * t637) * t647;
t1187 = t35 + t37;
t36 = (qJD(3) * t762 - t131) * t650 + (qJD(3) * t254 + (-t260 * t643 + t139) * t638 + (-t266 * t643 - t135) * t637) * t647;
t38 = (qJD(3) * t760 - t133) * t650 + (qJD(3) * t257 + (-t263 * t643 + t141) * t638 + (-t269 * t643 - t137) * t637) * t647;
t1186 = t36 + t38;
t1185 = t1129 * t1175 + t1211 * t435 - t1212 * t436;
t1178 = t635 * t636;
t1083 = pkin(5) * t643;
t1071 = pkin(4) * qJD(4);
t910 = t649 * t1071;
t549 = t1083 * t638 + t910;
t820 = t647 * t912;
t702 = (-t650 * t642 - t820) * qJD(3);
t911 = t646 * t1071;
t548 = -t1083 * t637 - t911;
t987 = t650 * t548;
t748 = -t922 - t987;
t1163 = t652 * t930 + t910;
t830 = t650 * t911;
t1164 = t652 * t926 + t830;
t888 = t635 * t927;
t794 = t1163 * t636 + t1164 * t635 + t634 * t888;
t798 = t226 * rSges(7,1) + t225 * rSges(7,2);
t990 = t647 * t642;
t983 = -rSges(7,3) * t714 - t798 - (-t549 + (-t749 - t990) * qJD(1)) * t636 - (-qJD(1) * t853 + t702 - t748) * t635 - t794;
t614 = pkin(4) * t1003;
t941 = t636 * t986 + t614;
t979 = t636 * t844 + t1205 - t941;
t537 = (-rSges(7,1) * t637 - rSges(7,2) * t638) * t647;
t921 = qJD(6) * t650;
t1174 = -t921 + (t548 + t911) * t647 + t537 * t643 + t1210 * qJD(3);
t1172 = (t1217 * t470 + t1231 + t453 + t456) * t436 + (t1217 * t472 - t1049 - t1052 - t1230) * t435 + (-t1224 + t1241) * t1129;
t1170 = t1209 * t647;
t1169 = t1211 * t636 + t1212 * t635;
t1168 = t1213 * t636 + t1214 * t635;
t1167 = t1215 * t636 + t1216 * t635;
t858 = t926 / 0.2e1;
t1166 = t636 * t858 - t893 / 0.2e1;
t931 = qJD(1) * t636;
t864 = t931 / 0.2e1;
t1165 = t635 * t858 + t647 * t864;
t586 = pkin(8) * t882;
t883 = t636 * t927;
t892 = t635 * t929;
t713 = -t883 - t892;
t375 = pkin(3) * t713 - pkin(8) * t893 + t586;
t584 = pkin(3) * t888;
t994 = t636 * t650;
t616 = pkin(3) * t994;
t376 = pkin(8) * t714 + qJD(1) * t616 - t584;
t612 = pkin(3) * t647 - pkin(8) * t650;
t527 = t612 * t635;
t1082 = pkin(8) * t647;
t1086 = pkin(3) * t650;
t613 = t1082 + t1086;
t528 = t613 * t635;
t529 = t612 * t636;
t1162 = t636 * t375 + t635 * t376 + t527 * t622 + t528 * t931 + t529 * t928;
t560 = t636 * pkin(2) + t635 * pkin(7);
t651 = cos(qJ(1));
t641 = t651 * pkin(1);
t1141 = t641 + t560;
t627 = t635 * rSges(4,3);
t444 = rSges(4,1) * t994 - rSges(4,2) * t996 + t627;
t1158 = t444 + t1141;
t1080 = pkin(8) + t652;
t1081 = pkin(3) - t634;
t869 = t1081 * t647;
t710 = -t1080 * t650 + t869;
t870 = t1081 * t650;
t1157 = t870 + t1082;
t1156 = t1129 * t1218 + t1175 * t563;
t803 = rSges(5,1) * t508 - rSges(5,2) * t507;
t313 = rSges(5,3) * t1002 + t803;
t802 = rSges(5,1) * t649 - rSges(5,2) * t646;
t525 = -rSges(5,3) * t650 + t647 * t802;
t1155 = -t313 * t623 - t525 * t543;
t800 = rSges(6,1) * t471 - rSges(6,2) * t470;
t272 = rSges(6,3) * t1002 + t800;
t799 = rSges(6,1) * t638 - rSges(6,2) * t637;
t487 = -rSges(6,3) * t650 + t647 * t799;
t1154 = -t1129 * t272 - t436 * t487;
t1153 = 0.2e1 * qJD(3);
t496 = t507 * pkin(4);
t497 = t736 * pkin(4);
t530 = pkin(8) * t996 + t616;
t917 = pkin(4) * t997;
t896 = qJD(1) * t917 + t1163 * t635;
t932 = qJD(1) * t635;
t167 = -t586 + t1157 * t932 + (-t830 + (-t650 * t652 + t869) * qJD(3)) * t636 + t896;
t168 = -pkin(8) * t887 + t584 + (-t1157 * t636 + t614) * qJD(1) - t794;
t1132 = t1080 * t647;
t322 = t917 - t557 + (t1086 + t1132) * t635;
t989 = t647 * t652;
t323 = -t636 * t989 - t530 + t941;
t429 = -t635 * t868 + t940;
t900 = t375 * t928 + t376 * t622 + t528 * t618;
t679 = t167 * t543 + t542 * t168 - t429 * t322 - t323 * t428 + t900;
t1137 = t224 * rSges(7,1) + t223 * rSges(7,2) + rSges(7,3) * t882 + t635 * t549 + t574 * t931 + t636 * t987 + t642 * t893 + t596;
t985 = t836 * t892 + (qJD(3) * t680 + t830) * t636 - t896 - rSges(7,3) * t893 + t1137;
t19 = t985 * t436 - t983 * t435 - t980 * t312 - t979 * t311 + (-t530 * t932 + t922) * qJD(3) + t679;
t1150 = t19 * t979;
t1149 = t623 * t322 + t543 * t710;
t1148 = t1002 * t1174 - t1171 * t714 + t436 * t537 - t650 * t983;
t1147 = t1224 * t635;
t1146 = t1224 * t636;
t1145 = t1202 * t635;
t1144 = t1202 * t636;
t845 = t636 * rSges(3,1) - rSges(3,2) * t635;
t1140 = t641 + t845;
t1139 = (t1129 * t1204 + t1206 * t435 - t1232 * t436) * t650 + t1208 * t647;
t373 = t710 * t635;
t374 = t710 * t636;
t879 = t636 * t924;
t880 = t635 * t924;
t1138 = t636 * t167 + t635 * t168 + t323 * t880 - t542 * t373 - t374 * t543 + t1162 + (t879 - t931) * t322;
t572 = qJD(3) * t613;
t469 = -t870 - t1132;
t397 = qJD(3) * t469 - t647 * t911;
t905 = -t397 - t1174;
t1136 = -t572 + t905 - t921;
t916 = t647 * t1085;
t1135 = t397 * t1002 + t650 * t168 - t496 * t623 - t543 * t916 - t710 * t714;
t801 = rSges(6,1) * t226 + rSges(6,2) * t225;
t146 = rSges(6,3) * t714 + t801;
t1075 = rSges(6,3) * t647;
t489 = t650 * t799 + t1075;
t538 = (-rSges(6,1) * t637 - rSges(6,2) * t638) * t647;
t333 = qJD(3) * t489 + t538 * t643;
t349 = -rSges(6,1) * t470 - rSges(6,2) * t471;
t1134 = t333 * t1002 + t349 * t1129 + t650 * t146 + t436 * t538 + t487 * t714;
t547 = t612 * t932;
t725 = qJD(1) * t527 - t613 * t928;
t878 = t322 * t925;
t1133 = t373 * t623 + t543 * t469 + t547 - t725 - t878 - (-t880 + t932) * t710;
t639 = Icges(4,4) * t650;
t789 = -Icges(4,2) * t647 + t639;
t594 = Icges(4,1) * t647 + t639;
t903 = t224 * rSges(6,1) + t223 * rSges(6,2) + rSges(6,3) * t882;
t144 = -rSges(6,3) * t893 + t903;
t1122 = qJD(1) * t272 + t144;
t1037 = Icges(4,3) * t636;
t591 = Icges(4,5) * t650 - Icges(4,6) * t647;
t590 = Icges(4,5) * t647 + Icges(4,6) * t650;
t719 = qJD(3) * t590;
t440 = Icges(4,6) * t635 + t636 * t789;
t1010 = t440 * t647;
t1056 = Icges(4,4) * t647;
t595 = Icges(4,1) * t650 - t1056;
t442 = Icges(4,5) * t635 + t595 * t636;
t755 = -t442 * t650 + t1010;
t1121 = -t636 * t719 + (-t591 * t635 + t1037 + t755) * qJD(1);
t438 = Icges(4,3) * t635 + t591 * t636;
t1000 = t635 * t650;
t1041 = Icges(4,6) * t636;
t439 = Icges(4,4) * t1000 - Icges(4,2) * t1002 - t1041;
t1011 = t439 * t647;
t1046 = Icges(4,5) * t636;
t606 = Icges(4,4) * t1002;
t441 = Icges(4,1) * t1000 - t1046 - t606;
t756 = -t441 * t650 + t1011;
t1120 = -t635 * t719 + (t438 + t756) * qJD(1);
t592 = Icges(4,2) * t650 + t1056;
t750 = t592 * t647 - t594 * t650;
t1119 = qJD(1) * t750 + t591 * qJD(3);
t437 = Icges(4,5) * t1000 - Icges(4,6) * t1002 - t1037;
t180 = -t437 * t636 - t635 * t756;
t952 = -Icges(4,2) * t1000 + t441 - t606;
t954 = t594 * t635 + t439;
t1118 = -t647 * t952 - t650 * t954;
t512 = Icges(5,3) * t647 + t650 * t785;
t751 = -t513 * t646 + t515 * t649;
t758 = -t307 * t646 + t310 * t649;
t1115 = t542 * (-t511 * t636 - t758) - t543 * (-t511 * t635 + t1194) + t623 * (t512 - t751);
t553 = (-Icges(5,2) * t649 - t1054) * t647;
t1114 = t542 * (-Icges(5,2) * t510 + t310 - t495) - t543 * (-Icges(5,2) * t508 - t309 - t493) + t623 * (t515 + t553);
t1113 = -m(7) / 0.2e1;
t1112 = m(7) / 0.2e1;
t1111 = t311 / 0.2e1;
t1110 = t312 / 0.2e1;
t1109 = t428 / 0.2e1;
t1108 = t429 / 0.2e1;
t1107 = -t435 / 0.2e1;
t1106 = t435 / 0.2e1;
t1105 = -t436 / 0.2e1;
t1104 = t436 / 0.2e1;
t1101 = -t542 / 0.2e1;
t1100 = t542 / 0.2e1;
t1099 = -t543 / 0.2e1;
t1098 = t543 / 0.2e1;
t1096 = -t1129 / 0.2e1;
t1095 = t1129 / 0.2e1;
t1093 = -t623 / 0.2e1;
t1092 = t623 / 0.2e1;
t1089 = -t650 / 0.2e1;
t1088 = -rSges(5,3) - pkin(8);
t648 = sin(qJ(1));
t1087 = pkin(1) * t648;
t1079 = rSges(4,1) * t650;
t1077 = rSges(5,3) * t647;
t1072 = pkin(1) * qJD(1);
t276 = t473 * rSges(6,1) + t472 * rSges(6,2) + rSges(6,3) * t996;
t818 = t530 * t617;
t22 = t144 * t436 + t146 * t435 + t272 * t312 - t276 * t311 + t679 - t818;
t1070 = t22 * t276;
t1069 = t35 * t436;
t1068 = t36 * t435;
t1067 = t37 * t436;
t1066 = t38 * t435;
t692 = t623 * t649 + t646 * t927;
t829 = -qJD(4) + t929;
t296 = t635 * t692 - t829 * t997;
t691 = t623 * t646 - t649 * t927;
t297 = t635 * t691 + t829 * t995;
t156 = Icges(5,5) * t297 + Icges(5,6) * t296 + Icges(5,3) * t714;
t158 = Icges(5,4) * t297 + Icges(5,2) * t296 + Icges(5,6) * t714;
t160 = Icges(5,1) * t297 + Icges(5,4) * t296 + Icges(5,5) * t714;
t47 = (-qJD(3) * t1194 - t156) * t650 + (qJD(3) * t302 - t158 * t646 + t160 * t649 + (-t305 * t649 + t309 * t646) * qJD(4)) * t647;
t1065 = t47 * t543;
t294 = t1003 * t829 + t636 * t692;
t295 = -t1001 * t829 + t636 * t691;
t155 = Icges(5,5) * t295 + Icges(5,6) * t294 + Icges(5,3) * t712;
t157 = Icges(5,4) * t295 + Icges(5,2) * t294 + Icges(5,6) * t712;
t159 = Icges(5,1) * t295 + Icges(5,4) * t294 + Icges(5,5) * t712;
t48 = (qJD(3) * t758 - t155) * t650 + (qJD(3) * t304 - t157 * t646 + t159 * t649 + (-t307 * t649 - t310 * t646) * qJD(4)) * t647;
t1064 = t48 * t542;
t711 = (t530 + t1141) * qJD(1) - t612 * t622;
t677 = t323 * t623 + t542 * t710 + t711;
t877 = t635 * t922;
t81 = t1129 * t979 + t1171 * t435 + t677 + t877;
t1063 = t636 * t81;
t552 = (-Icges(5,5) * t646 - Icges(5,6) * t649) * t647;
t380 = qJD(3) * t512 + qJD(4) * t552;
t514 = Icges(5,6) * t647 + t650 * t788;
t381 = qJD(3) * t514 + qJD(4) * t553;
t516 = Icges(5,5) * t647 + t650 * t792;
t554 = (-Icges(5,1) * t646 - t1053) * t647;
t382 = qJD(3) * t516 + qJD(4) * t554;
t100 = (qJD(3) * t751 - t380) * t650 + (qJD(3) * t511 - t381 * t646 + t382 * t649 + (-t513 * t649 - t515 * t646) * qJD(4)) * t647;
t210 = -t511 * t650 + t647 * t751;
t867 = t647 * t918;
t1059 = t100 * t623 + t210 * t867;
t1058 = t642 - rSges(7,3);
t631 = t636 * pkin(7);
t559 = pkin(2) * t635 - t631;
t856 = -t559 - t1087;
t884 = t612 * t928;
t678 = (-t528 + t856) * qJD(1) - t884;
t666 = t678 + t1149;
t80 = t1171 * t436 + t1197 + t666;
t1033 = qJD(3) * t80;
t1026 = t117 * t311;
t1025 = t118 * t312;
t1024 = t119 * t311;
t1023 = t120 * t312;
t1022 = t123 * t428;
t124 = -t304 * t650 + t647 * t758;
t1021 = t124 * t429;
t147 = t1155 + t678;
t1020 = t147 * t636;
t936 = rSges(4,2) * t1002 + t636 * rSges(4,3);
t443 = rSges(4,1) * t1000 - t936;
t600 = rSges(4,1) * t647 + rSges(4,2) * t650;
t885 = t600 * t928;
t238 = -t885 + (-t443 + t856) * qJD(1);
t1016 = t238 * t635;
t1015 = t238 * t636;
t889 = t600 * t622;
t239 = qJD(1) * t1158 - t889;
t524 = t600 * t636;
t1014 = t239 * t524;
t1008 = t572 * t636;
t1007 = t590 * t635;
t1006 = t590 * t636;
t984 = t146 * t996 + t272 * t882;
t982 = t168 * t996 - t322 * t882;
t981 = t980 * t996;
t977 = t276 * t927 + t487 * t893;
t972 = -t276 - t323;
t971 = t1171 * t635;
t970 = t1171 * t636;
t969 = t323 * t927 - t710 * t893;
t315 = t510 * rSges(5,1) - rSges(5,2) * t736 + rSges(5,3) * t996;
t966 = -t315 - t530;
t965 = -t333 - t397;
t545 = t560 * qJD(1);
t964 = -t376 - t545;
t526 = t650 * t802 + t1077;
t555 = (-rSges(5,1) * t646 - rSges(5,2) * t649) * t647;
t383 = qJD(3) * t526 + qJD(4) * t555;
t963 = -t383 - t572;
t960 = -t635 * t437 - t441 * t994;
t959 = t635 * t438 + t442 * t994;
t958 = -t1002 * t710 - t650 * t322;
t957 = t487 * t1002 + t650 * t272;
t955 = -t623 * t497 + t542 * t916;
t953 = -t594 * t636 - t440;
t951 = -t592 * t636 + t442;
t948 = t635 * t528 + t636 * t530;
t947 = t710 - t487;
t942 = -t525 - t612;
t939 = rSges(4,2) * t893 + rSges(4,3) * t931;
t938 = -t592 + t595;
t937 = t594 + t789;
t933 = qJD(1) * t591;
t330 = -t635 * t750 - t1006;
t920 = t330 * qJD(1);
t915 = t647 * t1084;
t914 = t653 * t1087;
t913 = t653 * t641;
t909 = t648 * t1072;
t908 = -t167 - t985;
t907 = t322 + t980;
t906 = -t323 - t979;
t902 = t295 * rSges(5,1) + t294 * rSges(5,2) + rSges(5,3) * t882;
t901 = -t572 + t965;
t898 = t710 + t1171;
t897 = -t612 + t947;
t894 = t612 * t931;
t875 = t1002 / 0.2e1;
t874 = t996 / 0.2e1;
t873 = t528 * t622 + t530 * t928 + qJD(2);
t872 = -pkin(2) - t1079;
t871 = -pkin(2) - t1075;
t863 = -t622 / 0.2e1;
t860 = t928 / 0.2e1;
t859 = t927 / 0.2e1;
t854 = t631 - t1087;
t722 = -t542 * t322 + t323 * t543 + t873;
t66 = -t435 * t980 + t436 * t979 + t722 - t921;
t852 = t66 * t980;
t851 = t66 * t979;
t850 = t80 * t980;
t849 = t80 * t1171;
t848 = t81 * t1171;
t86 = t272 * t435 + t276 * t436 + t722;
t847 = t86 * t972;
t351 = rSges(6,1) * t472 - rSges(6,2) * t473;
t843 = t435 * t349 + t351 * t436;
t842 = t1129 * t351 - t435 * t538;
t840 = -t542 * t496 - t497 * t543;
t404 = t442 * t1000;
t839 = t438 * t636 - t404;
t838 = -t437 + t1010;
t837 = t964 * qJD(1);
t832 = -t572 + t921;
t828 = -t882 * t980 - t983 * t996;
t825 = -t1171 * t893 + t927 * t979;
t823 = -t635 * t322 + t636 * t323 + t948;
t822 = -t1002 * t1171 - t650 * t980;
t821 = -t612 + t898;
t621 = pkin(7) * t931;
t813 = qJD(1) * (-pkin(2) * t932 + t621) - t914;
t812 = t612 * t617 - t913;
t811 = t66 * t906;
t808 = -qJD(1) * t529 - t613 * t622;
t556 = rSges(3,1) * t635 + rSges(3,2) * t636;
t805 = -rSges(4,2) * t647 + t1079;
t804 = rSges(5,1) * t297 + rSges(5,2) * t296;
t782 = t101 * t636 - t102 * t635;
t780 = t103 * t636 - t104 * t635;
t778 = t105 * t636 - t106 * t635;
t776 = t107 * t636 - t108 * t635;
t774 = t113 * t636 - t114 * t635;
t773 = t113 * t635 + t114 * t636;
t772 = t115 * t636 - t116 * t635;
t771 = t115 * t635 + t116 * t636;
t770 = t117 * t636 - t118 * t635;
t768 = t119 * t636 - t120 * t635;
t766 = t123 * t636 - t124 * t635;
t765 = t123 * t635 + t124 * t636;
t764 = -t239 * t635 - t1015;
t757 = t313 * t636 - t315 * t635;
t277 = t439 * t650 + t441 * t647;
t278 = t440 * t650 + t442 * t647;
t754 = t443 * t635 + t444 * t636;
t747 = qJD(1) * t375 + t813;
t746 = -pkin(2) - t819;
t744 = t871 - t986;
t726 = -t572 * t635 - t894;
t523 = t600 * t635;
t723 = t1088 * t647 - pkin(2) - t1086;
t721 = qJD(3) * t594;
t720 = qJD(3) * t592;
t181 = -t1002 * t440 - t839;
t718 = (-t180 * t636 + t181 * t635) * qJD(3);
t182 = -t439 * t996 - t960;
t183 = -t440 * t996 + t959;
t717 = (-t182 * t636 + t183 * t635) * qJD(3);
t709 = -t302 * t543 + t304 * t542 + t511 * t623;
t706 = (-Icges(5,5) * t507 - Icges(5,6) * t508) * t543 - (-Icges(5,5) * t736 - Icges(5,6) * t510) * t542 - t552 * t623;
t704 = -t647 * t951 + t650 * t953;
t703 = -t623 * t168 - t543 * t397 - t428 * t710 + t812;
t550 = qJD(1) * t559;
t701 = -qJD(1) * t528 - t550 - t884 - t909;
t698 = t647 * t706;
t681 = (-t647 * t937 + t650 * t938) * qJD(1);
t671 = (-Icges(5,1) * t736 - t1055 - t307) * t542 - (-Icges(5,1) * t507 - t305 - t494) * t543 + (-t513 + t554) * t623;
t669 = t623 * t167 + t323 * t867 - t542 * t397 + t429 * t710 + t747;
t346 = rSges(4,1) * t713 - rSges(4,2) * t882 + t939;
t347 = -qJD(3) * t523 + (t636 * t805 + t627) * qJD(1);
t668 = t346 * t636 + t347 * t635 + (t443 * t636 - t444 * t635) * qJD(1);
t667 = t701 + t1149;
t663 = t323 * t925 + t623 * t374 - t469 * t542 + t710 * t879 + t808;
t121 = t313 * t542 + t315 * t543 + t873;
t148 = t315 * t623 - t525 * t542 + t711;
t662 = t121 * t757 + (t147 * t635 - t148 * t636) * t525;
t565 = t789 * qJD(3);
t566 = t595 * qJD(3);
t658 = qJD(1) * t590 - t565 * t647 + t566 * t650 + (-t592 * t650 - t594 * t647) * qJD(3);
t657 = t1115 * t647;
t656 = t1192 * t875 + t1193 * t874 + t1185 * t859 + (t1167 * t647 - t1177 * t650) * t1111 + (t1168 * t647 - t1176 * t650) * t1110 + (-t1170 * t636 + t1172 * t473 - t1173 * t472) * t1107 + ((qJD(3) * t1168 - t1220) * t650 + (t1190 * t636 + t1191 * t635 + t1176 * qJD(3) + (t776 + t778) * qJD(1)) * t647) * t1106 + ((qJD(3) * t1167 - t1219) * t650 + (t1188 * t636 + t1189 * t635 + t1177 * qJD(3) + (t780 + t782) * qJD(1)) * t647) * t1105 + (-t1170 * t635 + t1172 * t471 + t1173 * t470) * t1104 + (t1169 * t647 - t1175 * t650) * t563 / 0.2e1 + (t1209 * t650 + (t1172 * t638 + t1173 * t637) * t647) * t1096 + ((qJD(3) * t1169 - t1218) * t650 + (t1186 * t636 + t1187 * t635 + t1175 * qJD(3) + (t768 + t770) * qJD(1)) * t647) * t1095 + (t1025 + t1026 + t1068 - t1069 + t1023 + t1024 + t1066 - t1067 + t1156) * t1089 + t1151 * t1166 + t1152 * t1165;
t567 = t805 * qJD(3);
t506 = t636 * t831;
t505 = t635 * t831;
t492 = t853 * t647;
t459 = t472 * pkin(5);
t458 = t470 * pkin(5);
t414 = t525 * t636;
t413 = t525 * t635;
t412 = t515 * t636;
t411 = t515 * t635;
t410 = t513 * t636;
t409 = t513 * t635;
t403 = t487 * t636;
t401 = t487 * t635;
t372 = -rSges(5,1) * t736 - rSges(5,2) * t510;
t371 = -rSges(5,1) * t507 - rSges(5,2) * t508;
t350 = rSges(7,1) * t472 - rSges(7,2) * t473;
t348 = -rSges(7,1) * t470 - rSges(7,2) * t471;
t331 = -t636 * t750 + t1007;
t301 = t331 * qJD(1);
t289 = -t574 * t994 + t575 * t635 + t497;
t288 = -t1000 * t574 - t575 * t636 + t496;
t270 = t322 * t996;
t244 = t1129 * t350;
t237 = qJD(3) * t754 + qJD(2);
t236 = t272 * t996;
t207 = t435 * t348;
t175 = -t913 - t567 * t928 + (-t347 - t545 + t889) * qJD(1);
t174 = -t567 * t622 + (t346 - t885) * qJD(1) + t813;
t163 = rSges(5,3) * t714 + t804;
t162 = -rSges(5,3) * t893 + t902;
t150 = -t1119 * t636 + t658 * t635;
t149 = t1119 * t635 + t658 * t636;
t128 = -qJD(3) * t755 + (-t636 * t720 + (-t635 * t789 + t1041) * qJD(1)) * t650 + (-t636 * t721 + (-t595 * t635 + t1046) * qJD(1)) * t647;
t127 = -qJD(3) * t756 + (qJD(1) * t440 - t635 * t720) * t650 + (qJD(1) * t442 - t635 * t721) * t647;
t122 = t668 * qJD(3);
t99 = t301 + t717;
t98 = t718 + t920;
t96 = t1129 * t276 - t435 * t487 + t677;
t95 = t1154 + t666;
t85 = t1002 * t380 + t296 * t513 + t297 * t515 - t381 * t507 + t382 * t508 + t511 * t714;
t84 = t294 * t513 + t295 * t515 + t380 * t996 - t381 * t736 + t382 * t510 + t511 * t712;
t83 = -t163 * t623 - t383 * t543 + t428 * t525 + (-t313 * t925 - t1008) * qJD(3) + t837 + t812;
t82 = t162 * t623 - t383 * t542 - t429 * t525 + (t315 * t925 + t726) * qJD(3) + t747;
t75 = -t123 * t543 + t124 * t542 + t210 * t623;
t65 = t162 * t543 + t163 * t542 + t313 * t429 - t315 * t428 - t818 + t900;
t46 = t1002 * t155 - t157 * t507 + t159 * t508 + t296 * t307 + t297 * t310 + t304 * t714;
t45 = t1002 * t156 - t158 * t507 + t160 * t508 + t296 * t305 - t297 * t309 + t302 * t714;
t44 = t155 * t996 - t157 * t736 + t159 * t510 + t294 * t307 + t295 * t310 + t304 * t712;
t43 = t156 * t996 - t158 * t736 + t160 * t510 + t294 * t305 - t295 * t309 + t302 * t712;
t42 = -t146 * t1129 - t272 * t563 + t311 * t487 - t333 * t436 + (t878 - t1008) * qJD(3) + t837 + t703;
t41 = qJD(3) * t726 + t1129 * t144 + t276 * t563 - t312 * t487 - t333 * t435 + t669;
t21 = t983 * t1129 + t980 * t563 - t1174 * t436 - t1171 * t311 + (t636 * t832 + t878) * qJD(3) + (-t877 + t964) * qJD(1) + t703;
t20 = qJD(1) * t596 + t985 * t1129 + t979 * t563 - t1174 * t435 + t1171 * t312 + (t635 * t832 - t894) * qJD(3) + t669;
t18 = t113 * t428 + t114 * t429 + t176 * t867 - t45 * t543 + t46 * t542 + t623 * t85;
t17 = t115 * t428 + t116 * t429 + t177 * t867 - t43 * t543 + t44 * t542 + t623 * t84;
t1 = [t1066 / 0.2e1 + (t98 - t920 + ((t636 * t838 + t183 - t959) * t636 + (t635 * t838 + t182 + t839) * t635) * qJD(3)) * t863 - t1067 / 0.2e1 + t1068 / 0.2e1 - t1065 / 0.2e1 + t1059 + t1064 / 0.2e1 + t1025 / 0.2e1 + t1026 / 0.2e1 + t1023 / 0.2e1 + t1024 / 0.2e1 + t1022 / 0.2e1 + (-qJD(3) * t750 + t565 * t650 + t566 * t647) * qJD(1) - (t127 + t150 + t99) * t928 / 0.2e1 + (t128 + t149) * t622 / 0.2e1 + (t83 * (-t803 + t854) + t147 * (t584 - t804) + t82 * (t1141 - t966) + t148 * (-pkin(3) * t883 + t586 + t621 + t902) + (t1088 * t147 * t926 + t723 * t83) * t635 + ((-t147 * t651 - t148 * t648) * pkin(1) + t723 * t1020 + (-t147 * pkin(7) + t148 * (-pkin(2) - t613 - t1077)) * t635) * qJD(1) - (t1155 - t147 + t701) * t148) * m(5) + t68 * t1098 + (t175 * (t635 * t872 + t854 + t936) + t174 * t1158 + t239 * (t621 + t939) + (t1016 * t600 - t1014) * qJD(3) + ((-t238 * t651 - t239 * t648) * pkin(1) + (-pkin(2) - t805) * t1015 + (t238 * (-rSges(4,3) - pkin(7)) + t239 * t872) * t635) * qJD(1) - (-t885 - t238 - t550 + (-t443 - t1087) * qJD(1)) * t239) * m(4) + (t301 + ((t181 - t404 + (t438 + t1011) * t636 + t960) * t636 + t959 * t635) * qJD(3)) * t860 + t177 * t1108 + t176 * t1109 + t84 * t1100 + t1220 * t1106 + (t1219 + t1151) * t1105 + t1151 * t1104 + t1021 / 0.2e1 + t1156 + m(3) * ((-t556 * t653 - t914) * t1140 + (-t913 + (-0.2e1 * t845 - t641 + t1140) * t653) * (-t556 - t1087)) + t1176 * t1110 + t1177 * t1111 - t1069 / 0.2e1 + (t1221 * t647 + t1182 + t1184 + t1212) * t435 * t1096 + (-(t1154 + t667 - t95) * t96 + t42 * (-t557 - t800 + t854) + t95 * (t794 - t801) + t41 * (t1141 + t276 + t941) + t96 * (t621 + t896 + t903) + (t42 * t1085 - t41 * t989 + t96 * (-t634 * t927 - t1164)) * t636 + (t42 * (t871 + t989) - t95 * rSges(6,3) * t926) * t635 + ((-t648 * t96 - t651 * t95) * pkin(1) + t95 * t744 * t636 + (t95 * (-pkin(7) - t1085) + t96 * t744) * t635) * qJD(1)) * m(6) + ((t277 + t330) * t635 + (t278 + t331) * t636) * t919 / 0.2e1 + (t21 * (t854 + t1201) + t80 * (-t651 * t1072 - t798) + t20 * (t1141 + t1205) + t81 * (t621 - t909 + t1137) + (-t20 * t990 + t21 * t574 + t81 * t702 + (t549 + (t1058 * t647 + t746) * qJD(1)) * t80) * t636 + (t21 * (-pkin(2) + t990 - t1073) + t80 * t748 + (t1058 * t650 + t820) * t1033 + (t80 * (-pkin(7) - t574) + t81 * (t746 - t1073)) * qJD(1)) * t635 - (t1197 + t667 - t80) * t81 - t848 * t436) * m(7) + (t85 + t68) * t1099; m(4) * t122 + m(5) * t65 + m(6) * t22 + m(7) * t19; ((-t1006 * t622 + t933) * t635 + (t681 + (-t1118 * t636 + (t1007 + t704) * t635) * qJD(3)) * t636) * t863 + ((-t1007 * t928 - t933) * t636 + (t681 + (t704 * t635 + (-t1118 + t1006) * t636) * qJD(3)) * t635) * t860 + (-t96 * (-t1129 * t403 + t276 * t588 - t435 * t489 - t487 * t506 + t663) + t22 * t823 + (t1070 + (qJD(1) * t96 + t42) * t897) * t636 + (t22 * t272 + t41 * t897 + t96 * t901) * t635 + (t272 * t588 - t401 * t1129 + t436 * t489 + t901 * t636 + (-t505 + t932) * t487 + t1133) * t95 + (-t272 * t506 + t276 * t505 + t401 * t435 + t403 * t436 + t1122 * t636 + (t146 + (-t530 + t972) * qJD(1)) * t635 + t1138) * t86) * m(6) + (((t410 * t646 - t412 * t649 + t304) * t542 - (t409 * t646 - t411 * t649 + t302) * t543 + (-t514 * t646 + t516 * t649 + t511) * t623 + t210 * qJD(4)) * t647 + (qJD(4) * t765 - t1115) * t650) * t1093 + ((t410 * t736 - t412 * t510) * t542 - (t409 * t736 - t411 * t510) * t543 + (t510 * t516 - t514 * t736) * t623 + (t1000 * t115 + t177 * t647) * qJD(4) + ((qJD(4) * t116 + t709) * t650 + t657) * t636) * t1101 + ((t410 * t507 - t412 * t508) * t542 - (t409 * t507 - t411 * t508) * t543 + (-t507 * t514 + t508 * t516) * t623 + (t114 * t994 + t176 * t647) * qJD(4) + ((qJD(4) * t113 + t709) * t650 + t657) * t635) * t1098 + (-(t238 * t523 - t1014) * qJD(1) - (t237 * (-t523 * t635 - t524 * t636) + t764 * t805) * qJD(3) + t122 * t754 + t237 * t668 + t764 * t567 + (-t174 * t635 - t175 * t636 + (-t239 * t636 + t1016) * qJD(1)) * t600) * m(4) + (t19 * t823 + (-qJD(1) * t852 + t21 * t821 + t1150) * t636 + (-qJD(1) * t849 - t19 * t980 + t20 * t821) * t635 - t850 * t588 - (t848 - t852) * t506 - (-t849 - t851) * t505 + (-t1129 * t970 + t1136 * t635 + t1210 * t435 - t979 * t588 + t821 * t931 - t663) * t81 + (t1129 * t971 + t1136 * t636 + t1210 * t436 + t1133) * t80 + (t985 * t636 + (-t983 + (-t530 + t906) * qJD(1)) * t635 - t922 - t970 * t436 - t971 * t435 + t1138) * t66) * m(7) + (-t776 / 0.2e1 - t778 / 0.2e1) * t312 - (t635 * t67 + t636 * t68) * t924 / 0.2e1 + (-t770 / 0.2e1 - t768 / 0.2e1) * t563 + (qJD(1) * t765 - t47 * t636 + t48 * t635) * t1092 - t1185 * t588 / 0.2e1 + (qJD(1) * t1169 + t1186 * t635 - t1187 * t636) * t1095 + (qJD(1) * t1167 + t1188 * t635 - t1189 * t636) * t1105 - t429 * t772 / 0.2e1 - t428 * t774 / 0.2e1 + (qJD(1) * t149 + t17 + (t1121 * t635 ^ 2 - t1120 * t1178 + (t182 * t635 + t183 * t636) * qJD(1)) * t1153 + t1193) * t635 / 0.2e1 + qJD(1) * (-t127 * t636 + t128 * t635 + (t277 * t635 + t278 * t636) * qJD(1)) / 0.2e1 + (-t1208 * t650 + ((t1142 * t638 - t1143 * t637 + t1204) * t1129 + (t1145 * t638 - t1147 * t637 - t1232) * t436 + (-t1144 * t638 + t1146 * t637 + t1206) * t435) * t647 + t1175 * t588 + t1211 * t506 + t1212 * t505) * t1096 + (t1139 * t636 + t1176 * t588 + (t1142 * t473 + t1143 * t472) * t1129 + t1213 * t506 + t1214 * t505 + (t1145 * t473 + t1147 * t472) * t436 + (-t1144 * t473 - t1146 * t472) * t435) * t1107 + (t1139 * t635 + t1177 * t588 + (t1142 * t471 - t1143 * t470) * t1129 + t1215 * t506 + t1216 * t505 + (t1145 * t471 - t1147 * t470) * t436 + (-t1144 * t471 + t1146 * t470) * t435) * t1104 + (qJD(1) * t1168 + t1190 * t635 - t1191 * t636) * t1106 - (qJD(1) * t150 + t18 + (-t1121 * t1178 + t636 ^ 2 * t1120 + (t180 * t635 + t181 * t636) * qJD(1)) * t1153 + t1192) * t636 / 0.2e1 + (t147 * t547 + t65 * t948 + (t147 * t963 + t65 * t315 + (qJD(1) * t148 + t83) * t942) * t636 + (qJD(1) * t147 * t525 + t148 * t963 + t65 * t313 + t82 * t942) * t635 - t147 * (t413 * t623 - t526 * t543 + t725) - t148 * (-t414 * t623 - t526 * t542 + t808) - ((-t147 * t313 + t148 * t315) * t647 + t662 * t650) * qJD(4) + ((qJD(1) * t313 + t162) * t636 + (qJD(1) * t966 + t163) * t635 + t413 * t542 + t414 * t543 + t1162) * t121) * m(5) - t766 * t867 / 0.2e1 + (qJD(1) * t773 - t45 * t636 + t46 * t635) * t1099 + (qJD(1) * t771 - t43 * t636 + t44 * t635) * t1100 + (-t780 / 0.2e1 - t782 / 0.2e1) * t311 + (t68 + t717 + t99 + t1151) * t864 - t1151 * t506 / 0.2e1 + (t718 + t98 + t67 + t1152) * t932 / 0.2e1 - t1152 * t505 / 0.2e1 - qJD(1) * ((t647 * t938 + t650 * t937) * qJD(1) + ((t635 * t951 - t636 * t952) * t650 + (t635 * t953 + t636 * t954) * t647) * qJD(3)) / 0.2e1 - t75 * t925 / 0.2e1; ((qJD(3) * t765 - t100) * t650 + (qJD(1) * t766 + qJD(3) * t210 + t47 * t635 + t48 * t636) * t647) * t1092 + ((qJD(3) * t771 - t84) * t650 + (qJD(1) * t772 + qJD(3) * t177 + t43 * t635 + t44 * t636) * t647) * t1100 + (t1021 + t1022 + t1059 + t1064 - t1065) * t1089 + t18 * t875 + t17 * t874 + (-t177 * t650 + t647 * t771) * t1108 + (t706 * t650 + (-t1114 * t646 + t649 * t671) * t647) * t1093 + (-t1114 * t736 + t671 * t510 - t636 * t698) * t1101 + (-t1114 * t507 + t508 * t671 - t635 * t698) * t1098 + ((qJD(3) * t773 - t85) * t650 + (qJD(1) * t774 + qJD(3) * t176 + t45 * t635 + t46 * t636) * t647) * t1099 + t656 + (-t176 * t650 + t647 * t773) * t1109 + (t75 + qJD(4) * (-t210 * t650 + t647 * t765)) * t859 + t1166 * t68 + t1165 * t67 + (-t81 * (t289 * t1129 + t244 + (-t492 - t537) * t435 + t955) + t21 * (t822 + t958) + t81 * (t825 + t969) + t19 * (-t270 - t981) + (t20 * t906 + t81 * t908 + (t1063 * t898 + t635 * t811) * qJD(3)) * t650 + (t907 * t1033 + (qJD(1) * t811 + t20 * t898 + t81 * t905) * t636 + t19 * t906 * t635) * t647 + (t436 * t492 - (-t288 - t348) * t1129 + t1135 + t1148) * t80 + (-t288 * t435 - t207 - (t289 + t350) * t436 - t840 + t828 + t982 + (qJD(1) * t907 + t908) * t1002) * t66) * m(7) + (t42 * (t957 + t958) + t96 * (t969 + t977) + t22 * (t236 - t270) + t86 * (t982 + t984) + (t41 * t972 + t96 * (-t144 - t167) + (t636 * t947 * t96 + t635 * t847) * qJD(3)) * t650 + ((qJD(1) * t847 + t41 * t947 + t96 * t965) * t636 + (t22 * t972 + t86 * (qJD(1) * t322 - t1122 - t167)) * t635) * t647 - t96 * (t842 + t955) - t86 * (t840 + t843) + ((-t272 + t322) * t927 + t1134 + t1135) * t95) * m(6) + ((qJD(3) * t662 + t147 * t163 - t148 * t162 + t83 * t313 - t82 * t315) * t650 + (t147 * (-qJD(3) * t313 + t383 * t635) + t148 * (qJD(3) * t315 - t383 * t636) + t65 * t757 + t121 * (-t162 * t635 + t163 * t636 - t313 * t932 - t315 * t931) + (t83 * t635 - t82 * t636 + (t148 * t635 + t1020) * qJD(1)) * t525) * t647 - t147 * (-t371 * t623 - t543 * t555) - t148 * (t372 * t623 - t542 * t555) - t121 * (t371 * t542 + t372 * t543)) * m(5); t656 + (-t81 * (t459 * t1129 + t244 + (-t537 + t915) * t435) - t66 * (-t435 * t458 + t207 + (t350 + t459) * t436) + t21 * t822 + t81 * t825 - t19 * t981 + t66 * t828 + (-t20 * t979 - t81 * t985 + (-t635 * t851 + t636 * t848) * qJD(3)) * t650 + (qJD(3) * t850 + (-qJD(1) * t851 + t1171 * t20 - t1174 * t81) * t636 + (-t1150 + t66 * (qJD(1) * t980 - t985)) * t635) * t647 + (-t436 * t915 - (-t348 + t458) * t1129 + t1148) * t80) * m(7) + (t42 * t957 + t41 * (-t276 * t650 - t487 * t996) + t22 * t236 - t1070 * t1002 + (-t144 * t650 + (-t333 * t647 - t487 * t926) * t636 + t977 - t842) * t96 + (-t272 * t927 + t1134) * t95 + (-t276 * t891 + t984 + (-t1122 * t647 - t276 * t926) * t635 - t843) * t86) * m(6); 0.2e1 * ((t622 * t81 + t80 * t928 - t19) * t1112 + (t435 * t81 + t436 * t80) * t1113) * t650 + 0.2e1 * ((qJD(3) * t66 + t20 * t635 + t21 * t636 - t80 * t932 + t81 * t931) * t1112 + (t66 * (t435 * t635 + t436 * t636) + (-t635 * t80 + t1063) * t1129) * t1113) * t647;];
tauc  = t1(:);
