% Calculate vector of centrifugal and Coriolis load on the joints for
% S5RRPRR3
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d2,d4,d5,theta3]';
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
% tauc [5x1]
%   joint torques required to compensate Coriolis and centrifugal load

% Quelle: HybrDyn-Toolbox
% Datum: 2020-01-03 12:00
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauc = S5RRPRR3_coriolisvecJ_fixb_slag_vp1(qJ, qJD, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(9,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRPRR3_coriolisvecJ_fixb_slag_vp1: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RRPRR3_coriolisvecJ_fixb_slag_vp1: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5RRPRR3_coriolisvecJ_fixb_slag_vp1: pkin has to be [9x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RRPRR3_coriolisvecJ_fixb_slag_vp1: m has to be [6x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [6,3]), ...
  'S5RRPRR3_coriolisvecJ_fixb_slag_vp1: rSges has to be [6x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [6 6]), ...
  'S5RRPRR3_coriolisvecJ_fixb_slag_vp1: Icges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From coriolisvec_joint_fixb_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2020-01-03 12:00:21
% EndTime: 2020-01-03 12:00:27
% DurationCPUTime: 3.44s
% Computational Cost: add. (10079->307), mult. (5719->396), div. (0->0), fcn. (4238->10), ass. (0->201)
t162 = qJD(1) + qJD(2);
t163 = qJ(1) + qJ(2);
t154 = pkin(9) + t163;
t146 = sin(t154);
t155 = sin(t163);
t148 = pkin(2) * t155;
t302 = pkin(3) * t146 + t148;
t202 = t302 * t162;
t165 = sin(qJ(1));
t278 = pkin(1) * qJD(1);
t240 = t165 * t278;
t181 = t240 + t202;
t153 = qJD(4) + t162;
t150 = qJ(4) + t154;
t142 = sin(t150);
t143 = cos(t150);
t303 = t142 * rSges(5,1) + t143 * rSges(5,2);
t312 = t303 * t153;
t51 = t181 + t312;
t164 = sin(qJ(5));
t166 = cos(qJ(5));
t157 = Icges(6,4) * t166;
t207 = -Icges(6,2) * t164 + t157;
t197 = t207 * t153;
t271 = Icges(6,4) * t164;
t123 = Icges(6,1) * t166 - t271;
t198 = t123 * t153;
t65 = -Icges(6,5) * t143 + t123 * t142;
t274 = t166 * t65;
t63 = -Icges(6,6) * t143 + t142 * t207;
t276 = t164 * t63;
t211 = t274 - t276;
t300 = Icges(6,1) * t164 + t157;
t293 = -Icges(6,5) * t153 + qJD(5) * t300;
t120 = Icges(6,2) * t166 + t271;
t296 = -Icges(6,6) * t153 + qJD(5) * t120;
t260 = t143 * t166;
t261 = t143 * t164;
t64 = Icges(6,4) * t260 - Icges(6,2) * t261 + Icges(6,6) * t142;
t112 = Icges(6,4) * t261;
t66 = Icges(6,1) * t260 + Icges(6,5) * t142 - t112;
t316 = -qJD(5) * t211 - t164 * (-t293 * t142 + t143 * t198) - t166 * (-t296 * t142 + t143 * t197) + t153 * (t164 * t66 + t166 * t64);
t147 = cos(t154);
t95 = t146 * rSges(4,1) + t147 * rSges(4,2);
t305 = t148 + t95;
t156 = cos(t163);
t301 = t155 * rSges(3,1) + t156 * rSges(3,2);
t87 = t301 * t162;
t315 = -t240 - t87;
t252 = t300 + t207;
t253 = t120 - t123;
t313 = (t164 * t252 + t166 * t253) * t153;
t238 = t153 * t261;
t262 = t143 * t153;
t265 = t142 * t153;
t255 = pkin(4) * t262 + pkin(8) * t265;
t242 = rSges(6,1) * t260;
t280 = rSges(6,3) * t265 + t153 * t242;
t68 = -rSges(6,2) * t261 + rSges(6,3) * t142 + t242;
t58 = t153 * t68;
t92 = t143 * pkin(4) + t142 * pkin(8);
t311 = -rSges(6,2) * t238 - t153 * t92 + t255 + t280 - t58;
t135 = t142 * pkin(4);
t127 = rSges(6,1) * t164 + rSges(6,2) * t166;
t246 = qJD(5) * t143;
t233 = t127 * t246;
t264 = t142 * t164;
t241 = rSges(6,2) * t264;
t263 = t142 * t166;
t220 = rSges(6,1) * t263 - t241;
t67 = -rSges(6,3) * t143 + t220;
t192 = -t153 * (-pkin(8) * t143 + t135 + t67) - t233;
t28 = t181 - t192;
t167 = cos(qJ(1));
t152 = t167 * t278;
t184 = t68 + t92;
t247 = qJD(5) * t142;
t234 = t127 * t247;
t258 = t147 * t162;
t117 = pkin(3) * t258;
t256 = t156 * t162;
t126 = pkin(2) * t256;
t254 = t117 + t126;
t203 = -t234 + t254;
t29 = t153 * t184 + t152 + t203;
t80 = t127 * t142;
t81 = t127 * t143;
t310 = -t28 * t80 - t29 * t81;
t309 = 0.2e1 * qJD(5);
t115 = rSges(4,1) * t258;
t279 = rSges(4,2) * t146;
t96 = t147 * rSges(4,1) - t279;
t85 = t162 * t96;
t306 = t115 - t85;
t250 = t126 + t152;
t90 = t143 * rSges(5,1) - rSges(5,2) * t142;
t79 = t153 * t90;
t52 = t117 + t250 + t79;
t304 = t162 * t305;
t299 = t254 - t203 + t311;
t118 = Icges(6,5) * t164 + Icges(6,6) * t166;
t297 = -Icges(6,3) * t153 + qJD(5) * t118;
t106 = t207 * qJD(5);
t107 = t123 * qJD(5);
t295 = t106 * t164 - t107 * t166 - t118 * t153 + (t120 * t166 + t164 * t300) * qJD(5);
t161 = t162 ^ 2;
t289 = t153 / 0.2e1;
t288 = pkin(2) * t161;
t141 = pkin(3) * t147;
t159 = t165 * pkin(1);
t160 = t167 * pkin(1);
t286 = t142 * t300 + t63;
t285 = t143 * t300 + t64;
t284 = -t120 * t142 + t65;
t283 = -Icges(6,2) * t260 - t112 + t66;
t281 = -rSges(6,3) * t262 - t153 * t241;
t275 = t164 * t64;
t266 = t118 * t142;
t50 = -t120 * t261 + t260 * t300 + t266;
t273 = t50 * t153;
t74 = t118 * t143;
t119 = Icges(6,5) * t166 - Icges(6,6) * t164;
t196 = t119 * t153;
t259 = t146 * t162;
t257 = t155 * t162;
t168 = qJD(1) ^ 2;
t151 = t168 * t160;
t251 = t156 * t288 + t151;
t149 = pkin(2) * t156;
t248 = t141 + t149;
t245 = qJD(5) * t164;
t244 = qJD(5) * t166;
t243 = t168 * t159;
t239 = t153 * t263;
t237 = t161 * t141 + t251;
t228 = t247 / 0.2e1;
t226 = t246 / 0.2e1;
t104 = rSges(3,1) * t156 - rSges(3,2) * t155;
t84 = t104 * t162 + t152;
t224 = t302 + t303;
t72 = rSges(5,1) * t262 - rSges(5,2) * t265;
t221 = -pkin(4) * t265 + pkin(8) * t262;
t88 = rSges(3,1) * t256 - rSges(3,2) * t257;
t219 = t149 + t96;
t214 = rSges(6,1) * t166 - rSges(6,2) * t164;
t213 = -t142 * t29 + t143 * t28;
t212 = t142 * t67 + t143 * t68;
t41 = t164 * t65 + t166 * t63;
t209 = -t166 * t66 + t275;
t205 = -t120 * t164 + t166 * t300;
t204 = t90 + t248;
t201 = t72 + t254;
t53 = t65 * t263;
t61 = -Icges(6,3) * t143 + t119 * t142;
t24 = -t143 * t61 - t264 * t63 + t53;
t54 = t66 * t263;
t62 = Icges(6,5) * t260 - Icges(6,6) * t261 + Icges(6,3) * t142;
t25 = t143 * t62 + t264 * t64 - t54;
t200 = (-t142 * t25 - t143 * t24) * qJD(5);
t55 = t63 * t261;
t26 = -t142 * t61 - t260 * t65 + t55;
t27 = t142 * t62 - t209 * t143;
t199 = (-t142 * t27 - t143 * t26) * qJD(5);
t191 = t164 * t284 + t166 * t286;
t190 = t164 * t283 + t166 * t285;
t188 = -t142 * t196 - t297 * t143 + t153 * t209;
t187 = t297 * t142 - t143 * t196 + t153 * t211;
t185 = -t119 * qJD(5) + t153 * t205;
t183 = t135 + (-rSges(6,3) - pkin(8)) * t143 + t220;
t182 = -t161 * t302 - t243;
t179 = -rSges(6,1) * t239 + t221 - t281;
t43 = rSges(6,2) * t143 * t244 + (t143 * t245 + t239) * rSges(6,1) + t281;
t44 = -rSges(6,1) * t142 * t245 + (-t142 * t244 - t238) * rSges(6,2) + t280;
t178 = (t153 * t67 - t43) * t143 + (t44 - t58) * t142;
t177 = t184 + t248;
t176 = t183 + t302;
t10 = t199 - t273;
t15 = qJD(5) * t209 + t164 * (t142 * t198 + t293 * t143) + t166 * (t142 * t197 + t296 * t143);
t20 = t185 * t142 + t295 * t143;
t21 = -t295 * t142 + t185 * t143;
t49 = t142 * t205 - t74;
t46 = t49 * t153;
t9 = t46 + t200;
t175 = (t46 + ((t26 + t54 - t55 + (t61 - t275) * t142) * t142 + (-t53 - t27 + (-t209 + t61) * t143 + (t274 + t276) * t142) * t143) * qJD(5)) * t228 + (qJD(5) * t205 + t106 * t166 + t107 * t164) * t153 + (t273 + ((-t25 + t55 + (t62 - t274) * t143) * t143 + (t24 - t53 + (t62 + t276) * t142) * t142) * qJD(5) + t10) * t226 - (t15 + t20 + t9) * t247 / 0.2e1 - (t21 - t316) * t246 / 0.2e1 + (t143 * t50 + (t41 + t49) * t142) * qJD(5) * t289;
t59 = t240 + t304;
t60 = t250 + t85;
t171 = (-t59 * t279 - t60 * t305) * t162;
t170 = -pkin(2) * t257 - pkin(3) * t259 + t179;
t169 = t310 * qJD(5);
t108 = t214 * qJD(5);
t70 = t162 * t88 + t151;
t69 = -t162 * t87 - t243;
t57 = t162 * (-rSges(4,2) * t259 + t115) + t251;
t56 = -t155 * t288 - t161 * t95 - t243;
t48 = t153 * t72 + t237;
t47 = -t153 * t312 + t182;
t32 = qJD(5) * t212 + qJD(3);
t19 = t108 * t246 + (t44 - t234 + t255) * t153 + t237;
t18 = -t108 * t247 + (t221 - t43 - t233) * t153 + t182;
t11 = t178 * qJD(5);
t1 = [m(3) * (t69 * (t104 + t160) + t70 * (t159 + t301) + (t84 - t152 - t88) * t315) + t175 + (t18 * (t160 + t177) + t29 * (t170 - t240) + t19 * (t159 + t176) + t169 + (t29 + t299) * t28) * m(6) + m(5) * (t47 * (t160 + t204) + t48 * (t159 + t224) + (-t52 + t152 + t201) * t51) + (t56 * (t160 + t219) - t60 * t240 + t57 * (t159 + t305) + t171 + (t60 + t306) * t59) * m(4); t175 + (t19 * t176 + t18 * t177 + t169 + (t170 + t202 - t192) * t29 + t299 * t28) * m(6) + (t47 * t204 + t48 * t224 + (t201 - t79 - t254) * t51) * m(5) + (t56 * t219 + t60 * t304 + t305 * t57 + t306 * t59 + t171) * m(4) + (t301 * t70 + t104 * t69 - t315 * t88 - t84 * t87 - (-t104 * t315 - t301 * t84) * t162) * m(3); m(6) * t11; t175 + (t18 * t184 + t19 * t183 + t169 + (t179 - t192) * t29 + (t234 + t311) * t28) * m(6) + (-(-t303 * t52 + t51 * t90) * t153 + t47 * t90 + t48 * t303 + t51 * t72 - t52 * t312) * m(5); (t316 * t143 + (t153 * t41 - t15) * t142) * t289 - t153 * ((-t253 * t164 + t252 * t166) * t153 + ((t142 * t283 - t143 * t284) * t166 + (-t142 * t285 + t143 * t286) * t164) * qJD(5)) / 0.2e1 + ((-t246 * t266 - t196) * t143 + (-t313 + (-t190 * t142 + (t74 + t191) * t143) * qJD(5)) * t142) * t226 + ((t74 * t247 - t196) * t142 + (t313 + (-t191 * t143 + (-t266 + t190) * t142) * qJD(5)) * t143) * t228 - (t153 * t20 + ((-t187 * t142 - t153 * t27) * t143 + (-t188 * t142 + t153 * t26) * t142) * t309) * t142 / 0.2e1 - (t153 * t21 + ((-t187 * t143 - t153 * t25) * t143 + (-t188 * t143 + t153 * t24) * t142) * t309) * t143 / 0.2e1 + (t9 + t200) * t265 / 0.2e1 - (t10 + t199) * t262 / 0.2e1 + (t11 * t212 + t32 * t178 + t213 * t108 + ((-t153 * t29 + t19) * t143 + (-t153 * t28 - t18) * t142) * t127 - t310 * t153 - (t32 * (-t142 * t80 - t143 * t81) + t213 * t214) * qJD(5)) * m(6);];
tauc = t1(:);
