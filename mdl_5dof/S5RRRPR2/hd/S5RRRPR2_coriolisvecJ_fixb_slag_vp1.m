% Calculate vector of centrifugal and Coriolis load on the joints for
% S5RRRPR2
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d2,d3,d5,theta4]';
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
% Datum: 2020-01-03 12:08
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauc = S5RRRPR2_coriolisvecJ_fixb_slag_vp1(qJ, qJD, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(9,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRRPR2_coriolisvecJ_fixb_slag_vp1: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RRRPR2_coriolisvecJ_fixb_slag_vp1: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5RRRPR2_coriolisvecJ_fixb_slag_vp1: pkin has to be [9x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RRRPR2_coriolisvecJ_fixb_slag_vp1: m has to be [6x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [6,3]), ...
  'S5RRRPR2_coriolisvecJ_fixb_slag_vp1: rSges has to be [6x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [6 6]), ...
  'S5RRRPR2_coriolisvecJ_fixb_slag_vp1: Icges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From coriolisvec_joint_fixb_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2020-01-03 12:07:22
% EndTime: 2020-01-03 12:07:28
% DurationCPUTime: 3.62s
% Computational Cost: add. (10345->310), mult. (5823->397), div. (0->0), fcn. (4286->10), ass. (0->206)
t166 = sin(qJ(1));
t281 = pkin(1) * qJD(1);
t241 = t166 * t281;
t164 = qJ(1) + qJ(2);
t155 = sin(t164);
t163 = qJD(1) + qJD(2);
t260 = t155 * t163;
t246 = pkin(2) * t260;
t192 = t241 + t246;
t154 = qJD(3) + t163;
t158 = qJ(3) + t164;
t149 = sin(t158);
t150 = cos(t158);
t309 = t149 * rSges(4,1) + t150 * rSges(4,2);
t316 = t309 * t154;
t59 = t192 + t316;
t148 = pkin(9) + t158;
t140 = sin(t148);
t141 = cos(t148);
t165 = sin(qJ(5));
t167 = cos(qJ(5));
t157 = Icges(6,4) * t167;
t206 = -Icges(6,2) * t165 + t157;
t194 = t206 * t154;
t273 = Icges(6,4) * t165;
t123 = Icges(6,1) * t167 - t273;
t195 = t123 * t154;
t65 = -Icges(6,5) * t141 + t123 * t140;
t276 = t167 * t65;
t63 = -Icges(6,6) * t141 + t140 * t206;
t278 = t165 * t63;
t210 = t276 - t278;
t307 = Icges(6,1) * t165 + t157;
t299 = -Icges(6,5) * t154 + qJD(5) * t307;
t120 = Icges(6,2) * t167 + t273;
t302 = -Icges(6,6) * t154 + qJD(5) * t120;
t262 = t141 * t167;
t263 = t141 * t165;
t64 = Icges(6,4) * t262 - Icges(6,2) * t263 + Icges(6,6) * t140;
t115 = Icges(6,4) * t263;
t66 = Icges(6,1) * t262 + Icges(6,5) * t140 - t115;
t322 = -qJD(5) * t210 - t165 * (-t140 * t299 + t141 * t195) - t167 * (-t140 * t302 + t141 * t194) + t154 * (t165 * t66 + t167 * t64);
t142 = pkin(3) * t149;
t89 = t140 * rSges(5,1) + t141 * rSges(5,2);
t311 = t142 + t89;
t265 = t140 * t167;
t116 = rSges(6,1) * t265;
t117 = rSges(6,2) * t263;
t127 = rSges(6,1) * t165 + rSges(6,2) * t167;
t249 = qJD(5) * t141;
t235 = t127 * t249;
t191 = -t235 - t246;
t310 = t140 * pkin(4) + t142;
t266 = t140 * t165;
t242 = rSges(6,2) * t266;
t222 = t116 - t242;
t67 = -rSges(6,3) * t141 + t222;
t225 = -pkin(8) * t141 + t310 + t67;
t28 = t154 * t225 - t191 + t241;
t261 = t150 * t154;
t112 = pkin(3) * t261;
t250 = qJD(5) * t140;
t236 = t127 * t250;
t217 = t112 - t236;
t156 = cos(t164);
t259 = t156 * t163;
t126 = pkin(2) * t259;
t168 = cos(qJ(1));
t152 = t168 * t281;
t251 = t126 + t152;
t243 = rSges(6,1) * t262;
t68 = rSges(6,3) * t140 - t117 + t243;
t92 = t141 * pkin(4) + t140 * pkin(8);
t286 = t68 + t92;
t29 = t154 * t286 + t217 + t251;
t80 = t127 * t140;
t81 = t127 * t141;
t315 = -t28 * t80 - t29 * t81;
t170 = (t29 * (-t116 - t310) - t28 * t117) * t154 + t315 * qJD(5);
t280 = t154 * t29;
t264 = t141 * t154;
t267 = t140 * t154;
t256 = pkin(4) * t264 + pkin(8) * t267;
t284 = rSges(6,3) * t267 + t154 * t243;
t58 = t154 * t68;
t305 = -t154 * t92 + t112 - t217 + t256 + t284 - t58;
t321 = t225 * t280 + t305 * t28 + t170;
t308 = t155 * rSges(3,1) + t156 * rSges(3,2);
t87 = t308 * t163;
t319 = -t241 - t87;
t253 = t307 + t206;
t254 = t120 - t123;
t317 = (t165 * t253 + t167 * t254) * t154;
t314 = 0.2e1 * qJD(5);
t282 = rSges(5,2) * t140;
t90 = t141 * rSges(5,1) - t282;
t79 = t154 * t90;
t99 = rSges(5,1) * t264;
t312 = -t79 + t99;
t283 = rSges(4,2) * t149;
t96 = t150 * rSges(4,1) - t283;
t85 = t154 * t96;
t60 = t251 + t85;
t219 = t311 * t154;
t204 = -t165 * t120 + t167 * t307;
t118 = Icges(6,5) * t165 + Icges(6,6) * t167;
t303 = -Icges(6,3) * t154 + qJD(5) * t118;
t106 = t206 * qJD(5);
t107 = t123 * qJD(5);
t301 = qJD(5) * (t120 * t167 + t165 * t307) + t106 * t165 - t107 * t167 - t118 * t154;
t51 = t192 + t219;
t52 = t112 + t251 + t79;
t171 = (-t51 * t282 - t52 * t311) * t154;
t297 = t52 * t219 + t312 * t51 + t171;
t153 = t154 ^ 2;
t294 = t154 / 0.2e1;
t293 = pkin(2) * t163 ^ 2;
t292 = pkin(3) * t153;
t160 = t166 * pkin(1);
t161 = t168 * pkin(1);
t290 = t140 * t307 + t63;
t289 = t141 * t307 + t64;
t288 = -t120 * t140 + t65;
t287 = -Icges(6,2) * t262 - t115 + t66;
t285 = rSges(6,3) * t264 + t154 * t242;
t277 = t165 * t64;
t268 = t118 * t140;
t50 = t141 * t204 + t268;
t275 = t50 * t154;
t74 = t118 * t141;
t119 = Icges(6,5) * t167 - Icges(6,6) * t165;
t193 = t119 * t154;
t169 = qJD(1) ^ 2;
t151 = t169 * t161;
t252 = t156 * t293 + t151;
t248 = qJD(5) * t165;
t247 = qJD(5) * t167;
t245 = t169 * t160;
t101 = pkin(8) * t264;
t244 = t101 + t285;
t240 = t150 * t292 + t252;
t146 = pkin(2) * t155;
t237 = t146 + t309;
t230 = t250 / 0.2e1;
t228 = t249 / 0.2e1;
t104 = rSges(3,1) * t156 - rSges(3,2) * t155;
t84 = t104 * t163 + t152;
t226 = t146 + t311;
t88 = rSges(3,1) * t259 - rSges(3,2) * t260;
t72 = rSges(4,1) * t261 - t154 * t283;
t147 = pkin(2) * t156;
t221 = t147 + t96;
t143 = pkin(3) * t150;
t220 = t143 + t90;
t213 = rSges(6,1) * t167 - rSges(6,2) * t165;
t212 = -t140 * t29 + t141 * t28;
t211 = t140 * t67 + t141 * t68;
t41 = t165 * t65 + t167 * t63;
t208 = -t167 * t66 + t277;
t203 = t126 + t72;
t202 = t147 + t220;
t55 = t65 * t265;
t61 = -Icges(6,3) * t141 + t119 * t140;
t24 = -t141 * t61 - t266 * t63 + t55;
t56 = t66 * t265;
t62 = Icges(6,5) * t262 - Icges(6,6) * t263 + Icges(6,3) * t140;
t25 = t141 * t62 + t266 * t64 - t56;
t198 = (-t140 * t25 - t141 * t24) * qJD(5);
t57 = t63 * t263;
t26 = -t140 * t61 - t262 * t65 + t57;
t27 = t140 * t62 - t208 * t141;
t197 = (-t140 * t27 - t141 * t26) * qJD(5);
t196 = -t155 * t293 - t245;
t187 = t165 * t288 + t167 * t290;
t186 = t165 * t287 + t167 * t289;
t184 = -t140 * t193 - t141 * t303 + t154 * t208;
t183 = t140 * t303 - t141 * t193 + t154 * t210;
t182 = -t119 * qJD(5) + t154 * t204;
t181 = -t149 * t292 + t196;
t180 = t143 + t286;
t179 = (-rSges(6,3) - pkin(8)) * t141 + t222 + t310;
t43 = rSges(6,2) * t141 * t247 + (t141 * t248 + t154 * t265) * rSges(6,1) - t285;
t44 = -rSges(6,1) * t140 * t248 + (-t140 * t247 - t154 * t263) * rSges(6,2) + t284;
t177 = (t154 * t67 - t43) * t141 + (t44 - t58) * t140;
t176 = t147 + t180;
t175 = t146 + t179;
t10 = t197 - t275;
t15 = qJD(5) * t208 + t165 * (t140 * t195 + t141 * t299) + t167 * (t140 * t194 + t141 * t302);
t20 = t182 * t140 + t141 * t301;
t21 = -t140 * t301 + t182 * t141;
t49 = t140 * t204 - t74;
t46 = t49 * t154;
t9 = t46 + t198;
t174 = (t46 + ((t26 + t56 - t57 + (t61 - t277) * t140) * t140 + (-t55 - t27 + (-t208 + t61) * t141 + (t276 + t278) * t140) * t141) * qJD(5)) * t230 + (qJD(5) * t204 + t106 * t167 + t107 * t165) * t154 + (t10 + t275 + ((-t25 + t57 + (t62 - t276) * t141) * t141 + (t24 - t55 + (t62 + t278) * t140) * t140) * qJD(5)) * t228 - (t15 + t20 + t9) * t250 / 0.2e1 - (t21 - t322) * t249 / 0.2e1 + (t141 * t50 + (t41 + t49) * t140) * qJD(5) * t294;
t109 = t213 * qJD(5);
t70 = t163 * t88 + t151;
t69 = -t163 * t87 - t245;
t54 = t154 * t72 + t252;
t53 = -t154 * t316 + t196;
t48 = t154 * (-rSges(5,2) * t267 + t99) + t240;
t47 = -t153 * t89 + t181;
t32 = qJD(5) * t211 + qJD(4);
t19 = t109 * t249 + (t44 - t236 + t256) * t154 + t240;
t18 = -t109 * t250 + (-pkin(4) * t267 + t101 - t235 - t43) * t154 + t181;
t11 = t177 * qJD(5);
t1 = [m(3) * (t69 * (t104 + t161) + t70 * (t160 + t308) + (t84 - t152 - t88) * t319) + t174 + (t18 * (t161 + t176) + t29 * (-t192 + t244) + t19 * (t160 + t175) + t170 + (t29 + t305) * t28) * m(6) + (t47 * (t161 + t202) - t52 * t192 + t48 * (t160 + t226) + t171 + (t52 + t312) * t51) * m(5) + m(4) * (t53 * (t161 + t221) + t54 * (t160 + t237) + (-t60 + t152 + t203) * t59); t174 + (t19 * t175 + t18 * t176 + (-t191 + t244 - t246) * t29 + t321) * m(6) + (t47 * t202 + t48 * t226 + t297) * m(5) + (t53 * t221 + t54 * t237 + (t203 - t126 - t85) * t59) * m(4) + (-(-t104 * t319 - t308 * t84) * t163 + t308 * t70 + t104 * t69 - t319 * t88 - t84 * t87) * m(3); t174 + (t19 * t179 + t18 * t180 + (t235 + t244) * t29 + t321) * m(6) + (t47 * t220 + t311 * t48 + t297) * m(5) + (-(-t309 * t60 + t59 * t96) * t154 + t53 * t96 + t54 * t309 + t59 * t72 - t60 * t316) * m(4); m(6) * t11; (t322 * t141 + (t154 * t41 - t15) * t140) * t294 - t154 * ((-t254 * t165 + t253 * t167) * t154 + ((t140 * t287 - t141 * t288) * t167 + (-t140 * t289 + t141 * t290) * t165) * qJD(5)) / 0.2e1 + ((-t249 * t268 - t193) * t141 + (-t317 + (-t186 * t140 + (t74 + t187) * t141) * qJD(5)) * t140) * t228 + ((t74 * t250 - t193) * t140 + (t317 + (-t187 * t141 + (-t268 + t186) * t140) * qJD(5)) * t141) * t230 - (t154 * t20 + ((-t183 * t140 - t154 * t27) * t141 + (-t184 * t140 + t154 * t26) * t140) * t314) * t140 / 0.2e1 - (t154 * t21 + ((-t183 * t141 - t154 * t25) * t141 + (-t184 * t141 + t154 * t24) * t140) * t314) * t141 / 0.2e1 + (t9 + t198) * t267 / 0.2e1 - (t10 + t197) * t264 / 0.2e1 + (t11 * t211 + t32 * t177 + t212 * t109 + ((t19 - t280) * t141 + (-t154 * t28 - t18) * t140) * t127 - t315 * t154 - (t32 * (-t140 * t80 - t141 * t81) + t212 * t213) * qJD(5)) * m(6);];
tauc = t1(:);
