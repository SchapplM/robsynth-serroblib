% Calculate minimal parameter regressor of coriolis joint torque vector for
% S6RPRPRR9
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [13x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,alpha3,d1,d3,d5,d6,theta2,theta4]';
% 
% Output:
% tauc_reg [6x30]
%   minimal parameter regressor of coriolis joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 04:06
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauc_reg = S6RPRPRR9_coriolisvecJ_fixb_regmin_slag_vp(qJ, qJD, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(13,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRPRR9_coriolisvecJ_fixb_regmin_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RPRPRR9_coriolisvecJ_fixb_regmin_slag_vp: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [13 1]), ...
  'S6RPRPRR9_coriolisvecJ_fixb_regmin_slag_vp: pkin has to be [13x1] (double)');

%% Symbolic Calculation
% From coriolisvec_joint_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 04:05:07
% EndTime: 2019-03-09 04:05:23
% DurationCPUTime: 5.74s
% Computational Cost: add. (14197->451), mult. (49705->662), div. (0->0), fcn. (43794->14), ass. (0->232)
t197 = sin(pkin(6));
t206 = cos(qJ(3));
t196 = sin(pkin(7));
t314 = cos(pkin(6));
t266 = t314 * t196;
t251 = t206 * t266;
t199 = cos(pkin(12));
t200 = cos(pkin(7));
t302 = t200 * t206;
t279 = t199 * t302;
t340 = t197 * t279 + t251;
t195 = sin(pkin(12));
t203 = sin(qJ(3));
t303 = t200 * t203;
t147 = (t195 * t206 + t199 * t303) * t197 + t203 * t266;
t141 = qJD(1) * t147;
t194 = sin(pkin(13));
t198 = cos(pkin(13));
t292 = qJD(1) * t197;
t276 = t195 * t292;
t336 = t340 * qJD(1) - t203 * t276;
t261 = -t141 * t194 + t198 * t336;
t330 = qJD(5) - t261;
t205 = cos(qJ(5));
t275 = t199 * t292;
t177 = t196 * t275;
t264 = qJD(1) * t314;
t226 = t200 * t264 - t177;
t220 = -qJD(3) - t226;
t158 = t205 * t220;
t202 = sin(qJ(5));
t236 = t198 * t141 + t194 * t336;
t91 = t202 * t236 + t158;
t339 = t330 * t91;
t163 = (-t194 * t203 + t198 * t206) * t196;
t90 = qJD(6) + t91;
t219 = t197 * (-t195 * t302 - t199 * t203);
t156 = qJD(1) * t219;
t304 = t199 * t206;
t157 = (-t195 * t303 + t304) * t292;
t338 = -qJD(3) * t163 + t156 * t194 + t157 * t198;
t128 = t336 * qJD(3);
t140 = t147 * qJD(3);
t129 = qJD(1) * t140;
t100 = t128 * t194 + t198 * t129;
t260 = t330 * t205;
t337 = -t100 * t202 - t330 * t260;
t201 = sin(qJ(6));
t204 = cos(qJ(6));
t93 = -t202 * t220 + t205 * t236;
t67 = t201 * t93 - t204 * t330;
t335 = t330 * t67;
t333 = t197 ^ 2 * (t195 ^ 2 + t199 ^ 2);
t331 = -qJD(5) + t330;
t256 = pkin(1) * t264;
t166 = qJ(2) * t275 + t195 * t256;
t306 = t197 * t199;
t217 = (t200 * t306 + t266) * pkin(9);
t130 = qJD(1) * t217 + t166;
t184 = t199 * t256;
t308 = t195 * t197;
t213 = t314 * pkin(2) + (-pkin(9) * t200 - qJ(2)) * t308;
t135 = t213 * qJD(1) + t184;
t159 = (-pkin(9) * t195 * t196 - pkin(2) * t199 - pkin(1)) * t197;
t152 = qJD(1) * t159 + qJD(2);
t310 = t152 * t196;
t215 = -(t135 * t200 + t310) * t203 - t206 * t130;
t84 = qJ(4) * t336 - t215;
t319 = t198 * t84;
t295 = t135 * t302 + t206 * t310;
t249 = -t203 * t130 + t295;
t83 = -t141 * qJ(4) + t249;
t71 = -t220 * pkin(3) + t83;
t36 = t194 * t71 + t319;
t34 = -t220 * pkin(10) + t36;
t116 = -t135 * t196 + t200 * t152;
t89 = -pkin(3) * t336 + qJD(4) + t116;
t46 = -pkin(4) * t261 - pkin(10) * t236 + t89;
t15 = t202 * t46 + t205 * t34;
t212 = t215 * qJD(3);
t232 = -qJ(4) * t128 - qJD(4) * t141;
t282 = qJD(1) * qJD(2);
t271 = t197 * t282;
t255 = t195 * t271;
t233 = t200 * t255;
t254 = t199 * t271;
t290 = qJD(3) * t206;
t272 = t196 * t290;
t273 = t200 * t290;
t278 = t135 * t273 + t152 * t272 + t206 * t254;
t210 = (-qJD(3) * t130 - t233) * t203 + t278;
t54 = -qJ(4) * t129 + qJD(4) * t336 + t210;
t27 = t198 * t54 + (-t203 * t254 - t206 * t233 + t212 + t232) * t194;
t101 = t128 * t198 - t129 * t194;
t171 = t196 * t255;
t124 = pkin(3) * t129 + t171;
t53 = pkin(4) * t100 - pkin(10) * t101 + t124;
t267 = t202 * t27 - t205 * t53;
t4 = -pkin(5) * t100 + t15 * qJD(5) + t267;
t329 = (pkin(5) * t93 + t90 * pkin(11)) * t90 + t4;
t288 = qJD(5) * t202;
t63 = -qJD(5) * t158 + t205 * t101 - t236 * t288;
t69 = t201 * t330 + t204 * t93;
t23 = t69 * qJD(6) - t204 * t100 + t201 * t63;
t277 = pkin(1) * t314;
t294 = qJ(2) * t306 + t195 * t277;
t144 = t217 + t294;
t187 = t199 * t277;
t148 = t187 + t213;
t235 = t148 * t200 + t159 * t196;
t328 = -t144 * t206 - t235 * t203;
t216 = qJD(2) * t219;
t214 = qJD(1) * t216;
t208 = t214 + t212;
t26 = t194 * t54 - t198 * (t208 + t232);
t301 = t202 * t101;
t64 = t93 * qJD(5) + t301;
t10 = pkin(5) * t64 - pkin(11) * t63 + t26;
t13 = pkin(11) * t330 + t15;
t79 = t194 * t84;
t35 = t198 * t71 - t79;
t33 = t220 * pkin(4) - t35;
t20 = t91 * pkin(5) - t93 * pkin(11) + t33;
t245 = t13 * t201 - t20 * t204;
t286 = qJD(5) * t205;
t224 = t202 * t53 + t205 * t27 + t46 * t286 - t34 * t288;
t3 = pkin(11) * t100 + t224;
t1 = -t245 * qJD(6) + t10 * t201 + t204 * t3;
t327 = t67 * t90;
t326 = t69 * t90;
t312 = t261 * t205;
t76 = t201 * t312 - t204 * t236;
t325 = t76 * t90;
t77 = t201 * t236 + t204 * t312;
t324 = t77 * t90;
t39 = t198 * t83 - t79;
t66 = pkin(3) * t141 + pkin(4) * t236 - pkin(10) * t261;
t323 = t202 * t66 + t205 * t39;
t265 = t314 * t200;
t167 = t196 * t306 - t265;
t85 = -pkin(3) * t167 - qJ(4) * t147 - t144 * t203 + t235 * t206;
t307 = t195 * t203;
t146 = t197 * t307 - t340;
t88 = -qJ(4) * t146 - t328;
t44 = t194 * t85 + t198 * t88;
t42 = -pkin(10) * t167 + t44;
t117 = t198 * t146 + t147 * t194;
t118 = -t146 * t194 + t147 * t198;
t119 = -t148 * t196 + t200 * t159;
t229 = pkin(3) * t146 + t119;
t49 = pkin(4) * t117 - pkin(10) * t118 + t229;
t241 = t202 * t49 + t205 * t42;
t322 = t236 * t91;
t321 = t236 * t93;
t320 = t261 * t69;
t284 = qJD(6) * t204;
t285 = qJD(6) * t201;
t22 = t201 * t100 + t204 * t63 + t284 * t330 - t93 * t285;
t318 = t201 * t22;
t317 = t201 * t64;
t316 = t204 * t64;
t263 = t204 * t90;
t38 = t194 * t83 + t319;
t315 = -t38 + t330 * (pkin(5) * t202 - pkin(11) * t205);
t164 = (t194 * t206 + t198 * t203) * t196;
t150 = t164 * t202 - t200 * t205;
t258 = t196 * t276;
t299 = t150 * qJD(5) + t202 * t258 + t338 * t205;
t298 = -qJD(3) * t164 - t198 * t156 + t157 * t194;
t151 = t164 * t205 + t200 * t202;
t296 = t151 * qJD(5) - t338 * t202 + t205 * t258;
t291 = qJD(2) * t197;
t289 = qJD(5) * t201;
t287 = qJD(5) * t204;
t281 = t90 * t289;
t280 = t90 * t287;
t190 = -pkin(3) * t198 - pkin(4);
t274 = t195 * t291;
t176 = t196 * t274;
t270 = pkin(3) * t140 + t176;
t211 = t148 * t273 + t159 * t272 + t291 * t304 + (-qJD(3) * t144 - t200 * t274) * t203;
t60 = -qJ(4) * t140 - qJD(4) * t146 + t211;
t139 = (t251 + (t279 - t307) * t197) * qJD(3);
t209 = qJD(3) * t328 + t216;
t61 = -qJ(4) * t139 - qJD(4) * t147 + t209;
t31 = t194 * t60 - t198 * t61;
t43 = -t194 * t88 + t198 * t85;
t169 = -pkin(5) * t205 - pkin(11) * t202 + t190;
t262 = pkin(11) * t236 - qJD(6) * t169 + t323;
t207 = qJD(1) ^ 2;
t253 = t197 * t207 * t314;
t247 = qJD(6) * t163 + t299;
t246 = qJD(6) * t151 + t298;
t6 = t13 * t204 + t20 * t201;
t17 = pkin(11) * t117 + t241;
t41 = pkin(4) * t167 - t43;
t94 = t118 * t202 + t167 * t205;
t95 = t118 * t205 - t167 * t202;
t21 = pkin(5) * t94 - pkin(11) * t95 + t41;
t244 = t17 * t204 + t201 * t21;
t243 = -t17 * t201 + t204 * t21;
t32 = t194 * t61 + t198 * t60;
t108 = t139 * t194 + t198 * t140;
t111 = t139 * t198 - t140 * t194;
t57 = pkin(4) * t108 - pkin(10) * t111 + t270;
t242 = -t202 * t32 + t205 * t57;
t14 = -t202 * t34 + t205 * t46;
t240 = -t202 * t42 + t205 * t49;
t239 = t204 * t117 - t201 * t95;
t75 = t117 * t201 + t204 * t95;
t234 = (-qJ(2) * t276 + t184) * t195 - t166 * t199;
t231 = t205 * t100 + (t202 * t261 - t288) * t330;
t228 = -t90 * t284 - t317;
t227 = t90 * t285 - t316;
t223 = t202 * t57 + t205 * t32 + t49 * t286 - t42 * t288;
t222 = -0.2e1 * t264 * t291;
t189 = pkin(3) * t194 + pkin(10);
t221 = -t189 * t100 + t33 * t330;
t12 = -pkin(5) * t330 - t14;
t218 = -pkin(11) * t64 + (t12 + t14) * t90;
t2 = -t6 * qJD(6) + t204 * t10 - t201 * t3;
t73 = t95 * qJD(5) + t111 * t202;
t72 = -t94 * qJD(5) + t111 * t205;
t65 = t69 * t288;
t29 = t75 * qJD(6) - t204 * t108 + t201 * t72;
t28 = t239 * qJD(6) + t108 * t201 + t204 * t72;
t18 = -pkin(5) * t236 + t202 * t39 - t205 * t66;
t16 = -pkin(5) * t117 - t240;
t11 = pkin(5) * t73 - pkin(11) * t72 + t31;
t8 = -pkin(5) * t108 + t241 * qJD(5) - t242;
t7 = pkin(11) * t108 + t223;
t5 = [0, 0, 0, t195 * t222, t199 * t222, 0.2e1 * t282 * t333 ((t199 * t294 + (qJ(2) * t308 - t187) * t195) * qJD(1) - t234) * t291, t128 * t147 + t139 * t141, -t128 * t146 - t129 * t147 + t139 * t336 - t140 * t141, -t128 * t167 - t139 * t220, t129 * t167 + t140 * t220, 0, t116 * t140 + t119 * t129 + t146 * t171 - t167 * t208 - t176 * t336 - t209 * t220, t116 * t139 + t119 * t128 + 0.2e1 * t141 * t176 + t167 * t210 + t211 * t220, -t100 * t44 - t101 * t43 - t108 * t36 - t111 * t35 - t117 * t27 + t118 * t26 + t236 * t31 + t261 * t32, t124 * t229 - t26 * t43 + t27 * t44 + t270 * t89 - t35 * t31 + t36 * t32, t63 * t95 + t72 * t93, -t63 * t94 - t64 * t95 - t72 * t91 - t73 * t93, t100 * t95 + t108 * t93 + t117 * t63 + t330 * t72, -t100 * t94 - t108 * t91 - t117 * t64 - t330 * t73, t100 * t117 + t108 * t330, t242 * t330 + t240 * t100 - t267 * t117 + t14 * t108 + t31 * t91 + t41 * t64 + t26 * t94 + t33 * t73 + (-t117 * t15 - t241 * t330) * qJD(5), -t241 * t100 - t15 * t108 - t224 * t117 - t223 * t330 + t26 * t95 + t31 * t93 + t33 * t72 + t41 * t63, t22 * t75 + t28 * t69, t22 * t239 - t23 * t75 - t28 * t67 - t29 * t69, t22 * t94 + t28 * t90 + t64 * t75 + t69 * t73, -t23 * t94 + t239 * t64 - t29 * t90 - t67 * t73, t64 * t94 + t73 * t90 (-qJD(6) * t244 + t11 * t204 - t201 * t7) * t90 + t243 * t64 + t2 * t94 - t245 * t73 + t8 * t67 + t16 * t23 - t4 * t239 + t12 * t29 -(qJD(6) * t243 + t11 * t201 + t204 * t7) * t90 - t244 * t64 - t1 * t94 - t6 * t73 + t8 * t69 + t16 * t22 + t4 * t75 + t12 * t28; 0, 0, 0, t195 * t253, t199 * t253, -t207 * t333, t234 * t292, 0, 0, 0, 0, 0, t200 * t129 + t156 * t220 + (qJD(3) * t203 * t220 + t276 * t336) * t196, t200 * t128 - t157 * t220 + (-t141 * t276 + t220 * t290) * t196, -t100 * t164 - t101 * t163 - t236 * t298 - t261 * t338, t124 * t200 - t163 * t26 + t164 * t27 - t258 * t89 + t298 * t35 - t338 * t36, 0, 0, 0, 0, 0, -t100 * t150 - t163 * t64 - t296 * t330 - t298 * t91, -t100 * t151 - t163 * t63 - t298 * t93 + t299 * t330, 0, 0, 0, 0, 0 (-t151 * t201 - t163 * t204) * t64 + t150 * t23 + (t201 * t247 - t204 * t246) * t90 + t296 * t67 -(t151 * t204 - t163 * t201) * t64 + t150 * t22 + (t201 * t246 + t204 * t247) * t90 + t296 * t69; 0, 0, 0, 0, 0, 0, 0, -t141 * t336, t141 ^ 2 - t336 ^ 2, t220 * t336 + t128, t141 * (qJD(3) - t177) + (t141 * t265 - t140) * qJD(1), 0, -t116 * t141 - t215 * t226 + t214, qJD(3) * t295 - t116 * t336 + t203 * t233 + t226 * t249 - t278 (-t100 * t194 - t101 * t198) * pkin(3) + (-t39 + t35) * t261 + (t36 - t38) * t236, t35 * t38 - t36 * t39 + (-t141 * t89 + t194 * t27 - t198 * t26) * pkin(3), t202 * t63 + t260 * t93 (t63 - t339) * t205 + (-t330 * t93 - t64) * t202, -t321 - t337, t231 + t322, -t330 * t236, -t14 * t236 + t190 * t64 - t38 * t91 + (-t26 + (-qJD(5) * t189 - t66) * t330) * t205 + (t330 * t39 + t221) * t202, t15 * t236 + t190 * t63 + t26 * t202 - t38 * t93 + (t189 * t288 + t323) * t330 + t221 * t205, t202 * t204 * t22 + (-t202 * t285 + t204 * t286 - t77) * t69, t67 * t77 + t69 * t76 + (-t201 * t69 - t204 * t67) * t286 + (-t318 - t204 * t23 + (t201 * t67 - t204 * t69) * qJD(6)) * t202, -t324 + t65 + (-t22 + t280) * t205 + (-t227 - t320) * t202, t325 + (t23 - t281) * t205 + (t228 - t335) * t202, t202 * t330 * t90 - t205 * t64, t169 * t316 - t12 * t76 - t18 * t67 + (t262 * t201 + t315 * t204) * t90 + (t12 * t289 - t2 + (qJD(5) * t67 + t228) * t189) * t205 + (t12 * t284 + t261 * t245 + t189 * t23 + t4 * t201 + (t189 * t201 * t90 - t245) * qJD(5)) * t202, -t169 * t317 - t12 * t77 - t18 * t69 + (-t315 * t201 + t262 * t204) * t90 + (t12 * t287 + t1 + (qJD(5) * t69 + t227) * t189) * t205 + (-t12 * t285 + t6 * t261 + t189 * t22 + t4 * t204 + (t189 * t263 - t6) * qJD(5)) * t202; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t236 ^ 2 - t261 ^ 2, t236 * t35 - t261 * t36 + t124, 0, 0, 0, 0, 0, t231 - t322, -t321 + t337, 0, 0, 0, 0, 0, t325 + (-t23 - t281) * t205 + (t228 + t335) * t202, t324 + t65 + (-t22 - t280) * t205 + (t227 - t320) * t202; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t93 * t91, -t91 ^ 2 + t93 ^ 2, t63 + t339, t331 * t93 - t301, t100, t331 * t15 - t33 * t93 - t267, t14 * t330 + t33 * t91 - t224, t263 * t69 + t318 (t22 - t327) * t204 + (-t23 - t326) * t201, t263 * t90 - t69 * t93 + t317, -t201 * t90 ^ 2 + t67 * t93 + t316, -t90 * t93, -pkin(5) * t23 - t15 * t67 + t218 * t201 - t204 * t329 + t245 * t93, -pkin(5) * t22 - t15 * t69 + t201 * t329 + t218 * t204 + t6 * t93; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t69 * t67, -t67 ^ 2 + t69 ^ 2, t22 + t327, -t23 + t326, t64, -t12 * t69 + t6 * t90 + t2, t12 * t67 - t245 * t90 - t1;];
tauc_reg  = t5;
