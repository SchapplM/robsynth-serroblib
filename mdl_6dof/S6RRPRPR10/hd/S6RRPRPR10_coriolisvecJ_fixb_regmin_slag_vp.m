% Calculate minimal parameter regressor of coriolis joint torque vector for
% S6RRPRPR10
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [11x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d1,d2,d4,d6,theta3]';
% 
% Output:
% tauc_reg [6x32]
%   minimal parameter regressor of coriolis joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 11:10
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauc_reg = S6RRPRPR10_coriolisvecJ_fixb_regmin_slag_vp(qJ, qJD, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(11,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPRPR10_coriolisvecJ_fixb_regmin_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRPRPR10_coriolisvecJ_fixb_regmin_slag_vp: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6RRPRPR10_coriolisvecJ_fixb_regmin_slag_vp: pkin has to be [11x1] (double)');

%% Symbolic Calculation
% From coriolisvec_joint_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 11:09:35
% EndTime: 2019-03-09 11:09:47
% DurationCPUTime: 5.68s
% Computational Cost: add. (8207->485), mult. (22175->662), div. (0->0), fcn. (17912->10), ass. (0->243)
t216 = cos(pkin(6));
t307 = qJD(1) * t216;
t205 = qJD(2) + t307;
t213 = sin(pkin(11));
t215 = cos(pkin(11));
t219 = sin(qJ(2));
t214 = sin(pkin(6));
t308 = qJD(1) * t214;
t290 = t219 * t308;
t158 = t205 * t215 - t213 * t290;
t159 = t205 * t213 + t215 * t290;
t218 = sin(qJ(4));
t345 = cos(qJ(4));
t245 = -t218 * t158 - t345 * t159;
t353 = qJD(6) - t245;
t360 = t353 ^ 2;
t109 = -t345 * t158 + t159 * t218;
t221 = cos(qJ(2));
t306 = qJD(1) * t221;
t289 = t214 * t306;
t195 = -qJD(4) + t289;
t217 = sin(qJ(6));
t220 = cos(qJ(6));
t83 = -t220 * t109 - t195 * t217;
t359 = t353 * t83;
t326 = t109 * t195;
t291 = t345 * t215;
t318 = t214 * t221;
t261 = t291 * t318;
t246 = qJD(1) * t261;
t297 = qJD(1) * qJD(2);
t284 = t214 * t297;
t274 = t221 * t284;
t258 = t213 * t274;
t286 = qJD(4) * t345;
t64 = t218 * (qJD(4) * t159 + t258) - qJD(2) * t246 - t158 * t286;
t358 = t64 - t326;
t276 = t213 * t289;
t302 = qJD(4) * t218;
t312 = t213 * t302 - t215 * t286 - t218 * t276 + t246;
t237 = t217 * t64 - t360 * t220;
t357 = t109 ^ 2;
t347 = t245 ^ 2;
t356 = t83 * t109;
t85 = t109 * t217 - t195 * t220;
t355 = t85 * t109;
t328 = t245 * t195;
t296 = pkin(1) * t307;
t172 = pkin(8) * t289 + t219 * t296;
t148 = qJ(3) * t205 + t172;
t167 = (-pkin(2) * t221 - qJ(3) * t219 - pkin(1)) * t214;
t153 = qJD(1) * t167;
t94 = -t148 * t213 + t215 * t153;
t53 = -pkin(3) * t289 - pkin(9) * t159 + t94;
t95 = t215 * t148 + t213 * t153;
t66 = pkin(9) * t158 + t95;
t27 = t218 * t66 - t345 * t53;
t315 = -qJD(5) - t27;
t266 = pkin(2) * t219 - qJ(3) * t221;
t170 = t266 * t308;
t171 = -pkin(8) * t290 + t221 * t296;
t120 = t213 * t170 + t215 * t171;
t101 = -pkin(9) * t276 + t120;
t341 = pkin(9) + qJ(3);
t190 = t341 * t213;
t191 = t341 * t215;
t119 = t215 * t170 - t171 * t213;
t316 = t215 * t221;
t240 = (pkin(3) * t219 - pkin(9) * t316) * t214;
t91 = qJD(1) * t240 + t119;
t354 = qJD(3) * t291 - t345 * t101 - t190 * t286 + (-qJD(3) * t213 - qJD(4) * t191 - t91) * t218;
t180 = t345 * t213 + t218 * t215;
t231 = t180 * t318;
t146 = qJD(1) * t231;
t176 = t180 * qJD(4);
t311 = t176 - t146;
t331 = qJ(5) * t290 - t354;
t145 = -t218 * t190 + t345 * t191;
t352 = qJD(3) * t180 + qJD(4) * t145 - t218 * t101 + t345 * t91;
t139 = pkin(3) * t276 + t172;
t351 = qJ(5) * t312 - t180 * qJD(5) - t139;
t350 = qJD(4) * t245;
t314 = -pkin(5) * t245 - t315;
t349 = -qJD(6) + t353;
t28 = t218 * t53 + t345 * t66;
t26 = qJ(5) * t195 - t28;
t342 = pkin(5) * t109;
t18 = -t26 - t342;
t346 = pkin(4) + pkin(10);
t348 = t346 * t64 + (t18 - t28 + t342) * t353;
t234 = qJD(2) * t240;
t151 = (qJD(2) * t266 - qJD(3) * t219) * t214;
t303 = qJD(2) * t219;
t288 = t214 * t303;
t295 = pkin(1) * qJD(2) * t216;
t255 = -pkin(8) * t288 + t221 * t295;
t157 = qJD(3) * t216 + t255;
t98 = t215 * t151 - t157 * t213;
t71 = t234 + t98;
t344 = pkin(1) * t219;
t166 = pkin(8) * t318 + (qJ(3) + t344) * t216;
t115 = -t166 * t213 + t215 * t167;
t319 = t214 * t219;
t174 = t213 * t216 + t215 * t319;
t73 = -pkin(3) * t318 - pkin(9) * t174 + t115;
t287 = qJD(2) * t318;
t275 = t213 * t287;
t99 = t213 * t151 + t215 * t157;
t86 = -pkin(9) * t275 + t99;
t116 = t215 * t166 + t213 * t167;
t293 = t213 * t319;
t317 = t215 * t216;
t248 = -t293 + t317;
t90 = pkin(9) * t248 + t116;
t247 = -t218 * t71 - t73 * t286 + t302 * t90 - t345 * t86;
t13 = -t214 * (qJ(5) * t303 - qJD(5) * t221) + t247;
t343 = pkin(1) * t221;
t340 = t218 * t73 + t345 * t90;
t338 = -t311 * pkin(5) - t331;
t141 = -pkin(2) * t205 + qJD(3) - t171;
t106 = -pkin(3) * t158 + t141;
t227 = qJ(5) * t245 + t106;
t35 = pkin(4) * t109 + t227;
t337 = t245 * t35;
t335 = t219 * t94;
t334 = t219 * t95;
t59 = t220 * t64;
t196 = t219 * t284;
t300 = qJD(6) * t220;
t301 = qJD(6) * t217;
t230 = qJD(2) * t231;
t65 = qJD(1) * t230 - t350;
t30 = t109 * t300 + t195 * t301 + t220 * t196 + t217 * t65;
t333 = t30 * t220;
t332 = t311 * pkin(4) + t351;
t330 = pkin(4) * t290 + t352;
t327 = t109 * qJ(5);
t325 = t245 * t109;
t277 = qJD(1) * t295;
t165 = pkin(8) * t274 + t219 * t277;
t324 = t165 * t213;
t179 = t213 * t218 - t291;
t323 = t179 * t217;
t322 = t179 * t220;
t210 = t214 ^ 2;
t321 = t210 * qJD(1) ^ 2;
t320 = t213 * t221;
t135 = qJD(1) * t151;
t242 = -pkin(8) * t196 + t221 * t277;
t136 = qJD(3) * t205 + t242;
t88 = t213 * t135 + t215 * t136;
t173 = pkin(8) * t287 + t219 * t295;
t309 = t219 ^ 2 - t221 ^ 2;
t144 = t345 * t190 + t218 * t191;
t305 = qJD(2) * t144;
t304 = qJD(2) * t145;
t294 = t221 * t321;
t292 = t220 * t318;
t131 = pkin(3) * t258 + t165;
t140 = pkin(3) * t275 + t173;
t209 = -pkin(3) * t215 - pkin(2);
t285 = t210 * t297;
t87 = t215 * t135 - t136 * t213;
t54 = qJD(1) * t234 + t87;
t67 = -pkin(9) * t258 + t88;
t283 = -t218 * t54 - t53 * t286 + t66 * t302 - t345 * t67;
t282 = t218 * t67 + t66 * t286 + t53 * t302 - t345 * t54;
t281 = t353 * t217;
t280 = t205 + t307;
t279 = 0.2e1 * t285;
t278 = t346 * t319;
t273 = -t218 * t90 + t345 * t73;
t271 = -qJD(5) * t195 - t283;
t254 = -qJ(5) * t180 + t209;
t114 = t346 * t179 + t254;
t270 = t312 * pkin(5) - qJD(1) * t278 + qJD(6) * t114 - t352;
t121 = t180 * pkin(5) + t144;
t269 = -qJD(6) * t121 - t311 * t346 - t351;
t268 = -0.2e1 * pkin(1) * t285;
t236 = qJ(5) * t64 + qJD(5) * t245 + t131;
t12 = t346 * t65 + t236;
t260 = qJD(2) * t278;
t7 = -pkin(5) * t64 - qJD(1) * t260 + t282;
t267 = t220 * t12 + t217 * t7;
t16 = t346 * t195 + t314;
t22 = t346 * t109 + t227;
t3 = t16 * t220 - t217 * t22;
t4 = t16 * t217 + t22 * t220;
t126 = t345 * t174 + t218 * t248;
t37 = pkin(4) * t318 - t273;
t23 = t126 * pkin(5) + pkin(10) * t318 + t37;
t232 = t345 * t248;
t125 = t174 * t218 - t232;
t206 = pkin(8) * t319;
t127 = pkin(3) * t293 + t206 + (t209 - t343) * t216;
t225 = -t126 * qJ(5) + t127;
t33 = t346 * t125 + t225;
t264 = -t217 * t33 + t220 * t23;
t263 = t217 * t23 + t220 * t33;
t262 = pkin(4) * t196;
t259 = qJ(5) * t196;
t36 = qJ(5) * t318 - t340;
t251 = -t195 * t28 - t282;
t250 = -t281 * t353 - t59;
t249 = t218 * t86 + t90 * t286 + t73 * t302 - t345 * t71;
t103 = t125 * t220 + t217 * t318;
t10 = -t259 - t271;
t6 = -pkin(5) * t65 - t10;
t244 = t6 + (t346 * t353 + t327) * t353;
t241 = t196 * t217 - t220 * t65;
t123 = -t220 * t146 + t217 * t290;
t239 = -t176 * t220 + t179 * t301 - t123;
t124 = t146 * t217 + t220 * t290;
t238 = t176 * t217 + t179 * t300 - t124;
t78 = -qJD(4) * t232 - qJD(2) * t261 + (qJD(4) * t174 + t275) * t218;
t235 = qJ(5) * t78 - qJD(5) * t126 + t140;
t11 = -t262 + t282;
t2 = -qJD(6) * t4 - t217 * t12 + t220 * t7;
t17 = pkin(4) * t65 + t236;
t229 = -qJ(3) * t303 + (-pkin(2) * qJD(2) + qJD(3) - t141) * t221;
t226 = -t64 - t326;
t224 = -t180 * t274 + t350;
t169 = t206 + (-pkin(2) - t343) * t216;
t128 = pkin(4) * t179 + t254;
t122 = -t179 * pkin(5) + t145;
t104 = t125 * t217 - t292;
t79 = qJD(4) * t126 + t230;
t46 = t64 * t180;
t45 = -pkin(4) * t245 + t327;
t43 = t125 * pkin(4) + t225;
t42 = t64 * t126;
t39 = qJD(6) * t103 + t217 * t79 + t220 * t288;
t38 = -qJD(6) * t292 - t220 * t79 + (qJD(6) * t125 + t288) * t217;
t31 = qJD(6) * t85 + t241;
t25 = pkin(4) * t195 - t315;
t24 = -pkin(5) * t125 - t36;
t21 = pkin(4) * t79 + t235;
t15 = t346 * t79 + t235;
t14 = -pkin(4) * t288 + t249;
t9 = -pkin(5) * t79 - t13;
t8 = -t78 * pkin(5) + t249 - t260;
t1 = qJD(6) * t3 + t267;
t5 = [0, 0, 0, t219 * t221 * t279, -t309 * t279, t280 * t287, -t280 * t288, 0, -t165 * t216 - t173 * t205 + t219 * t268, -t205 * t255 - t216 * t242 + t221 * t268, -t165 * t317 - t173 * t158 + (t219 * t324 + (-qJD(1) * t98 - t87) * t221 + (t141 * t320 + t335 + (t115 * t219 + t169 * t320) * qJD(1)) * qJD(2)) * t214, t159 * t173 + t165 * t174 + ((qJD(1) * t99 + t88) * t221 + (t141 * t316 - t334 + (-t116 * t219 + t169 * t316) * qJD(1)) * qJD(2)) * t214, t99 * t158 + t88 * t248 - t98 * t159 - t87 * t174 + (-t213 * t95 - t215 * t94 + (-t115 * t215 - t116 * t213) * qJD(1)) * t287, t115 * t87 + t116 * t88 + t141 * t173 + t165 * t169 + t94 * t98 + t95 * t99, t245 * t78 - t42, t109 * t78 + t125 * t64 - t126 * t65 + t245 * t79, t195 * t78 + (t221 * t64 + (qJD(1) * t126 - t245) * t303) * t214, t195 * t79 + (t221 * t65 + (-qJD(1) * t125 - t109) * t303) * t214 (-t195 * t214 - t210 * t306) * t303, t249 * t195 + t140 * t109 + t127 * t65 + t131 * t125 + t106 * t79 + (t282 * t221 + (qJD(1) * t273 - t27) * t303) * t214, -t247 * t195 - t140 * t245 - t127 * t64 + t131 * t126 - t106 * t78 + (-t283 * t221 + (-t340 * qJD(1) - t28) * t303) * t214, t10 * t125 + t109 * t13 + t11 * t126 - t14 * t245 - t25 * t78 + t26 * t79 + t36 * t65 - t37 * t64, -t109 * t21 - t125 * t17 - t14 * t195 - t35 * t79 - t43 * t65 + (-t11 * t221 + (qJD(1) * t37 + t25) * t303) * t214, t245 * t21 - t126 * t17 + t13 * t195 + t35 * t78 + t43 * t64 + (t10 * t221 + (-qJD(1) * t36 - t26) * t303) * t214, t10 * t36 + t11 * t37 + t13 * t26 + t14 * t25 + t17 * t43 + t21 * t35, t104 * t30 + t39 * t85, t103 * t30 - t104 * t31 - t38 * t85 - t39 * t83, -t104 * t64 + t126 * t30 + t353 * t39 - t78 * t85, -t103 * t64 - t126 * t31 - t353 * t38 + t78 * t83, -t353 * t78 - t42 (-qJD(6) * t263 - t217 * t15 + t220 * t8) * t353 - t264 * t64 + t2 * t126 - t3 * t78 + t9 * t83 + t24 * t31 - t6 * t103 + t18 * t38 -(qJD(6) * t264 + t220 * t15 + t217 * t8) * t353 + t263 * t64 - t1 * t126 + t4 * t78 + t9 * t85 + t24 * t30 + t6 * t104 + t18 * t39; 0, 0, 0, -t219 * t294, t309 * t321 (qJD(2) - t205) * t289, t205 * t290 - t196, 0, t172 * t205 + t321 * t344 - t165, pkin(1) * t294 + t171 * t205 - t242, t158 * t172 - t165 * t215 + (t119 * t221 + t213 * t229 - t335) * t308, -t159 * t172 + t324 + (-t120 * t221 + t215 * t229 + t334) * t308, t119 * t159 - t120 * t158 + (qJD(3) * t158 + t289 * t94 + t88) * t215 + (qJD(3) * t159 + t289 * t95 - t87) * t213, -pkin(2) * t165 - t119 * t94 - t120 * t95 - t141 * t172 + (-t213 * t94 + t215 * t95) * qJD(3) + (-t213 * t87 + t215 * t88) * qJ(3), t245 * t312 - t46, t109 * t312 + t179 * t64 - t180 * t65 + t245 * t311, t312 * t195 + (qJD(2) * t180 + t245) * t290, t311 * t195 + (-qJD(2) * t179 + t109) * t290, t195 * t290, -t139 * t109 + t131 * t179 + t209 * t65 + t352 * t195 + t311 * t106 + (t27 - t305) * t290, t139 * t245 + t131 * t180 - t209 * t64 + t354 * t195 - t312 * t106 + (t28 - t304) * t290, t10 * t179 + t331 * t109 + t11 * t180 - t144 * t64 - t145 * t65 - t245 * t330 - t312 * t25 + t311 * t26, -t128 * t65 - t17 * t179 - t311 * t35 - t330 * t195 - t332 * t109 + (-t25 + t305) * t290, t128 * t64 - t17 * t180 + t312 * t35 + t331 * t195 + t332 * t245 + (t26 + t304) * t290, -t10 * t145 + t11 * t144 + t128 * t17 + t330 * t25 + t331 * t26 + t332 * t35, t238 * t85 + t30 * t323, t123 * t85 + t124 * t83 + (-t217 * t83 + t220 * t85) * t176 + (-t217 * t31 + t333 + (-t217 * t85 - t220 * t83) * qJD(6)) * t179, t180 * t30 + t238 * t353 - t312 * t85 - t323 * t64, -t180 * t31 - t239 * t353 + t312 * t83 - t322 * t64, -t312 * t353 - t46 -(-t114 * t217 + t121 * t220) * t64 + t2 * t180 + t122 * t31 - t6 * t322 + t338 * t83 - t312 * t3 + (t217 * t269 - t220 * t270) * t353 + t239 * t18 (t114 * t220 + t121 * t217) * t64 - t1 * t180 + t122 * t30 + t6 * t323 + t338 * t85 + t312 * t4 + (t217 * t270 + t220 * t269) * t353 + t238 * t18; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 (qJD(2) * t213 - t159) * t289 (qJD(2) * t215 - t158) * t289, -t158 ^ 2 - t159 ^ 2, -t158 * t95 + t159 * t94 + t165, 0, 0, 0, 0, 0, t65 + t328, -t358, -t347 - t357, t224 - t328, t358, -t109 * t26 + t245 * t25 + t17, 0, 0, 0, 0, 0, t237 + t356, t217 * t360 + t355 + t59; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t325, t347 - t357, t226, t224 + t328, t196, t106 * t245 + t251, t106 * t109 + t195 * t27 + t283, pkin(4) * t64 - qJ(5) * t65 - (-t26 - t28) * t245 + (t25 + t315) * t109, t109 * t45 - t251 - 0.2e1 * t262 - t337, -t109 * t35 + t195 * t315 - t245 * t45 + 0.2e1 * t259 + t271, -pkin(4) * t11 - qJ(5) * t10 - t25 * t28 + t26 * t315 - t35 * t45, -t281 * t85 + t333 (-t353 * t85 - t31) * t220 + (-t30 + t359) * t217, t250 + t355, t237 - t356, t353 * t109, qJ(5) * t31 + t3 * t109 + t244 * t217 + t220 * t348 + t314 * t83, qJ(5) * t30 - t4 * t109 - t217 * t348 + t244 * t220 + t314 * t85; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t226, t196 + t325, -t195 ^ 2 - t347, -t195 * t26 + t11 - t337, 0, 0, 0, 0, 0, t195 * t83 + t250, t195 * t85 + t237; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t85 * t83, -t83 ^ 2 + t85 ^ 2, t30 + t359, t349 * t85 - t241, -t64, -t18 * t85 + t353 * t4 + t2, t18 * t83 + t3 * t349 - t267;];
tauc_reg  = t5;
