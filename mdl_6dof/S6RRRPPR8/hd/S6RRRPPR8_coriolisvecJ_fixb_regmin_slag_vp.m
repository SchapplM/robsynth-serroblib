% Calculate minimal parameter regressor of coriolis joint torque vector for
% S6RRRPPR8
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d1,d2,d3,d6]';
% 
% Output:
% tauc_reg [6x32]
%   minimal parameter regressor of coriolis joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 16:10
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauc_reg = S6RRRPPR8_coriolisvecJ_fixb_regmin_slag_vp(qJ, qJD, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(10,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRPPR8_coriolisvecJ_fixb_regmin_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRRPPR8_coriolisvecJ_fixb_regmin_slag_vp: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RRRPPR8_coriolisvecJ_fixb_regmin_slag_vp: pkin has to be [10x1] (double)');

%% Symbolic Calculation
% From coriolisvec_joint_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 16:09:30
% EndTime: 2019-03-09 16:09:43
% DurationCPUTime: 4.79s
% Computational Cost: add. (5163->499), mult. (13702->653), div. (0->0), fcn. (10153->8), ass. (0->230)
t192 = sin(qJ(3));
t195 = cos(qJ(3));
t305 = cos(pkin(6));
t251 = t305 * qJD(1);
t227 = t251 + qJD(2);
t193 = sin(qJ(2));
t189 = sin(pkin(6));
t287 = qJD(1) * t189;
t262 = t193 * t287;
t118 = t192 * t227 + t195 * t262;
t196 = cos(qJ(2));
t297 = t189 * t196;
t209 = t305 * t193 * pkin(1) + pkin(8) * t297;
t330 = t209 * qJD(2);
t128 = qJD(1) * t330;
t213 = t195 * t227;
t273 = qJD(1) * qJD(2);
t256 = t189 * t273;
t239 = t196 * t256;
t298 = t189 * t193;
t264 = t192 * t298;
t241 = qJD(3) * t264;
t72 = qJD(1) * t241 - qJD(3) * t213 - t195 * t239;
t282 = qJD(3) * t118;
t73 = t192 * t239 + t282;
t20 = t73 * pkin(3) + t72 * qJ(4) - t118 * qJD(4) + t128;
t286 = qJD(1) * t196;
t261 = t189 * t286;
t242 = t192 * t261;
t281 = qJD(3) * t192;
t342 = t242 - t281;
t145 = t305 * t192 + t195 * t298;
t285 = qJD(2) * t189;
t258 = t196 * t285;
t81 = qJD(3) * t145 + t192 * t258;
t82 = -t241 + (t305 * qJD(3) + t258) * t195;
t30 = t81 * pkin(3) - t82 * qJ(4) - t145 * qJD(4) + t330;
t159 = -qJD(3) + t261;
t131 = (-pkin(2) * t196 - pkin(9) * t193 - pkin(1)) * t189;
t100 = qJD(1) * t131;
t240 = pkin(1) * t251;
t136 = pkin(8) * t261 + t193 * t240;
t90 = pkin(9) * t227 + t136;
t324 = -t195 * t100 + t192 * t90;
t341 = qJD(4) + t324;
t325 = pkin(9) - qJ(5);
t133 = -pkin(8) * t262 + t196 * t240;
t220 = t189 * (pkin(2) * t193 - pkin(9) * t196);
t134 = qJD(1) * t220;
t290 = t195 * t133 + t192 * t134;
t59 = qJ(4) * t262 + t290;
t340 = -qJD(5) * t195 - t281 * t325 - t59;
t163 = t193 * t256;
t230 = pkin(3) * t163;
t135 = qJD(2) * t220;
t126 = qJD(1) * t135;
t255 = t196 * t305;
t332 = pkin(1) * t255 - pkin(8) * t298;
t137 = t332 * qJD(2);
t127 = qJD(1) * t137;
t279 = qJD(3) * t195;
t253 = t100 * t281 - t195 * t126 + t192 * t127 + t90 * t279;
t21 = -t230 + t253;
t152 = t159 * qJ(4);
t54 = t192 * t100 + t195 * t90;
t50 = -t152 + t54;
t339 = t159 * t50 + t21;
t116 = t192 * t262 - t213;
t303 = t116 * t159;
t338 = -t72 + t303;
t301 = t118 * t159;
t337 = t73 - t301;
t336 = t193 * t196;
t335 = -t118 * qJ(5) + t341;
t115 = t118 ^ 2;
t334 = -t159 ^ 2 - t115;
t247 = -t192 * t133 + t134 * t195;
t333 = -qJD(5) * t192 + t325 * t279 + t247;
t331 = t192 * qJD(4) + t136 + (-t195 * t261 + t279) * qJ(4);
t109 = qJD(6) + t118;
t326 = -pkin(4) - pkin(10);
t274 = pkin(3) - t326;
t62 = t118 * pkin(3) + t116 * qJ(4);
t157 = qJ(4) * t163;
t150 = qJD(4) * t159;
t254 = -t100 * t279 - t192 * t126 - t195 * t127 + t90 * t281;
t238 = -t150 - t254;
t217 = -t73 * qJ(5) - t116 * qJD(5) - t238;
t9 = -t157 + t217;
t7 = pkin(5) * t163 - t9;
t329 = t109 * (-pkin(5) * t116 - qJD(6) * t274 + t326 * t118 - t62) - t7;
t191 = sin(qJ(6));
t194 = cos(qJ(6));
t68 = t116 * t194 + t159 * t191;
t33 = qJD(6) * t68 + t163 * t194 + t191 * t73;
t311 = t194 * t72;
t210 = t109 ^ 2 * t191 + t311;
t198 = qJD(1) ^ 2;
t327 = pkin(3) + pkin(4);
t190 = qJ(4) + pkin(5);
t323 = qJ(4) * t73;
t66 = t116 * t191 - t194 * t159;
t322 = t109 * t66;
t321 = t109 * t68;
t320 = t116 * t66;
t319 = t116 * t68;
t89 = -t227 * pkin(2) - t133;
t43 = t116 * pkin(3) - t118 * qJ(4) + t89;
t318 = t118 * t43;
t36 = qJ(5) * t116 + t54;
t31 = t152 - t36;
t317 = t159 * t31;
t315 = t159 * t66;
t314 = t159 * t68;
t275 = qJD(6) * t194;
t32 = t159 * t275 + t194 * t73 + (-qJD(6) * t116 - t163) * t191;
t313 = t191 * t32;
t312 = t191 * t72;
t310 = t327 * t342 + t331;
t296 = t192 * t196;
t309 = -(pkin(5) * t193 + qJ(5) * t296) * t287 + t340;
t308 = qJ(5) * t242 - t340;
t307 = -pkin(3) * t342 - t331;
t294 = t195 * t196;
t265 = qJ(5) * t294;
t306 = -(-t327 * t193 - t265) * t287 + t333;
t144 = -t305 * t195 + t264;
t304 = qJ(5) * t144;
t302 = t118 * t116;
t300 = t159 * t195;
t185 = t189 ^ 2;
t299 = t185 * t198;
t295 = t194 * t196;
t130 = t305 * pkin(9) + t209;
t291 = t195 * t130 + t192 * t131;
t288 = t193 ^ 2 - t196 ^ 2;
t284 = qJD(2) * t193;
t283 = qJD(2) * t195;
t280 = qJD(3) * t194;
t278 = qJD(4) * t196;
t277 = qJD(5) * t118;
t276 = qJD(6) * t191;
t272 = pkin(9) * t159 * t192;
t271 = pkin(9) * t300;
t156 = -t195 * pkin(3) - t192 * qJ(4) - pkin(2);
t268 = pkin(9) * t284;
t267 = pkin(9) * t283;
t129 = -t305 * pkin(2) - t332;
t260 = t189 * t284;
t259 = t189 * t278;
t257 = t185 * t273;
t219 = qJD(5) - t43;
t29 = -pkin(4) * t116 + t219;
t252 = (-qJD(5) - t29) * t118;
t249 = -t192 * t130 + t131 * t195;
t248 = -t163 + t302;
t246 = t109 * t194;
t245 = t274 * t193;
t146 = t195 * pkin(4) - t156;
t244 = 0.2e1 * t257;
t161 = t325 * t192;
t237 = qJD(6) * t161 - t331 + t159 * (pkin(5) * t195 - t274 * t192);
t236 = -0.2e1 * pkin(1) * t257;
t95 = (-t191 * t193 + t192 * t295) * t287;
t235 = t192 * t280 - t95;
t6 = -t72 * pkin(5) + t326 * t73 - t20;
t223 = qJ(5) * t72 + t253;
t224 = t245 * t285;
t8 = -qJD(1) * t224 + t223 - t277;
t234 = t191 * t6 + t194 * t8;
t108 = pkin(5) * t192 + pkin(10) * t195 + t146;
t233 = -qJD(6) * t108 + (-t245 - t265) * t287 - t333;
t57 = pkin(3) * t297 - t249;
t17 = pkin(5) * t118 + t326 * t116 + t219;
t22 = t274 * t159 + t335;
t3 = t17 * t194 - t191 * t22;
t4 = t17 * t191 + t194 * t22;
t55 = t144 * pkin(3) - t145 * qJ(4) + t129;
t25 = pkin(5) * t145 + t326 * t144 - t55;
t39 = pkin(4) * t297 - qJ(5) * t145 + t57;
t34 = pkin(10) * t297 + t39;
t232 = -t191 * t34 + t194 * t25;
t231 = t191 * t25 + t194 * t34;
t229 = t327 * t260;
t226 = 0.2e1 * t251 + qJD(2);
t225 = -t130 * t279 - t131 * t281 + t135 * t195 - t192 * t137;
t56 = -qJ(4) * t297 + t291;
t222 = t159 * t324 + t254;
t221 = -t159 * t54 - t253;
t218 = -t144 * t191 + t189 * t295;
t84 = t144 * t194 + t191 * t297;
t214 = -t130 * t281 + t131 * t279 + t192 * t135 + t195 * t137;
t212 = qJ(4) * t260 + t214;
t211 = -t109 * t246 + t312;
t27 = -pkin(5) * t159 - t31;
t207 = -t274 * t72 + (-t27 + t36) * t109;
t206 = pkin(1) * (-qJD(2) * t251 + t299);
t2 = -qJD(6) * t4 - t191 * t8 + t194 * t6;
t204 = t72 + t303;
t203 = -qJ(5) * t82 - qJD(5) * t145 - t225;
t202 = -qJD(1) * t229 + t223;
t201 = qJ(5) * t81 + qJD(5) * t144 + t212;
t15 = -t73 * pkin(4) - t20;
t162 = t325 * t195;
t154 = 0.2e1 * t157;
t114 = t116 ^ 2;
t94 = (t191 * t296 + t193 * t194) * t287;
t65 = t72 * t192;
t61 = -pkin(3) * t262 - t247;
t58 = t72 * t145;
t49 = -pkin(4) * t118 - t62;
t48 = pkin(3) * t159 + t341;
t46 = -t56 - t304;
t44 = -pkin(4) * t144 - t55;
t41 = qJD(6) * t218 - t191 * t260 + t194 * t81;
t40 = qJD(6) * t84 + t191 * t81 + t194 * t260;
t38 = -t190 * t297 + t291 + t304;
t28 = -pkin(3) * t260 - t225;
t26 = t327 * t159 + t335;
t24 = t212 - t259;
t19 = -t81 * pkin(4) - t30;
t18 = t157 + t238;
t16 = -t201 + t259;
t14 = -t229 + t203;
t13 = (pkin(5) * t284 - t278) * t189 + t201;
t12 = -t224 + t203;
t11 = t82 * pkin(5) + t326 * t81 - t30;
t10 = t202 - t277;
t1 = t3 * qJD(6) + t234;
t5 = [0, 0, 0, t244 * t336, -t288 * t244, t226 * t258, -t226 * t260, 0, -t128 * t305 + t193 * t236 - t227 * t330, -t127 * t305 - t137 * t227 + t196 * t236, t118 * t82 - t58, -t116 * t82 - t118 * t81 + t144 * t72 - t145 * t73, -t159 * t82 + (t196 * t72 + (qJD(1) * t145 + t118) * t284) * t189, t159 * t81 + (t196 * t73 + (-qJD(1) * t144 - t116) * t284) * t189 (-t159 * t189 - t185 * t286) * t284, -t225 * t159 + t330 * t116 + t129 * t73 + t128 * t144 + t89 * t81 + (t253 * t196 + (qJD(1) * t249 - t324) * t284) * t189, t214 * t159 + t330 * t118 - t129 * t72 + t128 * t145 + t89 * t82 + (-t254 * t196 + (-t291 * qJD(1) - t54) * t284) * t189, t116 * t30 + t144 * t20 + t159 * t28 + t43 * t81 + t55 * t73 + (t196 * t21 + (-qJD(1) * t57 - t48) * t284) * t189, -t116 * t24 + t118 * t28 - t144 * t18 + t145 * t21 + t48 * t82 - t50 * t81 - t56 * t73 - t57 * t72, -t118 * t30 - t145 * t20 - t159 * t24 - t43 * t82 + t55 * t72 + (-t18 * t196 + (qJD(1) * t56 + t50) * t284) * t189, t18 * t56 + t20 * t55 + t21 * t57 + t24 * t50 + t28 * t48 + t30 * t43, t118 * t19 + t145 * t15 + t159 * t16 + t29 * t82 - t44 * t72 + (t196 * t9 + (-qJD(1) * t46 - t31) * t284) * t189, t116 * t19 - t14 * t159 + t144 * t15 + t29 * t81 + t44 * t73 + (-t10 * t196 + (qJD(1) * t39 + t26) * t284) * t189, -t10 * t145 - t116 * t16 - t118 * t14 - t144 * t9 - t26 * t82 - t31 * t81 + t39 * t72 - t46 * t73, t10 * t39 + t14 * t26 + t15 * t44 + t16 * t31 + t19 * t29 + t46 * t9, t32 * t84 + t41 * t68, t218 * t32 - t33 * t84 - t40 * t68 - t41 * t66, t109 * t41 + t145 * t32 + t68 * t82 - t72 * t84, -t109 * t40 - t145 * t33 - t218 * t72 - t66 * t82, t109 * t82 - t58 (-qJD(6) * t231 + t11 * t194 - t12 * t191) * t109 - t232 * t72 + t2 * t145 + t3 * t82 + t13 * t66 + t38 * t33 - t7 * t218 + t27 * t40 -(qJD(6) * t232 + t11 * t191 + t12 * t194) * t109 + t231 * t72 - t1 * t145 - t4 * t82 + t13 * t68 + t38 * t32 + t7 * t84 + t27 * t41; 0, 0, 0, -t299 * t336, t288 * t299, -t189 * t198 * t255, t227 * t262 - t163, 0, -pkin(8) * t239 + t136 * t227 + t193 * t206, pkin(8) * t163 + t133 * t227 + t196 * t206, -t118 * t300 - t65, -t192 * t337 + t195 * t338, -t159 * t279 + (t159 * t294 + (t192 * qJD(2) - t118) * t193) * t287, t159 * t281 + (-t159 * t296 + (t116 + t283) * t193) * t287, t159 * t262, -pkin(2) * t73 - t128 * t195 + t247 * t159 - t136 * t116 + (t192 * t89 + t271) * qJD(3) + (t324 * t193 + (-t196 * t89 - t268) * t192) * t287, pkin(2) * t72 + t128 * t192 - t290 * t159 - t136 * t118 + (t195 * t89 - t272) * qJD(3) + (-t89 * t294 + (t54 - t267) * t193) * t287, t156 * t73 - t159 * t61 - t195 * t20 + t307 * t116 + (t192 * t43 + t271) * qJD(3) + (t193 * t48 + (-t196 * t43 - t268) * t192) * t287, t116 * t59 - t118 * t61 + (t18 - t159 * t48 + (-t73 + t282) * pkin(9)) * t195 + ((qJD(3) * t116 - t72) * pkin(9) + t339) * t192, t156 * t72 + t159 * t59 - t192 * t20 - t307 * t118 + (-t195 * t43 + t272) * qJD(3) + (t43 * t294 + (-t50 + t267) * t193) * t287, t156 * t20 - t48 * t61 - t50 * t59 + t307 * t43 + (t18 * t195 + t192 * t21 + (-t192 * t50 + t195 * t48) * qJD(3)) * pkin(9), t29 * t279 - t146 * t72 + t15 * t192 + t308 * t159 + t310 * t118 + (-t29 * t294 + (qJD(2) * t162 + t31) * t193) * t287, t29 * t281 + t146 * t73 - t15 * t195 - t306 * t159 + t310 * t116 + (-t29 * t296 + (qJD(2) * t161 - t26) * t193) * t287, t161 * t72 + t162 * t73 - t306 * t118 - t308 * t116 + (t159 * t26 + t9) * t195 + (-t10 + t317) * t192, t10 * t161 + t146 * t15 - t162 * t9 + t306 * t26 + t310 * t29 + t308 * t31, -t194 * t195 * t32 + (t195 * t276 + t235) * t68, t66 * t95 + t68 * t94 + (-t191 * t68 - t194 * t66) * t281 + (t313 + t194 * t33 + (-t191 * t66 + t194 * t68) * qJD(6)) * t195, t192 * t32 + t235 * t109 + (t109 * t276 + t311 - t314) * t195, -t192 * t33 + (-t191 * t281 + t94) * t109 + (t109 * t275 - t312 + t315) * t195, -t109 * t300 - t65 -(t108 * t194 - t161 * t191) * t72 + t162 * t33 - t27 * t94 + t309 * t66 + (t27 * t191 * qJD(3) + t2) * t192 + (t191 * t233 - t194 * t237) * t109 + (-t159 * t3 - t7 * t191 - t27 * t275) * t195 (t108 * t191 + t161 * t194) * t72 + t162 * t32 - t27 * t95 + t309 * t68 + (t27 * t280 - t1) * t192 + (t191 * t237 + t194 * t233) * t109 + (t159 * t4 - t7 * t194 + t27 * t276) * t195; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t302, -t114 + t115, -t204, -t301 - t73, t163, -t118 * t89 + t221, t116 * t89 + t222, -t116 * t62 + t221 + 0.2e1 * t230 - t318, pkin(3) * t72 - t323 + (t50 - t54) * t118 + (t48 - t341) * t116, -t116 * t43 + t118 * t62 - 0.2e1 * t150 + t154 - t222, -pkin(3) * t21 + qJ(4) * t18 + t341 * t50 - t43 * t62 - t48 * t54, t116 * t29 - t118 * t49 - t159 * t335 + t154 - t217, -t116 * t49 + t159 * t36 - 0.2e1 * t163 * t327 + t223 + t252, t323 - t327 * t72 + (t31 + t36) * t118 + (-t26 + t335) * t116, -qJ(4) * t9 - t10 * t327 - t26 * t36 - t29 * t49 - t31 * t335, -t246 * t68 - t313 (-t32 + t322) * t194 + (t33 + t321) * t191, t211 + t319, t210 - t320, t109 * t116, t3 * t116 + t190 * t33 + t207 * t191 - t194 * t329 + t335 * t66, -t4 * t116 + t190 * t32 + t191 * t329 + t207 * t194 + t335 * t68; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t248, -t204, t334, t318 + t339, t334, -t248, t204, t202 + t252 - t317, 0, 0, 0, 0, 0, t211 + t315, t210 + t314; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t338, t337, -t115 - t114, t31 * t116 + t26 * t118 + t15, 0, 0, 0, 0, 0, -t210 - t320, t211 - t319; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t68 * t66, -t66 ^ 2 + t68 ^ 2, t32 + t322, t321 - t33, -t72, t109 * t4 - t27 * t68 + t2, t27 * t66 - t234 + (-qJD(6) + t109) * t3;];
tauc_reg  = t5;
