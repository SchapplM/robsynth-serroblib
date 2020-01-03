% Calculate minimal parameter regressor of inverse dynamics joint torque vector for
% S5RRRRR9
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
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d2,d3,d4,d5]';
% 
% Output:
% tau_reg [5x31]
%   minimal parameter regressor of inverse dynamics joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 22:31
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tau_reg = S5RRRRR9_invdynJ_fixb_regmin_slag_vp(qJ, qJD, qJDD, g, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(5,1),zeros(3,1),zeros(9,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRRRR9_invdynJ_fixb_regmin_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RRRRR9_invdynJ_fixb_regmin_slag_vp: qJD has to be [5x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [5 1]), ...
  'S5RRRRR9_invdynJ_fixb_regmin_slag_vp: qJDD has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RRRRR9_invdynJ_fixb_regmin_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5RRRRR9_invdynJ_fixb_regmin_slag_vp: pkin has to be [9x1] (double)');

%% Symbolic Calculation
% From invdyn_joint_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 22:29:56
% EndTime: 2019-12-31 22:30:10
% DurationCPUTime: 5.93s
% Computational Cost: add. (5735->468), mult. (12969->652), div. (0->0), fcn. (9569->14), ass. (0->233)
t215 = sin(qJ(5));
t220 = cos(qJ(5));
t217 = sin(qJ(3));
t218 = sin(qJ(2));
t302 = qJD(1) * t218;
t280 = t217 * t302;
t222 = cos(qJ(3));
t287 = t222 * qJD(2);
t157 = t280 - t287;
t298 = qJD(2) * t217;
t159 = t222 * t302 + t298;
t216 = sin(qJ(4));
t221 = cos(qJ(4));
t248 = t157 * t216 - t159 * t221;
t93 = t157 * t221 + t159 * t216;
t251 = t215 * t93 + t220 * t248;
t45 = t215 * t248 - t220 * t93;
t356 = t251 * t45;
t351 = t251 ^ 2 - t45 ^ 2;
t223 = cos(qJ(2));
t301 = qJD(1) * t223;
t191 = -qJD(3) + t301;
t183 = -qJD(4) + t191;
t175 = -qJD(5) + t183;
t286 = qJD(1) * qJD(2);
t271 = t223 * t286;
t285 = t218 * qJDD(1);
t294 = qJD(3) * t218;
t353 = -qJD(1) * t294 + qJDD(2);
t81 = qJD(3) * t287 + (t271 + t285) * t222 + t353 * t217;
t82 = (qJD(2) * (qJD(3) + t301) + t285) * t217 - t353 * t222;
t229 = qJD(4) * t248 - t216 * t81 - t221 * t82;
t291 = qJD(4) * t221;
t292 = qJD(4) * t216;
t27 = -t157 * t291 - t159 * t292 - t216 * t82 + t221 * t81;
t289 = qJD(5) * t220;
t290 = qJD(5) * t215;
t4 = t215 * t229 + t220 * t27 + t248 * t290 - t289 * t93;
t350 = t175 * t45 + t4;
t214 = qJ(3) + qJ(4);
t210 = qJ(5) + t214;
t196 = sin(t210);
t197 = cos(t210);
t224 = cos(qJ(1));
t219 = sin(qJ(1));
t314 = t219 * t223;
t115 = t196 * t224 - t197 * t314;
t311 = t223 * t224;
t117 = t196 * t219 + t197 * t311;
t170 = -pkin(2) * t223 - pkin(7) * t218 - pkin(1);
t150 = t170 * qJD(1);
t202 = pkin(6) * t301;
t177 = qJD(2) * pkin(7) + t202;
t103 = t150 * t222 - t177 * t217;
t67 = -pkin(8) * t159 + t103;
t59 = -pkin(3) * t191 + t67;
t104 = t150 * t217 + t177 * t222;
t68 = -pkin(8) * t157 + t104;
t64 = t221 * t68;
t30 = t216 * t59 + t64;
t359 = pkin(9) * t93;
t18 = t30 - t359;
t14 = t18 * t290;
t332 = g(3) * t218;
t176 = -qJD(2) * pkin(2) + pkin(6) * t302;
t112 = pkin(3) * t157 + t176;
t56 = pkin(4) * t93 + t112;
t348 = g(1) * t117 - g(2) * t115 + t197 * t332 - t56 * t45 + t14;
t114 = t196 * t314 + t197 * t224;
t116 = -t196 * t311 + t197 * t219;
t206 = t223 * qJDD(1);
t340 = -t218 * t286 + t206;
t154 = qJDD(3) - t340;
t149 = qJDD(4) + t154;
t132 = pkin(6) * t340 + qJDD(2) * pkin(7);
t255 = pkin(2) * t218 - pkin(7) * t223;
t168 = t255 * qJD(2);
t109 = qJD(1) * t168 + qJDD(1) * t170;
t98 = t222 * t109;
t19 = pkin(3) * t154 - pkin(8) * t81 - qJD(3) * t104 - t132 * t217 + t98;
t293 = qJD(3) * t222;
t295 = qJD(3) * t217;
t240 = t109 * t217 + t132 * t222 + t150 * t293 - t177 * t295;
t23 = -pkin(8) * t82 + t240;
t268 = t19 * t221 - t216 * t23;
t231 = -qJD(4) * t30 + t268;
t2 = pkin(4) * t149 - pkin(9) * t27 + t231;
t263 = -t19 * t216 - t221 * t23 - t291 * t59 + t292 * t68;
t3 = pkin(9) * t229 - t263;
t282 = t2 * t220 - t215 * t3;
t362 = -g(1) * t116 + g(2) * t114 + t196 * t332 + t56 * t251 + t282;
t230 = qJD(5) * t251 - t215 * t27 + t220 * t229;
t344 = t175 * t251 + t230;
t312 = t222 * t223;
t247 = pkin(3) * t218 - pkin(8) * t312;
t335 = pkin(7) + pkin(8);
t281 = qJD(3) * t335;
t165 = t255 * qJD(1);
t305 = pkin(6) * t280 + t165 * t222;
t361 = qJD(1) * t247 + t222 * t281 + t305;
t143 = t217 * t165;
t315 = t218 * t222;
t316 = t217 * t223;
t355 = -t143 - (-pkin(6) * t315 - pkin(8) * t316) * qJD(1) - t217 * t281;
t200 = pkin(6) * t285;
t133 = -qJDD(2) * pkin(2) + pkin(6) * t271 + t200;
t254 = g(1) * t224 + g(2) * t219;
t331 = g(3) * t223;
t235 = t218 * t254 - t331;
t360 = qJD(3) * pkin(7) * t191 - t133 + t235;
t358 = pkin(9) * t248;
t357 = t248 * t93;
t160 = t216 * t217 - t221 * t222;
t239 = t160 * t223;
t337 = qJD(3) + qJD(4);
t310 = qJD(1) * t239 - t160 * t337;
t161 = t216 * t222 + t217 * t221;
t309 = (-t301 + t337) * t161;
t296 = qJD(2) * t223;
t274 = t217 * t296;
t354 = t218 * t293 + t274;
t352 = t248 ^ 2 - t93 ^ 2;
t349 = -t183 * t93 + t27;
t208 = sin(t214);
t209 = cos(t214);
t125 = t208 * t224 - t209 * t314;
t127 = t208 * t219 + t209 * t311;
t347 = g(1) * t127 - g(2) * t125 + t112 * t93 + t209 * t332 + t263;
t62 = t216 * t68;
t29 = t221 * t59 - t62;
t17 = t29 + t358;
t13 = -pkin(4) * t183 + t17;
t325 = t18 * t220;
t7 = t13 * t215 + t325;
t346 = -qJD(5) * t7 + t362;
t124 = t208 * t314 + t209 * t224;
t126 = -t208 * t311 + t209 * t219;
t345 = -g(1) * t126 + g(2) * t124 + t112 * t248 + t208 * t332 + t231;
t343 = t183 * t248 + t229;
t342 = t361 * t221;
t129 = t161 * t218;
t256 = -t202 + (-t217 * t301 + t295) * pkin(3);
t178 = t335 * t217;
t179 = t335 * t222;
t306 = -t178 * t216 + t179 * t221;
t341 = -t178 * t291 - t179 * t292 - t216 * t361 + t221 * t355;
t339 = -t217 * t294 + t223 * t287;
t334 = pkin(3) * t216;
t333 = pkin(6) * t217;
t99 = t160 * t220 + t161 * t215;
t330 = -qJD(5) * t99 - t215 * t309 + t220 * t310;
t100 = -t160 * t215 + t161 * t220;
t329 = qJD(5) * t100 + t215 * t310 + t220 * t309;
t328 = t221 * t67 - t62;
t327 = pkin(4) * t309 + t256;
t324 = t217 * t81;
t156 = t222 * t170;
t102 = -pkin(8) * t315 + t156 + (-pkin(3) - t333) * t223;
t193 = pkin(6) * t312;
t304 = t170 * t217 + t193;
t317 = t217 * t218;
t110 = -pkin(8) * t317 + t304;
t323 = t102 * t216 + t110 * t221;
t322 = t157 * t191;
t321 = t159 * t191;
t320 = t159 * t222;
t134 = qJDD(5) + t149;
t319 = t215 * t134;
t318 = t216 * t220;
t313 = t220 * t134;
t308 = t168 * t217 + t170 * t293;
t297 = qJD(2) * t218;
t307 = t168 * t222 + t297 * t333;
t169 = pkin(3) * t317 + pkin(6) * t218;
t212 = t218 ^ 2;
t303 = -t223 ^ 2 + t212;
t300 = qJD(2) * t157;
t299 = qJD(2) * t159;
t288 = t176 * qJD(3);
t113 = pkin(3) * t354 + pkin(6) * t296;
t199 = -pkin(3) * t222 - pkin(2);
t279 = t191 * t287;
t278 = t191 * t295;
t277 = t191 * t293;
t269 = qJD(5) * t13 + t3;
t49 = t247 * qJD(2) + (-t193 + (pkin(8) * t218 - t170) * t217) * qJD(3) + t307;
t52 = -t354 * pkin(8) + (-t218 * t287 - t223 * t295) * pkin(6) + t308;
t266 = -t216 * t52 + t221 * t49;
t265 = -t216 * t67 - t64;
t262 = t102 * t221 - t110 * t216;
t261 = -t178 * t221 - t179 * t216;
t260 = -qJD(3) * t150 - t132;
t77 = -pkin(9) * t160 + t306;
t258 = pkin(4) * t302 + pkin(9) * t310 + t306 * qJD(4) + qJD(5) * t77 + t216 * t355 + t342;
t76 = -pkin(9) * t161 + t261;
t257 = -pkin(9) * t309 + qJD(5) * t76 + t341;
t253 = g(1) * t219 - g(2) * t224;
t252 = t177 * t293 - t98;
t250 = -pkin(7) * t154 + t288;
t130 = t160 * t218;
t69 = t129 * t220 - t130 * t215;
t70 = -t129 * t215 - t130 * t220;
t245 = -0.2e1 * pkin(1) * t286 - pkin(6) * qJDD(2);
t244 = t154 * t217 - t277;
t243 = t222 * t154 + t278;
t242 = t102 * t291 - t110 * t292 + t216 * t49 + t221 * t52;
t58 = pkin(3) * t82 + t133;
t226 = qJD(1) ^ 2;
t237 = pkin(1) * t226 + t254;
t225 = qJD(2) ^ 2;
t232 = 0.2e1 * qJDD(1) * pkin(1) - pkin(6) * t225 + t253;
t198 = pkin(3) * t221 + pkin(4);
t141 = t217 * t219 + t222 * t311;
t140 = -t217 * t311 + t219 * t222;
t139 = t217 * t224 - t219 * t312;
t138 = t217 * t314 + t222 * t224;
t121 = pkin(4) * t160 + t199;
t105 = pkin(4) * t129 + t169;
t61 = pkin(3) * t159 - pkin(4) * t248;
t54 = -t292 * t317 + (t315 * t337 + t274) * t221 + t339 * t216;
t53 = -qJD(2) * t239 - t129 * t337;
t39 = pkin(4) * t54 + t113;
t38 = -pkin(9) * t129 + t323;
t37 = -pkin(4) * t223 + pkin(9) * t130 + t262;
t21 = t328 + t358;
t20 = t265 + t359;
t12 = -pkin(4) * t229 + t58;
t11 = qJD(5) * t70 + t215 * t53 + t220 * t54;
t10 = -qJD(5) * t69 - t215 * t54 + t220 * t53;
t9 = -pkin(9) * t54 + t242;
t8 = pkin(4) * t297 - pkin(9) * t53 - qJD(4) * t323 + t266;
t6 = t13 * t220 - t18 * t215;
t1 = [qJDD(1), t253, t254, qJDD(1) * t212 + 0.2e1 * t218 * t271, 0.2e1 * t206 * t218 - 0.2e1 * t286 * t303, qJDD(2) * t218 + t223 * t225, qJDD(2) * t223 - t218 * t225, 0, t218 * t245 + t223 * t232, -t218 * t232 + t223 * t245, t159 * t339 + t315 * t81, (-t157 * t222 - t159 * t217) * t296 + (-t324 - t222 * t82 + (t157 * t217 - t320) * qJD(3)) * t218, (-t81 - t279) * t223 + (t243 + t299) * t218, (t191 * t298 + t82) * t223 + (-t244 - t300) * t218, -t154 * t223 - t191 * t297, -(-t170 * t295 + t307) * t191 + t156 * t154 - g(1) * t139 - g(2) * t141 + ((t277 + t300) * pkin(6) + (-pkin(6) * t154 + qJD(2) * t176 - t260) * t217 + t252) * t223 + (pkin(6) * t82 + qJD(2) * t103 + t133 * t217 + t222 * t288) * t218, t308 * t191 - t304 * t154 - g(1) * t138 - g(2) * t140 + (t176 * t287 + (-t278 + t299) * pkin(6) + t240) * t223 + (-t217 * t288 - t104 * qJD(2) + t133 * t222 + (t81 - t279) * pkin(6)) * t218, -t130 * t27 - t248 * t53, -t129 * t27 - t130 * t229 + t248 * t54 - t53 * t93, -t130 * t149 - t183 * t53 - t223 * t27 - t248 * t297, -t129 * t149 + t183 * t54 - t223 * t229 - t297 * t93, -t149 * t223 - t183 * t297, -t266 * t183 + t262 * t149 - t268 * t223 + t29 * t297 + t113 * t93 - t169 * t229 + t58 * t129 + t112 * t54 - g(1) * t125 - g(2) * t127 + (t183 * t323 + t223 * t30) * qJD(4), -g(1) * t124 - g(2) * t126 + t112 * t53 - t113 * t248 - t130 * t58 - t149 * t323 + t169 * t27 + t183 * t242 - t223 * t263 - t297 * t30, -t10 * t251 + t4 * t70, t10 * t45 + t11 * t251 + t230 * t70 - t4 * t69, -t10 * t175 + t134 * t70 - t223 * t4 - t251 * t297, t11 * t175 - t134 * t69 - t223 * t230 + t297 * t45, -t134 * t223 - t175 * t297, -(-t215 * t9 + t220 * t8) * t175 + (-t215 * t38 + t220 * t37) * t134 - t282 * t223 + t6 * t297 - t39 * t45 - t105 * t230 + t12 * t69 + t56 * t11 - g(1) * t115 - g(2) * t117 + (-(-t215 * t37 - t220 * t38) * t175 + t7 * t223) * qJD(5), -t7 * t297 - g(1) * t114 - g(2) * t116 + t56 * t10 + t105 * t4 + t12 * t70 - t14 * t223 - t39 * t251 + ((-qJD(5) * t38 + t8) * t175 - t37 * t134 + t2 * t223) * t215 + ((qJD(5) * t37 + t9) * t175 - t38 * t134 + t269 * t223) * t220; 0, 0, 0, -t218 * t226 * t223, t303 * t226, t285, t206, qJDD(2), t218 * t237 - t200 - t331, t332 + (-pkin(6) * qJDD(1) + t237) * t223, -t191 * t320 + t324, (t81 + t322) * t222 + (-t82 + t321) * t217, (-t159 * t218 + t191 * t312) * qJD(1) + t244, (t157 * t218 - t191 * t316) * qJD(1) + t243, t191 * t302, -pkin(2) * t82 + t305 * t191 + t250 * t217 + (-t103 * t218 + (-pkin(6) * t157 - t176 * t217) * t223) * qJD(1) + t360 * t222, -pkin(2) * t81 - t143 * t191 + t250 * t222 + (-t176 * t312 + t104 * t218 + (-t159 * t223 + t191 * t315) * pkin(6)) * qJD(1) - t360 * t217, t161 * t27 - t248 * t310, -t160 * t27 + t161 * t229 + t248 * t309 - t310 * t93, t149 * t161 - t183 * t310 + t248 * t302, -t149 * t160 + t183 * t309 + t302 * t93, t183 * t302, t261 * t149 - t199 * t229 + t58 * t160 - t29 * t302 + t256 * t93 + (t179 * t291 + (-qJD(4) * t178 + t355) * t216 + t342) * t183 + t309 * t112 + t235 * t209, t112 * t310 - t149 * t306 + t58 * t161 + t183 * t341 + t199 * t27 - t208 * t235 - t248 * t256 + t30 * t302, t100 * t4 - t251 * t330, t100 * t230 + t251 * t329 + t330 * t45 - t4 * t99, t100 * t134 - t175 * t330 + t251 * t302, -t134 * t99 + t175 * t329 - t302 * t45, t175 * t302, (-t215 * t77 + t220 * t76) * t134 - t121 * t230 + t12 * t99 - t6 * t302 + t329 * t56 - t327 * t45 + (t215 * t257 + t220 * t258) * t175 + t235 * t197, -(t215 * t76 + t220 * t77) * t134 + t121 * t4 + t12 * t100 + t7 * t302 + t330 * t56 - t327 * t251 + (-t215 * t258 + t220 * t257) * t175 - t235 * t196; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t159 * t157, -t157 ^ 2 + t159 ^ 2, t81 - t322, -t321 - t82, t154, -g(1) * t140 + g(2) * t138 - t104 * t191 - t159 * t176 + (t260 + t332) * t217 - t252, g(1) * t141 - g(2) * t139 + g(3) * t315 - t103 * t191 + t157 * t176 - t240, -t357, t352, t349, t343, t149, t265 * t183 + (t149 * t221 - t159 * t93 + t183 * t292) * pkin(3) + t345, -t328 * t183 + (-t149 * t216 + t159 * t248 + t183 * t291) * pkin(3) + t347, t356, t351, t350, t344, t134, t198 * t313 + (t20 * t220 - t21 * t215) * t175 + t61 * t45 + (-t216 * t319 - (-t215 * t221 - t318) * t175 * qJD(4)) * pkin(3) + (-(-pkin(3) * t318 - t198 * t215) * t175 - t7) * qJD(5) + t362, t61 * t251 + (-t198 * t134 - t2 + (-t20 + (-qJD(4) - qJD(5)) * t334) * t175) * t215 + (-t134 * t334 + (pkin(3) * t291 + qJD(5) * t198 - t21) * t175 - t269) * t220 + t348; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t357, t352, t349, t343, t149, -t183 * t30 + t345, -t183 * t29 + t347, t356, t351, t350, t344, t134, (-t17 * t215 - t325) * t175 + (t175 * t290 - t248 * t45 + t313) * pkin(4) + t346, (t175 * t18 - t2) * t215 + (-t17 * t175 - t269) * t220 + (t175 * t289 - t248 * t251 - t319) * pkin(4) + t348; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t356, t351, t350, t344, t134, -t175 * t7 + t346, -t175 * t6 - t215 * t2 - t220 * t269 + t348;];
tau_reg = t1;
