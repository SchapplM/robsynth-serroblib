% Calculate minimal parameter regressor of inverse dynamics joint torque vector for
% S6RPRRRP5
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% qJDD [6x1]
%   Generalized joint accelerations
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d3,d4,d5,theta2]';
% 
% Output:
% tau_reg [6x32]
%   minimal parameter regressor of inverse dynamics joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 06:13
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tau_reg = S6RPRRRP5_invdynJ_fixb_regmin_slag_vp(qJ, qJD, qJDD, g, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(3,1),zeros(10,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRRRP5_invdynJ_fixb_regmin_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RPRRRP5_invdynJ_fixb_regmin_slag_vp: qJD has to be [6x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [6 1]), ...
  'S6RPRRRP5_invdynJ_fixb_regmin_slag_vp: qJDD has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RPRRRP5_invdynJ_fixb_regmin_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RPRRRP5_invdynJ_fixb_regmin_slag_vp: pkin has to be [10x1] (double)');

%% Symbolic Calculation
% From invdyn_joint_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 06:12:29
% EndTime: 2019-03-09 06:12:43
% DurationCPUTime: 6.32s
% Computational Cost: add. (11621->507), mult. (28242->607), div. (0->0), fcn. (22473->14), ass. (0->261)
t200 = cos(qJ(5));
t289 = qJD(5) * t200;
t195 = cos(pkin(10));
t201 = cos(qJ(3));
t298 = t201 * t195;
t194 = sin(pkin(10));
t198 = sin(qJ(3));
t304 = t194 * t198;
t137 = -t298 + t304;
t128 = t137 * qJD(1);
t138 = t194 * t201 + t195 * t198;
t129 = t138 * qJD(1);
t197 = sin(qJ(4));
t347 = cos(qJ(4));
t222 = -t347 * t128 - t197 * t129;
t382 = t222 * t200;
t393 = t289 - t382;
t285 = qJD(5) - t222;
t253 = t285 ^ 2;
t196 = sin(qJ(5));
t280 = qJD(3) + qJD(4);
t337 = pkin(7) + qJ(2);
t159 = t337 * t194;
t139 = qJD(1) * t159;
t160 = t337 * t195;
t140 = qJD(1) * t160;
t235 = t139 * t198 - t140 * t201;
t93 = -pkin(8) * t128 - t235;
t90 = t347 * t93;
t363 = -t201 * t139 - t140 * t198;
t92 = -pkin(8) * t129 + t363;
t91 = qJD(3) * pkin(3) + t92;
t63 = t197 * t91 + t90;
t57 = pkin(9) * t280 + t63;
t269 = t195 * pkin(2) + pkin(1);
t147 = -t269 * qJD(1) + qJD(2);
t112 = pkin(3) * t128 + t147;
t221 = -t197 * t128 + t347 * t129;
t64 = -pkin(4) * t222 - pkin(9) * t221 + t112;
t26 = t196 * t64 + t200 * t57;
t19 = qJ(6) * t285 + t26;
t392 = t19 * t285;
t89 = t197 * t93;
t62 = t347 * t91 - t89;
t56 = -t280 * pkin(4) - t62;
t94 = t196 * t221 - t200 * t280;
t96 = t196 * t280 + t200 * t221;
t33 = t94 * pkin(5) - t96 * qJ(6) + t56;
t391 = t285 * t33;
t199 = sin(qJ(1));
t202 = cos(qJ(1));
t244 = g(1) * t202 + g(2) * t199;
t283 = qJD(1) * qJD(3);
t262 = t198 * t283;
t261 = t201 * t283;
t281 = t195 * qJDD(1);
t282 = t194 * qJDD(1);
t270 = t195 * t261 + t198 * t281 + t201 * t282;
t109 = -t194 * t262 + t270;
t271 = -t194 * t261 - t195 * t262 - t198 * t282;
t219 = t201 * t281 + t271;
t256 = t197 * t109 - t347 * t219;
t358 = qJD(4) * t221;
t55 = t256 + t358;
t53 = qJDD(5) + t55;
t45 = t196 * t53;
t233 = t285 * t393 + t45;
t326 = t221 * t96;
t390 = t233 - t326;
t210 = t347 * t109 + t197 * t219;
t287 = t222 * qJD(3);
t389 = -t287 + t210;
t388 = qJD(4) * t222;
t279 = qJDD(3) + qJDD(4);
t290 = qJD(5) * t196;
t206 = t210 + t388;
t379 = qJD(5) * t280 + t206;
t39 = -t196 * t279 - t200 * t379 + t221 * t290;
t37 = t39 * t196;
t387 = t393 * t96 - t37;
t384 = t222 * t33;
t383 = t56 * t222;
t381 = t221 * t222;
t193 = pkin(10) + qJ(3);
t188 = qJ(4) + t193;
t177 = sin(t188);
t307 = t177 * t202;
t308 = t177 * t199;
t380 = g(1) * t307 + g(2) * t308;
t313 = qJDD(1) * pkin(1);
t185 = qJDD(2) - t313;
t361 = g(1) * t199 - g(2) * t202;
t228 = -t185 + t361;
t241 = t200 * pkin(5) + t196 * qJ(6);
t178 = cos(t188);
t378 = t244 * t178;
t377 = t221 ^ 2 - t222 ^ 2;
t78 = pkin(4) * t221 - pkin(9) * t222;
t263 = qJD(4) * t347;
t291 = qJD(4) * t197;
t284 = qJD(1) * qJD(2);
t351 = t337 * qJDD(1) + t284;
t116 = t351 * t194;
t117 = t351 * t195;
t255 = -t201 * t116 - t198 * t117;
t51 = qJDD(3) * pkin(3) - pkin(8) * t109 + qJD(3) * t235 + t255;
t236 = -t198 * t116 + t201 * t117;
t54 = t219 * pkin(8) + qJD(3) * t363 + t236;
t213 = t197 * t51 + t91 * t263 - t93 * t291 + t347 * t54;
t171 = g(3) * t177;
t272 = t171 + t378;
t376 = -t112 * t222 - t213 + t272;
t240 = pkin(5) * t196 - qJ(6) * t200;
t375 = pkin(5) * t290 - qJ(6) * t289 - qJD(6) * t196 - t222 * t240;
t317 = t200 * t94;
t320 = t196 * t96;
t237 = t317 + t320;
t247 = -t39 * t200 - t96 * t290;
t40 = t196 * t379 - t200 * t279 + t221 * t289;
t336 = -t196 * t40 - t94 * t289;
t374 = t222 * t237 + t247 + t336;
t12 = pkin(9) * t279 + t213;
t16 = -pkin(2) * t281 - pkin(3) * t219 + t55 * pkin(4) - pkin(9) * t206 + t185;
t226 = t200 * t12 + t196 * t16 + t64 * t289 - t57 * t290;
t329 = qJ(6) * t53;
t2 = qJD(6) * t285 + t226 + t329;
t260 = t196 * t12 - t200 * t16 + t57 * t289 + t64 * t290;
t349 = pkin(5) * t53;
t4 = qJDD(6) + t260 - t349;
t371 = t4 * t196 + t2 * t200;
t369 = pkin(5) * t221;
t327 = t221 * t94;
t46 = t200 * t53;
t366 = t285 * t290 - t46;
t65 = t197 * t92 + t90;
t248 = pkin(3) * t291 - t65;
t365 = qJ(6) * t221;
t364 = t285 * t221;
t295 = -t198 * t159 + t201 * t160;
t362 = t178 * pkin(4) + t177 * pkin(9);
t288 = t221 * qJD(3);
t360 = t288 - t256;
t25 = -t196 * t57 + t200 * t64;
t296 = qJD(6) - t25;
t18 = -pkin(5) * t285 + t296;
t359 = t18 * t196 + t19 * t200;
t357 = qJ(2) * qJDD(1);
t294 = t380 * t200;
t32 = t33 * t290;
t356 = t18 * t221 + t294 + t32;
t341 = g(3) * t196;
t273 = -t178 * t341 + t196 * t380;
t258 = t197 * t54 + t93 * t263 + t91 * t291 - t347 * t51;
t13 = -t279 * pkin(4) + t258;
t5 = t40 * pkin(5) + t39 * qJ(6) - t96 * qJD(6) + t13;
t339 = t5 * t196;
t355 = -t19 * t221 + t273 - t339;
t354 = t13 * t196 + t26 * t221 + t56 * t289 - t273;
t47 = t56 * t290;
t353 = -t25 * t221 + t294 + t47;
t342 = g(3) * t178;
t352 = -t112 * t221 - t258 - t342 + t380;
t350 = t96 ^ 2;
t348 = pkin(9) * t53;
t131 = t138 * qJD(3);
t346 = pkin(3) * t131;
t338 = t96 * t94;
t334 = t196 * t78 + t200 * t62;
t66 = t347 * t92 - t89;
t69 = pkin(3) * t129 + t78;
t333 = t196 * t69 + t200 * t66;
t135 = t201 * t159;
t254 = -t160 * t198 - t135;
t100 = -pkin(8) * t138 + t254;
t101 = -pkin(8) * t137 + t295;
t74 = t197 * t100 + t347 * t101;
t111 = -t197 * t137 + t347 * t138;
t115 = pkin(3) * t137 - t269;
t220 = -t347 * t137 - t197 * t138;
t75 = -pkin(4) * t220 - pkin(9) * t111 + t115;
t332 = t196 * t75 + t200 * t74;
t330 = pkin(9) * qJD(5);
t328 = t285 * t26;
t182 = t197 * pkin(3) + pkin(9);
t324 = t182 * t53;
t130 = t137 * qJD(3);
t79 = qJD(4) * t220 - t130 * t347 - t197 * t131;
t322 = t196 * t79;
t321 = t196 * t94;
t319 = t200 * t40;
t318 = t200 * t79;
t316 = t200 * t96;
t315 = t248 + t375;
t314 = t375 - t63;
t312 = t285 * t196;
t302 = t196 * t199;
t300 = t199 * t200;
t299 = t200 * t202;
t297 = t202 * t196;
t293 = t194 ^ 2 + t195 ^ 2;
t292 = qJD(3) * t129;
t278 = t347 * pkin(3);
t268 = -t5 - t342;
t267 = qJD(1) * t304;
t265 = t182 * t290;
t264 = -t13 - t342;
t257 = t293 * qJD(1) ^ 2;
t251 = pkin(3) * t263;
t250 = t178 * t241 + t362;
t249 = 0.2e1 * t293;
t121 = t178 * t302 + t299;
t123 = t178 * t297 - t300;
t246 = -g(1) * t121 + g(2) * t123;
t122 = t178 * t300 - t297;
t124 = t178 * t299 + t302;
t245 = g(1) * t122 - g(2) * t124;
t239 = -t324 - t383;
t238 = t18 * t200 - t19 * t196;
t234 = t18 * t289 - t272 + t371;
t231 = t222 * t312 - t366;
t230 = pkin(4) + t241;
t187 = cos(t193);
t175 = pkin(3) * t187;
t229 = t175 + t269 + t362;
t186 = sin(t193);
t227 = t244 * t186;
t223 = t347 * t100 - t197 * t101;
t215 = -qJD(3) * t135 + qJD(2) * t298 + (-qJD(2) * t194 - qJD(3) * t160) * t198;
t81 = -pkin(8) * t131 + t215;
t207 = -t138 * qJD(2) - qJD(3) * t295;
t82 = pkin(8) * t130 + t207;
t28 = qJD(4) * t223 + t197 * t82 + t347 * t81;
t80 = qJD(4) * t111 - t197 * t130 + t131 * t347;
t42 = pkin(4) * t80 - pkin(9) * t79 + t346;
t225 = t196 * t42 + t200 * t28 + t75 * t289 - t74 * t290;
t217 = t228 + t313;
t216 = g(1) * t123 + g(2) * t121 + t177 * t341 - t260;
t214 = t249 * t284 - t244;
t212 = t238 * qJD(5) + t371;
t211 = t33 * t96 + qJDD(6) - t216;
t209 = -g(1) * t124 - g(2) * t122 - t200 * t171 + t226;
t29 = t74 * qJD(4) + t197 * t81 - t347 * t82;
t205 = t244 * t177 * t230 - pkin(9) * t378;
t190 = -pkin(8) - t337;
t183 = -t278 - pkin(4);
t146 = -t269 * qJDD(1) + qJDD(2);
t133 = -t278 - t230;
t98 = qJDD(2) - t271 * pkin(3) + (-pkin(1) + (-pkin(3) * t201 - pkin(2)) * t195) * qJDD(1);
t67 = pkin(5) * t96 + qJ(6) * t94;
t43 = t111 * t240 - t223;
t31 = pkin(5) * t220 + t196 * t74 - t200 * t75;
t30 = -qJ(6) * t220 + t332;
t24 = t196 * t62 - t200 * t78 - t369;
t23 = t334 + t365;
t22 = t196 * t66 - t200 * t69 - t369;
t21 = t333 + t365;
t20 = t285 * t94 - t39;
t8 = t240 * t79 + (qJD(5) * t241 - qJD(6) * t200) * t111 + t29;
t7 = -pkin(5) * t80 + qJD(5) * t332 + t196 * t28 - t200 * t42;
t6 = qJ(6) * t80 - qJD(6) * t220 + t225;
t1 = [qJDD(1), t361, t244, t217 * t195, -t217 * t194, t249 * t357 + t214, t228 * pkin(1) + (t293 * t357 + t214) * qJ(2), t109 * t138 - t129 * t130, -t109 * t137 + t130 * t128 - t129 * t131 + t138 * t219, -qJD(3) * t130 + qJDD(3) * t138, -qJD(3) * t131 - qJDD(3) * t137, 0, qJD(3) * t207 + qJDD(3) * t254 + t147 * t131 + t146 * t137 + t187 * t361 + t219 * t269, -qJD(3) * t215 - qJDD(3) * t295 - t109 * t269 - t147 * t130 + t146 * t138 - t186 * t361, t111 * t206 + t221 * t79, -t111 * t55 + t206 * t220 - t221 * t80 + t222 * t79, t111 * t279 + t280 * t79, t220 * t279 - t280 * t80, 0, t112 * t80 + t115 * t55 + t178 * t361 - t220 * t98 - t222 * t346 + t223 * t279 - t280 * t29, -g(1) * t308 + g(2) * t307 + t98 * t111 + t112 * t79 + t115 * t206 + t221 * t346 - t279 * t74 - t28 * t280, t111 * t247 + t79 * t316, -t237 * t79 + (t37 - t319 + (-t316 + t321) * qJD(5)) * t111, t111 * t46 + t220 * t39 + t80 * t96 + (-t111 * t290 + t318) * t285, -t111 * t45 + t220 * t40 - t80 * t94 + (-t111 * t289 - t322) * t285, -t220 * t53 + t285 * t80, t260 * t220 + t25 * t80 + t29 * t94 - t223 * t40 + ((-qJD(5) * t74 + t42) * t285 + t75 * t53 + t56 * qJD(5) * t111) * t200 + ((-qJD(5) * t75 - t28) * t285 - t74 * t53 + t13 * t111 + t56 * t79) * t196 + t245, -t225 * t285 - t332 * t53 + t226 * t220 - t26 * t80 + t29 * t96 + t223 * t39 + t56 * t318 + (t13 * t200 - t47) * t111 + t246, t33 * t322 - t285 * t7 + t220 * t4 - t18 * t80 - t31 * t53 + t40 * t43 + t8 * t94 + (t289 * t33 + t339) * t111 + t245, -t30 * t40 - t31 * t39 - t6 * t94 + t7 * t96 + t238 * t79 + t361 * t177 + (-qJD(5) * t359 - t196 * t2 + t200 * t4) * t111, -t33 * t318 + t285 * t6 - t220 * t2 + t19 * t80 + t30 * t53 + t39 * t43 - t8 * t96 + (-t5 * t200 + t32) * t111 - t246, t2 * t30 + t19 * t6 + t5 * t43 + t33 * t8 + t4 * t31 + t18 * t7 - g(1) * (-pkin(5) * t122 - qJ(6) * t121) - g(2) * (pkin(5) * t124 + qJ(6) * t123) + (g(1) * t190 - g(2) * t229) * t202 + (g(1) * t229 + g(2) * t190) * t199; 0, 0, 0, -t281, t282, -t257, -qJ(2) * t257 - t228, 0, 0, 0, 0, 0, -t219 + t292 (-t128 - t267) * qJD(3) + t270, 0, 0, 0, 0, 0, t256 + t288 + 0.2e1 * t358, t210 + t287 + 0.2e1 * t388, 0, 0, 0, 0, 0, t231 - t327, -t200 * t253 - t326 - t45, -t196 * t253 - t327 + t46 (t317 - t320) * t222 - t247 + t336, t233 + t326, -t221 * t33 + (-t4 + t392) * t200 + (t18 * t285 + t2) * t196 - t361; 0, 0, 0, 0, 0, 0, 0, t129 * t128, -t128 ^ 2 + t129 ^ 2 (t128 - t267) * qJD(3) + t270, t219 + t292, qJDD(3), -g(3) * t187 - t147 * t129 + t227 + t255, g(3) * t186 + t147 * t128 + t187 * t244 - t236, -t381, t377, t389, t360, t279, t65 * t280 + (t129 * t222 + t279 * t347 - t280 * t291) * pkin(3) + t352, t66 * t280 + (-t129 * t221 - t197 * t279 - t263 * t280) * pkin(3) + t376, t387, t374, t390, t231 + t327, -t364, t183 * t40 + t248 * t94 + t264 * t200 + t239 * t196 + ((-qJD(5) * t182 - t69) * t200 + (-t251 + t66) * t196) * t285 + t353, -t183 * t39 + t248 * t96 + t239 * t200 + (-t200 * t251 + t265 + t333) * t285 + t354, t133 * t40 + t315 * t94 + t268 * t200 + (-t324 - t384) * t196 + (-t182 * t289 - t196 * t251 + t22) * t285 + t356, t21 * t94 - t22 * t96 + (-t94 * t251 - t222 * t18 + (qJD(5) * t96 - t40) * t182) * t200 + (t96 * t251 + t222 * t19 - t182 * t39 + (t182 * t94 - t19) * qJD(5)) * t196 + t234, t133 * t39 - t315 * t96 + (-t21 - t265) * t285 + (t251 * t285 + t324 - t391) * t200 + t355, t5 * t133 - t19 * t21 - t18 * t22 - g(3) * (t175 + t250) + t315 * t33 + (t263 * t359 + t227) * pkin(3) + t212 * t182 + t205; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t381, t377, t389, t360, t279, t280 * t63 + t352, t280 * t62 + t376, t387, t374, t390, -t285 * t312 + t327 + t46, -t364, -pkin(4) * t40 - t63 * t94 + (t285 * t62 - t348 - t383) * t196 + ((-t78 - t330) * t285 + t264) * t200 + t353, pkin(4) * t39 + pkin(9) * t366 + t285 * t334 - t382 * t56 - t63 * t96 + t354, t285 * t24 - t230 * t40 + t314 * t94 + (-t348 - t384) * t196 + (-t285 * t330 + t268) * t200 + t356, -t19 * t290 + t23 * t94 - t24 * t96 - t238 * t222 + (-t37 - t319 + (t316 + t321) * qJD(5)) * pkin(9) + t234, -t230 * t39 - t314 * t96 + (-pkin(9) * t290 - t23) * t285 + (t348 - t391) * t200 + t355, pkin(9) * t212 - g(3) * t250 - t18 * t24 - t19 * t23 - t230 * t5 + t314 * t33 + t205; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t338, -t94 ^ 2 + t350, t20, t285 * t96 - t40, t53, -t56 * t96 + t216 + t328, t25 * t285 + t56 * t94 - t209, -t67 * t94 - t211 + t328 + 0.2e1 * t349, pkin(5) * t39 - qJ(6) * t40 + (t19 - t26) * t96 + (t18 - t296) * t94, 0.2e1 * t329 - t33 * t94 + t67 * t96 + (0.2e1 * qJD(6) - t25) * t285 + t209, t2 * qJ(6) - t4 * pkin(5) - t33 * t67 - t18 * t26 - g(1) * (-pkin(5) * t123 + qJ(6) * t124) - g(2) * (-pkin(5) * t121 + qJ(6) * t122) + t296 * t19 + t240 * t171; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t338 - t53, t20, -t350 - t253, t211 - t349 - t392;];
tau_reg  = t1;
