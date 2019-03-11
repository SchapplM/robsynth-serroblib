% Calculate inertial parameters regressor of inverse dynamics joint torque vector for
% S6RPRRPP6
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
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d3,d4,theta5]';
% 
% Output:
% tau_reg [6x(6*10)]
%   inertial parameter regressor of inverse dynamics joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 04:49
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tau_reg = S6RPRRPP6_invdynJ_fixb_reg2_slag_vp(qJ, qJD, qJDD, g, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(3,1),zeros(9,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRRPP6_invdynJ_fixb_reg2_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RPRRPP6_invdynJ_fixb_reg2_slag_vp: qJD has to be [6x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [6 1]), ...
  'S6RPRRPP6_invdynJ_fixb_reg2_slag_vp: qJDD has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RPRRPP6_invdynJ_fixb_reg2_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S6RPRRPP6_invdynJ_fixb_reg2_slag_vp: pkin has to be [9x1] (double)');

%% Symbolic Calculation
% From invdyn_joint_fixb_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 04:48:35
% EndTime: 2019-03-09 04:48:49
% DurationCPUTime: 7.37s
% Computational Cost: add. (7639->618), mult. (15240->751), div. (0->0), fcn. (9943->10), ass. (0->299)
t216 = sin(qJ(3));
t333 = qJD(1) * t216;
t182 = qJD(4) + t333;
t219 = cos(qJ(3));
t396 = g(3) * t216;
t220 = cos(qJ(1));
t207 = g(2) * t220;
t217 = sin(qJ(1));
t208 = g(1) * t217;
t406 = t208 - t207;
t237 = t219 * t406 - t396;
t221 = -pkin(1) - pkin(7);
t177 = qJD(1) * t221 + qJD(2);
t328 = qJD(3) * t216;
t285 = -qJDD(3) * pkin(3) + t177 * t328;
t173 = qJDD(1) * t221 + qJDD(2);
t366 = t173 * t219;
t96 = t285 - t366;
t231 = -t96 - t237;
t411 = -pkin(8) * qJD(4) * t182 + t231;
t215 = sin(qJ(4));
t218 = cos(qJ(4));
t320 = t218 * qJD(3);
t332 = qJD(1) * t219;
t146 = t215 * t332 - t320;
t329 = qJD(3) * t215;
t148 = t218 * t332 + t329;
t213 = sin(pkin(9));
t372 = cos(pkin(9));
t84 = t372 * t146 + t148 * t213;
t376 = t84 * t182;
t300 = t216 * t320;
t321 = qJD(4) * t219;
t242 = t215 * t321 + t300;
t315 = t219 * qJDD(1);
t77 = qJD(1) * t242 - qJD(4) * t320 - t215 * qJDD(3) - t218 * t315;
t301 = t215 * t328;
t325 = qJD(4) * t148;
t78 = -qJD(1) * t301 - t218 * qJDD(3) + t215 * t315 + t325;
t41 = -t213 * t78 - t372 * t77;
t22 = t41 + t376;
t288 = t372 * t218;
t363 = t213 * t215;
t247 = t288 - t363;
t318 = qJD(1) * qJD(3);
t293 = t219 * t318;
t316 = t216 * qJDD(1);
t140 = qJDD(4) + t293 + t316;
t209 = qJ(4) + pkin(9);
t197 = sin(t209);
t389 = qJ(5) + pkin(8);
t164 = t389 * t218;
t99 = t164 * t372 - t363 * t389;
t410 = -t99 * t140 + t197 * t237;
t248 = -t213 * t146 + t148 * t372;
t409 = t84 * t248;
t289 = t372 * t215;
t143 = t213 * t218 + t289;
t124 = t143 * qJD(1);
t323 = qJD(4) * t216;
t327 = qJD(3) * t219;
t373 = t143 * t323 - t247 * t327 + t124;
t123 = t143 * qJD(4);
t343 = t216 * t124 + t123;
t324 = qJD(4) * t215;
t125 = qJD(4) * t288 - t213 * t324;
t405 = t247 * qJD(1);
t342 = t405 * t216 + t125;
t297 = t218 * t321;
t408 = t297 - t301;
t355 = t217 * t218;
t357 = t216 * t220;
t129 = t215 * t357 + t355;
t393 = t248 ^ 2;
t290 = qJD(4) * t389;
t319 = t218 * qJD(5);
t117 = -t215 * t290 + t319;
t241 = -t215 * qJD(5) - t218 * t290;
t358 = t216 * t218;
t274 = pkin(3) * t219 + pkin(8) * t216;
t150 = t274 * qJD(1);
t361 = t215 * t219;
t93 = t218 * t150 - t177 * t361;
t63 = (pkin(4) * t219 + qJ(5) * t358) * qJD(1) + t93;
t302 = t215 * t333;
t352 = t218 * t219;
t94 = t215 * t150 + t177 * t352;
t72 = qJ(5) * t302 + t94;
t387 = (t241 - t63) * t372 + (-t117 + t72) * t213;
t273 = pkin(3) * t216 - pkin(8) * t219;
t156 = qJ(2) + t273;
t126 = t156 * qJD(1);
t155 = t216 * t177;
t134 = qJD(3) * pkin(8) + t155;
t322 = qJD(4) * t218;
t141 = qJD(3) * t274 + qJD(2);
t89 = qJD(1) * t141 + qJDD(1) * t156;
t97 = qJDD(3) * pkin(8) + t173 * t216 + t177 * t327;
t29 = t126 * t322 - t134 * t324 + t215 * t89 + t218 * t97;
t75 = t218 * t126 - t134 * t215;
t407 = -t182 * t75 + t29;
t211 = t216 ^ 2;
t212 = t219 ^ 2;
t335 = t211 + t212;
t286 = t335 * t173;
t278 = -t155 + (t302 + t324) * pkin(4);
t356 = t216 * t221;
t102 = t215 * t156 + t218 * t356;
t40 = -t213 * t77 + t372 * t78;
t404 = t40 * pkin(5) - t41 * qJ(6) - t248 * qJD(6);
t114 = t247 * t216;
t403 = -t114 * t140 + t182 * t373 - t219 * t41 + t248 * t328;
t402 = -t140 * t247 + t182 * t343 - t332 * t84;
t135 = -qJD(3) * pkin(3) - t177 * t219;
t91 = pkin(4) * t146 + qJD(5) + t135;
t31 = pkin(5) * t84 - qJ(6) * t248 + t91;
t401 = -t31 * t248 - qJDD(6);
t400 = t143 * t40 - t247 * t41 + t248 * t343 + t342 * t84;
t399 = 0.2e1 * qJ(2);
t76 = t126 * t215 + t134 * t218;
t80 = t218 * t89;
t30 = -qJD(4) * t76 - t215 * t97 + t80;
t11 = t140 * pkin(4) + t77 * qJ(5) - t148 * qJD(5) + t30;
t18 = -qJ(5) * t78 - qJD(5) * t146 + t29;
t3 = t372 * t11 - t213 * t18;
t4 = t213 * t11 + t372 * t18;
t397 = pkin(5) * t140;
t395 = g(3) * t219;
t59 = -qJ(5) * t146 + t76;
t55 = t372 * t59;
t58 = -qJ(5) * t148 + t75;
t26 = t213 * t58 + t55;
t394 = t26 * t248;
t391 = t84 ^ 2;
t235 = -qJD(4) * t102 + t218 * t141;
t292 = -t215 * t221 + pkin(4);
t39 = qJ(5) * t300 + (qJ(5) * t324 + qJD(3) * t292 - t319) * t219 + t235;
t326 = qJD(3) * t221;
t298 = t219 * t326;
t306 = t215 * t141 + t156 * t322 + t218 * t298;
t44 = -qJ(5) * t297 + (-qJD(5) * t219 + (qJ(5) * qJD(3) - qJD(4) * t221) * t216) * t215 + t306;
t15 = t213 * t39 + t372 * t44;
t388 = pkin(5) * t343 - qJ(6) * t342 - qJD(6) * t143 + t278;
t53 = pkin(4) * t182 + t58;
t24 = t213 * t53 + t55;
t36 = t213 * t63 + t372 * t72;
t386 = pkin(5) * t332 - t387;
t33 = qJ(6) * t332 + t36;
t70 = t117 * t372 + t213 * t241;
t385 = t70 - t33;
t384 = t70 - t36;
t139 = t218 * t156;
t82 = -qJ(5) * t352 + t216 * t292 + t139;
t92 = -qJ(5) * t361 + t102;
t49 = t213 * t82 + t372 * t92;
t382 = t182 * t248;
t381 = t213 * t59;
t379 = t76 * t182;
t378 = t77 * t215;
t377 = t78 * t218;
t287 = qJD(3) * t372;
t299 = t219 * t320;
t374 = -t125 * t216 - t213 * t299 - t287 * t361 - t405;
t371 = pkin(1) * qJDD(1);
t370 = t146 * t182;
t369 = t148 * t146;
t368 = t148 * t182;
t367 = t148 * t218;
t365 = t182 * t215;
t198 = cos(t209);
t364 = t198 * t217;
t362 = t215 * t140;
t360 = t215 * t220;
t359 = t216 * t217;
t354 = t217 * t219;
t353 = t218 * t140;
t351 = t218 * t220;
t350 = t219 * t389;
t349 = t219 * t220;
t348 = t220 * t198;
t222 = qJD(3) ^ 2;
t347 = t221 * t222;
t223 = qJD(1) ^ 2;
t346 = t223 * qJ(2);
t27 = t372 * t58 - t381;
t344 = qJD(6) - t27;
t341 = t129 * pkin(4);
t340 = g(1) * t349 + g(2) * t354;
t311 = 0.2e1 * qJD(1) * qJD(2);
t339 = (qJDD(1) * qJ(2) + t311) * qJ(2);
t338 = t220 * pkin(1) + t217 * qJ(2);
t336 = t211 - t212;
t334 = -t222 - t223;
t331 = qJD(3) * t146;
t330 = qJD(3) * t148;
t317 = qJDD(3) * t216;
t314 = g(1) * t364;
t310 = t215 * t359;
t308 = t215 * t356;
t307 = t219 * t223 * t216;
t195 = pkin(4) * t218 + pkin(3);
t204 = t220 * qJ(2);
t305 = t195 * t357 - t349 * t389 + t204;
t303 = t220 * pkin(7) + t338;
t296 = t182 * t332;
t291 = -g(1) * t354 + t396;
t284 = -t173 + t346;
t183 = pkin(4) * t361;
t144 = -t219 * t221 + t183;
t283 = t335 * qJDD(1);
t282 = qJDD(2) - t371;
t281 = qJD(1) + t323;
t280 = g(2) * t303;
t279 = t216 * t293;
t98 = t213 * t164 + t289 * t389;
t277 = g(2) * t219 * t348 - t98 * t140 + t198 * t396;
t108 = t197 * t359 - t348;
t110 = t197 * t357 + t364;
t272 = -g(1) * t110 - g(2) * t108;
t109 = t197 * t220 + t198 * t359;
t111 = -t197 * t217 + t216 * t348;
t271 = -g(1) * t111 - g(2) * t109;
t270 = g(1) * t220 + g(2) * t217;
t268 = -t391 - t393;
t267 = -t391 + t393;
t266 = -t406 - t346;
t113 = t143 * t219;
t112 = t143 * t216;
t66 = qJD(3) * t112 - t125 * t219;
t265 = t113 * t40 - t66 * t84;
t264 = pkin(5) * t198 + qJ(6) * t197;
t262 = t215 * t76 + t218 * t75;
t261 = t215 * t75 - t218 * t76;
t260 = g(1) * t359 + t395;
t259 = qJDD(1) * t399 + t311;
t258 = t40 + t382;
t257 = -t40 + t382;
t256 = t195 + t264;
t254 = (-pkin(4) * t215 + t221) * t208;
t253 = pkin(4) * t360 + t195 * t359 - t217 * t350 + t303;
t14 = -t213 * t44 + t372 * t39;
t23 = t372 * t53 - t381;
t48 = -t213 * t92 + t372 * t82;
t251 = t182 * t322 + t362;
t250 = -t182 * t324 + t353;
t249 = -t247 * t40 + t343 * t84;
t100 = pkin(4) * t408 + t216 * t326;
t188 = g(2) * t357;
t246 = t188 - t260;
t245 = qJDD(3) * t221 + t318 * t399;
t244 = t78 * pkin(4) + qJDD(5) + t285;
t243 = g(1) * t108 - g(2) * t110 + t197 * t395 + t3;
t240 = -pkin(8) * t140 + t135 * t182;
t238 = -t41 + t376;
t115 = -t213 * t361 + t219 * t288;
t68 = t123 * t219 - t213 * t301 + t287 * t358;
t234 = t113 * t41 + t115 * t40 - t248 * t66 - t68 * t84;
t233 = g(1) * t109 - g(2) * t111 + t198 * t395 - t4;
t50 = t244 - t366;
t232 = t26 * t182 + t243;
t230 = -t112 * t140 + t182 * t374 - t219 * t40 + t84 * t328;
t229 = t259 - t270;
t228 = t112 * t41 - t114 * t40 - t248 * t374 + t373 * t84;
t227 = -t99 * t40 + t98 * t41 - t70 * t84 + t246;
t226 = t113 * t140 - t182 * t66 + t216 * t40 + t327 * t84;
t225 = -qJD(4) * t262 - t30 * t215 + t29 * t218;
t224 = (-t173 - t207) * t219 + t244 - t291;
t201 = qJDD(3) * t219;
t193 = -pkin(4) * t372 - pkin(5);
t191 = pkin(4) * t213 + qJ(6);
t186 = pkin(4) * t351;
t171 = t389 * t357;
t153 = t195 * t354;
t133 = t140 * qJ(6);
t130 = -t215 * t217 + t216 * t351;
t128 = t216 * t355 + t360;
t127 = -t310 + t351;
t101 = t139 - t308;
t95 = t140 * t216 + t182 * t327;
t81 = -pkin(5) * t247 - qJ(6) * t143 - t195;
t61 = -t215 * t298 + t235;
t60 = -qJD(4) * t308 + t306;
t57 = pkin(5) * t113 - qJ(6) * t115 + t144;
t45 = -t216 * pkin(5) - t48;
t43 = qJ(6) * t216 + t49;
t37 = pkin(4) * t148 + pkin(5) * t248 + qJ(6) * t84;
t25 = t143 * t140 + t182 * t342 - t248 * t332;
t21 = qJ(6) * t182 + t24;
t20 = -t182 * pkin(5) + qJD(6) - t23;
t19 = -pkin(5) * t66 + qJ(6) * t68 - qJD(6) * t115 + t100;
t13 = -pkin(5) * t327 - t14;
t12 = qJ(6) * t327 + qJD(6) * t216 + t15;
t10 = t115 * t41 - t248 * t68;
t7 = t143 * t41 + t248 * t342;
t6 = t115 * t140 - t182 * t68 + t216 * t41 + t248 * t327;
t5 = t50 + t404;
t2 = qJDD(6) - t397 - t3;
t1 = qJD(6) * t182 + t133 + t4;
t8 = [0, 0, 0, 0, 0, qJDD(1), t406, t270, 0, 0, qJDD(1), 0, 0, 0, 0, 0, 0, qJDD(2) - t406 - 0.2e1 * t371, t229, -t282 * pkin(1) - g(1) * (-t217 * pkin(1) + t204) - g(2) * t338 + t339, qJDD(1) * t212 - 0.2e1 * t279, -0.2e1 * t216 * t315 + 0.2e1 * t318 * t336, -t216 * t222 + t201, qJDD(1) * t211 + 0.2e1 * t279, -t219 * t222 - t317, 0, t245 * t219 + (t229 - t347) * t216, -t245 * t216 + (t259 - t347) * t219 - t340, -t221 * t283 - t286 + t406, -g(1) * (t217 * t221 + t204) - t280 + t221 * t286 + t339, -t148 * t242 - t352 * t77 (t146 * t218 + t148 * t215) * t328 + (t378 - t377 + (t146 * t215 - t367) * qJD(4)) * t219 (-t182 * t320 - t77) * t216 + (t250 + t330) * t219, t146 * t408 + t78 * t361 (t182 * t329 - t78) * t216 + (-t251 - t331) * t219, t95, -g(1) * t130 - g(2) * t128 + t101 * t140 + t61 * t182 + (t30 + (-t135 * t215 + t146 * t221) * qJD(3)) * t216 + (qJD(3) * t75 + t135 * t322 + t96 * t215 - t221 * t78) * t219, g(1) * t129 - g(2) * t127 - t102 * t140 - t60 * t182 + (-t29 + (-t135 * t218 + t148 * t221) * qJD(3)) * t216 + (-qJD(3) * t76 - t135 * t324 + t96 * t218 + t221 * t77) * t219, t101 * t77 - t102 * t78 - t60 * t146 - t61 * t148 + t262 * t328 + (qJD(4) * t261 - t215 * t29 - t218 * t30) * t219 + t340, t29 * t102 + t76 * t60 + t30 * t101 + t75 * t61 - g(1) * (pkin(3) * t357 - pkin(8) * t349 + t204) - t280 + (t135 * t328 - t96 * t219) * t221 + (-g(1) * t221 - g(2) * t273) * t217, t10, -t234, t6, t265, -t226, t95, t100 * t84 + t113 * t50 + t14 * t182 + t140 * t48 + t144 * t40 + t216 * t3 + t23 * t327 - t66 * t91 + t271, t100 * t248 + t115 * t50 - t140 * t49 + t144 * t41 - t15 * t182 - t216 * t4 - t24 * t327 - t68 * t91 - t272, -t113 * t4 - t115 * t3 - t14 * t248 - t15 * t84 + t23 * t68 + t24 * t66 - t40 * t49 - t41 * t48 + t340, -g(1) * t305 - g(2) * t253 + t91 * t100 + t23 * t14 + t50 * t144 + t24 * t15 + t3 * t48 + t4 * t49 - t254, t10, t6, t234, t95, t226, t265, t113 * t5 - t13 * t182 - t140 * t45 + t19 * t84 - t2 * t216 - t20 * t327 - t31 * t66 + t40 * t57 + t271, -t1 * t113 + t115 * t2 - t12 * t84 + t13 * t248 - t20 * t68 + t21 * t66 - t40 * t43 + t41 * t45 + t340, t1 * t216 - t115 * t5 + t12 * t182 + t140 * t43 - t19 * t248 + t21 * t327 + t31 * t68 - t41 * t57 + t272, t1 * t43 + t21 * t12 + t5 * t57 + t31 * t19 + t2 * t45 + t20 * t13 - g(1) * (pkin(5) * t111 + qJ(6) * t110 + t305) - g(2) * (pkin(5) * t109 + qJ(6) * t108 + t253) - t254; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, qJDD(1), -t223, t266 + t282, 0, 0, 0, 0, 0, 0, t216 * t334 + t201, t219 * t334 - t317, -t283, t286 + t266, 0, 0, 0, 0, 0, 0, -t219 * t78 + (t331 - t362) * t216 + (-t215 * t327 - t218 * t281) * t182, t219 * t77 + (t330 - t353) * t216 + (t215 * t281 - t299) * t182 (-t146 * t327 + t148 * t281 - t78 * t216) * t218 + (t146 * t281 + t148 * t327 - t77 * t216) * t215, -t262 * qJD(1) + (-qJD(3) * t261 - t96) * t219 + (qJD(3) * t135 + t225) * t216 - t406, 0, 0, 0, 0, 0, 0, t230, t403, t228, -t3 * t112 + t4 * t114 - t50 * t219 + t23 * t374 - t24 * t373 + t328 * t91 - t406, 0, 0, 0, 0, 0, 0, t230, t228, -t403, t1 * t114 + t2 * t112 - t20 * t374 - t21 * t373 - t5 * t219 + t31 * t328 - t406; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t307, -t336 * t223, t315, -t307, -t316, qJDD(3) (-t284 + t207) * t219 + t291, t395 - t188 + (t284 + t208) * t216, 0, 0, t182 * t367 - t378 (-t77 - t370) * t218 + (-t78 - t368) * t215 (-t148 * t219 + t182 * t358) * qJD(1) + t251, t146 * t365 - t377 (t146 * t219 - t216 * t365) * qJD(1) + t250, -t296, -pkin(3) * t78 - t146 * t155 - t93 * t182 + t240 * t215 + t218 * t411 - t75 * t332, pkin(3) * t77 - t148 * t155 + t94 * t182 - t215 * t411 + t240 * t218 + t76 * t332, t94 * t146 + t93 * t148 + ((-t78 + t325) * pkin(8) + t407) * t218 + (-t30 - t379 + (qJD(4) * t146 - t77) * pkin(8)) * t215 + t246, -t135 * t155 - t75 * t93 - t76 * t94 + t231 * pkin(3) + (-t216 * t406 + t225 - t395) * pkin(8), t7, -t400, t25, t249, -t402, -t296, -t50 * t247 - t195 * t40 + t343 * t91 + t278 * t84 + (-qJD(1) * t23 - t314) * t219 + t387 * t182 + t277, t50 * t143 - t182 * t384 - t195 * t41 + t24 * t332 + t248 * t278 + t342 * t91 + t410, -t3 * t143 - t23 * t342 - t24 * t343 + t247 * t4 - t248 * t387 + t36 * t84 + t227, t4 * t99 - t3 * t98 - t50 * t195 - g(1) * (t359 * t389 + t153) - g(2) * (-t195 * t349 - t171) - g(3) * (-t195 * t216 + t350) + t278 * t91 + t384 * t24 + t387 * t23, t7, t25, t400, -t296, t402, t249, -t5 * t247 + t81 * t40 + t388 * t84 + t343 * t31 + (qJD(1) * t20 - t314) * t219 - t386 * t182 + t277, t1 * t247 + t2 * t143 + t20 * t342 - t21 * t343 + t248 * t386 + t33 * t84 + t227, -t5 * t143 + t182 * t385 - t21 * t332 - t248 * t388 - t31 * t342 - t81 * t41 - t410, -g(1) * t153 + g(2) * t171 + t1 * t99 + t2 * t98 + t5 * t81 + t388 * t31 + t385 * t21 + t386 * t20 + (g(3) * t256 - t208 * t389) * t216 + (-g(3) * t389 + t207 * t256 - t208 * t264) * t219; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t369, -t146 ^ 2 + t148 ^ 2, t370 - t77, -t369, t368 - t78, t140, -t134 * t322 - g(1) * t127 - g(2) * t129 - t135 * t148 + t379 + t80 + (-qJD(4) * t126 + t395 - t97) * t215, g(1) * t128 - g(2) * t130 + g(3) * t352 + t135 * t146 - t407, 0, 0, t409, t267, t22, -t409, t257, t140, -t91 * t248 + (t140 * t372 - t148 * t84) * pkin(4) + t232, t27 * t182 + t91 * t84 + (-t140 * t213 - t148 * t248) * pkin(4) + t233, t24 * t248 - t394 + (-t213 * t40 - t372 * t41) * pkin(4) + (-t23 + t27) * t84, -t24 * t27 + t23 * t26 - g(1) * t186 - g(2) * t341 + (-t91 * t148 + t4 * t213 + t215 * t260 + t3 * t372) * pkin(4), t409, t22, -t267, t140, -t257, -t409, -t37 * t84 + (pkin(5) - t193) * t140 + t232 + t401, -t191 * t40 + t193 * t41 + t21 * t248 - t394 + (t20 - t344) * t84, t191 * t140 - t31 * t84 + t37 * t248 + t133 + (0.2e1 * qJD(6) - t27) * t182 - t233, t1 * t191 + t2 * t193 - t31 * t37 - t20 * t26 - g(1) * (-pkin(4) * t310 - pkin(5) * t108 + qJ(6) * t109 + t186) - g(2) * (pkin(5) * t110 - qJ(6) * t111 + t341) - g(3) * (-t183 + (-pkin(5) * t197 + qJ(6) * t198) * t219) + t344 * t21; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t258, -t238, t268, t23 * t248 + t24 * t84 + t224, 0, 0, 0, 0, 0, 0, t258, t268, t238, -t20 * t248 + t21 * t84 + t224 + t404; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t140 + t409, t22, -t182 ^ 2 - t393, -t182 * t21 - t243 - t397 - t401;];
tau_reg  = t8;
