% Calculate inertial parameters regressor of inverse dynamics joint torque vector for
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
% tau_reg [5x(5*10)]
%   inertial parameter regressor of inverse dynamics joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 22:31
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tau_reg = S5RRRRR9_invdynJ_fixb_reg2_slag_vp(qJ, qJD, qJDD, g, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(5,1),zeros(3,1),zeros(9,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRRRR9_invdynJ_fixb_reg2_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RRRRR9_invdynJ_fixb_reg2_slag_vp: qJD has to be [5x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [5 1]), ...
  'S5RRRRR9_invdynJ_fixb_reg2_slag_vp: qJDD has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RRRRR9_invdynJ_fixb_reg2_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5RRRRR9_invdynJ_fixb_reg2_slag_vp: pkin has to be [9x1] (double)');

%% Symbolic Calculation
% From invdyn_joint_fixb_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 22:29:58
% EndTime: 2019-12-31 22:30:22
% DurationCPUTime: 10.75s
% Computational Cost: add. (13180->633), mult. (29806->849), div. (0->0), fcn. (20920->14), ass. (0->285)
t279 = sin(qJ(2));
t283 = cos(qJ(2));
t318 = pkin(2) * t279 - pkin(7) * t283;
t215 = t318 * qJD(1);
t282 = cos(qJ(3));
t278 = sin(qJ(3));
t355 = qJD(1) * t279;
t330 = t278 * t355;
t150 = pkin(6) * t330 + t282 * t215;
t371 = t282 * t283;
t307 = pkin(3) * t279 - pkin(8) * t371;
t414 = pkin(8) + pkin(7);
t336 = qJD(3) * t414;
t447 = -qJD(1) * t307 - t282 * t336 - t150;
t192 = t278 * t215;
t375 = t279 * t282;
t376 = t278 * t283;
t446 = t192 + (-pkin(6) * t375 - pkin(8) * t376) * qJD(1) + t278 * t336;
t345 = t282 * qJD(2);
t206 = t330 - t345;
t329 = t282 * t355;
t346 = t278 * qJD(2);
t208 = t329 + t346;
t277 = sin(qJ(4));
t281 = cos(qJ(4));
t130 = t281 * t206 + t277 * t208;
t276 = sin(qJ(5));
t312 = -t277 * t206 + t281 * t208;
t413 = cos(qJ(5));
t431 = -t276 * t130 + t312 * t413;
t67 = t413 * t130 + t276 * t312;
t400 = t67 * t431;
t378 = t277 * t278;
t209 = -t281 * t282 + t378;
t341 = qJD(3) + qJD(4);
t344 = t283 * qJD(1);
t348 = qJD(4) * t281;
t350 = qJD(3) * t282;
t361 = -t209 * t344 - t281 * t350 - t282 * t348 + t341 * t378;
t210 = t277 * t282 + t281 * t278;
t145 = t341 * t210;
t360 = -t210 * t344 + t145;
t246 = -qJD(3) + t344;
t342 = t279 * qJDD(1);
t257 = pkin(6) * t342;
t343 = qJD(1) * qJD(2);
t327 = t283 * t343;
t180 = -qJDD(2) * pkin(2) + pkin(6) * t327 + t257;
t280 = sin(qJ(1));
t284 = cos(qJ(1));
t317 = g(1) * t284 + g(2) * t280;
t403 = g(3) * t283;
t297 = t279 * t317 - t403;
t293 = -t180 + t297;
t445 = pkin(7) * qJD(3) * t246 + t293;
t437 = t431 ^ 2 - t67 ^ 2;
t351 = qJD(3) * t279;
t438 = qJD(1) * t351 - qJDD(2);
t114 = -qJD(3) * t345 + (-t327 - t342) * t282 + t438 * t278;
t115 = t278 * (qJD(2) * (qJD(3) + t344) + t342) + t438 * t282;
t349 = qJD(4) * t277;
t306 = -t277 * t114 + t281 * t115 - t206 * t349 + t208 * t348;
t328 = qJD(5) * t413;
t347 = qJD(5) * t276;
t45 = t281 * t114 + t277 * t115 + t206 * t348 + t208 * t349;
t10 = t130 * t328 + t276 * t306 + t312 * t347 + t413 * t45;
t236 = -qJD(4) + t246;
t228 = -qJD(5) + t236;
t435 = -t67 * t228 - t10;
t428 = pkin(9) * t312;
t319 = t283 * pkin(2) + t279 * pkin(7);
t222 = -pkin(1) - t319;
t199 = t222 * qJD(1);
t259 = pkin(6) * t344;
t230 = qJD(2) * pkin(7) + t259;
t140 = t282 * t199 - t278 * t230;
t100 = -t208 * pkin(8) + t140;
t90 = -t246 * pkin(3) + t100;
t141 = t278 * t199 + t282 * t230;
t101 = -t206 * pkin(8) + t141;
t93 = t277 * t101;
t46 = t281 * t90 - t93;
t35 = t46 - t428;
t31 = -t236 * pkin(4) + t35;
t441 = t130 * pkin(9);
t95 = t281 * t101;
t47 = t277 * t90 + t95;
t36 = t47 - t441;
t393 = t276 * t36;
t12 = t413 * t31 - t393;
t338 = t413 * t36;
t13 = t276 * t31 + t338;
t444 = -t12 * t67 + t13 * t431;
t263 = t283 * qJDD(1);
t325 = t279 * t343;
t424 = -t325 + t263;
t203 = qJDD(3) - t424;
t198 = qJDD(4) + t203;
t218 = t318 * qJD(2);
t146 = qJD(1) * t218 + qJDD(1) * t222;
t135 = t282 * t146;
t179 = t424 * pkin(6) + qJDD(2) * pkin(7);
t58 = -qJD(3) * t141 - t278 * t179 + t135;
t37 = t203 * pkin(3) + t114 * pkin(8) + t58;
t352 = qJD(3) * t278;
t57 = t278 * t146 + t282 * t179 + t199 * t350 - t230 * t352;
t41 = -t115 * pkin(8) + t57;
t9 = -qJD(4) * t47 - t277 * t41 + t281 * t37;
t6 = t198 * pkin(4) + t45 * pkin(9) + t9;
t8 = -t101 * t349 + t277 * t37 + t281 * t41 + t90 * t348;
t7 = -pkin(9) * t306 + t8;
t1 = t276 * t6 + t31 * t328 - t36 * t347 + t413 * t7;
t275 = qJ(3) + qJ(4);
t267 = qJ(5) + t275;
t254 = cos(t267);
t253 = sin(t267);
t368 = t284 * t253;
t372 = t280 * t283;
t157 = -t254 * t372 + t368;
t367 = t284 * t254;
t159 = t280 * t253 + t283 * t367;
t404 = g(3) * t279;
t229 = -qJD(2) * pkin(2) + pkin(6) * t355;
t153 = t206 * pkin(3) + t229;
t87 = t130 * pkin(4) + t153;
t434 = g(1) * t159 - g(2) * t157 + t254 * t404 + t87 * t67 - t1;
t231 = t414 * t278;
t232 = t414 * t282;
t395 = -t231 * t348 - t232 * t349 + t447 * t277 - t446 * t281;
t149 = -t277 * t231 + t281 * t232;
t394 = -t149 * qJD(4) + t446 * t277 + t447 * t281;
t11 = qJD(5) * t431 - t276 * t45 + t413 * t306;
t418 = -t228 * t431 - t11;
t156 = t253 * t372 + t367;
t158 = t280 * t254 - t283 * t368;
t2 = -qJD(5) * t13 - t276 * t7 + t413 * t6;
t420 = -g(1) * t158 + g(2) * t156 + t253 * t404 - t431 * t87 + t2;
t440 = -pkin(4) * t355 + t361 * pkin(9) + t394;
t439 = t360 * pkin(9) - t395;
t387 = t130 * t312;
t334 = t283 * t346;
t298 = t279 * t350 + t334;
t436 = -t130 ^ 2 + t312 ^ 2;
t433 = -t236 * t130 - t45;
t266 = cos(t275);
t265 = sin(t275);
t366 = t284 * t265;
t169 = -t266 * t372 + t366;
t365 = t284 * t266;
t171 = t280 * t265 + t283 * t365;
t432 = g(1) * t171 - g(2) * t169 + t130 * t153 + t266 * t404 - t8;
t427 = t140 * t246 + t57;
t205 = t282 * t222;
t139 = -pkin(8) * t375 + t205 + (-pkin(6) * t278 - pkin(3)) * t283;
t248 = pkin(6) * t371;
t162 = t278 * t222 + t248;
t377 = t278 * t279;
t147 = -pkin(8) * t377 + t162;
t76 = t277 * t139 + t281 * t147;
t402 = t278 * pkin(3);
t320 = pkin(3) * t352 - t344 * t402 - t259;
t333 = t283 * t345;
t423 = -t278 * t351 + t333;
t363 = t284 * t282;
t187 = t278 * t372 + t363;
t364 = t284 * t278;
t189 = t280 * t282 - t283 * t364;
t422 = -g(1) * t189 + g(2) * t187;
t168 = t265 * t372 + t365;
t170 = t280 * t266 - t283 * t366;
t421 = -g(1) * t170 + g(2) * t168 + t265 * t404;
t419 = -t153 * t312 + t421 + t9;
t417 = -t236 * t312 - t306;
t415 = -0.2e1 * pkin(1);
t410 = g(1) * t280;
t405 = g(2) * t284;
t268 = t279 * pkin(6);
t401 = t282 * pkin(3);
t148 = -t281 * t231 - t277 * t232;
t109 = -t210 * pkin(9) + t148;
t110 = -pkin(9) * t209 + t149;
t59 = t413 * t109 - t276 * t110;
t399 = t59 * qJD(5) + t440 * t276 - t439 * t413;
t60 = t276 * t109 + t413 * t110;
t398 = -t60 * qJD(5) + t439 * t276 + t440 * t413;
t397 = t209 * t328 + t210 * t347 + t360 * t276 + t361 * t413;
t137 = -t276 * t209 + t413 * t210;
t396 = t137 * qJD(5) - t361 * t276 + t360 * t413;
t53 = t281 * t100 - t93;
t255 = t281 * pkin(3) + pkin(4);
t379 = t276 * t277;
t52 = -t277 * t100 - t95;
t38 = t52 + t441;
t39 = t53 - t428;
t392 = -t276 * t38 - t413 * t39 + t255 * t328 + (-t277 * t347 + (t413 * t281 - t379) * qJD(4)) * pkin(3);
t335 = t413 * t277;
t391 = t276 * t39 - t413 * t38 - t255 * t347 + (-t277 * t328 + (-t276 * t281 - t335) * qJD(4)) * pkin(3);
t390 = pkin(6) * qJDD(1);
t389 = t114 * t278;
t388 = t115 * t282;
t385 = t141 * t246;
t384 = t206 * t246;
t383 = t206 * t278;
t382 = t208 * t206;
t381 = t208 * t246;
t380 = t208 * t282;
t374 = t279 * t284;
t373 = t279 * t414;
t370 = t283 * t284;
t220 = pkin(4) * t265 + t402;
t369 = t284 * t220;
t362 = t360 * pkin(4) + t320;
t359 = t282 * t218 + t346 * t268;
t219 = pkin(3) * t377 + t268;
t358 = t284 * pkin(1) + t280 * pkin(6);
t273 = t279 ^ 2;
t274 = t283 ^ 2;
t357 = t273 - t274;
t356 = t273 + t274;
t354 = qJD(2) * t279;
t353 = qJD(2) * t283;
t287 = qJD(1) ^ 2;
t337 = t279 * t287 * t283;
t154 = pkin(3) * t298 + pkin(6) * t353;
t256 = pkin(2) + t401;
t221 = pkin(4) * t266 + t401;
t75 = t281 * t139 - t277 * t147;
t322 = t283 * t325;
t249 = t279 * t410;
t321 = -g(2) * t374 + t249;
t316 = pkin(6) * t206 + t229 * t278;
t315 = pkin(6) * t208 + t229 * t282;
t314 = -pkin(7) * t203 + qJD(3) * t229;
t313 = -t140 * t282 - t141 * t278;
t214 = pkin(2) + t221;
t272 = -pkin(9) - t414;
t311 = t283 * t214 - t279 * t272;
t309 = t283 * t256 + t373;
t174 = t209 * t279;
t56 = -t283 * pkin(4) + t174 * pkin(9) + t75;
t173 = t210 * t279;
t61 = -pkin(9) * t173 + t76;
t26 = -t276 * t61 + t413 * t56;
t27 = t276 * t56 + t413 * t61;
t304 = -pkin(6) * qJDD(2) + t343 * t415;
t103 = -t276 * t173 - t413 * t174;
t302 = t278 * t203 - t246 * t350;
t301 = t282 * t203 + t246 * t352;
t74 = t307 * qJD(2) + (-t248 + (pkin(8) * t279 - t222) * t278) * qJD(3) + t359;
t98 = t278 * t218 + t222 * t350 + (-t279 * t345 - t283 * t352) * pkin(6);
t81 = -pkin(8) * t298 + t98;
t24 = t139 * t348 - t147 * t349 + t277 * t74 + t281 * t81;
t300 = pkin(1) * t287 + t317;
t89 = t115 * pkin(3) + t180;
t286 = qJD(2) ^ 2;
t299 = pkin(6) * t286 + qJDD(1) * t415 + t405;
t295 = -t283 * t317 - t404;
t25 = -t76 * qJD(4) - t277 * t81 + t281 * t74;
t270 = t284 * pkin(6);
t190 = t280 * t278 + t283 * t363;
t188 = -t280 * t371 + t364;
t183 = qJDD(5) + t198;
t182 = pkin(3) * t335 + t276 * t255;
t181 = -pkin(3) * t379 + t413 * t255;
t165 = t209 * pkin(4) - t256;
t161 = -pkin(6) * t376 + t205;
t151 = -pkin(6) * t329 + t192;
t142 = t173 * pkin(4) + t219;
t136 = t413 * t209 + t276 * t210;
t102 = t413 * t173 - t276 * t174;
t99 = -t162 * qJD(3) + t359;
t92 = pkin(3) * t208 + pkin(4) * t312;
t83 = -t349 * t377 + (t341 * t375 + t334) * t281 + t423 * t277;
t82 = t145 * t279 + t277 * t334 - t281 * t333;
t62 = t83 * pkin(4) + t154;
t30 = pkin(4) * t306 + t89;
t29 = t103 * qJD(5) - t276 * t82 + t413 * t83;
t28 = t173 * t328 - t174 * t347 + t276 * t83 + t413 * t82;
t21 = -t83 * pkin(9) + t24;
t18 = pkin(4) * t354 + t82 * pkin(9) + t25;
t15 = t413 * t35 - t393;
t14 = -t276 * t35 - t338;
t4 = -t27 * qJD(5) + t413 * t18 - t276 * t21;
t3 = t26 * qJD(5) + t276 * t18 + t413 * t21;
t5 = [0, 0, 0, 0, 0, qJDD(1), -t405 + t410, t317, 0, 0, t273 * qJDD(1) + 0.2e1 * t322, 0.2e1 * t263 * t279 - 0.2e1 * t343 * t357, qJDD(2) * t279 + t286 * t283, t274 * qJDD(1) - 0.2e1 * t322, qJDD(2) * t283 - t286 * t279, 0, t304 * t279 + (-t299 + t410) * t283, t279 * t299 + t283 * t304 - t249, 0.2e1 * t356 * t390 - t317, -g(1) * (-t280 * pkin(1) + t270) - g(2) * t358 + (pkin(6) ^ 2 * t356 + pkin(1) ^ 2) * qJDD(1), -t114 * t375 + t423 * t208, (-t206 * t282 - t208 * t278) * t353 + (t389 - t388 + (-t380 + t383) * qJD(3)) * t279, (-t246 * t345 + t114) * t283 + (qJD(2) * t208 + t301) * t279, t115 * t377 + t206 * t298, (t246 * t346 + t115) * t283 + (-qJD(2) * t206 - t302) * t279, -t203 * t283 - t246 * t354, -g(1) * t188 - g(2) * t190 + t161 * t203 - t99 * t246 + (qJD(2) * t316 - t58) * t283 + (pkin(6) * t115 + qJD(2) * t140 + t180 * t278 + t229 * t350) * t279, -g(1) * t187 - g(2) * t189 - t162 * t203 + t98 * t246 + (qJD(2) * t315 + t57) * t283 + (-pkin(6) * t114 - qJD(2) * t141 + t180 * t282 - t229 * t352) * t279, t161 * t114 - t162 * t115 - t98 * t206 - t99 * t208 + t249 + t313 * t353 + (-t405 - t278 * t57 - t282 * t58 + (t140 * t278 - t141 * t282) * qJD(3)) * t279, t57 * t162 + t141 * t98 + t58 * t161 + t140 * t99 - g(1) * t270 - g(2) * (t284 * t319 + t358) - t222 * t410 + (t180 * t279 + t229 * t353) * pkin(6), t174 * t45 - t312 * t82, t82 * t130 + t45 * t173 + t174 * t306 - t312 * t83, -t198 * t174 + t236 * t82 + t283 * t45 + t312 * t354, t130 * t83 + t173 * t306, -t130 * t354 - t198 * t173 + t236 * t83 + t283 * t306, -t198 * t283 - t236 * t354, -g(1) * t169 - g(2) * t171 + t154 * t130 + t153 * t83 + t89 * t173 + t75 * t198 + t219 * t306 - t25 * t236 - t9 * t283 + t354 * t46, -g(1) * t168 - g(2) * t170 - t153 * t82 + t154 * t312 - t89 * t174 - t76 * t198 - t219 * t45 + t24 * t236 + t8 * t283 - t354 * t47, -t24 * t130 - t8 * t173 + t9 * t174 - t25 * t312 - t306 * t76 + t75 * t45 + t46 * t82 - t47 * t83 + t321, t8 * t76 + t47 * t24 + t9 * t75 + t46 * t25 + t89 * t219 + t153 * t154 - g(1) * (pkin(3) * t364 + t270) - g(2) * (t256 * t370 + t284 * t373 + t358) + (-g(1) * (-pkin(1) - t309) - g(2) * t402) * t280, -t10 * t103 - t28 * t431, t10 * t102 - t103 * t11 + t28 * t67 - t29 * t431, t10 * t283 + t103 * t183 + t28 * t228 + t354 * t431, t102 * t11 + t29 * t67, -t102 * t183 + t11 * t283 + t29 * t228 - t354 * t67, -t183 * t283 - t228 * t354, -g(1) * t157 - g(2) * t159 + t30 * t102 + t142 * t11 + t12 * t354 + t26 * t183 - t2 * t283 - t4 * t228 + t87 * t29 + t62 * t67, -g(1) * t156 - g(2) * t158 + t1 * t283 - t142 * t10 + t30 * t103 - t13 * t354 - t27 * t183 + t3 * t228 - t87 * t28 + t431 * t62, -t1 * t102 + t26 * t10 - t2 * t103 - t27 * t11 + t12 * t28 - t13 * t29 - t3 * t67 - t4 * t431 + t321, t1 * t27 + t13 * t3 + t2 * t26 + t12 * t4 + t30 * t142 + t87 * t62 - g(1) * (t270 + t369) - g(2) * (t214 * t370 - t272 * t374 + t358) + (-g(1) * (-pkin(1) - t311) - g(2) * t220) * t280; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t337, t357 * t287, t342, t337, t263, qJDD(2), t279 * t300 - t257 - t403, t404 + (t300 - t390) * t283, 0, 0, -t246 * t380 - t389, (-t114 + t384) * t282 + (-t115 + t381) * t278, (-t208 * t279 + t246 * t371) * qJD(1) + t302, -t246 * t383 - t388, (t206 * t279 - t246 * t376) * qJD(1) + t301, t246 * t355, -pkin(2) * t115 + t150 * t246 + t314 * t278 + (-t140 * t279 - t283 * t316) * qJD(1) + t445 * t282, pkin(2) * t114 - t151 * t246 + t314 * t282 + (t141 * t279 - t283 * t315) * qJD(1) - t445 * t278, t150 * t208 + t151 * t206 + ((qJD(3) * t208 - t115) * pkin(7) + t427) * t282 + (-t58 + t385 + (qJD(3) * t206 - t114) * pkin(7)) * t278 + t295, -t229 * t259 - t140 * t150 - t141 * t151 + t293 * pkin(2) + (qJD(3) * t313 - t58 * t278 + t57 * t282 + t295) * pkin(7), -t210 * t45 - t312 * t361, t130 * t361 + t45 * t209 - t210 * t306 - t312 * t360, t198 * t210 + t236 * t361 - t312 * t355, t130 * t360 + t209 * t306, t130 * t355 - t198 * t209 + t236 * t360, t236 * t355, t130 * t320 + t148 * t198 + t153 * t360 + t89 * t209 - t236 * t394 - t256 * t306 + t266 * t297 - t355 * t46, -t149 * t198 - t153 * t361 + t89 * t210 + t236 * t395 + t256 * t45 - t265 * t297 + t312 * t320 + t355 * t47, -t130 * t395 + t148 * t45 - t149 * t306 - t8 * t209 - t9 * t210 - t312 * t394 - t360 * t47 + t361 * t46 + t295, -g(3) * t309 + t9 * t148 + t8 * t149 + t153 * t320 - t89 * t256 + t394 * t46 + t395 * t47 + t317 * (t256 * t279 - t283 * t414), -t10 * t137 - t397 * t431, t10 * t136 - t11 * t137 - t396 * t431 + t397 * t67, t137 * t183 + t228 * t397 - t355 * t431, t11 * t136 + t396 * t67, -t136 * t183 + t228 * t396 + t355 * t67, t228 * t355, t165 * t11 - t12 * t355 + t30 * t136 + t59 * t183 - t228 * t398 + t254 * t297 + t362 * t67 + t396 * t87, -t165 * t10 + t13 * t355 + t30 * t137 - t60 * t183 + t228 * t399 - t253 * t297 + t362 * t431 - t397 * t87, -t1 * t136 + t59 * t10 - t60 * t11 + t12 * t397 - t13 * t396 - t2 * t137 - t398 * t431 - t399 * t67 + t295, -g(3) * t311 + t1 * t60 + t12 * t398 + t13 * t399 + t30 * t165 + t2 * t59 + t362 * t87 + t317 * (t214 * t279 + t272 * t283); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t382, -t206 ^ 2 + t208 ^ 2, -t114 - t384, -t382, -t115 - t381, t203, -t230 * t350 - t385 - t229 * t208 + t135 + (-qJD(3) * t199 - t179 + t404) * t278 + t422, g(1) * t190 - g(2) * t188 + g(3) * t375 + t229 * t206 - t427, 0, 0, t387, t436, t433, -t387, t417, t198, t52 * t236 + (-t130 * t208 + t198 * t281 + t236 * t349) * pkin(3) + t419, -t53 * t236 + (-t198 * t277 - t208 * t312 + t236 * t348) * pkin(3) + t432, t47 * t312 + t53 * t130 - t46 * t130 + t52 * t312 + (-t277 * t306 + t281 * t45 + (-t281 * t130 + t277 * t312) * qJD(4)) * pkin(3), -t46 * t52 - t47 * t53 + (t8 * t277 + t9 * t281 - t153 * t208 + g(3) * t377 + (-t46 * t277 + t47 * t281) * qJD(4) + t422) * pkin(3), t400, t437, t435, -t400, t418, t183, t181 * t183 - t228 * t391 - t92 * t67 + t420, -t182 * t183 + t228 * t392 - t431 * t92 + t434, t10 * t181 - t11 * t182 - t391 * t431 - t392 * t67 + t444, t1 * t182 + t2 * t181 - t87 * t92 - g(1) * (t280 * t221 - t283 * t369) - g(2) * (-t220 * t372 - t284 * t221) + t220 * t404 + t392 * t13 + t391 * t12; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t387, t436, t433, -t387, t417, t198, -t47 * t236 + t419, -t236 * t46 + t432, 0, 0, t400, t437, t435, -t400, t418, t183, t14 * t228 + (t413 * t183 + t228 * t347 - t312 * t67) * pkin(4) + t420, -t15 * t228 + (-t183 * t276 + t228 * t328 - t312 * t431) * pkin(4) + t434, t14 * t431 + t15 * t67 + (t413 * t10 - t11 * t276 + (t276 * t431 - t413 * t67) * qJD(5)) * pkin(4) + t444, -t12 * t14 - t13 * t15 + (t1 * t276 + t2 * t413 - t87 * t312 + (-t12 * t276 + t13 * t413) * qJD(5) + t421) * pkin(4); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t400, t437, t435, -t400, t418, t183, -t13 * t228 + t420, -t12 * t228 + t434, 0, 0;];
tau_reg = t5;
