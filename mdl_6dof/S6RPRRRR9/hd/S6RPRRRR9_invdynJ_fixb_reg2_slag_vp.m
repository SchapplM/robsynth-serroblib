% Calculate inertial parameters regressor of inverse dynamics joint torque vector for
% S6RPRRRR9
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
%   pkin=[a2,a3,a4,a5,a6,d1,d3,d4,d5,d6]';
% 
% Output:
% tau_reg [6x(6*10)]
%   inertial parameter regressor of inverse dynamics joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 07:26
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tau_reg = S6RPRRRR9_invdynJ_fixb_reg2_slag_vp(qJ, qJD, qJDD, g, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(3,1),zeros(10,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRRRR9_invdynJ_fixb_reg2_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RPRRRR9_invdynJ_fixb_reg2_slag_vp: qJD has to be [6x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [6 1]), ...
  'S6RPRRRR9_invdynJ_fixb_reg2_slag_vp: qJDD has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RPRRRR9_invdynJ_fixb_reg2_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RPRRRR9_invdynJ_fixb_reg2_slag_vp: pkin has to be [10x1] (double)');

%% Symbolic Calculation
% From invdyn_joint_fixb_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 07:25:26
% EndTime: 2019-03-09 07:25:50
% DurationCPUTime: 12.88s
% Computational Cost: add. (17180->722), mult. (34363->942), div. (0->0), fcn. (23315->14), ass. (0->313)
t297 = sin(qJ(3));
t301 = cos(qJ(3));
t340 = pkin(3) * t301 + pkin(8) * t297;
t232 = t340 * qJD(1);
t304 = -pkin(1) - pkin(7);
t258 = qJD(1) * t304 + qJD(2);
t300 = cos(qJ(4));
t296 = sin(qJ(4));
t415 = t296 * t301;
t156 = t300 * t232 - t258 * t415;
t454 = pkin(9) + pkin(8);
t360 = qJD(4) * t454;
t412 = t297 * t300;
t369 = pkin(9) * t412;
t485 = -(pkin(4) * t301 + t369) * qJD(1) - t156 - t300 * t360;
t406 = t300 * t301;
t157 = t296 * t232 + t258 * t406;
t390 = qJD(1) * t297;
t358 = t296 * t390;
t484 = pkin(9) * t358 + t296 * t360 + t157;
t375 = t300 * qJD(3);
t389 = qJD(1) * t301;
t222 = t296 * t389 - t375;
t386 = qJD(3) * t296;
t224 = t300 * t389 + t386;
t295 = sin(qJ(5));
t299 = cos(qJ(5));
t145 = t222 * t299 + t224 * t295;
t294 = sin(qJ(6));
t333 = -t222 * t295 + t299 * t224;
t453 = cos(qJ(6));
t469 = -t294 * t145 + t333 * t453;
t72 = t453 * t145 + t294 * t333;
t443 = t72 * t469;
t370 = qJD(4) + qJD(5);
t377 = qJD(5) * t299;
t380 = qJD(4) * t300;
t408 = t299 * t300;
t417 = t295 * t296;
t399 = t295 * t358 - t299 * t380 - t300 * t377 + t370 * t417 - t390 * t408;
t227 = t295 * t300 + t296 * t299;
t155 = t370 * t227;
t200 = t227 * qJD(1);
t398 = t297 * t200 + t155;
t263 = qJD(4) + t390;
t254 = qJDD(1) * t304 + qJDD(2);
t385 = qJD(3) * t297;
t160 = -qJDD(3) * pkin(3) - t301 * t254 + t258 * t385;
t445 = g(3) * t297;
t302 = cos(qJ(1));
t287 = g(2) * t302;
t298 = sin(qJ(1));
t288 = g(1) * t298;
t462 = t288 - t287;
t316 = t301 * t462 - t445;
t313 = -t160 - t316;
t483 = -pkin(8) * qJD(4) * t263 + t313;
t475 = t469 ^ 2 - t72 ^ 2;
t256 = qJD(5) + t263;
t339 = pkin(3) * t297 - pkin(8) * t301;
t238 = qJ(2) + t339;
t201 = t238 * qJD(1);
t235 = t297 * t258;
t211 = qJD(3) * pkin(8) + t235;
t125 = t300 * t201 - t211 * t296;
t101 = -pkin(9) * t224 + t125;
t92 = pkin(4) * t263 + t101;
t126 = t201 * t296 + t211 * t300;
t102 = -pkin(9) * t222 + t126;
t96 = t295 * t102;
t44 = t299 * t92 - t96;
t466 = pkin(10) * t333;
t38 = t44 - t466;
t35 = pkin(5) * t256 + t38;
t98 = t299 * t102;
t45 = t295 * t92 + t98;
t479 = pkin(10) * t145;
t39 = t45 - t479;
t438 = t294 * t39;
t10 = t35 * t453 - t438;
t365 = t453 * t39;
t11 = t294 * t35 + t365;
t482 = -t10 * t72 + t11 * t469;
t356 = t297 * t375;
t379 = qJD(4) * t301;
t319 = t296 * t379 + t356;
t371 = t301 * qJDD(1);
t130 = qJD(1) * t319 - qJD(4) * t375 - t296 * qJDD(3) - t300 * t371;
t357 = t296 * t385;
t382 = qJD(4) * t224;
t131 = -qJD(1) * t357 - t300 * qJDD(3) + t296 * t371 + t382;
t378 = qJD(5) * t295;
t328 = -t295 * t130 + t131 * t299 - t222 * t378 + t224 * t377;
t352 = qJD(6) * t453;
t376 = qJD(6) * t294;
t50 = t299 * t130 + t295 * t131 + t222 * t377 + t224 * t378;
t14 = t145 * t352 + t294 * t328 + t333 * t376 + t453 * t50;
t248 = qJD(6) + t256;
t473 = t248 * t72 - t14;
t374 = qJD(1) * qJD(3);
t351 = t301 * t374;
t372 = t297 * qJDD(1);
t219 = qJDD(4) + t351 + t372;
t214 = qJDD(5) + t219;
t220 = qJD(3) * t340 + qJD(2);
t143 = qJD(1) * t220 + qJDD(1) * t238;
t133 = t300 * t143;
t384 = qJD(3) * t301;
t161 = qJDD(3) * pkin(8) + t254 * t297 + t258 * t384;
t55 = -qJD(4) * t126 - t296 * t161 + t133;
t36 = pkin(4) * t219 + pkin(9) * t130 + t55;
t381 = qJD(4) * t296;
t54 = t296 * t143 + t300 * t161 + t201 * t380 - t211 * t381;
t41 = -pkin(9) * t131 + t54;
t9 = -qJD(5) * t45 - t295 * t41 + t299 * t36;
t6 = pkin(5) * t214 + pkin(10) * t50 + t9;
t8 = -t102 * t378 + t295 * t36 + t299 * t41 + t92 * t377;
t7 = -pkin(10) * t328 + t8;
t1 = t294 * t6 + t35 * t352 - t39 * t376 + t453 * t7;
t293 = qJ(4) + qJ(5);
t284 = qJ(6) + t293;
t271 = sin(t284);
t272 = cos(t284);
t413 = t297 * t298;
t170 = t271 * t302 + t272 * t413;
t411 = t297 * t302;
t172 = -t271 * t298 + t272 * t411;
t444 = g(3) * t301;
t212 = -qJD(3) * pkin(3) - t258 * t301;
t159 = pkin(4) * t222 + t212;
t83 = pkin(5) * t145 + t159;
t472 = g(1) * t170 - g(2) * t172 + t272 * t444 + t83 * t72 - t1;
t249 = t454 * t296;
t250 = t454 * t300;
t440 = -t249 * t377 - t250 * t378 + t295 * t485 - t484 * t299;
t163 = -t295 * t249 + t299 * t250;
t439 = -qJD(5) * t163 + t484 * t295 + t299 * t485;
t15 = qJD(6) * t469 - t294 * t50 + t453 * t328;
t457 = t248 * t469 - t15;
t169 = -t271 * t413 + t272 * t302;
t171 = t271 * t411 + t272 * t298;
t2 = -qJD(6) * t11 - t294 * t7 + t453 * t6;
t459 = -g(1) * t169 - g(2) * t171 + t271 * t444 - t469 * t83 + t2;
t478 = -pkin(5) * t389 + pkin(10) * t399 + t439;
t477 = pkin(10) * t398 - t440;
t184 = t227 * t297;
t226 = -t408 + t417;
t187 = t226 * t301;
t431 = qJD(3) * t187 + t184 * t370 + t200;
t186 = t226 * t297;
t430 = t226 * qJD(1) + t186 * t370 - t227 * t384;
t424 = t145 * t333;
t476 = t300 * t379 - t357;
t474 = -t145 ^ 2 + t333 ^ 2;
t471 = t145 * t256 - t50;
t280 = sin(t293);
t281 = cos(t293);
t180 = t280 * t302 + t281 * t413;
t182 = -t280 * t298 + t281 * t411;
t470 = g(1) * t180 - g(2) * t182 + t145 * t159 + t281 * t444 - t8;
t465 = -t125 * t263 + t54;
t291 = t297 ^ 2;
t292 = t301 ^ 2;
t392 = t291 + t292;
t347 = t392 * t254;
t218 = t300 * t238;
t350 = -t296 * t304 + pkin(4);
t139 = -pkin(9) * t406 + t297 * t350 + t218;
t410 = t297 * t304;
t255 = t300 * t410;
t167 = t296 * t238 + t255;
t158 = -pkin(9) * t415 + t167;
t79 = t295 * t139 + t299 * t158;
t341 = -t235 + (t358 + t381) * pkin(4);
t405 = t300 * t302;
t202 = -t296 * t413 + t405;
t409 = t298 * t300;
t204 = t296 * t411 + t409;
t461 = -g(1) * t202 - g(2) * t204;
t179 = -t280 * t413 + t281 * t302;
t181 = t280 * t411 + t281 * t298;
t460 = -g(1) * t179 - g(2) * t181 + t280 * t444;
t458 = -t333 * t159 + t460 + t9;
t456 = t256 * t333 - t328;
t455 = 0.2e1 * qJ(2);
t452 = pkin(4) * t296;
t451 = pkin(4) * t300;
t446 = g(2) * t298;
t162 = -t299 * t249 - t250 * t295;
t118 = -pkin(10) * t227 + t162;
t119 = -pkin(10) * t226 + t163;
t64 = t118 * t453 - t294 * t119;
t442 = qJD(6) * t64 + t294 * t478 - t453 * t477;
t65 = t294 * t118 + t119 * t453;
t441 = -qJD(6) * t65 + t294 * t477 + t453 * t478;
t53 = t299 * t101 - t96;
t152 = -t294 * t226 + t227 * t453;
t437 = -qJD(6) * t152 + t294 * t399 - t398 * t453;
t436 = t226 * t352 + t227 * t376 + t294 * t398 + t399 * t453;
t111 = -t294 * t184 - t186 * t453;
t435 = qJD(6) * t111 - t294 * t431 - t430 * t453;
t109 = -t184 * t453 + t294 * t186;
t434 = -qJD(6) * t109 - t294 * t430 + t431 * t453;
t273 = pkin(4) * t299 + pkin(5);
t418 = t294 * t295;
t52 = -t101 * t295 - t98;
t42 = t52 + t479;
t43 = t53 - t466;
t433 = -t294 * t42 - t43 * t453 + t273 * t352 + (-t295 * t376 + (t299 * t453 - t418) * qJD(5)) * pkin(4);
t359 = t453 * t295;
t432 = t294 * t43 - t42 * t453 - t273 * t376 + (-t295 * t352 + (-t294 * t299 - t359) * qJD(5)) * pkin(4);
t429 = pkin(1) * qJDD(1);
t427 = t126 * t263;
t426 = t130 * t296;
t425 = t131 * t300;
t423 = t222 * t263;
t422 = t224 * t222;
t421 = t224 * t263;
t420 = t224 * t300;
t419 = t263 * t296;
t416 = t296 * t219;
t414 = t296 * t302;
t407 = t300 * t219;
t404 = t301 * t302;
t403 = t301 * t454;
t305 = qJD(3) ^ 2;
t402 = t304 * t305;
t306 = qJD(1) ^ 2;
t401 = t306 * qJ(2);
t400 = pkin(5) * t398 + t341;
t397 = g(1) * t404 + t301 * t446;
t366 = 0.2e1 * qJD(1) * qJD(2);
t396 = (qJDD(1) * qJ(2) + t366) * qJ(2);
t395 = t302 * pkin(1) + t298 * qJ(2);
t393 = t291 - t292;
t391 = -t305 - t306;
t388 = qJD(3) * t222;
t387 = qJD(3) * t224;
t383 = qJD(3) * t304;
t373 = qJDD(3) * t297;
t364 = t296 * t410;
t363 = t301 * t306 * t297;
t361 = t302 * pkin(7) + t395;
t274 = pkin(3) + t451;
t355 = t301 * t383;
t237 = pkin(5) * t281 + t451;
t349 = -g(2) * t411 + t444;
t78 = t299 * t139 - t158 * t295;
t221 = pkin(4) * t415 - t301 * t304;
t346 = t392 * qJDD(1);
t345 = qJDD(2) - t429;
t344 = qJD(4) * t297 + qJD(1);
t343 = g(2) * t361;
t342 = t297 * t351;
t338 = g(1) * t302 + t446;
t336 = -t462 - t401;
t335 = t125 * t300 + t126 * t296;
t334 = t125 * t296 - t126 * t300;
t231 = pkin(3) + t237;
t290 = -pkin(10) - t454;
t331 = t231 * t297 + t290 * t301;
t329 = t274 * t297 - t403;
t327 = qJDD(1) * t455 + t366;
t326 = -t254 + t401 + t288;
t58 = pkin(5) * t297 + pkin(10) * t187 + t78;
t185 = t227 * t301;
t61 = -pkin(10) * t185 + t79;
t26 = -t294 * t61 + t453 * t58;
t27 = t294 * t58 + t453 * t61;
t112 = -t294 * t185 - t187 * t453;
t323 = t263 * t380 + t416;
t322 = -t263 * t381 + t407;
t164 = pkin(4) * t476 + t297 * t383;
t195 = t300 * t220;
t69 = t195 + (-t255 + (pkin(9) * t301 - t238) * t296) * qJD(4) + (t301 * t350 + t369) * qJD(3);
t99 = -qJD(4) * t364 + t296 * t220 + t238 * t380 + t300 * t355;
t80 = -pkin(9) * t476 + t99;
t24 = t139 * t377 - t158 * t378 + t295 * t69 + t299 * t80;
t321 = -g(1) * t413 - t349;
t320 = qJDD(3) * t304 + t374 * t455;
t317 = -pkin(8) * t219 + t212 * t263;
t81 = t131 * pkin(4) + t160;
t312 = t327 - t338;
t25 = -qJD(5) * t79 - t295 * t80 + t299 * t69;
t310 = -qJD(4) * t335 - t55 * t296 + t54 * t300;
t283 = t302 * qJ(2);
t278 = qJDD(3) * t301;
t236 = pkin(5) * t280 + t452;
t205 = -t296 * t298 + t297 * t405;
t203 = t297 * t409 + t414;
t196 = qJDD(6) + t214;
t193 = pkin(4) * t359 + t294 * t273;
t192 = -pkin(4) * t418 + t273 * t453;
t177 = pkin(5) * t226 - t274;
t166 = t218 - t364;
t151 = t226 * t453 + t227 * t294;
t150 = pkin(5) * t185 + t221;
t110 = t185 * t453 - t187 * t294;
t103 = pkin(4) * t224 + pkin(5) * t333;
t100 = -qJD(4) * t167 - t296 * t355 + t195;
t91 = -t378 * t415 + (t370 * t406 - t357) * t299 - t319 * t295;
t89 = t155 * t301 - t295 * t357 + t299 * t356;
t66 = pkin(5) * t91 + t164;
t32 = qJD(6) * t112 - t294 * t89 + t453 * t91;
t30 = t185 * t352 - t187 * t376 + t294 * t91 + t453 * t89;
t28 = pkin(5) * t328 + t81;
t19 = -pkin(10) * t91 + t24;
t18 = pkin(5) * t384 + pkin(10) * t89 + t25;
t13 = t38 * t453 - t438;
t12 = -t294 * t38 - t365;
t4 = -qJD(6) * t27 + t18 * t453 - t294 * t19;
t3 = qJD(6) * t26 + t294 * t18 + t19 * t453;
t5 = [0, 0, 0, 0, 0, qJDD(1), t462, t338, 0, 0, qJDD(1), 0, 0, 0, 0, 0, 0, qJDD(2) - t462 - 0.2e1 * t429, t312, -t345 * pkin(1) - g(1) * (-t298 * pkin(1) + t283) - g(2) * t395 + t396, qJDD(1) * t292 - 0.2e1 * t342, -0.2e1 * t297 * t371 + 0.2e1 * t374 * t393, -t297 * t305 + t278, qJDD(1) * t291 + 0.2e1 * t342, -t301 * t305 - t373, 0, t320 * t301 + (t312 - t402) * t297, -t320 * t297 + (t327 - t402) * t301 - t397, -t304 * t346 - t347 + t462, -g(1) * (t298 * t304 + t283) - t343 + t304 * t347 + t396, -t130 * t406 - t224 * t319 (t222 * t300 + t224 * t296) * t385 + (t426 - t425 + (t222 * t296 - t420) * qJD(4)) * t301 (-t263 * t375 - t130) * t297 + (t322 + t387) * t301, t131 * t415 + t222 * t476 (t263 * t386 - t131) * t297 + (-t323 - t388) * t301, t219 * t297 + t263 * t384, -g(1) * t205 - g(2) * t203 + t100 * t263 + t166 * t219 + (t55 + (-t212 * t296 + t222 * t304) * qJD(3)) * t297 + (qJD(3) * t125 - t131 * t304 + t160 * t296 + t212 * t380) * t301, g(1) * t204 - g(2) * t202 - t167 * t219 - t263 * t99 + (-t54 + (-t212 * t300 + t224 * t304) * qJD(3)) * t297 + (-qJD(3) * t126 + t130 * t304 + t160 * t300 - t212 * t381) * t301, -t100 * t224 + t130 * t166 - t131 * t167 - t222 * t99 + t335 * t385 + (qJD(4) * t334 - t296 * t54 - t300 * t55) * t301 + t397, t54 * t167 + t126 * t99 + t55 * t166 + t125 * t100 - g(1) * (pkin(3) * t411 - pkin(8) * t404 + t283) - t343 + (-t160 * t301 + t212 * t385) * t304 + (-g(1) * t304 - g(2) * t339) * t298, t187 * t50 - t333 * t89, t89 * t145 + t50 * t185 + t187 * t328 - t333 * t91, -t187 * t214 - t256 * t89 - t297 * t50 + t333 * t384, t145 * t91 + t185 * t328, -t145 * t384 - t185 * t214 - t91 * t256 - t297 * t328, t214 * t297 + t256 * t384, -g(1) * t182 - g(2) * t180 + t164 * t145 + t159 * t91 + t81 * t185 + t78 * t214 + t221 * t328 + t25 * t256 + t9 * t297 + t384 * t44, g(1) * t181 - g(2) * t179 - t159 * t89 + t164 * t333 - t187 * t81 - t214 * t79 - t221 * t50 - t24 * t256 - t297 * t8 - t384 * t45, -t24 * t145 - t8 * t185 + t9 * t187 - t25 * t333 - t328 * t79 + t44 * t89 - t45 * t91 + t78 * t50 + t397, t8 * t79 + t45 * t24 + t9 * t78 + t44 * t25 + t81 * t221 + t159 * t164 - g(1) * (t274 * t411 - t302 * t403 + t283) - g(2) * (pkin(4) * t414 + t361) + (-g(1) * (-t452 + t304) - g(2) * t329) * t298, -t112 * t14 - t30 * t469, t110 * t14 - t112 * t15 + t30 * t72 - t32 * t469, t112 * t196 - t14 * t297 - t248 * t30 + t384 * t469, t110 * t15 + t32 * t72, -t110 * t196 - t15 * t297 - t248 * t32 - t384 * t72, t196 * t297 + t248 * t384, -g(1) * t172 - g(2) * t170 + t10 * t384 + t110 * t28 + t15 * t150 + t196 * t26 + t2 * t297 + t248 * t4 + t32 * t83 + t66 * t72, g(1) * t171 - g(2) * t169 - t1 * t297 - t11 * t384 + t112 * t28 - t14 * t150 - t196 * t27 - t248 * t3 - t30 * t83 + t469 * t66, -t1 * t110 + t10 * t30 - t11 * t32 - t112 * t2 + t14 * t26 - t15 * t27 - t3 * t72 - t4 * t469 + t397, t1 * t27 + t11 * t3 + t2 * t26 + t10 * t4 + t28 * t150 + t83 * t66 - g(1) * (t231 * t411 + t290 * t404 + t283) - g(2) * (t236 * t302 + t361) + (-g(1) * (-t236 + t304) - g(2) * t331) * t298; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, qJDD(1), -t306, t336 + t345, 0, 0, 0, 0, 0, 0, t297 * t391 + t278, t301 * t391 - t373, -t346, t347 + t336, 0, 0, 0, 0, 0, 0, -t131 * t301 + (t388 - t416) * t297 + (-t296 * t384 - t300 * t344) * t263, t130 * t301 + (t387 - t407) * t297 + (t296 * t344 - t301 * t375) * t263 (-t131 * t297 - t222 * t384 + t224 * t344) * t300 + (-t130 * t297 + t222 * t344 + t224 * t384) * t296, -t335 * qJD(1) + (-qJD(3) * t334 - t160) * t301 + (qJD(3) * t212 + t310) * t297 - t462, 0, 0, 0, 0, 0, 0, t145 * t385 - t184 * t214 + t256 * t430 - t301 * t328, t186 * t214 + t256 * t431 + t301 * t50 + t333 * t385, t145 * t431 - t184 * t50 + t186 * t328 - t333 * t430, t159 * t385 - t184 * t9 - t186 * t8 - t301 * t81 + t430 * t44 - t431 * t45 - t462, 0, 0, 0, 0, 0, 0, t109 * t196 - t301 * t15 - t248 * t435 + t385 * t72, -t111 * t196 + t301 * t14 + t248 * t434 + t385 * t469, t109 * t14 - t111 * t15 + t434 * t72 + t435 * t469, t1 * t111 - t10 * t435 + t109 * t2 - t11 * t434 - t28 * t301 + t385 * t83 - t462; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t363, -t393 * t306, t371, -t363, -t372, qJDD(3), t445 + (-t326 + t287) * t301, t297 * t326 + t349, 0, 0, t263 * t420 - t426 (-t130 - t423) * t300 + (-t131 - t421) * t296 (-t224 * t301 + t263 * t412) * qJD(1) + t323, t222 * t419 - t425 (t222 * t301 - t297 * t419) * qJD(1) + t322, -t263 * t389, -pkin(3) * t131 - t125 * t389 - t156 * t263 - t222 * t235 + t317 * t296 + t300 * t483, pkin(3) * t130 + t126 * t389 + t157 * t263 - t224 * t235 - t296 * t483 + t317 * t300, t156 * t224 + t157 * t222 + ((-t131 + t382) * pkin(8) + t465) * t300 + (-t55 - t427 + (qJD(4) * t222 - t130) * pkin(8)) * t296 + t321, -t212 * t235 - t125 * t156 - t126 * t157 + t313 * pkin(3) + (-t297 * t462 + t310 - t444) * pkin(8), -t50 * t227 - t333 * t399, t145 * t399 + t50 * t226 - t227 * t328 - t333 * t398, t227 * t214 - t256 * t399 - t333 * t389, t145 * t398 + t226 * t328, t145 * t389 - t226 * t214 - t256 * t398, -t256 * t389, t145 * t341 + t159 * t398 + t162 * t214 + t81 * t226 + t256 * t439 - t274 * t328 - t281 * t316 - t389 * t44, -t159 * t399 - t163 * t214 + t227 * t81 - t256 * t440 + t274 * t50 + t280 * t316 + t333 * t341 + t389 * t45, -t145 * t440 + t162 * t50 - t163 * t328 - t8 * t226 - t9 * t227 - t333 * t439 - t398 * t45 + t399 * t44 + t321, g(3) * t329 + t159 * t341 + t9 * t162 + t8 * t163 - t81 * t274 + t439 * t44 + t440 * t45 - t462 * (t274 * t301 + t297 * t454) -t14 * t152 - t436 * t469, t14 * t151 - t15 * t152 + t436 * t72 + t437 * t469, t152 * t196 - t248 * t436 - t389 * t469, t15 * t151 - t437 * t72, -t151 * t196 + t248 * t437 + t389 * t72, -t248 * t389, -t10 * t389 + t15 * t177 + t151 * t28 + t196 * t64 + t248 * t441 - t272 * t316 + t400 * t72 - t437 * t83, t11 * t389 - t14 * t177 + t152 * t28 - t196 * t65 - t248 * t442 + t271 * t316 + t400 * t469 - t436 * t83, -t1 * t151 + t10 * t436 + t11 * t437 + t14 * t64 - t15 * t65 - t152 * t2 - t441 * t469 - t442 * t72 + t321, g(3) * t331 + t1 * t65 + t10 * t441 + t11 * t442 + t28 * t177 + t2 * t64 + t400 * t83 - t462 * (t231 * t301 - t290 * t297); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t422, -t222 ^ 2 + t224 ^ 2, -t130 + t423, -t422, t421 - t131, t219, -t211 * t380 + t427 - t212 * t224 + t133 + (-qJD(4) * t201 - t161 + t444) * t296 + t461, g(1) * t203 - g(2) * t205 + g(3) * t406 + t212 * t222 - t465, 0, 0, t424, t474, t471, -t424, t456, t214, -t256 * t52 + (-t145 * t224 + t214 * t299 - t256 * t378) * pkin(4) + t458, t256 * t53 + (-t214 * t295 - t224 * t333 - t256 * t377) * pkin(4) + t470, t45 * t333 + t53 * t145 - t44 * t145 + t52 * t333 + (-t295 * t328 + t299 * t50 + (-t145 * t299 + t295 * t333) * qJD(5)) * pkin(4), -t44 * t52 - t45 * t53 + (t8 * t295 + t9 * t299 - t159 * t224 + g(3) * t415 + (-t295 * t44 + t299 * t45) * qJD(5) + t461) * pkin(4), t443, t475, t473, -t443, t457, t196, -t103 * t72 + t192 * t196 + t248 * t432 + t459, -t103 * t469 - t193 * t196 - t248 * t433 + t472, t14 * t192 - t15 * t193 - t432 * t469 - t433 * t72 + t482, t1 * t193 + t2 * t192 - t83 * t103 - g(1) * (-t236 * t413 + t237 * t302) - g(2) * (t236 * t411 + t237 * t298) + t236 * t444 + t433 * t11 + t432 * t10; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t424, t474, t471, -t424, t456, t214, t256 * t45 + t458, t256 * t44 + t470, 0, 0, t443, t475, t473, -t443, t457, t196, -t12 * t248 + (t196 * t453 - t248 * t376 - t333 * t72) * pkin(5) + t459, t13 * t248 + (-t196 * t294 - t248 * t352 - t333 * t469) * pkin(5) + t472, t12 * t469 + t13 * t72 + (t453 * t14 - t15 * t294 + (t294 * t469 - t453 * t72) * qJD(6)) * pkin(5) + t482, -t10 * t12 - t11 * t13 + (t1 * t294 + t2 * t453 - t83 * t333 + (-t10 * t294 + t11 * t453) * qJD(6) + t460) * pkin(5); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t443, t475, t473, -t443, t457, t196, t11 * t248 + t459, t10 * t248 + t472, 0, 0;];
tau_reg  = t5;
