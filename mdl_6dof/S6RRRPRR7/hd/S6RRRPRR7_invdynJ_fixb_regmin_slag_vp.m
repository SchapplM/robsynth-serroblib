% Calculate minimal parameter regressor of inverse dynamics joint torque vector for
% S6RRRPRR7
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
% pkin [12x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d1,d2,d3,d5,d6,theta4]';
% 
% Output:
% tau_reg [6x33]
%   minimal parameter regressor of inverse dynamics joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 18:43
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tau_reg = S6RRRPRR7_invdynJ_fixb_regmin_slag_vp(qJ, qJD, qJDD, g, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(3,1),zeros(12,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRPRR7_invdynJ_fixb_regmin_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRRPRR7_invdynJ_fixb_regmin_slag_vp: qJD has to be [6x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [6 1]), ...
  'S6RRRPRR7_invdynJ_fixb_regmin_slag_vp: qJDD has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RRRPRR7_invdynJ_fixb_regmin_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [12 1]), ...
  'S6RRRPRR7_invdynJ_fixb_regmin_slag_vp: pkin has to be [12x1] (double)');

%% Symbolic Calculation
% From invdyn_joint_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 18:39:31
% EndTime: 2019-03-09 18:40:03
% DurationCPUTime: 15.62s
% Computational Cost: add. (19651->639), mult. (48326->900), div. (0->0), fcn. (40381->16), ass. (0->336)
t320 = cos(qJ(6));
t416 = qJD(6) * t320;
t321 = cos(qJ(5));
t313 = cos(pkin(6));
t425 = qJD(1) * t313;
t295 = qJD(2) + t425;
t322 = cos(qJ(3));
t317 = sin(qJ(3));
t318 = sin(qJ(2));
t311 = sin(pkin(6));
t426 = qJD(1) * t311;
t399 = t318 * t426;
t377 = t317 * t399;
t215 = -t322 * t295 + t377;
t217 = t295 * t317 + t322 * t399;
t310 = sin(pkin(12));
t312 = cos(pkin(12));
t353 = t215 * t312 + t217 * t310;
t164 = t321 * t353;
t166 = t215 * t310 - t217 * t312;
t316 = sin(qJ(5));
t110 = -t166 * t316 + t164;
t522 = t110 * t320;
t529 = t416 + t522;
t323 = cos(qJ(2));
t369 = pkin(2) * t318 - pkin(9) * t323;
t239 = t369 * t426;
t221 = t322 * t239;
t408 = pkin(1) * t425;
t238 = -pkin(8) * t399 + t323 * t408;
t314 = -qJ(4) - pkin(9);
t388 = qJD(3) * t314;
t438 = t322 * t323;
t528 = -t221 - (pkin(3) * t318 - qJ(4) * t438) * t426 + t322 * t388 + (t238 - qJD(4)) * t317;
t424 = qJD(1) * t323;
t398 = t311 * t424;
t376 = t317 * t398;
t430 = t322 * t238 + t317 * t239;
t527 = -qJ(4) * t376 - qJD(4) * t322 - t317 * t388 + t430;
t410 = qJDD(1) * t323;
t293 = t311 * t410;
t413 = qJD(1) * qJD(2);
t393 = t318 * t413;
t375 = t311 * t393;
t235 = qJDD(3) - t293 + t375;
t224 = qJDD(5) + t235;
t277 = -qJD(3) + t398;
t265 = -qJD(5) + t277;
t315 = sin(qJ(6));
t417 = qJD(6) * t315;
t412 = qJDD(1) * t313;
t294 = qJDD(2) + t412;
t420 = qJD(3) * t322;
t392 = t323 * t413;
t411 = qJDD(1) * t318;
t516 = t392 + t411;
t501 = t516 * t311;
t159 = -qJD(3) * t377 + t317 * t294 + t295 * t420 + t322 * t501;
t422 = qJD(2) * t323;
t395 = t317 * t422;
t421 = qJD(3) * t317;
t160 = -t322 * t294 + t295 * t421 + t311 * (qJD(1) * (t318 * t420 + t395) + t317 * t411);
t104 = -t159 * t310 - t160 * t312;
t105 = t159 * t312 - t160 * t310;
t419 = qJD(5) * t316;
t43 = -qJD(5) * t164 + t316 * t104 + t321 * t105 + t166 * t419;
t497 = -t321 * t166 - t316 * t353;
t21 = t315 * t224 - t265 * t416 + t320 * t43 - t417 * t497;
t19 = t21 * t315;
t98 = -t265 * t315 + t320 * t497;
t526 = t529 * t98 + t19;
t44 = qJD(5) * t497 - t321 * t104 + t316 * t105;
t41 = qJDD(6) + t44;
t37 = t315 * t41;
t414 = -qJD(6) - t110;
t471 = t497 * t98;
t525 = -t529 * t414 + t37 - t471;
t504 = pkin(10) * t166;
t241 = pkin(8) * t398 + t318 * t408;
t204 = pkin(9) * t295 + t241;
t351 = -pkin(2) * t323 - pkin(9) * t318 - pkin(1);
t233 = t351 * t311;
t209 = qJD(1) * t233;
t157 = -t204 * t317 + t322 * t209;
t130 = -qJ(4) * t217 + t157;
t120 = -pkin(3) * t277 + t130;
t158 = t204 * t322 + t209 * t317;
t131 = -qJ(4) * t215 + t158;
t126 = t310 * t131;
t78 = t312 * t120 - t126;
t58 = -pkin(4) * t277 + t504 + t78;
t503 = pkin(10) * t353;
t445 = t312 * t131;
t79 = t310 * t120 + t445;
t62 = t79 - t503;
t29 = -t316 * t62 + t321 * t58;
t27 = pkin(5) * t265 - t29;
t524 = t110 * t27;
t463 = t110 * t265;
t523 = t43 - t463;
t521 = t497 * t110;
t258 = t310 * t322 + t312 * t317;
t197 = t258 * t398;
t247 = t258 * qJD(3);
t508 = t247 - t197;
t257 = -t310 * t317 + t312 * t322;
t198 = t257 * t398;
t248 = t257 * qJD(3);
t517 = -t248 + t198;
t520 = -t110 ^ 2 + t497 ^ 2;
t69 = pkin(5) * t497 + pkin(11) * t110;
t203 = -t295 * pkin(2) - t238;
t165 = t215 * pkin(3) + qJD(4) + t203;
t115 = pkin(4) * t353 + t165;
t319 = sin(qJ(1));
t440 = t319 * t323;
t324 = cos(qJ(1));
t441 = t318 * t324;
t253 = t313 * t441 + t440;
t306 = qJ(3) + pkin(12) + qJ(5);
t301 = sin(t306);
t302 = cos(t306);
t446 = t311 * t324;
t185 = t253 * t302 - t301 * t446;
t436 = t323 * t324;
t442 = t318 * t319;
t255 = -t313 * t442 + t436;
t449 = t311 * t319;
t188 = t255 * t302 + t301 * t449;
t450 = t311 * t318;
t211 = t301 * t313 + t302 * t450;
t407 = pkin(1) * qJD(2) * t313;
t379 = qJD(1) * t407;
t406 = pkin(1) * t412;
t402 = -pkin(8) * t293 - t318 * t406 - t323 * t379;
t335 = -pkin(8) * t375 - t402;
t172 = pkin(9) * t294 + t335;
t346 = t369 * qJD(2);
t174 = (qJD(1) * t346 + qJDD(1) * t351) * t311;
t330 = -qJD(3) * t158 - t317 * t172 + t322 * t174;
t55 = pkin(3) * t235 - qJ(4) * t159 - qJD(4) * t217 + t330;
t345 = -t322 * t172 - t317 * t174 + t204 * t421 - t209 * t420;
t60 = -qJ(4) * t160 - qJD(4) * t215 - t345;
t23 = -t310 * t60 + t312 * t55;
t14 = pkin(4) * t235 - pkin(10) * t105 + t23;
t24 = t310 * t55 + t312 * t60;
t16 = pkin(10) * t104 + t24;
t418 = qJD(5) * t321;
t387 = -t316 * t14 - t321 * t16 - t58 * t418 + t62 * t419;
t519 = g(1) * t188 + g(2) * t185 + g(3) * t211 + t110 * t115 + t387;
t22 = qJD(6) * t98 - t320 * t224 + t315 * t43;
t96 = t320 * t265 + t315 * t497;
t360 = t315 * t98 + t320 * t96;
t518 = -t110 * t360 + t21 * t320 - t315 * t22 - t96 * t416 - t417 * t98;
t433 = t527 * t310 + t528 * t312;
t493 = t528 * t310 - t527 * t312;
t472 = t497 * t96;
t512 = -pkin(4) * t399 + t517 * pkin(10) + t433;
t511 = -t508 * pkin(10) + t493;
t465 = t497 * t265;
t510 = -t44 - t465;
t509 = t414 * t497;
t30 = t316 * t58 + t321 * t62;
t28 = -pkin(11) * t265 + t30;
t49 = t110 * pkin(5) - pkin(11) * t497 + t115;
t11 = t28 * t320 + t315 * t49;
t187 = -t255 * t301 + t302 * t449;
t455 = t253 * t301;
t343 = -g(3) * (-t301 * t450 + t302 * t313) - g(2) * (-t302 * t446 - t455) - g(1) * t187;
t332 = -qJD(5) * t30 + t321 * t14 - t316 * t16;
t5 = -pkin(5) * t224 - t332;
t338 = t343 - t5;
t507 = t11 * t497 + t27 * t416 - t315 * t338;
t363 = t28 * t315 - t320 * t49;
t506 = t27 * t417 + t363 * t497;
t505 = -t115 * t497 + t332 + t343;
t352 = t321 * t257 - t258 * t316;
t136 = qJD(5) * t352 - t247 * t316 + t248 * t321;
t150 = -t197 * t316 + t198 * t321;
t435 = t136 - t150;
t190 = t257 * t316 + t258 * t321;
t434 = qJD(5) * t190 - t517 * t316 + t508 * t321;
t500 = -t241 + (-t376 + t421) * pkin(3);
t252 = -t313 * t436 + t442;
t499 = t185 * t315 - t252 * t320;
t498 = t185 * t320 + t252 * t315;
t307 = t311 ^ 2;
t409 = 0.2e1 * t307;
t278 = t314 * t317;
t279 = t314 * t322;
t201 = t312 * t278 + t279 * t310;
t178 = -pkin(10) * t258 + t201;
t202 = t310 * t278 - t312 * t279;
t179 = pkin(10) * t257 + t202;
t123 = t178 * t316 + t179 * t321;
t496 = qJD(5) * t123 + t511 * t316 - t512 * t321;
t38 = t320 * t41;
t495 = -t414 * t417 - t38;
t355 = t178 * t321 - t179 * t316;
t494 = qJD(5) * t355 + t512 * t316 + t511 * t321;
t432 = t508 * pkin(4) + t500;
t447 = t311 * t323;
t487 = pkin(1) * t318;
t428 = pkin(8) * t447 + t313 * t487;
t232 = pkin(9) * t313 + t428;
t431 = t322 * t232 + t317 * t233;
t303 = pkin(3) * t312 + pkin(4);
t484 = pkin(3) * t310;
t429 = t316 * t303 + t321 * t484;
t481 = g(1) * t324;
t492 = g(2) * t319 + t481;
t491 = pkin(3) * t160 + qJDD(4);
t448 = t311 * t322;
t194 = -t255 * t317 + t319 * t448;
t249 = -t313 * t322 + t317 * t450;
t437 = t322 * t324;
t454 = t253 * t317;
t490 = g(3) * t249 - g(2) * (-t311 * t437 - t454) - g(1) * t194;
t250 = t313 * t317 + t318 * t448;
t181 = -t249 * t310 + t250 * t312;
t383 = -t232 * t317 + t322 * t233;
t140 = -pkin(3) * t447 - qJ(4) * t250 + t383;
t147 = -qJ(4) * t249 + t431;
t91 = t312 * t140 - t147 * t310;
t74 = -pkin(4) * t447 - pkin(10) * t181 + t91;
t180 = -t249 * t312 - t250 * t310;
t92 = t310 * t140 + t312 * t147;
t76 = pkin(10) * t180 + t92;
t357 = t316 * t74 + t321 * t76;
t192 = qJD(3) * t250 + t311 * t395;
t396 = t311 * t422;
t193 = -qJD(3) * t249 + t322 * t396;
t146 = -t192 * t310 + t193 * t312;
t423 = qJD(2) * t318;
t397 = t311 * t423;
t240 = t311 * t346;
t296 = pkin(8) * t450;
t444 = t313 * t323;
t242 = (pkin(1) * t444 - t296) * qJD(2);
t331 = -qJD(3) * t431 + t322 * t240 - t242 * t317;
t89 = pkin(3) * t397 - qJ(4) * t193 - qJD(4) * t250 + t331;
t344 = -t232 * t421 + t233 * t420 + t317 * t240 + t322 * t242;
t93 = -qJ(4) * t192 - qJD(4) * t249 + t344;
t51 = -t310 * t93 + t312 * t89;
t36 = pkin(4) * t397 - pkin(10) * t146 + t51;
t145 = -t192 * t312 - t193 * t310;
t52 = t310 * t89 + t312 * t93;
t42 = pkin(10) * t145 + t52;
t488 = -qJD(5) * t357 - t316 * t42 + t321 * t36;
t4 = pkin(11) * t224 - t387;
t378 = pkin(8) * t501 + t318 * t379 - t323 * t406;
t486 = pkin(2) * t294;
t173 = t378 - t486;
t114 = t173 + t491;
t63 = -pkin(4) * t104 + t114;
t9 = pkin(5) * t44 - pkin(11) * t43 + t63;
t1 = -t363 * qJD(6) + t315 * t9 + t320 * t4;
t476 = g(3) * t311;
t474 = pkin(5) * t399 + t496;
t349 = t303 * t321 - t316 * t484;
t81 = -t130 * t310 - t445;
t67 = t81 + t503;
t82 = t312 * t130 - t126;
t68 = t82 + t504;
t469 = t349 * qJD(5) - t316 * t67 - t321 * t68;
t468 = qJD(5) * t429 - t316 * t68 + t321 * t67;
t466 = t414 * t315;
t461 = t190 * t315;
t460 = t190 * t320;
t459 = t215 * t277;
t458 = t217 * t277;
t453 = t302 * t315;
t452 = t302 * t320;
t451 = t307 * qJD(1) ^ 2;
t443 = t315 * t323;
t439 = t320 * t323;
t243 = pkin(8) * t396 + t318 * t407;
t308 = t318 ^ 2;
t427 = -t323 ^ 2 + t308;
t415 = qJD(2) - t295;
t405 = t323 * t451;
t404 = t311 * t443;
t403 = t311 * t439;
t304 = pkin(3) * t322 + pkin(2);
t133 = pkin(3) * t217 - pkin(4) * t166;
t237 = pkin(11) + t429;
t385 = qJD(6) * t237 + t133 + t69;
t382 = t253 * t322 - t317 * t446;
t381 = t295 + t425;
t380 = t294 + t412;
t373 = pkin(3) * t192 + t243;
t225 = -pkin(4) * t257 - t304;
t121 = -pkin(5) * t352 - pkin(11) * t190 + t225;
t372 = pkin(11) * t399 - qJD(6) * t121 - t494;
t371 = -pkin(5) * t434 + pkin(11) * t435 + qJD(6) * t123 - t432;
t254 = t313 * t440 + t441;
t368 = -g(1) * t252 + g(2) * t254;
t367 = g(1) * t255 + g(2) * t253;
t364 = -t237 * t41 + t524;
t34 = -pkin(11) * t447 + t357;
t125 = t180 * t316 + t181 * t321;
t231 = t296 + (-pkin(1) * t323 - pkin(2)) * t313;
t336 = pkin(3) * t249 + t231;
t132 = -pkin(4) * t180 + t336;
t354 = t321 * t180 - t181 * t316;
t56 = -pkin(5) * t354 - pkin(11) * t125 + t132;
t362 = t315 * t56 + t320 * t34;
t361 = -t315 * t34 + t320 * t56;
t358 = -t316 * t76 + t321 * t74;
t350 = t110 * t466 - t495;
t116 = t125 * t315 + t403;
t348 = t316 * t36 + t321 * t42 + t74 * t418 - t419 * t76;
t106 = -pkin(4) * t145 + t373;
t134 = t150 * t315 - t320 * t399;
t342 = t136 * t315 + t190 * t416 - t134;
t135 = t150 * t320 + t315 * t399;
t341 = t136 * t320 - t190 * t417 - t135;
t337 = -g(1) * t254 - g(2) * t252 + g(3) * t447;
t334 = -pkin(9) * t235 - t203 * t277;
t2 = -qJD(6) * t11 - t315 * t4 + t320 * t9;
t329 = -t337 - t378;
t328 = pkin(9) * qJD(3) * t277 - t173 - t337;
t236 = -pkin(5) - t349;
t195 = t255 * t322 + t317 * t449;
t154 = t188 * t320 + t254 * t315;
t153 = -t188 * t315 + t254 * t320;
t117 = t125 * t320 - t404;
t65 = qJD(5) * t125 - t321 * t145 + t146 * t316;
t64 = qJD(5) * t354 + t145 * t316 + t146 * t321;
t46 = -qJD(6) * t404 + t125 * t416 + t315 * t64 - t320 * t397;
t45 = -qJD(6) * t116 + t315 * t397 + t320 * t64;
t33 = pkin(5) * t447 - t358;
t17 = pkin(5) * t65 - pkin(11) * t64 + t106;
t7 = -pkin(5) * t397 - t488;
t6 = pkin(11) * t397 + t348;
t3 = [qJDD(1), g(1) * t319 - g(2) * t324, t492 (qJDD(1) * t308 + 0.2e1 * t318 * t392) * t307 (t318 * t410 - t413 * t427) * t409 (t318 * t380 + t381 * t422) * t311 (t323 * t380 - t381 * t423) * t311, t294 * t313, -t243 * t295 - t296 * t294 - t378 * t313 + g(1) * t253 - g(2) * t255 + (t294 * t444 + (-t393 + t410) * t409) * pkin(1), -t516 * pkin(1) * t409 - t242 * t295 - t428 * t294 - t335 * t313 + t368, t159 * t250 + t193 * t217, -t159 * t249 - t160 * t250 - t192 * t217 - t193 * t215, -t193 * t277 + t235 * t250 + (-t159 * t323 + t217 * t423) * t311, t192 * t277 - t235 * t249 + (t160 * t323 - t215 * t423) * t311 (-t235 * t323 - t277 * t423) * t311, -t331 * t277 + t383 * t235 + t243 * t215 + t231 * t160 + t173 * t249 + t203 * t192 + g(1) * t382 - g(2) * t195 + (t157 * t423 - t323 * t330) * t311, t344 * t277 - t431 * t235 + t243 * t217 + t231 * t159 + t173 * t250 + t203 * t193 - g(1) * t454 - g(2) * t194 + (-g(1) * t437 - t158 * t423 - t323 * t345) * t311, t92 * t104 - t91 * t105 + t79 * t145 - t78 * t146 + t166 * t51 + t24 * t180 - t23 * t181 - t353 * t52 - t368, t24 * t92 + t79 * t52 + t23 * t91 + t78 * t51 + t114 * t336 + t165 * t373 - g(1) * (-pkin(1) * t319 + t252 * t314 - t253 * t304) - g(2) * (pkin(1) * t324 - t254 * t314 + t255 * t304) - t492 * t311 * (pkin(3) * t317 + pkin(8)) t125 * t43 + t497 * t64, -t110 * t64 - t125 * t44 + t354 * t43 - t497 * t65, t125 * t224 - t265 * t64 + (-t323 * t43 + t423 * t497) * t311, t354 * t224 + t265 * t65 + (-t110 * t423 + t323 * t44) * t311 (-t224 * t323 - t265 * t423) * t311, -t488 * t265 + t358 * t224 + t106 * t110 + t132 * t44 - t63 * t354 + t115 * t65 + g(1) * t185 - g(2) * t188 + (t29 * t423 - t323 * t332) * t311, t348 * t265 - t357 * t224 + t106 * t497 + t132 * t43 + t63 * t125 + t115 * t64 - g(1) * t455 - g(2) * t187 + (-t30 * t423 - t302 * t481 - t323 * t387) * t311, t117 * t21 + t45 * t98, -t116 * t21 - t117 * t22 - t45 * t96 - t46 * t98, t117 * t41 - t21 * t354 - t414 * t45 + t65 * t98, -t116 * t41 + t22 * t354 + t414 * t46 - t65 * t96, -t354 * t41 - t414 * t65 -(-qJD(6) * t362 + t17 * t320 - t315 * t6) * t414 + t361 * t41 - t2 * t354 - t363 * t65 + t7 * t96 + t33 * t22 + t5 * t116 + t27 * t46 + g(1) * t498 - g(2) * t154 (qJD(6) * t361 + t17 * t315 + t320 * t6) * t414 - t362 * t41 + t1 * t354 - t11 * t65 + t7 * t98 + t33 * t21 + t5 * t117 + t27 * t45 - g(1) * t499 - g(2) * t153; 0, 0, 0, -t318 * t405, t427 * t451 (t415 * t424 + t411) * t311, -t399 * t415 + t293, t294, t241 * t295 + t451 * t487 + t329, pkin(1) * t405 + t238 * t295 + (pkin(8) * t413 + g(3)) * t450 + t367 + t402, t159 * t317 - t322 * t458 (t159 + t459) * t322 + (-t160 + t458) * t317, -t277 * t420 + t235 * t317 + (-t217 * t318 + t277 * t438) * t426, t277 * t421 + t235 * t322 + (-t277 * t317 * t323 + t215 * t318) * t426, t277 * t399, -t157 * t399 - pkin(2) * t160 - t241 * t215 + t221 * t277 + (-t238 * t277 + t334) * t317 + t328 * t322, -pkin(2) * t159 + t158 * t399 - t241 * t217 - t277 * t430 - t317 * t328 + t322 * t334, -g(3) * t450 + t202 * t104 - t201 * t105 + t433 * t166 - t23 * t258 + t24 * t257 - t493 * t353 - t508 * t79 + t517 * t78 - t367, t24 * t202 + t23 * t201 - t114 * t304 - g(1) * (-t254 * t304 - t255 * t314) - g(2) * (-t252 * t304 - t253 * t314) + t493 * t79 + t433 * t78 - (t304 * t323 - t314 * t318) * t476 + t500 * t165, t190 * t43 + t435 * t497, -t110 * t435 - t190 * t44 + t352 * t43 - t434 * t497, t190 * t224 - t265 * t435 - t399 * t497, t110 * t399 + t224 * t352 + t265 * t434, t265 * t399, t432 * t110 + t434 * t115 + t224 * t355 + t225 * t44 + t265 * t496 - t29 * t399 - t337 * t302 - t352 * t63, t435 * t115 - t123 * t224 + t63 * t190 + t225 * t43 + t265 * t494 + t30 * t399 + t337 * t301 + t432 * t497, t21 * t460 + t341 * t98, t134 * t98 + t135 * t96 - t360 * t136 + (-t19 - t22 * t320 + (t315 * t96 - t320 * t98) * qJD(6)) * t190, -t21 * t352 - t341 * t414 + t41 * t460 + t434 * t98, t22 * t352 + t342 * t414 - t41 * t461 - t434 * t96, -t352 * t41 - t414 * t434 (t121 * t320 - t123 * t315) * t41 - t2 * t352 - t355 * t22 + t5 * t461 - g(1) * (-t254 * t452 + t255 * t315) - g(2) * (-t252 * t452 + t253 * t315) + t474 * t96 - (t302 * t439 + t315 * t318) * t476 - (t315 * t372 - t320 * t371) * t414 - t434 * t363 + t342 * t27 -(t121 * t315 + t123 * t320) * t41 + t1 * t352 - t355 * t21 + t5 * t460 - g(1) * (t254 * t453 + t255 * t320) - g(2) * (t252 * t453 + t253 * t320) + t474 * t98 - (-t302 * t443 + t318 * t320) * t476 - t434 * t11 - (t315 * t371 + t320 * t372) * t414 + t341 * t27; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t217 * t215, -t215 ^ 2 + t217 ^ 2, t159 - t459, -t160 - t458, t235, -t158 * t277 - t203 * t217 + t330 + t490, g(1) * t195 + g(2) * t382 + g(3) * t250 - t157 * t277 + t203 * t215 + t345 (t104 * t310 - t105 * t312) * pkin(3) + (t82 - t78) * t353 + (-t79 - t81) * t166, -t78 * t81 - t79 * t82 + (-t165 * t217 + t23 * t312 + t24 * t310 + t490) * pkin(3), t521, t520, t523, t510, t224, -t133 * t110 + t224 * t349 + t265 * t468 + t505, -t133 * t497 - t224 * t429 + t265 * t469 + t519, t526, t518, t525, t350 + t472, t509, t236 * t22 + t468 * t96 + (t414 * t469 + t364) * t315 + (t385 * t414 + t338) * t320 + t506, t236 * t21 + t468 * t98 + t364 * t320 - (t315 * t385 - t320 * t469) * t414 + t507; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t166 ^ 2 - t353 ^ 2, -t166 * t78 + t353 * t79 - t329 - t486 + t491, 0, 0, 0, 0, 0, t44 - t465, t43 + t463, 0, 0, 0, 0, 0, t350 - t472, -t320 * t414 ^ 2 - t37 - t471; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t521, t520, t523, t510, t224, -t265 * t30 + t505, -t265 * t29 + t519, t526, t518, t525, -t414 * t466 + t38 + t472, t509, -pkin(5) * t22 - t30 * t96 + (-pkin(11) * t41 - t29 * t414 + t524) * t315 + (-(-pkin(11) * qJD(6) - t69) * t414 + t338) * t320 + t506, -pkin(5) * t21 - (t29 * t320 + t315 * t69) * t414 - t30 * t98 + t27 * t522 + t495 * pkin(11) + t507; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t98 * t96, -t96 ^ 2 + t98 ^ 2, -t414 * t96 + t21, -t414 * t98 - t22, t41, -t11 * t414 - t27 * t98 - g(1) * t153 + g(2) * t499 - g(3) * (-t211 * t315 - t403) + t2, t363 * t414 + t27 * t96 + g(1) * t154 + g(2) * t498 - g(3) * (-t211 * t320 + t404) - t1;];
tau_reg  = t3;
