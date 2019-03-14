% Calculate minimal parameter regressor of inverse dynamics joint torque vector for
% S6RRPRPR9
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
%   pkin=[a2,a3,a4,a5,a6,alpha2,d1,d2,d4,d6,theta3,theta5]';
% 
% Output:
% tau_reg [6x32]
%   minimal parameter regressor of inverse dynamics joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 11:03
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tau_reg = S6RRPRPR9_invdynJ_fixb_regmin_slag_vp(qJ, qJD, qJDD, g, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(3,1),zeros(12,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPRPR9_invdynJ_fixb_regmin_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRPRPR9_invdynJ_fixb_regmin_slag_vp: qJD has to be [6x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [6 1]), ...
  'S6RRPRPR9_invdynJ_fixb_regmin_slag_vp: qJDD has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RRPRPR9_invdynJ_fixb_regmin_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [12 1]), ...
  'S6RRPRPR9_invdynJ_fixb_regmin_slag_vp: pkin has to be [12x1] (double)');

%% Symbolic Calculation
% From invdyn_joint_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 11:01:19
% EndTime: 2019-03-09 11:01:44
% DurationCPUTime: 13.26s
% Computational Cost: add. (16687->697), mult. (41731->960), div. (0->0), fcn. (34698->18), ass. (0->318)
t323 = cos(pkin(6));
t423 = qJD(1) * t323;
t295 = qJD(2) + t423;
t319 = sin(pkin(11));
t322 = cos(pkin(11));
t327 = sin(qJ(2));
t320 = sin(pkin(6));
t424 = qJD(1) * t320;
t398 = t327 * t424;
t230 = t295 * t322 - t319 * t398;
t231 = t295 * t319 + t322 * t398;
t326 = sin(qJ(4));
t330 = cos(qJ(4));
t172 = -t330 * t230 + t231 * t326;
t318 = sin(pkin(12));
t321 = cos(pkin(12));
t325 = sin(qJ(6));
t329 = cos(qJ(6));
t261 = t318 * t329 + t321 * t325;
t251 = t261 * qJD(6);
t454 = t261 * t172 + t251;
t331 = cos(qJ(2));
t412 = qJD(1) * qJD(2);
t393 = t331 * t412;
t410 = qJDD(1) * t327;
t350 = t393 + t410;
t484 = t350 * t320;
t260 = t319 * t326 - t330 * t322;
t438 = t320 * t331;
t343 = t260 * t438;
t218 = qJD(1) * t343;
t252 = t260 * qJD(4);
t496 = t252 - t218;
t262 = t319 * t330 + t322 * t326;
t344 = t262 * t438;
t429 = -qJD(1) * t344 + t262 * qJD(4);
t367 = pkin(2) * t327 - qJ(3) * t331;
t244 = t367 * t424;
t406 = pkin(1) * t423;
t245 = -pkin(8) * t398 + t331 * t406;
t184 = t322 * t244 - t245 * t319;
t345 = (-pkin(9) * t322 * t331 + pkin(3) * t327) * t320;
t148 = qJD(1) * t345 + t184;
t185 = t319 * t244 + t322 * t245;
t422 = qJD(1) * t331;
t397 = t320 * t422;
t381 = t319 * t397;
t162 = -pkin(9) * t381 + t185;
t465 = pkin(9) + qJ(3);
t273 = t465 * t319;
t275 = t465 * t322;
t361 = -t273 * t330 - t275 * t326;
t495 = -qJD(3) * t260 + qJD(4) * t361 - t326 * t148 - t330 * t162;
t168 = qJD(6) + t172;
t259 = t318 * t325 - t329 * t321;
t483 = t168 * t259;
t411 = qJDD(1) * t323;
t294 = qJDD(2) + t411;
t446 = t294 * t319;
t199 = t322 * t484 + t446;
t362 = t230 * t326 + t231 * t330;
t428 = t484 * t319;
t368 = t294 * t322 - t428;
t85 = qJD(4) * t362 + t199 * t326 - t330 * t368;
t83 = qJDD(6) + t85;
t494 = -t168 * t483 + t261 * t83;
t493 = -qJ(5) * t398 + t495;
t246 = pkin(8) * t397 + t327 * t406;
t210 = pkin(3) * t381 + t246;
t492 = pkin(4) * t429 + qJ(5) * t496 - qJD(5) * t262 - t210;
t491 = -t168 * t454 - t259 * t83;
t279 = -qJD(4) + t397;
t139 = t321 * t279 + t318 * t362;
t141 = -t279 * t318 + t321 * t362;
t75 = t329 * t139 + t141 * t325;
t490 = t168 * t75;
t363 = t139 * t325 - t141 * t329;
t488 = t168 * t363;
t487 = t172 * t279;
t216 = -t273 * t326 + t275 * t330;
t486 = qJD(3) * t262 + qJD(4) * t216 + t148 * t330 - t326 * t162;
t186 = -t218 * t318 - t321 * t398;
t485 = -t252 * t318 - t186;
t472 = cos(qJ(1));
t400 = t472 * t327;
t328 = sin(qJ(1));
t436 = t328 * t331;
t255 = t323 * t400 + t436;
t315 = pkin(11) + qJ(4);
t309 = sin(t315);
t311 = cos(t315);
t401 = t320 * t472;
t201 = t255 * t311 - t309 * t401;
t399 = t472 * t331;
t437 = t327 * t328;
t254 = -t323 * t399 + t437;
t314 = pkin(12) + qJ(6);
t308 = sin(t314);
t310 = cos(t314);
t482 = t201 * t308 - t254 * t310;
t481 = t201 * t310 + t254 * t308;
t313 = t320 ^ 2;
t408 = 0.2e1 * t313;
t480 = pkin(3) * t319;
t463 = -t318 * t493 + t492 * t321;
t479 = t362 * t75;
t478 = t362 * t363;
t459 = t492 * t318 + t321 * t493;
t456 = pkin(4) * t398 + t486;
t477 = t279 * t362;
t256 = t323 * t436 + t400;
t374 = g(1) * t256 + g(2) * t254;
t475 = g(3) * t438 - t374;
t476 = t475 * t309;
t466 = t294 * pkin(2);
t474 = qJDD(3) - t466;
t409 = qJDD(1) * t331;
t293 = t320 * t409;
t394 = t327 * t412;
t379 = t320 * t394;
t243 = qJDD(4) - t293 + t379;
t473 = -pkin(4) * t243 + qJDD(5);
t417 = qJD(4) * t330;
t418 = qJD(4) * t326;
t84 = t330 * t199 + t230 * t417 - t231 * t418 + t326 * t368;
t69 = -t321 * t243 + t318 * t84;
t70 = t243 * t318 + t321 * t84;
t20 = -qJD(6) * t363 + t325 * t70 + t329 * t69;
t471 = pkin(1) * t327;
t467 = g(3) * t320;
t464 = pkin(10) + qJ(5);
t219 = qJ(3) * t295 + t246;
t359 = -pkin(2) * t331 - qJ(3) * t327 - pkin(1);
t240 = t359 * t320;
t224 = qJD(1) * t240;
t152 = -t219 * t319 + t322 * t224;
t107 = -pkin(3) * t397 - pkin(9) * t231 + t152;
t153 = t322 * t219 + t319 * t224;
t121 = pkin(9) * t230 + t153;
t349 = t394 - t409;
t405 = pkin(1) * qJD(2) * t323;
t383 = qJD(1) * t405;
t404 = pkin(1) * t411;
t402 = -pkin(8) * t293 - t327 * t404 - t331 * t383;
t338 = -pkin(8) * t379 - t402;
t159 = qJ(3) * t294 + qJD(3) * t295 + t338;
t334 = qJD(2) * t367 - qJD(3) * t327;
t167 = (qJD(1) * t334 + qJDD(1) * t359) * t320;
t90 = -t159 * t319 + t322 * t167;
t60 = pkin(3) * t320 * t349 - pkin(9) * t199 + t90;
t91 = t322 * t159 + t319 * t167;
t72 = pkin(9) * t368 + t91;
t351 = t107 * t417 - t121 * t418 + t326 * t60 + t330 * t72;
t14 = qJ(5) * t243 - qJD(5) * t279 + t351;
t382 = pkin(8) * t484 + t327 * t383 - t331 * t404;
t177 = t382 + t474;
t118 = -pkin(3) * t368 + t177;
t25 = t85 * pkin(4) - t84 * qJ(5) - qJD(5) * t362 + t118;
t7 = t321 * t14 + t318 * t25;
t222 = t334 * t320;
t420 = qJD(2) * t327;
t396 = t320 * t420;
t360 = -pkin(8) * t396 + t331 * t405;
t228 = qJD(3) * t323 + t360;
t160 = t322 * t222 - t228 * t319;
t125 = qJD(2) * t345 + t160;
t426 = pkin(8) * t438 + t323 * t471;
t239 = qJ(3) * t323 + t426;
t178 = -t239 * t319 + t322 * t240;
t440 = t320 * t327;
t249 = t319 * t323 + t322 * t440;
t128 = -pkin(3) * t438 - pkin(9) * t249 + t178;
t161 = t319 * t222 + t322 * t228;
t419 = qJD(2) * t331;
t395 = t320 * t419;
t380 = t319 * t395;
t144 = -pkin(9) * t380 + t161;
t179 = t322 * t239 + t319 * t240;
t248 = t319 * t440 - t323 * t322;
t147 = -pkin(9) * t248 + t179;
t348 = t326 * t125 + t128 * t417 + t330 * t144 - t147 * t418;
t33 = (qJ(5) * t420 - qJD(5) * t331) * t320 + t348;
t188 = t330 * t248 + t249 * t326;
t135 = -qJD(2) * t343 - qJD(4) * t188;
t189 = -t248 * t326 + t249 * t330;
t136 = qJD(2) * t344 + qJD(4) * t189;
t247 = pkin(8) * t395 + t327 * t405;
t211 = pkin(3) * t380 + t247;
t49 = pkin(4) * t136 - qJ(5) * t135 - qJD(5) * t189 + t211;
t17 = t318 * t49 + t321 * t33;
t56 = t326 * t107 + t330 * t121;
t53 = -qJ(5) * t279 + t56;
t212 = -pkin(2) * t295 + qJD(3) - t245;
t169 = -pkin(3) * t230 + t212;
t66 = pkin(4) * t172 - qJ(5) * t362 + t169;
t29 = t318 * t66 + t321 * t53;
t55 = t107 * t330 - t326 * t121;
t98 = pkin(4) * t362 + qJ(5) * t172;
t37 = t318 * t98 + t321 * t55;
t432 = t326 * t128 + t330 * t147;
t73 = -qJ(5) * t438 + t432;
t296 = pkin(8) * t440;
t242 = t296 + (-pkin(1) * t331 - pkin(2)) * t323;
t197 = pkin(3) * t248 + t242;
t89 = pkin(4) * t188 - qJ(5) * t189 + t197;
t41 = t318 * t89 + t321 * t73;
t461 = t318 * t85;
t460 = t321 * t85;
t457 = pkin(5) * t485 + t456;
t453 = t172 * t318;
t452 = t172 * t321;
t448 = t262 * t318;
t447 = t262 * t321;
t445 = t294 * t323;
t444 = t308 * t311;
t443 = t310 * t311;
t442 = t311 * t331;
t441 = t313 * qJD(1) ^ 2;
t439 = t320 * t328;
t52 = pkin(4) * t279 + qJD(5) - t55;
t435 = -qJD(5) + t52;
t187 = -t218 * t321 + t318 * t398;
t434 = t186 * t325 - t187 * t329 - t251 * t262 + t252 * t259;
t415 = qJD(6) * t329;
t416 = qJD(6) * t325;
t433 = -t329 * t186 - t187 * t325 - t252 * t261 + t415 * t447 - t416 * t448;
t306 = pkin(3) * t322 + pkin(2);
t198 = pkin(4) * t260 - qJ(5) * t262 - t306;
t131 = t318 * t198 + t321 * t216;
t427 = pkin(1) * t472 + pkin(8) * t439;
t316 = t327 ^ 2;
t425 = -t331 ^ 2 + t316;
t421 = qJD(2) * t322;
t414 = qJD(2) - t295;
t413 = qJ(3) * qJDD(1);
t407 = g(3) * t440;
t403 = t331 * t441;
t391 = -t328 * pkin(1) + pkin(8) * t401;
t6 = -t14 * t318 + t321 * t25;
t16 = -t318 * t33 + t321 * t49;
t28 = -t318 * t53 + t321 * t66;
t36 = -t318 * t55 + t321 * t98;
t40 = -t318 * t73 + t321 * t89;
t389 = t128 * t330 - t326 * t147;
t130 = t321 * t198 - t216 * t318;
t386 = -t107 * t418 - t121 * t417 - t326 * t72 + t330 * t60;
t385 = t295 + t423;
t384 = t294 + t411;
t101 = pkin(5) * t260 - pkin(10) * t447 + t130;
t377 = pkin(10) * t485 - qJD(6) * t101 - t459;
t110 = -pkin(10) * t448 + t131;
t376 = qJD(6) * t110 - t463 + (-t252 * t321 - t187) * pkin(10) - t429 * pkin(5);
t200 = t255 * t309 + t311 * t401;
t257 = -t323 * t437 + t399;
t204 = t257 * t309 - t311 * t439;
t375 = -g(1) * t200 + g(2) * t204;
t373 = -g(1) * t254 + g(2) * t256;
t372 = g(1) * t257 + g(2) * t255;
t371 = g(1) * t255 - g(2) * t257;
t370 = -t318 * t7 - t321 * t6;
t4 = pkin(5) * t85 - pkin(10) * t70 + t6;
t5 = -pkin(10) * t69 + t7;
t369 = t325 * t4 + t329 * t5;
t74 = pkin(4) * t438 - t389;
t18 = pkin(5) * t172 - pkin(10) * t141 + t28;
t21 = -pkin(10) * t139 + t29;
t8 = t18 * t329 - t21 * t325;
t9 = t18 * t325 + t21 * t329;
t166 = t189 * t321 - t318 * t438;
t27 = pkin(5) * t188 - pkin(10) * t166 + t40;
t165 = t189 * t318 + t321 * t438;
t31 = -pkin(10) * t165 + t41;
t366 = t27 * t329 - t31 * t325;
t365 = t27 * t325 + t31 * t329;
t364 = t28 * t318 - t29 * t321;
t94 = t329 * t165 + t166 * t325;
t95 = -t165 * t325 + t166 * t329;
t358 = -t177 + t374;
t356 = t125 * t330 - t128 * t418 - t326 * t144 - t147 * t417;
t355 = g(1) * t472 + g(2) * t328;
t274 = t464 * t321;
t353 = pkin(5) * t362 + pkin(10) * t452 + qJD(5) * t318 + qJD(6) * t274 + t36;
t272 = t464 * t318;
t352 = pkin(10) * t453 - qJD(5) * t321 + qJD(6) * t272 + t37;
t19 = -t139 * t415 - t141 * t416 - t325 * t69 + t329 * t70;
t347 = g(1) * t204 + g(2) * t200 + g(3) * (t309 * t440 - t323 * t311);
t205 = t257 * t311 + t309 * t439;
t238 = t309 * t323 + t311 * t440;
t346 = -g(1) * t205 - g(2) * t201 - g(3) * t238;
t15 = -t386 + t473;
t342 = -t15 + t347;
t340 = -t372 - t407;
t339 = -qJ(3) * t420 + (qJD(3) - t212) * t331;
t337 = t475 * t311;
t336 = t15 * t262 - t252 * t52 - t372;
t335 = -t475 - t382;
t34 = -pkin(4) * t396 - t356;
t2 = -qJD(6) * t9 - t325 * t5 + t329 * t4;
t333 = t347 + t386;
t307 = -pkin(5) * t321 - pkin(4);
t191 = t259 * t262;
t190 = t261 * t262;
t183 = pkin(5) * t448 - t361;
t150 = t205 * t310 + t256 * t308;
t149 = -t205 * t308 + t256 * t310;
t120 = t135 * t321 + t318 * t396;
t119 = t135 * t318 - t321 * t396;
t50 = pkin(5) * t165 + t74;
t44 = -pkin(5) * t453 + t56;
t43 = pkin(5) * t139 + t52;
t39 = qJD(6) * t95 + t329 * t119 + t120 * t325;
t38 = -qJD(6) * t94 - t119 * t325 + t120 * t329;
t26 = pkin(5) * t119 + t34;
t12 = -pkin(10) * t119 + t17;
t11 = pkin(5) * t136 - pkin(10) * t120 + t16;
t10 = pkin(5) * t69 + t15;
t1 = t8 * qJD(6) + t369;
t3 = [qJDD(1), g(1) * t328 - g(2) * t472, t355 (qJDD(1) * t316 + 0.2e1 * t327 * t393) * t313 (t327 * t409 - t412 * t425) * t408 (t327 * t384 + t385 * t419) * t320 (t331 * t384 - t385 * t420) * t320, t445, -t247 * t295 - t296 * t294 - t382 * t323 + (t331 * t445 - t349 * t408) * pkin(1) + t371, -pkin(1) * t350 * t408 - t294 * t426 - t295 * t360 - t323 * t338 + t373, -t247 * t230 + t242 * t428 + t177 * t248 + (-t242 * t294 + t371) * t322 + (-t355 * t319 + (qJD(1) * t178 + t152) * t420 + (qJD(2) * t212 * t319 - qJD(1) * t160 - qJDD(1) * t178 - t90) * t331) * t320, t177 * t249 + t242 * t199 + t247 * t231 - t371 * t319 + (-t355 * t322 + (-qJD(1) * t179 - t153) * t420 + (qJD(1) * t161 + qJDD(1) * t179 + t212 * t421 + t91) * t331) * t320, t161 * t230 + t179 * t368 - t91 * t248 - t160 * t231 - t178 * t199 - t90 * t249 + (-t152 * t322 - t153 * t319) * t395 - t373, t91 * t179 + t153 * t161 + t90 * t178 + t152 * t160 + t177 * t242 + t212 * t247 - g(1) * (-pkin(2) * t255 - qJ(3) * t254 + t391) - g(2) * (pkin(2) * t257 + qJ(3) * t256 + t427) t135 * t362 + t189 * t84, -t135 * t172 - t136 * t362 - t188 * t84 - t189 * t85, -t135 * t279 + t189 * t243 + (-t331 * t84 + t362 * t420) * t320, t136 * t279 - t188 * t243 + (-t172 * t420 + t331 * t85) * t320 (-t243 * t331 - t279 * t420) * t320, -t356 * t279 + t389 * t243 + t211 * t172 + t197 * t85 + t118 * t188 + t169 * t136 + g(1) * t201 - g(2) * t205 + (-t331 * t386 + t420 * t55) * t320, t348 * t279 - t432 * t243 + t211 * t362 + t197 * t84 + t118 * t189 + t169 * t135 + (t331 * t351 - t420 * t56) * t320 + t375, t16 * t172 + t40 * t85 + t6 * t188 + t28 * t136 + t34 * t139 + t74 * t69 + t15 * t165 + t52 * t119 - g(1) * (-t201 * t321 - t254 * t318) - g(2) * (t205 * t321 + t256 * t318) -t17 * t172 - t41 * t85 - t7 * t188 - t29 * t136 + t34 * t141 + t74 * t70 + t15 * t166 + t52 * t120 - g(1) * (t201 * t318 - t254 * t321) - g(2) * (-t205 * t318 + t256 * t321) -t119 * t29 - t120 * t28 - t139 * t17 - t141 * t16 - t165 * t7 - t166 * t6 - t40 * t70 - t41 * t69 - t375, t7 * t41 + t29 * t17 + t6 * t40 + t28 * t16 + t15 * t74 + t52 * t34 - g(1) * (-pkin(4) * t201 - qJ(5) * t200 - t254 * t465 - t255 * t306 + t401 * t480 + t391) - g(2) * (pkin(4) * t205 + qJ(5) * t204 + t256 * t465 + t257 * t306 + t439 * t480 + t427) t19 * t95 - t363 * t38, -t19 * t94 - t20 * t95 + t363 * t39 - t38 * t75, -t136 * t363 + t168 * t38 + t188 * t19 + t83 * t95, -t136 * t75 - t168 * t39 - t188 * t20 - t83 * t94, t136 * t168 + t188 * t83 (-qJD(6) * t365 + t11 * t329 - t12 * t325) * t168 + t366 * t83 + t2 * t188 + t8 * t136 + t26 * t75 + t50 * t20 + t10 * t94 + t43 * t39 + g(1) * t481 - g(2) * t150 -(qJD(6) * t366 + t11 * t325 + t12 * t329) * t168 - t365 * t83 - t1 * t188 - t9 * t136 - t26 * t363 + t50 * t19 + t10 * t95 + t43 * t38 - g(1) * t482 - g(2) * t149; 0, 0, 0, -t327 * t403, t425 * t441 (t414 * t422 + t410) * t320, -t398 * t414 + t293, t294, t246 * t295 + t441 * t471 + t335, pkin(1) * t403 + t245 * t295 + (pkin(8) * t412 + g(3)) * t440 + t372 + t402, -pkin(2) * t428 + t246 * t230 + (t358 + t466) * t322 + ((-g(3) * t322 + t319 * t413) * t331 + (-t152 * t327 + t184 * t331 + t319 * t339) * qJD(1)) * t320, -pkin(2) * t199 - t231 * t246 - t358 * t319 + ((g(3) * t319 + t322 * t413) * t331 + (t153 * t327 - t185 * t331 + t322 * t339) * qJD(1)) * t320, t184 * t231 - t185 * t230 + (qJ(3) * t368 + qJD(3) * t230 + t152 * t397 + t91) * t322 + (qJ(3) * t199 + qJD(3) * t231 + t153 * t397 - t90) * t319 + t340, -t152 * t184 - t153 * t185 - t212 * t246 + (-t152 * t319 + t153 * t322) * qJD(3) + (-t177 - t475) * pkin(2) + (-t90 * t319 + t91 * t322 + t340) * qJ(3), t262 * t84 - t362 * t496, t172 * t496 - t260 * t84 - t262 * t85 - t362 * t429, t243 * t262 + t279 * t496 - t362 * t398, t172 * t398 - t243 * t260 + t279 * t429, t279 * t398, t118 * t260 + t429 * t169 - t210 * t172 + t361 * t243 + t279 * t486 - t306 * t85 - t55 * t398 - t337, t118 * t262 - t169 * t496 - t210 * t362 - t216 * t243 + t279 * t495 - t306 * t84 + t56 * t398 + t476, t130 * t85 - t52 * t186 - t361 * t69 + t6 * t260 - t321 * t337 + (t336 - t407) * t318 + t429 * t28 + t463 * t172 + t456 * t139, -t131 * t85 - t52 * t187 - t361 * t70 - t7 * t260 - t374 * t318 * t311 + t336 * t321 - (-t318 * t442 + t321 * t327) * t467 - t429 * t29 - t459 * t172 + t456 * t141, -t130 * t70 - t131 * t69 + t186 * t29 + t187 * t28 + t370 * t262 - (-t28 * t321 - t29 * t318) * t252 - t463 * t141 - t459 * t139 - t476, t6 * t130 + t7 * t131 - t15 * t361 + t463 * t28 + t459 * t29 + t456 * t52 - (t327 * t467 + t372) * t465 + (-t331 * t467 + t374) * (pkin(4) * t311 + qJ(5) * t309 + t306) -t19 * t191 - t363 * t434, -t19 * t190 + t191 * t20 + t363 * t433 - t434 * t75, t168 * t434 + t19 * t260 - t191 * t83 - t363 * t429, -t168 * t433 - t190 * t83 - t20 * t260 - t429 * t75, t168 * t429 + t260 * t83 (t101 * t329 - t110 * t325) * t83 + t2 * t260 + t183 * t20 + t10 * t190 - g(1) * (-t256 * t443 + t257 * t308) - g(2) * (-t254 * t443 + t255 * t308) + t429 * t8 + t457 * t75 + t433 * t43 - (t308 * t327 + t310 * t442) * t467 + (t325 * t377 - t329 * t376) * t168 -(t101 * t325 + t110 * t329) * t83 - t1 * t260 + t183 * t19 - t10 * t191 - g(1) * (t256 * t444 + t257 * t310) - g(2) * (t254 * t444 + t255 * t310) - t429 * t9 - t457 * t363 + t434 * t43 - (-t308 * t442 + t310 * t327) * t467 + (t325 * t376 + t329 * t377) * t168; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t231 * t397 - t368, t446 + (t322 * t410 + (-t230 + t421) * t422) * t320, -t230 ^ 2 - t231 ^ 2, t152 * t231 - t153 * t230 - t335 + t474, 0, 0, 0, 0, 0, t85 - t477, t84 + t487, -t139 * t362 - t172 * t453 + t460, -t141 * t362 - t172 * t452 - t461, -t318 * t69 - t321 * t70 - (t139 * t321 - t141 * t318) * t172, -t172 * t364 - t362 * t52 - t370 + t475, 0, 0, 0, 0, 0, -t479 + t491, t478 - t494; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t362 * t172, -t172 ^ 2 + t362 ^ 2, t84 - t487, -t85 - t477, t243, -t169 * t362 - t279 * t56 + t333, t169 * t172 - t279 * t55 - t346 - t351, -qJ(5) * t461 - pkin(4) * t69 - t139 * t56 - t362 * t28 + (t318 * t435 - t36) * t172 + t342 * t321, -qJ(5) * t460 - pkin(4) * t70 - t141 * t56 + t362 * t29 + (t321 * t435 + t37) * t172 - t342 * t318, t139 * t37 + t141 * t36 + (-qJ(5) * t69 - qJD(5) * t139 - t172 * t28 + t7) * t321 + (qJ(5) * t70 + qJD(5) * t141 - t172 * t29 - t6) * t318 + t346, -t28 * t36 - t29 * t37 - t52 * t56 - t364 * qJD(5) + t342 * pkin(4) + (-t318 * t6 + t321 * t7 + t346) * qJ(5), t19 * t261 + t363 * t483, -t19 * t259 - t20 * t261 + t363 * t454 + t483 * t75, t478 + t494, t479 + t491, -t168 * t362 (-t272 * t329 - t274 * t325) * t83 + t307 * t20 + t10 * t259 - t8 * t362 - t44 * t75 + t454 * t43 + (t325 * t352 - t329 * t353) * t168 + t347 * t310 -(-t272 * t325 + t274 * t329) * t83 + t307 * t19 + t10 * t261 + t9 * t362 + t44 * t363 - t483 * t43 + (t325 * t353 + t329 * t352) * t168 - t347 * t308; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t141 * t172 + t69, -t139 * t172 + t70, -t139 ^ 2 - t141 ^ 2, t139 * t29 + t141 * t28 - t333 + t473, 0, 0, 0, 0, 0, t20 - t488, t19 - t490; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t363 * t75, t363 ^ 2 - t75 ^ 2, t19 + t490, -t20 - t488, t83, t9 * t168 + t43 * t363 - g(1) * t149 + g(2) * t482 - g(3) * (-t238 * t308 - t310 * t438) + t2, t43 * t75 + g(1) * t150 + g(2) * t481 - g(3) * (-t238 * t310 + t308 * t438) - t369 + (t168 - qJD(6)) * t8;];
tau_reg  = t3;