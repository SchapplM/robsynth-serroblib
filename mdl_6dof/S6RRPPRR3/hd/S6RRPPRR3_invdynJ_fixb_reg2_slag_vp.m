% Calculate inertial parameters regressor of inverse dynamics joint torque vector for
% S6RRPPRR3
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
%   pkin=[a2,a3,a4,a5,a6,alpha2,d1,d2,d5,d6,theta3,theta4]';
% 
% Output:
% tau_reg [6x(6*10)]
%   inertial parameter regressor of inverse dynamics joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 09:00
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tau_reg = S6RRPPRR3_invdynJ_fixb_reg2_slag_vp(qJ, qJD, qJDD, g, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(3,1),zeros(12,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPPRR3_invdynJ_fixb_reg2_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRPPRR3_invdynJ_fixb_reg2_slag_vp: qJD has to be [6x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [6 1]), ...
  'S6RRPPRR3_invdynJ_fixb_reg2_slag_vp: qJDD has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RRPPRR3_invdynJ_fixb_reg2_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [12 1]), ...
  'S6RRPPRR3_invdynJ_fixb_reg2_slag_vp: pkin has to be [12x1] (double)');

%% Symbolic Calculation
% From invdyn_joint_fixb_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 08:58:47
% EndTime: 2019-03-09 08:59:22
% DurationCPUTime: 22.07s
% Computational Cost: add. (27386->817), mult. (77875->1077), div. (0->0), fcn. (64995->16), ass. (0->386)
t331 = sin(qJ(6));
t335 = cos(qJ(6));
t327 = sin(pkin(6));
t336 = cos(qJ(2));
t510 = cos(pkin(11));
t420 = t336 * t510;
t397 = t327 * t420;
t289 = qJD(1) * t397;
t326 = sin(pkin(11));
t333 = sin(qJ(2));
t461 = qJD(1) * t333;
t434 = t327 * t461;
t250 = t326 * t434 - t289;
t363 = qJD(5) + t250;
t278 = -t336 * t326 - t333 * t510;
t463 = qJD(1) * t327;
t254 = t278 * t463;
t329 = cos(pkin(6));
t462 = qJD(1) * t329;
t306 = qJD(2) + t462;
t325 = sin(pkin(12));
t328 = cos(pkin(12));
t210 = t254 * t328 - t306 * t325;
t332 = sin(qJ(5));
t534 = cos(qJ(5));
t483 = t325 * t254;
t554 = t306 * t328 + t483;
t561 = -t210 * t534 + t332 * t554;
t104 = t331 * t363 + t335 * t561;
t192 = t534 * t554;
t125 = t332 * t210 + t192;
t569 = qJD(6) - t125;
t508 = t569 * t331;
t574 = t104 * t508;
t310 = t329 * t336 * pkin(1);
t302 = qJD(1) * t310;
t524 = pkin(8) + qJ(3);
t430 = t524 * t333;
t403 = t327 * t430;
t231 = -qJD(1) * t403 + t302;
t214 = pkin(2) * t306 + t231;
t478 = t329 * t333;
t309 = pkin(1) * t478;
t480 = t327 * t336;
t232 = (t480 * t524 + t309) * qJD(1);
t419 = t510 * t232;
t137 = t326 * t214 + t419;
t129 = qJ(4) * t306 + t137;
t318 = pkin(2) * t336 + pkin(1);
t264 = -t318 * t463 + qJD(3);
t151 = pkin(3) * t250 + qJ(4) * t254 + t264;
t84 = -t129 * t325 + t328 * t151;
t62 = pkin(4) * t250 + pkin(9) * t210 + t84;
t85 = t328 * t129 + t325 * t151;
t66 = pkin(9) * t554 + t85;
t29 = t332 * t62 + t534 * t66;
t26 = pkin(10) * t363 + t29;
t222 = t326 * t232;
t136 = t214 * t510 - t222;
t128 = -t306 * pkin(3) + qJD(4) - t136;
t100 = -pkin(4) * t554 + t128;
t48 = -pkin(5) * t125 - pkin(10) * t561 + t100;
t14 = t26 * t335 + t331 * t48;
t395 = qJD(2) * t420;
t450 = qJD(1) * qJD(2);
t427 = t333 * t450;
t402 = t327 * t427;
t186 = -t326 * t402 + (qJD(1) * t395 - qJDD(1) * t278) * t327;
t449 = qJDD(1) * t329;
t305 = qJDD(2) + t449;
t154 = t186 * t325 - t328 * t305;
t532 = pkin(4) * t154;
t442 = pkin(1) * t449;
t301 = t336 * t442;
t428 = t329 * t450;
t405 = pkin(1) * t428;
t371 = -t333 * t405 + t301;
t422 = qJD(2) * t524;
t459 = qJD(3) * t333;
t135 = pkin(2) * t305 + (-qJDD(1) * t430 + (-t336 * t422 - t459) * qJD(1)) * t327 + t371;
t358 = qJD(3) * t336 - t333 * t422;
t308 = pkin(8) * t480;
t436 = qJDD(1) * t308 + t333 * t442 + t336 * t405;
t447 = qJDD(1) * t336;
t144 = (qJ(3) * t447 + qJD(1) * t358) * t327 + t436;
t88 = t510 * t135 - t326 * t144;
t79 = -pkin(3) * t305 + qJDD(4) - t88;
t53 = t79 + t532;
t155 = t186 * t328 + t305 * t325;
t458 = qJD(5) * t332;
t60 = -qJD(5) * t192 + t332 * t154 - t534 * t155 - t210 * t458;
t415 = t534 * t154 + t332 * t155;
t564 = qJD(5) * t561;
t61 = t415 + t564;
t16 = pkin(5) * t61 + pkin(10) * t60 + t53;
t286 = qJDD(1) * t397;
t448 = qJDD(1) * t333;
t425 = t326 * t448;
t542 = qJD(2) * t278;
t342 = (qJD(1) * t542 - t425) * t327;
t341 = t286 + t342;
t89 = t326 * t135 + t510 * t144;
t75 = qJ(4) * t305 + qJD(4) * t306 + t89;
t388 = t318 * qJDD(1);
t446 = pkin(2) * t402 + qJDD(3);
t230 = -t327 * t388 + t446;
t90 = -pkin(3) * t341 - qJ(4) * t186 + qJD(4) * t254 + t230;
t46 = -t325 * t75 + t328 * t90;
t27 = -pkin(4) * t341 - pkin(9) * t155 + t46;
t47 = t325 * t90 + t328 * t75;
t34 = -pkin(9) * t154 + t47;
t432 = qJD(5) * t534;
t378 = -t332 * t27 - t534 * t34 - t62 * t432 + t458 * t66;
t445 = qJDD(5) - t286;
t543 = -t342 + t445;
t5 = pkin(10) * t543 - t378;
t2 = -qJD(6) * t14 + t335 * t16 - t331 * t5;
t573 = t14 * t569 + t2;
t387 = t26 * t331 - t335 * t48;
t1 = -t387 * qJD(6) + t331 * t16 + t335 * t5;
t559 = t387 * t569 + t1;
t572 = t125 ^ 2;
t351 = t335 * t363;
t102 = t331 * t561 - t351;
t571 = t102 * t125;
t570 = t125 * t363;
t339 = qJDD(5) - t341;
t435 = t534 * t328;
t375 = -t332 * t325 + t435;
t468 = -t375 * t250 + t325 * t458 - t328 * t432;
t568 = t561 ^ 2;
t567 = t102 * t561;
t566 = t104 * t561;
t565 = t250 * t561;
t277 = t325 * t534 + t332 * t328;
t467 = t363 * t277;
t452 = qJD(5) - t289;
t477 = t333 * t326;
t563 = t277 * t445 - t327 * (qJD(1) * (t277 * t542 + t468 * t477) - t277 * t425) - t468 * t452;
t313 = pkin(2) * t326 + qJ(4);
t523 = pkin(9) + t313;
t271 = t523 * t325;
t272 = t523 * t328;
t376 = -t271 * t534 - t332 * t272;
t496 = t250 * t328;
t158 = t231 * t510 - t222;
t173 = pkin(2) * t434 - pkin(3) * t254 + qJ(4) * t250;
t96 = -t158 * t325 + t328 * t173;
t73 = -pkin(4) * t254 + pkin(9) * t496 + t96;
t497 = t250 * t325;
t97 = t328 * t158 + t325 * t173;
t83 = pkin(9) * t497 + t97;
t513 = qJD(4) * t375 + qJD(5) * t376 - t332 * t73 - t534 * t83;
t322 = pkin(12) + qJ(5);
t319 = sin(t322);
t334 = sin(qJ(1));
t337 = cos(qJ(1));
t373 = t420 - t477;
t356 = t329 * t373;
t197 = t334 * t278 + t337 * t356;
t200 = t278 * t337 - t334 * t356;
t256 = t327 * t477 - t397;
t541 = -g(1) * t200 - g(2) * t197 + g(3) * t256;
t562 = t319 * t541;
t416 = -t254 * t331 + t335 * t468;
t558 = pkin(10) * t254 + t513;
t157 = t231 * t326 + t419;
t112 = -pkin(4) * t497 + t157;
t557 = pkin(5) * t467 + pkin(10) * t468 - t112;
t456 = qJD(6) * t331;
t361 = t277 * t456 + t416;
t466 = t278 * t329;
t196 = -t334 * t373 + t337 * t466;
t320 = cos(t322);
t479 = t327 * t337;
t168 = -t196 * t320 - t319 * t479;
t553 = t168 * t331 + t197 * t335;
t552 = t168 * t335 - t197 * t331;
t321 = t327 ^ 2;
t548 = 0.2e1 * t321;
t206 = -t332 * t271 + t272 * t534;
t511 = qJD(4) * t277 + qJD(5) * t206 - t332 * t83 + t534 * t73;
t344 = -t541 + t79;
t546 = t254 * t554;
t457 = qJD(6) * t104;
t36 = -t331 * t60 - t335 * t339 + t457;
t465 = t308 + t309;
t262 = t465 * qJD(2);
t423 = -t534 * t27 + t332 * t34;
t8 = -qJD(5) * t29 - t423;
t540 = -t102 * t467 + t36 * t375;
t35 = -qJD(6) * t351 - t331 * t339 + t335 * t60 + t456 * t561;
t417 = t335 * t254 + t331 * t468;
t455 = qJD(6) * t335;
t362 = t277 * t455 - t417;
t489 = t277 * t331;
t539 = -t104 * t362 + t35 * t489;
t488 = t277 * t335;
t59 = qJDD(6) + t61;
t538 = -t361 * t569 + t59 * t488;
t537 = -t375 * t60 - t467 * t561;
t245 = t250 ^ 2;
t504 = t341 * t325;
t536 = -t245 * t328 + t504;
t201 = t334 * t466 + t337 * t373;
t257 = t278 * t327;
t226 = -t257 * t328 + t325 * t329;
t229 = pkin(2) * t329 + t310 - t403;
t246 = qJ(3) * t480 + t465;
t164 = t326 * t229 + t510 * t246;
t152 = qJ(4) * t329 + t164;
t307 = pkin(2) * t480;
t381 = -pkin(3) * t256 - qJ(4) * t257 + t307;
t533 = pkin(1) * t327;
t176 = -t381 - t533;
t98 = -t152 * t325 + t328 * t176;
t72 = pkin(4) * t256 - pkin(9) * t226 + t98;
t225 = -t257 * t325 - t329 * t328;
t99 = t328 * t152 + t325 * t176;
t81 = -pkin(9) * t225 + t99;
t520 = t332 * t72 + t534 * t81;
t253 = t327 * t542;
t460 = qJD(2) * t333;
t433 = t327 * t460;
t252 = t326 * t433 - t327 * t395;
t494 = t252 * t328;
t303 = qJD(2) * t310;
t216 = t327 * t358 + t303;
t431 = t524 * t327;
t217 = -t327 * t459 + (-t336 * t431 - t309) * qJD(2);
t134 = t510 * t216 + t326 * t217;
t120 = qJD(4) * t329 + t134;
t407 = pkin(2) * t433;
t132 = -pkin(3) * t253 + qJ(4) * t252 + qJD(4) * t257 + t407;
t76 = -t120 * t325 + t328 * t132;
t58 = -pkin(4) * t253 + pkin(9) * t494 + t76;
t495 = t252 * t325;
t77 = t328 * t120 + t325 * t132;
t65 = pkin(9) * t495 + t77;
t12 = -qJD(5) * t520 - t332 * t65 + t534 * t58;
t535 = t254 ^ 2;
t473 = t334 * t336;
t475 = t333 * t337;
t269 = -t329 * t473 - t475;
t530 = g(1) * t269;
t529 = g(1) * t334;
t526 = g(3) * t336;
t525 = t328 * pkin(4);
t317 = -pkin(2) * t510 - pkin(3);
t287 = t317 - t525;
t194 = -pkin(5) * t375 - t277 * pkin(10) + t287;
t115 = t194 * t335 - t206 * t331;
t522 = qJD(6) * t115 + t331 * t557 + t335 * t558;
t116 = t194 * t331 + t206 * t335;
t521 = -qJD(6) * t116 - t331 * t558 + t335 * t557;
t515 = t331 * t59;
t514 = -t102 * t455 - t331 * t36;
t512 = -t254 * pkin(5) + t511;
t509 = t104 * t102;
t507 = t561 * t125;
t506 = t561 * t254;
t505 = t155 * t328;
t501 = t210 * t254;
t500 = t210 * t325;
t499 = t250 * t254;
t498 = t250 * t306;
t493 = t254 * t125;
t492 = t254 * t306;
t491 = t256 * t331;
t486 = t320 * t331;
t485 = t320 * t335;
t484 = t321 * qJD(1) ^ 2;
t482 = t327 * t333;
t481 = t327 * t334;
t179 = t328 * t341;
t476 = t333 * t334;
t471 = t336 * t337;
t469 = -t245 * t325 - t179;
t323 = t333 ^ 2;
t324 = t336 ^ 2;
t464 = t323 - t324;
t454 = qJD(2) - t306;
t453 = -qJD(4) + t128;
t441 = t336 * t484;
t440 = t325 * t481;
t439 = t325 * t479;
t438 = t329 * t471;
t316 = pkin(3) + t525;
t330 = -pkin(9) - qJ(4);
t437 = -t256 * t316 + t257 * t330 + t307;
t429 = pkin(1) * t548;
t426 = t336 * t450;
t424 = g(2) * t479 - g(3) * t329;
t414 = t196 * t319 - t320 * t479;
t133 = t216 * t326 - t510 * t217;
t263 = pkin(2) * t478 - t431;
t412 = -t263 * t334 + t337 * t318;
t411 = t569 ^ 2;
t410 = t569 * t335;
t408 = t306 + t462;
t406 = t305 + t449;
t404 = t333 * t441;
t401 = t333 * t426;
t294 = pkin(2) * t438;
t399 = -pkin(2) * t476 + t294;
t398 = t104 * t467 + t35 * t375;
t396 = -t125 * t468 - t277 * t61;
t394 = pkin(5) * t320 + pkin(10) * t319;
t171 = t201 * t319 - t320 * t481;
t393 = g(1) * t414 + g(2) * t171;
t392 = g(1) * t197 - g(2) * t200;
t391 = g(1) * t337 + g(2) * t334;
t107 = -pkin(4) * t495 + t133;
t386 = -t325 * t84 + t328 * t85;
t38 = pkin(10) * t256 + t520;
t163 = t229 * t510 - t326 * t246;
t153 = -t329 * pkin(3) - t163;
t111 = t225 * pkin(4) + t153;
t140 = t225 * t534 + t226 * t332;
t141 = -t332 * t225 + t226 * t534;
t51 = t140 * pkin(5) - t141 * pkin(10) + t111;
t18 = t331 * t51 + t335 * t38;
t17 = -t331 * t38 + t335 * t51;
t385 = t325 * t554;
t384 = t328 * t554;
t110 = t141 * t335 + t491;
t383 = -t250 * t253 - t256 * t341;
t382 = -t263 * t337 - t318 * t334;
t28 = -t332 * t66 + t534 * t62;
t42 = -t332 * t81 + t534 * t72;
t11 = t332 * t58 + t72 * t432 - t458 * t81 + t534 * t65;
t25 = -pkin(5) * t363 - t28;
t377 = -pkin(10) * t59 + t25 * t569;
t374 = t196 * t330 + t197 * t316 + t399;
t370 = -g(1) * t481 + t424;
t369 = pkin(4) * t440 + t200 * t330 + t201 * t316 + t412;
t368 = -t426 - t448;
t367 = g(1) * t171 - g(2) * t414 - g(3) * (t257 * t319 + t320 * t329);
t172 = t201 * t320 + t319 * t481;
t219 = -t257 * t320 + t319 * t329;
t366 = -g(1) * t172 - g(2) * t168 - g(3) * t219;
t365 = -g(1) * t201 + g(2) * t196 + g(3) * t257;
t360 = t269 * pkin(2);
t6 = -pkin(5) * t543 - t8;
t359 = t367 - t6;
t357 = t368 * pkin(8);
t355 = pkin(4) * t439 + t196 * t316 - t197 * t330 + t382;
t352 = t384 + t500;
t350 = t102 * t361 - t36 * t488;
t349 = t157 * t306 + t541;
t348 = t200 * t316 - t201 * t330 + t360;
t347 = t375 * t339 - t363 * t467;
t346 = pkin(10) * qJD(6) * t569 - t359;
t345 = t1 * t335 - t2 * t331 + (-t14 * t331 + t335 * t387) * qJD(6);
t343 = -t362 * t569 - t489 * t59;
t285 = t305 * t329;
t279 = -t307 - t533;
t273 = -pkin(8) * t482 + t310;
t270 = -t329 * t476 + t471;
t268 = -t329 * t475 - t473;
t267 = -t438 + t476;
t261 = -pkin(8) * t433 + t303;
t260 = t465 * qJD(1);
t259 = -pkin(8) * t434 + t302;
t249 = t256 * t335;
t204 = t327 * t357 + t371;
t203 = -pkin(8) * t402 + t436;
t146 = t325 * t154;
t109 = t141 * t331 - t249;
t106 = t172 * t335 - t200 * t331;
t105 = -t172 * t331 - t200 * t335;
t95 = qJD(5) * t141 - t252 * t277;
t94 = t225 * t432 + t252 * t435 + (qJD(5) * t226 - t495) * t332;
t78 = pkin(5) * t561 - pkin(10) * t125;
t55 = t335 * t59;
t50 = qJD(6) * t110 + t253 * t335 - t331 * t94;
t49 = t141 * t456 + t253 * t331 - t256 * t455 + t335 * t94;
t39 = pkin(5) * t95 + pkin(10) * t94 + t107;
t37 = -t256 * pkin(5) - t42;
t20 = t28 * t335 + t331 * t78;
t19 = -t28 * t331 + t335 * t78;
t10 = t253 * pkin(5) - t12;
t9 = -pkin(10) * t253 + t11;
t4 = -qJD(6) * t18 - t331 * t9 + t335 * t39;
t3 = qJD(6) * t17 + t331 * t39 + t335 * t9;
t7 = [0, 0, 0, 0, 0, qJDD(1), -g(2) * t337 + t529, t391, 0, 0 (qJDD(1) * t323 + 0.2e1 * t401) * t321 (t333 * t447 - t450 * t464) * t548 (qJD(2) * t336 * t408 + t333 * t406) * t327 (qJDD(1) * t324 - 0.2e1 * t401) * t321 (t336 * t406 - t408 * t460) * t327, t285, -g(1) * t268 - g(2) * t270 + t204 * t329 - t262 * t306 + t273 * t305 + (-t427 + t447) * t429, -g(1) * t267 - g(2) * t269 - t203 * t329 - t261 * t306 - t305 * t465 + t368 * t429 ((-qJD(2) * t259 + qJDD(1) * t465 + t203 + (-qJD(2) * t273 + t261) * qJD(1)) * t336 + (-qJD(2) * t260 - qJDD(1) * t273 - t204) * t333 - t391) * t327, t203 * t465 + t260 * t261 + t204 * t273 - t259 * t262 + t321 * qJDD(1) * pkin(1) ^ 2 - g(1) * (-pkin(1) * t334 + pkin(8) * t479) - g(2) * (pkin(1) * t337 + pkin(8) * t481) -t186 * t257 + t252 * t254, -t186 * t256 + t250 * t252 - t253 * t254 - t257 * t341, t186 * t329 - t252 * t306 - t257 * t305, t383, t253 * t306 - t256 * t305 + t329 * t341, t285, -g(1) * t196 - g(2) * t201 - t133 * t306 + t163 * t305 + t230 * t256 + t250 * t407 - t253 * t264 - t279 * t341 + t329 * t88, -t134 * t306 - t164 * t305 + t186 * t279 - t230 * t257 - t252 * t264 - t254 * t407 - t329 * t89 + t392, -t133 * t254 - t134 * t250 + t136 * t252 + t137 * t253 - t163 * t186 + t164 * t341 - t256 * t89 + t257 * t88 - t327 * t391, -g(1) * t382 - g(2) * t412 - t136 * t133 + t137 * t134 + t88 * t163 + t89 * t164 + t230 * t279 + t264 * t407, t155 * t226 + t210 * t494, -t226 * t154 - t155 * t225 - t252 * t352, t155 * t256 + t210 * t253 - t226 * t341 - t250 * t494, t154 * t225 + t252 * t385, -t154 * t256 + t225 * t341 + t250 * t495 - t253 * t554, t383, t76 * t250 - t98 * t341 + t46 * t256 - t84 * t253 - t133 * t554 + t153 * t154 + t79 * t225 - t128 * t495 - g(1) * (t196 * t328 + t439) - g(2) * (t201 * t328 + t440) -t77 * t250 + t99 * t341 - t47 * t256 + t85 * t253 - t133 * t210 + t153 * t155 + t79 * t226 - t128 * t494 - g(1) * (-t196 * t325 + t328 * t479) - g(2) * (-t201 * t325 + t328 * t481) t77 * t554 - t99 * t154 - t47 * t225 + t76 * t210 - t98 * t155 - t46 * t226 + (t325 * t85 + t328 * t84) * t252 - t392, t47 * t99 + t85 * t77 + t46 * t98 + t84 * t76 + t79 * t153 + t128 * t133 - g(1) * (pkin(3) * t196 + qJ(4) * t197 + t382) - g(2) * (pkin(3) * t201 - qJ(4) * t200 + t412) -t141 * t60 - t561 * t94, -t125 * t94 + t140 * t60 - t141 * t61 - t561 * t95, -t94 * t452 + t141 * t445 - t60 * t256 - t561 * t253 + (t141 * t425 + (-t141 * t542 - t477 * t94) * qJD(1)) * t327, -t125 * t95 + t140 * t61, -t95 * t452 - t140 * t445 - t61 * t256 - t125 * t253 + (-t140 * t425 + (t140 * t542 - t477 * t95) * qJD(1)) * t327, t445 * t256 - t452 * t253 + (t256 * t425 + (-t253 * t477 - t256 * t542) * qJD(1)) * t327, g(1) * t168 - g(2) * t172 + t100 * t95 - t107 * t125 + t111 * t61 + t12 * t363 + t53 * t140 - t28 * t253 + t8 * t256 + t339 * t42, -t100 * t94 + t107 * t561 - t11 * t363 - t111 * t60 + t53 * t141 + t29 * t253 + t256 * t378 - t339 * t520 + t393, t11 * t125 - t12 * t561 + t140 * t378 - t141 * t8 + t28 * t94 - t29 * t95 + t42 * t60 - t520 * t61 - t392, -g(1) * t355 - g(2) * t369 + t100 * t107 + t29 * t11 + t53 * t111 + t28 * t12 - t378 * t520 + t8 * t42, -t104 * t49 - t110 * t35, t102 * t49 - t104 * t50 + t109 * t35 - t110 * t36, t104 * t95 + t110 * t59 - t140 * t35 - t49 * t569, t102 * t50 + t109 * t36, -t102 * t95 - t109 * t59 - t140 * t36 - t50 * t569, t140 * t59 + t569 * t95, g(1) * t552 - g(2) * t106 + t10 * t102 + t6 * t109 + t2 * t140 + t17 * t59 + t25 * t50 + t37 * t36 - t387 * t95 + t4 * t569, -g(1) * t553 - g(2) * t105 - t1 * t140 + t10 * t104 + t6 * t110 - t14 * t95 - t18 * t59 - t25 * t49 - t3 * t569 - t37 * t35, -t1 * t109 - t102 * t3 - t104 * t4 - t110 * t2 - t14 * t50 + t17 * t35 - t18 * t36 - t387 * t49 - t393, t1 * t18 + t14 * t3 + t2 * t17 - t387 * t4 + t6 * t37 + t25 * t10 - g(1) * (-pkin(5) * t168 + pkin(10) * t414 + t355) - g(2) * (pkin(5) * t172 + pkin(10) * t171 + t369); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t404, t464 * t484 (qJD(1) * t336 * t454 + t448) * t327, t404 (-t454 * t461 + t447) * t327, t305, -t530 + g(2) * t267 + t260 * t306 + t301 + (-t428 + t484) * t333 * pkin(1) + (t357 - t526) * t327, pkin(1) * t441 + g(1) * t270 - g(2) * t268 + t259 * t306 + (pkin(8) * t450 + g(3)) * t482 - t436, 0, 0, -t499, -t245 + t535, t186 + t498, t499, t341 - t492, t305, t264 * t254 + (-t250 * t434 + t305 * t510) * pkin(2) + t349 + t88, t158 * t306 + t250 * t264 + (t254 * t434 - t305 * t326) * pkin(2) - t365 - t89 -(t137 - t157) * t254 + (-t136 + t158) * t250 + (-t186 * t510 + t326 * t341) * pkin(2), -g(2) * t294 + t136 * t157 - t137 * t158 + (t89 * t326 + t88 * t510 - t530 + g(2) * t476 + (-t264 * t461 - t526) * t327) * pkin(2), t155 * t325 - t210 * t496, t250 * t352 - t146 + t505, -t501 - t536, -t154 * t328 - t250 * t385, t469 + t546, t499, t313 * t504 + t317 * t154 + t157 * t483 + t84 * t254 + (t325 * t453 - t96) * t250 + (t349 - t79) * t328, t313 * t179 + t155 * t317 + t157 * t210 - t254 * t85 + (t328 * t453 + t97) * t250 + t344 * t325, -t96 * t210 - t97 * t483 + (qJD(4) * t554 - t313 * t154 - t84 * t250 - t97 * t306 + t47) * t328 + (-qJD(4) * t210 + t155 * t313 - t250 * t85 - t46) * t325 + t365, t79 * t317 - t85 * t97 - t84 * t96 - t128 * t157 - g(1) * (pkin(3) * t200 + qJ(4) * t201 + t360) - g(2) * (pkin(3) * t197 - qJ(4) * t196 + t399) - g(3) * t381 + (-t46 * t325 + t47 * t328) * t313 + t386 * qJD(4), -t60 * t277 - t468 * t561, t396 + t537, t506 + t563, -t125 * t467 - t375 * t61, t347 + t493, t363 * t254, t100 * t467 + t112 * t125 + t28 * t254 + t287 * t61 + t320 * t541 + t339 * t376 - t363 * t511 - t375 * t53, -t100 * t468 - t112 * t561 - t206 * t339 - t29 * t254 + t53 * t277 - t287 * t60 - t363 * t513 - t562, t125 * t513 - t206 * t61 - t277 * t8 + t28 * t468 - t29 * t467 - t375 * t378 + t376 * t60 + t511 * t561 + t365, -g(1) * t348 - g(2) * t374 - g(3) * t437 - t100 * t112 - t206 * t378 - t28 * t511 + t53 * t287 + t29 * t513 + t376 * t8, -t104 * t361 - t35 * t488, t350 + t539, t398 + t538, t102 * t362 + t36 * t489, t343 + t540, -t375 * t59 + t467 * t569, t115 * t59 - t2 * t375 - t376 * t36 + t6 * t489 - g(1) * (t200 * t485 + t201 * t331) - g(2) * (-t196 * t331 + t197 * t485) - g(3) * (-t256 * t485 - t257 * t331) - t467 * t387 + t521 * t569 + t512 * t102 + t362 * t25, -t116 * t59 + t1 * t375 + t376 * t35 + t6 * t488 - g(1) * (-t200 * t486 + t201 * t335) - g(2) * (-t196 * t335 - t197 * t486) - g(3) * (t256 * t486 - t257 * t335) - t467 * t14 - t522 * t569 + t512 * t104 - t361 * t25, t115 * t35 - t116 * t36 + t417 * t14 - t416 * t387 - t521 * t104 - t522 * t102 + t562 + (-t1 * t331 - t2 * t335 + (-t14 * t335 - t331 * t387) * qJD(6)) * t277, t1 * t116 + t2 * t115 - t6 * t376 - g(1) * (t200 * t394 + t348) - g(2) * (t197 * t394 + t374) - g(3) * (-t256 * t394 + t437) + t512 * t25 + t522 * t14 - t521 * t387; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t341 - t492, t186 - t498, -t245 - t535, -t136 * t254 + t137 * t250 + (-t388 - t529) * t327 + t424 + t446, 0, 0, 0, 0, 0, 0, t469 - t546, -t501 + t536, -t505 - t146 + (t384 - t500) * t250, t128 * t254 + t250 * t386 + t325 * t47 + t328 * t46 + t370, 0, 0, 0, 0, 0, 0, t347 - t493, t506 - t563, t396 - t537, t100 * t254 - t277 * t378 - t28 * t467 - t29 * t468 + t375 * t8 + t370, 0, 0, 0, 0, 0, 0, t343 - t540, t398 - t538, t350 - t539, -t14 * t416 + t25 * t467 + t277 * t345 - t375 * t6 - t387 * t417 + t370; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t210 * t250 + t154, t250 * t554 + t155, -t210 ^ 2 - t554 ^ 2, -t210 * t84 - t554 * t85 + t344, 0, 0, 0, 0, 0, 0, t415 + 0.2e1 * t564 + t565, -t60 + t570, -t568 - t572, -t125 * t29 + t28 * t561 + t344 + t532, 0, 0, 0, 0, 0, 0, -t331 * t411 + t55 - t567, -t335 * t411 - t515 - t566 (t35 + t571) * t335 + t574 + t514, -t561 * t25 + t331 * t559 + t335 * t573 - t541; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t507, t568 - t572, -t60 - t570, t507, -t415 + t565, t339, -t100 * t561 + t250 * t29 + t367 - t423, -t100 * t125 + t28 * t363 - t366 + t378, 0, 0, t104 * t410 - t35 * t331 (-t35 + t571) * t335 - t574 + t514, t410 * t569 + t515 - t566, t102 * t508 - t36 * t335, -t508 * t569 + t55 + t567, -t569 * t561, -pkin(5) * t36 - t102 * t29 - t19 * t569 + t331 * t377 - t335 * t346 + t387 * t561, pkin(5) * t35 - t104 * t29 + t14 * t561 + t20 * t569 + t331 * t346 + t335 * t377, t102 * t20 + t104 * t19 + ((-t36 + t457) * pkin(10) + t559) * t335 + ((qJD(6) * t102 - t35) * pkin(10) - t573) * t331 + t366, t387 * t19 - t14 * t20 - t25 * t29 + t359 * pkin(5) + (t345 + t366) * pkin(10); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t509, -t102 ^ 2 + t104 ^ 2, t102 * t569 - t35, -t509, t104 * t569 - t36, t59, -t25 * t104 - g(1) * t105 + g(2) * t553 - g(3) * (-t219 * t331 + t249) + t573, t25 * t102 + g(1) * t106 + g(2) * t552 - g(3) * (-t219 * t335 - t491) - t559, 0, 0;];
tau_reg  = t7;