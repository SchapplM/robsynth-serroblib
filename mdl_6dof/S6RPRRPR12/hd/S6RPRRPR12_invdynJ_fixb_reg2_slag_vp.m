% Calculate inertial parameters regressor of inverse dynamics joint torque vector for
% S6RPRRPR12
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
%   pkin=[a2,a3,a4,a5,a6,alpha2,alpha3,d1,d3,d4,d6,theta2]';
% 
% Output:
% tau_reg [6x(6*10)]
%   inertial parameter regressor of inverse dynamics joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 05:54
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tau_reg = S6RPRRPR12_invdynJ_fixb_reg2_slag_vp(qJ, qJD, qJDD, g, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(3,1),zeros(12,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRRPR12_invdynJ_fixb_reg2_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RPRRPR12_invdynJ_fixb_reg2_slag_vp: qJD has to be [6x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [6 1]), ...
  'S6RPRRPR12_invdynJ_fixb_reg2_slag_vp: qJDD has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RPRRPR12_invdynJ_fixb_reg2_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [12 1]), ...
  'S6RPRRPR12_invdynJ_fixb_reg2_slag_vp: pkin has to be [12x1] (double)');

%% Symbolic Calculation
% From invdyn_joint_fixb_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 05:53:06
% EndTime: 2019-03-09 05:53:38
% DurationCPUTime: 18.02s
% Computational Cost: add. (28458->806), mult. (90669->1032), div. (0->0), fcn. (78451->14), ass. (0->370)
t297 = cos(pkin(6));
t296 = cos(pkin(12));
t525 = sin(qJ(1));
t447 = t525 * t296;
t294 = sin(pkin(12));
t527 = cos(qJ(1));
t450 = t527 * t294;
t242 = t297 * t450 + t447;
t300 = sin(qJ(3));
t526 = cos(qJ(3));
t287 = t525 * t294;
t449 = t527 * t296;
t379 = -t297 * t449 + t287;
t295 = sin(pkin(6));
t503 = sin(pkin(7));
t440 = t295 * t503;
t504 = cos(pkin(7));
t555 = t379 * t504 + t527 * t440;
t166 = -t242 * t526 + t300 * t555;
t299 = sin(qJ(4));
t302 = cos(qJ(4));
t441 = t295 * t504;
t547 = t379 * t503 - t527 * t441;
t120 = t166 * t299 + t302 * t547;
t163 = t242 * t300 + t526 * t555;
t298 = sin(qJ(6));
t301 = cos(qJ(6));
t559 = t120 * t298 - t163 * t301;
t558 = t120 * t301 + t163 * t298;
t438 = t300 * t504;
t344 = t295 * (-t294 * t438 + t296 * t526);
t226 = qJD(1) * t344;
t410 = t503 * t526;
t380 = qJD(3) * t410;
t539 = t226 - t380;
t485 = t294 * t295;
t446 = qJD(1) * t485;
t411 = t504 * t526;
t484 = t295 * t296;
t372 = t411 * t484;
t329 = t297 * t410 + t372;
t538 = t329 * qJD(1);
t197 = t300 * t446 - t538;
t188 = qJD(4) + t197;
t121 = t166 * t302 - t299 * t547;
t456 = t300 * t485;
t417 = qJD(3) * t456;
t254 = qJD(1) * t417;
t432 = t503 * t297;
t210 = t295 * (t294 * t526 + t296 * t438) + t300 * t432;
t317 = t210 * qJDD(1);
t313 = t317 - t254;
t324 = t329 * qJD(3);
t307 = qJD(1) * t324 + t313;
t230 = (-pkin(9) * t294 * t503 - pkin(2) * t296 - pkin(1)) * t295;
t218 = qJD(1) * t230 + qJD(2);
t524 = pkin(1) * t297;
t458 = qJD(1) * t524;
t277 = t296 * t458;
t347 = t297 * pkin(2) + (-pkin(9) * t504 - qJ(2)) * t485;
t190 = qJD(1) * t347 + t277;
t435 = t504 * t190;
t279 = qJ(2) * t484;
t237 = qJD(1) * t279 + t294 * t458;
t342 = (t296 * t441 + t432) * pkin(9);
t183 = qJD(1) * t342 + t237;
t448 = t526 * t183;
t103 = t448 + (t218 * t503 + t435) * t300;
t472 = qJD(4) * t299;
t488 = t197 * t299;
t552 = -qJD(5) * t299 - t103 + (t472 + t488) * pkin(4);
t200 = t210 * qJD(3);
t499 = qJ(5) * t302;
t551 = t552 + t188 * (pkin(11) * t299 - t499);
t201 = t210 * qJD(1);
t135 = pkin(3) * t201 + pkin(10) * t197;
t528 = pkin(5) + pkin(10);
t273 = t528 * t302;
t529 = pkin(4) + pkin(11);
t387 = t190 * t411;
t102 = -t300 * t183 + t218 * t410 + t387;
t98 = t299 * t102;
t550 = qJD(4) * t273 - t98 - (-pkin(5) * t197 - t135) * t302 + t529 * t201;
t142 = -t190 * t503 + t218 * t504;
t331 = t197 * pkin(3) - t201 * pkin(10) + t142;
t465 = qJD(1) * qJD(2);
t444 = t295 * t465;
t462 = qJDD(1) * t297;
t457 = pkin(1) * t462;
t222 = qJDD(1) * t279 + t294 * t457 + t296 * t444;
t173 = qJDD(1) * t342 + t222;
t275 = t296 * t457;
t416 = t294 * t444;
t174 = qJDD(1) * t347 + t275 - t416;
t214 = qJDD(1) * t230 + qJDD(2);
t437 = t300 * t503;
t474 = qJD(3) * t300;
t351 = -qJD(3) * t387 - t526 * t173 - t174 * t438 + t183 * t474 - t214 * t437 - t218 * t380;
t439 = t297 * t504;
t415 = t296 * t440;
t461 = -qJDD(1) * t415 + qJDD(3);
t364 = -qJDD(1) * t439 - t461;
t549 = -pkin(10) * t364 + qJD(4) * t331 - t351;
t428 = qJD(1) * t503;
t409 = t295 * t428;
t466 = -t296 * t409 + qJD(3);
t365 = -qJD(1) * t439 - t466;
t85 = -pkin(10) * t365 + t103;
t42 = t299 * t85 - t302 * t331;
t481 = qJD(5) + t42;
t247 = t299 * t504 + t302 * t437;
t383 = t294 * t409;
t478 = -qJD(4) * t247 + t299 * t539 - t302 * t383;
t209 = -t329 + t456;
t359 = t297 * t447 + t450;
t548 = t359 * t503 + t525 * t441;
t343 = t295 * (t294 * t411 + t296 * t300);
t225 = qJD(1) * t343;
t408 = qJD(3) * t437;
t376 = t408 - t225;
t245 = t294 * t524 + t279;
t203 = t342 + t245;
t284 = t296 * t524;
t211 = t284 + t347;
t108 = t526 * t203 + (t504 * t211 + t503 * t230) * t300;
t148 = t302 * t201 - t299 * t365;
t378 = qJD(3) * t448 + t300 * t173 - t174 * t411 - t214 * t410 + t218 * t408 + t435 * t474;
t50 = pkin(3) * t364 + t378;
t350 = t302 * t365;
t75 = qJD(4) * t350 + t201 * t472 + t299 * t364 - t302 * t307;
t312 = t75 * qJ(5) - t148 * qJD(5) + t50;
t473 = qJD(4) * t148;
t314 = qJD(3) * t538 - t254;
t537 = t314 + t317;
t76 = t299 * t537 + t302 * t364 + t473;
t10 = t529 * t76 + t312;
t482 = pkin(5) * t148 + t481;
t29 = -t188 * t529 + t482;
t146 = t201 * t299 + t350;
t84 = pkin(3) * t365 - t102;
t309 = -t148 * qJ(5) + t84;
t38 = t146 * t529 + t309;
t13 = t29 * t301 - t298 * t38;
t366 = t209 * qJDD(1);
t134 = qJD(1) * t200 + t366;
t133 = qJDD(4) + t134;
t470 = qJD(4) * t302;
t126 = -t174 * t503 + t504 * t214;
t63 = t134 * pkin(3) - pkin(10) * t537 + t126;
t430 = t299 * t549 - t302 * t63 + t85 * t470;
t403 = qJDD(5) + t430;
t6 = -pkin(5) * t75 - t133 * t529 + t403;
t1 = qJD(6) * t13 + t301 * t10 + t298 * t6;
t145 = qJD(6) + t148;
t546 = -t13 * t145 + t1;
t14 = t29 * t298 + t301 * t38;
t2 = -qJD(6) * t14 - t298 * t10 + t301 * t6;
t545 = t14 * t145 + t2;
t468 = qJD(6) * t301;
t469 = qJD(6) * t298;
t35 = -t301 * t133 - t146 * t468 + t188 * t469 - t298 * t76;
t104 = -t301 * t146 + t188 * t298;
t420 = t104 * t145;
t544 = t35 - t420;
t106 = t146 * t298 + t188 * t301;
t36 = qJD(6) * t106 + t133 * t298 - t301 * t76;
t496 = t106 * t145;
t543 = -t36 + t496;
t492 = t148 * t188;
t542 = t76 - t492;
t291 = t294 ^ 2;
t292 = t295 ^ 2;
t293 = t296 ^ 2;
t541 = t292 * (t291 + t293);
t107 = -t300 * t203 + t211 * t411 + t230 * t410;
t536 = t359 * t504 - t525 * t440;
t246 = t299 * t437 - t302 * t504;
t43 = t299 * t331 + t302 * t85;
t40 = -t188 * qJ(5) - t43;
t521 = t146 * pkin(5);
t30 = -t40 - t521;
t74 = -qJDD(6) + t75;
t534 = t145 * t30 + t529 * t74;
t533 = -t246 * t133 + t146 * t376 + t188 * t478 - t410 * t76;
t477 = qJD(4) * t246 + t299 * t383 + t302 * t539;
t532 = -t247 * t133 + t148 * t376 + t188 * t477 + t410 * t75;
t531 = t148 ^ 2;
t530 = t188 ^ 2;
t523 = pkin(4) * t133;
t522 = pkin(10) * t133;
t241 = t415 - t439;
t161 = t210 * t299 + t241 * t302;
t520 = t161 * pkin(11);
t519 = t163 * pkin(10);
t243 = -t287 * t297 + t449;
t167 = t243 * t300 + t526 * t536;
t518 = t167 * pkin(10);
t101 = -pkin(10) * t241 + t108;
t150 = -t211 * t503 + t504 * t230;
t208 = t209 * pkin(3);
t94 = -t210 * pkin(10) + t150 + t208;
t56 = t302 * t101 + t299 * t94;
t515 = t188 * t40;
t514 = t188 * t42;
t513 = t188 * t43;
t512 = t298 * t36;
t511 = t298 * t74;
t510 = t301 * t35;
t71 = t301 * t74;
t500 = qJ(5) * t299;
t443 = -pkin(3) - t500;
t248 = -t302 * t529 + t443;
t272 = t528 * t299;
t212 = -t248 * t298 + t272 * t301;
t508 = qJD(6) * t212 + t550 * t298 + t301 * t551;
t213 = t248 * t301 + t272 * t298;
t507 = -qJD(6) * t213 - t298 * t551 + t550 * t301;
t506 = -qJ(5) * t470 - t197 * t499 + t552;
t67 = t302 * t102 + t299 * t135;
t57 = -qJ(5) * t201 - t67;
t505 = -pkin(5) * t488 - t528 * t472 + t57;
t502 = pkin(1) * qJDD(1);
t501 = qJ(5) * t146;
t498 = t104 * t188;
t497 = t106 * t104;
t495 = t106 * t188;
t130 = t133 * qJ(5);
t494 = t146 * t148;
t493 = t146 * t188;
t491 = t163 * t302;
t490 = t167 * t302;
t489 = t188 * t201;
t487 = t201 * t197;
t486 = t209 * t302;
t304 = qJD(1) ^ 2;
t483 = t297 * t304;
t353 = -t298 * t246 + t301 * t410;
t480 = -qJD(6) * t353 + t298 * t376 + t301 * t478;
t219 = t301 * t246 + t298 * t410;
t479 = -qJD(6) * t219 + t298 * t478 - t301 * t376;
t451 = t295 * t525;
t476 = t527 * pkin(1) + qJ(2) * t451;
t471 = qJD(4) * t301;
t467 = qJD(6) * t302;
t464 = qJDD(1) * t292;
t463 = qJDD(1) * t296;
t460 = g(1) * t525;
t459 = g(2) * t527;
t157 = t163 * pkin(3);
t455 = -pkin(4) * t491 - t163 * t500 - t157;
t159 = t167 * pkin(3);
t454 = -pkin(4) * t490 - t167 * t500 - t159;
t453 = -pkin(4) * t486 - t209 * t500 - t208;
t452 = t295 * t527;
t445 = -t146 ^ 2 + t531;
t442 = g(2) * t452 - g(3) * t297;
t55 = -t299 * t101 + t302 * t94;
t431 = -t299 * t63 - t302 * t549 + t85 * t472;
t66 = t135 * t302 - t98;
t427 = pkin(4) * t120 - qJ(5) * t121;
t168 = t243 * t526 - t300 * t536;
t122 = t168 * t299 - t302 * t548;
t123 = t168 * t302 + t299 * t548;
t426 = -t122 * pkin(4) + qJ(5) * t123;
t162 = t210 * t302 - t241 * t299;
t425 = -t161 * pkin(4) + t162 * qJ(5);
t423 = t188 * t299;
t422 = t188 * t302;
t421 = t145 * t298;
t418 = 0.2e1 * t295 * t462;
t47 = -qJ(5) * t209 - t56;
t412 = -pkin(1) * t525 + qJ(2) * t452;
t407 = g(1) * t120 + g(2) * t122;
t406 = -g(1) * t121 - g(2) * t123;
t405 = -g(1) * t163 + g(2) * t167;
t127 = t201 * t298 + t301 * t488;
t402 = t299 * t471 + t127;
t128 = t201 * t301 - t298 * t488;
t401 = t298 * t472 - t128;
t399 = (qJD(4) * t146 - t75) * pkin(10);
t398 = (-t76 + t473) * pkin(10);
t396 = -t13 * t298 + t14 * t301;
t32 = pkin(5) * t162 - t209 * t529 - t55;
t100 = t241 * pkin(3) - t107;
t62 = t100 - t425;
t44 = t62 + t520;
t16 = -t298 * t44 + t301 * t32;
t17 = t298 * t32 + t301 * t44;
t199 = -qJD(3) * t372 - t297 * t380 + t417;
t110 = qJD(4) * t162 - t199 * t299;
t395 = t110 * t146 + t161 * t76;
t111 = -t210 * t472 + (-qJD(4) * t241 - t199) * t302;
t394 = t111 * t148 - t162 * t75;
t393 = -t133 * t302 - t146 * t201;
t392 = -t133 * t299 + t148 * t201;
t391 = t161 * t301 - t209 * t298;
t113 = t161 * t298 + t209 * t301;
t390 = (-qJ(2) * t446 + t277) * t294 - t237 * t296;
t382 = t294 * qJD(2) * t440;
t125 = t200 * pkin(3) + t199 * pkin(10) + t382;
t88 = qJD(2) * t344 + qJD(3) * t107;
t28 = -t101 * t470 + t125 * t302 - t299 * t88 - t94 * t472;
t184 = t188 * qJD(5);
t8 = -t130 - t184 + t431;
t377 = -t145 * t421 - t71;
t278 = -t295 * t502 + qJDD(2);
t375 = pkin(1) * t464 - t278 * t295;
t374 = t188 * t84 - t522;
t52 = t146 * pkin(4) + t309;
t373 = -t188 * t52 + t522;
t27 = -t101 * t472 + t299 * t125 + t302 * t88 + t94 * t470;
t371 = g(1) * t527 + g(2) * t525;
t370 = g(1) * t122 - g(2) * t120 + g(3) * t161;
t369 = -g(1) * t123 + g(2) * t121 - g(3) * t162;
t368 = g(1) * t167 + g(2) * t163 + g(3) * t209;
t367 = g(1) * t168 - g(2) * t166 + g(3) * t210;
t363 = -t145 ^ 2 * t301 + t511;
t362 = t148 * t422 - t299 * t75;
t361 = t146 * t423 - t302 * t76;
t360 = -g(1) * t451 + t442;
t349 = t75 - t493;
t348 = -t110 * t148 - t111 * t146 + t161 * t75 - t162 * t76;
t346 = t110 * t188 + t133 * t161 + t146 * t200 + t209 * t76;
t345 = t111 * t188 + t133 * t162 + t148 * t200 - t209 * t75;
t341 = pkin(10) * qJD(4) * t188 - t368;
t340 = t370 - t430;
t339 = t369 - t431;
t15 = t76 * pkin(4) + t312;
t337 = t15 + t341;
t336 = t341 + t50;
t7 = -pkin(5) * t76 - t8;
t335 = qJD(6) * t145 * t529 + t369 + t7;
t24 = -qJ(5) * t200 - qJD(5) * t209 - t27;
t330 = t146 * t477 - t148 * t478 - t246 * t75 - t247 * t76;
t328 = t148 * t52 + qJDD(5) - t340;
t326 = (-t75 - t493) * t302 + (-t76 - t492) * t299;
t320 = -t242 * pkin(2) - pkin(9) * t547 + t412;
t319 = t243 * pkin(2) + pkin(9) * t548 + t476;
t318 = t166 * pkin(3) + t320;
t316 = t168 * pkin(3) + t319;
t311 = t121 * pkin(4) + t120 * qJ(5) + t318;
t310 = t123 * pkin(4) + t122 * qJ(5) + t316;
t89 = qJD(2) * t343 + qJD(3) * t108;
t308 = -t111 * qJ(5) - t162 * qJD(5) + t89;
t264 = -pkin(4) * t302 + t443;
t244 = -qJ(2) * t485 + t284;
t221 = t275 + (-qJ(2) * qJDD(1) - t465) * t485;
t92 = pkin(4) * t148 + t501;
t78 = t122 * t298 + t167 * t301;
t77 = t122 * t301 - t167 * t298;
t70 = t133 * t209 + t188 * t200;
t69 = t148 * t529 + t501;
t65 = qJD(6) * t391 + t110 * t298 + t200 * t301;
t64 = qJD(6) * t113 - t110 * t301 + t200 * t298;
t59 = -pkin(4) * t201 - t66;
t48 = -pkin(4) * t209 - t55;
t39 = -pkin(4) * t188 + t481;
t37 = -pkin(5) * t161 - t47;
t34 = t43 - t521;
t31 = t110 * pkin(4) + t308;
t26 = t110 * t529 + t308;
t25 = -pkin(4) * t200 - t28;
t21 = t298 * t34 + t301 * t69;
t20 = -t298 * t69 + t301 * t34;
t19 = -pkin(5) * t110 - t24;
t18 = pkin(5) * t111 - t200 * t529 - t28;
t9 = t403 - t523;
t4 = -qJD(6) * t17 + t18 * t301 - t26 * t298;
t3 = qJD(6) * t16 + t18 * t298 + t26 * t301;
t5 = [0, 0, 0, 0, 0, qJDD(1), t460 - t459, t371, 0, 0, t291 * t464, 0.2e1 * t292 * t294 * t463, t294 * t418, t293 * t464, t296 * t418, t297 ^ 2 * qJDD(1), g(1) * t242 - g(2) * t243 + t375 * t296 + (qJDD(1) * t244 + t221 - t416) * t297, -g(1) * t287 + (t459 - t375) * t294 + (-t245 * qJDD(1) - t222 + (t371 - t444) * t296) * t297, t465 * t541 + (-t221 * t294 + t222 * t296 + (-t244 * t294 + t245 * t296) * qJDD(1) - t371) * t295, t222 * t245 + t221 * t244 - g(1) * t412 - g(2) * t476 + (-t278 * pkin(1) - qJD(2) * t390) * t295, -t201 * t199 + t210 * t537, -t210 * t134 + t199 * t197 - t201 * t200 - t209 * t537, t199 * t365 - t314 * t241 + (t461 + (-t241 + t439) * qJDD(1)) * t210, t134 * t209 + t197 * t200, t134 * t241 + t200 * t365 + t209 * t364, t364 * t241, -g(1) * t166 - g(2) * t168 - t107 * t364 + t126 * t209 + t150 * t134 + t142 * t200 + t197 * t382 + t241 * t378 + t365 * t89, t108 * t364 + t126 * t210 - t142 * t199 + t150 * t307 + t201 * t382 - t241 * t351 + t365 * t88 + t405, g(1) * t547 - g(2) * t548 + t102 * t199 - t103 * t200 - t107 * t307 - t108 * t134 - t88 * t197 + t89 * t201 + t209 * t351 + t210 * t378, -g(1) * t320 - g(2) * t319 - t102 * t89 + t103 * t88 - t107 * t378 - t108 * t351 + t126 * t150 + t142 * t382, t394, t348, t345, t395, -t346, t70, t100 * t76 + t110 * t84 + t133 * t55 + t146 * t89 + t161 * t50 + t188 * t28 - t200 * t42 - t209 * t430 + t406, -t100 * t75 + t111 * t84 - t133 * t56 + t148 * t89 + t162 * t50 - t188 * t27 - t200 * t43 + t209 * t431 + t407, -t110 * t43 + t111 * t42 - t146 * t27 - t148 * t28 + t161 * t431 + t162 * t430 + t55 * t75 - t56 * t76 - t405, -t431 * t56 + t43 * t27 - t430 * t55 - t42 * t28 + t50 * t100 + t84 * t89 - g(1) * (t318 - t519) - g(2) * (t316 + t518) t70, -t345, t346, t394, t348, t395, t110 * t40 + t111 * t39 + t146 * t24 + t148 * t25 + t161 * t8 + t162 * t9 + t47 * t76 - t48 * t75 - t405, -t110 * t52 + t133 * t48 - t146 * t31 - t15 * t161 + t188 * t25 + t200 * t39 + t209 * t9 - t62 * t76 - t406, -t111 * t52 - t133 * t47 - t148 * t31 - t15 * t162 - t188 * t24 - t200 * t40 - t209 * t8 + t62 * t75 - t407, t15 * t62 + t52 * t31 + t8 * t47 + t40 * t24 + t9 * t48 + t39 * t25 - g(1) * (t311 - t519) - g(2) * (t310 + t518) t106 * t65 - t113 * t35, -t104 * t65 - t106 * t64 - t113 * t36 - t35 * t391, t106 * t111 - t113 * t74 + t145 * t65 - t162 * t35, t104 * t64 - t36 * t391, -t104 * t111 - t145 * t64 - t162 * t36 - t391 * t74, t111 * t145 - t162 * t74, -g(1) * t559 - g(2) * t78 + t19 * t104 + t13 * t111 + t4 * t145 - t16 * t74 + t2 * t162 + t30 * t64 + t37 * t36 - t7 * t391, -g(1) * t558 - g(2) * t77 - t1 * t162 + t19 * t106 - t14 * t111 + t7 * t113 - t3 * t145 + t17 * t74 + t30 * t65 - t37 * t35, t1 * t391 - t104 * t3 - t106 * t4 - t113 * t2 - t13 * t65 - t14 * t64 + t16 * t35 - t17 * t36 + t406, t1 * t17 + t14 * t3 + t2 * t16 + t13 * t4 + t7 * t37 + t30 * t19 - g(1) * (t121 * pkin(11) - t163 * t528 + t311) - g(2) * (t123 * pkin(11) + t167 * t528 + t310); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 (t294 * t483 - t463) * t295 (qJDD(1) * t294 + t296 * t483) * t295, -t304 * t541, qJDD(2) + (qJD(1) * t390 - t460 - t502) * t295 + t442, 0, 0, 0, 0, 0, 0, t504 * t134 - t197 * t383 - t364 * t410 + t365 * t376, -t201 * t383 + t504 * t307 + t364 * t437 - t365 * t539, -t134 * t437 + t197 * t539 + t201 * t376 - t307 * t410, -t378 * t410 - t351 * t437 + t126 * t504 + t102 * t225 - t103 * t226 + (-t142 * t294 * t428 - t460) * t295 + (-t102 * t437 + t103 * t410) * qJD(3) + t442, 0, 0, 0, 0, 0, 0, t533, t532, t330, t246 * t430 - t247 * t431 + t376 * t84 - t410 * t50 - t42 * t478 - t43 * t477 + t360, 0, 0, 0, 0, 0, 0, t330, -t533, -t532, -t15 * t410 + t9 * t246 - t8 * t247 + t376 * t52 - t39 * t478 + t40 * t477 + t360, 0, 0, 0, 0, 0, 0, -t104 * t477 - t145 * t480 - t219 * t74 + t247 * t36, -t106 * t477 + t145 * t479 - t247 * t35 - t353 * t74, t104 * t479 + t106 * t480 + t219 * t35 + t353 * t36, -t1 * t353 - t13 * t480 - t14 * t479 + t2 * t219 + t7 * t247 - t30 * t477 + t360; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t487, -t197 ^ 2 + t201 ^ 2, t197 * t466 + (t197 * t439 + t324) * qJD(1) + t313, -t487, t201 * t466 + (t201 * t439 - t200) * qJD(1) - t366, -t364, -t103 * t365 - t142 * t201 + t368 - t378, -t102 * t365 + t142 * t197 + t351 + t367, 0, 0, t362, t326, t188 * t422 - t392, t361, -t299 * t530 - t393, -t489, -pkin(3) * t76 - t103 * t146 - t188 * t66 + t201 * t42 + t299 * t374 - t302 * t336, pkin(3) * t75 - t103 * t148 + t188 * t67 + t201 * t43 + t299 * t336 + t302 * t374, t146 * t67 + t148 * t66 + (-t431 + t398 + t514) * t302 + (t430 + t399 - t513) * t299 - t367, -t50 * pkin(3) + g(1) * t159 + g(2) * t157 + g(3) * t208 - t84 * t103 + t42 * t66 - t43 * t67 + (-t431 * t302 + t430 * t299 + (-t299 * t43 + t302 * t42) * qJD(4) - t367) * pkin(10), -t489, -t302 * t530 + t392, t188 * t423 + t393, t362, t326, t361, -t146 * t57 - t148 * t59 + (t188 * t39 + t398 - t8) * t302 + (t399 + t9 + t515) * t299 - t367, -t146 * t506 - t188 * t59 - t201 * t39 - t264 * t76 + t299 * t373 + t302 * t337, -t148 * t506 + t188 * t57 + t201 * t40 + t264 * t75 - t299 * t337 + t302 * t373, t15 * t264 - t40 * t57 - t39 * t59 - g(1) * t454 - g(2) * t455 - g(3) * t453 + t506 * t52 + (t9 * t299 - t8 * t302 + (t299 * t40 + t302 * t39) * qJD(4) - t367) * pkin(10), t298 * t302 * t35 + (-t301 * t467 + t401) * t106, t104 * t128 + t106 * t127 + (-t104 * t298 + t106 * t301) * t472 + (t512 + t510 + (t104 * t301 + t106 * t298) * qJD(6)) * t302, -t299 * t35 + t401 * t145 + (-t145 * t468 + t495 + t511) * t302, t301 * t302 * t36 + (-t298 * t467 - t402) * t104, -t299 * t36 + t402 * t145 + (t145 * t469 - t498 + t71) * t302, t145 * t422 - t299 * t74, -t30 * t127 - t212 * t74 + t273 * t36 - t367 * t301 + (t298 * t368 - t30 * t471 + t2) * t299 + t507 * t145 + t505 * t104 + (t13 * t188 - t30 * t469 + t7 * t301) * t302, -t30 * t128 + t213 * t74 - t273 * t35 + t367 * t298 + (t30 * t298 * qJD(4) + t301 * t368 - t1) * t299 - t508 * t145 + t505 * t106 + (-t14 * t188 - t7 * t298 - t30 * t468) * t302, t127 * t14 + t128 * t13 + t212 * t35 - t213 * t36 - t507 * t106 - t508 * t104 + t396 * t472 + (-t1 * t301 + t2 * t298 + (t13 * t301 + t14 * t298) * qJD(6) + t368) * t302, t1 * t213 + t2 * t212 + t7 * t273 - g(1) * (-pkin(11) * t490 + t168 * t528 + t454) - g(2) * (-pkin(11) * t491 - t166 * t528 + t455) - g(3) * (-pkin(11) * t486 + t210 * t528 + t453) + t505 * t30 + t508 * t14 + t507 * t13; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t494, t445, -t349, -t494, -t542, t133, -t148 * t84 + t340 + t513, t146 * t84 - t339 - t514, 0, 0, t133, t349, t542, t494, t445, -t494, pkin(4) * t75 - qJ(5) * t76 + (-t40 - t43) * t148 + (t39 - t481) * t146, t146 * t92 + t328 - t513 - 0.2e1 * t523, -t146 * t52 + t148 * t92 + t188 * t481 + 0.2e1 * t130 + t184 + t339, -t9 * pkin(4) - g(1) * t426 - g(2) * t427 - g(3) * t425 - t8 * qJ(5) - t39 * t43 - t40 * t481 - t52 * t92, -t106 * t421 - t510 (-t36 - t496) * t301 + (t35 + t420) * t298, t106 * t146 + t377, t301 * t420 + t512, -t104 * t146 + t363, t145 * t146, qJ(5) * t36 + t482 * t104 + t13 * t146 - t145 * t20 + t335 * t298 + t301 * t534, -qJ(5) * t35 + t482 * t106 - t14 * t146 + t145 * t21 - t298 * t534 + t335 * t301, t104 * t21 + t106 * t20 + (-t14 * t148 - t529 * t35 - t2 + (t104 * t529 - t14) * qJD(6)) * t301 + (t13 * t148 + t529 * t36 - t1 + (-t106 * t529 + t13) * qJD(6)) * t298 + t370, t7 * qJ(5) - t14 * t21 - t13 * t20 - g(1) * (-pkin(11) * t122 + t426) - g(2) * (pkin(11) * t120 + t427) - g(3) * (t425 - t520) + t482 * t30 - (qJD(6) * t396 + t1 * t298 + t2 * t301) * t529; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t349, t133 - t494, -t530 - t531, t328 + t515 - t523, 0, 0, 0, 0, 0, 0, t377 - t498, t363 - t495, t298 * t543 + t301 * t544, -t188 * t30 + t298 * t546 + t545 * t301 - t370; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t497, -t104 ^ 2 + t106 ^ 2, -t544, -t497, t543, -t74, -g(1) * t77 + g(2) * t558 - g(3) * t391 - t30 * t106 + t545, g(1) * t78 - g(2) * t559 + g(3) * t113 + t30 * t104 - t546, 0, 0;];
tau_reg  = t5;