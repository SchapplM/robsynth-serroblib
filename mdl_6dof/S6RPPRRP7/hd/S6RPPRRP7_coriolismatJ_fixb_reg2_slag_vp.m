% Calculate inertial parameters regressor of coriolis matrix for
% S6RPPRRP7
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d4,d5,theta3]';
% 
% Output:
% cmat_reg [(6*6)x(6*10)]
%   inertial parameter regressor of coriolis matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 02:14
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function cmat_reg = S6RPPRRP7_coriolismatJ_fixb_reg2_slag_vp(qJ, qJD, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(9,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPPRRP7_coriolismatJ_fixb_reg2_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RPPRRP7_coriolismatJ_fixb_reg2_slag_vp: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S6RPPRRP7_coriolismatJ_fixb_reg2_slag_vp: pkin has to be [9x1] (double)');

%% Symbolic Calculation
% From coriolismat_joint_fixb_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 02:14:01
% EndTime: 2019-03-09 02:14:11
% DurationCPUTime: 7.85s
% Computational Cost: add. (9296->464), mult. (16656->570), div. (0->0), fcn. (17951->6), ass. (0->379)
t404 = sin(pkin(9));
t588 = sin(qJ(4));
t388 = t588 * t404;
t405 = cos(pkin(9));
t589 = cos(qJ(4));
t470 = t589 * t405;
t370 = t470 - t388;
t406 = sin(qJ(5));
t580 = pkin(8) + qJ(6);
t375 = t580 * t406;
t407 = cos(qJ(5));
t524 = t407 * t375;
t376 = t580 * t407;
t527 = t406 * t376;
t420 = -t527 / 0.2e1 + t524 / 0.2e1;
t387 = pkin(3) * t404 + qJ(2);
t469 = t588 * t405;
t471 = t589 * t404;
t368 = t469 + t471;
t586 = t368 * pkin(4);
t440 = -pkin(8) * t370 + t586;
t417 = t387 + t440;
t581 = -pkin(1) - qJ(3);
t478 = -pkin(7) + t581;
t445 = t478 * t405;
t423 = t588 * t445;
t372 = t478 * t404;
t472 = t589 * t372;
t277 = t472 + t423;
t525 = t407 * t277;
t149 = t406 * t417 + t525;
t245 = t406 * t370;
t118 = -qJ(6) * t245 + t149;
t558 = t118 * t407;
t110 = t558 / 0.2e1;
t528 = t406 * t277;
t148 = -t407 * t417 + t528;
t250 = t407 * t370;
t117 = qJ(6) * t250 + t148;
t585 = t368 * pkin(5);
t86 = -t117 + t585;
t574 = t86 * t406;
t570 = t110 - t574 / 0.2e1;
t603 = t370 * t420 + t570;
t402 = t406 ^ 2;
t403 = t407 ^ 2;
t592 = -t403 / 0.2e1;
t448 = t592 - t402 / 0.2e1;
t474 = t406 * t250;
t447 = 0.2e1 * t474;
t602 = t117 + t86;
t487 = t370 * qJD(1);
t340 = t407 * t487;
t508 = qJD(4) * t406;
t307 = t340 + t508;
t365 = t368 ^ 2;
t366 = t370 ^ 2;
t600 = -t366 - t365;
t477 = t366 - t365;
t398 = t404 ^ 2;
t399 = t405 ^ 2;
t379 = t398 + t399;
t442 = t365 / 0.2e1 + t366 / 0.2e1;
t276 = t588 * t372 - t589 * t445;
t253 = t276 * t407;
t583 = t370 * pkin(4);
t584 = t368 * pkin(8);
t278 = t583 + t584;
t259 = t406 * t278;
t159 = -t253 + t259;
t548 = t159 * t407;
t260 = t407 * t278;
t537 = t276 * t406;
t158 = t260 + t537;
t551 = t158 * t406;
t553 = t149 * t407;
t554 = t148 * t406;
t14 = (t553 / 0.2e1 + t554 / 0.2e1 - t277 / 0.2e1) * t370 + (t548 / 0.2e1 - t551 / 0.2e1 + t276 / 0.2e1) * t368;
t384 = t402 + t403;
t131 = (-0.1e1 + t384) * t370 * t368;
t500 = t131 * qJD(2);
t599 = -t14 * qJD(1) - t500;
t243 = t406 * t368;
t194 = -pkin(5) * t243 + t277;
t139 = qJ(6) * t243 + t159;
t556 = t139 * t407;
t248 = t407 * t368;
t358 = t370 * pkin(5);
t98 = qJ(6) * t248 + t158 + t358;
t572 = t98 * t406;
t124 = pkin(5) * t245;
t193 = t276 + t124;
t596 = t193 / 0.2e1;
t71 = t86 * t245;
t408 = (t110 - t194 / 0.2e1) * t370 + (t556 / 0.2e1 - t572 / 0.2e1 + t596) * t368 - t71 / 0.2e1;
t2 = t408 + t420;
t598 = t2 * qJD(1) + t500;
t359 = t470 / 0.2e1 - t388 / 0.2e1;
t597 = t98 / 0.2e1;
t239 = t253 / 0.2e1;
t595 = -t260 / 0.2e1;
t594 = -t376 / 0.2e1;
t591 = -t406 / 0.2e1;
t590 = -t407 / 0.2e1;
t587 = pkin(5) * t366;
t582 = t407 * pkin(5);
t579 = t14 * qJD(4);
t473 = t117 / 0.2e1 + t86 / 0.2e1;
t441 = t473 * t368;
t480 = -t366 / 0.2e1;
t6 = (-t441 + (t480 - 0.1e1 / 0.2e1) * pkin(5)) * t407;
t577 = t6 * qJD(1);
t7 = t118 * t139 + t193 * t194 + t86 * t98;
t576 = t7 * qJD(1);
t475 = t585 / 0.2e1;
t426 = t475 + t473;
t8 = t426 * t407;
t575 = t8 * qJD(1);
t573 = t86 * t407;
t571 = t98 * t407;
t569 = t245 * qJD(3);
t12 = -t117 * t245 - t71;
t567 = qJD(1) * t12;
t545 = t193 * t370;
t24 = t545 + (-t558 + t574) * t368;
t566 = qJD(1) * t24;
t389 = -pkin(4) - t582;
t456 = t370 * t389 / 0.2e1;
t533 = t376 * t407;
t534 = t375 * t406;
t410 = (-t533 / 0.2e1 - t534 / 0.2e1) * t368 + t456;
t557 = t139 * t406;
t422 = t557 / 0.2e1 + t571 / 0.2e1;
t25 = t410 - t422;
t565 = qJD(1) * t25;
t538 = t276 * t370;
t32 = t538 + (-t553 - t554) * t368;
t564 = qJD(1) * t32;
t37 = t118 * t406 + t573;
t563 = qJD(1) * t37;
t431 = t148 * t407 - t149 * t406;
t562 = qJD(1) * t431;
t64 = -t149 * t368 + t276 * t250;
t561 = qJD(1) * t64;
t10 = (t557 + t571) * t370 - t37 * t368;
t560 = t10 * qJD(1);
t11 = t193 * pkin(5) * t250 - t602 * t118;
t559 = t11 * qJD(1);
t15 = t426 * t406;
t552 = t15 * qJD(1);
t550 = t158 * t407;
t549 = t159 * t406;
t542 = t194 * t406;
t544 = t193 * t406;
t17 = (t86 + t542) * t370 + (t98 - t544) * t368;
t547 = t17 * qJD(1);
t18 = (t549 + t550) * t370 + t431 * t368;
t546 = t18 * qJD(1);
t543 = t193 * t407;
t541 = t194 * t407;
t20 = (-t118 + t541) * t370 + (-t139 - t543) * t368;
t540 = t20 * qJD(1);
t27 = (-t148 + t528) * t370 + (t158 - t537) * t368;
t539 = t27 * qJD(1);
t28 = (-t149 + t525) * t370 + (-t159 - t253) * t368;
t536 = t28 * qJD(1);
t31 = t37 * t370;
t535 = t31 * qJD(1);
t439 = t448 * t368;
t413 = pkin(8) * t439 - t583 / 0.2e1;
t421 = t550 / 0.2e1 + t549 / 0.2e1;
t38 = t413 - t421;
t532 = t38 * qJD(1);
t531 = t389 * t406;
t530 = t389 * t407;
t40 = -t118 * t368 + (t406 * t587 + t545) * t407;
t529 = t40 * qJD(1);
t346 = t402 * t368;
t347 = t403 * t368;
t526 = t406 * t407;
t41 = t117 * t368 - t193 * t245 + t403 * t587;
t523 = t41 * qJD(1);
t63 = t148 * t368 - t276 * t245;
t522 = t63 * qJD(1);
t173 = -0.1e1 / 0.2e1 - t442;
t185 = t173 * t407;
t521 = -t185 * qJD(2) + t243 * qJD(3);
t186 = t365 * t590 + (t480 + 0.1e1 / 0.2e1) * t407;
t520 = t186 * qJD(2);
t258 = t600 * t407;
t510 = qJD(2) * t368;
t519 = -t258 * qJD(3) - t406 * t510;
t518 = t239 - t259 / 0.2e1;
t196 = t600 * t406;
t517 = -t196 * qJD(3) + t407 * t510;
t385 = t403 - t402;
t125 = -t277 * t368 + t538;
t516 = qJD(1) * t125;
t195 = t477 * t406;
t515 = qJD(1) * t195;
t514 = qJD(1) * t196;
t197 = t477 * t407;
t513 = qJD(1) * t197;
t367 = t379 * t581;
t512 = qJD(1) * t367;
t511 = qJD(1) * t387;
t509 = qJD(3) * t407;
t507 = qJD(4) * t407;
t506 = qJD(5) * t118;
t505 = qJD(5) * t376;
t504 = qJD(5) * t406;
t395 = qJD(5) * t407;
t503 = qJD(6) * t406;
t502 = qJD(6) * t407;
t418 = t448 * t365 + t480;
t106 = t418 + t448;
t501 = t106 * qJD(1);
t162 = t173 * t406;
t91 = t162 * qJD(1);
t499 = t173 * qJD(1);
t120 = t185 * qJD(1);
t443 = t448 * t370;
t189 = t443 - t359;
t498 = t189 * qJD(1);
t497 = t477 * qJD(1);
t449 = t402 / 0.2e1 + t592;
t242 = t449 * t370;
t496 = t242 * qJD(5);
t495 = t243 * qJD(1);
t494 = t245 * qJD(1);
t222 = t248 * qJD(1);
t254 = t346 + t347;
t493 = t254 * qJD(1);
t255 = t384 * t366;
t492 = t255 * qJD(1);
t257 = t384 * t370;
t491 = t257 * qJD(1);
t150 = t258 * qJD(1);
t490 = t600 * qJD(1);
t489 = t359 * qJD(1);
t488 = t368 * qJD(1);
t352 = t368 * qJD(4);
t355 = t370 * qJD(4);
t450 = -t398 / 0.2e1 - t399 / 0.2e1;
t374 = -0.1e1 / 0.2e1 + t450;
t486 = t374 * qJD(1);
t485 = t379 * qJD(1);
t484 = t384 * qJD(4);
t483 = t404 * qJD(1);
t482 = t405 * qJD(1);
t476 = pkin(5) * t504;
t468 = t406 * t488;
t467 = t406 * t507;
t466 = t368 * t395;
t465 = t370 * t504;
t464 = t370 * t395;
t463 = t370 * t502;
t462 = t368 * t487;
t287 = t368 * t355;
t386 = t406 * t395;
t461 = t406 * t487;
t460 = t370 * t503;
t459 = t407 * t488;
t458 = t368 * t502;
t457 = t544 / 0.2e1;
t238 = t537 / 0.2e1;
t455 = t530 / 0.2e1;
t454 = -t243 / 0.2e1;
t452 = t248 / 0.2e1;
t273 = t368 * t340;
t157 = qJD(4) * t243 + t273;
t446 = qJD(3) + t511;
t444 = t366 * t386;
t438 = qJD(4) * t447;
t437 = t162 * qJD(2) + t368 * t509;
t436 = t370 * t509;
t435 = -t245 * qJD(4) - t466;
t434 = t273 + t464;
t19 = -t148 * t158 + t149 * t159 + t276 * t277;
t433 = t19 * qJD(1) + t14 * qJD(2);
t187 = pkin(5) * t531;
t3 = t473 * t376 + (t597 - t389 * t250 / 0.2e1 - t544 / 0.2e1) * pkin(5);
t432 = -qJD(1) * t3 + qJD(4) * t187;
t430 = t548 - t551;
t298 = t533 + t534;
t412 = (t449 * pkin(5) + t455) * t370;
t30 = -t358 + t595 + (t596 - t276 / 0.2e1) * t406 + (qJ(6) * t590 + t594) * t368 + t412;
t345 = pkin(5) * t526 - t531;
t429 = qJD(1) * t30 - qJD(4) * t345;
t415 = -pkin(5) * t474 - t543 / 0.2e1 + t406 * t456;
t35 = (qJ(6) * t591 - t375 / 0.2e1) * t368 + t415 + t518;
t364 = t402 * pkin(5) + t530;
t428 = -t35 * qJD(1) + t364 * qJD(4);
t427 = t584 / 0.2e1 + t583 / 0.2e1;
t419 = t427 * t407;
t59 = t595 - t419;
t425 = pkin(4) * t508 - t59 * qJD(1);
t411 = t427 * t406 + t239;
t57 = -t253 / 0.2e1 + t259 / 0.2e1 + t411;
t424 = pkin(4) * t507 - t57 * qJD(1);
t132 = (-qJD(5) - t488) * t245;
t207 = -t242 * qJD(1) + t467;
t305 = t461 - t507;
t198 = qJD(5) * t359 + t462;
t170 = qJD(1) * t366 * t526 + t242 * qJD(4);
t256 = t385 * t366;
t176 = t256 * qJD(1) + t438;
t274 = qJD(1) * t447 - t385 * qJD(4);
t414 = t471 / 0.2e1 + t469 / 0.2e1;
t190 = t439 + t414;
t409 = t472 / 0.2e1 + pkin(5) * t454 + t423 / 0.2e1;
t22 = -t409 + t603;
t416 = qJD(1) * t22 - qJD(2) * t190 + qJD(4) * t298;
t401 = qJ(2) * qJD(2);
t400 = qJD(1) * qJ(2);
t377 = t385 * qJD(5);
t373 = 0.1e1 / 0.2e1 + t450;
t342 = t359 * qJD(4);
t339 = t407 * t355;
t306 = -t395 - t459;
t299 = t307 * pkin(5);
t252 = t260 / 0.2e1;
t235 = t254 * qJD(3);
t233 = t257 * qJD(2);
t231 = t257 * qJD(4);
t229 = t254 * qJD(4);
t214 = t243 * qJD(5);
t209 = t395 + t222;
t208 = -t495 - t504;
t191 = t347 / 0.2e1 + t346 / 0.2e1 + t414;
t188 = t443 + t359;
t184 = -t403 * t287 - t444;
t183 = -t402 * t287 + t444;
t174 = 0.2e1 * t407 * t132;
t172 = 0.1e1 / 0.2e1 - t442;
t163 = t442 * t406 + t591;
t161 = t163 * qJD(2);
t156 = t403 * t462 - t496;
t155 = -qJD(4) * t248 + t368 * t461;
t154 = t402 * t462 + t496;
t153 = -t256 * qJD(5) + t368 * t438;
t141 = qJD(4) * t197 - t368 * t465;
t140 = -t195 * qJD(4) - t368 * t464;
t136 = qJD(5) * t248 - t513;
t134 = -t214 + t515;
t128 = t131 * qJD(4);
t127 = -t496 + (-t403 * t487 - t467) * t368;
t126 = t496 + (-t402 * t487 + t467) * t368;
t121 = (t346 - t347) * qJD(4) + (-qJD(5) + t488) * t447;
t119 = t435 + t150;
t105 = t418 - t448;
t103 = t106 * qJD(2);
t101 = t106 * qJD(3);
t100 = t105 * qJD(2);
t99 = t105 * qJD(3);
t90 = t406 * t355 + t513;
t89 = t339 - t515;
t88 = -t214 + t339 + t514;
t81 = t435 + t120;
t74 = -t245 * qJD(5) - t407 * t352;
t72 = -qJD(4) * t250 + t368 * t504 - t91;
t62 = t185 * qJD(5) - t459;
t61 = t186 * qJD(5) + t459;
t60 = 0.2e1 * t238 + t252 - t419;
t58 = t411 + t518;
t47 = -qJD(5) * t250 + t406 * t352;
t39 = t413 + t421;
t36 = t375 * t368 / 0.2e1 + qJ(6) * t454 - t415 + t518;
t34 = -qJD(5) * t162 + t468;
t33 = qJD(5) * t163 - t468;
t29 = qJ(6) * t452 + t368 * t594 + t238 + t252 + t358 + t412 + t457;
t26 = t410 + t422;
t21 = t409 + t603;
t16 = t117 * t591 - t558 / 0.2e1 + t406 * t475 + t570;
t9 = t117 * t590 - t573 / 0.2e1 + pkin(5) * t452;
t5 = t582 / 0.2e1 + (pkin(5) * t480 - t441) * t407;
t4 = t455 * t358 + t602 * t594 + (t457 + t597) * pkin(5);
t1 = t408 - t420;
t13 = [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, qJD(2), t401, 0, 0, 0, 0, 0, 0, qJD(2) * t404, qJD(2) * t405, t379 * qJD(3), -qJD(3) * t367 + t401, -t287, -t477 * qJD(4), 0, t287, 0, 0, t387 * t355 + t510, qJD(2) * t370 - t387 * t352, -t600 * qJD(3), qJD(2) * t387 + qJD(3) * t125, t184, t153, t141, t183, t140, t287, qJD(4) * t27 + qJD(5) * t64 + t517, qJD(4) * t28 + qJD(5) * t63 + t519, -qJD(4) * t18 - t233, -qJD(2) * t431 + qJD(3) * t32 + qJD(4) * t19, t184, t153, t141, t183, t140, t287, t17 * qJD(4) + t40 * qJD(5) - t370 * t458 + t517, qJD(4) * t20 + qJD(5) * t41 + t368 * t460 + t519, -qJD(4) * t10 - qJD(5) * t12 + qJD(6) * t255 - t233, qJD(2) * t37 + qJD(3) * t24 + qJD(4) * t7 + qJD(5) * t11 - qJD(6) * t31; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, qJD(1), t400, 0, 0, 0, 0, 0, 0, t483, t482, 0, qJD(3) * t373 + t400, 0, 0, 0, 0, 0, 0, t488, t487, 0, qJD(3) * t172 + t511, 0, 0, 0, 0, 0, 0, t61, t33, -t491, t99 - t562 + t579, 0, 0, 0, 0, 0, 0, t61, t33, -t491, qJD(4) * t1 + qJD(5) * t5 + t563 + t99; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t485, qJD(2) * t373 - t512, 0, 0, 0, 0, 0, 0, 0, 0, -t490, qJD(2) * t172 + t516, 0, 0, 0, 0, 0, 0, -t514, -t150, 0, qJD(4) * t39 + t100 + t564, 0, 0, 0, 0, 0, 0, -t514, -t150, 0, qJD(4) * t26 + qJD(5) * t16 + qJD(6) * t188 + t100 + t566; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t462, -t497, -t352, t462, -t355, 0, -qJD(4) * t277 + t387 * t487, qJD(4) * t276 - t387 * t488, 0, 0, t127, t121, t90, t126, t89, t198, t539 + (t406 * t440 - t525) * qJD(4) + t60 * qJD(5), t536 + (t407 * t440 + t528) * qJD(4) + t58 * qJD(5), qJD(4) * t430 - t546, t39 * qJD(3) + (-t277 * pkin(4) + pkin(8) * t430) * qJD(4) + t433, t127, t121, t90, t126, t89, t198, t547 + (-t243 * t389 - t375 * t370 - t541) * qJD(4) + t29 * qJD(5) - t243 * qJD(6), t540 + (-t248 * t389 - t376 * t370 + t542) * qJD(4) + t36 * qJD(5) - t458, -t560 + (t556 - t572 + (-t524 + t527) * t368) * qJD(4) + t9 * qJD(5), t576 + t1 * qJD(2) + t26 * qJD(3) + (t139 * t376 + t194 * t389 - t375 * t98) * qJD(4) + t4 * qJD(5) + t21 * qJD(6); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t170, -t176, t132, t170, -t434, t342, qJD(4) * t60 - qJD(5) * t149 + t520 + t561, qJD(4) * t58 + qJD(5) * t148 + t161 + t522, 0, 0, -t170, -t176, t132, t170, -t434, t342, qJD(4) * t29 - t506 + t520 + t529, qJD(4) * t36 + qJD(5) * t117 + t161 + t523, pkin(5) * t465 + qJD(4) * t9 - t567, -pkin(5) * t506 + qJD(2) * t5 + qJD(3) * t16 + qJD(4) * t4 + t559; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t157, t305 * t368, t492, qJD(3) * t188 + qJD(4) * t21 - t535; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -qJD(1), -t400, 0, 0, 0, 0, 0, 0, -t483, -t482, 0, qJD(3) * t374 - t400, 0, 0, 0, 0, 0, 0, -t488, -t487, 0, qJD(3) * t173 - t511, 0, 0, 0, 0, 0, 0, t62, t34, t491, t101 + t562 + t579, 0, 0, 0, 0, 0, 0, t62, t34, t491, qJD(4) * t2 + qJD(5) * t6 + t101 - t563; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t128, 0, 0, 0, 0, 0, 0, 0, 0, 0, t128; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t486, 0, 0, 0, 0, 0, 0, 0, 0, 0, t499, 0, 0, 0, 0, 0, 0, 0, 0, 0, t501, 0, 0, 0, 0, 0, 0, 0, 0, 0, t501; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t352, -t355, 0, 0, 0, 0, 0, 0, 0, 0, t74, t47, t231 (pkin(8) * t257 - t586) * qJD(4) - t599, 0, 0, 0, 0, 0, 0, t74, t47, t231 (t298 * t370 + t368 * t389) * qJD(4) - t124 * qJD(5) + t191 * qJD(6) + t598; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t81, t72, 0, 0, 0, 0, 0, 0, 0, 0, t81, t72, 0, -pkin(5) * t466 - t124 * qJD(4) + t577; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t191 * qJD(4); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t485, -qJD(2) * t374 + t512, 0, 0, 0, 0, 0, 0, t355, -t352, t490, -qJD(2) * t173 - t516, 0, 0, 0, 0, 0, 0, t88, t119, t229, -qJD(4) * t38 - t103 - t564, 0, 0, 0, 0, 0, 0, t88, t119, t229, -qJD(4) * t25 - qJD(5) * t15 + qJD(6) * t189 - t103 - t566; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t486, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t499, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t501, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t501; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t487, -t488, 0, 0, 0, 0, 0, 0, 0, 0, t340, -t494, t493, -t532, 0, 0, 0, 0, 0, 0, t340, -t494, t493, -t565; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t208, t306, 0, 0, 0, 0, 0, 0, 0, 0, t208, t306, 0, -t476 - t552; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t498; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t462, t497, 0, -t462, 0, 0, -t446 * t370, t446 * t368, 0, 0, t156, t174, t136, t154, t134, -t198, t59 * qJD(5) - t436 - t539, qJD(5) * t57 - t536 + t569, -t235 + t546, qJD(3) * t38 - t433, t156, t174, t136, t154, t134, -t198, t30 * qJD(5) - t436 - t547, -qJD(5) * t35 - t540 + t569, -qJD(5) * t8 - t235 + t560, -qJD(2) * t2 + qJD(3) * t25 - qJD(5) * t3 + qJD(6) * t22 - t576; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t599, 0, 0, 0, 0, 0, 0, 0, 0, 0, -qJD(6) * t190 - t598; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t487, t488, 0, 0, 0, 0, 0, 0, 0, 0, -t340, t494, -t493, t532, 0, 0, 0, 0, 0, 0, -t340, t494, -t493, t565; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t386, t377, 0, -t386, 0, 0, -pkin(4) * t504, -pkin(4) * t395, 0, 0, t386, t377, 0, -t386, 0, 0, -t345 * qJD(5), t364 * qJD(5), t384 * qJD(6), qJD(5) * t187 + qJD(6) * t298; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t207, -t274, t209, -t207, t208, -t489, -pkin(8) * t395 - t425, pkin(8) * t504 - t424, 0, 0, t207, -t274, t209, -t207, t208, -t489, t429 - t505, qJD(5) * t375 + t428, -pkin(5) * t395 - t575, -pkin(5) * t505 + t432; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t484, t416; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t170, t176, t155, -t170, t157, t342, -qJD(4) * t59 + t521 - t561, -t57 * qJD(4) + t437 - t522, 0, 0, t170, t176, t155, -t170, t157, t342, -t30 * qJD(4) - t463 + t521 - t529, t35 * qJD(4) + t437 + t460 - t523, qJD(4) * t8 + t567, -pkin(5) * t463 - t6 * qJD(2) + t15 * qJD(3) + t3 * qJD(4) - t559; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t120, t91, 0, 0, 0, 0, 0, 0, 0, 0, -t120, t91, 0, -t577; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t495, t459, 0, 0, 0, 0, 0, 0, 0, 0, t495, t459, 0, t552; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t207, t274, -t222, t207, t495, t489, t425, t424, 0, 0, -t207, t274, -t222, t207, t495, t489, -t429 - t503, -t428 - t502, t575, -pkin(5) * t503 - t432; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t307, t305, 0, -t299; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t434, t132, -t492, pkin(5) * t464 - t189 * qJD(3) - t22 * qJD(4) + t535; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t190 * qJD(4); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t498; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t504, t395, -t484, -t416 + t476; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t307, -t305, 0, t299; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0;];
cmat_reg  = t13;
