% Calculate minimal parameter regressor of coriolis matrix for
% S6PRRPRP1
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [11x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d2,d3,d5,theta1,theta4]';
% 
% Output:
% cmat_reg [(6*%NQJ)%x26]
%   minimal parameter regressor of coriolis matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2021-01-16 02:40
% Revision: 24b2e7d74a0c1a3b64fa2f8f5ad758691ad61af3 (2021-01-15)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function cmat_reg = S6PRRPRP1_coriolismatJ_fixb_regmin_slag_vp(qJ, qJD, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(11,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRRPRP1_coriolismatJ_fixb_regmin_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6PRRPRP1_coriolismatJ_fixb_regmin_slag_vp: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6PRRPRP1_coriolismatJ_fixb_regmin_slag_vp: pkin has to be [11x1] (double)');

%% Symbolic Calculation
% From coriolismat_joint_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2021-01-16 02:38:55
% EndTime: 2021-01-16 02:39:19
% DurationCPUTime: 9.12s
% Computational Cost: add. (8422->512), mult. (18941->731), div. (0->0), fcn. (21647->10), ass. (0->418)
t374 = sin(qJ(3));
t604 = qJ(4) + pkin(8);
t354 = t604 * t374;
t377 = cos(qJ(3));
t355 = t604 * t377;
t371 = sin(pkin(11));
t597 = cos(pkin(11));
t281 = t597 * t354 + t371 * t355;
t376 = cos(qJ(5));
t635 = t281 * t376;
t638 = t635 / 0.2e1;
t372 = sin(pkin(6));
t375 = sin(qJ(2));
t553 = t372 * t375;
t598 = cos(pkin(6));
t328 = t374 * t553 - t598 * t377;
t329 = t598 * t374 + t377 * t553;
t217 = t597 * t328 + t371 * t329;
t637 = t217 / 0.2e1;
t448 = t637 - t217 / 0.2e1;
t378 = cos(qJ(2));
t552 = t372 * t378;
t456 = t552 / 0.2e1;
t554 = t371 * t377;
t636 = t554 / 0.2e1;
t373 = sin(qJ(5));
t349 = t597 * t355;
t556 = t371 * t354;
t631 = t349 - t556;
t570 = t631 * t373;
t569 = t631 * t376;
t446 = t597 * t374;
t340 = t446 + t554;
t549 = t373 * t376;
t483 = t340 * t549;
t444 = 0.2e1 * t483;
t366 = -t377 * pkin(3) - pkin(2);
t445 = t597 * t377;
t555 = t371 * t374;
t403 = t445 - t555;
t398 = -pkin(4) * t403 - t340 * pkin(9) + t366;
t147 = -t376 * t398 + t570;
t548 = t376 * qJ(6);
t130 = t340 * t548 + t147;
t607 = t403 * pkin(5);
t112 = -t130 - t607;
t634 = t112 + t130;
t302 = t597 * t329;
t557 = t371 * t328;
t633 = t302 - t557;
t497 = t340 * qJD(2);
t323 = t376 * t497;
t492 = t373 * qJD(3);
t294 = t323 + t492;
t335 = t403 ^ 2;
t336 = t340 ^ 2;
t632 = -t336 - t335;
t365 = -t597 * pkin(3) - pkin(4);
t351 = -t376 * pkin(5) + t365;
t560 = t351 * t373;
t295 = pkin(5) * t549 - t560;
t126 = t448 * t373;
t510 = t126 * qJD(1);
t364 = t371 * pkin(3) + pkin(9);
t532 = qJ(6) + t364;
t332 = t532 * t376;
t369 = t373 ^ 2;
t370 = t376 ^ 2;
t612 = -t370 / 0.2e1;
t447 = t369 / 0.2e1 + t612;
t559 = t351 * t376;
t460 = -t559 / 0.2e1;
t248 = t373 * t340;
t206 = pkin(5) * t248 + t281;
t580 = t206 * t373;
t466 = -t580 / 0.2e1;
t605 = t374 * pkin(3);
t251 = t340 * pkin(4) - pkin(9) * t403 + t605;
t227 = t376 * t251;
t551 = t281 * t373;
t529 = t227 / 0.2e1 + t551 / 0.2e1;
t52 = t466 - (t548 / 0.2e1 + t332 / 0.2e1) * t403 + (t460 + (0.1e1 - t447) * pkin(5)) * t340 + t529;
t630 = -t52 * qJD(2) - t295 * qJD(3) - t510;
t610 = t373 / 0.2e1;
t611 = -t373 / 0.2e1;
t128 = (t610 - t611) * t217;
t490 = t376 * qJD(3);
t127 = t448 * t376;
t509 = t127 * qJD(1);
t614 = -t340 / 0.2e1;
t615 = -t403 / 0.2e1;
t410 = t364 * t615 + t365 * t614;
t383 = t410 * t373 + t638;
t226 = t373 * t251;
t530 = -t226 / 0.2e1 + t638;
t80 = t383 - t530;
t629 = -t80 * qJD(2) - t365 * t490 + t509;
t313 = t369 * pkin(5) + t559;
t330 = t532 * t373;
t613 = t340 / 0.2e1;
t462 = t351 * t613;
t579 = t206 * t376;
t390 = -pkin(5) * t483 - t579 / 0.2e1 + t373 * t462;
t455 = qJ(6) * t611;
t59 = -(t455 - t330 / 0.2e1) * t403 + t390 + t530;
t628 = t59 * qJD(2) - t313 * qJD(3) + t509;
t392 = t410 * t376;
t82 = -t227 / 0.2e1 - t392;
t627 = -t82 * qJD(2) - t365 * t492 + t510;
t526 = t446 * t456 + t552 * t636;
t331 = t446 / 0.2e1 + t636;
t626 = pkin(5) / 0.2e1;
t250 = t376 * t403;
t606 = t340 * pkin(5);
t113 = -qJ(6) * t250 + t227 + t551 + t606;
t625 = t113 / 0.2e1;
t623 = t633 / 0.2e1;
t287 = t403 * t552;
t550 = t373 * t287;
t228 = t376 * t553 - t550;
t618 = t228 / 0.2e1;
t617 = -t332 / 0.2e1;
t616 = t403 / 0.2e1;
t609 = -t376 / 0.2e1;
t608 = pkin(5) * t336;
t603 = qJD(3) * pkin(3);
t185 = -t373 * t552 + t376 * t633;
t436 = t448 * t403;
t389 = t613 * t633 - t436;
t380 = t185 * t614 + t389 * t376;
t286 = t340 * t552;
t567 = t286 * t373;
t46 = -t567 / 0.2e1 + t380;
t544 = t46 * qJD(1);
t600 = t248 * qJD(4) - t544;
t245 = t373 * t403;
t458 = t553 / 0.2e1;
t576 = t217 * t340;
t465 = -t576 / 0.2e1;
t409 = t458 + t465;
t453 = -t550 / 0.2e1;
t585 = t185 * t403;
t71 = t453 - t585 / 0.2e1 + t409 * t376;
t537 = t71 * qJD(1);
t599 = -t245 * qJD(4) + t537;
t562 = t332 * t376;
t565 = t330 * t373;
t384 = -(-t562 / 0.2e1 - t565 / 0.2e1) * t403 + t462;
t528 = -t226 + t635;
t133 = -qJ(6) * t245 - t528;
t589 = t133 * t373;
t591 = t113 * t376;
t415 = -t591 / 0.2e1 - t589 / 0.2e1;
t34 = t384 + t415;
t596 = qJD(2) * t34;
t148 = t373 * t398 + t569;
t561 = t340 * t376;
t94 = t148 * t403 + t281 * t561;
t595 = qJD(2) * t94;
t449 = t130 / 0.2e1 + t112 / 0.2e1;
t484 = -t607 / 0.2e1;
t416 = t484 + t449;
t10 = t416 * t376;
t594 = t10 * qJD(2);
t593 = t112 * t373;
t592 = t112 * t376;
t131 = -qJ(6) * t248 + t148;
t590 = t131 * t376;
t15 = t634 * t248;
t588 = t15 * qJD(2);
t184 = t373 * t633 + t376 * t552;
t587 = t184 * t403;
t586 = t184 * t373;
t584 = t185 * t376;
t205 = pkin(5) * t245 + t631;
t583 = t205 * t373;
t582 = t205 * t376;
t581 = t206 * t340;
t577 = t217 * t286;
t575 = t633 * t403;
t574 = t228 * t376;
t482 = t373 * t553;
t547 = t376 * t287;
t229 = t482 + t547;
t573 = t229 * t373;
t25 = (t112 + t583) * t340 - (t113 - t580) * t403;
t572 = t25 * qJD(2);
t28 = (-t131 + t582) * t340 - (-t133 - t579) * t403;
t571 = t28 * qJD(2);
t566 = t286 * t376;
t564 = t330 * t376;
t563 = t332 * t373;
t558 = t369 * t340;
t38 = (t633 - t584 - t586) * t217;
t546 = t38 * qJD(1);
t381 = t184 * t614 + t389 * t373;
t43 = t566 / 0.2e1 + t381;
t545 = t43 * qJD(1);
t49 = (-t147 + t570) * t340 - t227 * t403;
t543 = t49 * qJD(2);
t50 = (-t148 + t569) * t340 - (t528 - t635) * t403;
t542 = t50 * qJD(2);
t51 = -t184 * t228 + t185 * t229 + t577;
t541 = t51 * qJD(1);
t58 = (-t633 / 0.2e1 + t623) * t340 - t436;
t540 = t58 * qJD(2);
t62 = t131 * t403 + (t373 * t608 + t581) * t376;
t539 = t62 * qJD(2);
t63 = -t130 * t403 - t206 * t248 + t370 * t608;
t538 = t63 * qJD(2);
t451 = -t547 / 0.2e1;
t72 = t451 + t587 / 0.2e1 - t409 * t373;
t536 = t72 * qJD(1);
t91 = -t372 ^ 2 * t378 * t375 + t287 * t633 + t577;
t534 = t91 * qJD(1);
t93 = -t147 * t403 - t248 * t281;
t533 = t93 * qJD(2);
t531 = -t593 / 0.2e1 + t590 / 0.2e1;
t359 = t369 + t370;
t360 = t370 - t369;
t488 = t336 - t335;
t192 = t488 * t373;
t524 = qJD(2) * t192;
t193 = t632 * t373;
t523 = qJD(2) * t193;
t522 = qJD(2) * t375;
t521 = qJD(2) * t377;
t520 = qJD(3) * t217;
t519 = qJD(3) * t374;
t518 = qJD(3) * t377;
t517 = qJD(3) * t378;
t172 = t584 / 0.2e1;
t467 = -t584 / 0.2e1;
t103 = t172 + t467;
t516 = qJD(4) * t103;
t515 = qJD(4) * t403;
t514 = qJD(5) * t131;
t513 = qJD(5) * t185;
t512 = qJD(5) * t332;
t511 = qJD(5) * t373;
t368 = qJD(5) * t376;
t459 = t340 * t612;
t439 = -t558 / 0.2e1 + t459;
t188 = t439 - t331;
t508 = t188 * qJD(2);
t194 = t488 * t376;
t507 = t194 * qJD(2);
t481 = -t597 / 0.2e1;
t396 = t340 * t481 + t371 * t616;
t209 = (-t374 / 0.2e1 + t396) * pkin(3);
t506 = t209 * qJD(2);
t505 = t245 * qJD(2);
t504 = t248 * qJD(2);
t503 = t250 * qJD(2);
t326 = t369 * t403;
t327 = t370 * t403;
t252 = -t326 - t327;
t502 = t252 * qJD(2);
t253 = t359 * t336;
t501 = t253 * qJD(2);
t255 = t632 * t376;
t138 = t255 * qJD(2);
t500 = t632 * qJD(2);
t499 = t331 * qJD(2);
t498 = t403 * qJD(2);
t496 = t340 * qJD(3);
t495 = t340 * qJD(4);
t494 = t359 * qJD(3);
t361 = -t374 ^ 2 + t377 ^ 2;
t493 = t361 * qJD(2);
t491 = t373 * qJD(6);
t489 = t376 * qJD(6);
t487 = pkin(2) * t374 * qJD(2);
t486 = pkin(2) * t521;
t485 = pkin(5) * t511;
t480 = t340 * t511;
t479 = t340 * t368;
t478 = t403 * t497;
t477 = t403 * t496;
t476 = t372 * t522;
t475 = qJD(2) * t552;
t474 = t373 * t368;
t473 = t373 * t497;
t472 = t340 * t491;
t471 = t373 * t490;
t470 = t374 * t521;
t469 = t376 * t498;
t468 = t340 * t489;
t464 = t576 / 0.2e1;
t461 = t561 / 0.2e1;
t457 = -t552 / 0.2e1;
t452 = -t248 / 0.2e1;
t450 = -t250 / 0.2e1;
t266 = t403 * t323;
t443 = -qJD(3) * t245 - t266;
t442 = t373 * t484;
t440 = t217 * t461;
t438 = -t376 * t515 + t536;
t437 = -t376 * t495 - t545;
t435 = qJD(3) * t444;
t434 = t302 / 0.2e1 - t557 / 0.2e1;
t433 = -t266 + t479;
t379 = -t531 * t217 - t184 * t113 / 0.2e1 + t185 * t133 / 0.2e1 + t206 * t623 + t205 * t637;
t387 = t330 * t618 + t229 * t617 - t286 * t351 / 0.2e1;
t2 = t379 + t387;
t9 = t112 * t113 + t131 * t133 + t205 * t206;
t432 = t2 * qJD(1) + t9 * qJD(2);
t16 = t206 * pkin(5) * t561 - t634 * t131;
t3 = t449 * t185 + (t376 * t465 + t618) * pkin(5);
t431 = -t3 * qJD(1) + t16 * qJD(2);
t425 = t131 * t373 + t592;
t12 = (t589 + t591) * t340 + t425 * t403;
t414 = t184 * t609 + t185 * t610;
t393 = t414 * t403;
t411 = t228 * t610 + t229 * t609;
t56 = -t393 + t411;
t430 = t56 * qJD(1) - t12 * qJD(2);
t382 = t448 * t631;
t397 = t287 * t371 / 0.2e1 + t286 * t481;
t26 = (t374 * t456 + t397) * pkin(3) + t382;
t95 = t366 * t605;
t429 = -t26 * qJD(1) + t95 * qJD(2);
t35 = t581 - (-t590 + t593) * t403;
t413 = t467 - t586 / 0.2e1;
t385 = -t403 * t413 + t464;
t412 = -t574 / 0.2e1 - t573 / 0.2e1;
t37 = t385 + t412;
t428 = qJD(1) * t37 + qJD(2) * t35;
t54 = t425 * t340;
t64 = t340 * t414 + t526;
t427 = -t64 * qJD(1) - t54 * qJD(2);
t426 = t58 * qJD(1);
t242 = t562 + t565;
t424 = -t340 * t364 + t365 * t403;
t19 = t416 * t373;
t423 = -t103 * qJD(1) + t19 * qJD(2);
t136 = t281 * t340 + t403 * t631;
t98 = -t575 / 0.2e1 + t409;
t422 = -qJD(1) * t98 + qJD(2) * t136;
t195 = t340 * t456 - t526;
t224 = t366 * t340 - t403 * t605;
t418 = -t195 * qJD(1) + t224 * qJD(2);
t395 = -t445 / 0.2e1 + t555 / 0.2e1;
t196 = (t616 + t395) * t552;
t225 = t340 * t605 + t366 * t403;
t417 = -t196 * qJD(1) + t225 * qJD(2);
t406 = (-qJD(5) + t498) * t248;
t244 = t447 * t340;
t405 = -t244 * qJD(2) + t471;
t292 = t473 - t490;
t404 = t331 * qJD(5) - t478;
t402 = t340 * t460 + t466;
t401 = t336 * qJD(2) * t549 + t244 * qJD(3);
t254 = t360 * t336;
t400 = t254 * qJD(2) + t435;
t399 = qJD(2) * t444 - t360 * qJD(3);
t137 = pkin(5) * t560;
t21 = pkin(5) * t126;
t5 = t449 * t332 + (t625 + t402) * pkin(5);
t394 = -t21 * qJD(1) - t5 * qJD(2) + t137 * qJD(3);
t386 = (-t563 / 0.2e1 + t564 / 0.2e1) * t340 + t531;
t388 = -t349 / 0.2e1 + t556 / 0.2e1 + t442;
t32 = t386 + t388;
t78 = t413 + t434;
t391 = -t78 * qJD(1) + t32 * qJD(2) + t242 * qJD(3);
t325 = t331 * qJD(3);
t322 = t340 * t490;
t293 = -t368 + t469;
t285 = t294 * pkin(5);
t243 = t255 * qJD(4);
t231 = t245 * qJD(5);
t230 = t244 * qJD(5);
t208 = t605 / 0.2e1 + t396 * pkin(3);
t207 = t505 - t511;
t198 = t340 * t457 - t526;
t197 = t395 * t552 - t403 * t456;
t189 = t193 * qJD(4);
t187 = t439 + t331;
t129 = -t217 * t609 + t376 * t637;
t120 = -t248 * qJD(3) + t368 * t403 + t138;
t107 = t231 + t322 + t523;
t101 = t103 * qJD(5);
t99 = t575 / 0.2e1 + t464 + t458;
t83 = t281 * t610 - t392 + t529;
t81 = t383 + t530;
t79 = t172 + t586 / 0.2e1 + t434;
t74 = t585 / 0.2e1 + t440 + t453 + t376 * t458;
t73 = -t587 / 0.2e1 + t217 * t452 + t451 - t482 / 0.2e1;
t65 = t184 * t461 + t185 * t452 + t526;
t60 = t330 * t615 - t403 * t455 - t390 + t530;
t57 = t58 * qJD(3);
t55 = -t393 - t411;
t53 = pkin(5) * t459 + qJ(6) * t450 + t332 * t616 + t558 * t626 - t402 + t529 + t606;
t48 = qJD(2) * t72 + qJD(3) * t127;
t47 = qJD(2) * t71 + qJD(3) * t126;
t45 = t567 / 0.2e1 + t380;
t44 = -t566 / 0.2e1 + t381;
t36 = t385 - t412;
t33 = t384 - t415;
t31 = t386 - t388;
t30 = qJD(2) * t74 + qJD(3) * t128 - t513;
t29 = qJD(2) * t73 + qJD(3) * t129 + qJD(5) * t184;
t27 = pkin(3) * t397 + t457 * t605 - t382;
t24 = -qJD(2) * t43 - qJD(5) * t126;
t23 = -qJD(2) * t46 - qJD(5) * t127;
t22 = t128 * pkin(5);
t20 = t130 * t611 - t590 / 0.2e1 + t442 + t531;
t18 = t45 * qJD(2) + t129 * qJD(5) + t492 * t633;
t17 = t44 * qJD(2) + t128 * qJD(5) - t490 * t633;
t14 = qJD(3) * t46 - qJD(5) * t72;
t13 = qJD(3) * t43 - qJD(5) * t71;
t11 = t130 * t609 - t592 / 0.2e1 + pkin(5) * t450;
t8 = (t229 * t403 + t286 * t561) * qJD(2) + t45 * qJD(3) + t73 * qJD(5);
t7 = (-t228 * t403 + t248 * t286) * qJD(2) + t44 * qJD(3) + t74 * qJD(5);
t6 = t580 * t626 + t634 * t617 + (t351 * t461 + t625) * pkin(5);
t4 = -t634 * t185 / 0.2e1 + (t440 + t618) * pkin(5);
t1 = t379 - t387;
t39 = [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, qJD(2) * t91, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, qJD(2) * t51 + qJD(3) * t38; 0, 0, -t476, -t475, 0, 0, 0, 0, 0, (-t374 * t517 - t375 * t521) * t372, (t374 * t522 - t377 * t517) * t372, t198 * qJD(3) - t403 * t476, t197 * qJD(3) + t340 * t476, (t286 * t340 + t287 * t403) * qJD(2) + t57, t534 + (t286 * t281 + t287 * t631 + t366 * t553) * qJD(2) + t27 * qJD(3) + t99 * qJD(4), 0, 0, 0, 0, 0, t7, t8, t7, t8, t55 * qJD(3) + (-t573 - t574) * t497, t541 + (t112 * t228 + t131 * t229 + t206 * t286) * qJD(2) + t1 * qJD(3) + t36 * qJD(4) + t4 * qJD(5) + t65 * qJD(6); 0, 0, 0, 0, 0, 0, 0, 0, 0, -t329 * qJD(3) - t374 * t475, t328 * qJD(3) - t377 * t475, qJD(2) * t198 - qJD(3) * t633, qJD(2) * t197 + t520, t540, t27 * qJD(2) + (-t217 * t371 - t597 * t633) * t603, 0, 0, 0, 0, 0, t17, t18, t17, t18, t55 * qJD(2) - t359 * t520, t546 + t1 * qJD(2) + (-t217 * t242 + t351 * t633) * qJD(3) + t22 * qJD(5) + t79 * qJD(6); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t99 * qJD(2), 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, qJD(2) * t36 + t101; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t30, t29, t30, t29, 0, -pkin(5) * t513 + qJD(2) * t4 + qJD(3) * t22 + t516; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, qJD(2) * t65 + qJD(3) * t79; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t195 * qJD(3), -t196 * qJD(3), t57, -qJD(3) * t26 - qJD(4) * t98 - t534, 0, 0, 0, 0, 0, t13, t14, t13, t14, qJD(3) * t56, qJD(3) * t2 + qJD(4) * t37 - qJD(5) * t3 - qJD(6) * t64 - t541; 0, 0, 0, 0, t374 * t518, t361 * qJD(3), 0, 0, 0, -pkin(2) * t519, -pkin(2) * t518, t224 * qJD(3), t225 * qJD(3), -qJD(4) * t632, qJD(3) * t95 + qJD(4) * t136, -t336 * t474 + t370 * t477, -t254 * qJD(5) - t403 * t435, t194 * qJD(3) + t403 * t480, -t192 * qJD(3) + t403 * t479, -t477, qJD(3) * t49 + qJD(5) * t94 - t189, qJD(3) * t50 + qJD(5) * t93 - t243, t25 * qJD(3) + t62 * qJD(5) + t403 * t468 - t189, t28 * qJD(3) + t63 * qJD(5) - t403 * t472 - t243, -qJD(3) * t12 + qJD(5) * t15 + qJD(6) * t253, qJD(3) * t9 + qJD(4) * t35 + qJD(5) * t16 - qJD(6) * t54; 0, 0, 0, 0, t470, t493, t518, -t519, 0, -pkin(8) * t518 - t487, pkin(8) * t519 - t486, -qJD(3) * t631 + t418, qJD(3) * t281 + t417, (-t340 * t371 - t403 * t597) * t603 + t426, (-t281 * t371 - t597 * t631) * t603 + t208 * qJD(4) + t429, -t230 - (-t370 * t497 - t471) * t403, (-t326 + t327) * qJD(3) + (-qJD(5) - t498) * t444, t340 * t492 + t507, t322 - t524, t404, t545 + t543 + (t373 * t424 - t569) * qJD(3) + t83 * qJD(5), t542 + (t376 * t424 + t570) * qJD(3) + t81 * qJD(5) + t544, t545 + t572 + (t245 * t351 - t330 * t340 - t582) * qJD(3) + t53 * qJD(5) + t245 * qJD(6), t571 + (t250 * t351 - t332 * t340 + t583) * qJD(3) + t60 * qJD(5) + t403 * t489 + t544, (-t113 * t373 + t133 * t376 - (t563 - t564) * t403) * qJD(3) + t11 * qJD(5) + t430, (-t113 * t330 + t133 * t332 + t205 * t351) * qJD(3) + t33 * qJD(4) + t6 * qJD(5) + t31 * qJD(6) + t432; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t500, qJD(3) * t208 + t422, 0, 0, 0, 0, 0, -t523, -t138, -t523, -t138, 0, qJD(3) * t33 + qJD(5) * t20 + qJD(6) * t187 + t428; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t401, -t400, t406, -t433, t325, qJD(3) * t83 - qJD(5) * t148 - t537 + t595, qJD(3) * t81 + qJD(5) * t147 + t533 - t536, qJD(3) * t53 - t514 - t537 + t539, qJD(3) * t60 + qJD(5) * t130 - t536 + t538, pkin(5) * t480 + t11 * qJD(3) + t588, -pkin(5) * t514 + qJD(3) * t6 + qJD(4) * t20 + t431; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t443, -t292 * t403, t501, qJD(3) * t31 + qJD(4) * t187 + t427; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t195 * qJD(2), t196 * qJD(2), -t540, qJD(2) * t26, 0, 0, 0, 0, 0, t24, t23, t24, t23, -qJD(2) * t56, -qJD(2) * t2 - qJD(5) * t21 - qJD(6) * t78 - t546; 0, 0, 0, 0, -t470, -t493, 0, 0, 0, t487, t486, -t418 - t495, -t417 - t515, -t426, qJD(4) * t209 - t429, -t370 * t478 - t230, 0.2e1 * t376 * t406, -qJD(5) * t250 - t507, t231 + t524, -t404, t82 * qJD(5) + t437 - t543, qJD(5) * t80 - t542 + t600, -t52 * qJD(5) + t437 - t572, -qJD(5) * t59 - t571 + t600, -qJD(4) * t252 - qJD(5) * t10 - t430, qJD(4) * t34 - qJD(5) * t5 + qJD(6) * t32 - t432; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t474, t360 * qJD(5), 0, 0, 0, t365 * t511, t365 * t368, -t295 * qJD(5), t313 * qJD(5), t359 * qJD(6), qJD(5) * t137 + qJD(6) * t242; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t497, -t498, 0, t506, 0, 0, 0, 0, 0, -t323, t504, -t323, t504, -t502, t596; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t405, -t399, t368 - t503, t207, -t499, -t364 * t368 - t627, t364 * t511 - t629, -t512 + t630, qJD(5) * t330 - t628, -pkin(5) * t368 - t594, -pkin(5) * t512 + t394; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t494, t391; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t98 * qJD(2), 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -qJD(2) * t37 + t101; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t496, t403 * qJD(3), t500, -qJD(3) * t209 - t422, 0, 0, 0, 0, 0, t107, t120, t107, t120, t252 * qJD(3), -qJD(3) * t34 - qJD(5) * t19 + qJD(6) * t188 - t428; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t497, t498, 0, -t506, 0, 0, 0, 0, 0, t323, -t504, t323, -t504, t502, -t596; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t207, t293, t207, t293, 0, -t423 - t485; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t508; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t47, t48, t47, t48, 0, qJD(2) * t3 + qJD(3) * t21 - t516; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t401, t400, t250 * qJD(3) - t403 * t473, t443, t325, -qJD(3) * t82 - t595 + t599, -t80 * qJD(3) + t438 - t533, t52 * qJD(3) - t468 - t539 + t599, t59 * qJD(3) + t438 + t472 - t538, qJD(3) * t10 - t588, -pkin(5) * t468 + t5 * qJD(3) + t19 * qJD(4) - t431; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t405, t399, t503, -t505, t499, t627, t629, -t491 - t630, -t489 + t628, t594, -pkin(5) * t491 - t394; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t505, -t469, -t505, -t469, 0, t423; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t294, t292, 0, -t285; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, qJD(2) * t64 + qJD(3) * t78; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t433, t406, -t501, pkin(5) * t479 - t32 * qJD(3) - t188 * qJD(4) - t427; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t511, t368, -t494, -t391 + t485; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t508; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t294, -t292, 0, t285; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0;];
cmat_reg = t39;
