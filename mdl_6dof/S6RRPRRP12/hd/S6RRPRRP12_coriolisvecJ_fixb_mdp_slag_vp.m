% Calculate minimal parameter regressor of Coriolis joint torque vector for
% S6RRPRRP12
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d4,d5]';
% MDP [32x1]
%   Minimal dynamic parameter vector (fixed base model)
%   see S6RRPRRP12_convert_par2_MPV_fixb.m
% 
% Output:
% tauc [6x1]
%   joint torques required to compensate Coriolis and centrifugal load

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 12:54
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauc = S6RRPRRP12_coriolisvecJ_fixb_mdp_slag_vp(qJ, qJD, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(9,1),zeros(32,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPRRP12_coriolisvecJ_fixb_mdp_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRPRRP12_coriolisvecJ_fixb_mdp_slag_vp: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S6RRPRRP12_coriolisvecJ_fixb_mdp_slag_vp: pkin has to be [9x1] (double)');
assert(isreal(MDP) && all(size(MDP) == [32 1]), ...
  'S6RRPRRP12_coriolisvecJ_fixb_mdp_slag_vp: MDP has to be [32x1] (double)'); 

%% Symbolic Calculation
% From coriolisvec_joint_fixb_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 12:53:46
% EndTime: 2019-03-09 12:53:55
% DurationCPUTime: 6.34s
% Computational Cost: add. (5828->517), mult. (12720->657), div. (0->0), fcn. (7884->6), ass. (0->227)
t502 = cos(qJ(4));
t503 = cos(qJ(2));
t581 = qJD(1) * t503;
t555 = t502 * t581;
t500 = sin(qJ(4));
t570 = t500 * qJD(2);
t436 = -t555 - t570;
t556 = t500 * t581;
t577 = qJD(2) * t502;
t437 = -t556 + t577;
t499 = sin(qJ(5));
t618 = cos(qJ(5));
t384 = -t618 * t436 + t437 * t499;
t523 = t499 * t436 + t437 * t618;
t648 = t384 * t523;
t501 = sin(qJ(2));
t602 = t500 * t501;
t529 = pkin(4) * t503 - pkin(9) * t602;
t582 = qJD(1) * t501;
t490 = pkin(2) * t582;
t534 = pkin(8) * t501 - qJ(3) * t503;
t416 = qJD(1) * t534 + t490;
t486 = pkin(7) * t581;
t445 = pkin(3) * t581 + t486;
t543 = -t416 * t500 + t502 * t445;
t574 = qJD(4) * t500;
t504 = -pkin(2) - pkin(8);
t616 = pkin(9) - t504;
t647 = qJD(1) * t529 - t616 * t574 + t543;
t450 = t616 * t502;
t561 = t502 * t582;
t587 = t502 * t416 + t500 * t445;
t646 = pkin(9) * t561 + qJD(4) * t450 + t587;
t572 = qJD(4) * t503;
t557 = t502 * t572;
t560 = t501 * t570;
t513 = -t557 + t560;
t566 = qJD(2) * qJD(4);
t510 = -qJD(1) * t513 + t500 * t566;
t482 = t502 * t566;
t558 = t500 * t572;
t559 = t501 * t577;
t514 = t558 + t559;
t511 = qJD(1) * t514 - t482;
t553 = t618 * qJD(5);
t571 = qJD(5) * t499;
t351 = -t436 * t553 + t437 * t571 - t499 * t511 + t618 * t510;
t478 = qJD(4) + t582;
t463 = qJD(5) + t478;
t335 = t384 * t463 - t351;
t352 = qJD(5) * t523 - t499 * t510 - t618 * t511;
t567 = qJD(1) * qJD(2);
t551 = t503 * t567;
t621 = t523 ^ 2;
t645 = MDP(26) * t551 + (t463 * t523 - t352) * MDP(25) + MDP(22) * t648 + (-t384 ^ 2 + t621) * MDP(23) + t335 * MDP(24);
t619 = pkin(3) + pkin(7);
t496 = qJD(2) * qJ(3);
t425 = t496 + t445;
t392 = -pkin(4) * t436 + t425;
t345 = pkin(5) * t384 - qJ(6) * t523 + t392;
t644 = t345 * t384;
t643 = t384 * t392;
t573 = qJD(4) * t502;
t628 = t618 * qJD(4) + t553;
t390 = -t499 * t573 - t500 * t628 - t502 * t571;
t537 = t618 * t582;
t408 = -t499 * t561 - t500 * t537;
t589 = t390 + t408;
t603 = t499 * t500;
t588 = -t499 * t574 - t500 * t571 - t582 * t603 + (t537 + t628) * t502;
t485 = pkin(7) * t582;
t641 = qJD(3) + t485;
t452 = t463 * qJD(6);
t467 = qJ(6) * t551;
t614 = qJ(3) * t501;
t549 = -pkin(1) - t614;
t632 = t503 * t504;
t432 = t549 + t632;
t406 = t432 * qJD(1);
t569 = pkin(3) * t582 + t641;
t411 = qJD(2) * t504 + t569;
t373 = t502 * t406 + t500 * t411;
t365 = pkin(9) * t436 + t373;
t515 = t529 * qJD(2);
t552 = t501 * t567;
t477 = pkin(2) * t552;
t575 = qJD(3) * t501;
t512 = qJD(2) * t534 - t575;
t396 = qJD(1) * t512 + t477;
t476 = pkin(7) * t551;
t431 = pkin(3) * t551 + t476;
t545 = -t500 * t396 + t502 * t431;
t336 = qJD(1) * t515 - qJD(4) * t365 + t545;
t518 = t502 * t396 - t406 * t574 + t411 * t573 + t500 * t431;
t341 = pkin(9) * t511 + t518;
t372 = -t406 * t500 + t502 * t411;
t364 = -pkin(9) * t437 + t372;
t357 = pkin(4) * t478 + t364;
t539 = -t499 * t336 - t618 * t341 - t357 * t553 + t365 * t571;
t323 = t467 + t452 - t539;
t538 = t618 * t336 - t499 * t341 - t357 * t571 - t365 * t553;
t540 = pkin(5) * t551;
t324 = -t538 - t540;
t604 = t499 * t365;
t331 = t357 * t618 - t604;
t568 = qJD(6) - t331;
t329 = -t463 * pkin(5) + t568;
t563 = t618 * t365;
t332 = t499 * t357 + t563;
t330 = t463 * qJ(6) + t332;
t562 = t618 * t502;
t438 = -t562 + t603;
t439 = t499 * t502 + t500 * t618;
t640 = t323 * t439 + t324 * t438 - t329 * t589 + t330 * t588;
t358 = pkin(5) * t523 + qJ(6) * t384;
t638 = -0.2e1 * t567;
t497 = t501 ^ 2;
t498 = t503 ^ 2;
t584 = t497 - t498;
t636 = MDP(5) * t584;
t613 = t345 * t523;
t633 = t392 * t523;
t449 = t616 * t500;
t522 = t499 * t449 - t450 * t618;
t631 = -qJD(5) * t522 + t647 * t499 + t646 * t618;
t394 = -t449 * t618 - t499 * t450;
t630 = -qJD(5) * t394 + t646 * t499 - t647 * t618;
t458 = t619 * t501;
t441 = t502 * t458;
t548 = pkin(9) * t503 - t432;
t378 = pkin(4) * t501 + t500 * t548 + t441;
t440 = t500 * t458;
t586 = t502 * t432 + t440;
t599 = t502 * t503;
t381 = -pkin(9) * t599 + t586;
t629 = t499 * t378 + t618 * t381;
t564 = -pkin(4) * t502 - pkin(3);
t585 = pkin(4) * t573 - t564 * t582 + t641;
t627 = -MDP(28) + MDP(31);
t626 = t351 * t438 + t523 * t589;
t609 = t425 * t501;
t625 = qJD(2) * (t614 - t632) - t609;
t578 = qJD(2) * t501;
t489 = pkin(2) * t578;
t402 = t489 + t512;
t576 = qJD(2) * t503;
t446 = t619 * t576;
t544 = -t402 * t500 + t502 * t446;
t346 = t515 + (t502 * t548 - t440) * qJD(4) + t544;
t517 = t502 * t402 - t432 * t574 + t500 * t446 + t458 * t573;
t349 = pkin(9) * t514 + t517;
t624 = -qJD(5) * t629 + t346 * t618 - t499 * t349;
t623 = t503 * (qJD(4) + qJD(5));
t620 = t502 ^ 2;
t615 = qJD(2) * pkin(2);
t444 = t619 * t578;
t495 = qJD(2) * qJD(3);
t413 = -qJD(1) * t444 + t495;
t611 = t413 * t500;
t610 = t413 * t502;
t608 = t437 * t478;
t607 = t437 * t503;
t606 = t478 * t502;
t605 = t478 * t504;
t601 = t500 * t503;
t505 = qJD(2) ^ 2;
t600 = t501 * t505;
t598 = t503 * t505;
t506 = qJD(1) ^ 2;
t597 = t503 * t506;
t481 = t500 * pkin(4) + qJ(3);
t338 = t364 * t618 - t604;
t596 = -pkin(4) * t553 - qJD(6) + t338;
t595 = qJ(6) * t581 + t631;
t594 = -pkin(5) * t581 + t630;
t593 = pkin(5) * t588 - qJ(6) * t589 + qJD(6) * t438 + t585;
t590 = t463 * t390 - t438 * t551;
t459 = t619 * t503;
t451 = -pkin(2) * t503 + t549;
t426 = qJD(1) * t451;
t580 = qJD(2) * t522;
t579 = qJD(2) * t394;
t565 = t501 * t597;
t424 = pkin(4) * t599 + t459;
t550 = qJD(1) * t572;
t547 = pkin(1) * t638;
t546 = qJD(3) - t615;
t542 = -t436 - t570;
t541 = -qJD(4) + t582;
t460 = t501 * t551;
t337 = t499 * t364 + t563;
t536 = pkin(4) * t571 - t337;
t531 = t502 * (-qJD(1) * t459 - t425);
t530 = -0.2e1 * qJD(2) * t426;
t521 = -qJ(3) * t576 - t575;
t404 = qJD(1) * t521 + t477;
t421 = t489 + t521;
t528 = pkin(7) * t505 + qJD(1) * t421 + t404;
t526 = t378 * t618 - t499 * t381;
t520 = t331 * t463 + t539;
t519 = t332 * t463 + t538;
t516 = t499 * t346 + t618 * t349 + t378 * t553 - t381 * t571;
t395 = -pkin(4) * t558 + (-pkin(7) + t564) * t578;
t447 = pkin(7) * t552 - t495;
t448 = t485 + t546;
t457 = -t486 - t496;
t509 = -t447 * t503 + (t448 * t503 + (t457 + t486) * t501) * qJD(2);
t377 = t482 * pkin(4) + qJD(1) * t395 + t495;
t484 = -pkin(4) * t618 - pkin(5);
t480 = pkin(4) * t499 + qJ(6);
t461 = t502 * t551;
t442 = -qJ(3) * t581 + t490;
t418 = t439 * t503;
t417 = t499 * t601 - t503 * t562;
t410 = t426 * t582;
t382 = pkin(5) * t439 + qJ(6) * t438 + t481;
t368 = -pkin(5) * t417 + qJ(6) * t418 + t424;
t367 = t439 * t623 - t499 * t560 + t618 * t559;
t366 = t438 * t623 + t439 * t578;
t354 = pkin(4) * t437 + t358;
t350 = -t501 * pkin(5) - t526;
t348 = qJ(6) * t501 + t629;
t328 = -pkin(5) * t367 - qJ(6) * t366 + qJD(6) * t418 + t395;
t327 = t352 * pkin(5) + t351 * qJ(6) - qJD(6) * t523 + t377;
t326 = -pkin(5) * t576 - t624;
t325 = qJ(6) * t576 + qJD(6) * t501 + t516;
t1 = [(t544 * t478 + t444 * t436 + t459 * t482 + (qJD(2) * t531 + t545) * t501 + (-t373 * t501 - t478 * t586) * qJD(4) + (-t425 * t574 + t372 * qJD(2) + t610 + ((-t432 * t500 + t441) * qJD(2) - t459 * t574) * qJD(1)) * t503) * MDP(20) + 0.2e1 * MDP(4) * t460 + (t331 * t576 + t424 * t352 - t392 * t367 - t377 * t417 + t395 * t384 + t463 * t624 + t538 * t501 + t526 * t551) * MDP(27) + (t463 * t576 + t460) * MDP(26) + (t437 * t513 + t510 * t601) * MDP(15) + (-t501 * t528 + t503 * t530) * MDP(13) + (t501 * t530 + t503 * t528) * MDP(12) + ((t500 * t482 + (-t502 * t436 + (t437 + t577) * t500 + (-t500 ^ 2 + t620) * t581) * qJD(4)) * t503 + (t437 * t502 + (t436 - 0.2e1 * t555) * t500) * t578) * MDP(16) + MDP(6) * t598 + t636 * t638 + (t478 * t558 + (t500 * t550 - t482) * t501 + (t436 * t503 + (qJD(1) * t584 + t478 * t501) * t502) * qJD(2)) * MDP(18) + (-t516 * t463 + t539 * t501 + t395 * t523 - t424 * t351 - t377 * t418 + t392 * t366 + (-qJD(1) * t629 - t332) * t576) * MDP(28) + ((-t478 - t582) * t557 + (t607 + (-qJD(1) * t498 + (t478 + t541) * t501) * t500) * qJD(2)) * MDP(17) + t509 * MDP(11) + (pkin(7) * t509 + t404 * t451 + t421 * t426) * MDP(14) + (-t352 * t501 + t367 * t463 + (qJD(1) * t417 - t384) * t576) * MDP(25) + (-t324 * t501 - t326 * t463 - t327 * t417 + t328 * t384 - t345 * t367 + t352 * t368 + (-qJD(1) * t350 - t329) * t576) * MDP(29) + (t478 * t576 + t460) * MDP(19) + (-pkin(7) * t598 + t501 * t547) * MDP(9) - MDP(7) * t600 + (pkin(7) * t600 + t503 * t547) * MDP(10) + (-t351 * t501 + t366 * t463 + (-qJD(1) * t418 + t523) * t576) * MDP(24) + (t323 * t501 + t325 * t463 + t327 * t418 - t328 * t523 - t345 * t366 + t351 * t368 + (qJD(1) * t348 + t330) * t576) * MDP(31) + (-t351 * t417 + t352 * t418 - t366 * t384 + t367 * t523) * MDP(23) + (t323 * t417 - t324 * t418 - t325 * t384 + t326 * t523 + t329 * t366 + t330 * t367 - t348 * t352 - t350 * t351) * MDP(30) + (t351 * t418 + t366 * t523) * MDP(22) + (-t517 * t478 - t518 * t501 - t444 * t437 + (qJD(4) * t531 - t611) * t503 + ((-qJD(1) * t586 - t373) * t503 + (t459 * t541 + t609) * t500) * qJD(2)) * MDP(21) + (t323 * t348 + t324 * t350 + t325 * t330 + t326 * t329 + t327 * t368 + t328 * t345) * MDP(32); t410 * MDP(12) + (-t620 * t550 + (t541 * t577 - t608) * t500) * MDP(15) + (qJ(3) * t482 + t611 - t543 * t478 - t569 * t436 + (t425 * t502 - t500 * t605) * qJD(4) + ((-qJ(3) * t574 - t372) * t503 - t625 * t502) * qJD(1)) * MDP(20) + t626 * MDP(22) + (t481 * t352 + t377 * t439 + t585 * t384 + t588 * t392 + t630 * t463) * MDP(27) - t588 * t463 * MDP(25) + (t351 * t439 + t352 * t438 - t384 * t589 - t523 * t588) * MDP(23) + (t408 * t463 + t590) * MDP(24) + (-t481 * t351 - t377 * t438 + t589 * t392 + t631 * t463 + t585 * t523) * MDP(28) + (t327 * t439 + t588 * t345 + t352 * t382 + t593 * t384 + t594 * t463) * MDP(29) + (t351 * t522 - t352 * t394 + t384 * t595 - t523 * t594 - t640) * MDP(30) + (t327 * t438 - t589 * t345 + t351 * t382 - t595 * t463 - t523 * t593) * MDP(31) + (t323 * t394 - t324 * t522 + t327 * t382 - t329 * t594 - t330 * t595 + t345 * t593) * MDP(32) + (-t478 * t574 + t461 + (-t478 * t602 - t607) * qJD(1)) * MDP(17) + (t610 + t587 * t478 + t569 * t437 + (-t502 * t605 + (-t425 - t496) * t500) * qJD(4) + ((-qJ(3) * t573 + t373) * t503 + t625 * t500) * qJD(1)) * MDP(21) + (-t478 * t573 + (-t501 * t606 + t503 * t542) * qJD(1)) * MDP(18) - MDP(4) * t565 + (-qJ(3) * t447 - qJD(3) * t457 - t426 * t442 + (-t457 * t501 + (-t448 - t615) * t503) * qJD(1) * pkin(7)) * MDP(14) + (0.2e1 * t495 + (t426 * t503 + t442 * t501) * qJD(1)) * MDP(13) + ((-t437 * qJD(4) - t482 + (-t437 + t577) * t582) * t502 + ((-t436 + t570) * qJD(4) + (t501 * t542 + 0.2e1 * t557) * qJD(1)) * t500) * MDP(16) + t506 * t636 + ((-t457 - t496) * t501 + (-t448 + t546) * t503) * qJD(1) * MDP(11) + (-t442 * MDP(12) + (-t331 + t580) * MDP(27) + (-qJD(2) * t439 + t384) * MDP(25) - t523 * MDP(24) + (t332 - t579) * MDP(28) + (t329 + t580) * MDP(29) + (-t330 + t579) * MDP(31) - t478 * MDP(19) - t463 * MDP(26)) * t581 + (MDP(9) * t501 * t506 + MDP(10) * t597) * pkin(1); MDP(12) * t565 + (-t497 * t506 - t505) * MDP(13) + (t410 + t476) * MDP(14) + t461 * MDP(20) + t590 * MDP(27) + (-t352 * t439 - t384 * t588 - t626) * MDP(30) + t640 * MDP(32) + (-MDP(20) * t478 * t500 - MDP(21) * t606) * t478 + (t408 * MDP(27) + MDP(29) * t589 + t588 * t627) * t463 + (t457 * MDP(14) + t436 * MDP(20) + (-t437 - t556) * MDP(21) - t384 * MDP(27) + (-t438 * t581 - t384) * MDP(29) - t345 * MDP(32) + t627 * (t439 * t581 + t523)) * qJD(2); -t437 * t436 * MDP(15) + (-t436 ^ 2 + t437 ^ 2) * MDP(16) + (-t436 * t478 - t510) * MDP(17) + (t511 + t608) * MDP(18) + MDP(19) * t551 + (-t425 * t437 + t545 + (-qJD(4) + t478) * t373) * MDP(20) + (t372 * t478 - t425 * t436 - t518) * MDP(21) + (t337 * t463 - t633 + (-t384 * t437 - t463 * t571 + t551 * t618) * pkin(4) + t538) * MDP(27) + (t338 * t463 + t643 + (-t437 * t523 - t463 * t553 - t499 * t551) * pkin(4) + t539) * MDP(28) + (-t613 - t354 * t384 - t536 * t463 + (pkin(5) - t484) * t551 + t538) * MDP(29) + (-t351 * t484 - t352 * t480 + (t330 + t536) * t523 + (t329 + t596) * t384) * MDP(30) + (t354 * t523 - t463 * t596 + t480 * t551 + t323 - t644) * MDP(31) + (t323 * t480 + t324 * t484 + t329 * t536 - t330 * t596 - t345 * t354) * MDP(32) + t645; (t519 - t633) * MDP(27) + (t520 + t643) * MDP(28) + (-t358 * t384 + t519 + 0.2e1 * t540 - t613) * MDP(29) + (pkin(5) * t351 - qJ(6) * t352 + (t330 - t332) * t523 + (t329 - t568) * t384) * MDP(30) + (t358 * t523 + 0.2e1 * t452 + 0.2e1 * t467 - t520 - t644) * MDP(31) + (-pkin(5) * t324 + qJ(6) * t323 - t329 * t332 + t330 * t568 - t345 * t358) * MDP(32) + t645; (-t551 + t648) * MDP(29) + t335 * MDP(30) + (-t463 ^ 2 - t621) * MDP(31) + (-t330 * t463 + t324 + t613) * MDP(32);];
tauc  = t1;
