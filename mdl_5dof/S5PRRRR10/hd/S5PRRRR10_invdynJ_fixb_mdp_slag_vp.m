% Calculate vector of inverse dynamics joint torques for
% S5PRRRR10
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
% pkin [11x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,alpha2,alpha3,d2,d3,d4,d5,theta1]';
% MDP [25x1]
%   Minimal dynamic parameter vector (fixed base model)
%   see S5PRRRR10_convert_par2_MPV_fixb.m
% 
% Output:
% tau [5x1]
%   joint torques of inverse dynamics (contains inertial, gravitational Coriolis and centrifugal forces)

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-05 17:27
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tau = S5PRRRR10_invdynJ_fixb_mdp_slag_vp(qJ, qJD, qJDD, g, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(5,1),zeros(3,1),zeros(11,1),zeros(25,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PRRRR10_invdynJ_fixb_mdp_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5PRRRR10_invdynJ_fixb_mdp_slag_vp: qJD has to be [5x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [5 1]), ...
  'S5PRRRR10_invdynJ_fixb_mdp_slag_vp: qJDD has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5PRRRR10_invdynJ_fixb_mdp_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S5PRRRR10_invdynJ_fixb_mdp_slag_vp: pkin has to be [11x1] (double)');
assert(isreal(MDP) && all(size(MDP) == [25 1]), ...
  'S5PRRRR10_invdynJ_fixb_mdp_slag_vp: MDP has to be [25x1] (double)'); 

%% Symbolic Calculation
% From invdyn_joint_fixb_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 17:26:24
% EndTime: 2019-12-05 17:26:37
% DurationCPUTime: 8.36s
% Computational Cost: add. (3775->507), mult. (9769->741), div. (0->0), fcn. (8480->14), ass. (0->231)
t491 = sin(pkin(5));
t502 = cos(qJ(2));
t593 = qJD(1) * t502;
t463 = qJD(2) * pkin(2) + t491 * t593;
t493 = cos(pkin(6));
t490 = sin(pkin(6));
t494 = cos(pkin(5));
t594 = qJD(1) * t494;
t572 = t490 * t594;
t644 = t463 * t493 + t572;
t497 = sin(qJ(3));
t501 = cos(qJ(3));
t527 = (pkin(3) * t497 - pkin(9) * t501) * t490;
t498 = sin(qJ(2));
t615 = t491 * t498;
t571 = qJD(1) * t615;
t643 = qJD(3) * t527 - t490 * t571;
t603 = t501 * t502;
t609 = t497 * t498;
t524 = -t493 * t609 + t603;
t426 = t524 * t491;
t618 = t490 * t497;
t483 = pkin(8) * t618;
t604 = t501 * t493;
t638 = pkin(2) * t604 - t483;
t642 = qJD(1) * t426 - t638 * qJD(3);
t592 = qJD(2) * t490;
t569 = t501 * t592;
t475 = -qJD(4) + t569;
t500 = cos(qJ(4));
t496 = sin(qJ(4));
t570 = t497 * t592;
t543 = t496 * t570;
t591 = qJD(2) * t493;
t549 = qJD(3) + t591;
t422 = -t500 * t549 + t543;
t415 = qJD(5) + t422;
t641 = -pkin(8) * qJDD(2) * t490 - (qJD(2) * t593 + qJDD(1) * t498) * t491 - t644 * qJD(3);
t580 = qJDD(2) * t493;
t482 = qJDD(3) + t580;
t529 = qJD(4) * t549;
t589 = qJD(3) * t501;
t564 = t496 * t589;
t579 = qJDD(2) * t497;
t585 = qJD(4) * t500;
t367 = -t500 * t482 + t490 * (qJD(2) * (t497 * t585 + t564) + t496 * t579) + t496 * t529;
t610 = t497 * t493;
t617 = t490 * t501;
t596 = pkin(2) * t610 + pkin(8) * t617;
t433 = pkin(9) * t493 + t596;
t538 = -pkin(3) * t501 - pkin(9) * t497;
t434 = (-pkin(2) + t538) * t490;
t598 = t500 * t433 + t496 * t434;
t640 = qJD(4) * t598 - t642 * t496 - t643 * t500;
t587 = qJD(4) * t496;
t639 = -t433 * t587 + t434 * t585 + t643 * t496 - t642 * t500;
t608 = t497 * t502;
t526 = t498 * t604 + t608;
t425 = t526 * t491;
t597 = -qJD(1) * t425 + t596 * qJD(3);
t637 = t493 * t603 - t609;
t578 = qJDD(2) * t501;
t480 = t490 * t578;
t582 = qJD(2) * qJD(3);
t562 = t497 * t582;
t435 = t490 * t562 + qJDD(4) - t480;
t614 = t491 * t502;
t481 = qJDD(1) * t614;
t568 = qJD(2) * t615;
t541 = qJD(1) * t568;
t430 = qJDD(2) * pkin(2) + t481 - t541;
t454 = pkin(8) * t592 + t571;
t581 = qJDD(1) * t494;
t560 = t490 * t581;
t590 = qJD(3) * t497;
t510 = -t430 * t610 + t454 * t590 - t497 * t560 + t641 * t501;
t334 = pkin(9) * t482 - t510;
t382 = t501 * t454 + t463 * t610 + t497 * t572;
t369 = pkin(9) * t549 + t382;
t479 = t493 * t594;
t390 = t479 + (qJD(2) * t538 - t463) * t490;
t341 = t369 * t500 + t390 * t496;
t478 = t493 * t581;
t521 = t562 - t578;
t561 = t501 * t582;
t522 = t561 + t579;
t358 = t478 + (pkin(3) * t521 - pkin(9) * t522 - t430) * t490;
t506 = -qJD(4) * t341 - t334 * t496 + t500 * t358;
t326 = -pkin(4) * t435 - t506;
t424 = t496 * t549 + t500 * t570;
t632 = sin(pkin(11));
t556 = t632 * t498;
t492 = cos(pkin(11));
t612 = t492 * t502;
t447 = t494 * t612 - t556;
t555 = t632 * t502;
t613 = t492 * t498;
t448 = t494 * t613 + t555;
t616 = t491 * t492;
t576 = t490 * t616;
t371 = t448 * t501 + (t447 * t493 - t576) * t497;
t449 = -t494 * t555 - t613;
t450 = -t494 * t556 + t612;
t557 = t491 * t632;
t539 = t490 * t557;
t373 = t450 * t501 + (t449 * t493 + t539) * t497;
t525 = t493 * t608 + t498 * t501;
t397 = t491 * t525 + t494 * t618;
t446 = -t490 * t614 + t493 * t494;
t374 = t397 * t496 - t446 * t500;
t398 = -t447 * t490 - t493 * t616;
t399 = -t449 * t490 + t493 * t557;
t518 = g(1) * (-t373 * t496 + t399 * t500) + g(2) * (-t371 * t496 + t398 * t500) - g(3) * t374;
t636 = t415 * (pkin(4) * t424 + t415 * pkin(10)) + t326 + t518;
t363 = qJDD(5) + t367;
t474 = -pkin(4) * t500 - pkin(10) * t496 - pkin(3);
t635 = t415 * (t382 + t475 * (pkin(4) * t496 - pkin(10) * t500)) - t474 * t363;
t381 = -t497 * t454 + t501 * t644;
t503 = qJD(2) ^ 2;
t633 = pkin(9) * qJD(4);
t565 = t490 * t589;
t542 = t500 * t565;
t559 = t490 * t579;
t366 = qJD(2) * t542 - qJD(4) * t543 + t496 * t482 + (t529 + t559) * t500;
t495 = sin(qJ(5));
t499 = cos(qJ(5));
t583 = qJD(5) * t499;
t573 = t499 * t366 + t495 * t435 - t475 * t583;
t584 = qJD(5) * t495;
t336 = -t424 * t584 + t573;
t631 = t336 * t495;
t626 = t424 * t495;
t391 = t499 * t475 + t626;
t630 = t391 * t415;
t393 = t424 * t499 - t475 * t495;
t629 = t393 * t415;
t628 = t422 * t475;
t627 = t424 * t475;
t625 = t446 * t490;
t622 = t475 * t496;
t621 = t482 * MDP(9);
t619 = t490 * t496;
t611 = t495 * t363;
t607 = t498 * t503;
t606 = t499 * t363;
t605 = t500 * t501;
t602 = qJDD(1) - g(3);
t566 = t490 * t590;
t601 = -pkin(4) * t566 + t640;
t439 = qJD(2) * t527;
t599 = t500 * t381 + t496 * t439;
t488 = t497 ^ 2;
t595 = -t501 ^ 2 + t488;
t588 = qJD(4) * t495;
t586 = qJD(4) * t499;
t575 = t495 * t617;
t563 = t490 * t493 * t503;
t520 = t500 * t334 + t496 * t358 - t369 * t587 + t390 * t585;
t325 = pkin(10) * t435 + t520;
t528 = t430 * t604 - t454 * t589 + t641 * t497 + t501 * t560;
t335 = -pkin(3) * t482 - t528;
t328 = pkin(4) * t367 - pkin(10) * t366 + t335;
t553 = -t495 * t325 + t499 * t328;
t551 = t366 * t495 - t499 * t435;
t550 = t415 * t499;
t548 = qJD(3) + 0.2e1 * t591;
t547 = t482 + t580;
t487 = t490 ^ 2;
t546 = t487 * t491 * t607;
t544 = t490 * t568;
t536 = g(1) * t450 + g(2) * t448;
t412 = (t495 * t497 + t499 * t605) * t592;
t535 = t499 * t585 - t412;
t432 = t483 + (-pkin(2) * t501 - pkin(3)) * t493;
t451 = -t500 * t493 + t496 * t618;
t452 = t493 * t496 + t500 * t618;
t376 = pkin(4) * t451 - pkin(10) * t452 + t432;
t534 = -pkin(10) * t566 - qJD(5) * t376 - t639;
t378 = -pkin(10) * t617 + t598;
t400 = -qJD(4) * t451 + t542;
t401 = qJD(4) * t452 + t490 * t564;
t533 = -pkin(4) * t401 + pkin(10) * t400 + qJD(5) * t378 - t597;
t532 = t499 * t325 + t495 * t328;
t339 = -pkin(10) * t475 + t341;
t368 = -pkin(3) * t549 - t381;
t342 = t422 * pkin(4) - t424 * pkin(10) + t368;
t330 = t339 * t499 + t342 * t495;
t531 = t339 * t495 - t342 * t499;
t340 = -t369 * t496 + t390 * t500;
t375 = t397 * t500 + t446 * t496;
t396 = -t491 * t637 - t494 * t617;
t346 = t375 * t499 + t396 * t495;
t345 = -t375 * t495 + t396 * t499;
t530 = -t433 * t496 + t434 * t500;
t402 = t452 * t495 + t499 * t617;
t370 = -t447 * t604 + t448 * t497 + t501 * t576;
t372 = -t449 * t604 + t450 * t497 - t501 * t539;
t517 = g(1) * t372 + g(2) * t370 + g(3) * t396;
t516 = -g(1) * t373 - g(2) * t371 - g(3) * t397;
t385 = t447 * t501 - t448 * t610;
t387 = t449 * t501 - t450 * t610;
t515 = g(1) * t387 + g(2) * t385 + g(3) * t426;
t512 = -t335 + t517;
t338 = pkin(4) * t475 - t340;
t509 = -pkin(10) * t363 + (t338 + t340) * t415;
t508 = -pkin(9) * t435 - t368 * t475;
t507 = pkin(9) * qJD(5) * t415 - t517;
t504 = (pkin(10) * t570 - qJD(5) * t474 + t599) * t415 + t516;
t414 = -t463 * t490 + t479;
t411 = t495 * t500 * t569 - t499 * t570;
t403 = t452 * t499 - t575;
t395 = -t430 * t490 + t478;
t394 = t426 * t500 + t615 * t619;
t386 = t449 * t497 + t450 * t604;
t384 = t447 * t497 + t448 * t604;
t377 = pkin(4) * t617 - t530;
t365 = t494 * t565 + (t524 * qJD(2) + qJD(3) * t637) * t491;
t364 = t494 * t566 + (qJD(2) * t526 + qJD(3) * t525) * t491;
t360 = t387 * t500 + t450 * t619;
t359 = t385 * t500 + t448 * t619;
t355 = -qJD(5) * t575 + t400 * t495 + t452 * t583 - t499 * t566;
t354 = -qJD(5) * t402 + t400 * t499 + t495 * t566;
t351 = -pkin(4) * t570 + t381 * t496 - t439 * t500;
t350 = t373 * t500 + t399 * t496;
t348 = t371 * t500 + t398 * t496;
t337 = qJD(5) * t393 + t551;
t332 = -qJD(4) * t374 + t365 * t500 + t496 * t544;
t331 = qJD(4) * t375 + t365 * t496 - t500 * t544;
t324 = -qJD(5) * t330 + t553;
t323 = -qJD(5) * t531 + t532;
t1 = [t602 * MDP(1) + (-t364 * t549 - t396 * t482 - t501 * t546 + t521 * t625) * MDP(10) + (-t365 * t549 - t397 * t482 + t497 * t546 + t522 * t625) * MDP(11) + (t331 * t475 + t364 * t422 + t367 * t396 - t374 * t435) * MDP(17) + (t332 * t475 + t364 * t424 + t366 * t396 - t375 * t435) * MDP(18) + ((-qJD(5) * t346 - t332 * t495 + t364 * t499) * t415 + t345 * t363 + t331 * t391 + t374 * t337) * MDP(24) + (-(qJD(5) * t345 + t332 * t499 + t364 * t495) * t415 - t346 * t363 + t331 * t393 + t374 * t336) * MDP(25) + ((qJDD(2) * t502 - t607) * MDP(3) + (-qJDD(2) * t498 - t502 * t503) * MDP(4)) * t491; qJDD(2) * MDP(2) + (-g(1) * t449 - g(2) * t447 - g(3) * t614 + t481) * MDP(3) + (-t602 * t615 + t536) * MDP(4) + (t638 * t482 + t528 * t493 - t597 * t549 - t515) * MDP(10) + (g(1) * t386 + g(2) * t384 + g(3) * t425 - t596 * t482 + t510 * t493 + t642 * t549) * MDP(11) + (t366 * t452 + t400 * t424) * MDP(12) + (-t366 * t451 - t367 * t452 - t400 * t422 - t401 * t424) * MDP(13) + (-t400 * t475 + t435 * t452) * MDP(14) + (t401 * t475 - t435 * t451) * MDP(15) + (-g(1) * t360 - g(2) * t359 - g(3) * t394 + t335 * t451 + t432 * t367 + t368 * t401 + t597 * t422 + t530 * t435 + t640 * t475) * MDP(17) + (t335 * t452 + t432 * t366 + t368 * t400 + t597 * t424 - t598 * t435 + t639 * t475 + t515 * t496) * MDP(18) + (t336 * t403 + t354 * t393) * MDP(19) + (-t336 * t402 - t337 * t403 - t354 * t391 - t355 * t393) * MDP(20) + (t336 * t451 + t354 * t415 + t363 * t403 + t393 * t401) * MDP(21) + (-t337 * t451 - t355 * t415 - t363 * t402 - t391 * t401) * MDP(22) + (t363 * t451 + t401 * t415) * MDP(23) + ((t376 * t499 - t378 * t495) * t363 + t324 * t451 - t531 * t401 + t377 * t337 + t326 * t402 + t338 * t355 - g(1) * (t360 * t499 + t386 * t495) - g(2) * (t359 * t499 + t384 * t495) - g(3) * (t394 * t499 + t425 * t495) + (t495 * t534 - t499 * t533) * t415 + t601 * t391) * MDP(24) + (-(t376 * t495 + t378 * t499) * t363 - t323 * t451 - t330 * t401 + t377 * t336 + t326 * t403 + t338 * t354 - g(1) * (-t360 * t495 + t386 * t499) - g(2) * (-t359 * t495 + t384 * t499) - g(3) * (-t394 * t495 + t425 * t499) + (t495 * t533 + t499 * t534) * t415 + t601 * t393) * MDP(25) + t493 * t621 + ((qJDD(2) * t488 + 0.2e1 * t497 * t561) * MDP(5) + 0.2e1 * (t497 * t578 - t582 * t595) * MDP(6) + (-pkin(2) * t521 + t501 * t541) * MDP(10) + (-pkin(2) * t522 - t497 * t541) * MDP(11)) * t487 + ((t497 * t547 + t548 * t589) * MDP(7) + (t501 * t547 - t548 * t590) * MDP(8) + (-t395 * t501 + t414 * t590) * MDP(10) + (t395 * t497 + t414 * t589) * MDP(11) + (-t366 * t501 + t424 * t590) * MDP(14) + (t367 * t501 - t422 * t590) * MDP(15) + (-t435 * t501 - t475 * t590) * MDP(16) + (t340 * t590 - t501 * t506) * MDP(17) + (t520 * t501 - t341 * t590 + (-g(3) * t615 - t536) * t500) * MDP(18)) * t490; (-t501 * t563 + t559) * MDP(7) + (t497 * t563 + t480) * MDP(8) + t621 + (t382 * t549 - t414 * t570 + t517 + t528) * MDP(10) + (t381 * t549 - t414 * t569 + t510 - t516) * MDP(11) + (t366 * t496 - t500 * t627) * MDP(12) + ((t366 + t628) * t500 + (-t367 + t627) * t496) * MDP(13) + (-t475 * t585 + t435 * t496 + (-t424 * t497 + t475 * t605) * t592) * MDP(14) + (t475 * t587 + t435 * t500 + (t422 * t497 - t501 * t622) * t592) * MDP(15) + t475 * MDP(16) * t570 + (-t340 * t570 - pkin(3) * t367 - t382 * t422 + (-t381 * t475 + t508) * t496 + ((t439 + t633) * t475 + t512) * t500) * MDP(17) + (-pkin(3) * t366 - t599 * t475 + t341 * t570 - t382 * t424 + t508 * t500 + (-t475 * t633 - t512) * t496) * MDP(18) + (t336 * t496 * t499 + (-t496 * t584 + t535) * t393) * MDP(19) + (t391 * t412 + t393 * t411 + (-t391 * t499 - t393 * t495) * t585 + (-t631 - t337 * t499 + (t391 * t495 - t393 * t499) * qJD(5)) * t496) * MDP(20) + (-t336 * t500 + t535 * t415 + (-t393 * t475 - t415 * t584 + t606) * t496) * MDP(21) + (t337 * t500 + (-t495 * t585 + t411) * t415 + (t391 * t475 - t415 * t583 - t611) * t496) * MDP(22) + (-t363 * t500 - t415 * t622) * MDP(23) + (-t338 * t411 - t351 * t391 - t635 * t499 + t504 * t495 + (t338 * t588 - t324 + (qJD(4) * t391 - t611) * pkin(9) - t507 * t499) * t500 + (t338 * t583 + t326 * t495 + t475 * t531 + (t415 * t588 + t337) * pkin(9)) * t496) * MDP(24) + (-t338 * t412 - t351 * t393 + t635 * t495 + t504 * t499 + (t338 * t586 + t323 + (qJD(4) * t393 - t606) * pkin(9) + t507 * t495) * t500 + (-t338 * t584 + t326 * t499 + t475 * t330 + (t415 * t586 + t336) * pkin(9)) * t496) * MDP(25) + (-MDP(5) * t497 * t501 + MDP(6) * t595) * t487 * t503; -t422 ^ 2 * MDP(13) + (t366 - t628) * MDP(14) + (-t367 - t627) * MDP(15) + t435 * MDP(16) + (-t341 * t475 + t506 - t518) * MDP(17) + (g(1) * t350 + g(2) * t348 + g(3) * t375 - t340 * t475 + t368 * t422 - t520) * MDP(18) + (t393 * t550 + t631) * MDP(19) + ((t336 - t630) * t499 + (-t337 - t629) * t495) * MDP(20) + (t415 * t550 + t611) * MDP(21) + (-t415 ^ 2 * t495 + t606) * MDP(22) + (-pkin(4) * t337 - t341 * t391 + t509 * t495 - t499 * t636) * MDP(24) + (-pkin(4) * t336 - t341 * t393 + t495 * t636 + t509 * t499) * MDP(25) + (MDP(12) * t422 + t424 * MDP(13) - t368 * MDP(17) - t393 * MDP(21) + t391 * MDP(22) - t415 * MDP(23) + MDP(24) * t531 + t330 * MDP(25)) * t424; t393 * t391 * MDP(19) + (-t391 ^ 2 + t393 ^ 2) * MDP(20) + (t573 + t630) * MDP(21) + (-t551 + t629) * MDP(22) + t363 * MDP(23) + (t330 * t415 - t338 * t393 - g(1) * (-t350 * t495 + t372 * t499) - g(2) * (-t348 * t495 + t370 * t499) - g(3) * t345 + t553) * MDP(24) + (-t531 * t415 + t338 * t391 - g(1) * (-t350 * t499 - t372 * t495) - g(2) * (-t348 * t499 - t370 * t495) + g(3) * t346 - t532) * MDP(25) + (-MDP(21) * t626 - MDP(22) * t393 - MDP(24) * t330 + MDP(25) * t531) * qJD(5);];
tau = t1;
