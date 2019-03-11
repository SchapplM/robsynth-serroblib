% Calculate vector of inverse dynamics joint torques for
% S6PRPRRR3
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
%   pkin=[a2,a3,a4,a5,a6,alpha2,d2,d4,d5,d6,theta1,theta3]';
% MDP [29x1]
%   Minimal dynamic parameter vector (fixed base model)
%   see S6PRPRRR3_convert_par2_MPV_fixb.m
% 
% Output:
% tau [6x1]
%   joint torques of inverse dynamics (contains inertial, gravitational Coriolis and centrifugal forces)

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-08 20:34
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tau = S6PRPRRR3_invdynJ_fixb_mdp_slag_vp(qJ, qJD, qJDD, g, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(3,1),zeros(12,1),zeros(29,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRPRRR3_invdynJ_fixb_mdp_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6PRPRRR3_invdynJ_fixb_mdp_slag_vp: qJD has to be [6x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [6 1]), ...
  'S6PRPRRR3_invdynJ_fixb_mdp_slag_vp: qJDD has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6PRPRRR3_invdynJ_fixb_mdp_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [12 1]), ...
  'S6PRPRRR3_invdynJ_fixb_mdp_slag_vp: pkin has to be [12x1] (double)');
assert(isreal(MDP) && all(size(MDP) == [29 1]), ...
  'S6PRPRRR3_invdynJ_fixb_mdp_slag_vp: MDP has to be [29x1] (double)'); 

%% Symbolic Calculation
% From invdyn_joint_fixb_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-08 20:34:27
% EndTime: 2019-03-08 20:34:36
% DurationCPUTime: 7.62s
% Computational Cost: add. (5200->476), mult. (12402->630), div. (0->0), fcn. (10655->18), ass. (0->218)
t507 = cos(qJ(6));
t581 = qJD(6) * t507;
t498 = sin(pkin(12));
t505 = sin(qJ(4));
t586 = qJD(2) * t505;
t565 = t498 * t586;
t501 = cos(pkin(12));
t509 = cos(qJ(4));
t585 = qJD(2) * t509;
t567 = t501 * t585;
t450 = -t565 + t567;
t460 = t498 * t509 + t501 * t505;
t451 = t460 * qJD(2);
t504 = sin(qJ(5));
t508 = cos(qJ(5));
t405 = -t508 * t450 + t451 * t504;
t651 = t405 * t507;
t657 = t581 + t651;
t500 = sin(pkin(6));
t510 = cos(qJ(2));
t602 = t500 * t510;
t523 = t460 * t602;
t625 = pkin(8) + qJ(3);
t464 = t625 * t498;
t465 = t625 * t501;
t593 = -t505 * t464 + t509 * t465;
t656 = qJD(1) * t523 - t460 * qJD(3) - qJD(4) * t593;
t457 = t509 * t464;
t600 = t509 * t501;
t459 = t498 * t505 - t600;
t522 = t459 * t602;
t655 = -qJD(1) * t522 + t505 * (qJD(3) * t498 + qJD(4) * t465) - qJD(3) * t600 + qJD(4) * t457;
t535 = t450 * t504 + t508 * t451;
t574 = qJDD(2) * t509;
t575 = qJDD(2) * t505;
t570 = qJD(4) * t567 + t498 * t574 + t501 * t575;
t419 = -qJD(4) * t565 + t570;
t453 = t460 * qJD(4);
t477 = t501 * t574;
t544 = -t498 * t575 + t477;
t420 = qJD(2) * t453 - t544;
t583 = qJD(5) * t508;
t584 = qJD(5) * t504;
t365 = t508 * t419 - t504 * t420 + t450 * t583 - t451 * t584;
t493 = qJDD(4) + qJDD(5);
t497 = qJD(4) + qJD(5);
t503 = sin(qJ(6));
t571 = t507 * t365 + t503 * t493 + t497 * t581;
t582 = qJD(6) * t503;
t353 = -t535 * t582 + t571;
t352 = t353 * t507;
t392 = t497 * t503 + t507 * t535;
t558 = t365 * t503 - t507 * t493;
t354 = qJD(6) * t392 + t558;
t609 = t535 * t503;
t390 = -t507 * t497 + t609;
t654 = -t503 * t354 - t657 * t390 + t352;
t351 = t353 * t503;
t366 = qJD(5) * t535 + t419 * t504 + t508 * t420;
t363 = qJDD(6) + t366;
t360 = t503 * t363;
t645 = -qJD(6) - t405;
t595 = -t581 * t645 + t360;
t611 = t405 * t497;
t613 = t535 * t497;
t615 = t392 * t535;
t653 = t493 * MDP(20) + (-t366 + t613) * MDP(19) - t405 ^ 2 * MDP(17) + (MDP(16) * t405 + MDP(17) * t535 + MDP(27) * t645) * t535 + (t365 + t611) * MDP(18) + (t657 * t392 + t351) * MDP(23) + (-t645 * t651 + t595 - t615) * MDP(25);
t506 = sin(qJ(2));
t588 = qJD(1) * t506;
t569 = t500 * t588;
t463 = qJD(2) * qJ(3) + t569;
t502 = cos(pkin(6));
t589 = qJD(1) * t502;
t474 = t501 * t589;
t424 = t474 + (-pkin(8) * qJD(2) - t463) * t498;
t433 = t501 * t463 + t498 * t589;
t587 = qJD(2) * t501;
t425 = pkin(8) * t587 + t433;
t639 = t509 * t424 - t425 * t505;
t375 = -pkin(9) * t451 + t639;
t374 = qJD(4) * pkin(4) + t375;
t537 = -t424 * t505 - t425 * t509;
t376 = pkin(9) * t450 - t537;
t618 = t376 * t504;
t355 = t374 * t508 - t618;
t348 = -pkin(5) * t497 - t355;
t652 = t348 * t405;
t554 = t645 * t503;
t568 = qJD(1) * t602;
t576 = qJDD(2) * qJ(3);
t578 = qJDD(1) * t500;
t431 = t506 * t578 + t576 + (qJD(3) + t568) * qJD(2);
t577 = qJDD(1) * t502;
t471 = t501 * t577;
t622 = pkin(8) * qJDD(2);
t395 = t471 + (-t431 - t622) * t498;
t413 = t501 * t431 + t498 * t577;
t396 = t501 * t622 + t413;
t557 = t509 * t395 - t505 * t396;
t344 = qJDD(4) * pkin(4) - pkin(9) * t419 + qJD(4) * t537 + t557;
t540 = t505 * t395 + t509 * t396;
t345 = -pkin(9) * t420 + qJD(4) * t639 + t540;
t617 = t376 * t508;
t356 = t374 * t504 + t617;
t630 = -qJD(5) * t356 + t508 * t344 - t345 * t504;
t334 = -pkin(5) * t493 - t630;
t624 = cos(pkin(11));
t562 = t624 * t506;
t499 = sin(pkin(11));
t604 = t499 * t510;
t447 = t502 * t562 + t604;
t561 = t624 * t510;
t605 = t499 * t506;
t449 = -t502 * t605 + t561;
t496 = pkin(12) + qJ(4);
t492 = qJ(5) + t496;
t484 = sin(t492);
t485 = cos(t492);
t563 = t500 * t624;
t603 = t500 * t506;
t606 = t499 * t500;
t525 = -g(3) * (-t484 * t603 + t485 * t502) - g(2) * (-t447 * t484 - t485 * t563) - g(1) * (-t449 * t484 + t485 * t606);
t520 = -t334 + t525;
t650 = -pkin(9) * t453 - t655;
t452 = t459 * qJD(4);
t649 = -pkin(9) * t452 - t656;
t601 = t501 * MDP(5);
t647 = -MDP(6) * t498 + t601;
t421 = t508 * t459 + t460 * t504;
t422 = -t459 * t504 + t460 * t508;
t486 = -pkin(3) * t501 - pkin(2);
t437 = pkin(4) * t459 + t486;
t370 = pkin(5) * t421 - pkin(10) * t422 + t437;
t446 = -t502 * t561 + t605;
t448 = t502 * t604 + t562;
t548 = g(1) * t448 + g(2) * t446;
t521 = g(3) * t602 - t548;
t518 = t521 * t485;
t646 = t370 * t363 - t518;
t378 = pkin(5) * t535 + pkin(10) * t405;
t415 = t447 * t485 - t484 * t563;
t417 = t449 * t485 + t484 * t606;
t545 = qJD(3) - t568;
t442 = qJD(2) * t486 + t545;
t418 = -pkin(4) * t450 + t442;
t435 = t484 * t502 + t485 * t603;
t631 = (qJD(5) * t374 + t345) * t508 + t344 * t504 - t376 * t584;
t644 = g(1) * t417 + g(2) * t415 + g(3) * t435 + t405 * t418 - t631;
t592 = t498 ^ 2 + t501 ^ 2;
t640 = MDP(7) * t592;
t616 = t390 * t535;
t361 = t507 * t363;
t527 = -t582 * t645 - t361;
t531 = pkin(4) * t453 - t569;
t349 = pkin(10) * t497 + t356;
t364 = pkin(5) * t405 - pkin(10) * t535 + t418;
t541 = t349 * t503 - t364 * t507;
t638 = t348 * t582 + t535 * t541;
t336 = t349 * t507 + t364 * t503;
t637 = t336 * t535 + t348 * t581 - t503 * t520;
t636 = -t535 * t418 + t525 + t630;
t333 = pkin(10) * t493 + t631;
t555 = -t465 * t505 - t457;
t397 = -pkin(9) * t460 + t555;
t398 = -pkin(9) * t459 + t593;
t369 = t397 * t504 + t398 * t508;
t379 = -qJD(5) * t421 - t452 * t508 - t453 * t504;
t547 = g(1) * t449 + g(2) * t447;
t573 = g(3) * t603;
t539 = t397 * t508 - t398 * t504;
t598 = -qJD(5) * t539 + t504 * t649 - t508 * t650;
t635 = -(qJD(6) * t364 + t333) * t421 + t334 * t422 + t348 * t379 - (-qJD(6) * t370 + t598) * t645 - t369 * t363 - t573 - t547;
t621 = qJDD(2) * pkin(2);
t620 = t348 * t422;
t614 = t392 * t503;
t607 = (-t463 * t498 + t474) * t498;
t599 = qJDD(1) - g(3);
t597 = qJD(5) * t369 + t504 * t650 + t508 * t649;
t564 = qJD(2) * t588;
t488 = pkin(4) * t504 + pkin(10);
t551 = pkin(4) * t451 + qJD(6) * t488 + t378;
t357 = t375 * t504 + t617;
t550 = pkin(4) * t584 - t357;
t358 = t375 * t508 - t618;
t549 = -pkin(4) * t583 + t358;
t380 = qJD(5) * t422 - t452 * t504 + t508 * t453;
t546 = pkin(5) * t380 - pkin(10) * t379 + t531;
t542 = -t488 * t363 + t652;
t443 = -t498 * t603 + t501 * t502;
t444 = t498 * t502 + t501 * t603;
t402 = t443 * t509 - t444 * t505;
t403 = t443 * t505 + t444 * t509;
t538 = t402 * t508 - t403 * t504;
t373 = t402 * t504 + t403 * t508;
t536 = -t433 * t501 + t607;
t528 = t500 * t564 - t510 * t578 + qJDD(3);
t436 = t528 - t621;
t533 = -t436 + t548;
t532 = t405 * t554 - t527;
t529 = MDP(3) + t647;
t526 = t379 * t507 - t422 * t582;
t524 = -qJD(2) * t607 + t433 * t587;
t412 = -t431 * t498 + t471;
t517 = -t412 * t498 + t413 * t501 - t547;
t426 = qJDD(2) * t486 + t528;
t383 = pkin(4) * t420 + t426;
t511 = qJD(2) ^ 2;
t491 = cos(t496);
t490 = sin(t496);
t489 = -pkin(4) * t508 - pkin(5);
t461 = -qJD(2) * pkin(2) + t545;
t382 = -qJD(2) * t523 - qJD(4) * t403;
t381 = -qJD(2) * t522 + qJD(4) * t402;
t340 = qJD(5) * t373 + t381 * t504 - t382 * t508;
t339 = qJD(5) * t538 + t381 * t508 + t382 * t504;
t338 = pkin(5) * t366 - pkin(10) * t365 + t383;
t337 = t507 * t338;
t1 = [t599 * MDP(1) + (t412 * t443 + t413 * t444 - g(3)) * MDP(8) + (qJD(4) * t382 + qJDD(4) * t402) * MDP(14) + (-qJD(4) * t381 - qJDD(4) * t403) * MDP(15) + (-t340 * t497 + t493 * t538) * MDP(21) + (-t339 * t497 - t373 * t493) * MDP(22) + (-(-t339 * t503 - t373 * t581) * t645 - t373 * t360 + t340 * t390 - t538 * t354) * MDP(28) + ((t339 * t507 - t373 * t582) * t645 - t373 * t361 + t340 * t392 - t538 * t353) * MDP(29) + (-t443 * t498 + t444 * t501) * MDP(7) * qJDD(2) + ((-qJDD(2) * MDP(4) - t529 * t511 + (-MDP(14) * t450 + MDP(15) * t451 + MDP(21) * t405 + MDP(22) * t535 + MDP(8) * t461 - (MDP(28) * t507 - MDP(29) * t503) * t645) * qJD(2)) * t506 + ((-t436 + t524) * MDP(8) - t420 * MDP(14) - t419 * MDP(15) - t366 * MDP(21) - t365 * MDP(22) + t527 * MDP(28) + t595 * MDP(29) + (-MDP(4) + t640) * t511 + t529 * qJDD(2)) * t510) * t500; (qJD(4) * t655 - t593 * qJDD(4) + t486 * t419 + t426 * t460 - t442 * t452 - t451 * t569 + t521 * t490) * MDP(15) + (qJD(4) * t656 + t555 * qJDD(4) + t486 * t420 + t426 * t459 + t442 * t453 + t450 * t569 - t521 * t491) * MDP(14) + t647 * (t500 * (-g(3) * t510 + t564) + t533 + t621) + (-t573 + t517 + (qJD(2) * t545 + t576) * t592) * MDP(7) + (-t336 * t380 - t539 * t353 + t597 * t392 + (-(-qJD(6) * t349 + t338) * t421 - qJD(6) * t620 - (qJD(6) * t369 - t546) * t645 - t646) * t503 + t635 * t507) * MDP(29) + (-t541 * t380 + t337 * t421 - t539 * t354 + t597 * t390 + (-t546 * t645 + (-t349 * t421 + t369 * t645 + t620) * qJD(6) + t646) * t507 + t635 * t503) * MDP(28) + (t353 * t421 + t361 * t422 + t380 * t392 - t526 * t645) * MDP(25) + (-t422 * t360 - t354 * t421 - t380 * t390 - (-t379 * t503 - t422 * t581) * t645) * MDP(26) + (t363 * t421 - t380 * t645) * MDP(27) + qJDD(2) * MDP(2) + (-t536 * qJD(3) + t533 * pkin(2) + t517 * qJ(3) + (-g(3) * (pkin(2) * t510 + qJ(3) * t506) + (-t461 * t506 + t510 * t536) * qJD(1)) * t500) * MDP(8) + (t599 * t602 + t548) * MDP(3) + (-t599 * t603 + t547) * MDP(4) + (t365 * t437 - t369 * t493 + t379 * t418 + t383 * t422 + t484 * t521 + t497 * t598 + t531 * t535) * MDP(22) + (t365 * t422 + t379 * t535) * MDP(16) + (-t365 * t421 - t366 * t422 - t379 * t405 - t380 * t535) * MDP(17) + (-qJD(4) * t452 + qJDD(4) * t460) * MDP(11) + (t419 * t460 - t451 * t452) * MDP(9) + (-t419 * t459 - t420 * t460 - t450 * t452 - t451 * t453) * MDP(10) + (t352 * t422 + t392 * t526) * MDP(23) + ((-t390 * t507 - t614) * t379 + (-t351 - t354 * t507 + (t390 * t503 - t392 * t507) * qJD(6)) * t422) * MDP(24) + (t366 * t437 + t380 * t418 + t383 * t421 + t405 * t531 + t493 * t539 - t497 * t597 - t518) * MDP(21) + (-qJD(4) * t453 - qJDD(4) * t459) * MDP(12) + (-t380 * t497 - t421 * t493) * MDP(19) + (t379 * t497 + t422 * t493) * MDP(18); (t521 - t524 + t528) * MDP(8) - t477 * MDP(14) + t570 * MDP(15) + (t366 + t613) * MDP(21) + (t365 - t611) * MDP(22) + (t532 - t616) * MDP(28) + (-t507 * t645 ^ 2 - t360 - t615) * MDP(29) - t511 * t640 + (-t601 - pkin(2) * MDP(8) + (MDP(14) * t505 + MDP(6)) * t498) * qJDD(2) + ((t498 * t585 + t501 * t586 + t451) * MDP(14) + (t450 - t565) * MDP(15)) * qJD(4); (t357 * t497 + (-t405 * t451 + t493 * t508 - t497 * t584) * pkin(4) + t636) * MDP(21) + qJDD(4) * MDP(13) + (-t442 * t450 - g(1) * (-t449 * t491 - t490 * t606) - g(2) * (-t447 * t491 + t490 * t563) - g(3) * (-t490 * t502 - t491 * t603) - t540) * MDP(15) + (-t442 * t451 - g(1) * (-t449 * t490 + t491 * t606) - g(2) * (-t447 * t490 - t491 * t563) - g(3) * (-t490 * t603 + t491 * t502) + t557) * MDP(14) + (t532 + t616) * MDP(26) + t544 * MDP(12) - t451 * t450 * MDP(9) + (t614 * t645 + t654) * MDP(24) + (t489 * t354 + t550 * t390 + (-t549 * t645 + t542) * t503 + (t551 * t645 + t520) * t507 + t638) * MDP(28) + (t489 * t353 + t542 * t507 + t550 * t392 - (t503 * t551 + t507 * t549) * t645 + t637) * MDP(29) + (t358 * t497 + (-t451 * t535 - t493 * t504 - t497 * t583) * pkin(4) + t644) * MDP(22) + ((-t450 - t565) * qJD(4) + t570) * MDP(11) + (-t450 ^ 2 + t451 ^ 2) * MDP(10) + t653; (t356 * t497 + t636) * MDP(21) + (t355 * t497 + t644) * MDP(22) + (t392 * t554 + t654) * MDP(24) + (-t554 * t645 + t361 + t616) * MDP(26) + (-pkin(5) * t354 - t356 * t390 + (-pkin(10) * t363 - t355 * t645 + t652) * t503 + (-(-pkin(10) * qJD(6) - t378) * t645 + t520) * t507 + t638) * MDP(28) + (-pkin(5) * t353 - (t355 * t507 + t378 * t503) * t645 - t356 * t392 + t348 * t651 + t527 * pkin(10) + t637) * MDP(29) + t653; t392 * t390 * MDP(23) + (-t390 ^ 2 + t392 ^ 2) * MDP(24) + (-t390 * t645 + t571) * MDP(25) + (-t392 * t645 - t558) * MDP(26) + t363 * MDP(27) + (-t503 * t333 + t337 - t336 * t645 - t348 * t392 - g(1) * (-t417 * t503 + t448 * t507) - g(2) * (-t415 * t503 + t446 * t507) - g(3) * (-t435 * t503 - t507 * t602)) * MDP(28) + (-t507 * t333 - t503 * t338 + t541 * t645 + t348 * t390 - g(1) * (-t417 * t507 - t448 * t503) - g(2) * (-t415 * t507 - t446 * t503) - g(3) * (-t435 * t507 + t503 * t602)) * MDP(29) + (-MDP(25) * t609 - MDP(26) * t392 - MDP(28) * t336 + MDP(29) * t541) * qJD(6);];
tau  = t1;
