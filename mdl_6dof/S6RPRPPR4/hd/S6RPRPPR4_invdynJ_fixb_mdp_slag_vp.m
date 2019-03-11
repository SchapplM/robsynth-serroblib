% Calculate vector of inverse dynamics joint torques for
% S6RPRPPR4
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
%   pkin=[a2,a3,a4,a5,a6,d1,d3,d6,theta2,theta4]';
% MDP [29x1]
%   Minimal dynamic parameter vector (fixed base model)
%   see S6RPRPPR4_convert_par2_MPV_fixb.m
% 
% Output:
% tau [6x1]
%   joint torques of inverse dynamics (contains inertial, gravitational Coriolis and centrifugal forces)

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 02:49
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tau = S6RPRPPR4_invdynJ_fixb_mdp_slag_vp(qJ, qJD, qJDD, g, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(3,1),zeros(10,1),zeros(29,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRPPR4_invdynJ_fixb_mdp_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RPRPPR4_invdynJ_fixb_mdp_slag_vp: qJD has to be [6x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [6 1]), ...
  'S6RPRPPR4_invdynJ_fixb_mdp_slag_vp: qJDD has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RPRPPR4_invdynJ_fixb_mdp_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RPRPPR4_invdynJ_fixb_mdp_slag_vp: pkin has to be [10x1] (double)');
assert(isreal(MDP) && all(size(MDP) == [29 1]), ...
  'S6RPRPPR4_invdynJ_fixb_mdp_slag_vp: MDP has to be [29x1] (double)'); 

%% Symbolic Calculation
% From invdyn_joint_fixb_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 02:48:35
% EndTime: 2019-03-09 02:48:47
% DurationCPUTime: 8.61s
% Computational Cost: add. (4916->528), mult. (11688->646), div. (0->0), fcn. (9036->12), ass. (0->226)
t526 = sin(pkin(9));
t528 = cos(pkin(9));
t531 = sin(qJ(3));
t653 = cos(qJ(3));
t485 = t526 * t653 + t531 * t528;
t525 = sin(pkin(10));
t527 = cos(pkin(10));
t530 = sin(qJ(6));
t533 = cos(qJ(6));
t667 = -t525 * t533 + t527 * t530;
t419 = t667 * t485;
t645 = pkin(7) + qJ(2);
t494 = t645 * t526;
t486 = qJD(1) * t494;
t496 = t645 * t528;
t487 = qJD(1) * t496;
t434 = -t531 * t486 + t653 * t487;
t684 = qJD(3) * t434;
t665 = t485 * qJD(1);
t442 = -t527 * qJD(3) + t525 * t665;
t561 = qJD(3) * t525 + t527 * t665;
t631 = t561 * t530;
t393 = -t533 * t442 + t631;
t600 = t653 * t528;
t507 = qJD(1) * t600;
t621 = t531 * t526;
t599 = qJD(1) * t621;
t471 = -t507 + t599;
t607 = qJD(6) - t471;
t683 = t393 * t607;
t395 = t442 * t530 + t533 * t561;
t682 = t395 * t607;
t681 = t442 * t471;
t610 = qJD(6) * t533;
t611 = qJD(6) * t530;
t617 = t667 * t471 + t525 * t610 - t527 * t611;
t638 = qJDD(1) * pkin(1);
t532 = sin(qJ(1));
t534 = cos(qJ(1));
t669 = g(1) * t532 - g(2) * t534;
t559 = -qJDD(2) + t638 + t669;
t441 = -t531 * t494 + t496 * t653;
t409 = qJD(2) * t485 + qJD(3) * t441;
t680 = t527 * t485 * qJD(5) - t409;
t524 = pkin(9) + qJ(3);
t520 = cos(t524);
t511 = g(3) * t520;
t519 = sin(t524);
t582 = g(1) * t534 + g(2) * t532;
t673 = t582 * t519;
t679 = t673 - t511;
t678 = qJD(3) * t665;
t677 = MDP(15) + MDP(19);
t676 = MDP(16) - MDP(21);
t675 = t561 ^ 2;
t477 = t485 * qJD(3);
t593 = qJDD(1) * t653;
t605 = qJDD(1) * t531;
t577 = t526 * t605 - t528 * t593;
t431 = qJD(1) * t477 + t577;
t601 = qJD(3) * t507 + t526 * t593 + t528 * t605;
t430 = -t599 * qJD(3) + t601;
t515 = pkin(2) * t528 + pkin(1);
t490 = -qJDD(1) * t515 + qJDD(2);
t371 = pkin(3) * t431 - qJ(4) * t430 - qJD(4) * t665 + t490;
t552 = t653 * t486 + t531 * t487;
t606 = qJD(1) * qJD(2);
t656 = qJDD(1) * t645 + t606;
t451 = t656 * t526;
t452 = t656 * t528;
t553 = -t531 * t451 + t452 * t653;
t380 = qJDD(3) * qJ(4) + (qJD(4) - t552) * qJD(3) + t553;
t349 = t371 * t527 - t525 * t380;
t576 = qJDD(5) - t349;
t347 = -pkin(4) * t431 + t576;
t491 = -qJD(1) * t515 + qJD(2);
t407 = pkin(3) * t471 - qJ(4) * t665 + t491;
t425 = qJD(3) * qJ(4) + t434;
t379 = t525 * t407 + t527 * t425;
t365 = t471 * qJ(5) + t379;
t672 = -t365 * t471 + t347;
t427 = pkin(3) * t665 + qJ(4) * t471;
t388 = t525 * t427 - t527 * t552;
t372 = qJ(5) * t665 + t388;
t613 = qJD(4) * t527;
t671 = -t372 + t613;
t670 = -t653 * t494 - t531 * t496;
t615 = t520 * pkin(3) + t519 * qJ(4);
t413 = qJDD(3) * t525 + t430 * t527;
t612 = qJD(5) * t561;
t668 = qJ(5) * t413 + t612;
t666 = MDP(4) * t528 - MDP(5) * t526;
t664 = qJ(2) * qJDD(1);
t663 = MDP(17) + MDP(20);
t428 = -qJDD(6) + t431;
t482 = t525 * t530 + t527 * t533;
t618 = t607 * t482;
t661 = -t428 * t667 + t607 * t618;
t646 = g(3) * t519;
t658 = t582 * t520 + t646;
t639 = qJ(5) * t527;
t654 = pkin(4) + pkin(5);
t657 = t525 * t654 - t639;
t412 = -t527 * qJDD(3) + t430 * t525;
t590 = -t533 * t412 + t530 * t413;
t354 = t395 * qJD(6) + t590;
t587 = -t653 * t451 - t531 * t452 - t684;
t547 = qJDD(3) * pkin(3) - qJDD(4) + t587;
t538 = t547 + t668;
t346 = -t412 * t654 + t538;
t655 = t346 + t679;
t466 = t471 ^ 2;
t652 = pkin(4) * t412;
t651 = pkin(4) * t527;
t650 = pkin(8) * t525;
t644 = -pkin(8) + qJ(4);
t641 = qJ(4) * t527;
t636 = t393 * t665;
t635 = t395 * t665;
t633 = t431 * t525;
t632 = t431 * t527;
t630 = t519 * t532;
t629 = t519 * t534;
t628 = t520 * t534;
t627 = t525 * t532;
t624 = t527 * t534;
t623 = t645 * t534;
t620 = t532 * t527;
t619 = t534 * t525;
t350 = t525 * t371 + t527 * t380;
t551 = t600 - t621;
t476 = t551 * qJD(3);
t397 = pkin(3) * t477 - qJ(4) * t476 - qJD(4) * t485;
t408 = t551 * qJD(2) + qJD(3) * t670;
t369 = t525 * t397 + t527 * t408;
t429 = -pkin(3) * t551 - qJ(4) * t485 - t515;
t391 = t525 * t429 + t527 * t441;
t616 = (g(1) * t624 + g(2) * t620) * t519;
t614 = t526 ^ 2 + t528 ^ 2;
t558 = qJD(3) * pkin(3) - qJD(4) - t552;
t609 = -qJD(4) - t558;
t543 = qJ(5) * t561 + t558;
t377 = pkin(4) * t442 - t543;
t608 = qJD(4) - t377;
t604 = qJ(4) * t632;
t603 = t530 * t412 + t533 * t413 + t442 * t610;
t383 = -qJ(5) * t551 + t391;
t602 = pkin(3) * t628 + qJ(4) * t629 + t534 * t515;
t597 = qJ(5) * t525 + pkin(3);
t348 = -t538 + t652;
t595 = -t348 - t511;
t594 = t547 - t511;
t592 = t614 * qJD(1) ^ 2;
t341 = -pkin(8) * t413 - t431 * t654 + t576;
t345 = t431 * qJ(5) + t471 * qJD(5) + t350;
t342 = pkin(8) * t412 + t345;
t591 = t533 * t341 - t530 * t342;
t400 = t525 * t408;
t368 = t397 * t527 - t400;
t378 = t407 * t527 - t525 * t425;
t423 = t525 * t552;
t387 = t427 * t527 + t423;
t435 = t525 * t441;
t390 = t429 * t527 - t435;
t589 = t607 ^ 2;
t588 = qJD(5) * t525 - t471 * t657 + t434;
t356 = t477 * qJ(5) - qJD(5) * t551 + t369;
t586 = 0.2e1 * t614;
t585 = g(1) * t630 - g(2) * t629;
t455 = t520 * t627 + t624;
t457 = t520 * t619 - t620;
t584 = -g(1) * t455 + g(2) * t457;
t456 = t520 * t620 - t619;
t458 = t520 * t624 + t627;
t583 = g(1) * t456 - g(2) * t458;
t580 = t482 * t428 - t607 * t617;
t579 = qJD(5) - t378;
t578 = pkin(4) * t525 - t639;
t575 = t530 * t341 + t533 * t342;
t574 = -t345 * t525 + t347 * t527;
t573 = t348 * t485 + t377 * t476;
t572 = -t349 * t527 - t350 * t525;
t355 = -pkin(8) * t561 - t471 * t654 + t579;
t357 = pkin(8) * t442 + t365;
t343 = t355 * t533 - t357 * t530;
t344 = t355 * t530 + t357 * t533;
t361 = t435 + (-pkin(8) * t485 - t429) * t527 + t654 * t551;
t370 = t485 * t650 + t383;
t571 = t361 * t533 - t370 * t530;
t570 = t361 * t530 + t370 * t533;
t569 = -t378 * t525 + t379 * t527;
t568 = -t476 * t558 - t485 * t547;
t564 = t455 * t533 - t456 * t530;
t563 = t455 * t530 + t456 * t533;
t562 = qJ(4) * t413 + qJD(4) * t561;
t560 = t597 + t651;
t495 = t644 * t527;
t550 = qJD(4) * t525 - qJD(6) * t495 + t423 - (pkin(8) * t471 - t427) * t527 + t654 * t665;
t493 = t644 * t525;
t549 = qJD(6) * t493 + t471 * t650 + t671;
t420 = t482 * t485;
t548 = t561 * t611 - t603;
t541 = -t412 * t641 - t442 * t613 - t658;
t540 = (-g(1) * (-t515 - t615) - g(2) * t645) * t532;
t539 = t586 * t606 - t582;
t537 = -t547 - t679;
t499 = qJ(4) * t628;
t497 = t532 * t520 * qJ(4);
t478 = t527 * t654 + t597;
t415 = t457 * t530 + t458 * t533;
t414 = t457 * t533 - t458 * t530;
t392 = t485 * t578 - t670;
t389 = -t471 * t578 + t434;
t386 = -t485 * t657 + t670;
t385 = pkin(4) * t551 - t390;
t382 = qJD(6) * t420 + t476 * t667;
t381 = -qJD(6) * t419 + t476 * t482;
t374 = -pkin(4) * t665 - t387;
t367 = t476 * t578 - t680;
t364 = -pkin(4) * t471 + t579;
t362 = -t442 * t654 + t543;
t360 = -pkin(4) * t477 - t368;
t359 = -t476 * t657 + t680;
t352 = t476 * t650 + t356;
t351 = t400 + (-pkin(8) * t476 - t397) * t527 - t654 * t477;
t1 = [(-g(1) * t623 - g(2) * t602 + t349 * t390 + t350 * t391 + t378 * t368 + t379 * t369 - t409 * t558 + t547 * t670 + t540) * MDP(18) + (t350 * t551 - t369 * t471 - t379 * t477 - t391 * t431 + t409 * t561 - t413 * t670 + t527 * t568 + t584) * MDP(16) + (-t349 * t551 + t368 * t471 + t378 * t477 + t390 * t431 + t409 * t442 - t412 * t670 + t525 * t568 + t583) * MDP(15) + (-qJD(3) * t409 + qJDD(3) * t670 - t431 * t515 + t477 * t491 - t490 * t551 + t520 * t669) * MDP(13) + (t345 * t383 + t365 * t356 + t348 * t392 + t377 * t367 + t347 * t385 + t364 * t360 - g(1) * (-pkin(4) * t456 - qJ(5) * t455 + t623) - g(2) * (pkin(4) * t458 + qJ(5) * t457 + t602) + t540) * MDP(22) + (-t368 * t561 - t369 * t442 - t390 * t413 - t391 * t412 + t572 * t485 + (-t378 * t527 - t379 * t525) * t476 + t585) * MDP(17) + (-t356 * t442 + t360 * t561 - t383 * t412 + t385 * t413 + t574 * t485 + (t364 * t527 - t365 * t525) * t476 + t585) * MDP(20) + (t381 * t395 - t420 * t548) * MDP(23) + (-t354 * t420 - t381 * t393 - t382 * t395 + t419 * t548) * MDP(24) + t666 * (t559 + t638) + (pkin(1) * t559 + (t614 * t664 + t539) * qJ(2)) * MDP(7) + (t586 * t664 + t539) * MDP(6) + qJDD(1) * MDP(1) + (t430 * t485 + t476 * t665) * MDP(8) + (-t345 * t551 + t356 * t471 + t365 * t477 - t367 * t561 + t383 * t431 - t392 * t413 - t527 * t573 - t584) * MDP(21) + (-qJD(3) * t477 + qJDD(3) * t551) * MDP(11) + (t347 * t551 - t360 * t471 - t364 * t477 + t367 * t442 - t385 * t431 + t392 * t412 + t525 * t573 + t583) * MDP(19) + ((t351 * t533 - t352 * t530) * t607 - t571 * t428 + t591 * t551 - t343 * t477 + t359 * t393 + t386 * t354 + t346 * t419 + t362 * t382 + g(1) * t563 - g(2) * t415 + (-t344 * t551 - t570 * t607) * qJD(6)) * MDP(28) + (-t428 * t551 - t477 * t607) * MDP(27) + (-t354 * t551 - t382 * t607 + t393 * t477 + t419 * t428) * MDP(26) + (t381 * t607 - t395 * t477 - t420 * t428 - t548 * t551) * MDP(25) + (-(t351 * t530 + t352 * t533) * t607 + t570 * t428 - t575 * t551 + t344 * t477 + t359 * t395 - t386 * t548 + t346 * t420 + t362 * t381 + g(1) * t564 - g(2) * t414 + (-t343 * t551 - t571 * t607) * qJD(6)) * MDP(29) + (t430 * t551 - t431 * t485 - t471 * t476 - t477 * t665) * MDP(9) + t669 * MDP(2) + (qJD(3) * t476 + qJDD(3) * t485) * MDP(10) + t582 * MDP(3) + (-qJD(3) * t408 - qJDD(3) * t441 - t430 * t515 + t476 * t491 + t485 * t490 - t585) * MDP(14); -MDP(6) * t592 + (-qJ(2) * t592 - t559) * MDP(7) + (t577 + 0.2e1 * t678) * MDP(13) + ((-t471 - t599) * qJD(3) + t601) * MDP(14) + (t471 * t569 + t558 * t665 - t572 - t669) * MDP(18) + (-t377 * t665 + (t364 * t525 + t365 * t527) * t471 - t574 - t669) * MDP(22) + (t580 + t636) * MDP(28) + (t635 + t661) * MDP(29) + t663 * (-t413 * t527 - t525 * t412 + (-t442 * t527 + t525 * t561) * t471) - t666 * qJDD(1) + t677 * (-t442 * t665 - t466 * t525 + t632) - t676 * (t466 * t527 + t561 * t665 + t633); -t466 * MDP(9) + ((t471 - t599) * qJD(3) + t601) * MDP(10) - t577 * MDP(11) + qJDD(3) * MDP(12) + (t587 + t679 + t684) * MDP(13) + (t491 * t471 - t553 + t658) * MDP(14) + (-qJ(4) * t633 - pkin(3) * t412 - t434 * t442 + t594 * t527 + (t525 * t609 - t387) * t471 + t616) * MDP(15) + (-t604 - pkin(3) * t413 - t434 * t561 + (t527 * t609 + t388) * t471 + (-t673 - t594) * t525) * MDP(16) + (t387 * t561 + t388 * t442 + (-t378 * t471 + t350) * t527 + (-t379 * t471 - t349 + t562) * t525 + t541) * MDP(17) + (t547 * pkin(3) - t379 * t388 - t378 * t387 + t558 * t434 - g(1) * (-pkin(3) * t629 + t499) - g(2) * (-pkin(3) * t630 + t497) - g(3) * t615 + t569 * qJD(4) + (-t349 * t525 + t350 * t527) * qJ(4)) * MDP(18) + (t374 * t471 - t389 * t442 - t412 * t560 + t595 * t527 + (-qJ(4) * t431 - qJD(5) * t442 - t471 * t608) * t525 + t616) * MDP(19) + (t372 * t442 - t374 * t561 + (t364 * t471 + t345) * t527 + (t562 + t672) * t525 + t541) * MDP(20) + (t604 + t389 * t561 + t413 * t560 + (t527 * t608 - t372) * t471 + (t673 + t595 + t612) * t525) * MDP(21) + (t345 * t641 - t377 * t389 - t364 * t374 - g(1) * t499 - g(2) * t497 - g(3) * (t520 * t651 + t615) + t671 * t365 + (qJ(4) * t347 - qJ(5) * t511 + qJD(4) * t364 - qJD(5) * t377) * t525 + (-t348 + t673) * t560) * MDP(22) + (-t395 * t618 + t548 * t667) * MDP(23) + (t354 * t667 + t393 * t618 - t395 * t617 + t482 * t548) * MDP(24) + (t635 - t661) * MDP(25) + (t580 - t636) * MDP(26) + (-(t493 * t533 - t495 * t530) * t428 + t478 * t354 - (t530 * t549 - t533 * t550) * t607 + t588 * t393 + t617 * t362 + t655 * t482) * MDP(28) + ((t493 * t530 + t495 * t533) * t428 - t478 * t548 - (t530 * t550 + t533 * t549) * t607 + t588 * t395 - t618 * t362 - t655 * t667) * MDP(29) + (-t491 * MDP(13) - t378 * MDP(15) + t379 * MDP(16) + t364 * MDP(19) - t365 * MDP(21) + t607 * MDP(27) + t343 * MDP(28) - t344 * MDP(29) + t471 * MDP(8) + MDP(9) * t665) * t665; (t378 * t561 + t379 * t442 + t537) * MDP(18) + (-t364 * t561 + t365 * t442 + t537 + t652 - t668) * MDP(22) + (-t354 - t682) * MDP(28) + (t548 + t683) * MDP(29) + t663 * (-t442 ^ 2 - t675) + t677 * (t471 * t561 + t412) + t676 * (t413 - t681); (t442 * t561 - t577 - t678) * MDP(19) + (t413 + t681) * MDP(20) + (-t466 - t675) * MDP(21) + (-g(1) * t457 - g(2) * t455 + t377 * t561 - t525 * t646 + t672) * MDP(22) + (-t393 * t561 - t533 * t428 - t530 * t589) * MDP(28) + (-t395 * t561 + t530 * t428 - t533 * t589) * MDP(29); t395 * t393 * MDP(23) + (-t393 ^ 2 + t395 ^ 2) * MDP(24) + (t603 + t683) * MDP(25) + (-t590 + t682) * MDP(26) - t428 * MDP(27) + (-g(1) * t414 - g(2) * t564 + t344 * t607 - t362 * t395 + t591) * MDP(28) + (g(1) * t415 + g(2) * t563 + t343 * t607 + t362 * t393 - t575) * MDP(29) + (MDP(28) * t667 + MDP(29) * t482) * t646 + (-MDP(25) * t631 - MDP(26) * t395 - MDP(28) * t344 - MDP(29) * t343) * qJD(6);];
tau  = t1;
