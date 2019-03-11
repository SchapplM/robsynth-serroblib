% Calculate vector of inverse dynamics joint torques for
% S6RRPPRR1
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
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d5,d6,theta3]';
% MDP [30x1]
%   Minimal dynamic parameter vector (fixed base model)
%   see S6RRPPRR1_convert_par2_MPV_fixb.m
% 
% Output:
% tau [6x1]
%   joint torques of inverse dynamics (contains inertial, gravitational Coriolis and centrifugal forces)

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 08:48
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tau = S6RRPPRR1_invdynJ_fixb_mdp_slag_vp(qJ, qJD, qJDD, g, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(3,1),zeros(10,1),zeros(30,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPPRR1_invdynJ_fixb_mdp_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRPPRR1_invdynJ_fixb_mdp_slag_vp: qJD has to be [6x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [6 1]), ...
  'S6RRPPRR1_invdynJ_fixb_mdp_slag_vp: qJDD has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RRPPRR1_invdynJ_fixb_mdp_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RRPPRR1_invdynJ_fixb_mdp_slag_vp: pkin has to be [10x1] (double)');
assert(isreal(MDP) && all(size(MDP) == [30 1]), ...
  'S6RRPPRR1_invdynJ_fixb_mdp_slag_vp: MDP has to be [30x1] (double)'); 

%% Symbolic Calculation
% From invdyn_joint_fixb_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 08:48:16
% EndTime: 2019-03-09 08:48:26
% DurationCPUTime: 7.95s
% Computational Cost: add. (5158->490), mult. (11882->608), div. (0->0), fcn. (9070->12), ass. (0->211)
t530 = sin(pkin(10));
t531 = cos(pkin(10));
t536 = sin(qJ(2));
t540 = cos(qJ(2));
t485 = t530 * t540 + t531 * t536;
t471 = t485 * qJD(2);
t598 = qJDD(1) * t540;
t599 = qJDD(1) * t536;
t570 = t530 * t599 - t531 * t598;
t427 = qJD(1) * t471 + t570;
t600 = qJD(1) * qJD(2);
t593 = t536 * t600;
t498 = t530 * t593;
t553 = qJDD(1) * t485 - t498;
t592 = t540 * t600;
t428 = t531 * t592 + t553;
t612 = qJD(1) * t540;
t594 = t531 * t612;
t613 = qJD(1) * t536;
t469 = t530 * t613 - t594;
t472 = t485 * qJD(1);
t535 = sin(qJ(5));
t539 = cos(qJ(5));
t609 = qJD(5) * t539;
t610 = qJD(5) * t535;
t365 = t535 * t427 + t539 * t428 + t469 * t609 - t472 * t610;
t521 = qJDD(2) - qJDD(5);
t522 = qJD(2) - qJD(5);
t534 = sin(qJ(6));
t538 = cos(qJ(6));
t605 = qJD(6) * t538;
t596 = t538 * t365 - t534 * t521 - t522 * t605;
t606 = qJD(6) * t534;
t653 = t469 * t535 + t539 * t472;
t348 = -t606 * t653 + t596;
t396 = -t522 * t534 + t538 * t653;
t507 = t538 * t521;
t349 = t396 * qJD(6) + t365 * t534 + t507;
t366 = qJD(5) * t653 - t539 * t427 + t428 * t535;
t364 = qJDD(6) + t366;
t620 = t538 * t364;
t656 = -t539 * t469 + t472 * t535;
t661 = qJD(6) + t656;
t587 = t606 * t661 - t620;
t574 = t534 * t656 * t661 + t587;
t629 = t656 * t522;
t630 = t653 * t522;
t632 = t396 * t653;
t633 = t396 * t661;
t394 = t538 * t522 + t534 * t653;
t634 = t394 * t653;
t635 = t394 * t661;
t637 = t348 * t534;
t621 = t534 * t364;
t671 = t661 * t538;
t667 = t661 * t671 + t621;
t679 = -((t349 + t633) * t534 - (t348 - t635) * t538) * MDP(25) + (t396 * t671 + t637) * MDP(24) + (-t632 + t667) * MDP(26) + (t365 - t629) * MDP(19) - (t366 + t630) * MDP(20) - (t574 - t634) * MDP(27) + (t653 ^ 2 - t656 ^ 2) * MDP(18) - t521 * MDP(21) + (MDP(17) * t656 - MDP(28) * t661) * t653;
t533 = -qJ(3) - pkin(7);
t496 = t533 * t536;
t489 = qJD(1) * t496;
t497 = t533 * t540;
t490 = qJD(1) * t497;
t624 = t530 * t490;
t433 = t531 * t489 + t624;
t602 = qJD(4) - t433;
t665 = pkin(5) * t653;
t493 = -qJD(1) * pkin(1) - pkin(2) * t612 + qJD(3);
t399 = t469 * pkin(3) - t472 * qJ(4) + t493;
t383 = -pkin(4) * t469 - t399;
t350 = pkin(5) * t656 - pkin(9) * t653 + t383;
t481 = qJD(2) * pkin(2) + t489;
t425 = t481 * t531 + t624;
t573 = qJD(4) - t425;
t644 = pkin(8) * t472;
t646 = -pkin(3) - pkin(4);
t386 = qJD(2) * t646 + t573 - t644;
t623 = t531 * t490;
t426 = t530 * t481 - t623;
t415 = qJD(2) * qJ(4) + t426;
t645 = pkin(8) * t469;
t392 = t415 + t645;
t357 = t386 * t535 + t392 * t539;
t354 = -pkin(9) * t522 + t357;
t343 = t350 * t538 - t354 * t534;
t664 = t343 * t653;
t344 = t350 * t534 + t354 * t538;
t663 = t344 * t653;
t662 = -t644 + t602;
t541 = cos(qJ(1));
t523 = qJ(2) + pkin(10);
t517 = sin(t523);
t518 = cos(t523);
t564 = t517 * t535 + t518 * t539;
t443 = t564 * t541;
t589 = qJD(2) * t533;
t465 = -qJD(3) * t536 + t540 * t589;
t422 = qJDD(2) * pkin(2) + qJD(1) * t465 + qJDD(1) * t496;
t464 = qJD(3) * t540 + t536 * t589;
t431 = qJD(1) * t464 - qJDD(1) * t497;
t376 = t531 * t422 - t530 * t431;
t591 = -qJDD(4) + t376;
t362 = -pkin(8) * t428 + qJDD(2) * t646 - t591;
t525 = qJD(2) * qJD(4);
t377 = t530 * t422 + t531 * t431;
t595 = qJDD(2) * qJ(4) + t377;
t372 = t525 + t595;
t363 = pkin(8) * t427 + t372;
t560 = t535 * t362 + t539 * t363 + t386 * t609 - t392 * t610;
t537 = sin(qJ(1));
t441 = t564 * t537;
t652 = -t517 * t539 + t518 * t535;
t575 = g(2) * t441 - g(3) * t652;
t659 = g(1) * t443 + t383 * t656 - t560 + t575;
t440 = t652 * t537;
t442 = t652 * t541;
t559 = g(1) * t442 + g(2) * t440 + g(3) * t564;
t563 = -t539 * t362 + t535 * t363 + t386 * t610 + t392 * t609;
t658 = -t383 * t653 + t559 - t563;
t513 = -pkin(2) * t531 - pkin(3);
t508 = -pkin(4) + t513;
t511 = pkin(2) * t530 + qJ(4);
t616 = t535 * t508 + t539 * t511;
t520 = g(2) * t541;
t654 = g(1) * t537 - t520;
t641 = g(1) * t541;
t577 = g(2) * t537 + t641;
t651 = t427 * pkin(3) - t428 * qJ(4) - t472 * qJD(4);
t401 = pkin(2) * t613 + t472 * pkin(3) + t469 * qJ(4);
t387 = -pkin(4) * t472 - t401;
t445 = -pkin(9) + t616;
t342 = pkin(5) * t521 + t563;
t556 = -t342 + t559;
t650 = t661 * (-pkin(9) * t656 + qJD(6) * t445 + t387 - t665) + t556;
t649 = t661 * (t661 * pkin(9) + t665) - t556;
t648 = -g(3) * t518 - t399 * t472 + t577 * t517 + t591;
t402 = t464 * t530 - t531 * t465;
t611 = qJD(2) * t536;
t622 = t531 * t540;
t474 = qJD(2) * t622 - t530 * t611;
t388 = -pkin(8) * t474 + t402;
t403 = t531 * t464 + t530 * t465;
t389 = pkin(8) * t471 + t403;
t438 = -t531 * t496 - t497 * t530;
t404 = -pkin(8) * t485 + t438;
t439 = t530 * t496 - t531 * t497;
t484 = t530 * t536 - t622;
t405 = pkin(8) * t484 + t439;
t568 = t404 * t539 - t405 * t535;
t346 = qJD(5) * t568 + t388 * t535 + t389 * t539;
t356 = t386 * t539 - t392 * t535;
t353 = pkin(5) * t522 - t356;
t519 = t540 * pkin(2);
t514 = t519 + pkin(1);
t424 = t484 * pkin(3) - t485 * qJ(4) - t514;
t393 = -pkin(4) * t484 - t424;
t430 = t484 * t535 + t485 * t539;
t566 = t539 * t484 - t485 * t535;
t360 = -pkin(5) * t566 - pkin(9) * t430 + t393;
t371 = t404 * t535 + t405 * t539;
t378 = qJD(5) * t566 + t471 * t535 + t474 * t539;
t583 = -pkin(9) * t521 + qJD(6) * t350 + t560;
t647 = t342 * t430 + t353 * t378 - t371 * t364 - (qJD(6) * t360 + t346) * t661 + t566 * t583 + t641;
t467 = t472 ^ 2;
t643 = g(1) * t441;
t639 = g(3) * t540;
t529 = qJDD(1) * pkin(1);
t638 = qJDD(2) * pkin(3);
t636 = t353 * t430;
t432 = t489 * t530 - t623;
t397 = t432 + t645;
t565 = t508 * t539 - t511 * t535;
t618 = qJD(5) * t565 - t397 * t535 + t662 * t539;
t617 = t616 * qJD(5) + t397 * t539 + t662 * t535;
t527 = t536 ^ 2;
t615 = -t540 ^ 2 + t527;
t608 = qJD(6) * t354;
t607 = qJD(6) * t653;
t597 = pkin(2) * t611;
t580 = t522 ^ 2;
t579 = t539 * t522;
t578 = qJDD(3) - t529 + (t593 - t598) * pkin(2);
t572 = pkin(3) * t518 + qJ(4) * t517;
t571 = t360 * t364 + t643;
t569 = -t608 - t520;
t562 = -0.2e1 * pkin(1) * t600 - pkin(7) * qJDD(2);
t561 = t378 * t538 - t430 * t606;
t391 = t471 * pkin(3) - t474 * qJ(4) - t485 * qJD(4) + t597;
t367 = t578 + t651;
t557 = t578 - t654;
t554 = t575 - t583;
t380 = -pkin(4) * t471 - t391;
t552 = -pkin(9) * t364 + (t353 + t356) * t661;
t542 = qJD(2) ^ 2;
t551 = -pkin(7) * t542 + 0.2e1 * t529 + t654;
t543 = qJD(1) ^ 2;
t550 = pkin(1) * t543 - pkin(7) * qJDD(1) + t577;
t351 = -pkin(4) * t427 - t367;
t547 = -t445 * t364 + (-t353 - t618) * t661;
t546 = t402 * t472 - t403 * t469 - t439 * t427 + t428 * t438 - t577;
t499 = t541 * t514;
t444 = pkin(5) - t565;
t435 = t443 * t538 - t534 * t537;
t434 = -t443 * t534 - t537 * t538;
t407 = -qJD(2) * pkin(3) + t573;
t379 = qJD(5) * t430 - t539 * t471 + t474 * t535;
t375 = -t591 - t638;
t347 = qJD(5) * t371 - t388 * t539 + t389 * t535;
t345 = pkin(5) * t379 - pkin(9) * t378 + t380;
t340 = pkin(5) * t366 - pkin(9) * t365 + t351;
t339 = t538 * t340;
t1 = [(-g(2) * t434 - t344 * t379 + t347 * t396 - t568 * t348 + (-(-qJD(6) * t371 + t345) * t661 + (t340 - t608) * t566 - qJD(6) * t636 - t571) * t534 + t647 * t538) * MDP(30) + (-g(2) * t435 - t339 * t566 + t343 * t379 + t347 * t394 - t568 * t349 + (t345 * t661 + (t354 * t566 - t371 * t661 + t636) * qJD(6) + t571) * t538 + t647 * t534) * MDP(29) + (-t348 * t566 + t379 * t396 + t430 * t620 + t561 * t661) * MDP(26) + (-t430 * t621 + t349 * t566 - t379 * t394 + (-t378 * t534 - t430 * t605) * t661) * MDP(27) + (-t364 * t566 + t379 * t661) * MDP(28) + (-t376 * t485 - t377 * t484 - t425 * t474 - t426 * t471 + t546) * MDP(11) + (-t372 * t484 + t375 * t485 + t407 * t474 - t415 * t471 + t546) * MDP(14) + (t365 * t430 + t378 * t653) * MDP(17) + (-g(1) * t440 + g(2) * t442 + t346 * t522 + t351 * t430 + t365 * t393 + t371 * t521 + t378 * t383 + t380 * t653) * MDP(23) + (-g(2) * t443 + t347 * t522 - t351 * t566 + t366 * t393 + t379 * t383 + t380 * t656 - t521 * t568 + t643) * MDP(22) + (t365 * t566 - t366 * t430 - t378 * t656 - t379 * t653) * MDP(18) + (t379 * t522 - t521 * t566) * MDP(20) + (qJD(2) * t403 + qJDD(2) * t439 - t367 * t485 - t391 * t472 - t399 * t474 - t424 * t428 + t517 * t654) * MDP(15) + (-qJD(2) * t402 - qJDD(2) * t438 + t367 * t484 + t391 * t469 + t399 * t471 + t424 * t427 + t518 * t654) * MDP(13) + t654 * MDP(2) + qJDD(1) * MDP(1) + (t348 * t430 * t538 + t396 * t561) * MDP(24) + (t536 * t562 + t540 * t551) * MDP(9) + (-t536 * t551 + t540 * t562) * MDP(10) + (-g(2) * t499 + t367 * t424 + t372 * t439 + t375 * t438 + t399 * t391 + t407 * t402 + t415 * t403 + (g(1) * t533 - g(2) * t572) * t541 + (-g(1) * (-t514 - t572) + g(2) * t533) * t537) * MDP(16) + (qJDD(2) * t536 + t540 * t542) * MDP(6) + (qJDD(2) * t540 - t536 * t542) * MDP(7) + t577 * MDP(3) + (qJDD(1) * t527 + 0.2e1 * t536 * t592) * MDP(4) + (t377 * t439 + t426 * t403 - t376 * t438 - t425 * t402 - t578 * t514 + t493 * t597 - g(1) * (-t514 * t537 - t533 * t541) - g(2) * (-t533 * t537 + t499)) * MDP(12) + 0.2e1 * (t536 * t598 - t600 * t615) * MDP(5) + ((-t394 * t538 - t396 * t534) * t378 + (-t637 - t349 * t538 + (t394 * t534 - t396 * t538) * qJD(6)) * t430) * MDP(25) + (-t378 * t522 - t430 * t521) * MDP(19); (g(3) * t536 + t540 * t550) * MDP(10) + (t444 * t348 + t617 * t396 + t534 * t650 + t547 * t538 - t663) * MDP(30) + (t444 * t349 + t617 * t394 + t547 * t534 - t538 * t650 + t664) * MDP(29) + (-t387 * t656 - t521 * t565 + t522 * t617 - t658) * MDP(22) + (-t387 * t653 + t521 * t616 + t522 * t618 - t659) * MDP(23) + ((t426 - t432) * t472 + (-t425 + t433) * t469 + (-t427 * t530 - t428 * t531) * pkin(2)) * MDP(11) + (-g(3) * t517 - qJD(2) * t433 + qJDD(2) * t511 - t399 * t469 + t401 * t472 - t518 * t577 + 0.2e1 * t525 + t595) * MDP(15) + (-t427 * t511 + t428 * t513 + (t415 - t432) * t472 + (t407 - t602) * t469) * MDP(14) + qJDD(2) * MDP(8) + MDP(7) * t598 + MDP(6) * t599 + (-MDP(4) * t536 * t540 + MDP(5) * t615) * t543 + (qJD(2) * t432 - t401 * t469 + (pkin(3) - t513) * qJDD(2) + t648) * MDP(13) + (t372 * t511 + t375 * t513 - t399 * t401 - t407 * t432 - g(3) * (t519 + t572) + t602 * t415 + t577 * (pkin(2) * t536 + pkin(3) * t517 - qJ(4) * t518)) * MDP(16) + (t536 * t550 - t639) * MDP(9) + (t425 * t432 - t426 * t433 + (-t639 + t376 * t531 + t377 * t530 + (-qJD(1) * t493 + t577) * t536) * pkin(2)) * MDP(12) - t679; (t425 * t472 + t426 * t469 + t557) * MDP(12) + t570 * MDP(13) + (-t530 * t598 - t531 * t599 + t498) * MDP(15) + (-t407 * t472 + t415 * t469 + t557 + t651) * MDP(16) + (-t366 + t630) * MDP(22) + (-t365 - t629) * MDP(23) + (t574 + t634) * MDP(29) + (t632 + t667) * MDP(30) + ((t530 * t612 + t531 * t613 + t472) * MDP(13) + (t469 - t594) * MDP(15)) * qJD(2) + (MDP(14) + MDP(11)) * (-t469 ^ 2 - t467); (t469 * t472 - qJDD(2)) * MDP(13) + ((t469 + t594) * qJD(2) + t553) * MDP(14) + (-t467 - t542) * MDP(15) + (-qJD(2) * t415 - t638 - t648) * MDP(16) + (-t472 * t656 - t521 * t539 - t535 * t580) * MDP(22) + (-t472 * t653 + t521 * t535 - t539 * t580) * MDP(23) + (-t539 * t349 + (-t538 * t472 + t534 * t579) * t661 + (-t394 * t522 - t605 * t661 - t621) * t535) * MDP(29) + (-t539 * t348 + (t534 * t472 + t538 * t579) * t661 + (-t396 * t522 + t587) * t535) * MDP(30); (-t357 * t522 + t658) * MDP(22) + (-t356 * t522 + t659) * MDP(23) + (-pkin(5) * t349 - t357 * t394 + t552 * t534 - t538 * t649 - t664) * MDP(29) + (-pkin(5) * t348 - t357 * t396 + t534 * t649 + t552 * t538 + t663) * MDP(30) + t679; t396 * t394 * MDP(24) + (-t394 ^ 2 + t396 ^ 2) * MDP(25) + (t596 + t635) * MDP(26) + (-t507 + t633) * MDP(27) + t364 * MDP(28) + (-g(1) * t434 + t344 * t661 - t353 * t396 + t339) * MDP(29) + (g(1) * t435 + t343 * t661 + t353 * t394) * MDP(30) + (-MDP(27) * t607 + MDP(29) * t569 + MDP(30) * t554) * t538 + (-MDP(26) * t607 + (qJD(6) * t522 - t365) * MDP(27) + t554 * MDP(29) + (-t340 - t569) * MDP(30)) * t534;];
tau  = t1;
