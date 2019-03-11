% Calculate vector of inverse dynamics joint torques for
% S6RPRPRR3
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
% pkin [11x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d3,d5,d6,theta2,theta4]';
% MDP [29x1]
%   Minimal dynamic parameter vector (fixed base model)
%   see S6RPRPRR3_convert_par2_MPV_fixb.m
% 
% Output:
% tau [6x1]
%   joint torques of inverse dynamics (contains inertial, gravitational Coriolis and centrifugal forces)

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 03:43
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tau = S6RPRPRR3_invdynJ_fixb_mdp_slag_vp(qJ, qJD, qJDD, g, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(3,1),zeros(11,1),zeros(29,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRPRR3_invdynJ_fixb_mdp_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RPRPRR3_invdynJ_fixb_mdp_slag_vp: qJD has to be [6x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [6 1]), ...
  'S6RPRPRR3_invdynJ_fixb_mdp_slag_vp: qJDD has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RPRPRR3_invdynJ_fixb_mdp_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6RPRPRR3_invdynJ_fixb_mdp_slag_vp: pkin has to be [11x1] (double)');
assert(isreal(MDP) && all(size(MDP) == [29 1]), ...
  'S6RPRPRR3_invdynJ_fixb_mdp_slag_vp: MDP has to be [29x1] (double)'); 

%% Symbolic Calculation
% From invdyn_joint_fixb_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 03:42:42
% EndTime: 2019-03-09 03:42:55
% DurationCPUTime: 9.55s
% Computational Cost: add. (5469->550), mult. (11793->728), div. (0->0), fcn. (8729->18), ass. (0->229)
t565 = cos(qJ(3));
t638 = qJD(1) * t565;
t534 = -qJD(5) + t638;
t531 = -qJD(6) + t534;
t555 = sin(pkin(11));
t557 = cos(pkin(11));
t627 = t557 * qJD(3);
t561 = sin(qJ(3));
t639 = qJD(1) * t561;
t501 = t555 * t639 - t627;
t636 = qJD(3) * t555;
t503 = t557 * t639 + t636;
t560 = sin(qJ(5));
t564 = cos(qJ(5));
t441 = t564 * t501 + t503 * t560;
t563 = cos(qJ(6));
t440 = t501 * t560 - t503 * t564;
t559 = sin(qJ(6));
t659 = t440 * t559;
t678 = -t563 * t441 + t659;
t677 = t531 * t678;
t556 = sin(pkin(10));
t538 = pkin(1) * t556 + pkin(7);
t526 = t538 * qJD(1);
t634 = qJD(3) * t565;
t524 = t538 * qJDD(1);
t672 = -qJD(2) * qJD(3) - t524;
t621 = -t526 * t634 + t561 * t672;
t581 = -qJDD(3) * pkin(3) + qJDD(4) - t621;
t420 = -qJDD(2) * t565 + t581;
t552 = qJ(1) + pkin(10);
t545 = sin(t552);
t547 = cos(t552);
t604 = g(1) * t547 + g(2) * t545;
t588 = t604 * t561;
t577 = -g(3) * t565 + t588;
t680 = t420 - t577;
t484 = qJD(2) * t565 - t561 * t526;
t599 = pkin(3) * t561 - qJ(4) * t565;
t515 = t599 * qJD(1);
t433 = t557 * t484 + t555 * t515;
t620 = t555 * t638;
t410 = -pkin(8) * t620 + t433;
t679 = qJD(4) * t557 - t410;
t676 = t440 * t534;
t675 = t441 * t534;
t593 = t440 * t563 + t441 * t559;
t674 = t531 * t593;
t656 = t555 * t560;
t509 = -t564 * t557 + t656;
t586 = t565 * t509;
t644 = qJD(1) * t586 - t509 * qJD(5);
t510 = t555 * t564 + t557 * t560;
t585 = t510 * t565;
t643 = -qJD(1) * t585 + t510 * qJD(5);
t485 = t561 * qJD(2) + t565 * t526;
t471 = qJD(3) * qJ(4) + t485;
t590 = pkin(3) * t565 + qJ(4) * t561 + pkin(2);
t558 = cos(pkin(10));
t665 = pkin(1) * t558;
t497 = -t590 - t665;
t474 = t497 * qJD(1);
t400 = -t471 * t555 + t557 * t474;
t388 = -pkin(4) * t638 - pkin(8) * t503 + t400;
t401 = t557 * t471 + t555 * t474;
t390 = -pkin(8) * t501 + t401;
t358 = t388 * t560 + t390 * t564;
t354 = -pkin(9) * t441 + t358;
t629 = qJD(6) * t559;
t352 = t354 * t629;
t468 = -qJD(3) * pkin(3) + qJD(4) - t484;
t434 = pkin(4) * t501 + t468;
t387 = pkin(5) * t441 + t434;
t551 = pkin(11) + qJ(5);
t549 = qJ(6) + t551;
t536 = sin(t549);
t537 = cos(t549);
t658 = t545 * t565;
t452 = t536 * t547 - t537 * t658;
t657 = t547 * t565;
t454 = t536 * t545 + t537 * t657;
t664 = g(3) * t561;
t671 = g(1) * t454 - g(2) * t452 - t387 * t678 + t537 * t664 + t352;
t451 = t536 * t658 + t537 * t547;
t453 = -t536 * t657 + t537 * t545;
t542 = t557 * qJDD(3);
t626 = qJD(1) * qJD(3);
t616 = t565 * t626;
t624 = qJDD(1) * t561;
t584 = t616 + t624;
t466 = t555 * t584 - t542;
t623 = qJDD(3) * t555;
t467 = t557 * t584 + t623;
t630 = qJD(5) * t564;
t632 = qJD(5) * t560;
t375 = -t560 * t466 + t564 * t467 - t501 * t630 - t503 * t632;
t548 = t565 * qJDD(1);
t666 = t561 * t626 - t548;
t508 = qJDD(5) + t666;
t412 = qJDD(3) * qJ(4) + qJDD(2) * t561 + t524 * t565 + (qJD(4) + t484) * qJD(3);
t488 = qJD(3) * t599 - qJD(4) * t561;
t425 = qJD(1) * t488 + qJDD(1) * t497;
t370 = -t412 * t555 + t557 * t425;
t364 = pkin(4) * t666 - pkin(8) * t467 + t370;
t371 = t557 * t412 + t555 * t425;
t367 = -pkin(8) * t466 + t371;
t611 = t564 * t364 - t367 * t560;
t571 = -qJD(5) * t358 + t611;
t343 = pkin(5) * t508 - pkin(9) * t375 + t571;
t376 = -qJD(5) * t440 + t564 * t466 + t467 * t560;
t583 = t560 * t364 + t564 * t367 + t388 * t630 - t390 * t632;
t344 = -pkin(9) * t376 + t583;
t612 = t563 * t343 - t559 * t344;
t670 = -g(1) * t453 + g(2) * t451 + t387 * t593 + t536 * t664 + t612;
t496 = qJDD(6) + t508;
t669 = t496 * MDP(27) + (t593 ^ 2 - t678 ^ 2) * MDP(24) + t678 * MDP(23) * t593;
t651 = qJDD(2) - g(3);
t668 = t565 * t651;
t483 = t557 * t497;
t653 = t557 * t561;
t422 = -pkin(8) * t653 + t483 + (-t538 * t555 - pkin(4)) * t565;
t652 = t557 * t565;
t449 = t555 * t497 + t538 * t652;
t655 = t555 * t561;
t431 = -pkin(8) * t655 + t449;
t645 = t560 * t422 + t564 * t431;
t662 = pkin(8) + qJ(4);
t522 = t662 * t555;
t523 = t662 * t557;
t642 = -t560 * t522 + t564 * t523;
t432 = -t484 * t555 + t557 * t515;
t591 = pkin(4) * t561 - pkin(8) * t652;
t397 = qJD(1) * t591 + t432;
t592 = qJD(4) * t555 + qJD(5) * t523;
t667 = -t522 * t630 + t679 * t564 + (-t397 - t592) * t560;
t610 = t559 * t375 + t563 * t376;
t348 = -qJD(6) * t593 + t610;
t357 = t564 * t388 - t390 * t560;
t353 = pkin(9) * t440 + t357;
t351 = -pkin(5) * t534 + t353;
t661 = t351 * t563;
t660 = t354 * t563;
t654 = t555 * t565;
t528 = t561 * t538;
t480 = t510 * t561;
t481 = t509 * t561;
t419 = -t480 * t559 - t481 * t563;
t631 = qJD(5) * t561;
t427 = -qJD(3) * t586 - t510 * t631;
t428 = qJD(3) * t585 + t630 * t653 - t631 * t656;
t360 = qJD(6) * t419 + t427 * t559 + t563 * t428;
t418 = t563 * t480 - t481 * t559;
t650 = t360 * t531 - t418 * t496;
t445 = t563 * t509 + t510 * t559;
t649 = -qJD(6) * t445 - t559 * t643 + t563 * t644;
t446 = -t509 * t559 + t510 * t563;
t648 = qJD(6) * t446 + t559 * t644 + t563 * t643;
t646 = t428 * t534 - t480 * t508;
t635 = qJD(3) * t561;
t619 = t538 * t635;
t438 = t557 * t488 + t555 * t619;
t478 = (pkin(4) * t555 + t538) * t634;
t487 = pkin(4) * t655 + t528;
t553 = t561 ^ 2;
t641 = -t565 ^ 2 + t553;
t539 = -pkin(2) - t665;
t527 = qJD(1) * t539;
t628 = qJD(6) * t563;
t622 = t563 * t375 - t559 * t376 - t441 * t628;
t450 = pkin(4) * t620 + t485;
t540 = -pkin(4) * t557 - pkin(3);
t618 = qJ(4) * t548;
t615 = t555 * t624;
t614 = t557 * t624;
t613 = pkin(5) * t643 - t450;
t408 = qJD(3) * t591 + t438;
t476 = t555 * t488;
t421 = t476 + (-pkin(8) * t654 - t528 * t557) * qJD(3);
t609 = t564 * t408 - t421 * t560;
t608 = t564 * t422 - t431 * t560;
t606 = -t564 * t522 - t523 * t560;
t605 = qJD(6) * t351 + t344;
t603 = g(1) * t545 - g(2) * t547;
t562 = sin(qJ(1));
t566 = cos(qJ(1));
t602 = g(1) * t562 - g(2) * t566;
t396 = t564 * t397;
t424 = -pkin(9) * t509 + t642;
t601 = pkin(5) * t639 + pkin(9) * t644 + t510 * qJD(4) + t642 * qJD(5) + qJD(6) * t424 - t410 * t560 + t396;
t423 = -pkin(9) * t510 + t606;
t600 = -pkin(9) * t643 + qJD(6) * t423 + t667;
t346 = t351 * t559 + t660;
t359 = -qJD(6) * t418 + t427 * t563 - t428 * t559;
t598 = t359 * t531 - t419 * t496;
t597 = -t370 * t555 + t371 * t557;
t596 = -t400 * t555 + t401 * t557;
t594 = t427 * t534 + t481 * t508;
t589 = t602 * pkin(1);
t582 = t560 * t408 + t564 * t421 + t422 * t630 - t431 * t632;
t347 = t440 * t629 + t622;
t580 = -qJD(1) * t527 + t604;
t579 = -qJ(4) * t635 + (qJD(4) - t468) * t565;
t578 = 0.2e1 * qJD(3) * t527 - qJDD(3) * t538;
t575 = -t565 * t604 - t664;
t567 = qJD(3) ^ 2;
t570 = -0.2e1 * qJDD(1) * t539 - t538 * t567 + t603;
t389 = pkin(4) * t466 + t420;
t546 = cos(t551);
t544 = sin(t551);
t518 = qJDD(3) * t565 - t561 * t567;
t517 = qJDD(3) * t561 + t565 * t567;
t475 = pkin(5) * t509 + t540;
t465 = t544 * t545 + t546 * t657;
t464 = -t544 * t657 + t545 * t546;
t463 = t544 * t547 - t546 * t658;
t462 = t544 * t658 + t546 * t547;
t448 = -t538 * t654 + t483;
t439 = -t557 * t619 + t476;
t437 = pkin(5) * t480 + t487;
t429 = t440 * t635;
t391 = pkin(5) * t428 + t478;
t377 = t593 * t635;
t369 = -pkin(9) * t480 + t645;
t368 = -pkin(5) * t565 + pkin(9) * t481 + t608;
t355 = pkin(5) * t376 + t389;
t350 = -pkin(9) * t428 + t582;
t349 = pkin(5) * t635 - pkin(9) * t427 - qJD(5) * t645 + t609;
t345 = -t354 * t559 + t661;
t1 = [(t561 * t578 + t565 * t570) * MDP(10) + (-t561 * t570 + t565 * t578) * MDP(11) + (-t609 * t534 + t608 * t508 - t611 * t565 + t357 * t635 + t478 * t441 + t487 * t376 + t389 * t480 + t434 * t428 - g(1) * t463 - g(2) * t465 + (t358 * t565 + t534 * t645) * qJD(5)) * MDP(21) + (-t347 * t418 - t348 * t419 + t359 * t678 + t360 * t593) * MDP(24) + (-(t349 * t563 - t350 * t559) * t531 + (t368 * t563 - t369 * t559) * t496 - t612 * t565 + t345 * t635 - t391 * t678 + t437 * t348 + t355 * t418 + t387 * t360 - g(1) * t452 - g(2) * t454 + (-(-t368 * t559 - t369 * t563) * t531 + t346 * t565) * qJD(6)) * MDP(28) + (t348 * t565 + t635 * t678 + t650) * MDP(26) + t517 * MDP(7) + t518 * MDP(8) + ((t556 ^ 2 + t558 ^ 2) * pkin(1) ^ 2 * qJDD(1) + t589) * MDP(4) + (g(1) * t566 + g(2) * t562) * MDP(3) + (t347 * t419 - t359 * t593) * MDP(23) + (-t346 * t635 - g(1) * t451 - g(2) * t453 + t437 * t347 - t352 * t565 + t355 * t419 + t387 * t359 - t391 * t593 + ((-qJD(6) * t369 + t349) * t531 - t368 * t496 + t343 * t565) * t559 + ((qJD(6) * t368 + t350) * t531 - t369 * t496 + t605 * t565) * t563) * MDP(29) + (-t375 * t480 + t376 * t481 - t427 * t441 + t428 * t440) * MDP(17) + (-t375 * t481 - t427 * t440) * MDP(16) + (-g(1) * t462 - g(2) * t464 - t358 * t635 + t487 * t375 - t389 * t481 + t434 * t427 - t440 * t478 - t508 * t645 + t534 * t582 + t565 * t583) * MDP(22) + (-t347 * t565 - t377 - t598) * MDP(25) + (-t375 * t565 - t429 - t594) * MDP(18) + qJDD(1) * MDP(1) + t602 * MDP(2) + (-t604 * t555 + (t420 * t555 + t538 * t466 + (qJD(1) * t448 + t400) * qJD(3)) * t561 + (-t438 * qJD(1) - t448 * qJDD(1) - t370 + t603 * t557 + (t468 * t555 + t501 * t538) * qJD(3)) * t565) * MDP(12) + (-t604 * t557 + (t420 * t557 + t538 * t467 + (-qJD(1) * t449 - t401) * qJD(3)) * t561 + (t439 * qJD(1) + t449 * qJDD(1) + t371 - t603 * t555 + (t468 * t557 + t503 * t538) * qJD(3)) * t565) * MDP(13) + (qJDD(1) * t553 + 0.2e1 * t561 * t616) * MDP(5) + (-t438 * t503 - t439 * t501 - t448 * t467 - t449 * t466 + (-t400 * t557 - t401 * t555) * t634 + (-t370 * t557 - t371 * t555 + t603) * t561) * MDP(14) + (t370 * t448 + t371 * t449 + t400 * t438 + t401 * t439 + t589 + (-g(1) * pkin(7) - g(2) * t590) * t547 + (-g(2) * pkin(7) + g(1) * t590) * t545 + (t420 * t561 + t468 * t634) * t538) * MDP(15) + (-t508 * t565 - t534 * t635) * MDP(20) + (-t496 * t565 - t531 * t635) * MDP(27) + 0.2e1 * (t548 * t561 - t626 * t641) * MDP(6) + (t376 * t565 - t441 * t635 + t646) * MDP(19); t651 * MDP(4) + t518 * MDP(10) - t517 * MDP(11) - g(3) * MDP(15) + t646 * MDP(21) + (-t429 + t594) * MDP(22) + t650 * MDP(28) + (-t377 + t598) * MDP(29) + ((-t466 * t557 + t467 * t555) * MDP(14) + t597 * MDP(15)) * t561 + ((-t466 + t615) * MDP(12) + (-t467 + t614) * MDP(13) - t420 * MDP(15) - t376 * MDP(21) - t375 * MDP(22) - t348 * MDP(28) - t347 * MDP(29)) * t565 + ((t501 * MDP(12) + t503 * MDP(13) + t468 * MDP(15) + MDP(21) * t441 - MDP(28) * t678) * t561 + ((-t501 * t557 + t503 * t555) * MDP(14) + t596 * MDP(15)) * t565 + (-MDP(12) * t555 - MDP(13) * t557) * qJD(1) * t641) * qJD(3); qJDD(3) * MDP(9) + (-t400 * t432 - t401 * t433 - t468 * t485 + t596 * qJD(4) - t680 * pkin(3) + (t575 + t597) * qJ(4)) * MDP(15) + MDP(8) * t548 + MDP(7) * t624 + (t557 * t618 - pkin(3) * t467 - t485 * t503 + t680 * t555 + (t401 * t561 - t433 * t565 + t557 * t579) * qJD(1)) * MDP(13) + (t555 * t618 - pkin(3) * t466 - t485 * t501 - t680 * t557 + (-t400 * t561 + t432 * t565 + t555 * t579) * qJD(1)) * MDP(12) + (t432 * t503 + t433 * t501 + (-qJ(4) * t466 - qJD(4) * t501 + t400 * t638 + t371) * t557 + (qJ(4) * t467 + qJD(4) * t503 + t401 * t638 - t370) * t555 + t575) * MDP(14) + (t606 * t508 + t540 * t376 + t389 * t509 - t450 * t441 + (t396 + t592 * t564 + (-qJD(5) * t522 + t679) * t560) * t534 + t643 * t434 + t577 * t546) * MDP(21) + (-t508 * t509 + t534 * t643) * MDP(19) + (t375 * t510 - t440 * t644) * MDP(16) + (t508 * t510 - t534 * t644) * MDP(18) + (-t375 * t509 - t376 * t510 + t440 * t643 - t441 * t644) * MDP(17) + (t540 * t375 + t389 * t510 + t644 * t434 + t440 * t450 - t642 * t508 + t534 * t667 - t544 * t577) * MDP(22) + ((t423 * t563 - t424 * t559) * t496 + t475 * t348 + t355 * t445 + (t559 * t600 + t563 * t601) * t531 + t648 * t387 - t613 * t678 + t577 * t537) * MDP(28) + (-t445 * t496 + t531 * t648) * MDP(26) + (t347 * t446 - t593 * t649) * MDP(23) + (-(t423 * t559 + t424 * t563) * t496 + t475 * t347 + t355 * t446 + (-t559 * t601 + t563 * t600) * t531 + t649 * t387 - t613 * t593 - t577 * t536) * MDP(29) + (t446 * t496 - t531 * t649) * MDP(25) + (-t347 * t445 - t348 * t446 + t593 * t648 + t649 * t678) * MDP(24) + (qJD(3) * t484 + (qJD(3) * t526 - t651) * t561 + (t580 + t672) * t565) * MDP(11) + (qJD(3) * t485 + t561 * t580 + t621 + t668) * MDP(10) + (MDP(18) * t440 + t441 * MDP(19) + t534 * MDP(20) - t357 * MDP(21) + t358 * MDP(22) + MDP(25) * t593 - MDP(26) * t678 + t531 * MDP(27) - t345 * MDP(28) + t346 * MDP(29)) * t639 + (-MDP(5) * t561 * t565 + MDP(6) * t641) * qJD(1) ^ 2; (t615 - t542 + (-t503 + t636) * t638) * MDP(12) + (t614 + t623 + (t501 + t627) * t638) * MDP(13) + (-t501 ^ 2 - t503 ^ 2) * MDP(14) + (t400 * t503 + t401 * t501 + t581 - t588 - t668) * MDP(15) + (t376 + t676) * MDP(21) + (t375 + t675) * MDP(22) + (t348 + t674) * MDP(28) + (t347 - t677) * MDP(29); -t440 * t441 * MDP(16) + (t440 ^ 2 - t441 ^ 2) * MDP(17) + (t375 - t675) * MDP(18) + (-t376 + t676) * MDP(19) + t508 * MDP(20) + (-g(1) * t464 + g(2) * t462 - t358 * t534 + t434 * t440 + t544 * t664 + t571) * MDP(21) + (g(1) * t465 - g(2) * t463 - t357 * t534 + t434 * t441 + t546 * t664 - t583) * MDP(22) + (t347 + t677) * MDP(25) + (-t348 + t674) * MDP(26) + ((-t353 * t559 - t660) * t531 - t346 * qJD(6) + (-t440 * t678 + t496 * t563 + t531 * t629) * pkin(5) + t670) * MDP(28) + ((t354 * t531 - t343) * t559 + (-t353 * t531 - t605) * t563 + (-t440 * t593 - t496 * t559 + t531 * t628) * pkin(5) + t671) * MDP(29) + t669; (t622 + t677) * MDP(25) + (-t610 + t674) * MDP(26) + (-t346 * t531 + t670) * MDP(28) + (-t559 * t343 - t563 * t344 - t345 * t531 + t671) * MDP(29) + (MDP(25) * t659 + MDP(26) * t593 - MDP(28) * t346 - MDP(29) * t661) * qJD(6) + t669;];
tau  = t1;
