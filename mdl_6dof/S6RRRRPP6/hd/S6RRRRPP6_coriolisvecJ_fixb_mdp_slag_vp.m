% Calculate minimal parameter regressor of Coriolis joint torque vector for
% S6RRRRPP6
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d3,d4]';
% MDP [32x1]
%   Minimal dynamic parameter vector (fixed base model)
%   see S6RRRRPP6_convert_par2_MPV_fixb.m
% 
% Output:
% tauc [6x1]
%   joint torques required to compensate Coriolis and centrifugal load

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 21:16
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauc = S6RRRRPP6_coriolisvecJ_fixb_mdp_slag_vp(qJ, qJD, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(9,1),zeros(32,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRRPP6_coriolisvecJ_fixb_mdp_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRRRPP6_coriolisvecJ_fixb_mdp_slag_vp: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S6RRRRPP6_coriolisvecJ_fixb_mdp_slag_vp: pkin has to be [9x1] (double)');
assert(isreal(MDP) && all(size(MDP) == [32 1]), ...
  'S6RRRRPP6_coriolisvecJ_fixb_mdp_slag_vp: MDP has to be [32x1] (double)'); 

%% Symbolic Calculation
% From coriolisvec_joint_fixb_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 21:15:21
% EndTime: 2019-03-09 21:15:34
% DurationCPUTime: 8.21s
% Computational Cost: add. (7748->587), mult. (18478->704), div. (0->0), fcn. (12390->6), ass. (0->234)
t548 = cos(qJ(3));
t547 = sin(qJ(2));
t549 = cos(qJ(2));
t653 = t548 * t549;
t574 = pkin(3) * t547 - pkin(9) * t653;
t675 = -pkin(9) - pkin(8);
t607 = qJD(3) * t675;
t581 = pkin(2) * t547 - pkin(8) * t549;
t487 = t581 * qJD(1);
t546 = sin(qJ(3));
t628 = qJD(1) * t547;
t600 = t546 * t628;
t632 = pkin(7) * t600 + t548 * t487;
t701 = -qJD(1) * t574 + t548 * t607 - t632;
t471 = t546 * t487;
t655 = t547 * t548;
t700 = t471 + (-pkin(9) * t546 * t549 - pkin(7) * t655) * qJD(1) - t546 * t607;
t622 = qJD(2) * t548;
t483 = -t600 + t622;
t624 = qJD(2) * t546;
t484 = t548 * t628 + t624;
t545 = sin(qJ(4));
t674 = cos(qJ(4));
t437 = -t674 * t483 + t484 * t545;
t571 = t545 * t483 + t484 * t674;
t662 = t437 * t571;
t627 = qJD(1) * t549;
t605 = t546 * t627;
t606 = t674 * t548;
t658 = t545 * t546;
t678 = qJD(3) + qJD(4);
t597 = t674 * qJD(4);
t681 = t674 * qJD(3) + t597;
t638 = -t545 * t605 - t548 * t681 + t606 * t627 + t658 * t678;
t619 = qJD(3) * t547;
t603 = t546 * t619;
t621 = qJD(2) * t549;
t563 = t548 * t621 - t603;
t611 = qJD(2) * qJD(3);
t555 = qJD(1) * t563 + t548 * t611;
t595 = qJD(1) * t619;
t612 = qJD(1) * qJD(2);
t596 = t549 * t612;
t608 = t548 * t595 + (t596 + t611) * t546;
t617 = qJD(4) * t545;
t387 = -t483 * t597 + t484 * t617 + t545 * t608 - t674 * t555;
t527 = -qJD(3) + t627;
t508 = -qJD(4) + t527;
t366 = -t437 * t508 - t387;
t435 = t571 ^ 2;
t531 = t547 * t612;
t388 = qJD(4) * t571 + t545 * t555 + t674 * t608;
t552 = -t508 * t571 - t388;
t676 = t437 ^ 2;
t699 = MDP(18) * t662 + t552 * MDP(21) + MDP(22) * t531 + t366 * MDP(20) + (t435 - t676) * MDP(19);
t665 = qJ(5) * t437;
t667 = qJD(2) * pkin(2);
t500 = pkin(7) * t628 - t667;
t457 = -pkin(3) * t483 + t500;
t562 = -qJ(5) * t571 + t457;
t668 = pkin(4) + qJ(6);
t368 = t437 * t668 + t562;
t698 = t368 * t437;
t383 = pkin(4) * t437 + t562;
t697 = t383 * t437;
t502 = t675 * t546;
t503 = t675 * t548;
t696 = t502 * t597 + t503 * t617 + t545 * t701 - t700 * t674;
t486 = t545 * t548 + t546 * t674;
t450 = t678 * t486;
t637 = -t486 * t627 + t450;
t537 = pkin(7) * t627;
t620 = qJD(3) * t546;
t583 = -t537 + (-t605 + t620) * pkin(3);
t604 = t546 * t621;
t618 = qJD(3) * t548;
t695 = t547 * t618 + t604;
t494 = -pkin(2) * t549 - pkin(8) * t547 - pkin(1);
t476 = t494 * qJD(1);
t501 = qJD(2) * pkin(8) + t537;
t445 = t548 * t476 - t501 * t546;
t422 = -pkin(9) * t484 + t445;
t416 = -pkin(3) * t527 + t422;
t657 = t546 * t476;
t446 = t548 * t501 + t657;
t423 = pkin(9) * t483 + t446;
t659 = t545 * t423;
t373 = -t674 * t416 + t659;
t689 = pkin(5) * t571;
t572 = t373 + t689;
t614 = qJD(5) + t572;
t679 = pkin(5) * t437 - qJD(6);
t564 = t574 * qJD(2);
t490 = t581 * qJD(2);
t477 = qJD(1) * t490;
t587 = pkin(7) * t531;
t635 = t548 * t477 + t546 * t587;
t382 = qJD(1) * t564 - qJD(3) * t423 + t635;
t636 = t476 * t618 + t546 * t477;
t556 = -t501 * t620 - t548 * t587 + t636;
t392 = -pkin(9) * t608 + t556;
t586 = t545 * t382 + t674 * t392 + t416 * t597 - t423 * t617;
t693 = t437 * t457 - t586;
t691 = -0.2e1 * t612;
t690 = pkin(4) * t571;
t542 = t547 ^ 2;
t688 = MDP(5) * (-t549 ^ 2 + t542);
t687 = t383 * t571;
t686 = t571 * t668;
t420 = t674 * t423;
t377 = t545 * t422 + t420;
t582 = pkin(3) * t617 - t377;
t378 = t422 * t674 - t659;
t644 = -pkin(3) * t597 - qJD(5) + t378;
t642 = qJ(5) * t628 - t696;
t493 = qJD(5) * t508;
t520 = qJ(5) * t531;
t684 = t520 - t493;
t683 = qJ(5) * t638 - qJD(5) * t486 + t583;
t456 = t545 * t502 - t503 * t674;
t682 = -qJD(4) * t456 + t700 * t545 + t674 * t701;
t616 = t542 * qJD(1);
t525 = pkin(4) * t531;
t585 = -t674 * t382 + t545 * t392 + t416 * t617 + t423 * t597;
t352 = -t525 + t585;
t348 = -t387 * pkin(5) - qJ(6) * t531 + qJD(6) * t508 + t352;
t559 = -t368 * t571 - t348;
t677 = -t457 * t571 - t585;
t504 = t508 ^ 2;
t673 = pkin(5) * t388;
t671 = pkin(7) * t546;
t670 = pkin(8) * t527;
t666 = qJ(5) * t388;
t374 = t545 * t416 + t420;
t664 = t374 * t508;
t530 = pkin(3) * t545 + qJ(5);
t663 = t388 * t530;
t661 = t484 * t527;
t660 = t500 * t546;
t656 = t546 * t547;
t550 = qJD(2) ^ 2;
t654 = t547 * t550;
t652 = t549 * t527;
t651 = t549 * t550;
t551 = qJD(1) ^ 2;
t650 = t549 * t551;
t649 = t582 + t679;
t648 = t644 - t689;
t485 = -t606 + t658;
t647 = qJD(6) * t485 + t637 * t668 + t683;
t599 = t547 * t668;
t646 = -pkin(5) * t638 + qJD(1) * t599 - t682;
t645 = pkin(5) * t637 + t642;
t643 = pkin(4) * t637 + t683;
t641 = -pkin(4) * t628 + t682;
t482 = t548 * t494;
t444 = -pkin(9) * t655 + t482 + (-pkin(3) - t671) * t549;
t529 = pkin(7) * t653;
t631 = t546 * t494 + t529;
t451 = -pkin(9) * t656 + t631;
t639 = t545 * t444 + t674 * t451;
t634 = t546 * t490 + t494 * t618;
t623 = qJD(2) * t547;
t633 = t548 * t490 + t623 * t671;
t491 = pkin(3) * t656 + t547 * pkin(7);
t455 = -t502 * t674 - t545 * t503;
t626 = qJD(2) * t455;
t625 = qJD(2) * t456;
t615 = -qJD(5) - t373;
t613 = -t374 + t679;
t609 = t545 * t656;
t458 = pkin(3) * t695 + pkin(7) * t621;
t535 = -pkin(3) * t548 - pkin(2);
t601 = t527 * t618;
t594 = MDP(15) * t623;
t593 = pkin(1) * t691;
t591 = t527 + t627;
t590 = -t483 + t622;
t589 = -t484 + t624;
t588 = qJD(3) + t627;
t534 = -pkin(3) * t674 - pkin(4);
t584 = t674 * t621;
t371 = qJ(5) * t508 - t374;
t399 = qJ(5) * t549 - t639;
t466 = t547 * t606 - t609;
t579 = -qJ(5) * t466 + t491;
t577 = t444 * t674 - t545 * t451;
t576 = pkin(3) * t484 + t665;
t575 = t493 - t586;
t573 = -qJ(5) * t486 + t535;
t400 = t549 * pkin(4) - t577;
t434 = pkin(3) * t608 + pkin(7) * t596;
t351 = -t520 + t575;
t569 = t586 - t673;
t568 = t588 * t624;
t567 = t585 + t687;
t403 = t564 + (-t529 + (pkin(9) * t547 - t494) * t546) * qJD(3) + t633;
t406 = -t695 * pkin(9) + (-t547 * t622 - t549 * t620) * pkin(7) + t634;
t566 = -t403 * t674 + t545 * t406 + t444 * t617 + t451 * t597;
t565 = t545 * t403 + t674 * t406 + t444 * t597 - t451 * t617;
t407 = t450 * t547 + t545 * t604 - t548 * t584;
t561 = qJ(5) * t407 - qJD(5) * t466 + t458;
t560 = t569 - t698;
t558 = t530 * t531 - t351;
t557 = t387 * qJ(5) - qJD(5) * t571 + t434;
t357 = -qJ(5) * t623 + qJD(5) * t549 - t565;
t524 = -qJ(6) + t534;
t515 = 0.2e1 * t520;
t465 = t486 * t547;
t433 = pkin(4) * t485 + t573;
t426 = -t485 * pkin(5) + t456;
t425 = t486 * pkin(5) + t455;
t421 = t485 * t668 + t573;
t417 = pkin(4) * t465 + t579;
t408 = t546 * t584 - t545 * t603 - qJD(4) * t609 + (t545 * t621 + t547 * t681) * t548;
t402 = t665 + t690;
t401 = t465 * t668 + t579;
t394 = t576 + t690;
t389 = -pkin(5) * t465 - t399;
t381 = t466 * pkin(5) + t549 * qJ(6) + t400;
t376 = t665 + t686;
t370 = pkin(4) * t508 - t615;
t369 = t576 + t686;
t361 = -t371 - t679;
t360 = t508 * t668 + t614;
t359 = pkin(4) * t408 + t561;
t358 = -pkin(4) * t623 + t566;
t356 = qJD(6) * t465 + t408 * t668 + t561;
t355 = t388 * pkin(4) + t557;
t354 = -pkin(5) * t408 - t357;
t353 = -t407 * pkin(5) - qJD(2) * t599 + t549 * qJD(6) + t566;
t350 = t437 * qJD(6) + t388 * t668 + t557;
t349 = t569 + t684;
t1 = [(t547 * t601 + t608 * t549 + (t483 * t547 + (-t616 + t652) * t546) * qJD(2)) * MDP(14) - MDP(7) * t654 + (pkin(7) * t654 + t549 * t593) * MDP(10) + (t484 * t563 + t555 * t655) * MDP(11) + (t591 * t603 + (t484 * t547 + (t616 + (-t527 - t588) * t549) * t548) * qJD(2)) * MDP(13) + ((t483 * t548 - t484 * t546) * t621 + ((-t483 + t600) * t620 + (-t484 * qJD(3) - t568 - t608) * t548) * t547) * MDP(12) + (-(-t494 * t620 + t633) * t527 + (pkin(7) * t608 + t500 * t618 + (qJD(1) * t482 + t445) * qJD(2)) * t547 + ((-pkin(7) * t483 + t660) * qJD(2) + (t657 + (pkin(7) * t527 + t501) * t548) * qJD(3) - t635) * t549) * MDP(16) + (-t508 - t627) * MDP(22) * t623 + (t388 * t549 + t408 * t508 + (-qJD(1) * t465 - t437) * t623) * MDP(21) + (t566 * t508 + t585 * t549 + t458 * t437 + t491 * t388 + t434 * t465 + t457 * t408 + (qJD(1) * t577 - t373) * t623) * MDP(23) + (-t352 * t549 - t355 * t465 - t358 * t508 - t359 * t437 - t383 * t408 - t388 * t417 + (qJD(1) * t400 + t370) * t623) * MDP(26) + (t348 * t549 + t350 * t465 + t353 * t508 + t356 * t437 + t368 * t408 + t388 * t401 + (-qJD(1) * t381 - t360) * t623) * MDP(31) + 0.2e1 * t549 * MDP(4) * t531 - t591 * t594 + (t565 * t508 + t586 * t549 + t458 * t571 - t491 * t387 + t434 * t466 - t457 * t407 + (-qJD(1) * t639 - t374) * t623) * MDP(24) + (t387 * t549 + t407 * t508 + (qJD(1) * t466 + t571) * t623) * MDP(20) + (t351 * t549 - t355 * t466 + t357 * t508 - t359 * t571 + t383 * t407 + t387 * t417 + (-qJD(1) * t399 - t371) * t623) * MDP(27) + (-t349 * t549 - t350 * t466 - t354 * t508 - t356 * t571 + t368 * t407 + t387 * t401 + (qJD(1) * t389 + t361) * t623) * MDP(30) + (t351 * t465 + t352 * t466 + t357 * t437 + t358 * t571 - t370 * t407 + t371 * t408 - t387 * t400 + t388 * t399) * MDP(25) + (t387 * t465 - t388 * t466 + t407 * t437 - t408 * t571) * MDP(19) + (t348 * t466 - t349 * t465 + t353 * t571 - t354 * t437 - t360 * t407 - t361 * t408 - t381 * t387 - t388 * t389) * MDP(29) + (-t387 * t466 - t407 * t571) * MDP(18) + (t351 * t399 + t352 * t400 + t355 * t417 + t357 * t371 + t358 * t370 + t359 * t383) * MDP(28) + (t348 * t381 + t349 * t389 + t350 * t401 + t353 * t360 + t354 * t361 + t356 * t368) * MDP(32) + (-pkin(7) * t651 + t547 * t593) * MDP(9) + (t634 * t527 + t636 * t549 + (-t500 * t547 - t501 * t549 + (-t652 - t616) * pkin(7)) * t620 + ((pkin(7) * t484 + t500 * t548) * t549 + (-t631 * qJD(1) - t446 + (-t527 + t588) * pkin(7) * t548) * t547) * qJD(2)) * MDP(17) + t688 * t691 + MDP(6) * t651; (-t387 * t486 - t571 * t638) * MDP(18) + (-t535 * t387 + t434 * t486 - t638 * t457 + t571 * t583) * MDP(24) + (t351 * t485 + t352 * t486 - t370 * t638 + t371 * t637 - t387 * t455 - t388 * t456 + t437 * t642 - t571 * t641) * MDP(25) + (-t351 * t456 + t352 * t455 + t355 * t433 - t370 * t641 + t371 * t642 + t383 * t643) * MDP(28) + (-t355 * t485 - t637 * t383 - t388 * t433 - t643 * t437) * MDP(26) + (-t355 * t486 + t638 * t383 + t387 * t433 - t571 * t643) * MDP(27) + (-t546 ^ 2 * t595 + (t568 - t661) * t548) * MDP(11) + ((-t608 + t661) * t546 + ((t483 + t622) * qJD(3) + (t549 * t590 - t603) * qJD(1)) * t548) * MDP(12) + t551 * t688 + (t535 * t388 + t434 * t485 + t583 * t437 + t637 * t457) * MDP(23) + (t387 * t485 - t388 * t486 + t437 * t638 - t571 * t637) * MDP(19) + (-t471 * t527 + (-t546 * t670 + (t500 - t667) * t548) * qJD(3) + ((-t500 - t667) * t653 + (pkin(2) * t620 - pkin(8) * t622 + t446) * t547 + (t527 * t655 + t549 * t589) * pkin(7)) * qJD(1)) * MDP(17) + (-pkin(2) * t608 + t632 * t527 + (t548 * t670 + t660) * qJD(3) + ((-pkin(8) * t624 - t445) * t547 + (-pkin(7) * t590 - t660) * t549) * qJD(1)) * MDP(16) - t547 * MDP(4) * t650 + (t348 * t486 - t349 * t485 - t360 * t638 - t361 * t637 - t387 * t425 - t388 * t426 + t437 * t645 + t571 * t646) * MDP(29) + (t348 * t425 + t349 * t426 + t350 * t421 + t360 * t646 - t361 * t645 + t368 * t647) * MDP(32) + (t350 * t485 + t637 * t368 + t388 * t421 + t647 * t437) * MDP(31) + (-t350 * t486 + t638 * t368 + t387 * t421 - t571 * t647) * MDP(30) + (t527 * t620 + (-t546 * t652 + t547 * t590) * qJD(1)) * MDP(14) + (-t601 + (t547 * t589 + t548 * t652) * qJD(1)) * MDP(13) + (t638 * MDP(20) + t637 * MDP(21) - t682 * MDP(23) + MDP(24) * t696 + t641 * MDP(26) + t642 * MDP(27) + t645 * MDP(30) + t646 * MDP(31)) * t508 + ((t374 - t625) * MDP(24) + (-t370 + t626) * MDP(26) + (t371 + t625) * MDP(27) + (-qJD(2) * t485 + t437) * MDP(21) + (t373 - t626) * MDP(23) + (qJD(2) * t486 - t571) * MDP(20) + (-qJD(2) * t425 + t360) * MDP(31) + (qJD(2) * t426 - t361) * MDP(30) + t527 * MDP(15) + t508 * MDP(22)) * t628 + (MDP(9) * t547 * t551 + MDP(10) * t650) * pkin(1); (-t445 * t527 - t483 * t500 - t556) * MDP(17) + (-t351 * t530 + t352 * t534 + t370 * t582 + t371 * t644 - t383 * t394) * MDP(28) + (-t608 - t661) * MDP(14) + (t483 * t527 + t555) * MDP(13) - t484 * t483 * MDP(11) + (-t484 * t500 + t635 + (-qJD(3) - t527) * t446) * MDP(16) + (-t483 ^ 2 + t484 ^ 2) * MDP(12) + (t394 * t437 - t508 * t582 + t531 * t534 + t352 + t687) * MDP(26) + (-t369 * t437 + t508 * t649 - t524 * t531 + t559) * MDP(31) + (-t377 * t508 + (-t437 * t484 + t508 * t617 + t531 * t674) * pkin(3) + t677) * MDP(23) + (-t378 * t508 + (-t484 * t571 + t508 * t597 - t531 * t545) * pkin(3) + t693) * MDP(24) + (-t387 * t534 - t663 + (-t371 + t582) * t571 + (t370 + t644) * t437) * MDP(25) + (t394 * t571 + t508 * t644 + t558 - t697) * MDP(27) + (-t387 * t524 - t663 + (t361 + t649) * t571 + (t360 + t648) * t437) * MDP(29) + (t369 * t571 + t508 * t648 + t558 - t673 - t698) * MDP(30) + (t348 * t524 + t349 * t530 + t360 * t649 - t361 * t648 - t368 * t369) * MDP(32) + qJD(1) * t594 + t699; (-t664 + t677) * MDP(23) + (t373 * t508 + t693) * MDP(24) + (pkin(4) * t387 - t666 + (-t371 - t374) * t571 + (t370 + t615) * t437) * MDP(25) + (t402 * t437 - 0.2e1 * t525 + t567 + t664) * MDP(26) + (t402 * t571 + t508 * t615 + t515 - t575 - t697) * MDP(27) + (-pkin(4) * t352 - qJ(5) * t351 - t370 * t374 + t371 * t615 - t383 * t402) * MDP(28) + (-t666 + t387 * t668 + (t361 + t613) * t571 + (t360 - t614) * t437) * MDP(29) + (t376 * t571 - t508 * t572 - 0.2e1 * t493 + t515 + t560) * MDP(30) + (-t376 * t437 + t508 * t613 + t531 * t668 + t559) * MDP(31) + (qJ(5) * t349 - t348 * t668 + t360 * t613 + t361 * t614 - t368 * t376) * MDP(32) + t699; (-t371 * t508 - t525 + t567) * MDP(28) + (t361 * t508 - t559) * MDP(32) + (-MDP(26) + MDP(31)) * (-t531 + t662) + (MDP(25) + MDP(29)) * t366 + (MDP(27) + MDP(30)) * (-t435 - t504); t552 * MDP(29) + (t531 + t662) * MDP(30) + (-t504 - t676) * MDP(31) + (-t360 * t508 + t560 + t684) * MDP(32);];
tauc  = t1;
