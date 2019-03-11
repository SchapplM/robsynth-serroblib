% Calculate minimal parameter regressor of Coriolis joint torque vector for
% S6RRPRRR7
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d4,d5,d6]';
% MDP [35x1]
%   Minimal dynamic parameter vector (fixed base model)
%   see S6RRPRRR7_convert_par2_MPV_fixb.m
% 
% Output:
% tauc [6x1]
%   joint torques required to compensate Coriolis and centrifugal load

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 14:01
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauc = S6RRPRRR7_coriolisvecJ_fixb_mdp_slag_vp(qJ, qJD, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(10,1),zeros(35,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPRRR7_coriolisvecJ_fixb_mdp_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRPRRR7_coriolisvecJ_fixb_mdp_slag_vp: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RRPRRR7_coriolisvecJ_fixb_mdp_slag_vp: pkin has to be [10x1] (double)');
assert(isreal(MDP) && all(size(MDP) == [35 1]), ...
  'S6RRPRRR7_coriolisvecJ_fixb_mdp_slag_vp: MDP has to be [35x1] (double)'); 

%% Symbolic Calculation
% From coriolisvec_joint_fixb_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 13:59:58
% EndTime: 2019-03-09 14:00:10
% DurationCPUTime: 7.81s
% Computational Cost: add. (5459->509), mult. (12236->675), div. (0->0), fcn. (8676->8), ass. (0->230)
t545 = sin(qJ(4));
t549 = cos(qJ(4));
t550 = cos(qJ(2));
t628 = qJD(1) * t550;
t546 = sin(qJ(2));
t629 = qJD(1) * t546;
t482 = -t545 * t628 + t549 * t629;
t537 = qJD(2) - qJD(4);
t544 = sin(qJ(5));
t548 = cos(qJ(5));
t449 = t482 * t548 - t537 * t544;
t488 = t545 * t546 + t549 * t550;
t563 = t488 * qJD(4);
t615 = qJD(1) * qJD(2);
t600 = t550 * t615;
t601 = t546 * t615;
t635 = t545 * t601 + t549 * t600;
t432 = -qJD(1) * t563 + t635;
t620 = qJD(5) * t548;
t621 = qJD(5) * t544;
t391 = t548 * t432 - t482 * t621 - t537 * t620;
t665 = t432 * t544;
t392 = qJD(5) * t449 + t665;
t447 = t482 * t544 + t548 * t537;
t543 = sin(qJ(6));
t547 = cos(qJ(6));
t618 = qJD(6) * t547;
t609 = t547 * t391 - t543 * t392 - t447 * t618;
t619 = qJD(6) * t543;
t360 = -t449 * t619 + t609;
t571 = t447 * t543 - t547 * t449;
t595 = t391 * t543 + t547 * t392;
t361 = -qJD(6) * t571 + t595;
t661 = t449 * t543;
t396 = t547 * t447 + t661;
t622 = qJD(4) * t549;
t623 = qJD(4) * t545;
t625 = qJD(2) * t550;
t709 = t545 * t625 + t546 * t622 - t550 * t623;
t433 = qJD(1) * t709 - t549 * t601;
t677 = t488 * qJD(1);
t695 = qJD(5) + t677;
t471 = qJD(6) + t695;
t654 = t543 * t544;
t487 = -t547 * t548 + t654;
t653 = t543 * t548;
t490 = t544 * t547 + t653;
t614 = qJD(5) + qJD(6);
t639 = (t614 + t677) * t490;
t679 = -qJD(6) * t548 - t620;
t560 = t679 * t547;
t641 = t487 * t677 + t614 * t654 + t560;
t649 = t548 * t433;
t652 = t544 * t433;
t663 = t433 * t490;
t664 = t433 * t487;
t666 = t391 * t544;
t683 = t695 * t449;
t684 = t695 ^ 2;
t710 = t447 * t695;
t711 = -(t544 * (t392 + t683) + (-t391 + t710) * t548) * MDP(23) + (t360 * t490 + t571 * t641) * MDP(29) + (-t360 * t487 - t361 * t490 + t396 * t641 + t571 * t639) * MDP(30) - (t471 * t641 - t482 * t571 - t663) * MDP(31) - (-t396 * t482 + t471 * t639 + t664) * MDP(32) + (t548 * t683 + t666) * MDP(22) + (-t449 * t482 + t548 * t684 + t652) * MDP(24) - (-t447 * t482 + t544 * t684 - t649) * MDP(25) - (t482 * t537 + t433) * MDP(18) - (-t482 ^ 2 + t677 ^ 2) * MDP(16) + (MDP(15) * t677 - MDP(26) * t695 - MDP(33) * t471) * t482;
t528 = pkin(7) * t629;
t495 = pkin(8) * t629 - t528;
t708 = qJD(3) - t495;
t484 = -qJD(1) * pkin(1) - pkin(2) * t628 - qJ(3) * t629;
t460 = pkin(3) * t628 - t484;
t406 = pkin(4) * t677 - pkin(9) * t482 + t460;
t551 = -pkin(2) - pkin(3);
t607 = t551 * qJD(2);
t466 = t607 + t708;
t529 = pkin(7) * t628;
t497 = -pkin(8) * t628 + t529;
t539 = qJD(2) * qJ(3);
t483 = t497 + t539;
t425 = t545 * t466 + t549 * t483;
t412 = -pkin(9) * t537 + t425;
t373 = t548 * t406 - t412 * t544;
t367 = -pkin(10) * t449 + t373;
t362 = pkin(5) * t695 + t367;
t374 = t406 * t544 + t412 * t548;
t368 = -pkin(10) * t447 + t374;
t667 = t368 * t547;
t358 = t362 * t543 + t667;
t627 = qJD(2) * t546;
t674 = pkin(7) - pkin(8);
t496 = t674 * t627;
t538 = qJD(2) * qJD(3);
t472 = -qJD(1) * t496 + t538;
t520 = pkin(7) * t600;
t485 = -pkin(8) * t600 + t520;
t386 = t466 * t623 + t545 * t472 + t483 * t622 - t485 * t549;
t363 = pkin(5) * t392 + t386;
t424 = t466 * t549 - t545 * t483;
t411 = pkin(4) * t537 - t424;
t387 = pkin(5) * t447 + t411;
t707 = -t358 * t482 - t363 * t490 + t641 * t387;
t668 = t362 * t547;
t357 = -t368 * t543 + t668;
t706 = t357 * t482 - t363 * t487 - t639 * t387;
t699 = t396 * t471;
t698 = t471 * t571;
t437 = t495 * t549 + t497 * t545;
t678 = -t545 * qJ(3) + t549 * t551;
t467 = qJD(3) * t549 + qJD(4) * t678;
t643 = t437 - t467;
t366 = t368 * t619;
t694 = t387 * t396 + t366;
t587 = t546 * t607;
t532 = t546 * qJD(3);
t633 = qJ(3) * t600 + qJD(1) * t532;
t446 = qJD(1) * t587 + t633;
t379 = pkin(4) * t433 - pkin(9) * t432 + t446;
t377 = t548 * t379;
t385 = t466 * t622 + t549 * t472 - t483 * t623 + t545 * t485;
t556 = -qJD(5) * t374 - t385 * t544 + t377;
t354 = pkin(5) * t433 - pkin(10) * t391 + t556;
t565 = t544 * t379 + t548 * t385 + t406 * t620 - t412 * t621;
t355 = -pkin(10) * t392 + t565;
t596 = t547 * t354 - t543 * t355;
t693 = t387 * t571 + t596;
t688 = t433 * MDP(33) + (-t396 ^ 2 + t571 ^ 2) * MDP(30) - t396 * MDP(29) * t571;
t687 = -0.2e1 * t615;
t540 = t546 ^ 2;
t686 = MDP(5) * (-t550 ^ 2 + t540);
t685 = t411 * t695;
t489 = -t545 * t550 + t546 * t549;
t429 = t490 * t489;
t581 = qJ(3) * t549 + t545 * t551;
t636 = qJD(4) * t581 + t497 * t549 + t545 * t708;
t506 = t674 * t546;
t508 = t674 * t550;
t450 = -t506 * t549 + t545 * t508;
t681 = -t543 * t621 - t544 * t619;
t659 = t677 * t544;
t680 = (t621 + t659) * pkin(5);
t673 = pkin(9) + pkin(10);
t672 = pkin(5) * t548;
t494 = -pkin(9) + t581;
t671 = pkin(10) - t494;
t670 = qJD(2) * pkin(2);
t442 = qJD(2) * t488 - t563;
t662 = t442 * t548;
t660 = t677 * t537;
t657 = t489 * t544;
t656 = t489 * t548;
t651 = t544 * t442;
t552 = qJD(2) ^ 2;
t650 = t546 * t552;
t451 = t506 * t545 + t508 * t549;
t445 = t548 * t451;
t647 = t550 * t552;
t553 = qJD(1) ^ 2;
t646 = t550 * t553;
t434 = pkin(4) * t482 + pkin(9) * t677;
t645 = t548 * t424 + t544 * t434;
t524 = qJ(3) * t628;
t474 = t551 * t629 + t524;
t407 = -t434 + t474;
t644 = t544 * t407 + t548 * t437;
t503 = -t550 * pkin(2) - t546 * qJ(3) - pkin(1);
t486 = t550 * pkin(3) - t503;
t422 = pkin(4) * t488 - pkin(9) * t489 + t486;
t638 = t544 * t422 + t445;
t637 = -t680 + t636;
t632 = qJ(3) * t625 + t532;
t626 = qJD(2) * t549;
t624 = qJD(4) * t471;
t612 = pkin(10) * t659;
t610 = t546 * t646;
t608 = qJD(5) * t673;
t606 = t489 * t620;
t599 = qJD(5) * t671;
t598 = pkin(1) * t687;
t597 = qJD(3) - t670;
t590 = qJD(1) * t503 + t484;
t589 = qJD(6) * t362 + t355;
t588 = t537 ^ 2;
t493 = pkin(4) - t678;
t586 = -t425 + t680;
t405 = t548 * t407;
t462 = t671 * t548;
t570 = -pkin(10) * t548 * t677 - pkin(5) * t482;
t585 = -qJD(6) * t462 - t544 * t643 - t548 * t599 + t405 + t570;
t428 = t548 * t434;
t507 = t673 * t548;
t584 = qJD(6) * t507 - t424 * t544 + t548 * t608 + t428 - t570;
t461 = t671 * t544;
t583 = -qJD(6) * t461 - t467 * t548 - t544 * t599 - t612 + t644;
t505 = t673 * t544;
t582 = qJD(6) * t505 + t544 * t608 + t612 + t645;
t577 = -t373 * t482 - t386 * t548;
t576 = t374 * t482 + t386 * t544;
t459 = pkin(2) * t601 - t633;
t475 = pkin(2) * t627 - t632;
t569 = -pkin(7) * t552 - qJD(1) * t475 - t459;
t568 = t606 + t651;
t567 = -t489 * t621 + t662;
t566 = -pkin(9) * t433 + t685;
t455 = t587 + t632;
t441 = -t546 * t626 + t709;
t384 = pkin(4) * t441 - pkin(9) * t442 + t455;
t498 = qJD(2) * t508;
t401 = -qJD(4) * t450 - t496 * t549 + t498 * t545;
t564 = t544 * t384 + t548 * t401 + t422 * t620 - t451 * t621;
t561 = -t494 * t433 - t685;
t558 = t460 * t482 + t386;
t557 = -t460 * t677 + t385;
t402 = qJD(4) * t451 - t496 * t545 - t498 * t549;
t501 = -pkin(7) * t601 + t538;
t502 = t528 + t597;
t504 = t529 + t539;
t554 = t501 * t550 + (t502 * t550 + (-t504 + t529) * t546) * qJD(2);
t527 = -pkin(4) - t672;
t492 = pkin(2) * t629 - t524;
t481 = t544 * t629 + t548 * t626;
t478 = -t544 * t626 + t548 * t629;
t477 = t493 + t672;
t430 = t487 * t489;
t423 = pkin(5) * t657 + t450;
t414 = t548 * t422;
t409 = t433 * t488;
t382 = t548 * t384;
t380 = -pkin(10) * t657 + t638;
t375 = pkin(5) * t488 - pkin(10) * t656 - t451 * t544 + t414;
t371 = pkin(5) * t568 + t402;
t365 = t442 * t653 + (t614 * t656 + t651) * t547 + t681 * t489;
t364 = -t429 * t614 - t487 * t442;
t359 = -pkin(10) * t568 + t564;
t356 = -pkin(10) * t662 + pkin(5) * t441 - t401 * t544 + t382 + (-t445 + (pkin(10) * t489 - t422) * t544) * qJD(5);
t1 = [(-t360 * t430 - t364 * t571) * MDP(29) + (-t358 * t441 + t423 * t360 - t363 * t430 + t387 * t364 + t366 * t488 - t371 * t571 + (-(-qJD(6) * t380 + t356) * t471 - t375 * t433 - t354 * t488) * t543 + (-(qJD(6) * t375 + t359) * t471 - t380 * t433 - t589 * t488) * t547) * MDP(35) + (t360 * t488 + t364 * t471 - t430 * t433 - t441 * t571) * MDP(31) + (-t360 * t429 + t361 * t430 - t364 * t396 + t365 * t571) * MDP(30) + (t433 * t486 + t441 * t460 + t446 * t488 + t455 * t677) * MDP(20) + (-t432 * t488 - t433 * t489 - t441 * t482 - t442 * t677) * MDP(16) + ((-t451 * t620 + t382) * t695 + t414 * t433 + (-t412 * t620 + t377) * t488 + t373 * t441 + t402 * t447 + t450 * t392 + t411 * t606 + ((-qJD(5) * t422 - t401) * t695 - t451 * t433 + (-qJD(5) * t406 - t385) * t488 + t386 * t489 + t411 * t442) * t544) * MDP(27) + (-t374 * t441 + t386 * t656 + t450 * t391 + t402 * t449 + t411 * t567 - t433 * t638 - t488 * t565 - t564 * t695) * MDP(28) + (-t392 * t488 - t441 * t447 - t489 * t652 - t568 * t695) * MDP(25) + (t391 * t488 + t441 * t449 + t489 * t649 + t567 * t695) * MDP(24) + (t441 * t695 + t409) * MDP(26) + t554 * MDP(12) + (pkin(7) * t554 + t459 * t503 + t475 * t484) * MDP(14) + ((t356 * t547 - t359 * t543) * t471 + (t375 * t547 - t380 * t543) * t433 + t596 * t488 + t357 * t441 + t371 * t396 + t423 * t361 + t363 * t429 + t387 * t365 + ((-t375 * t543 - t380 * t547) * t471 - t358 * t488) * qJD(6)) * MDP(34) + (-MDP(17) * t442 + MDP(18) * t441 + MDP(20) * t402 + MDP(21) * t401) * t537 + (t391 * t656 + t449 * t567) * MDP(22) + ((-t447 * t548 - t449 * t544) * t442 + (-t666 - t392 * t548 + (t447 * t544 - t449 * t548) * qJD(5)) * t489) * MDP(23) - MDP(7) * t650 + (pkin(7) * t650 + t550 * t598) * MDP(10) + (-pkin(7) * t647 + t546 * t598) * MDP(9) + (t432 * t486 + t442 * t460 + t446 * t489 + t455 * t482) * MDP(21) + (t546 * t569 - t590 * t625) * MDP(13) + (t550 * t569 + t590 * t627) * MDP(11) + 0.2e1 * t546 * MDP(4) * t600 + MDP(6) * t647 + (t432 * t489 + t442 * t482) * MDP(15) + (-t361 * t488 - t365 * t471 - t396 * t441 - t429 * t433) * MDP(32) + (t441 * t471 + t409) * MDP(33) + t686 * t687; -t711 + 0.2e1 * t538 * MDP(13) + (qJ(3) * t501 + qJD(3) * t504 - t484 * t492) * MDP(14) + (MDP(9) * t546 * t553 + MDP(10) * t646) * pkin(1) - MDP(4) * t610 + (-t474 * t677 + t537 * t636 + t558) * MDP(20) + (qJD(4) * t677 - t635 + t660) * MDP(17) + (t493 * t392 + (-t494 * t620 - t405) * t695 + t636 * t447 + (t643 * t695 + t561) * t544 - t577) * MDP(27) + (t493 * t391 + (t494 * t621 + t644) * t695 + t636 * t449 + (-t467 * t695 + t561) * t548 - t576) * MDP(28) + ((-t484 * t546 + t492 * t550) * MDP(11) + (t484 * t550 + t492 * t546) * MDP(13) + (t504 * t546 + (-t502 - t670) * t550) * pkin(7) * MDP(14) + ((t504 - t539) * t546 + (-t502 + t597) * t550) * MDP(12)) * qJD(1) + (-t474 * t482 - t537 * t643 + t557) * MDP(21) + ((t461 * t547 + t462 * t543) * t433 + t477 * t361 + (t543 * t583 - t547 * t585) * t471 + t637 * t396 + t706) * MDP(34) + (-(t461 * t543 - t462 * t547) * t433 + t477 * t360 + (t543 * t585 + t547 * t583) * t471 - t637 * t571 + t707) * MDP(35) + t553 * t686; -MDP(11) * t610 + (-t540 * t553 - t552) * MDP(13) + (-qJD(2) * t504 + t484 * t629 + t520) * MDP(14) + (-t545 * t588 - t629 * t677) * MDP(20) + (-t482 * t629 - t549 * t588) * MDP(21) + (-t392 * t549 + (-t544 * t622 - t478) * t695 + (-t447 * t537 - t620 * t695 - t652) * t545) * MDP(27) + (-t391 * t549 + (-t548 * t622 + t481) * t695 + (-t449 * t537 + t621 * t695 - t649) * t545) * MDP(28) + (-(t478 * t547 - t481 * t543) * t471 + (-t490 * t624 - t361) * t549 + ((t560 - t681) * t471 - t663 - t537 * t396) * t545) * MDP(34) + ((t478 * t543 + t481 * t547) * t471 + (t487 * t624 - t360) * t549 + (-(t543 * t679 - t544 * t618 - t547 * t621) * t471 + t664 + t537 * t571) * t545) * MDP(35); (t432 - t660) * MDP(17) + (-t425 * t537 - t558) * MDP(20) + (-t424 * t537 - t557) * MDP(21) + (-pkin(4) * t392 - t425 * t447 + (-pkin(9) * t620 - t428) * t695 + (t424 * t695 + t566) * t544 + t577) * MDP(27) + (-pkin(4) * t391 - t425 * t449 + (pkin(9) * t621 + t645) * t695 + t566 * t548 + t576) * MDP(28) + ((-t505 * t547 - t507 * t543) * t433 + t527 * t361 + (t543 * t582 - t547 * t584) * t471 + t586 * t396 - t706) * MDP(34) + (-(-t505 * t543 + t507 * t547) * t433 + t527 * t360 + (t543 * t584 + t547 * t582) * t471 - t586 * t571 - t707) * MDP(35) + t711; t449 * t447 * MDP(22) + (-t447 ^ 2 + t449 ^ 2) * MDP(23) + (t391 + t710) * MDP(24) + (-t665 + (-qJD(5) + t695) * t449) * MDP(25) + t433 * MDP(26) + (t374 * t695 - t411 * t449 + t556) * MDP(27) + (t373 * t695 + t411 * t447 - t565) * MDP(28) + (t360 + t699) * MDP(31) + (-t361 - t698) * MDP(32) + (-(-t367 * t543 - t667) * t471 - t358 * qJD(6) + (-t396 * t449 + t433 * t547 - t471 * t619) * pkin(5) + t693) * MDP(34) + ((-t368 * t471 - t354) * t543 + (t367 * t471 - t589) * t547 + (-t433 * t543 + t449 * t571 - t471 * t618) * pkin(5) + t694) * MDP(35) + t688; (t609 + t699) * MDP(31) + (-t595 - t698) * MDP(32) + (t358 * t471 + t693) * MDP(34) + (-t543 * t354 - t547 * t355 + t357 * t471 + t694) * MDP(35) + (-MDP(31) * t661 + MDP(32) * t571 - MDP(34) * t358 - MDP(35) * t668) * qJD(6) + t688;];
tauc  = t1;
