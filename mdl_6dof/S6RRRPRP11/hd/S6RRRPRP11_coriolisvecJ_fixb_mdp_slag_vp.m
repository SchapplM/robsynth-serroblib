% Calculate minimal parameter regressor of Coriolis joint torque vector for
% S6RRRPRP11
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d1,d2,d3,d5]';
% MDP [30x1]
%   Minimal dynamic parameter vector (fixed base model)
%   see S6RRRPRP11_convert_par2_MPV_fixb.m
% 
% Output:
% tauc [6x1]
%   joint torques required to compensate Coriolis and centrifugal load

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 17:50
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauc = S6RRRPRP11_coriolisvecJ_fixb_mdp_slag_vp(qJ, qJD, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(10,1),zeros(30,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRPRP11_coriolisvecJ_fixb_mdp_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRRPRP11_coriolisvecJ_fixb_mdp_slag_vp: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RRRPRP11_coriolisvecJ_fixb_mdp_slag_vp: pkin has to be [10x1] (double)');
assert(isreal(MDP) && all(size(MDP) == [30 1]), ...
  'S6RRRPRP11_coriolisvecJ_fixb_mdp_slag_vp: MDP has to be [30x1] (double)'); 

%% Symbolic Calculation
% From coriolisvec_joint_fixb_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 17:49:22
% EndTime: 2019-03-09 17:49:40
% DurationCPUTime: 10.87s
% Computational Cost: add. (6905->576), mult. (17732->747), div. (0->0), fcn. (13231->8), ass. (0->229)
t547 = sin(qJ(3));
t545 = sin(pkin(6));
t551 = cos(qJ(2));
t638 = qJD(1) * t551;
t614 = t545 * t638;
t590 = t547 * t614;
t635 = qJD(3) * t547;
t681 = t590 - t635;
t526 = -qJD(3) + t614;
t548 = sin(qJ(2));
t671 = cos(pkin(6));
t605 = t671 * qJD(1);
t587 = pkin(1) * t605;
t533 = t548 * t587;
t495 = pkin(8) * t614 + t533;
t696 = pkin(3) * t681 + qJD(4) * t547 + t495;
t639 = qJD(1) * t545;
t615 = t548 * t639;
t492 = -pkin(8) * t615 + t551 * t587;
t570 = (pkin(2) * t548 - pkin(9) * t551) * t545;
t493 = qJD(1) * t570;
t550 = cos(qJ(3));
t600 = -t547 * t492 + t493 * t550;
t634 = qJD(3) * t550;
t656 = t550 * t551;
t674 = pkin(4) + pkin(9);
t675 = pkin(3) + pkin(10);
t695 = (pkin(4) * t656 - t548 * t675) * t639 - t600 - t674 * t634;
t694 = t696 + t526 * (pkin(10) * t547 - qJ(4) * t550);
t682 = t550 * t614 - t634;
t574 = t605 + qJD(2);
t475 = t547 * t574 + t550 * t615;
t454 = pkin(9) * t574 + t495;
t488 = (-pkin(2) * t551 - pkin(9) * t548 - pkin(1)) * t545;
t467 = qJD(1) * t488;
t417 = t547 * t454 - t550 * t467;
t626 = -qJD(4) - t417;
t679 = t475 * pkin(4) - t626;
t387 = t526 * t675 + t679;
t564 = t550 * t574;
t473 = t547 * t615 - t564;
t453 = -pkin(2) * t574 - t492;
t556 = -t475 * qJ(4) + t453;
t391 = t473 * t675 + t556;
t546 = sin(qJ(5));
t549 = cos(qJ(5));
t360 = t387 * t546 + t391 * t549;
t624 = qJD(1) * qJD(2);
t607 = t545 * t624;
t585 = t551 * t607;
t660 = t545 * t548;
t621 = t547 * t660;
t588 = qJD(3) * t621;
t440 = qJD(1) * t588 - qJD(3) * t564 - t550 * t585;
t637 = qJD(2) * t548;
t613 = t545 * t637;
t577 = t675 * t613;
t494 = qJD(2) * t570;
t483 = qJD(1) * t494;
t616 = pkin(1) * t671;
t680 = -pkin(8) * t660 + t551 * t616;
t496 = t680 * qJD(2);
t484 = qJD(1) * t496;
t593 = t454 * t634 + t467 * t635 - t550 * t483 + t547 * t484;
t369 = -pkin(4) * t440 - qJD(1) * t577 + t593;
t485 = pkin(8) * t585 + qJD(2) * t533;
t627 = t475 * qJD(3);
t683 = t547 * t585 + t627;
t381 = pkin(3) * t683 + t440 * qJ(4) - t475 * qJD(4) + t485;
t371 = pkin(10) * t683 + t381;
t603 = t549 * t369 - t371 * t546;
t354 = -qJD(5) * t360 + t603;
t586 = t548 * t607;
t631 = qJD(5) * t549;
t618 = t473 * t631 + t546 * t683 + t549 * t586;
t632 = qJD(5) * t546;
t392 = -t526 * t632 - t618;
t437 = t473 * t546 - t526 * t549;
t349 = -pkin(5) * t440 + qJ(6) * t392 - qJD(6) * t437 + t354;
t435 = -t549 * t473 - t526 * t546;
t358 = -qJ(6) * t435 + t360;
t468 = qJD(5) + t475;
t693 = t358 * t468 + t349;
t619 = -t546 * t369 - t549 * t371 - t387 * t631;
t353 = -t391 * t632 - t619;
t563 = t546 * t586 - t549 * t683;
t633 = qJD(5) * t437;
t393 = t563 + t633;
t350 = -qJ(6) * t393 - qJD(6) * t435 + t353;
t359 = t549 * t387 - t391 * t546;
t357 = -qJ(6) * t437 + t359;
t355 = pkin(5) * t468 + t357;
t692 = -t355 * t468 + t350;
t527 = t674 * t547;
t691 = -qJD(5) * t527 + t694;
t461 = t546 * t615 - t549 * t590;
t582 = t549 * t635 + t461;
t542 = t545 ^ 2;
t690 = -0.2e1 * t542 * t624;
t689 = MDP(5) * (t548 ^ 2 - t551 ^ 2);
t686 = t435 * t468;
t592 = t548 * t616;
t659 = t545 * t551;
t487 = pkin(8) * t659 + pkin(9) * t671 + t592;
t601 = -t547 * t487 + t488 * t550;
t422 = pkin(3) * t659 - t601;
t501 = t547 * t671 + t550 * t660;
t402 = pkin(4) * t501 + pkin(10) * t659 + t422;
t500 = -t550 * t671 + t621;
t486 = -pkin(2) * t671 - t680;
t557 = -t501 * qJ(4) + t486;
t410 = t500 * t675 + t557;
t650 = t546 * t402 + t549 * t410;
t642 = -t682 * qJ(4) + t696;
t685 = t695 * t549;
t606 = -qJ(4) * t547 - pkin(2);
t504 = -t550 * t675 + t606;
t684 = t504 * t632 - t527 * t631 + t695 * t546 + t694 * t549;
t418 = t550 * t454 + t547 * t467;
t401 = -pkin(4) * t473 + t418;
t515 = t526 * qJ(4);
t394 = t401 - t515;
t678 = t394 * t468 + t440 * t675;
t677 = t437 ^ 2;
t676 = t475 ^ 2;
t553 = qJD(1) ^ 2;
t673 = pkin(5) * t549;
t672 = pkin(9) * t547;
t540 = t550 * pkin(9);
t670 = qJ(4) * t473;
t669 = qJ(6) * t475;
t521 = qJ(4) * t586;
t594 = t454 * t635 - t467 * t634 - t547 * t483 - t550 * t484;
t578 = -qJD(4) * t526 - t594;
t379 = -t521 - t578;
t367 = -pkin(4) * t683 - t379;
t668 = t367 * t546;
t667 = t367 * t549;
t666 = t392 * t549;
t409 = t473 * pkin(3) + t556;
t665 = t409 * t475;
t664 = t437 * t468;
t663 = t440 * t546;
t433 = t440 * t547;
t662 = t473 * t475;
t661 = t542 * t553;
t658 = t547 * t551;
t431 = t549 * t440;
t657 = t549 * t550;
t654 = qJ(6) + t675;
t653 = t355 - t357;
t462 = (t546 * t658 + t548 * t549) * t639;
t604 = qJ(6) * t550 - t504;
t628 = qJD(6) * t550;
t652 = qJ(6) * t462 + t604 * t631 - t685 + (-qJ(6) * t635 + t628 + t691) * t546 - t682 * pkin(5);
t630 = qJD(5) * t550;
t610 = t546 * t630;
t651 = -t549 * t628 - t684 + (t610 + t582) * qJ(6);
t414 = t475 * t675 + t670;
t649 = t546 * t401 + t549 * t414;
t647 = t550 * t487 + t547 * t488;
t646 = t550 * t492 + t547 * t493;
t396 = t549 * t401;
t645 = -qJD(6) * t549 + t632 * t654 + pkin(5) * t473 - t396 - (-t414 - t669) * t546;
t518 = t654 * t549;
t644 = -qJD(5) * t518 - qJD(6) * t546 - t549 * t669 - t649;
t643 = t549 * t504 + t546 * t527;
t424 = -qJ(4) * t615 - t646;
t560 = pkin(4) * t590 + t424;
t641 = -t674 * t635 + t560;
t609 = qJD(2) * t659;
t497 = pkin(8) * t609 + qJD(2) * t592;
t528 = t550 * pkin(4) + t540;
t636 = qJD(2) * t550;
t629 = qJD(5) * t675;
t623 = t526 * t672;
t622 = pkin(9) * t636;
t620 = t549 * t659;
t612 = t526 * t634;
t602 = t549 * t402 - t410 * t546;
t599 = t468 ^ 2;
t598 = t468 * t546;
t595 = t542 * t548 * t551 * MDP(4);
t583 = pkin(1) * t690;
t581 = t546 * t635 - t462;
t579 = pkin(3) * t586;
t382 = -t579 + t593;
t576 = -t379 * t550 + t382 * t547;
t421 = qJ(4) * t659 - t647;
t572 = -t487 * t634 - t488 * t635 + t494 * t550 - t547 * t496;
t447 = t500 * t549 + t546 * t659;
t569 = t526 * t550;
t568 = -t418 * t526 - t593;
t567 = -t487 * t635 + t488 * t634 + t547 * t494 + t550 * t496;
t446 = -t588 + (qJD(3) * t671 + t609) * t550;
t375 = pkin(4) * t446 - t572 - t577;
t445 = qJD(3) * t501 + t547 * t609;
t561 = -qJ(4) * t446 - qJD(4) * t501 + t497;
t380 = t445 * t675 + t561;
t566 = t546 * t375 + t549 * t380 + t402 * t631 - t410 * t632;
t411 = -pkin(4) * t500 - t421;
t559 = -qJD(5) * t650 + t549 * t375 - t380 * t546;
t558 = pkin(9) * t612 - t586 * t672;
t384 = -qJ(4) * t613 + qJD(4) * t659 - t567;
t376 = -pkin(4) * t445 - t384;
t356 = t393 * pkin(5) + t367;
t520 = -pkin(3) * t550 + t606;
t517 = t654 * t546;
t508 = t549 * t527;
t448 = -t500 * t546 + t620;
t439 = -qJ(6) * t657 + t643;
t434 = t435 ^ 2;
t428 = pkin(5) * t547 + t546 * t604 + t508;
t427 = pkin(3) * t475 + t670;
t426 = -pkin(3) * t615 - t600;
t423 = t440 * t501;
t420 = t500 * pkin(3) + t557;
t415 = t515 - t418;
t413 = pkin(3) * t526 - t626;
t407 = qJD(5) * t447 + t445 * t546 + t549 * t613;
t406 = -t445 * t549 - qJD(5) * t620 + (qJD(5) * t500 + t613) * t546;
t390 = pkin(3) * t445 + t561;
t388 = -pkin(3) * t613 - t572;
t377 = pkin(5) * t435 + qJD(6) + t394;
t363 = qJ(6) * t447 + t650;
t362 = pkin(5) * t501 + qJ(6) * t448 + t602;
t352 = -qJ(6) * t406 + qJD(6) * t447 + t566;
t351 = pkin(5) * t446 - qJ(6) * t407 + qJD(6) * t448 + t559;
t1 = [(-t417 * t613 + t453 * t445 + t497 * t473 + t485 * t500 + t486 * t683 - t526 * t572 + t586 * t601 + t593 * t659) * MDP(16) + (t379 * t500 + t382 * t501 + t384 * t473 + t388 * t475 + t413 * t446 + t415 * t445 + t421 * t683 - t422 * t440) * MDP(18) + (-t381 * t500 - t382 * t659 - t388 * t526 - t390 * t473 - t409 * t445 + t413 * t613 - t420 * t683 + t422 * t586) * MDP(19) + (t445 * t526 - t473 * t613 - t500 * t586 + t659 * t683) * MDP(14) + (t440 * t500 - t475 * t445 - t446 * t473 - t501 * t683) * MDP(12) + (t379 * t421 + t381 * t420 + t382 * t422 + t384 * t415 + t388 * t413 + t390 * t409) * MDP(21) + (-t484 * t671 - t496 * t574 + t551 * t583) * MDP(10) + (-t485 * t671 - t497 * t574 + t548 * t583) * MDP(9) + (-t526 * t545 - t542 * t638) * MDP(15) * t637 + (t567 * t526 + t497 * t475 - t486 * t440 + t485 * t501 + t453 * t446 + (-t594 * t551 + (-qJD(1) * t647 - t418) * t637) * t545) * MDP(17) + (-t353 * t501 - t360 * t446 - t367 * t448 + t376 * t437 - t411 * t392 + t394 * t407 + t440 * t650 - t468 * t566) * MDP(28) + (t446 * t468 - t423) * MDP(26) + (t446 * t475 - t423) * MDP(11) + (-t446 * t526 + (t440 * t551 + (qJD(1) * t501 + t475) * t637) * t545) * MDP(13) + (-t381 * t501 + t384 * t526 - t390 * t475 - t409 * t446 + t420 * t440 + (t379 * t551 + (-qJD(1) * t421 - t415) * t637) * t545) * MDP(20) + (MDP(6) * t609 - MDP(7) * t613) * (0.2e1 * t605 + qJD(2)) + 0.2e1 * t595 * t624 + (t354 * t501 + t359 * t446 - t367 * t447 + t376 * t435 + t411 * t393 + t394 * t406 - t440 * t602 + t468 * t559) * MDP(27) + t689 * t690 + (t350 * t363 + t358 * t352 + t349 * t362 + t355 * t351 + t356 * (-pkin(5) * t447 + t411) + t377 * (pkin(5) * t406 + t376)) * MDP(30) + (t349 * t448 + t350 * t447 - t351 * t437 - t352 * t435 - t355 * t407 - t358 * t406 + t362 * t392 - t363 * t393) * MDP(29) + (t392 * t448 + t407 * t437) * MDP(22) + (-t392 * t447 + t393 * t448 - t406 * t437 - t407 * t435) * MDP(23) + (-t393 * t501 - t406 * t468 - t435 * t446 - t440 * t447) * MDP(25) + (-t392 * t501 + t407 * t468 + t437 * t446 + t440 * t448) * MDP(24); (-(-t504 * t546 + t508) * t440 + t354 * t547 + t528 * t393 + (-t504 * t631 + t546 * t691 - t685) * t468 + t641 * t435 - t582 * t394 + (-t359 * t526 - t394 * t632 + t667) * t550) * MDP(27) + (-pkin(2) * t683 + t417 * t615 - t453 * t681 - t495 * t473 - t485 * t550 + t526 * t600 + t558) * MDP(16) + (-t424 * t473 - t426 * t475 + t576 + (-t683 + t627) * t540 - t681 * t415 - t682 * t413 + (t473 * t635 - t433) * pkin(9)) * MDP(18) + (t381 * t550 + t409 * t681 - t413 * t615 + t426 * t526 + t473 * t642 - t520 * t683 - t558) * MDP(19) + (-t440 * t550 + t475 * t590 + (-t683 - t627) * t547 + t682 * t473) * MDP(12) + (t643 * t440 - t353 * t547 - t528 * t392 + t684 * t468 + t641 * t437 + t581 * t394 + (t360 * t526 - t394 * t631 - t668) * t550) * MDP(28) + (t435 * t462 + t437 * t461 + (-t435 * t546 + t437 * t549) * t635 + (t666 + t393 * t546 + (t435 * t549 + t437 * t546) * qJD(5)) * t550) * MDP(23) + (-t392 * t547 + t581 * t468 + (-t437 * t526 - t468 * t631 + t663) * t550) * MDP(24) + (pkin(1) * t548 * t661 + t495 * t574 - t485) * MDP(9) + (pkin(8) * t586 + t492 * t574 + (-qJD(2) * t605 + t661) * t551 * pkin(1)) * MDP(10) + (-t393 * t547 + t582 * t468 + (t435 * t526 + t468 * t632 + t431) * t550) * MDP(25) + (t526 * t635 + (-t526 * t658 + (t473 + t636) * t548) * t639) * MDP(14) + (pkin(2) * t440 + t485 * t547 - t646 * t526 - t495 * t475 + (t453 * t550 - t623) * qJD(3) + (-t453 * t656 + (t418 - t622) * t548) * t639) * MDP(17) + (-t612 + (t526 * t656 + (t547 * qJD(2) - t475) * t548) * t639) * MDP(13) + (-t381 * t547 - t424 * t526 + t440 * t520 + t642 * t475 + (-t409 * t550 + t623) * qJD(3) + (t409 * t656 + (t415 + t622) * t548) * t639) * MDP(20) + (t355 * t462 + t358 * t461 + t392 * t428 - t393 * t439 - t652 * t437 - t651 * t435 + (-t355 * t546 + t358 * t549) * t635 + (t349 * t546 - t350 * t549 + (t355 * t549 + t358 * t546) * qJD(5)) * t550) * MDP(29) + (t381 * t520 - t413 * t426 - t415 * t424 - t642 * t409 + ((t413 * t550 + t415 * t547) * qJD(3) + t576) * pkin(9)) * MDP(21) + (t392 * t546 * t550 + (-t549 * t630 + t581) * t437) * MDP(22) + (t350 * t439 + t349 * t428 + t356 * (pkin(5) * t657 + t528) + t651 * t358 + t652 * t355 + ((-t461 - t610) * pkin(5) + t560 + (-t673 - t674) * t635) * t377) * MDP(30) + (-t475 * t569 - t433) * MDP(11) + (-t468 * t569 - t433) * MDP(26) + t661 * t689 + t526 * MDP(15) * t615 + (-t595 + (-MDP(6) * t551 + MDP(7) * t548) * t545 * t671) * t553; MDP(11) * t662 + t676 * MDP(12) - t440 * MDP(13) + (-t475 * t526 - t683) * MDP(14) + MDP(15) * t586 + (-t453 * t475 + t568) * MDP(16) + (t417 * t526 + t594) * MDP(17) + (pkin(3) * t440 - qJ(4) * t683 + (-t415 - t418) * t475) * MDP(18) + (-t568 - 0.2e1 * t579 + t665) * MDP(19) + (t427 * t475 + t526 * t626 + 0.2e1 * t521 + t578) * MDP(20) + (-pkin(3) * t382 - qJ(4) * t379 - t409 * t427 - t413 * t418 + t415 * t626) * MDP(21) + (-t437 * t598 - t666) * MDP(22) + ((-t393 - t664) * t549 + (t392 + t686) * t546) * MDP(23) + (-t468 * t598 - t431) * MDP(24) + (-t549 * t599 + t663) * MDP(25) + (qJ(4) * t393 + t668 + (-t396 + (t414 + t629) * t546) * t468 + t679 * t435 + t678 * t549) * MDP(27) + (-qJ(4) * t392 + t667 + (t549 * t629 + t649) * t468 + t679 * t437 - t678 * t546) * MDP(28) + (-t392 * t518 + t393 * t517 - t644 * t435 - t645 * t437 - t692 * t546 - t549 * t693) * MDP(29) + (-t350 * t517 - t349 * t518 + t356 * (pkin(5) * t546 + qJ(4)) + (t468 * t673 + t679) * t377 + t644 * t358 + t645 * t355) * MDP(30) + (-t526 * MDP(13) + t453 * MDP(17) + (t413 + t626) * MDP(18) + t427 * MDP(19) - t409 * MDP(20) + t437 * MDP(24) - t435 * MDP(25) + t468 * MDP(26) + t359 * MDP(27) - t360 * MDP(28) - t473 * MDP(12)) * t473; -t440 * MDP(18) + (t586 - t662) * MDP(19) - t676 * MDP(20) + (t382 + t665) * MDP(21) - t431 * MDP(27) + (-MDP(18) * t473 - MDP(20) * t526 - MDP(21) * t415 + MDP(27) * t435 + MDP(28) * t437 + MDP(30) * t377) * t526 + ((t392 - t686) * MDP(29) + t693 * MDP(30) - MDP(28) * t599) * t549 + (t440 * MDP(28) + (t437 * t475 - t393 + t633) * MDP(29) + t692 * MDP(30) - MDP(27) * t599) * t546; (-t434 + t677) * MDP(23) + t618 * MDP(24) + (-t563 + t664) * MDP(25) - t440 * MDP(26) + (t360 * t468 - t394 * t437 + t603) * MDP(27) + (t359 * t468 + t619) * MDP(28) + t653 * MDP(30) * t358 + (t392 * MDP(29) + (-t377 * t437 + t349) * MDP(30)) * pkin(5) + (t437 * MDP(22) + t468 * MDP(24) + t394 * MDP(28) - MDP(29) * t653) * t435 + ((MDP(25) * t526 - MDP(27) * t391) * t549 + (MDP(24) * t526 - MDP(25) * t473 - MDP(27) * t387 + MDP(28) * t391) * t546) * qJD(5); (-t434 - t677) * MDP(29) + (t355 * t437 + t358 * t435 + t356) * MDP(30);];
tauc  = t1;
