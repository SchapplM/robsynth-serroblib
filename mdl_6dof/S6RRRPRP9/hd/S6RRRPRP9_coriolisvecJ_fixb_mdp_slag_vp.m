% Calculate minimal parameter regressor of Coriolis joint torque vector for
% S6RRRPRP9
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d3,d5]';
% MDP [32x1]
%   Minimal dynamic parameter vector (fixed base model)
%   see S6RRRPRP9_convert_par2_MPV_fixb.m
% 
% Output:
% tauc [6x1]
%   joint torques required to compensate Coriolis and centrifugal load

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 17:28
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauc = S6RRRPRP9_coriolisvecJ_fixb_mdp_slag_vp(qJ, qJD, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(9,1),zeros(32,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRPRP9_coriolisvecJ_fixb_mdp_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRRPRP9_coriolisvecJ_fixb_mdp_slag_vp: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S6RRRPRP9_coriolisvecJ_fixb_mdp_slag_vp: pkin has to be [9x1] (double)');
assert(isreal(MDP) && all(size(MDP) == [32 1]), ...
  'S6RRRPRP9_coriolisvecJ_fixb_mdp_slag_vp: MDP has to be [32x1] (double)'); 

%% Symbolic Calculation
% From coriolisvec_joint_fixb_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 17:27:30
% EndTime: 2019-03-09 17:27:44
% DurationCPUTime: 8.85s
% Computational Cost: add. (5878->601), mult. (13376->752), div. (0->0), fcn. (8518->6), ass. (0->245)
t561 = cos(qJ(3));
t651 = qJD(3) * t561;
t562 = cos(qJ(2));
t660 = qJD(1) * t562;
t732 = t561 * t660 - t651;
t558 = sin(qJ(3));
t626 = t558 * t660;
t653 = qJD(3) * t558;
t731 = t626 - t653;
t639 = qJD(1) * qJD(2);
t614 = t562 * t639;
t559 = sin(qJ(2));
t652 = qJD(3) * t559;
t620 = t558 * t652;
t638 = qJD(2) * qJD(3);
t460 = qJD(1) * t620 + (-t614 - t638) * t561;
t661 = qJD(1) * t559;
t627 = t558 * t661;
t642 = t561 * qJD(2);
t493 = t627 - t642;
t625 = t561 * t661;
t657 = qJD(2) * t558;
t495 = t625 + t657;
t557 = sin(qJ(5));
t560 = cos(qJ(5));
t618 = t559 * t651;
t655 = qJD(2) * t562;
t623 = t558 * t655;
t572 = t618 + t623;
t613 = t558 * t638;
t567 = qJD(1) * t572 + t613;
t648 = qJD(5) * t560;
t649 = qJD(5) * t557;
t390 = t560 * t460 - t493 * t648 + t495 * t649 - t557 * t567;
t535 = -qJD(3) + t660;
t524 = qJD(5) + t535;
t584 = -t560 * t493 + t495 * t557;
t724 = -t524 * t584 + t390;
t730 = MDP(30) * t724;
t701 = pkin(8) * t559;
t593 = -pkin(2) * t562 - t701;
t508 = -pkin(1) + t593;
t481 = t508 * qJD(1);
t545 = pkin(7) * t660;
t515 = qJD(2) * pkin(8) + t545;
t442 = t561 * t481 - t558 * t515;
t641 = qJD(4) - t442;
t702 = pkin(8) - pkin(9);
t517 = t702 * t561;
t628 = -pkin(7) * t558 - pkin(3);
t685 = t561 * t562;
t570 = -pkin(9) * t685 + (-pkin(4) + t628) * t559;
t592 = pkin(2) * t559 - pkin(8) * t562;
t503 = t592 * qJD(1);
t694 = t503 * t561;
t729 = qJD(1) * t570 - qJD(3) * t517 - t694;
t479 = t558 * t503;
t669 = qJ(4) * t661 + t479;
t687 = t559 * t561;
t688 = t558 * t562;
t728 = (-pkin(7) * t687 + pkin(9) * t688) * qJD(1) + t669 + t702 * t653;
t721 = pkin(9) * t495 - t641;
t439 = t493 * t557 + t495 * t560;
t727 = t439 ^ 2;
t556 = qJD(2) * pkin(2);
t514 = pkin(7) * t661 - t556;
t432 = t493 * pkin(3) - t495 * qJ(4) + t514;
t417 = -pkin(4) * t493 - t432;
t376 = pkin(5) * t584 - qJ(6) * t439 + t417;
t726 = t376 * t439;
t615 = t559 * t639;
t534 = pkin(5) * t615;
t622 = t559 * t642;
t506 = t592 * qJD(2);
t483 = qJD(1) * t506;
t672 = t481 * t651 + t558 * t483;
t700 = pkin(9) * qJD(2);
t518 = t535 * qJD(4);
t526 = qJ(4) * t615;
t706 = t526 - t518;
t382 = (-t515 + t700) * t653 + (-pkin(7) * t622 + pkin(9) * t572) * qJD(1) + t672 + t706;
t601 = pkin(7) * t615;
t598 = t481 * t653 - t561 * t483 + t515 * t651 - t558 * t601;
t703 = pkin(3) + pkin(4);
t383 = pkin(9) * t460 - t615 * t703 + t598;
t406 = t535 * t703 - t721;
t443 = t558 * t481 + t561 * t515;
t421 = pkin(9) * t493 + t443;
t520 = t535 * qJ(4);
t415 = t421 - t520;
t599 = t557 * t382 - t560 * t383 + t406 * t649 + t415 * t648;
t366 = t534 + t599;
t375 = t406 * t557 + t415 * t560;
t372 = qJ(6) * t524 + t375;
t725 = t372 * t524 - t366;
t497 = t557 * t558 + t560 * t561;
t445 = t497 * qJD(5) - t557 * t653 - t560 * t651;
t575 = t497 * t562;
t675 = qJD(1) * t575 + t445;
t674 = -t557 * t732 + t558 * t648 + t731 * t560 - t561 * t649;
t723 = qJ(4) * t732 - t558 * qJD(4) - t545;
t643 = t559 * MDP(26);
t611 = qJD(2) * t643;
t722 = -MDP(24) * t724 - qJD(1) * t611;
t401 = pkin(5) * t439 + qJ(6) * t584;
t720 = t439 * MDP(22) - MDP(23) * t584 + t417 * MDP(28) - t376 * MDP(31);
t719 = t524 ^ 2;
t718 = -0.2e1 * t639;
t717 = MDP(4) * t559;
t554 = t559 ^ 2;
t716 = MDP(5) * (-t562 ^ 2 + t554);
t604 = -t495 + t657;
t712 = t562 * t604;
t711 = qJ(4) * t649 + t421 * t557 + t721 * t560 + t648 * t703;
t602 = pkin(3) * t615;
t402 = t598 - t602;
t429 = -t520 + t443;
t710 = t429 * t535 + t402;
t516 = t702 * t558;
t583 = t516 * t560 - t517 * t557;
t709 = -qJD(5) * t583 + t729 * t557 + t728 * t560;
t457 = t516 * t557 + t517 * t560;
t708 = -qJD(5) * t457 + t728 * t557 - t729 * t560;
t537 = pkin(7) * t688;
t553 = t562 * pkin(3);
t433 = pkin(4) * t562 + t537 + t553 + (-pkin(9) * t559 - t508) * t561;
t538 = pkin(7) * t685;
t665 = t558 * t508 + t538;
t458 = -qJ(4) * t562 + t665;
t690 = t558 * t559;
t441 = pkin(9) * t690 + t458;
t707 = t557 * t433 + t560 * t441;
t631 = t703 * t558;
t673 = -qJD(3) * t631 + t626 * t703 - t723;
t663 = t560 * qJ(4) - t557 * t703;
t637 = -MDP(27) - MDP(29);
t636 = MDP(28) - MDP(31);
t591 = -qJD(3) * t538 + t506 * t561 - t508 * t653;
t398 = pkin(9) * t620 + qJD(2) * t570 - t591;
t656 = qJD(2) * t559;
t668 = t558 * t506 + t508 * t651;
t629 = qJ(4) * t656 + t668;
t400 = (-pkin(7) * qJD(2) + pkin(9) * qJD(3)) * t687 + (-qJD(4) + (-pkin(7) * qJD(3) + t700) * t558) * t562 + t629;
t704 = -qJD(5) * t707 + t398 * t560 - t400 * t557;
t630 = -pkin(3) * t567 - t460 * qJ(4) + t495 * qJD(4);
t399 = pkin(7) * t614 - t630;
t699 = t399 * t558;
t698 = t399 * t561;
t696 = t432 * t562;
t695 = t493 * t535;
t693 = t514 * t561;
t692 = t535 * t558;
t691 = t535 * t561;
t549 = t558 * qJ(4);
t689 = t558 * t560;
t564 = qJD(2) ^ 2;
t686 = t559 * t564;
t684 = t562 * t535;
t683 = t562 * t564;
t565 = qJD(1) ^ 2;
t682 = t562 * t565;
t582 = t557 * t561 - t689;
t681 = pkin(5) * t674 + qJ(6) * t675 + qJD(6) * t582 + t673;
t680 = qJD(6) + t711;
t679 = -qJ(6) * t661 + t709;
t678 = pkin(5) * t661 + t708;
t671 = -qJD(5) * t663 - t421 * t560 + t557 * t721;
t670 = -pkin(3) * t731 + t723;
t449 = t495 * pkin(3) + t493 * qJ(4);
t621 = t562 * t642;
t666 = qJ(4) * t621 + qJD(4) * t687;
t659 = qJD(2) * t583;
t658 = qJD(2) * t457;
t654 = qJD(3) * t493;
t650 = qJD(5) * t439;
t647 = qJD(6) * t524;
t645 = t554 * qJD(1);
t644 = t559 * MDP(15);
t374 = t406 * t560 - t415 * t557;
t640 = qJD(6) - t374;
t635 = pkin(8) * t692;
t634 = pkin(8) * t691;
t507 = -t561 * pkin(3) - pkin(2) - t549;
t633 = pkin(8) * t642;
t632 = pkin(7) * t655;
t619 = t562 * t653;
t612 = qJD(2) * t644;
t610 = pkin(1) * t718;
t609 = -t460 * t557 - t560 * t567;
t608 = t508 * t561 - t537;
t488 = t561 * pkin(4) - t507;
t607 = qJD(1) * t458 + t429;
t606 = t535 + t660;
t605 = t493 + t642;
t603 = qJD(3) + t660;
t600 = t560 * t382 + t557 * t383 + t406 * t648 - t415 * t649;
t597 = qJ(6) * t615;
t596 = -pkin(7) - t631;
t423 = -pkin(4) * t495 - t449;
t594 = t559 * t628;
t590 = pkin(7) * t495 + t693;
t589 = -qJ(4) * t557 - t560 * t703;
t586 = t433 * t560 - t441 * t557;
t581 = t561 * t603;
t579 = -t515 * t653 + t672;
t578 = t374 * t524 - t600;
t577 = t375 * t524 - t599;
t576 = -t443 * t535 - t598;
t574 = t557 * t398 + t560 * t400 + t433 * t648 - t441 * t649;
t532 = qJ(4) * t687;
t455 = t559 * t596 + t532;
t571 = t524 * t671 + t599;
t391 = t609 + t650;
t569 = -t561 * t601 + t579;
t566 = -t442 * t535 - t569;
t407 = (-t561 * t703 - t549) * t652 + t596 * t655 + t666;
t385 = -pkin(4) * t613 + (-pkin(4) * t572 - t632) * qJD(1) + t630;
t502 = pkin(5) - t589;
t501 = -qJ(6) + t663;
t475 = t497 * t559;
t474 = t557 * t687 - t559 * t689;
t468 = -t532 + (pkin(3) * t558 + pkin(7)) * t559;
t459 = t553 - t608;
t448 = qJD(1) * t594 - t694;
t447 = -pkin(7) * t625 + t669;
t427 = pkin(3) * t535 + t641;
t425 = pkin(5) * t497 + qJ(6) * t582 + t488;
t422 = -t460 - t695;
t418 = pkin(3) * t572 + qJ(4) * t620 + t632 - t666;
t416 = qJD(2) * t594 - t591;
t414 = -qJD(4) * t562 + (-t619 - t622) * pkin(7) + t629;
t409 = qJD(2) * t575 + (qJD(3) - qJD(5)) * t559 * t582;
t408 = t445 * t559 + t557 * t621 - t560 * t623;
t405 = pkin(5) * t474 - qJ(6) * t475 + t455;
t394 = t569 + t706;
t393 = -pkin(5) * t562 - t586;
t392 = qJ(6) * t562 + t707;
t377 = -t401 + t423;
t371 = -pkin(5) * t524 + t640;
t370 = pkin(5) * t408 - qJ(6) * t409 - qJD(6) * t475 + t407;
t369 = pkin(5) * t656 - t704;
t368 = -qJ(6) * t656 + qJD(6) * t562 + t574;
t367 = t391 * pkin(5) + t390 * qJ(6) - t439 * qJD(6) + t385;
t365 = -t597 + t600 + t647;
t1 = [(-t574 * t524 - t600 * t562 + t407 * t439 - t455 * t390 + t385 * t475 + t417 * t409 + (qJD(1) * t707 + t375) * t656) * MDP(28) + ((-t493 * t655 + (-t495 - t625) * t652) * t561 + ((t460 + t654) * t559 + (-t495 * t562 - t559 * t581) * qJD(2)) * t558) * MDP(12) + 0.2e1 * t614 * t717 + MDP(6) * t683 + (-t390 * t475 + t409 * t439) * MDP(22) + (t394 * t458 + t399 * t468 + t402 * t459 + t414 * t429 + t416 * t427 + t418 * t432) * MDP(21) + (t365 * t392 + t366 * t393 + t367 * t405 + t368 * t372 + t369 * t371 + t370 * t376) * MDP(32) + (t704 * t524 - t599 * t562 + t407 * t584 + t455 * t391 + t385 * t474 + t417 * t408 + (-qJD(1) * t586 - t374) * t656) * MDP(27) + (t402 * t562 + t416 * t535 + t418 * t493 + (t699 + (qJD(1) * t468 + t432) * t651) * t559 + ((-qJD(1) * t459 - t427) * t559 + (t468 * t603 + t696) * t558) * qJD(2)) * MDP(18) + (-t414 * t535 - t418 * t495 + t460 * t468 + (-t432 * t642 - t394) * t562 + (qJD(2) * t607 + t432 * t653 - t698) * t559) * MDP(20) + ((-pkin(7) * t619 + t668) * t535 + t579 * t562 + (-pkin(7) * t460 - t514 * t653) * t559 + (t590 * t562 + (-pkin(7) * t691 - qJD(1) * t665 - t443) * t559) * qJD(2)) * MDP(17) + (-t591 * t535 + (-pkin(7) * t692 + qJD(1) * t608 + t442) * t656 + t590 * t652 + ((t514 * t558 + (t493 + 0.2e1 * t627) * pkin(7)) * qJD(2) + t598) * t562) * MDP(16) - MDP(7) * t686 + (pkin(7) * t686 + t562 * t610) * MDP(10) + (-t460 * t687 + (-t620 + t621) * t495) * MDP(11) + (-pkin(7) * t683 + t559 * t610) * MDP(9) + (t535 * t620 + t460 * t562 + (t495 * t559 + (t645 - t684) * t561) * qJD(2)) * MDP(13) + (-t524 - t660) * t611 + (-t414 * t493 + t416 * t495 - t459 * t460 + (t427 * t655 + (-qJD(3) * t607 + t402) * t559) * t561 + ((-qJD(3) * t427 - t394) * t559 + (-t429 * t562 - t458 * t603) * qJD(2)) * t558) * MDP(19) + (-t390 * t562 + t409 * t524 + (-qJD(1) * t475 - t439) * t656) * MDP(24) + (t365 * t562 - t367 * t475 + t368 * t524 - t370 * t439 - t376 * t409 + t390 * t405 + (-qJD(1) * t392 - t372) * t656) * MDP(31) + (t606 * t618 + (-t493 * t559 + (-t645 + (t535 + t603) * t562) * t558) * qJD(2)) * MDP(14) - t606 * t612 + t716 * t718 + (t390 * t474 - t391 * t475 - t408 * t439 - t409 * t584) * MDP(23) + (-t365 * t474 + t366 * t475 - t368 * t584 + t369 * t439 + t371 * t409 - t372 * t408 - t390 * t393 - t391 * t392) * MDP(30) + (-t391 * t562 - t408 * t524 + (qJD(1) * t474 + t584) * t656) * MDP(25) + (-t366 * t562 + t367 * t474 - t369 * t524 + t370 * t584 + t376 * t408 + t391 * t405 + (qJD(1) * t393 + t371) * t656) * MDP(29); (MDP(9) * t559 * t565 + MDP(10) * t682) * pkin(1) + (-t385 * t582 - t488 * t390 + t709 * t524 + t673 * t439 - t675 * t417 + (-t375 + t658) * t661) * MDP(28) + (t365 * t457 - t366 * t583 + t367 * t425 - t371 * t678 - t372 * t679 + t376 * t681) * MDP(32) + (-t365 * t497 - t366 * t582 - t371 * t675 - t372 * t674 + t390 * t583 - t391 * t457 - t439 * t678 + t584 * t679) * MDP(30) + (t367 * t582 + t390 * t425 - t679 * t524 - t681 * t439 + t675 * t376 + (t372 - t658) * t661) * MDP(31) + (t390 * t582 - t439 * t675) * MDP(22) + (-t675 * t524 + (qJD(2) * t582 + t439) * t661) * MDP(24) + (t390 * t497 + t391 * t582 - t439 * t674 + t584 * t675) * MDP(23) + (t367 * t497 + t391 * t425 + t678 * t524 + t681 * t584 + t674 * t376 + (-t371 - t659) * t661) * MDP(29) - t682 * t717 + (t385 * t497 + t488 * t391 + t708 * t524 + t673 * t584 + t674 * t417 + (t374 - t659) * t661) * MDP(27) + (t447 * t493 - t448 * t495 + (-t427 * t660 + t394 + (t427 + (t495 - t625) * pkin(8)) * qJD(3)) * t561 + ((-qJD(2) * t581 - t460 + t654) * pkin(8) + t710) * t558) * MDP(19) + ((-t460 + t695) * t561 + ((-t495 - t657) * qJD(3) + (-t618 - t712) * qJD(1)) * t558) * MDP(12) + (pkin(2) * t460 - t479 * t535 + (-t635 + t693) * qJD(3) + (-t514 * t685 + (t443 - t633) * t559 + (t535 * t687 + t712) * pkin(7)) * qJD(1)) * MDP(17) + (t503 * t691 + (t634 + (t514 - t556) * t558) * qJD(3) + ((-pkin(2) * t651 - t442) * t559 + (qJD(2) * t593 - t514 * t562) * t558 + (t535 * t690 - t562 * t605) * pkin(7)) * qJD(1)) * MDP(16) + (-t698 - t448 * t535 + t670 * t493 + (t634 + (qJD(2) * t507 + t432) * t558) * qJD(3) + ((t507 * t651 + t427) * t559 + (-t696 + (t507 * t562 - t701) * qJD(2)) * t558) * qJD(1)) * MDP(18) + (-t699 + t447 * t535 + t460 * t507 - t670 * t495 + (-t432 * t561 + t635) * qJD(3) + (t432 * t685 + (-t429 + t633) * t559) * qJD(1)) * MDP(20) + (-t460 * t558 - t495 * t691) * MDP(11) + (t535 * t653 + (-t558 * t684 + t559 * t605) * qJD(1)) * MDP(14) + (-t535 * t651 + (t559 * t604 + t561 * t684) * qJD(1)) * MDP(13) + (t399 * t507 - t427 * t448 - t429 * t447 + t670 * t432 + (t394 * t561 + t402 * t558 + (t427 * t561 - t429 * t558) * qJD(3)) * pkin(8)) * MDP(21) + t524 * qJD(1) * t643 + t535 * qJD(1) * t644 + t565 * t716 + (-t674 * t524 + (qJD(2) * t497 - t584) * t661) * MDP(25); -t722 + (-t495 * t535 - t567) * MDP(14) + (-t567 * qJ(4) + pkin(3) * t460 + (t429 - t443) * t495 + (t427 - t641) * t493) * MDP(19) + (-t439 * t524 + t391) * MDP(25) + (t417 * t439 - t589 * t615 + t571) * MDP(27) + (t502 * t615 + t534 + t571 + t726) * MDP(29) + (-t390 * t502 - t391 * t501 + (-t372 - t671) * t439) * MDP(30) + (-t432 * t493 + t449 * t495 - 0.2e1 * t518 + 0.2e1 * t526 - t566) * MDP(20) + (t493 * t514 + t566) * MDP(17) + t495 * t493 * MDP(11) + qJD(1) * t612 - t727 * MDP(23) + (t377 * t439 + (qJ(6) - t501) * t615 + (-qJD(6) - t680) * t524 - t600) * MDP(31) + (-t423 * MDP(27) - t377 * MDP(29) + (-t371 + t680) * MDP(30) - t720) * t584 + (-t493 ^ 2 + t495 ^ 2) * MDP(12) + t422 * MDP(13) + (-t495 * t514 + t576) * MDP(16) + (t365 * t501 + t366 * t502 - t371 * t671 - t372 * t680 - t376 * t377) * MDP(32) + (-pkin(3) * t402 + qJ(4) * t394 - t427 * t443 + t429 * t641 - t432 * t449) * MDP(21) + (-t423 * t439 + t524 * t711 + t663 * t615 + t600) * MDP(28) + (-t432 * t495 - t449 * t493 + t576 + 0.2e1 * t602) * MDP(18); -MDP(18) * t615 + t422 * MDP(19) - t535 ^ 2 * MDP(20) + t710 * MDP(21) + (t493 * MDP(18) - MDP(20) * t495 + t432 * MDP(21) - t376 * MDP(32) - t439 * t636 + t584 * t637) * t495 + (MDP(32) * t725 + t637 * t615 - t636 * t719 + t730) * t560 + ((t439 * t535 - t391 + t650) * MDP(30) + (t371 * t524 + t365) * MDP(32) + t636 * t615 + t637 * t719) * t557; (-t493 * t649 - t495 * t648 - t609) * MDP(25) + t577 * MDP(27) + t578 * MDP(28) + (-0.2e1 * t534 + t577) * MDP(29) + (pkin(5) * t390 - qJ(6) * t391) * MDP(30) + (-t578 - 0.2e1 * t597 + 0.2e1 * t647) * MDP(31) + (-pkin(5) * t366 + qJ(6) * t365 - t371 * t375 + t372 * t640 - t376 * t401) * MDP(32) + (t524 * MDP(25) - t417 * MDP(27) - t376 * MDP(29) + (t372 - t375) * MDP(30) + t401 * MDP(31) + MDP(23) * t439) * t439 + (-t401 * MDP(29) + (t371 - t640) * MDP(30) + t720) * t584 + t722; (t439 * t584 + t615) * MDP(29) - t730 + (-t719 - t727) * MDP(31) + (-t725 + t726) * MDP(32);];
tauc  = t1;
