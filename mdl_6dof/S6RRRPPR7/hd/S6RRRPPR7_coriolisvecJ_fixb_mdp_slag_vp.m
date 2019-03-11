% Calculate minimal parameter regressor of Coriolis joint torque vector for
% S6RRRPPR7
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d3,d6,theta5]';
% MDP [32x1]
%   Minimal dynamic parameter vector (fixed base model)
%   see S6RRRPPR7_convert_par2_MPV_fixb.m
% 
% Output:
% tauc [6x1]
%   joint torques required to compensate Coriolis and centrifugal load

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 16:02
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauc = S6RRRPPR7_coriolisvecJ_fixb_mdp_slag_vp(qJ, qJD, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(10,1),zeros(32,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRPPR7_coriolisvecJ_fixb_mdp_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRRPPR7_coriolisvecJ_fixb_mdp_slag_vp: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RRRPPR7_coriolisvecJ_fixb_mdp_slag_vp: pkin has to be [10x1] (double)');
assert(isreal(MDP) && all(size(MDP) == [32 1]), ...
  'S6RRRPPR7_coriolisvecJ_fixb_mdp_slag_vp: MDP has to be [32x1] (double)'); 

%% Symbolic Calculation
% From coriolisvec_joint_fixb_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 16:01:23
% EndTime: 2019-03-09 16:01:40
% DurationCPUTime: 11.49s
% Computational Cost: add. (5815->614), mult. (13911->806), div. (0->0), fcn. (9343->8), ass. (0->249)
t628 = sin(qJ(3));
t631 = cos(qJ(3));
t696 = t631 * qJD(2);
t629 = sin(qJ(2));
t707 = qJD(1) * t629;
t570 = t628 * t707 - t696;
t683 = t631 * t707;
t705 = qJD(2) * t628;
t572 = t683 + t705;
t624 = sin(pkin(10));
t625 = cos(pkin(10));
t506 = t570 * t624 + t572 * t625;
t627 = sin(qJ(6));
t630 = cos(qJ(6));
t753 = -t625 * t570 + t572 * t624;
t444 = t506 * t630 - t627 * t753;
t694 = qJD(1) * qJD(2);
t676 = t629 * t694;
t758 = t506 * t627 + t630 * t753;
t773 = MDP(26) * t444 * t758 + (t444 ^ 2 - t758 ^ 2) * MDP(27) - MDP(30) * t676;
t632 = cos(qJ(2));
t706 = qJD(1) * t632;
t605 = -qJD(3) + t706;
t598 = qJD(6) + t605;
t760 = t444 * t598;
t700 = qJD(3) * t631;
t772 = t631 * t706 - t700;
t684 = t628 * t706;
t702 = qJD(3) * t628;
t771 = t684 - t702;
t590 = -qJD(2) * pkin(2) + pkin(7) * t707;
t493 = t570 * pkin(3) - t572 * qJ(4) + t590;
t461 = -pkin(4) * t570 + qJD(5) - t493;
t427 = pkin(5) * t753 + t461;
t675 = t632 * t694;
t701 = qJD(3) * t629;
t679 = t628 * t701;
t693 = qJD(2) * qJD(3);
t524 = qJD(1) * t679 + (-t675 - t693) * t631;
t677 = t629 * t700;
t703 = qJD(2) * t632;
t638 = t628 * t703 + t677;
t525 = qJD(1) * t638 + t628 * t693;
t458 = -t524 * t625 + t525 * t624;
t593 = t605 * qJD(4);
t602 = qJ(4) * t676;
t586 = -pkin(2) * t632 - pkin(8) * t629 - pkin(1);
t553 = t586 * qJD(1);
t659 = pkin(2) * t629 - pkin(8) * t632;
t578 = t659 * qJD(2);
t555 = qJD(1) * t578;
t615 = pkin(7) * t706;
t591 = qJD(2) * pkin(8) + t615;
t644 = t553 * t700 + t628 * t555 - t591 * t702;
t664 = pkin(7) * t676;
t637 = t631 * t664 - t644;
t434 = -t593 + t602 - t637;
t415 = qJ(5) * t525 + qJD(5) * t570 + t434;
t663 = -t553 * t702 + t631 * t555 - t591 * t700 + t628 * t664;
t750 = pkin(3) + pkin(4);
t417 = qJ(5) * t524 - qJD(5) * t572 - t676 * t750 - t663;
t671 = t415 * t624 - t625 * t417;
t391 = -pkin(5) * t676 - pkin(9) * t458 - t671;
t395 = t625 * t415 + t624 * t417;
t457 = -t524 * t624 - t625 * t525;
t393 = -pkin(9) * t457 + t395;
t672 = -t630 * t391 + t627 * t393;
t769 = t427 * t444 + t672;
t767 = t598 * t758;
t748 = pkin(8) - qJ(5);
t589 = t748 * t631;
t749 = pkin(7) * t628;
t685 = -pkin(3) - t749;
t666 = -pkin(4) + t685;
t728 = t631 * t632;
t577 = t659 * qJD(1);
t738 = t577 * t631;
t766 = -t738 + (-qJ(5) * t728 + t629 * t666) * qJD(1) - qJD(3) * t589 + qJD(5) * t628;
t699 = qJD(5) * t631;
t551 = t628 * t577;
t716 = qJ(4) * t707 + t551;
t730 = t629 * t631;
t731 = t628 * t632;
t765 = (-pkin(7) * t730 + qJ(5) * t731) * qJD(1) + t716 + t702 * t748 + t699;
t508 = t631 * t553 - t628 * t591;
t471 = qJ(5) * t572 + t508;
t448 = t605 * t750 + qJD(4) - t471;
t509 = t628 * t553 + t631 * t591;
t472 = qJ(5) * t570 + t509;
t595 = t605 * qJ(4);
t459 = t472 - t595;
t410 = t625 * t448 - t459 * t624;
t762 = pkin(9) * t506;
t403 = pkin(5) * t605 + t410 - t762;
t697 = qJD(6) * t630;
t687 = t627 * t391 + t630 * t393 + t403 * t697;
t764 = -t427 * t758 + t687;
t761 = pkin(9) * t753;
t718 = t772 * t624 - t771 * t625;
t564 = t624 * t628 + t625 * t631;
t639 = t632 * t564;
t717 = qJD(1) * t639 - t564 * qJD(3);
t759 = -t772 * qJ(4) + t628 * qJD(4) + t615;
t757 = -0.2e1 * t694;
t756 = MDP(4) * t629;
t622 = t629 ^ 2;
t755 = MDP(5) * (-t632 ^ 2 + t622);
t725 = -t624 * t765 + t625 * t766;
t724 = t624 * t766 + t625 * t765;
t665 = pkin(3) * t676;
t445 = -t663 - t665;
t489 = -t595 + t509;
t754 = t489 * t605 + t445;
t688 = t750 * t628;
t720 = -qJD(3) * t688 + t684 * t750 + t759;
t680 = t632 * t696;
t752 = -t679 + t680;
t695 = qJD(4) - t508;
t751 = t572 ^ 2;
t601 = t605 ^ 2;
t439 = t525 * pkin(3) + pkin(7) * t675 + t524 * qJ(4) - t572 * qJD(4);
t747 = t439 * t628;
t746 = t439 * t631;
t744 = t493 * t572;
t743 = t524 * t628;
t742 = t570 * t572;
t741 = t570 * t605;
t740 = t572 * t605;
t737 = t590 * t628;
t736 = t590 * t631;
t735 = t598 * t605;
t734 = t605 * t631;
t411 = t624 * t448 + t625 * t459;
t404 = t411 - t761;
t733 = t627 * t404;
t618 = t628 * qJ(4);
t732 = t628 * t629;
t634 = qJD(2) ^ 2;
t729 = t629 * t634;
t727 = t632 * t634;
t635 = qJD(1) ^ 2;
t726 = t632 * t635;
t608 = pkin(7) * t728;
t714 = qJD(3) * t608 + t586 * t702;
t431 = (-qJ(5) * t703 - t578) * t631 + (qJ(5) * t702 + qJD(2) * t666 - t699) * t629 + t714;
t704 = qJD(2) * t629;
t715 = t628 * t578 + t586 * t700;
t643 = qJ(4) * t704 - qJD(4) * t632 + t715;
t432 = (-pkin(7) * qJD(2) + qJ(5) * qJD(3)) * t730 + (qJD(5) * t629 + (-pkin(7) * qJD(3) + qJ(5) * qJD(2)) * t632) * t628 + t643;
t402 = t624 * t431 + t625 * t432;
t646 = t624 * t631 - t625 * t628;
t499 = t630 * t564 - t627 * t646;
t723 = -qJD(6) * t499 + t627 * t718 - t630 * t717;
t500 = -t564 * t627 - t630 * t646;
t722 = qJD(6) * t500 - t627 * t717 - t630 * t718;
t423 = t625 * t471 + t624 * t472;
t721 = -pkin(5) * t718 + t720;
t607 = pkin(7) * t731;
t621 = t632 * pkin(3);
t497 = pkin(4) * t632 + t607 + t621 + (-qJ(5) * t629 - t586) * t631;
t711 = t628 * t586 + t608;
t521 = -qJ(4) * t632 + t711;
t507 = qJ(5) * t732 + t521;
t436 = t624 * t497 + t625 * t507;
t719 = t771 * pkin(3) + t759;
t513 = t572 * pkin(3) + t570 * qJ(4);
t588 = t748 * t628;
t517 = t624 * t588 + t625 * t589;
t712 = qJ(4) * t680 + qJD(4) * t730;
t698 = qJD(6) * t627;
t692 = pkin(8) * t605 * t628;
t691 = pkin(8) * t734;
t581 = -t631 * pkin(3) - pkin(2) - t618;
t690 = pkin(8) * t704;
t689 = pkin(8) * t696;
t686 = -t627 * t457 + t630 * t458 - t697 * t753;
t678 = t632 * t702;
t673 = pkin(1) * t757;
t579 = -qJ(4) * t624 - t625 * t750;
t401 = t625 * t431 - t432 * t624;
t670 = t630 * t457 + t627 * t458;
t422 = -t471 * t624 + t625 * t472;
t435 = t625 * t497 - t507 * t624;
t516 = t625 * t588 - t589 * t624;
t669 = t586 * t631 - t607;
t560 = t631 * pkin(4) - t581;
t668 = t570 + t696;
t667 = -t572 + t705;
t662 = -pkin(7) - t688;
t478 = -pkin(4) * t572 - t513;
t660 = t685 * t629;
t658 = t578 * t631 - t714;
t486 = -pkin(9) * t564 + t517;
t657 = -pkin(5) * t707 - pkin(9) * t717 + qJD(6) * t486 + t725;
t485 = pkin(9) * t646 + t516;
t656 = -pkin(9) * t718 - qJD(6) * t485 + t724;
t388 = t627 * t403 + t630 * t404;
t655 = -t410 * t624 + t411 * t625;
t539 = t564 * t629;
t420 = pkin(5) * t632 - pkin(9) * t539 + t435;
t538 = t624 * t730 - t625 * t732;
t424 = -pkin(9) * t538 + t436;
t654 = t420 * t630 - t424 * t627;
t653 = t420 * t627 + t424 * t630;
t487 = pkin(3) * t605 + t695;
t652 = t487 * t631 - t489 * t628;
t651 = t506 * t624 - t625 * t753;
t476 = t630 * t538 + t539 * t627;
t477 = -t538 * t627 + t539 * t630;
t574 = -pkin(5) + t579;
t580 = qJ(4) * t625 - t624 * t750;
t650 = t574 * t630 - t580 * t627;
t649 = t574 * t627 + t580 * t630;
t648 = t624 * t630 + t625 * t627;
t647 = t624 * t627 - t625 * t630;
t645 = qJD(1) * t622 - t605 * t632;
t642 = -t509 * t605 + t663;
t641 = t598 * t648;
t640 = t598 * t647;
t398 = -t506 * t698 + t686;
t604 = qJ(4) * t730;
t520 = t629 * t662 + t604;
t421 = -pkin(4) * t525 - t439;
t399 = qJD(6) * t444 + t670;
t636 = -t508 * t605 + t637;
t449 = (-t631 * t750 - t618) * t701 + t662 * t703 + t712;
t537 = -t604 + (pkin(3) * t628 + pkin(7)) * t629;
t522 = t621 - t669;
t514 = pkin(5) * t564 + t560;
t512 = qJD(1) * t660 - t738;
t511 = -pkin(7) * t683 + t716;
t482 = qJD(2) * t639 + t646 * t701;
t481 = t624 * t752 - t625 * t638;
t475 = -t524 - t741;
t473 = pkin(5) * t538 + t520;
t470 = pkin(3) * t638 + pkin(7) * t703 + qJ(4) * t679 - t712;
t462 = qJD(2) * t660 - t658;
t456 = (-t629 * t696 - t678) * pkin(7) + t643;
t433 = -pkin(5) * t506 + t478;
t419 = pkin(5) * t481 + t449;
t409 = t423 + t762;
t408 = t422 - t761;
t407 = qJD(6) * t477 + t630 * t481 + t482 * t627;
t406 = -qJD(6) * t476 - t481 * t627 + t482 * t630;
t405 = pkin(5) * t457 + t421;
t397 = -pkin(9) * t481 + t402;
t396 = -pkin(5) * t704 - pkin(9) * t482 + t401;
t387 = t403 * t630 - t733;
t1 = [((-t570 * t631 - t572 * t628) * t703 + (t743 - t525 * t631 + (t570 * t628 - t572 * t631) * qJD(3)) * t629) * MDP(12) + ((-pkin(7) * t678 + t715) * t605 + t644 * t632 + (-pkin(7) * t524 - t590 * t702) * t629 + ((pkin(7) * t572 + t736) * t632 + (-pkin(7) * t734 - qJD(1) * t711 - t509) * t629) * qJD(2)) * MDP(17) - MDP(7) * t729 + (pkin(7) * t729 + t632 * t673) * MDP(10) + (-pkin(7) * t727 + t629 * t673) * MDP(9) + (-t456 * t570 + t462 * t572 - t521 * t525 - t522 * t524 + t652 * t703 + (-t434 * t628 + t445 * t631 + (-t487 * t628 - t489 * t631) * qJD(3)) * t629) * MDP(19) + (t605 * t679 + t524 * t632 + (t572 * t629 + t631 * t645) * qJD(2)) * MDP(13) + (-t395 * t632 - t402 * t605 + t421 * t539 + t449 * t506 + t458 * t520 + t461 * t482) * MDP(23) + (-t399 * t632 - t407 * t598) * MDP(29) + (t398 * t632 + t406 * t598) * MDP(28) + (-t658 * t605 - t663 * t632 + (pkin(7) * t525 + t590 * t700) * t629 + ((pkin(7) * t570 + t737) * t632 + (t669 * qJD(1) + t508 + (-t605 + t706) * t749) * t629) * qJD(2)) * MDP(16) + MDP(6) * t727 + t755 * t757 + (t401 * t605 + t421 * t538 + t449 * t753 + t457 * t520 + t461 * t481 - t632 * t671) * MDP(22) + (-t395 * t538 - t401 * t506 - t402 * t753 - t410 * t482 - t411 * t481 - t435 * t458 - t436 * t457 + t539 * t671) * MDP(24) + (t605 * t677 + t525 * t632 + (-t570 * t629 - t628 * t645) * qJD(2)) * MDP(14) + (-t456 * t605 - t470 * t572 + t524 * t537 + (-t493 * t696 - t434) * t632 + (t493 * t702 - t746 + (qJD(1) * t521 + t489) * qJD(2)) * t629) * MDP(20) + (t462 * t605 + t470 * t570 + t525 * t537 + (t493 * t705 + t445) * t632 + (t493 * t700 + t747 + (-qJD(1) * t522 - t487) * qJD(2)) * t629) * MDP(18) + (-(qJD(6) * t654 + t396 * t627 + t397 * t630) * t598 - (-t404 * t698 + t687) * t632 + t419 * t444 + t473 * t398 + t405 * t477 + t427 * t406) * MDP(32) + (t395 * t436 + t401 * t410 + t402 * t411 + t421 * t520 - t435 * t671 + t449 * t461) * MDP(25) + ((t396 * t630 - t397 * t627) * t598 - t672 * t632 + t419 * t758 + t473 * t399 + t405 * t476 + t427 * t407 + (-t388 * t632 - t598 * t653) * qJD(6)) * MDP(31) + ((-qJD(1) * t477 - t444) * MDP(28) + (qJD(1) * t476 + t758) * MDP(29) + (qJD(1) * t436 + t411) * MDP(23) + (-qJD(1) * t435 - t410) * MDP(22) + (qJD(1) * t653 + t388) * MDP(32) + (-qJD(1) * t654 - t387) * MDP(31) + (-t605 - t706) * MDP(15) + (-t598 - t706) * MDP(30)) * t704 + (-t398 * t476 - t399 * t477 - t406 * t758 - t407 * t444) * MDP(27) + (t434 * t521 + t439 * t537 + t445 * t522 + t456 * t489 + t462 * t487 + t470 * t493) * MDP(21) + (-t524 * t730 + t572 * t752) * MDP(11) + 0.2e1 * t675 * t756 + (t398 * t477 + t406 * t444) * MDP(26); ((-t524 + t741) * t631 + (-t525 + t740) * t628) * MDP(12) + (-t572 * t734 - t743) * MDP(11) + (pkin(2) * t524 - t551 * t605 + (-t692 + t736) * qJD(3) + (-t590 * t728 + (t509 - t689) * t629 + (t605 * t730 + t632 * t667) * pkin(7)) * qJD(1)) * MDP(17) + (t577 * t734 - pkin(2) * t525 + (t691 + t737) * qJD(3) + (-t508 * t629 + (-t590 * t632 - t690) * t628 + (t605 * t732 - t632 * t668) * pkin(7)) * qJD(1)) * MDP(16) + (t605 * t702 + (-t605 * t731 + t629 * t668) * qJD(1)) * MDP(14) + (-t605 * t700 + (t605 * t728 + t629 * t667) * qJD(1)) * MDP(13) + (t398 * t500 + t444 * t723) * MDP(26) + (t514 * t398 + t405 * t500 + t723 * t427 + t721 * t444) * MDP(32) + (-t398 * t499 - t399 * t500 - t444 * t722 - t723 * t758) * MDP(27) + (-t421 * t646 + t458 * t560 - t717 * t461 + t720 * t506 + t724 * t605) * MDP(23) + (-t395 * t564 + t410 * t717 + t411 * t718 - t457 * t517 - t458 * t516 + t506 * t725 - t646 * t671 + t724 * t753) * MDP(24) + (t395 * t517 - t410 * t725 - t411 * t724 + t421 * t560 + t461 * t720 - t516 * t671) * MDP(25) + (t421 * t564 + t457 * t560 - t718 * t461 - t725 * t605 + t720 * t753) * MDP(22) + (t439 * t581 - t487 * t512 - t489 * t511 - t719 * t493 + (qJD(3) * t652 + t434 * t631 + t445 * t628) * pkin(8)) * MDP(21) + (t514 * t399 + t405 * t499 + t722 * t427 + t721 * t758) * MDP(31) + (t511 * t570 - t512 * t572 + (t434 - t605 * t487 + (qJD(3) * t572 - t525) * pkin(8)) * t631 + ((qJD(3) * t570 - t524) * pkin(8) + t754) * t628) * MDP(19) + (-t746 - t512 * t605 + t525 * t581 - t719 * t570 + (t493 * t628 + t691) * qJD(3) + (t487 * t629 + (-t493 * t632 - t690) * t628) * qJD(1)) * MDP(18) + (-t747 + t511 * t605 + t524 * t581 + t719 * t572 + (-t493 * t631 + t692) * qJD(3) + (t493 * t728 + (-t489 + t689) * t629) * qJD(1)) * MDP(20) - t726 * t756 + t635 * t755 + ((t627 * t657 + t630 * t656) * MDP(32) + t723 * MDP(28) + (t627 * t656 - t630 * t657) * MDP(31) - t722 * MDP(29)) * t598 + (((t485 * t627 + t486 * t630) * qJD(2) - t388) * MDP(32) + (-qJD(2) * t500 + t444) * MDP(28) + (qJD(2) * t517 - t411) * MDP(23) + (-qJD(2) * t516 + t410) * MDP(22) + (-(t485 * t630 - t486 * t627) * qJD(2) + t387) * MDP(31) + (qJD(2) * t499 - t758) * MDP(29) + t605 * MDP(15) + t598 * MDP(30)) * t707 + (MDP(9) * t629 * t635 + MDP(10) * t726) * pkin(1); (-t513 * t570 + t642 + 0.2e1 * t665 - t744) * MDP(18) + (-t572 * t590 + t642) * MDP(16) + MDP(11) * t742 + (-t579 * t676 + t461 * t506 - t478 * t753 + (-qJD(4) * t624 - t422) * t605 + t671) * MDP(22) + (t580 * t676 - t461 * t753 - t478 * t506 + (-qJD(4) * t625 + t423) * t605 + t395) * MDP(23) + (t399 - t760) * MDP(29) + (-t650 * t676 - (t408 * t630 - t409 * t627) * t598 - t433 * t758 - qJD(4) * t641 + (-t598 * t649 + t388) * qJD(6) + t769) * MDP(31) + MDP(15) * t676 + (-t493 * t570 + t513 * t572 - 0.2e1 * t593 + 0.2e1 * t602 - t636) * MDP(20) + (t570 * t590 + t636) * MDP(17) + (-t570 ^ 2 + t751) * MDP(12) + (qJD(4) * t655 + t395 * t580 - t410 * t422 - t411 * t423 - t461 * t478 - t579 * t671) * MDP(25) + (-t398 - t767) * MDP(28) + (t649 * t676 + (t408 * t627 + t409 * t630) * t598 - t433 * t444 + qJD(4) * t640 + (-t598 * t650 - t733) * qJD(6) + t764) * MDP(32) + (qJD(4) * t651 - t457 * t580 - t458 * t579 + (t410 + t423) * t753 + (-t411 + t422) * t506) * MDP(24) + (pkin(3) * t524 - qJ(4) * t525 + (t489 - t509) * t572 + (t487 - t695) * t570) * MDP(19) + (-pkin(3) * t445 + qJ(4) * t434 - t487 * t509 + t489 * t695 - t493 * t513) * MDP(21) + (-t525 - t740) * MDP(14) + t475 * MDP(13) - t773; (-t676 + t742) * MDP(18) + t475 * MDP(19) + (-t601 - t751) * MDP(20) + (t744 + t754) * MDP(21) + (-t572 * t753 - t601 * t624 - t625 * t676) * MDP(22) + (-t506 * t572 - t601 * t625 + t624 * t676) * MDP(23) + (-t457 * t624 - t458 * t625 + t605 * t651) * MDP(24) + (t395 * t624 - t461 * t572 + t605 * t655 - t625 * t671) * MDP(25) + (-qJD(6) * t641 - t572 * t758 + t647 * t676 - t648 * t735) * MDP(31) + (qJD(6) * t640 - t572 * t444 + t647 * t735 + t648 * t676) * MDP(32); (t506 * t605 + t457) * MDP(22) + (-t605 * t753 + t458) * MDP(23) + (-t506 ^ 2 - t753 ^ 2) * MDP(24) + (t410 * t506 + t411 * t753 + t421) * MDP(25) + (t399 + t760) * MDP(31) + (t398 - t767) * MDP(32); (t686 + t767) * MDP(28) + (-t670 + t760) * MDP(29) + (t388 * t598 - t769) * MDP(31) + (t387 * t598 - t764) * MDP(32) + ((-MDP(29) * t506 - MDP(31) * t404) * t630 + (-MDP(28) * t506 + MDP(29) * t753 - MDP(31) * t403 + MDP(32) * t404) * t627) * qJD(6) + t773;];
tauc  = t1;
