% Calculate minimal parameter regressor of Coriolis joint torque vector for
% S6RRPRRR2
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [11x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d4,d5,d6,theta3]';
% MDP [33x1]
%   Minimal dynamic parameter vector (fixed base model)
%   see S6RRPRRR2_convert_par2_MPV_fixb.m
% 
% Output:
% tauc [6x1]
%   joint torques required to compensate Coriolis and centrifugal load

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 13:20
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauc = S6RRPRRR2_coriolisvecJ_fixb_mdp_slag_vp(qJ, qJD, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(11,1),zeros(33,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPRRR2_coriolisvecJ_fixb_mdp_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRPRRR2_coriolisvecJ_fixb_mdp_slag_vp: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6RRPRRR2_coriolisvecJ_fixb_mdp_slag_vp: pkin has to be [11x1] (double)');
assert(isreal(MDP) && all(size(MDP) == [33 1]), ...
  'S6RRPRRR2_coriolisvecJ_fixb_mdp_slag_vp: MDP has to be [33x1] (double)'); 

%% Symbolic Calculation
% From coriolisvec_joint_fixb_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 13:19:27
% EndTime: 2019-03-09 13:19:42
% DurationCPUTime: 10.14s
% Computational Cost: add. (9655->506), mult. (24687->664), div. (0->0), fcn. (19899->10), ass. (0->240)
t636 = sin(pkin(11));
t637 = cos(pkin(11));
t641 = sin(qJ(2));
t644 = cos(qJ(2));
t606 = -t636 * t641 + t637 * t644;
t596 = t606 * qJD(1);
t607 = t636 * t644 + t637 * t641;
t598 = t607 * qJD(1);
t640 = sin(qJ(4));
t750 = cos(qJ(4));
t550 = t750 * t596 - t598 * t640;
t754 = qJD(5) + qJD(6);
t786 = t550 - t754;
t638 = sin(qJ(6));
t639 = sin(qJ(5));
t642 = cos(qJ(6));
t643 = cos(qJ(5));
t611 = t638 * t639 - t642 * t643;
t783 = t786 * t611;
t692 = -qJD(5) + t550;
t539 = qJD(6) - t692;
t721 = t638 * t643;
t612 = t639 * t642 + t721;
t691 = qJD(1) * qJD(2);
t679 = t644 * t691;
t680 = t641 * t691;
t588 = -t636 * t680 + t637 * t679;
t597 = t607 * qJD(2);
t650 = qJD(1) * t597;
t657 = -t640 * t596 - t598 * t750;
t648 = qJD(4) * t657 - t640 * t588 - t750 * t650;
t775 = t783 * t539 - t612 * t648;
t785 = t786 * t612;
t696 = qJD(5) * t643;
t768 = t550 * t643;
t784 = t696 - t768;
t769 = t550 * t639;
t782 = pkin(10) * t769;
t502 = -pkin(4) * t657 - pkin(9) * t550;
t699 = qJD(1) * t641;
t570 = pkin(2) * t699 + pkin(3) * t598;
t479 = t502 + t570;
t744 = -qJ(3) - pkin(7);
t618 = t744 * t641;
t613 = qJD(1) * t618;
t619 = t744 * t644;
t614 = qJD(1) * t619;
t722 = t637 * t614;
t556 = -t613 * t636 + t722;
t748 = pkin(8) * t596;
t535 = t556 - t748;
t601 = t636 * t614;
t557 = t637 * t613 + t601;
t747 = pkin(8) * t598;
t536 = t557 - t747;
t627 = pkin(2) * t637 + pkin(3);
t749 = pkin(2) * t636;
t757 = t750 * t627 - t640 * t749;
t710 = -t757 * qJD(4) + t640 * t535 + t536 * t750;
t781 = -t643 * t479 + t639 * t710;
t697 = qJD(5) * t639;
t780 = (t697 - t769) * pkin(5);
t779 = -pkin(5) * t657 - pkin(10) * t768;
t678 = qJD(2) * t744;
t593 = qJD(3) * t644 + t641 * t678;
t576 = t593 * qJD(1);
t594 = -qJD(3) * t641 + t644 * t678;
t577 = t594 * qJD(1);
t530 = -t576 * t636 + t637 * t577;
t511 = -pkin(8) * t588 + t530;
t533 = t637 * t576 + t636 * t577;
t512 = -pkin(8) * t650 + t533;
t743 = qJD(2) * pkin(2);
t605 = t613 + t743;
t552 = t637 * t605 + t601;
t521 = qJD(2) * pkin(3) + t552 - t747;
t553 = t636 * t605 - t722;
t527 = t553 + t748;
t681 = qJD(4) * t750;
t698 = qJD(4) * t640;
t425 = -t750 * t511 + t640 * t512 + t521 * t698 + t527 * t681;
t633 = qJD(2) + qJD(4);
t534 = t633 * t639 - t643 * t657;
t497 = t750 * t588 + t596 * t681 - t598 * t698 - t640 * t650;
t738 = t497 * t639;
t453 = qJD(5) * t534 + t738;
t410 = pkin(5) * t453 + t425;
t464 = t521 * t750 - t640 * t527;
t460 = -t633 * pkin(4) - t464;
t531 = -t643 * t633 - t639 * t657;
t443 = t531 * pkin(5) + t460;
t778 = t410 * t612 + t783 * t443;
t777 = t410 * t611 - t785 * t443;
t452 = t643 * t497 + t633 * t696 + t657 * t697;
t776 = t452 * t643 - t639 * t453 - t784 * t531;
t666 = t785 * t539 + t611 * t648;
t493 = t639 * t648;
t706 = -t692 * t696 - t493;
t774 = t692 * t768 + t706;
t694 = qJD(6) * t642;
t689 = t642 * t452 - t638 * t453 - t531 * t694;
t695 = qJD(6) * t638;
t411 = -t534 * t695 + t689;
t659 = t531 * t638 - t642 * t534;
t673 = t452 * t638 + t642 * t453;
t412 = -qJD(6) * t659 + t673;
t450 = t452 * t639;
t733 = t534 * t638;
t473 = t642 * t531 + t733;
t773 = (t784 * t534 + t450) * MDP(20) + t411 * t612 * MDP(27) + (-t411 * t611 - t612 * t412 - t783 * t473) * MDP(28) + (-t783 * MDP(27) - t785 * MDP(28)) * t659;
t771 = t473 * t539;
t770 = t539 * t659;
t465 = t640 * t521 + t750 * t527;
t461 = pkin(9) * t633 + t465;
t687 = -pkin(2) * t644 - pkin(1);
t665 = t687 * qJD(1);
t617 = qJD(3) + t665;
t560 = -pkin(3) * t596 + t617;
t469 = -pkin(4) * t550 + pkin(9) * t657 + t560;
t430 = t461 * t643 + t469 * t639;
t421 = -pkin(10) * t531 + t430;
t418 = t421 * t695;
t764 = t443 * t473 + t418;
t625 = pkin(2) * t680;
t561 = pkin(3) * t650 + t625;
t437 = -pkin(4) * t648 - t497 * pkin(9) + t561;
t436 = t643 * t437;
t647 = -t640 * t511 - t512 * t750 - t521 * t681 + t527 * t698;
t649 = -qJD(5) * t430 + t639 * t647 + t436;
t399 = -pkin(5) * t648 - pkin(10) * t452 + t649;
t653 = t639 * t437 - t461 * t697 + t469 * t696 - t643 * t647;
t400 = -pkin(10) * t453 + t653;
t675 = t642 * t399 - t638 * t400;
t763 = t443 * t659 + t675;
t762 = -t648 * MDP(31) + (-t473 ^ 2 + t659 ^ 2) * MDP(28) - t473 * MDP(27) * t659;
t761 = -0.2e1 * t691;
t760 = MDP(4) * t641;
t759 = MDP(5) * (t641 ^ 2 - t644 ^ 2);
t555 = t640 * t606 + t607 * t750;
t503 = t612 * t555;
t651 = t640 * t627 + t749 * t750;
t704 = t651 * qJD(4) + t535 * t750 - t640 * t536;
t495 = t643 * t648;
t758 = -t692 * t697 + t495;
t562 = t637 * t618 + t619 * t636;
t542 = -pkin(8) * t607 + t562;
t563 = t636 * t618 - t637 * t619;
t543 = pkin(8) * t606 + t563;
t756 = t750 * t542 - t640 * t543;
t755 = t639 * t479 + t643 * t710;
t753 = MDP(13) * t657 - MDP(14) * t550 - t560 * MDP(19);
t429 = -t461 * t639 + t643 * t469;
t420 = -pkin(10) * t534 + t429;
t414 = -pkin(5) * t692 + t420;
t742 = t414 * t642;
t403 = -t421 * t638 + t742;
t741 = t421 * t642;
t404 = t414 * t638 + t741;
t752 = -MDP(14) * t657 - MDP(18) * t560 + MDP(24) * t692 - MDP(25) * t429 + MDP(26) * t430 - MDP(31) * t539 - MDP(32) * t403 + MDP(33) * t404;
t751 = -pkin(9) - pkin(10);
t746 = t643 * pkin(5);
t592 = pkin(9) + t651;
t745 = -pkin(10) - t592;
t740 = t473 * t657;
t739 = t659 * t657;
t600 = t606 * qJD(2);
t656 = t606 * t750 - t640 * t607;
t505 = qJD(4) * t656 - t640 * t597 + t600 * t750;
t737 = t505 * t639;
t736 = t505 * t643;
t735 = t531 * t657;
t734 = t534 * t657;
t732 = t534 * t639;
t731 = t657 * t633;
t728 = t550 * t633;
t725 = t555 * t639;
t724 = t555 * t643;
t645 = qJD(2) ^ 2;
t720 = t641 * t645;
t492 = t640 * t542 + t543 * t750;
t485 = t643 * t492;
t719 = t644 * t645;
t646 = qJD(1) ^ 2;
t718 = t644 * t646;
t459 = t460 * t696;
t714 = t425 * t639 + t459;
t712 = t643 * t464 + t639 * t502;
t579 = -pkin(3) * t606 + t687;
t486 = -pkin(4) * t656 - pkin(9) * t555 + t579;
t707 = t639 * t486 + t485;
t705 = t704 + t780;
t541 = t637 * t593 + t636 * t594;
t631 = t641 * t743;
t686 = qJD(5) * t751;
t683 = t555 * t697;
t571 = pkin(3) * t597 + t631;
t677 = qJD(5) * t745;
t676 = pkin(1) * t761;
t674 = -t425 * t643 + t460 * t697;
t672 = -t464 * t639 + t643 * t502;
t540 = -t593 * t636 + t637 * t594;
t670 = t692 * t639;
t669 = qJD(6) * t414 + t400;
t667 = -t465 + t780;
t632 = t643 * pkin(10);
t565 = t592 * t643 + t632;
t664 = qJD(6) * t565 - t643 * t677 + t779 - t781;
t621 = pkin(9) * t643 + t632;
t663 = qJD(6) * t621 - t643 * t686 + t672 + t779;
t564 = t745 * t639;
t662 = -qJD(6) * t564 - t639 * t677 + t755 - t782;
t620 = t751 * t639;
t661 = -qJD(6) * t620 - t639 * t686 + t712 - t782;
t660 = -t460 * t550 + t592 * t648;
t591 = -pkin(4) - t757;
t658 = -t692 * t769 - t758;
t655 = t555 * t696 + t737;
t654 = -t683 + t736;
t514 = -pkin(8) * t600 + t540;
t515 = -pkin(8) * t597 + t541;
t433 = qJD(4) * t756 + t640 * t514 + t750 * t515;
t506 = qJD(4) * t555 + t597 * t750 + t640 * t600;
t442 = pkin(4) * t506 - pkin(9) * t505 + t571;
t652 = t643 * t433 + t639 * t442 + t486 * t696 - t492 * t697;
t434 = qJD(4) * t492 - t514 * t750 + t640 * t515;
t629 = -pkin(4) - t746;
t572 = t591 - t746;
t504 = t611 * t555;
t484 = t643 * t486;
t462 = t648 * t656;
t457 = pkin(5) * t725 - t756;
t441 = t643 * t442;
t431 = -pkin(10) * t725 + t707;
t427 = -pkin(5) * t656 - pkin(10) * t724 - t492 * t639 + t484;
t416 = t505 * t721 - t638 * t683 - t695 * t725 + (t724 * t754 + t737) * t642;
t415 = -t503 * t754 - t611 * t505;
t413 = pkin(5) * t655 + t434;
t402 = -pkin(10) * t655 + t652;
t401 = -pkin(10) * t736 + pkin(5) * t506 - t433 * t639 + t441 + (-t485 + (pkin(10) * t555 - t486) * t639) * qJD(5);
t1 = [((-t531 * t643 - t732) * t505 + (-t450 - t453 * t643 + (t531 * t639 - t534 * t643) * qJD(5)) * t555) * MDP(21) + (t452 * t724 + t534 * t654) * MDP(20) - MDP(7) * t720 + (pkin(7) * t720 + t644 * t676) * MDP(10) + (-pkin(7) * t719 + t641 * t676) * MDP(9) + (t530 * t562 + t533 * t563 + t552 * t540 + t553 * t541 + (t617 + t665) * t631) * MDP(12) + 0.2e1 * t679 * t760 + (-t411 * t503 + t412 * t504 - t415 * t473 + t416 * t659) * MDP(28) + (-t411 * t504 - t415 * t659) * MDP(27) + (t497 * t555 - t505 * t657) * MDP(13) + (t497 * t579 + t505 * t560 + t555 * t561 - t571 * t657) * MDP(19) + (MDP(15) * t505 - MDP(16) * t506 - MDP(18) * t434 - MDP(19) * t433) * t633 + (-t530 * t607 + t533 * t606 - t540 * t598 + t541 * t596 - t552 * t600 - t553 * t597 - t562 * t588 - t563 * t650) * MDP(11) + (-t506 * t692 + t462) * MDP(24) + (t453 * t656 + t493 * t555 - t506 * t531 + t655 * t692) * MDP(23) + (-t452 * t656 - t495 * t555 + t506 * t534 - t654 * t692) * MDP(22) + (t506 * t539 + t462) * MDP(31) + t759 * t761 + ((t401 * t642 - t402 * t638) * t539 - (t427 * t642 - t431 * t638) * t648 - t675 * t656 + t403 * t506 + t413 * t473 + t457 * t412 + t410 * t503 + t443 * t416 + ((-t427 * t638 - t431 * t642) * t539 + t404 * t656) * qJD(6)) * MDP(32) + (-(-t492 * t696 + t441) * t692 - t484 * t648 - (-t461 * t696 + t436) * t656 + t429 * t506 + t434 * t531 - t756 * t453 + t555 * t459 + (-(-qJD(5) * t486 - t433) * t692 + t492 * t648 - (-qJD(5) * t469 + t647) * t656 + t425 * t555 + t460 * t505) * t639) * MDP(25) + (t425 * t724 - t430 * t506 + t434 * t534 - t452 * t756 + t460 * t654 + t648 * t707 + t652 * t692 + t653 * t656) * MDP(26) + (-t404 * t506 - t410 * t504 + t457 * t411 - t413 * t659 + t443 * t415 - t418 * t656 + (-(-qJD(6) * t431 + t401) * t539 + t427 * t648 + t399 * t656) * t638 + (-(qJD(6) * t427 + t402) * t539 + t431 * t648 + t669 * t656) * t642) * MDP(33) + (-t411 * t656 + t415 * t539 + t504 * t648 - t506 * t659) * MDP(29) + (t506 * t560 - t550 * t571 - t561 * t656 - t579 * t648) * MDP(18) + (t497 * t656 + t505 * t550 + t506 * t657 + t555 * t648) * MDP(14) + (t412 * t656 - t416 * t539 - t473 * t506 + t503 * t648) * MDP(30) + MDP(6) * t719; t773 + (t497 - t728) * MDP(15) + (-t552 * t556 - t553 * t557 + (t530 * t637 + t533 * t636 - t617 * t699) * pkin(2)) * MDP(12) - t718 * t760 + t753 * t550 + (t591 * t452 + t660 * t643 + t704 * t534 - (t592 * t697 + t755) * t692 + t714) * MDP(26) - t752 * t657 + (t550 * t570 - t633 * t704 - t425) * MDP(18) + (t570 * t657 + t633 * t710 + t647) * MDP(19) + (MDP(9) * t641 * t646 + MDP(10) * t718) * pkin(1) + ((t553 + t556) * t598 + (-t557 + t552) * t596 + (-t637 * t588 - t636 * t650) * pkin(2)) * MDP(11) + t646 * t759 + (-t739 + t775) * MDP(29) + ((t564 * t638 + t565 * t642) * t648 + t572 * t411 + (t638 * t664 + t642 * t662) * t539 - t705 * t659 + t778) * MDP(33) + (-(t564 * t642 - t565 * t638) * t648 + t572 * t412 + (t638 * t662 - t642 * t664) * t539 + t705 * t473 + t777) * MDP(32) + (t692 * t732 + t776) * MDP(21) + (t734 + t774) * MDP(22) + (t591 * t453 + t660 * t639 + t704 * t531 - (-t592 * t696 + t781) * t692 + t674) * MDP(25) + (t648 - t731) * MDP(16) + (t658 - t735) * MDP(23) + (t666 - t740) * MDP(30); (-t596 ^ 2 - t598 ^ 2) * MDP(11) + (t552 * t598 - t553 * t596 + t625) * MDP(12) + (-t648 - t731) * MDP(18) + (t497 + t728) * MDP(19) + (t658 + t735) * MDP(25) + (-t643 * t692 ^ 2 + t493 + t734) * MDP(26) + (t666 + t740) * MDP(32) + (-t739 - t775) * MDP(33); t497 * MDP(15) + t648 * MDP(16) + (t465 * t633 - t425) * MDP(18) + (t464 * t633 + t647) * MDP(19) + (t534 * t670 + t776) * MDP(21) + t774 * MDP(22) + (-t670 * t692 - t495) * MDP(23) + (-pkin(4) * t453 - pkin(9) * t706 - t460 * t769 - t465 * t531 + t672 * t692 + t674) * MDP(25) + (-pkin(4) * t452 + pkin(9) * t758 - t460 * t768 - t465 * t534 - t692 * t712 + t714) * MDP(26) + t775 * MDP(29) + t666 * MDP(30) + (-(t620 * t642 - t621 * t638) * t648 + t629 * t412 + (t638 * t661 - t642 * t663) * t539 + t667 * t473 + t777) * MDP(32) + ((t620 * t638 + t621 * t642) * t648 + t629 * t411 + (t638 * t663 + t642 * t661) * t539 - t667 * t659 + t778) * MDP(33) - (t633 * MDP(15) - t753) * t550 - (MDP(16) * t633 - MDP(22) * t534 + MDP(23) * t531 + MDP(29) * t659 + MDP(30) * t473 + t752) * t657 + t773; t534 * t531 * MDP(20) + (-t531 ^ 2 + t534 ^ 2) * MDP(21) + (-t531 * t692 + t452) * MDP(22) + (-t738 + (-qJD(5) - t692) * t534) * MDP(23) - t648 * MDP(24) + (-t430 * t692 - t460 * t534 + t649) * MDP(25) + (-t429 * t692 + t460 * t531 - t653) * MDP(26) + (t411 + t771) * MDP(29) + (-t412 - t770) * MDP(30) + (-(-t420 * t638 - t741) * t539 - t404 * qJD(6) + (-t473 * t534 - t539 * t695 - t642 * t648) * pkin(5) + t763) * MDP(32) + ((-t421 * t539 - t399) * t638 + (t420 * t539 - t669) * t642 + (t534 * t659 - t539 * t694 + t638 * t648) * pkin(5) + t764) * MDP(33) + t762; (t689 + t771) * MDP(29) + (-t673 - t770) * MDP(30) + (t404 * t539 + t763) * MDP(32) + (-t638 * t399 - t642 * t400 + t403 * t539 + t764) * MDP(33) + (-MDP(29) * t733 + MDP(30) * t659 - MDP(32) * t404 - MDP(33) * t742) * qJD(6) + t762;];
tauc  = t1;
