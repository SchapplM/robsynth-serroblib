% Calculate vector of inverse dynamics joint torques for
% S6RPRRPR6
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
%   pkin=[a2,a3,a4,a5,a6,d1,d3,d4,d6,theta2,theta5]';
% MDP [30x1]
%   Minimal dynamic parameter vector (fixed base model)
%   see S6RPRRPR6_convert_par2_MPV_fixb.m
% 
% Output:
% tau [6x1]
%   joint torques of inverse dynamics (contains inertial, gravitational Coriolis and centrifugal forces)

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 05:18
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tau = S6RPRRPR6_invdynJ_fixb_mdp_slag_vp(qJ, qJD, qJDD, g, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(3,1),zeros(11,1),zeros(30,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRRPR6_invdynJ_fixb_mdp_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RPRRPR6_invdynJ_fixb_mdp_slag_vp: qJD has to be [6x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [6 1]), ...
  'S6RPRRPR6_invdynJ_fixb_mdp_slag_vp: qJDD has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RPRRPR6_invdynJ_fixb_mdp_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6RPRRPR6_invdynJ_fixb_mdp_slag_vp: pkin has to be [11x1] (double)');
assert(isreal(MDP) && all(size(MDP) == [30 1]), ...
  'S6RPRRPR6_invdynJ_fixb_mdp_slag_vp: MDP has to be [30x1] (double)'); 

%% Symbolic Calculation
% From invdyn_joint_fixb_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 05:18:01
% EndTime: 2019-03-09 05:18:15
% DurationCPUTime: 11.39s
% Computational Cost: add. (9477->555), mult. (22728->717), div. (0->0), fcn. (18095->16), ass. (0->250)
t644 = cos(pkin(10));
t653 = cos(qJ(3));
t717 = qJD(1) * t653;
t622 = t644 * t717;
t642 = sin(pkin(10));
t649 = sin(qJ(3));
t718 = qJD(1) * t649;
t700 = t642 * t718;
t591 = t622 - t700;
t578 = qJD(4) - t591;
t651 = cos(qJ(6));
t602 = t642 * t653 + t644 * t649;
t593 = t602 * qJD(1);
t648 = sin(qJ(4));
t652 = cos(qJ(4));
t710 = t652 * qJD(3);
t555 = t593 * t648 - t710;
t557 = qJD(3) * t648 + t593 * t652;
t641 = sin(pkin(11));
t643 = cos(pkin(11));
t680 = -t555 * t643 - t557 * t641;
t731 = t651 * t680;
t505 = t555 * t641 - t557 * t643;
t647 = sin(qJ(6));
t756 = t505 * t647;
t457 = t731 + t756;
t574 = qJD(6) + t578;
t759 = t457 * t574;
t627 = pkin(2) * t644 + pkin(1);
t611 = -qJD(1) * t627 + qJD(2);
t513 = -pkin(3) * t591 - pkin(8) * t593 + t611;
t764 = pkin(7) + qJ(2);
t612 = t764 * t642;
t606 = qJD(1) * t612;
t613 = t764 * t644;
t607 = qJD(1) * t613;
t549 = -t606 * t649 + t607 * t653;
t540 = qJD(3) * pkin(8) + t549;
t480 = t513 * t648 + t540 * t652;
t707 = qJDD(1) * t653;
t708 = qJDD(1) * t649;
t701 = qJD(3) * t622 + t642 * t707 + t644 * t708;
t544 = -qJD(3) * t700 + t701;
t595 = t602 * qJD(3);
t619 = t644 * t707;
t682 = -t642 * t708 + t619;
t545 = qJD(1) * t595 - t682;
t610 = -qJDD(1) * t627 + qJDD(2);
t493 = pkin(3) * t545 - pkin(8) * t544 + t610;
t486 = t652 * t493;
t709 = qJD(1) * qJD(2);
t773 = qJDD(1) * t764 + t709;
t571 = t773 * t642;
t572 = t773 * t644;
t679 = -t571 * t649 + t572 * t653;
t780 = -t606 * t653 - t607 * t649;
t491 = qJDD(3) * pkin(8) + qJD(3) * t780 + t679;
t713 = qJD(4) * t648;
t497 = qJD(4) * t710 + qJDD(3) * t648 + t544 * t652 - t593 * t713;
t538 = qJDD(4) + t545;
t415 = pkin(4) * t538 - qJ(5) * t497 - qJD(4) * t480 - qJD(5) * t557 - t491 * t648 + t486;
t498 = qJD(4) * t557 - qJDD(3) * t652 + t544 * t648;
t712 = qJD(4) * t652;
t663 = t491 * t652 + t493 * t648 + t513 * t712 - t540 * t713;
t417 = -qJ(5) * t498 - qJD(5) * t555 + t663;
t403 = t415 * t643 - t417 * t641;
t449 = t497 * t643 - t498 * t641;
t401 = pkin(5) * t538 - pkin(9) * t449 + t403;
t404 = t415 * t641 + t417 * t643;
t448 = -t497 * t641 - t498 * t643;
t402 = pkin(9) * t448 + t404;
t479 = t513 * t652 - t540 * t648;
t465 = -qJ(5) * t557 + t479;
t452 = pkin(4) * t578 + t465;
t466 = -qJ(5) * t555 + t480;
t739 = t643 * t466;
t429 = t452 * t641 + t739;
t784 = pkin(9) * t680;
t420 = t429 + t784;
t711 = qJD(6) * t647;
t419 = t420 * t711;
t539 = -qJD(3) * pkin(3) - t780;
t499 = pkin(4) * t555 + qJD(5) + t539;
t453 = -pkin(5) * t680 + t499;
t635 = qJ(4) + pkin(11) + qJ(6);
t624 = sin(t635);
t625 = cos(t635);
t654 = cos(qJ(1));
t640 = pkin(10) + qJ(3);
t632 = cos(t640);
t650 = sin(qJ(1));
t743 = t632 * t650;
t559 = t624 * t654 - t625 * t743;
t742 = t632 * t654;
t561 = t624 * t650 + t625 * t742;
t631 = sin(t640);
t766 = g(3) * t631;
t798 = g(1) * t561 - g(2) * t559 - t647 * t401 - t651 * t402 - t453 * t457 + t625 * t766 + t419;
t535 = qJDD(6) + t538;
t785 = -t505 * t651 + t647 * t680;
t797 = t535 * MDP(28) + (-t457 ^ 2 + t785 ^ 2) * MDP(25) - t457 * MDP(24) * t785;
t601 = t641 * t652 + t643 * t648;
t789 = t578 * t601;
t674 = t641 * t648 - t643 * t652;
t788 = t578 * t674;
t761 = t785 * t574;
t541 = pkin(3) * t593 - pkin(8) * t591;
t530 = t652 * t541;
t763 = qJ(5) + pkin(8);
t696 = qJD(4) * t763;
t795 = -pkin(4) * t593 - t530 + (qJ(5) * t591 - t696) * t652 + (-qJD(5) + t780) * t648;
t723 = t541 * t648 + t652 * t780;
t749 = t591 * t648;
t794 = -qJ(5) * t749 - qJD(5) * t652 + t648 * t696 + t723;
t715 = qJD(3) * t653;
t716 = qJD(3) * t649;
t671 = -t571 * t653 - t572 * t649 + t606 * t716 - t607 * t715;
t492 = -qJDD(3) * pkin(3) - t671;
t689 = g(1) * t654 + g(2) * t650;
t661 = -g(3) * t632 + t631 * t689;
t793 = -pkin(8) * qJD(4) * t578 - t492 + t661;
t792 = t713 - t749;
t558 = t624 * t743 + t625 * t654;
t560 = -t624 * t742 + t625 * t650;
t695 = t401 * t651 - t647 * t402;
t791 = -g(1) * t560 + g(2) * t558 - t453 * t785 + t624 * t766 + t695;
t790 = pkin(9) * t505;
t547 = t601 * t651 - t647 * t674;
t728 = qJD(6) * t547 - t647 * t788 + t651 * t789;
t600 = t642 * t649 - t644 * t653;
t594 = t600 * qJD(3);
t699 = t602 * t712;
t787 = -t594 * t648 + t699;
t779 = g(1) * t650 - g(2) * t654;
t786 = qJDD(2) - t779;
t783 = t631 * t779;
t725 = t641 * t794 + t643 * t795;
t724 = t641 * t795 - t643 * t794;
t551 = t612 * t653 + t613 * t649;
t778 = pkin(4) * t792 - t549;
t730 = t652 * t654;
t735 = t648 * t650;
t579 = t632 * t735 + t730;
t733 = t650 * t652;
t734 = t648 * t654;
t581 = -t632 * t734 + t733;
t777 = -g(1) * t581 + g(2) * t579;
t776 = qJ(2) * qJDD(1);
t678 = -t601 * t647 - t651 * t674;
t729 = qJD(6) * t678 - t647 * t789 - t651 * t788;
t775 = -t535 * t547 - t574 * t729;
t774 = t689 * t632 + t766;
t694 = -t448 * t651 + t449 * t647;
t412 = qJD(6) * t785 + t694;
t772 = pkin(4) * t641;
t762 = qJDD(1) * pkin(1);
t760 = t457 * t593;
t758 = t785 * t593;
t757 = t497 * t648;
t754 = t538 * t648;
t753 = t555 * t578;
t752 = t555 * t593;
t751 = t557 * t578;
t750 = t557 * t593;
t747 = t602 * t648;
t746 = t602 * t652;
t459 = t641 * t466;
t737 = t644 * MDP(4);
t428 = t452 * t643 - t459;
t418 = pkin(5) * t578 + t428 + t790;
t732 = t651 * t418;
t527 = t652 * t538;
t552 = -t612 * t649 + t613 * t653;
t550 = t652 * t552;
t514 = -qJD(2) * t600 - qJD(3) * t551;
t542 = pkin(3) * t595 + pkin(8) * t594;
t531 = t652 * t542;
t543 = pkin(3) * t600 - pkin(8) * t602 - t627;
t673 = qJ(5) * t594 - qJD(5) * t602;
t436 = pkin(4) * t595 - t514 * t648 + t531 + t673 * t652 + (-t550 + (qJ(5) * t602 - t543) * t648) * qJD(4);
t703 = t514 * t652 + t542 * t648 + t543 * t712;
t440 = -qJ(5) * t699 + (-qJD(4) * t552 + t673) * t648 + t703;
t410 = t436 * t641 + t440 * t643;
t433 = t465 * t643 - t459;
t533 = t652 * t543;
t470 = pkin(4) * t600 - qJ(5) * t746 - t552 * t648 + t533;
t722 = t543 * t648 + t550;
t481 = -qJ(5) * t747 + t722;
t442 = t470 * t641 + t481 * t643;
t721 = pkin(5) * t789 + t778;
t614 = t763 * t648;
t615 = t763 * t652;
t554 = -t614 * t641 + t615 * t643;
t720 = t642 ^ 2 + t644 ^ 2;
t714 = qJD(4) * t602;
t704 = qJD(6) * t731 + t448 * t647 + t449 * t651;
t629 = pkin(4) * t652 + pkin(3);
t698 = pkin(4) * t648 + t764;
t409 = t436 * t643 - t440 * t641;
t432 = -t465 * t641 - t739;
t441 = t470 * t643 - t481 * t641;
t553 = -t614 * t643 - t615 * t641;
t692 = t578 * t652;
t691 = -qJD(4) * t513 - t491;
t515 = qJD(2) * t602 - t612 * t716 + t613 * t715;
t690 = 0.2e1 * t720;
t687 = t535 * t678 - t574 * t728;
t686 = pkin(4) * t747 + t551;
t685 = -t540 * t712 + t486;
t522 = -pkin(9) * t674 + t554;
t684 = pkin(5) * t593 - pkin(9) * t788 + qJD(6) * t522 - t725;
t521 = -pkin(9) * t601 + t553;
t683 = pkin(9) * t789 - qJD(6) * t521 - t724;
t406 = t647 * t418 + t651 * t420;
t524 = t601 * t602;
t525 = t674 * t602;
t681 = -t524 * t651 + t525 * t647;
t488 = -t524 * t647 - t525 * t651;
t676 = t629 * t632 + t631 * t763;
t672 = t762 - t786;
t670 = -t578 * t792 + t527;
t626 = pkin(4) * t643 + pkin(5);
t669 = t626 * t647 + t651 * t772;
t668 = t626 * t651 - t647 * t772;
t667 = pkin(4) * t787 + t515;
t666 = t627 + t676;
t665 = -t594 * t652 - t602 * t713;
t664 = -pkin(8) * t538 + t539 * t578;
t411 = t505 * t711 + t704;
t658 = t690 * t709 - t689;
t444 = pkin(4) * t498 + qJDD(5) + t492;
t582 = t632 * t730 + t735;
t580 = -t632 * t733 + t734;
t565 = pkin(5) * t674 - t629;
t490 = pkin(4) * t557 - pkin(5) * t505;
t489 = pkin(5) * t524 + t686;
t484 = -t594 * t674 + t601 * t714;
t483 = t594 * t601 + t674 * t714;
t443 = -pkin(5) * t483 + t667;
t431 = -pkin(9) * t524 + t442;
t430 = pkin(5) * t600 + pkin(9) * t525 + t441;
t426 = qJD(6) * t488 - t483 * t651 - t484 * t647;
t425 = qJD(6) * t681 + t483 * t647 - t484 * t651;
t423 = t433 + t790;
t422 = t432 - t784;
t421 = -pkin(5) * t448 + t444;
t408 = pkin(9) * t483 + t410;
t407 = pkin(5) * t595 + pkin(9) * t484 + t409;
t405 = -t420 * t647 + t732;
t1 = [(-t412 * t600 - t426 * t574 + t457 * t595 + t535 * t681) * MDP(27) + ((t407 * t651 - t408 * t647) * t574 + (t430 * t651 - t431 * t647) * t535 + t695 * t600 + t405 * t595 - t443 * t457 + t489 * t412 - t421 * t681 + t453 * t426 - g(1) * t559 - g(2) * t561 + ((-t430 * t647 - t431 * t651) * t574 - t406 * t600) * qJD(6)) * MDP(29) + (t411 * t681 - t412 * t488 + t425 * t457 - t426 * t785) * MDP(25) + (-MDP(5) * t642 + t737) * (t672 + t762) + (-t544 * t600 - t545 * t602 - t591 * t594 - t593 * t595) * MDP(9) + (-qJD(3) * t594 + qJDD(3) * t602) * MDP(10) + (t544 * t602 - t593 * t594) * MDP(8) + (-(-t555 * t652 - t557 * t648) * t594 + (-t757 - t498 * t652 + (t555 * t648 - t557 * t652) * qJD(4)) * t602) * MDP(16) + ((-t552 * t712 + t531) * t578 + t533 * t538 + t685 * t600 + t479 * t595 + t515 * t555 + t551 * t498 + t539 * t699 - g(1) * t580 - g(2) * t582 + ((-qJD(4) * t543 - t514) * t578 - t552 * t538 + t691 * t600 + t492 * t602 - t539 * t594) * t648) * MDP(20) + (-t498 * t600 - t538 * t747 - t555 * t595 - t578 * t787) * MDP(18) + (-qJD(3) * t515 - qJDD(3) * t551 - t545 * t627 + t595 * t611 + t600 * t610 + t632 * t779) * MDP(13) + t779 * MDP(2) + (t403 * t525 - t404 * t524 + t409 * t505 + t410 * t680 + t428 * t484 + t429 * t483 - t441 * t449 + t442 * t448 + t783) * MDP(22) + qJDD(1) * MDP(1) + (-qJD(3) * t514 - qJDD(3) * t552 - t544 * t627 - t594 * t611 + t602 * t610 - t783) * MDP(14) + (pkin(1) * t672 + (t720 * t776 + t658) * qJ(2)) * MDP(7) + (t690 * t776 + t658) * MDP(6) + (t538 * t600 + t578 * t595) * MDP(19) + (t535 * t600 + t574 * t595) * MDP(28) + (-qJD(3) * t595 - qJDD(3) * t600) * MDP(11) + (-g(1) * t558 - g(2) * t560 - t406 * t595 + t489 * t411 + t419 * t600 + t421 * t488 + t453 * t425 + t443 * t785 + (-(-qJD(6) * t431 + t407) * t574 - t430 * t535 - t401 * t600) * t647 + (-(qJD(6) * t430 + t408) * t574 - t431 * t535 - (qJD(6) * t418 + t402) * t600) * t651) * MDP(30) + (t411 * t600 + t425 * t574 + t488 * t535 + t595 * t785) * MDP(26) + (t411 * t488 + t425 * t785) * MDP(24) + t689 * MDP(3) + (t404 * t442 + t429 * t410 + t403 * t441 + t428 * t409 + t444 * t686 + t499 * t667 + (-g(1) * t698 - g(2) * t666) * t654 + (g(1) * t666 - g(2) * t698) * t650) * MDP(23) + (t497 * t600 + t527 * t602 + t557 * t595 + t578 * t665) * MDP(17) + (t497 * t746 + t557 * t665) * MDP(15) + (-(-t552 * t713 + t703) * t578 - t722 * t538 - t663 * t600 - t480 * t595 + t515 * t557 + t551 * t497 + t492 * t746 - g(1) * t579 - g(2) * t581 + t665 * t539) * MDP(21); t786 * MDP(7) - t619 * MDP(13) + t701 * MDP(14) + (t670 - t752) * MDP(20) + (-t578 ^ 2 * t652 - t750 - t754) * MDP(21) + (t448 * t601 + t449 * t674 - t505 * t789 - t680 * t788) * MDP(22) + (-t403 * t674 + t404 * t601 - t428 * t789 - t429 * t788 - t499 * t593 - t779) * MDP(23) + (t687 + t760) * MDP(29) + (-t758 + t775) * MDP(30) + (-t737 - pkin(1) * MDP(7) + (MDP(13) * t649 + MDP(5)) * t642) * qJDD(1) + ((t642 * t717 + t644 * t718 + t593) * MDP(13) + (t591 - t700) * MDP(14)) * qJD(3) + (-MDP(7) * qJ(2) - MDP(6)) * qJD(1) ^ 2 * t720; -t591 ^ 2 * MDP(9) + ((-t591 - t700) * qJD(3) + t701) * MDP(10) + t682 * MDP(11) + qJDD(3) * MDP(12) + (qJD(3) * t549 + t661 + t671) * MDP(13) + (-t591 * t611 - t679 + t774) * MDP(14) + (t557 * t692 + t757) * MDP(15) + ((t497 - t753) * t652 + (-t498 - t751) * t648) * MDP(16) + (t578 * t692 - t750 + t754) * MDP(17) + (t670 + t752) * MDP(18) + (-pkin(3) * t498 - t530 * t578 - t549 * t555 + (t578 * t780 + t664) * t648 + t793 * t652) * MDP(20) + (-pkin(3) * t497 - t549 * t557 + t723 * t578 - t648 * t793 + t664 * t652) * MDP(21) + (-t403 * t601 - t404 * t674 + t428 * t788 - t429 * t789 + t448 * t554 - t449 * t553 + t505 * t725 + t680 * t724 - t774) * MDP(22) + (t404 * t554 + t403 * t553 - t444 * t629 - g(3) * t676 + t778 * t499 + t724 * t429 + t725 * t428 + t689 * (t629 * t631 - t632 * t763)) * MDP(23) + (t411 * t547 + t729 * t785) * MDP(24) + (t411 * t678 - t412 * t547 + t457 * t729 - t728 * t785) * MDP(25) + (-t758 - t775) * MDP(26) + (t687 - t760) * MDP(27) + ((t521 * t651 - t522 * t647) * t535 + t565 * t412 - t421 * t678 + (t647 * t683 - t651 * t684) * t574 - t721 * t457 + t728 * t453 + t661 * t625) * MDP(29) + (-(t521 * t647 + t522 * t651) * t535 + t565 * t411 + t421 * t547 + (t647 * t684 + t651 * t683) * t574 + t721 * t785 + t729 * t453 - t661 * t624) * MDP(30) + (-MDP(13) * t611 - MDP(19) * t578 - MDP(20) * t479 + MDP(21) * t480 - MDP(28) * t574 - MDP(29) * t405 + MDP(30) * t406 - MDP(8) * t591 + MDP(9) * t593) * t593; t557 * t555 * MDP(15) + (-t555 ^ 2 + t557 ^ 2) * MDP(16) + (t497 + t753) * MDP(17) + (-t498 + t751) * MDP(18) + t538 * MDP(19) + (t480 * t578 - t539 * t557 + (t691 + t766) * t648 + t685 + t777) * MDP(20) + (g(1) * t582 - g(2) * t580 + t479 * t578 + t539 * t555 + t652 * t766 - t663) * MDP(21) + ((t448 * t641 - t449 * t643) * pkin(4) + (t428 - t433) * t680 + (-t429 - t432) * t505) * MDP(22) + (-t428 * t432 - t429 * t433 + (t403 * t643 + t404 * t641 - t499 * t557 + t648 * t766 + t777) * pkin(4)) * MDP(23) + (t411 - t759) * MDP(26) + (-t412 + t761) * MDP(27) + (t668 * t535 - (t422 * t651 - t423 * t647) * t574 + t490 * t457 + (-t574 * t669 - t406) * qJD(6) + t791) * MDP(29) + (-t669 * t535 + (t422 * t647 + t423 * t651) * t574 - t490 * t785 + (-t574 * t668 - t732) * qJD(6) + t798) * MDP(30) + t797; (-t505 ^ 2 - t680 ^ 2) * MDP(22) + (-t428 * t505 - t429 * t680 + t444 - t661) * MDP(23) + (t412 + t761) * MDP(29) + (t411 + t759) * MDP(30); (t704 - t759) * MDP(26) + (-t694 + t761) * MDP(27) + (t406 * t574 + t791) * MDP(29) + (t405 * t574 + t798) * MDP(30) + (MDP(26) * t756 - MDP(27) * t785 - MDP(29) * t406 - MDP(30) * t732) * qJD(6) + t797;];
tau  = t1;
