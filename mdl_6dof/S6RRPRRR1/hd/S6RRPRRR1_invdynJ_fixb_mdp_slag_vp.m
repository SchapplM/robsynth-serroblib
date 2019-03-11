% Calculate vector of inverse dynamics joint torques for
% S6RRPRRR1
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
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d4,d5,d6,theta3]';
% MDP [33x1]
%   Minimal dynamic parameter vector (fixed base model)
%   see S6RRPRRR1_convert_par2_MPV_fixb.m
% 
% Output:
% tau [6x1]
%   joint torques of inverse dynamics (contains inertial, gravitational Coriolis and centrifugal forces)

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 13:15
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tau = S6RRPRRR1_invdynJ_fixb_mdp_slag_vp(qJ, qJD, qJDD, g, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(3,1),zeros(11,1),zeros(33,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPRRR1_invdynJ_fixb_mdp_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRPRRR1_invdynJ_fixb_mdp_slag_vp: qJD has to be [6x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [6 1]), ...
  'S6RRPRRR1_invdynJ_fixb_mdp_slag_vp: qJDD has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RRPRRR1_invdynJ_fixb_mdp_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6RRPRRR1_invdynJ_fixb_mdp_slag_vp: pkin has to be [11x1] (double)');
assert(isreal(MDP) && all(size(MDP) == [33 1]), ...
  'S6RRPRRR1_invdynJ_fixb_mdp_slag_vp: MDP has to be [33x1] (double)'); 

%% Symbolic Calculation
% From invdyn_joint_fixb_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 13:14:46
% EndTime: 2019-03-09 13:14:59
% DurationCPUTime: 9.70s
% Computational Cost: add. (12810->510), mult. (31585->651), div. (0->0), fcn. (25747->16), ass. (0->241)
t629 = sin(qJ(5));
t634 = cos(qJ(5));
t625 = sin(pkin(11));
t626 = cos(pkin(11));
t631 = sin(qJ(2));
t636 = cos(qJ(2));
t577 = -t625 * t631 + t626 * t636;
t565 = t577 * qJD(1);
t578 = t625 * t636 + t626 * t631;
t567 = t578 * qJD(1);
t630 = sin(qJ(4));
t635 = cos(qJ(4));
t661 = t565 * t630 + t635 * t567;
t662 = t565 * t635 - t630 * t567;
t472 = t629 * t661 - t634 * t662;
t633 = cos(qJ(6));
t706 = qJD(6) * t633;
t796 = t472 * t633 + t706;
t566 = t578 * qJD(2);
t530 = -qJD(1) * t566 + qJDD(1) * t577;
t702 = qJD(1) * qJD(2);
t694 = t636 * t702;
t695 = t631 * t702;
t531 = qJDD(1) * t578 - t625 * t695 + t626 * t694;
t710 = qJD(4) * t635;
t711 = qJD(4) * t630;
t456 = t630 * t530 + t635 * t531 + t565 * t710 - t567 * t711;
t457 = qJD(4) * t661 - t635 * t530 + t531 * t630;
t754 = t629 * t662 + t634 * t661;
t424 = qJD(5) * t754 + t456 * t629 + t634 * t457;
t422 = qJDD(6) + t424;
t418 = t633 * t422;
t622 = qJD(2) + qJD(4);
t619 = qJD(5) + t622;
t628 = sin(qJ(6));
t465 = t619 * t628 + t633 * t754;
t463 = -t633 * t619 + t628 * t754;
t737 = t463 * t754;
t774 = qJD(6) + t472;
t708 = qJD(5) * t634;
t709 = qJD(5) * t629;
t423 = t634 * t456 - t629 * t457 - t661 * t709 + t662 * t708;
t621 = qJDD(2) + qJDD(4);
t618 = qJDD(5) + t621;
t696 = t633 * t423 + t628 * t618 + t619 * t706;
t707 = qJD(6) * t628;
t413 = -t707 * t754 + t696;
t411 = t413 * t628;
t417 = t628 * t422;
t731 = t472 * t619;
t734 = t754 * t619;
t736 = t465 * t754;
t783 = (t423 + t731) * MDP(22) + (-t424 + t734) * MDP(23) - t472 ^ 2 * MDP(21) + (MDP(20) * t472 + MDP(21) * t754 - MDP(31) * t774) * t754 + (t796 * t465 + t411) * MDP(27) + (t796 * t774 + t417 - t736) * MDP(29);
t412 = t413 * t633;
t601 = t633 * t618;
t414 = t465 * qJD(6) + t423 * t628 - t601;
t784 = -t628 * t414 - t796 * t463 + t412;
t793 = t774 * t628;
t795 = t783 + (-t465 * t793 + t784) * MDP(28) + (-t774 * t793 + t418 + t737) * MDP(30);
t620 = qJ(2) + pkin(11) + qJ(4);
t611 = qJ(5) + t620;
t605 = sin(t611);
t632 = sin(qJ(1));
t637 = cos(qJ(1));
t670 = g(1) * t637 + g(2) * t632;
t794 = t670 * t605;
t743 = qJ(3) + pkin(7);
t596 = t743 * t636;
t583 = qJD(1) * t596;
t570 = t625 * t583;
t595 = t743 * t631;
t582 = qJD(1) * t595;
t742 = qJD(2) * pkin(2);
t574 = -t582 + t742;
t528 = t626 * t574 - t570;
t748 = pkin(8) * t567;
t493 = qJD(2) * pkin(3) + t528 - t748;
t727 = t626 * t583;
t529 = t625 * t574 + t727;
t749 = pkin(8) * t565;
t499 = t529 + t749;
t686 = t635 * t493 - t499 * t630;
t767 = pkin(9) * t661;
t446 = t686 - t767;
t444 = pkin(4) * t622 + t446;
t665 = -t493 * t630 - t499 * t635;
t766 = pkin(9) * t662;
t447 = -t665 + t766;
t738 = t447 * t634;
t426 = t444 * t629 + t738;
t421 = pkin(10) * t619 + t426;
t615 = pkin(2) * t636 + pkin(1);
t589 = -qJD(1) * t615 + qJD(3);
t537 = -pkin(3) * t565 + t589;
t484 = -pkin(4) * t662 + t537;
t435 = pkin(5) * t472 - pkin(10) * t754 + t484;
t404 = -t421 * t628 + t435 * t633;
t739 = t447 * t629;
t425 = t444 * t634 - t739;
t420 = -pkin(5) * t619 - t425;
t792 = -t404 * t754 + t420 * t707 + t633 * t794;
t405 = t421 * t633 + t435 * t628;
t692 = qJD(2) * t743;
t563 = -qJD(3) * t631 - t636 * t692;
t527 = qJDD(2) * pkin(2) + qJD(1) * t563 - qJDD(1) * t595;
t562 = qJD(3) * t636 - t631 * t692;
t534 = qJD(1) * t562 + qJDD(1) * t596;
t478 = t626 * t527 - t534 * t625;
t458 = qJDD(2) * pkin(3) - pkin(8) * t531 + t478;
t479 = t625 * t527 + t626 * t534;
t459 = pkin(8) * t530 + t479;
t647 = qJD(4) * t665 + t635 * t458 - t630 * t459;
t408 = pkin(4) * t621 - pkin(9) * t456 + t647;
t756 = t635 * (qJD(4) * t493 + t459) + t630 * t458 - t499 * t711;
t409 = -pkin(9) * t457 + t756;
t755 = qJD(5) * t426 - t634 * t408 + t629 * t409;
t399 = -pkin(5) * t618 + t755;
t606 = cos(t611);
t745 = g(3) * t606;
t778 = t399 + t745;
t791 = t405 * t754 + t420 * t706 + t778 * t628;
t781 = t420 * t472;
t535 = t582 * t625 - t727;
t500 = t535 - t749;
t536 = -t626 * t582 - t570;
t501 = t536 - t748;
t610 = pkin(2) * t626 + pkin(3);
t750 = pkin(2) * t625;
t673 = t635 * t610 - t630 * t750;
t789 = -t673 * qJD(4) + t630 * t500 + t635 * t501;
t561 = t610 * t630 + t635 * t750;
t788 = -t561 * qJD(4) - t635 * t500 + t501 * t630;
t787 = -t628 * t794 + t791;
t786 = -t633 * t778 + t792;
t768 = pkin(5) * t754;
t785 = pkin(10) * t472 + t768;
t599 = g(3) * t605;
t757 = (qJD(5) * t444 + t409) * t634 + t629 * t408 - t447 * t709;
t773 = t472 * t484 + t670 * t606 + t599 - t757;
t777 = -t766 - t788;
t776 = -t767 + t789;
t760 = -t484 * t754 - t745 - t755 + t794;
t604 = t618 * MDP(24);
t728 = t662 * t622;
t729 = t661 * t622;
t772 = t621 * MDP(17) + t604 + (t661 ^ 2 - t662 ^ 2) * MDP(14) - t662 * MDP(13) * t661 + (t456 - t728) * MDP(15) + (-t457 + t729) * MDP(16);
t769 = pkin(4) * t661;
t763 = t774 * (t774 * pkin(10) + t768);
t538 = -t626 * t595 - t596 * t625;
t514 = -pkin(8) * t578 + t538;
t539 = -t625 * t595 + t626 * t596;
t515 = pkin(8) * t577 + t539;
t717 = t630 * t514 + t635 * t515;
t560 = pkin(4) + t673;
t716 = t629 * t560 + t634 * t561;
t608 = sin(t620);
t609 = cos(t620);
t759 = -g(3) * t609 - t537 * t661 + t670 * t608 + t647;
t758 = g(3) * t608 - t537 * t662 + t670 * t609 - t756;
t533 = t577 * t630 + t578 * t635;
t569 = t577 * qJD(2);
t481 = qJD(4) * t533 + t635 * t566 + t569 * t630;
t512 = -t562 * t625 + t626 * t563;
t489 = -pkin(8) * t569 + t512;
t513 = t626 * t562 + t625 * t563;
t490 = -pkin(8) * t566 + t513;
t655 = t630 * t489 + t635 * t490 + t514 * t710 - t515 * t711;
t431 = -pkin(9) * t481 + t655;
t660 = t635 * t577 - t578 * t630;
t480 = qJD(4) * t660 - t566 * t630 + t569 * t635;
t645 = -qJD(4) * t717 + t635 * t489 - t490 * t630;
t432 = -pkin(9) * t480 + t645;
t684 = t635 * t514 - t515 * t630;
t450 = -pkin(9) * t533 + t684;
t451 = pkin(9) * t660 + t717;
t666 = t450 * t634 - t451 * t629;
t400 = qJD(5) * t666 + t431 * t634 + t432 * t629;
t434 = t450 * t629 + t451 * t634;
t482 = t533 * t629 - t634 * t660;
t437 = -qJD(5) * t482 + t480 * t634 - t481 * t629;
t483 = t533 * t634 + t629 * t660;
t544 = -pkin(3) * t577 - t615;
t494 = -pkin(4) * t660 + t544;
t440 = pkin(5) * t482 - pkin(10) * t483 + t494;
t398 = pkin(10) * t618 + t757;
t678 = qJD(6) * t435 + t398;
t753 = t399 * t483 + t420 * t437 - t434 * t422 - (qJD(6) * t440 + t400) * t774 - t482 * t678;
t744 = g(3) * t636;
t741 = t420 * t483;
t740 = t440 * t422;
t735 = t465 * t628;
t726 = t628 * t632;
t725 = t628 * t637;
t724 = t632 * t633;
t723 = t633 * t637;
t663 = t560 * t634 - t561 * t629;
t720 = -qJD(5) * t663 + t777 * t629 + t776 * t634;
t719 = t716 * qJD(5) - t776 * t629 + t777 * t634;
t623 = t631 ^ 2;
t714 = -t636 ^ 2 + t623;
t701 = qJDD(1) * t636;
t617 = t631 * t742;
t541 = pkin(3) * t566 + t617;
t540 = t631 * qJD(1) * pkin(2) + pkin(3) * t567;
t653 = pkin(2) * t695 - qJDD(1) * t615 + qJDD(3);
t495 = -pkin(3) * t530 + t653;
t443 = pkin(4) * t457 + t495;
t403 = pkin(5) * t424 - pkin(10) * t423 + t443;
t676 = qJD(6) * t421 - t403;
t485 = t540 + t769;
t508 = pkin(10) + t716;
t675 = qJD(6) * t508 + t485 + t785;
t613 = pkin(4) * t629 + pkin(10);
t674 = qJD(6) * t613 + t769 + t785;
t427 = t446 * t629 + t738;
t672 = pkin(4) * t709 - t427;
t428 = t446 * t634 - t739;
t671 = -pkin(4) * t708 + t428;
t669 = g(1) * t632 - g(2) * t637;
t668 = -t422 * t613 + t781;
t667 = -t422 * t508 + t781;
t462 = pkin(4) * t481 + t541;
t659 = t418 - (t472 * t628 + t707) * t774;
t657 = -0.2e1 * pkin(1) * t702 - pkin(7) * qJDD(2);
t656 = t437 * t633 - t483 * t707;
t654 = -pkin(10) * t422 + t425 * t774 + t781;
t638 = qJD(2) ^ 2;
t650 = 0.2e1 * qJDD(1) * pkin(1) - pkin(7) * t638 + t669;
t639 = qJD(1) ^ 2;
t649 = pkin(1) * t639 - pkin(7) * qJDD(1) + t670;
t614 = -pkin(4) * t634 - pkin(5);
t554 = t606 * t723 + t726;
t553 = -t606 * t725 + t724;
t552 = -t606 * t724 + t725;
t551 = t606 * t726 + t723;
t507 = -pkin(5) - t663;
t438 = qJD(5) * t483 + t480 * t629 + t634 * t481;
t406 = pkin(5) * t438 - pkin(10) * t437 + t462;
t402 = t633 * t403;
t401 = qJD(5) * t434 + t431 * t629 - t432 * t634;
t1 = [(-t481 * t622 + t621 * t660) * MDP(16) + (t544 * t457 + t537 * t481 - t495 * t660 - t541 * t662 + t609 * t669 + t621 * t684 + t622 * t645) * MDP(18) + (t456 * t660 - t457 * t533 + t480 * t662 - t481 * t661) * MDP(14) + (t479 * t539 + t529 * t513 + t478 * t538 + t528 * t512 - t653 * t615 + t589 * t617 - g(1) * (-t615 * t632 + t637 * t743) - g(2) * (t615 * t637 + t632 * t743)) * MDP(12) + (-t400 * t619 + t423 * t494 - t434 * t618 + t437 * t484 + t443 * t483 + t462 * t754 - t605 * t669) * MDP(26) + (t423 * t483 + t437 * t754) * MDP(20) + (-t423 * t482 - t424 * t483 - t437 * t472 - t438 * t754) * MDP(21) + (-t401 * t619 + t424 * t494 + t438 * t484 + t443 * t482 + t462 * t472 + t606 * t669 + t618 * t666) * MDP(25) + (-g(1) * t551 - g(2) * t553 + t401 * t465 - t405 * t438 - t666 * t413 + (-(-qJD(6) * t434 + t406) * t774 - t740 + t676 * t482 - qJD(6) * t741) * t628 + t753 * t633) * MDP(33) + (-g(1) * t552 - g(2) * t554 + t401 * t463 + t402 * t482 + t404 * t438 - t666 * t414 + (t406 * t774 + t740 + (-t421 * t482 - t434 * t774 + t741) * qJD(6)) * t633 + t753 * t628) * MDP(32) + (t422 * t482 + t438 * t774) * MDP(31) + (t413 * t482 + t418 * t483 + t438 * t465 + t656 * t774) * MDP(29) + (-t483 * t417 - t414 * t482 - t438 * t463 + (-t437 * t628 - t483 * t706) * t774) * MDP(30) + t669 * MDP(2) + t670 * MDP(3) + (t631 * t657 + t636 * t650) * MDP(9) + (-t631 * t650 + t636 * t657) * MDP(10) + (-t478 * t578 + t479 * t577 - t512 * t567 + t513 * t565 - t528 * t569 - t529 * t566 + t530 * t539 - t531 * t538 - t670) * MDP(11) + (qJDD(2) * t636 - t631 * t638) * MDP(7) + (qJDD(2) * t631 + t636 * t638) * MDP(6) + (t480 * t622 + t533 * t621) * MDP(15) + (t437 * t619 + t483 * t618) * MDP(22) + (-t438 * t619 - t482 * t618) * MDP(23) + (t456 * t533 + t480 * t661) * MDP(13) + (t544 * t456 + t537 * t480 + t495 * t533 + t541 * t661 - t608 * t669 - t621 * t717 - t622 * t655) * MDP(19) + (qJDD(1) * t623 + 0.2e1 * t631 * t694) * MDP(4) + qJDD(1) * MDP(1) + 0.2e1 * (t631 * t701 - t702 * t714) * MDP(5) + (t412 * t483 + t465 * t656) * MDP(27) + ((-t463 * t633 - t735) * t437 + (-t411 - t414 * t633 + (t463 * t628 - t465 * t633) * qJD(6)) * t483) * MDP(28); t772 + (-t485 * t472 + t618 * t663 - t619 * t719 + t760) * MDP(25) + t783 + (t507 * t413 + t667 * t633 + t719 * t465 + (t628 * t675 + t633 * t720) * t774 + t787) * MDP(33) + (t507 * t414 + t667 * t628 + t719 * t463 + (t628 * t720 - t633 * t675) * t774 + t786) * MDP(32) + (-t735 * t774 + t784) * MDP(28) + (-MDP(4) * t631 * t636 + MDP(5) * t714) * t639 + (g(3) * t631 + t636 * t649) * MDP(10) + ((t529 + t535) * t567 + (t528 - t536) * t565 + (t530 * t625 - t531 * t626) * pkin(2)) * MDP(11) + MDP(7) * t701 + (-t485 * t754 - t618 * t716 + t619 * t720 + t773) * MDP(26) + qJDD(2) * MDP(8) + (t540 * t662 + t673 * t621 + t622 * t788 + t759) * MDP(18) + (-t540 * t661 - t561 * t621 + t622 * t789 + t758) * MDP(19) + (t659 + t737) * MDP(30) + (-t528 * t535 - t529 * t536 + (-t744 + t478 * t626 + t479 * t625 + (-qJD(1) * t589 + t670) * t631) * pkin(2)) * MDP(12) + (t631 * t649 - t744) * MDP(9) + t631 * qJDD(1) * MDP(6); (-t565 ^ 2 - t567 ^ 2) * MDP(11) + (t528 * t567 - t529 * t565 + t653 - t669) * MDP(12) + (t457 + t729) * MDP(18) + (t456 + t728) * MDP(19) + (t424 + t734) * MDP(25) + (t423 - t731) * MDP(26) + (t659 - t737) * MDP(32) + (-t633 * t774 ^ 2 - t417 - t736) * MDP(33); (-t622 * t665 + t759) * MDP(18) + (t622 * t686 + t758) * MDP(19) + (t427 * t619 + (-t472 * t661 + t618 * t634 - t619 * t709) * pkin(4) + t760) * MDP(25) + (t428 * t619 + (-t618 * t629 - t619 * t708 - t661 * t754) * pkin(4) + t773) * MDP(26) + (t614 * t414 + t668 * t628 + t672 * t463 + (t628 * t671 - t633 * t674) * t774 + t786) * MDP(32) + (t614 * t413 + t668 * t633 + t672 * t465 + (t628 * t674 + t633 * t671) * t774 + t787) * MDP(33) + t772 + t795; t604 + (t426 * t619 + t760) * MDP(25) + (t425 * t619 + t773) * MDP(26) + (-pkin(5) * t414 - t426 * t463 + t654 * t628 + (-t778 - t763) * t633 + t792) * MDP(32) + (-pkin(5) * t413 - t426 * t465 + t654 * t633 + (-t794 + t763) * t628 + t791) * MDP(33) + t795; t465 * t463 * MDP(27) + (-t463 ^ 2 + t465 ^ 2) * MDP(28) + (t463 * t774 + t696) * MDP(29) + (t465 * t774 + t601) * MDP(30) + t422 * MDP(31) + (-g(1) * t553 + g(2) * t551 + t405 * t774 - t420 * t465 + t402) * MDP(32) + (g(1) * t554 - g(2) * t552 + t404 * t774 + t420 * t463) * MDP(33) + ((-t398 + t599) * MDP(33) + (-MDP(30) * t754 - MDP(32) * t421 - MDP(33) * t435) * qJD(6)) * t633 + (-qJD(6) * t754 * MDP(29) + (-qJD(6) * t619 - t423) * MDP(30) + (-t678 + t599) * MDP(32) + t676 * MDP(33)) * t628;];
tau  = t1;
