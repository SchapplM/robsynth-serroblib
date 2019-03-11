% Calculate vector of inverse dynamics joint torques for
% S6RRRRPR3
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
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d3,d4,d6]';
% MDP [35x1]
%   Minimal dynamic parameter vector (fixed base model)
%   see S6RRRRPR3_convert_par2_MPV_fixb.m
% 
% Output:
% tau [6x1]
%   joint torques of inverse dynamics (contains inertial, gravitational Coriolis and centrifugal forces)

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 22:05
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tau = S6RRRRPR3_invdynJ_fixb_mdp_slag_vp(qJ, qJD, qJDD, g, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(3,1),zeros(10,1),zeros(35,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRRPR3_invdynJ_fixb_mdp_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRRRPR3_invdynJ_fixb_mdp_slag_vp: qJD has to be [6x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [6 1]), ...
  'S6RRRRPR3_invdynJ_fixb_mdp_slag_vp: qJDD has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RRRRPR3_invdynJ_fixb_mdp_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RRRRPR3_invdynJ_fixb_mdp_slag_vp: pkin has to be [10x1] (double)');
assert(isreal(MDP) && all(size(MDP) == [35 1]), ...
  'S6RRRRPR3_invdynJ_fixb_mdp_slag_vp: MDP has to be [35x1] (double)'); 

%% Symbolic Calculation
% From invdyn_joint_fixb_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 22:04:16
% EndTime: 2019-03-09 22:04:30
% DurationCPUTime: 10.21s
% Computational Cost: add. (9663->595), mult. (22321->723), div. (0->0), fcn. (16729->14), ass. (0->265)
t628 = sin(qJ(2));
t631 = cos(qJ(2));
t785 = sin(qJ(3));
t720 = qJD(1) * t785;
t787 = cos(qJ(3));
t721 = qJD(1) * t787;
t552 = -t628 * t720 + t631 * t721;
t553 = -t628 * t721 - t631 * t720;
t627 = sin(qJ(4));
t562 = t628 * t787 + t631 * t785;
t621 = qJD(2) + qJD(3);
t804 = t621 * t562;
t643 = t804 * qJD(1);
t675 = t628 * t785 - t631 * t787;
t638 = -t675 * qJDD(1) - t643;
t644 = t621 * t675;
t639 = -t644 * qJD(1) + t562 * qJDD(1);
t786 = cos(qJ(4));
t718 = qJD(4) * t786;
t738 = qJD(4) * t627;
t679 = t552 * t718 + t553 * t738 + t627 * t638 + t786 * t639;
t451 = -qJDD(6) - t679;
t630 = cos(qJ(6));
t447 = t630 * t451;
t683 = -t627 * t552 + t553 * t786;
t801 = qJD(6) - t683;
t626 = sin(qJ(6));
t811 = t801 * t626;
t815 = -t801 * t811 - t447;
t771 = t451 * t626;
t810 = t801 * t630;
t814 = -t801 * t810 + t771;
t614 = qJD(4) + t621;
t519 = t786 * t552 + t553 * t627;
t620 = qJDD(2) + qJDD(3);
t613 = qJDD(4) + t620;
t636 = qJD(4) * t683 - t627 * t639 + t786 * t638;
t734 = qJD(6) * t630;
t725 = -t519 * t734 + t630 * t613 - t626 * t636;
t735 = qJD(6) * t626;
t431 = -t614 * t735 + t725;
t430 = t431 * t630;
t449 = t630 * t636;
t505 = -t519 * t626 + t614 * t630;
t432 = t505 * qJD(6) + t613 * t626 + t449;
t664 = -t519 * t614 + t679;
t503 = t519 * t630 + t614 * t626;
t812 = t503 * t801;
t813 = t613 * MDP(22) + (-t614 * t683 + t636) * MDP(21) + t664 * MDP(20) + (-t505 * t811 + t430) * MDP(29) + (-t505 * t810 + (-t431 + t812) * t626 - t630 * t432) * MDP(30);
t548 = t553 * pkin(9);
t788 = pkin(7) + pkin(8);
t581 = t788 * t631;
t569 = qJD(1) * t581;
t554 = t785 * t569;
t580 = t788 * t628;
t567 = qJD(1) * t580;
t774 = qJD(2) * pkin(2);
t560 = -t567 + t774;
t694 = t787 * t560 - t554;
t501 = t548 + t694;
t493 = t621 * pkin(3) + t501;
t558 = t787 * t569;
t676 = -t560 * t785 - t558;
t776 = t552 * pkin(9);
t502 = -t676 + t776;
t497 = t786 * t502;
t468 = t627 * t493 + t497;
t463 = -qJ(5) * t614 - t468;
t806 = t519 * pkin(5);
t438 = -t463 + t806;
t712 = t801 * t438;
t773 = qJ(5) * t519;
t485 = -pkin(4) * t683 - t773;
t790 = t683 ^ 2;
t778 = t683 * pkin(5);
t619 = t631 * pkin(2);
t775 = pkin(1) + t619;
t805 = t438 * t683;
t695 = t567 * t785 - t558;
t507 = t695 - t776;
t748 = -t787 * t567 - t554;
t508 = t548 + t748;
t609 = pkin(2) * t787 + pkin(3);
t717 = t785 * qJD(3);
t719 = qJD(3) * t787;
t784 = pkin(2) * t627;
t803 = t627 * t507 + (qJD(4) * t785 + t717) * t784 - t609 * t718 + (-pkin(2) * t719 + t508) * t786;
t625 = qJ(2) + qJ(3);
t617 = cos(t625);
t618 = qJ(4) + t625;
t606 = sin(t618);
t607 = cos(t618);
t743 = t607 * pkin(4) + t606 * qJ(5);
t724 = pkin(3) * t617 + t743;
t496 = t627 * t502;
t467 = -t786 * t493 + t496;
t732 = qJD(5) + t467;
t579 = t775 * qJD(1);
t533 = -pkin(3) * t552 - t579;
t657 = qJ(5) * t683 + t533;
t471 = -pkin(4) * t519 + t657;
t597 = g(3) * t607;
t729 = qJD(1) * qJD(2);
t715 = t631 * t729;
t728 = qJDD(1) * t628;
t531 = qJDD(2) * pkin(2) + t788 * (-t715 - t728);
t716 = t628 * t729;
t727 = qJDD(1) * t631;
t532 = t788 * (-t716 + t727);
t651 = qJD(3) * t676 + t787 * t531 - t532 * t785;
t441 = t620 * pkin(3) - pkin(9) * t639 + t651;
t653 = t531 * t785 + t532 * t787 + t560 * t719 - t569 * t717;
t450 = pkin(9) * t638 + t653;
t706 = -t786 * t441 + t627 * t450 + t493 * t738 + t502 * t718;
t632 = cos(qJ(1));
t764 = t606 * t632;
t629 = sin(qJ(1));
t765 = t606 * t629;
t678 = -g(1) * t764 - g(2) * t765 + t597 + t706;
t800 = -t471 * t683 + qJDD(5) + t678;
t660 = t533 * t683 - t678;
t704 = t785 * t786;
t749 = t507 * t786 - t627 * t508 + t609 * t738 + (qJD(4) * t704 + (t627 * t787 + t704) * qJD(3)) * pkin(2);
t798 = t749 * t614;
t469 = t627 * t501 + t497;
t701 = pkin(3) * t738 - t469;
t750 = -qJD(5) + t803;
t470 = t501 * t786 - t496;
t745 = -pkin(3) * t718 - qJD(5) + t470;
t595 = t613 * qJ(5);
t797 = -t614 * qJD(5) - t595;
t536 = pkin(3) * t675 - t775;
t796 = MDP(18) * t683 - MDP(33) * t801;
t616 = sin(t625);
t700 = -pkin(3) * t616 - pkin(4) * t606;
t733 = t732 - t778;
t608 = -pkin(3) * t786 - pkin(4);
t601 = -pkin(10) + t608;
t795 = t801 * (-t806 + t701) - t601 * t451;
t696 = t609 * t786 - t785 * t784;
t546 = -pkin(4) - t696;
t539 = -pkin(10) + t546;
t794 = -t801 * (t806 - t749) - t539 * t451;
t722 = qJD(2) * t788;
t568 = t628 * t722;
t570 = t631 * t722;
t672 = -t787 * t568 - t785 * t570 - t580 * t719 - t581 * t717;
t480 = -pkin(9) * t804 + t672;
t697 = t568 * t785 - t787 * t570;
t481 = pkin(9) * t644 + t580 * t717 - t581 * t719 + t697;
t693 = -t787 * t580 - t581 * t785;
t512 = -t562 * pkin(9) + t693;
t747 = -t785 * t580 + t787 * t581;
t513 = -pkin(9) * t675 + t747;
t427 = -t786 * t480 - t627 * t481 - t512 * t718 + t513 * t738;
t484 = t627 * t512 + t513 * t786;
t698 = g(1) * t629 - g(2) * t632;
t793 = -t427 * t614 + t484 * t613 + t606 * t698;
t428 = qJD(4) * t484 + t627 * t480 - t481 * t786;
t483 = -t512 * t786 + t627 * t513;
t792 = t428 * t614 + t483 * t613 - t607 * t698;
t789 = pkin(4) + pkin(10);
t783 = pkin(3) * t553;
t780 = pkin(4) * t613;
t596 = g(3) * t606;
t772 = qJDD(1) * pkin(1);
t658 = t786 * t675;
t529 = t562 * t627 + t658;
t669 = t627 * t675;
t530 = t562 * t786 - t669;
t652 = -t530 * qJ(5) + t536;
t461 = t529 * t789 + t652;
t770 = t461 * t451;
t769 = t468 * t614;
t768 = t529 * t626;
t763 = t607 * t629;
t762 = t607 * t632;
t761 = t616 * t629;
t760 = t616 * t632;
t759 = t617 * t629;
t758 = t617 * t632;
t757 = t626 * t629;
t756 = t626 * t632;
t755 = t629 * t630;
t754 = t630 * t632;
t705 = -t627 * t441 - t786 * t450 - t493 * t718 + t502 * t738;
t415 = t705 + t797;
t413 = pkin(5) * t636 - t415;
t753 = t413 * t626 + t438 * t734;
t751 = -t778 - t750;
t744 = -t778 - t745;
t623 = t628 ^ 2;
t742 = -t631 ^ 2 + t623;
t739 = qJD(1) * t628;
t456 = -t519 * t789 + t657;
t737 = qJD(6) * t456;
t736 = qJD(6) * t614;
t612 = t628 * t774;
t714 = -pkin(2) * t628 + t700;
t411 = t413 * t630;
t436 = -t614 * t789 + t733;
t425 = t436 * t626 + t456 * t630;
t713 = t425 * t519 + t411;
t598 = pkin(2) * t716;
t417 = -pkin(2) * t727 - pkin(3) * t638 - pkin(4) * t636 - qJ(5) * t679 + qJD(5) * t683 + t598 - t772;
t414 = -pkin(10) * t636 + t417;
t707 = qJD(6) * t436 + t414;
t699 = g(1) * t632 + g(2) * t629;
t692 = -t737 + t597;
t689 = (t468 + t806) * t801 - t789 * t451;
t462 = -pkin(4) * t614 + t732;
t688 = -t462 * t519 + t463 * t683;
t687 = qJDD(5) + t706;
t424 = t436 * t630 - t456 * t626;
t686 = -t424 * t519 - t630 * t805 + t753;
t685 = t775 + t724;
t684 = -0.2e1 * pkin(1) * t729 - pkin(7) * qJDD(2);
t476 = -qJD(4) * t669 + t562 * t718 - t627 * t644 + t786 * t804;
t682 = t476 * t626 + t529 * t734;
t677 = g(1) * t762 + g(2) * t763 + t596 + t705;
t477 = -t783 + t485;
t464 = t530 * pkin(5) + t483;
t673 = t413 * t529 + t438 * t476 + t464 * t451;
t671 = pkin(2) * t704 + t627 * t609;
t611 = pkin(2) * t739;
t474 = t477 + t611;
t668 = -t607 * t699 - t596;
t667 = -t677 - t797;
t634 = qJD(2) ^ 2;
t663 = -pkin(7) * t634 + t698 + 0.2e1 * t772;
t635 = qJD(1) ^ 2;
t662 = pkin(1) * t635 - pkin(7) * qJDD(1) + t699;
t661 = -t533 * t519 + t677;
t654 = t471 * t519 + t667;
t514 = t683 * pkin(10);
t650 = (-qJD(6) * t539 + t474 - t514) * t801 + t668;
t649 = (-qJD(6) * t601 + t477 - t514) * t801 + t668;
t648 = (t789 * t801 - t773) * t801 + t668;
t647 = g(1) * t758 + g(2) * t759 + g(3) * t616 + t579 * t552 - t653;
t641 = g(1) * t760 + g(2) * t761 - g(3) * t617 - t579 * t553 + t651;
t522 = pkin(3) * t804 + t612;
t640 = t553 * t552 * MDP(11) + (-t505 * t519 + t815) * MDP(31) + (t503 * t519 + t814) * MDP(32) + (-t519 ^ 2 + t790) * MDP(19) + (-t552 * t621 + t639) * MDP(13) + (-t553 * t621 + t638) * MDP(14) + (-t552 ^ 2 + t553 ^ 2) * MDP(12) + t620 * MDP(15) + t796 * t519 + t813;
t475 = qJD(4) * t658 + t562 * t738 + t627 * t804 + t644 * t786;
t426 = t476 * pkin(4) + t475 * qJ(5) - t530 * qJD(5) + t522;
t622 = -pkin(9) - t788;
t602 = pkin(3) * t627 + qJ(5);
t574 = qJ(5) * t762;
t573 = qJ(5) * t763;
t549 = -qJDD(1) * t775 + t598;
t545 = qJ(5) + t671;
t543 = -t606 * t757 + t754;
t542 = t606 * t755 + t756;
t541 = t606 * t756 + t755;
t540 = t606 * t754 - t757;
t534 = t611 - t783;
t488 = pkin(3) * t643 + qJDD(1) * t536 + t598;
t482 = t529 * pkin(4) + t652;
t465 = -t529 * pkin(5) + t484;
t420 = -t475 * pkin(5) + t428;
t419 = -pkin(5) * t476 - t427;
t418 = t476 * pkin(10) + t426;
t416 = t687 - t780;
t412 = pkin(5) * t679 - t613 * t789 + t687;
t409 = t630 * t412;
t1 = [(g(1) * t542 - g(2) * t540 + t419 * t505 + t425 * t475 + t465 * t431 + (-(qJD(6) * t464 + t418) * t801 + t770 - t707 * t530 + t438 * qJD(6) * t529) * t630 + (-(-qJD(6) * t461 + t420) * t801 - (t412 - t737) * t530 + t673) * t626) * MDP(35) + (-g(1) * t543 - g(2) * t541 + t409 * t530 + t419 * t503 - t424 * t475 + t465 * t432 + (-t414 * t530 - t418 * t801 + t770) * t626 + (t420 * t801 - t673) * t630 + ((-t461 * t630 - t464 * t626) * t801 - t425 * t530 + t438 * t768) * qJD(6)) * MDP(34) + (-t620 * t675 - t621 * t804) * MDP(14) + (-t552 * t644 + t553 * t804 + t562 * t638 - t639 * t675) * MDP(12) + (-g(1) * t761 + g(2) * t760 + t549 * t562 - t553 * t612 + t579 * t644 - t620 * t747 - t621 * t672 - t639 * t775) * MDP(17) + (t476 * t533 + t488 * t529 - t519 * t522 - t536 * t636 - t792) * MDP(23) + (-t417 * t529 + t426 * t519 - t471 * t476 + t482 * t636 + t792) * MDP(26) + (-t475 * t519 + t476 * t683 - t529 * t679 + t530 * t636) * MDP(19) + (t415 * t529 + t416 * t530 - t427 * t519 - t428 * t683 - t462 * t475 + t463 * t476 + t483 * t679 + t484 * t636 - t699) * MDP(25) + (-t451 * t530 - t475 * t801) * MDP(33) + (-t529 * t447 - t432 * t530 + t475 * t503 + (t476 * t630 - t529 * t735) * t801) * MDP(32) + (t431 * t530 - t451 * t768 - t475 * t505 + t682 * t801) * MDP(31) + (-t552 * t612 + t775 * t638 + t549 * t675 - t579 * t804 + (-qJD(3) * t747 + t697) * t621 + t693 * t620 + g(1) * t759 - g(2) * t758) * MDP(16) + (-t628 * t663 + t684 * t631) * MDP(10) + (t628 * t684 + t631 * t663) * MDP(9) + (-t475 * t533 + t488 * t530 - t522 * t683 + t536 * t679 - t793) * MDP(24) + (-t417 * t530 + t426 * t683 + t471 * t475 - t482 * t679 + t793) * MDP(27) + (t475 * t683 + t530 * t679) * MDP(18) + (t562 * t620 - t621 * t644) * MDP(13) + (t553 * t644 + t562 * t639) * MDP(11) + qJDD(1) * MDP(1) + (qJDD(2) * t628 + t631 * t634) * MDP(6) + (qJDD(2) * t631 - t628 * t634) * MDP(7) + (-t475 * t614 + t530 * t613) * MDP(20) + (-t476 * t614 - t529 * t613) * MDP(21) + (-t415 * t484 + t416 * t483 + t417 * t482 + t471 * t426 + t463 * t427 + t462 * t428 + (g(1) * t622 - g(2) * t685) * t632 + (g(1) * t685 + g(2) * t622) * t629) * MDP(28) + t698 * MDP(2) + t699 * MDP(3) + (qJDD(1) * t623 + 0.2e1 * t628 * t715) * MDP(4) + 0.2e1 * (t628 * t727 - t729 * t742) * MDP(5) + (t431 * t768 + t505 * t682) * MDP(29) + ((-t503 * t626 + t505 * t630) * t476 + (t430 - t432 * t626 + (-t503 * t630 - t505 * t626) * qJD(6)) * t529) * MDP(30); (t519 * t534 + t613 * t696 + t660 - t798) * MDP(23) + (-t415 * t545 + t416 * t546 - t471 * t474 - g(1) * (t632 * t714 + t574) - g(2) * (t629 * t714 + t573) - g(3) * (t619 + t724) + t750 * t463 + t749 * t462) * MDP(28) + (t748 * t621 + (t553 * t739 - t620 * t785 - t621 * t719) * pkin(2) + t647) * MDP(17) + (-t695 * t621 + (t552 * t739 + t620 * t787 - t621 * t717) * pkin(2) + t641) * MDP(16) + (-t474 * t683 + t545 * t613 - t614 * t750 + t654) * MDP(27) + (-t474 * t519 + t798 + (-pkin(4) + t546) * t613 + t800) * MDP(26) + (-g(3) * t631 + t628 * t662) * MDP(9) + (t545 * t431 + t751 * t505 + t650 * t630 + (-t712 - t794) * t626 + t713) * MDP(35) + (g(3) * t628 + t631 * t662) * MDP(10) + (t545 * t432 + t751 * t503 + t650 * t626 + t630 * t794 + t686) * MDP(34) + t640 + (-t519 * t750 + t545 * t636 + t546 * t679 - t683 * t749 + t688) * MDP(25) + qJDD(2) * MDP(8) + MDP(7) * t727 + MDP(6) * t728 + (t534 * t683 - t671 * t613 + t614 * t803 + t661) * MDP(24) + (-MDP(4) * t628 * t631 + MDP(5) * t742) * t635; (-t415 * t602 + t416 * t608 - t471 * t477 - g(1) * (t632 * t700 + t574) - g(2) * (t629 * t700 + t573) - g(3) * t724 + t745 * t463 + t701 * t462) * MDP(28) + (-t621 * t676 + t641) * MDP(16) + (-t477 * t519 + t701 * t614 + (-pkin(4) + t608) * t613 + t800) * MDP(26) + (-t477 * t683 + t602 * t613 - t614 * t745 + t654) * MDP(27) + (t469 * t614 + (-t519 * t553 + t613 * t786 - t614 * t738) * pkin(3) + t660) * MDP(23) + (t602 * t431 + t744 * t505 + t649 * t630 + (-t712 - t795) * t626 + t713) * MDP(35) + (t602 * t432 + t744 * t503 + t649 * t626 + t630 * t795 + t686) * MDP(34) + t640 + (-t519 * t745 + t602 * t636 + t608 * t679 - t683 * t701 + t688) * MDP(25) + (t470 * t614 + (-t553 * t683 - t613 * t627 - t614 * t718) * pkin(3) + t661) * MDP(24) + (t621 * t694 + t647) * MDP(17); t790 * MDP(19) + (t660 + t769) * MDP(23) + (-t467 * t614 + t677) * MDP(24) + (-pkin(4) * t679 + qJ(5) * t636 - (-t463 - t468) * t683) * MDP(25) + (t800 - t769 - 0.2e1 * t780) * MDP(26) + (-t485 * t683 + t614 * t732 + t595 + t667) * MDP(27) + (-t415 * qJ(5) - t416 * pkin(4) - t471 * t485 - t462 * t468 - g(1) * (-pkin(4) * t764 + t574) - g(2) * (-pkin(4) * t765 + t573) - g(3) * t743 - t732 * t463) * MDP(28) + t815 * MDP(31) + t814 * MDP(32) + (qJ(5) * t432 + t733 * t503 + (-t689 - t805) * t630 + t648 * t626 + t753) * MDP(34) + (qJ(5) * t431 + t411 + t733 * t505 + (t689 - t712) * t626 + t648 * t630) * MDP(35) - (t533 * MDP(24) + (t462 - t732) * MDP(25) + t485 * MDP(26) - t471 * MDP(27) + t505 * MDP(31) - t503 * MDP(32) + t424 * MDP(34) - t425 * MDP(35) + MDP(19) * t519 - t796) * t519 + t813; t664 * MDP(25) + (-t519 * t683 + t613) * MDP(26) + (-t614 ^ 2 - t790) * MDP(27) + (t463 * t614 - t780 + t800) * MDP(28) + (-t503 * t614 - t447) * MDP(34) + (-t505 * t614 + t771) * MDP(35) + (-MDP(34) * t811 - MDP(35) * t810) * t801; t505 * t503 * MDP(29) + (-t503 ^ 2 + t505 ^ 2) * MDP(30) + (t725 + t812) * MDP(31) + (t505 * t801 - t449) * MDP(32) - t451 * MDP(33) + (-g(1) * t540 - g(2) * t542 + t425 * t801 - t438 * t505 + t409) * MDP(34) + (g(1) * t541 - g(2) * t543 + t424 * t801 + t438 * t503) * MDP(35) + (-MDP(32) * t736 + MDP(34) * t692 - MDP(35) * t707) * t630 + (-MDP(31) * t736 + (qJD(6) * t519 - t613) * MDP(32) - t707 * MDP(34) + (-t412 - t692) * MDP(35)) * t626;];
tau  = t1;
