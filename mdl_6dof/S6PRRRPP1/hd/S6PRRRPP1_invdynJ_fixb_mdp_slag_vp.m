% Calculate vector of inverse dynamics joint torques for
% S6PRRRPP1
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
%   pkin=[a2,a3,a4,a5,a6,alpha2,d2,d3,d4,theta1,theta5]';
% MDP [24x1]
%   Minimal dynamic parameter vector (fixed base model)
%   see S6PRRRPP1_convert_par2_MPV_fixb.m
% 
% Output:
% tau [6x1]
%   joint torques of inverse dynamics (contains inertial, gravitational Coriolis and centrifugal forces)

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-08 22:47
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tau = S6PRRRPP1_invdynJ_fixb_mdp_slag_vp(qJ, qJD, qJDD, g, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(3,1),zeros(11,1),zeros(24,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRRRPP1_invdynJ_fixb_mdp_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6PRRRPP1_invdynJ_fixb_mdp_slag_vp: qJD has to be [6x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [6 1]), ...
  'S6PRRRPP1_invdynJ_fixb_mdp_slag_vp: qJDD has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6PRRRPP1_invdynJ_fixb_mdp_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6PRRRPP1_invdynJ_fixb_mdp_slag_vp: pkin has to be [11x1] (double)');
assert(isreal(MDP) && all(size(MDP) == [24 1]), ...
  'S6PRRRPP1_invdynJ_fixb_mdp_slag_vp: MDP has to be [24x1] (double)'); 

%% Symbolic Calculation
% From invdyn_joint_fixb_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-08 22:47:10
% EndTime: 2019-03-08 22:47:22
% DurationCPUTime: 9.76s
% Computational Cost: add. (6404->594), mult. (14621->781), div. (0->0), fcn. (11086->14), ass. (0->265)
t629 = sin(qJ(3));
t632 = cos(qJ(3));
t670 = pkin(3) * t629 - pkin(9) * t632;
t577 = t670 * qJD(3);
t628 = sin(qJ(4));
t630 = sin(qJ(2));
t631 = cos(qJ(4));
t714 = qJD(3) * t629;
t624 = sin(pkin(6));
t721 = qJD(1) * t624;
t633 = cos(qJ(2));
t742 = t632 * t633;
t770 = pkin(8) * t628;
t796 = (-t628 * t742 + t630 * t631) * t721 - t631 * t577 - t714 * t770;
t582 = -pkin(3) * t632 - pkin(9) * t629 - pkin(2);
t709 = qJD(4) * t631;
t745 = t628 * t630;
t795 = -(t631 * t742 + t745) * t721 + t628 * t577 + t582 * t709;
t717 = qJD(2) * t632;
t691 = t628 * t717;
t711 = qJD(4) * t628;
t794 = -t691 + t711;
t743 = t631 * t632;
t605 = pkin(8) * t743;
t772 = pkin(4) * t629;
t662 = -qJ(5) * t743 + t772;
t708 = qJD(5) * t631;
t793 = -t629 * t708 + t662 * qJD(3) + (-t605 + (qJ(5) * t629 - t582) * t628) * qJD(4) - t796;
t744 = t629 * t631;
t792 = (-pkin(8) * qJD(3) - qJ(5) * qJD(4)) * t744 + (-qJD(5) * t629 + (-pkin(8) * qJD(4) - qJ(5) * qJD(3)) * t632) * t628 + t795;
t769 = qJ(5) + pkin(9);
t791 = t769 * t629 + pkin(2);
t790 = pkin(4) * t628 + pkin(8);
t707 = t631 * qJD(3);
t718 = qJD(2) * t629;
t573 = t628 * t718 - t707;
t715 = qJD(3) * t628;
t575 = t631 * t718 + t715;
t623 = sin(pkin(11));
t625 = cos(pkin(11));
t491 = t625 * t573 + t575 * t623;
t602 = -qJD(4) + t717;
t789 = t491 * t602;
t569 = t623 * t631 + t625 * t628;
t553 = t569 * qJD(4);
t728 = t569 * t717 - t553;
t748 = t625 * t631;
t727 = t623 * t794 - t625 * t709 + t717 * t748;
t704 = qJD(2) * qJD(3);
t684 = t632 * t704;
t702 = qJDD(2) * t629;
t788 = t684 + t702;
t710 = qJD(4) * t629;
t787 = -qJD(2) * t710 + qJDD(3);
t611 = pkin(4) * t631 + pkin(3);
t786 = -t611 * t632 - t791;
t485 = t628 * (qJD(3) * (qJD(4) + t717) + t702) - t787 * t631;
t664 = -t573 * t623 + t625 * t575;
t785 = t664 ^ 2;
t738 = -t792 * t623 + t625 * t793;
t737 = t623 * t793 + t792 * t625;
t576 = t670 * qJD(2);
t560 = t631 * t576;
t578 = qJD(2) * pkin(8) + t630 * t721;
t626 = cos(pkin(6));
t720 = qJD(1) * t626;
t784 = -t629 * t578 + t632 * t720;
t463 = qJD(2) * t662 - t628 * t784 + t560;
t730 = t628 * t576 + t631 * t784;
t470 = -qJ(5) * t691 + t730;
t679 = qJD(4) * t769;
t546 = -t628 * t679 + t708;
t650 = -qJD(5) * t628 - t631 * t679;
t735 = (-t463 + t650) * t625 + (t470 - t546) * t623;
t713 = qJD(3) * t632;
t689 = t628 * t713;
t695 = t633 * t721;
t783 = pkin(4) * t689 + pkin(8) * t713 - t629 * t695 + t709 * t772;
t617 = t632 * qJDD(2);
t782 = -t629 * t704 + t617;
t598 = t629 * t720;
t527 = t632 * t578 + t598;
t781 = pkin(4) * t794 - t527;
t780 = MDP(19) + MDP(22);
t779 = pkin(4) * t485 + qJDD(5);
t516 = -qJD(3) * pkin(3) - t784;
t482 = pkin(4) * t573 + qJD(5) + t516;
t427 = pkin(5) * t491 - qJ(6) * t664 + t482;
t766 = sin(pkin(10));
t673 = t766 * t633;
t767 = cos(pkin(10));
t676 = t767 * t630;
t550 = t626 * t676 + t673;
t678 = t624 * t767;
t503 = t550 * t632 - t629 * t678;
t674 = t766 * t630;
t675 = t767 * t633;
t549 = -t626 * t675 + t674;
t620 = qJ(4) + pkin(11);
t613 = sin(t620);
t614 = cos(t620);
t458 = t503 * t613 - t549 * t614;
t552 = -t626 * t674 + t675;
t677 = t624 * t766;
t505 = t552 * t632 + t629 * t677;
t551 = t626 * t673 + t676;
t460 = t505 * t613 - t551 * t614;
t750 = t624 * t630;
t557 = t626 * t629 + t632 * t750;
t749 = t624 * t633;
t495 = t557 * t613 + t614 * t749;
t484 = qJD(4) * t707 + t628 * t787 + t631 * t788;
t567 = qJDD(4) - t782;
t705 = qJD(1) * qJD(2);
t535 = qJDD(2) * pkin(8) + (qJDD(1) * t630 + t633 * t705) * t624;
t703 = qJDD(1) * t626;
t682 = t629 * t703;
t455 = qJDD(3) * pkin(9) + qJD(3) * t784 + t535 * t632 + t682;
t517 = qJD(3) * pkin(9) + t527;
t528 = qJD(2) * t582 - t695;
t467 = t517 * t631 + t528 * t628;
t686 = t630 * t705;
t666 = -qJDD(1) * t749 + t624 * t686;
t477 = qJD(2) * t577 + qJDD(2) * t582 + t666;
t472 = t631 * t477;
t640 = -qJD(4) * t467 - t628 * t455 + t472;
t404 = pkin(4) * t567 - qJ(5) * t484 - qJD(5) * t575 + t640;
t698 = t631 * t455 + t628 * t477 + t528 * t709;
t655 = t517 * t711 - t698;
t407 = -qJ(5) * t485 - qJD(5) * t573 - t655;
t399 = t625 * t404 - t623 * t407;
t681 = -qJDD(6) + t399;
t778 = g(1) * t460 + g(2) * t458 + g(3) * t495 - t427 * t664 + t681;
t634 = qJD(3) ^ 2;
t669 = g(1) * t551 + g(2) * t549;
t777 = 0.2e1 * qJDD(2) * pkin(2) - pkin(8) * t634 + t624 * (-g(3) * t633 + t686) - t666 + t669;
t649 = g(3) * t749 - t669;
t775 = qJD(4) * (pkin(8) * t602 + t517) - t649;
t771 = pkin(5) * t567;
t768 = qJD(2) * pkin(2);
t764 = qJDD(3) * pkin(3);
t449 = -qJ(5) * t573 + t467;
t443 = t625 * t449;
t465 = -t517 * t628 + t631 * t528;
t448 = -qJ(5) * t575 + t465;
t419 = t448 * t623 + t443;
t763 = t419 * t664;
t762 = t449 * t623;
t761 = t484 * t628;
t760 = t503 * t628;
t759 = t505 * t628;
t758 = t549 * t631;
t757 = t551 * t631;
t756 = t573 * t602;
t755 = t575 * t602;
t754 = t575 * t631;
t752 = t613 * t632;
t751 = t614 * t632;
t746 = t628 * t629;
t741 = qJDD(1) - g(3);
t400 = t623 * t404 + t625 * t407;
t740 = qJ(6) * t714 - qJD(6) * t632 + t737;
t739 = -pkin(5) * t714 - t738;
t435 = -pkin(4) * t602 + t448;
t414 = t623 * t435 + t443;
t736 = -pkin(5) * t728 + qJ(6) * t727 - qJD(6) * t569 + t781;
t426 = t623 * t463 + t625 * t470;
t734 = pkin(5) * t718 - t735;
t423 = qJ(6) * t718 + t426;
t481 = t625 * t546 + t623 * t650;
t733 = t481 - t423;
t502 = t550 * t629 + t632 * t678;
t732 = -t502 * t611 + t503 * t769;
t504 = t552 * t629 - t632 * t677;
t731 = -t504 * t611 + t505 * t769;
t571 = t631 * t582;
t500 = -qJ(5) * t744 + t571 + (-pkin(4) - t770) * t632;
t724 = t628 * t582 + t605;
t510 = -qJ(5) * t746 + t724;
t451 = t623 * t500 + t625 * t510;
t556 = -t626 * t632 + t629 * t750;
t729 = -t556 * t611 + t557 * t769;
t723 = pkin(4) * t746 + t629 * pkin(8);
t621 = t629 ^ 2;
t722 = -t632 ^ 2 + t621;
t719 = qJD(2) * t624;
t716 = qJD(3) * t573;
t712 = qJD(4) * t602;
t420 = t448 * t625 - t762;
t706 = qJD(6) - t420;
t700 = t631 * t749;
t699 = t567 * qJ(6) + t400;
t693 = t630 * t719;
t692 = t633 * t719;
t690 = t602 * t707;
t688 = t632 * t707;
t687 = t769 * t628;
t432 = t484 * t623 + t625 * t485;
t668 = g(1) * t552 + g(2) * t550;
t667 = -pkin(5) * t614 - qJ(6) * t613;
t413 = t435 * t625 - t762;
t433 = t484 * t625 - t485 * t623;
t450 = t500 * t625 - t510 * t623;
t568 = t623 * t628 - t748;
t478 = t568 * t710 - t569 * t713;
t479 = t553 * t629 + t623 * t689 - t625 * t688;
t543 = -t623 * t746 + t625 * t744;
t663 = pkin(5) * t478 - qJ(6) * t479 + qJD(6) * t543 - t783;
t508 = -t557 * t628 - t700;
t659 = -t557 * t631 + t628 * t749;
t657 = t567 * t628 - t602 * t709;
t656 = t567 * t631 + t602 * t711;
t654 = g(1) * t504 + g(2) * t502 + g(3) * t556;
t653 = g(1) * t505 + g(2) * t503 + g(3) * t557;
t652 = -qJD(3) * t598 - t629 * t535 - t578 * t713 + t632 * t703;
t648 = -g(3) * t750 - t668;
t647 = -pkin(9) * t567 - t516 * t602;
t646 = pkin(8) * t750 + t791 * t749 + (pkin(4) * t745 + t611 * t742) * t624;
t645 = t549 * t786 + t550 * t790;
t644 = t551 * t786 + t552 * t790;
t643 = t649 * t629;
t456 = -t652 - t764;
t639 = pkin(9) * t712 - t456 + t654;
t584 = t769 * t631;
t513 = t584 * t623 + t625 * t687;
t514 = t625 * t584 - t623 * t687;
t638 = -t514 * t432 + t433 * t513 - t481 * t491 - t653;
t579 = -t695 - t768;
t637 = -pkin(8) * qJDD(3) + (t579 + t695 - t768) * qJD(3);
t429 = t456 + t779;
t636 = t652 + t654;
t401 = pkin(5) * t432 - qJ(6) * t433 - qJD(6) * t664 + t429;
t635 = qJD(2) ^ 2;
t609 = -pkin(4) * t625 - pkin(5);
t606 = pkin(4) * t623 + qJ(6);
t542 = t569 * t629;
t539 = pkin(4) * t757;
t537 = pkin(4) * t758;
t520 = (t613 * t630 + t614 * t742) * t624;
t519 = (t613 * t742 - t614 * t630) * t624;
t507 = qJD(3) * t557 + t629 * t692;
t506 = -qJD(3) * t556 + t632 * t692;
t496 = t557 * t614 - t613 * t749;
t488 = pkin(5) * t568 - qJ(6) * t569 - t611;
t476 = -t551 * t751 + t552 * t613;
t475 = -t551 * t752 - t552 * t614;
t474 = -t549 * t751 + t550 * t613;
t473 = -t549 * t752 - t550 * t614;
t469 = pkin(5) * t542 - qJ(6) * t543 + t723;
t461 = t505 * t614 + t551 * t613;
t459 = t503 * t614 + t549 * t613;
t453 = t508 * t623 - t625 * t659;
t452 = -t625 * t508 - t623 * t659;
t446 = pkin(5) * t632 - t450;
t445 = -qJ(6) * t632 + t451;
t442 = qJD(4) * t508 + t506 * t631 + t628 * t693;
t441 = qJD(4) * t659 - t506 * t628 + t631 * t693;
t430 = pkin(4) * t575 + pkin(5) * t664 + qJ(6) * t491;
t418 = t441 * t623 + t442 * t625;
t416 = -t625 * t441 + t442 * t623;
t411 = -qJ(6) * t602 + t414;
t409 = pkin(5) * t602 + qJD(6) - t413;
t398 = -t681 - t771;
t397 = -qJD(6) * t602 + t699;
t1 = [t741 * MDP(1) + (-qJD(3) * t507 - qJDD(3) * t556) * MDP(10) + (-qJD(3) * t506 - qJDD(3) * t557) * MDP(11) + (-t441 * t602 + t485 * t556 + t507 * t573 + t508 * t567) * MDP(17) + (t442 * t602 + t484 * t556 + t507 * t575 + t567 * t659) * MDP(18) + (-t399 * t452 + t400 * t453 - t413 * t416 + t414 * t418 + t429 * t556 + t482 * t507 - g(3)) * MDP(20) + (t416 * t602 + t432 * t556 - t452 * t567 + t491 * t507) * MDP(21) + (-t418 * t602 - t433 * t556 + t453 * t567 - t507 * t664) * MDP(23) + (t397 * t453 + t398 * t452 + t401 * t556 + t409 * t416 + t411 * t418 + t427 * t507 - g(3)) * MDP(24) + ((-qJDD(2) * MDP(4) + (-MDP(10) * t632 + MDP(11) * t629 - MDP(3)) * t635) * t630 + (t782 * MDP(10) - MDP(11) * t788 + qJDD(2) * MDP(3) - t635 * MDP(4)) * t633) * t624 + t780 * (t416 * t664 - t418 * t491 - t453 * t432 + t433 * t452); qJDD(2) * MDP(2) + (t741 * t749 + t669) * MDP(3) + (-t741 * t750 + t668) * MDP(4) + (qJDD(2) * t621 + 0.2e1 * t629 * t684) * MDP(5) + 0.2e1 * (t617 * t629 - t704 * t722) * MDP(6) + (qJDD(3) * t629 + t632 * t634) * MDP(7) + (qJDD(3) * t632 - t629 * t634) * MDP(8) + (t637 * t629 + t632 * t777) * MDP(10) + (-t629 * t777 + t637 * t632) * MDP(11) + (t484 * t744 + (-t628 * t710 + t688) * t575) * MDP(12) + ((-t573 * t631 - t575 * t628) * t713 + (-t761 - t485 * t631 + (t573 * t628 - t754) * qJD(4)) * t629) * MDP(13) + ((-t484 - t690) * t632 + (qJD(3) * t575 + t656) * t629) * MDP(14) + ((t602 * t715 + t485) * t632 + (-t657 - t716) * t629) * MDP(15) + (-t567 * t632 - t602 * t714) * MDP(16) + (t571 * t567 + t796 * t602 + (t582 * t712 + t648) * t628 + (pkin(8) * t716 - t472 + (-pkin(8) * t567 + qJD(3) * t516 + qJD(4) * t528 + t455) * t628 + t775 * t631) * t632 + (pkin(8) * t485 + qJD(3) * t465 + t456 * t628 + t516 * t709 - t573 * t695) * t629) * MDP(17) + (-t724 * t567 + t795 * t602 + t648 * t631 + ((pkin(8) * t575 + t516 * t631) * qJD(3) - t775 * t628 + t698) * t632 + (-t575 * t695 - t516 * t711 - t467 * qJD(3) + t456 * t631 + (t484 - t690) * pkin(8)) * t629) * MDP(18) + (-t399 * t543 - t400 * t542 + t413 * t479 + t414 * t478 - t432 * t451 - t433 * t450 - t491 * t737 - t664 * t738 - t643) * MDP(19) + (-g(1) * t644 - g(2) * t645 - g(3) * t646 + t399 * t450 + t400 * t451 + t738 * t413 + t737 * t414 + t429 * t723 + t482 * t783) * MDP(20) + (-g(1) * t476 - g(2) * t474 - g(3) * t520 + t398 * t632 + t401 * t542 - t409 * t714 - t427 * t478 + t432 * t469 - t446 * t567 - t491 * t663 + t602 * t739) * MDP(21) + (-t397 * t542 + t398 * t543 - t409 * t479 + t411 * t478 - t432 * t445 + t433 * t446 - t491 * t740 + t664 * t739 - t643) * MDP(22) + (-g(1) * t475 - g(2) * t473 - g(3) * t519 - t397 * t632 - t401 * t543 + t411 * t714 + t427 * t479 - t433 * t469 + t445 * t567 - t602 * t740 + t663 * t664) * MDP(23) + (t397 * t445 + t401 * t469 + t398 * t446 - g(1) * (pkin(5) * t476 + qJ(6) * t475 + t644) - g(2) * (pkin(5) * t474 + qJ(6) * t473 + t645) - g(3) * (pkin(5) * t520 + qJ(6) * t519 + t646) - t663 * t427 + t740 * t411 + t739 * t409) * MDP(24); MDP(7) * t702 + MDP(8) * t617 + qJDD(3) * MDP(9) + (qJD(3) * t527 - t579 * t718 + t636) * MDP(10) + (-t682 + (-qJD(2) * t579 - t535) * t632 + t653) * MDP(11) + (-t602 * t754 + t761) * MDP(12) + ((t484 + t756) * t631 + (-t485 + t755) * t628) * MDP(13) + ((-t575 * t629 + t602 * t743) * qJD(2) + t657) * MDP(14) + ((-t602 * t628 * t632 + t573 * t629) * qJD(2) + t656) * MDP(15) + t602 * MDP(16) * t718 + (-t465 * t718 - pkin(3) * t485 - t527 * t573 + t560 * t602 + (-t602 * t784 + t647) * t628 + t639 * t631) * MDP(17) + (-pkin(3) * t484 + t467 * t718 - t527 * t575 - t602 * t730 - t628 * t639 + t631 * t647) * MDP(18) + (-t399 * t569 - t400 * t568 + t413 * t727 + t414 * t728 + t426 * t491 - t664 * t735 + t638) * MDP(19) + (t400 * t514 - t399 * t513 - t429 * t611 - g(1) * t731 - g(2) * t732 - g(3) * t729 + t781 * t482 + (t481 - t426) * t414 + t735 * t413) * MDP(20) + (t401 * t568 + t409 * t718 - t427 * t728 + t432 * t488 + t491 * t736 - t513 * t567 + t602 * t734 + t614 * t654) * MDP(21) + (-t397 * t568 + t398 * t569 - t409 * t727 + t411 * t728 + t423 * t491 + t664 * t734 + t638) * MDP(22) + (-t401 * t569 - t411 * t718 + t427 * t727 - t433 * t488 + t514 * t567 - t602 * t733 + t613 * t654 - t664 * t736) * MDP(23) + (t397 * t514 + t401 * t488 + t398 * t513 - g(1) * (t504 * t667 + t731) - g(2) * (t502 * t667 + t732) - g(3) * (t556 * t667 + t729) + t736 * t427 + t733 * t411 + t734 * t409) * MDP(24) + (-MDP(5) * t629 * t632 + MDP(6) * t722) * t635; t575 * t573 * MDP(12) + (-t573 ^ 2 + t575 ^ 2) * MDP(13) + (t484 - t756) * MDP(14) + (-t485 - t755) * MDP(15) + t567 * MDP(16) + (-t467 * t602 - t516 * t575 - g(1) * (t757 - t759) - g(2) * (t758 - t760) - g(3) * t508 + t640) * MDP(17) + (-t465 * t602 + t516 * t573 - g(1) * (-t505 * t631 - t551 * t628) - g(2) * (-t503 * t631 - t549 * t628) - g(3) * t659 + t655) * MDP(18) + (t414 * t664 - t763 + (-t432 * t623 - t433 * t625) * pkin(4) + (-t413 + t420) * t491) * MDP(19) + (-g(1) * t539 - g(2) * t537 + t413 * t419 - t414 * t420 + (g(3) * t700 + t399 * t625 + t400 * t623 - t482 * t575 + t628 * t653) * pkin(4)) * MDP(20) + (-t419 * t602 - t430 * t491 + (pkin(5) - t609) * t567 + t778) * MDP(21) + (t411 * t664 - t432 * t606 + t433 * t609 - t763 + (t409 - t706) * t491) * MDP(22) + (-g(1) * t461 - g(2) * t459 - g(3) * t496 - t427 * t491 + t430 * t664 + t567 * t606 + (-0.2e1 * qJD(6) + t420) * t602 + t699) * MDP(23) + (t397 * t606 + t398 * t609 - t427 * t430 - t409 * t419 - g(1) * (-pkin(4) * t759 - pkin(5) * t460 + qJ(6) * t461 + t539) - g(2) * (-pkin(4) * t760 - pkin(5) * t458 + qJ(6) * t459 + t537) - g(3) * (pkin(4) * t508 - pkin(5) * t495 + qJ(6) * t496) + t706 * t411) * MDP(24); (t413 * t664 + t414 * t491 - t636 - t764 + t779) * MDP(20) + (-t602 * t664 + t432) * MDP(21) + (-t433 - t789) * MDP(23) + (-t409 * t664 + t411 * t491 + t401 - t654) * MDP(24) + t780 * (-t491 ^ 2 - t785); (t491 * t664 - t567) * MDP(21) + (t433 - t789) * MDP(22) + (-t602 ^ 2 - t785) * MDP(23) + (t411 * t602 - t771 - t778) * MDP(24);];
tau  = t1;
