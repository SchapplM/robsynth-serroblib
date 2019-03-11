% Calculate vector of inverse dynamics joint torques for
% S6RRPRPR1
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
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d4,d6,theta3,theta5]';
% MDP [30x1]
%   Minimal dynamic parameter vector (fixed base model)
%   see S6RRPRPR1_convert_par2_MPV_fixb.m
% 
% Output:
% tau [6x1]
%   joint torques of inverse dynamics (contains inertial, gravitational Coriolis and centrifugal forces)

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 10:11
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tau = S6RRPRPR1_invdynJ_fixb_mdp_slag_vp(qJ, qJD, qJDD, g, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(3,1),zeros(11,1),zeros(30,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPRPR1_invdynJ_fixb_mdp_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRPRPR1_invdynJ_fixb_mdp_slag_vp: qJD has to be [6x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [6 1]), ...
  'S6RRPRPR1_invdynJ_fixb_mdp_slag_vp: qJDD has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RRPRPR1_invdynJ_fixb_mdp_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6RRPRPR1_invdynJ_fixb_mdp_slag_vp: pkin has to be [11x1] (double)');
assert(isreal(MDP) && all(size(MDP) == [30 1]), ...
  'S6RRPRPR1_invdynJ_fixb_mdp_slag_vp: MDP has to be [30x1] (double)'); 

%% Symbolic Calculation
% From invdyn_joint_fixb_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 10:10:16
% EndTime: 2019-03-09 10:10:33
% DurationCPUTime: 13.01s
% Computational Cost: add. (11293->573), mult. (27033->736), div. (0->0), fcn. (21228->18), ass. (0->268)
t683 = sin(pkin(10));
t685 = cos(pkin(10));
t689 = sin(qJ(2));
t692 = cos(qJ(2));
t622 = -t683 * t689 + t685 * t692;
t608 = t622 * qJD(1);
t624 = t683 * t692 + t685 * t689;
t610 = t624 * qJD(1);
t688 = sin(qJ(4));
t807 = cos(qJ(4));
t564 = t807 * t608 - t610 * t688;
t558 = qJD(6) - t564;
t682 = sin(pkin(11));
t684 = cos(pkin(11));
t687 = sin(qJ(6));
t691 = cos(qJ(6));
t623 = t682 * t687 - t691 * t684;
t839 = t558 * t623;
t752 = qJD(1) * qJD(2);
t743 = t692 * t752;
t744 = t689 * t752;
t573 = qJDD(1) * t624 - t683 * t744 + t685 * t743;
t609 = t624 * qJD(2);
t706 = qJD(1) * t609;
t696 = qJDD(1) * t622 - t706;
t716 = -t688 * t608 - t610 * t807;
t699 = qJD(4) * t716 - t688 * t573 + t807 * t696;
t493 = qJDD(6) - t699;
t625 = t682 * t691 + t684 * t687;
t842 = -t625 * t493 + t558 * t839;
t678 = qJD(2) + qJD(4);
t547 = -t684 * t678 - t682 * t716;
t549 = t678 * t682 - t684 * t716;
t787 = t549 * t687;
t834 = -t691 * t547 - t787;
t841 = t558 * t834;
t614 = t625 * qJD(6);
t840 = -t625 * t564 + t614;
t745 = qJD(4) * t807;
t757 = qJD(4) * t688;
t494 = t807 * t573 + t608 * t745 - t610 * t757 + t688 * t696;
t675 = qJDD(2) + qJDD(4);
t482 = t494 * t682 - t684 * t675;
t483 = t494 * t684 + t675 * t682;
t755 = qJD(6) * t691;
t748 = -t687 * t482 + t691 * t483 - t547 * t755;
t756 = qJD(6) * t687;
t438 = -t549 * t756 + t748;
t723 = t547 * t687 - t549 * t691;
t738 = t691 * t482 + t687 * t483;
t439 = -qJD(6) * t723 + t738;
t731 = -t623 * t493 - t558 * t840;
t784 = t564 * t678;
t786 = t716 * t678;
t789 = t723 * t716;
t790 = t834 * t716;
t838 = (t731 + t790) * MDP(27) - t564 ^ 2 * MDP(14) + t675 * MDP(17) + (MDP(13) * t564 + MDP(14) * t716 + MDP(28) * t558) * t716 + (t494 - t784) * MDP(15) + (t699 - t786) * MDP(16) + t438 * t625 * MDP(24) + (-t789 - t842) * MDP(26) + (-t438 * t623 - t625 * t439 - t839 * t834) * MDP(25) + (MDP(24) * t839 + MDP(25) * t840) * t723;
t827 = t564 * t682;
t837 = pkin(5) * t827;
t836 = pkin(9) * t827;
t673 = t692 * pkin(2);
t666 = t673 + pkin(1);
t593 = -pkin(3) * t622 - t666;
t794 = qJ(3) + pkin(7);
t648 = t794 * t692;
t634 = qJD(1) * t648;
t615 = t683 * t634;
t647 = t794 * t689;
t633 = qJD(1) * t647;
t793 = qJD(2) * pkin(2);
t619 = -t633 + t793;
t571 = t685 * t619 - t615;
t802 = pkin(8) * t610;
t542 = qJD(2) * pkin(3) + t571 - t802;
t768 = t685 * t634;
t572 = t683 * t619 + t768;
t803 = pkin(8) * t608;
t546 = t572 + t803;
t499 = t688 * t542 + t807 * t546;
t490 = qJ(5) * t678 + t499;
t637 = -qJD(1) * t666 + qJD(3);
t579 = -pkin(3) * t608 + t637;
t504 = -pkin(4) * t564 + qJ(5) * t716 + t579;
t457 = -t490 * t682 + t684 * t504;
t447 = -pkin(5) * t564 - pkin(9) * t549 + t457;
t458 = t684 * t490 + t682 * t504;
t450 = -pkin(9) * t547 + t458;
t427 = t447 * t687 + t450 * t691;
t740 = qJD(2) * t794;
t606 = -qJD(3) * t689 - t692 * t740;
t570 = qJDD(2) * pkin(2) + qJD(1) * t606 - qJDD(1) * t647;
t605 = qJD(3) * t692 - t689 * t740;
t576 = qJD(1) * t605 + qJDD(1) * t648;
t528 = t685 * t570 - t576 * t683;
t497 = qJDD(2) * pkin(3) - pkin(8) * t573 + t528;
t529 = t683 * t570 + t685 * t576;
t502 = pkin(8) * t696 + t529;
t734 = -t807 * t497 + t688 * t502 + t542 * t757 + t546 * t745;
t808 = -pkin(4) * t675 + qJDD(5);
t443 = t734 + t808;
t432 = pkin(5) * t482 + t443;
t498 = t542 * t807 - t688 * t546;
t489 = -t678 * pkin(4) + qJD(5) - t498;
t474 = t547 * pkin(5) + t489;
t679 = qJ(2) + pkin(10);
t671 = qJ(4) + t679;
t662 = cos(t671);
t677 = pkin(11) + qJ(6);
t669 = sin(t677);
t661 = sin(t671);
t690 = sin(qJ(1));
t693 = cos(qJ(1));
t733 = g(1) * t693 + g(2) * t690;
t720 = t661 * t733;
t799 = g(3) * t669;
t833 = -t427 * t716 + t432 * t625 - t474 * t839 + t662 * t799 - t669 * t720;
t426 = t447 * t691 - t450 * t687;
t670 = cos(t677);
t798 = g(3) * t670;
t776 = t661 * t693;
t777 = t661 * t690;
t823 = g(1) * t776 + g(2) * t777;
t832 = t426 * t716 + t432 * t623 + t474 * t840 - t662 * t798 + t670 * t823;
t830 = pkin(5) * t716;
t829 = t457 * t564;
t828 = t558 * t723;
t762 = t662 * pkin(4) + t661 * qJ(5);
t822 = pkin(3) * cos(t679) + t673 + t762;
t821 = -g(3) * t662 + t823;
t698 = t688 * t497 + t502 * t807 + t542 * t745 - t546 * t757;
t441 = t675 * qJ(5) + t678 * qJD(5) + t698;
t749 = pkin(2) * t744 + qJDD(3);
t544 = pkin(3) * t706 + qJDD(1) * t593 + t749;
t446 = -pkin(4) * t699 - t494 * qJ(5) + qJD(5) * t716 + t544;
t424 = -t441 * t682 + t684 * t446;
t820 = t458 * t564 - t424;
t817 = t457 * t716 + (-t443 + t821) * t684;
t775 = t662 * t682;
t816 = -t458 * t716 + g(3) * t775 + (t443 - t720) * t682;
t773 = t662 * t693;
t774 = t662 * t690;
t747 = -g(1) * t773 - g(2) * t774 - g(3) * t661;
t815 = -t579 * t564 - t698 - t747;
t705 = -t734 + t821;
t814 = t579 * t716 + t705;
t525 = -pkin(4) * t716 - qJ(5) * t564;
t809 = g(1) * t690 - g(2) * t693;
t812 = t809 * t661;
t577 = t633 * t683 - t768;
t551 = t577 - t803;
t578 = -t685 * t633 - t615;
t552 = t578 - t802;
t511 = t688 * t551 + t552 * t807;
t795 = t689 * pkin(2);
t584 = pkin(3) * t610 + qJD(1) * t795;
t512 = t525 + t584;
t461 = -t511 * t682 + t684 * t512;
t806 = pkin(2) * t683;
t655 = t688 * t806;
t663 = pkin(2) * t685 + pkin(3);
t722 = -qJD(4) * t655 + t663 * t745;
t591 = qJD(5) + t722;
t811 = -t591 * t682 - t461;
t462 = t684 * t511 + t682 * t512;
t810 = t591 * t684 - t462;
t761 = t688 * t663 + t807 * t806;
t763 = t761 * qJD(4) + t551 * t807 - t688 * t552;
t797 = g(3) * t692;
t796 = t684 * pkin(5);
t672 = t684 * pkin(9);
t425 = t684 * t441 + t682 * t446;
t423 = t425 * t684;
t792 = t699 * t682;
t791 = t699 * t684;
t612 = t622 * qJD(2);
t715 = t622 * t807 - t688 * t624;
t530 = t715 * qJD(4) - t688 * t609 + t612 * t807;
t788 = t530 * t682;
t782 = t564 * t684;
t575 = t688 * t622 + t624 * t807;
t781 = t575 * t682;
t780 = t575 * t684;
t772 = t669 * t690;
t771 = t669 * t693;
t770 = t670 * t690;
t769 = t670 * t693;
t553 = -t605 * t683 + t685 * t606;
t533 = -pkin(8) * t612 + t553;
t554 = t685 * t605 + t683 * t606;
t534 = -pkin(8) * t609 + t554;
t582 = -t685 * t647 - t648 * t683;
t555 = -pkin(8) * t624 + t582;
t583 = -t683 * t647 + t685 * t648;
t556 = pkin(8) * t622 + t583;
t717 = t555 * t807 - t688 * t556;
t463 = t717 * qJD(4) + t688 * t533 + t534 * t807;
t531 = t575 * qJD(4) + t609 * t807 + t688 * t612;
t668 = t689 * t793;
t585 = pkin(3) * t609 + t668;
t467 = pkin(4) * t531 - qJ(5) * t530 - qJD(5) * t575 + t585;
t434 = t684 * t463 + t682 * t467;
t469 = t684 * t498 + t682 * t525;
t520 = -pkin(4) * t715 - qJ(5) * t575 + t593;
t522 = t688 * t555 + t556 * t807;
t473 = t682 * t520 + t684 * t522;
t764 = t763 - t837;
t680 = t689 ^ 2;
t759 = -t692 ^ 2 + t680;
t753 = -qJD(5) + t489;
t751 = qJDD(1) * t689;
t750 = qJDD(1) * t692;
t742 = -pkin(4) * t661 - pkin(3) * sin(t679) - t795;
t421 = -pkin(5) * t699 - pkin(9) * t483 + t424;
t422 = -pkin(9) * t482 + t425;
t739 = t691 * t421 - t422 * t687;
t433 = -t463 * t682 + t684 * t467;
t468 = -t498 * t682 + t684 * t525;
t472 = t684 * t520 - t522 * t682;
t736 = t423 + t747;
t730 = t663 * t807 - t655;
t729 = t421 * t687 + t422 * t691;
t728 = -t424 * t684 - t425 * t682;
t727 = -t424 * t682 + t423;
t453 = -pkin(5) * t715 - pkin(9) * t780 + t472;
t459 = -pkin(9) * t781 + t473;
t726 = t453 * t691 - t459 * t687;
t725 = t453 * t687 + t459 * t691;
t724 = t457 * t682 - t458 * t684;
t604 = -pkin(4) - t730;
t721 = pkin(1) + t822;
t719 = t809 * t662;
t718 = -0.2e1 * pkin(1) * t752 - pkin(7) * qJDD(2);
t603 = qJ(5) + t761;
t581 = t603 * t684 + t672;
t714 = -pkin(9) * t782 + qJD(6) * t581 - t811 - t830;
t580 = (-pkin(9) - t603) * t682;
t713 = -qJD(6) * t580 - t810 - t836;
t646 = qJ(5) * t684 + t672;
t712 = qJD(5) * t682 + qJD(6) * t646 - t564 * t672 + t468 - t830;
t645 = (-pkin(9) - qJ(5)) * t682;
t711 = -qJD(5) * t684 - qJD(6) * t645 + t469 - t836;
t708 = -qJDD(1) * t666 + t749;
t707 = t699 * t603 + (-t489 + t591) * t564;
t694 = qJD(2) ^ 2;
t704 = 0.2e1 * qJDD(1) * pkin(1) - pkin(7) * t694 + t809;
t695 = qJD(1) ^ 2;
t703 = pkin(1) * t695 - pkin(7) * qJDD(1) + t733;
t701 = t443 * t575 + t489 * t530 - t733;
t464 = t522 * qJD(4) - t533 * t807 + t688 * t534;
t676 = -pkin(8) - t794;
t664 = -pkin(4) - t796;
t639 = qJ(5) * t773;
t638 = qJ(5) * t774;
t590 = t604 - t796;
t589 = t662 * t769 + t772;
t588 = -t662 * t771 + t770;
t587 = -t662 * t770 + t771;
t586 = t662 * t772 + t769;
t527 = t623 * t575;
t526 = t625 * t575;
t484 = pkin(5) * t781 - t717;
t475 = t499 + t837;
t456 = t530 * t625 + t755 * t780 - t756 * t781;
t455 = -t530 * t623 - t575 * t614;
t449 = pkin(5) * t788 + t464;
t429 = -pkin(9) * t788 + t434;
t428 = pkin(5) * t531 - t530 * t672 + t433;
t1 = [(t424 * t472 + t425 * t473 + t457 * t433 + t458 * t434 - t443 * t717 + t489 * t464 + (g(1) * t676 - g(2) * t721) * t693 + (g(1) * t721 + g(2) * t676) * t690) * MDP(23) + ((t428 * t691 - t429 * t687) * t558 + t726 * t493 - t739 * t715 + t426 * t531 - t449 * t834 + t484 * t439 + t432 * t526 + t474 * t456 - g(1) * t587 - g(2) * t589 + (t427 * t715 - t558 * t725) * qJD(6)) * MDP(29) + (t439 * t715 - t456 * t558 - t493 * t526 + t531 * t834) * MDP(27) + (-t438 * t526 + t439 * t527 + t455 * t834 + t456 * t723) * MDP(25) + (-t424 * t715 - t433 * t564 + t457 * t531 + t464 * t547 - t472 * t699 - t482 * t717 + t682 * t701 + t684 * t719) * MDP(20) + (-t464 * t678 + t531 * t579 - t544 * t715 - t564 * t585 - t593 * t699 + t675 * t717 + t719) * MDP(18) + (t494 * t715 + t530 * t564 + t531 * t716 + t575 * t699) * MDP(14) + (t425 * t715 + t434 * t564 - t458 * t531 + t464 * t549 + t473 * t699 - t483 * t717 + t684 * t701 - t775 * t809) * MDP(21) + (t689 * t718 + t692 * t704) * MDP(9) + (-t689 * t704 + t692 * t718) * MDP(10) + (-t463 * t678 + t494 * t593 - t522 * t675 + t530 * t579 + t544 * t575 - t585 * t716 - t812) * MDP(19) + (t494 * t575 - t530 * t716) * MDP(13) + (-(t428 * t687 + t429 * t691) * t558 - t725 * t493 + t729 * t715 - t427 * t531 - t449 * t723 + t484 * t438 - t432 * t527 + t474 * t455 - g(1) * t586 - g(2) * t588 + (t426 * t715 - t558 * t726) * qJD(6)) * MDP(30) + (-t438 * t715 + t455 * t558 - t493 * t527 - t531 * t723) * MDP(26) + (-t531 * t678 + t675 * t715) * MDP(16) + (-t493 * t715 + t531 * t558) * MDP(28) + (t529 * t583 + t572 * t554 + t528 * t582 + t571 * t553 - t708 * t666 + t637 * t668 - g(1) * (-t666 * t690 + t693 * t794) - g(2) * (t666 * t693 + t690 * t794)) * MDP(12) + t809 * MDP(2) + (qJDD(1) * t680 + 0.2e1 * t689 * t743) * MDP(4) + 0.2e1 * (t689 * t750 - t752 * t759) * MDP(5) + (-t433 * t549 - t434 * t547 - t472 * t483 - t473 * t482 + t812 + t728 * t575 + (-t457 * t684 - t458 * t682) * t530) * MDP(22) + qJDD(1) * MDP(1) + t733 * MDP(3) + (-t438 * t527 - t455 * t723) * MDP(24) + (qJDD(2) * t689 + t692 * t694) * MDP(6) + (qJDD(2) * t692 - t689 * t694) * MDP(7) + (t530 * t678 + t575 * t675) * MDP(15) + (-t528 * t624 + t529 * t622 - t553 * t610 + t554 * t608 - t571 * t612 - t572 * t609 - t582 * t573 + t583 * t696 - t733) * MDP(11); (t689 * t703 - t797) * MDP(9) + (-t571 * t577 - t572 * t578 + (-t797 + t528 * t685 + t529 * t683 + (-qJD(1) * t637 + t733) * t689) * pkin(2)) * MDP(12) + ((t580 * t691 - t581 * t687) * t493 + t590 * t439 + (t687 * t713 - t691 * t714) * t558 - t764 * t834 + t832) * MDP(29) + (-(t580 * t687 + t581 * t691) * t493 + t590 * t438 + (t687 * t714 + t691 * t713) * t558 - t764 * t723 + t833) * MDP(30) + (g(3) * t689 + t692 * t703) * MDP(10) + (-MDP(4) * t689 * t692 + MDP(5) * t759) * t695 + (-t462 * t564 + t483 * t604 + t549 * t763 + t684 * t707 + t816) * MDP(21) + (t461 * t564 + t482 * t604 + t547 * t763 + t682 * t707 + t817) * MDP(20) + (t564 * t584 + t675 * t730 - t678 * t763 + t814) * MDP(18) + (-t761 * t675 + t584 * t716 + (t511 - t722) * t678 + t815) * MDP(19) + (t443 * t604 - g(1) * (t693 * t742 + t639) - g(2) * (t690 * t742 + t638) - g(3) * t822 + t727 * t603 + t763 * t489 + t810 * t458 + t811 * t457) * MDP(23) + (t461 * t549 + t462 * t547 + (-t482 * t603 - t547 * t591 + t829) * t684 + (t483 * t603 + t549 * t591 + t820) * t682 + t736) * MDP(22) + qJDD(2) * MDP(8) + ((t572 + t577) * t610 + (-t578 + t571) * t608 + (-t685 * t573 + ((-t743 - t751) * t683 + (-t744 + t750) * t685) * t683) * pkin(2)) * MDP(11) + MDP(7) * t750 + MDP(6) * t751 + t838; (-t608 ^ 2 - t610 ^ 2) * MDP(11) + (t571 * t610 - t572 * t608 + t708 - t809) * MDP(12) + (-t699 - t786) * MDP(18) + (t494 + t784) * MDP(19) + (t547 * t716 - t564 * t827 - t791) * MDP(20) + (t549 * t716 - t564 * t782 + t792) * MDP(21) + (-t482 * t682 - t483 * t684 + (t547 * t684 - t549 * t682) * t564) * MDP(22) + (t489 * t716 + t564 * t724 - t728 - t809) * MDP(23) + (t731 - t790) * MDP(29) + (-t789 + t842) * MDP(30); (t499 * t678 + t814) * MDP(18) + (t498 * t678 + t815) * MDP(19) + (qJ(5) * t792 - pkin(4) * t482 - t499 * t547 - (t682 * t753 - t468) * t564 + t817) * MDP(20) + (qJ(5) * t791 - pkin(4) * t483 - t499 * t549 - (t684 * t753 + t469) * t564 + t816) * MDP(21) + (t468 * t549 + t469 * t547 + (-qJ(5) * t482 - qJD(5) * t547 + t829) * t684 + (qJ(5) * t483 + qJD(5) * t549 + t820) * t682 + t736) * MDP(22) + (-t443 * pkin(4) - t458 * t469 - t457 * t468 - t489 * t499 - g(1) * (-pkin(4) * t776 + t639) - g(2) * (-pkin(4) * t777 + t638) - g(3) * t762 - t724 * qJD(5) + t727 * qJ(5)) * MDP(23) + ((t645 * t691 - t646 * t687) * t493 + t664 * t439 + t475 * t834 + (t687 * t711 - t691 * t712) * t558 + t832) * MDP(29) + (-(t645 * t687 + t646 * t691) * t493 + t664 * t438 + t475 * t723 + (t687 * t712 + t691 * t711) * t558 + t833) * MDP(30) + t838; (-t549 * t564 + t482) * MDP(20) + (t547 * t564 + t483) * MDP(21) + (-t547 ^ 2 - t549 ^ 2) * MDP(22) + (t457 * t549 + t458 * t547 - t705 + t808) * MDP(23) + (t439 - t828) * MDP(29) + (t438 + t841) * MDP(30); t723 * t834 * MDP(24) + (t723 ^ 2 - t834 ^ 2) * MDP(25) + (t748 - t841) * MDP(26) + (-t738 - t828) * MDP(27) + t493 * MDP(28) + (-g(1) * t588 + g(2) * t586 + t427 * t558 + t474 * t723 + t661 * t799 + t739) * MDP(29) + (g(1) * t589 - g(2) * t587 + t426 * t558 - t474 * t834 + t661 * t798 - t729) * MDP(30) + (-MDP(26) * t787 + MDP(27) * t723 - MDP(29) * t427 - MDP(30) * t426) * qJD(6);];
tau  = t1;
