% Calculate minimal parameter regressor of Coriolis joint torque vector for
% S6RRPRRR8
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
% MDP [35x1]
%   Minimal dynamic parameter vector (fixed base model)
%   see S6RRPRRR8_convert_par2_MPV_fixb.m
% 
% Output:
% tauc [6x1]
%   joint torques required to compensate Coriolis and centrifugal load

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 14:07
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauc = S6RRPRRR8_coriolisvecJ_fixb_mdp_slag_vp(qJ, qJD, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(11,1),zeros(35,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPRRR8_coriolisvecJ_fixb_mdp_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRPRRR8_coriolisvecJ_fixb_mdp_slag_vp: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6RRPRRR8_coriolisvecJ_fixb_mdp_slag_vp: pkin has to be [11x1] (double)');
assert(isreal(MDP) && all(size(MDP) == [35 1]), ...
  'S6RRPRRR8_coriolisvecJ_fixb_mdp_slag_vp: MDP has to be [35x1] (double)'); 

%% Symbolic Calculation
% From coriolisvec_joint_fixb_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 14:05:49
% EndTime: 2019-03-09 14:06:09
% DurationCPUTime: 13.71s
% Computational Cost: add. (11232->583), mult. (28205->781), div. (0->0), fcn. (22244->10), ass. (0->256)
t669 = sin(qJ(2));
t673 = cos(qJ(2));
t635 = -pkin(2) * t673 - qJ(3) * t669 - pkin(1);
t613 = t635 * qJD(1);
t750 = qJD(1) * t673;
t658 = pkin(7) * t750;
t641 = qJD(2) * qJ(3) + t658;
t664 = sin(pkin(11));
t665 = cos(pkin(11));
t577 = t665 * t613 - t641 * t664;
t751 = qJD(1) * t669;
t728 = t665 * t751;
t749 = qJD(2) * t664;
t622 = t728 + t749;
t733 = pkin(3) * t750;
t531 = -pkin(8) * t622 + t577 - t733;
t578 = t664 * t613 + t665 * t641;
t729 = t664 * t751;
t738 = t665 * qJD(2);
t620 = t729 - t738;
t536 = -pkin(8) * t620 + t578;
t668 = sin(qJ(4));
t672 = cos(qJ(4));
t478 = t672 * t531 - t536 * t668;
t568 = t620 * t668 - t622 * t672;
t466 = pkin(9) * t568 + t478;
t653 = -qJD(4) + t750;
t461 = -pkin(4) * t653 + t466;
t479 = t531 * t668 + t536 * t672;
t569 = t672 * t620 + t622 * t668;
t467 = -pkin(9) * t569 + t479;
t671 = cos(qJ(5));
t465 = t671 * t467;
t667 = sin(qJ(5));
t428 = t461 * t667 + t465;
t506 = t568 * t667 - t671 * t569;
t807 = pkin(10) * t506;
t420 = t428 + t807;
t666 = sin(qJ(6));
t740 = qJD(6) * t666;
t418 = t420 * t740;
t670 = cos(qJ(6));
t795 = t568 * t671 + t569 * t667;
t781 = t795 * t666;
t457 = t506 * t670 + t781;
t657 = pkin(7) * t751;
t784 = qJD(2) * pkin(2);
t721 = qJD(3) - t784;
t634 = t657 + t721;
t584 = pkin(3) * t620 + t634;
t520 = pkin(4) * t569 + t584;
t469 = -pkin(5) * t506 + t520;
t821 = -t469 * t457 + t418;
t697 = -t506 * t666 + t670 * t795;
t748 = qJD(2) * t669;
t725 = MDP(33) * t748;
t820 = qJD(1) * t725 + (-t457 ^ 2 + t697 ^ 2) * MDP(30) + t457 * MDP(29) * t697;
t700 = pkin(2) * t669 - qJ(3) * t673;
t630 = t700 * qJD(1);
t611 = t664 * t630;
t777 = t665 * t669;
t778 = t664 * t673;
t689 = -pkin(7) * t777 - pkin(8) * t778;
t575 = qJD(1) * t689 + t611;
t819 = qJD(3) * t665 - t575;
t735 = qJD(1) * qJD(2);
t726 = t673 * t735;
t704 = t664 * t726;
t743 = qJD(4) * t672;
t772 = t672 * t665;
t523 = -t620 * t743 + t726 * t772 + (-qJD(4) * t622 - t704) * t668;
t628 = t664 * t672 + t665 * t668;
t688 = t673 * t628;
t681 = qJD(2) * t688;
t787 = qJD(4) * t568;
t524 = qJD(1) * t681 - t787;
t741 = qJD(5) * t671;
t742 = qJD(5) * t667;
t448 = t671 * t523 - t667 * t524 + t568 * t742 - t569 * t741;
t606 = qJD(2) * t700 - qJD(3) * t669;
t595 = t606 * qJD(1);
t633 = (qJD(3) - t657) * qJD(2);
t555 = t665 * t595 - t633 * t664;
t776 = t665 * t673;
t693 = pkin(3) * t669 - pkin(8) * t776;
t682 = t693 * qJD(2);
t527 = qJD(1) * t682 + t555;
t556 = t664 * t595 + t665 * t633;
t533 = -pkin(8) * t704 + t556;
t713 = t672 * t527 - t533 * t668;
t679 = -qJD(4) * t479 + t713;
t727 = t669 * t735;
t431 = pkin(4) * t727 - pkin(9) * t523 + t679;
t745 = qJD(4) * t668;
t686 = t668 * t527 + t531 * t743 + t672 * t533 - t536 * t745;
t433 = -pkin(9) * t524 + t686;
t719 = t671 * t431 - t667 * t433;
t678 = -qJD(5) * t428 + t719;
t409 = pkin(5) * t727 - pkin(10) * t448 + t678;
t676 = qJD(5) * t795 - t523 * t667 - t671 * t524;
t705 = -t667 * t431 - t671 * t433 - t461 * t741 + t467 * t742;
t410 = pkin(10) * t676 - t705;
t818 = -t666 * t409 - t670 * t410 + t821;
t739 = qJD(6) * t670;
t730 = t670 * t448 + t506 * t739 + t666 * t676;
t415 = t740 * t795 + t730;
t717 = t666 * t448 - t670 * t676;
t677 = qJD(6) * t697 - t717;
t646 = -qJD(5) + t653;
t640 = -qJD(6) + t646;
t801 = t640 * t697;
t810 = t646 * t795;
t812 = t506 * t646;
t816 = t457 * t640;
t817 = MDP(26) * t727 + (-t506 ^ 2 + t795 ^ 2) * MDP(23) + (t448 + t812) * MDP(24) + (t676 + t810) * MDP(25) + t506 * MDP(22) * t795 + (t677 + t801) * MDP(32) + (t415 + t816) * MDP(31) + t820;
t780 = t664 * t668;
t627 = -t772 + t780;
t687 = t627 * t673;
t758 = qJD(1) * t687 - t627 * qJD(4);
t757 = -qJD(1) * t688 + t628 * qJD(4);
t585 = pkin(7) * t729 + t665 * t630;
t557 = qJD(1) * t693 + t585;
t785 = pkin(8) + qJ(3);
t637 = t785 * t664;
t638 = t785 * t665;
t694 = qJD(3) * t664 + qJD(4) * t638;
t811 = -t637 * t743 + t819 * t672 + (-t557 - t694) * t668;
t720 = t670 * t409 - t666 * t410;
t798 = t469 * t697 + t720;
t549 = t672 * t557;
t756 = -t668 * t637 + t672 * t638;
t808 = pkin(4) * t751 + pkin(9) * t758 + t628 * qJD(3) + qJD(4) * t756 - t575 * t668 + t549;
t799 = -pkin(9) * t757 + t811;
t806 = pkin(10) * t795;
t803 = t568 * t653;
t802 = t569 * t653;
t573 = t671 * t627 + t628 * t667;
t763 = -qJD(5) * t573 - t667 * t757 + t671 * t758;
t574 = -t627 * t667 + t628 * t671;
t762 = qJD(5) * t574 + t667 * t758 + t671 * t757;
t797 = -t506 * t520 + t705;
t796 = t520 * t795 + t678;
t792 = -0.2e1 * t735;
t791 = MDP(4) * t669;
t790 = MDP(5) * (t669 ^ 2 - t673 ^ 2);
t789 = t808 * t671;
t603 = t627 * t669;
t619 = t665 * t635;
t576 = -pkin(8) * t777 + t619 + (-pkin(7) * t664 - pkin(3)) * t673;
t651 = pkin(7) * t776;
t591 = t664 * t635 + t651;
t779 = t664 * t669;
t583 = -pkin(8) * t779 + t591;
t710 = t672 * t576 - t583 * t668;
t490 = -pkin(4) * t673 + pkin(9) * t603 + t710;
t602 = t628 * t669;
t759 = t668 * t576 + t672 * t583;
t492 = -pkin(9) * t602 + t759;
t764 = t667 * t490 + t671 * t492;
t709 = -t672 * t637 - t638 * t668;
t550 = -pkin(9) * t628 + t709;
t551 = -pkin(9) * t627 + t756;
t761 = t667 * t550 + t671 * t551;
t614 = t664 * t733 + t658;
t723 = pkin(4) * t757 - t614;
t788 = t550 * t741 - t551 * t742 - t667 * t808 + t671 * t799;
t463 = t667 * t467;
t427 = t671 * t461 - t463;
t419 = t427 + t806;
t417 = -pkin(5) * t646 + t419;
t783 = t417 * t670;
t782 = t420 * t670;
t775 = t666 * t667;
t774 = t667 * t670;
t674 = qJD(2) ^ 2;
t773 = t669 * t674;
t771 = t673 * t674;
t675 = qJD(1) ^ 2;
t770 = t673 * t675;
t510 = t670 * t573 + t574 * t666;
t769 = -qJD(6) * t510 - t666 * t762 + t670 * t763;
t511 = -t573 * t666 + t574 * t670;
t768 = qJD(6) * t511 + t666 * t763 + t670 * t762;
t767 = t671 * t466 - t463;
t765 = pkin(5) * t762 + t723;
t731 = pkin(7) * t748;
t581 = t665 * t606 + t664 * t731;
t652 = pkin(7) * t726;
t605 = pkin(3) * t704 + t652;
t747 = qJD(2) * t673;
t659 = pkin(7) * t747;
t615 = t664 * pkin(3) * t747 + t659;
t631 = pkin(3) * t779 + t669 * pkin(7);
t744 = qJD(4) * t669;
t737 = MDP(11) * qJD(1);
t736 = MDP(12) * qJD(1);
t734 = pkin(7) * t778;
t732 = pkin(4) * qJD(5) * t640;
t655 = -pkin(3) * t665 - pkin(2);
t724 = MDP(19) * t748;
t722 = pkin(1) * t792;
t552 = -qJD(2) * t687 - t628 * t744;
t547 = t682 + t581;
t597 = t664 * t606;
t559 = qJD(2) * t689 + t597;
t712 = t672 * t547 - t559 * t668;
t443 = pkin(4) * t748 - pkin(9) * t552 - qJD(4) * t759 + t712;
t553 = t743 * t777 - t744 * t780 + t681;
t685 = t668 * t547 + t672 * t559 + t576 * t743 - t583 * t745;
t450 = -pkin(9) * t553 + t685;
t718 = t671 * t443 - t450 * t667;
t716 = -t466 * t667 - t465;
t715 = t671 * t490 - t492 * t667;
t711 = t671 * t550 - t551 * t667;
t708 = t620 + t738;
t707 = -t622 + t749;
t706 = qJD(6) * t417 + t410;
t495 = pkin(4) * t524 + t605;
t525 = pkin(4) * t553 + t615;
t580 = pkin(4) * t602 + t631;
t703 = -t634 + t721;
t472 = -pkin(10) * t574 + t711;
t702 = -t762 * pkin(10) + qJD(6) * t472 + t788;
t473 = -pkin(10) * t573 + t761;
t701 = pkin(5) * t751 + t763 * pkin(10) + qJD(5) * t761 + qJD(6) * t473 + t667 * t799 + t789;
t412 = t417 * t666 + t782;
t546 = -t602 * t667 - t603 * t671;
t437 = -pkin(5) * t673 - pkin(10) * t546 + t715;
t545 = t671 * t602 - t603 * t667;
t438 = -pkin(10) * t545 + t764;
t699 = t437 * t666 + t438 * t670;
t487 = t670 * t545 + t546 * t666;
t488 = -t545 * t666 + t546 * t670;
t596 = pkin(4) * t627 + t655;
t656 = pkin(4) * t671 + pkin(5);
t692 = pkin(4) * t774 + t656 * t666;
t691 = -pkin(4) * t775 + t656 * t670;
t684 = t667 * t443 + t671 * t450 + t490 * t741 - t492 * t742;
t590 = t619 - t734;
t586 = -pkin(7) * t728 + t611;
t582 = -t665 * t731 + t597;
t529 = pkin(5) * t573 + t596;
t500 = pkin(5) * t545 + t580;
t474 = -pkin(4) * t568 - pkin(5) * t795;
t471 = qJD(5) * t546 + t552 * t667 + t671 * t553;
t470 = -qJD(5) * t545 + t552 * t671 - t553 * t667;
t451 = pkin(5) * t471 + t525;
t434 = -pkin(5) * t676 + t495;
t424 = t767 + t806;
t423 = t716 - t807;
t422 = qJD(6) * t488 + t470 * t666 + t670 * t471;
t421 = -qJD(6) * t487 + t470 * t670 - t471 * t666;
t414 = -pkin(10) * t471 + t684;
t413 = pkin(5) * t748 - pkin(10) * t470 - qJD(5) * t764 + t718;
t411 = -t420 * t666 + t783;
t1 = [(-t523 * t603 - t552 * t568) * MDP(15) + (-t523 * t602 + t524 * t603 - t552 * t569 + t553 * t568) * MDP(16) + (t685 * t653 + t686 * t673 - t615 * t568 + t631 * t523 - t605 * t603 + t584 * t552 + (-qJD(1) * t759 - t479) * t748) * MDP(21) + (-t523 * t673 - t552 * t653 + (-qJD(1) * t603 - t568) * t748) * MDP(17) + (-t712 * t653 - t713 * t673 + t615 * t569 + t631 * t524 + t605 * t602 + t584 * t553 + (t479 * t673 + t653 * t759) * qJD(4) + (qJD(1) * t710 + t478) * t748) * MDP(20) + (t448 * t546 - t470 * t795) * MDP(22) + (t684 * t646 - t705 * t673 - t525 * t795 + t580 * t448 + t495 * t546 + t520 * t470 + (-qJD(1) * t764 - t428) * t748) * MDP(28) + (-t448 * t673 - t470 * t646 + (qJD(1) * t546 - t795) * t748) * MDP(24) + (t415 * t488 - t421 * t697) * MDP(29) + (-t415 * t673 - t421 * t640 + (qJD(1) * t488 - t697) * t748) * MDP(31) + (t500 * t415 - t418 * t673 + t469 * t421 + t434 * t488 - t451 * t697 + ((-qJD(6) * t438 + t413) * t640 + t409 * t673) * t666 + ((qJD(6) * t437 + t414) * t640 + t706 * t673) * t670 + (-qJD(1) * t699 - t412) * t748) * MDP(35) + (-t677 * t673 + t422 * t640 + (-qJD(1) * t487 + t457) * t748) * MDP(32) + (-t415 * t487 + t421 * t457 + t422 * t697 + t488 * t677) * MDP(30) + (-(t413 * t670 - t414 * t666) * t640 - t720 * t673 - t451 * t457 - t500 * t677 + t434 * t487 + t469 * t422 + (t412 * t673 + t640 * t699) * qJD(6) + ((t437 * t670 - t438 * t666) * qJD(1) + t411) * t748) * MDP(34) + MDP(6) * t771 + t790 * t792 - MDP(7) * t773 + (pkin(7) * t773 + t673 * t722) * MDP(10) + (-pkin(7) * t771 + t669 * t722) * MDP(9) + (t524 * t673 + t553 * t653 + (-qJD(1) * t602 - t569) * t748) * MDP(18) + (-t653 - t750) * t724 + (-t640 - t750) * t725 + (-t581 * t622 - t582 * t620 + (-t555 * t665 - t556 * t664) * t669 + (-t577 * t665 - t578 * t664 + (-t590 * t665 - t591 * t664) * qJD(1)) * t747) * MDP(13) + ((-qJD(1) * t581 - t555) * t673 + ((pkin(7) * t620 + t634 * t664) * t673 + (t577 + (t590 + 0.2e1 * t734) * qJD(1)) * t669) * qJD(2)) * MDP(11) + (t555 * t590 + t556 * t591 + t577 * t581 + t578 * t582 + (t634 + t657) * t659) * MDP(14) + ((qJD(1) * t582 + t556) * t673 + ((pkin(7) * t622 + t634 * t665) * t673 + (-t578 + (-t591 + 0.2e1 * t651) * qJD(1)) * t669) * qJD(2)) * MDP(12) + (-t646 - t750) * MDP(26) * t748 + 0.2e1 * t726 * t791 + (-t718 * t646 - t719 * t673 - t525 * t506 - t580 * t676 + t495 * t545 + t520 * t471 + (t428 * t673 + t646 * t764) * qJD(5) + (qJD(1) * t715 + t427) * t748) * MDP(27) + (-t676 * t673 + t471 * t646 + (-qJD(1) * t545 + t506) * t748) * MDP(25) + (-t448 * t545 + t470 * t506 + t471 * t795 + t546 * t676) * MDP(23); (-t577 * t585 - t578 * t586 + (-t577 * t664 + t578 * t665) * qJD(3) + (-t555 * t664 + t556 * t665) * qJ(3) + (-t634 - t784) * t658) * MDP(14) + (t448 * t574 - t763 * t795) * MDP(22) + (-t448 * t573 + t506 * t763 + t574 * t676 + t762 * t795) * MDP(23) + (t415 * t511 - t697 * t769) * MDP(29) + (-t415 * t510 + t457 * t769 + t511 * t677 + t697 * t768) * MDP(30) + (t523 * t628 - t568 * t758) * MDP(15) + (-t523 * t627 - t524 * t628 + t568 * t757 - t569 * t758) * MDP(16) + ((-qJ(3) * t749 - t577) * t669 + (-pkin(7) * t708 + t664 * t703 + t585) * t673) * t737 + (t585 * t622 + t586 * t620 + (-qJD(3) * t620 + t577 * t750 + t556) * t665 + (qJD(3) * t622 + t578 * t750 - t555) * t664) * MDP(13) + ((-qJ(3) * t738 + t578) * t669 + (pkin(7) * t707 + t665 * t703 - t586) * t673) * t736 + t675 * t790 - t770 * t791 + (t596 * t448 + t495 * t574 + t763 * t520 - t723 * t795) * MDP(28) + (t434 * t510 - t457 * t765 + t768 * t469 - t529 * t677) * MDP(34) + (t529 * t415 + t434 * t511 + t769 * t469 - t697 * t765) * MDP(35) + (t655 * t524 - t614 * t569 + t757 * t584 + t605 * t627) * MDP(20) + (t655 * t523 + t568 * t614 + t758 * t584 + t605 * t628) * MDP(21) + (t495 * t573 - t506 * t723 + t762 * t520 - t596 * t676) * MDP(27) + (t757 * MDP(18) - t758 * MDP(17) + (t549 + t694 * t672 + (-qJD(4) * t637 + t819) * t668) * MDP(20) + t811 * MDP(21)) * t653 + (-t763 * MDP(24) + t762 * MDP(25) + t788 * MDP(28) + (t551 * t741 + (qJD(5) * t550 + t799) * t667 + t789) * MDP(27)) * t646 + (-t769 * MDP(31) + t768 * MDP(32) + (t666 * t702 + t670 * t701) * MDP(34) + (-t666 * t701 + t670 * t702) * MDP(35)) * t640 + (MDP(9) * t669 * t675 + MDP(10) * t770) * pkin(1) + ((qJD(2) * t511 + t697) * MDP(31) + (qJD(2) * t574 + t795) * MDP(24) + (-qJD(2) * t761 + t428) * MDP(28) + ((t472 * t670 - t473 * t666) * qJD(2) - t411) * MDP(34) + (-qJD(2) * t510 - t457) * MDP(32) + (-(t472 * t666 + t473 * t670) * qJD(2) + t412) * MDP(35) + (-qJD(2) * t627 + t569) * MDP(18) + (qJD(2) * t709 - t478) * MDP(20) + (qJD(2) * t628 + t568) * MDP(17) + (-qJD(2) * t756 + t479) * MDP(21) + (-qJD(2) * t573 - t506) * MDP(25) + (qJD(2) * t711 - t427) * MDP(27) + t653 * MDP(19) + t646 * MDP(26) + t640 * MDP(33)) * t751; (-t620 ^ 2 - t622 ^ 2) * MDP(13) + (t577 * t622 + t578 * t620 + t652) * MDP(14) + (t524 + t803) * MDP(20) + (t523 + t802) * MDP(21) + (-t676 + t810) * MDP(27) + (t448 - t812) * MDP(28) + (-t677 + t801) * MDP(34) + (t415 - t816) * MDP(35) + (t707 * t737 + t708 * t736) * t673; (t523 - t802) * MDP(17) + (t568 ^ 2 - t569 ^ 2) * MDP(16) + (-t628 * t726 + t787 + t803) * MDP(18) + (-t479 * t653 + t568 * t584 + t679) * MDP(20) + (t716 * t646 + (-t506 * t568 + t646 * t742 + t671 * t727) * pkin(4) + t796) * MDP(27) + (-t767 * t646 + (-t568 * t795 + t646 * t741 - t667 * t727) * pkin(4) + t797) * MDP(28) + qJD(1) * t724 + (-t478 * t653 + t569 * t584 - t686) * MDP(21) - t568 * t569 * MDP(15) + (-t692 * t727 - (t423 * t666 + t424 * t670) * t640 + t474 * t697 + (t670 * t671 - t775) * t732 + (t640 * t691 - t783) * qJD(6) + t818) * MDP(35) + (t691 * t727 + (t423 * t670 - t424 * t666) * t640 + t474 * t457 - (-t666 * t671 - t774) * t732 + (t640 * t692 - t412) * qJD(6) + t798) * MDP(34) + t817; (-t428 * t646 + t796) * MDP(27) + (-t427 * t646 + t797) * MDP(28) + ((-t419 * t666 - t782) * t640 - t412 * qJD(6) + (-t457 * t795 + t640 * t740 + t670 * t727) * pkin(5) + t798) * MDP(34) + ((t420 * t640 - t409) * t666 + (-t419 * t640 - t706) * t670 + (t640 * t739 - t666 * t727 - t697 * t795) * pkin(5) + t821) * MDP(35) + t817; (t730 + t816) * MDP(31) + (-t717 + t801) * MDP(32) + (-t412 * t640 + t798) * MDP(34) + (-t411 * t640 + t818) * MDP(35) + (MDP(31) * t781 + MDP(32) * t697 - MDP(34) * t412 - MDP(35) * t783) * qJD(6) + t820;];
tauc  = t1;
