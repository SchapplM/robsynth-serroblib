% Calculate minimal parameter regressor of Coriolis joint torque vector for
% S6RRPRRR13
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [11x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d1,d2,d4,d5,d6]';
% MDP [35x1]
%   Minimal dynamic parameter vector (fixed base model)
%   see S6RRPRRR13_convert_par2_MPV_fixb.m
% 
% Output:
% tauc [6x1]
%   joint torques required to compensate Coriolis and centrifugal load

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 14:56
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauc = S6RRPRRR13_coriolisvecJ_fixb_mdp_slag_vp(qJ, qJD, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(11,1),zeros(35,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPRRR13_coriolisvecJ_fixb_mdp_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRPRRR13_coriolisvecJ_fixb_mdp_slag_vp: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6RRPRRR13_coriolisvecJ_fixb_mdp_slag_vp: pkin has to be [11x1] (double)');
assert(isreal(MDP) && all(size(MDP) == [35 1]), ...
  'S6RRPRRR13_coriolisvecJ_fixb_mdp_slag_vp: MDP has to be [35x1] (double)'); 

%% Symbolic Calculation
% From coriolisvec_joint_fixb_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 14:53:38
% EndTime: 2019-03-09 14:53:59
% DurationCPUTime: 12.35s
% Computational Cost: add. (8239->667), mult. (20746->899), div. (0->0), fcn. (15949->10), ass. (0->276)
t659 = cos(qJ(2));
t650 = sin(pkin(6));
t655 = sin(qJ(2));
t752 = qJD(1) * t655;
t727 = t650 * t752;
t651 = cos(pkin(6));
t753 = qJD(1) * t651;
t736 = pkin(1) * t753;
t579 = pkin(8) * t727 - t659 * t736;
t738 = qJD(3) + t579;
t654 = sin(qJ(4));
t658 = cos(qJ(4));
t693 = pkin(4) * t658 + pkin(10) * t654;
t813 = (-pkin(3) - t693) * t727 - qJD(4) * t693 - t738;
t653 = sin(qJ(5));
t657 = cos(qJ(5));
t754 = qJD(1) * t650;
t726 = t659 * t754;
t567 = t653 * t654 * t727 - t657 * t726;
t747 = qJD(4) * t654;
t806 = t653 * t747 + t567;
t770 = t655 * t657;
t568 = (t653 * t659 + t654 * t770) * t754;
t691 = -t657 * t747 - t568;
t743 = qJD(5) * t658;
t812 = -t653 * t743 + t691;
t637 = qJD(2) + t753;
t660 = -pkin(2) - pkin(9);
t739 = pkin(3) * t727 + t738;
t501 = t637 * t660 + t739;
t714 = -qJ(3) * t655 - pkin(1);
t555 = (t659 * t660 + t714) * t650;
t530 = qJD(1) * t555;
t459 = t501 * t658 - t654 * t530;
t620 = qJD(4) + t727;
t446 = -pkin(4) * t620 - t459;
t697 = t654 * t726;
t573 = t637 * t658 - t697;
t511 = t573 * t653 - t657 * t620;
t440 = pkin(5) * t511 + t446;
t656 = cos(qJ(6));
t513 = t573 * t657 + t620 * t653;
t652 = sin(qJ(6));
t784 = t513 * t652;
t463 = t656 * t511 + t784;
t811 = t440 * t463;
t683 = t511 * t652 - t656 * t513;
t810 = t440 * t683;
t571 = t637 * t654 + t658 * t726;
t570 = qJD(5) + t571;
t562 = qJD(6) + t570;
t809 = t463 * t562;
t808 = t562 * t683;
t807 = t657 * t743 - t806;
t797 = qJD(5) + qJD(6);
t805 = t571 + t797;
t737 = qJD(1) * qJD(2);
t715 = t650 * t737;
t696 = t655 * t715;
t746 = qJD(4) * t658;
t519 = -qJD(4) * t697 + t637 * t746 - t658 * t696;
t804 = t519 * MDP(33) + (-t463 ^ 2 + t683 ^ 2) * MDP(30) - t463 * MDP(29) * t683;
t648 = t655 ^ 2;
t803 = MDP(5) * (-t659 ^ 2 + t648);
t802 = t511 * t620;
t596 = t652 * t657 + t653 * t656;
t582 = t596 * t658;
t774 = t650 * t655;
t639 = pkin(8) * t774;
t729 = -pkin(1) * t659 - pkin(2);
t535 = pkin(3) * t774 + t639 + (-pkin(9) + t729) * t651;
t763 = t654 * t535 + t658 * t555;
t474 = pkin(10) * t774 + t763;
t794 = pkin(1) * t655;
t642 = t651 * t794;
t773 = t650 * t659;
t800 = pkin(8) * t773 + t642;
t575 = -t651 * qJ(3) - t800;
t554 = pkin(3) * t773 - t575;
t586 = t651 * t654 + t658 * t773;
t732 = t654 * t773;
t587 = t651 * t658 - t732;
t484 = pkin(4) * t586 - pkin(10) * t587 + t554;
t767 = t657 * t474 + t653 * t484;
t801 = t813 * t657;
t635 = pkin(2) * t727;
t686 = pkin(9) * t655 - qJ(3) * t659;
t557 = t686 * t754 + t635;
t580 = pkin(8) * t726 + t655 * t736;
t559 = pkin(3) * t726 + t580;
t760 = t658 * t557 + t654 * t559;
t477 = pkin(10) * t726 + t760;
t605 = pkin(4) * t654 - pkin(10) * t658 + qJ(3);
t720 = t657 * t746;
t744 = qJD(5) * t657;
t799 = t657 * t477 - t605 * t744 + t813 * t653 - t660 * t720;
t798 = -qJD(6) * t657 - t744;
t518 = -qJD(4) * t571 + t654 * t696;
t625 = t659 * t715;
t456 = qJD(5) * t513 + t518 * t653 - t657 * t625;
t745 = qJD(5) * t653;
t455 = t657 * t518 - t573 * t745 + t620 * t744 + t653 * t625;
t710 = t455 * t652 + t656 * t456;
t417 = -qJD(6) * t683 + t710;
t796 = pkin(3) + pkin(8);
t795 = pkin(10) + pkin(11);
t793 = pkin(11) * t657;
t792 = pkin(11) * t658;
t460 = t654 * t501 + t658 * t530;
t447 = pkin(10) * t620 + t460;
t626 = t637 * qJ(3);
t524 = t626 + t559;
t467 = pkin(4) * t571 - pkin(10) * t573 + t524;
t428 = -t447 * t653 + t657 * t467;
t421 = -pkin(11) * t513 + t428;
t419 = pkin(5) * t570 + t421;
t791 = t419 * t656;
t429 = t447 * t657 + t467 * t653;
t422 = -pkin(11) * t511 + t429;
t790 = t422 * t656;
t611 = pkin(2) * t696;
t749 = qJD(3) * t655;
t662 = (qJD(2) * t686 - t749) * t650;
t507 = qJD(1) * t662 + t611;
t735 = pkin(1) * qJD(2) * t651;
t701 = qJD(1) * t735;
t574 = pkin(8) * t625 + t655 * t701;
t536 = pkin(3) * t625 + t574;
t699 = t501 * t747 + t654 * t507 + t530 * t746 - t658 * t536;
t432 = -pkin(4) * t625 + t699;
t789 = t432 * t653;
t788 = t432 * t657;
t787 = t455 * t653;
t786 = t511 * t570;
t785 = t513 * t570;
t595 = t652 * t653 - t656 * t657;
t783 = t519 * t595;
t782 = t519 * t596;
t781 = t519 * t653;
t780 = t519 * t657;
t779 = t571 * t620;
t778 = t571 * t653;
t777 = t573 * t620;
t776 = t620 * t654;
t671 = t620 * t658;
t647 = t650 ^ 2;
t775 = t647 * qJD(1) ^ 2;
t772 = t653 * t660;
t771 = t654 * t660;
t769 = t657 * t658;
t500 = pkin(4) * t573 + pkin(10) * t571;
t768 = t657 * t459 + t653 * t500;
t765 = t567 * t652 - t568 * t656 - t582 * t797 + t595 * t747;
t742 = qJD(6) * t652;
t717 = t653 * t742;
t764 = -t658 * t717 + (t769 * t797 - t806) * t656 + t812 * t652;
t762 = t805 * t595;
t761 = t805 * t596;
t706 = -t654 * t557 + t559 * t658;
t476 = -pkin(4) * t726 - t706;
t759 = pkin(5) * t807 + t660 * t747 - t476;
t758 = pkin(8) * t696 - t659 * t701;
t636 = t657 * t771;
t757 = t653 * t605 + t636;
t751 = qJD(2) * t655;
t750 = qJD(2) * t659;
t748 = qJD(4) * t562;
t741 = qJD(6) * t656;
t734 = t655 * t671;
t733 = t659 * t775;
t731 = t656 * t455 - t652 * t456 - t511 * t741;
t544 = -t637 * qJD(3) + t758;
t728 = qJD(5) * t795;
t725 = t650 * t751;
t724 = t650 * t750;
t723 = t620 * t746;
t718 = t660 * t745;
t716 = t647 * t737;
t713 = pkin(5) - t772;
t668 = t501 * t746 + t658 * t507 - t530 * t747 + t654 * t536;
t431 = pkin(10) * t625 + t668;
t508 = -pkin(3) * t696 - t544;
t443 = pkin(4) * t519 - pkin(10) * t518 + t508;
t412 = -qJD(5) * t429 - t431 * t653 + t657 * t443;
t408 = pkin(5) * t519 - pkin(11) * t455 + t412;
t411 = t657 * t431 + t653 * t443 - t447 * t745 + t467 * t744;
t409 = -pkin(11) * t456 + t411;
t712 = t656 * t408 - t652 * t409;
t420 = t422 * t742;
t711 = t652 * t408 - t420;
t709 = -t474 * t653 + t657 * t484;
t707 = t535 * t658 - t654 * t555;
t705 = t570 * t657;
t703 = qJD(6) * t419 + t409;
t700 = t655 * t733;
t698 = t658 * t727;
t695 = -t460 + (t745 + t778) * pkin(5);
t694 = -0.2e1 * pkin(1) * t716;
t692 = t580 * t637 - t574;
t495 = t657 * t500;
t624 = t795 * t657;
t690 = pkin(5) * t573 + qJD(6) * t624 - t459 * t653 + t571 * t793 + t657 * t728 + t495;
t623 = t795 * t653;
t689 = pkin(11) * t778 + qJD(6) * t623 + t653 * t728 + t768;
t547 = -t653 * t792 + t757;
t688 = -pkin(5) * t698 - pkin(11) * t568 + qJD(6) * t547 - t477 * t653 - (-t636 + (-t605 + t792) * t653) * qJD(5) - (t654 * t793 + t658 * t713) * qJD(4) + t801;
t591 = t657 * t605;
t529 = -pkin(11) * t769 + t654 * t713 + t591;
t687 = pkin(11) * t807 - qJD(6) * t529 + t654 * t718 + t799;
t414 = t419 * t652 + t790;
t551 = t587 * t657 + t653 * t774;
t426 = pkin(5) * t586 - pkin(11) * t551 + t709;
t678 = -t587 * t653 + t650 * t770;
t433 = pkin(11) * t678 + t767;
t685 = t426 * t656 - t433 * t652;
t684 = t426 * t652 + t433 * t656;
t482 = t551 * t652 - t656 * t678;
t483 = t551 * t656 + t652 * t678;
t581 = t800 * qJD(2);
t682 = t574 * t651 + t581 * t637;
t634 = t659 * t735;
t681 = -pkin(8) * t725 + t634;
t630 = pkin(2) * t725;
t526 = t630 + t662;
t560 = (t773 * t796 + t642) * qJD(2);
t680 = -t654 * t526 - t535 * t747 - t555 * t746 + t560 * t658;
t679 = -t637 * t726 + t625;
t677 = t524 * t655 + t660 * t750;
t676 = -t570 * t744 - t781;
t675 = -t570 * t745 + t780;
t672 = t620 * t513;
t669 = -pkin(10) * t519 + t446 * t570;
t667 = t658 * t526 + t535 * t746 - t555 * t747 + t654 * t560;
t438 = pkin(10) * t724 + t667;
t645 = t651 * qJD(3);
t534 = -t725 * t796 + t634 + t645;
t548 = -qJD(4) * t586 + t654 * t725;
t549 = -qJD(4) * t732 + t651 * t746 - t658 * t725;
t454 = pkin(4) * t549 - pkin(10) * t548 + t534;
t666 = t657 * t438 + t653 * t454 - t474 * t745 + t484 * t744;
t416 = -t513 * t742 + t731;
t473 = -pkin(4) * t774 - t707;
t576 = (-pkin(2) * t659 + t714) * t650;
t664 = (-qJ(3) * t750 - t749) * t650;
t439 = -pkin(4) * t724 - t680;
t663 = -qJD(5) * t767 - t438 * t653 + t657 * t454;
t406 = -qJD(6) * t414 + t712;
t644 = -pkin(5) * t657 - pkin(4);
t603 = t658 * t625;
t593 = (pkin(5) * t653 - t660) * t658;
t583 = t595 * t658;
t578 = -qJ(3) * t726 + t635;
t577 = t651 * t729 + t639;
t569 = -t645 - t681;
t566 = t637 * t653 - t657 * t698;
t565 = t637 * t657 + t653 * t698;
t564 = qJD(1) * t576;
t561 = t630 + t664;
t556 = -t626 - t580;
t552 = -pkin(2) * t637 + t738;
t541 = qJD(1) * t664 + t611;
t532 = t564 * t727;
t510 = t519 * t654;
t491 = t519 * t586;
t471 = qJD(5) * t678 + t548 * t657 + t653 * t724;
t470 = qJD(5) * t551 + t548 * t653 - t657 * t724;
t445 = -pkin(5) * t678 + t473;
t425 = pkin(5) * t470 + t439;
t424 = qJD(6) * t483 + t656 * t470 + t471 * t652;
t423 = -qJD(6) * t482 - t470 * t652 + t471 * t656;
t418 = pkin(5) * t456 + t432;
t415 = -pkin(11) * t470 + t666;
t413 = -t422 * t652 + t791;
t410 = pkin(5) * t549 - pkin(11) * t471 + t663;
t405 = t703 * t656 + t711;
t1 = [0.2e1 * (t655 * t659 * MDP(4) - t803) * t716 + (-(qJD(6) * t685 + t410 * t652 + t415 * t656) * t562 - t684 * t519 - t405 * t586 - t414 * t549 - t425 * t683 + t445 * t416 + t418 * t483 + t440 * t423) * MDP(35) + (t416 * t586 + t423 * t562 + t483 * t519 - t549 * t683) * MDP(31) + (-t416 * t482 - t417 * t483 - t423 * t463 + t424 * t683) * MDP(30) + (t416 * t483 - t423 * t683) * MDP(29) + (-t456 * t586 - t470 * t570 - t511 * t549 + t519 * t678) * MDP(25) + (t412 * t586 + t428 * t549 - t432 * t678 + t439 * t511 + t446 * t470 + t473 * t456 + t519 * t709 + t570 * t663) * MDP(27) + (t455 * t678 - t456 * t551 - t470 * t513 - t471 * t511) * MDP(23) + (MDP(6) * t724 - MDP(7) * t725) * (t637 + t753) + (-t544 * t651 - t569 * t637 + (-t564 * t750 - t541 * t655 + (-t561 * t655 - t576 * t750) * qJD(1)) * t650) * MDP(13) + ((-t564 * t751 + t541 * t659 + (t561 * t659 - t576 * t751) * qJD(1)) * t650 + t682) * MDP(12) + ((-qJD(6) * t684 + t410 * t656 - t415 * t652) * t562 + t685 * t519 + t406 * t586 + t413 * t549 + t425 * t463 + t445 * t417 + t418 * t482 + t440 * t424) * MDP(34) + (t548 * t620 + (t518 * t655 + (qJD(1) * t587 + t573) * t750) * t650) * MDP(17) + (-t549 * t620 + (-t519 * t655 + (-qJD(1) * t586 - t571) * t750) * t650) * MDP(18) + (t680 * t620 + t534 * t571 + t554 * t519 + t508 * t586 + t524 * t549 + (-t699 * t655 + (qJD(1) * t707 + t459) * t750) * t650) * MDP(20) + (t620 * t650 + t647 * t752) * MDP(19) * t750 + (-t667 * t620 + t534 * t573 + t554 * t518 + t508 * t587 + t524 * t548 + (-t668 * t655 + (-qJD(1) * t763 - t460) * t750) * t650) * MDP(21) + (-t411 * t586 - t429 * t549 + t432 * t551 + t439 * t513 + t446 * t471 + t473 * t455 - t519 * t767 - t570 * t666) * MDP(28) + (t455 * t586 + t471 * t570 + t513 * t549 + t519 * t551) * MDP(24) + (-t417 * t586 - t424 * t562 - t463 * t549 - t482 * t519) * MDP(32) + (-t518 * t586 - t519 * t587 - t548 * t571 - t549 * t573) * MDP(16) + (t518 * t587 + t548 * t573) * MDP(15) + (t541 * t576 + t544 * t575 + t552 * t581 + t556 * t569 + t561 * t564 + t574 * t577) * MDP(14) + (t549 * t570 + t491) * MDP(26) + (t549 * t562 + t491) * MDP(33) + (t455 * t551 + t471 * t513) * MDP(22) + (-t637 * t681 + t651 * t758 + t659 * t694) * MDP(10) + (t655 * t694 - t682) * MDP(9) + (-t544 * t659 + t574 * t655 + (t552 * t659 + t556 * t655) * qJD(2) + (-t569 * t659 + t581 * t655 + (t575 * t655 + t577 * t659) * qJD(2)) * qJD(1)) * t650 * MDP(11); (-t757 * t519 - t476 * t513 - t446 * t568 + t799 * t570 + (t570 * t718 - t411 + (-t446 * t657 + t513 * t660) * qJD(4)) * t654 + (-t429 * t620 - t446 * t745 - t660 * t455 + t788) * t658) * MDP(28) + (t455 * t769 + t812 * t513) * MDP(22) + (-pkin(2) * t574 - qJ(3) * t544 - t552 * t580 - t556 * t738 - t564 * t578) * MDP(14) + (-t456 * t654 + t806 * t570 + (t676 - t802) * t658) * MDP(25) + (-t578 * t726 + t532 - t692) * MDP(12) + (t416 * t654 - t519 * t583 + t562 * t765 - t671 * t683) * MDP(31) + (-t416 * t583 - t683 * t765) * MDP(29) + (-t416 * t582 + t417 * t583 - t463 * t765 + t683 * t764) * MDP(30) + (t518 * t658 - t573 * t776) * MDP(15) + (-(t529 * t652 + t547 * t656) * t519 - t405 * t654 + t593 * t416 - t418 * t583 + (t652 * t688 + t656 * t687) * t562 - t759 * t683 + t765 * t440 - t414 * t671) * MDP(35) + (-t417 * t654 - t463 * t671 - t519 * t582 - t562 * t764) * MDP(32) + (qJ(3) * t519 + t508 * t654 - t706 * t620 + t739 * t571 + (t524 * t658 - t620 * t771) * qJD(4) + (-t459 * t659 + t658 * t677) * t754) * MDP(20) + (-t446 * t567 - t476 * t511 + t591 * t519 + ((-qJD(5) * t605 + t477) * t653 - t801) * t570 + (-t446 * t653 * qJD(4) + t412 + (qJD(4) * t511 + t676) * t660) * t654 + (t428 * t727 + t446 * t744 + t789 - t660 * t456 + (-t570 * t772 + t428) * qJD(4)) * t658) * MDP(27) + (-qJD(2) + t637) * MDP(7) * t727 + (t570 * t671 + t510) * MDP(26) + (t562 * t671 + t510) * MDP(33) + (-t723 + (-t734 + (-qJD(2) * t654 + t571) * t659) * t754) * MDP(18) + (t738 * t637 + (t564 * t659 + t578 * t655) * t754 - t544) * MDP(13) + t679 * MDP(6) + (t455 * t654 + t691 * t570 + (t672 + t675) * t658) * MDP(24) - MDP(4) * t700 - t620 * MDP(19) * t726 + ((t529 * t656 - t547 * t652) * t519 + t406 * t654 + t593 * t417 + t418 * t582 + (t652 * t687 - t656 * t688) * t562 + t759 * t463 + t764 * t440 + t413 * t671) * MDP(34) + ((-qJ(3) * qJD(2) - t556 - t580) * t655 + (-pkin(2) * qJD(2) - t552 + t738) * t659) * MDP(11) * t754 + (pkin(1) * t733 - t579 * t637 + t758) * MDP(10) + t775 * t803 + (qJ(3) * t518 + t508 * t658 + t760 * t620 + t739 * t573 + (-t524 * t654 - t660 * t671) * qJD(4) + (t460 * t659 - t654 * t677) * t754) * MDP(21) + (-t620 * t747 + t603 + (-t573 * t659 - t655 * t776) * t754) * MDP(17) + ((-t519 - t777) * t658 + (-t518 + t779) * t654) * MDP(16) + (t511 * t568 + t513 * t567 + (t511 * t657 + t513 * t653) * t747 + (-t787 - t456 * t657 + (t511 * t653 - t513 * t657) * qJD(5)) * t658) * MDP(23) + (t775 * t794 + t692) * MDP(9); t679 * MDP(11) + MDP(12) * t700 + (-t637 ^ 2 - t648 * t775) * MDP(13) + (t556 * t637 + t532 + t574) * MDP(14) + (-t571 * t637 - t620 * t776 + t603) * MDP(20) + (-t723 - t573 * t637 + (-t654 * t750 - t734) * t754) * MDP(21) + (-t456 * t658 + (-t653 * t746 - t565) * t570 + (t676 + t802) * t654) * MDP(27) + (-t455 * t658 + (t566 - t720) * t570 + (t672 - t675) * t654) * MDP(28) + (-(t565 * t656 - t566 * t652) * t562 + (-t596 * t748 - t417) * t658 + ((t652 * t745 + t656 * t798 + t717) * t562 - t782 + t620 * t463) * t654) * MDP(34) + ((t565 * t652 + t566 * t656) * t562 + (t595 * t748 - t416) * t658 + (-(t652 * t798 - t653 * t741 - t656 * t745) * t562 + t783 - t620 * t683) * t654) * MDP(35); -t571 ^ 2 * MDP(16) + (t518 + t779) * MDP(17) + (-t519 + t777) * MDP(18) + MDP(19) * t625 + (t460 * t620 - t699) * MDP(20) + (t459 * t620 + t524 * t571 - t668) * MDP(21) + (t513 * t705 + t787) * MDP(22) + ((t455 - t786) * t657 + (-t456 - t785) * t653) * MDP(23) + (t570 * t705 + t781) * MDP(24) + (-t570 ^ 2 * t653 + t780) * MDP(25) + (-pkin(4) * t456 - t788 - t460 * t511 + (-pkin(10) * t744 - t495) * t570 + (t459 * t570 + t669) * t653) * MDP(27) + (-pkin(4) * t455 + t789 - t460 * t513 + (pkin(10) * t745 + t768) * t570 + t669 * t657) * MDP(28) + (t416 * t596 + t683 * t762) * MDP(29) + (-t416 * t595 - t417 * t596 + t463 * t762 + t683 * t761) * MDP(30) + (-t562 * t762 + t782) * MDP(31) + (-t562 * t761 - t783) * MDP(32) + ((-t623 * t656 - t624 * t652) * t519 + t644 * t417 + t418 * t595 + (t652 * t689 - t656 * t690) * t562 + t695 * t463 + t761 * t440) * MDP(34) + (-(-t623 * t652 + t624 * t656) * t519 + t644 * t416 + t418 * t596 + (t652 * t690 + t656 * t689) * t562 - t695 * t683 - t762 * t440) * MDP(35) + (MDP(15) * t571 + t573 * MDP(16) - t524 * MDP(20) - t513 * MDP(24) + t511 * MDP(25) - t570 * MDP(26) - t428 * MDP(27) + t429 * MDP(28) + MDP(31) * t683 + t463 * MDP(32) - t562 * MDP(33) - t413 * MDP(34) + t414 * MDP(35)) * t573; t513 * t511 * MDP(22) + (-t511 ^ 2 + t513 ^ 2) * MDP(23) + (t455 + t786) * MDP(24) + (-t456 + t785) * MDP(25) + t519 * MDP(26) + (t429 * t570 - t446 * t513 + t412) * MDP(27) + (t428 * t570 + t446 * t511 - t411) * MDP(28) + (t416 + t809) * MDP(31) + (-t417 - t808) * MDP(32) + (-(-t421 * t652 - t790) * t562 + t810 + (-t463 * t513 + t519 * t656 - t562 * t742) * pkin(5) + t406) * MDP(34) + (t811 + t420 + (-t422 * t562 - t408) * t652 + (t421 * t562 - t703) * t656 + (t513 * t683 - t519 * t652 - t562 * t741) * pkin(5)) * MDP(35) + t804; (t731 + t809) * MDP(31) + (-t710 - t808) * MDP(32) + (t414 * t562 + t712 + t810) * MDP(34) + (-t656 * t409 + t413 * t562 - t711 + t811) * MDP(35) + (-MDP(31) * t784 + MDP(32) * t683 - MDP(34) * t414 - MDP(35) * t791) * qJD(6) + t804;];
tauc  = t1;
