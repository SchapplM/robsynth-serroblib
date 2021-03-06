% Calculate minimal parameter regressor of Coriolis joint torque vector for
% S6RRRRPR10
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [11x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d1,d2,d3,d4,d6]';
% MDP [35x1]
%   Minimal dynamic parameter vector (fixed base model)
%   see S6RRRRPR10_convert_par2_MPV_fixb.m
% 
% Output:
% tauc [6x1]
%   joint torques required to compensate Coriolis and centrifugal load

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 23:11
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauc = S6RRRRPR10_coriolisvecJ_fixb_mdp_slag_vp(qJ, qJD, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(11,1),zeros(35,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRRPR10_coriolisvecJ_fixb_mdp_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRRRPR10_coriolisvecJ_fixb_mdp_slag_vp: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6RRRRPR10_coriolisvecJ_fixb_mdp_slag_vp: pkin has to be [11x1] (double)');
assert(isreal(MDP) && all(size(MDP) == [35 1]), ...
  'S6RRRRPR10_coriolisvecJ_fixb_mdp_slag_vp: MDP has to be [35x1] (double)'); 

%% Symbolic Calculation
% From coriolisvec_joint_fixb_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 23:08:53
% EndTime: 2019-03-09 23:09:15
% DurationCPUTime: 13.75s
% Computational Cost: add. (11748->652), mult. (30360->843), div. (0->0), fcn. (23959->10), ass. (0->263)
t650 = cos(qJ(2));
t644 = sin(pkin(6));
t795 = sin(qJ(3));
t742 = t644 * t795;
t708 = qJD(1) * t742;
t686 = t650 * t708;
t731 = qJD(3) * t795;
t822 = t686 - t731;
t797 = cos(qJ(3));
t748 = t797 * pkin(9);
t624 = pkin(10) * t797 + t748;
t648 = sin(qJ(2));
t762 = qJD(1) * t644;
t738 = t648 * t762;
t645 = cos(pkin(6));
t761 = qJD(1) * t645;
t747 = pkin(1) * t761;
t594 = -pkin(8) * t738 + t650 * t747;
t678 = (pkin(2) * t648 - pkin(9) * t650) * t644;
t595 = qJD(1) * t678;
t696 = -t594 * t795 + t797 * t595;
t741 = t650 * t797;
t824 = -(pkin(3) * t648 - pkin(10) * t741) * t762 - t696 - t624 * qJD(3);
t623 = (-pkin(10) - pkin(9)) * t795;
t767 = t797 * t594 + t795 * t595;
t823 = -pkin(10) * t686 - t623 * qJD(3) + t767;
t647 = sin(qJ(4));
t733 = qJD(1) * t797;
t709 = t644 * t733;
t683 = t650 * t709;
t796 = cos(qJ(4));
t710 = t796 * t797;
t802 = qJD(3) + qJD(4);
t770 = t683 * t796 - t802 * t710 + (qJD(4) * t795 - t822) * t647;
t760 = qJD(1) * t650;
t737 = t644 * t760;
t622 = -qJD(3) + t737;
t614 = -qJD(4) + t622;
t717 = qJD(2) + t761;
t578 = -t648 * t708 + t717 * t797;
t579 = t648 * t709 + t795 * t717;
t729 = t797 * qJD(2);
t707 = t644 * t729;
t682 = qJD(1) * t707;
t653 = qJD(3) * t578 + t650 * t682;
t730 = qJD(4) * t796;
t728 = t795 * qJD(2);
t782 = t644 * t650;
t685 = t728 * t782;
t743 = qJD(1) * t685 + qJD(3) * t579;
t756 = qJD(4) * t647;
t479 = t578 * t756 + t579 * t730 + t647 * t653 + t796 * t743;
t533 = -t796 * t578 + t579 * t647;
t646 = sin(qJ(6));
t649 = cos(qJ(6));
t757 = qJD(2) * t648;
t736 = t644 * t757;
t705 = qJD(1) * t736;
t754 = qJD(6) * t649;
t744 = t646 * t479 + t533 * t754 + t649 * t705;
t755 = qJD(6) * t646;
t443 = t614 * t755 + t744;
t440 = t443 * t649;
t512 = t533 * t646 - t614 * t649;
t669 = -t649 * t479 + t646 * t705;
t444 = qJD(6) * t512 + t669;
t478 = -t578 * t730 + t579 * t756 + t647 * t743 - t796 * t653;
t473 = t649 * t478;
t677 = t647 * t578 + t579 * t796;
t781 = t646 * t478;
t799 = t677 ^ 2;
t803 = qJD(6) + t677;
t817 = t803 * t649;
t818 = t803 * t646;
t786 = t614 * t646;
t510 = -t649 * t533 - t786;
t819 = t510 * t803;
t821 = t799 * MDP(19) + MDP(22) * t705 + (-t614 * t677 - t479) * MDP(21) + (-t512 * t818 + t440) * MDP(29) + (-t803 * t818 - t473) * MDP(31) + (-t512 * t817 + (-t443 + t819) * t646 - t649 * t444) * MDP(30) + (-t803 * t817 + t781) * MDP(32);
t793 = pkin(5) * t533;
t790 = qJ(5) * t533;
t597 = pkin(8) * t737 + t648 * t747;
t562 = pkin(9) * t717 + t597;
t792 = pkin(9) * t648;
t592 = (-pkin(2) * t650 - pkin(1) - t792) * t644;
t574 = qJD(1) * t592;
t522 = -t562 * t795 + t797 * t574;
t504 = -t579 * pkin(10) + t522;
t496 = -t622 * pkin(3) + t504;
t666 = -t562 * t797 - t574 * t795;
t505 = t578 * pkin(10) - t666;
t503 = t796 * t505;
t451 = t647 * t496 + t503;
t448 = qJ(5) * t614 - t451;
t436 = -t448 - t793;
t820 = t436 * t803;
t816 = t623 * t730 - t624 * t756 + t647 * t824 - t823 * t796;
t607 = t647 * t797 + t795 * t796;
t769 = (-t737 + t802) * t607;
t681 = -pkin(3) * t822 - t597;
t732 = qJD(3) * t797;
t805 = -t683 + t732;
t780 = t647 * t505;
t450 = -t796 * t496 + t780;
t753 = -qJD(5) - t450;
t641 = t644 ^ 2;
t749 = qJD(1) * qJD(2);
t813 = -0.2e1 * t641 * t749;
t812 = pkin(4) * t677;
t811 = pkin(5) * t677;
t810 = MDP(5) * (t648 ^ 2 - t650 ^ 2);
t809 = t436 * t677;
t561 = -pkin(2) * t717 - t594;
t538 = -t578 * pkin(3) + t561;
t657 = -qJ(5) * t677 + t538;
t465 = t533 * pkin(4) + t657;
t789 = t465 * t677;
t808 = t538 * t677;
t798 = pkin(4) + pkin(11);
t807 = t677 * t798;
t457 = t647 * t504 + t503;
t703 = pkin(3) * t756 - t457;
t458 = t504 * t796 - t780;
t777 = -pkin(3) * t730 - qJD(5) + t458;
t776 = qJ(5) * t738 - t816;
t806 = -qJ(5) * t770 + qJD(5) * t607 - t681;
t567 = t647 * t623 + t624 * t796;
t804 = -qJD(4) * t567 + t823 * t647 + t796 * t824;
t752 = t811 - t753;
t794 = pkin(1) * t648;
t591 = pkin(8) * t782 + (pkin(9) + t794) * t645;
t711 = t648 * t742;
t663 = -t645 * t797 + t711;
t656 = -qJD(3) * t663 + t650 * t707;
t596 = qJD(2) * t678;
t783 = t644 * t648;
t633 = pkin(8) * t783;
t791 = t650 * pkin(1);
t598 = (t645 * t791 - t633) * qJD(2);
t698 = t797 * t596 - t598 * t795;
t480 = pkin(3) * t736 - pkin(10) * t656 - t591 * t732 - t592 * t731 + t698;
t600 = t645 * t795 + t783 * t797;
t655 = qJD(3) * t600 + t685;
t664 = -t591 * t731 + t592 * t732 + t795 * t596 + t797 * t598;
t484 = -pkin(10) * t655 + t664;
t697 = -t591 * t795 + t797 * t592;
t509 = -pkin(3) * t782 - t600 * pkin(10) + t697;
t768 = t797 * t591 + t795 * t592;
t517 = -pkin(10) * t663 + t768;
t670 = -t647 * t480 - t796 * t484 - t509 * t730 + t517 * t756;
t429 = -t644 * (qJ(5) * t757 - qJD(5) * t650) + t670;
t434 = t614 * t798 + t752;
t445 = t533 * t798 + t657;
t427 = t434 * t649 - t445 * t646;
t428 = t434 * t646 + t445 * t649;
t800 = MDP(18) * t677 - MDP(19) * t533 + t538 * MDP(24) - t465 * MDP(27) + t512 * MDP(31) - t510 * MDP(32) + MDP(33) * t803 + t427 * MDP(34) - t428 * MDP(35);
t652 = qJD(1) ^ 2;
t617 = qJ(5) * t705;
t586 = qJD(1) * t596;
t587 = qJD(1) * t598;
t699 = t797 * t586 - t587 * t795;
t464 = pkin(3) * t705 - pkin(10) * t653 - t562 * t732 - t574 * t731 + t699;
t665 = -t562 * t731 + t574 * t732 + t795 * t586 + t797 * t587;
t470 = -pkin(10) * t743 + t665;
t713 = -t647 * t464 - t796 * t470 - t496 * t730 + t505 * t756;
t691 = qJD(5) * t614 + t713;
t424 = -t617 + t691;
t419 = -pkin(5) * t479 - t424;
t418 = t419 * t649;
t606 = t647 * t795 - t710;
t787 = t606 * t646;
t638 = -pkin(3) * t796 - pkin(4);
t635 = -pkin(11) + t638;
t785 = t635 * t478;
t784 = t641 * t652;
t779 = t798 * t478;
t778 = t419 * t646 + t436 * t754;
t775 = -pkin(4) * t738 + t804;
t774 = -pkin(5) * t769 - t776;
t773 = -pkin(4) * t769 + t806;
t772 = t647 * t509 + t796 * t517;
t588 = t597 * qJD(2);
t735 = qJD(2) * t782;
t599 = t645 * pkin(1) * t757 + pkin(8) * t735;
t766 = -t777 + t811;
t566 = -t623 * t796 + t647 * t624;
t759 = qJD(2) * t566;
t758 = qJD(2) * t567;
t745 = t649 * t782;
t740 = t650 * t795;
t739 = t641 * t760;
t715 = t798 * t783;
t690 = qJD(2) * t715;
t712 = -t796 * t464 + t647 * t470 + t496 * t756 + t505 * t730;
t421 = -pkin(5) * t478 - qJD(1) * t690 + t712;
t521 = pkin(3) * t743 + t588;
t659 = t478 * qJ(5) - qJD(5) * t677 + t521;
t425 = t479 * t798 + t659;
t726 = t649 * t421 - t425 * t646;
t725 = t646 * t738 + t649 * t769;
t724 = -t646 * t769 + t649 * t738;
t693 = pkin(3) * t579 + t790;
t719 = -qJD(6) * t635 + t693 + t807;
t718 = qJD(6) * t798 + t790 + t807;
t714 = t641 * t650 * t648 * MDP(4);
t639 = -pkin(3) * t797 - pkin(2);
t706 = t622 * t731;
t704 = t703 + t793;
t702 = pkin(1) * t813;
t701 = t509 * t796 - t647 * t517;
t672 = -t607 * qJ(5) + t639;
t539 = t606 * t798 + t672;
t695 = pkin(5) * t770 - qJD(1) * t715 + qJD(6) * t539 + t804;
t543 = t607 * pkin(5) + t566;
t694 = -qJD(6) * t543 - t769 * t798 + t806;
t692 = pkin(4) * t705;
t689 = t421 * t646 + t425 * t649;
t467 = pkin(4) * t782 - t701;
t661 = t647 * t663;
t549 = t600 * t796 - t661;
t446 = t549 * pkin(5) + pkin(11) * t782 + t467;
t660 = t796 * t663;
t548 = t600 * t647 + t660;
t550 = pkin(3) * t711 + t633 + (t639 - t791) * t645;
t654 = -t549 * qJ(5) + t550;
t460 = t548 * t798 + t654;
t688 = t446 * t649 - t460 * t646;
t687 = t446 * t646 + t460 * t649;
t466 = qJ(5) * t782 - t772;
t529 = t548 * t649 + t646 * t782;
t673 = -t451 * t614 - t712;
t671 = -t480 * t796 + t647 * t484 + t509 * t756 + t517 * t730;
t668 = t606 * t755 - t725;
t667 = t606 * t754 - t724;
t662 = -t533 * t614 - t478;
t426 = -t692 + t712;
t658 = qJD(3) * t666 + t699;
t540 = pkin(3) * t655 + t599;
t489 = qJD(4) * t660 + t600 * t756 + t647 * t655 - t656 * t796;
t490 = -qJD(4) * t661 + t600 * t730 + t647 * t656 + t655 * t796;
t433 = t490 * pkin(4) + t489 * qJ(5) - t549 * qJD(5) + t540;
t636 = pkin(3) * t647 + qJ(5);
t590 = t633 + (-pkin(2) - t791) * t645;
t551 = t606 * pkin(4) + t672;
t544 = -t606 * pkin(5) + t567;
t530 = t548 * t646 - t745;
t491 = t790 + t812;
t485 = t548 * pkin(4) + t654;
t481 = t693 + t812;
t471 = t478 * t607;
t455 = t478 * t549;
t454 = qJD(6) * t529 + t490 * t646 + t649 * t736;
t453 = -t649 * t490 - qJD(6) * t745 + (qJD(6) * t548 + t736) * t646;
t449 = -pkin(5) * t548 - t466;
t447 = pkin(4) * t614 - t753;
t438 = t451 - t793;
t432 = t479 * pkin(4) + t659;
t431 = -pkin(4) * t736 + t671;
t430 = t490 * pkin(11) + t433;
t423 = -pkin(5) * t490 - t429;
t422 = -t489 * pkin(5) + t671 - t690;
t416 = -qJD(6) * t428 + t726;
t415 = qJD(6) * t427 + t689;
t1 = [0.2e1 * t714 * t749 + (t561 * t656 + t599 * t579 + t588 * t600 + t590 * t653 + t622 * t664 + t665 * t782 + t666 * t736 - t705 * t768) * MDP(17) + (-(qJD(6) * t688 + t422 * t646 + t430 * t649) * t803 + t687 * t478 - t415 * t549 + t428 * t489 + t423 * t512 + t449 * t443 + t419 * t530 + t436 * t454) * MDP(35) + ((-qJD(6) * t687 + t422 * t649 - t430 * t646) * t803 - t688 * t478 + t416 * t549 - t427 * t489 + t423 * t510 + t449 * t444 - t419 * t529 + t436 * t453) * MDP(34) + (-t489 * t803 - t455) * MDP(33) + (t443 * t549 + t454 * t803 - t478 * t530 - t489 * t512) * MDP(31) + (-t444 * t549 - t453 * t803 - t478 * t529 + t489 * t510) * MDP(32) + (-(-qJD(3) * t768 + t698) * t622 + t697 * t705 - t658 * t782 + t522 * t736 - t599 * t578 + t590 * t743 + t588 * t663 + t561 * t655) * MDP(16) + (t578 * t656 - t579 * t655 - t600 * t743 - t653 * t663) * MDP(12) + (-t622 * t644 - t739) * MDP(15) * t757 + (-t614 * t644 - t739) * MDP(22) * t757 + (t490 * t614 + (t479 * t650 + (-qJD(1) * t548 - t533) * t757) * t644) * MDP(21) + (t671 * t614 + t540 * t533 + t550 * t479 + t521 * t548 + t538 * t490 + (t712 * t650 + (qJD(1) * t701 - t450) * t757) * t644) * MDP(23) + (-t431 * t614 - t432 * t548 - t433 * t533 - t465 * t490 - t479 * t485 + (-t426 * t650 + (qJD(1) * t467 + t447) * t757) * t644) * MDP(26) + (t424 * t466 + t426 * t467 + t429 * t448 + t431 * t447 + t432 * t485 + t433 * t465) * MDP(28) + (-t588 * t645 - t599 * t717 + t648 * t702) * MDP(9) + (-t587 * t645 - t598 * t717 + t650 * t702) * MDP(10) + ((t579 * t797 + t600 * t733) * t735 + (t578 * t600 - t579 * t663) * qJD(3)) * MDP(11) + (MDP(6) * t735 - MDP(7) * t736) * (qJD(2) + 0.2e1 * t761) + (-t489 * t677 - t455) * MDP(18) + (t489 * t614 + (t478 * t650 + (qJD(1) * t549 + t677) * t757) * t644) * MDP(20) + (t429 * t614 - t432 * t549 - t433 * t677 + t465 * t489 + t478 * t485 + (t424 * t650 + (-qJD(1) * t466 - t448) * t757) * t644) * MDP(27) + (-t670 * t614 + t540 * t677 - t550 * t478 + t521 * t549 - t538 * t489 + (-t713 * t650 + (-qJD(1) * t772 - t451) * t757) * t644) * MDP(24) + (t424 * t548 + t426 * t549 + t429 * t533 + t431 * t677 - t447 * t489 + t448 * t490 + t466 * t479 - t467 * t478) * MDP(25) + (t478 * t548 - t479 * t549 + t489 * t533 - t490 * t677) * MDP(19) + (t443 * t529 - t444 * t530 - t453 * t512 - t454 * t510) * MDP(30) + (t443 * t530 + t454 * t512) * MDP(29) + t810 * t813 + (t579 * t736 + t600 * t705 - t622 * t656 - t653 * t782) * MDP(13) + (t578 * t736 + t622 * t655 - t663 * t705 + t743 * t782) * MDP(14); (t622 * MDP(15) + t614 * MDP(22) - t666 * MDP(17) + (-qJD(2) * t606 + t533) * MDP(21) + (t450 - t759) * MDP(23) + (qJD(2) * t607 - t677) * MDP(20) + (t451 - t758) * MDP(24) + (-t447 + t759) * MDP(26) + (t448 + t758) * MDP(27)) * t738 + ((t539 * t649 + t543 * t646) * t478 - t415 * t607 + t544 * t443 + t419 * t787 + (t646 * t695 + t649 * t694) * t803 + t774 * t512 + t770 * t428 + t667 * t436) * MDP(35) + (-t770 * t803 - t471) * MDP(33) + (-t444 * t607 - t473 * t606 + t510 * t770 - t668 * t803) * MDP(32) + (-(-t539 * t646 + t543 * t649) * t478 + t416 * t607 + t544 * t444 - t606 * t418 + (t646 * t694 - t649 * t695) * t803 + t774 * t510 - t770 * t427 + t668 * t436) * MDP(34) + (t443 * t607 - t512 * t770 - t606 * t781 + t667 * t803) * MDP(31) + (-t714 + (-MDP(6) * t650 + MDP(7) * t648) * t644 * t645) * t652 + (t805 * t578 + t579 * t822 + t653 * t797 - t743 * t795) * MDP(12) + (t579 * t805 + t653 * t795) * MDP(11) + (t725 * t512 + t724 * t510 + (t440 - t444 * t646 + (-t510 * t649 - t512 * t646) * qJD(6)) * t606) * MDP(30) + (t443 * t787 + t512 * t667) * MDP(29) + (-t622 * t732 + (t622 * t741 + (t728 - t579) * t648) * t762) * MDP(13) + (t706 + (-t622 * t740 + (t729 - t578) * t648) * t762) * MDP(14) + (-pkin(2) * t743 - t588 * t797 + t696 * t622 + t597 * t578 + (t561 * t795 + t622 * t748) * qJD(3) + (-t561 * t740 + (-pkin(9) * t728 - t522) * t648) * t762) * MDP(16) + (t478 * t606 - t479 * t607 + t533 * t770 - t677 * t769) * MDP(19) + (-t677 * t770 - t471) * MDP(18) + (t424 * t606 + t426 * t607 - t447 * t770 + t448 * t769 - t478 * t566 - t479 * t567 + t533 * t776 - t677 * t775) * MDP(25) + (-t639 * t478 + t521 * t607 - t770 * t538 + t677 * t681) * MDP(24) + (-t432 * t607 + t770 * t465 + t478 * t551 + t677 * t773) * MDP(27) + (-t424 * t567 + t426 * t566 + t432 * t551 - t447 * t775 + t448 * t776 - t465 * t773) * MDP(28) + (pkin(8) * t705 + t594 * t717 + (-t645 * t749 + t784) * t791) * MDP(10) + (t597 * t717 + t784 * t794 - t588) * MDP(9) + (t770 * MDP(20) + t769 * MDP(21) - t804 * MDP(23) + MDP(24) * t816 + t775 * MDP(26) + t776 * MDP(27)) * t614 + (-pkin(2) * t653 - pkin(9) * t706 + t561 * t805 - t597 * t579 + t588 * t795 - t622 * t767 - t682 * t792) * MDP(17) + t784 * t810 + (t639 * t479 + t521 * t606 + t681 * t533 + t769 * t538) * MDP(23) + (-t432 * t606 - t769 * t465 - t479 * t551 + t773 * t533) * MDP(26); (-t561 * t579 + t622 * t666 + t658) * MDP(16) + (t636 * t444 + (-t785 + t809) * t649 + t766 * t510 + (t646 * t719 + t649 * t704) * t803 + t778) * MDP(34) - (MDP(20) * t614 - MDP(25) * t447 - t800) * t533 + t821 + (-t457 * t614 - t808 + (-t533 * t579 + t614 * t756 + t705 * t796) * pkin(3) - t712) * MDP(23) + (-t579 * t622 - t743) * MDP(14) + MDP(15) * t705 + (t636 * t443 + t418 + t719 * t817 + t766 * t512 + (-t704 * t803 + t785 - t820) * t646) * MDP(35) + (-t458 * t614 + (-t579 * t677 + t614 * t730 - t647 * t705) * pkin(3) + t713) * MDP(24) + (t481 * t677 + t614 * t777 + t636 * t705 - t424) * MDP(27) + (-t424 * t636 + t426 * t638 + t447 * t703 + t448 * t777 - t465 * t481) * MDP(28) + (t578 * t622 + t653) * MDP(13) + (t789 + t481 * t533 - t703 * t614 + (-pkin(4) + t638) * t705 + t712) * MDP(26) + (-t578 ^ 2 + t579 ^ 2) * MDP(12) - t478 * MDP(20) + (-t522 * t622 - t561 * t578 - t665) * MDP(17) + (-t478 * t638 - t479 * t636 + t533 * t777 + (-t448 + t703) * t677) * MDP(25) - t579 * t578 * MDP(11); t662 * MDP(20) + (t673 - t808) * MDP(23) + (t450 * t614 + t713) * MDP(24) + (pkin(4) * t478 - qJ(5) * t479 + (-t448 - t451) * t677) * MDP(25) + (-t673 - 0.2e1 * t692 + t789) * MDP(26) + (t491 * t677 + t614 * t753 + 0.2e1 * t617 - t691) * MDP(27) + (-pkin(4) * t426 - qJ(5) * t424 - t447 * t451 + t448 * t753 - t465 * t491) * MDP(28) + (qJ(5) * t444 + (t779 + t809) * t649 + (-t438 * t649 + t646 * t718) * t803 + t752 * t510 + t778) * MDP(34) + (qJ(5) * t443 + t418 + t718 * t817 + t752 * t512 + (t438 * t803 - t779 - t820) * t646) * MDP(35) + ((t447 + t753) * MDP(25) + t491 * MDP(26) + t800) * t533 + t821; t662 * MDP(25) + (-t533 * t677 + t705) * MDP(26) + (-t614 ^ 2 - t799) * MDP(27) + (-t448 * t614 + t426 + t789) * MDP(28) + (t510 * t614 - t473) * MDP(34) + (t512 * t614 + t781) * MDP(35) + (-MDP(34) * t818 - MDP(35) * t817) * t803; t512 * t510 * MDP(29) + (-t510 ^ 2 + t512 ^ 2) * MDP(30) + (t744 + t819) * MDP(31) + (t512 * t803 - t669) * MDP(32) - t478 * MDP(33) + (t428 * t803 - t436 * t512 + t726) * MDP(34) + (t427 * t803 + t436 * t510 - t689) * MDP(35) + (MDP(31) * t786 - MDP(32) * t512 - MDP(34) * t428 - MDP(35) * t427) * qJD(6);];
tauc  = t1;
