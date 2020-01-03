% Calculate inertial parameters regressor of inverse dynamics cutting torque vector with Newton-Euler for
% S4RRPR10
% Use Code from Maple symbolic Code Generation
%
% Input:
% qJ [4x1]
%   Generalized joint coordinates (joint angles)
% qJD [4x1]
%   Generalized joint velocities
% qJDD [4x1]
%   Generalized joint accelerations
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [6x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,d1,d2,d4]';
%
% Output:
% m_new_reg [(3*5)x(%Nl%*10)]
%   inertial parameter regressor of inverse dynamics cutting torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 17:12
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function m_new_reg = S4RRPR10_invdynm_fixb_reg2_snew_vp(qJ, qJD, qJDD, g, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(4,1),zeros(4,1),zeros(3,1),zeros(6,1)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4RRPR10_invdynm_fixb_reg2_snew_vp: qJ has to be [4x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [4 1]), ...
  'S4RRPR10_invdynm_fixb_reg2_snew_vp: qJD has to be [4x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [4 1]), ...
  'S4RRPR10_invdynm_fixb_reg2_snew_vp: qJDD has to be [4x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S4RRPR10_invdynm_fixb_reg2_snew_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [6 1]), ...
  'S4RRPR10_invdynm_fixb_reg2_snew_vp: pkin has to be [6x1] (double)');

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_m_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 17:12:03
% EndTime: 2019-12-31 17:12:08
% DurationCPUTime: 4.76s
% Computational Cost: add. (9831->416), mult. (20991->463), div. (0->0), fcn. (11528->6), ass. (0->270)
t697 = sin(qJ(2));
t693 = t697 ^ 2;
t702 = qJD(1) ^ 2;
t687 = t693 * t702;
t789 = qJD(2) ^ 2;
t673 = -t687 - t789;
t700 = cos(qJ(2));
t755 = t700 * t702;
t678 = t697 * t755;
t668 = qJDD(2) - t678;
t756 = t700 * t668;
t622 = t697 * t673 + t756;
t748 = qJD(1) * qJD(2);
t684 = t700 * t748;
t745 = t697 * qJDD(1);
t657 = 0.2e1 * t684 + t745;
t698 = sin(qJ(1));
t701 = cos(qJ(1));
t823 = pkin(4) * (t701 * t622 - t698 * t657);
t822 = pkin(4) * (t698 * t622 + t701 * t657);
t790 = t700 ^ 2;
t688 = t790 * t702;
t675 = -t688 - t789;
t667 = qJDD(2) + t678;
t757 = t700 * t667;
t610 = t697 * t675 + t757;
t821 = pkin(1) * t610;
t768 = t697 * t667;
t619 = -t700 * t675 + t768;
t686 = t700 * qJDD(1);
t735 = t697 * t748;
t660 = t686 - 0.2e1 * t735;
t820 = pkin(4) * (t701 * t619 + t698 * t660);
t819 = pkin(4) * (t698 * t619 - t701 * t660);
t818 = pkin(5) * t610;
t817 = pkin(5) * t622;
t672 = -t687 + t789;
t621 = -t697 * t672 + t757;
t743 = t701 * qJDD(1);
t816 = t698 * t621 - t697 * t743;
t744 = t698 * qJDD(1);
t815 = t701 * t621 + t697 * t744;
t812 = 2 * qJD(3);
t767 = t697 * t668;
t614 = -t700 * t673 + t767;
t811 = pkin(1) * t614;
t810 = pkin(5) * t614;
t809 = pkin(5) * t619;
t696 = sin(qJ(4));
t699 = cos(qJ(4));
t750 = qJD(1) * t700;
t650 = t696 * qJD(2) + t699 * t750;
t652 = t699 * qJD(2) - t696 * t750;
t605 = t652 * t650;
t658 = t684 + t745;
t643 = qJDD(4) + t658;
t798 = -t605 + t643;
t808 = t696 * t798;
t807 = t699 * t798;
t690 = t697 * g(3);
t671 = t701 * g(1) + t698 * g(2);
t639 = -t702 * pkin(1) + qJDD(1) * pkin(5) - t671;
t776 = qJ(3) * t697;
t785 = pkin(2) * t700;
t721 = -t776 - t785;
t727 = t702 * t721 + t639;
t706 = -t789 * pkin(2) + t727 * t700 - t690;
t705 = qJD(2) * t812 + t706;
t806 = t705 + (qJDD(2) + t668) * qJ(3) - pkin(2) * t673;
t674 = t688 - t789;
t620 = -t700 * t674 + t767;
t803 = t698 * t620 + t700 * t743;
t802 = t701 * t620 - t698 * t686;
t717 = t658 + t684;
t801 = t717 * qJ(3);
t659 = t686 - t735;
t749 = t697 * qJD(1);
t669 = pkin(3) * t749 - qJD(2) * pkin(6);
t670 = t698 * g(1) - t701 * g(2);
t713 = -qJDD(1) * pkin(1) - t670;
t797 = -pkin(2) * t735 + t749 * t812;
t704 = t713 - t797 - t801;
t788 = pkin(2) + pkin(6);
t531 = -t669 * t749 + (-pkin(3) * t790 - pkin(5)) * t702 - t788 * t659 + t704;
t780 = t700 * g(3);
t714 = -qJDD(2) * pkin(2) - t789 * qJ(3) + qJDD(3) + t780;
t542 = -qJDD(2) * pkin(6) + (t658 - t684) * pkin(3) + (-pkin(6) * t755 + t727) * t697 + t714;
t504 = t696 * t531 - t699 * t542;
t505 = t699 * t531 + t696 * t542;
t800 = -t699 * t504 + t696 * t505;
t596 = -t650 * qJD(4) + t699 * qJDD(2) - t696 * t659;
t680 = qJD(4) + t749;
t635 = t680 * t650;
t799 = t596 - t635;
t609 = t697 * t674 + t756;
t578 = t727 * t697 + t714;
t796 = -pkin(2) * t667 - qJ(3) * t675 + t578;
t746 = qJDD(2) * qJ(3);
t541 = t746 + t659 * pkin(3) - pkin(6) * t688 + (t812 + t669) * qJD(2) + t706;
t794 = qJ(3) * t541 - t788 * t800;
t726 = -t696 * qJDD(2) - t699 * t659;
t568 = (-qJD(4) + t680) * t652 + t726;
t571 = t596 + t635;
t525 = t696 * t568 - t699 * t571;
t641 = t650 ^ 2;
t642 = t652 ^ 2;
t587 = -t641 - t642;
t793 = qJ(3) * t587 - t788 * t525 - t800;
t676 = t680 ^ 2;
t597 = -t676 - t641;
t539 = t696 * t597 + t807;
t636 = t680 * t652;
t710 = t652 * qJD(4) - t726;
t566 = t636 + t710;
t772 = t696 * t541;
t792 = qJ(3) * t566 - t788 * t539 + t772;
t537 = t699 * t541;
t741 = -t642 - t676;
t589 = t605 + t643;
t771 = t696 * t589;
t546 = t699 * t741 - t771;
t791 = qJ(3) * t799 - t788 * t546 + t537;
t784 = pkin(3) * t800;
t783 = pkin(3) * t525;
t782 = pkin(3) * t541;
t740 = t693 + t790;
t662 = t740 * qJDD(1);
t665 = t687 + t688;
t781 = pkin(4) * (t698 * t662 + t701 * t665);
t779 = t702 * pkin(5);
t775 = t680 * t696;
t774 = t680 * t699;
t638 = -t713 + t779;
t770 = t697 * t638;
t769 = t697 * t657;
t760 = t699 * t589;
t759 = t700 * t638;
t758 = t700 * t660;
t752 = pkin(1) * t665 + pkin(5) * t662;
t739 = t652 * t775;
t738 = t697 * t605;
t737 = t700 * t605;
t625 = t697 * t639 + t780;
t626 = t700 * t639 - t690;
t574 = t697 * t625 + t700 * t626;
t728 = -t698 * t670 - t701 * t671;
t725 = t698 * t678;
t724 = t701 * t678;
t664 = -t698 * t702 + t743;
t723 = -pkin(4) * t664 - t698 * g(3);
t575 = t705 + t746;
t722 = -pkin(2) * t578 + qJ(3) * t575;
t720 = pkin(2) * t697 - qJ(3) * t700;
t719 = pkin(3) * t566 + t537;
t718 = pkin(3) * t799 - t772;
t493 = t696 * t504 + t699 * t505;
t573 = t700 * t625 - t697 * t626;
t613 = t700 * t672 + t768;
t715 = t701 * t670 - t698 * t671;
t712 = -pkin(3) * t539 + t504;
t709 = pkin(3) * t587 - t493;
t708 = -pkin(3) * t546 + t505;
t703 = t659 * pkin(2) + t638 + t797;
t666 = t687 - t688;
t663 = t701 * t702 + t744;
t653 = t720 * qJDD(1);
t645 = t740 * t748;
t637 = -pkin(4) * t663 + t701 * g(3);
t633 = -t642 + t676;
t632 = t641 - t676;
t631 = t698 * qJDD(2) + t701 * t645;
t630 = t700 * t658 - t693 * t748;
t629 = -t701 * qJDD(2) + t698 * t645;
t628 = -t697 * t659 - t790 * t748;
t627 = t650 * t774;
t608 = t717 * t697;
t607 = (t659 - t735) * t700;
t603 = t642 - t641;
t601 = pkin(4) * (t701 * t662 - t698 * t665);
t600 = t758 - t769;
t599 = t700 * t657 + t697 * t660;
t595 = t701 * t630 - t725;
t594 = t701 * t628 + t725;
t593 = t698 * t630 + t724;
t592 = t698 * t628 - t724;
t584 = -t759 + t810;
t583 = -t770 - t818;
t582 = t627 - t739;
t581 = (-t650 * t696 - t652 * t699) * t680;
t580 = t701 * t600 + t698 * t666;
t579 = t698 * t600 - t701 * t666;
t577 = t626 + t811;
t576 = t625 - t821;
t567 = -t636 + t710;
t565 = pkin(1) * t660 + t759 - t809;
t564 = -pkin(1) * t657 - t770 - t817;
t561 = qJ(3) * t665 + t578;
t560 = pkin(2) * t665 + t575;
t559 = t696 * t596 + t652 * t774;
t558 = t696 * t710 + t627;
t557 = -t699 * t596 + t739;
t556 = -t650 * t775 + t699 * t710;
t555 = t703 + t801;
t554 = t697 * t581 + t700 * t643;
t553 = -t700 * t581 + t697 * t643;
t552 = t696 * t633 - t807;
t551 = t696 * t632 + t760;
t550 = -t699 * t632 + t771;
t549 = t699 * t633 + t808;
t548 = pkin(1) * t638 + pkin(5) * t574;
t547 = -t696 * t741 - t760;
t545 = -t779 + (-t659 - t660) * pkin(2) + t704;
t544 = (t657 + t717) * qJ(3) + t703;
t543 = t574 + t752;
t540 = t699 * t597 - t808;
t536 = -t796 + t821;
t535 = t697 * t559 + t737;
t534 = -t697 * t556 - t737;
t533 = -t700 * t559 + t738;
t532 = t700 * t556 - t738;
t530 = -t806 - t811;
t529 = t700 * t575 + t697 * t578;
t528 = t697 * t575 - t700 * t578;
t527 = t699 * t568 + t696 * t571;
t526 = t699 * t566 + t696 * t799;
t524 = -t696 * t566 + t699 * t799;
t523 = -pkin(2) * t769 + t700 * t544 - t810;
t522 = -qJ(3) * t758 - t697 * t545 + t818;
t521 = -t697 * t560 + t700 * t561;
t520 = t697 * t549 + t700 * t571;
t519 = t697 * t551 - t700 * t567;
t518 = -t700 * t549 + t697 * t571;
t517 = -t700 * t551 - t697 * t567;
t516 = t817 + t697 * t544 + (pkin(1) + t785) * t657;
t515 = t809 + t700 * t545 + (-pkin(1) - t776) * t660;
t514 = t697 * t546 + t700 * t799;
t513 = -t700 * t546 + t697 * t799;
t512 = t697 * t539 + t700 * t566;
t511 = -t700 * t539 + t697 * t566;
t510 = t700 * t560 + t697 * t561 + t752;
t509 = t697 * t524 + t700 * t603;
t508 = -t700 * t524 + t697 * t603;
t507 = t697 * t525 + t700 * t587;
t506 = -t700 * t525 + t697 * t587;
t501 = -pkin(1) * t528 - t722;
t500 = -qJ(3) * t527 + t783;
t499 = -pkin(5) * t528 - t720 * t555;
t498 = -t788 * t547 + t718;
t497 = pkin(5) * t529 + (pkin(1) - t721) * t555;
t496 = -t788 * t540 + t719;
t495 = -qJ(3) * t547 - t708;
t494 = -qJ(3) * t540 - t712;
t491 = -pkin(1) * t513 - t791;
t490 = -pkin(1) * t511 - t792;
t489 = t700 * t541 + t697 * t800;
t488 = t697 * t541 - t700 * t800;
t487 = -t788 * t527 + t709;
t486 = -qJ(3) * t493 + t784;
t485 = -pkin(5) * t513 + t700 * t495 - t697 * t498;
t484 = -pkin(5) * t511 + t700 * t494 - t697 * t496;
t483 = -pkin(1) * t506 - t793;
t482 = -t788 * t493 + t782;
t481 = -pkin(1) * t547 + pkin(5) * t514 + t697 * t495 + t700 * t498;
t480 = -pkin(1) * t540 + pkin(5) * t512 + t697 * t494 + t700 * t496;
t479 = -pkin(5) * t506 - t697 * t487 + t700 * t500;
t478 = -pkin(1) * t527 + pkin(5) * t507 + t700 * t487 + t697 * t500;
t477 = -pkin(1) * t488 - t794;
t476 = -pkin(5) * t488 - t697 * t482 + t700 * t486;
t475 = -pkin(1) * t493 + pkin(5) * t489 + t700 * t482 + t697 * t486;
t1 = [0, 0, 0, 0, 0, 0, 0, -g(3), g(2), 0, 0, 0, t664, 0, -t663, 0, t723, -t637, -t715, -pkin(4) * t715, t595, t580, t815, t594, -t802, t631, -t698 * t576 + t701 * t583 + t819, -t698 * t577 + t701 * t584 + t822, t701 * t573 - t781, -pkin(4) * (t698 * t574 + t701 * t638) - (t698 * pkin(1) - t701 * pkin(5)) * t573, t631, -t815, t802, t595, t580, t594, t701 * t521 - t698 * t653 - t781, t701 * t522 - t698 * t536 - t819, t701 * t523 - t698 * t530 - t822, t701 * t499 - t698 * t501 - pkin(4) * (t698 * t529 + t701 * t555), t701 * t535 - t698 * t557, t701 * t509 - t698 * t526, t701 * t520 - t698 * t552, t701 * t534 + t698 * t558, t701 * t519 - t698 * t550, t701 * t554 - t698 * t582, t701 * t484 - t698 * t490 - pkin(4) * (t698 * t512 - t701 * t540), t701 * t485 - t698 * t491 - pkin(4) * (t698 * t514 - t701 * t547), t701 * t479 - t698 * t483 - pkin(4) * (t698 * t507 - t701 * t527), t701 * t476 - t698 * t477 - pkin(4) * (t698 * t489 - t701 * t493); 0, 0, 0, 0, 0, 0, g(3), 0, -g(1), 0, 0, 0, t663, 0, t664, 0, t637, t723, t728, pkin(4) * t728, t593, t579, t816, t592, -t803, t629, t701 * t576 + t698 * t583 - t820, t701 * t577 + t698 * t584 - t823, t698 * t573 + t601, pkin(4) * (t701 * t574 - t698 * t638) - (-t701 * pkin(1) - t698 * pkin(5)) * t573, t629, -t816, t803, t593, t579, t592, t698 * t521 + t701 * t653 + t601, t698 * t522 + t701 * t536 + t820, t698 * t523 + t701 * t530 + t823, t698 * t499 + t701 * t501 + pkin(4) * (t701 * t529 - t698 * t555), t698 * t535 + t701 * t557, t698 * t509 + t701 * t526, t698 * t520 + t701 * t552, t698 * t534 - t701 * t558, t698 * t519 + t701 * t550, t698 * t554 + t701 * t582, t698 * t484 + t701 * t490 + pkin(4) * (t701 * t512 + t698 * t540), t698 * t485 + t701 * t491 + pkin(4) * (t701 * t514 + t698 * t547), t698 * t479 + t701 * t483 + pkin(4) * (t701 * t507 + t698 * t527), t698 * t476 + t701 * t477 + pkin(4) * (t701 * t489 + t698 * t493); 0, 0, 0, 0, 0, 0, -g(2), g(1), 0, 0, 0, 0, 0, 0, 0, qJDD(1), t670, t671, 0, 0, t608, t599, t613, t607, t609, 0, t565, t564, t543, t548, 0, -t613, -t609, t608, t599, t607, t510, t515, t516, t497, t533, t508, t518, t532, t517, t553, t480, t481, t478, t475; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, qJDD(1), 0, -t702, 0, 0, -g(3), -t670, 0, t630, t600, t621, t628, -t620, t645, t583, t584, t573, pkin(5) * t573, t645, -t621, t620, t630, t600, t628, t521, t522, t523, t499, t535, t509, t520, t534, t519, t554, t484, t485, t479, t476; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t702, 0, qJDD(1), 0, g(3), 0, -t671, 0, t678, -t666, -t745, -t678, -t686, -qJDD(2), t576, t577, 0, pkin(1) * t573, -qJDD(2), t745, t686, t678, -t666, -t678, t653, t536, t530, t501, t557, t526, t552, -t558, t550, t582, t490, t491, t483, t477; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, qJDD(1), t670, t671, 0, 0, t608, t599, t613, t607, t609, 0, t565, t564, t543, t548, 0, -t613, -t609, t608, t599, t607, t510, t515, t516, t497, t533, t508, t518, t532, t517, t553, t480, t481, t478, t475; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t658, t660, t667, -t684, t674, t684, 0, -t638, t625, 0, t684, -t667, -t674, t658, t660, -t684, t561, -qJ(3) * t660, t544, qJ(3) * t555, t605, t603, t571, -t605, -t567, t643, t494, t495, t500, t486; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t735, t657, t672, t659, t668, -t735, t638, 0, t626, 0, -t735, -t672, -t668, t735, t657, t659, t560, t545, pkin(2) * t657, pkin(2) * t555, -t559, -t524, -t549, t556, -t551, -t581, t496, t498, t487, t482; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t678, t666, t745, t678, t686, qJDD(2), -t625, -t626, 0, 0, qJDD(2), -t745, -t686, -t678, t666, t678, -t653, t796, t806, t722, -t557, -t526, -t552, t558, -t550, -t582, t792, t791, t793, t794; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, qJDD(2), -t745, -t686, -t678, t666, t678, 0, t578, t575, 0, -t557, -t526, -t552, t558, -t550, -t582, -pkin(6) * t539 + t772, -pkin(6) * t546 + t537, -pkin(6) * t525 - t800, -pkin(6) * t800; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t684, t667, t674, -t658, -t660, t684, -t578, 0, -t555, 0, -t605, -t603, -t571, t605, t567, -t643, t712, t708, -t783, -t784; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t735, t672, t668, -t735, -t657, -t659, -t575, t555, 0, 0, t559, t524, t549, -t556, t551, t581, pkin(6) * t540 - t719, pkin(6) * t547 - t718, pkin(6) * t527 - t709, pkin(6) * t493 - t782; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t596, -t566, t798, t635, t632, -t635, 0, t541, t504, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t636, t799, t633, -t710, t589, -t636, -t541, 0, t505, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t605, t603, t571, -t605, -t567, t643, -t504, -t505, 0, 0;];
m_new_reg = t1;
