% Calculate inertial parameters regressor of inverse dynamics cutting torque vector with Newton-Euler for
% S5RPPRP3
% Use Code from Maple symbolic Code Generation
%
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% qJDD [5x1]
%   Generalized joint accelerations
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [7x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d4,theta2]';
%
% Output:
% m_new_reg [(3*6)x(%Nl%*10)]
%   inertial parameter regressor of inverse dynamics cutting torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 17:51
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function m_new_reg = S5RPPRP3_invdynm_fixb_reg2_snew_vp(qJ, qJD, qJDD, g, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(5,1),zeros(3,1),zeros(7,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPPRP3_invdynm_fixb_reg2_snew_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RPPRP3_invdynm_fixb_reg2_snew_vp: qJD has to be [5x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [5 1]), ...
  'S5RPPRP3_invdynm_fixb_reg2_snew_vp: qJDD has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RPPRP3_invdynm_fixb_reg2_snew_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [7 1]), ...
  'S5RPPRP3_invdynm_fixb_reg2_snew_vp: pkin has to be [7x1] (double)');

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_m_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 17:51:09
% EndTime: 2019-12-31 17:51:15
% DurationCPUTime: 5.64s
% Computational Cost: add. (9451->401), mult. (18320->426), div. (0->0), fcn. (9047->6), ass. (0->276)
t673 = sin(pkin(7));
t674 = cos(pkin(7));
t682 = qJD(1) ^ 2;
t632 = -t674 * qJDD(1) + t673 * t682;
t670 = g(3) - qJDD(2);
t600 = qJ(2) * t632 - t673 * t670;
t678 = sin(qJ(1));
t680 = cos(qJ(1));
t631 = t673 * qJDD(1) + t674 * t682;
t723 = t680 * t631 - t678 * t632;
t727 = -qJ(2) * t631 + t674 * t670;
t818 = -pkin(5) * t723 + t678 * t600 + t680 * t727;
t681 = qJD(4) ^ 2;
t677 = sin(qJ(4));
t668 = t677 ^ 2;
t772 = t668 * t682;
t645 = -t681 - t772;
t679 = cos(qJ(4));
t754 = t679 * t682;
t648 = t677 * t754;
t640 = qJDD(4) - t648;
t755 = t679 * t640;
t587 = t677 * t645 + t755;
t580 = pkin(3) * t587;
t741 = qJD(1) * qJD(4);
t657 = t677 * t741;
t662 = t679 * qJDD(1);
t627 = t662 - t657;
t672 = qJDD(1) * pkin(2);
t641 = t678 * g(1) - t680 * g(2);
t620 = qJDD(1) * pkin(1) + t641;
t642 = t680 * g(1) + t678 * g(2);
t621 = -t682 * pkin(1) - t642;
t562 = -t674 * t620 + t673 * t621;
t709 = qJDD(3) + t562;
t552 = -t682 * qJ(3) - t672 + t709;
t687 = -qJDD(1) * pkin(6) + t552;
t536 = t679 * t687;
t740 = qJD(1) * qJD(5);
t702 = -qJDD(4) * pkin(4) + 0.2e1 * t679 * t740 - t536;
t716 = pkin(4) * t754 - t670;
t499 = t702 + t677 * (qJ(5) * t741 + t716) + t627 * qJ(5);
t683 = pkin(4) * t640 - t499;
t819 = -t580 - t683;
t669 = t679 ^ 2;
t771 = t669 * t682;
t647 = -t681 - t771;
t639 = qJDD(4) + t648;
t764 = t677 * t639;
t589 = t679 * t647 - t764;
t581 = pkin(3) * t589;
t658 = t679 * t741;
t737 = t677 * qJDD(1);
t626 = -t658 - t737;
t742 = qJD(1) * t679;
t638 = qJD(4) * pkin(4) - qJ(5) * t742;
t744 = t679 * t670 - t677 * t687;
t500 = -pkin(4) * t772 + t626 * qJ(5) - qJD(4) * t638 - 0.2e1 * t677 * t740 - t744;
t699 = pkin(4) * t647 - t500;
t815 = -t581 - t699;
t563 = t673 * t620 + t674 * t621;
t724 = t673 * t562 + t674 * t563;
t520 = t674 * t562 - t673 * t563;
t753 = t680 * t520;
t814 = -t678 * t724 + t753;
t762 = t678 * t520;
t813 = t680 * t724 + t762;
t721 = t678 * t631 + t680 * t632;
t795 = pkin(5) * t721 + t680 * t600 - t678 * t727;
t738 = qJDD(1) * qJ(3);
t689 = -t682 * pkin(2) + t563 + t738;
t739 = qJD(3) * qJD(1);
t736 = 0.2e1 * t739;
t550 = t689 + t736;
t502 = t673 * t550 - t674 * t552;
t725 = t674 * t550 + t673 * t552;
t810 = -t678 * t502 + t680 * t725;
t809 = t680 * t502 + t678 * t725;
t525 = -t677 * t670 - t536;
t479 = -t679 * t525 - t677 * t744;
t625 = 0.2e1 * t658 + t737;
t802 = -pkin(2) * t587 + qJ(3) * t625;
t628 = t662 - 0.2e1 * t657;
t801 = -pkin(2) * t589 + qJ(3) * t628;
t743 = t668 + t669;
t633 = t743 * qJDD(1);
t636 = t743 * t682;
t800 = pkin(2) * t633 - qJ(3) * t636;
t794 = -pkin(2) - pkin(6);
t793 = pkin(1) * t631;
t792 = pkin(1) * t632;
t763 = t677 * t640;
t591 = t679 * t645 - t763;
t791 = pkin(2) * t591;
t756 = t679 * t639;
t594 = -t677 * t647 - t756;
t790 = pkin(2) * t594;
t789 = pkin(3) * t479;
t675 = t682 * pkin(6);
t542 = t550 - t675;
t788 = pkin(3) * t542;
t787 = pkin(3) * t636;
t545 = -t674 * t587 + t673 * t625;
t547 = t673 * t587 + t674 * t625;
t786 = pkin(5) * (t680 * t545 + t678 * t547);
t546 = -t674 * t589 + t673 * t628;
t548 = t673 * t589 + t674 * t628;
t785 = pkin(5) * (t680 * t546 + t678 * t548);
t769 = t674 * t633;
t576 = -t673 * t636 + t769;
t770 = t673 * t633;
t577 = -t674 * t636 - t770;
t784 = pkin(5) * (t680 * t576 + t678 * t577);
t783 = pkin(6) * t479;
t782 = pkin(6) * t587;
t781 = pkin(6) * t589;
t780 = pkin(6) * t633;
t779 = qJ(2) * t545;
t778 = qJ(2) * t546;
t777 = qJ(2) * t576;
t775 = qJ(3) * t591;
t774 = qJ(3) * t594;
t767 = t677 * t499;
t765 = t677 * t542;
t757 = t679 * t499;
t533 = t679 * t542;
t747 = -pkin(2) * t552 + qJ(3) * t550;
t746 = pkin(3) * t625 - pkin(6) * t591;
t745 = pkin(3) * t628 - pkin(6) * t594;
t735 = pkin(4) * t662;
t734 = -t581 - t744;
t731 = -pkin(1) * t591 + qJ(2) * t547;
t730 = -pkin(1) * t594 + qJ(2) * t548;
t462 = t677 * t500 - t757;
t497 = pkin(4) * t499;
t729 = -pkin(3) * t462 + t497;
t728 = t533 - t781;
t719 = -t678 * t641 - t680 * t642;
t718 = t673 * t648;
t717 = t674 * t648;
t715 = -pkin(2) * t479 + qJ(3) * t542 - t783;
t635 = t680 * qJDD(1) - t678 * t682;
t714 = -pkin(5) * t635 - t678 * g(3);
t713 = t525 - t580;
t712 = -t746 - t533;
t711 = -t745 + t765;
t710 = t765 - t782;
t480 = t677 * t525 - t679 * t744;
t707 = t680 * t641 - t678 * t642;
t706 = t728 + t801;
t485 = pkin(4) * t636 - qJ(5) * t737 + t500;
t492 = t716 * t677 + (t627 + t657 + t662) * qJ(5) + t702;
t705 = -t677 * t485 + t679 * t492 + t780;
t684 = -t626 * pkin(4) - qJ(5) * t772 + qJDD(5) - t675 + t689;
t511 = (t638 * t679 + 0.2e1 * qJD(3)) * qJD(1) + t684;
t505 = -qJ(5) * t647 + t511;
t584 = -pkin(4) * t628 - qJ(5) * t639;
t704 = t679 * t505 - t677 * t584 - t781;
t703 = t780 - t479;
t701 = -0.2e1 * t672 + t709;
t700 = t710 + t802;
t698 = t677 * t505 + t679 * t584 - t745;
t697 = t705 + t800;
t696 = t704 + t801;
t695 = t703 + t800;
t490 = -pkin(4) * t625 + qJ(5) * t645 - t638 * t742 - t684 - 0.2e1 * t739;
t693 = -qJ(5) * t763 + t679 * t490 - t746;
t470 = -pkin(4) * t511 + qJ(5) * t500;
t692 = pkin(3) * t511 - qJ(5) * t767 - t679 * t470;
t691 = -pkin(6) * t462 + qJ(5) * t757 - t677 * t470;
t690 = -qJ(5) * t755 - t677 * t490 - t782;
t688 = t563 + 0.2e1 * t738 + t736;
t686 = -pkin(2) * t462 + qJ(3) * t511 + t691;
t685 = t690 + t802;
t646 = t681 - t771;
t644 = -t681 + t772;
t637 = (-t668 + t669) * t682;
t634 = t678 * qJDD(1) + t680 * t682;
t623 = pkin(3) * t633;
t618 = t743 * t741;
t606 = -pkin(5) * t634 + t680 * g(3);
t605 = t623 + t735;
t598 = t677 * t627 + t669 * t741;
t597 = t679 * t626 + t668 * t741;
t596 = t674 * qJDD(4) - t673 * t618;
t595 = t673 * qJDD(4) + t674 * t618;
t593 = -t677 * t646 + t755;
t592 = (t627 - t657) * t679;
t590 = t679 * t644 - t764;
t588 = t679 * t646 + t763;
t586 = t677 * t644 + t756;
t585 = (-t626 + t658) * t677;
t567 = pkin(1) * t576;
t566 = qJ(2) * t577;
t565 = -t679 * t625 - t677 * t628;
t564 = -t677 * t625 + t679 * t628;
t560 = t673 * t597 - t717;
t559 = t673 * t598 + t717;
t558 = -t674 * t597 - t718;
t557 = -t674 * t598 + t718;
t556 = t673 * t588 + t674 * t662;
t555 = t673 * t586 - t674 * t737;
t554 = -t674 * t588 + t673 * t662;
t553 = -t674 * t586 - t673 * t737;
t544 = -t562 - t792;
t543 = -t563 - t793;
t541 = pkin(1) * t546;
t540 = pkin(1) * t545;
t532 = t673 * t564 + t674 * t637;
t531 = -t674 * t564 + t673 * t637;
t530 = t701 + t792;
t529 = t688 + t793;
t528 = -t678 * t595 + t680 * t596;
t527 = t680 * t595 + t678 * t596;
t522 = pkin(5) * (-t678 * t576 + t680 * t577);
t517 = pkin(1) * t520;
t516 = pkin(1) * t670 + qJ(2) * t724;
t515 = -t678 * t558 + t680 * t560;
t514 = -t678 * t557 + t680 * t559;
t513 = t680 * t558 + t678 * t560;
t512 = t680 * t557 + t678 * t559;
t510 = -t678 * t554 + t680 * t556;
t509 = -t678 * t553 + t680 * t555;
t508 = t680 * t554 + t678 * t556;
t507 = t680 * t553 + t678 * t555;
t494 = pkin(5) * (-t678 * t546 + t680 * t548);
t493 = pkin(5) * (-t678 * t545 + t680 * t547);
t489 = -t678 * t531 + t680 * t532;
t488 = t680 * t531 + t678 * t532;
t487 = -t734 - t774;
t486 = -t713 - t775;
t484 = -qJ(2) * t502 + (-pkin(2) * t673 + qJ(3) * t674) * t670;
t483 = qJ(2) * t725 + (pkin(2) * t674 + qJ(3) * t673 + pkin(1)) * t670;
t482 = -t712 - t791;
t481 = -t711 - t790;
t477 = t480 + t787;
t476 = -t774 - t815;
t475 = -t775 - t819;
t474 = t541 + t706;
t473 = t540 + t700;
t472 = t673 * t479 + t674 * t542;
t471 = -t674 * t479 + t673 * t542;
t469 = -pkin(3) * t769 + t673 * t477 - t777;
t468 = -pkin(3) * t770 - t674 * t477 + t566;
t467 = pkin(1) * t502 + t747;
t466 = t567 + t695;
t465 = -t693 - t791;
t464 = -t698 - t790;
t463 = t679 * t500 + t767;
t460 = t679 * t485 + t677 * t492 + t787;
t459 = t540 + t685;
t458 = t541 + t696;
t457 = -qJ(3) * t480 + t789;
t456 = -t673 * t481 + t674 * t487 - t778;
t455 = -t673 * t482 + t674 * t486 - t779;
t454 = t673 * t462 + t674 * t511;
t453 = -t674 * t462 + t673 * t511;
t452 = t567 + t697;
t451 = t794 * t480 + t788;
t450 = t673 * t460 - t674 * t605 - t777;
t449 = -t674 * t460 - t673 * t605 + t566;
t448 = t674 * t481 + t673 * t487 + t730;
t447 = t674 * t482 + t673 * t486 + t731;
t446 = -t673 * t464 + t674 * t476 - t778;
t445 = -t673 * t465 + t674 * t475 - t779;
t444 = t674 * t464 + t673 * t476 + t730;
t443 = t674 * t465 + t673 * t475 + t731;
t442 = pkin(1) * t471 + t715;
t441 = -qJ(3) * t463 - t729;
t440 = -qJ(2) * t471 - t673 * t451 + t674 * t457;
t439 = t794 * t463 + t692;
t438 = -pkin(1) * t480 + qJ(2) * t472 + t674 * t451 + t673 * t457;
t437 = pkin(1) * t453 + t686;
t436 = -qJ(2) * t453 - t673 * t439 + t674 * t441;
t435 = -pkin(1) * t463 + qJ(2) * t454 + t674 * t439 + t673 * t441;
t1 = [0, 0, 0, 0, 0, 0, 0, -g(3), g(2), 0, 0, 0, t635, 0, -t634, 0, t714, -t606, -t707, -pkin(5) * t707, 0, 0, -t721, 0, -t723, 0, t795, -t818, t814, pkin(5) * t814 + qJ(2) * t753 - t678 * t516, 0, t721, t723, 0, 0, 0, -t809, -t795, t818, -pkin(5) * t809 - t678 * t483 + t680 * t484, t514, t489, t510, t515, t509, t528, -t678 * t447 + t680 * t455 - t786, -t678 * t448 + t680 * t456 - t785, -t678 * t468 + t680 * t469 - t784, t680 * t440 - t678 * t438 - pkin(5) * (t680 * t471 + t678 * t472), t514, t489, t510, t515, t509, t528, -t678 * t443 + t680 * t445 - t786, -t678 * t444 + t680 * t446 - t785, -t678 * t449 + t680 * t450 - t784, t680 * t436 - t678 * t435 - pkin(5) * (t680 * t453 + t678 * t454); 0, 0, 0, 0, 0, 0, g(3), 0, -g(1), 0, 0, 0, t634, 0, t635, 0, t606, t714, t719, pkin(5) * t719, 0, 0, t723, 0, -t721, 0, t818, t795, t813, pkin(5) * t813 + qJ(2) * t762 + t680 * t516, 0, -t723, t721, 0, 0, 0, t810, -t818, -t795, pkin(5) * t810 + t680 * t483 + t678 * t484, t512, t488, t508, t513, t507, t527, t680 * t447 + t678 * t455 + t493, t680 * t448 + t678 * t456 + t494, t680 * t468 + t678 * t469 + t522, t678 * t440 + t680 * t438 + pkin(5) * (-t678 * t471 + t680 * t472), t512, t488, t508, t513, t507, t527, t680 * t443 + t678 * t445 + t493, t680 * t444 + t678 * t446 + t494, t680 * t449 + t678 * t450 + t522, t678 * t436 + t680 * t435 + pkin(5) * (-t678 * t453 + t680 * t454); 0, 0, 0, 0, 0, 0, -g(2), g(1), 0, 0, 0, 0, 0, 0, 0, qJDD(1), t641, t642, 0, 0, 0, 0, 0, 0, 0, qJDD(1), t544, t543, 0, -t517, qJDD(1), 0, 0, 0, 0, 0, 0, t530, t529, t467, t592, t565, t593, t585, t590, 0, t473, t474, t466, t442, t592, t565, t593, t585, t590, 0, t459, t458, t452, t437; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, qJDD(1), 0, -t682, 0, 0, -g(3), -t641, 0, 0, 0, -t632, 0, -t631, 0, t600, -t727, t520, qJ(2) * t520, 0, t632, t631, 0, 0, 0, -t502, -t600, t727, t484, t559, t532, t556, t560, t555, t596, t455, t456, t469, t440, t559, t532, t556, t560, t555, t596, t445, t446, t450, t436; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t682, 0, qJDD(1), 0, g(3), 0, -t642, 0, 0, 0, t631, 0, -t632, 0, t727, t600, t724, t516, 0, -t631, t632, 0, 0, 0, t725, -t727, -t600, t483, t557, t531, t554, t558, t553, t595, t447, t448, t468, t438, t557, t531, t554, t558, t553, t595, t443, t444, t449, t435; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, qJDD(1), t641, t642, 0, 0, 0, 0, 0, 0, 0, qJDD(1), t544, t543, 0, -t517, qJDD(1), 0, 0, 0, 0, 0, 0, t530, t529, t467, t592, t565, t593, t585, t590, 0, t473, t474, t466, t442, t592, t565, t593, t585, t590, 0, t459, t458, t452, t437; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, qJDD(1), 0, -t682, 0, 0, -t670, t562, 0, 0, -qJDD(1), t682, 0, 0, 0, t552, 0, t670, qJ(3) * t670, t648, t637, t662, -t648, -t737, qJDD(4), t486, t487, -t623, t457, t648, t637, t662, -t648, -t737, qJDD(4), t475, t476, -t605, t441; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t682, 0, qJDD(1), 0, t670, 0, t563, 0, 0, -t682, -qJDD(1), 0, 0, 0, t550, -t670, 0, pkin(2) * t670, -t598, -t564, -t588, -t597, -t586, t618, t482, t481, -t477, t451, -t598, -t564, -t588, -t597, -t586, t618, t465, t464, -t460, t439; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, qJDD(1), -t562, -t563, 0, 0, qJDD(1), 0, 0, 0, 0, 0, 0, t701, t688, t747, t592, t565, t593, t585, t590, 0, t700, t706, t695, t715, t592, t565, t593, t585, t590, 0, t685, t696, t697, t686; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, qJDD(1), 0, 0, 0, 0, 0, 0, t552, t550, 0, t592, t565, t593, t585, t590, 0, t710, t728, t703, -t783, t592, t565, t593, t585, t590, 0, t690, t704, t705, t691; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, qJDD(1), -t682, 0, 0, 0, -t552, 0, -t670, 0, -t648, -t637, -t662, t648, t737, -qJDD(4), t713, t734, t623, -t789, -t648, -t637, -t662, t648, t737, -qJDD(4), t819, t815, t605, t729; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t682, qJDD(1), 0, 0, 0, -t550, t670, 0, 0, t598, t564, t588, t597, t586, -t618, t712, t711, t477, pkin(6) * t480 - t788, t598, t564, t588, t597, t586, -t618, t693, t698, t460, pkin(6) * t463 - t692; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t627, -t625, t640, t657, t644, -t657, 0, t542, t525, 0, t627, -t625, t640, t657, t644, -t657, -qJ(5) * t640, t505, t492, qJ(5) * t499; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t658, t628, t646, t626, t639, -t658, -t542, 0, -t744, 0, t658, t628, t646, t626, t639, -t658, t490, t584, t485, t470; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t648, t637, t662, -t648, -t737, qJDD(4), -t525, t744, 0, 0, t648, t637, t662, -t648, -t737, qJDD(4), t683, t699, -t735, -t497; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t627, -t625, t640, t657, t644, -t657, 0, t511, t499, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t658, t628, t646, t626, t639, -t658, -t511, 0, t500, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t648, t637, t662, -t648, -t737, qJDD(4), -t499, -t500, 0, 0;];
m_new_reg = t1;
