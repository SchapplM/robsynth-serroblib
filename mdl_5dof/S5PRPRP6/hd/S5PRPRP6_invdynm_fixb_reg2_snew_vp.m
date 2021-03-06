% Calculate inertial parameters regressor of inverse dynamics cutting torque vector with Newton-Euler for
% S5PRPRP6
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
%   pkin=[a2,a3,a4,a5,d2,d4,theta1]';
%
% Output:
% m_new_reg [(3*6)x(%Nl%*10)]
%   inertial parameter regressor of inverse dynamics cutting torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-05 15:41
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function m_new_reg = S5PRPRP6_invdynm_fixb_reg2_snew_vp(qJ, qJD, qJDD, g, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(5,1),zeros(3,1),zeros(7,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PRPRP6_invdynm_fixb_reg2_snew_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5PRPRP6_invdynm_fixb_reg2_snew_vp: qJD has to be [5x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [5 1]), ...
  'S5PRPRP6_invdynm_fixb_reg2_snew_vp: qJDD has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5PRPRP6_invdynm_fixb_reg2_snew_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [7 1]), ...
  'S5PRPRP6_invdynm_fixb_reg2_snew_vp: pkin has to be [7x1] (double)');

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_m_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 15:41:30
% EndTime: 2019-12-05 15:41:36
% DurationCPUTime: 5.77s
% Computational Cost: add. (7714->431), mult. (14819->422), div. (0->0), fcn. (7833->6), ass. (0->283)
t653 = qJD(4) ^ 2;
t651 = cos(qJ(4));
t640 = t651 ^ 2;
t654 = qJD(2) ^ 2;
t721 = t654 * t640;
t621 = t653 + t721;
t649 = sin(qJ(4));
t623 = t651 * t654 * t649;
t615 = qJDD(4) + t623;
t734 = t649 * t615;
t562 = t621 * t651 + t734;
t712 = qJD(2) * qJD(4);
t628 = t649 * t712;
t705 = t651 * qJDD(2);
t601 = -0.2e1 * t628 + t705;
t650 = sin(qJ(2));
t652 = cos(qJ(2));
t513 = t562 * t650 - t601 * t652;
t728 = t651 * t615;
t569 = -t621 * t649 + t728;
t644 = sin(pkin(7));
t645 = cos(pkin(7));
t807 = qJ(1) * (t513 * t645 + t644 * t569);
t806 = qJ(1) * (t513 * t644 - t645 * t569);
t509 = t562 * t652 + t601 * t650;
t805 = pkin(1) * t509;
t804 = pkin(5) * t509;
t803 = -pkin(1) * t569 + pkin(5) * t513;
t629 = t651 * t712;
t707 = t649 * qJDD(2);
t599 = 0.2e1 * t629 + t707;
t736 = t649 * t599;
t535 = -t601 * t651 + t736;
t639 = t649 ^ 2;
t612 = (-t639 + t640) * t654;
t502 = t535 * t650 - t612 * t652;
t729 = t651 * t599;
t735 = t649 * t601;
t537 = t729 + t735;
t802 = t502 * t644 - t645 * t537;
t801 = t502 * t645 + t644 * t537;
t722 = t654 * t639;
t619 = -t653 + t722;
t559 = t619 * t649 + t728;
t704 = t652 * qJDD(2);
t521 = t559 * t650 - t649 * t704;
t568 = -t619 * t651 + t734;
t798 = t521 * t644 + t645 * t568;
t797 = t521 * t645 - t644 * t568;
t706 = t650 * qJDD(2);
t609 = t652 * t654 + t706;
t613 = g(1) * t644 - g(2) * t645;
t699 = -pkin(5) * t609 + t613 * t652;
t796 = t644 * t699;
t795 = t645 * t699;
t793 = pkin(3) * t562;
t792 = pkin(6) * t562;
t791 = pkin(6) * t569;
t790 = qJ(3) * t569;
t610 = t650 * t654 - t704;
t550 = pkin(5) * t610 - t613 * t650;
t788 = t644 * t550;
t784 = t645 * t550;
t769 = pkin(6) + pkin(2);
t781 = t769 * t569;
t499 = t535 * t652 + t612 * t650;
t517 = t559 * t652 + t649 * t706;
t780 = qJ(3) * t601 + t562 * t769;
t600 = -t628 + t705;
t647 = t654 * pkin(6);
t614 = g(1) * t645 + g(2) * t644;
t641 = g(3) - qJDD(1);
t671 = t652 * t614 + t650 * t641;
t709 = qJDD(2) * qJ(3);
t661 = -pkin(2) * t654 - t671 + t709;
t686 = t629 + t707;
t658 = pkin(4) * t686 - qJ(5) * t600 - t647 + t661;
t759 = pkin(4) * t651;
t674 = qJ(5) * t649 + t759;
t713 = qJD(5) * t651;
t476 = (qJD(4) * t674 + (2 * qJD(3)) - 0.2e1 * t713) * qJD(2) + t658;
t779 = (pkin(3) + t674) * t476;
t643 = qJDD(2) * pkin(2);
t577 = -t614 * t650 + t641 * t652;
t673 = -qJDD(3) - t577;
t541 = -qJ(3) * t654 - t643 - t673;
t534 = -qJDD(2) * pkin(6) + t541;
t691 = t534 * t651 + t613 * t649;
t716 = -t534 * t649 + t613 * t651;
t462 = -t649 * t716 + t651 * t691;
t620 = -t653 - t722;
t616 = qJDD(4) - t623;
t727 = t651 * t616;
t560 = t620 * t649 + t727;
t778 = -pkin(2) * t560 + qJ(3) * t599;
t586 = t645 * t613;
t777 = -t614 * t644 + t586;
t715 = t639 + t640;
t608 = t715 * qJDD(2);
t611 = t715 * t654;
t776 = pkin(2) * t608 - qJ(3) * t611;
t775 = t674 * qJDD(2);
t711 = qJD(3) * qJD(2);
t703 = 0.2e1 * t711;
t539 = t661 + t703;
t530 = t539 - t647;
t774 = qJ(3) * t530 - t462 * t769;
t635 = -0.2e1 * t711;
t655 = 0.2e1 * qJD(2) * t713 + t635 - t658;
t470 = -pkin(4) * t629 + (t601 - t628) * qJ(5) + t655;
t718 = pkin(4) * t735 - t470 * t651;
t773 = -t718 - t780;
t516 = t651 * t530;
t772 = t516 + t780;
t771 = (pkin(3) + t759) * t601 + t649 * t470;
t710 = qJD(5) * qJD(4);
t633 = 0.2e1 * t710;
t748 = qJ(5) * t651;
t760 = pkin(4) * t649;
t714 = t654 * (-t748 + t760);
t672 = -pkin(4) * t653 + qJDD(4) * qJ(5) - t649 * t714 - t716;
t481 = t633 + t672;
t483 = -qJDD(4) * pkin(4) - qJ(5) * t653 + t651 * t714 + qJDD(5) - t691;
t451 = t481 * t649 - t483 * t651;
t475 = t476 * t760;
t770 = -t769 * t451 - (-qJ(3) + t748) * t476 + t475;
t508 = -t560 * t652 + t599 * t650;
t768 = pkin(1) * t508;
t725 = t652 * t608;
t542 = -t611 * t650 + t725;
t767 = pkin(1) * t542;
t766 = pkin(1) * t609;
t765 = pkin(1) * t610;
t733 = t649 * t616;
t566 = t620 * t651 - t733;
t764 = pkin(2) * t566;
t763 = pkin(3) * t462;
t762 = pkin(3) * t530;
t761 = pkin(3) * t611;
t758 = pkin(5) * t508;
t757 = pkin(5) * t542;
t511 = t560 * t650 + t599 * t652;
t755 = qJ(1) * (t511 * t644 - t566 * t645);
t731 = t650 * t608;
t543 = -t611 * t652 - t731;
t754 = qJ(1) * t543;
t753 = qJ(1) * t609;
t752 = qJ(1) * t610;
t750 = qJ(3) * t566;
t747 = t644 * t609;
t746 = t644 * t610;
t745 = t644 * t613;
t743 = t644 * t641;
t742 = t645 * t609;
t741 = t645 * t610;
t740 = t645 * t641;
t737 = t649 * t530;
t720 = -pkin(4) * t483 + qJ(5) * t481;
t719 = -pkin(1) * t566 + pkin(5) * t511;
t717 = pkin(3) * t599 - pkin(6) * t566;
t708 = t645 * qJDD(2);
t663 = -t671 + 0.2e1 * t709;
t514 = t635 - t663 - t766;
t697 = -t514 + t752;
t668 = -0.2e1 * t643 - t673;
t515 = -t668 - t765;
t696 = t515 + t753;
t532 = -t671 + t766;
t695 = t532 + t752;
t533 = t577 + t765;
t694 = -t533 + t753;
t690 = t539 * t652 + t541 * t650;
t544 = pkin(6) * t560;
t689 = -t544 + t737;
t688 = t577 * t650 - t652 * t671;
t687 = -t614 * t645 - t745;
t685 = t650 * t623;
t684 = t652 * t623;
t683 = -pkin(3) * t451 - t720;
t682 = -t716 + t793;
t473 = pkin(4) * t611 + t481;
t474 = qJ(5) * t611 + t483;
t596 = pkin(6) * t608;
t680 = -t473 * t649 + t474 * t651 + t596;
t679 = t596 - t462;
t546 = pkin(3) * t560;
t678 = t546 + t691;
t677 = -t717 - t516;
t676 = -pkin(2) * t541 + qJ(3) * t539;
t675 = pkin(3) * t601 - t737;
t463 = -t649 * t691 - t651 * t716;
t491 = t539 * t650 - t541 * t652;
t505 = t577 * t652 + t650 * t671;
t669 = t689 + t778;
t667 = t680 + t776;
t666 = t679 + t776;
t665 = pkin(4) * t621 + qJ(5) * t615 + t672;
t471 = -qJ(5) * t628 + (-t599 - t629) * pkin(4) + t655;
t664 = -qJ(5) * t729 - t471 * t649 - t544;
t662 = -qJ(5) * t736 + t471 * t651 - t717;
t660 = t664 + t778;
t659 = -t665 - t793;
t657 = pkin(4) * t616 + qJ(5) * t620 - t483;
t656 = -t546 - t657;
t630 = t644 * qJDD(2);
t622 = -t653 + t721;
t603 = pkin(1) * t613;
t597 = pkin(3) * t608;
t593 = t715 * t712;
t576 = qJDD(4) * t652 - t593 * t650;
t575 = qJDD(4) * t650 + t593 * t652;
t574 = t600 * t649 + t640 * t712;
t573 = t639 * t712 - t651 * t686;
t570 = t622 * t649 + t727;
t567 = (t600 - t628) * t651;
t563 = -t622 * t651 + t733;
t558 = (t629 + t686) * t649;
t557 = t597 + t775;
t556 = t645 * t576;
t555 = t644 * t576;
t540 = pkin(5) * t543;
t531 = t645 * t754;
t526 = t573 * t650 - t684;
t525 = t574 * t650 + t684;
t524 = -t573 * t652 - t685;
t523 = -t574 * t652 + t685;
t522 = t563 * t650 + t651 * t704;
t519 = -t563 * t652 + t650 * t705;
t494 = pkin(5) * t688 + t603;
t489 = t526 * t645 + t558 * t644;
t488 = t525 * t645 + t567 * t644;
t487 = t526 * t644 - t558 * t645;
t486 = t525 * t644 - t567 * t645;
t485 = t522 * t645 + t570 * t644;
t484 = t522 * t644 - t570 * t645;
t477 = qJ(1) * (t511 * t645 + t566 * t644);
t469 = -pkin(5) * t491 + (-pkin(2) * t650 + qJ(3) * t652) * t613;
t467 = -t682 + t790;
t466 = t678 - t750;
t465 = -t677 - t764;
t464 = t675 + t781;
t461 = pkin(5) * t690 + t603 + (pkin(2) * t652 + qJ(3) * t650) * t613;
t460 = t463 + t761;
t459 = -pkin(1) * t491 - t676;
t458 = t462 * t650 + t530 * t652;
t457 = -t462 * t652 + t530 * t650;
t456 = -t656 - t750;
t455 = -t772 - t805;
t454 = -t669 - t768;
t453 = t633 - t659 - t790;
t452 = t481 * t651 + t483 * t649;
t450 = -pkin(3) * t725 + t460 * t650 - t757;
t449 = -pkin(3) * t731 - t460 * t652 + t540;
t448 = -t666 - t767;
t447 = -t771 - t781;
t446 = -t662 - t764;
t445 = t473 * t651 + t474 * t649 + t761;
t444 = -qJ(3) * t463 + t763;
t443 = -t773 + t805;
t442 = -t660 - t768;
t441 = t451 * t650 + t476 * t652;
t440 = -t451 * t652 + t476 * t650;
t439 = -t667 - t767;
t438 = -t463 * t769 + t762;
t437 = t445 * t650 - t557 * t652 - t757;
t436 = -t445 * t652 - t557 * t650 + t540;
t435 = -t464 * t650 + t467 * t652 - t804;
t434 = -t465 * t650 + t466 * t652 - t758;
t433 = t464 * t652 + t467 * t650 - t803;
t432 = t465 * t652 + t466 * t650 + t719;
t431 = -t446 * t650 + t456 * t652 - t758;
t430 = -t447 * t650 + t453 * t652 + t804;
t429 = -pkin(1) * t457 - t774;
t428 = t446 * t652 + t456 * t650 + t719;
t427 = t447 * t652 + t453 * t650 + t803;
t426 = -qJ(3) * t452 - t683;
t425 = -t452 * t769 + t779;
t424 = -pkin(5) * t457 - t438 * t650 + t444 * t652;
t423 = -pkin(1) * t463 + pkin(5) * t458 + t438 * t652 + t444 * t650;
t422 = -pkin(1) * t440 - t770;
t421 = -pkin(5) * t440 - t425 * t650 + t426 * t652;
t420 = -pkin(1) * t452 + pkin(5) * t441 + t425 * t652 + t426 * t650;
t1 = [0, 0, 0, 0, 0, 0, 0, -g(3), g(2), 0, 0, 0, 0, 0, 0, 0, -t743, -t740, -t777, -qJ(1) * t777, 0, 0, -t741, 0, -t742, t630, t644 * t694 + t784, -t644 * t695 - t795, t645 * t505, -qJ(1) * (t644 * t688 + t586) - (pkin(1) * t644 - pkin(5) * t645) * t505, t630, t741, t742, 0, 0, 0, -t645 * t491, -t644 * t696 - t784, t644 * t697 + t795, t645 * t469 - t644 * t459 - qJ(1) * (t644 * t690 + t586), t488, -t801, t485, t489, t797, t556, t434 * t645 - t454 * t644 - t755, t645 * t435 - t644 * t455 + t806, t645 * t450 + (-t448 - t754) * t644, t645 * t424 - t644 * t429 - qJ(1) * (t458 * t644 - t463 * t645), t488, t485, t801, t556, -t797, t489, t431 * t645 - t442 * t644 - t755, t645 * t437 + (-t439 - t754) * t644, t645 * t430 - t644 * t443 - t806, t645 * t421 - t644 * t422 - qJ(1) * (t441 * t644 - t452 * t645); 0, 0, 0, 0, 0, 0, g(3), 0, -g(1), 0, 0, 0, 0, 0, 0, 0, t740, -t743, t687, qJ(1) * t687, 0, 0, -t746, 0, -t747, -t708, -t645 * t694 + t788, t645 * t695 - t796, t644 * t505, qJ(1) * (t645 * t688 - t745) - (-pkin(1) * t645 - pkin(5) * t644) * t505, -t708, t746, t747, 0, 0, 0, -t644 * t491, t645 * t696 - t788, -t645 * t697 + t796, t644 * t469 + t645 * t459 + qJ(1) * (t645 * t690 - t745), t486, -t802, t484, t487, t798, t555, t434 * t644 + t454 * t645 + t477, t644 * t435 + t645 * t455 - t807, t448 * t645 + t450 * t644 + t531, t644 * t424 + t645 * t429 + qJ(1) * (t458 * t645 + t463 * t644), t486, t484, t802, t555, -t798, t487, t431 * t644 + t442 * t645 + t477, t437 * t644 + t439 * t645 + t531, t644 * t430 + t645 * t443 + t807, t644 * t421 + t645 * t422 + qJ(1) * (t441 * t645 + t452 * t644); 0, 0, 0, 0, 0, 0, -g(2), g(1), 0, 0, 0, 0, 0, 0, 0, 0, t613, t614, 0, 0, 0, 0, t609, 0, -t610, 0, t699, t550, t688, t494, 0, -t609, t610, 0, 0, 0, t690, -t699, -t550, t461, t523, t499, t519, t524, -t517, t575, t432, t433, t449, t423, t523, t519, -t499, t575, t517, t524, t428, t436, t427, t420; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t641, -t613, 0, 0, 0, -t610, 0, -t609, 0, t550, -t699, t505, pkin(5) * t505, 0, t610, t609, 0, 0, 0, -t491, -t550, t699, t469, t525, -t502, t522, t526, t521, t576, t434, t435, t450, t424, t525, t522, t502, t576, -t521, t526, t431, t437, t430, t421; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t641, 0, -t614, 0, 0, 0, 0, 0, 0, -qJDD(2), t533, t532, 0, pkin(1) * t505, -qJDD(2), 0, 0, 0, 0, 0, 0, t515, t514, t459, -t567, t537, -t570, -t558, t568, 0, t454, t455, t448, t429, -t567, -t570, -t537, 0, -t568, -t558, t442, t439, t443, t422; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t613, t614, 0, 0, 0, 0, t609, 0, -t610, 0, t699, t550, t688, t494, 0, -t609, t610, 0, 0, 0, t690, -t699, -t550, t461, t523, t499, t519, t524, -t517, t575, t432, t433, t449, t423, t523, t519, -t499, t575, t517, t524, t428, t436, t427, t420; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, qJDD(2), 0, -t654, 0, 0, -t613, t577, 0, 0, -qJDD(2), t654, 0, 0, 0, t541, 0, t613, qJ(3) * t613, t623, t612, t705, -t623, -t707, qJDD(4), t466, t467, -t597, t444, t623, t705, -t612, qJDD(4), t707, -t623, t456, -t557, t453, t426; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t654, 0, qJDD(2), 0, t613, 0, -t671, 0, 0, -t654, -qJDD(2), 0, 0, 0, t539, -t613, 0, pkin(2) * t613, -t574, t535, -t563, -t573, -t559, t593, t465, t464, -t460, t438, -t574, -t563, -t535, t593, t559, -t573, t446, -t445, t447, t425; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, qJDD(2), -t577, t671, 0, 0, qJDD(2), 0, 0, 0, 0, 0, 0, t668, t663 + t703, t676, t567, -t537, t570, t558, -t568, 0, t669, t772, t666, t774, t567, t570, t537, 0, t568, t558, t660, t667, t773, t770; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, qJDD(2), 0, 0, 0, 0, 0, 0, t541, t539, 0, t567, -t537, t570, t558, -t568, 0, t689, t516 + t792, t679, -pkin(6) * t462, t567, t570, t537, 0, t568, t558, t664, t680, -t718 - t792, -pkin(6) * t451 - t476 * t748 + t475; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, qJDD(2), -t654, 0, 0, 0, -t541, 0, -t613, 0, -t623, -t612, -t705, t623, t707, -qJDD(4), -t678, t682, t597, -t763, -t623, -t705, t612, -qJDD(4), -t707, t623, t656, t557, t659 - 0.2e1 * t710, t683; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t654, qJDD(2), 0, 0, 0, -t539, t613, 0, 0, t574, -t535, t563, t573, t559, -t593, t677, -t675 - t791, t460, pkin(6) * t463 - t762, t574, t563, t535, -t593, -t559, t573, t662, t445, t771 + t791, pkin(6) * t452 - t779; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t600, -t599, t616, t628, t619, -t628, 0, t530, -t691, 0, t600, t616, t599, -t628, -t619, t628, -qJ(5) * t599, t474, t470, -qJ(5) * t476; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t629, t601, -t622, -t686, t615, -t629, -t530, 0, -t716, 0, t629, -t622, -t601, -t629, -t615, -t686, t471, t473, pkin(4) * t601, -pkin(4) * t476; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t623, t612, t705, -t623, -t707, qJDD(4), t691, t716, 0, 0, t623, t705, -t612, qJDD(4), t707, -t623, t657, -t775, t633 + t665, t720; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t600, t616, t599, -t628, -t619, t628, 0, t483, -t476, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t623, t705, -t612, qJDD(4), t707, -t623, -t483, 0, t481, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t629, t622, t601, t629, t615, t686, t476, -t481, 0, 0;];
m_new_reg = t1;
