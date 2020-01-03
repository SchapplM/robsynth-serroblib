% Calculate inertial parameters regressor of inverse dynamics base forces vector with Newton-Euler for
% S5RRRPP7
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
%   pkin=[a2,a3,a4,a5,d1,d2,d3]';
% 
% Output:
% tauB_reg [6x(6*10)]
%   inertial parameter regressor of inverse dynamics base forces vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 21:06
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauB_reg = S5RRRPP7_invdynB_fixb_reg2_snew_vp(qJ, qJD, qJDD, g, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(5,1),zeros(3,1),zeros(7,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRRPP7_invdynB_fixb_reg2_snew_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RRRPP7_invdynB_fixb_reg2_snew_vp: qJD has to be [5x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [5 1]), ...
  'S5RRRPP7_invdynB_fixb_reg2_snew_vp: qJDD has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RRRPP7_invdynB_fixb_reg2_snew_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [7 1]), ...
  'S5RRRPP7_invdynB_fixb_reg2_snew_vp: pkin has to be [7x1] (double)');

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_tauB_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 21:06:04
% EndTime: 2019-12-31 21:06:20
% DurationCPUTime: 13.07s
% Computational Cost: add. (15254->492), mult. (30609->615), div. (0->0), fcn. (19560->6), ass. (0->358)
t665 = sin(qJ(3));
t668 = cos(qJ(3));
t666 = sin(qJ(2));
t711 = qJD(1) * t666;
t625 = t665 * qJD(2) + t668 * t711;
t669 = cos(qJ(2));
t708 = t669 * qJD(1);
t652 = -qJD(3) + t708;
t705 = qJD(1) * qJD(2);
t698 = t669 * t705;
t704 = t666 * qJDD(1);
t629 = t698 + t704;
t695 = -t668 * qJDD(2) + t665 * t629;
t679 = (qJD(3) + t652) * t625 + t695;
t623 = -t668 * qJD(2) + t665 * t711;
t683 = -t665 * qJDD(2) - t668 * t629;
t678 = t623 * qJD(3) + t683;
t751 = t652 * t623;
t774 = t751 + t678;
t452 = t665 * t774 + t668 * t679;
t619 = t625 ^ 2;
t761 = t623 ^ 2;
t547 = -t761 - t619;
t418 = t669 * t452 - t666 * t547;
t446 = t665 * t679 - t668 * t774;
t667 = sin(qJ(1));
t670 = cos(qJ(1));
t381 = t667 * t418 - t670 * t446;
t891 = pkin(5) * t381;
t384 = t670 * t418 + t667 * t446;
t890 = pkin(5) * t384;
t760 = t652 ^ 2;
t772 = t761 - t760;
t577 = t625 * t623;
t655 = t666 * t705;
t702 = t669 * qJDD(1);
t630 = -t655 + t702;
t620 = -qJDD(3) + t630;
t780 = -t577 + t620;
t801 = t665 * t780;
t808 = t668 * t772 + t801;
t836 = t669 * t808;
t443 = t666 * t679 - t836;
t796 = t668 * t780;
t811 = -t665 * t772 + t796;
t834 = t670 * t811;
t889 = t667 * t443 - t834;
t839 = t667 * t811;
t888 = t670 * t443 + t839;
t415 = t666 * t452 + t669 * t547;
t887 = pkin(6) * t415;
t886 = pkin(6) * t418;
t885 = pkin(1) * t415 + pkin(2) * t547 + pkin(7) * t452;
t884 = -pkin(1) * t446 + t886;
t595 = t619 - t760;
t779 = t577 + t620;
t797 = t668 * t779;
t486 = -t665 * t595 + t797;
t438 = t669 * t486 + t666 * t774;
t802 = t665 * t779;
t480 = -t668 * t595 - t802;
t883 = t667 * t438 + t670 * t480;
t882 = t670 * t438 - t667 * t480;
t777 = t619 + t760;
t812 = t668 * t777 - t801;
t838 = t667 * t812;
t775 = t751 - t678;
t477 = t665 * t777 + t796;
t852 = t669 * t477;
t846 = -t666 * t775 - t852;
t856 = t670 * t846 + t838;
t881 = pkin(5) * t856;
t773 = -t760 - t761;
t809 = t665 * t773 - t797;
t819 = t667 * t809;
t560 = t625 * qJD(3) + t695;
t601 = t625 * t652;
t514 = t560 - t601;
t807 = t668 * t773 + t802;
t816 = t669 * t807;
t850 = t666 * t514 + t816;
t857 = t670 * t850 + t819;
t880 = pkin(5) * t857;
t833 = t670 * t812;
t858 = t667 * t846 - t833;
t879 = pkin(5) * t858;
t815 = t670 * t809;
t859 = t667 * t850 - t815;
t878 = pkin(5) * t859;
t841 = t666 * t808;
t877 = t669 * t679 + t841;
t707 = qJD(3) - t652;
t516 = t707 * t625 + t695;
t746 = t665 * t775;
t453 = t668 * t516 + t746;
t573 = -t761 + t619;
t734 = t666 * t573;
t423 = t669 * t453 - t734;
t818 = t668 * t775;
t447 = -t665 * t516 + t818;
t876 = t667 * t423 + t670 * t447;
t875 = t670 * t423 - t667 * t447;
t853 = t666 * t477;
t847 = t669 * t775 - t853;
t873 = pkin(6) * t847;
t820 = t666 * t807;
t851 = -t669 * t514 + t820;
t872 = pkin(6) * t851;
t871 = pkin(7) * t446;
t825 = pkin(7) * t807;
t866 = -pkin(1) * t851 - t825;
t854 = pkin(7) * t477;
t865 = -pkin(1) * t847 + t854;
t864 = pkin(2) * t446 + qJ(4) * t679;
t828 = pkin(1) * t809;
t863 = pkin(6) * t850 - t828;
t843 = pkin(1) * t812;
t862 = pkin(6) * t846 - t843;
t860 = t666 * t486 - t669 * t774;
t855 = pkin(2) * t812;
t842 = pkin(7) * t812;
t826 = pkin(2) * t809;
t849 = -qJ(4) * t773 - t826;
t564 = t669 * t573;
t848 = t666 * t453 + t564;
t824 = pkin(7) * t809;
t831 = qJ(4) * t780 - t855;
t823 = qJ(4) * t547;
t822 = qJ(4) * t775;
t814 = t774 * qJ(5);
t642 = t670 * g(1) + t667 * g(2);
t671 = qJD(1) ^ 2;
t611 = -t671 * pkin(1) + qJDD(1) * pkin(6) - t642;
t755 = pkin(2) * t669;
t690 = -pkin(7) * t666 - t755;
t627 = t690 * qJD(1);
t752 = t669 * g(3);
t758 = qJD(2) ^ 2;
t538 = -qJDD(2) * pkin(2) - t758 * pkin(7) + (qJD(1) * t627 + t611) * t666 + t752;
t794 = t560 * pkin(3) + t538 - t822;
t749 = t652 * t668;
t590 = t625 * t749;
t750 = t652 * t665;
t700 = t623 * t750;
t687 = -t590 - t700;
t589 = t625 * t750;
t699 = t623 * t749;
t688 = -t589 + t699;
t732 = t666 * t620;
t763 = t669 * t688 - t732;
t792 = t667 * t763 + t670 * t687;
t681 = t668 * t560 + t700;
t567 = t666 * t577;
t682 = t665 * t560 - t699;
t767 = t669 * t682 - t567;
t791 = t667 * t767 + t670 * t681;
t789 = -t667 * t687 + t670 * t763;
t788 = -t667 * t681 + t670 * t767;
t713 = -t668 * t678 + t589;
t766 = t669 * t713 + t567;
t776 = t665 * t678 + t590;
t762 = -t667 * t776 + t670 * t766;
t768 = t667 * t766 + t670 * t776;
t710 = qJD(4) * t652;
t640 = 0.2e1 * t710;
t709 = qJD(5) * t623;
t778 = -0.2e1 * t709 + t640;
t586 = t652 * pkin(4) - t625 * qJ(5);
t771 = t625 * t586 + qJDD(5);
t769 = -t560 * pkin(4) + t771;
t568 = t669 * t577;
t765 = t666 * t682 + t568;
t715 = t666 * t713 - t568;
t606 = t669 * t620;
t764 = t666 * t688 + t606;
t759 = 0.2e1 * t625;
t757 = pkin(3) + pkin(4);
t756 = pkin(2) * t666;
t754 = pkin(3) * t668;
t661 = t666 ^ 2;
t748 = t661 * t671;
t744 = t665 * t538;
t641 = t667 * g(1) - t670 * g(2);
t610 = qJDD(1) * pkin(1) + t671 * pkin(6) + t641;
t733 = t666 * t610;
t651 = t669 * t671 * t666;
t637 = -t651 + qJDD(2);
t731 = t666 * t637;
t638 = qJDD(2) + t651;
t730 = t666 * t638;
t728 = t668 * t514;
t725 = t668 * t538;
t717 = t669 * t610;
t716 = t669 * t637;
t685 = -t630 + t655;
t686 = t629 + t698;
t512 = pkin(2) * t685 - pkin(7) * t686 - t610;
t588 = -t666 * g(3) + t669 * t611;
t539 = -t758 * pkin(2) + qJDD(2) * pkin(7) + t627 * t708 + t588;
t458 = t665 * t512 + t668 * t539;
t662 = t669 ^ 2;
t712 = t661 + t662;
t703 = t667 * qJDD(1);
t701 = t670 * qJDD(1);
t697 = qJ(4) * t665 + pkin(2);
t566 = t623 * pkin(3) - t625 * qJ(4);
t696 = -pkin(4) * t623 - t566;
t457 = -t668 * t512 + t665 * t539;
t587 = t666 * t611 + t752;
t528 = t666 * t587 + t669 * t588;
t579 = -t667 * t641 - t670 * t642;
t693 = t667 * t651;
t692 = t670 * t651;
t634 = -t667 * t671 + t701;
t689 = -pkin(5) * t634 - t667 * g(3);
t684 = pkin(3) * t760 + t620 * qJ(4) + t623 * t566 - t458;
t401 = -t668 * t457 + t665 * t458;
t402 = t665 * t457 + t668 * t458;
t527 = t669 * t587 - t666 * t588;
t578 = t670 * t641 - t667 * t642;
t639 = -0.2e1 * t710;
t412 = t639 - t684;
t680 = t620 * pkin(3) - qJ(4) * t760 + qJDD(4) + t457;
t677 = pkin(4) * t761 + t652 * t586 + t684;
t676 = t625 * t566 + t680;
t674 = t620 * pkin(4) + t680 + t814;
t673 = -t560 * qJ(5) + t677;
t672 = qJD(4) * t759 - t794;
t380 = (-0.2e1 * qJD(5) - t696) * t625 + t674;
t413 = (-pkin(3) * t652 - 0.2e1 * qJD(4)) * t625 + t794;
t393 = (-t514 + t601) * pkin(3) + t672;
t392 = pkin(3) * t601 + t672 + t822;
t659 = t662 * t671;
t648 = -t659 - t758;
t647 = t659 - t758;
t646 = -t748 - t758;
t645 = -t748 + t758;
t636 = t659 - t748;
t635 = t659 + t748;
t633 = t670 * t671 + t703;
t632 = t712 * qJDD(1);
t631 = -0.2e1 * t655 + t702;
t628 = 0.2e1 * t698 + t704;
t622 = t669 * t638;
t621 = t712 * t705;
t605 = -pkin(5) * t633 + t670 * g(3);
t592 = t669 * t629 - t661 * t705;
t591 = -t666 * t630 - t662 * t705;
t585 = -t666 * t646 - t716;
t584 = -t666 * t645 + t622;
t583 = t669 * t648 - t730;
t582 = t669 * t647 - t731;
t581 = t669 * t646 - t731;
t580 = t666 * t648 + t622;
t572 = t670 * t632 - t667 * t635;
t571 = t667 * t632 + t670 * t635;
t565 = -t666 * t628 + t669 * t631;
t545 = t670 * t585 + t667 * t628;
t544 = t670 * t583 - t667 * t631;
t543 = t667 * t585 - t670 * t628;
t542 = t667 * t583 + t670 * t631;
t541 = -pkin(6) * t581 - t717;
t540 = -pkin(6) * t580 - t733;
t536 = (t623 * t668 - t625 * t665) * t652;
t533 = (-t623 * t665 - t625 * t668) * t652;
t530 = -pkin(1) * t581 + t588;
t529 = -pkin(1) * t580 + t587;
t523 = t707 * t623 + t683;
t513 = t560 + t601;
t493 = t669 * t536 - t732;
t473 = t670 * t528 - t667 * t610;
t472 = t667 * t528 + t670 * t610;
t459 = -qJ(4) * t514 - qJ(5) * t779;
t451 = t728 + t746;
t445 = -t665 * t514 + t818;
t444 = t725 + t842;
t441 = -t666 * t513 + t836;
t437 = t744 - t824;
t434 = -t666 * t523 + t852;
t431 = t669 * t523 + t853;
t430 = t666 * t516 + t816;
t427 = -t669 * t516 + t820;
t422 = t669 * t451 - t734;
t421 = qJ(5) * t780 + t757 * t775;
t411 = t458 + t855;
t410 = t457 - t826;
t409 = t676 - t823;
t406 = t670 * t434 - t838;
t403 = t667 * t434 + t833;
t400 = -pkin(3) * t547 + t412;
t399 = t670 * t430 + t819;
t396 = t667 * t430 - t815;
t391 = -pkin(3) * t774 + t864;
t390 = qJ(5) * t761 + t413 - t769;
t389 = t669 * t402 + t666 * t538;
t388 = t666 * t402 - t669 * t538;
t387 = t639 - t673 + 0.2e1 * t709;
t379 = t757 * t774 - t864;
t378 = -pkin(1) * t431 - pkin(2) * t523 - t744 - t854;
t377 = pkin(3) * t779 + t676 + t849;
t376 = -pkin(1) * t427 + pkin(2) * t516 + t725 - t825;
t375 = -pkin(3) * t777 + t640 + t684 + t831;
t374 = -t401 + t871;
t373 = t392 + (t777 - t761) * qJ(5) + t769;
t372 = t668 * t412 + t665 * t676;
t371 = t665 * t412 - t668 * t676;
t370 = -qJ(4) * t728 - t665 * t393 - t824;
t369 = -pkin(3) * t746 + t668 * t392 - t842;
t368 = qJD(5) * t759 + t625 * t696 - t674 - t814 + t823;
t367 = (-t773 - t761) * qJ(5) + (-t514 - t560) * pkin(4) + t393 + t771;
t366 = t757 * t547 + (-t679 - t560) * qJ(5) + t677 + t778;
t365 = -pkin(6) * t431 - t666 * t411 + t669 * t444;
t364 = -t757 * t777 + t673 + t778 + t831;
t363 = -pkin(6) * t427 - t666 * t410 + t669 * t437;
t362 = t757 * t779 + t380 + t849;
t361 = t670 * t389 + t667 * t401;
t360 = t667 * t389 - t670 * t401;
t359 = t669 * t372 + t666 * t413;
t358 = t666 * t372 - t669 * t413;
t357 = -t402 + t885;
t356 = -pkin(1) * t388 + pkin(2) * t538 - pkin(7) * t402;
t355 = -qJ(4) * t390 - qJ(5) * t380;
t354 = -t668 * t393 + t514 * t697 + t866;
t353 = -t665 * t400 + t668 * t409 + t871;
t352 = t665 * t380 + t668 * t387;
t351 = -t668 * t380 + t665 * t387;
t350 = -t665 * t392 + (-pkin(2) - t754) * t775 + t865;
t349 = t668 * t373 - t665 * t421 - t842;
t348 = -t665 * t367 + t668 * t459 - t824;
t347 = t669 * t374 - t446 * t756 + t887;
t346 = -pkin(7) * t371 + (pkin(3) * t665 - qJ(4) * t668) * t413;
t345 = -pkin(2) * t371 + pkin(3) * t676 - qJ(4) * t412;
t344 = -pkin(6) * t388 + (-pkin(7) * t669 + t756) * t401;
t343 = -t668 * t400 - t665 * t409 + t885;
t342 = -pkin(2) * t775 - t665 * t373 - t668 * t421 + t865;
t341 = pkin(2) * t514 - t668 * t367 - t665 * t459 + t866;
t340 = -qJ(5) * t387 - t757 * t390;
t339 = t669 * t370 - t666 * t377 - t872;
t338 = t669 * t352 + t666 * t390;
t337 = t666 * t352 - t669 * t390;
t336 = t669 * t369 - t666 * t375 - t873;
t335 = t670 * t359 + t667 * t371;
t334 = t667 * t359 - t670 * t371;
t333 = -t665 * t366 + t668 * t368 - t871;
t332 = t669 * t353 - t666 * t391 + t887;
t331 = -t668 * t366 - t665 * t368 - t885;
t330 = t669 * t349 - t666 * t364 - t873;
t329 = t669 * t348 - t666 * t362 - t872;
t328 = t669 * t333 - t666 * t379 - t887;
t327 = t670 * t338 + t667 * t351;
t326 = t667 * t338 - t670 * t351;
t325 = -pkin(2) * t351 - qJ(4) * t387 + t757 * t380;
t324 = -pkin(1) * t358 - pkin(7) * t372 + (t697 + t754) * t413;
t323 = -pkin(7) * t351 - t665 * t340 + t668 * t355;
t322 = -pkin(6) * t358 - t666 * t345 + t669 * t346;
t321 = -pkin(1) * t337 + pkin(2) * t390 - pkin(7) * t352 - t668 * t340 - t665 * t355;
t320 = -pkin(6) * t337 + t669 * t323 - t666 * t325;
t1 = [0, 0, 0, 0, 0, 0, 0, 0, 0, -g(1), 0, 0, 0, 0, 0, 0, -t633, -t634, 0, t579, 0, 0, 0, 0, 0, 0, t544, t545, t572, t473, 0, 0, 0, 0, 0, 0, t399, t406, -t384, t361, 0, 0, 0, 0, 0, 0, t857, -t384, t856, t335, 0, 0, 0, 0, 0, 0, t857, t856, t384, t327; 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(2), 0, 0, 0, 0, 0, 0, t634, -t633, 0, t578, 0, 0, 0, 0, 0, 0, t542, t543, t571, t472, 0, 0, 0, 0, 0, 0, t396, t403, -t381, t360, 0, 0, 0, 0, 0, 0, t859, -t381, t858, t334, 0, 0, 0, 0, 0, 0, t859, t858, t381, t326; 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(3), 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(3), 0, 0, 0, 0, 0, 0, t580, t581, 0, -t527, 0, 0, 0, 0, 0, 0, t427, t431, -t415, t388, 0, 0, 0, 0, 0, 0, t851, -t415, t847, t358, 0, 0, 0, 0, 0, 0, t851, t847, t415, t337; 0, 0, 0, 0, 0, 0, 0, -g(3), g(2), 0, 0, 0, t634, 0, -t633, 0, t689, -t605, -t578, -pkin(5) * t578, t670 * t592 - t693, t670 * t565 - t667 * t636, t670 * t584 + t666 * t703, t670 * t591 + t693, t670 * t582 + t667 * t702, t667 * qJDD(2) + t670 * t621, -pkin(5) * t542 - t667 * t529 + t670 * t540, -pkin(5) * t543 - t667 * t530 + t670 * t541, -pkin(5) * t571 + t670 * t527, -pkin(5) * t472 - (pkin(1) * t667 - pkin(6) * t670) * t527, t762, -t875, -t882, t788, -t888, t789, -pkin(5) * t396 + t670 * t363 - t667 * t376, -pkin(5) * t403 + t670 * t365 - t667 * t378, t670 * t347 - t667 * t357 + t891, -pkin(5) * t360 + t670 * t344 - t667 * t356, t762, -t882, t875, t789, t888, t788, t670 * t339 - t667 * t354 - t878, t670 * t332 - t667 * t343 + t891, t670 * t336 - t667 * t350 - t879, -pkin(5) * t334 + t670 * t322 - t667 * t324, t762, t670 * t422 - t667 * t445, t882, t788, t670 * t441 - t839, t670 * t493 - t667 * t533, t670 * t329 - t667 * t341 - t878, t670 * t330 - t667 * t342 - t879, t670 * t328 - t667 * t331 - t891, -pkin(5) * t326 + t670 * t320 - t667 * t321; 0, 0, 0, 0, 0, 0, g(3), 0, -g(1), 0, 0, 0, t633, 0, t634, 0, t605, t689, t579, pkin(5) * t579, t667 * t592 + t692, t667 * t565 + t670 * t636, t667 * t584 - t666 * t701, t667 * t591 - t692, t667 * t582 - t669 * t701, -t670 * qJDD(2) + t667 * t621, pkin(5) * t544 + t670 * t529 + t667 * t540, pkin(5) * t545 + t670 * t530 + t667 * t541, pkin(5) * t572 + t667 * t527, pkin(5) * t473 - (-pkin(1) * t670 - pkin(6) * t667) * t527, t768, -t876, -t883, t791, -t889, t792, pkin(5) * t399 + t667 * t363 + t670 * t376, pkin(5) * t406 + t667 * t365 + t670 * t378, t667 * t347 + t670 * t357 - t890, pkin(5) * t361 + t667 * t344 + t670 * t356, t768, -t883, t876, t792, t889, t791, t667 * t339 + t670 * t354 + t880, t667 * t332 + t670 * t343 - t890, t667 * t336 + t670 * t350 + t881, pkin(5) * t335 + t667 * t322 + t670 * t324, t768, t667 * t422 + t670 * t445, t883, t791, t667 * t441 + t834, t667 * t493 + t670 * t533, t667 * t329 + t670 * t341 + t880, t667 * t330 + t670 * t342 + t881, t667 * t328 + t670 * t331 + t890, pkin(5) * t327 + t667 * t320 + t670 * t321; 0, 0, 0, 0, 0, 0, -g(2), g(1), 0, 0, 0, 0, 0, 0, 0, qJDD(1), t641, t642, 0, 0, t686 * t666, t669 * t628 + t666 * t631, t669 * t645 + t730, -t685 * t669, t666 * t647 + t716, 0, pkin(1) * t631 + pkin(6) * t583 + t717, -pkin(1) * t628 + pkin(6) * t585 - t733, pkin(1) * t635 + pkin(6) * t632 + t528, pkin(1) * t610 + pkin(6) * t528, t715, -t848, -t860, t765, t877, t764, pkin(6) * t430 + t669 * t410 + t666 * t437 - t828, pkin(6) * t434 + t669 * t411 + t666 * t444 + t843, -t886 + t666 * t374 - (-pkin(1) - t755) * t446, pkin(6) * t389 + (-pkin(1) + t690) * t401, t715, -t860, t848, t764, -t877, t765, t666 * t370 + t669 * t377 + t863, t666 * t353 + t669 * t391 - t884, t666 * t369 + t669 * t375 + t862, -pkin(1) * t371 + pkin(6) * t359 + t669 * t345 + t666 * t346, t715, t666 * t451 + t564, t860, t765, t669 * t513 + t841, t666 * t536 + t606, t666 * t348 + t669 * t362 + t863, t666 * t349 + t669 * t364 + t862, t666 * t333 + t669 * t379 + t884, -pkin(1) * t351 + pkin(6) * t338 + t666 * t323 + t669 * t325;];
tauB_reg = t1;
