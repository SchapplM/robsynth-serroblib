% Calculate inertial parameters regressor of inverse dynamics base forces vector with Newton-Euler for
% S5RRRPP1
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
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d2,d3,theta4]';
% 
% Output:
% tauB_reg [6x(6*10)]
%   inertial parameter regressor of inverse dynamics base forces vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 20:50
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauB_reg = S5RRRPP1_invdynB_fixb_reg2_snew_vp(qJ, qJD, qJDD, g, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(5,1),zeros(3,1),zeros(8,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRRPP1_invdynB_fixb_reg2_snew_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RRRPP1_invdynB_fixb_reg2_snew_vp: qJD has to be [5x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [5 1]), ...
  'S5RRRPP1_invdynB_fixb_reg2_snew_vp: qJDD has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RRRPP1_invdynB_fixb_reg2_snew_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RRRPP1_invdynB_fixb_reg2_snew_vp: pkin has to be [8x1] (double)');

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_tauB_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 20:50:05
% EndTime: 2019-12-31 20:50:19
% DurationCPUTime: 13.59s
% Computational Cost: add. (33131->493), mult. (45559->678), div. (0->0), fcn. (28838->8), ass. (0->366)
t657 = qJD(1) + qJD(2);
t660 = sin(pkin(8));
t661 = cos(pkin(8));
t662 = sin(qJ(3));
t665 = cos(qJ(3));
t612 = (t660 * t665 + t661 * t662) * t657;
t609 = t612 ^ 2;
t668 = qJD(3) ^ 2;
t548 = t668 + t609;
t740 = t657 * t662;
t610 = -t661 * t665 * t657 + t660 * t740;
t566 = t612 * t610;
t753 = qJDD(3) + t566;
t734 = t660 * t753;
t474 = t661 * t548 + t734;
t729 = t661 * t753;
t501 = t660 * t548 - t729;
t439 = t662 * t474 + t665 * t501;
t663 = sin(qJ(2));
t666 = cos(qJ(2));
t602 = qJD(3) * t610;
t701 = qJDD(1) + qJDD(2);
t689 = t662 * t701;
t708 = qJD(3) * t657;
t696 = t665 * t708;
t619 = t689 + t696;
t687 = t665 * t701;
t697 = t662 * t708;
t673 = t687 - t697;
t672 = t661 * t619 + t660 * t673;
t770 = t602 - t672;
t406 = t663 * t439 + t666 * t770;
t412 = t666 * t439 - t663 * t770;
t664 = sin(qJ(1));
t667 = cos(qJ(1));
t352 = t664 * t406 - t667 * t412;
t855 = pkin(5) * t352;
t361 = t667 * t406 + t664 * t412;
t854 = pkin(5) * t361;
t853 = pkin(6) * t406;
t852 = -pkin(1) * t406 - pkin(2) * t770 - pkin(7) * t439;
t419 = t665 * t474 - t662 * t501;
t851 = -pkin(1) * t419 - pkin(6) * t412;
t693 = t660 * t619 - t661 * t673;
t709 = qJD(3) * t612;
t520 = t693 + t709;
t454 = -t660 * t520 - t661 * t770;
t730 = t661 * t520;
t736 = t660 * t770;
t456 = t730 - t736;
t388 = t662 * t454 + t665 * t456;
t746 = t610 ^ 2;
t555 = t746 - t609;
t376 = t663 * t388 - t666 * t555;
t378 = t666 * t388 + t663 * t555;
t850 = t667 * t376 + t664 * t378;
t849 = t664 * t376 - t667 * t378;
t592 = t746 - t668;
t490 = t660 * t592 + t729;
t496 = t661 * t592 - t734;
t435 = t662 * t490 - t665 * t496;
t522 = -t693 + t709;
t402 = t663 * t435 + t666 * t522;
t408 = t666 * t435 - t663 * t522;
t848 = t667 * t402 + t664 * t408;
t847 = t664 * t402 - t667 * t408;
t845 = pkin(7) * t419;
t844 = -pkin(2) * t419 - pkin(3) * t474;
t754 = qJDD(3) - t566;
t733 = t660 * t754;
t755 = -t746 - t668;
t765 = t661 * t755 - t733;
t540 = t661 * t754;
t767 = t660 * t755 + t540;
t784 = -t662 * t767 + t665 * t765;
t802 = t663 * t520 + t666 * t784;
t804 = -t666 * t520 + t663 * t784;
t816 = t664 * t802 + t667 * t804;
t842 = pkin(5) * t816;
t817 = -t664 * t804 + t667 * t802;
t841 = pkin(5) * t817;
t518 = -t746 - t609;
t759 = -t602 - t672;
t768 = t661 * t522 - t660 * t759;
t769 = t660 * t522 + t661 * t759;
t782 = -t662 * t769 + t665 * t768;
t801 = t663 * t518 + t666 * t782;
t803 = -t666 * t518 + t663 * t782;
t818 = t664 * t801 + t667 * t803;
t840 = pkin(5) * t818;
t819 = -t664 * t803 + t667 * t801;
t839 = pkin(5) * t819;
t593 = -t609 + t668;
t787 = t661 * t593 + t733;
t788 = -t660 * t593 + t540;
t799 = -t662 * t787 + t665 * t788;
t821 = -t663 * t759 + t666 * t799;
t822 = t663 * t799 + t666 * t759;
t838 = -t664 * t822 + t667 * t821;
t837 = t664 * t821 + t667 * t822;
t656 = t657 ^ 2;
t688 = t663 * t701;
t627 = -t666 * t656 - t688;
t686 = t666 * t701;
t629 = t663 * t656 - t686;
t570 = t667 * t627 + t664 * t629;
t605 = pkin(6) * t627 + t666 * g(3);
t807 = pkin(6) * t629 - t663 * g(3);
t836 = pkin(5) * t570 + t667 * t605 + t664 * t807;
t571 = t664 * t627 - t667 * t629;
t835 = -pkin(5) * t571 - t664 * t605 + t667 * t807;
t833 = pkin(6) * t803;
t832 = pkin(6) * t804;
t831 = qJ(4) * t474;
t830 = qJ(4) * t501;
t827 = pkin(1) * t803 - pkin(2) * t518 + pkin(7) * t782;
t826 = pkin(1) * t804 - pkin(2) * t520 + pkin(7) * t784;
t781 = t662 * t768 + t665 * t769;
t825 = -pkin(1) * t781 + pkin(6) * t801;
t783 = t662 * t765 + t665 * t767;
t824 = -pkin(1) * t783 + pkin(6) * t802;
t823 = -t665 * t454 + t662 * t456;
t820 = t665 * t490 + t662 * t496;
t813 = pkin(7) * t781;
t812 = pkin(7) * t783;
t809 = t770 * qJ(5);
t368 = -pkin(2) * t781 - pkin(3) * t769;
t808 = -pkin(2) * t783 - pkin(3) * t767;
t800 = t662 * t788 + t665 * t787;
t795 = qJ(4) * t765;
t794 = qJ(4) * t767;
t793 = qJ(4) * t769;
t646 = t667 * g(1) + t664 * g(2);
t669 = qJD(1) ^ 2;
t633 = -t669 * pkin(1) - t646;
t645 = t664 * g(1) - t667 * g(2);
t678 = qJDD(1) * pkin(1) + t645;
t576 = t663 * t633 - t666 * t678;
t577 = t666 * t633 + t663 * t678;
t694 = t663 * t576 + t666 * t577;
t505 = t666 * t576 - t663 * t577;
t713 = t667 * t505;
t786 = -t664 * t694 + t713;
t720 = t664 * t505;
t447 = t667 * t694 + t720;
t785 = -pkin(3) * t518 + qJ(4) * t768;
t700 = t663 * t566;
t706 = qJD(3) * t661;
t698 = t610 * t706;
t676 = t660 * t693 + t698;
t707 = qJD(3) * t660;
t681 = t610 * t707 - t661 * t693;
t750 = -t662 * t681 + t665 * t676;
t762 = t666 * t750 - t700;
t699 = t666 * t566;
t764 = t663 * t750 + t699;
t780 = -t664 * t764 + t667 * t762;
t779 = t664 * t762 + t667 * t764;
t702 = t666 * qJDD(3);
t675 = (-t610 * t660 - t661 * t612) * qJD(3);
t591 = t612 * t707;
t680 = t591 - t698;
t751 = -t662 * t675 + t665 * t680;
t763 = t663 * t751 - t702;
t652 = t663 * qJDD(3);
t766 = t666 * t751 + t652;
t778 = t664 * t766 + t667 * t763;
t777 = -t664 * t763 + t667 * t766;
t776 = 2 * qJD(5);
t659 = t662 ^ 2;
t745 = t665 ^ 2;
t757 = t659 + t745;
t542 = t610 * pkin(4) - t612 * qJ(5);
t562 = -t656 * pkin(2) + t701 * pkin(7) + t577;
t726 = t662 * t562;
t741 = t656 * t662;
t481 = qJDD(3) * pkin(3) - t619 * qJ(4) - t726 + (pkin(3) * t741 + qJ(4) * t708 - g(3)) * t665;
t539 = -t662 * g(3) + t665 * t562;
t634 = qJD(3) * pkin(3) - qJ(4) * t740;
t649 = t745 * t656;
t482 = -pkin(3) * t649 + qJ(4) * t673 - qJD(3) * t634 + t539;
t712 = t660 * t481 + t661 * t482;
t756 = qJDD(3) * qJ(5) + qJD(3) * t776 - t610 * t542 + t712;
t752 = t662 * t680 + t665 * t675;
t749 = t662 * t676 + t665 * t681;
t511 = t612 * t706 + t660 * t672;
t512 = t661 * t672 - t591;
t451 = -t662 * t511 + t665 * t512;
t682 = t666 * t451 + t700;
t683 = t663 * t451 - t699;
t748 = -t664 * t683 + t667 * t682;
t747 = t664 * t682 + t667 * t683;
t744 = pkin(4) * t661;
t739 = t659 * t656;
t561 = -t701 * pkin(2) - t656 * pkin(7) + t576;
t483 = -t673 * pkin(3) - qJ(4) * t649 + t634 * t740 + qJDD(4) + t561;
t738 = t660 * t483;
t731 = t661 * t483;
t704 = qJD(4) * t612;
t601 = 0.2e1 * t704;
t711 = -t661 * t481 + t660 * t482;
t416 = t601 + t711;
t705 = qJD(4) * t610;
t599 = -0.2e1 * t705;
t417 = t599 + t712;
t366 = -t661 * t416 + t660 * t417;
t728 = t662 * t366;
t727 = t662 * t561;
t644 = t665 * t741;
t635 = qJDD(3) + t644;
t725 = t662 * t635;
t636 = qJDD(3) - t644;
t724 = t662 * t636;
t719 = t665 * t366;
t718 = t665 * t561;
t620 = t687 - 0.2e1 * t697;
t717 = t665 * t620;
t716 = t665 * t636;
t710 = -t518 - t668;
t695 = -qJ(5) * t660 - pkin(3);
t367 = t660 * t416 + t661 * t417;
t538 = t665 * g(3) + t726;
t472 = t662 * t538 + t665 * t539;
t590 = -t664 * t645 - t667 * t646;
t692 = t663 * t644;
t691 = t666 * t644;
t685 = t612 * t542 + qJDD(5) + t711;
t638 = t667 * qJDD(1) - t664 * t669;
t684 = -pkin(5) * t638 - t664 * g(3);
t471 = t665 * t538 - t662 * t539;
t589 = t667 * t645 - t664 * t646;
t677 = t599 + t756;
t674 = -qJDD(3) * pkin(4) + t685;
t671 = t693 * pkin(4) + t483 + t809;
t670 = t612 * t776 - t671;
t643 = -t649 - t668;
t642 = t649 - t668;
t641 = -t668 - t739;
t640 = t668 - t739;
t637 = t664 * qJDD(1) + t667 * t669;
t631 = t649 - t739;
t630 = t649 + t739;
t625 = t665 * t635;
t624 = t757 * t701;
t618 = t689 + 0.2e1 * t696;
t617 = -pkin(5) * t637 + t667 * g(3);
t616 = t757 * t708;
t588 = t666 * t616 + t652;
t587 = t663 * t616 - t702;
t585 = t665 * t619 - t659 * t708;
t584 = -t662 * t673 - t745 * t708;
t583 = -t662 * t641 - t716;
t582 = -t662 * t640 + t625;
t581 = t665 * t643 - t725;
t580 = t665 * t642 - t724;
t579 = t665 * t641 - t724;
t578 = t662 * t643 + t625;
t572 = t666 * t624 - t663 * t630;
t569 = t663 * t624 + t666 * t630;
t568 = -t662 * t618 + t717;
t560 = t666 * t582 + t662 * t688;
t559 = t666 * t580 + t663 * t687;
t558 = t663 * t582 - t662 * t686;
t557 = t663 * t580 - t665 * t686;
t546 = t666 * t585 - t692;
t545 = t666 * t584 + t692;
t544 = t663 * t585 + t691;
t543 = t663 * t584 - t691;
t533 = t666 * t583 + t663 * t618;
t532 = t666 * t581 - t663 * t620;
t531 = t663 * t583 - t666 * t618;
t530 = t663 * t581 + t666 * t620;
t528 = t666 * t568 - t663 * t631;
t527 = t663 * t568 + t666 * t631;
t502 = pkin(1) * g(3) + pkin(6) * t694;
t489 = -pkin(7) * t579 + t718;
t488 = -pkin(7) * t578 + t727;
t487 = -pkin(2) * t579 + t539;
t486 = -pkin(2) * t578 + t538;
t485 = -t664 * t569 + t667 * t572;
t484 = t667 * t569 + t664 * t572;
t467 = -t664 * t531 + t667 * t533;
t466 = -t664 * t530 + t667 * t532;
t465 = t667 * t531 + t664 * t533;
t464 = t667 * t530 + t664 * t532;
t448 = t665 * t511 + t662 * t512;
t445 = -pkin(6) * t569 + t666 * t471;
t444 = pkin(6) * t572 + t663 * t471;
t441 = t666 * t472 + t663 * t561;
t440 = t663 * t472 - t666 * t561;
t431 = t731 + t831;
t430 = t738 - t794;
t414 = -pkin(6) * t531 - t663 * t487 + t666 * t489;
t413 = -pkin(6) * t530 - t663 * t486 + t666 * t488;
t400 = pkin(3) * t770 + t738 + t830;
t399 = (pkin(4) * qJD(3) - (2 * qJD(5))) * t612 + t671;
t398 = -pkin(1) * t579 + pkin(6) * t533 + t666 * t487 + t663 * t489;
t397 = -pkin(1) * t578 + pkin(6) * t532 + t666 * t486 + t663 * t488;
t396 = -pkin(3) * t520 - t731 + t795;
t385 = t668 * qJ(5) - t674 - 0.2e1 * t704;
t384 = -t668 * pkin(4) + t677;
t383 = t670 + (-t520 - t709) * pkin(4);
t382 = -pkin(4) * t709 + t670 - t809;
t381 = -t664 * t440 + t667 * t441;
t380 = t667 * t440 + t664 * t441;
t371 = t710 * qJ(5) + t601 + t674;
t370 = t710 * pkin(4) + t677;
t369 = -pkin(6) * t440 - (pkin(2) * t663 - pkin(7) * t666) * t471;
t365 = t417 - t844;
t364 = -qJ(5) * t730 - t660 * t383 - t794;
t359 = pkin(4) * t736 + t661 * t382 - t831;
t358 = pkin(6) * t441 - (-pkin(2) * t666 - pkin(7) * t663 - pkin(1)) * t471;
t357 = t416 + t808;
t356 = t661 * t383 + t695 * t520 + t795;
t355 = -pkin(3) * t483 + qJ(4) * t367;
t354 = -t830 + t660 * t382 - (pkin(3) + t744) * t770;
t353 = -pkin(4) * t759 - qJ(5) * t522 + t368;
t348 = -t366 - t793;
t347 = -t662 * t400 + t665 * t431 + t845;
t346 = t661 * t384 - t660 * t385;
t345 = t660 * t384 + t661 * t385;
t344 = t601 + (-t755 - t668) * qJ(5) + (-qJDD(3) - t754) * pkin(4) + t685 + t808;
t343 = t367 + t785;
t342 = -t662 * t396 + t665 * t430 - t812;
t341 = -qJ(5) * t753 + 0.2e1 * t705 + (-t548 + t668) * pkin(4) - t756 + t844;
t336 = -t660 * t370 + t661 * t371 - t793;
t335 = t665 * t367 - t728;
t334 = t662 * t367 + t719;
t333 = t661 * t370 + t660 * t371 + t785;
t332 = t666 * t335 + t663 * t483;
t331 = t663 * t335 - t666 * t483;
t330 = -qJ(4) * t345 + (pkin(4) * t660 - qJ(5) * t661) * t399;
t329 = -t662 * t356 + t665 * t364 - t812;
t328 = -t662 * t345 + t665 * t346;
t327 = t665 * t345 + t662 * t346;
t326 = -t662 * t354 + t665 * t359 - t845;
t325 = t666 * t347 - t663 * t365 - t853;
t324 = -pkin(2) * t334 - pkin(3) * t366;
t323 = qJ(4) * t346 + (t695 - t744) * t399;
t322 = t663 * t347 + t666 * t365 - t851;
t321 = t666 * t328 + t663 * t399;
t320 = t663 * t328 - t666 * t399;
t319 = t666 * t342 - t663 * t357 - t832;
t318 = -t662 * t343 + t665 * t348 - t813;
t317 = t663 * t342 + t666 * t357 + t824;
t316 = -pkin(7) * t334 - qJ(4) * t719 - t662 * t355;
t315 = -t664 * t331 + t667 * t332;
t314 = t667 * t331 + t664 * t332;
t313 = t666 * t329 - t663 * t344 - t832;
t312 = -t662 * t333 + t665 * t336 - t813;
t311 = t663 * t329 + t666 * t344 + t824;
t310 = t666 * t326 - t663 * t341 + t853;
t309 = t666 * t318 - t663 * t368 - t833;
t308 = -pkin(2) * t327 - pkin(3) * t345 - pkin(4) * t385 - qJ(5) * t384;
t307 = t663 * t326 + t666 * t341 + t851;
t306 = t663 * t318 + t666 * t368 + t825;
t305 = -t664 * t320 + t667 * t321;
t304 = t667 * t320 + t664 * t321;
t303 = t666 * t312 - t663 * t353 - t833;
t302 = t663 * t312 + t666 * t353 + t825;
t301 = -pkin(7) * t327 - t662 * t323 + t665 * t330;
t300 = -pkin(6) * t331 + t666 * t316 - t663 * t324;
t299 = -pkin(1) * t334 + pkin(6) * t332 + t663 * t316 + t666 * t324;
t298 = -pkin(6) * t320 + t666 * t301 - t663 * t308;
t297 = -pkin(1) * t327 + pkin(6) * t321 + t663 * t301 + t666 * t308;
t1 = [0, 0, 0, 0, 0, 0, 0, 0, 0, -g(1), 0, 0, 0, 0, 0, 0, -t637, -t638, 0, t590, 0, 0, 0, 0, 0, 0, t570, -t571, 0, t447, 0, 0, 0, 0, 0, 0, t466, t467, t485, t381, 0, 0, 0, 0, 0, 0, t817, -t352, t819, t315, 0, 0, 0, 0, 0, 0, t817, t819, t352, t305; 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(2), 0, 0, 0, 0, 0, 0, t638, -t637, 0, t589, 0, 0, 0, 0, 0, 0, t571, t570, 0, -t786, 0, 0, 0, 0, 0, 0, t464, t465, t484, t380, 0, 0, 0, 0, 0, 0, t816, t361, t818, t314, 0, 0, 0, 0, 0, 0, t816, t818, -t361, t304; 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(3), 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(3), 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(3), 0, 0, 0, 0, 0, 0, t578, t579, 0, -t471, 0, 0, 0, 0, 0, 0, t783, -t419, t781, t334, 0, 0, 0, 0, 0, 0, t783, t781, t419, t327; 0, 0, 0, 0, 0, 0, 0, -g(3), g(2), 0, 0, 0, t638, 0, -t637, 0, t684, -t617, -t589, -pkin(5) * t589, 0, 0, t571, 0, t570, 0, t835, -t836, t786, pkin(5) * t786 + pkin(6) * t713 - t664 * t502, -t664 * t544 + t667 * t546, -t664 * t527 + t667 * t528, -t664 * t558 + t667 * t560, -t664 * t543 + t667 * t545, -t664 * t557 + t667 * t559, -t664 * t587 + t667 * t588, -pkin(5) * t464 - t664 * t397 + t667 * t413, -pkin(5) * t465 - t664 * t398 + t667 * t414, -pkin(5) * t484 - t664 * t444 + t667 * t445, -pkin(5) * t380 - t664 * t358 + t667 * t369, t748, t849, t838, t780, t847, t777, -t664 * t317 + t667 * t319 - t842, -t664 * t322 + t667 * t325 - t854, -t664 * t306 + t667 * t309 - t840, -pkin(5) * t314 - t664 * t299 + t667 * t300, t748, t838, -t849, t777, -t847, t780, -t664 * t311 + t667 * t313 - t842, -t664 * t302 + t667 * t303 - t840, -t664 * t307 + t667 * t310 + t854, -pkin(5) * t304 - t664 * t297 + t667 * t298; 0, 0, 0, 0, 0, 0, g(3), 0, -g(1), 0, 0, 0, t637, 0, t638, 0, t617, t684, t590, pkin(5) * t590, 0, 0, -t570, 0, t571, 0, t836, t835, t447, pkin(5) * t447 + pkin(6) * t720 + t667 * t502, t667 * t544 + t664 * t546, t667 * t527 + t664 * t528, t667 * t558 + t664 * t560, t667 * t543 + t664 * t545, t667 * t557 + t664 * t559, t667 * t587 + t664 * t588, pkin(5) * t466 + t667 * t397 + t664 * t413, pkin(5) * t467 + t667 * t398 + t664 * t414, pkin(5) * t485 + t667 * t444 + t664 * t445, pkin(5) * t381 + t667 * t358 + t664 * t369, t747, -t850, t837, t779, -t848, t778, t667 * t317 + t664 * t319 + t841, t667 * t322 + t664 * t325 - t855, t667 * t306 + t664 * t309 + t839, pkin(5) * t315 + t667 * t299 + t664 * t300, t747, t837, t850, t778, t848, t779, t667 * t311 + t664 * t313 + t841, t667 * t302 + t664 * t303 + t839, t667 * t307 + t664 * t310 + t855, pkin(5) * t305 + t667 * t297 + t664 * t298; 0, 0, 0, 0, 0, 0, -g(2), g(1), 0, 0, 0, 0, 0, 0, 0, qJDD(1), t645, t646, 0, 0, 0, 0, 0, 0, 0, t701, -pkin(1) * t629 - t576, pkin(1) * t627 - t577, 0, -pkin(1) * t505, (t619 + t696) * t662, t665 * t618 + t662 * t620, t665 * t640 + t725, t717, t662 * t642 + t716, 0, pkin(1) * t530 + pkin(2) * t620 + pkin(7) * t581 - t718, pkin(1) * t531 - pkin(2) * t618 + pkin(7) * t583 + t727, pkin(1) * t569 + pkin(2) * t630 + pkin(7) * t624 + t472, pkin(1) * t440 - pkin(2) * t561 + pkin(7) * t472, t448, -t823, t800, t749, t820, t752, t665 * t396 + t662 * t430 + t826, t665 * t400 + t662 * t431 - t852, t665 * t343 + t662 * t348 + t827, pkin(1) * t331 - pkin(2) * t483 + pkin(7) * t335 - qJ(4) * t728 + t665 * t355, t448, t800, t823, t752, -t820, t749, t665 * t356 + t662 * t364 + t826, t665 * t333 + t662 * t336 + t827, t665 * t354 + t662 * t359 + t852, pkin(1) * t320 - pkin(2) * t399 + pkin(7) * t328 + t665 * t323 + t662 * t330;];
tauB_reg = t1;
