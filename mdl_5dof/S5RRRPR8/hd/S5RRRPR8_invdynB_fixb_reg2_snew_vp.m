% Calculate inertial parameters regressor of inverse dynamics base forces vector with Newton-Euler for
% S5RRRPR8
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
%   pkin=[a2,a3,a4,a5,d1,d2,d3,d5]';
% 
% Output:
% tauB_reg [6x(6*10)]
%   inertial parameter regressor of inverse dynamics base forces vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 21:21
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauB_reg = S5RRRPR8_invdynB_fixb_reg2_snew_vp(qJ, qJD, qJDD, g, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(5,1),zeros(3,1),zeros(8,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRRPR8_invdynB_fixb_reg2_snew_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RRRPR8_invdynB_fixb_reg2_snew_vp: qJD has to be [5x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [5 1]), ...
  'S5RRRPR8_invdynB_fixb_reg2_snew_vp: qJDD has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RRRPR8_invdynB_fixb_reg2_snew_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RRRPR8_invdynB_fixb_reg2_snew_vp: pkin has to be [8x1] (double)');

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_tauB_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 21:20:39
% EndTime: 2019-12-31 21:20:55
% DurationCPUTime: 16.13s
% Computational Cost: add. (32924->574), mult. (69387->820), div. (0->0), fcn. (47382->8), ass. (0->409)
t704 = qJD(2) + qJD(3);
t702 = t704 ^ 2;
t708 = sin(qJ(3));
t709 = sin(qJ(2));
t712 = cos(qJ(3));
t713 = cos(qJ(2));
t663 = (t708 * t713 + t709 * t712) * qJD(1);
t797 = t663 ^ 2;
t635 = t797 + t702;
t703 = qJDD(2) + qJDD(3);
t760 = qJD(1) * t709;
t661 = -t712 * t713 * qJD(1) + t708 * t760;
t790 = t663 * t661;
t828 = t790 + t703;
t836 = t708 * t828;
t554 = t712 * t635 + t836;
t835 = t712 * t828;
t556 = -t708 * t635 + t835;
t495 = t713 * t554 + t709 * t556;
t860 = pkin(6) * t495;
t498 = t709 * t554 - t713 * t556;
t859 = pkin(6) * t498;
t710 = sin(qJ(1));
t858 = t710 * t498;
t643 = -t797 + t702;
t610 = t790 - t703;
t841 = t708 * t610;
t560 = -t712 * t643 + t841;
t840 = t712 * t610;
t564 = t708 * t643 + t840;
t504 = t709 * t560 - t713 * t564;
t857 = t710 * t504;
t798 = t661 ^ 2;
t642 = t798 - t702;
t561 = t708 * t642 + t835;
t566 = -t712 * t642 + t836;
t507 = t709 * t561 + t713 * t566;
t856 = t710 * t507;
t714 = cos(qJ(1));
t855 = t714 * t498;
t854 = t714 * t504;
t853 = t714 * t507;
t852 = -pkin(1) * t495 - pkin(2) * t554;
t850 = pkin(7) * t554;
t849 = pkin(7) * t556;
t848 = t713 * t560 + t709 * t564;
t847 = t713 * t561 - t709 * t566;
t605 = -t702 - t798;
t529 = t708 * t605 - t840;
t532 = -t712 * t605 - t841;
t473 = t713 * t529 - t709 * t532;
t846 = pkin(6) * t473;
t476 = t709 * t529 + t713 * t532;
t845 = pkin(6) * t476;
t844 = t710 * t476;
t843 = t714 * t476;
t842 = -pkin(1) * t473 - pkin(2) * t529;
t838 = pkin(7) * t529;
t837 = pkin(7) * t532;
t827 = 2 * qJD(4);
t805 = -t798 - t797;
t826 = pkin(1) * t805;
t825 = pkin(2) * t805;
t707 = sin(qJ(5));
t711 = cos(qJ(5));
t624 = -t711 * t661 + t707 * t704;
t626 = t707 * t661 + t711 * t704;
t577 = t626 * t624;
t756 = qJD(1) * qJD(2);
t745 = t713 * t756;
t755 = t709 * qJDD(1);
t673 = t745 + t755;
t698 = t713 * qJDD(1);
t746 = t709 * t756;
t674 = t698 - t746;
t734 = t712 * t673 + t708 * t674;
t723 = t661 * qJD(3) - t734;
t584 = qJDD(5) - t723;
t809 = -t577 + t584;
t824 = t707 * t809;
t823 = t710 * t805;
t806 = -t797 + t798;
t822 = t710 * t806;
t821 = t711 * t809;
t820 = t714 * t805;
t819 = t714 * t806;
t706 = t713 ^ 2;
t716 = qJD(1) ^ 2;
t683 = t710 * g(1) - t714 * g(2);
t729 = qJDD(1) * pkin(1) + t683;
t730 = qJD(2) * pkin(2) - pkin(7) * t760;
t592 = t674 * pkin(2) + (pkin(7) * t706 + pkin(6)) * t716 - t730 * t760 + t729;
t788 = t704 * t663;
t818 = pkin(3) * t788 - t663 * t827 - t592;
t648 = t661 * t704;
t817 = -t723 + t648;
t742 = t723 + t648;
t725 = (-t661 * t708 - t663 * t712) * t704;
t726 = (-t661 * t712 + t663 * t708) * t704;
t803 = -t709 * t725 + t713 * t726;
t816 = -t714 * t703 + t710 * t803;
t748 = t714 * t790;
t787 = t704 * t708;
t735 = -t663 * t787 - t712 * t723;
t786 = t704 * t712;
t737 = t663 * t786 - t708 * t723;
t800 = -t709 * t737 + t713 * t735;
t815 = t710 * t800 - t748;
t741 = t708 * t673 - t712 * t674;
t588 = t663 * qJD(3) + t741;
t728 = t708 * t588 + t661 * t786;
t736 = -t712 * t588 + t661 * t787;
t802 = -t709 * t736 + t713 * t728;
t814 = t710 * t802 + t748;
t750 = t710 * t790;
t813 = t714 * t800 + t750;
t812 = t714 * t802 - t750;
t811 = t710 * t703 + t714 * t803;
t810 = t742 * qJ(4);
t804 = t709 * t726 + t713 * t725;
t801 = t709 * t728 + t713 * t736;
t799 = t709 * t735 + t713 * t737;
t621 = t624 ^ 2;
t622 = t626 ^ 2;
t655 = qJD(5) + t663;
t653 = t655 ^ 2;
t796 = pkin(3) + pkin(8);
t795 = pkin(3) * t708;
t794 = pkin(3) * t712;
t793 = t655 * t624;
t792 = t655 * t707;
t791 = t655 * t711;
t789 = t703 * qJ(4);
t705 = t709 ^ 2;
t785 = t705 * t716;
t700 = t706 * t716;
t639 = t663 * pkin(4) - t704 * pkin(8);
t773 = t709 * t716;
t684 = t714 * g(1) + t710 * g(2);
t666 = -t716 * pkin(1) + qJDD(1) * pkin(6) - t684;
t776 = t709 * t666;
t571 = qJDD(2) * pkin(2) - t673 * pkin(7) - t776 + (pkin(2) * t773 + pkin(7) * t756 - g(3)) * t713;
t638 = -t709 * g(3) + t713 * t666;
t576 = -pkin(2) * t700 + t674 * pkin(7) - qJD(2) * t730 + t638;
t518 = t708 * t571 + t712 * t576;
t612 = t661 * pkin(3) - t663 * qJ(4);
t732 = -t702 * pkin(3) - t661 * t612 + t518;
t441 = t789 - t588 * pkin(4) - t798 * pkin(8) + (t827 + t639) * t704 + t732;
t784 = t707 * t441;
t521 = t577 + t584;
t783 = t707 * t521;
t782 = t708 * t592;
t517 = -t712 * t571 + t708 * t576;
t452 = -t712 * t517 + t708 * t518;
t778 = t709 * t452;
t665 = t716 * pkin(6) + t729;
t777 = t709 * t665;
t690 = t713 * t773;
t681 = qJDD(2) + t690;
t775 = t709 * t681;
t682 = qJDD(2) - t690;
t774 = t709 * t682;
t771 = t711 * t441;
t770 = t711 * t521;
t546 = t588 + t788;
t769 = t712 * t546;
t768 = t712 * t592;
t764 = t713 * t452;
t763 = t713 * t665;
t762 = t713 * t682;
t761 = t705 + t706;
t758 = qJD(3) - t704;
t757 = qJD(3) + t704;
t754 = t710 * qJDD(1);
t753 = t714 * qJDD(1);
t752 = -t622 - t653;
t751 = t708 * t577;
t749 = t712 * t577;
t527 = -t624 * qJD(5) + t707 * t588 + t711 * t703;
t744 = qJ(4) * t708 + pkin(2);
t453 = t708 * t517 + t712 * t518;
t743 = -t711 * t588 + t707 * t703;
t637 = t713 * g(3) + t776;
t575 = t709 * t637 + t713 * t638;
t628 = -t710 * t683 - t714 * t684;
t740 = t710 * t690;
t739 = t714 * t690;
t678 = -t710 * t716 + t753;
t738 = -pkin(5) * t678 - t710 * g(3);
t717 = t810 + t818;
t425 = -pkin(4) * t798 + t796 * t588 - t663 * t639 + t717;
t727 = -t703 * pkin(3) - t702 * qJ(4) + t663 * t612 + qJDD(4) + t517;
t440 = t817 * pkin(4) + pkin(8) * t610 + t727;
t388 = t707 * t425 - t711 * t440;
t389 = t711 * t425 + t707 * t440;
t355 = -t711 * t388 + t707 * t389;
t356 = t707 * t388 + t711 * t389;
t574 = t713 * t637 - t709 * t638;
t627 = t714 * t683 - t710 * t684;
t731 = t527 - t793;
t724 = (-qJD(5) + t655) * t626 - t743;
t722 = -0.2e1 * qJD(4) * t704 - t732;
t720 = -t757 * t661 + t734;
t467 = -t722 + t789;
t718 = -t588 * pkin(3) - t818;
t715 = qJD(2) ^ 2;
t688 = -t700 - t715;
t687 = t700 - t715;
t686 = -t715 - t785;
t685 = t715 - t785;
t680 = t700 - t785;
t679 = t700 + t785;
t677 = t714 * t716 + t754;
t676 = t761 * qJDD(1);
t675 = t698 - 0.2e1 * t746;
t672 = 0.2e1 * t745 + t755;
t670 = t713 * t681;
t669 = t761 * t756;
t656 = -pkin(5) * t677 + t714 * g(3);
t641 = t713 * t673 - t705 * t756;
t640 = -t709 * t674 - t706 * t756;
t634 = -t709 * t686 - t762;
t633 = -t709 * t685 + t670;
t632 = t713 * t688 - t775;
t631 = t713 * t687 - t774;
t630 = t713 * t686 - t774;
t629 = t709 * t688 + t670;
t620 = t714 * t676 - t710 * t679;
t619 = t710 * t676 + t714 * t679;
t617 = -t709 * t672 + t713 * t675;
t604 = t714 * t634 + t710 * t672;
t603 = t714 * t632 - t710 * t675;
t602 = t710 * t634 - t714 * t672;
t601 = t710 * t632 + t714 * t675;
t600 = -t622 + t653;
t599 = t621 - t653;
t594 = -pkin(6) * t630 - t763;
t593 = -pkin(6) * t629 - t777;
t587 = -pkin(1) * t630 + t638;
t586 = -pkin(1) * t629 + t637;
t572 = t622 - t621;
t548 = -t758 * t661 + t734;
t547 = -t588 + t788;
t545 = t758 * t663 + t741;
t544 = t757 * t663 + t741;
t543 = -t653 - t621;
t534 = t714 * t575 - t710 * t665;
t533 = t710 * t575 + t714 * t665;
t528 = -t621 - t622;
t526 = -t626 * qJD(5) - t743;
t525 = (t624 * t711 - t626 * t707) * t655;
t524 = (t624 * t707 + t626 * t711) * t655;
t515 = -t768 + t850;
t513 = t527 + t793;
t509 = (qJD(5) + t655) * t626 + t743;
t508 = -t782 - t838;
t502 = -t711 * t527 + t626 * t792;
t501 = -t707 * t527 - t626 * t791;
t500 = t707 * t526 - t624 * t791;
t499 = -t711 * t526 - t624 * t792;
t494 = t712 * t547 + t708 * t817;
t493 = -t708 * t720 - t769;
t492 = -t712 * t544 + t708 * t742;
t491 = -t712 * t545 + t708 * t548;
t490 = t708 * t547 - t712 * t817;
t489 = -t708 * t546 + t712 * t720;
t488 = -t708 * t544 - t712 * t742;
t487 = -t708 * t545 - t712 * t548;
t486 = -t708 * t524 + t712 * t584;
t485 = t712 * t524 + t708 * t584;
t484 = -t711 * t599 + t783;
t483 = t707 * t600 - t821;
t482 = -t707 * t599 - t770;
t481 = -t711 * t600 - t824;
t472 = -t707 * t752 - t770;
t471 = t711 * t752 - t783;
t469 = t711 * t543 - t824;
t468 = t707 * t543 + t821;
t466 = pkin(2) * t742 - t782 - t849;
t465 = -t708 * t501 + t749;
t464 = -t708 * t499 - t749;
t463 = t712 * t501 + t751;
t462 = t712 * t499 - t751;
t461 = -pkin(2) * t544 + t768 - t837;
t460 = t718 - t810;
t459 = -t710 * t742 + t855;
t458 = -t710 * t720 - t855;
t457 = t714 * t742 + t858;
t456 = t714 * t720 - t858;
t455 = -qJ(4) * t805 + t727;
t454 = -pkin(3) * t805 + t467;
t451 = -t710 * t546 + t843;
t450 = t710 * t544 - t843;
t449 = t714 * t546 + t844;
t448 = -t714 * t544 - t844;
t447 = t707 * t513 + t711 * t724;
t446 = t711 * t509 + t707 * t731;
t445 = -t711 * t513 + t707 * t724;
t444 = t707 * t509 - t711 * t731;
t443 = (t546 + t588) * pkin(3) + t717;
t442 = (t720 - t742) * qJ(4) + t718;
t439 = pkin(2) * t592 + pkin(7) * t453;
t438 = -t709 * t490 + t713 * t494;
t437 = -t709 * t489 + t713 * t493;
t436 = -t709 * t488 + t713 * t492;
t435 = -t709 * t487 + t713 * t491;
t434 = t713 * t490 + t709 * t494;
t433 = t713 * t487 + t709 * t491;
t432 = -t708 * t481 + t712 * t513;
t431 = -t708 * t482 + t712 * t724;
t430 = t712 * t481 + t708 * t513;
t429 = t712 * t482 + t708 * t724;
t428 = -t709 * t485 + t713 * t486;
t427 = t708 * t471 + t712 * t731;
t426 = -t712 * t471 + t708 * t731;
t424 = t708 * t468 + t712 * t509;
t423 = -t712 * t468 + t708 * t509;
t422 = t518 - t852;
t421 = -t708 * t444 + t712 * t572;
t420 = t712 * t444 + t708 * t572;
t419 = t714 * t438 + t823;
t418 = t714 * t435 + t823;
t417 = t710 * t438 - t820;
t416 = t710 * t435 - t820;
t415 = t708 * t445 + t712 * t528;
t414 = -t712 * t445 + t708 * t528;
t413 = t517 + t842;
t412 = t712 * t467 + t708 * t727;
t411 = t708 * t467 - t712 * t727;
t410 = -pkin(7) * t490 - t452;
t409 = -t709 * t463 + t713 * t465;
t408 = -t709 * t462 + t713 * t464;
t407 = t712 * t442 - t720 * t795 - t850;
t406 = qJ(4) * t769 - t708 * t443 + t838;
t405 = pkin(7) * t494 + t453 - t825;
t404 = -pkin(1) * t434 - pkin(2) * t490;
t403 = -pkin(3) * t635 + (-t828 - t703) * qJ(4) + t722 + t852;
t402 = -t709 * t466 + t713 * t515 + t860;
t401 = t713 * t453 - t778;
t400 = t709 * t453 + t764;
t399 = t849 + t708 * t442 + (pkin(2) + t794) * t720;
t398 = -pkin(3) * t610 + qJ(4) * t605 - t727 - t842;
t397 = t712 * t443 + t744 * t546 + t837;
t396 = -t709 * t461 + t713 * t508 - t846;
t395 = t714 * t401 - t710 * t592;
t394 = t710 * t401 + t714 * t592;
t393 = pkin(4) * t445 - qJ(4) * t447;
t392 = -pkin(1) * t433 - pkin(2) * t487 + pkin(3) * t548 + qJ(4) * t545;
t391 = -t709 * t430 + t713 * t432;
t390 = -t709 * t429 + t713 * t431;
t386 = -pkin(7) * t487 - t708 * t454 + t712 * t455;
t385 = -t709 * t426 + t713 * t427;
t384 = t713 * t426 + t709 * t427;
t383 = -t709 * t423 + t713 * t424;
t382 = t713 * t423 + t709 * t424;
t381 = pkin(7) * t491 + t712 * t454 + t708 * t455 - t825;
t380 = -t709 * t420 + t713 * t421;
t379 = -t709 * t414 + t713 * t415;
t378 = t713 * t414 + t709 * t415;
t377 = -pkin(1) * t400 - pkin(2) * t452;
t376 = -t709 * t411 + t713 * t412;
t375 = t713 * t411 + t709 * t412;
t374 = pkin(4) * t731 - t796 * t472 - t784;
t373 = pkin(4) * t509 - t796 * t469 + t771;
t372 = -pkin(7) * t411 + (qJ(4) * t712 - t795) * t460;
t371 = t714 * t385 + t710 * t472;
t370 = t710 * t385 - t714 * t472;
t369 = t714 * t383 + t710 * t469;
t368 = t710 * t383 - t714 * t469;
t367 = pkin(4) * t471 - qJ(4) * t472 - t389;
t366 = t714 * t376 - t710 * t460;
t365 = t710 * t376 + t714 * t460;
t364 = pkin(4) * t468 - qJ(4) * t469 - t388;
t363 = t714 * t379 + t710 * t447;
t362 = t710 * t379 - t714 * t447;
t361 = pkin(7) * t412 + (t744 + t794) * t460;
t360 = -pkin(6) * t400 - pkin(7) * t764 - t709 * t439;
t359 = -t709 * t399 + t713 * t407 - t860;
t358 = -t709 * t397 + t713 * t406 + t846;
t357 = -pkin(6) * t434 - t709 * t405 + t713 * t410;
t354 = -pkin(1) * t375 - pkin(2) * t411 + pkin(3) * t727 - qJ(4) * t467;
t353 = t708 * t355 + t712 * t441;
t352 = -t712 * t355 + t708 * t441;
t351 = -pkin(6) * t433 - t709 * t381 + t713 * t386;
t350 = -pkin(1) * t384 - pkin(2) * t426 - qJ(4) * t731 + t796 * t471 - t771;
t349 = -pkin(1) * t382 - pkin(2) * t423 - qJ(4) * t509 + t796 * t468 - t784;
t348 = pkin(4) * t528 - t796 * t447 - t356;
t347 = -pkin(7) * t426 + t712 * t367 - t708 * t374;
t346 = -pkin(7) * t423 + t712 * t364 - t708 * t373;
t345 = -pkin(2) * t472 + pkin(7) * t427 + t708 * t367 + t712 * t374;
t344 = -pkin(2) * t469 + pkin(7) * t424 + t708 * t364 + t712 * t373;
t343 = pkin(4) * t355 - qJ(4) * t356;
t342 = -pkin(6) * t375 - t709 * t361 + t713 * t372;
t341 = pkin(4) * t441 - t796 * t356;
t340 = -pkin(7) * t414 - t708 * t348 + t712 * t393;
t339 = -pkin(1) * t378 - pkin(2) * t414 - qJ(4) * t528 + t796 * t445 + t355;
t338 = -t709 * t352 + t713 * t353;
t337 = t713 * t352 + t709 * t353;
t336 = -pkin(2) * t447 + pkin(7) * t415 + t712 * t348 + t708 * t393;
t335 = t714 * t338 + t710 * t356;
t334 = t710 * t338 - t714 * t356;
t333 = -pkin(6) * t384 - t709 * t345 + t713 * t347;
t332 = -pkin(6) * t382 - t709 * t344 + t713 * t346;
t331 = -pkin(6) * t378 - t709 * t336 + t713 * t340;
t330 = -pkin(7) * t352 - t708 * t341 + t712 * t343;
t329 = -pkin(1) * t337 - pkin(2) * t352 - qJ(4) * t441 + t796 * t355;
t328 = -pkin(2) * t356 + pkin(7) * t353 + t712 * t341 + t708 * t343;
t327 = -pkin(6) * t337 - t709 * t328 + t713 * t330;
t1 = [0, 0, 0, 0, 0, 0, 0, 0, 0, -g(1), 0, 0, 0, 0, 0, 0, -t677, -t678, 0, t628, 0, 0, 0, 0, 0, 0, t603, t604, t620, t534, 0, 0, 0, 0, 0, 0, t450, t459, t419, t395, 0, 0, 0, 0, 0, 0, t418, t451, t458, t366, 0, 0, 0, 0, 0, 0, t369, t371, t363, t335; 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(2), 0, 0, 0, 0, 0, 0, t678, -t677, 0, t627, 0, 0, 0, 0, 0, 0, t601, t602, t619, t533, 0, 0, 0, 0, 0, 0, t448, t457, t417, t394, 0, 0, 0, 0, 0, 0, t416, t449, t456, t365, 0, 0, 0, 0, 0, 0, t368, t370, t362, t334; 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(3), 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(3), 0, 0, 0, 0, 0, 0, t629, t630, 0, -t574, 0, 0, 0, 0, 0, 0, t473, -t495, t434, t400, 0, 0, 0, 0, 0, 0, t433, -t473, t495, t375, 0, 0, 0, 0, 0, 0, t382, t384, t378, t337; 0, 0, 0, 0, 0, 0, 0, -g(3), g(2), 0, 0, 0, t678, 0, -t677, 0, t738, -t656, -t627, -pkin(5) * t627, t714 * t641 - t740, t714 * t617 - t710 * t680, t714 * t633 + t709 * t754, t714 * t640 + t740, t714 * t631 + t698 * t710, t710 * qJDD(2) + t714 * t669, -pkin(5) * t601 - t710 * t586 + t714 * t593, -pkin(5) * t602 - t710 * t587 + t714 * t594, -pkin(5) * t619 + t714 * t574, -pkin(5) * t533 - (pkin(1) * t710 - pkin(6) * t714) * t574, t813, t714 * t436 - t822, t710 * t817 + t854, t812, -t710 * t545 - t853, t811, -pkin(5) * t448 + t714 * t396 - t710 * t413, -pkin(5) * t457 + t714 * t402 - t710 * t422, -pkin(5) * t417 + t714 * t357 - t710 * t404, -pkin(5) * t394 + t714 * t360 - t710 * t377, t811, -t710 * t548 - t854, -t710 * t547 + t853, t813, t714 * t437 - t822, t812, -pkin(5) * t416 + t714 * t351 - t710 * t392, -pkin(5) * t449 + t714 * t358 - t710 * t398, -pkin(5) * t456 + t714 * t359 - t710 * t403, -pkin(5) * t365 + t714 * t342 - t710 * t354, t714 * t409 - t710 * t502, t714 * t380 - t710 * t446, t714 * t391 - t710 * t483, t714 * t408 - t710 * t500, t714 * t390 - t710 * t484, t714 * t428 - t710 * t525, -pkin(5) * t368 + t714 * t332 - t710 * t349, -pkin(5) * t370 + t714 * t333 - t710 * t350, -pkin(5) * t362 + t714 * t331 - t710 * t339, -pkin(5) * t334 + t714 * t327 - t710 * t329; 0, 0, 0, 0, 0, 0, g(3), 0, -g(1), 0, 0, 0, t677, 0, t678, 0, t656, t738, t628, pkin(5) * t628, t710 * t641 + t739, t710 * t617 + t714 * t680, t710 * t633 - t709 * t753, t710 * t640 - t739, t710 * t631 - t713 * t753, -t714 * qJDD(2) + t710 * t669, pkin(5) * t603 + t714 * t586 + t710 * t593, pkin(5) * t604 + t714 * t587 + t710 * t594, pkin(5) * t620 + t710 * t574, pkin(5) * t534 - (-pkin(1) * t714 - pkin(6) * t710) * t574, t815, t710 * t436 + t819, -t714 * t817 + t857, t814, t714 * t545 - t856, t816, pkin(5) * t450 + t710 * t396 + t714 * t413, pkin(5) * t459 + t710 * t402 + t714 * t422, pkin(5) * t419 + t710 * t357 + t714 * t404, pkin(5) * t395 + t710 * t360 + t714 * t377, t816, t714 * t548 - t857, t714 * t547 + t856, t815, t710 * t437 + t819, t814, pkin(5) * t418 + t710 * t351 + t714 * t392, pkin(5) * t451 + t710 * t358 + t714 * t398, pkin(5) * t458 + t710 * t359 + t714 * t403, pkin(5) * t366 + t710 * t342 + t714 * t354, t710 * t409 + t714 * t502, t710 * t380 + t714 * t446, t710 * t391 + t714 * t483, t710 * t408 + t714 * t500, t710 * t390 + t714 * t484, t710 * t428 + t714 * t525, pkin(5) * t369 + t710 * t332 + t714 * t349, pkin(5) * t371 + t710 * t333 + t714 * t350, pkin(5) * t363 + t710 * t331 + t714 * t339, pkin(5) * t335 + t710 * t327 + t714 * t329; 0, 0, 0, 0, 0, 0, -g(2), g(1), 0, 0, 0, 0, 0, 0, 0, qJDD(1), t683, t684, 0, 0, (t673 + t745) * t709, t713 * t672 + t709 * t675, t713 * t685 + t775, (t674 - t746) * t713, t709 * t687 + t762, 0, pkin(1) * t675 + pkin(6) * t632 + t763, -pkin(1) * t672 + pkin(6) * t634 - t777, pkin(1) * t679 + pkin(6) * t676 + t575, pkin(1) * t665 + pkin(6) * t575, t799, t713 * t488 + t709 * t492, -t848, t801, t847, t804, -pkin(1) * t544 + t713 * t461 + t709 * t508 - t845, pkin(1) * t742 + t713 * t466 + t709 * t515 + t859, pkin(6) * t438 + t713 * t405 + t709 * t410 - t826, pkin(1) * t592 + pkin(6) * t401 - pkin(7) * t778 + t713 * t439, t804, t848, -t847, t799, t713 * t489 + t709 * t493, t801, pkin(6) * t435 + t713 * t381 + t709 * t386 - t826, pkin(1) * t546 + t713 * t397 + t709 * t406 + t845, pkin(1) * t720 + t713 * t399 + t709 * t407 - t859, pkin(1) * t460 + pkin(6) * t376 + t713 * t361 + t709 * t372, t713 * t463 + t709 * t465, t713 * t420 + t709 * t421, t713 * t430 + t709 * t432, t713 * t462 + t709 * t464, t713 * t429 + t709 * t431, t713 * t485 + t709 * t486, -pkin(1) * t469 + pkin(6) * t383 + t713 * t344 + t709 * t346, -pkin(1) * t472 + pkin(6) * t385 + t713 * t345 + t709 * t347, -pkin(1) * t447 + pkin(6) * t379 + t713 * t336 + t709 * t340, -pkin(1) * t356 + pkin(6) * t338 + t713 * t328 + t709 * t330;];
tauB_reg = t1;
