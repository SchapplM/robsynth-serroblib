% Calculate minimal parameter regressor of Coriolis joint torque vector for
% S6RRRPRP5
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d3,d5,theta4]';
% MDP [30x1]
%   Minimal dynamic parameter vector (fixed base model)
%   see S6RRRPRP5_convert_par2_MPV_fixb.m
% 
% Output:
% tauc [6x1]
%   joint torques required to compensate Coriolis and centrifugal load

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 16:52
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauc = S6RRRPRP5_coriolisvecJ_fixb_mdp_slag_vp(qJ, qJD, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(10,1),zeros(30,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRPRP5_coriolisvecJ_fixb_mdp_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRRPRP5_coriolisvecJ_fixb_mdp_slag_vp: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RRRPRP5_coriolisvecJ_fixb_mdp_slag_vp: pkin has to be [10x1] (double)');
assert(isreal(MDP) && all(size(MDP) == [30 1]), ...
  'S6RRRPRP5_coriolisvecJ_fixb_mdp_slag_vp: MDP has to be [30x1] (double)'); 

%% Symbolic Calculation
% From coriolisvec_joint_fixb_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 16:51:50
% EndTime: 2019-03-09 16:52:06
% DurationCPUTime: 10.90s
% Computational Cost: add. (10591->551), mult. (25961->710), div. (0->0), fcn. (18551->8), ass. (0->245)
t605 = sin(qJ(3));
t606 = sin(qJ(2));
t694 = qJD(1) * t606;
t673 = t605 * t694;
t607 = cos(qJ(3));
t681 = t607 * qJD(2);
t566 = t673 - t681;
t690 = qJD(2) * t605;
t568 = t607 * t694 + t690;
t603 = sin(pkin(10));
t732 = cos(pkin(10));
t515 = -t603 * t566 + t568 * t732;
t604 = sin(qJ(5));
t678 = qJD(2) * qJD(3);
t685 = qJD(3) * t607;
t670 = t606 * t685;
t608 = cos(qJ(2));
t688 = qJD(2) * t608;
t752 = t605 * t688 + t670;
t527 = qJD(1) * t752 + t605 * t678;
t686 = qJD(3) * t606;
t671 = t605 * t686;
t637 = t608 * t681 - t671;
t618 = qJD(1) * t637 + t607 * t678;
t616 = -t603 * t527 + t618 * t732;
t617 = -t527 * t732 - t603 * t618;
t643 = t566 * t732 + t603 * t568;
t738 = cos(qJ(5));
t624 = t738 * t643;
t683 = qJD(5) * t604;
t414 = qJD(5) * t624 + t515 * t683 - t604 * t617 - t738 * t616;
t464 = t515 * t604 + t624;
t693 = qJD(1) * t608;
t591 = -qJD(3) + t693;
t580 = -qJD(5) + t591;
t726 = t464 * t580;
t397 = -t414 - t726;
t751 = t515 * t738 - t604 * t643;
t415 = qJD(5) * t751 + t604 * t616 - t738 * t617;
t664 = MDP(24) * t694;
t727 = t464 ^ 2;
t728 = t751 * t580;
t762 = t751 ^ 2;
t771 = t464 * t751;
t772 = t397 * MDP(22) + MDP(20) * t771 + qJD(2) * t664 + (-t415 - t728) * MDP(23) + (-t727 + t762) * MDP(21);
t719 = t607 * t608;
t650 = pkin(3) * t606 - qJ(4) * t719;
t734 = -qJ(4) - pkin(8);
t663 = qJD(3) * t734;
t652 = pkin(2) * t606 - pkin(8) * t608;
t569 = t652 * qJD(1);
t701 = pkin(7) * t673 + t607 * t569;
t770 = -qJD(1) * t650 - qJD(4) * t605 + t607 * t663 - t701;
t553 = t605 * t569;
t684 = qJD(4) * t607;
t721 = t606 * t607;
t769 = t553 + (-qJ(4) * t605 * t608 - pkin(7) * t721) * qJD(1) - t605 * t663 - t684;
t733 = qJD(2) * pkin(2);
t578 = pkin(7) * t694 - t733;
t525 = t566 * pkin(3) + qJD(4) + t578;
t473 = pkin(4) * t643 + t525;
t408 = t464 * pkin(5) - qJ(6) * t751 + t473;
t767 = t408 * t464;
t766 = t464 * t473;
t560 = t603 * t607 + t605 * t732;
t631 = t608 * t560;
t532 = qJD(1) * t631;
t628 = qJD(3) * t560;
t757 = t532 - t628;
t642 = t603 * t605 - t607 * t732;
t632 = t608 * t642;
t533 = qJD(1) * t632;
t741 = qJD(3) * t642;
t764 = t533 - t741;
t422 = pkin(5) * t751 + qJ(6) * t464;
t743 = t769 * t603 + t770 * t732;
t742 = t770 * t603 - t769 * t732;
t730 = t408 * t751;
t760 = t473 * t751;
t759 = pkin(4) * t694 + pkin(9) * t764 - t743;
t758 = pkin(9) * t757 + t742;
t756 = pkin(9) * t515;
t755 = pkin(9) * t643;
t623 = t738 * t642;
t709 = t560 * t683 - (-qJD(3) - qJD(5)) * t623 - t533 * t738 - t757 * t604;
t634 = t604 * t642;
t669 = qJD(5) * t738;
t674 = t738 * t560;
t708 = qJD(3) * t674 - qJD(5) * t634 - t532 * t738 + t560 * t669 + t604 * t764;
t596 = pkin(7) * t693;
t687 = qJD(3) * t605;
t753 = -t596 + (-t605 * t693 + t687) * pkin(3);
t679 = qJD(1) * qJD(2);
t750 = -0.2e1 * t679;
t749 = MDP(4) * t606;
t601 = t606 ^ 2;
t748 = MDP(5) * (-t608 ^ 2 + t601);
t573 = -pkin(2) * t608 - pkin(8) * t606 - pkin(1);
t557 = t573 * qJD(1);
t579 = qJD(2) * pkin(8) + t596;
t520 = t607 * t557 - t579 * t605;
t489 = -qJ(4) * t568 + t520;
t723 = t605 * t557;
t521 = t607 * t579 + t723;
t490 = -qJ(4) * t566 + t521;
t662 = t732 * t490;
t444 = -t489 * t603 - t662;
t430 = t444 + t755;
t485 = t603 * t490;
t445 = t732 * t489 - t485;
t431 = t445 - t756;
t594 = pkin(3) * t732 + pkin(4);
t737 = pkin(3) * t603;
t677 = t604 * t737;
t747 = -qJD(5) * t677 - t604 * t430 - t431 * t738 + t594 * t669;
t575 = t734 * t605;
t576 = t734 * t607;
t523 = t732 * t575 + t576 * t603;
t500 = -pkin(9) * t560 + t523;
t524 = t603 * t575 - t732 * t576;
t501 = -pkin(9) * t642 + t524;
t645 = t500 * t738 - t604 * t501;
t746 = qJD(5) * t645 - t604 * t759 + t738 * t758;
t453 = t604 * t500 + t501 * t738;
t745 = qJD(5) * t453 + t604 * t758 + t738 * t759;
t562 = t607 * t573;
t736 = pkin(7) * t605;
t517 = -qJ(4) * t721 + t562 + (-pkin(3) - t736) * t608;
t593 = pkin(7) * t719;
t699 = t605 * t573 + t593;
t722 = t605 * t606;
t522 = -qJ(4) * t722 + t699;
t470 = t732 * t517 - t522 * t603;
t544 = t642 * t606;
t448 = -pkin(4) * t608 + pkin(9) * t544 + t470;
t471 = t603 * t517 + t732 * t522;
t633 = t606 * t560;
t454 = -pkin(9) * t633 + t471;
t744 = t604 * t448 + t738 * t454;
t707 = -pkin(4) * t757 + t753;
t700 = t604 * t594 + t738 * t737;
t682 = t601 * qJD(1);
t635 = t650 * qJD(2);
t570 = t652 * qJD(2);
t689 = qJD(2) * t606;
t702 = t607 * t570 + t689 * t736;
t459 = -t606 * t684 + t635 + (-t593 + (qJ(4) * t606 - t573) * t605) * qJD(3) + t702;
t703 = t605 * t570 + t573 * t685;
t469 = (-pkin(7) * qJD(2) - qJ(4) * qJD(3)) * t721 + (-qJD(4) * t606 + (-pkin(7) * qJD(3) - qJ(4) * qJD(2)) * t608) * t605 + t703;
t420 = t732 * t459 - t603 * t469;
t622 = t606 * t628;
t410 = pkin(9) * t622 + (t606 * pkin(4) + pkin(9) * t632) * qJD(2) + t420;
t421 = t603 * t459 + t732 * t469;
t614 = -qJD(2) * t631 + t606 * t741;
t413 = pkin(9) * t614 + t421;
t740 = -qJD(5) * t744 + t410 * t738 - t604 * t413;
t735 = pkin(8) * t591;
t725 = t568 * t591;
t724 = t578 * t605;
t609 = qJD(2) ^ 2;
t720 = t606 * t609;
t718 = t608 * t591;
t717 = t608 * t609;
t610 = qJD(1) ^ 2;
t716 = t608 * t610;
t715 = -qJD(6) - t747;
t714 = qJ(6) * t694 - t746;
t713 = pkin(5) * t694 + t745;
t511 = t674 - t634;
t712 = pkin(5) * t708 + qJ(6) * t709 - t511 * qJD(6) + t707;
t558 = qJD(1) * t570;
t668 = t606 * t679;
t655 = pkin(7) * t668;
t704 = t607 * t558 + t605 * t655;
t436 = qJD(1) * t635 - qJD(3) * t490 - t568 * qJD(4) + t704;
t705 = t557 * t685 + t605 * t558;
t619 = -t579 * t687 - t607 * t655 + t705;
t442 = -qJ(4) * t527 - qJD(4) * t566 + t619;
t404 = t603 * t436 + t732 * t442;
t481 = -pkin(3) * t591 + t489;
t438 = t603 * t481 + t662;
t706 = qJD(5) * t700 + t430 * t738 - t604 * t431;
t697 = pkin(3) * t722 + t606 * pkin(7);
t692 = qJD(2) * t645;
t691 = qJD(2) * t453;
t437 = t732 * t481 - t485;
t423 = -pkin(4) * t591 + t437 - t756;
t426 = t438 - t755;
t390 = t423 * t738 - t604 * t426;
t680 = qJD(6) - t390;
t676 = pkin(3) * t752 + pkin(7) * t688;
t675 = -t607 * pkin(3) - pkin(2);
t667 = t608 * t679;
t666 = MDP(15) * t694;
t512 = t527 * pkin(3) + pkin(7) * t667;
t661 = pkin(1) * t750;
t403 = t732 * t436 - t603 * t442;
t660 = t591 + t693;
t659 = t566 + t681;
t658 = -t568 + t690;
t657 = qJD(3) + t693;
t656 = pkin(5) * t668;
t394 = pkin(4) * t668 - pkin(9) * t616 + t403;
t398 = pkin(9) * t617 + t404;
t654 = -t604 * t394 - t738 * t398 - t423 * t669 + t426 * t683;
t653 = -t738 * t394 + t604 * t398 + t423 * t683 + t426 * t669;
t482 = pkin(3) * t568 + pkin(4) * t515;
t572 = t580 * qJD(6);
t588 = qJ(6) * t668;
t382 = t588 - t572 - t654;
t648 = t448 * t738 - t604 * t454;
t391 = t604 * t423 + t426 * t738;
t644 = t657 * t690;
t641 = -t390 * t580 + t654;
t640 = -t391 * t580 - t653;
t639 = t604 * t410 + t738 * t413 + t448 * t669 - t454 * t683;
t638 = t594 * t738 - t677;
t630 = t580 * t706 - t653;
t383 = t653 - t656;
t626 = t604 * t633;
t620 = t606 * t674;
t534 = pkin(4) * t642 + t675;
t518 = pkin(4) * t633 + t697;
t615 = -qJD(2) * t632 - t622;
t472 = -pkin(4) * t614 + t676;
t450 = -pkin(4) * t617 + t512;
t386 = t415 * pkin(5) + t414 * qJ(6) - qJD(6) * t751 + t450;
t545 = -pkin(5) - t638;
t543 = qJ(6) + t700;
t510 = t560 * t604 + t623;
t495 = -t544 * t738 - t626;
t494 = -t544 * t604 + t620;
t449 = t510 * pkin(5) - t511 * qJ(6) + t534;
t432 = t494 * pkin(5) - t495 * qJ(6) + t518;
t429 = -qJD(5) * t626 - t544 * t669 + t604 * t615 - t614 * t738;
t428 = qJD(5) * t620 - t544 * t683 - t604 * t614 - t615 * t738;
t417 = t608 * pkin(5) - t648;
t416 = -qJ(6) * t608 + t744;
t412 = t422 + t482;
t389 = -t580 * qJ(6) + t391;
t388 = t580 * pkin(5) + t680;
t387 = t429 * pkin(5) + t428 * qJ(6) - t495 * qJD(6) + t472;
t385 = -pkin(5) * t689 - t740;
t384 = qJ(6) * t689 - qJD(6) * t608 + t639;
t1 = [(-t421 * t643 + t471 * t617 - t404 * t633 + t438 * t614 - t420 * t515 - t470 * t616 + t403 * t544 + t437 * (t560 * t686 + t642 * t688)) * MDP(18) + (t415 * t608 + t429 * t580) * MDP(23) + (t414 * t608 + t428 * t580) * MDP(22) + (t383 * t608 + t385 * t580 + t386 * t494 + t387 * t464 + t408 * t429 + t415 * t432) * MDP(27) + (-t518 * t414 - t473 * t428 + t450 * t495 + t472 * t751 + t639 * t580 - t654 * t608) * MDP(26) + (-t382 * t608 - t384 * t580 - t386 * t495 - t387 * t751 + t408 * t428 + t414 * t432) * MDP(29) + (t414 * t494 - t415 * t495 + t428 * t464 - t429 * t751) * MDP(21) + (-t382 * t494 + t383 * t495 - t384 * t464 + t385 * t751 - t388 * t428 - t389 * t429 - t414 * t417 - t415 * t416) * MDP(28) + (-t414 * t495 - t428 * t751) * MDP(20) + ((-qJD(1) * t744 - t391) * MDP(26) + (qJD(1) * t495 + t751) * MDP(22) + (-qJD(1) * t494 - t464) * MDP(23) + (qJD(1) * t416 + t389) * MDP(29) + (-qJD(1) * t417 - t388) * MDP(27) + t390 * MDP(25) - t660 * MDP(15) + (-t580 - t693) * MDP(24)) * t689 + (t703 * t591 + t705 * t608 + (-t578 * t606 - t579 * t608 + (-t718 - t682) * pkin(7)) * t687 + ((pkin(7) * t568 + t578 * t607) * t608 + (-t699 * qJD(1) - t521 + (-t591 + t657) * pkin(7) * t607) * t606) * qJD(2)) * MDP(17) + (-(-t573 * t687 + t702) * t591 + (t578 * t685 + pkin(7) * t527 + (t562 * qJD(1) + t520) * qJD(2)) * t606 + ((pkin(7) * t566 + t724) * qJD(2) + (t723 + (pkin(7) * t591 + t579) * t607) * qJD(3) - t704) * t608) * MDP(16) - MDP(7) * t720 + (pkin(7) * t720 + t608 * t661) * MDP(10) + (t568 * t637 + t618 * t721) * MDP(11) + (-pkin(7) * t717 + t606 * t661) * MDP(9) + (t591 * t670 + t527 * t608 + (-t566 * t606 + (-t682 + t718) * t605) * qJD(2)) * MDP(14) + (t518 * t415 + t473 * t429 + t450 * t494 + t472 * t464 - t580 * t740 + t653 * t608 + t648 * t668) * MDP(25) + (t403 * t470 + t404 * t471 + t437 * t420 + t438 * t421 + t512 * t697 + t525 * t676) * MDP(19) + ((-t607 * t566 - t568 * t605) * t688 + ((t566 + t673) * t687 + (-t568 * qJD(3) - t527 - t644) * t607) * t606) * MDP(12) + (t660 * t671 + (t568 * t606 + (t682 + (-t591 - t657) * t608) * t607) * qJD(2)) * MDP(13) + 0.2e1 * t667 * t749 + MDP(6) * t717 + (t382 * t416 + t383 * t417 + t384 * t389 + t385 * t388 + t386 * t432 + t387 * t408) * MDP(30) + t748 * t750; t610 * t748 + (-t605 ^ 2 * qJD(1) * t686 + (t644 - t725) * t607) * MDP(11) + ((-t527 + t725) * t605 + ((-t566 + t681) * qJD(3) + (t608 * t659 - t671) * qJD(1)) * t607) * MDP(12) + (-t591 * t685 + (t606 * t658 + t607 * t718) * qJD(1)) * MDP(13) + (t591 * t687 + (-t718 * t605 + t606 * t659) * qJD(1)) * MDP(14) + (-pkin(2) * t527 + t701 * t591 + (t607 * t735 + t724) * qJD(3) + ((-pkin(8) * t690 - t520) * t606 + (-pkin(7) * t659 - t724) * t608) * qJD(1)) * MDP(16) + (-t553 * t591 + (-t605 * t735 + (t578 - t733) * t607) * qJD(3) + ((-t578 - t733) * t719 + (pkin(2) * t687 - pkin(8) * t681 + t521) * t606 + (t591 * t721 + t608 * t658) * pkin(7)) * qJD(1)) * MDP(17) + (-t403 * t560 - t404 * t642 - t437 * t764 + t757 * t438 - t743 * t515 - t523 * t616 + t524 * t617 - t742 * t643) * MDP(18) + (t403 * t523 + t404 * t524 + t743 * t437 + t742 * t438 + t512 * t675 + t525 * t753) * MDP(19) + (-t414 * t511 - t709 * t751) * MDP(20) + (t414 * t510 - t415 * t511 + t464 * t709 - t708 * t751) * MDP(21) + (t709 * t580 + (qJD(2) * t511 - t751) * t694) * MDP(22) + (t708 * t580 + (-qJD(2) * t510 + t464) * t694) * MDP(23) + (t534 * t415 + t450 * t510 + t745 * t580 + t708 * t473 + t707 * t464 + (-t390 + t692) * t694) * MDP(25) + (-t534 * t414 + t450 * t511 + t746 * t580 - t709 * t473 + t707 * t751 + (t391 - t691) * t694) * MDP(26) + (t386 * t510 + t415 * t449 + t713 * t580 + t712 * t464 + t708 * t408 + (t388 + t692) * t694) * MDP(27) + (-t382 * t510 + t383 * t511 - t388 * t709 - t389 * t708 + t414 * t645 - t415 * t453 + t464 * t714 + t713 * t751) * MDP(28) + (-t386 * t511 + t414 * t449 + t714 * t580 - t712 * t751 + t709 * t408 + (-t389 + t691) * t694) * MDP(29) + (t382 * t453 - t383 * t645 + t386 * t449 + t388 * t713 - t389 * t714 + t408 * t712) * MDP(30) - t716 * t749 + t591 * t666 + t580 * t664 + (MDP(9) * t606 * t610 + MDP(10) * t716) * pkin(1); t568 * t566 * MDP(11) + (-t566 ^ 2 + t568 ^ 2) * MDP(12) + (-t566 * t591 + t618) * MDP(13) + (-t527 - t725) * MDP(14) + qJD(2) * t666 + (-t568 * t578 + t704 + (-t591 - qJD(3)) * t521) * MDP(16) + (-t520 * t591 + t566 * t578 - t619) * MDP(17) + ((-t603 ^ 2 - t732 ^ 2) * pkin(3) * t618 + (-t437 + t445) * t643 + (t438 + t444) * t515) * MDP(18) + (-t437 * t444 - t438 * t445 + (t403 * t732 + t404 * t603 - t525 * t568) * pkin(3)) * MDP(19) + (-t482 * t464 + t638 * t668 + t630 - t760) * MDP(25) + (-t482 * t751 + t580 * t747 - t700 * t668 + t654 + t766) * MDP(26) + (-t730 - t412 * t464 + (pkin(5) - t545) * t668 + t630) * MDP(27) + (-t414 * t545 - t415 * t543 + (t389 + t706) * t751 + (t388 + t715) * t464) * MDP(28) + (t412 * t751 + t543 * t668 + t580 * t715 + t382 - t767) * MDP(29) + (t382 * t543 + t383 * t545 + t388 * t706 - t389 * t715 - t408 * t412) * MDP(30) + t772; (-t515 ^ 2 - t643 ^ 2) * MDP(18) + (t437 * t515 + t438 * t643 + t512) * MDP(19) + (-t727 - t762) * MDP(28) + (-t388 * t751 + t389 * t464 + t386) * MDP(30) + (-MDP(26) + MDP(29)) * (t414 - t726) + (MDP(25) + MDP(27)) * (t415 - t728); (t640 - t760) * MDP(25) + (t641 + t766) * MDP(26) + (-t422 * t464 + t640 + 0.2e1 * t656 - t730) * MDP(27) + (pkin(5) * t414 - qJ(6) * t415 + (t389 - t391) * t751 + (t388 - t680) * t464) * MDP(28) + (t422 * t751 - 0.2e1 * t572 + 0.2e1 * t588 - t641 - t767) * MDP(29) + (-pkin(5) * t383 + qJ(6) * t382 - t388 * t391 + t389 * t680 - t408 * t422) * MDP(30) + t772; (-t668 + t771) * MDP(27) + t397 * MDP(28) + (-t580 ^ 2 - t762) * MDP(29) + (t389 * t580 + t383 + t730) * MDP(30);];
tauc  = t1;
