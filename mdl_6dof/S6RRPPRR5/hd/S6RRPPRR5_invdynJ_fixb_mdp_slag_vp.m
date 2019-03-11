% Calculate vector of inverse dynamics joint torques for
% S6RRPPRR5
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% qJDD [6x1]
%   Generalized joint accelerations
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d1,d2,d5,d6]';
% MDP [32x1]
%   Minimal dynamic parameter vector (fixed base model)
%   see S6RRPPRR5_convert_par2_MPV_fixb.m
% 
% Output:
% tau [6x1]
%   joint torques of inverse dynamics (contains inertial, gravitational Coriolis and centrifugal forces)

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 09:12
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tau = S6RRPPRR5_invdynJ_fixb_mdp_slag_vp(qJ, qJD, qJDD, g, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(3,1),zeros(10,1),zeros(32,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPPRR5_invdynJ_fixb_mdp_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRPPRR5_invdynJ_fixb_mdp_slag_vp: qJD has to be [6x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [6 1]), ...
  'S6RRPPRR5_invdynJ_fixb_mdp_slag_vp: qJDD has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RRPPRR5_invdynJ_fixb_mdp_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RRPPRR5_invdynJ_fixb_mdp_slag_vp: pkin has to be [10x1] (double)');
assert(isreal(MDP) && all(size(MDP) == [32 1]), ...
  'S6RRPPRR5_invdynJ_fixb_mdp_slag_vp: MDP has to be [32x1] (double)'); 

%% Symbolic Calculation
% From invdyn_joint_fixb_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 09:11:19
% EndTime: 2019-03-09 09:11:33
% DurationCPUTime: 11.65s
% Computational Cost: add. (4583->665), mult. (11089->843), div. (0->0), fcn. (8249->10), ass. (0->272)
t601 = sin(pkin(6));
t610 = cos(qJ(2));
t707 = qJD(1) * qJD(2);
t689 = t610 * t707;
t606 = sin(qJ(2));
t705 = qJDD(1) * t606;
t639 = t689 + t705;
t779 = t601 * t639;
t602 = cos(pkin(6));
t711 = t602 * qJD(1);
t577 = qJD(2) + t711;
t609 = cos(qJ(5));
t605 = sin(qJ(5));
t721 = qJD(1) * t601;
t697 = t606 * t721;
t670 = t605 * t697;
t775 = t609 * t577 + t670;
t488 = qJD(6) + t775;
t780 = t488 ^ 2;
t494 = -t577 * t605 + t609 * t697;
t696 = t610 * t721;
t536 = qJD(5) + t696;
t755 = t494 * t536;
t756 = t775 * t536;
t706 = qJDD(1) * t602;
t574 = qJDD(2) + t706;
t773 = t574 * qJ(3) + t577 * qJD(3);
t607 = sin(qJ(1));
t737 = t607 * t610;
t611 = cos(qJ(1));
t738 = t606 * t611;
t514 = t602 * t738 + t737;
t741 = t601 * t611;
t473 = t514 * t609 + t605 * t741;
t733 = t610 * t611;
t739 = t606 * t607;
t513 = -t602 * t733 + t739;
t604 = sin(qJ(6));
t608 = cos(qJ(6));
t778 = t473 * t604 + t513 * t608;
t777 = t473 * t608 - t513 * t604;
t719 = qJD(2) * t610;
t693 = t605 * t719;
t714 = qJD(5) * t609;
t632 = t606 * t714 + t693;
t715 = qJD(5) * t605;
t657 = -t609 * t574 + t577 * t715;
t687 = t605 * t705;
t433 = t601 * (qJD(1) * t632 + t687) - t657;
t592 = t601 * pkin(1);
t742 = t601 * t610;
t745 = t601 * t606;
t726 = pkin(2) * t742 + qJ(3) * t745;
t498 = -t592 - t726;
t580 = pkin(3) * t742;
t479 = t580 - t498;
t667 = pkin(4) * t610 - pkin(9) * t606;
t645 = t667 * t601;
t453 = t645 + t479;
t765 = pkin(1) * t606;
t642 = pkin(8) * t742 + t602 * t765;
t497 = t602 * qJ(3) + t642;
t478 = -qJ(4) * t742 + t497;
t465 = -pkin(9) * t602 + t478;
t776 = t605 * t453 + t609 * t465;
t505 = t642 * qJD(2);
t571 = qJD(4) * t745;
t694 = t601 * t719;
t462 = -qJ(4) * t694 + t505 - t571;
t774 = 0.2e1 * t773;
t515 = t602 * t737 + t738;
t630 = -g(1) * t515 - g(2) * t513 + g(3) * t742;
t703 = pkin(1) * t711;
t503 = -pkin(8) * t697 + t610 * t703;
t709 = qJD(3) - t503;
t544 = t574 * pkin(2);
t772 = t544 - qJDD(3);
t704 = qJDD(1) * t610;
t570 = t601 * t704;
t690 = t606 * t707;
t669 = t601 * t690;
t500 = -qJDD(5) - t570 + t669;
t612 = -pkin(2) - pkin(3);
t598 = pkin(4) - t612;
t643 = -pkin(9) * t610 - t598 * t606;
t627 = qJD(2) * t643;
t572 = qJD(3) * t745;
t656 = pkin(2) * t570 + qJ(3) * t779 + qJD(1) * t572 + qJDD(1) * t592;
t633 = pkin(3) * t570 + qJDD(4) + t656;
t412 = (qJD(1) * t627 + qJDD(1) * t667) * t601 + t633;
t485 = -pkin(1) * t721 - pkin(2) * t696 - qJ(3) * t697;
t463 = pkin(3) * t696 + qJD(4) - t485;
t442 = qJD(1) * t645 + t463;
t504 = pkin(8) * t696 + t606 * t703;
t483 = -qJ(4) * t696 + t504;
t545 = t577 * qJ(3);
t459 = t545 + t483;
t450 = -pkin(9) * t577 + t459;
t414 = t442 * t605 + t450 * t609;
t528 = qJ(4) * t669;
t691 = qJ(4) * t704;
t766 = pkin(1) * t602;
t702 = qJD(2) * t766;
t674 = qJD(1) * t702;
t701 = pkin(1) * t706;
t698 = pkin(8) * t570 + t606 * t701 + t610 * t674;
t718 = qJD(4) * t610;
t720 = qJD(2) * t606;
t425 = t528 + (-t691 + (-pkin(8) * t720 - t718) * qJD(1)) * t601 + t698 + t773;
t419 = -pkin(9) * t574 + t425;
t619 = -qJD(5) * t414 + t609 * t412 - t419 * t605;
t399 = pkin(5) * t500 - t619;
t516 = -t602 * t739 + t733;
t743 = t601 * t609;
t475 = -t516 * t605 - t607 * t743;
t511 = t602 * t609 + t605 * t745;
t734 = t609 * t611;
t751 = t514 * t605;
t636 = g(1) * t475 + g(2) * (t601 * t734 - t751) - g(3) * t511;
t771 = t488 * (pkin(5) * t494 + pkin(10) * t488) + t399 + t636;
t431 = qJDD(6) + t433;
t522 = pkin(5) * t609 + pkin(10) * t605 + t598;
t770 = t488 * (-t483 + t536 * (pkin(5) * t605 - pkin(10) * t609)) - t522 * t431;
t728 = qJ(3) * t694 + t572;
t439 = t601 * t627 + t728;
t569 = t610 * t702;
t589 = t602 * qJD(3);
t454 = t569 + t589 + (-t718 + (-pkin(8) + qJ(4)) * t720) * t601;
t768 = -qJD(5) * t776 + t439 * t609 - t454 * t605;
t767 = -t574 * pkin(3) - t779 * qJ(4) - qJD(1) * t571;
t764 = g(3) * t606;
t763 = pkin(2) * qJD(2);
t762 = pkin(8) * qJD(2);
t761 = qJ(3) * t610;
t672 = t605 * t574 + t577 * t714 - t609 * t779;
t432 = -qJD(5) * t670 - t672;
t712 = qJD(6) * t608;
t699 = t608 * t432 - t604 * t500 + t536 * t712;
t713 = qJD(6) * t604;
t406 = -t494 * t713 + t699;
t760 = t406 * t604;
t754 = t494 * t604;
t455 = -t608 * t536 + t754;
t759 = t455 * t488;
t457 = t494 * t608 + t536 * t604;
t758 = t457 * t488;
t603 = qJ(3) - pkin(9);
t757 = t488 * t603;
t680 = t488 * t608;
t749 = t536 * t605;
t748 = t536 * t609;
t747 = t574 * MDP(8);
t597 = t601 ^ 2;
t746 = t597 * qJD(1) ^ 2;
t744 = t601 * t607;
t740 = t604 * t431;
t736 = t608 * t431;
t735 = t609 * t610;
t549 = qJ(3) * t696;
t451 = t643 * t721 + t549;
t546 = qJ(4) * t697;
t481 = t546 + t503;
t730 = t605 * t451 + t609 * t481;
t727 = -pkin(8) * t745 + t610 * t766;
t599 = t606 ^ 2;
t600 = t610 ^ 2;
t725 = t599 - t600;
t724 = MDP(12) * t601;
t723 = MDP(17) * t601;
t722 = qJ(3) * qJD(2);
t717 = qJD(5) * t536;
t716 = qJD(5) * t603;
t710 = qJD(3) - t481;
t700 = t610 * t746;
t499 = -t602 * pkin(2) - t727;
t695 = t601 * t720;
t692 = 0.2e1 * pkin(1) * t597;
t685 = -t513 * pkin(2) + qJ(3) * t514;
t684 = -t515 * pkin(2) + qJ(3) * t516;
t638 = t605 * t412 + t609 * t419 + t442 * t714 - t450 * t715;
t398 = -pkin(10) * t500 + t638;
t671 = pkin(8) * t779 + t606 * t674 - t610 * t701;
t444 = t671 - t772;
t424 = t444 + t767;
t418 = pkin(4) * t574 - t424;
t403 = pkin(5) * t433 - pkin(10) * t432 + t418;
t683 = -t604 * t398 + t608 * t403;
t681 = t432 * t604 + t608 * t500;
t679 = qJD(1) * t479 + t463;
t678 = qJD(1) * t498 + t485;
t677 = t577 + t711;
t676 = t574 + t706;
t675 = t612 * t745;
t673 = t606 * t700;
t477 = -t577 * pkin(2) + t709;
t664 = g(1) * t513 - g(2) * t515;
t663 = g(1) * t516 + g(2) * t514;
t662 = g(1) * t514 - g(2) * t516;
t661 = g(1) * t611 + g(2) * t607;
t660 = g(1) * t607 - g(2) * t611;
t486 = (t604 * t735 + t606 * t608) * t721;
t659 = t604 * t714 + t486;
t487 = (-t604 * t606 + t608 * t735) * t721;
t658 = -t608 * t714 - t487;
t466 = -t602 * pkin(3) - qJ(4) * t745 + t499;
t655 = qJD(2) * t675;
t654 = t608 * t398 + t604 * t403;
t411 = pkin(10) * t536 + t414;
t448 = -t577 * pkin(3) + t477 - t546;
t437 = pkin(4) * t577 - t448;
t415 = pkin(5) * t775 - pkin(10) * t494 + t437;
t401 = t411 * t608 + t415 * t604;
t653 = t411 * t604 - t415 * t608;
t421 = pkin(10) * t742 + t776;
t458 = t602 * pkin(4) - t466;
t512 = -t602 * t605 + t606 * t743;
t429 = pkin(5) * t511 - pkin(10) * t512 + t458;
t652 = t421 * t608 + t429 * t604;
t651 = -t421 * t604 + t429 * t608;
t413 = t442 * t609 - t450 * t605;
t649 = t453 * t609 - t465 * t605;
t647 = -pkin(8) * t695 + t569;
t646 = t611 * pkin(1) + t516 * pkin(2) + pkin(8) * t744 + qJ(3) * t515;
t644 = -t512 * t604 + t608 * t742;
t471 = t512 * t608 + t604 * t742;
t426 = qJD(1) * t655 + t633;
t461 = t655 + t728;
t641 = qJD(1) * t461 + qJDD(1) * t479 + t426;
t438 = pkin(2) * t669 - t656;
t484 = pkin(2) * t695 - t728;
t640 = -qJD(1) * t484 - qJDD(1) * t498 - t438;
t637 = t605 * t439 + t453 * t714 + t609 * t454 - t465 * t715;
t635 = -pkin(1) * t607 - t514 * pkin(2) + pkin(8) * t741 - qJ(3) * t513;
t634 = -t505 * t577 + t662;
t631 = t663 - t698;
t410 = -pkin(5) * t536 - t413;
t628 = -qJD(3) * t488 - qJD(5) * t410 - t431 * t603;
t626 = t605 * t713 + t658;
t625 = -pkin(8) * t669 + t698;
t624 = t418 - t630;
t623 = -pkin(10) * t431 + (t410 + t413) * t488;
t622 = -t437 * t536 + t603 * t500;
t621 = -t671 - t630;
t620 = t503 * t577 + t631;
t618 = qJD(6) * t757 + t630;
t617 = t621 + t772;
t616 = t504 * t577 + t621;
t615 = g(3) * t745 + (-pkin(10) * t697 - qJD(6) * t522 + t730) * t488 + t663;
t614 = t617 - t767;
t517 = t577 * t696;
t501 = pkin(2) * t697 - t549;
t489 = t589 + t647;
t482 = qJD(1) * t675 + t549;
t480 = t545 + t504;
t476 = t516 * t609 - t605 * t744;
t469 = -qJD(5) * t511 + t609 * t694;
t468 = t601 * t632 - t602 * t715;
t464 = -t517 + t779;
t441 = t476 * t608 - t515 * t604;
t440 = -t476 * t604 - t515 * t608;
t434 = t625 + t773;
t428 = qJD(6) * t644 + t469 * t608 - t604 * t695;
t427 = qJD(6) * t471 + t469 * t604 + t608 * t695;
t422 = pkin(5) * t697 - t451 * t609 + t481 * t605;
t420 = -pkin(5) * t742 - t649;
t417 = pkin(5) * t468 - pkin(10) * t469 - t462;
t407 = qJD(6) * t457 + t681;
t405 = pkin(5) * t695 - t768;
t404 = -pkin(10) * t695 + t637;
t397 = -t401 * qJD(6) + t683;
t396 = -qJD(6) * t653 + t654;
t1 = [(-t468 * t536 + t500 * t511) * MDP(22) + (t469 * t536 - t500 * t512) * MDP(21) + ((-qJD(6) * t652 - t404 * t604 + t417 * t608) * t488 + t651 * t431 + t397 * t511 - t653 * t468 + t405 * t455 + t420 * t407 - t399 * t644 + t410 * t427 + g(1) * t777 - g(2) * t441) * MDP(31) + (-(qJD(6) * t651 + t404 * t608 + t417 * t604) * t488 - t652 * t431 - t396 * t511 - t401 * t468 + t405 * t457 + t420 * t406 + t399 * t471 + t410 * t428 - g(1) * t778 - g(2) * t440) * MDP(32) + (t406 * t471 + t428 * t457) * MDP(26) + (t431 * t511 + t468 * t488) * MDP(30) + (t406 * t511 + t428 * t488 + t431 * t471 + t457 * t468) * MDP(28) + (t432 * t512 + t469 * t494) * MDP(19) + (t406 * t644 - t407 * t471 - t427 * t457 - t428 * t455) * MDP(27) + (-t407 * t511 - t427 * t488 + t431 * t644 - t455 * t468) * MDP(29) + (-t432 * t511 - t433 * t512 - t468 * t494 - t469 * t775) * MDP(20) + (g(1) * t473 - g(2) * t476 + t418 * t511 + t458 * t433 + t437 * t468 - t462 * t775 - t649 * t500 + t536 * t768) * MDP(24) + ((t606 * t676 + t677 * t719) * MDP(6) + (t610 * t676 - t677 * t720) * MDP(7) + (-t500 * t610 - t536 * t720) * MDP(23) + (g(1) * t734 + t414 * t720 - t610 * t638) * MDP(25) + (t610 * t641 - t679 * t720) * MDP(15) + (-t413 * t720 + t610 * t619) * MDP(24) + (t432 * t610 - t494 * t720) * MDP(21) + (t610 * t640 + t678 * t720) * MDP(11) + (-t433 * t610 + t720 * t775) * MDP(22) + (t606 * t640 - t678 * t719) * MDP(13) + (t606 * t641 + t679 * t719) * MDP(16)) * t601 + (-g(1) * t751 - g(2) * t475 + t418 * t512 + t458 * t432 + t437 * t469 - t462 * t494 + t500 * t776 - t637 * t536) * MDP(25) + (-t574 * t642 - t577 * t647 - t602 * t625 - t639 * t692 - t664) * MDP(10) + (t727 * t574 - t671 * t602 + (-t690 + t704) * t692 + t634) * MDP(9) + t660 * MDP(2) + t661 * MDP(3) + qJDD(1) * MDP(1) + (t425 * t478 + t459 * t454 + t424 * t466 + t448 * t462 + t426 * t479 + t463 * t461 - g(1) * (-pkin(3) * t514 - qJ(4) * t741 + t635) - g(2) * (pkin(3) * t516 - qJ(4) * t744 + t646)) * MDP(18) + (-g(1) * t635 - g(2) * t646 + t434 * t497 + t438 * t498 + t444 * t499 + t477 * t505 + t480 * t489 + t485 * t484) * MDP(14) + (-t424 * t602 - t462 * t577 - t466 * t574 + t662) * MDP(15) + (-t444 * t602 - t499 * t574 + t634) * MDP(11) + (t434 * t602 + t489 * t577 + t497 * t574 + t664) * MDP(13) + (t425 * t602 + t454 * t577 + t478 * t574 + t664) * MDP(16) + ((qJDD(1) * t599 + 0.2e1 * t606 * t689) * MDP(4) + 0.2e1 * (t606 * t704 - t707 * t725) * MDP(5)) * t597 + ((-qJD(2) * t448 - qJDD(1) * t478 - t425 + (-qJD(2) * t466 - t454) * qJD(1)) * t610 + (qJD(2) * t459 - qJDD(1) * t466 - t424 + (qJD(2) * t478 - t462) * qJD(1)) * t606 + t661) * t723 + ((qJD(2) * t477 + qJDD(1) * t497 + t434 + (qJD(2) * t499 + t489) * qJD(1)) * t610 + (-qJD(2) * t480 + qJDD(1) * t499 + t444 + (-qJD(2) * t497 + t505) * qJD(1)) * t606 - t661) * t724 + t602 * t747; t464 * MDP(6) + t747 + (-t410 * t487 - t422 * t457 + t770 * t604 + t615 * t608 + (t457 * t716 + t604 * t618 + t608 * t628 - t396) * t609 + (t401 * t696 + t410 * t713 + qJD(3) * t457 - t399 * t608 + t603 * t406 + (t603 * t680 + t401) * qJD(5)) * t605) * MDP(32) + (t536 * t715 + t609 * t500 + (-t606 * t775 + t610 * t749) * t721) * MDP(22) + (t413 * t697 + t598 * t433 + t483 * t775 + (-t536 * t710 + t622) * t605 + ((-t451 - t716) * t536 + t624) * t609) * MDP(24) + (-t410 * t486 - t422 * t455 - t770 * t608 + t615 * t604 + (t455 * t716 + t604 * t628 - t608 * t618 + t397) * t609 + (t653 * t696 - t410 * t712 + qJD(3) * t455 - t399 * t604 + t603 * t407 + (t604 * t757 + t653) * qJD(5)) * t605) * MDP(31) + (t598 * t432 + t730 * t536 - t414 * t697 + t483 * t494 + (-qJD(3) * t536 + t622) * t609 + (t536 * t716 - t624) * t605) * MDP(25) + (t425 * qJ(3) + t424 * t612 - t448 * t483 - t463 * t482 - g(1) * (-pkin(3) * t515 + t684) - g(2) * (-pkin(3) * t513 + t685) - g(3) * (t580 + t726) + t710 * t459) * MDP(18) + (-t444 * pkin(2) - g(1) * t684 - g(2) * t685 - g(3) * t726 + t434 * qJ(3) - t477 * t504 + t480 * t709 - t485 * t501) * MDP(14) - MDP(4) * t673 + (t455 * t487 + t457 * t486 + (t455 * t608 + t457 * t604) * t714 + (t760 + t407 * t608 + (-t455 * t604 + t457 * t608) * qJD(6)) * t605) * MDP(27) + ((-t606 * t612 - t761) * qJDD(1) + ((-t459 + t483 + t722) * t606 + (-qJD(2) * t612 + t448 - t710) * t610) * qJD(1)) * t723 + ((-t432 + t756) * t609 + (t433 + t755) * t605) * MDP(20) + (-t432 * t605 - t494 * t748) * MDP(19) + (t431 * t609 - t488 * t749) * MDP(30) + (t406 * t609 + t658 * t488 + (-t457 * t536 + t488 * t713 - t736) * t605) * MDP(28) + (t746 * t765 + t616) * MDP(9) + (0.2e1 * t544 - qJDD(3) + (-t485 * t606 + t501 * t610) * t721 + t616) * MDP(11) + ((-pkin(2) * t606 + t761) * qJDD(1) + ((t480 - t504 - t722) * t606 + (-t477 + t709 - t763) * t610) * qJD(1)) * t724 + (t570 + (-qJD(2) + t577) * t697) * MDP(7) + (-t536 * t714 + t605 * t500 + (t494 * t606 - t536 * t735) * t721) * MDP(21) + (t614 + t483 * t577 - t574 * t612 + (t463 * t606 - t482 * t610) * t721) * MDP(15) + (pkin(1) * t700 + (pkin(8) * t707 + g(3)) * t745 + t620) * MDP(10) + (-t407 * t609 + t659 * t488 + (t455 * t536 + t488 * t712 + t740) * t605) * MDP(29) + (-t406 * t605 * t608 + t457 * t626) * MDP(26) + t536 * MDP(23) * t697 + (-t481 * t577 + t528 + (-t691 - t764 + ((-qJD(4) - t463) * t610 + (-t482 - t762) * t606) * qJD(1)) * t601 - t631 + t774) * MDP(16) + ((-t764 + (t485 * t610 + (t501 - t762) * t606) * qJD(1)) * t601 - t620 + t774) * MDP(13) + t725 * MDP(5) * t746; (-t480 * t577 - t617) * MDP(14) + (-t459 * t577 - t614) * MDP(18) + (t657 - t755) * MDP(24) + (t672 + t756) * MDP(25) + (t455 * t494 - t736) * MDP(31) + (t457 * t494 + t740) * MDP(32) + (-MDP(24) * t687 + (-MDP(24) * t693 + (t485 * MDP(14) - t463 * MDP(18) + (-MDP(24) * t609 + MDP(25) * t605) * qJD(5)) * t606) * qJD(1)) * t601 + (MDP(13) + MDP(16)) * (-t577 ^ 2 - t599 * t746) + (MDP(11) + MDP(15)) * (-t574 - t673) + (MDP(12) - MDP(17)) * t464 + (MDP(31) * t604 + MDP(32) * t608) * t780; t570 * MDP(15) + t517 * MDP(16) + (g(3) * t602 + t633) * MDP(18) + (-MDP(24) * t500 - MDP(25) * t717 - MDP(31) * t407 - MDP(32) * t406) * t609 + (-MDP(24) * t717 + t500 * MDP(25) + (qJD(5) * t455 - t740) * MDP(31) + (qJD(5) * t457 - t736) * MDP(32)) * t605 + (-t599 - t600) * MDP(17) * t746 + ((-t605 * t712 - t659) * MDP(31) + t626 * MDP(32)) * t488 + (MDP(16) * t705 + t660 * MDP(18) + ((-MDP(25) * t748 + qJD(2) * MDP(16) + t459 * MDP(18) + (-MDP(24) * t536 + MDP(31) * t455 + MDP(32) * t457) * t605) * t610 + ((-qJD(2) - t577) * MDP(15) + (-pkin(3) * qJD(2) + t448 - t763) * MDP(18) - t775 * MDP(24) - t494 * MDP(25)) * t606) * qJD(1)) * t601; -t775 ^ 2 * MDP(20) + (t432 + t756) * MDP(21) + (-t433 + t755) * MDP(22) - t500 * MDP(23) + (t414 * t536 + t619 - t636) * MDP(24) + (g(1) * t476 + g(2) * t473 + g(3) * t512 + t413 * t536 + t437 * t775 - t638) * MDP(25) + (t457 * t680 + t760) * MDP(26) + ((t406 - t759) * t608 + (-t407 - t758) * t604) * MDP(27) + (t488 * t680 + t740) * MDP(28) + (-t604 * t780 + t736) * MDP(29) + (-pkin(5) * t407 - t414 * t455 + t623 * t604 - t608 * t771) * MDP(31) + (-pkin(5) * t406 - t414 * t457 + t604 * t771 + t623 * t608) * MDP(32) + (MDP(19) * t775 + t494 * MDP(20) - t437 * MDP(24) - t457 * MDP(28) + t455 * MDP(29) - t488 * MDP(30) + MDP(31) * t653 + MDP(32) * t401) * t494; t457 * t455 * MDP(26) + (-t455 ^ 2 + t457 ^ 2) * MDP(27) + (t699 + t759) * MDP(28) + (-t681 + t758) * MDP(29) + t431 * MDP(30) + (-g(1) * t440 + g(2) * t778 - g(3) * t644 + t401 * t488 - t410 * t457 + t683) * MDP(31) + (g(1) * t441 + g(2) * t777 + g(3) * t471 + t410 * t455 - t653 * t488 - t654) * MDP(32) + (-MDP(28) * t754 - MDP(29) * t457 - MDP(31) * t401 + MDP(32) * t653) * qJD(6);];
tau  = t1;
