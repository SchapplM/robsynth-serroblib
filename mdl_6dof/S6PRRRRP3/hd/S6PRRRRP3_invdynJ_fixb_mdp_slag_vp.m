% Calculate vector of inverse dynamics joint torques for
% S6PRRRRP3
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
% pkin [11x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d2,d3,d4,d5,theta1]';
% MDP [27x1]
%   Minimal dynamic parameter vector (fixed base model)
%   see S6PRRRRP3_convert_par2_MPV_fixb.m
% 
% Output:
% tau [6x1]
%   joint torques of inverse dynamics (contains inertial, gravitational Coriolis and centrifugal forces)

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 00:12
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tau = S6PRRRRP3_invdynJ_fixb_mdp_slag_vp(qJ, qJD, qJDD, g, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(3,1),zeros(11,1),zeros(27,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRRRRP3_invdynJ_fixb_mdp_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6PRRRRP3_invdynJ_fixb_mdp_slag_vp: qJD has to be [6x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [6 1]), ...
  'S6PRRRRP3_invdynJ_fixb_mdp_slag_vp: qJDD has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6PRRRRP3_invdynJ_fixb_mdp_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6PRRRRP3_invdynJ_fixb_mdp_slag_vp: pkin has to be [11x1] (double)');
assert(isreal(MDP) && all(size(MDP) == [27 1]), ...
  'S6PRRRRP3_invdynJ_fixb_mdp_slag_vp: MDP has to be [27x1] (double)'); 

%% Symbolic Calculation
% From invdyn_joint_fixb_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 00:11:26
% EndTime: 2019-03-09 00:11:38
% DurationCPUTime: 9.90s
% Computational Cost: add. (5717->569), mult. (12991->767), div. (0->0), fcn. (9970->14), ass. (0->242)
t602 = sin(qJ(3));
t606 = cos(qJ(3));
t640 = pkin(3) * t602 - pkin(9) * t606;
t545 = t640 * qJD(3);
t601 = sin(qJ(4));
t603 = sin(qJ(2));
t605 = cos(qJ(4));
t690 = qJD(3) * t602;
t598 = sin(pkin(6));
t698 = qJD(1) * t598;
t607 = cos(qJ(2));
t717 = t606 * t607;
t745 = pkin(8) * t601;
t777 = (-t601 * t717 + t603 * t605) * t698 - t605 * t545 - t690 * t745;
t551 = -pkin(3) * t606 - pkin(9) * t602 - pkin(2);
t684 = qJD(4) * t605;
t776 = -(t601 * t603 + t605 * t717) * t698 + t601 * t545 + t551 * t684;
t718 = t605 * t606;
t579 = pkin(8) * t718;
t634 = pkin(4) * t602 - pkin(10) * t718;
t775 = t634 * qJD(3) + (-t579 + (pkin(10) * t602 - t551) * t601) * qJD(4) - t777;
t686 = qJD(4) * t601;
t689 = qJD(3) * t605;
t688 = qJD(3) * t606;
t665 = t601 * t688;
t766 = t602 * t684 + t665;
t774 = -t766 * pkin(10) + (-t602 * t689 - t606 * t686) * pkin(8) + t776;
t542 = t640 * qJD(2);
t524 = t605 * t542;
t746 = pkin(9) + pkin(10);
t673 = qJD(4) * t746;
t547 = qJD(2) * pkin(8) + t603 * t698;
t599 = cos(pkin(6));
t722 = t599 * t606;
t756 = qJD(1) * t722 - t602 * t547;
t773 = qJD(2) * t634 - t601 * t756 + t605 * t673 + t524;
t693 = qJD(2) * t606;
t667 = t601 * t693;
t705 = t601 * t542 + t605 * t756;
t768 = -pkin(10) * t667 + t601 * t673 + t705;
t733 = sin(pkin(11));
t651 = t733 * t603;
t734 = cos(pkin(11));
t652 = t734 * t607;
t510 = -t599 * t652 + t651;
t650 = t733 * t607;
t653 = t734 * t603;
t512 = t599 * t650 + t653;
t639 = g(1) * t512 + g(2) * t510;
t723 = t598 * t607;
t772 = -g(3) * t723 + t639;
t685 = qJD(4) * t602;
t658 = qJD(2) * t685;
t680 = qJD(2) * qJD(3);
t659 = t606 * t680;
t677 = qJDD(2) * t602;
t767 = -t659 - t677;
t765 = qJD(3) * qJD(4) - t767;
t463 = (qJDD(3) - t658) * t601 + t765 * t605;
t694 = qJD(2) * t602;
t534 = -t601 * t694 + t689;
t691 = qJD(3) * t601;
t535 = t605 * t694 + t691;
t600 = sin(qJ(5));
t604 = cos(qJ(5));
t642 = t601 * t765 + t605 * t658;
t621 = qJDD(3) * t605 - t642;
t682 = qJD(5) * t604;
t683 = qJD(5) * t600;
t410 = -t604 * t463 - t534 * t682 + t535 * t683 - t600 * t621;
t635 = t534 * t600 + t604 * t535;
t411 = qJD(5) * t635 + t463 * t600 - t604 * t621;
t469 = -t604 * t534 + t535 * t600;
t467 = t469 ^ 2;
t589 = t606 * qJDD(2);
t754 = -t602 * t680 + t589;
t531 = qJDD(4) - t754;
t526 = qJDD(5) + t531;
t577 = -qJD(4) + t693;
t566 = -qJD(5) + t577;
t747 = t635 ^ 2;
t771 = t526 * MDP(23) + (-t566 * t635 - t411) * MDP(22) + t469 * MDP(19) * t635 + (-t469 * t566 - t410) * MDP(21) + (-t467 + t747) * MDP(20);
t769 = qJ(6) * t469;
t721 = t600 * t601;
t536 = -t604 * t605 + t721;
t676 = qJD(4) + qJD(5);
t707 = -t536 * t693 - t604 * t684 - t605 * t682 + t676 * t721;
t537 = t600 * t605 + t601 * t604;
t480 = t676 * t537;
t706 = -t537 * t693 + t480;
t511 = t599 * t653 + t650;
t513 = -t599 * t651 + t652;
t638 = g(1) * t513 + g(2) * t511;
t724 = t598 * t603;
t619 = -g(3) * t724 - t638;
t492 = -qJD(3) * pkin(3) - t756;
t457 = -pkin(4) * t534 + t492;
t655 = t598 * t734;
t482 = t511 * t606 - t602 * t655;
t654 = t598 * t733;
t484 = t513 * t606 + t602 * t654;
t518 = t599 * t602 + t606 * t724;
t597 = qJ(4) + qJ(5);
t590 = sin(t597);
t591 = cos(t597);
t681 = qJD(1) * qJD(2);
t503 = qJDD(2) * pkin(8) + (qJDD(1) * t603 + t607 * t681) * t598;
t678 = qJDD(1) * t599;
t657 = t602 * t678;
t437 = qJDD(3) * pkin(9) + qJD(3) * t756 + t503 * t606 + t657;
t697 = qJD(1) * t602;
t571 = t599 * t697;
t499 = t606 * t547 + t571;
t493 = qJD(3) * pkin(9) + t499;
t696 = qJD(1) * t607;
t672 = t598 * t696;
t501 = qJD(2) * t551 - t672;
t449 = t493 * t605 + t501 * t601;
t661 = t603 * t681;
t637 = -qJDD(1) * t723 + t598 * t661;
t456 = qJD(2) * t545 + qJDD(2) * t551 + t637;
t455 = t605 * t456;
t615 = -qJD(4) * t449 - t601 * t437 + t455;
t395 = pkin(4) * t531 - pkin(10) * t463 + t615;
t674 = t605 * t437 + t601 * t456 + t501 * t684;
t627 = -t493 * t686 + t674;
t400 = pkin(10) * t621 + t627;
t448 = -t493 * t601 + t605 * t501;
t429 = -pkin(10) * t535 + t448;
t419 = -pkin(4) * t577 + t429;
t430 = pkin(10) * t534 + t449;
t643 = -t600 * t395 - t604 * t400 - t419 * t682 + t430 * t683;
t764 = t457 * t469 - g(1) * (-t484 * t591 - t512 * t590) - g(2) * (-t482 * t591 - t510 * t590) - g(3) * (-t518 * t591 + t590 * t723) + t643;
t761 = qJ(6) * t635;
t420 = pkin(5) * t469 + qJD(6) + t457;
t760 = t420 * t635;
t759 = t773 * t604;
t758 = t774 * t600 - t604 * t775;
t533 = t605 * t551;
t719 = t602 * t605;
t476 = -pkin(10) * t719 + t533 + (-pkin(4) - t745) * t606;
t701 = t601 * t551 + t579;
t720 = t601 * t602;
t489 = -pkin(10) * t720 + t701;
t757 = t476 * t682 - t489 * t683 + t600 * t775 + t774 * t604;
t641 = -t499 + (-t667 + t686) * pkin(4);
t708 = t600 * t476 + t604 * t489;
t559 = t746 * t601;
t560 = t746 * t605;
t702 = -t600 * t559 + t604 * t560;
t755 = -t559 * t682 - t560 * t683 - t600 * t773 - t604 * t768;
t664 = t605 * t688;
t753 = -t601 * t685 + t664;
t752 = -g(3) * (-t518 * t590 - t591 * t723) - g(2) * (-t482 * t590 + t510 * t591) - g(1) * (-t484 * t590 + t512 * t591);
t426 = t604 * t430;
t407 = t419 * t600 + t426;
t649 = t604 * t395 - t600 * t400;
t613 = -qJD(5) * t407 + t649;
t751 = -t457 * t635 + t613 + t752;
t608 = qJD(3) ^ 2;
t750 = 0.2e1 * qJDD(2) * pkin(2) - pkin(8) * t608 + t598 * (-g(3) * t607 + t661) - t637 + t639;
t749 = qJD(4) * (pkin(8) * t577 + t493) + t772;
t593 = t605 * pkin(4);
t737 = pkin(3) + t593;
t735 = qJD(2) * pkin(2);
t731 = qJDD(3) * pkin(3);
t730 = t463 * t601;
t729 = t534 * t577;
t728 = t535 * t577;
t727 = t535 * t605;
t726 = t590 * t606;
t725 = t591 * t606;
t424 = t600 * t430;
t716 = qJDD(1) - g(3);
t439 = t480 * t602 + t600 * t665 - t604 * t664;
t507 = t536 * t602;
t715 = pkin(5) * t690 + qJ(6) * t439 - qJD(5) * t708 + qJD(6) * t507 - t758;
t440 = -t683 * t720 + (t676 * t719 + t665) * t604 + t753 * t600;
t506 = t537 * t602;
t714 = -qJ(6) * t440 - qJD(6) * t506 + t757;
t406 = t604 * t419 - t424;
t402 = t406 - t761;
t401 = -pkin(5) * t566 + t402;
t713 = -t402 + t401;
t712 = -qJ(6) * t706 - qJD(6) * t536 + t755;
t711 = -pkin(5) * t694 + qJ(6) * t707 - t702 * qJD(5) - qJD(6) * t537 + t600 * t768 - t759;
t710 = t604 * t429 - t424;
t550 = pkin(5) * t591 + t593;
t546 = pkin(4) * t720 + t602 * pkin(8);
t595 = t602 ^ 2;
t700 = -t606 ^ 2 + t595;
t695 = qJD(2) * t598;
t692 = qJD(3) * t534;
t687 = qJD(4) * t577;
t500 = pkin(4) * t766 + pkin(8) * t688;
t670 = t602 * t696;
t669 = t603 * t695;
t668 = t607 * t695;
t666 = t577 * t689;
t648 = -t429 * t600 - t426;
t645 = t604 * t476 - t489 * t600;
t644 = -t604 * t559 - t560 * t600;
t487 = -t518 * t601 - t605 * t723;
t631 = -t518 * t605 + t601 * t723;
t433 = t487 * t604 + t600 * t631;
t434 = t487 * t600 - t604 * t631;
t517 = t602 * t724 - t722;
t629 = t531 * t601 - t577 * t684;
t628 = t531 * t605 + t577 * t686;
t481 = t511 * t602 + t606 * t655;
t483 = t513 * t602 - t606 * t654;
t624 = g(1) * t483 + g(2) * t481 + g(3) * t517;
t623 = g(1) * t484 + g(2) * t482 + g(3) * t518;
t622 = qJD(3) * t571 + t602 * t503 + t547 * t688 - t606 * t678;
t618 = -pkin(9) * t531 - t492 * t577;
t438 = t622 - t731;
t612 = -pkin(9) * t687 + t438 - t624;
t548 = -t672 - t735;
t611 = -pkin(8) * qJDD(3) + (t548 + t672 - t735) * qJD(3);
t414 = -pkin(4) * t621 + t438;
t398 = t411 * pkin(5) + qJDD(6) + t414;
t609 = qJD(2) ^ 2;
t594 = -qJ(6) - t746;
t585 = pkin(4) * t604 + pkin(5);
t549 = pkin(4) * t601 + pkin(5) * t590;
t541 = pkin(3) + t550;
t486 = qJD(3) * t518 + t602 * t668;
t485 = -qJD(3) * t517 + t606 * t668;
t459 = -qJ(6) * t536 + t702;
t458 = -qJ(6) * t537 + t644;
t423 = qJD(4) * t487 + t485 * t605 + t601 * t669;
t422 = qJD(4) * t631 - t485 * t601 + t605 * t669;
t416 = -qJ(6) * t506 + t708;
t415 = -pkin(5) * t606 + qJ(6) * t507 + t645;
t405 = t710 - t761;
t404 = t648 + t769;
t403 = t407 - t769;
t397 = -qJD(5) * t434 + t422 * t604 - t423 * t600;
t396 = qJD(5) * t433 + t422 * t600 + t423 * t604;
t390 = -qJ(6) * t411 - qJD(6) * t469 - t643;
t389 = pkin(5) * t526 + qJ(6) * t410 - qJD(6) * t635 + t613;
t1 = [t716 * MDP(1) + (-qJD(3) * t486 - qJDD(3) * t517) * MDP(10) + (-qJD(3) * t485 - qJDD(3) * t518) * MDP(11) + (-t422 * t577 - t486 * t534 + t487 * t531 - t517 * t621) * MDP(17) + (t423 * t577 + t463 * t517 + t486 * t535 + t531 * t631) * MDP(18) + (-t397 * t566 + t411 * t517 + t433 * t526 + t469 * t486) * MDP(24) + (t396 * t566 - t410 * t517 - t434 * t526 + t486 * t635) * MDP(25) + (-t396 * t469 - t397 * t635 + t410 * t433 - t411 * t434) * MDP(26) + (t389 * t433 + t390 * t434 + t396 * t403 + t397 * t401 + t398 * t517 + t420 * t486 - g(3)) * MDP(27) + ((-qJDD(2) * MDP(4) + (-MDP(10) * t606 + MDP(11) * t602 - MDP(3)) * t609) * t603 + (t754 * MDP(10) + MDP(11) * t767 + qJDD(2) * MDP(3) - t609 * MDP(4)) * t607) * t598; (t611 * t602 + t606 * t750) * MDP(10) + (-t602 * t750 + t611 * t606) * MDP(11) + (t390 * t416 + t389 * t415 + t398 * (pkin(5) * t506 + t546) + t714 * t403 + t715 * t401 + t619 * (pkin(8) + t549) + (pkin(5) * t440 - t697 * t723 + t500) * t420 + t772 * (t541 * t606 - t594 * t602 + pkin(2))) * MDP(27) + ((t534 * t605 - t535 * t601) * t688 + (t605 * t621 - t730 + (-t601 * t534 - t727) * qJD(4)) * t602) * MDP(13) + (t716 * t723 + t639) * MDP(3) + (-t716 * t724 + t638) * MDP(4) + 0.2e1 * (t589 * t602 - t680 * t700) * MDP(6) + (-t531 * t606 - t577 * t690) * MDP(16) + (-t526 * t606 - t566 * t690) * MDP(23) + (t411 * t606 + t440 * t566 - t469 * t690 - t506 * t526) * MDP(22) + ((t577 * t691 - t621) * t606 + (-t629 + t692) * t602) * MDP(15) + ((-t463 - t666) * t606 + (qJD(3) * t535 + t628) * t602) * MDP(14) + (t463 * t719 + t535 * t753) * MDP(12) + (qJDD(2) * t595 + 0.2e1 * t602 * t659) * MDP(5) + (t389 * t507 - t390 * t506 + t401 * t439 - t403 * t440 + t410 * t415 - t411 * t416 - t469 * t714 + t602 * t772 - t635 * t715) * MDP(26) + (t410 * t606 + t439 * t566 - t507 * t526 + t635 * t690) * MDP(21) + (t410 * t506 + t411 * t507 + t439 * t469 - t440 * t635) * MDP(20) + (t410 * t507 - t439 * t635) * MDP(19) + (-t708 * t526 - t643 * t606 - t407 * t690 + t500 * t635 - t546 * t410 - t414 * t507 - t457 * t439 - g(1) * (t512 * t726 + t513 * t591) - g(2) * (t510 * t726 + t511 * t591) + t757 * t566 + (-t635 * t670 - g(3) * (-t590 * t717 + t591 * t603)) * t598) * MDP(25) + (qJDD(3) * t602 + t606 * t608) * MDP(7) + (qJDD(3) * t606 - t602 * t608) * MDP(8) + (t645 * t526 - t649 * t606 + t406 * t690 + t500 * t469 + t546 * t411 + t414 * t506 + t457 * t440 - g(1) * (-t512 * t725 + t513 * t590) - g(2) * (-t510 * t725 + t511 * t590) + t758 * t566 + (t407 * t606 + t566 * t708) * qJD(5) + (-t469 * t670 - g(3) * (t590 * t603 + t591 * t717)) * t598) * MDP(24) + qJDD(2) * MDP(2) + (-t701 * t531 + t776 * t577 + t619 * t605 + ((pkin(8) * t535 + t492 * t605) * qJD(3) - t749 * t601 + t674) * t606 + (-t535 * t672 - t492 * t686 - t449 * qJD(3) + t438 * t605 + (t463 - t666) * pkin(8)) * t602) * MDP(18) + (t533 * t531 + t777 * t577 + (t551 * t687 + t619) * t601 + (-pkin(8) * t692 - t455 + (-pkin(8) * t531 + qJD(3) * t492 + qJD(4) * t501 + t437) * t601 + t749 * t605) * t606 + (-pkin(8) * t621 + t448 * qJD(3) + t438 * t601 + t492 * t684 + t534 * t672) * t602) * MDP(17); MDP(7) * t677 + MDP(8) * t589 + qJDD(3) * MDP(9) + (qJD(3) * t499 - t622 + t624) * MDP(10) + (-t657 + (-qJD(2) * t548 - t503) * t606 + t623) * MDP(11) + (-t577 * t727 + t730) * MDP(12) + ((t463 - t729) * t605 + (t621 + t728) * t601) * MDP(13) + ((-t535 * t602 + t577 * t718) * qJD(2) + t629) * MDP(14) + ((-t577 * t601 * t606 - t534 * t602) * qJD(2) + t628) * MDP(15) + (-pkin(3) * t642 + t524 * t577 + t499 * t534 + (-t577 * t756 + t618) * t601 + (-t612 + t731) * t605) * MDP(17) + (-pkin(3) * t463 - t499 * t535 - t577 * t705 + t601 * t612 + t605 * t618) * MDP(18) + (-t410 * t537 - t635 * t707) * MDP(19) + (t410 * t536 - t411 * t537 + t469 * t707 - t635 * t706) * MDP(20) + (t526 * t537 + t566 * t707) * MDP(21) + (-t526 * t536 + t566 * t706) * MDP(22) + (t644 * t526 - t737 * t411 + t414 * t536 + (t560 * t682 + (-qJD(5) * t559 - t768) * t600 + t759) * t566 + t641 * t469 + t706 * t457 + t624 * t591) * MDP(24) + (t410 * t737 + t414 * t537 - t707 * t457 - t702 * t526 + t566 * t755 - t624 * t590 + t641 * t635) * MDP(25) + (-t389 * t537 - t390 * t536 + t401 * t707 - t403 * t706 + t410 * t458 - t411 * t459 - t469 * t712 - t635 * t711 - t623) * MDP(26) + (t390 * t459 + t389 * t458 + t398 * (pkin(5) * t536 - t737) - g(1) * (-t483 * t541 - t484 * t594) - g(2) * (-t481 * t541 - t482 * t594) - g(3) * (-t517 * t541 - t518 * t594) + (pkin(5) * t706 + t641) * t420 + t712 * t403 + t711 * t401) * MDP(27) + (-MDP(10) * t548 + t577 * MDP(16) - t448 * MDP(17) + MDP(18) * t449 - MDP(21) * t635 + t469 * MDP(22) + t566 * MDP(23) - t406 * MDP(24) + t407 * MDP(25)) * t694 + (-MDP(5) * t602 * t606 + MDP(6) * t700) * t609; -t535 * t534 * MDP(12) + (-t534 ^ 2 + t535 ^ 2) * MDP(13) + (t463 + t729) * MDP(14) + (t621 - t728) * MDP(15) + t531 * MDP(16) + (-t449 * t577 - t492 * t535 - g(1) * (-t484 * t601 + t512 * t605) - g(2) * (-t482 * t601 + t510 * t605) - g(3) * t487 + t615) * MDP(17) + (-t448 * t577 - t492 * t534 - g(1) * (-t484 * t605 - t512 * t601) - g(2) * (-t482 * t605 - t510 * t601) - g(3) * t631 - t627) * MDP(18) + (t648 * t566 + (-t469 * t535 + t526 * t604 + t566 * t683) * pkin(4) + t751) * MDP(24) + (-t710 * t566 + (-t526 * t600 - t535 * t635 + t566 * t682) * pkin(4) + t764) * MDP(25) + (-t401 * t469 + t403 * t635 + t404 * t635 + t405 * t469 + t410 * t585 + (-t411 * t600 + (-t469 * t604 + t600 * t635) * qJD(5)) * pkin(4)) * MDP(26) + (t389 * t585 - t403 * t405 - t401 * t404 - pkin(5) * t760 - g(1) * (-t484 * t549 + t512 * t550) - g(2) * (-t482 * t549 + t510 * t550) - g(3) * (-t518 * t549 - t550 * t723) + (t390 * t600 - t420 * t535 + (-t401 * t600 + t403 * t604) * qJD(5)) * pkin(4)) * MDP(27) + t771; (-t407 * t566 + t751) * MDP(24) + (-t406 * t566 + t764) * MDP(25) + (pkin(5) * t410 - t469 * t713) * MDP(26) + (t713 * t403 + (t389 + t752 - t760) * pkin(5)) * MDP(27) + t771; (-t467 - t747) * MDP(26) + (t401 * t635 + t403 * t469 + t398 - t624) * MDP(27);];
tau  = t1;
