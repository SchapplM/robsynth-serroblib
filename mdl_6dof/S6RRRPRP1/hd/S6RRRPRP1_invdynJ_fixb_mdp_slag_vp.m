% Calculate vector of inverse dynamics joint torques for
% S6RRRPRP1
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
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d3,d5,theta4]';
% MDP [28x1]
%   Minimal dynamic parameter vector (fixed base model)
%   see S6RRRPRP1_convert_par2_MPV_fixb.m
% 
% Output:
% tau [6x1]
%   joint torques of inverse dynamics (contains inertial, gravitational Coriolis and centrifugal forces)

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 16:34
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tau = S6RRRPRP1_invdynJ_fixb_mdp_slag_vp(qJ, qJD, qJDD, g, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(3,1),zeros(10,1),zeros(28,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRPRP1_invdynJ_fixb_mdp_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRRPRP1_invdynJ_fixb_mdp_slag_vp: qJD has to be [6x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [6 1]), ...
  'S6RRRPRP1_invdynJ_fixb_mdp_slag_vp: qJDD has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RRRPRP1_invdynJ_fixb_mdp_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RRRPRP1_invdynJ_fixb_mdp_slag_vp: pkin has to be [10x1] (double)');
assert(isreal(MDP) && all(size(MDP) == [28 1]), ...
  'S6RRRPRP1_invdynJ_fixb_mdp_slag_vp: MDP has to be [28x1] (double)'); 

%% Symbolic Calculation
% From invdyn_joint_fixb_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 16:33:16
% EndTime: 2019-03-09 16:33:28
% DurationCPUTime: 7.49s
% Computational Cost: add. (9869->538), mult. (23472->671), div. (0->0), fcn. (17480->14), ass. (0->257)
t620 = cos(qJ(2));
t756 = pkin(7) + pkin(8);
t567 = t756 * t620;
t560 = qJD(1) * t567;
t619 = cos(qJ(3));
t542 = t619 * t560;
t616 = sin(qJ(2));
t566 = t756 * t616;
t558 = qJD(1) * t566;
t615 = sin(qJ(3));
t698 = qJD(1) * t620;
t679 = t619 * t698;
t699 = qJD(1) * t616;
t682 = t615 * t699;
t535 = -t679 + t682;
t739 = qJ(4) * t535;
t493 = t558 * t615 - t542 + t739;
t537 = -t615 * t698 - t619 * t699;
t531 = t537 * qJ(4);
t538 = t615 * t560;
t703 = -t619 * t558 - t538;
t494 = t531 + t703;
t611 = sin(pkin(10));
t612 = cos(pkin(10));
t726 = t611 * t615;
t741 = pkin(2) * qJD(3);
t704 = -t493 * t611 - t494 * t612 + (t612 * t619 - t726) * t741;
t610 = qJ(2) + qJ(3);
t602 = sin(t610);
t603 = cos(t610);
t617 = sin(qJ(1));
t621 = cos(qJ(1));
t660 = g(1) * t621 + g(2) * t617;
t766 = -g(3) * t603 + t602 * t660;
t691 = qJD(1) * qJD(2);
t677 = t620 * t691;
t690 = qJDD(1) * t616;
t764 = t677 + t690;
t605 = t620 * pkin(2);
t742 = pkin(1) + t605;
t740 = qJD(2) * pkin(2);
t545 = -t558 + t740;
t670 = t619 * t545 - t538;
t483 = t531 + t670;
t600 = pkin(10) + t610;
t587 = sin(t600);
t588 = cos(t600);
t618 = cos(qJ(5));
t752 = pkin(5) * t618;
t595 = pkin(4) + t752;
t613 = -qJ(6) - pkin(9);
t647 = -t587 * t613 + t588 * t595;
t702 = -t615 * t566 + t619 * t567;
t614 = sin(qJ(5));
t671 = -t612 * t535 + t537 * t611;
t731 = t671 * t614;
t763 = qJ(6) * t731 + t618 * qJD(6);
t762 = g(1) * t617 - g(2) * t621;
t649 = -t535 * t611 - t612 * t537;
t755 = pkin(3) * t537;
t460 = pkin(4) * t649 - pkin(9) * t671 - t755;
t598 = pkin(2) * t699;
t457 = t460 + t598;
t761 = t614 * t457 - t618 * t704;
t718 = t618 * t621;
t722 = t614 * t617;
t520 = t588 * t722 + t718;
t719 = t617 * t618;
t721 = t614 * t621;
t522 = -t588 * t721 + t719;
t760 = -g(1) * t522 + g(2) * t520;
t759 = qJD(5) - t671;
t688 = qJD(2) + qJD(3);
t511 = qJDD(2) * pkin(2) - t756 * t764;
t678 = t616 * t691;
t689 = qJDD(1) * t620;
t515 = t756 * (-t678 + t689);
t665 = qJD(3) * t545 + t515;
t697 = qJD(3) * t615;
t758 = t615 * t511 - t560 * t697 + t619 * t665;
t491 = t614 * t688 + t618 * t649;
t757 = t491 ^ 2;
t754 = pkin(3) * t602;
t753 = pkin(3) * t612;
t746 = g(3) * t587;
t745 = g(3) * t588;
t743 = g(3) * t614;
t474 = pkin(3) * t688 + t483;
t648 = -t545 * t615 - t542;
t484 = -t648 - t739;
t725 = t612 * t484;
t438 = t611 * t474 + t725;
t436 = pkin(9) * t688 + t438;
t565 = t742 * qJD(1);
t512 = pkin(3) * t535 + qJD(4) - t565;
t449 = -pkin(4) * t671 - pkin(9) * t649 + t512;
t415 = -t436 * t614 + t618 * t449;
t410 = -qJ(6) * t491 + t415;
t403 = pkin(5) * t759 + t410;
t738 = t403 * t618;
t485 = qJD(3) * t679 + t615 * t689 + t619 * t764 - t682 * t688;
t552 = t615 * t620 + t616 * t619;
t510 = t688 * t552;
t655 = t615 * t690 - t619 * t689;
t486 = qJD(1) * t510 + t655;
t650 = t612 * t485 - t611 * t486;
t663 = t618 * t688;
t687 = qJDD(2) + qJDD(3);
t694 = qJD(5) * t614;
t430 = -qJD(5) * t663 - t614 * t687 - t618 * t650 + t649 * t694;
t737 = t430 * t614;
t477 = t611 * t484;
t437 = t612 * t474 - t477;
t435 = -pkin(4) * t688 - t437;
t736 = t435 * t671;
t489 = t614 * t649 - t663;
t735 = t489 * t649;
t734 = t489 * t671;
t733 = t491 * t649;
t732 = t491 * t614;
t551 = t615 * t616 - t619 * t620;
t507 = -t551 * t611 + t552 * t612;
t730 = t507 * t614;
t729 = t507 * t618;
t724 = t612 * t615;
t447 = -t485 * t611 - t612 * t486;
t444 = qJDD(5) - t447;
t723 = t614 * t444;
t720 = t615 * t515;
t604 = t618 * qJ(6);
t439 = t618 * t444;
t669 = -t619 * t566 - t567 * t615;
t499 = -qJ(4) * t552 + t669;
t500 = -qJ(4) * t551 + t702;
t467 = t499 * t611 + t500 * t612;
t464 = t618 * t467;
t596 = pkin(2) * t619 + pkin(3);
t530 = pkin(2) * t724 + t611 * t596;
t525 = pkin(9) + t530;
t717 = -qJ(6) - t525;
t589 = pkin(3) * t611 + pkin(9);
t716 = -qJ(6) - t589;
t715 = -t410 + t403;
t631 = t614 * t650 - t618 * t687;
t431 = qJD(5) * t491 + t631;
t693 = qJD(5) * t618;
t714 = -t614 * t431 - t489 * t693;
t508 = t619 * t511;
t426 = pkin(3) * t687 - t485 * qJ(4) + qJD(3) * t648 + t537 * qJD(4) + t508 - t720;
t432 = -qJ(4) * t486 - qJD(4) * t535 + t758;
t399 = t611 * t426 + t612 * t432;
t713 = t731 * t759 + t439;
t446 = t483 * t612 - t477;
t712 = t618 * t446 + t614 * t460;
t506 = t612 * t551 + t552 * t611;
t662 = pkin(3) * t551 - t742;
t465 = pkin(4) * t506 - pkin(9) * t507 + t662;
t710 = t614 * t465 + t464;
t668 = qJD(5) * t717;
t709 = t614 * t668 - t761 + t763;
t452 = t618 * t457;
t658 = pkin(5) * t649 - t604 * t671;
t708 = t618 * t668 - t452 - t658 + (-qJD(6) - t704) * t614;
t667 = qJD(5) * t716;
t707 = t614 * t667 - t712 + t763;
t459 = t618 * t460;
t706 = t618 * t667 - t459 - t658 + (-qJD(6) + t446) * t614;
t705 = t612 * t493 - t611 * t494 + (t611 * t619 + t724) * t741;
t592 = pkin(3) * t603;
t701 = t592 + t605;
t608 = t616 ^ 2;
t700 = -t620 ^ 2 + t608;
t696 = qJD(3) * t619;
t695 = qJD(5) * t759;
t599 = t616 * t740;
t683 = qJD(2) * t756;
t559 = t616 * t683;
t561 = t620 * t683;
t635 = -t619 * t559 - t615 * t561 - t566 * t696 - t567 * t697;
t453 = -qJ(4) * t510 - qJD(4) * t551 + t635;
t509 = t688 * t551;
t627 = -qJD(3) * t702 + t559 * t615 - t619 * t561;
t454 = qJ(4) * t509 - qJD(4) * t552 + t627;
t420 = t453 * t612 + t454 * t611;
t471 = -t509 * t611 + t612 * t510;
t472 = -t509 * t612 - t510 * t611;
t676 = pkin(3) * t510 + t599;
t424 = pkin(4) * t471 - pkin(9) * t472 + t676;
t685 = t618 * t420 + t614 * t424 + t465 * t693;
t681 = t507 * t693;
t680 = t589 * t695;
t434 = t435 * t693;
t607 = -qJ(4) - t756;
t675 = pkin(5) * t614 - t607;
t398 = t612 * t426 - t611 * t432;
t396 = -pkin(4) * t687 - t398;
t673 = -t396 - t745;
t419 = t453 * t611 - t612 * t454;
t445 = t483 * t611 + t725;
t466 = -t612 * t499 + t500 * t611;
t529 = -pkin(2) * t726 + t596 * t612;
t666 = t759 * t618;
t397 = pkin(9) * t687 + t399;
t664 = -qJD(5) * t449 - t397;
t524 = -pkin(4) - t529;
t661 = (t694 - t731) * pkin(5);
t657 = t592 + t647;
t532 = pkin(2) * t678 - qJDD(1) * t742;
t628 = t486 * pkin(3) + qJDD(4) + t532;
t409 = -t447 * pkin(4) - pkin(9) * t650 + t628;
t408 = t618 * t409;
t656 = -t436 * t693 + t408;
t416 = t436 * t618 + t449 * t614;
t411 = -qJ(6) * t489 + t416;
t654 = -t411 * t614 - t738;
t653 = -t444 * t525 - t736;
t652 = -t444 * t589 - t736;
t651 = t437 * t671 + t438 * t649;
t646 = -t587 * t595 - t588 * t613;
t645 = -qJ(6) * t472 - qJD(6) * t507;
t644 = t396 * t614 + t416 * t649 + t588 * t743 + t434;
t643 = -t415 * t649 + t435 * t694 + (g(1) * t718 + g(2) * t719) * t587;
t642 = t660 * t587;
t641 = -0.2e1 * pkin(1) * t691 - pkin(7) * qJDD(2);
t639 = t472 * t614 + t681;
t638 = t472 * t618 - t507 * t694;
t634 = t618 * t397 + t614 * t409 - t436 * t694 + t449 * t693;
t632 = -t565 * t537 + t508 + t766;
t622 = qJD(2) ^ 2;
t630 = 0.2e1 * qJDD(1) * pkin(1) - pkin(7) * t622 + t762;
t623 = qJD(1) ^ 2;
t629 = pkin(1) * t623 - pkin(7) * qJDD(1) + t660;
t392 = t431 * pkin(5) + qJDD(6) + t396;
t626 = g(3) * t602 - t565 * t535 + t603 * t660 - t758;
t625 = -t537 * t535 * MDP(11) - t759 * t649 * MDP(24) + ((-t430 + t734) * t618 - t759 * t732 + t714) * MDP(21) + (-t694 * t759 + t713 + t735) * MDP(23) + (t666 * t759 + t723 - t733) * MDP(22) + (t491 * t666 - t737) * MDP(20) + (t535 * t688 + t485) * MDP(13) + (-t655 + (-qJD(1) * t552 - t537) * t688) * MDP(14) + (-t535 ^ 2 + t537 ^ 2) * MDP(12) + t687 * MDP(15);
t386 = pkin(5) * t444 + qJ(6) * t430 - qJD(5) * t416 - qJD(6) * t491 - t397 * t614 + t408;
t388 = -qJ(6) * t431 - qJD(6) * t489 + t634;
t624 = qJD(5) * t654 - t386 * t614 + t388 * t618 + t411 * t731 - t588 * t660 + t671 * t738 - t746;
t590 = -pkin(4) - t753;
t562 = -pkin(2) * t616 - t754;
t557 = pkin(1) + t701;
t547 = t589 * t618 + t604;
t546 = t716 * t614;
t544 = t621 * t557;
t523 = t588 * t718 + t722;
t521 = -t588 * t719 + t721;
t514 = t525 * t618 + t604;
t513 = t717 * t614;
t488 = t489 ^ 2;
t463 = t618 * t465;
t427 = t489 * pkin(5) + qJD(6) + t435;
t422 = t618 * t424;
t417 = -qJ(6) * t730 + t710;
t414 = pkin(5) * t506 - t467 * t614 - t507 * t604 + t463;
t390 = -qJ(6) * t681 + (-qJD(5) * t467 + t645) * t614 + t685;
t389 = pkin(5) * t471 - t420 * t614 + t422 + t645 * t618 + (-t464 + (qJ(6) * t507 - t465) * t614) * qJD(5);
t1 = [(-t398 * t507 - t399 * t506 + t419 * t649 + t420 * t671 - t437 * t472 - t438 * t471 + t467 * t447 + t466 * t650 - t660) * MDP(18) + (t444 * t506 + t471 * t759) * MDP(24) + (-(-t467 * t694 + t685) * t759 - t710 * t444 - t634 * t506 - t416 * t471 + t419 * t491 - t466 * t430 + t396 * t729 - g(1) * t520 - g(2) * t522 + t638 * t435) * MDP(26) + (-t431 * t506 - t471 * t489 - t507 * t723 - t639 * t759) * MDP(23) + (-t430 * t506 + t439 * t507 + t471 * t491 + t638 * t759) * MDP(22) + ((-t467 * t693 + t422) * t759 + t463 * t444 + t656 * t506 + t415 * t471 + t419 * t489 + t466 * t431 + t507 * t434 - g(1) * t521 - g(2) * t523 + ((-qJD(5) * t465 - t420) * t759 - t467 * t444 + t664 * t506 + t396 * t507 + t435 * t472) * t614) * MDP(25) + (-t389 * t491 - t390 * t489 + t414 * t430 - t417 * t431 + t762 * t587 + t654 * t472 + (-t386 * t618 - t388 * t614 + (t403 * t614 - t411 * t618) * qJD(5)) * t507) * MDP(27) + t762 * MDP(2) + (-t485 * t551 - t486 * t552 + t509 * t535 + t510 * t537) * MDP(12) + (t485 * t552 + t509 * t537) * MDP(11) + (-t509 * t688 + t552 * t687) * MDP(13) + (-t486 * t742 - t565 * t510 + t532 * t551 + t535 * t599 + t603 * t762 + t627 * t688 + t669 * t687) * MDP(16) + (-t485 * t742 + t565 * t509 + t532 * t552 - t537 * t599 - t602 * t762 - t635 * t688 - t687 * t702) * MDP(17) + qJDD(1) * MDP(1) + (qJDD(2) * t616 + t620 * t622) * MDP(6) + (qJDD(2) * t620 - t616 * t622) * MDP(7) + ((-t489 * t618 - t732) * t472 + (t737 - t431 * t618 + (t489 * t614 - t491 * t618) * qJD(5)) * t507) * MDP(21) + (-t430 * t729 + t491 * t638) * MDP(20) + (t388 * t417 + t411 * t390 + t386 * t414 + t403 * t389 + t392 * (pkin(5) * t730 + t466) + t427 * (pkin(5) * t639 + t419) - g(2) * t544 + (-g(1) * t675 - g(2) * t647) * t621 + (-g(1) * (-t557 - t647) - g(2) * t675) * t617) * MDP(28) + (t616 * t641 + t620 * t630) * MDP(9) + (-t616 * t630 + t620 * t641) * MDP(10) + t660 * MDP(3) + (t399 * t467 + t438 * t420 - t398 * t466 - t437 * t419 + t628 * t662 + t512 * t676 - g(1) * (-t557 * t617 - t607 * t621) - g(2) * (-t607 * t617 + t544)) * MDP(19) + (qJDD(1) * t608 + 0.2e1 * t616 * t677) * MDP(4) + (-t510 * t688 - t551 * t687) * MDP(14) + 0.2e1 * (t616 * t689 - t691 * t700) * MDP(5); (t530 * t447 - t529 * t650 + t649 * t705 + t671 * t704 + t651) * MDP(18) + (t430 * t513 - t431 * t514 - t489 * t709 - t491 * t708 + t624) * MDP(27) + (-g(3) * t620 + t616 * t629) * MDP(9) + (-t524 * t430 + t653 * t618 - t614 * t642 + t705 * t491 + (t525 * t694 + t761) * t759 + t644) * MDP(26) + (g(3) * t616 + t620 * t629) * MDP(10) + MDP(7) * t689 + MDP(6) * t690 + (t388 * t514 + t386 * t513 + t392 * (t524 - t752) - g(3) * (t605 + t657) + (t661 + t705) * t427 + t709 * t411 + t708 * t403 + t660 * (-t562 - t646)) * MDP(28) + (t399 * t530 + t398 * t529 - t512 * (t598 - t755) - g(3) * t701 - t660 * t562 + t704 * t438 - t705 * t437) * MDP(19) + (t524 * t431 + t673 * t618 + t653 * t614 + t705 * t489 + (-t525 * t693 - t614 * t704 - t452) * t759 + t643) * MDP(25) + qJDD(2) * MDP(8) + (t703 * t688 + (t537 * t699 - t615 * t687 - t688 * t696) * pkin(2) + t626) * MDP(17) + (-t560 * t696 + t542 * t688 + (-t558 * t688 - t665) * t615 + (-t535 * t699 + t619 * t687 - t688 * t697) * pkin(2) + t632) * MDP(16) + t625 + (-MDP(4) * t616 * t620 + MDP(5) * t700) * t623; (t430 * t546 - t431 * t547 - t489 * t707 - t491 * t706 + t624) * MDP(27) + (-t445 * t649 - t446 * t671 + (t611 * t447 - t612 * t650) * pkin(3) + t651) * MDP(18) + (t670 * t688 + t626) * MDP(17) + (t590 * t431 - t445 * t489 - t459 * t759 + (t446 * t759 + t652) * t614 + (t673 - t680) * t618 + t643) * MDP(25) + (-qJD(2) * t648 + t632 - t720) * MDP(16) + (t388 * t547 + t386 * t546 + t392 * (-t595 - t753) - g(3) * t657 + (-t445 + t661) * t427 + t707 * t411 + t706 * t403 + t660 * (-t646 + t754)) * MDP(28) + (-t590 * t430 + t712 * t759 - t445 * t491 + t652 * t618 + (-t642 + t680) * t614 + t644) * MDP(26) + t625 + (t437 * t445 - t438 * t446 + (t398 * t612 + t399 * t611 + t512 * t537 + t766) * pkin(3)) * MDP(19); (-t649 ^ 2 - t671 ^ 2) * MDP(18) + (t437 * t649 - t438 * t671 + t628 - t762) * MDP(19) + (t713 - t735) * MDP(25) - MDP(26) * t733 + t714 * MDP(27) + (-t427 * t649 - t762) * MDP(28) + (-MDP(25) * t695 - t444 * MDP(26) + (-t403 * t759 + t388) * MDP(28) + t759 * MDP(27) * t491) * t614 + ((t430 + t734) * MDP(27) + (t411 * t759 + t386) * MDP(28) - t759 ^ 2 * MDP(26)) * t618; t491 * t489 * MDP(20) + (-t488 + t757) * MDP(21) + (t489 * t759 - t430) * MDP(22) + (-t631 + (-qJD(5) + t759) * t491) * MDP(23) + t444 * MDP(24) + (t416 * t759 - t435 * t491 + (t664 + t746) * t614 + t656 + t760) * MDP(25) + (g(1) * t523 - g(2) * t521 + t415 * t759 + t435 * t489 + t618 * t746 - t634) * MDP(26) + (pkin(5) * t430 - t489 * t715) * MDP(27) + (t715 * t411 + (-t427 * t491 + t587 * t743 + t386 + t760) * pkin(5)) * MDP(28); (-t488 - t757) * MDP(27) + (t403 * t491 + t411 * t489 + t392 - t642 + t745) * MDP(28);];
tau  = t1;
