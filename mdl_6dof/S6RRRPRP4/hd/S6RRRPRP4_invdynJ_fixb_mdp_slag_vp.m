% Calculate vector of inverse dynamics joint torques for
% S6RRRPRP4
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
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d3,d5]';
% MDP [32x1]
%   Minimal dynamic parameter vector (fixed base model)
%   see S6RRRPRP4_convert_par2_MPV_fixb.m
% 
% Output:
% tau [6x1]
%   joint torques of inverse dynamics (contains inertial, gravitational Coriolis and centrifugal forces)

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 16:46
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tau = S6RRRPRP4_invdynJ_fixb_mdp_slag_vp(qJ, qJD, qJDD, g, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(3,1),zeros(9,1),zeros(32,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRPRP4_invdynJ_fixb_mdp_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRRPRP4_invdynJ_fixb_mdp_slag_vp: qJD has to be [6x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [6 1]), ...
  'S6RRRPRP4_invdynJ_fixb_mdp_slag_vp: qJDD has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RRRPRP4_invdynJ_fixb_mdp_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S6RRRPRP4_invdynJ_fixb_mdp_slag_vp: pkin has to be [9x1] (double)');
assert(isreal(MDP) && all(size(MDP) == [32 1]), ...
  'S6RRRPRP4_invdynJ_fixb_mdp_slag_vp: MDP has to be [32x1] (double)'); 

%% Symbolic Calculation
% From invdyn_joint_fixb_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 16:46:14
% EndTime: 2019-03-09 16:46:30
% DurationCPUTime: 12.46s
% Computational Cost: add. (7919->650), mult. (17058->756), div. (0->0), fcn. (11921->10), ass. (0->272)
t598 = cos(qJ(2));
t752 = cos(qJ(3));
t677 = t752 * t598;
t657 = qJD(1) * t677;
t594 = sin(qJ(3));
t595 = sin(qJ(2));
t701 = qJD(1) * t595;
t676 = t594 * t701;
t516 = -t657 + t676;
t587 = qJDD(2) + qJDD(3);
t588 = qJD(2) + qJD(3);
t593 = sin(qJ(5));
t597 = cos(qJ(5));
t528 = t594 * t598 + t595 * t752;
t487 = t588 * t528;
t670 = qJDD(1) * t752;
t560 = t598 * t670;
t689 = qJDD(1) * t595;
t647 = t594 * t689 - t560;
t604 = qJD(1) * t487 + t647;
t696 = qJD(5) * t597;
t697 = qJD(5) * t593;
t434 = -t516 * t696 - t597 * t587 + t588 * t697 - t593 * t604;
t492 = -t597 * t516 + t588 * t593;
t518 = t528 * qJD(1);
t766 = qJD(5) + t518;
t765 = t492 * t766;
t769 = t434 - t765;
t770 = t769 * MDP(30);
t667 = t766 ^ 2;
t494 = t516 * t593 + t588 * t597;
t665 = t766 * t494;
t600 = -pkin(8) - pkin(7);
t541 = t600 * t598;
t531 = qJD(1) * t541;
t520 = t594 * t531;
t540 = t600 * t595;
t529 = qJD(1) * t540;
t485 = t529 * t752 + t520;
t673 = qJD(3) * t752;
t768 = -pkin(2) * t673 + t485;
t767 = t518 * t597 + t696;
t764 = t518 * t588;
t717 = t594 * t595;
t527 = -t677 + t717;
t749 = pkin(2) * t598;
t577 = pkin(1) + t749;
t640 = -qJ(4) * t528 - t577;
t753 = pkin(3) + pkin(9);
t467 = t527 * t753 + t640;
t496 = -t540 * t752 - t594 * t541;
t475 = t528 * pkin(4) + t496;
t710 = t597 * t467 + t593 * t475;
t523 = t752 * t531;
t484 = t594 * t529 - t523;
t700 = qJD(3) * t594;
t685 = pkin(2) * t700;
t763 = t484 - t685;
t648 = -pkin(5) * t597 - qJ(6) * t593;
t639 = -pkin(4) + t648;
t762 = pkin(5) * t696 + qJ(6) * t697 - t597 * qJD(6) - t639 * t518 + qJD(4);
t705 = -qJD(4) + t768;
t592 = qJ(2) + qJ(3);
t584 = sin(t592);
t585 = cos(t592);
t703 = t585 * pkin(3) + t584 * qJ(4);
t578 = t587 * qJ(4);
t761 = -t588 * qJD(4) - t578;
t605 = -t764 - t647;
t699 = qJD(5) * t494;
t435 = t587 * t593 - t597 * t604 + t699;
t586 = t593 * pkin(5);
t738 = qJ(6) * t597;
t760 = -t738 + t586;
t596 = sin(qJ(1));
t599 = cos(qJ(1));
t651 = g(1) * t599 + g(2) * t596;
t740 = qJD(2) * pkin(2);
t525 = t529 + t740;
t480 = -t752 * t525 - t520;
t693 = qJD(4) + t480;
t759 = t577 * qJDD(1);
t758 = -MDP(27) - MDP(29);
t687 = MDP(28) - MDP(31);
t678 = qJD(2) * t600;
t530 = t595 * t678;
t532 = t598 * t678;
t449 = -t752 * t530 - t594 * t532 - t540 * t673 - t541 * t700;
t497 = t594 * t540 - t541 * t752;
t650 = g(1) * t596 - g(2) * t599;
t757 = -t449 * t588 + t497 * t587 + t584 * t650;
t539 = t577 * qJD(1);
t621 = -qJ(4) * t518 - t539;
t448 = t516 * t753 + t621;
t716 = t596 * t597;
t718 = t593 * t599;
t508 = t584 * t718 + t716;
t715 = t597 * t599;
t719 = t593 * t596;
t510 = -t584 * t719 + t715;
t572 = g(3) * t585;
t649 = t588 * t717;
t688 = qJDD(1) * t598;
t659 = t588 * t657 + t594 * t688 + t595 * t670;
t465 = qJD(1) * t649 - t659;
t690 = qJD(1) * qJD(2);
t672 = t595 * t690;
t565 = pkin(2) * t672;
t632 = t465 * qJ(4) - t518 * qJD(4) + t565;
t737 = qJDD(1) * pkin(1);
t399 = -pkin(2) * t688 - pkin(3) * t605 + pkin(9) * t604 + t632 - t737;
t671 = t598 * t690;
t490 = qJDD(2) * pkin(2) - t600 * (-t671 - t689);
t491 = t600 * (-t672 + t688);
t661 = -t752 * t490 - t594 * t491 + t525 * t700 - t531 * t673;
t641 = qJDD(4) + t661;
t405 = -pkin(4) * t465 - t587 * t753 + t641;
t741 = t518 * pkin(4);
t694 = t741 + t693;
t445 = -t588 * t753 + t694;
t682 = t597 * t399 + t593 * t405 + t445 * t696;
t756 = -g(1) * t508 + g(2) * t510 + (-qJD(5) * t448 + t572) * t593 + t682;
t755 = t494 ^ 2;
t754 = t518 ^ 2;
t751 = pkin(2) * t594;
t750 = pkin(2) * t595;
t748 = pkin(3) * t587;
t459 = -qJDD(5) + t465;
t747 = pkin(5) * t459;
t746 = pkin(5) * t516;
t571 = g(3) * t584;
t743 = g(3) * t598;
t742 = t516 * pkin(4);
t569 = t585 * pkin(9);
t739 = qJ(6) * t459;
t662 = t593 * t399 - t597 * t405 + t445 * t697 + t448 * t696;
t390 = qJDD(6) + t662 + t747;
t389 = t390 * t597;
t660 = -t594 * t490 + t752 * t491 - t525 * t673 - t531 * t700;
t426 = t660 + t761;
t406 = -pkin(4) * t604 - t426;
t392 = t435 * pkin(5) + t434 * qJ(6) - t494 * qJD(6) + t406;
t736 = t392 * t597;
t416 = t445 * t593 + t448 * t597;
t410 = qJ(6) * t766 + t416;
t735 = t410 * t518;
t734 = t416 * t766;
t733 = t434 * t597;
t481 = t594 * t525 - t523;
t732 = t481 * t588;
t731 = t492 * t494;
t730 = t516 * t518;
t728 = t527 * t593;
t576 = -pkin(2) * t752 - pkin(3);
t566 = -pkin(9) + t576;
t727 = t566 * t459;
t726 = t584 * t596;
t725 = t584 * t597;
t724 = t584 * t599;
t723 = t585 * t596;
t722 = t585 * t597;
t721 = t585 * t599;
t720 = t593 * t459;
t454 = t597 * t459;
t714 = t599 * t600;
t713 = t753 * t459;
t478 = pkin(3) * t518 + qJ(4) * t516;
t474 = pkin(2) * t701 + t478;
t511 = t518 * pkin(9);
t451 = t474 + t511;
t468 = t484 - t742;
t712 = t597 * t451 + t593 * t468;
t460 = t478 + t511;
t464 = t481 - t742;
t711 = t597 * t460 + t593 * t464;
t709 = t762 - t768;
t708 = t480 + t762;
t542 = qJ(4) * t723;
t686 = t585 * t586;
t707 = t596 * t686 + t542;
t544 = qJ(4) * t721;
t706 = t599 * t686 + t544;
t704 = t741 - t705;
t590 = t595 ^ 2;
t702 = -t598 ^ 2 + t590;
t698 = qJD(5) * t566;
t695 = qJD(5) * t753;
t415 = t445 * t597 - t448 * t593;
t691 = qJD(6) - t415;
t581 = t595 * t740;
t684 = t585 * t716;
t683 = t585 * t715;
t681 = -g(1) * t683 - g(2) * t684 - g(3) * t725;
t680 = pkin(3) * t721 + qJ(4) * t724 + t599 * t577;
t679 = -g(1) * t724 - g(2) * t726 + t572;
t675 = t566 * t696;
t674 = t597 * t695;
t627 = -t448 * t697 + t682;
t388 = qJD(6) * t766 + t627 - t739;
t409 = -pkin(5) * t766 + t691;
t669 = -t409 * t518 - t388;
t583 = t588 * qJ(4);
t453 = t464 + t583;
t666 = t766 * t453;
t663 = t389 - t679;
t658 = t584 * t586 + t569 + t703;
t656 = -g(1) * t723 + g(2) * t721;
t654 = -pkin(3) * t584 - t750;
t507 = -t584 * t715 + t719;
t509 = t584 * t716 + t718;
t653 = g(1) * t509 + g(2) * t507;
t652 = -g(1) * t510 - g(2) * t508;
t537 = qJ(4) + t760;
t645 = t409 * t597 - t410 * t593;
t644 = t409 * t593 + t410 * t597;
t473 = -pkin(3) * t588 + t693;
t477 = -t583 - t481;
t642 = t473 * t516 - t477 * t518;
t425 = pkin(5) * t492 - qJ(6) * t494 + t453;
t638 = t392 * t593 - t409 * t516 + t425 * t767;
t637 = t406 * t593 + t415 * t516 + t453 * t767;
t636 = -t577 - t703;
t635 = -0.2e1 * pkin(1) * t690 - pkin(7) * qJDD(2);
t634 = t487 * t593 + t527 * t696;
t633 = -t487 * t597 + t527 * t697;
t486 = -qJD(2) * t677 - t598 * t673 + t649;
t631 = qJ(4) * t486 - qJD(4) * t528 + t581;
t630 = t406 * t597 - t416 * t516 + t681;
t424 = t487 * t753 + t631;
t450 = qJD(3) * t497 + t594 * t530 - t532 * t752;
t433 = -t486 * pkin(4) + t450;
t626 = t597 * t424 + t593 * t433 - t467 * t697 + t475 * t696;
t625 = t661 + t679;
t624 = g(1) * t721 + g(2) * t723 + t571 + t660;
t620 = -t585 * t651 - t571;
t602 = qJD(2) ^ 2;
t618 = -pkin(7) * t602 + t650 + 0.2e1 * t737;
t603 = qJD(1) ^ 2;
t617 = pkin(1) * t603 - pkin(7) * qJDD(1) + t651;
t616 = t450 * t588 + t496 * t587 + t656;
t615 = g(1) * t507 - g(2) * t509 + g(3) * t722 - t662;
t614 = -t539 * t516 + t624;
t613 = t539 * t518 - t625;
t612 = t410 * t516 - t681 - t736 + (t518 * t593 + t697) * t425;
t471 = pkin(3) * t516 + t621;
t611 = t471 * t518 + qJDD(4) + t625;
t610 = -t471 * t516 - t624 - t761;
t609 = qJD(5) * t644 + t388 * t593 - t389;
t608 = t425 * t494 + qJDD(6) - t615;
t607 = t766 * t516 * MDP(26) + MDP(11) * t730 + ((-t435 - t665) * t597 + (t434 + t765) * t593) * MDP(23) + (-t593 * t665 - t733) * MDP(22) + (t494 * t516 - t593 * t667 - t454) * MDP(24) + (-t492 * t516 - t597 * t667 + t720) * MDP(25) + (t659 + (t516 - t676) * t588) * MDP(13) + (t605 + t764) * MDP(14) + (-t516 ^ 2 + t754) * MDP(12) + t587 * MDP(15);
t568 = qJ(4) + t751;
t524 = t537 + t751;
t513 = t565 - t759;
t506 = t516 * qJ(6);
t479 = pkin(3) * t527 + t640;
t476 = -t527 * pkin(4) + t497;
t452 = pkin(5) * t494 + qJ(6) * t492;
t439 = t527 * t639 + t497;
t437 = pkin(3) * t487 + t631;
t432 = -pkin(4) * t487 - t449;
t429 = t641 - t748;
t428 = -pkin(5) * t528 + t467 * t593 - t475 * t597;
t427 = qJ(6) * t528 + t710;
t420 = t460 * t593 - t464 * t597 + t746;
t419 = -t506 + t711;
t418 = t451 * t593 - t468 * t597 + t746;
t417 = -t506 + t712;
t411 = -t759 + t632 + (qJDD(1) * t717 - t560 + t764) * pkin(3);
t396 = (qJD(5) * t760 - qJD(6) * t593) * t527 + t639 * t487 - t449;
t394 = pkin(5) * t486 + qJD(5) * t710 + t424 * t593 - t433 * t597;
t393 = -qJ(6) * t486 + qJD(6) * t528 + t626;
t1 = [(t426 * t527 + t429 * t528 + t449 * t516 + t450 * t518 - t496 * t465 - t473 * t486 + t477 * t487 - t497 * t604 - t651) * MDP(18) + (-t539 * t487 + t513 * t527 + t516 * t581 - t577 * t604 - t616) * MDP(16) + (t465 * t527 + t486 * t516 - t518 * t487 + t528 * t605) * MDP(12) + (-t459 * t528 - t486 * t766) * MDP(26) + (-t390 * t528 - t394 * t766 + t396 * t492 + t409 * t486 + t425 * t633 + t428 * t459 + t435 * t439 - t527 * t736 + t652) * MDP(29) + (t388 * t528 - t392 * t728 + t393 * t766 - t396 * t494 - t410 * t486 - t425 * t634 - t427 * t459 + t434 * t439 - t653) * MDP(31) + (t406 * t728 + t416 * t486 + t432 * t494 - t476 * t434 + t453 * t634 + t459 * t710 - t528 * t627 - t626 * t766 + t653) * MDP(28) + (-t434 * t528 - t486 * t494 - t527 * t720 + t634 * t766) * MDP(24) + (-t435 * t528 - t454 * t527 + t486 * t492 - t633 * t766) * MDP(25) + (-t662 * t528 - t415 * t486 + t432 * t492 + t476 * t435 + ((-qJD(5) * t475 - t424) * t766 + t467 * t459 + t453 * qJD(5) * t527) * t593 + ((-qJD(5) * t467 + t433) * t766 - t475 * t459 - t406 * t527 - t453 * t487) * t597 + t652) * MDP(27) + (t465 * t577 + t486 * t539 + t513 * t528 + t518 * t581 - t757) * MDP(17) + (-t411 * t528 - t437 * t518 + t465 * t479 + t471 * t486 + t757) * MDP(20) + (t595 * t635 + t598 * t618) * MDP(9) + (-t595 * t618 + t598 * t635) * MDP(10) + qJDD(1) * MDP(1) + (qJDD(2) * t595 + t598 * t602) * MDP(6) + (qJDD(2) * t598 - t595 * t602) * MDP(7) + (-t486 * t588 + t528 * t587) * MDP(13) + (-t487 * t588 - t527 * t587) * MDP(14) + (-t465 * t528 - t486 * t518) * MDP(11) + (t388 * t427 + t410 * t393 + t392 * t439 + t425 * t396 + t390 * t428 + t409 * t394 - g(1) * (pkin(4) * t599 + pkin(5) * t510 + qJ(6) * t509 - t714) - g(2) * (pkin(5) * t508 + pkin(9) * t721 + qJ(6) * t507 + t680) + (-g(1) * (t636 - t569) - g(2) * (pkin(4) - t600)) * t596) * MDP(32) + ((-t492 * t593 + t494 * t597) * t487 + (-t733 - t435 * t593 + (-t492 * t597 - t494 * t593) * qJD(5)) * t527) * MDP(23) + (-t434 * t728 + t494 * t634) * MDP(22) + (t411 * t479 + t471 * t437 - t426 * t497 + t477 * t449 + t429 * t496 + t473 * t450 + g(1) * t714 - g(2) * t680 + (-g(1) * t636 + g(2) * t600) * t596) * MDP(21) + 0.2e1 * (t595 * t688 - t690 * t702) * MDP(5) + t650 * MDP(2) + t651 * MDP(3) + (-t393 * t492 + t394 * t494 - t427 * t435 - t428 * t434 + t644 * t487 + (qJD(5) * t645 + t388 * t597 + t390 * t593) * t527 - t656) * MDP(30) + (qJDD(1) * t590 + 0.2e1 * t595 * t671) * MDP(4) + (-t411 * t527 - t437 * t516 - t471 * t487 + t479 * t605 + t616) * MDP(19); (-t576 * t465 + t516 * t705 - t518 * t763 - t568 * t604 + t642) * MDP(18) + (-t566 * t454 + t435 * t524 + (t597 * t685 + t418) * t766 + t709 * t492 + (-t698 * t766 + t620) * t593 + t638) * MDP(29) + (-t426 * t568 + t429 * t576 - t471 * t474 - g(1) * (t599 * t654 + t544) - g(2) * (t596 * t654 + t542) - g(3) * (t703 + t749) + t705 * t477 - t763 * t473) * MDP(21) + (-t566 * t720 + t434 * t524 - t709 * t494 + (t593 * t685 - t417 + t675) * t766 + t612) * MDP(31) + qJDD(2) * MDP(8) + t607 + (g(3) * t595 + t598 * t617) * MDP(10) + (t568 * t435 + t704 * t492 + (-t727 + (-t468 + t685) * t766) * t597 + ((t451 - t698) * t766 + t620) * t593 + t637) * MDP(27) + (t417 * t492 - t418 * t494 + (-t494 * t685 - t735 + t434 * t566 + (-t492 * t566 - t410) * qJD(5)) * t597 + (-t492 * t685 - t435 * t566 + (t494 * t566 - t409) * qJD(5) + t669) * t593 + t663) * MDP(30) + (t474 * t518 + t568 * t587 - t588 * t705 + t610) * MDP(20) + (t474 * t516 - t763 * t588 + (-pkin(3) + t576) * t587 + t611) * MDP(19) + (t595 * t617 - t743) * MDP(9) + (-t568 * t434 + (-t675 + t712) * t766 + t704 * t494 + (-t685 * t766 - t666 + t727) * t593 + t630) * MDP(28) + (t485 * t588 + (-t518 * t701 - t587 * t594 - t588 * t673) * pkin(2) + t614) * MDP(17) + (t484 * t588 + (-t516 * t701 + t587 * t752 - t588 * t700) * pkin(2) + t613) * MDP(16) + (t392 * t524 - t410 * t417 - t409 * t418 - g(1) * t706 - g(2) * t707 - g(3) * (-qJ(6) * t725 + t658) + t709 * t425 + (-t645 * t700 - t743) * pkin(2) + t609 * t566 + t651 * (qJ(6) * t722 + t584 * t753 + t750)) * MDP(32) + MDP(6) * t689 + MDP(7) * t688 + (-t595 * t598 * MDP(4) + MDP(5) * t702) * t603; (pkin(3) * t465 - qJ(4) * t604 - t481 * t518 - t516 * t693 + t642) * MDP(18) + (t597 * t713 + t420 * t766 + t435 * t537 + t708 * t492 + (t695 * t766 + t620) * t593 + t638) * MDP(29) + (t593 * t713 + t434 * t537 + (-t419 - t674) * t766 - t708 * t494 + t612) * MDP(31) + t607 + (t613 + t732) * MDP(16) + (t478 * t516 + t611 - t732 - 0.2e1 * t748) * MDP(19) + (t478 * t518 + t588 * t693 + t578 + t610) * MDP(20) + (-qJ(4) * t434 + (t674 + t711) * t766 + t694 * t494 + (-t666 - t713) * t593 + t630) * MDP(28) + (t419 * t492 - t420 * t494 + (-t735 - t434 * t753 + (t492 * t753 - t410) * qJD(5)) * t597 + (t435 * t753 + (-t494 * t753 - t409) * qJD(5) + t669) * t593 + t663) * MDP(30) + (qJ(4) * t435 + (-t464 * t766 + t713) * t597 + t694 * t492 + ((t460 + t695) * t766 + t620) * t593 + t637) * MDP(27) + (-t480 * t588 + t614) * MDP(17) + (-t426 * qJ(4) - t429 * pkin(3) - t471 * t478 - t473 * t481 - g(1) * (-pkin(3) * t724 + t544) - g(2) * (-pkin(3) * t726 + t542) - g(3) * t703 - t693 * t477) * MDP(21) + (t392 * t537 - t410 * t419 - t409 * t420 - g(1) * (-qJ(6) * t683 + t706) - g(2) * (-qJ(6) * t684 + t707) - g(3) * t658 + (g(3) * t738 + t651 * t753) * t584 + t708 * t425 - t609 * t753) * MDP(32); (-qJD(3) * t676 - t594 * t672 + t659) * MDP(18) + (t587 - t730) * MDP(19) - t754 * MDP(20) + (t611 - t748) * MDP(21) - t454 * MDP(27) - t663 * MDP(32) + (t516 * MDP(18) - MDP(20) * t588 + t477 * MDP(21) - t425 * MDP(32) + t492 * t758 - t687 * t494) * t588 + (-t459 * MDP(29) + t770 + (MDP(32) * t410 - t687 * t766) * t766) * t597 + ((t494 * t518 - t435 + t699) * MDP(30) + (qJD(5) * t409 - t669) * MDP(32) + t687 * t459 + t758 * t667) * t593; MDP(22) * t731 + (-t492 ^ 2 + t755) * MDP(23) - t769 * MDP(24) + (-t435 + t665) * MDP(25) - t459 * MDP(26) + (-t453 * t494 + t615 + t734) * MDP(27) + (t415 * t766 + t453 * t492 - t756) * MDP(28) + (-t452 * t492 - t608 + t734 - 0.2e1 * t747) * MDP(29) + (pkin(5) * t434 - qJ(6) * t435 + (t410 - t416) * t494 + (t409 - t691) * t492) * MDP(30) + (-0.2e1 * t739 - t425 * t492 + t452 * t494 + (0.2e1 * qJD(6) - t415) * t766 + t756) * MDP(31) + (t388 * qJ(6) - t390 * pkin(5) - t425 * t452 - t409 * t416 - g(1) * (-pkin(5) * t507 + qJ(6) * t508) - g(2) * (pkin(5) * t509 - qJ(6) * t510) - t648 * t572 + t691 * t410) * MDP(32); (t459 + t731) * MDP(29) - t770 + (-t755 - t667) * MDP(31) + (-t410 * t766 + t608 + t747) * MDP(32);];
tau  = t1;
