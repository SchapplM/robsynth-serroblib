% Calculate vector of inverse dynamics joint torques for
% S6PRRPRR5
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
% pkin [12x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d2,d3,d5,d6,theta1,theta4]';
% MDP [29x1]
%   Minimal dynamic parameter vector (fixed base model)
%   see S6PRRPRR5_convert_par2_MPV_fixb.m
% 
% Output:
% tau [6x1]
%   joint torques of inverse dynamics (contains inertial, gravitational Coriolis and centrifugal forces)

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-08 22:21
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tau = S6PRRPRR5_invdynJ_fixb_mdp_slag_vp(qJ, qJD, qJDD, g, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(3,1),zeros(12,1),zeros(29,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRRPRR5_invdynJ_fixb_mdp_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6PRRPRR5_invdynJ_fixb_mdp_slag_vp: qJD has to be [6x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [6 1]), ...
  'S6PRRPRR5_invdynJ_fixb_mdp_slag_vp: qJDD has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6PRRPRR5_invdynJ_fixb_mdp_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [12 1]), ...
  'S6PRRPRR5_invdynJ_fixb_mdp_slag_vp: pkin has to be [12x1] (double)');
assert(isreal(MDP) && all(size(MDP) == [29 1]), ...
  'S6PRRPRR5_invdynJ_fixb_mdp_slag_vp: MDP has to be [29x1] (double)'); 

%% Symbolic Calculation
% From invdyn_joint_fixb_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-08 22:21:06
% EndTime: 2019-03-08 22:21:26
% DurationCPUTime: 14.18s
% Computational Cost: add. (5826->632), mult. (13557->873), div. (0->0), fcn. (11063->18), ass. (0->251)
t621 = cos(qJ(3));
t704 = qJD(2) * t621;
t594 = -qJD(5) + t704;
t589 = -qJD(6) + t594;
t611 = sin(pkin(12));
t613 = cos(pkin(12));
t694 = t613 * qJD(3);
t617 = sin(qJ(3));
t705 = qJD(2) * t617;
t562 = t611 * t705 - t694;
t703 = qJD(3) * t611;
t564 = t613 * t705 + t703;
t616 = sin(qJ(5));
t620 = cos(qJ(5));
t487 = t620 * t562 + t564 * t616;
t619 = cos(qJ(6));
t486 = t562 * t616 - t564 * t620;
t615 = sin(qJ(6));
t738 = t486 * t615;
t763 = -t619 * t487 + t738;
t762 = t589 * t763;
t656 = pkin(3) * t617 - qJ(4) * t621;
t540 = qJD(3) * t656 - qJD(4) * t617;
t618 = sin(qJ(2));
t702 = qJD(3) * t617;
t688 = pkin(8) * t702;
t612 = sin(pkin(6));
t709 = qJD(1) * t612;
t622 = cos(qJ(2));
t724 = t621 * t622;
t716 = t613 * t540 + t611 * t688 - (-t611 * t724 + t613 * t618) * t709;
t767 = t611 * t540 - (t611 * t618 + t613 * t724) * t709;
t572 = t656 * qJD(2);
t574 = qJD(2) * pkin(8) + t618 * t709;
t614 = cos(pkin(6));
t726 = t614 * t621;
t752 = qJD(1) * t726 - t617 * t574;
t470 = t611 * t572 + t613 * t752;
t683 = t611 * t704;
t458 = -pkin(9) * t683 + t470;
t766 = qJD(4) * t613 - t458;
t727 = t613 * t621;
t648 = pkin(4) * t617 - pkin(9) * t727;
t765 = qJD(3) * t648 + t716;
t728 = t613 * t617;
t731 = t611 * t621;
t764 = (-pkin(8) * t728 - pkin(9) * t731) * qJD(3) + t767;
t761 = t486 * t594;
t760 = t487 * t594;
t651 = t486 * t619 + t487 * t615;
t759 = t589 * t651;
t733 = t611 * t616;
t570 = -t620 * t613 + t733;
t644 = t621 * t570;
t714 = qJD(2) * t644 - t570 * qJD(5);
t571 = t611 * t620 + t613 * t616;
t643 = t571 * t621;
t713 = -qJD(2) * t643 + t571 * qJD(5);
t708 = qJD(1) * t617;
t590 = t614 * t708;
t522 = t621 * t574 + t590;
t513 = qJD(3) * qJ(4) + t522;
t647 = pkin(3) * t621 + qJ(4) * t617 + pkin(2);
t707 = qJD(1) * t622;
t686 = t612 * t707;
t523 = -qJD(2) * t647 - t686;
t449 = -t513 * t611 + t613 * t523;
t425 = -pkin(4) * t704 - pkin(9) * t564 + t449;
t450 = t613 * t513 + t611 * t523;
t433 = -pkin(9) * t562 + t450;
t405 = t425 * t616 + t433 * t620;
t397 = -pkin(10) * t487 + t405;
t696 = qJD(6) * t615;
t395 = t397 * t696;
t509 = -qJD(3) * pkin(3) + qJD(4) - t752;
t471 = pkin(4) * t562 + t509;
t424 = pkin(5) * t487 + t471;
t743 = sin(pkin(11));
t669 = t743 * t622;
t744 = cos(pkin(11));
t672 = t744 * t618;
t545 = t614 * t672 + t669;
t674 = t612 * t744;
t503 = t545 * t621 - t617 * t674;
t670 = t743 * t618;
t671 = t744 * t622;
t547 = -t614 * t670 + t671;
t673 = t612 * t743;
t505 = t547 * t621 + t617 * t673;
t544 = -t614 * t671 + t670;
t546 = t614 * t669 + t672;
t730 = t612 * t618;
t552 = t614 * t617 + t621 * t730;
t608 = pkin(12) + qJ(5);
t606 = qJ(6) + t608;
t597 = sin(t606);
t598 = cos(t606);
t729 = t612 * t622;
t757 = -t424 * t763 - g(1) * (-t505 * t598 - t546 * t597) - g(2) * (-t503 * t598 - t544 * t597) - g(3) * (-t552 * t598 + t597 * t729) + t395;
t601 = t613 * qJDD(3);
t692 = qJD(2) * qJD(3);
t677 = t621 * t692;
t690 = qJDD(2) * t617;
t642 = t677 + t690;
t524 = t611 * t642 - t601;
t689 = qJDD(3) * t611;
t525 = t613 * t642 + t689;
t697 = qJD(5) * t620;
t699 = qJD(5) * t616;
t419 = -t616 * t524 + t620 * t525 - t562 * t697 - t564 * t699;
t605 = t621 * qJDD(2);
t641 = t617 * t692 - t605;
t569 = qJDD(5) + t641;
t693 = qJD(1) * qJD(2);
t533 = qJDD(2) * pkin(8) + (qJDD(1) * t618 + t622 * t693) * t612;
t691 = qJDD(1) * t614;
t676 = t617 * t691;
t437 = t676 + qJDD(3) * qJ(4) + t533 * t621 + (qJD(4) + t752) * qJD(3);
t679 = t618 * t693;
t655 = -qJDD(1) * t729 + t612 * t679;
t457 = qJD(2) * t540 - qJDD(2) * t647 + t655;
t412 = -t437 * t611 + t613 * t457;
t403 = pkin(4) * t641 - pkin(9) * t525 + t412;
t413 = t613 * t437 + t611 * t457;
t407 = -pkin(9) * t524 + t413;
t667 = t620 * t403 - t616 * t407;
t627 = -qJD(5) * t405 + t667;
t386 = pkin(5) * t569 - pkin(10) * t419 + t627;
t420 = -qJD(5) * t486 + t620 * t524 + t525 * t616;
t640 = -t616 * t403 - t620 * t407 - t425 * t697 + t433 * t699;
t387 = -pkin(10) * t420 - t640;
t668 = t619 * t386 - t615 * t387;
t756 = t424 * t651 - g(1) * (-t505 * t597 + t546 * t598) - g(2) * (-t503 * t597 + t544 * t598) - g(3) * (-t552 * t597 - t598 * t729) + t668;
t555 = qJDD(6) + t569;
t755 = t555 * MDP(27) + (t651 ^ 2 - t763 ^ 2) * MDP(24) + t763 * MDP(23) * t651;
t754 = t616 * t764 - t620 * t765;
t561 = t613 * t647;
t493 = -pkin(9) * t728 - t561 + (-pkin(8) * t611 - pkin(4)) * t621;
t527 = pkin(8) * t727 - t611 * t647;
t732 = t611 * t617;
t508 = -pkin(9) * t732 + t527;
t753 = t493 * t697 - t508 * t699 + t616 * t765 + t620 * t764;
t717 = t616 * t493 + t620 * t508;
t746 = pkin(9) + qJ(4);
t578 = t746 * t611;
t579 = t746 * t613;
t712 = -t616 * t578 + t620 * t579;
t469 = t613 * t572 - t611 * t752;
t448 = qJD(2) * t648 + t469;
t649 = qJD(4) * t611 + qJD(5) * t579;
t751 = -t578 * t697 + t766 * t620 + (-t448 - t649) * t616;
t750 = -qJDD(3) * pkin(3) + qJDD(4);
t623 = qJD(3) ^ 2;
t747 = g(2) * t544;
t748 = g(1) * t546;
t660 = t747 + t748;
t749 = 0.2e1 * qJDD(2) * pkin(2) - pkin(8) * t623 + t612 * (-g(3) * t622 + t679) - t655 + t660;
t666 = t615 * t419 + t619 * t420;
t391 = -qJD(6) * t651 + t666;
t745 = qJD(2) * pkin(2);
t404 = t620 * t425 - t433 * t616;
t396 = pkin(10) * t486 + t404;
t394 = -pkin(5) * t594 + t396;
t740 = t394 * t619;
t739 = t397 * t619;
t737 = t597 * t621;
t736 = t598 * t621;
t603 = sin(t608);
t735 = t603 * t621;
t604 = cos(t608);
t734 = t604 * t621;
t725 = t615 * t386;
t723 = qJDD(1) - g(3);
t698 = qJD(5) * t617;
t474 = -qJD(3) * t644 - t571 * t698;
t722 = -pkin(5) * t702 + pkin(10) * t474 + qJD(5) * t717 + t754;
t475 = qJD(3) * t643 + t697 * t728 - t698 * t733;
t721 = -pkin(10) * t475 + t753;
t491 = t619 * t570 + t571 * t615;
t720 = -qJD(6) * t491 - t615 * t713 + t619 * t714;
t492 = -t570 * t615 + t571 * t619;
t719 = qJD(6) * t492 + t615 * t714 + t619 * t713;
t715 = -t613 * t688 + t767;
t701 = qJD(3) * t621;
t556 = (pkin(4) * t611 + pkin(8)) * t701;
t573 = pkin(4) * t732 + t617 * pkin(8);
t609 = t617 ^ 2;
t711 = -t621 ^ 2 + t609;
t706 = qJD(2) * t612;
t695 = qJD(6) * t619;
t687 = t619 * t419 - t615 * t420 - t487 * t695;
t496 = pkin(4) * t683 + t522;
t599 = -pkin(4) * t613 - pkin(3);
t684 = t617 * t707;
t682 = t618 * t706;
t681 = t622 * t706;
t680 = qJ(4) * t605;
t675 = pkin(5) * t713 - t496;
t664 = t620 * t493 - t508 * t616;
t662 = -t620 * t578 - t579 * t616;
t661 = qJD(6) * t394 + t387;
t659 = g(1) * t547 + g(2) * t545;
t442 = t620 * t448;
t473 = -pkin(10) * t570 + t712;
t658 = pkin(5) * t705 + pkin(10) * t714 + t571 * qJD(4) + qJD(5) * t712 + qJD(6) * t473 - t458 * t616 + t442;
t472 = -pkin(10) * t571 + t662;
t657 = -pkin(10) * t713 + qJD(6) * t472 + t751;
t389 = t394 * t615 + t739;
t538 = t570 * t617;
t414 = -pkin(5) * t621 + pkin(10) * t538 + t664;
t537 = t571 * t617;
t418 = -pkin(10) * t537 + t717;
t654 = t414 * t615 + t418 * t619;
t500 = -t552 * t611 - t613 * t729;
t501 = t552 * t613 - t611 * t729;
t435 = t500 * t620 - t501 * t616;
t436 = t500 * t616 + t501 * t620;
t653 = t435 * t619 - t436 * t615;
t652 = t435 * t615 + t436 * t619;
t464 = t619 * t537 - t538 * t615;
t465 = -t537 * t615 - t538 * t619;
t551 = t617 * t730 - t726;
t390 = t486 * t696 + t687;
t638 = g(1) * (t547 * t617 - t621 * t673) + g(2) * (t545 * t617 + t621 * t674) + g(3) * t551;
t637 = g(1) * t505 + g(2) * t503 + g(3) * t552;
t636 = qJD(3) * t590 + t617 * t533 + t574 * t701 - t621 * t691;
t635 = g(3) * t729 - t660;
t634 = -g(3) * t730 - t659;
t440 = t636 + t750;
t632 = -t440 + t638;
t631 = -qJ(4) * t702 + (qJD(4) - t509) * t621;
t575 = -t686 - t745;
t626 = -pkin(8) * qJDD(3) + (t575 + t686 - t745) * qJD(3);
t421 = pkin(4) * t524 + t440;
t625 = -t636 + t638;
t624 = qJD(2) ^ 2;
t531 = pkin(5) * t570 + t599;
t526 = -pkin(8) * t731 - t561;
t507 = qJD(3) * t552 + t617 * t681;
t506 = -qJD(3) * t551 + t621 * t681;
t495 = pkin(5) * t537 + t573;
t468 = t506 * t613 + t611 * t682;
t467 = -t506 * t611 + t613 * t682;
t447 = pkin(5) * t475 + t556;
t411 = qJD(6) * t465 + t474 * t615 + t619 * t475;
t410 = -qJD(6) * t464 + t474 * t619 - t475 * t615;
t402 = -qJD(5) * t436 + t467 * t620 - t468 * t616;
t401 = qJD(5) * t435 + t467 * t616 + t468 * t620;
t398 = pkin(5) * t420 + t421;
t388 = -t397 * t615 + t740;
t1 = [t723 * MDP(1) + (-qJD(3) * t507 - qJDD(3) * t551) * MDP(10) + (-qJD(3) * t506 - qJDD(3) * t552) * MDP(11) + (-t500 * t605 + t507 * t562 + t524 * t551) * MDP(12) + (t501 * t605 + t507 * t564 + t525 * t551) * MDP(13) + (-t467 * t564 - t468 * t562 - t500 * t525 - t501 * t524) * MDP(14) + (t412 * t500 + t413 * t501 + t440 * t551 + t449 * t467 + t450 * t468 + t507 * t509 - g(3)) * MDP(15) + (-t402 * t594 + t420 * t551 + t435 * t569 + t487 * t507) * MDP(21) + (t401 * t594 + t419 * t551 - t436 * t569 - t486 * t507) * MDP(22) + (-(-qJD(6) * t652 - t401 * t615 + t402 * t619) * t589 + t653 * t555 - t507 * t763 + t551 * t391) * MDP(28) + ((qJD(6) * t653 + t401 * t619 + t402 * t615) * t589 - t652 * t555 - t507 * t651 + t551 * t390) * MDP(29) + ((-t467 * t621 + t500 * t702) * MDP(12) + (t468 * t621 - t501 * t702) * MDP(13)) * qJD(2) + ((-qJDD(2) * MDP(4) + (-MDP(10) * t621 + MDP(11) * t617 - MDP(3)) * t624) * t618 + (-MDP(10) * t641 - MDP(11) * t642 + qJDD(2) * MDP(3) - t624 * MDP(4)) * t622) * t612; (qJDD(2) * t609 + 0.2e1 * t617 * t677) * MDP(5) + (t626 * t617 + t621 * t749) * MDP(10) + (-t617 * t749 + t626 * t621) * MDP(11) + (t412 * t526 + t413 * t527 + t715 * t450 + t716 * t449 + t647 * t748 + t647 * t747 + (t440 * t617 + t509 * t701 - t659) * pkin(8) + (-g(3) * pkin(8) * t618 + (-g(3) * t647 - t509 * t708) * t622) * t612) * MDP(15) + (t664 * t569 - t667 * t621 + t404 * t702 + t556 * t487 + t573 * t420 + t421 * t537 + t471 * t475 - g(1) * (-t546 * t734 + t547 * t603) - g(2) * (-t544 * t734 + t545 * t603) + t754 * t594 + (t405 * t621 + t594 * t717) * qJD(5) + (-t487 * t684 - g(3) * (t603 * t618 + t604 * t724)) * t612) * MDP(21) + (t634 * t613 + (-t564 * t686 + pkin(8) * t525 + t440 * t613 + (-qJD(2) * t527 - t450) * qJD(3)) * t617 + (t527 * qJDD(2) + t413 + (pkin(8) * t564 + t509 * t613) * qJD(3) + t715 * qJD(2) + t635 * t611) * t621) * MDP(13) + (-t524 * t527 - t525 * t526 - t716 * t564 - t715 * t562 + (-t449 * t613 - t450 * t611) * t701 + (-t412 * t613 - t413 * t611 - t635) * t617) * MDP(14) + (t634 * t611 + (-t562 * t686 + pkin(8) * t524 + t440 * t611 + (qJD(2) * t526 + t449) * qJD(3)) * t617 + (-t526 * qJDD(2) - t412 + (pkin(8) * t562 + t509 * t611) * qJD(3) - t716 * qJD(2) - t635 * t613) * t621) * MDP(12) + ((t414 * t619 - t418 * t615) * t555 - t668 * t621 + t388 * t702 - t447 * t763 + t495 * t391 + t398 * t464 + t424 * t411 - g(1) * (-t546 * t736 + t547 * t597) - g(2) * (-t544 * t736 + t545 * t597) + (t615 * t721 + t619 * t722) * t589 + (t389 * t621 + t589 * t654) * qJD(6) + (t763 * t684 - g(3) * (t597 * t618 + t598 * t724)) * t612) * MDP(28) + (t391 * t621 + t411 * t589 - t464 * t555 + t702 * t763) * MDP(26) + (-t390 * t464 - t391 * t465 + t410 * t763 + t411 * t651) * MDP(24) + 0.2e1 * (t605 * t617 - t692 * t711) * MDP(6) + (-t569 * t621 - t594 * t702) * MDP(20) + (-t555 * t621 - t589 * t702) * MDP(27) + (t420 * t621 + t475 * t594 - t487 * t702 - t537 * t569) * MDP(19) + (t723 * t729 + t660) * MDP(3) + (-t723 * t730 + t659) * MDP(4) + qJDD(2) * MDP(2) + (-t654 * t555 + (t661 * t619 - t395 + t725) * t621 - t389 * t702 - t447 * t651 + t495 * t390 + t398 * t465 + t424 * t410 - g(1) * (t546 * t737 + t547 * t598) - g(2) * (t544 * t737 + t545 * t598) + ((qJD(6) * t414 + t721) * t619 + (-qJD(6) * t418 - t722) * t615) * t589 + (t651 * t684 - g(3) * (-t597 * t724 + t598 * t618)) * t612) * MDP(29) + (-t390 * t621 - t410 * t589 + t465 * t555 - t651 * t702) * MDP(25) + (t390 * t465 - t410 * t651) * MDP(23) + (-t717 * t569 - t640 * t621 - t405 * t702 - t556 * t486 + t573 * t419 - t421 * t538 + t471 * t474 - g(1) * (t546 * t735 + t547 * t604) - g(2) * (t544 * t735 + t545 * t604) + t753 * t594 + (t486 * t684 - g(3) * (-t603 * t724 + t604 * t618)) * t612) * MDP(22) + (-t419 * t621 - t474 * t594 - t486 * t702 - t538 * t569) * MDP(18) + (-t419 * t538 - t474 * t486) * MDP(16) + (-t419 * t537 + t420 * t538 - t474 * t487 + t475 * t486) * MDP(17) + (qJDD(3) * t617 + t621 * t623) * MDP(7) + (qJDD(3) * t621 - t617 * t623) * MDP(8); (-t676 + (-qJD(2) * t575 - t533) * t621 + t637) * MDP(11) + (t599 * t419 + t421 * t571 + t714 * t471 + t486 * t496 - t712 * t569 + t594 * t751 - t638 * t603) * MDP(22) + ((t472 * t619 - t473 * t615) * t555 + t531 * t391 + t398 * t491 + (t615 * t657 + t619 * t658) * t589 - t675 * t763 + t719 * t424 + t638 * t598) * MDP(28) + (-t491 * t555 + t589 * t719) * MDP(26) + (t390 * t492 - t651 * t720) * MDP(23) + (-(t472 * t615 + t473 * t619) * t555 + t531 * t390 + t398 * t492 + (-t615 * t658 + t619 * t657) * t589 - t675 * t651 + t720 * t424 - t638 * t597) * MDP(29) + (-t390 * t491 - t391 * t492 + t651 * t719 + t720 * t763) * MDP(24) + (t492 * t555 - t589 * t720) * MDP(25) + (t662 * t569 + t599 * t420 + t421 * t570 - t496 * t487 + (t442 + t649 * t620 + (-qJD(5) * t578 + t766) * t616) * t594 + t713 * t471 + t638 * t604) * MDP(21) + (-t569 * t570 + t594 * t713) * MDP(19) + (t419 * t571 - t486 * t714) * MDP(16) + (-t419 * t570 - t420 * t571 + t486 * t713 - t487 * t714) * MDP(17) + (t569 * t571 - t594 * t714) * MDP(18) + (qJD(3) * t522 + t625) * MDP(10) + (t469 * t564 + t470 * t562 + (-qJ(4) * t524 - qJD(4) * t562 + t449 * t704 + t413) * t613 + (qJ(4) * t525 + qJD(4) * t564 + t450 * t704 - t412) * t611 - t637) * MDP(14) + (t613 * t680 - pkin(3) * t525 - t522 * t564 - t632 * t611 + (t450 * t617 - t470 * t621 + t613 * t631) * qJD(2)) * MDP(13) + (t611 * t680 - pkin(3) * t524 - t522 * t562 + t632 * t613 + (-t449 * t617 + t469 * t621 + t611 * t631) * qJD(2)) * MDP(12) + (-t449 * t469 - t450 * t470 - t509 * t522 + (-t449 * t611 + t450 * t613) * qJD(4) + t632 * pkin(3) + (-t412 * t611 + t413 * t613 - t637) * qJ(4)) * MDP(15) + qJDD(3) * MDP(9) + MDP(8) * t605 + MDP(7) * t690 + (-t575 * MDP(10) + MDP(18) * t486 + t487 * MDP(19) + t594 * MDP(20) - t404 * MDP(21) + t405 * MDP(22) + MDP(25) * t651 - MDP(26) * t763 + t589 * MDP(27) - t388 * MDP(28) + t389 * MDP(29)) * t705 + (-MDP(5) * t617 * t621 + MDP(6) * t711) * t624; (t611 * t690 - t601 + (-t564 + t703) * t704) * MDP(12) + (t613 * t690 + t689 + (t562 + t694) * t704) * MDP(13) + (-t562 ^ 2 - t564 ^ 2) * MDP(14) + (t449 * t564 + t450 * t562 - t625 + t750) * MDP(15) + (t420 + t761) * MDP(21) + (t419 + t760) * MDP(22) + (t391 + t759) * MDP(28) + (t390 - t762) * MDP(29); -t486 * t487 * MDP(16) + (t486 ^ 2 - t487 ^ 2) * MDP(17) + (t419 - t760) * MDP(18) + (-t420 + t761) * MDP(19) + t569 * MDP(20) + (-t405 * t594 + t471 * t486 - g(1) * (-t505 * t603 + t546 * t604) - g(2) * (-t503 * t603 + t544 * t604) - g(3) * (-t552 * t603 - t604 * t729) + t627) * MDP(21) + (-t404 * t594 + t471 * t487 - g(1) * (-t505 * t604 - t546 * t603) - g(2) * (-t503 * t604 - t544 * t603) - g(3) * (-t552 * t604 + t603 * t729) + t640) * MDP(22) + (t390 + t762) * MDP(25) + (-t391 + t759) * MDP(26) + ((-t396 * t615 - t739) * t589 - t389 * qJD(6) + (-t486 * t763 + t555 * t619 + t589 * t696) * pkin(5) + t756) * MDP(28) + ((t397 * t589 - t386) * t615 + (-t396 * t589 - t661) * t619 + (-t486 * t651 - t555 * t615 + t589 * t695) * pkin(5) + t757) * MDP(29) + t755; (t687 + t762) * MDP(25) + (-t666 + t759) * MDP(26) + (-t389 * t589 + t756) * MDP(28) + (-t619 * t387 - t388 * t589 - t725 + t757) * MDP(29) + (MDP(25) * t738 + MDP(26) * t651 - MDP(28) * t389 - MDP(29) * t740) * qJD(6) + t755;];
tau  = t1;
