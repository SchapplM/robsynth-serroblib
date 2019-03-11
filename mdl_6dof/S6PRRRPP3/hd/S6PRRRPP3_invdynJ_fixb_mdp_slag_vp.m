% Calculate vector of inverse dynamics joint torques for
% S6PRRRPP3
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
%   pkin=[a2,a3,a4,a5,a6,alpha2,d2,d3,d4,theta1]';
% MDP [26x1]
%   Minimal dynamic parameter vector (fixed base model)
%   see S6PRRRPP3_convert_par2_MPV_fixb.m
% 
% Output:
% tau [6x1]
%   joint torques of inverse dynamics (contains inertial, gravitational Coriolis and centrifugal forces)

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-08 22:59
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tau = S6PRRRPP3_invdynJ_fixb_mdp_slag_vp(qJ, qJD, qJDD, g, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(3,1),zeros(10,1),zeros(26,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRRRPP3_invdynJ_fixb_mdp_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6PRRRPP3_invdynJ_fixb_mdp_slag_vp: qJD has to be [6x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [6 1]), ...
  'S6PRRRPP3_invdynJ_fixb_mdp_slag_vp: qJDD has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6PRRRPP3_invdynJ_fixb_mdp_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6PRRRPP3_invdynJ_fixb_mdp_slag_vp: pkin has to be [10x1] (double)');
assert(isreal(MDP) && all(size(MDP) == [26 1]), ...
  'S6PRRRPP3_invdynJ_fixb_mdp_slag_vp: MDP has to be [26x1] (double)'); 

%% Symbolic Calculation
% From invdyn_joint_fixb_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-08 22:58:46
% EndTime: 2019-03-08 22:58:57
% DurationCPUTime: 8.55s
% Computational Cost: add. (4574->616), mult. (10249->752), div. (0->0), fcn. (7789->10), ass. (0->255)
t578 = cos(qJ(3));
t769 = -pkin(3) * t578 - pkin(2);
t574 = sin(qJ(4));
t768 = qJ(5) * t574 + pkin(3);
t565 = t578 * qJDD(2);
t575 = sin(qJ(3));
t676 = qJD(2) * qJD(3);
t755 = -t575 * t676 + t565;
t516 = qJDD(4) - t755;
t738 = pkin(4) + qJ(6);
t654 = t738 * t516;
t572 = sin(pkin(6));
t576 = sin(qJ(2));
t577 = cos(qJ(4));
t579 = cos(qJ(2));
t707 = t578 * t579;
t481 = (t574 * t576 + t577 * t707) * t572;
t630 = pkin(3) * t575 - pkin(9) * t578;
t525 = t630 * qJD(3);
t530 = -pkin(9) * t575 + t769;
t682 = qJD(4) * t577;
t767 = -qJD(1) * t481 + t574 * t525 + t530 * t682;
t651 = t578 * t676;
t675 = qJDD(2) * t575;
t766 = -t651 - t675;
t692 = qJD(1) * t572;
t526 = qJD(2) * pkin(8) + t576 * t692;
t735 = cos(pkin(6));
t646 = qJD(1) * t735;
t549 = t575 * t646;
t477 = t578 * t526 + t549;
t468 = qJD(3) * pkin(9) + t477;
t663 = t579 * t692;
t479 = qJD(2) * t530 - t663;
t409 = t468 * t574 - t577 * t479;
t687 = qJD(3) * t574;
t690 = qJD(2) * t575;
t520 = t577 * t690 + t687;
t617 = pkin(5) * t520 + t409;
t679 = qJD(5) + t617;
t765 = MDP(20) - MDP(25);
t764 = MDP(21) + MDP(24);
t683 = qJD(4) * t575;
t763 = qJD(2) * t683 - qJDD(3);
t689 = qJD(2) * t578;
t439 = (qJD(3) * (qJD(4) + t689) + t675) * t574 + t763 * t577;
t708 = t577 * t578;
t558 = pkin(8) * t708;
t684 = qJD(4) * t574;
t754 = t574 * t707 - t576 * t577;
t762 = qJD(4) * t558 - t525 * t577 + t530 * t684 - t754 * t692;
t761 = qJD(5) * t578 - t767;
t760 = -t575 * t526 + t578 * t646;
t509 = t516 * qJ(5);
t554 = -qJD(4) + t689;
t542 = qJD(5) * t554;
t758 = t542 - t509;
t717 = t520 * t554;
t757 = -t439 - t717;
t711 = t574 * t578;
t674 = pkin(4) * t711;
t756 = pkin(4) * t684 - qJD(2) * t674 - qJD(5) * t574 - t477;
t753 = MDP(19) + MDP(23);
t752 = -pkin(5) * t439 + qJDD(6);
t580 = qJD(3) ^ 2;
t677 = qJD(1) * qJD(2);
t653 = t576 * t677;
t713 = t572 * t579;
t624 = -qJDD(1) * t713 + t572 * t653;
t571 = sin(pkin(10));
t734 = cos(pkin(10));
t625 = t735 * t734;
t500 = t571 * t576 - t579 * t625;
t648 = t571 * t735;
t502 = t576 * t734 + t579 * t648;
t629 = g(1) * t502 + g(2) * t500;
t749 = 0.2e1 * qJDD(2) * pkin(2) - pkin(8) * t580 + t572 * (-g(3) * t579 + t653) - t624 + t629;
t681 = t577 * qJD(3);
t518 = t574 * t690 - t681;
t747 = t518 ^ 2;
t515 = t520 ^ 2;
t551 = t554 ^ 2;
t746 = 0.2e1 * t509;
t745 = pkin(5) + pkin(9);
t743 = pkin(4) * t516;
t438 = -qJD(4) * t681 + t574 * t763 + t577 * t766;
t742 = pkin(5) * t438;
t740 = pkin(5) * t518;
t739 = pkin(9) * t516;
t737 = pkin(9) * qJD(4);
t736 = qJD(2) * pkin(2);
t733 = qJ(5) * t439;
t732 = qJ(5) * t518;
t730 = qJ(5) * t577;
t387 = t554 * t738 + t679;
t728 = t387 * t554;
t410 = t577 * t468 + t574 * t479;
t400 = qJ(5) * t554 - t410;
t727 = t400 * t554;
t726 = t410 * t554;
t725 = t438 * t574;
t501 = t571 * t579 + t576 * t625;
t647 = t572 * t734;
t452 = -t501 * t575 - t578 * t647;
t724 = t452 * t577;
t503 = -t576 * t648 + t579 * t734;
t714 = t572 * t578;
t454 = -t503 * t575 + t571 * t714;
t723 = t454 * t577;
t722 = t500 * t575;
t721 = t502 * t575;
t715 = t572 * t576;
t505 = t575 * t715 - t578 * t735;
t720 = t505 * t577;
t719 = t518 * t520;
t718 = t518 * t554;
t716 = t520 * t577;
t712 = t574 * t575;
t710 = t575 * t577;
t706 = qJDD(1) - g(3);
t655 = t574 * t683;
t664 = -pkin(8) * t574 - pkin(4);
t672 = pkin(5) * t708;
t705 = -pkin(5) * t655 + qJD(6) * t578 + (t672 + (-qJ(6) + t664) * t575) * qJD(3) + t762;
t557 = pkin(8) * t711;
t673 = pkin(5) * t711;
t704 = (-pkin(5) * t710 - t557) * qJD(4) + (-t673 + (-pkin(8) * t577 + qJ(5)) * t575) * qJD(3) - t761;
t686 = qJD(3) * t575;
t703 = -qJ(5) * t686 + (t575 * t681 + t578 * t684) * pkin(8) + t761;
t702 = t664 * t686 + t762;
t623 = qJ(6) * t574 - t730;
t613 = t623 * t578;
t701 = -qJD(2) * t613 + qJD(4) * t623 - qJD(6) * t577 + t756;
t522 = t630 * qJD(2);
t700 = t574 * t522 + t577 * t760;
t699 = -qJ(5) * t682 + t689 * t730 + t756;
t697 = -t745 * t684 - (qJ(5) * t575 - t673) * qJD(2) - t700;
t537 = t745 * t577;
t464 = t574 * t760;
t642 = -t522 * t577 + t464;
t696 = qJD(4) * t537 - (-t575 * t738 + t672) * qJD(2) - t642;
t695 = pkin(4) * t712 + t575 * pkin(8);
t694 = t574 * t530 + t558;
t569 = t575 ^ 2;
t693 = -t578 ^ 2 + t569;
t691 = qJD(2) * t572;
t688 = qJD(3) * t520;
t685 = qJD(3) * t578;
t680 = -qJD(5) - t409;
t396 = t410 - t740;
t678 = -qJD(6) - t396;
t671 = t575 * t713;
t670 = t577 * t713;
t668 = pkin(4) * t724 + t452 * t768;
t667 = pkin(4) * t723 + t454 * t768;
t666 = -pkin(4) * t720 - t505 * t768;
t662 = t576 * t691;
t661 = t579 * t691;
t660 = t554 * t687;
t659 = t554 * t681;
t467 = -qJD(3) * pkin(3) - t760;
t589 = -t520 * qJ(5) + t467;
t397 = t518 * t738 + t589;
t658 = t397 * t684;
t657 = t397 * t682;
t656 = t554 * t684;
t453 = t501 * t578 - t575 * t647;
t416 = t453 * t574 - t500 * t577;
t417 = t453 * t577 + t500 * t574;
t645 = -t416 * pkin(4) + qJ(5) * t417;
t455 = t571 * t572 * t575 + t503 * t578;
t418 = t455 * t574 - t502 * t577;
t419 = t455 * t577 + t502 * t574;
t644 = -t418 * pkin(4) + qJ(5) * t419;
t506 = t575 * t735 + t576 * t714;
t459 = t506 * t574 + t670;
t460 = t506 * t577 - t574 * t713;
t643 = -t459 * pkin(4) + qJ(5) * t460;
t641 = t530 * t577 - t557;
t640 = qJDD(1) * t735;
t483 = qJDD(2) * pkin(8) + (qJDD(1) * t576 + t579 * t677) * t572;
t627 = t575 * t640;
t404 = qJDD(3) * pkin(9) + qJD(3) * t760 + t578 * t483 + t627;
t427 = qJD(2) * t525 + qJDD(2) * t530 + t624;
t637 = t574 * t404 - t577 * t427 + t468 * t682 + t479 * t684;
t636 = -t577 * t404 - t574 * t427 + t468 * t684 - t479 * t682;
t635 = t575 * pkin(4) * t682 + pkin(8) * t685 + qJ(5) * t655 + qJD(3) * t674;
t634 = t520 * t663;
t633 = t575 * t663;
t632 = t518 * t663;
t471 = qJ(5) * t578 - t694;
t399 = pkin(4) * t554 - t680;
t620 = t399 * t577 + t400 * t574;
t619 = -qJDD(5) - t637;
t384 = t636 + t758;
t615 = t516 * t574 - t554 * t682;
t614 = t516 * t577 + t656;
t431 = -t500 * t711 - t501 * t577;
t433 = -t502 * t711 - t503 * t577;
t480 = t754 * t572;
t612 = g(1) * t433 + g(2) * t431 + g(3) * t480;
t432 = -t500 * t708 + t501 * t574;
t434 = -t502 * t708 + t503 * t574;
t611 = -g(1) * t434 - g(2) * t432 - g(3) * t481;
t610 = -g(1) * t454 - g(2) * t452 + g(3) * t505;
t609 = g(1) * t455 + g(2) * t453 + g(3) * t506;
t605 = t572 * pkin(3) * t707 + pkin(2) * t713 + t481 * pkin(4) + pkin(8) * t715 + pkin(9) * t671 + qJ(5) * t480;
t603 = -g(3) * t713 + t629;
t599 = -qJD(3) * t549 - t575 * t483 - t526 * t685 + t578 * t640;
t405 = -qJDD(3) * pkin(3) - t599;
t584 = t438 * qJ(5) - t520 * qJD(5) + t405;
t383 = t518 * qJD(6) + t439 * t738 + t584;
t602 = -t383 + t610;
t601 = -t467 * t554 - t739;
t411 = t518 * pkin(4) + t589;
t600 = t411 * t554 + t739;
t596 = t432 * pkin(4) + pkin(8) * t501 - pkin(9) * t722 + qJ(5) * t431 + t500 * t769;
t595 = t434 * pkin(4) + pkin(8) * t503 - pkin(9) * t721 + qJ(5) * t433 + t502 * t769;
t594 = t554 * t737 + t610;
t591 = g(1) * t418 + g(2) * t416 + g(3) * t459 - t637;
t590 = g(1) * t419 + g(2) * t417 + g(3) * t460 + t636;
t386 = t439 * pkin(4) + t584;
t588 = -t386 + t594;
t407 = -t438 - t718;
t527 = -t663 - t736;
t587 = -pkin(8) * qJDD(3) + (t527 + t663 - t736) * qJD(3);
t586 = -qJDD(5) + t591;
t585 = t411 * t520 - t586;
t583 = t397 * t520 - t586 - t742;
t582 = -t397 * t518 - t590 + t752;
t581 = qJD(2) ^ 2;
t568 = t578 * pkin(4);
t536 = t745 * t574;
t528 = -pkin(4) * t577 - t768;
t510 = -t577 * t738 - t768;
t495 = -qJ(5) * t710 + t695;
t472 = t568 - t641;
t466 = t575 * t623 + t695;
t458 = qJD(3) * t506 + t575 * t661;
t457 = -qJD(3) * t505 + t578 * t661;
t451 = pkin(4) * t520 + t732;
t445 = -pkin(5) * t712 - t471;
t440 = qJ(6) * t578 + t557 + t568 + (pkin(5) * t575 - t530) * t577;
t428 = t520 * t738 + t732;
t426 = (-qJ(5) * t685 - qJD(5) * t575) * t577 + t635;
t425 = -pkin(4) * t690 + t642;
t424 = -qJ(5) * t690 - t700;
t394 = qJD(3) * t613 + (qJD(6) * t574 + (qJ(6) * qJD(4) - qJD(5)) * t577) * t575 + t635;
t392 = -qJD(4) * t670 + t457 * t577 - t506 * t684 + t574 * t662;
t391 = qJD(4) * t460 + t457 * t574 - t577 * t662;
t389 = qJD(6) - t400 - t740;
t385 = -t619 - t743;
t382 = -t384 + t752;
t381 = qJD(6) * t554 - t619 - t654 - t742;
t1 = [t706 * MDP(1) + (-qJD(3) * t458 - qJDD(3) * t505) * MDP(10) + (-qJD(3) * t457 - qJDD(3) * t506) * MDP(11) + (-t384 * t460 + t385 * t459 + t386 * t505 + t391 * t399 - t392 * t400 + t411 * t458 - g(3)) * MDP(22) + (t381 * t459 + t382 * t460 + t383 * t505 + t387 * t391 + t389 * t392 + t397 * t458 - g(3)) * MDP(26) + ((-qJDD(2) * MDP(4) + (-MDP(10) * t578 + MDP(11) * t575 - MDP(3)) * t581) * t576 + (t755 * MDP(10) + MDP(11) * t766 + qJDD(2) * MDP(3) - t581 * MDP(4)) * t579) * t572 + t753 * (t391 * t520 - t392 * t518 - t459 * t438 - t439 * t460) + (-MDP(18) + t764) * (-t392 * t554 + t505 * t438 - t458 * t520 + t460 * t516) + (MDP(17) - t765) * (t391 * t554 + t505 * t439 + t458 * t518 - t459 * t516); qJDD(2) * MDP(2) + (t706 * t713 + t629) * MDP(3) + (g(1) * t503 + g(2) * t501 - t706 * t715) * MDP(4) + (qJDD(2) * t569 + 0.2e1 * t575 * t651) * MDP(5) + 0.2e1 * (t565 * t575 - t676 * t693) * MDP(6) + (qJDD(3) * t575 + t578 * t580) * MDP(7) + (qJDD(3) * t578 - t575 * t580) * MDP(8) + (t587 * t575 + t749 * t578) * MDP(10) + (-t749 * t575 + t587 * t578) * MDP(11) + (-t438 * t710 + (t578 * t681 - t655) * t520) * MDP(12) + ((-t518 * t577 - t520 * t574) * t685 + (t725 - t439 * t577 + (t518 * t574 - t716) * qJD(4)) * t575) * MDP(13) + ((t438 - t659) * t578 + (t614 + t688) * t575) * MDP(14) + ((t439 + t660) * t578 + (-qJD(3) * t518 - t615) * t575) * MDP(15) + (-t516 * t578 - t554 * t686) * MDP(16) + (t641 * t516 + t762 * t554 + ((pkin(8) * t518 + t467 * t574) * qJD(3) + t637) * t578 + (-t632 + t467 * t682 - t409 * qJD(3) + t405 * t574 + (t439 - t660) * pkin(8)) * t575 + t611) * MDP(17) + (-t694 * t516 + t767 * t554 + (t467 * t681 + (-t656 + t688) * pkin(8) - t636) * t578 + (-t634 - t467 * t684 - qJD(3) * t410 + t405 * t577 + (-t438 - t659) * pkin(8)) * t575 + t612) * MDP(18) + (-t438 * t472 + t439 * t471 + t702 * t520 + t703 * t518 + t620 * t685 + (t384 * t574 + t385 * t577 + (-t399 * t574 + t400 * t577) * qJD(4) + t603) * t575) * MDP(19) + (-t426 * t518 - t439 * t495 + t472 * t516 + (-t411 * t687 - t385) * t578 - t702 * t554 + (qJD(3) * t399 - t386 * t574 - t411 * t682 + t632) * t575 - t611) * MDP(20) + (-t426 * t520 + t438 * t495 - t471 * t516 + (-t411 * t681 + t384) * t578 + t703 * t554 + (-qJD(3) * t400 - t386 * t577 + t411 * t684 + t634) * t575 - t612) * MDP(21) + (t386 * t495 + t384 * t471 + t385 * t472 - g(1) * t595 - g(2) * t596 - g(3) * t605 + (t426 - t633) * t411 + t703 * t400 + t702 * t399) * MDP(22) + (-t438 * t440 - t439 * t445 + t705 * t520 - t704 * t518 + (t387 * t577 - t389 * t574) * t685 + (t381 * t577 - t382 * t574 + (-t387 * t574 - t389 * t577) * qJD(4) + t603) * t575) * MDP(23) + (-t394 * t520 + t438 * t466 + t445 * t516 + (-t397 * t681 - t382) * t578 - t704 * t554 + (qJD(3) * t389 - t383 * t577 + t634 + t658) * t575 - t612) * MDP(24) + (t394 * t518 + t439 * t466 - t440 * t516 + (t397 * t687 + t381) * t578 + t705 * t554 + (-qJD(3) * t387 + t383 * t574 - t632 + t657) * t575 + t611) * MDP(25) + (t383 * t466 + t381 * t440 + t382 * t445 - g(1) * (-pkin(5) * t721 + qJ(6) * t434 + t595) - g(2) * (-pkin(5) * t722 + qJ(6) * t432 + t596) - g(3) * (pkin(5) * t671 + qJ(6) * t481 + t605) + (t394 - t633) * t397 + t704 * t389 + t705 * t387) * MDP(26); MDP(7) * t675 + MDP(8) * t565 + qJDD(3) * MDP(9) + (t477 * qJD(3) - t527 * t690 + t599 + t610) * MDP(10) + (-t627 + (-qJD(2) * t527 - t483) * t578 + t609) * MDP(11) + (-t554 * t716 - t725) * MDP(12) + ((-t438 + t718) * t577 + (-t439 + t717) * t574) * MDP(13) + ((-t520 * t575 + t554 * t708) * qJD(2) + t615) * MDP(14) + ((t518 * t575 - t554 * t711) * qJD(2) + t614) * MDP(15) + t554 * MDP(16) * t690 + (t409 * t690 - pkin(3) * t439 - t464 * t554 - t477 * t518 + t601 * t574 + (-t405 + (t522 + t737) * t554 + t610) * t577) * MDP(17) + (pkin(3) * t438 - t700 * t554 + t410 * t690 - t477 * t520 + t601 * t577 + (t405 - t594) * t574) * MDP(18) + (-t424 * t518 - t425 * t520 + (-t384 - t554 * t399 + (qJD(4) * t520 - t439) * pkin(9)) * t577 + (t385 - t727 + (qJD(4) * t518 - t438) * pkin(9)) * t574 - t609) * MDP(19) + (-t399 * t690 + t425 * t554 - t439 * t528 - t518 * t699 + t574 * t600 - t577 * t588) * MDP(20) + (t400 * t690 - t424 * t554 + t438 * t528 - t520 * t699 + t574 * t588 + t577 * t600) * MDP(21) + (t386 * t528 - t400 * t424 - t399 * t425 - g(1) * t667 - g(2) * t668 - g(3) * t666 + t699 * t411 + (qJD(4) * t620 - t384 * t577 + t385 * t574 - t609) * pkin(9)) * MDP(22) + (-t438 * t536 - t439 * t537 + t696 * t520 - t697 * t518 + (t382 - t728) * t577 + (t389 * t554 + t381) * t574 - t609) * MDP(23) + (-t657 + t438 * t510 + t516 * t537 - t697 * t554 - t701 * t520 + (-t389 * t575 + t397 * t708) * qJD(2) + t602 * t574) * MDP(24) + (t658 + t439 * t510 - t516 * t536 + t696 * t554 + t701 * t518 + (t387 * t575 - t397 * t711) * qJD(2) + t602 * t577) * MDP(25) + (t383 * t510 + t381 * t536 + t382 * t537 - g(1) * (qJ(6) * t723 + t455 * t745 + t667) - g(2) * (qJ(6) * t724 + t453 * t745 + t668) - g(3) * (-qJ(6) * t720 + t506 * t745 + t666) + t701 * t397 + t697 * t389 + t696 * t387) * MDP(26) + (-MDP(5) * t575 * t578 + MDP(6) * t693) * t581; MDP(12) * t719 + (t515 - t747) * MDP(13) + t407 * MDP(14) + t757 * MDP(15) + t516 * MDP(16) + (-t467 * t520 + t591 - t726) * MDP(17) + (t409 * t554 + t467 * t518 + t590) * MDP(18) + (pkin(4) * t438 - t733 + (-t400 - t410) * t520 + (t399 + t680) * t518) * MDP(19) + (t451 * t518 + t585 + t726 - 0.2e1 * t743) * MDP(20) + (-t411 * t518 + t451 * t520 + t554 * t680 - t542 - t590 + t746) * MDP(21) + (-t385 * pkin(4) - g(1) * t644 - g(2) * t645 - g(3) * t643 - t384 * qJ(5) - t399 * t410 + t400 * t680 - t411 * t451) * MDP(22) + (-t733 + t438 * t738 + (t389 + t678) * t520 + (t387 - t679) * t518) * MDP(23) + (t428 * t520 - t554 * t617 - 0.2e1 * t542 + t582 + t746) * MDP(24) + (-t428 * t518 + (-0.2e1 * qJD(6) - t396) * t554 + 0.2e1 * t654 - t583) * MDP(25) + (-t381 * t738 + t382 * qJ(5) - t397 * t428 - g(1) * (-qJ(6) * t418 + t644) - g(2) * (-qJ(6) * t416 + t645) - g(3) * (-qJ(6) * t459 + t643) + t679 * t389 + t678 * t387) * MDP(26); (t585 - t727 - t743) * MDP(22) + ((qJD(6) + t389) * t554 - t654 + t583) * MDP(26) + t765 * (t516 - t719) + t753 * t407 + t764 * (-t515 - t551); (t516 + t719) * MDP(24) + (-t551 - t747) * MDP(25) + (t582 - t728 - t758) * MDP(26) + t757 * MDP(23);];
tau  = t1;
