% Calculate vector of inverse dynamics joint torques for
% S6RRPRPP2
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
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d4,theta3]';
% MDP [27x1]
%   Minimal dynamic parameter vector (fixed base model)
%   see S6RRPRPP2_convert_par2_MPV_fixb.m
% 
% Output:
% tau [6x1]
%   joint torques of inverse dynamics (contains inertial, gravitational Coriolis and centrifugal forces)

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 09:53
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tau = S6RRPRPP2_invdynJ_fixb_mdp_slag_vp(qJ, qJD, qJDD, g, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(3,1),zeros(9,1),zeros(27,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPRPP2_invdynJ_fixb_mdp_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRPRPP2_invdynJ_fixb_mdp_slag_vp: qJD has to be [6x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [6 1]), ...
  'S6RRPRPP2_invdynJ_fixb_mdp_slag_vp: qJDD has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RRPRPP2_invdynJ_fixb_mdp_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S6RRPRPP2_invdynJ_fixb_mdp_slag_vp: pkin has to be [9x1] (double)');
assert(isreal(MDP) && all(size(MDP) == [27 1]), ...
  'S6RRPRPP2_invdynJ_fixb_mdp_slag_vp: MDP has to be [27x1] (double)'); 

%% Symbolic Calculation
% From invdyn_joint_fixb_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 09:52:42
% EndTime: 2019-03-09 09:52:54
% DurationCPUTime: 8.10s
% Computational Cost: add. (6877->609), mult. (15574->707), div. (0->0), fcn. (11135->10), ass. (0->260)
t566 = sin(pkin(9));
t569 = sin(qJ(2));
t572 = cos(qJ(2));
t710 = cos(pkin(9));
t522 = t566 * t572 + t569 * t710;
t509 = t522 * qJD(1);
t636 = t710 * t572;
t655 = qJDD(1) * t569;
t615 = qJDD(1) * t636 - t566 * t655;
t463 = -qJD(2) * t509 + t615;
t544 = qJD(1) * t636;
t657 = qJD(1) * qJD(2);
t641 = t569 * t657;
t581 = qJD(2) * t544 + qJDD(1) * t522 - t566 * t641;
t652 = pkin(2) * t641 + qJDD(3);
t654 = qJDD(1) * t572;
t705 = qJDD(1) * pkin(1);
t404 = -pkin(2) * t654 - t463 * pkin(3) - pkin(8) * t581 + t652 - t705;
t567 = -qJ(3) - pkin(7);
t638 = qJD(2) * t567;
t505 = -qJD(3) * t569 + t572 * t638;
t531 = t567 * t569;
t456 = qJDD(2) * pkin(2) + qJD(1) * t505 + qJDD(1) * t531;
t504 = qJD(3) * t572 + t569 * t638;
t532 = t567 * t572;
t464 = qJD(1) * t504 - qJDD(1) * t532;
t414 = t566 * t456 + t710 * t464;
t410 = qJDD(2) * pkin(8) + t414;
t666 = qJD(1) * t569;
t507 = -t566 * t666 + t544;
t560 = t572 * pkin(2);
t554 = t560 + pkin(1);
t530 = -qJD(1) * t554 + qJD(3);
t428 = -pkin(3) * t507 - pkin(8) * t509 + t530;
t524 = qJD(1) * t531;
t711 = qJD(2) * pkin(2);
t516 = t524 + t711;
t525 = qJD(1) * t532;
t637 = t710 * t525;
t462 = t566 * t516 - t637;
t449 = qJD(2) * pkin(8) + t462;
t568 = sin(qJ(4));
t571 = cos(qJ(4));
t663 = qJD(4) * t571;
t664 = qJD(4) * t568;
t626 = -t571 * t404 + t568 * t410 + t428 * t664 + t449 * t663;
t508 = t522 * qJD(2);
t457 = qJD(1) * t508 + qJDD(4) - t615;
t452 = t457 * pkin(4);
t728 = t452 - qJDD(5);
t376 = t626 - t728;
t478 = qJD(2) * t568 + t509 * t571;
t656 = qJD(2) * qJD(4);
t597 = t568 * qJDD(2) - t509 * t664 + (t581 + t656) * t571;
t706 = qJ(6) * t597;
t720 = pkin(5) * t457;
t372 = -qJD(6) * t478 + t376 - t706 - t720;
t402 = t568 * t428 + t571 * t449;
t476 = -t571 * qJD(2) + t509 * t568;
t390 = qJ(6) * t476 + t402;
t499 = qJD(4) - t507;
t488 = t499 * qJ(5);
t385 = t390 + t488;
t704 = t385 * t499;
t745 = -t372 + t704;
t599 = t568 * t404 + t571 * t410 + t428 * t663 - t449 * t664;
t743 = t457 * qJ(5) + t499 * qJD(5);
t374 = t599 + t743;
t578 = -t571 * qJDD(2) + t568 * t581;
t421 = qJD(4) * t478 + t578;
t736 = t421 * qJ(6) + t476 * qJD(6);
t373 = t374 + t736;
t401 = t571 * t428 - t568 * t449;
t389 = qJ(6) * t478 + t401;
t660 = qJD(5) - t389;
t723 = pkin(4) + pkin(5);
t383 = -t499 * t723 + t660;
t744 = t383 * t499 + t373;
t562 = qJ(2) + pkin(9);
t555 = sin(t562);
t573 = cos(qJ(1));
t689 = t555 * t573;
t570 = sin(qJ(1));
t716 = g(2) * t570;
t742 = g(1) * t689 + t555 * t716;
t653 = MDP(20) + MDP(24);
t741 = MDP(22) + MDP(25);
t620 = g(1) * t573 + t716;
t739 = t555 * t620;
t738 = t571 * t723;
t730 = g(1) * t570 - g(2) * t573;
t737 = t730 * t555;
t735 = 0.2e1 * t743;
t465 = t524 * t566 - t637;
t734 = qJD(5) * t568 + t465;
t601 = -t566 * t569 + t636;
t732 = t508 * qJ(5) - qJD(5) * t601;
t556 = cos(t562);
t731 = -t556 * pkin(3) - t555 * pkin(8);
t659 = qJD(5) - t401;
t729 = MDP(21) - MDP(26);
t512 = t566 * t525;
t461 = t710 * t516 + t512;
t623 = qJD(2) * pkin(3) + t461;
t592 = qJ(5) * t478 + t623;
t403 = pkin(4) * t476 - t592;
t722 = pkin(2) * t566;
t549 = pkin(8) + t722;
t692 = t549 * t457;
t727 = t403 * t499 - t692;
t388 = -t476 * t723 + qJD(6) + t592;
t681 = t571 * t573;
t683 = t568 * t570;
t500 = t556 * t683 + t681;
t680 = t573 * t568;
t682 = t570 * t571;
t502 = t556 * t680 - t682;
t691 = t555 * t568;
t587 = g(1) * t502 + g(2) * t500 + g(3) * t691 - t626;
t586 = t587 + t728;
t726 = (qJD(6) + t388) * t478 + t586 + t706;
t725 = -g(3) * t555 - t620 * t556;
t724 = t476 ^ 2;
t475 = t478 ^ 2;
t721 = pkin(2) * t569;
t717 = g(2) * t567;
t714 = g(3) * t556;
t713 = g(3) * t572;
t712 = t571 * pkin(4);
t709 = qJ(5) * t421;
t708 = qJ(5) * t476;
t707 = qJ(5) * t571;
t703 = t402 * t499;
t702 = t597 * t568;
t701 = t476 * t499;
t700 = t476 * t507;
t699 = t476 * t509;
t698 = t478 * t476;
t697 = t478 * t499;
t696 = t478 * t509;
t628 = t499 * t571;
t695 = t507 * t568;
t694 = t522 * t568;
t693 = t522 * t571;
t690 = t555 * t571;
t688 = t556 * t570;
t687 = t556 * t571;
t686 = t556 * t573;
t685 = t567 * t573;
t684 = t568 * qJ(5);
t444 = t568 * t457;
t445 = t571 * t457;
t679 = qJ(6) - t549;
t678 = -t568 * t421 - t476 * t663;
t677 = t499 * t663 + t444;
t676 = t499 * t695 + t445;
t413 = t710 * t456 - t566 * t464;
t438 = pkin(2) * t666 + pkin(3) * t509 - pkin(8) * t507;
t466 = t524 * t710 + t512;
t675 = t568 * t438 + t571 * t466;
t460 = -pkin(3) * t601 - pkin(8) * t522 - t554;
t472 = t566 * t531 - t532 * t710;
t674 = t568 * t460 + t571 * t472;
t648 = t723 * t568;
t607 = -t648 + t707;
t673 = t499 * t607 + t734;
t396 = t509 * qJ(5) + t675;
t672 = -qJ(6) * t695 - qJD(6) * t571 + t664 * t679 - t396;
t458 = t568 * t466;
t520 = t679 * t571;
t671 = -qJD(4) * t520 - qJD(6) * t568 - t458 - (-qJ(6) * t507 - t438) * t571 + t723 * t509;
t616 = pkin(4) * t568 - t707;
t670 = t499 * t616 - t734;
t669 = t742 * t568;
t668 = t742 * t571;
t563 = t569 ^ 2;
t667 = -t572 ^ 2 + t563;
t665 = qJD(4) * t549;
t661 = qJD(5) * t571;
t649 = t569 * t711;
t435 = t504 * t710 + t566 * t505;
t647 = t568 * t435 + t460 * t664 + t472 * t663;
t511 = t601 * qJD(2);
t439 = pkin(3) * t508 - pkin(8) * t511 + t649;
t646 = t571 * t435 + t568 * t439 + t460 * t663;
t409 = -qJDD(2) * pkin(3) - t413;
t405 = -qJ(5) * t601 + t674;
t644 = t710 * pkin(2);
t643 = t522 * t664;
t642 = t522 * t663;
t483 = t499 * t664;
t379 = t421 * pkin(4) - qJ(5) * t597 - t478 * qJD(5) + t409;
t375 = -pkin(5) * t421 + qJDD(6) - t379;
t639 = t375 - t714;
t501 = t556 * t682 - t680;
t635 = -t500 * pkin(4) + qJ(5) * t501;
t503 = t556 * t681 + t683;
t634 = -t502 * pkin(4) + qJ(5) * t503;
t393 = -pkin(4) * t499 + t659;
t633 = -t393 * t507 + t374;
t394 = t488 + t402;
t632 = t394 * t507 + t376;
t631 = -t597 + t700;
t468 = t568 * t472;
t630 = t460 * t571 - t468;
t434 = t504 * t566 - t710 * t505;
t471 = -t710 * t531 - t532 * t566;
t627 = t499 * t568;
t625 = pkin(8) * t688 - t570 * t721;
t624 = pkin(8) * t686 - t573 * t721;
t550 = -t644 - pkin(3);
t622 = g(1) * t500 - g(2) * t502;
t621 = g(1) * t501 - g(2) * t503;
t617 = t684 + t712;
t614 = pkin(4) * t687 + t556 * t684 + t560 - t731;
t613 = t393 * t571 - t394 * t568;
t612 = -qJ(6) * t511 - qJD(6) * t522;
t610 = -t554 + t731;
t609 = t439 * t571 - t647;
t608 = -t499 * t665 - t714;
t606 = -t684 - t738;
t605 = -0.2e1 * pkin(1) * t657 - pkin(7) * qJDD(2);
t604 = -t501 * pkin(4) - qJ(5) * t500 - t685;
t603 = t511 * t568 + t642;
t602 = -t511 * t571 + t643;
t542 = t573 * t554;
t600 = pkin(3) * t686 + t503 * pkin(4) + pkin(8) * t689 + qJ(5) * t502 + t542;
t598 = -t472 * t664 + t646;
t596 = -t499 * t623 - t692;
t594 = t550 - t684;
t593 = -t379 + t608;
t590 = -qJDD(1) * t554 + t652;
t392 = t597 + t701;
t575 = qJD(2) ^ 2;
t589 = -pkin(7) * t575 + 0.2e1 * t705 + t730;
t576 = qJD(1) ^ 2;
t588 = pkin(1) * t576 - pkin(7) * qJDD(1) + t620;
t584 = t403 * t478 - t586;
t583 = g(1) * t503 + g(2) * t501 + g(3) * t690 - t599;
t580 = t401 * t499 + t583;
t577 = -t509 * t663 - t568 * t656 - t578;
t533 = qJ(5) * t690;
t519 = t679 * t568;
t518 = t594 - t712;
t489 = -t594 + t738;
t423 = pkin(4) * t478 + t708;
t422 = t522 * t616 + t471;
t412 = t522 * t607 - t471;
t411 = -t478 * t723 - t708;
t406 = pkin(4) * t601 - t630;
t397 = -pkin(4) * t509 - t438 * t571 + t458;
t395 = qJ(6) * t694 + t405;
t391 = t468 + (-qJ(6) * t522 - t460) * t571 + t723 * t601;
t384 = t616 * t511 + (qJD(4) * t617 - t661) * t522 + t434;
t382 = t607 * t511 + (qJD(4) * t606 + t661) * t522 - t434;
t381 = -pkin(4) * t508 - t609;
t380 = t598 + t732;
t378 = qJ(6) * t642 + (-qJD(4) * t472 - t612) * t568 + t646 + t732;
t377 = qJ(6) * t643 - t723 * t508 + (-t439 + t612) * t571 + t647;
t1 = [t730 * MDP(2) + (-t377 * t478 + t378 * t476 - t391 * t597 + t395 * t421 - t737 + (-t383 * t571 + t385 * t568) * t511 + (-t372 * t571 + t373 * t568 + (t383 * t568 + t385 * t571) * qJD(4)) * t522) * MDP(26) + (-t380 * t476 + t381 * t478 - t405 * t421 + t406 * t597 + t737 + t613 * t511 + (-t374 * t568 + t376 * t571 + (-t393 * t568 - t394 * t571) * qJD(4)) * t522) * MDP(21) + (-t478 * t602 + t597 * t693) * MDP(13) + (-t374 * t601 - t379 * t693 + t380 * t499 - t384 * t478 + t394 * t508 + t403 * t602 + t405 * t457 - t422 * t597 + t622) * MDP(22) + (-t373 * t601 + t375 * t693 + t378 * t499 + t382 * t478 + t385 * t508 - t388 * t602 + t395 * t457 + t412 * t597 + t622) * MDP(25) + (t445 * t522 + t478 * t508 - t499 * t602 - t597 * t601) * MDP(15) + (-t402 * t508 + t409 * t693 + t434 * t478 - t457 * t674 + t471 * t597 - t499 * t598 + t599 * t601 + t602 * t623 - t622) * MDP(19) + (t401 * t508 + t409 * t694 + t471 * t421 + t434 * t476 + t457 * t630 + t499 * t609 + t601 * t626 - t603 * t623 + t621) * MDP(18) + (-t457 * t601 + t499 * t508) * MDP(17) + (t421 * t601 - t444 * t522 - t476 * t508 - t499 * t603) * MDP(16) + (t372 * t601 - t375 * t694 - t377 * t499 - t382 * t476 - t383 * t508 - t388 * t603 - t391 * t457 - t412 * t421 + t621) * MDP(24) + (t376 * t601 + t379 * t694 - t381 * t499 + t384 * t476 - t393 * t508 + t403 * t603 - t406 * t457 + t421 * t422 + t621) * MDP(20) + t620 * MDP(3) + (qJDD(2) * t569 + t572 * t575) * MDP(6) + (qJDD(2) * t572 - t569 * t575) * MDP(7) + (t569 * t605 + t572 * t589) * MDP(9) + (-t569 * t589 + t572 * t605) * MDP(10) + (-t413 * t522 + t414 * t601 + t434 * t509 + t435 * t507 - t461 * t511 - t462 * t508 + t472 * t463 + t471 * t581 - t620) * MDP(11) + ((-t476 * t571 - t478 * t568) * t511 + (-t702 - t421 * t571 + (t476 * t568 - t478 * t571) * qJD(4)) * t522) * MDP(14) + qJDD(1) * MDP(1) + (qJDD(1) * t563 + 0.2e1 * t572 * t641) * MDP(4) + 0.2e1 * (t569 * t654 - t657 * t667) * MDP(5) + (t414 * t472 + t462 * t435 - t413 * t471 - t461 * t434 - t590 * t554 + t530 * t649 - g(1) * (-t554 * t570 - t685) - g(2) * (-t567 * t570 + t542)) * MDP(12) + (t374 * t405 + t394 * t380 + t379 * t422 + t403 * t384 + t376 * t406 + t393 * t381 - g(1) * t604 - g(2) * t600 + (-g(1) * t610 + t717) * t570) * MDP(23) + (t373 * t395 + t385 * t378 + t372 * t391 + t383 * t377 + t375 * t412 + t388 * t382 - g(1) * (-pkin(5) * t501 + t604) - g(2) * (pkin(5) * t503 - qJ(6) * t689 + t600) + (-g(1) * (qJ(6) * t555 + t610) + t717) * t570) * MDP(27); -t499 * t509 * MDP(17) + MDP(6) * t655 + MDP(7) * t654 + qJDD(2) * MDP(8) + (t569 * t588 - t713) * MDP(9) + (g(3) * t569 + t572 * t588) * MDP(10) + (t463 * t722 - t581 * t644 - (-t462 + t465) * t509 + (t461 - t466) * t507) * MDP(11) + (t461 * t465 - t462 * t466 + (t710 * t413 - t713 + t414 * t566 + (-qJD(1) * t530 + t620) * t569) * pkin(2)) * MDP(12) + (t478 * t628 + t702) * MDP(13) + ((t597 + t700) * t571 - t478 * t627 + t678) * MDP(14) + (-t507 * t628 + t677 - t696) * MDP(15) + (-t483 + t676 + t699) * MDP(16) + (-t401 * t509 + t550 * t421 + t458 * t499 - t465 * t476 + (-t714 - t409 + (-t438 - t665) * t499) * t571 + t596 * t568 + t668) * MDP(18) + (t550 * t597 + t675 * t499 + t402 * t509 - t465 * t478 + t596 * t571 + (t409 - t608) * t568 - t669) * MDP(19) + (t393 * t509 + t397 * t499 + t421 * t518 + t670 * t476 + t568 * t727 + t593 * t571 + t668) * MDP(20) + (t396 * t476 - t397 * t478 + (-t421 * t549 + (t478 * t549 + t393) * qJD(4) + t633) * t571 + (t597 * t549 + (t476 * t549 - t394) * qJD(4) + t632) * t568 + t725) * MDP(21) + (-t394 * t509 - t396 * t499 - t670 * t478 - t518 * t597 + t593 * t568 - t571 * t727 + t669) * MDP(22) + (t379 * t518 - t394 * t396 - t393 * t397 - g(1) * t624 - g(2) * t625 - g(3) * t614 + t670 * t403 + (qJD(4) * t613 + t374 * t571 + t376 * t568) * t549 + (pkin(3) + t617) * t739) * MDP(23) + (t383 * t509 - t388 * t627 - t421 * t489 + t457 * t519 - t476 * t673 - t499 * t671 + t571 * t639 + t668) * MDP(24) + (-t385 * t509 + t388 * t628 - t457 * t520 + t478 * t673 + t489 * t597 + t499 * t672 + t568 * t639 + t669) * MDP(25) + (-t421 * t520 + t672 * t476 - t671 * t478 + t597 * t519 + t568 * t745 - t744 * t571 - t725) * MDP(26) + (-t373 * t520 - t372 * t519 + t375 * t489 - g(1) * (-qJ(6) * t686 + t624) - g(2) * (-qJ(6) * t688 + t625) - g(3) * (pkin(5) * t687 + t614) + t673 * t388 + t672 * t385 + t671 * t383 + (g(3) * qJ(6) + t620 * (pkin(3) - t606)) * t555) * MDP(27) + (-MDP(4) * t569 * t572 + MDP(5) * t667) * t576; (-t507 ^ 2 - t509 ^ 2) * MDP(11) + (t461 * t509 - t462 * t507 + t590 - t730) * MDP(12) + (t676 - t699) * MDP(18) - MDP(19) * t696 + t678 * MDP(21) + (-t403 * t509 - t730) * MDP(23) + (t388 * t509 - t730) * MDP(27) + (-t457 * MDP(19) + (qJD(4) * t393 + t633) * MDP(23) + t421 * MDP(26) + t744 * MDP(27) + (-qJD(4) * MDP(18) + t478 * t729 + t653 * t507) * t499) * t568 + (t631 * MDP(21) + (qJD(4) * t394 - t632) * MDP(23) + (qJD(4) * t476 - t631) * MDP(26) + t745 * MDP(27) + t653 * t457 + (-qJD(4) * MDP(19) + (MDP(19) - t741) * t507) * t499) * t571 + t653 * (-t483 - t699) + t741 * (t677 + t696); MDP(13) * t698 + (t475 - t724) * MDP(14) + t392 * MDP(15) + (t577 + t697) * MDP(16) + t457 * MDP(17) + (t478 * t623 + t587 + t703) * MDP(18) + (-t476 * t623 + t580) * MDP(19) + (-t423 * t476 + t452 - t584 + t703) * MDP(20) + (-pkin(4) * t597 - t709 + (t394 - t402) * t478 + (t393 - t659) * t476) * MDP(21) + (-t403 * t476 + t423 * t478 - t580 + t735) * MDP(22) + (t374 * qJ(5) - t376 * pkin(4) - t403 * t423 - t393 * t402 - g(1) * t634 - g(2) * t635 - g(3) * (-pkin(4) * t691 + t533) + t659 * t394) * MDP(23) + ((pkin(5) + t723) * t457 + t390 * t499 + t411 * t476 + t726) * MDP(24) + (t388 * t476 - t389 * t499 - t411 * t478 - t583 + t735 + t736) * MDP(25) + (t709 + t597 * t723 + (-t385 + t390) * t478 + (-t383 + t660) * t476) * MDP(26) + (t373 * qJ(5) - t372 * t723 - t383 * t390 - t388 * t411 - g(1) * (-pkin(5) * t502 + t634) - g(2) * (-pkin(5) * t500 + t635) - g(3) * (-t555 * t648 + t533) + t660 * t385) * MDP(27); (-t394 * t499 + t584) * MDP(23) + (-t704 - t720 - t726) * MDP(27) + t653 * (-qJDD(4) + t463 + t698) + t729 * t392 + t741 * (-t499 ^ 2 - t475); (t577 - t697) * MDP(24) + (t597 - t701) * MDP(25) + (-t475 - t724) * MDP(26) + (t383 * t478 - t385 * t476 + t639 + t739) * MDP(27);];
tau  = t1;
