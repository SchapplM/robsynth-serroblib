% Calculate vector of inverse dynamics joint torques for
% S6PRRPRR7
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
%   pkin=[a2,a3,a4,a5,a6,alpha2,d2,d3,d5,d6,theta1]';
% MDP [29x1]
%   Minimal dynamic parameter vector (fixed base model)
%   see S6PRRPRR7_convert_par2_MPV_fixb.m
% 
% Output:
% tau [6x1]
%   joint torques of inverse dynamics (contains inertial, gravitational Coriolis and centrifugal forces)

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-08 22:34
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tau = S6PRRPRR7_invdynJ_fixb_mdp_slag_vp(qJ, qJD, qJDD, g, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(3,1),zeros(11,1),zeros(29,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRRPRR7_invdynJ_fixb_mdp_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6PRRPRR7_invdynJ_fixb_mdp_slag_vp: qJD has to be [6x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [6 1]), ...
  'S6PRRPRR7_invdynJ_fixb_mdp_slag_vp: qJDD has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6PRRPRR7_invdynJ_fixb_mdp_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6PRRPRR7_invdynJ_fixb_mdp_slag_vp: pkin has to be [11x1] (double)');
assert(isreal(MDP) && all(size(MDP) == [29 1]), ...
  'S6PRRPRR7_invdynJ_fixb_mdp_slag_vp: MDP has to be [29x1] (double)'); 

%% Symbolic Calculation
% From invdyn_joint_fixb_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-08 22:34:28
% EndTime: 2019-03-08 22:34:39
% DurationCPUTime: 9.65s
% Computational Cost: add. (3310->558), mult. (7247->752), div. (0->0), fcn. (5448->14), ass. (0->259)
t553 = sin(qJ(3));
t654 = qJD(4) * t553;
t557 = cos(qJ(3));
t655 = qJD(3) * t557;
t586 = -qJ(4) * t655 - t654;
t549 = sin(pkin(6));
t554 = sin(qJ(2));
t646 = qJD(1) * qJD(2);
t626 = t554 * t646;
t558 = cos(qJ(2));
t687 = t549 * t558;
t603 = -qJDD(1) * t687 + t549 * t626;
t645 = qJD(2) * qJD(3);
t623 = t553 * t645;
t589 = pkin(3) * t623 + t603;
t621 = -qJ(4) * t553 - pkin(2);
t591 = pkin(3) * t557 - t621;
t733 = qJDD(2) * t591;
t406 = qJD(2) * t586 + t589 - t733;
t741 = -g(3) * t687 - t406;
t625 = t558 * t646;
t550 = cos(pkin(6));
t659 = qJD(3) * t550;
t703 = qJDD(2) * pkin(8);
t739 = t703 + (qJDD(1) * t554 + t625) * t549 + qJD(1) * t659;
t642 = qJDD(2) * t557;
t736 = t623 - t642;
t707 = cos(pkin(11));
t617 = t707 * t558;
t548 = sin(pkin(11));
t691 = t548 * t554;
t471 = -t550 * t617 + t691;
t618 = t707 * t554;
t690 = t548 * t558;
t473 = t550 * t690 + t618;
t610 = g(1) * t473 + g(2) * t471;
t666 = qJD(1) * t558;
t715 = pkin(4) + pkin(8);
t556 = cos(qJ(5));
t552 = sin(qJ(5));
t658 = qJD(3) * t552;
t663 = qJD(2) * t557;
t492 = t556 * t663 + t658;
t555 = cos(qJ(6));
t627 = t552 * t663;
t656 = qJD(3) * t556;
t494 = -t627 + t656;
t551 = sin(qJ(6));
t695 = t494 * t551;
t421 = t555 * t492 + t695;
t664 = qJD(2) * t553;
t530 = qJD(5) + t664;
t521 = qJD(6) + t530;
t738 = t421 * t521;
t598 = t492 * t551 - t555 * t494;
t737 = t521 * t598;
t495 = t551 * t556 + t552 * t555;
t585 = t553 * t495;
t723 = qJD(5) + qJD(6);
t675 = -qJD(2) * t585 - t723 * t495;
t735 = qJD(2) * t591;
t689 = t549 * t554;
t635 = qJD(1) * t689;
t502 = qJD(2) * pkin(8) + t635;
t667 = qJD(1) * t550;
t453 = t553 * t502 - t557 * t667;
t734 = qJD(4) + t453;
t641 = MDP(12) * qJDD(2);
t559 = -pkin(3) - pkin(9);
t647 = pkin(4) * t664 + t734;
t419 = qJD(3) * t559 + t647;
t487 = t557 * t559 + t621;
t634 = t549 * t666;
t442 = qJD(2) * t487 - t634;
t393 = t419 * t552 + t442 * t556;
t385 = -pkin(10) * t492 + t393;
t649 = qJD(6) * t551;
t383 = t385 * t649;
t454 = t557 * t502 + t553 * t667;
t441 = pkin(4) * t663 + t454;
t544 = qJD(3) * qJ(4);
t427 = t544 + t441;
t407 = pkin(5) * t492 + t427;
t472 = t550 * t618 + t690;
t619 = t549 * t707;
t432 = t472 * t553 + t557 * t619;
t474 = -t550 * t691 + t617;
t688 = t549 * t557;
t434 = t474 * t553 - t548 * t688;
t640 = t553 * t689;
t480 = -t550 * t557 + t640;
t547 = qJ(5) + qJ(6);
t538 = sin(t547);
t539 = cos(t547);
t732 = t407 * t421 - g(1) * (-t434 * t538 - t473 * t539) - g(2) * (-t432 * t538 - t471 * t539) - g(3) * (-t480 * t538 + t539 * t687) + t383;
t417 = -qJD(5) * t492 + t556 * qJDD(3) + t552 * t736;
t622 = t557 * t645;
t643 = qJDD(2) * t553;
t583 = t622 + t643;
t491 = qJDD(5) + t583;
t644 = qJDD(1) * t550;
t721 = t502 * t655 + t553 * t739;
t577 = -t557 * t644 + t721;
t572 = qJDD(4) + t577;
t388 = t583 * pkin(4) + qJDD(3) * t559 + t572;
t387 = t556 * t388;
t705 = qJ(4) * t557;
t604 = pkin(9) * t553 - t705;
t571 = qJD(3) * t604 - t654;
t396 = qJD(2) * t571 + qJDD(2) * t487 + t589;
t568 = -t393 * qJD(5) - t552 * t396 + t387;
t373 = pkin(5) * t491 - pkin(10) * t417 + t568;
t418 = -qJD(5) * t627 + qJDD(3) * t552 + (qJD(3) * qJD(5) - t736) * t556;
t652 = qJD(5) * t556;
t638 = -t552 * t388 - t556 * t396 - t419 * t652;
t653 = qJD(5) * t552;
t581 = t442 * t653 + t638;
t374 = -pkin(10) * t418 - t581;
t615 = t555 * t373 - t551 * t374;
t731 = t407 * t598 - g(1) * (t434 * t539 - t473 * t538) - g(2) * (t432 * t539 - t471 * t538) - g(3) * (t480 * t539 + t538 * t687) + t615;
t484 = qJDD(6) + t491;
t730 = t484 * MDP(27) + (-t421 ^ 2 + t598 ^ 2) * MDP(24) - t421 * MDP(23) * t598;
t657 = qJD(3) * t553;
t535 = pkin(3) * t657;
t456 = t535 + t571;
t501 = t715 * t655;
t680 = t556 * t558;
t576 = t549 * (-t552 * t554 + t553 * t680);
t729 = -qJD(1) * t576 - t456 * t552 + t556 * t501;
t513 = t715 * t553;
t683 = t553 * t558;
t575 = t549 * (t552 * t683 + t554 * t556);
t728 = qJD(1) * t575 - t556 * t456 + t487 * t653 - t552 * t501 - t513 * t652;
t497 = t552 * t513;
t672 = t556 * t487 + t497;
t447 = -t544 - t454;
t446 = -qJD(3) * pkin(3) + t734;
t470 = t556 * t491;
t727 = -t530 * t653 + t470;
t725 = -MDP(10) + MDP(13);
t724 = MDP(14) - MDP(11);
t455 = -t634 - t735;
t722 = t455 * t664 + qJDD(4);
t720 = g(1) * t434 + g(2) * t432 + g(3) * t480;
t719 = t427 * t530 + t559 * t491;
t560 = qJD(3) ^ 2;
t579 = pkin(8) * t560 - t610;
t718 = 0.2e1 * qJDD(2) * pkin(2) + t549 * (-g(3) * t558 + t626) - t579 - t603;
t468 = t535 + t586;
t717 = qJD(2) * (-t468 + t635) + t733 - t579 + t741;
t614 = t417 * t551 + t555 * t418;
t379 = -qJD(6) * t598 + t614;
t716 = t557 * t723;
t709 = pkin(10) - t559;
t708 = qJD(2) * pkin(2);
t706 = pkin(8) * qJDD(3);
t702 = qJDD(3) * pkin(3);
t392 = t556 * t419 - t442 * t552;
t384 = -pkin(10) * t494 + t392;
t381 = pkin(5) * t530 + t384;
t701 = t381 * t555;
t700 = t385 * t555;
t699 = t417 * t556;
t698 = t491 * t552;
t697 = t492 * t530;
t696 = t494 * t530;
t694 = t530 * t556;
t693 = t538 * t553;
t692 = t539 * t553;
t686 = t551 * t373;
t685 = t551 * t552;
t684 = t552 * t553;
t682 = t555 * t556;
t681 = t556 * t557;
t678 = qJDD(1) - g(3);
t592 = pkin(5) * t557 - pkin(10) * t684;
t620 = pkin(10) * t557 - t487;
t677 = t592 * qJD(3) + (t556 * t620 - t497) * qJD(5) + t729;
t651 = qJD(5) * t557;
t628 = t552 * t651;
t676 = -(t553 * t656 + t628) * pkin(10) + t728;
t630 = t556 * t664;
t674 = -t551 * t653 - t552 * t649 + t555 * t630 - t664 * t685 + t682 * t723;
t536 = pkin(3) * t664;
t463 = qJD(2) * t604 + t536;
t673 = t552 * t441 + t556 * t463;
t636 = -pkin(5) * t556 - pkin(4);
t671 = pkin(5) * t652 - t636 * t664 + t734;
t514 = t715 * t557;
t545 = t553 ^ 2;
t546 = t557 ^ 2;
t670 = t545 - t546;
t669 = t545 + t546;
t665 = qJD(2) * t549;
t662 = qJD(3) * t454;
t661 = qJD(3) * t492;
t660 = qJD(3) * t494;
t650 = qJD(5) * t559;
t648 = qJD(6) * t555;
t561 = qJD(2) ^ 2;
t639 = t553 * t557 * t561;
t637 = t555 * t417 - t551 * t418 - t492 * t648;
t633 = t557 * t666;
t632 = t554 * t665;
t631 = t558 * t665;
t505 = t709 * t556;
t612 = qJD(6) * t381 + t374;
t611 = t502 * t657 - t553 * t644 - t557 * t739;
t609 = g(1) * t474 + g(2) * t472;
t597 = -t682 + t685;
t608 = -t597 * t484 + t521 * t675;
t429 = t556 * t441;
t504 = t709 * t552;
t606 = qJD(2) * t592 - qJD(6) * t504 - t463 * t552 - t709 * t653 + t429;
t605 = pkin(10) * t630 + t505 * t723 + t673;
t376 = t381 * t551 + t700;
t498 = t556 * t513;
t408 = pkin(5) * t553 + t552 * t620 + t498;
t415 = -pkin(10) * t681 + t672;
t601 = t408 * t551 + t415 * t555;
t438 = t480 * t556 + t552 * t687;
t439 = -t480 * t552 + t549 * t680;
t600 = t438 * t555 + t439 * t551;
t599 = t438 * t551 - t439 * t555;
t542 = qJDD(3) * qJ(4);
t543 = qJD(3) * qJD(4);
t397 = -t542 - t543 + t611;
t481 = t550 * t553 + t554 * t688;
t587 = -t530 * t652 - t698;
t584 = -t484 * t495 - t521 * t674;
t378 = -t494 * t649 + t637;
t433 = t472 * t557 - t553 * t619;
t435 = t548 * t549 * t553 + t474 * t557;
t578 = -g(1) * t435 - g(2) * t433 - g(3) * t481;
t389 = -pkin(4) * t736 - t397;
t573 = t389 + t578;
t503 = -t634 - t708;
t567 = -t706 + (t503 + t634 - t708) * qJD(3);
t566 = t706 + (-t455 - t634 + t735) * qJD(3);
t565 = qJD(3) * t453 + t578 - t611;
t564 = -t577 + t720;
t398 = t572 - t702;
t562 = -t397 * t557 + t398 * t553 + (t446 * t557 + t447 * t553) * qJD(3) - t609;
t532 = pkin(5) * t552 + qJ(4);
t500 = t715 * t657;
t499 = -qJ(4) * t663 + t536;
t479 = pkin(5) * t681 + t514;
t465 = t495 * t557;
t464 = t597 * t557;
t443 = -pkin(5) * t628 + (-pkin(8) + t636) * t657;
t437 = qJD(3) * t481 + t553 * t631;
t436 = -qJD(3) * t640 + (t631 + t659) * t557;
t401 = t495 * t716 - t597 * t657;
t400 = qJD(3) * t585 + t597 * t716;
t391 = qJD(5) * t438 + t437 * t552 + t556 * t632;
t390 = qJD(5) * t439 + t437 * t556 - t552 * t632;
t377 = pkin(5) * t418 + t389;
t375 = -t385 * t551 + t701;
t1 = [t678 * MDP(1) + (-t397 * t481 + t398 * t480 - t436 * t447 + t437 * t446 - g(3)) * MDP(15) + (t390 * t530 + t418 * t481 + t436 * t492 + t438 * t491) * MDP(21) + (-t391 * t530 + t417 * t481 + t436 * t494 + t439 * t491) * MDP(22) + ((-qJD(6) * t599 + t390 * t555 - t391 * t551) * t521 + t600 * t484 + t436 * t421 + t481 * t379) * MDP(28) + (-(qJD(6) * t600 + t390 * t551 + t391 * t555) * t521 - t599 * t484 - t436 * t598 + t481 * t378) * MDP(29) + (t480 * t553 + t481 * t557) * t641 + (t436 * t557 + t437 * t553 + (t480 * t557 - t481 * t553) * qJD(3)) * MDP(12) * qJD(2) + ((qJD(2) * t455 * MDP(15) - qJDD(2) * MDP(4) + (-t553 * t724 + t557 * t725 - MDP(3)) * t561) * t554 + (-t406 * MDP(15) + qJDD(2) * MDP(3) - t561 * MDP(4) + t724 * t583 + t725 * t736) * t558) * t549 + t724 * (qJD(3) * t436 + qJDD(3) * t481) + t725 * (qJD(3) * t437 + qJDD(3) * t480); (t567 * t553 + t557 * t718) * MDP(10) + (-t553 * t718 + t567 * t557) * MDP(11) + (t553 * t717 + t566 * t557) * MDP(14) + (t566 * t553 - t557 * t717) * MDP(13) + (qJDD(3) * t553 + t557 * t560) * MDP(7) + (qJDD(3) * t557 - t553 * t560) * MDP(8) + (t378 * t464 + t379 * t465 - t400 * t421 - t401 * t598) * MDP(24) + (-t378 * t465 - t400 * t598) * MDP(23) + (-t601 * t484 - (t555 * t612 - t383 + t686) * t553 - t376 * t655 - t443 * t598 + t479 * t378 - t377 * t465 + t407 * t400 - g(1) * (-t473 * t692 - t474 * t538) - g(2) * (-t471 * t692 - t472 * t538) + ((-qJD(6) * t408 + t676) * t555 + (qJD(6) * t415 - t677) * t551) * t521 + (t598 * t633 - g(3) * (-t538 * t554 + t539 * t683)) * t549) * MDP(29) + (t378 * t553 + t400 * t521 - t465 * t484 - t598 * t655) * MDP(25) + 0.2e1 * (t553 * t642 - t645 * t670) * MDP(6) + (t678 * t687 + t610) * MDP(3) + (-t678 * t689 + t609) * MDP(4) + ((t408 * t555 - t415 * t551) * t484 + t615 * t553 + t375 * t655 + t443 * t421 + t479 * t379 - t377 * t464 - t407 * t401 - g(1) * (-t473 * t693 + t474 * t539) - g(2) * (-t471 * t693 + t472 * t539) + (t551 * t676 + t555 * t677) * t521 + (-t376 * t553 - t521 * t601) * qJD(6) + (-t421 * t633 - g(3) * (t538 * t683 + t539 * t554)) * t549) * MDP(28) + (t491 * t553 + t530 * t655) * MDP(20) + (t484 * t553 + t521 * t655) * MDP(27) + (-t417 * t552 * t557 + (t552 * t657 - t556 * t651) * t494) * MDP(16) + ((t530 * t658 + t417) * t553 + (t587 + t660) * t557) * MDP(18) + ((t530 * t656 - t418) * t553 + (-t661 - t727) * t557) * MDP(19) + (-t672 * t491 - t500 * t494 + t514 * t417 + t609 * t552 + (t610 * t556 + (qJD(3) * t427 + qJD(5) * t442) * t552 + t638) * t553 - g(3) * t576 + t728 * t530 + (-qJD(3) * t393 - t389 * t552 - t427 * t652 - t494 * t634) * t557) * MDP(22) + ((-t492 * t552 + t494 * t556) * t657 + (-t699 + t418 * t552 + (t492 * t556 + t494 * t552) * qJD(5)) * t557) * MDP(17) + (qJDD(2) * t545 + 0.2e1 * t553 * t622) * MDP(5) + (-t379 * t553 + t401 * t521 - t421 * t655 + t464 * t484) * MDP(26) + (t669 * t703 + (-g(3) * t554 - t625 * t669) * t549 + t562) * MDP(12) + ((-t487 * t552 + t498) * t491 - t500 * t492 + t514 * t418 - t609 * t556 + (-t427 * t656 + t387 + (-t396 + t610) * t552) * t553 - g(3) * t575 + t729 * t530 + (-t393 * t553 - t530 * t672) * qJD(5) + (qJD(3) * t392 + t389 * t556 - t427 * t653 - t492 * t634) * t557) * MDP(21) + (t455 * t468 + t562 * pkin(8) + ((-pkin(8) * g(3) - qJD(1) * t455) * t554 + (-t446 * t553 + t447 * t557) * t666) * t549 + (t610 + t741) * t591) * MDP(15) + qJDD(2) * MDP(2); (-t397 * qJ(4) - t398 * pkin(3) - t455 * t499 - t446 * t454 - g(1) * (-pkin(3) * t434 + qJ(4) * t435) - g(2) * (-pkin(3) * t432 + qJ(4) * t433) - g(3) * (-pkin(3) * t480 + qJ(4) * t481) - t734 * t447) * MDP(15) - MDP(5) * t639 + ((t504 * t551 - t505 * t555) * t484 + t532 * t379 + t377 * t495 + (t551 * t605 - t555 * t606) * t521 + t671 * t421 + t674 * t407 + t578 * t538) * MDP(28) + ((t492 * t557 - t553 * t694) * qJD(2) + t587) * MDP(19) + ((-t418 - t696) * t556 + (-t417 + t697) * t552) * MDP(17) - t565 * MDP(11) + t608 * MDP(25) + t584 * MDP(26) + (-t503 * t664 + t564 + t662) * MDP(10) + (-0.2e1 * t702 - t662 + (-qJD(2) * t499 - t644) * t557 - t720 + t721 + t722) * MDP(13) + ((-t494 * t557 - t530 * t684) * qJD(2) + t727) * MDP(18) + (-t552 * t696 + t699) * MDP(16) + (-(-t504 * t555 - t505 * t551) * t484 + t532 * t378 - t377 * t597 + (t551 * t606 + t555 * t605) * t521 - t671 * t598 + t675 * t407 + t578 * t539) * MDP(29) + (-t378 * t597 - t598 * t675) * MDP(23) + (-t378 * t495 + t379 * t597 - t421 * t675 + t598 * t674) * MDP(24) + (qJ(4) * t417 + t673 * t530 + t647 * t494 - t719 * t552 + (-t530 * t650 + t573) * t556) * MDP(22) + (qJ(4) * t418 - t429 * t530 + t647 * t492 + t719 * t556 + ((t463 - t650) * t530 + t573) * t552) * MDP(21) + MDP(8) * t642 + MDP(7) * t643 + (-pkin(3) * t553 + t705) * t641 + qJDD(3) * MDP(9) + (0.2e1 * t542 + 0.2e1 * t543 + (t455 * t557 + t499 * t553) * qJD(2) + t565) * MDP(14) + t670 * t561 * MDP(6) + (-MDP(11) * t503 - t530 * MDP(20) - t392 * MDP(21) + t393 * MDP(22) + MDP(25) * t598 + t421 * MDP(26) - t521 * MDP(27) - t375 * MDP(28) + t376 * MDP(29)) * t663; t553 * t641 + (qJDD(3) + t639) * MDP(13) + (-t545 * t561 - t560) * MDP(14) + (qJD(3) * t447 - t564 - t702 + t722) * MDP(15) + (t470 - t661) * MDP(21) + (-t660 - t698) * MDP(22) + (-qJD(3) * t421 + t608) * MDP(28) + (qJD(3) * t598 + t584) * MDP(29) + (-MDP(21) * t530 * t552 - MDP(22) * t694) * t530; t494 * t492 * MDP(16) + (-t492 ^ 2 + t494 ^ 2) * MDP(17) + (t417 + t697) * MDP(18) + (-t418 + t696) * MDP(19) + t491 * MDP(20) + (t393 * t530 - t427 * t494 - g(1) * (t434 * t556 - t473 * t552) - g(2) * (t432 * t556 - t471 * t552) - g(3) * t438 + t568) * MDP(21) + (t392 * t530 + t427 * t492 - g(1) * (-t434 * t552 - t473 * t556) - g(2) * (-t432 * t552 - t471 * t556) - g(3) * t439 + t581) * MDP(22) + (t378 + t738) * MDP(25) + (-t379 - t737) * MDP(26) + (-(-t384 * t551 - t700) * t521 - t376 * qJD(6) + (-t421 * t494 + t484 * t555 - t521 * t649) * pkin(5) + t731) * MDP(28) + ((-t385 * t521 - t373) * t551 + (t384 * t521 - t612) * t555 + (-t484 * t551 + t494 * t598 - t521 * t648) * pkin(5) + t732) * MDP(29) + t730; (t637 + t738) * MDP(25) + (-t614 - t737) * MDP(26) + (t376 * t521 + t731) * MDP(28) + (-t555 * t374 + t375 * t521 - t686 + t732) * MDP(29) + (-MDP(25) * t695 + MDP(26) * t598 - MDP(28) * t376 - MDP(29) * t701) * qJD(6) + t730;];
tau  = t1;
