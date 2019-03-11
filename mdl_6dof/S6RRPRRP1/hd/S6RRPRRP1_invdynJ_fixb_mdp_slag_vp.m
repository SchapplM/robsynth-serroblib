% Calculate vector of inverse dynamics joint torques for
% S6RRPRRP1
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
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d4,d5,theta3]';
% MDP [28x1]
%   Minimal dynamic parameter vector (fixed base model)
%   see S6RRPRRP1_convert_par2_MPV_fixb.m
% 
% Output:
% tau [6x1]
%   joint torques of inverse dynamics (contains inertial, gravitational Coriolis and centrifugal forces)

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 11:42
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tau = S6RRPRRP1_invdynJ_fixb_mdp_slag_vp(qJ, qJD, qJDD, g, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(3,1),zeros(10,1),zeros(28,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPRRP1_invdynJ_fixb_mdp_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRPRRP1_invdynJ_fixb_mdp_slag_vp: qJD has to be [6x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [6 1]), ...
  'S6RRPRRP1_invdynJ_fixb_mdp_slag_vp: qJDD has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RRPRRP1_invdynJ_fixb_mdp_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RRPRRP1_invdynJ_fixb_mdp_slag_vp: pkin has to be [10x1] (double)');
assert(isreal(MDP) && all(size(MDP) == [28 1]), ...
  'S6RRPRRP1_invdynJ_fixb_mdp_slag_vp: MDP has to be [28x1] (double)'); 

%% Symbolic Calculation
% From invdyn_joint_fixb_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 11:41:31
% EndTime: 2019-03-09 11:41:42
% DurationCPUTime: 8.40s
% Computational Cost: add. (9295->509), mult. (22090->630), div. (0->0), fcn. (16951->14), ass. (0->238)
t604 = cos(qJ(5));
t671 = qJD(5) * t604;
t596 = sin(pkin(10));
t597 = cos(pkin(10));
t602 = sin(qJ(2));
t606 = cos(qJ(2));
t546 = -t596 * t602 + t597 * t606;
t536 = t546 * qJD(1);
t547 = t596 * t606 + t597 * t602;
t538 = t547 * qJD(1);
t601 = sin(qJ(4));
t605 = cos(qJ(4));
t490 = t605 * t536 - t538 * t601;
t749 = t490 * t604;
t758 = t671 - t749;
t600 = sin(qJ(5));
t537 = t547 * qJD(2);
t499 = -qJD(1) * t537 + qJDD(1) * t546;
t666 = qJD(1) * qJD(2);
t653 = t606 * t666;
t654 = t602 * t666;
t500 = qJDD(1) * t547 - t596 * t654 + t597 * t653;
t673 = qJD(4) * t605;
t674 = qJD(4) * t601;
t617 = t601 * t499 + t605 * t500 + t536 * t673 - t538 * t674;
t628 = t536 * t601 + t605 * t538;
t664 = qJD(2) + qJD(4);
t637 = t604 * t664;
t663 = qJDD(2) + qJDD(4);
t672 = qJD(5) * t600;
t411 = -qJD(5) * t637 - t600 * t663 - t604 * t617 + t628 * t672;
t409 = t411 * t600;
t410 = t604 * t411;
t646 = -t605 * t499 + t601 * t500;
t730 = t628 * qJD(4);
t427 = t646 + t730;
t426 = qJDD(5) + t427;
t419 = t600 * t426;
t473 = t600 * t628 - t637;
t475 = t600 * t664 + t604 * t628;
t669 = t628 * qJD(2);
t611 = t600 * t617 - t604 * t663;
t412 = qJD(5) * t475 + t611;
t688 = -t600 * t412 - t473 * t671;
t703 = t475 * t628;
t737 = t490 * t664;
t745 = qJD(5) - t490;
t753 = t745 * t600;
t738 = t475 * t753;
t757 = t663 * MDP(17) + (t669 - t646) * MDP(16) - t490 ^ 2 * MDP(14) + (-MDP(13) * t490 + t628 * MDP(14) - t745 * MDP(24)) * t628 + (t617 - t737) * MDP(15) + (t475 * t758 - t409) * MDP(20) + (t745 * t758 + t419 - t703) * MDP(22) + (t473 * t749 - t410 + t688 - t738) * MDP(21);
t593 = qJ(2) + pkin(10);
t587 = qJ(4) + t593;
t579 = sin(t587);
t603 = sin(qJ(1));
t607 = cos(qJ(1));
t634 = g(1) * t607 + g(2) * t603;
t756 = t634 * t579;
t599 = -qJ(3) - pkin(7);
t649 = qJD(2) * t599;
t532 = -qJD(3) * t602 + t606 * t649;
t565 = t599 * t602;
t496 = qJDD(2) * pkin(2) + qJD(1) * t532 + qJDD(1) * t565;
t531 = qJD(3) * t606 + t602 * t649;
t567 = t599 * t606;
t503 = qJD(1) * t531 - qJDD(1) * t567;
t453 = t597 * t496 - t503 * t596;
t429 = qJDD(2) * pkin(3) - pkin(8) * t500 + t453;
t454 = t596 * t496 + t597 * t503;
t434 = pkin(8) * t499 + t454;
t556 = qJD(1) * t567;
t541 = t596 * t556;
t555 = qJD(1) * t565;
t707 = qJD(2) * pkin(2);
t545 = t555 + t707;
t497 = t597 * t545 + t541;
t717 = pkin(8) * t538;
t467 = qJD(2) * pkin(3) + t497 - t717;
t695 = t597 * t556;
t498 = t596 * t545 - t695;
t718 = pkin(8) * t536;
t471 = t498 + t718;
t723 = -t605 * (qJD(4) * t467 + t434) - t601 * t429 + t471 * t674;
t389 = pkin(9) * t663 - t723;
t590 = t606 * pkin(2);
t584 = t590 + pkin(1);
t615 = pkin(2) * t654 - qJDD(1) * t584 + qJDD(3);
t470 = -t499 * pkin(3) + t615;
t395 = t427 * pkin(4) - pkin(9) * t617 + t470;
t431 = t601 * t467 + t605 * t471;
t422 = pkin(9) * t664 + t431;
t558 = -qJD(1) * t584 + qJD(3);
t506 = -pkin(3) * t536 + t558;
t436 = -pkin(4) * t490 - pkin(9) * t628 + t506;
t616 = t604 * t389 + t600 * t395 - t422 * t672 + t436 * t671;
t383 = -qJ(6) * t412 - qJD(6) * t473 + t616;
t571 = g(3) * t579;
t580 = cos(t587);
t657 = t580 * t634 + t571;
t394 = t604 * t395;
t403 = t422 * t604 + t436 * t600;
t381 = pkin(5) * t426 + qJ(6) * t411 - qJD(5) * t403 - qJD(6) * t475 - t389 * t600 + t394;
t398 = -qJ(6) * t473 + t403;
t728 = t398 * t745 + t381;
t755 = t383 * t604 - t600 * t728 - t657;
t504 = -t555 * t596 + t695;
t477 = t504 - t718;
t505 = t597 * t555 + t541;
t478 = t505 - t717;
t581 = pkin(2) * t597 + pkin(3);
t721 = pkin(2) * t596;
t644 = t581 * t605 - t601 * t721;
t735 = -t644 * qJD(4) + t477 * t601 + t478 * t605;
t430 = t605 * t467 - t601 * t471;
t421 = -pkin(4) * t664 - t430;
t750 = t421 * t490;
t719 = pkin(5) * t604;
t583 = pkin(4) + t719;
t598 = -qJ(6) - pkin(9);
t645 = -t579 * t598 + t580 * t583;
t747 = t645 + pkin(3) * cos(t593) + t590;
t452 = pkin(4) * t628 - pkin(9) * t490;
t744 = -t506 * t490 + t657 + t723;
t589 = t604 * qJ(6);
t743 = -pkin(5) * t628 + t490 * t589;
t704 = t473 * t628;
t734 = g(1) * t603 - g(2) * t607;
t736 = t734 * t579;
t677 = t601 * t581 + t605 * t721;
t680 = t677 * qJD(4) + t605 * t477 - t601 * t478;
t708 = t602 * pkin(2);
t511 = pkin(3) * t538 + qJD(1) * t708;
t441 = t452 + t511;
t732 = t600 * t441 + t604 * t735;
t691 = t604 * t607;
t694 = t600 * t603;
t525 = t580 * t694 + t691;
t692 = t603 * t604;
t693 = t600 * t607;
t527 = -t580 * t693 + t692;
t731 = -g(1) * t527 + g(2) * t525;
t711 = g(3) * t580;
t729 = t711 - t756;
t402 = -t422 * t600 + t604 * t436;
t727 = -t402 * t628 + t421 * t672 + t604 * t756;
t635 = -t605 * t429 + t601 * t434 + t467 * t674 + t471 * t673;
t390 = -pkin(4) * t663 + t635;
t418 = t421 * t671;
t710 = g(3) * t600;
t726 = t390 * t600 + t403 * t628 + t580 * t710 + t418;
t725 = -t506 * t628 - t635 - t729;
t722 = t475 ^ 2;
t720 = pkin(5) * t600;
t709 = g(3) * t606;
t502 = t546 * t601 + t547 * t605;
t706 = t390 * t502;
t397 = -qJ(6) * t475 + t402;
t391 = pkin(5) * t745 + t397;
t705 = t391 * t604;
t701 = t490 * t600;
t699 = t502 * t600;
t420 = t604 * t426;
t507 = t597 * t565 + t567 * t596;
t482 = -pkin(8) * t547 + t507;
t508 = t596 * t565 - t597 * t567;
t483 = pkin(8) * t546 + t508;
t448 = t482 * t601 + t483 * t605;
t445 = t604 * t448;
t530 = pkin(9) + t677;
t690 = -qJ(6) - t530;
t689 = -t397 + t391;
t686 = t604 * t430 + t600 * t452;
t514 = -pkin(3) * t546 - t584;
t627 = t605 * t546 - t547 * t601;
t446 = -pkin(4) * t627 - pkin(9) * t502 + t514;
t684 = t600 * t446 + t445;
t588 = t604 * qJD(6);
t643 = qJD(5) * t690;
t683 = qJ(6) * t701 + t600 * t643 + t588 - t732;
t438 = t604 * t441;
t682 = t604 * t643 - t438 + (-qJD(6) + t735) * t600 + t743;
t481 = t597 * t531 + t596 * t532;
t648 = qJD(5) * t598;
t679 = t588 - t686 + (qJ(6) * t490 + t648) * t600;
t450 = t604 * t452;
t678 = t604 * t648 - t450 + (-qJD(6) + t430) * t600 + t743;
t594 = t602 ^ 2;
t675 = -t606 ^ 2 + t594;
t665 = qJDD(1) * t606;
t586 = t602 * t707;
t662 = qJD(5) * pkin(9) * t745;
t480 = -t531 * t596 + t597 * t532;
t540 = t546 * qJD(2);
t458 = -pkin(8) * t540 + t480;
t459 = -pkin(8) * t537 + t481;
t629 = t482 * t605 - t483 * t601;
t406 = qJD(4) * t629 + t458 * t601 + t459 * t605;
t455 = qJD(4) * t627 - t537 * t601 + t540 * t605;
t456 = qJD(4) * t502 + t605 * t537 + t540 * t601;
t512 = pkin(3) * t537 + t586;
t415 = pkin(4) * t456 - pkin(9) * t455 + t512;
t660 = t604 * t406 + t600 * t415 + t446 * t671;
t655 = t502 * t671;
t652 = pkin(8) - t599 + t720;
t650 = -t390 - t711;
t641 = t745 * t604;
t638 = -qJD(5) * t436 - t389;
t529 = -pkin(4) - t644;
t632 = -t422 * t671 + t394;
t631 = -pkin(9) * t426 - t750;
t630 = -t426 * t530 - t750;
t626 = t579 * t583 + t580 * t598;
t625 = -qJ(6) * t455 - qJD(6) * t502;
t624 = t420 + (-t672 + t701) * t745;
t622 = -0.2e1 * pkin(1) * t666 - pkin(7) * qJDD(2);
t621 = -pkin(1) - t747;
t619 = t455 * t600 + t655;
t618 = t455 * t604 - t502 * t672;
t608 = qJD(2) ^ 2;
t613 = 0.2e1 * qJDD(1) * pkin(1) - pkin(7) * t608 + t734;
t609 = qJD(1) ^ 2;
t612 = pkin(1) * t609 - pkin(7) * qJDD(1) + t634;
t386 = t412 * pkin(5) + qJDD(6) + t390;
t407 = qJD(4) * t448 - t458 * t605 + t459 * t601;
t566 = pkin(9) * t604 + t589;
t564 = t598 * t600;
t528 = t580 * t691 + t694;
t526 = -t580 * t692 + t693;
t510 = t530 * t604 + t589;
t509 = t690 * t600;
t472 = t473 ^ 2;
t444 = t604 * t446;
t416 = t473 * pkin(5) + qJD(6) + t421;
t414 = t604 * t415;
t404 = -qJ(6) * t699 + t684;
t400 = -pkin(5) * t627 - t448 * t600 - t502 * t589 + t444;
t385 = -qJ(6) * t655 + (-qJD(5) * t448 + t625) * t600 + t660;
t384 = pkin(5) * t456 - t406 * t600 + t414 + t625 * t604 + (-t445 + (qJ(6) * t502 - t446) * t600) * qJD(5);
t1 = [(-t406 * t664 - t448 * t663 + t506 * t455 + t470 * t502 + t512 * t628 + t514 * t617 - t736) * MDP(19) + (-t384 * t475 - t385 * t473 + t400 * t411 - t404 * t412 + t736 + (-t398 * t600 - t705) * t455 + (-t381 * t604 - t383 * t600 + (t391 * t600 - t398 * t604) * qJD(5)) * t502) * MDP(27) + (-t453 * t547 + t454 * t546 - t480 * t538 + t481 * t536 - t497 * t540 - t498 * t537 + t499 * t508 - t500 * t507 - t634) * MDP(11) + (t455 * t664 + t502 * t663) * MDP(15) + (t454 * t508 + t498 * t481 + t453 * t507 + t497 * t480 - t615 * t584 + t558 * t586 - g(1) * (-t584 * t603 - t599 * t607) - g(2) * (t584 * t607 - t599 * t603)) * MDP(12) + (t602 * t622 + t606 * t613) * MDP(9) + (-t602 * t613 + t606 * t622) * MDP(10) + t634 * MDP(3) + 0.2e1 * (t602 * t665 - t666 * t675) * MDP(5) + (-t456 * t664 + t627 * t663) * MDP(16) + (-t502 * t427 + t455 * t490 - t456 * t628 + t617 * t627) * MDP(14) + (-t407 * t664 + t514 * t427 + t506 * t456 - t470 * t627 - t490 * t512 + t580 * t734 + t629 * t663) * MDP(18) + t734 * MDP(2) + (t455 * t628 + t502 * t617) * MDP(13) + (qJDD(1) * t594 + 0.2e1 * t602 * t653) * MDP(4) + (-t426 * t627 + t456 * t745) * MDP(24) + (t412 * t627 - t419 * t502 - t456 * t473 - t619 * t745) * MDP(23) + (t411 * t627 + t420 * t502 + t456 * t475 + t618 * t745) * MDP(22) + (-(-t448 * t672 + t660) * t745 - t684 * t426 + t616 * t627 - t403 * t456 + t407 * t475 + t629 * t411 + t604 * t706 - g(1) * t525 - g(2) * t527 + t618 * t421) * MDP(26) + ((-t448 * t671 + t414) * t745 + t444 * t426 - t632 * t627 + t402 * t456 + t407 * t473 - t629 * t412 + t502 * t418 - g(1) * t526 - g(2) * t528 + ((-qJD(5) * t446 - t406) * t745 - t448 * t426 - t638 * t627 + t706 + t421 * t455) * t600) * MDP(25) + (-t410 * t502 + t475 * t618) * MDP(20) + (qJDD(2) * t602 + t606 * t608) * MDP(6) + (qJDD(2) * t606 - t602 * t608) * MDP(7) + (t383 * t404 + t398 * t385 + t381 * t400 + t391 * t384 + t386 * (pkin(5) * t699 - t629) + t416 * (pkin(5) * t619 + t407) + (-g(1) * t652 + g(2) * t621) * t607 + (-g(1) * t621 - g(2) * t652) * t603) * MDP(28) + qJDD(1) * MDP(1) + ((-t473 * t604 - t475 * t600) * t455 + (t409 - t412 * t604 + (t473 * t600 - t475 * t604) * qJD(5)) * t502) * MDP(21); (g(3) * t602 + t606 * t612) * MDP(10) + ((t498 + t504) * t538 + (t497 - t505) * t536 + (t499 * t596 - t500 * t597) * pkin(2)) * MDP(11) + (-t391 * t641 + t411 * t509 - t412 * t510 - t683 * t473 - t682 * t475 + t755) * MDP(27) + (-t511 * t628 - t677 * t663 + t664 * t735 + t744) * MDP(19) + (t490 * t511 + t644 * t663 - t664 * t680 + t725) * MDP(18) + (t602 * t612 - t709) * MDP(9) + (-t497 * t504 - t498 * t505 + (-t709 + t453 * t597 + t454 * t596 + (-qJD(1) * t558 + t634) * t602) * pkin(2)) * MDP(12) + (t529 * t412 + t650 * t604 + t630 * t600 + t680 * t473 + (-t530 * t671 + t600 * t735 - t438) * t745 + t727) * MDP(25) + (-t529 * t411 + t630 * t604 - t600 * t756 + t680 * t475 + (t530 * t672 + t732) * t745 + t726) * MDP(26) + (t383 * t510 + t381 * t509 + t386 * (t529 - t719) - g(3) * t747 + (pkin(5) * t753 + t680) * t416 + t683 * t398 + t682 * t391 + t634 * (pkin(3) * sin(t593) + t708 + t626)) * MDP(28) + (-MDP(4) * t602 * t606 + MDP(5) * t675) * t609 + t602 * qJDD(1) * MDP(6) + qJDD(2) * MDP(8) + (t624 + t704) * MDP(23) + MDP(7) * t665 + t757; (-t536 ^ 2 - t538 ^ 2) * MDP(11) + (t497 * t538 - t498 * t536 + t615 - t734) * MDP(12) + (t646 + t669 + 0.2e1 * t730) * MDP(18) + (t617 + t737) * MDP(19) + (t624 - t704) * MDP(25) + (-t641 * t745 - t419 - t703) * MDP(26) + ((t473 * t490 + t411) * t604 + t738 + t688) * MDP(27) + (-t416 * t628 + t728 * t604 + (-t391 * t745 + t383) * t600 - t734) * MDP(28); (t431 * t664 + t725) * MDP(18) + (t430 * t664 + t744) * MDP(19) + (-t745 * t753 + t420 + t704) * MDP(23) + (-pkin(4) * t412 - t431 * t473 - t450 * t745 + (t430 * t745 + t631) * t600 + (t650 - t662) * t604 + t727) * MDP(25) + (pkin(4) * t411 + t686 * t745 - t431 * t475 + t631 * t604 + (-t756 + t662) * t600 + t726) * MDP(26) + (t411 * t564 - t412 * t566 - t679 * t473 - t678 * t475 - t745 * t705 + t755) * MDP(27) + (t383 * t566 + t381 * t564 - t386 * t583 - g(3) * t645 + (t720 * t745 - t431) * t416 + t679 * t398 + t678 * t391 + t634 * t626) * MDP(28) + t757; t475 * t473 * MDP(20) + (-t472 + t722) * MDP(21) + (t473 * t745 - t411) * MDP(22) + (-t611 + (-qJD(5) + t745) * t475) * MDP(23) + t426 * MDP(24) + (t403 * t745 - t421 * t475 + (t638 + t571) * t600 + t632 + t731) * MDP(25) + (g(1) * t528 - g(2) * t526 + t402 * t745 + t421 * t473 + t571 * t604 - t616) * MDP(26) + (pkin(5) * t411 - t473 * t689) * MDP(27) + (t689 * t398 + (-t416 * t475 + t579 * t710 + t381 + t731) * pkin(5)) * MDP(28); (-t472 - t722) * MDP(27) + (t391 * t475 + t398 * t473 + t386 + t729) * MDP(28);];
tau  = t1;
