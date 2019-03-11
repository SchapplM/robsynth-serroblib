% Calculate vector of inverse dynamics joint torques for
% S6RPRRPR3
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
%   pkin=[a2,a3,a4,a5,a6,d1,d3,d4,d6,theta2]';
% MDP [29x1]
%   Minimal dynamic parameter vector (fixed base model)
%   see S6RPRRPR3_convert_par2_MPV_fixb.m
% 
% Output:
% tau [6x1]
%   joint torques of inverse dynamics (contains inertial, gravitational Coriolis and centrifugal forces)

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 05:08
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tau = S6RPRRPR3_invdynJ_fixb_mdp_slag_vp(qJ, qJD, qJDD, g, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(3,1),zeros(10,1),zeros(29,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRRPR3_invdynJ_fixb_mdp_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RPRRPR3_invdynJ_fixb_mdp_slag_vp: qJD has to be [6x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [6 1]), ...
  'S6RPRRPR3_invdynJ_fixb_mdp_slag_vp: qJDD has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RPRRPR3_invdynJ_fixb_mdp_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RPRRPR3_invdynJ_fixb_mdp_slag_vp: pkin has to be [10x1] (double)');
assert(isreal(MDP) && all(size(MDP) == [29 1]), ...
  'S6RPRRPR3_invdynJ_fixb_mdp_slag_vp: MDP has to be [29x1] (double)'); 

%% Symbolic Calculation
% From invdyn_joint_fixb_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 05:07:31
% EndTime: 2019-03-09 05:07:41
% DurationCPUTime: 8.35s
% Computational Cost: add. (4096->577), mult. (8234->730), div. (0->0), fcn. (5461->12), ass. (0->239)
t542 = sin(qJ(4));
t546 = cos(qJ(4));
t628 = t546 * qJD(3);
t543 = sin(qJ(3));
t643 = qJD(1) * t543;
t479 = t542 * t643 - t628;
t641 = qJD(3) * t542;
t481 = t546 * t643 + t641;
t539 = sin(pkin(10));
t520 = pkin(1) * t539 + pkin(7);
t502 = t520 * qJD(1);
t547 = cos(qJ(3));
t456 = t547 * qJD(2) - t543 * t502;
t598 = qJD(3) * pkin(3) + t456;
t562 = qJ(5) * t481 + t598;
t688 = pkin(4) + pkin(5);
t383 = -t479 * t688 + t562;
t536 = qJ(1) + pkin(10);
t526 = sin(t536);
t527 = cos(t536);
t661 = t542 * t547;
t451 = -t526 * t546 + t527 * t661;
t659 = t546 * t547;
t452 = t526 * t542 + t527 * t659;
t541 = sin(qJ(6));
t545 = cos(qJ(6));
t397 = t451 * t545 - t452 * t541;
t426 = t479 * t541 + t481 * t545;
t660 = t543 * t546;
t662 = t542 * t545;
t458 = t541 * t660 - t543 * t662;
t624 = qJD(1) * qJD(3);
t609 = t547 * t624;
t622 = qJDD(1) * t543;
t636 = qJD(4) * t543;
t706 = qJD(1) * t636 - qJDD(3);
t415 = -qJD(4) * t628 + (-t609 - t622) * t546 + t706 * t542;
t531 = t547 * qJDD(1);
t475 = t543 * t624 + qJDD(4) - t531;
t496 = t520 * qJDD(1);
t709 = qJD(3) * t456;
t407 = qJDD(3) * pkin(8) + qJDD(2) * t543 + t496 * t547 + t709;
t574 = pkin(3) * t547 + pkin(8) * t543 + pkin(2);
t540 = cos(pkin(10));
t686 = pkin(1) * t540;
t470 = -t574 - t686;
t597 = pkin(3) * t543 - pkin(8) * t547;
t490 = t597 * qJD(3);
t417 = qJD(1) * t490 + qJDD(1) * t470;
t457 = t543 * qJD(2) + t547 * t502;
t444 = qJD(3) * pkin(8) + t457;
t445 = t470 * qJD(1);
t635 = qJD(4) * t546;
t637 = qJD(4) * t542;
t600 = t542 * t407 - t546 * t417 + t444 * t635 + t445 * t637;
t575 = qJDD(5) + t600;
t362 = pkin(9) * t415 - t475 * t688 + t575;
t468 = t475 * qJ(5);
t642 = qJD(1) * t547;
t708 = qJD(4) - t642;
t506 = t708 * qJD(5);
t616 = t546 * t407 + t542 * t417 + t445 * t635;
t563 = -t444 * t637 + t616;
t368 = t468 + t506 + t563;
t416 = t542 * (qJD(3) * (qJD(4) + t642) + t622) + t706 * t546;
t364 = pkin(9) * t416 + t368;
t606 = t545 * t362 - t541 * t364;
t449 = t526 * t661 + t527 * t546;
t450 = t526 * t659 - t527 * t542;
t698 = t449 * t545 - t450 * t541;
t714 = g(1) * t397 + g(2) * t698 - g(3) * t458 + t383 * t426 - t606;
t469 = -qJDD(6) + t475;
t577 = -t545 * t479 + t481 * t541;
t713 = MDP(23) * t426 * t577 + (t426 ^ 2 - t577 ^ 2) * MDP(24) - t469 * MDP(27);
t625 = -qJD(6) + t708;
t712 = t426 * t625;
t711 = -qJD(2) * qJD(3) - t496;
t678 = g(3) * t547;
t594 = g(1) * t527 + g(2) * t526;
t696 = t543 * t594;
t710 = t678 - t696;
t707 = qJD(6) - qJD(4);
t702 = t625 * t577;
t701 = qJD(5) * t542 + t457;
t611 = t542 * t636;
t612 = t547 * t628;
t699 = -t611 + t612;
t395 = -t542 * t444 + t546 * t445;
t626 = qJD(5) - t395;
t394 = pkin(4) * t479 - t562;
t683 = pkin(8) * t475;
t695 = -t394 * t708 + t683;
t621 = MDP(17) + MDP(19);
t694 = t546 * MDP(18) + t542 * t621;
t675 = qJ(5) * t542;
t693 = -t546 * t688 - t675;
t398 = t451 * t541 + t452 * t545;
t483 = t541 * t542 + t545 * t546;
t459 = t483 * t543;
t578 = t449 * t541 + t450 * t545;
t627 = pkin(9) * t481 - t626;
t377 = -t688 * t708 - t627;
t631 = qJD(6) * t545;
t617 = t541 * t362 + t545 * t364 + t377 * t631;
t690 = -g(1) * t398 - g(2) * t578 - g(3) * t459 - t383 * t577 + t617;
t689 = t481 ^ 2;
t687 = pkin(8) - pkin(9);
t685 = pkin(4) * t475;
t684 = pkin(4) * t542;
t677 = pkin(8) * qJD(4);
t676 = qJ(5) * t479;
t674 = qJ(5) * t546;
t396 = t546 * t444 + t542 * t445;
t508 = t708 * qJ(5);
t389 = t508 + t396;
t673 = t389 * t708;
t672 = t396 * t708;
t670 = t479 * t481;
t669 = t479 * t708;
t668 = t481 * t708;
t667 = t708 * t546;
t666 = t520 * t542;
t385 = pkin(9) * t479 + t396;
t379 = t385 + t508;
t665 = t541 * t379;
t664 = t541 * t546;
t663 = t542 * t543;
t658 = qJDD(2) - g(3);
t610 = t543 * t635;
t639 = qJD(3) * t547;
t386 = qJD(6) * t459 + t541 * t699 - t545 * t610 - t639 * t662;
t657 = t386 * t625 + t458 * t469;
t656 = -t416 * t660 - t479 * t612;
t565 = t547 * t483;
t655 = -qJD(1) * t565 - t483 * t707;
t613 = t542 * t642;
t632 = qJD(6) * t541;
t654 = t541 * t635 + t542 * t631 - t546 * t632 - t642 * t664 + (t613 - t637) * t545;
t618 = t688 * t542;
t568 = -t618 + t674;
t653 = t568 * t708 + t701;
t487 = t597 * qJD(1);
t652 = t546 * t456 + t542 * t487;
t651 = t475 * t660 + t612 * t708;
t650 = t470 * t635 + t542 * t490;
t587 = -t674 + t684;
t649 = t587 * t708 - t701;
t486 = t520 * t659;
t648 = t542 * t470 + t486;
t647 = t594 * t660;
t537 = t543 ^ 2;
t646 = -t547 ^ 2 + t537;
t645 = MDP(20) * t542;
t521 = -pkin(2) - t686;
t503 = qJD(1) * t521;
t640 = qJD(3) * t543;
t638 = qJD(4) * t479;
t633 = qJD(5) * t546;
t620 = MDP(18) - MDP(21);
t619 = pkin(9) * t659;
t505 = t687 * t546;
t615 = -t545 * t415 + t541 * t416 + t479 * t631;
t403 = qJ(5) * t643 + t652;
t614 = qJ(5) * t640 + t650;
t607 = -pkin(4) - t666;
t605 = -t415 * t541 - t545 * t416;
t441 = t542 * t456;
t603 = -t487 * t546 + t441;
t485 = t520 * t661;
t602 = t470 * t546 - t485;
t599 = t547 * qJDD(2) - t502 * t639 + t543 * t711;
t596 = -g(1) * t449 + g(2) * t451;
t595 = g(1) * t450 - g(2) * t452;
t593 = g(1) * t526 - g(2) * t527;
t544 = sin(qJ(1));
t548 = cos(qJ(1));
t592 = g(1) * t544 - g(2) * t548;
t419 = -qJ(5) * t547 + t648;
t591 = -qJD(4) * t486 - t470 * t637 + t490 * t546;
t590 = (-t543 * t688 - t619) * qJD(1) + t603 + t707 * t505;
t504 = t687 * t542;
t589 = pkin(9) * t613 - qJD(6) * t504 + t687 * t637 + t403;
t588 = pkin(4) * t546 + t675;
t586 = qJ(5) * t545 - t541 * t688;
t585 = qJ(5) * t541 + t545 * t688;
t369 = t575 - t685;
t584 = t368 * t546 + t369 * t542;
t367 = t541 * t377 + t545 * t379;
t576 = -t662 + t664;
t387 = -t543 * t576 * t707 + qJD(3) * t565;
t583 = t387 * t625 + t459 * t469;
t388 = -pkin(4) * t708 + t626;
t582 = t388 * t546 - t389 * t542;
t581 = t388 * t542 + t389 * t546;
t535 = t547 * pkin(4);
t399 = pkin(5) * t547 + t485 + t535 + (-pkin(9) * t543 - t470) * t546;
t409 = pkin(9) * t663 + t419;
t580 = t399 * t545 - t409 * t541;
t579 = t399 * t541 + t409 * t545;
t572 = pkin(3) + t588;
t571 = t677 * t708 + t678;
t567 = t475 * t542 + t635 * t708;
t566 = qJDD(3) * pkin(3) + t599;
t555 = -qJ(5) * t415 + qJD(5) * t481 + t566;
t370 = pkin(4) * t416 - t555;
t564 = -t370 - t571;
t372 = -t481 * t632 + t615;
t561 = -qJD(1) * t503 + t594;
t560 = -t598 * t708 - t683;
t559 = 0.2e1 * qJD(3) * t503 - qJDD(3) * t520;
t558 = -g(3) * t543 - t547 * t594;
t557 = g(1) * t451 + g(2) * t449 + g(3) * t663 - t600;
t373 = qJD(6) * t426 + t605;
t550 = qJD(3) ^ 2;
t556 = -0.2e1 * qJDD(1) * t521 - t520 * t550 + t593;
t554 = t394 * t481 + qJDD(5) - t557;
t553 = g(1) * t452 + g(2) * t450 + g(3) * t660 + t395 * t708 - t563;
t516 = qJ(5) * t660;
t495 = qJDD(3) * t547 - t543 * t550;
t494 = qJDD(3) * t543 + t547 * t550;
t471 = pkin(3) - t693;
t461 = t481 * t640;
t436 = -t516 + (t520 + t684) * t543;
t430 = pkin(4) * t481 + t676;
t427 = t516 + (-t520 - t618) * t543;
t420 = t535 - t602;
t405 = -pkin(4) * t643 + t603;
t402 = -t481 * t688 - t676;
t392 = (qJD(4) * t588 - t633) * t543 + (t520 + t587) * t639;
t390 = -t415 + t669;
t382 = t607 * t640 - t591;
t381 = (qJD(4) * t693 + t633) * t543 + (-t520 + t568) * t639;
t378 = -qJD(5) * t547 + (-t543 * t628 - t547 * t637) * t520 + t614;
t375 = (pkin(9) * qJD(4) - qJD(3) * t520) * t660 + (-qJD(5) + (pkin(9) * qJD(3) - qJD(4) * t520) * t542) * t547 + t614;
t374 = pkin(9) * t611 + (-t619 + (-pkin(5) + t607) * t543) * qJD(3) - t591;
t371 = t372 * t547;
t366 = t377 * t545 - t665;
t365 = -t416 * t688 + t555;
t1 = [(g(1) * t548 + g(2) * t544) * MDP(3) + t494 * MDP(7) + t495 * MDP(8) + (t372 * t459 + t387 * t426) * MDP(23) + (t591 * t708 + t602 * t475 + ((t479 * t520 - t542 * t598) * qJD(3) + t600) * t547 + (-t598 * t635 - t566 * t542 + t520 * t416 + (t666 * t708 + t395) * qJD(3)) * t543 + t595) * MDP(17) + (t415 * t547 - t611 * t708 + t461 + t651) * MDP(14) + (-t475 * t547 + t640 * t708) * MDP(16) + ((-t641 * t708 + t416) * t547 + (-qJD(3) * t479 - t567) * t543) * MDP(15) + (-t382 * t708 + t392 * t479 + t416 * t436 - t420 * t475 + (t394 * t641 + t369) * t547 + (-qJD(3) * t388 + t370 * t542 + t394 * t635) * t543 + t595) * MDP(19) + (t378 * t708 - t392 * t481 + t415 * t436 + t419 * t475 + (-t394 * t628 - t368) * t547 + (qJD(3) * t389 - t370 * t546 + t394 * t637) * t543 - t596) * MDP(21) + (-t650 * t708 - t648 * t475 + ((t520 * t708 - t444) * t637 + (t481 * t520 - t546 * t598) * qJD(3) + t616) * t547 + (t598 * t637 - t566 * t546 - t520 * t415 + (t520 * t667 - t396) * qJD(3)) * t543 + t596) * MDP(18) + (-(t374 * t545 - t375 * t541) * t625 - t580 * t469 + t606 * t547 - t366 * t640 + t381 * t577 + t427 * t373 + t365 * t458 + t383 * t386 + g(1) * t578 - g(2) * t398 + (-t367 * t547 + t579 * t625) * qJD(6)) * MDP(28) + (-t469 * t547 + t625 * t640) * MDP(27) + (-t372 * t458 - t373 * t459 - t386 * t426 - t387 * t577) * MDP(24) + (-t373 * t547 + t577 * t640 + t657) * MDP(26) + (t543 * t559 + t547 * t556) * MDP(10) + (-t543 * t556 + t547 * t559) * MDP(11) + qJDD(1) * MDP(1) + ((qJD(6) * t580 + t374 * t541 + t375 * t545) * t625 + t579 * t469 - (-t379 * t632 + t617) * t547 + t367 * t640 + t381 * t426 + t427 * t372 + t365 * t459 + t383 * t387 + g(1) * t698 - g(2) * t397) * MDP(29) + (-t415 * t660 + t481 * t699) * MDP(12) + (-t481 * t610 + (-t481 * t639 + (t415 + t638) * t543) * t542 + t656) * MDP(13) + 0.2e1 * (t531 * t543 - t624 * t646) * MDP(6) + (-t378 * t479 + t382 * t481 - t415 * t420 - t416 * t419 + t582 * t639 + (-qJD(4) * t581 - t368 * t542 + t369 * t546 + t593) * t543) * MDP(20) + (-t426 * t640 + t371 - t583) * MDP(25) + (t592 + (t539 ^ 2 + t540 ^ 2) * qJDD(1) * pkin(1)) * pkin(1) * MDP(4) + (qJDD(1) * t537 + 0.2e1 * t543 * t609) * MDP(5) + (t368 * t419 + t389 * t378 + t370 * t436 + t394 * t392 + t369 * t420 + t388 * t382 - g(1) * (-pkin(1) * t544 - pkin(4) * t450 - qJ(5) * t449) - g(2) * (pkin(1) * t548 + pkin(4) * t452 + qJ(5) * t451) + (-g(1) * pkin(7) - g(2) * t574) * t527 + (-g(2) * pkin(7) + g(1) * t574) * t526) * MDP(22) + t592 * MDP(2); t658 * MDP(4) + t495 * MDP(10) - t494 * MDP(11) + t461 * MDP(18) + t656 * MDP(20) + t651 * MDP(21) - g(3) * MDP(22) + t657 * MDP(28) + (t371 + t583) * MDP(29) + t621 * t479 * t640 + (-t370 * MDP(22) + t373 * MDP(28) - t621 * t416 + t620 * t415 + (t581 * MDP(22) + t481 * t645 - t694 * t708) * qJD(3)) * t547 + (-t415 * t645 + t584 * MDP(22) - t694 * t475 + ((t479 * t542 + t481 * t546) * MDP(20) + t582 * MDP(22) - (-t542 * t620 + t546 * t621) * t708) * qJD(4) + (-MDP(21) * t481 + MDP(22) * t394 - MDP(28) * t577 - MDP(29) * t426) * qJD(3)) * t543; MDP(8) * t531 + MDP(7) * t622 + (-pkin(3) * t416 + t441 * t708 - t457 * t479 + (-t678 + t566 - (t487 + t677) * t708) * t546 + t560 * t542 + t647) * MDP(17) + qJDD(3) * MDP(9) + (t469 * t483 + t625 * t654) * MDP(26) + (t469 * t576 - t625 * t655) * MDP(25) + (-t372 * t483 + t373 * t576 - t426 * t654 - t577 * t655) * MDP(24) + (-t372 * t576 + t426 * t655) * MDP(23) + (pkin(3) * t415 + t652 * t708 - t457 * t481 + t560 * t546 + (-t566 + t571 - t696) * t542) * MDP(18) + ((t504 * t541 + t505 * t545) * t469 + t471 * t372 - (t541 * t590 + t545 * t589) * t625 + t653 * t426 + t655 * t383 + (-t365 + t710) * t576) * MDP(29) + (-t403 * t708 - t415 * t572 - t649 * t481 + t695 * t546 + (t564 + t696) * t542) * MDP(21) + (-(t504 * t545 - t505 * t541) * t469 + t471 * t373 - g(3) * t565 - (t541 * t589 - t545 * t590) * t625 + t653 * t577 + t654 * t383 + (t365 + t696) * t483) * MDP(28) + (t405 * t708 - t416 * t572 + t649 * t479 - t542 * t695 + t564 * t546 + t647) * MDP(19) + (-t388 * t405 - t389 * t403 + t649 * t394 + (qJD(4) * t582 + t558 + t584) * pkin(8) + (-t370 - t710) * t572) * MDP(22) + (qJD(3) * t457 + t543 * t561 + t599 - t678) * MDP(10) + (t403 * t479 - t405 * t481 + (t368 + t708 * t388 + (qJD(4) * t481 - t416) * pkin(8)) * t546 + (t369 - t673 + (-t415 + t638) * pkin(8)) * t542 + t558) * MDP(20) + ((-t415 - t669) * t546 + (-t416 - t668) * t542) * MDP(13) + (-t415 * t542 + t481 * t667) * MDP(12) + (-t708 * t637 + t475 * t546 + (t479 * t543 + t661 * t708) * qJD(1)) * MDP(15) + (t709 + (qJD(3) * t502 - t658) * t543 + (t561 + t711) * t547) * MDP(11) + ((-t481 * t543 - t659 * t708) * qJD(1) + t567) * MDP(14) + (-MDP(16) * t708 - MDP(17) * t395 + MDP(18) * t396 + t388 * MDP(19) - t389 * MDP(21) + t426 * MDP(25) - MDP(26) * t577 - MDP(27) * t625 + t366 * MDP(28) - t367 * MDP(29)) * t643 + (-MDP(5) * t543 * t547 + MDP(6) * t646) * qJD(1) ^ 2; MDP(12) * t670 + (-t479 ^ 2 + t689) * MDP(13) + t390 * MDP(14) + (-t416 + t668) * MDP(15) + t475 * MDP(16) + (t481 * t598 + t557 + t672) * MDP(17) + (-t479 * t598 + t553) * MDP(18) + (-t430 * t479 - t554 + t672 + 0.2e1 * t685) * MDP(19) + (pkin(4) * t415 - qJ(5) * t416 + (t389 - t396) * t481 + (t388 - t626) * t479) * MDP(20) + (-t394 * t479 + t430 * t481 + 0.2e1 * t468 + 0.2e1 * t506 - t553) * MDP(21) + (t368 * qJ(5) - t369 * pkin(4) - t394 * t430 - t388 * t396 - g(1) * (-pkin(4) * t451 + qJ(5) * t452) - g(2) * (-pkin(4) * t449 + qJ(5) * t450) - g(3) * (-pkin(4) * t663 + t516) + t626 * t389) * MDP(22) + (-t372 + t702) * MDP(25) + (t373 + t712) * MDP(26) + (t585 * t469 - t402 * t577 - (-t545 * t385 + t541 * t627) * t625 + (t586 * t625 + t367) * qJD(6) + t714) * MDP(28) + (t586 * t469 - t402 * t426 - (t541 * t385 + t545 * t627) * t625 + (-t585 * t625 - t665) * qJD(6) + t690) * MDP(29) - t713; (-t475 + t670) * MDP(19) + t390 * MDP(20) + (-t708 ^ 2 - t689) * MDP(21) + (t554 - t673 - t685) * MDP(22) + (-t545 * t469 - t481 * t577) * MDP(28) + (-t426 * t481 + t541 * t469) * MDP(29) - (MDP(28) * t541 + MDP(29) * t545) * t625 ^ 2; (t615 - t702) * MDP(25) + (-t605 - t712) * MDP(26) + (-t367 * t625 - t714) * MDP(28) + (-t366 * t625 - t690) * MDP(29) + ((-MDP(26) * t481 - MDP(28) * t379) * t545 + (-MDP(25) * t481 - MDP(26) * t479 - MDP(28) * t377 + MDP(29) * t379) * t541) * qJD(6) + t713;];
tau  = t1;
