% Calculate vector of inverse dynamics joint torques for
% S5RRRRP11
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% qJDD [5x1]
%   Generalized joint accelerations
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,alpha2,d1,d2,d3,d4]';
% MDP [28x1]
%   Minimal dynamic parameter vector (fixed base model)
%   see S5RRRRP11_convert_par2_MPV_fixb.m
% 
% Output:
% tau [5x1]
%   joint torques of inverse dynamics (contains inertial, gravitational Coriolis and centrifugal forces)

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 22:20
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tau = S5RRRRP11_invdynJ_fixb_mdp_slag_vp(qJ, qJD, qJDD, g, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(5,1),zeros(3,1),zeros(9,1),zeros(28,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRRRP11_invdynJ_fixb_mdp_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RRRRP11_invdynJ_fixb_mdp_slag_vp: qJD has to be [5x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [5 1]), ...
  'S5RRRRP11_invdynJ_fixb_mdp_slag_vp: qJDD has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RRRRP11_invdynJ_fixb_mdp_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5RRRRP11_invdynJ_fixb_mdp_slag_vp: pkin has to be [9x1] (double)');
assert(isreal(MDP) && all(size(MDP) == [28 1]), ...
  'S5RRRRP11_invdynJ_fixb_mdp_slag_vp: MDP has to be [28x1] (double)'); 

%% Symbolic Calculation
% From invdyn_joint_fixb_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 22:18:55
% EndTime: 2019-12-31 22:19:16
% DurationCPUTime: 13.48s
% Computational Cost: add. (7455->646), mult. (18467->849), div. (0->0), fcn. (14404->10), ass. (0->245)
t571 = cos(qJ(2));
t565 = sin(pkin(5));
t683 = qJD(1) * t565;
t555 = t571 * t683;
t617 = t555 - qJD(3);
t567 = sin(qJ(3));
t568 = sin(qJ(2));
t570 = cos(qJ(3));
t714 = cos(pkin(5));
t648 = t714 * qJD(1);
t611 = t648 + qJD(2);
t590 = qJD(3) * t611;
t642 = t714 * qJDD(1);
t606 = t642 + qJDD(2);
t681 = qJD(2) * t571;
t656 = t567 * t681;
t670 = qJDD(1) * t568;
t677 = qJD(3) * t570;
t424 = (qJD(1) * (t568 * t677 + t656) + t567 * t670) * t565 + t590 * t567 - t570 * t606;
t419 = qJDD(4) + t424;
t634 = pkin(1) * t648;
t514 = pkin(7) * t555 + t568 * t634;
t741 = t514 + t617 * (pkin(3) * t567 - pkin(9) * t570);
t607 = qJD(2) * t634;
t631 = pkin(1) * t642;
t671 = qJD(1) * qJD(2);
t651 = t571 * t671;
t632 = t565 * t651;
t650 = t565 * t670;
t638 = t568 * t607 - t571 * t631 + (t632 + t650) * pkin(7);
t446 = -pkin(2) * t606 + t638;
t702 = t565 * t571;
t725 = cos(qJ(1));
t630 = t714 * t725;
t724 = sin(qJ(1));
t520 = t568 * t724 - t571 * t630;
t629 = t714 * t724;
t522 = t568 * t725 + t571 * t629;
t731 = g(1) * t522 + g(2) * t520;
t585 = -g(3) * t702 + t731;
t740 = t446 - t585;
t679 = qJD(3) * t567;
t738 = t567 * t555 - t679;
t658 = t568 * t683;
t491 = t567 * t658 - t570 * t611;
t485 = qJD(4) + t491;
t521 = t568 * t630 + t571 * t724;
t660 = t565 * t725;
t467 = t521 * t570 - t567 * t660;
t566 = sin(qJ(4));
t569 = cos(qJ(4));
t433 = t467 * t566 - t520 * t569;
t434 = t467 * t569 + t520 * t566;
t493 = t567 * t611 + t570 * t658;
t737 = t617 * t493;
t661 = pkin(1) * t714;
t703 = t565 * t568;
t602 = -pkin(7) * t703 + t571 * t661;
t502 = -pkin(2) * t714 - t602;
t518 = t567 * t703 - t570 * t714;
t519 = t567 * t714 + t570 * t703;
t427 = t518 * pkin(3) - t519 * pkin(9) + t502;
t686 = pkin(7) * t702 + t568 * t661;
t503 = pkin(8) * t714 + t686;
t722 = pkin(8) * t568;
t609 = -pkin(2) * t571 - pkin(1) - t722;
t504 = t609 * t565;
t690 = t570 * t503 + t567 * t504;
t429 = -pkin(9) * t702 + t690;
t736 = t566 * t427 + t569 * t429;
t589 = t567 * t606;
t655 = t570 * t681;
t669 = qJDD(1) * t570;
t574 = t570 * t590 + (t568 * t669 + (-t568 * t679 + t655) * qJD(1)) * t565 + t589;
t652 = t568 * t671;
t633 = t565 * t652;
t668 = qJDD(1) * t571;
t554 = t565 * t668;
t666 = qJDD(3) - t554;
t587 = t633 + t666;
t457 = t569 * t493 - t566 * t617;
t676 = qJD(4) * t457;
t396 = t566 * t574 - t569 * t587 + t676;
t734 = (qJDD(2) + 0.2e1 * t642) * t565;
t635 = t570 * t555;
t733 = t635 - t677;
t511 = -pkin(7) * t658 + t571 * t634;
t627 = pkin(2) * t568 - pkin(8) * t571;
t512 = t627 * t683;
t689 = t570 * t511 + t567 * t512;
t432 = pkin(9) * t658 + t689;
t608 = pkin(3) * t570 + pkin(9) * t567 + pkin(2);
t674 = qJD(4) * t569;
t732 = t569 * t432 + t741 * t566 + t608 * t674;
t523 = -t568 * t629 + t571 * t725;
t659 = t565 * t724;
t470 = t523 * t567 - t570 * t659;
t646 = -t521 * t567 - t570 * t660;
t592 = g(1) * t470 - g(2) * t646 + g(3) * t518;
t476 = pkin(8) * t611 + t514;
t484 = qJD(1) * t504;
t422 = -t567 * t476 + t570 * t484;
t414 = pkin(3) * t617 - t422;
t645 = t569 * t617;
t455 = t493 * t566 + t645;
t388 = t455 * pkin(4) - t457 * qJ(5) + t414;
t721 = pkin(9) * t419;
t730 = t388 * t485 - t721;
t601 = t627 * qJD(2);
t513 = t565 * t601;
t515 = t602 * qJD(2);
t596 = -t503 * t679 + t504 * t677 + t567 * t513 + t570 * t515;
t682 = qJD(2) * t568;
t657 = t565 * t682;
t403 = pkin(9) * t657 + t596;
t462 = qJD(3) * t519 + t565 * t656;
t463 = -qJD(3) * t518 + t565 * t655;
t516 = t686 * qJD(2);
t407 = t462 * pkin(3) - t463 * pkin(9) + t516;
t728 = -qJD(4) * t736 - t403 * t566 + t407 * t569;
t727 = t457 ^ 2;
t726 = t485 ^ 2;
t572 = qJD(1) ^ 2;
t723 = pkin(4) * t419;
t715 = pkin(9) * qJD(4);
t713 = qJ(5) * t419;
t475 = -pkin(2) * t611 - t511;
t413 = t491 * pkin(3) - t493 * pkin(9) + t475;
t423 = t570 * t476 + t567 * t484;
t415 = -pkin(9) * t617 + t423;
t390 = t413 * t566 + t415 * t569;
t384 = qJ(5) * t485 + t390;
t712 = t384 * t485;
t711 = t390 * t485;
t675 = qJD(4) * t566;
t395 = qJD(4) * t645 + t493 * t675 - t566 * t587 - t569 * t574;
t710 = t395 * t566;
t709 = t455 * t485;
t708 = t457 * t455;
t707 = t457 * t485;
t562 = t565 ^ 2;
t704 = t562 * t572;
t701 = t566 * t419;
t700 = t566 * t570;
t699 = t569 * t419;
t698 = t569 * t570;
t697 = t569 * t571;
t447 = pkin(3) * t493 + pkin(9) * t491;
t696 = t569 * t422 + t566 * t447;
t673 = qJD(4) * t570;
t678 = qJD(3) * t569;
t695 = -qJD(5) * t570 + (-t566 * t673 - t567 * t678) * pkin(8) - t732 - t738 * qJ(5);
t692 = -t608 * t675 + (-t679 * pkin(8) - t432) * t566 + (pkin(8) * t673 + t741) * t569 + t738 * pkin(4);
t496 = t567 * t511;
t431 = -pkin(3) * t658 - t512 * t570 + t496;
t482 = t566 * t635 - t569 * t658;
t495 = (t566 * t568 + t570 * t697) * t565;
t483 = qJD(1) * t495;
t618 = pkin(4) * t566 - qJ(5) * t569;
t604 = pkin(8) + t618;
t619 = pkin(4) * t569 + qJ(5) * t566;
t691 = -pkin(4) * t482 + qJ(5) * t483 - t431 + (qJD(4) * t619 - qJD(5) * t569) * t567 + t604 * t677;
t688 = -qJD(5) * t566 + t485 * t618 - t423;
t685 = pkin(8) * t698 - t566 * t608;
t563 = t568 ^ 2;
t684 = -t571 ^ 2 + t563;
t680 = qJD(3) * t566;
t389 = t413 * t569 - t415 * t566;
t672 = qJD(5) - t389;
t664 = t571 * t704;
t663 = t566 * t702;
t662 = pkin(7) * t554 + t568 * t631 + t571 * t607;
t654 = t485 * t675;
t653 = 0.2e1 * pkin(1) * t562;
t647 = -t567 * t503 + t504 * t570;
t644 = t571 * t617;
t643 = t485 * t569;
t641 = qJD(3) * t617;
t584 = -pkin(7) * t633 + t662;
t445 = pkin(8) * t606 + t584;
t448 = (qJD(1) * t601 + qJDD(1) * t609) * t565;
t597 = -t570 * t445 - t567 * t448 + t476 * t679 - t484 * t677;
t381 = pkin(9) * t587 - t597;
t387 = t424 * pkin(3) - pkin(9) * t574 + t446;
t640 = t566 * t381 - t569 * t387 + t413 * t675 + t415 * t674;
t639 = t567 * t445 - t570 * t448 + t476 * t677 + t484 * t679;
t628 = t565 * t572 * t714;
t471 = t523 * t570 + t567 * t659;
t437 = t471 * t566 - t522 * t569;
t625 = -g(1) * t433 + g(2) * t437;
t438 = t471 * t569 + t522 * t566;
t624 = g(1) * t434 - g(2) * t438;
t623 = g(1) * t646 + g(2) * t470;
t428 = pkin(3) * t702 - t647;
t622 = t566 * t677 - t482;
t621 = t569 * t677 - t483;
t464 = t519 * t566 + t565 * t697;
t465 = t519 * t569 - t663;
t620 = -pkin(4) * t464 + qJ(5) * t465;
t383 = -pkin(4) * t485 + t672;
t616 = t383 * t569 - t384 * t566;
t614 = t427 * t569 - t429 * t566;
t610 = 0.2e1 * t648 + qJD(2);
t605 = pkin(3) + t619;
t603 = -t503 * t677 - t504 * t679 + t513 * t570 - t567 * t515;
t600 = -t485 * t674 - t701;
t599 = t414 * t485 - t721;
t598 = t651 + t670;
t372 = t569 * t381 + t566 * t387 + t413 * t674 - t415 * t675;
t595 = t569 * t403 + t566 * t407 + t427 * t674 - t429 * t675;
t450 = -t520 * t700 - t521 * t569;
t452 = -t522 * t700 - t523 * t569;
t494 = -t569 * t703 + t570 * t663;
t594 = g(1) * t452 + g(2) * t450 + g(3) * t494;
t451 = -t520 * t698 + t521 * t566;
t453 = -t522 * t698 + t523 * t566;
t593 = -g(1) * t453 - g(2) * t451 - g(3) * t495;
t591 = -g(1) * t471 - g(2) * t467 - g(3) * t519;
t588 = t606 * MDP(8);
t404 = -pkin(3) * t657 - t603;
t582 = -t485 * t715 + t592;
t581 = g(1) * t437 + g(2) * t433 + g(3) * t464 - t640;
t382 = -pkin(3) * t587 + t639;
t374 = t396 * pkin(4) + t395 * qJ(5) - t457 * qJD(5) + t382;
t580 = -t374 + t582;
t579 = -g(1) * t438 - g(2) * t434 - g(3) * t465 + t372;
t578 = t388 * t457 + qJDD(5) - t581;
t505 = t604 * t567;
t478 = t608 * t569 + (pkin(8) * t566 + pkin(4)) * t570;
t477 = -qJ(5) * t570 + t685;
t411 = -qJD(4) * t464 + t463 * t569 + t566 * t657;
t410 = -qJD(4) * t663 + t463 * t566 + t519 * t674 - t569 * t657;
t406 = pkin(4) * t457 + qJ(5) * t455;
t397 = t428 - t620;
t394 = -pkin(4) * t518 - t614;
t393 = qJ(5) * t518 + t736;
t392 = -pkin(4) * t493 + t422 * t566 - t447 * t569;
t391 = qJ(5) * t493 + t696;
t378 = -t395 + t709;
t377 = pkin(4) * t410 - qJ(5) * t411 - qJD(5) * t465 + t404;
t376 = -pkin(4) * t462 - t728;
t375 = qJ(5) * t462 + qJD(5) * t518 + t595;
t371 = qJDD(5) + t640 - t723;
t370 = qJD(5) * t485 + t372 + t713;
t1 = [(t382 * t464 + t389 * t462 + t428 * t396 + t404 * t455 + t414 * t410 + t614 * t419 + t485 * t728 - t518 * t640 + t624) * MDP(23) + (g(1) * t467 - g(2) * t471 + t422 * t657 + t502 * t424 + t446 * t518 + t475 * t462 + t516 * t491 + t587 * t647 - t603 * t617 + t639 * t702) * MDP(16) + (-t372 * t518 + t382 * t465 - t390 * t462 - t428 * t395 + t404 * t457 + t414 * t411 - t419 * t736 - t485 * t595 + t625) * MDP(24) + (t571 * t734 - t610 * t657) * MDP(7) + (t565 * t610 * t681 + t568 * t734) * MDP(6) + (t493 * t463 + t519 * t574) * MDP(11) + (-t519 * t424 - t493 * t462 - t463 * t491 - t518 * t574) * MDP(12) + ((qJDD(1) * t563 + 0.2e1 * t568 * t651) * MDP(4) + 0.2e1 * (t568 * t668 - t671 * t684) * MDP(5)) * t562 + (t370 * t393 + t384 * t375 + t374 * t397 + t388 * t377 + t371 * t394 + t383 * t376 - g(1) * (-pkin(1) * t724 - t521 * pkin(2) - pkin(3) * t467 - pkin(4) * t434 + pkin(7) * t660 - t520 * pkin(8) + pkin(9) * t646 - qJ(5) * t433) - g(2) * (pkin(1) * t725 + t523 * pkin(2) + t471 * pkin(3) + t438 * pkin(4) + pkin(7) * t659 + t522 * pkin(8) + t470 * pkin(9) + t437 * qJ(5))) * MDP(28) + (g(1) * t724 - g(2) * t725) * MDP(2) + (g(1) * t725 + g(2) * t724) * MDP(3) + t714 * t588 + (-g(1) * t520 + g(2) * t522 - t515 * t611 - t584 * t714 - t598 * t653 - t606 * t686) * MDP(10) + (-t516 * t611 + t602 * t606 - t638 * t714 + g(1) * t521 - g(2) * t523 + (-t652 + t668) * t653) * MDP(9) + (-t463 * t617 + t519 * t666 + ((-t589 + (-t590 - t632) * t570) * t571 + (-(-qJD(1) * t679 + t669) * t702 + (qJD(1) * t519 + t493) * qJD(2)) * t568) * t565) * MDP(13) + (-t423 * t657 + t446 * t519 + t475 * t463 + t516 * t493 + t502 * t574 - t587 * t690 + t596 * t617 - t597 * t702 + t623) * MDP(17) + (t462 * t617 - t518 * t666 + (t424 * t571 + (-qJD(1) * t518 - t491) * t682) * t565) * MDP(14) + (-t396 * t518 - t410 * t485 - t419 * t464 - t455 * t462) * MDP(21) + (-t395 * t518 + t411 * t485 + t419 * t465 + t457 * t462) * MDP(20) + (t419 * t518 + t462 * t485) * MDP(22) + (t395 * t464 - t396 * t465 - t410 * t457 - t411 * t455) * MDP(19) + (-t395 * t465 + t411 * t457) * MDP(18) + (-t666 * t571 + (-t555 - t617) * t682) * t565 * MDP(15) + (-t370 * t464 + t371 * t465 - t375 * t455 + t376 * t457 + t383 * t411 - t384 * t410 - t393 * t396 - t394 * t395 - t623) * MDP(26) + (-t371 * t518 + t374 * t464 - t376 * t485 + t377 * t455 - t383 * t462 + t388 * t410 - t394 * t419 + t396 * t397 + t624) * MDP(25) + (t370 * t518 - t374 * t465 + t375 * t485 - t377 * t457 + t384 * t462 - t388 * t411 + t393 * t419 + t395 * t397 - t625) * MDP(27) + qJDD(1) * MDP(1); (-t567 * t424 + t733 * t491 + t493 * t738 + t570 * t574) * MDP(12) + (t370 * t477 + t374 * t505 + t371 * t478 - g(1) * (pkin(4) * t453 + pkin(8) * t523 + qJ(5) * t452) - g(2) * (pkin(4) * t451 + pkin(8) * t521 + qJ(5) * t450) + t691 * t388 + t695 * t384 + t692 * t383 + t731 * t608 + (-pkin(4) * t495 - qJ(5) * t494 - (t571 * t608 + t722) * t565) * g(3)) * MDP(28) + (-t685 * t419 - t431 * t457 - t414 * t483 + t732 * t485 + (t414 * t678 + t372 + (qJD(3) * t457 + t654) * pkin(8)) * t570 + (-t414 * t675 + t382 * t569 + t617 * t390 + (t485 * t678 - t395) * pkin(8)) * t567 + t594) * MDP(24) + ((-qJD(3) * t658 + t606) * t567 ^ 2 + ((t565 * t598 + t590) * t567 - t737) * t570) * MDP(11) + (-pkin(2) * t424 - t496 * t617 - t422 * t658 - t514 * t491 + (-pkin(8) * t587 - t475 * t617) * t567 + (pkin(8) * t641 + t512 * t617 - t740) * t570) * MDP(16) + (-pkin(2) * t574 + t423 * t658 - t514 * t493 - t617 * t689 - t733 * t475 + (-t570 * t587 - t617 * t679) * pkin(8) + t740 * t567) * MDP(17) + (-t570 * t641 + t567 * t666 + (t570 * t644 + (qJD(2) * t567 - t493) * t568) * t683) * MDP(13) + (-t571 * t628 + t650) * MDP(6) + t588 + (t455 * t483 + t457 * t482 + (-t455 * t569 - t457 * t566) * t677 + (t710 - t396 * t569 + (t455 * t566 - t457 * t569) * qJD(4)) * t567) * MDP(19) + (pkin(1) * t664 + t511 * t611 + g(1) * t523 + g(2) * t521 + (pkin(7) * t671 + g(3)) * t703 - t662) * MDP(10) + (pkin(1) * t568 * t704 + t514 * t611 + t585 - t638) * MDP(9) + (t395 * t570 + t621 * t485 + (-t457 * t617 - t654 + t699) * t567) * MDP(20) + (t371 * t570 + t396 * t505 - t419 * t478 - t692 * t485 + t691 * t455 + t622 * t388 + (t374 * t566 + t383 * t617 + t388 * t674) * t567 + t593) * MDP(25) + (-t370 * t570 + t395 * t505 + t419 * t477 + t695 * t485 - t691 * t457 - t621 * t388 + (-t374 * t569 - t384 * t617 + t388 * t675) * t567 - t594) * MDP(27) + (-t383 * t483 + t384 * t482 - t395 * t478 - t396 * t477 + t692 * t457 - t695 * t455 + t616 * t677 + (-t370 * t566 + t371 * t569 + (-t383 * t566 - t384 * t569) * qJD(4) + t585) * t567) * MDP(26) + (-t395 * t567 * t569 + (-t567 * t675 + t621) * t457) * MDP(18) + (t567 * t641 + t570 * t666 + (-t567 * t644 + (qJD(2) * t570 + t491) * t568) * t683) * MDP(14) + (-t608 * t699 - t414 * t482 - t431 * t455 + (-t741 * t569 + (qJD(4) * t608 + t432) * t566) * t485 + (t414 * t680 + t640 + (qJD(3) * t455 + t600) * pkin(8)) * t570 + (t414 * t674 + t382 * t566 - t617 * t389 + (t485 * t680 + t396) * pkin(8)) * t567 + t593) * MDP(23) - t568 * MDP(4) * t664 + t617 * MDP(15) * t658 + t684 * MDP(5) * t704 + (-t485 * t567 * t617 - t419 * t570) * MDP(22) + (t396 * t570 - t622 * t485 + (t455 * t617 + t600) * t567) * MDP(21) + (t568 * t628 + t554) * MDP(7); -t491 ^ 2 * MDP(12) + (-t491 * t617 + t574) * MDP(13) + (-t424 - t737) * MDP(14) + t587 * MDP(15) + (-t423 * t617 + t592 - t639) * MDP(16) + (-t422 * t617 + t475 * t491 - t591 + t597) * MDP(17) + (t457 * t643 - t710) * MDP(18) + ((-t395 - t709) * t569 + (-t396 - t707) * t566) * MDP(19) + (t485 * t643 + t701) * MDP(20) + (-t566 * t726 + t699) * MDP(21) + (-pkin(3) * t396 - t423 * t455 + (t422 * t485 + t599) * t566 + (-t382 + (-t447 - t715) * t485 + t592) * t569) * MDP(23) + (pkin(3) * t395 + t696 * t485 - t423 * t457 + t599 * t569 + (t382 - t582) * t566) * MDP(24) + (t392 * t485 - t396 * t605 + t688 * t455 + t566 * t730 + t580 * t569) * MDP(25) + (t391 * t455 - t392 * t457 + (t370 + t485 * t383 + (-t396 + t676) * pkin(9)) * t569 + (t371 - t712 + (qJD(4) * t455 - t395) * pkin(9)) * t566 + t591) * MDP(26) + (-t391 * t485 - t395 * t605 - t688 * t457 + t580 * t566 - t730 * t569) * MDP(27) + (-t383 * t392 - t384 * t391 + t688 * t388 + (qJD(4) * t616 + t370 * t569 + t371 * t566 + t591) * pkin(9) + (-t374 + t592) * t605) * MDP(28) + (MDP(11) * t491 + MDP(12) * t493 - t475 * MDP(16) - t457 * MDP(20) + t455 * MDP(21) - t485 * MDP(22) - t389 * MDP(23) + t390 * MDP(24) + MDP(25) * t383 - MDP(27) * t384) * t493; MDP(18) * t708 + (-t455 ^ 2 + t727) * MDP(19) + t378 * MDP(20) + (t707 - t396) * MDP(21) + t419 * MDP(22) + (-t414 * t457 + t581 + t711) * MDP(23) + (t389 * t485 + t414 * t455 - t579) * MDP(24) + (-t406 * t455 - t578 + t711 + 0.2e1 * t723) * MDP(25) + (pkin(4) * t395 - qJ(5) * t396 + (t384 - t390) * t457 + (t383 - t672) * t455) * MDP(26) + (0.2e1 * t713 - t388 * t455 + t406 * t457 + (0.2e1 * qJD(5) - t389) * t485 + t579) * MDP(27) + (t370 * qJ(5) - t371 * pkin(4) - t388 * t406 - t383 * t390 - g(1) * (-pkin(4) * t437 + qJ(5) * t438) - g(2) * (-pkin(4) * t433 + qJ(5) * t434) - g(3) * t620 + t672 * t384) * MDP(28); t378 * MDP(26) + (-t726 - t727) * MDP(27) + (t578 - t712 - t723) * MDP(28) + (t708 - t419) * MDP(25);];
tau = t1;
