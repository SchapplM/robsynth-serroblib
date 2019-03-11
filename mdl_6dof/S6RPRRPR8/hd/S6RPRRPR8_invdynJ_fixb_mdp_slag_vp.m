% Calculate vector of inverse dynamics joint torques for
% S6RPRRPR8
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
%   pkin=[a2,a3,a4,a5,a6,d1,d3,d4,d6,theta5]';
% MDP [29x1]
%   Minimal dynamic parameter vector (fixed base model)
%   see S6RPRRPR8_convert_par2_MPV_fixb.m
% 
% Output:
% tau [6x1]
%   joint torques of inverse dynamics (contains inertial, gravitational Coriolis and centrifugal forces)

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 05:25
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tau = S6RPRRPR8_invdynJ_fixb_mdp_slag_vp(qJ, qJD, qJDD, g, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(3,1),zeros(10,1),zeros(29,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRRPR8_invdynJ_fixb_mdp_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RPRRPR8_invdynJ_fixb_mdp_slag_vp: qJD has to be [6x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [6 1]), ...
  'S6RPRRPR8_invdynJ_fixb_mdp_slag_vp: qJDD has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RPRRPR8_invdynJ_fixb_mdp_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RPRRPR8_invdynJ_fixb_mdp_slag_vp: pkin has to be [10x1] (double)');
assert(isreal(MDP) && all(size(MDP) == [29 1]), ...
  'S6RPRRPR8_invdynJ_fixb_mdp_slag_vp: MDP has to be [29x1] (double)'); 

%% Symbolic Calculation
% From invdyn_joint_fixb_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 05:25:07
% EndTime: 2019-03-09 05:25:21
% DurationCPUTime: 9.42s
% Computational Cost: add. (5827->552), mult. (11858->731), div. (0->0), fcn. (8231->12), ass. (0->243)
t575 = cos(qJ(6));
t572 = sin(qJ(4));
t576 = cos(qJ(4));
t644 = t576 * qJD(3);
t577 = cos(qJ(3));
t656 = qJD(1) * t577;
t528 = t572 * t656 - t644;
t653 = qJD(3) * t572;
t530 = t576 * t656 + t653;
t568 = sin(pkin(10));
t569 = cos(pkin(10));
t605 = -t528 * t569 - t530 * t568;
t672 = t575 * t605;
t463 = t528 * t568 - t530 * t569;
t571 = sin(qJ(6));
t689 = t463 * t571;
t411 = t672 + t689;
t573 = sin(qJ(3));
t657 = qJD(1) * t573;
t550 = qJD(4) + t657;
t545 = qJD(6) + t550;
t691 = t411 * t545;
t534 = pkin(3) * t573 - pkin(8) * t577 + qJ(2);
t507 = t534 * qJD(1);
t579 = -pkin(1) - pkin(7);
t547 = qJD(1) * t579 + qJD(2);
t533 = t573 * t547;
t514 = qJD(3) * pkin(8) + t533;
t453 = t507 * t572 + t514 * t576;
t630 = t573 * t644;
t647 = qJD(4) * t577;
t590 = -t572 * t647 - t630;
t640 = qJDD(1) * t577;
t456 = qJD(1) * t590 + qJD(4) * t644 + t572 * qJDD(3) + t576 * t640;
t615 = pkin(3) * t577 + pkin(8) * t573;
t522 = qJD(3) * t615 + qJD(2);
t469 = qJD(1) * t522 + qJDD(1) * t534;
t459 = t576 * t469;
t543 = qJDD(1) * t579 + qJDD(2);
t651 = qJD(3) * t577;
t474 = qJDD(3) * pkin(8) + t543 * t573 + t547 * t651;
t643 = qJD(1) * qJD(3);
t626 = t577 * t643;
t641 = qJDD(1) * t573;
t521 = qJDD(4) + t626 + t641;
t373 = pkin(4) * t521 - qJ(5) * t456 - qJD(4) * t453 - qJD(5) * t530 - t474 * t572 + t459;
t652 = qJD(3) * t573;
t631 = t572 * t652;
t457 = -qJD(1) * t631 + qJD(4) * t530 - t576 * qJDD(3) + t572 * t640;
t648 = qJD(4) * t576;
t635 = -t572 * t469 - t576 * t474 - t507 * t648;
t649 = qJD(4) * t572;
t593 = -t514 * t649 - t635;
t378 = -qJ(5) * t457 - qJD(5) * t528 + t593;
t363 = t569 * t373 - t378 * t568;
t403 = t456 * t569 - t457 * t568;
t361 = pkin(5) * t521 - pkin(9) * t403 + t363;
t364 = t568 * t373 + t569 * t378;
t402 = -t456 * t568 - t457 * t569;
t362 = pkin(9) * t402 + t364;
t452 = t576 * t507 - t514 * t572;
t427 = -qJ(5) * t530 + t452;
t420 = pkin(4) * t550 + t427;
t428 = -qJ(5) * t528 + t453;
t680 = t569 * t428;
t387 = t568 * t420 + t680;
t706 = pkin(9) * t605;
t380 = t387 + t706;
t645 = qJD(6) * t571;
t379 = t380 * t645;
t682 = t547 * t577;
t515 = -qJD(3) * pkin(3) - t682;
t471 = pkin(4) * t528 + qJD(5) + t515;
t417 = -pkin(5) * t605 + t471;
t560 = qJ(4) + pkin(10) + qJ(6);
t552 = sin(t560);
t553 = cos(t560);
t578 = cos(qJ(1));
t574 = sin(qJ(1));
t677 = t573 * t574;
t480 = t552 * t578 + t553 * t677;
t675 = t573 * t578;
t482 = -t552 * t574 + t553 * t675;
t695 = g(3) * t577;
t717 = g(1) * t480 - g(2) * t482 - t571 * t361 - t575 * t362 - t417 * t411 + t553 * t695 + t379;
t517 = qJDD(6) + t521;
t707 = -t575 * t463 + t571 * t605;
t716 = t517 * MDP(27) + (-t411 ^ 2 + t707 ^ 2) * MDP(24) - t411 * MDP(23) * t707;
t692 = t707 * t545;
t532 = t615 * qJD(1);
t513 = t576 * t532;
t570 = -qJ(5) - pkin(8);
t623 = qJD(4) * t570;
t676 = t573 * t576;
t678 = t572 * t577;
t714 = t547 * t678 - t513 - (pkin(4) * t577 + qJ(5) * t676) * qJD(1) - qJD(5) * t572 + t576 * t623;
t632 = t572 * t657;
t646 = qJD(5) * t576;
t671 = t576 * t577;
t662 = t572 * t532 + t547 * t671;
t713 = qJ(5) * t632 - t572 * t623 - t646 + t662;
t622 = -qJDD(3) * pkin(3) + t547 * t652;
t683 = t543 * t577;
t473 = t622 - t683;
t696 = g(3) * t573;
t703 = -g(1) * t574 + g(2) * t578;
t587 = -t577 * t703 - t696;
t712 = qJD(4) * pkin(8) * t550 + t473 + t587;
t479 = -t552 * t677 + t553 * t578;
t481 = t552 * t675 + t553 * t574;
t621 = t575 * t361 - t571 * t362;
t711 = -g(1) * t479 - g(2) * t481 - t417 * t707 + t552 * t695 + t621;
t710 = pkin(9) * t463;
t524 = t568 * t576 + t569 * t572;
t504 = t524 * qJD(4);
t505 = t524 * qJD(1);
t709 = t573 * t505 + t504;
t679 = t569 * t576;
t602 = t568 * t572 - t679;
t701 = qJD(4) * t602;
t708 = -t568 * t632 + t657 * t679 - t701;
t665 = t568 * t713 + t569 * t714;
t664 = t568 * t714 - t569 * t713;
t661 = t572 * t534 + t579 * t676;
t670 = t576 * t578;
t508 = -t572 * t677 + t670;
t674 = t574 * t576;
t510 = t572 * t675 + t674;
t702 = -g(1) * t508 - g(2) * t510;
t620 = -t575 * t402 + t403 * t571;
t368 = qJD(6) * t707 + t620;
t700 = pkin(4) * t568;
t694 = pkin(1) * qJDD(1);
t581 = qJD(1) ^ 2;
t693 = qJ(2) * t581;
t690 = t456 * t572;
t688 = t521 * t572;
t687 = t521 * t576;
t686 = t528 * t550;
t685 = t530 * t550;
t684 = t530 * t576;
t681 = t550 * t572;
t423 = t568 * t428;
t386 = t569 * t420 - t423;
t377 = pkin(5) * t550 + t386 + t710;
t673 = t575 * t377;
t502 = t576 * t522;
t625 = -t572 * t579 + pkin(4);
t401 = qJ(5) * t630 + t502 - t661 * qJD(4) + (qJ(5) * t649 + qJD(3) * t625 - t646) * t577;
t628 = t576 * t647;
t650 = qJD(3) * t579;
t629 = t577 * t650;
t634 = t572 * t522 + t534 * t648 + t576 * t629;
t405 = -qJ(5) * t628 + (-qJD(5) * t577 + (qJ(5) * qJD(3) - qJD(4) * t579) * t573) * t572 + t634;
t375 = t568 * t401 + t569 * t405;
t606 = -t524 * t571 - t575 * t602;
t669 = qJD(6) * t606 - t571 * t709 + t575 * t708;
t462 = t524 * t575 - t571 * t602;
t668 = qJD(6) * t462 + t571 * t708 + t575 * t709;
t391 = t569 * t427 - t423;
t493 = t524 * t577;
t667 = t602 * qJD(1) - qJD(3) * t493 + t573 * t701;
t495 = t602 * t577;
t666 = -qJD(3) * t495 - t504 * t573 - t505;
t520 = t576 * t534;
t460 = -qJ(5) * t671 + t573 * t625 + t520;
t472 = -qJ(5) * t678 + t661;
t414 = t568 * t460 + t569 * t472;
t663 = -t533 + t709 * pkin(5) + (t632 + t649) * pkin(4);
t540 = t570 * t572;
t541 = t570 * t576;
t476 = t568 * t540 - t569 * t541;
t567 = t577 ^ 2;
t660 = t573 ^ 2 - t567;
t580 = qJD(3) ^ 2;
t659 = -t580 - t581;
t655 = qJD(3) * t528;
t654 = qJD(3) * t530;
t642 = qJDD(1) * qJ(2);
t639 = qJDD(3) * t573;
t637 = 0.2e1 * qJD(1) * qJD(2);
t636 = qJD(6) * t672 + t571 * t402 + t575 * t403;
t556 = pkin(4) * t576 + pkin(3);
t633 = pkin(4) * t572 + pkin(7);
t627 = g(2) * (t578 * pkin(1) + t574 * qJ(2));
t374 = t569 * t401 - t405 * t568;
t390 = -t427 * t568 - t680;
t413 = t569 * t460 - t472 * t568;
t619 = t550 * t579 + t514;
t475 = t569 * t540 + t541 * t568;
t618 = pkin(4) * t678 - t577 * t579;
t617 = -qJD(4) * t507 - t474;
t616 = qJD(4) * t573 + qJD(1);
t614 = g(1) * t578 + g(2) * t574;
t612 = qJDD(2) + t703;
t447 = -pkin(9) * t602 + t476;
t611 = pkin(5) * t656 + pkin(9) * t708 + qJD(6) * t447 - t665;
t446 = -pkin(9) * t524 + t475;
t610 = pkin(9) * t709 - qJD(6) * t446 - t664;
t494 = t602 * t573;
t609 = -qJD(6) * t494 - t667;
t492 = t524 * t573;
t608 = qJD(6) * t492 - t666;
t366 = t571 * t377 + t575 * t380;
t607 = -t575 * t493 + t495 * t571;
t436 = -t493 * t571 - t495 * t575;
t603 = t556 * t573 + t570 * t577;
t601 = -t543 - t703;
t554 = pkin(4) * t569 + pkin(5);
t600 = t554 * t571 + t575 * t700;
t599 = t554 * t575 - t571 * t700;
t597 = t550 * t648 + t688;
t596 = -t550 * t649 + t687;
t595 = t573 * t650 + (t628 - t631) * pkin(4);
t594 = 0.2e1 * qJ(2) * t643 + qJDD(3) * t579;
t367 = t463 * t645 + t636;
t591 = pkin(4) * t457 + qJDD(5) + t622;
t589 = t601 + t693;
t588 = -pkin(8) * t521 + t515 * t550;
t415 = t591 - t683;
t584 = -t614 + t637 + 0.2e1 * t642;
t582 = -t579 * t580 + t584;
t563 = t578 * qJ(2);
t559 = qJDD(3) * t577;
t511 = -t572 * t574 + t573 * t670;
t509 = t572 * t578 + t573 * t674;
t486 = pkin(5) * t602 - t556;
t468 = pkin(5) * t493 + t618;
t443 = t524 * t647 - t568 * t631 + t569 * t630;
t441 = t524 * t652 + t577 * t701;
t429 = pkin(4) * t530 - pkin(5) * t463;
t416 = -pkin(5) * t441 + t595;
t393 = -pkin(9) * t493 + t414;
t392 = pkin(5) * t573 + pkin(9) * t495 + t413;
t385 = qJD(6) * t436 - t575 * t441 - t443 * t571;
t384 = qJD(6) * t607 + t441 * t571 - t443 * t575;
t383 = t391 + t710;
t382 = t390 - t706;
t381 = -pkin(5) * t402 + t415;
t370 = pkin(9) * t441 + t375;
t369 = pkin(5) * t651 + pkin(9) * t443 + t374;
t365 = -t380 * t571 + t673;
t1 = [(-t368 * t573 - t385 * t545 + t411 * t651 + t517 * t607) * MDP(26) + ((t369 * t575 - t370 * t571) * t545 + (t392 * t575 - t393 * t571) * t517 + t621 * t573 + t365 * t651 - t416 * t411 + t468 * t368 - t381 * t607 + t417 * t385 - g(1) * t482 - g(2) * t480 + ((-t392 * t571 - t393 * t575) * t545 - t366 * t573) * qJD(6)) * MDP(28) + (t367 * t607 - t368 * t436 + t384 * t411 - t385 * t707) * MDP(24) + t584 * MDP(5) + (t363 * t495 - t364 * t493 + t374 * t463 + t375 * t605 + t386 * t443 + t387 * t441 + t402 * t414 - t403 * t413 + t577 * t614) * MDP(21) + qJDD(1) * MDP(1) + (-t573 * t580 + t559) * MDP(9) - t703 * MDP(2) + (t367 * t436 + t384 * t707) * MDP(23) + (-t366 * t651 + g(1) * t481 - g(2) * t479 + t468 * t367 + t379 * t573 + t381 * t436 + t417 * t384 + t416 * t707 + (-(-qJD(6) * t393 + t369) * t545 - t392 * t517 - t361 * t573) * t571 + (-(qJD(6) * t392 + t370) * t545 - t393 * t517 - (qJD(6) * t377 + t362) * t573) * t575) * MDP(29) + (t367 * t573 + t384 * t545 + t436 * t517 + t651 * t707) * MDP(25) + (-(qJDD(2) - t694) * pkin(1) - g(1) * (-pkin(1) * t574 + t563) - t627 + (t637 + t642) * qJ(2)) * MDP(6) + (t612 - 0.2e1 * t694) * MDP(4) + ((t528 * t576 + t530 * t572) * t652 + (-t690 - t457 * t576 + (t528 * t572 - t684) * qJD(4)) * t577) * MDP(15) + (t456 * t671 + t530 * t590) * MDP(14) + 0.2e1 * (-t573 * t640 + t643 * t660) * MDP(8) + (-t634 * t550 - t661 * t521 + g(1) * t510 - g(2) * t508 + (t619 * t649 + (-t515 * t576 + t530 * t579) * qJD(3) + t635) * t573 + (-qJD(3) * t453 - t456 * t579 + t473 * t576 - t515 * t649) * t577) * MDP(20) + (t521 * t573 + t550 * t651) * MDP(18) + (t517 * t573 + t545 * t651) * MDP(27) + ((-t550 * t644 + t456) * t573 + (t596 + t654) * t577) * MDP(16) + ((t550 * t653 - t457) * t573 + (-t597 - t655) * t577) * MDP(17) + (-g(1) * t511 - g(2) * t509 + t502 * t550 + t520 * t521 + (t528 * t650 - t619 * t648 + t459) * t573 + (qJD(3) * t452 - t457 * t579 + t515 * t648) * t577 + ((-qJD(4) * t534 - t629) * t550 + t473 * t577 + (-qJD(3) * t515 - t521 * t579 + t617) * t573) * t572) * MDP(19) + (-t577 * t580 - t639) * MDP(10) + (t364 * t414 + t387 * t375 + t363 * t413 + t386 * t374 + t415 * t618 + t471 * t595 - g(1) * t563 - t627 + (-g(1) * t603 - g(2) * t633) * t578 + (-g(1) * (-pkin(1) - t633) - g(2) * t603) * t574) * MDP(22) + (qJDD(1) * t567 - 0.2e1 * t573 * t626) * MDP(7) + (t573 * t582 + t577 * t594) * MDP(12) + (-t573 * t594 + t577 * t582) * MDP(13) + t614 * MDP(3); qJDD(1) * MDP(4) - t581 * MDP(5) + (t612 - t693 - t694) * MDP(6) + (t573 * t659 + t559) * MDP(12) + (t577 * t659 - t639) * MDP(13) + (-t457 * t577 + (t655 - t688) * t573 + (-t572 * t651 - t576 * t616) * t550) * MDP(19) + (-t456 * t577 + (t654 - t687) * t573 + (t572 * t616 - t577 * t644) * t550) * MDP(20) + (-t402 * t494 + t403 * t492 + t463 * t667 + t605 * t666) * MDP(21) + (-t363 * t492 - t364 * t494 + t386 * t667 + t387 * t666 - t415 * t577 + t471 * t652 + t703) * MDP(22) + ((-t492 * t575 + t494 * t571) * t517 - t411 * t652 - t577 * t368 + (t571 * t608 - t575 * t609) * t545) * MDP(28) + (-(-t492 * t571 - t494 * t575) * t517 + t707 * t652 - t577 * t367 + (t571 * t609 + t575 * t608) * t545) * MDP(29); MDP(9) * t640 - MDP(10) * t641 + qJDD(3) * MDP(11) + (-t577 * t589 + t696) * MDP(12) + (t573 * t589 + t695) * MDP(13) + (t550 * t684 + t690) * MDP(14) + ((t456 - t686) * t576 + (-t457 - t685) * t572) * MDP(15) + ((-t530 * t577 + t550 * t676) * qJD(1) + t597) * MDP(16) + ((t528 * t577 - t573 * t681) * qJD(1) + t596) * MDP(17) + (-t528 * t533 - pkin(3) * t457 - t513 * t550 + (t550 * t682 + t588) * t572 - t712 * t576) * MDP(19) + (-pkin(3) * t456 - t530 * t533 + t662 * t550 + t572 * t712 + t588 * t576) * MDP(20) + (-t363 * t524 - t364 * t602 - t386 * t708 - t387 * t709 + t402 * t476 - t403 * t475 + t665 * t463 + t573 * t703 + t664 * t605 - t695) * MDP(21) + (t364 * t476 + t363 * t475 - t415 * t556 + g(3) * t603 + (pkin(4) * t681 - t533) * t471 + t664 * t387 + t665 * t386 + t703 * (t556 * t577 - t570 * t573)) * MDP(22) + (t367 * t462 + t669 * t707) * MDP(23) + (t367 * t606 - t368 * t462 + t411 * t669 - t668 * t707) * MDP(24) + (t462 * t517 + t545 * t669) * MDP(25) + (t517 * t606 - t545 * t668) * MDP(26) + ((t446 * t575 - t447 * t571) * t517 + t486 * t368 - t381 * t606 + (t571 * t610 - t575 * t611) * t545 + t668 * t417 - t663 * t411 - t587 * t553) * MDP(28) + (-(t446 * t571 + t447 * t575) * t517 + t486 * t367 + t381 * t462 + (t571 * t611 + t575 * t610) * t545 + t669 * t417 + t663 * t707 + t587 * t552) * MDP(29) + (-t550 * MDP(18) - MDP(19) * t452 + t453 * MDP(20) - MDP(25) * t707 - MDP(26) * t411 - t545 * MDP(27) - t365 * MDP(28) + t366 * MDP(29)) * t656 + (MDP(7) * t573 * t577 - MDP(8) * t660) * t581; t530 * t528 * MDP(14) + (-t528 ^ 2 + t530 ^ 2) * MDP(15) + (t456 + t686) * MDP(16) + (-t457 + t685) * MDP(17) + t521 * MDP(18) + (-t514 * t648 + t453 * t550 - t515 * t530 + t459 + (t617 + t695) * t572 + t702) * MDP(19) + (g(1) * t509 - g(2) * t511 + g(3) * t671 + t452 * t550 + t515 * t528 - t593) * MDP(20) + ((t402 * t568 - t403 * t569) * pkin(4) + (t386 - t391) * t605 + (-t387 - t390) * t463) * MDP(21) + (-t386 * t390 - t387 * t391 + (g(3) * t678 + t363 * t569 + t364 * t568 - t471 * t530 + t702) * pkin(4)) * MDP(22) + (t367 - t691) * MDP(25) + (-t368 + t692) * MDP(26) + (t599 * t517 - (t382 * t575 - t383 * t571) * t545 + t429 * t411 + (-t545 * t600 - t366) * qJD(6) + t711) * MDP(28) + (-t600 * t517 + (t382 * t571 + t383 * t575) * t545 - t429 * t707 + (-t545 * t599 - t673) * qJD(6) + t717) * MDP(29) + t716; (-t463 ^ 2 - t605 ^ 2) * MDP(21) + (t368 + t692) * MDP(28) + (t367 + t691) * MDP(29) + (-t386 * t463 - t387 * t605 + t577 * t601 + t591 - t696) * MDP(22); (t636 - t691) * MDP(25) + (-t620 + t692) * MDP(26) + (t366 * t545 + t711) * MDP(28) + (t365 * t545 + t717) * MDP(29) + (MDP(25) * t689 - MDP(26) * t707 - MDP(28) * t366 - MDP(29) * t673) * qJD(6) + t716;];
tau  = t1;
