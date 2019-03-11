% Calculate vector of inverse dynamics joint torques for
% S6PRRPPR4
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
%   pkin=[a2,a3,a4,a5,a6,alpha2,d2,d3,d6,theta1,theta4]';
% MDP [26x1]
%   Minimal dynamic parameter vector (fixed base model)
%   see S6PRRPPR4_convert_par2_MPV_fixb.m
% 
% Output:
% tau [6x1]
%   joint torques of inverse dynamics (contains inertial, gravitational Coriolis and centrifugal forces)

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-08 21:17
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tau = S6PRRPPR4_invdynJ_fixb_mdp_slag_vp(qJ, qJD, qJDD, g, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(3,1),zeros(11,1),zeros(26,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRRPPR4_invdynJ_fixb_mdp_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6PRRPPR4_invdynJ_fixb_mdp_slag_vp: qJD has to be [6x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [6 1]), ...
  'S6PRRPPR4_invdynJ_fixb_mdp_slag_vp: qJDD has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6PRRPPR4_invdynJ_fixb_mdp_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6PRRPPR4_invdynJ_fixb_mdp_slag_vp: pkin has to be [11x1] (double)');
assert(isreal(MDP) && all(size(MDP) == [26 1]), ...
  'S6PRRPPR4_invdynJ_fixb_mdp_slag_vp: MDP has to be [26x1] (double)'); 

%% Symbolic Calculation
% From invdyn_joint_fixb_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-08 21:16:41
% EndTime: 2019-03-08 21:16:51
% DurationCPUTime: 8.74s
% Computational Cost: add. (3429->594), mult. (8147->780), div. (0->0), fcn. (6473->12), ass. (0->249)
t547 = sin(qJ(3));
t704 = qJ(4) * t547 + pkin(2);
t703 = MDP(14) + MDP(17);
t550 = cos(qJ(3));
t505 = -pkin(3) * t550 - t704;
t596 = pkin(3) * t547 - qJ(4) * t550;
t499 = t596 * qJD(2);
t543 = sin(pkin(11));
t545 = cos(pkin(11));
t548 = sin(qJ(2));
t544 = sin(pkin(6));
t652 = qJD(1) * t544;
t626 = t548 * t652;
t500 = qJD(2) * pkin(8) + t626;
t683 = cos(pkin(6));
t613 = qJD(1) * t683;
t567 = t547 * t500 - t550 * t613;
t396 = t543 * t499 - t545 * t567;
t650 = qJD(2) * t547;
t389 = qJ(5) * t650 + t396;
t645 = qJD(4) * t545;
t702 = -t389 + t645;
t474 = qJD(3) * t596 - qJD(4) * t547;
t551 = cos(qJ(2));
t664 = t550 * t551;
t633 = t544 * t664;
t608 = t543 * t633;
t701 = qJD(1) * t608 + (t474 - t626) * t545;
t539 = t550 * qJDD(2);
t637 = qJD(2) * qJD(3);
t618 = t547 * t637;
t575 = -t618 + t539;
t699 = MDP(12) + MDP(16);
t698 = MDP(13) - MDP(18);
t697 = qJDD(3) * pkin(3) - qJDD(4);
t451 = (t543 * t548 + t545 * t664) * t544;
t443 = qJD(1) * t451;
t647 = qJD(3) * t547;
t696 = qJ(5) * t647 - qJD(5) * t550 - t443;
t682 = cos(pkin(10));
t598 = t683 * t682;
t681 = sin(pkin(10));
t477 = t548 * t598 + t551 * t681;
t615 = t544 * t682;
t424 = t477 * t547 + t550 * t615;
t597 = t683 * t681;
t479 = -t548 * t597 + t551 * t682;
t614 = t544 * t681;
t426 = t479 * t547 - t550 * t614;
t668 = t544 * t548;
t483 = t547 * t668 - t550 * t683;
t570 = g(1) * t426 + g(2) * t424 + g(3) * t483;
t552 = qJD(3) ^ 2;
t638 = qJD(1) * qJD(2);
t619 = t548 * t638;
t667 = t544 * t551;
t593 = -qJDD(1) * t667 + t544 * t619;
t476 = t548 * t681 - t551 * t598;
t478 = t548 * t682 + t551 * t597;
t601 = g(1) * t478 + g(2) * t476;
t695 = 0.2e1 * qJDD(2) * pkin(2) - pkin(8) * t552 + t544 * (-g(3) * t551 + t619) - t593 + t601;
t639 = t545 * qJD(3);
t492 = t543 * t650 - t639;
t648 = qJD(3) * t543;
t494 = t545 * t650 + t648;
t546 = sin(qJ(6));
t549 = cos(qJ(6));
t413 = t492 * t546 + t494 * t549;
t537 = t545 * qJDD(3);
t617 = t550 * t637;
t636 = qJDD(2) * t547;
t576 = t617 + t636;
t448 = t543 * t576 - t537;
t449 = qJDD(3) * t543 + t545 * t576;
t611 = -t549 * t448 + t546 * t449;
t363 = qJD(6) * t413 + t611;
t457 = qJDD(2) * pkin(8) + (qJDD(1) * t548 + t551 * t638) * t544;
t523 = t547 * t613;
t610 = qJDD(1) * t683;
t646 = qJD(3) * t550;
t606 = -qJD(3) * t523 - t547 * t457 - t500 * t646 + t550 * t610;
t563 = t606 + t697;
t644 = qJD(5) * t494;
t678 = qJ(5) * t449;
t555 = t563 + t644 + t678;
t693 = pkin(4) + pkin(5);
t356 = -t448 * t693 + t555;
t694 = t356 + t570;
t490 = t494 ^ 2;
t691 = pkin(4) * t448;
t690 = pkin(4) * t545;
t689 = pkin(9) * t547;
t685 = -pkin(9) + qJ(4);
t684 = qJD(2) * pkin(2);
t680 = qJ(4) * t545;
t677 = qJ(5) * t545;
t671 = t494 * t546;
t411 = -t549 * t492 + t671;
t649 = qJD(2) * t550;
t531 = qJD(6) + t649;
t674 = t411 * t531;
t673 = t413 * t531;
t670 = t543 * t549;
t669 = t543 * t550;
t666 = t545 * t547;
t665 = t545 * t550;
t663 = qJDD(1) - g(3);
t599 = t547 * t610;
t371 = t599 + qJDD(3) * qJ(4) + t550 * t457 + (qJD(4) - t567) * qJD(3);
t386 = qJD(2) * t474 + qJDD(2) * t505 + t593;
t361 = t545 * t371 + t543 * t386;
t628 = -pkin(8) * t543 - pkin(4);
t635 = pkin(9) * t665;
t662 = (-t635 + (-pkin(5) + t628) * t547) * qJD(3) - t701;
t460 = t543 * t474;
t661 = -t460 - (-pkin(8) * t666 + pkin(9) * t669) * qJD(3) - t696;
t634 = pkin(8) * t647;
t421 = -t545 * t634 + t460;
t660 = t421 + t696;
t659 = t628 * t647 - t701;
t446 = t550 * t500 + t523;
t437 = qJD(3) * qJ(4) + t446;
t625 = t551 * t652;
t447 = qJD(2) * t505 - t625;
t380 = t545 * t437 + t543 * t447;
t658 = t543 * t634 + t701;
t657 = t421 - t443;
t497 = t543 * t546 + t545 * t549;
t577 = t497 * t550;
t656 = -qJD(2) * t577 - t497 * qJD(6);
t624 = t543 * t649;
t632 = t546 * t665;
t640 = qJD(6) * t549;
t641 = qJD(6) * t546;
t655 = -qJD(2) * t632 + t543 * t640 - t545 * t641 + t549 * t624;
t620 = qJ(4) * t539;
t654 = qJD(4) * t624 + t543 * t620;
t453 = pkin(8) * t665 + t543 * t505;
t541 = t547 ^ 2;
t542 = t550 ^ 2;
t653 = t541 - t542;
t651 = qJD(2) * t544;
t643 = qJD(5) * t543;
t631 = t546 * t448 + t549 * t449 + t492 * t640;
t630 = t505 * t476;
t629 = t505 * t478;
t627 = qJ(4) * t647;
t623 = t548 * t651;
t622 = t551 * t651;
t621 = qJD(5) * t666;
t616 = qJ(5) * t543 + pkin(3);
t360 = -t543 * t371 + t386 * t545;
t584 = pkin(4) * t539 + qJDD(5) - t360;
t358 = -pkin(4) * t618 + t584;
t352 = pkin(5) * t575 - pkin(9) * t449 + t358;
t357 = qJ(5) * t618 + (-qJ(5) * qJDD(2) - qJD(2) * qJD(5)) * t550 + t361;
t355 = pkin(9) * t448 + t357;
t612 = t549 * t352 - t546 * t355;
t379 = -t543 * t437 + t447 * t545;
t395 = t499 * t545 + t543 * t567;
t527 = pkin(8) * t669;
t452 = t505 * t545 - t527;
t583 = -t543 * t693 + t677;
t609 = -t583 * t649 + t446 + t643;
t607 = pkin(3) * t633 + pkin(8) * t668 + t667 * t704;
t605 = t492 * t625;
t604 = t494 * t625;
t603 = t547 * t625;
t602 = t545 * t620;
t600 = g(1) * t479 + g(2) * t477;
t436 = -qJ(5) * t550 + t453;
t595 = pkin(4) * t543 - t677;
t592 = t546 * t352 + t549 * t355;
t372 = pkin(4) * t649 + qJD(5) - t379;
t364 = pkin(5) * t649 - pkin(9) * t494 + t372;
t373 = -qJ(5) * t649 + t380;
t365 = pkin(9) * t492 + t373;
t353 = t364 * t549 - t365 * t546;
t354 = t364 * t546 + t365 * t549;
t540 = t550 * pkin(4);
t409 = pkin(5) * t550 + t527 + t540 + (-t505 - t689) * t545;
t414 = t543 * t689 + t436;
t591 = t409 * t549 - t414 * t546;
t590 = t409 * t546 + t414 * t549;
t484 = t547 * t683 + t550 * t668;
t422 = t484 * t543 + t545 * t667;
t423 = t484 * t545 - t543 * t667;
t369 = t422 * t549 - t423 * t546;
t370 = t422 * t546 + t423 * t549;
t498 = -t545 * t546 + t670;
t589 = qJ(4) * t449 + qJD(4) * t494;
t588 = pkin(8) + t595;
t587 = qJD(3) * pkin(3) - qJD(4) - t567;
t573 = -pkin(8) + t583;
t586 = t573 * t646 + t603 + t621;
t509 = t685 * t545;
t579 = qJD(4) * t543 - qJD(6) * t509 - (-t547 * t693 - t635) * qJD(2) + t395;
t508 = t685 * t543;
t578 = pkin(9) * t624 - qJD(6) * t508 - t702;
t467 = t497 * t547;
t574 = t494 * t641 - t631;
t397 = -t476 * t669 - t477 * t545;
t399 = -t478 * t669 - t479 * t545;
t450 = -t545 * t668 + t608;
t572 = g(1) * t399 + g(2) * t397 + g(3) * t450;
t398 = -t476 * t665 + t477 * t543;
t400 = -t478 * t665 + t479 * t543;
t571 = -g(1) * t400 - g(2) * t398 - g(3) * t451;
t425 = t477 * t550 - t547 * t615;
t427 = t479 * t550 + t547 * t614;
t569 = g(1) * t427 + g(2) * t425 + g(3) * t484;
t566 = -g(3) * t667 + t601;
t359 = -t555 + t691;
t565 = -t359 + t570;
t564 = t563 + t570;
t562 = qJ(5) * t494 + t587;
t561 = -t448 * t680 - t492 * t645 - t569;
t558 = t570 + t606;
t501 = -t625 - t684;
t556 = -pkin(8) * qJDD(3) + (t501 + t625 - t684) * qJD(3);
t554 = -t558 - t697;
t553 = qJD(2) ^ 2;
t502 = -t616 - t690;
t496 = -qJDD(6) - t575;
t486 = t545 * t693 + t616;
t473 = t483 * pkin(3);
t469 = t492 * t649;
t466 = t546 * t666 - t547 * t670;
t465 = t588 * t547;
t441 = -t452 + t540;
t432 = t573 * t547;
t429 = qJD(3) * t484 + t547 * t622;
t428 = -qJD(3) * t483 + t550 * t622;
t417 = t426 * pkin(3);
t416 = t424 * pkin(3);
t415 = t588 * t646 - t621;
t406 = t595 * t649 + t446;
t403 = qJD(3) * t632 + qJD(6) * t467 - t646 * t670;
t402 = qJD(6) * t498 * t547 + qJD(3) * t577;
t394 = t428 * t545 + t543 * t623;
t393 = t428 * t543 - t545 * t623;
t392 = -pkin(4) * t650 - t395;
t385 = t427 * t545 + t478 * t543;
t384 = t427 * t543 - t478 * t545;
t383 = t425 * t545 + t476 * t543;
t382 = t425 * t543 - t476 * t545;
t377 = pkin(4) * t492 - t562;
t366 = -t492 * t693 + t562;
t1 = [t663 * MDP(1) + (-qJD(3) * t429 - qJDD(3) * t483) * MDP(10) + (-qJD(3) * t428 - qJDD(3) * t484) * MDP(11) + (-t360 * t422 + t361 * t423 - t379 * t393 + t380 * t394 - t429 * t587 - t483 * t563 - g(3)) * MDP(15) + (t357 * t423 + t358 * t422 + t359 * t483 + t372 * t393 + t373 * t394 + t377 * t429 - g(3)) * MDP(19) + ((-qJD(6) * t370 + t393 * t549 - t394 * t546) * t531 - t369 * t496 - t429 * t411 - t483 * t363) * MDP(25) + (-(qJD(6) * t369 + t393 * t546 + t394 * t549) * t531 + t370 * t496 - t429 * t413 + t483 * t574) * MDP(26) + (t699 * (t393 * t550 - t422 * t647) - t698 * (-t394 * t550 + t423 * t647)) * qJD(2) + ((-qJDD(2) * MDP(4) + (-MDP(10) * t550 + MDP(11) * t547 - MDP(3)) * t553) * t548 + (MDP(10) * t575 - MDP(11) * t576 + qJDD(2) * MDP(3) - t553 * MDP(4)) * t551) * t544 + t699 * (t422 * t539 + t429 * t492 + t483 * t448) + t703 * (t393 * t494 - t394 * t492 + t422 * t449 - t423 * t448) + t698 * (t423 * t539 + t429 * t494 + t449 * t483); (t357 * t436 + t359 * t465 + t358 * t441 - g(1) * (pkin(4) * t400 + pkin(8) * t479 + qJ(5) * t399 + t629) - g(2) * (pkin(4) * t398 + pkin(8) * t477 + qJ(5) * t397 + t630) - g(3) * (pkin(4) * t451 + qJ(5) * t450 + t607) + (t415 - t603) * t377 + t660 * t373 + t659 * t372) * MDP(19) + (t402 * t413 - t467 * t574) * MDP(20) + (-t363 * t467 - t402 * t411 - t403 * t413 + t466 * t574) * MDP(21) + (t402 * t531 - t413 * t647 - t467 * t496 - t550 * t574) * MDP(22) + (-t363 * t550 - t403 * t531 + t411 * t647 + t466 * t496) * MDP(23) + (-t496 * t550 - t531 * t647) * MDP(24) + (-t591 * t496 + t612 * t550 - t353 * t647 + t432 * t363 + t356 * t466 + t366 * t403 - g(1) * (t399 * t546 + t400 * t549) - g(2) * (t397 * t546 + t398 * t549) - g(3) * (t450 * t546 + t451 * t549) + (t546 * t661 + t549 * t662) * t531 + t586 * t411 + (-t354 * t550 - t531 * t590) * qJD(6)) * MDP(25) + (t590 * t496 - t592 * t550 + t354 * t647 - t432 * t574 + t356 * t467 + t366 * t402 - g(1) * (t399 * t549 - t400 * t546) - g(2) * (t397 * t549 - t398 * t546) - g(3) * (t450 * t549 - t451 * t546) + (-t546 * t662 + t549 * t661) * t531 + t586 * t413 + (-t353 * t550 - t531 * t591) * qJD(6)) * MDP(26) + qJDD(2) * MDP(2) + (t663 * t667 + t601) * MDP(3) + (-t663 * t668 + t600) * MDP(4) + (qJDD(2) * t541 + 0.2e1 * t547 * t617) * MDP(5) + 0.2e1 * (t539 * t547 - t637 * t653) * MDP(6) + (qJDD(3) * t547 + t550 * t552) * MDP(7) + (qJDD(3) * t550 - t547 * t552) * MDP(8) + (t556 * t547 + t550 * t695) * MDP(10) + (-t547 * t695 + t556 * t550) * MDP(11) + ((-t605 + pkin(8) * t448 - t563 * t543 + (qJD(2) * t452 + t379) * qJD(3)) * t547 + (-qJDD(2) * t452 - t360 + (pkin(8) * t492 - t543 * t587) * qJD(3) - t658 * qJD(2)) * t550 + t571) * MDP(12) + ((-t604 + pkin(8) * t449 - t563 * t545 + (-qJD(2) * t453 - t380) * qJD(3)) * t547 + (qJDD(2) * t453 + t361 + (pkin(8) * t494 - t545 * t587) * qJD(3) + t657 * qJD(2)) * t550 + t572) * MDP(13) + (-t448 * t453 - t449 * t452 - t658 * t494 - t657 * t492 + (-t379 * t545 - t380 * t543) * t646 + (-t360 * t545 - t361 * t543 + t566) * t547) * MDP(14) + (t361 * t453 + t360 * t452 + t587 * t603 - g(1) * t629 - g(2) * t630 - g(3) * t607 + t657 * t380 + t658 * t379 + (-t547 * t563 - t587 * t646 - t600) * pkin(8)) * MDP(15) + (t415 * t492 + t448 * t465 + (-t605 + t359 * t543 + (-qJD(2) * t441 - t372) * qJD(3)) * t547 + (qJD(2) * t659 + qJDD(2) * t441 + t377 * t648 + t358) * t550 + t571) * MDP(16) + (-t436 * t448 + t441 * t449 + t659 * t494 - t660 * t492 + (t372 * t545 - t373 * t543) * t646 + (-t357 * t543 + t358 * t545 + t566) * t547) * MDP(17) + (-t415 * t494 - t449 * t465 + (t604 - t359 * t545 + (qJD(2) * t436 + t373) * qJD(3)) * t547 + (-qJD(2) * t660 - qJDD(2) * t436 - t377 * t639 - t357) * t550 - t572) * MDP(18); MDP(7) * t636 + MDP(8) * t539 + qJDD(3) * MDP(9) + (qJD(3) * t446 - t501 * t650 + t558) * MDP(10) + (-t599 + (-qJD(2) * t501 - t457) * t550 + t569) * MDP(11) + (-pkin(3) * t448 - t446 * t492 + t564 * t545 + (-t379 * t547 + t395 * t550 + (t550 * t587 - t627) * t543) * qJD(2) + t654) * MDP(12) + (t602 - pkin(3) * t449 - t446 * t494 - t564 * t543 + (t380 * t547 - t396 * t550 + (-t627 + (qJD(4) + t587) * t550) * t545) * qJD(2)) * MDP(13) + (t395 * t494 + t396 * t492 + (t379 * t649 + t361) * t545 + (t380 * t649 - t360 + t589) * t543 + t561) * MDP(14) + (t563 * pkin(3) + g(1) * t417 + g(2) * t416 + g(3) * t473 - t379 * t395 - t380 * t396 + t587 * t446 + (-t379 * t543 + t380 * t545) * qJD(4) + (-t360 * t543 + t361 * t545 - t569) * qJ(4)) * MDP(15) + (t448 * t502 + (-t406 - t643) * t492 + t565 * t545 + (t372 * t547 - t392 * t550 + (-t377 * t550 - t627) * t543) * qJD(2) + t654) * MDP(16) + (t389 * t492 - t392 * t494 + (-t372 * t649 + t357) * t545 + (t373 * t649 + t358 + t589) * t543 + t561) * MDP(17) + (-t602 + t406 * t494 - t449 * t502 + (t565 + t644) * t543 + (-t373 * t547 + t389 * t550 + (t627 + (-qJD(4) + t377) * t550) * t545) * qJD(2)) * MDP(18) + (t357 * t680 + t359 * t502 - t377 * t406 - t372 * t392 - g(1) * (qJ(4) * t427 - t426 * t690 - t417) - g(2) * (qJ(4) * t425 - t424 * t690 - t416) - g(3) * (qJ(4) * t484 - t483 * t690 - t473) + t702 * t373 + (t358 * qJ(4) + qJ(5) * t570 + t372 * qJD(4) - t377 * qJD(5)) * t543) * MDP(19) + (t413 * t656 - t498 * t574) * MDP(20) + (-t363 * t498 - t411 * t656 - t413 * t655 + t497 * t574) * MDP(21) + (t413 * t650 - t496 * t498 + t531 * t656) * MDP(22) + (-t411 * t650 + t496 * t497 - t531 * t655) * MDP(23) + t531 * MDP(24) * t650 + (-(t508 * t549 - t509 * t546) * t496 + t486 * t363 + t353 * t650 + (t546 * t578 + t549 * t579) * t531 + t609 * t411 + t655 * t366 + t694 * t497) * MDP(25) + ((t508 * t546 + t509 * t549) * t496 - t486 * t574 - t354 * t650 + (-t546 * t579 + t549 * t578) * t531 + t609 * t413 + t656 * t366 + t694 * t498) * MDP(26) + (-MDP(5) * t547 * t550 + MDP(6) * t653) * t553; (t379 * t494 + t380 * t492 + t554) * MDP(15) + (t691 - t678 + t373 * t492 + (-qJD(5) - t372) * t494 + t554) * MDP(19) + (-t363 - t673) * MDP(25) + (t574 + t674) * MDP(26) + t699 * (t543 * t636 - t537 + (-t494 + t648) * t649) + t698 * (t469 + t449) + t703 * (-t492 ^ 2 - t490); (t492 * t494 + t575) * MDP(16) + (-t469 + t449) * MDP(17) + (-t542 * t553 - t490) * MDP(18) + (-g(1) * t384 - g(2) * t382 - g(3) * t422 + t377 * t494 + (-pkin(4) * t647 + t373 * t550) * qJD(2) + t584) * MDP(19) + (-t411 * t494 - t549 * t496) * MDP(25) + (-t413 * t494 + t546 * t496) * MDP(26) + (-MDP(25) * t546 - MDP(26) * t549) * t531 ^ 2; t413 * t411 * MDP(20) + (-t411 ^ 2 + t413 ^ 2) * MDP(21) + (t631 + t674) * MDP(22) + (-t611 + t673) * MDP(23) - t496 * MDP(24) + (t354 * t531 - t366 * t413 - g(1) * (t384 * t549 - t385 * t546) - g(2) * (t382 * t549 - t383 * t546) - g(3) * t369 + t612) * MDP(25) + (t353 * t531 + t366 * t411 - g(1) * (-t384 * t546 - t385 * t549) - g(2) * (-t382 * t546 - t383 * t549) + g(3) * t370 - t592) * MDP(26) + (-MDP(22) * t671 - MDP(23) * t413 - MDP(25) * t354 - MDP(26) * t353) * qJD(6);];
tau  = t1;
