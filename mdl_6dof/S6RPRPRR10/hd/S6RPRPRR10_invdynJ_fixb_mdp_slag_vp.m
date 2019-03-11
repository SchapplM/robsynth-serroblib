% Calculate vector of inverse dynamics joint torques for
% S6RPRPRR10
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
%   pkin=[a2,a3,a4,a5,a6,d1,d3,d5,d6,theta4]';
% MDP [31x1]
%   Minimal dynamic parameter vector (fixed base model)
%   see S6RPRPRR10_convert_par2_MPV_fixb.m
% 
% Output:
% tau [6x1]
%   joint torques of inverse dynamics (contains inertial, gravitational Coriolis and centrifugal forces)

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 04:10
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tau = S6RPRPRR10_invdynJ_fixb_mdp_slag_vp(qJ, qJD, qJDD, g, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(3,1),zeros(10,1),zeros(31,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRPRR10_invdynJ_fixb_mdp_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RPRPRR10_invdynJ_fixb_mdp_slag_vp: qJD has to be [6x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [6 1]), ...
  'S6RPRPRR10_invdynJ_fixb_mdp_slag_vp: qJDD has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RPRPRR10_invdynJ_fixb_mdp_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RPRPRR10_invdynJ_fixb_mdp_slag_vp: pkin has to be [10x1] (double)');
assert(isreal(MDP) && all(size(MDP) == [31 1]), ...
  'S6RPRPRR10_invdynJ_fixb_mdp_slag_vp: MDP has to be [31x1] (double)'); 

%% Symbolic Calculation
% From invdyn_joint_fixb_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 04:09:36
% EndTime: 2019-03-09 04:09:50
% DurationCPUTime: 10.24s
% Computational Cost: add. (5275->580), mult. (10814->759), div. (0->0), fcn. (7909->14), ass. (0->234)
t558 = sin(qJ(3));
t642 = qJD(1) * t558;
t536 = qJD(5) + t642;
t529 = qJD(6) + t536;
t554 = sin(pkin(10));
t562 = cos(qJ(3));
t641 = qJD(1) * t562;
t617 = t554 * t641;
t555 = cos(pkin(10));
t630 = t555 * qJD(3);
t505 = t617 - t630;
t616 = t555 * t641;
t640 = qJD(3) * t554;
t507 = t616 + t640;
t557 = sin(qJ(5));
t561 = cos(qJ(5));
t437 = t505 * t561 + t507 * t557;
t560 = cos(qJ(6));
t436 = t505 * t557 - t507 * t561;
t556 = sin(qJ(6));
t669 = t436 * t556;
t694 = -t437 * t560 + t669;
t692 = t529 * t694;
t592 = pkin(3) * t562 + qJ(4) * t558;
t516 = t592 * qJD(1);
t564 = -pkin(1) - pkin(7);
t530 = qJD(1) * t564 + qJD(2);
t662 = t555 * t562;
t451 = t516 * t554 + t530 * t662;
t618 = t554 * t642;
t426 = pkin(8) * t618 + t451;
t695 = -qJD(4) * t555 + t426;
t681 = -t554 * t557 + t555 * t561;
t693 = qJD(5) * t681;
t691 = t436 * t536;
t690 = t437 * t536;
t585 = t436 * t560 + t437 * t556;
t689 = t529 * t585;
t684 = t558 * t681;
t651 = qJD(1) * t684 + t693;
t514 = t554 * t561 + t555 * t557;
t492 = t514 * qJD(1);
t680 = t514 * qJD(5);
t650 = t492 * t558 + t680;
t672 = qJ(4) * t562;
t591 = pkin(3) * t558 - t672;
t518 = qJ(2) + t591;
t497 = t518 * qJD(1);
t517 = t558 * t530;
t498 = qJD(3) * qJ(4) + t517;
t427 = t497 * t555 - t498 * t554;
t398 = pkin(4) * t642 - pkin(8) * t507 + t427;
t428 = t497 * t554 + t498 * t555;
t400 = -pkin(8) * t505 + t428;
t368 = t398 * t557 + t400 * t561;
t363 = -pkin(9) * t437 + t368;
t632 = qJD(6) * t556;
t361 = t363 * t632;
t665 = t530 * t562;
t489 = -qJD(3) * pkin(3) + qJD(4) - t665;
t448 = pkin(4) * t505 + t489;
t392 = pkin(5) * t437 + t448;
t551 = pkin(10) + qJ(5);
t546 = qJ(6) + t551;
t538 = sin(t546);
t539 = cos(t546);
t563 = cos(qJ(1));
t559 = sin(qJ(1));
t661 = t558 * t559;
t462 = t538 * t563 + t539 * t661;
t660 = t558 * t563;
t464 = -t538 * t559 + t539 * t660;
t676 = g(3) * t562;
t687 = g(1) * t462 - g(2) * t464 - t392 * t694 + t539 * t676 + t361;
t461 = -t538 * t661 + t539 * t563;
t463 = t538 * t660 + t539 * t559;
t623 = qJDD(1) * t562;
t600 = -qJDD(3) * t555 + t554 * t623;
t628 = qJD(1) * qJD(3);
t613 = t558 * t628;
t459 = t554 * t613 - t600;
t647 = qJDD(3) * t554 + t555 * t623;
t460 = t555 * t613 - t647;
t634 = qJD(5) * t561;
t635 = qJD(5) * t557;
t381 = t459 * t557 - t460 * t561 - t505 * t634 - t507 * t635;
t612 = t562 * t628;
t624 = qJDD(1) * t558;
t578 = t612 + t624;
t512 = qJDD(5) + t578;
t484 = qJD(3) * t592 - qJD(4) * t562 + qJD(2);
t424 = qJD(1) * t484 + qJDD(1) * t518;
t528 = qJDD(1) * t564 + qJDD(2);
t447 = qJDD(3) * qJ(4) + t528 * t558 + (qJD(4) + t665) * qJD(3);
t383 = t424 * t555 - t447 * t554;
t372 = pkin(4) * t578 + pkin(8) * t460 + t383;
t384 = t424 * t554 + t447 * t555;
t375 = pkin(8) * t459 + t384;
t607 = t372 * t561 - t375 * t557;
t568 = -qJD(5) * t368 + t607;
t352 = pkin(5) * t512 - pkin(9) * t381 + t568;
t382 = -qJD(5) * t436 - t459 * t561 - t460 * t557;
t577 = t372 * t557 + t375 * t561 + t398 * t634 - t400 * t635;
t353 = -pkin(9) * t382 + t577;
t608 = t352 * t560 - t556 * t353;
t686 = -g(1) * t461 - g(2) * t463 + t392 * t585 + t538 * t676 + t608;
t499 = qJDD(6) + t512;
t685 = t499 * MDP(29) + (t585 ^ 2 - t694 ^ 2) * MDP(26) + t694 * MDP(25) * t585;
t682 = g(1) * t559 - g(2) * t563;
t683 = t682 * t555;
t503 = t555 * t518;
t611 = -t554 * t564 + pkin(4);
t435 = -pkin(8) * t662 + t558 * t611 + t503;
t659 = t558 * t564;
t458 = t518 * t554 + t555 * t659;
t663 = t554 * t562;
t449 = -pkin(8) * t663 + t458;
t652 = t435 * t557 + t449 * t561;
t675 = pkin(8) + qJ(4);
t524 = t675 * t554;
t525 = t675 * t555;
t648 = -t524 * t557 + t525 * t561;
t677 = g(3) * t558;
t573 = t562 * t682 - t677;
t606 = t381 * t556 + t382 * t560;
t357 = -qJD(6) * t585 + t606;
t450 = t516 * t555 - t530 * t663;
t621 = pkin(8) * t555 * t558;
t409 = (pkin(4) * t562 + t621) * qJD(1) + t450;
t584 = qJD(4) * t554 + qJD(5) * t525;
t679 = t524 * t634 + t695 * t561 + (t409 + t584) * t557;
t674 = pkin(1) * qJDD(1);
t566 = qJD(1) ^ 2;
t673 = qJ(2) * t566;
t367 = t398 * t561 - t400 * t557;
t362 = pkin(9) * t436 + t367;
t360 = pkin(5) * t536 + t362;
t671 = t360 * t560;
t670 = t363 * t560;
t639 = qJD(3) * t558;
t590 = -qJDD(3) * pkin(3) + t530 * t639 + qJDD(4);
t452 = -t528 * t562 + t590;
t668 = t452 * t554;
t667 = t452 * t555;
t666 = t452 * t562;
t442 = t514 * t556 - t560 * t681;
t657 = -qJD(6) * t442 - t556 * t650 + t560 * t651;
t443 = t514 * t560 + t556 * t681;
t656 = qJD(6) * t443 + t556 * t651 + t560 * t650;
t483 = t681 * t562;
t654 = -qJD(3) * t483 + t558 * t680 + t492;
t638 = qJD(3) * t562;
t653 = -qJD(1) * t681 - qJD(5) * t684 - t514 * t638;
t637 = qJD(3) * t564;
t614 = t562 * t637;
t446 = t484 * t554 + t555 * t614;
t646 = pkin(1) * t563 + qJ(2) * t559;
t552 = t558 ^ 2;
t553 = t562 ^ 2;
t645 = t552 - t553;
t565 = qJD(3) ^ 2;
t644 = -t565 - t566;
t631 = qJD(6) * t560;
t629 = -qJD(4) + t489;
t627 = qJDD(1) * qJ(2);
t626 = qJDD(1) * t554;
t625 = qJDD(1) * t555;
t622 = qJDD(3) * t558;
t620 = 0.2e1 * qJD(1) * qJD(2);
t619 = t381 * t560 - t382 * t556 - t437 * t631;
t540 = -pkin(4) * t555 - pkin(3);
t615 = t554 * t639;
t471 = -pkin(4) * t618 + t517;
t610 = pkin(5) * t650 - t471;
t466 = t555 * t484;
t407 = t466 + (t562 * t611 + t621) * qJD(3);
t423 = pkin(8) * t615 + t446;
t605 = t407 * t561 - t423 * t557;
t604 = t435 * t561 - t449 * t557;
t602 = -t524 * t561 - t525 * t557;
t504 = pkin(4) * t663 - t562 * t564;
t601 = qJD(6) * t360 + t353;
t599 = g(1) * t563 + g(2) * t559;
t597 = qJDD(2) - t682;
t406 = t561 * t409;
t418 = pkin(9) * t681 + t648;
t596 = pkin(5) * t641 + pkin(9) * t651 + qJD(4) * t514 + qJD(5) * t648 + qJD(6) * t418 - t426 * t557 + t406;
t417 = -pkin(9) * t514 + t602;
t595 = pkin(9) * t650 - qJD(6) * t417 + t679;
t480 = t514 * t558;
t594 = qJD(6) * t480 + t654;
t593 = qJD(6) * t684 - t653;
t355 = t360 * t556 + t670;
t589 = -t383 * t554 + t384 * t555;
t588 = t427 * t555 + t428 * t554;
t587 = -t427 * t554 + t428 * t555;
t481 = t514 * t562;
t411 = t481 * t560 + t483 * t556;
t412 = -t481 * t556 + t483 * t560;
t490 = -pkin(4) * t615 + t558 * t637;
t583 = -t528 + t682;
t582 = t682 * t554;
t579 = 0.2e1 * qJ(2) * t628 + qJDD(3) * t564;
t576 = t407 * t557 + t423 * t561 + t435 * t634 - t449 * t635;
t356 = t436 * t632 + t619;
t574 = t583 + t673;
t571 = -t558 * t682 - t676;
t399 = -pkin(4) * t459 + t452;
t569 = -t599 + t620 + 0.2e1 * t627;
t567 = -t564 * t565 + t569;
t548 = t563 * qJ(2);
t545 = qJDD(3) * t562;
t544 = cos(t551);
t543 = sin(t551);
t475 = -t543 * t559 + t544 * t660;
t474 = t543 * t660 + t544 * t559;
t473 = t543 * t563 + t544 * t661;
t472 = -t543 * t661 + t544 * t563;
t470 = -pkin(5) * t681 + t540;
t457 = -t554 * t659 + t503;
t445 = -t554 * t614 + t466;
t444 = pkin(5) * t481 + t504;
t422 = -t557 * t558 * t630 - t561 * t615 + t562 * t693;
t420 = -qJD(3) * t684 - t562 * t680;
t397 = pkin(5) * t422 + t490;
t377 = -pkin(9) * t481 + t652;
t376 = pkin(5) * t558 - pkin(9) * t483 + t604;
t366 = qJD(6) * t412 + t420 * t556 + t422 * t560;
t365 = -qJD(6) * t411 + t420 * t560 - t422 * t556;
t364 = pkin(5) * t382 + t399;
t359 = -pkin(9) * t422 + t576;
t358 = pkin(5) * t638 - pkin(9) * t420 - qJD(5) * t652 + t605;
t354 = -t363 * t556 + t671;
t1 = [(t683 + (t667 + t564 * t460 + (-qJD(1) * t458 - t428) * qJD(3)) * t562 + (-t446 * qJD(1) - t458 * qJDD(1) - t384 + t599 * t554 + (-t489 * t555 + t507 * t564) * qJD(3)) * t558) * MDP(15) + (t356 * t412 - t365 * t585) * MDP(25) + (-t355 * t638 + g(1) * t463 - g(2) * t461 + t444 * t356 + t361 * t558 + t364 * t412 + t392 * t365 - t397 * t585 + (-(-qJD(6) * t377 + t358) * t529 - t376 * t499 - t352 * t558) * t556 + (-(qJD(6) * t376 + t359) * t529 - t377 * t499 - t601 * t558) * t560) * MDP(31) + (t356 * t558 + t365 * t529 + t412 * t499 - t585 * t638) * MDP(27) + t682 * MDP(2) + (-t558 * t565 + t545) * MDP(9) + (t582 + (t668 + t564 * t459 + (qJD(1) * t457 + t427) * qJD(3)) * t562 + (t445 * qJD(1) + t457 * qJDD(1) + t383 - t599 * t555 + (-t489 * t554 + t505 * t564) * qJD(3)) * t558) * MDP(14) + t599 * MDP(3) + (-t382 * t558 - t422 * t536 - t437 * t638 - t481 * t512) * MDP(21) + (t512 * t558 + t536 * t638) * MDP(22) + (t499 * t558 + t529 * t638) * MDP(29) + (-t445 * t507 - t446 * t505 + t457 * t460 + t458 * t459 + t588 * t639 + (-t383 * t555 - t384 * t554 + t599) * t562) * MDP(16) + (t558 * t567 + t562 * t579) * MDP(12) + (-t558 * t579 + t562 * t567) * MDP(13) + (t384 * t458 + t428 * t446 + t383 * t457 + t427 * t445 - g(1) * (pkin(3) * t660 - t563 * t672 + t548) - g(2) * (pkin(7) * t563 + t646) + (t489 * t639 - t666) * t564 + (-g(1) * t564 - g(2) * t591) * t559) * MDP(17) + (-(qJDD(2) - t674) * pkin(1) - g(1) * (-pkin(1) * t559 + t548) - g(2) * t646 + (t620 + t627) * qJ(2)) * MDP(6) + (t597 - 0.2e1 * t674) * MDP(4) + t569 * MDP(5) + (t605 * t536 + t604 * t512 + t607 * t558 + t367 * t638 + t490 * t437 + t504 * t382 + t399 * t481 + t448 * t422 - g(1) * t475 - g(2) * t473 + (-t368 * t558 - t536 * t652) * qJD(5)) * MDP(23) + (-t381 * t481 - t382 * t483 - t420 * t437 + t422 * t436) * MDP(19) + (t381 * t483 - t420 * t436) * MDP(18) + (t381 * t558 + t420 * t536 - t436 * t638 + t483 * t512) * MDP(20) + (g(1) * t474 - g(2) * t472 - t368 * t638 + t504 * t381 + t399 * t483 + t448 * t420 - t436 * t490 - t512 * t652 - t536 * t576 - t558 * t577) * MDP(24) + 0.2e1 * (-t558 * t623 + t628 * t645) * MDP(8) + qJDD(1) * MDP(1) + (-t562 * t565 - t622) * MDP(10) + (qJDD(1) * t553 - 0.2e1 * t558 * t612) * MDP(7) + (-t356 * t411 - t357 * t412 + t365 * t694 + t366 * t585) * MDP(26) + (-t357 * t558 - t366 * t529 - t411 * t499 + t638 * t694) * MDP(28) + ((t358 * t560 - t359 * t556) * t529 + (t376 * t560 - t377 * t556) * t499 + t608 * t558 + t354 * t638 - t397 * t694 + t444 * t357 + t364 * t411 + t392 * t366 - g(1) * t464 - g(2) * t462 + ((-t376 * t556 - t377 * t560) * t529 - t355 * t558) * qJD(6)) * MDP(30); qJDD(1) * MDP(4) - t566 * MDP(5) + (t597 - t673 - t674) * MDP(6) + (t558 * t644 + t545) * MDP(12) + (t562 * t644 - t622) * MDP(13) + (-t552 * t626 + t459 * t562 + (-t555 * t566 + (t505 - 0.2e1 * t617) * qJD(3)) * t558) * MDP(14) + (-t552 * t625 + t460 * t562 + (t554 * t566 + (t507 - 0.2e1 * t616) * qJD(3)) * t558) * MDP(15) + ((qJD(1) * t507 + t459 * t558 - t505 * t638) * t555 + (qJD(1) * t505 - t460 * t558 + t507 * t638) * t554) * MDP(16) + (-t666 + t589 * t558 - t588 * qJD(1) + (t489 * t558 + t562 * t587) * qJD(3) - t682) * MDP(17) + (-t382 * t562 + t437 * t639 - t480 * t512 + t536 * t653) * MDP(23) + (-t381 * t562 - t436 * t639 - t512 * t684 + t536 * t654) * MDP(24) + ((-t480 * t560 - t556 * t684) * t499 - t694 * t639 - t562 * t357 + (t556 * t594 - t560 * t593) * t529) * MDP(30) + (-(-t480 * t556 + t560 * t684) * t499 - t585 * t639 - t562 * t356 + (t556 * t593 + t560 * t594) * t529) * MDP(31); MDP(9) * t623 + (pkin(3) * t460 + t668 + (t582 + (-qJ(4) * t630 + t428) * qJD(1)) * t562 + (-qJ(4) * t625 - g(3) * t554 - t507 * t530 + (t555 * t629 + t451) * qJD(1)) * t558) * MDP(15) + (pkin(3) * t459 - t667 + (-t683 + (-qJ(4) * t640 - t427) * qJD(1)) * t562 + (-qJ(4) * t626 + g(3) * t555 - t505 * t530 + (t554 * t629 - t450) * qJD(1)) * t558) * MDP(14) + (-t489 * t517 - t427 * t450 - t428 * t451 + t587 * qJD(4) + (-t452 - t573) * pkin(3) + (t571 + t589) * qJ(4)) * MDP(17) + (t450 * t507 + t451 * t505 + (qJ(4) * t459 - qJD(4) * t505 - t427 * t642 + t384) * t555 + (-qJ(4) * t460 + qJD(4) * t507 - t428 * t642 - t383) * t554 + t571) * MDP(16) + (t558 * t574 + t676) * MDP(13) + (-t562 * t574 + t677) * MDP(12) + qJDD(3) * MDP(11) + (t540 * t381 + t399 * t514 + t436 * t471 + t651 * t448 - t648 * t512 + t536 * t679 + t573 * t543) * MDP(24) + ((t417 * t560 - t418 * t556) * t499 + t470 * t357 + t364 * t442 + (t556 * t595 - t560 * t596) * t529 + t656 * t392 - t610 * t694 - t573 * t539) * MDP(30) + (-t442 * t499 - t529 * t656) * MDP(28) + (t443 * t499 + t529 * t657) * MDP(27) + (-(t417 * t556 + t418 * t560) * t499 + t470 * t356 + t364 * t443 + (t556 * t596 + t560 * t595) * t529 + t657 * t392 - t610 * t585 + t573 * t538) * MDP(31) + (t356 * t443 - t585 * t657) * MDP(25) + (-t356 * t442 - t357 * t443 + t585 * t656 + t657 * t694) * MDP(26) + (t381 * t514 - t436 * t651) * MDP(18) + (t512 * t514 + t536 * t651) * MDP(20) + (t381 * t681 - t382 * t514 + t436 * t650 - t437 * t651) * MDP(19) + (t602 * t512 + t540 * t382 - t399 * t681 - t471 * t437 + (-t406 - t584 * t561 + (qJD(5) * t524 + t695) * t557) * t536 + t650 * t448 - t573 * t544) * MDP(23) + (t512 * t681 - t536 * t650) * MDP(21) - MDP(10) * t624 + (MDP(20) * t436 + t437 * MDP(21) - t536 * MDP(22) - t367 * MDP(23) + t368 * MDP(24) + MDP(27) * t585 - MDP(28) * t694 - t529 * MDP(29) - t354 * MDP(30) + t355 * MDP(31)) * t641 + (MDP(7) * t558 * t562 - MDP(8) * t645) * t566; ((t507 - t640) * t642 + t600) * MDP(14) + ((-t505 - t630) * t642 + t647) * MDP(15) + (-t505 ^ 2 - t507 ^ 2) * MDP(16) + (t427 * t507 + t428 * t505 + t562 * t583 + t590 - t677) * MDP(17) + (t382 - t691) * MDP(23) + (t381 - t690) * MDP(24) + (t357 - t689) * MDP(30) + (t356 + t692) * MDP(31); -t436 * t437 * MDP(18) + (t436 ^ 2 - t437 ^ 2) * MDP(19) + (t381 + t690) * MDP(20) + (-t382 - t691) * MDP(21) + t512 * MDP(22) + (-g(1) * t472 - g(2) * t474 + t368 * t536 + t436 * t448 + t543 * t676 + t568) * MDP(23) + (g(1) * t473 - g(2) * t475 + t367 * t536 + t437 * t448 + t544 * t676 - t577) * MDP(24) + (t356 - t692) * MDP(27) + (-t357 - t689) * MDP(28) + (-(-t362 * t556 - t670) * t529 - t355 * qJD(6) + (-t436 * t694 + t499 * t560 - t529 * t632) * pkin(5) + t686) * MDP(30) + ((-t363 * t529 - t352) * t556 + (t362 * t529 - t601) * t560 + (-t436 * t585 - t499 * t556 - t529 * t631) * pkin(5) + t687) * MDP(31) + t685; (t619 - t692) * MDP(27) + (-t606 - t689) * MDP(28) + (t355 * t529 + t686) * MDP(30) + (-t556 * t352 - t560 * t353 + t354 * t529 + t687) * MDP(31) + (MDP(27) * t669 + MDP(28) * t585 - MDP(30) * t355 - MDP(31) * t671) * qJD(6) + t685;];
tau  = t1;
