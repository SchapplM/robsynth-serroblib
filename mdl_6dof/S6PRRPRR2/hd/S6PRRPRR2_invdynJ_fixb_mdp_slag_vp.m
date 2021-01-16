% Calculate vector of inverse dynamics joint torques for
% S6PRRPRR2
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
%   see S6PRRPRR2_convert_par2_MPV_fixb.m
% 
% Output:
% tau [6x1]
%   joint torques of inverse dynamics (contains inertial, gravitational Coriolis and centrifugal forces)

% Quelle: HybrDyn-Toolbox
% Datum: 2021-01-16 03:49
% Revision: 24b2e7d74a0c1a3b64fa2f8f5ad758691ad61af3 (2021-01-15)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tau = S6PRRPRR2_invdynJ_fixb_mdp_slag_vp(qJ, qJD, qJDD, g, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(3,1),zeros(12,1),zeros(29,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRRPRR2_invdynJ_fixb_mdp_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6PRRPRR2_invdynJ_fixb_mdp_slag_vp: qJD has to be [6x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [6 1]), ...
  'S6PRRPRR2_invdynJ_fixb_mdp_slag_vp: qJDD has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6PRRPRR2_invdynJ_fixb_mdp_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [12 1]), ...
  'S6PRRPRR2_invdynJ_fixb_mdp_slag_vp: pkin has to be [12x1] (double)');
assert(isreal(MDP) && all(size(MDP) == [29 1]), ...
  'S6PRRPRR2_invdynJ_fixb_mdp_slag_vp: MDP has to be [29x1] (double)'); 

%% Symbolic Calculation
% From invdyn_joint_fixb_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2021-01-16 03:47:00
% EndTime: 2021-01-16 03:47:25
% DurationCPUTime: 11.60s
% Computational Cost: add. (5403->597), mult. (12723->816), div. (0->0), fcn. (10376->18), ass. (0->266)
t562 = cos(qJ(3));
t696 = cos(pkin(12));
t615 = t696 * t562;
t530 = qJD(2) * t615;
t551 = sin(pkin(12));
t558 = sin(qJ(3));
t649 = qJD(2) * t558;
t506 = t551 * t649 - t530;
t707 = qJD(5) + qJD(6);
t717 = t506 + t707;
t555 = qJ(4) + pkin(8);
t622 = qJD(3) * t555;
t501 = qJD(4) * t562 - t558 * t622;
t581 = -qJD(4) * t558 - t562 * t622;
t438 = t501 * t696 + t551 * t581;
t673 = t551 * t558;
t587 = t615 - t673;
t553 = sin(pkin(6));
t563 = cos(qJ(2));
t667 = t553 * t563;
t630 = qJD(1) * t667;
t482 = t587 * t630;
t657 = t438 - t482;
t616 = t696 * t558;
t517 = t551 * t562 + t616;
t508 = t517 * qJD(3);
t511 = t587 * qJD(3);
t559 = sin(qJ(2));
t650 = qJD(1) * t559;
t631 = t553 * t650;
t698 = qJD(3) * pkin(3);
t636 = t558 * t698;
t716 = pkin(4) * t508 - pkin(9) * t511 - t631 + t636;
t509 = t517 * qJD(2);
t557 = sin(qJ(5));
t561 = cos(qJ(5));
t643 = t561 * qJD(3);
t485 = t509 * t557 - t643;
t560 = cos(qJ(6));
t487 = qJD(3) * t557 + t509 * t561;
t556 = sin(qJ(6));
t683 = t487 * t556;
t419 = t560 * t485 + t683;
t498 = qJD(5) + t506;
t496 = qJD(6) + t498;
t715 = t419 * t496;
t596 = t485 * t556 - t560 * t487;
t714 = t496 * t596;
t665 = t556 * t561;
t520 = t557 * t560 + t665;
t654 = t717 * t520;
t647 = qJD(5) * t557;
t682 = t506 * t557;
t713 = t647 + t682;
t609 = qJD(2) * t555 + t631;
t554 = cos(pkin(6));
t651 = qJD(1) * t554;
t477 = -t558 * t609 + t562 * t651;
t471 = t477 + t698;
t478 = t558 * t651 + t562 * t609;
t617 = t696 * t478;
t408 = t551 * t471 + t617;
t403 = qJD(3) * pkin(9) + t408;
t539 = pkin(3) * t562 + pkin(2);
t497 = -qJD(2) * t539 + qJD(4) - t630;
t424 = pkin(4) * t506 - pkin(9) * t509 + t497;
t382 = t403 * t561 + t424 * t557;
t376 = -pkin(10) * t485 + t382;
t645 = qJD(6) * t556;
t374 = t376 * t645;
t467 = t551 * t478;
t407 = t471 * t696 - t467;
t402 = -qJD(3) * pkin(4) - t407;
t389 = t485 * pkin(5) + t402;
t697 = cos(pkin(11));
t619 = t697 * t559;
t552 = sin(pkin(11));
t670 = t552 * t563;
t504 = t554 * t619 + t670;
t547 = qJ(3) + pkin(12);
t540 = sin(t547);
t541 = cos(t547);
t620 = t553 * t697;
t464 = t504 * t541 - t540 * t620;
t618 = t697 * t563;
t671 = t552 * t559;
t502 = t554 * t671 - t618;
t672 = t552 * t553;
t466 = -t502 * t541 + t540 * t672;
t669 = t553 * t559;
t494 = t540 * t554 + t541 * t669;
t503 = -t554 * t618 + t671;
t505 = t554 * t670 + t619;
t550 = qJ(5) + qJ(6);
t545 = sin(t550);
t546 = cos(t550);
t712 = t389 * t419 - g(1) * (-t466 * t546 - t505 * t545) - g(2) * (-t464 * t546 - t503 * t545) - g(3) * (-t494 * t546 + t545 * t667) + t374;
t641 = qJD(2) * qJD(3);
t625 = t558 * t641;
t573 = qJDD(2) * t517 - t551 * t625;
t462 = qJD(3) * t530 + t573;
t404 = qJD(5) * t643 + t557 * qJDD(3) + t561 * t462 - t509 * t647;
t640 = qJDD(2) * t558;
t600 = -qJDD(2) * t615 + t551 * t640;
t461 = qJD(2) * t508 + t600;
t457 = qJDD(5) + t461;
t638 = t554 * qJDD(1);
t529 = t562 * t638;
t642 = qJD(1) * qJD(2);
t489 = qJDD(2) * pkin(8) + (qJDD(1) * t559 + t563 * t642) * t553;
t570 = qJ(4) * qJDD(2) + qJD(2) * qJD(4) + qJD(3) * t651 + t489;
t595 = t609 * qJD(3);
t397 = qJDD(3) * pkin(3) - t558 * t570 - t562 * t595 + t529;
t398 = (-t595 + t638) * t558 + t570 * t562;
t373 = t551 * t397 + t696 * t398;
t371 = qJDD(3) * pkin(9) + t373;
t626 = t559 * t642;
t527 = t553 * t626;
t575 = pkin(3) * t625 - qJDD(2) * t539 + qJDD(4) + t527;
t623 = qJDD(1) * t667;
t460 = t575 - t623;
t387 = pkin(4) * t461 - pkin(9) * t462 + t460;
t386 = t561 * t387;
t569 = -qJD(5) * t382 - t557 * t371 + t386;
t359 = pkin(5) * t457 - pkin(10) * t404 + t569;
t405 = qJD(5) * t487 - t561 * qJDD(3) + t462 * t557;
t646 = qJD(5) * t561;
t586 = -t561 * t371 - t557 * t387 + t403 * t647 - t424 * t646;
t360 = -pkin(10) * t405 - t586;
t613 = t560 * t359 - t556 * t360;
t711 = t389 * t596 - g(1) * (-t466 * t545 + t505 * t546) - g(2) * (-t464 * t545 + t503 * t546) - g(3) * (-t494 * t545 - t546 * t667) + t613;
t453 = qJDD(6) + t457;
t710 = t453 * MDP(27) + (-t419 ^ 2 + t596 ^ 2) * MDP(24) - t419 * MDP(23) * t596;
t448 = t520 * t517;
t709 = t482 * t557 + t716 * t561;
t658 = t501 * t551 - t517 * t630 - t696 * t581;
t459 = -pkin(4) * t587 - pkin(9) * t517 - t539;
t524 = t555 * t562;
t484 = t524 * t696 - t555 * t673;
t708 = -t459 * t646 + t484 * t647 - t716 * t557 - t657 * t561;
t519 = t556 * t557 - t560 * t561;
t655 = t717 * t519;
t706 = -t453 * t520 + t496 * t655;
t564 = qJD(3) ^ 2;
t604 = g(1) * t505 + g(2) * t503;
t705 = 0.2e1 * qJDD(2) * pkin(2) - pkin(8) * t564 + t553 * (-g(3) * t563 + t626) - t527 + t604 + t623;
t612 = t404 * t556 + t560 * t405;
t366 = -qJD(6) * t596 + t612;
t704 = pkin(3) * t551;
t703 = pkin(3) * t558;
t512 = t554 * t562 - t558 * t669;
t702 = g(3) * t512;
t701 = g(3) * t553;
t535 = pkin(9) + t704;
t700 = pkin(10) + t535;
t699 = qJD(2) * pkin(2);
t381 = -t403 * t557 + t561 * t424;
t375 = -pkin(10) * t487 + t381;
t368 = pkin(5) * t498 + t375;
t694 = t368 * t560;
t693 = t376 * t560;
t692 = t404 * t557;
t691 = t419 * t509;
t690 = t596 * t509;
t688 = t457 * t557;
t687 = t485 * t498;
t686 = t485 * t509;
t685 = t487 * t498;
t684 = t487 * t509;
t681 = t511 * t557;
t680 = t511 * t561;
t679 = t517 * t557;
t678 = t517 * t561;
t677 = t541 * t545;
t676 = t541 * t546;
t675 = t541 * t557;
t674 = t541 * t563;
t668 = t553 * t562;
t666 = t556 * t359;
t664 = t557 * t563;
t447 = t561 * t457;
t472 = t561 * t484;
t663 = qJDD(1) - g(3);
t662 = -pkin(10) * t680 + pkin(5) * t508 - t438 * t557 + (-t472 + (pkin(10) * t517 - t459) * t557) * qJD(5) + t709;
t589 = t517 * t646 + t681;
t661 = pkin(10) * t589 + t708;
t660 = pkin(5) * t589 + t658;
t372 = t696 * t397 - t551 * t398;
t412 = t477 * t696 - t467;
t637 = pkin(3) * t649;
t439 = pkin(4) * t509 + pkin(9) * t506 + t637;
t659 = t561 * t412 + t557 * t439;
t656 = t557 * t459 + t472;
t548 = t558 ^ 2;
t653 = -t562 ^ 2 + t548;
t648 = qJD(2) * t559;
t644 = qJD(6) * t560;
t639 = qJDD(2) * t562;
t635 = t553 * t664;
t634 = t561 * t667;
t633 = t560 * t404 - t556 * t405 - t485 * t644;
t632 = t696 * pkin(3);
t629 = t553 * t648;
t628 = qJD(2) * t667;
t627 = t517 * t647;
t624 = t562 * t641;
t621 = qJD(5) * t700;
t614 = t553 * t663;
t410 = t477 * t551 + t617;
t483 = t524 * t551 + t555 * t616;
t610 = t498 * t561;
t608 = qJD(6) * t368 + t360;
t607 = t558 * t628;
t536 = -t632 - pkin(4);
t606 = pkin(5) * t713 - t410;
t605 = g(1) * t502 - g(2) * t504;
t370 = -qJDD(3) * pkin(4) - t372;
t603 = -t519 * t453 - t496 * t654;
t432 = t561 * t439;
t515 = t700 * t561;
t602 = pkin(5) * t509 + qJD(6) * t515 - t412 * t557 + t432 + (pkin(10) * t506 + t621) * t561;
t514 = t700 * t557;
t601 = pkin(10) * t682 + qJD(6) * t514 + t557 * t621 + t659;
t362 = t368 * t556 + t693;
t451 = t561 * t459;
t388 = -pkin(5) * t587 - pkin(10) * t678 - t484 * t557 + t451;
t390 = -pkin(10) * t679 + t656;
t599 = t388 * t556 + t390 * t560;
t513 = t554 * t558 + t559 * t668;
t445 = t551 * t512 + t513 * t696;
t425 = -t445 * t557 - t634;
t590 = -t445 * t561 + t635;
t598 = t425 * t560 + t556 * t590;
t597 = t425 * t556 - t560 * t590;
t565 = qJD(2) ^ 2;
t594 = qJDD(2) * t563 - t559 * t565;
t593 = -t498 * t713 + t447;
t592 = -g(1) * t552 + g(2) * t697;
t588 = -t627 + t680;
t365 = -t487 * t645 + t633;
t584 = t402 * t498 - t535 * t457;
t583 = g(1) * (t502 * t540 + t541 * t672) + g(2) * (-t504 * t540 - t541 * t620) + g(3) * (-t540 * t669 + t541 * t554);
t582 = t594 * t553;
t580 = t513 * qJD(3);
t579 = -g(3) * t669 + t605;
t578 = g(3) * t667 - t604;
t522 = -t630 - t699;
t576 = -qJD(2) * t522 - t489 - t605;
t574 = t578 * t541;
t568 = qJD(5) * t498 * t535 + t370 + t583;
t567 = -pkin(8) * qJDD(3) + (t522 + t630 - t699) * qJD(3);
t566 = -t580 - t607;
t523 = -t561 * pkin(5) + t536;
t476 = qJD(3) * t512 + t562 * t628;
t449 = t519 * t517;
t444 = -t512 * t696 + t513 * t551;
t436 = pkin(5) * t679 + t483;
t411 = t476 * t696 + t551 * t566;
t409 = t476 * t551 - t566 * t696;
t384 = t511 * t665 - t556 * t627 - t645 * t679 + (t678 * t707 + t681) * t560;
t383 = -t448 * t707 - t519 * t511;
t379 = qJD(5) * t590 - t411 * t557 + t561 * t629;
t378 = qJD(5) * t425 + t411 * t561 + t557 * t629;
t363 = pkin(5) * t405 + t370;
t361 = -t376 * t556 + t694;
t1 = [t663 * MDP(1) + MDP(3) * t582 + (-qJDD(2) * t559 - t563 * t565) * t553 * MDP(4) + (t512 * qJDD(3) + t562 * t582 + (-t580 - 0.2e1 * t607) * qJD(3)) * MDP(10) + (-qJD(3) * t476 - qJDD(3) * t513 + (-t558 * t594 - t563 * t624) * t553) * MDP(11) + (-qJD(3) * t409 - qJDD(3) * t444 + (-t461 * t563 + t506 * t648) * t553) * MDP(12) + (-qJD(3) * t411 - qJDD(3) * t445 + (-t462 * t563 + t509 * t648) * t553) * MDP(13) + (t409 * t509 - t411 * t506 + t444 * t462 - t445 * t461) * MDP(14) + (-t372 * t444 + t373 * t445 - t407 * t409 + t408 * t411 - g(3) + (-t460 * t563 + t497 * t648) * t553) * MDP(15) + (t379 * t498 + t405 * t444 + t409 * t485 + t425 * t457) * MDP(21) + (-t378 * t498 + t404 * t444 + t409 * t487 + t457 * t590) * MDP(22) + ((-qJD(6) * t597 - t378 * t556 + t379 * t560) * t496 + t598 * t453 + t409 * t419 + t444 * t366) * MDP(28) + (-(qJD(6) * t598 + t378 * t560 + t379 * t556) * t496 - t597 * t453 - t409 * t596 + t444 * t365) * MDP(29); (t381 * t508 - t386 * t587 + t483 * t405 + t451 * t457 + t709 * t498 + t658 * t485 + (-t574 + (t402 * t517 + t403 * t587 - t484 * t498) * qJD(5)) * t561 + ((-qJD(5) * t459 - t438) * t498 - t484 * t457 - (-qJD(5) * t424 - t371) * t587 + t370 * t517 + t402 * t511 + t579) * t557) * MDP(21) + (-t365 * t448 + t366 * t449 - t383 * t419 + t384 * t596) * MDP(24) + (-t365 * t449 - t383 * t596) * MDP(23) + (-t599 * t453 + (t608 * t560 - t374 + t666) * t587 - t362 * t508 + t436 * t365 - t363 * t449 + t389 * t383 - g(1) * (-t502 * t546 + t505 * t677) - g(2) * (t503 * t677 + t504 * t546) - (-t545 * t674 + t546 * t559) * t701 + ((-qJD(6) * t388 + t661) * t560 + (qJD(6) * t390 - t662) * t556) * t496 - t660 * t596) * MDP(29) + (-t365 * t587 + t383 * t496 - t449 * t453 - t508 * t596) * MDP(25) + (t405 * t587 - t457 * t679 - t485 * t508 - t498 * t589) * MDP(19) + (-t404 * t587 + t447 * t517 + t487 * t508 + t498 * t588) * MDP(18) + (t366 * t587 - t384 * t496 - t419 * t508 - t448 * t453) * MDP(26) + (-t453 * t587 + t496 * t508) * MDP(27) + (-t457 * t587 + t498 * t508) * MDP(20) + (-t372 * t517 + t373 * t587 - t407 * t511 - t408 * t508 - t461 * t484 + t462 * t483 - t506 * t657 + t509 * t658 + t579) * MDP(14) + (-t506 * t631 - qJDD(3) * t483 - t460 * t587 - t461 * t539 + t497 * t508 - t574 + (t506 * t703 - t658) * qJD(3)) * MDP(12) + ((t388 * t560 - t390 * t556) * t453 - t613 * t587 + t361 * t508 + t436 * t366 + t363 * t448 + t389 * t384 - g(1) * (-t502 * t545 - t505 * t676) - g(2) * (-t503 * t676 + t504 * t545) - (t545 * t559 + t546 * t674) * t701 + (t556 * t661 + t560 * t662) * t496 + t660 * t419 + (t362 * t587 - t496 * t599) * qJD(6)) * MDP(28) + (-t656 * t457 - t586 * t587 - t382 * t508 + t483 * t404 + t370 * t678 - g(1) * (-t502 * t561 + t505 * t675) - g(2) * (t503 * t675 + t504 * t561) - (-t541 * t664 + t559 * t561) * t701 + t708 * t498 + t658 * t487 + t588 * t402) * MDP(22) + ((-t485 * t561 - t487 * t557) * t511 + (-t692 - t405 * t561 + (t485 * t557 - t487 * t561) * qJD(5)) * t517) * MDP(17) + (t404 * t678 + t487 * t588) * MDP(16) + (t663 * t667 + t604) * MDP(3) + (t373 * t484 - t372 * t483 - t460 * t539 + t497 * t636 - g(1) * (-t502 * t555 - t505 * t539) - g(2) * (-t503 * t539 + t504 * t555) + t657 * t408 - t658 * t407 + (-t497 * t650 - g(3) * (t539 * t563 + t555 * t559)) * t553) * MDP(15) + 0.2e1 * (t558 * t639 - t641 * t653) * MDP(6) + qJDD(2) * MDP(2) + (-t559 * t614 - t605) * MDP(4) + (-t509 * t631 - qJDD(3) * t484 + t460 * t517 - t462 * t539 + t497 * t511 + t578 * t540 + (t509 * t703 - t657) * qJD(3)) * MDP(13) + (qJDD(3) * t558 + t562 * t564) * MDP(7) + (qJDD(3) * t562 - t558 * t564) * MDP(8) + (qJDD(2) * t548 + 0.2e1 * t558 * t624) * MDP(5) + (t567 * t558 + t562 * t705) * MDP(10) + (-t558 * t705 + t567 * t562) * MDP(11); MDP(8) * t639 + MDP(7) * t640 + (t558 * t576 + t592 * t668 + t529 - t702) * MDP(10) + (t372 * t632 + t373 * t704 + t407 * t410 - t408 * t412 - t497 * t637 + (-g(1) * (t502 * t558 + t552 * t668) - g(2) * (-t504 * t558 - t562 * t620) - t702) * pkin(3)) * MDP(15) + (t487 * t610 + t692) * MDP(16) + ((t404 - t687) * t561 + (-t405 - t685) * t557) * MDP(17) - t496 * t509 * MDP(27) - t498 * t509 * MDP(20) + qJDD(3) * MDP(9) + (g(3) * t513 + (-t553 * t592 - t638) * t558 + t576 * t562) * MDP(11) + (t498 * t610 - t684 + t688) * MDP(18) + (t593 + t686) * MDP(19) + (t603 + t691) * MDP(26) + (-t365 * t519 - t366 * t520 + t419 * t655 + t596 * t654) * MDP(24) + (t365 * t520 + t596 * t655) * MDP(23) + ((t408 - t410) * t509 + (-t407 + t412) * t506 + (-t461 * t551 - t462 * t696) * pkin(3)) * MDP(14) + (t410 * qJD(3) - t497 * t509 + (qJDD(3) * t696 - t506 * t649) * pkin(3) - t583 + t372) * MDP(12) + (t382 * t509 + t536 * t404 - t410 * t487 + t498 * t659 + t557 * t568 + t561 * t584) * MDP(22) + (t412 * qJD(3) + t497 * t506 + g(1) * t466 + g(2) * t464 + g(3) * t494 + (-qJDD(3) * t551 - t509 * t649) * pkin(3) - t373) * MDP(13) + (-t381 * t509 + t536 * t405 - t410 * t485 - t432 * t498 + (t412 * t498 + t584) * t557 - t568 * t561) * MDP(21) + ((-t514 * t560 - t515 * t556) * t453 + t523 * t366 + t363 * t519 - t361 * t509 + (t556 * t601 - t560 * t602) * t496 + t606 * t419 + t654 * t389 - t583 * t546) * MDP(28) + (-(-t514 * t556 + t515 * t560) * t453 + t523 * t365 + t363 * t520 + t362 * t509 + (t556 * t602 + t560 * t601) * t496 - t606 * t596 - t655 * t389 + t583 * t545) * MDP(29) + (t690 - t706) * MDP(25) + (-MDP(5) * t558 * t562 + MDP(6) * t653) * t565; (0.2e1 * t509 * qJD(3) + t600) * MDP(12) + ((t530 - t506) * qJD(3) + t573) * MDP(13) + (-t506 ^ 2 - t509 ^ 2) * MDP(14) + (t407 * t509 + t408 * t506 - t563 * t614 + t575 - t604) * MDP(15) + (t593 - t686) * MDP(21) + (-t498 ^ 2 * t561 - t684 - t688) * MDP(22) + (t603 - t691) * MDP(28) + (t690 + t706) * MDP(29); t487 * t485 * MDP(16) + (-t485 ^ 2 + t487 ^ 2) * MDP(17) + (t404 + t687) * MDP(18) + (-t405 + t685) * MDP(19) + t457 * MDP(20) + (t382 * t498 - t402 * t487 - g(1) * (-t466 * t557 + t505 * t561) - g(2) * (-t464 * t557 + t503 * t561) - g(3) * (-t494 * t557 - t634) + t569) * MDP(21) + (t381 * t498 + t402 * t485 - g(1) * (-t466 * t561 - t505 * t557) - g(2) * (-t464 * t561 - t503 * t557) - g(3) * (-t494 * t561 + t635) + t586) * MDP(22) + (t365 + t715) * MDP(25) + (-t366 - t714) * MDP(26) + (-(-t375 * t556 - t693) * t496 - t362 * qJD(6) + (-t419 * t487 + t453 * t560 - t496 * t645) * pkin(5) + t711) * MDP(28) + ((-t376 * t496 - t359) * t556 + (t375 * t496 - t608) * t560 + (-t453 * t556 + t487 * t596 - t496 * t644) * pkin(5) + t712) * MDP(29) + t710; (t633 + t715) * MDP(25) + (-t612 - t714) * MDP(26) + (t362 * t496 + t711) * MDP(28) + (-t560 * t360 + t361 * t496 - t666 + t712) * MDP(29) + (-MDP(25) * t683 + MDP(26) * t596 - MDP(28) * t362 - MDP(29) * t694) * qJD(6) + t710;];
tau = t1;
