% Calculate vector of inverse dynamics joint torques for
% S6RRPPRP1
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
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d5,theta3,theta4]';
% MDP [27x1]
%   Minimal dynamic parameter vector (fixed base model)
%   see S6RRPPRP1_convert_par2_MPV_fixb.m
% 
% Output:
% tau [6x1]
%   joint torques of inverse dynamics (contains inertial, gravitational Coriolis and centrifugal forces)

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 08:28
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tau = S6RRPPRP1_invdynJ_fixb_mdp_slag_vp(qJ, qJD, qJDD, g, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(3,1),zeros(10,1),zeros(27,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPPRP1_invdynJ_fixb_mdp_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRPPRP1_invdynJ_fixb_mdp_slag_vp: qJD has to be [6x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [6 1]), ...
  'S6RRPPRP1_invdynJ_fixb_mdp_slag_vp: qJDD has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RRPPRP1_invdynJ_fixb_mdp_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RRPPRP1_invdynJ_fixb_mdp_slag_vp: pkin has to be [10x1] (double)');
assert(isreal(MDP) && all(size(MDP) == [27 1]), ...
  'S6RRPPRP1_invdynJ_fixb_mdp_slag_vp: MDP has to be [27x1] (double)'); 

%% Symbolic Calculation
% From invdyn_joint_fixb_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 08:28:06
% EndTime: 2019-03-09 08:28:19
% DurationCPUTime: 11.94s
% Computational Cost: add. (8855->579), mult. (20558->713), div. (0->0), fcn. (15452->14), ass. (0->239)
t564 = sin(pkin(9));
t569 = sin(qJ(2));
t571 = cos(qJ(2));
t681 = cos(pkin(9));
t527 = t564 * t571 + t569 * t681;
t512 = t527 * qJD(1);
t563 = sin(pkin(10));
t565 = cos(pkin(10));
t491 = qJD(2) * t563 + t512 * t565;
t568 = sin(qJ(5));
t694 = cos(qJ(5));
t660 = t563 * t512;
t711 = qJD(2) * t565 - t660;
t596 = t694 * t711;
t436 = -t491 * t568 + t596;
t721 = t436 ^ 2;
t625 = t681 * t571;
t543 = qJD(1) * t625;
t647 = qJD(1) * t569;
t509 = t564 * t647 - t543;
t505 = qJD(5) + t509;
t720 = t436 * t505;
t606 = t568 * t711;
t718 = t491 * t694 + t606;
t695 = t718 ^ 2;
t719 = t505 * t718;
t717 = MDP(22) + MDP(24);
t716 = -MDP(23) + MDP(26);
t683 = qJ(3) + pkin(7);
t627 = qJD(2) * t683;
t507 = -qJD(3) * t569 - t571 * t627;
t535 = t683 * t569;
t473 = qJDD(2) * pkin(2) + qJD(1) * t507 - qJDD(1) * t535;
t506 = qJD(3) * t571 - t569 * t627;
t536 = t683 * t571;
t482 = qJD(1) * t506 + qJDD(1) * t536;
t422 = t473 * t681 - t564 * t482;
t421 = -qJDD(2) * pkin(3) + qJDD(4) - t422;
t560 = qJ(2) + pkin(9);
t556 = cos(t560);
t547 = g(3) * t556;
t554 = sin(t560);
t570 = sin(qJ(1));
t572 = cos(qJ(1));
t616 = g(1) * t572 + g(2) * t570;
t588 = t616 * t554 - t547;
t715 = t588 - t421;
t511 = t527 * qJD(2);
t624 = qJDD(1) * t681;
t641 = qJDD(1) * t569;
t610 = t564 * t641 - t571 * t624;
t479 = qJD(1) * t511 + t610;
t474 = qJDD(5) + t479;
t528 = t563 * t694 + t568 * t565;
t597 = -t563 * t568 + t694 * t565;
t631 = qJD(5) * t694;
t644 = qJD(5) * t568;
t698 = -t563 * t644 + t565 * t631;
t706 = -t597 * t509 - t698;
t612 = t528 * t474 - t505 * t706;
t693 = pkin(2) * t564;
t546 = qJ(4) + t693;
t684 = pkin(8) + t546;
t520 = t684 * t563;
t521 = t684 * t565;
t470 = -t568 * t520 + t521 * t694;
t559 = pkin(10) + qJ(5);
t553 = sin(t559);
t712 = t470 * t474 + t553 * t588;
t516 = t528 * qJD(5);
t649 = t528 * t509 + t516;
t613 = t597 * t474 - t505 * t649;
t710 = t491 * t509;
t663 = t554 * t572;
t664 = t554 * t570;
t705 = g(1) * t663 + g(2) * t664 - t547;
t557 = t571 * pkin(2);
t552 = t557 + pkin(1);
t595 = -t564 * t569 + t625;
t476 = -pkin(3) * t595 - qJ(4) * t527 - t552;
t489 = -t564 * t535 + t536 * t681;
t425 = t565 * t476 - t489 * t563;
t667 = t527 * t565;
t409 = -pkin(4) * t595 - pkin(8) * t667 + t425;
t426 = t563 * t476 + t565 * t489;
t668 = t527 * t563;
t414 = -pkin(8) * t668 + t426;
t704 = t568 * t409 + t694 * t414;
t457 = pkin(2) * t647 + pkin(3) * t512 + qJ(4) * t509;
t531 = qJD(1) * t536;
t517 = t564 * t531;
t530 = qJD(1) * t535;
t484 = -t530 * t681 - t517;
t420 = t563 * t457 + t565 * t484;
t703 = -qJD(4) * t565 + t420;
t419 = t565 * t457 - t484 * t563;
t671 = t509 * t565;
t399 = pkin(4) * t512 + pkin(8) * t671 + t419;
t672 = t509 * t563;
t410 = pkin(8) * t672 + t420;
t598 = -t520 * t694 - t568 * t521;
t702 = -qJD(4) * t597 - qJD(5) * t598 + t568 * t399 + t694 * t410;
t701 = -qJD(4) * t528 - qJD(5) * t470 - t399 * t694 + t568 * t410;
t700 = g(1) * t570 - g(2) * t572;
t642 = qJD(1) * qJD(2);
t630 = t569 * t642;
t591 = qJD(2) * t543 - t564 * t630;
t579 = qJDD(1) * t527 + t591;
t578 = t563 * qJDD(2) + t565 * t579;
t468 = t563 * t579;
t622 = qJDD(2) * t565 - t468;
t384 = -qJD(5) * t596 + t491 * t644 - t568 * t622 - t694 * t578;
t697 = -t384 * t597 - t649 * t718;
t514 = t595 * qJD(2);
t682 = qJD(2) * pkin(2);
t636 = t569 * t682;
t437 = pkin(3) * t511 - qJ(4) * t514 - qJD(4) * t527 + t636;
t459 = t506 * t681 + t564 * t507;
t405 = t565 * t437 - t459 * t563;
t669 = t514 * t565;
t388 = pkin(4) * t511 - pkin(8) * t669 + t405;
t406 = t563 * t437 + t565 * t459;
t670 = t514 * t563;
t394 = -pkin(8) * t670 + t406;
t696 = -qJD(5) * t704 + t388 * t694 - t568 * t394;
t508 = t509 ^ 2;
t692 = pkin(2) * t569;
t691 = pkin(5) * t474;
t687 = g(3) * t554;
t686 = g(3) * t571;
t685 = t565 * pkin(4);
t680 = qJ(6) * t474;
t679 = qJDD(1) * pkin(1);
t534 = -qJD(1) * t552 + qJD(3);
t449 = pkin(3) * t509 - qJ(4) * t512 + t534;
t522 = -t530 + t682;
t626 = t681 * t531;
t478 = t564 * t522 + t626;
t467 = qJD(2) * qJ(4) + t478;
t411 = t565 * t449 - t467 * t563;
t391 = pkin(4) * t509 - pkin(8) * t491 + t411;
t412 = t563 * t449 + t565 * t467;
t400 = pkin(8) * t711 + t412;
t371 = t568 * t391 + t400 * t694;
t678 = t371 * t505;
t676 = t436 * t512;
t675 = t718 * t436;
t674 = t718 * t512;
t666 = t546 * t563;
t567 = -pkin(8) - qJ(4);
t665 = t554 * t567;
t662 = t556 * t570;
t661 = t556 * t572;
t658 = t565 * t479;
t555 = cos(t559);
t657 = t570 * t555;
t656 = t572 * t553;
t639 = pkin(2) * t630 + qJDD(3);
t640 = qJDD(1) * t571;
t403 = -pkin(2) * t640 + t479 * pkin(3) - qJ(4) * t579 - t512 * qJD(4) + t639 - t679;
t423 = t564 * t473 + t681 * t482;
t418 = qJDD(2) * qJ(4) + qJD(2) * qJD(4) + t423;
t381 = t563 * t403 + t565 * t418;
t653 = -qJ(6) * t512 - t702;
t652 = t512 * pkin(5) - t701;
t483 = -t530 * t564 + t626;
t440 = -pkin(4) * t672 + t483;
t651 = pkin(5) * t649 + qJ(6) * t706 - qJD(6) * t528 - t440;
t561 = t569 ^ 2;
t648 = -t571 ^ 2 + t561;
t370 = t391 * t694 - t568 * t400;
t643 = qJD(6) - t370;
t633 = t681 * pkin(2);
t629 = pkin(4) * t563 + t683;
t380 = t565 * t403 - t563 * t418;
t458 = t506 * t564 - t681 * t507;
t488 = t681 * t535 + t536 * t564;
t368 = t479 * pkin(4) - pkin(8) * t578 + t380;
t376 = pkin(8) * t622 + t381;
t621 = -t694 * t368 + t568 * t376 + t391 * t644 + t400 * t631;
t619 = g(1) * t664 - g(2) * t663;
t551 = -t633 - pkin(3);
t494 = t553 * t662 + t555 * t572;
t496 = t556 * t656 - t657;
t618 = -g(1) * t494 + g(2) * t496;
t495 = t556 * t657 - t656;
t497 = t553 * t570 + t555 * t661;
t617 = g(1) * t495 - g(2) * t497;
t577 = t568 * t578 - t694 * t622;
t385 = qJD(5) * t718 + t577;
t614 = -t528 * t385 - t436 * t706;
t428 = pkin(4) * t670 + t458;
t455 = pkin(4) * t668 + t488;
t611 = pkin(3) * t556 + qJ(4) * t554;
t477 = t522 * t681 - t517;
t609 = -t380 * t563 + t381 * t565;
t608 = -t411 * t563 + t412 * t565;
t550 = pkin(3) + t685;
t607 = t550 * t556 - t665;
t605 = pkin(5) * t555 + qJ(6) * t553 + t550;
t603 = -0.2e1 * pkin(1) * t642 - pkin(7) * qJDD(2);
t601 = t409 * t694 - t568 * t414;
t533 = t551 - t685;
t594 = t568 * t368 + t694 * t376 + t391 * t631 - t400 * t644;
t593 = t568 * t388 + t694 * t394 + t409 * t631 - t414 * t644;
t590 = -qJDD(1) * t552 + t639;
t589 = t598 * t474 + t555 * t705;
t461 = -qJD(2) * pkin(3) + qJD(4) - t477;
t573 = qJD(2) ^ 2;
t586 = -pkin(7) * t573 + 0.2e1 * t679 + t700;
t574 = qJD(1) ^ 2;
t585 = pkin(1) * t574 - pkin(7) * qJDD(1) + t616;
t584 = g(1) * t496 + g(2) * t494 + t553 * t687 - t621;
t427 = -pkin(4) * t711 + t461;
t382 = -pkin(5) * t436 - qJ(6) * t718 + t427;
t582 = t382 * t718 + qJDD(6) - t584;
t395 = -pkin(4) * t622 + t421;
t581 = -g(1) * t497 - g(2) * t495 - t555 * t687 + t594;
t362 = t385 * pkin(5) + t384 * qJ(6) - qJD(6) * t718 + t395;
t576 = -qJD(5) * t606 - t491 * t631 - t577;
t539 = t572 * t552;
t463 = t597 * t527;
t462 = t528 * t527;
t460 = -pkin(5) * t597 - t528 * qJ(6) + t533;
t417 = t514 * t528 + t527 * t698;
t416 = -t514 * t597 + t516 * t527;
t398 = pkin(5) * t718 - qJ(6) * t436;
t393 = pkin(5) * t462 - qJ(6) * t463 + t455;
t379 = pkin(5) * t595 - t601;
t378 = -qJ(6) * t595 + t704;
t377 = -t384 - t720;
t369 = pkin(5) * t417 + qJ(6) * t416 - qJD(6) * t463 + t428;
t367 = t505 * qJ(6) + t371;
t366 = -t505 * pkin(5) + t643;
t363 = -t511 * pkin(5) - t696;
t361 = qJ(6) * t511 - qJD(6) * t595 + t593;
t360 = qJDD(6) + t621 - t691;
t359 = qJD(6) * t505 + t594 + t680;
t1 = [0.2e1 * (t569 * t640 - t642 * t648) * MDP(5) + (t359 * t378 + t367 * t361 + t362 * t393 + t382 * t369 + t360 * t379 + t366 * t363 - g(1) * (-pkin(5) * t495 - qJ(6) * t494) - g(2) * (pkin(5) * t497 + qJ(6) * t496 + t539) + (-g(1) * t629 - g(2) * t607) * t572 + (-g(1) * (-t552 - t607) - g(2) * t629) * t570) * MDP(27) + (qJDD(1) * t561 + 0.2e1 * t571 * t630) * MDP(4) + t616 * MDP(3) + (t569 * t603 + t571 * t586) * MDP(9) + (-t569 * t586 + t571 * t603) * MDP(10) + (t370 * t511 + t455 * t385 + t395 * t462 + t427 * t417 - t428 * t436 + t601 * t474 + t505 * t696 + t595 * t621 + t617) * MDP(22) + (t360 * t595 + t362 * t462 - t363 * t505 - t366 * t511 - t369 * t436 - t379 * t474 + t382 * t417 + t385 * t393 + t617) * MDP(24) + (t385 * t595 - t417 * t505 + t436 * t511 - t462 * t474) * MDP(20) + (-t359 * t462 + t360 * t463 + t361 * t436 + t363 * t718 - t366 * t416 - t367 * t417 - t378 * t385 - t379 * t384 + t619) * MDP(25) + (t384 * t462 - t385 * t463 - t416 * t436 - t417 * t718) * MDP(18) + (-t474 * t595 + t505 * t511) * MDP(21) + (-t406 * t509 - t426 * t479 + t381 * t595 - t412 * t511 + t458 * t491 + t488 * t578 + t421 * t667 + t461 * t669 - g(1) * (t563 * t662 + t565 * t572) - g(2) * (-t563 * t661 + t565 * t570)) * MDP(14) + (t423 * t489 + t478 * t459 - t422 * t488 - t477 * t458 - t590 * t552 + t534 * t636 - g(1) * (-t552 * t570 + t572 * t683) - g(2) * (t570 * t683 + t539)) * MDP(12) + (-g(2) * t539 + t380 * t425 + t381 * t426 + t411 * t405 + t412 * t406 + t421 * t488 + t461 * t458 + (-g(1) * t683 - g(2) * t611) * t572 + (-g(1) * (-t552 - t611) - g(2) * t683) * t570) * MDP(16) + t700 * MDP(2) + (qJDD(2) * t569 + t571 * t573) * MDP(6) + (qJDD(2) * t571 - t569 * t573) * MDP(7) + (t405 * t509 + t425 * t479 - t380 * t595 + t411 * t511 - t458 * t711 - t488 * t622 + t700 * t565 * t556 + (t421 * t527 + t461 * t514 - t616) * t563) * MDP(13) + (-t380 * t667 - t381 * t668 - t405 * t491 + t406 * t711 - t411 * t669 - t412 * t670 - t425 * t578 + t426 * t622 + t619) * MDP(15) + (-t422 * t527 + t423 * t595 + t458 * t512 - t459 * t509 - t477 * t514 - t478 * t511 - t489 * t479 + t488 * t579 - t616) * MDP(11) + (-t384 * t463 - t416 * t718) * MDP(17) + (-t359 * t595 + t361 * t505 - t362 * t463 + t367 * t511 - t369 * t718 + t378 * t474 + t382 * t416 + t384 * t393 - t618) * MDP(26) + (t384 * t595 - t416 * t505 + t463 * t474 + t511 * t718) * MDP(19) + (-t371 * t511 - t455 * t384 + t395 * t463 - t427 * t416 + t428 * t718 - t474 * t704 - t505 * t593 + t594 * t595 + t618) * MDP(23) + qJDD(1) * MDP(1); -t505 * t512 * MDP(21) + MDP(6) * t641 + MDP(7) * t640 + qJDD(2) * MDP(8) + (t569 * t585 - t686) * MDP(9) + (g(3) * t569 + t571 * t585) * MDP(10) + (-t479 * t693 - t579 * t633 - (-t478 + t483) * t512 + (t484 - t477) * t509) * MDP(11) + (t477 * t483 - t478 * t484 + (t681 * t422 - t686 + t423 * t564 + (-qJD(1) * t534 + t616) * t569) * pkin(2)) * MDP(12) + (-t479 * t666 - t411 * t512 + t551 * t468 - t483 * t660 + (-t419 + (-qJD(4) + t461) * t563) * t509 + (t483 * qJD(2) - t551 * qJDD(2) + t715) * t565) * MDP(13) + (t412 * t512 + t461 * t671 - t483 * t491 - t546 * t658 + t551 * t578 + t703 * t509 + (t421 - t705) * t563) * MDP(14) + (t546 * t565 * t622 - g(1) * t661 - g(2) * t662 - t411 * t671 - t412 * t672 + t578 * t666 + t609 - t687 - t703 * t711 + (qJD(4) * t563 + t419) * t491) * MDP(15) + (t421 * t551 - t412 * t420 - t411 * t419 - t461 * t483 - g(3) * (t557 + t611) + t609 * t546 + t608 * qJD(4) + t616 * (pkin(3) * t554 - qJ(4) * t556 + t692)) * MDP(16) + (-t384 * t528 - t706 * t718) * MDP(17) + (t614 + t697) * MDP(18) + (t612 - t674) * MDP(19) + (t613 - t676) * MDP(20) + (-t370 * t512 + t533 * t385 - t395 * t597 + t649 * t427 + t436 * t440 + t505 * t701 + t589) * MDP(22) + (t371 * t512 - t533 * t384 + t395 * t528 - t427 * t706 - t440 * t718 + t505 * t702 - t712) * MDP(23) + (-t362 * t597 + t366 * t512 + t382 * t649 + t385 * t460 - t436 * t651 - t505 * t652 + t589) * MDP(24) + (t359 * t597 + t360 * t528 - t366 * t706 - t367 * t649 + t384 * t598 - t385 * t470 + t436 * t653 - t556 * t616 + t652 * t718 - t687) * MDP(25) + (-t362 * t528 - t367 * t512 + t382 * t706 + t384 * t460 + t505 * t653 - t651 * t718 + t712) * MDP(26) + (t359 * t470 + t362 * t460 - t360 * t598 - g(3) * (t557 - t665) - t605 * t547 + t651 * t382 + t653 * t367 + t652 * t366 + t616 * (t554 * t605 + t556 * t567 + t692)) * MDP(27) + (-MDP(4) * t569 * t571 + MDP(5) * t648) * t574; (-t512 ^ 2 - t508) * MDP(11) + (t477 * t512 + t478 * t509 + t590 - t700) * MDP(12) + (-t508 * t563 + t512 * t711 + t658) * MDP(13) + (-t479 * t563 - t491 * t512 - t508 * t565) * MDP(14) + ((-t509 * t660 + (t509 * qJD(2) - t564 * t640 - t569 * t624 - t591) * t565) * t565 + (-t468 + t710) * t563) * MDP(15) + (t380 * t565 + t381 * t563 - t461 * t512 + t509 * t608 - t700) * MDP(16) + (t614 - t697) * MDP(25) + (t359 * t528 - t360 * t597 + t366 * t649 - t367 * t706 - t382 * t512 - t700) * MDP(27) + t716 * (t612 + t674) + t717 * (t613 + t676); (-t622 + t710) * MDP(13) + (t509 * t711 + t578) * MDP(14) + (-t491 ^ 2 - t711 ^ 2) * MDP(15) + (t411 * t491 - t412 * t711 - t715) * MDP(16) + (-t695 - t721) * MDP(25) + (-t366 * t718 - t367 * t436 + t362 - t588) * MDP(27) + t716 * (t384 - t720) + t717 * (-t576 + t719); -MDP(17) * t675 + (t695 - t721) * MDP(18) + t377 * MDP(19) + (t576 + t719) * MDP(20) + t474 * MDP(21) + (-t427 * t718 + t584 + t678) * MDP(22) + (t370 * t505 - t427 * t436 - t581) * MDP(23) + (t398 * t436 - t582 + t678 + 0.2e1 * t691) * MDP(24) + (pkin(5) * t384 - qJ(6) * t385 + (t367 - t371) * t718 - (t366 - t643) * t436) * MDP(25) + (0.2e1 * t680 + t382 * t436 + t398 * t718 + (0.2e1 * qJD(6) - t370) * t505 + t581) * MDP(26) + (t359 * qJ(6) - t360 * pkin(5) - t382 * t398 - t366 * t371 - g(1) * (-pkin(5) * t496 + qJ(6) * t497) - g(2) * (-pkin(5) * t494 + qJ(6) * t495) - (-pkin(5) * t553 + qJ(6) * t555) * t687 + t643 * t367) * MDP(27); (-qJD(2) * t512 - qJDD(5) - t610 - t675) * MDP(24) + t377 * MDP(25) + (-t505 ^ 2 - t695) * MDP(26) + (-t367 * t505 + t582 - t691) * MDP(27);];
tau  = t1;
