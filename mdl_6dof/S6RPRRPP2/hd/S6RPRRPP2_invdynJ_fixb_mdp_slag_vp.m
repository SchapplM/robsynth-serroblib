% Calculate vector of inverse dynamics joint torques for
% S6RPRRPP2
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
%   pkin=[a2,a3,a4,a5,a6,d1,d3,d4,theta2]';
% MDP [26x1]
%   Minimal dynamic parameter vector (fixed base model)
%   see S6RPRRPP2_convert_par2_MPV_fixb.m
% 
% Output:
% tau [6x1]
%   joint torques of inverse dynamics (contains inertial, gravitational Coriolis and centrifugal forces)

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 04:34
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tau = S6RPRRPP2_invdynJ_fixb_mdp_slag_vp(qJ, qJD, qJDD, g, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(3,1),zeros(9,1),zeros(26,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRRPP2_invdynJ_fixb_mdp_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RPRRPP2_invdynJ_fixb_mdp_slag_vp: qJD has to be [6x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [6 1]), ...
  'S6RPRRPP2_invdynJ_fixb_mdp_slag_vp: qJDD has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RPRRPP2_invdynJ_fixb_mdp_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S6RPRRPP2_invdynJ_fixb_mdp_slag_vp: pkin has to be [9x1] (double)');
assert(isreal(MDP) && all(size(MDP) == [26 1]), ...
  'S6RPRRPP2_invdynJ_fixb_mdp_slag_vp: MDP has to be [26x1] (double)'); 

%% Symbolic Calculation
% From invdyn_joint_fixb_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 04:33:37
% EndTime: 2019-03-09 04:33:45
% DurationCPUTime: 7.13s
% Computational Cost: add. (4194->581), mult. (8242->677), div. (0->0), fcn. (5111->10), ass. (0->232)
t684 = MDP(20) - MDP(25);
t526 = cos(qJ(3));
t507 = t526 * qJDD(1);
t523 = sin(qJ(3));
t603 = qJD(1) * qJD(3);
t445 = t523 * t603 + qJDD(4) - t507;
t618 = qJD(1) * t526;
t682 = qJD(4) - t618;
t683 = t445 * qJ(5) + qJD(5) * t682;
t516 = qJ(1) + pkin(9);
t503 = cos(t516);
t641 = t503 * t523;
t502 = sin(t516);
t659 = g(2) * t502;
t675 = g(1) * t641 + t523 * t659;
t681 = MDP(21) + MDP(24);
t587 = t526 * t603;
t601 = qJDD(1) * t523;
t680 = qJD(3) * qJD(4) + t587 + t601;
t522 = sin(qJ(4));
t525 = cos(qJ(4));
t612 = qJD(4) * t523;
t586 = qJD(1) * t612;
t553 = (-qJDD(3) + t586) * t525;
t388 = t522 * (qJD(3) * (qJD(4) + t618) + t601) + t553;
t639 = t522 * qJ(5);
t666 = pkin(4) + pkin(5);
t551 = -t525 * t666 - t639;
t438 = pkin(3) - t551;
t607 = t525 * qJD(3);
t619 = qJD(1) * t523;
t449 = t522 * t619 - t607;
t679 = t388 * qJ(6) + t449 * qJD(6);
t520 = sin(pkin(9));
t496 = pkin(1) * t520 + pkin(7);
t466 = t496 * qJD(1);
t508 = t523 * qJD(2);
t423 = t526 * t466 + t508;
t678 = qJD(5) * t522 + t423;
t677 = 0.2e1 * t683;
t676 = t526 * qJD(2) - t523 * t466;
t673 = g(1) * t503 + t659;
t407 = qJD(3) * pkin(8) + t423;
t511 = t523 * pkin(8);
t513 = t526 * pkin(3);
t594 = -pkin(2) - t513;
t558 = t594 - t511;
t521 = cos(pkin(9));
t665 = pkin(1) * t521;
t544 = t558 - t665;
t412 = t544 * qJD(1);
t369 = -t522 * t407 + t525 * t412;
t605 = qJD(5) - t369;
t672 = MDP(19) + MDP(23);
t434 = t445 * pkin(4);
t671 = t434 - qJDD(5);
t617 = qJD(3) * t522;
t451 = t525 * t619 + t617;
t573 = qJD(3) * pkin(3) + t676;
t545 = qJ(5) * t451 + t573;
t368 = pkin(4) * t449 - t545;
t662 = pkin(8) * t445;
t670 = -t368 * t682 + t662;
t359 = -t449 * t666 + qJD(6) + t545;
t637 = t522 * t526;
t414 = t502 * t637 + t503 * t525;
t416 = -t502 * t525 + t503 * t637;
t464 = t496 * qJDD(1);
t600 = qJDD(2) * t523;
t379 = qJDD(3) * pkin(8) + qJD(3) * t676 + t464 * t526 + t600;
t572 = pkin(3) * t523 - pkin(8) * t526;
t456 = t572 * qJD(3);
t389 = qJD(1) * t456 + qJDD(1) * t544;
t611 = qJD(4) * t525;
t613 = qJD(4) * t522;
t577 = t522 * t379 - t525 * t389 + t407 * t611 + t412 * t613;
t638 = t522 * t523;
t538 = g(1) * t416 + g(2) * t414 + g(3) * t638 - t577;
t536 = t538 + t671;
t576 = t522 * qJDD(3) + t525 * t680;
t387 = t522 * t586 - t576;
t652 = qJ(6) * t387;
t669 = (qJD(6) + t359) * t451 + t536 - t652;
t667 = t449 ^ 2;
t444 = t451 ^ 2;
t664 = pkin(4) * t522;
t663 = pkin(5) * t445;
t661 = g(1) * t502;
t658 = g(3) * t526;
t657 = pkin(8) - qJ(6);
t656 = pkin(8) * qJD(4);
t655 = qJ(5) * t388;
t654 = qJ(5) * t449;
t653 = qJ(5) * t525;
t370 = t525 * t407 + t522 * t412;
t361 = qJ(6) * t449 + t370;
t482 = t682 * qJ(5);
t355 = t361 + t482;
t651 = t355 * t682;
t363 = t482 + t370;
t650 = t363 * t682;
t649 = t370 * t682;
t648 = t387 * t522;
t647 = t449 * t451;
t646 = t449 * t682;
t645 = t451 * t682;
t644 = t682 * t525;
t643 = t496 * t522;
t642 = t502 * t526;
t640 = t503 * t526;
t636 = t523 * t525;
t635 = t525 * t526;
t592 = t526 * t607;
t634 = -t388 * t636 - t449 * t592;
t599 = t666 * t522;
t552 = -t599 + t653;
t633 = t552 * t682 + t678;
t455 = t572 * qJD(1);
t632 = t522 * t455 + t525 * t676;
t631 = t445 * t636 + t592 * t682;
t497 = -pkin(2) - t665;
t622 = t511 + t513;
t437 = t497 - t622;
t454 = t496 * t635;
t630 = qJD(4) * t454 + t437 * t613;
t629 = t437 * t611 + t522 * t456;
t375 = qJ(5) * t619 + t632;
t608 = qJD(6) * t525;
t628 = -qJ(6) * t522 * t618 - t613 * t657 - t375 - t608;
t566 = -t653 + t664;
t627 = t566 * t682 - t678;
t470 = t657 * t525;
t404 = t522 * t676;
t580 = -t455 * t525 + t404;
t598 = t666 * t523;
t626 = qJD(4) * t470 - qJD(6) * t522 - (-qJ(6) * t635 - t598) * qJD(1) - t580;
t625 = t522 * t437 + t454;
t624 = t675 * t522;
t623 = t675 * t525;
t517 = t523 ^ 2;
t621 = -t526 ^ 2 + t517;
t620 = MDP(20) * t522;
t467 = qJD(1) * t497;
t616 = qJD(3) * t523;
t615 = qJD(3) * t526;
t614 = qJD(4) * t449;
t609 = qJD(5) * t525;
t360 = qJ(6) * t451 + t369;
t606 = qJD(5) - t360;
t597 = t525 * t379 + t522 * t389 + t412 * t611;
t588 = t523 * t611;
t593 = t682 * t617;
t596 = -t445 * t638 - t526 * t593 - t588 * t682;
t595 = g(1) * t640 + g(2) * t642 + g(3) * t523;
t591 = t359 * t613;
t590 = t359 * t611;
t589 = t522 * t612;
t585 = -pkin(4) - t643;
t574 = -qJD(3) * t508 + t526 * qJDD(2) - t523 * t464 - t466 * t615;
t380 = -qJDD(3) * pkin(3) - t574;
t351 = t388 * pkin(4) + t387 * qJ(5) - t451 * qJD(5) + t380;
t348 = -pkin(5) * t388 + qJDD(6) - t351;
t584 = t348 - t658;
t583 = MDP(18) - t681;
t415 = t502 * t635 - t503 * t522;
t582 = -t414 * pkin(4) + qJ(5) * t415;
t417 = t502 * t522 + t503 * t635;
t581 = -t416 * pkin(4) + qJ(5) * t417;
t453 = t496 * t637;
t579 = t437 * t525 - t453;
t575 = pkin(4) * t635 + qJ(5) * t637 + t622;
t571 = g(1) * t414 - g(2) * t416;
t570 = g(1) * t415 - g(2) * t417;
t524 = sin(qJ(1));
t527 = cos(qJ(1));
t569 = g(1) * t524 - g(2) * t527;
t391 = -qJ(5) * t526 + t625;
t568 = t456 * t525 - t630;
t567 = pkin(4) * t525 + t639;
t548 = -t407 * t613 + t597;
t349 = t548 + t683;
t350 = t577 - t671;
t565 = t349 * t525 + t350 * t522;
t354 = -t666 * t682 + t606;
t564 = t354 * t525 - t355 * t522;
t563 = t354 * t522 + t355 * t525;
t362 = -pkin(4) * t682 + t605;
t562 = t362 * t525 - t363 * t522;
t561 = t362 * t522 + t363 * t525;
t557 = MDP(17) * t522 + MDP(18) * t525;
t556 = pkin(3) + t567;
t555 = -t656 * t682 - t658;
t554 = qJ(5) * t616 - qJD(5) * t526 + t629;
t550 = t445 * t522 + t611 * t682;
t549 = -t351 + t555;
t546 = -pkin(1) * t524 - t415 * pkin(4) + t503 * pkin(7) - qJ(5) * t414;
t543 = -t573 * t682 - t662;
t542 = t527 * pkin(1) + t503 * pkin(2) + pkin(3) * t640 + t417 * pkin(4) + t502 * pkin(7) + pkin(8) * t641 + qJ(5) * t416;
t529 = qJD(3) ^ 2;
t541 = 0.2e1 * qJDD(1) * t497 + t496 * t529 - t661;
t539 = 0.2e1 * qJD(3) * t467 - qJDD(3) * t496;
t537 = t387 - t646;
t534 = t368 * t451 - t536;
t533 = g(1) * t417 + g(2) * t415 + g(3) * t636 - t548;
t532 = t369 * t682 + t533;
t512 = t526 * pkin(4);
t487 = qJ(5) * t636;
t477 = g(2) * t641;
t473 = pkin(8) * t640;
t471 = pkin(8) * t642;
t469 = t657 * t522;
t463 = qJDD(3) * t526 - t523 * t529;
t462 = qJDD(3) * t523 + t526 * t529;
t426 = t449 * t616;
t425 = t451 * t616;
t400 = -t487 + (t496 + t664) * t523;
t394 = pkin(4) * t451 + t654;
t393 = t487 + (-t496 - t599) * t523;
t392 = t512 - t579;
t381 = qJ(6) * t638 + t391;
t377 = -pkin(4) * t619 + t580;
t374 = -t451 * t666 - t654;
t371 = pkin(5) * t526 + t453 + t512 + (-qJ(6) * t523 - t437) * t525;
t366 = (qJD(4) * t567 - t609) * t523 + (t496 + t566) * t615;
t358 = t585 * t616 - t568;
t357 = (qJD(4) * t551 + t609) * t523 + (-t496 + t552) * t615;
t356 = (-t523 * t607 - t526 * t613) * t496 + t554;
t353 = (qJ(6) * qJD(4) - qJD(3) * t496) * t636 + (qJD(6) * t523 + (qJ(6) * qJD(3) - qJD(4) * t496) * t526) * t522 + t554;
t352 = (-qJ(6) * t615 - t456) * t525 + (qJ(6) * t613 - t608 + (-pkin(5) + t585) * qJD(3)) * t523 + t630;
t347 = t349 + t679;
t346 = -qJD(6) * t451 + t350 + t652 - t663;
t1 = [(g(1) * t527 + g(2) * t524) * MDP(3) + t462 * MDP(7) + t463 * MDP(8) + (t353 * t682 + t357 * t451 + t381 * t445 - t387 * t393 + (t359 * t607 - t347) * t526 + (qJD(3) * t355 + t348 * t525 - t591) * t523 + t571) * MDP(24) + (t356 * t682 - t366 * t451 + t387 * t400 + t391 * t445 + (-t368 * t607 - t349) * t526 + (qJD(3) * t363 - t351 * t525 + t368 * t613) * t523 + t571) * MDP(21) + (-t629 * t682 - t625 * t445 + ((t496 * t682 - t407) * t613 + (t451 * t496 - t525 * t573) * qJD(3) + t597) * t526 + (t573 * t613 + t380 * t525 - t496 * t387 + (t496 * t644 - t370) * qJD(3)) * t523 - t571) * MDP(18) + 0.2e1 * (t507 * t523 - t603 * t621) * MDP(6) + t569 * MDP(2) + (t347 * t381 + t355 * t353 + t346 * t371 + t354 * t352 + t348 * t393 + t359 * t357 - g(1) * (-pkin(5) * t415 + t546) - g(2) * (pkin(5) * t417 - qJ(6) * t641 + t542) - (-t523 * t657 + t594) * t661) * MDP(26) + (-g(1) * t546 - g(2) * t542 + t349 * t391 + t350 * t392 + t351 * t400 + t363 * t356 + t362 * t358 + t368 * t366 - t558 * t661) * MDP(22) + (-t352 * t451 + t353 * t449 + t371 * t387 + t381 * t388 + t477 - t564 * t615 + (qJD(4) * t563 - t346 * t525 + t347 * t522 - t661) * t523) * MDP(25) + (-t356 * t449 + t358 * t451 - t387 * t392 - t388 * t391 - t477 + t562 * t615 + (-qJD(4) * t561 - t349 * t522 + t350 * t525 + t661) * t523) * MDP(20) + (qJDD(1) * t517 + 0.2e1 * t523 * t587) * MDP(5) + (-t445 * t526 + t616 * t682) * MDP(16) + (-t352 * t682 - t357 * t449 - t371 * t445 - t388 * t393 + (-t359 * t617 + t346) * t526 + (-qJD(3) * t354 - t348 * t522 - t590) * t523 + t570) * MDP(23) + (-t358 * t682 + t366 * t449 + t388 * t400 - t392 * t445 + (t368 * t617 + t350) * t526 + (-qJD(3) * t362 + t351 * t522 + t368 * t611) * t523 + t570) * MDP(19) + (t568 * t682 + t579 * t445 + ((t449 * t496 - t522 * t573) * qJD(3) + t577) * t526 + (-t573 * t611 + t380 * t522 + t496 * t388 + (t643 * t682 + t369) * qJD(3)) * t523 + t570) * MDP(17) + (t387 * t526 - t589 * t682 + t425 + t631) * MDP(14) + ((t388 - t593) * t526 + (-qJD(3) * t449 - t550) * t523) * MDP(15) + (-t387 * t636 + (-t589 + t592) * t451) * MDP(12) + (t539 * t523 + (-g(2) * t503 - t541) * t526) * MDP(10) + (t523 * t541 + t526 * t539 + t477) * MDP(11) + qJDD(1) * MDP(1) + (t569 + (t520 ^ 2 + t521 ^ 2) * qJDD(1) * pkin(1)) * pkin(1) * MDP(4) + (-t451 * t588 + (-t451 * t615 + (t387 + t614) * t523) * t522 + t634) * MDP(13); qJDD(2) * MDP(4) + t463 * MDP(10) - t462 * MDP(11) + t426 * MDP(17) + t425 * MDP(18) + (t426 + t596) * MDP(19) + t634 * MDP(20) + t631 * MDP(21) + t596 * MDP(23) + (-t425 + t631) * MDP(24) + (-MDP(4) - MDP(22) - MDP(26)) * g(3) + (-t351 * MDP(22) + t348 * MDP(26) + (-MDP(17) - t672) * t388 + t583 * t387 + (t451 * t620 + t561 * MDP(22) + (t449 * t525 - t451 * t522) * MDP(25) + t563 * MDP(26) - t557 * t682) * qJD(3)) * t526 + (-t387 * t620 + t565 * MDP(22) + (t388 * t525 + t648) * MDP(25) + (t346 * t522 + t347 * t525) * MDP(26) - t557 * t445 + (-MDP(21) * t451 + MDP(22) * t368 + MDP(23) * t449 - MDP(26) * t359) * qJD(3) + (t562 * MDP(22) + t564 * MDP(26) - (t525 * MDP(17) - t522 * t583) * t682 + t684 * (t449 * t522 + t451 * t525)) * qJD(4)) * t523; MDP(7) * t601 + MDP(8) * t507 + qJDD(3) * MDP(9) + (qJD(3) * t423 - t467 * t619 + t574 - t658 + t675) * MDP(10) + (-t600 + (-qJD(1) * t467 - t464) * t526 + t595) * MDP(11) + (t451 * t644 - t648) * MDP(12) + ((-t387 - t646) * t525 + (-t388 - t645) * t522) * MDP(13) + ((-t451 * t523 - t635 * t682) * qJD(1) + t550) * MDP(14) + (-t682 * t613 + t445 * t525 + (t449 * t523 + t637 * t682) * qJD(1)) * MDP(15) - t682 * MDP(16) * t619 + (-t369 * t619 - pkin(3) * t388 + t404 * t682 - t423 * t449 + (-t658 - t380 - (t455 + t656) * t682) * t525 + t543 * t522 + t623) * MDP(17) + (pkin(3) * t387 + t632 * t682 + t370 * t619 - t423 * t451 + t543 * t525 + (t380 - t555) * t522 - t624) * MDP(18) + (t362 * t619 + t377 * t682 - t388 * t556 + t627 * t449 - t522 * t670 + t549 * t525 + t623) * MDP(19) + (t375 * t449 - t377 * t451 + (t349 + t682 * t362 + (qJD(4) * t451 - t388) * pkin(8)) * t525 + (t350 - t650 + (-t387 + t614) * pkin(8)) * t522 - t595) * MDP(20) + (-t363 * t619 - t375 * t682 - t387 * t556 - t627 * t451 + t549 * t522 + t525 * t670 + t624) * MDP(21) + (-t363 * t375 - t362 * t377 - g(1) * t473 - g(2) * t471 - g(3) * t575 + t627 * t368 + (qJD(4) * t562 + t565) * pkin(8) + (t673 * t523 - t351) * t556) * MDP(22) + (-t591 - t388 * t438 - t445 * t469 + t584 * t525 - t626 * t682 - t633 * t449 + (t354 * t523 + t359 * t637) * qJD(1) + t623) * MDP(23) + (t590 - t387 * t438 + t445 * t470 + t584 * t522 + t628 * t682 + t633 * t451 + (-t355 * t523 - t359 * t635) * qJD(1) + t624) * MDP(24) + (t387 * t469 + t388 * t470 - t626 * t451 + t628 * t449 + (-t354 * t682 - t347) * t525 + (-t346 + t651) * t522 + t595) * MDP(25) + (t347 * t470 + t346 * t469 + t348 * t438 - g(1) * (-qJ(6) * t640 + t473) - g(2) * (-qJ(6) * t642 + t471) - g(3) * (pkin(5) * t635 + t575) + t633 * t359 + t628 * t355 + t626 * t354 + (g(3) * qJ(6) + t438 * t673) * t523) * MDP(26) + (-MDP(5) * t523 * t526 + MDP(6) * t621) * qJD(1) ^ 2; MDP(12) * t647 + (t444 - t667) * MDP(13) - t537 * MDP(14) + (-t388 + t645) * MDP(15) + t445 * MDP(16) + (t451 * t573 + t538 + t649) * MDP(17) + (-t449 * t573 + t532) * MDP(18) + (-t394 * t449 + t434 - t534 + t649) * MDP(19) + (pkin(4) * t387 - t655 + (t363 - t370) * t451 + (t362 - t605) * t449) * MDP(20) + (-t368 * t449 + t394 * t451 - t532 + t677) * MDP(21) + (t349 * qJ(5) - t350 * pkin(4) - t368 * t394 - t362 * t370 - g(1) * t581 - g(2) * t582 - g(3) * (-pkin(4) * t638 + t487) + t605 * t363) * MDP(22) + ((pkin(5) + t666) * t445 + t361 * t682 + t374 * t449 + t669) * MDP(23) + (t359 * t449 - t360 * t682 - t374 * t451 - t533 + t677 + t679) * MDP(24) + (t655 - t387 * t666 + (-t355 + t361) * t451 + (-t354 + t606) * t449) * MDP(25) + (t347 * qJ(5) - t346 * t666 - t354 * t361 - t359 * t374 - g(1) * (-pkin(5) * t416 + t581) - g(2) * (-pkin(5) * t414 + t582) - g(3) * (-t522 * t598 + t487) + t606 * t355) * MDP(26); (t534 - t650) * MDP(22) + (-t651 - t663 - t669) * MDP(26) + t672 * (-t445 + t647) + t681 * (-t682 ^ 2 - t444) - t684 * t537; (-t645 - t553) * MDP(23) + (t576 - t646) * MDP(24) + (-t444 - t667) * MDP(25) + (t354 * t451 - t355 * t449 + t584 + t675) * MDP(26) + (-MDP(23) * t680 - MDP(24) * t586) * t522;];
tau  = t1;
