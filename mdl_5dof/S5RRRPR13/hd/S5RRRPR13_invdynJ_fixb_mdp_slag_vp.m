% Calculate vector of inverse dynamics joint torques for
% S5RRRPR13
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
%   pkin=[a2,a3,a4,a5,alpha2,d1,d2,d3,d5]';
% MDP [28x1]
%   Minimal dynamic parameter vector (fixed base model)
%   see S5RRRPR13_convert_par2_MPV_fixb.m
% 
% Output:
% tau [5x1]
%   joint torques of inverse dynamics (contains inertial, gravitational Coriolis and centrifugal forces)

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 21:48
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tau = S5RRRPR13_invdynJ_fixb_mdp_slag_vp(qJ, qJD, qJDD, g, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(5,1),zeros(3,1),zeros(9,1),zeros(28,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRRPR13_invdynJ_fixb_mdp_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RRRPR13_invdynJ_fixb_mdp_slag_vp: qJD has to be [5x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [5 1]), ...
  'S5RRRPR13_invdynJ_fixb_mdp_slag_vp: qJDD has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RRRPR13_invdynJ_fixb_mdp_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5RRRPR13_invdynJ_fixb_mdp_slag_vp: pkin has to be [9x1] (double)');
assert(isreal(MDP) && all(size(MDP) == [28 1]), ...
  'S5RRRPR13_invdynJ_fixb_mdp_slag_vp: MDP has to be [28x1] (double)'); 

%% Symbolic Calculation
% From invdyn_joint_fixb_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 21:46:51
% EndTime: 2019-12-31 21:47:05
% DurationCPUTime: 8.78s
% Computational Cost: add. (4434->565), mult. (11076->753), div. (0->0), fcn. (8574->10), ass. (0->239)
t512 = sin(qJ(3));
t516 = cos(qJ(3));
t659 = cos(pkin(5));
t592 = t659 * qJD(1);
t555 = t592 + qJD(2);
t513 = sin(qJ(2));
t510 = sin(pkin(5));
t625 = qJD(1) * t510;
t604 = t513 * t625;
t439 = t512 * t555 + t516 * t604;
t431 = qJD(5) + t439;
t615 = qJDD(1) * t513;
t597 = t510 * t615;
t517 = cos(qJ(2));
t616 = qJD(1) * qJD(2);
t598 = t517 * t616;
t677 = t510 * t598 + t597;
t603 = t517 * t625;
t484 = -qJD(3) + t603;
t577 = pkin(1) * t592;
t456 = pkin(7) * t603 + t513 * t577;
t422 = pkin(8) * t555 + t456;
t552 = -pkin(2) * t517 - pkin(8) * t513 - pkin(1);
t450 = t552 * t510;
t430 = qJD(1) * t450;
t379 = t422 * t512 - t516 * t430;
t618 = -qJD(4) - t379;
t514 = sin(qJ(1));
t667 = cos(qJ(1));
t572 = t659 * t667;
t463 = t513 * t572 + t514 * t517;
t605 = t510 * t667;
t411 = t463 * t512 + t516 * t605;
t462 = t513 * t514 - t517 * t572;
t511 = sin(qJ(5));
t515 = cos(qJ(5));
t676 = t411 * t511 + t462 * t515;
t675 = t411 * t515 - t462 * t511;
t540 = qJD(3) * t555;
t585 = t659 * qJDD(1);
t550 = t585 + qJDD(2);
t621 = qJD(3) * t516;
t623 = qJD(2) * t517;
t382 = t510 * (qJD(1) * (t512 * t623 + t513 * t621) + t512 * t615) + t512 * t540 - t516 * t550;
t593 = t517 * t659;
t641 = t510 * t513;
t547 = pkin(1) * t593 - pkin(7) * t641;
t674 = (qJDD(2) + 0.2e1 * t585) * t510;
t579 = t512 * t603;
t622 = qJD(3) * t512;
t673 = -qJD(4) * t512 - t456 + (-t579 + t622) * pkin(3);
t617 = pkin(4) * t439 - t618;
t464 = t513 * t667 + t514 * t593;
t638 = t510 * t517;
t535 = g(1) * t464 + g(2) * t462 - g(3) * t638;
t380 = t516 * t422 + t512 * t430;
t373 = qJ(4) * t484 - t380;
t437 = t512 * t604 - t516 * t555;
t665 = pkin(4) * t437;
t361 = -t373 - t665;
t610 = t512 * t641;
t578 = qJD(3) * t610;
t381 = qJD(1) * t578 - t512 * t550 + (-t540 - t677) * t516;
t378 = -qJDD(5) + t381;
t669 = pkin(3) + pkin(9);
t672 = t669 * t378 + (t361 - t380 + t665) * t431;
t594 = t513 * t659;
t627 = pkin(1) * t594 + pkin(7) * t638;
t449 = pkin(8) * t659 + t627;
t570 = pkin(2) * t513 - pkin(8) * t517;
t545 = t570 * qJD(2);
t455 = t510 * t545;
t457 = t547 * qJD(2);
t542 = t449 * t622 - t450 * t621 - t512 * t455 - t516 * t457;
t624 = qJD(2) * t513;
t356 = -t510 * (qJ(4) * t624 - qJD(4) * t517) + t542;
t670 = t439 ^ 2;
t519 = qJD(1) ^ 2;
t668 = pkin(4) + pkin(8);
t614 = qJDD(1) * t517;
t500 = t510 * t614;
t599 = t513 * t616;
t575 = t510 * t599;
t451 = qJDD(3) - t500 + t575;
t666 = pkin(3) * t451;
t664 = pkin(8) * t451;
t460 = -t516 * t659 + t610;
t661 = t460 * pkin(3);
t660 = pkin(8) * qJD(3);
t658 = qJ(4) * t437;
t657 = qJ(4) * t516;
t619 = qJD(5) * t515;
t608 = t511 * t382 + t437 * t619 + t515 * t451;
t620 = qJD(5) * t511;
t354 = t484 * t620 + t608;
t656 = t354 * t515;
t655 = t373 * t484;
t654 = t380 * t484;
t643 = t484 * t511;
t400 = -t515 * t437 - t643;
t653 = t400 * t431;
t652 = t400 * t484;
t402 = t437 * t511 - t484 * t515;
t651 = t402 * t431;
t650 = t402 * t484;
t649 = t437 * t439;
t648 = t437 * t484;
t647 = t439 * t484;
t447 = t451 * qJ(4);
t639 = t510 * t516;
t461 = t512 * t659 + t513 * t639;
t646 = t461 * qJ(4);
t507 = t510 ^ 2;
t642 = t507 * t519;
t640 = t510 * t514;
t637 = t511 * t378;
t636 = t511 * t512;
t635 = t512 * t515;
t634 = t512 * t517;
t375 = t515 * t378;
t633 = t515 * t517;
t632 = t516 * t517;
t631 = t516 * t449 + t512 * t450;
t453 = -pkin(7) * t604 + t517 * t577;
t454 = t570 * t625;
t630 = t516 * t453 + t512 * t454;
t629 = -qJ(4) * t621 + t603 * t657 + t673;
t628 = -t668 * t622 - (-pkin(4) * t634 + qJ(4) * t513) * t625 - t630;
t576 = pkin(1) * qJD(2) * t659;
t601 = t510 * t623;
t458 = pkin(7) * t601 + t513 * t576;
t508 = t513 ^ 2;
t626 = -t517 ^ 2 + t508;
t486 = t668 * t516;
t611 = t517 * t642;
t609 = t510 * t633;
t551 = qJD(1) * t576;
t573 = pkin(1) * t585;
t606 = pkin(7) * t500 + t513 * t573 + t517 * t551;
t602 = t510 * t624;
t600 = 0.2e1 * pkin(1) * t507;
t595 = -qJ(4) * t512 - pkin(2);
t533 = -pkin(7) * t575 + t606;
t395 = pkin(8) * t550 + t533;
t398 = (qJD(1) * t545 + qJDD(1) * t552) * t510;
t582 = t512 * t395 - t516 * t398 + t422 * t621 + t430 * t622;
t553 = qJDD(4) + t582;
t343 = -pkin(4) * t381 - t451 * t669 + t553;
t581 = pkin(7) * t677 + t513 * t551 - t517 * t573;
t396 = -pkin(2) * t550 + t581;
t520 = t381 * qJ(4) - t439 * qJD(4) + t396;
t345 = t382 * t669 + t520;
t591 = t515 * t343 - t511 * t345;
t590 = -t515 * t382 + t451 * t511;
t589 = -t512 * t449 + t450 * t516;
t441 = t512 * t453;
t588 = -t454 * t516 + t441;
t412 = t463 * t516 - t512 * t605;
t587 = t431 * t511;
t586 = t431 * t515;
t583 = -t516 * t395 - t512 * t398 + t422 * t622 - t430 * t621;
t571 = t510 * t519 * t659;
t465 = -t514 * t594 + t517 * t667;
t415 = t465 * t512 - t514 * t639;
t569 = g(1) * t411 - g(2) * t415;
t416 = t465 * t516 + t512 * t640;
t568 = -g(1) * t412 + g(2) * t416;
t567 = -g(1) * t462 + g(2) * t464;
t566 = -g(1) * t465 - g(2) * t463;
t385 = pkin(3) * t638 - t589;
t428 = t511 * t604 - t515 * t579;
t565 = t515 * t622 + t428;
t536 = t510 * (t511 * t634 + t513 * t515);
t429 = qJD(1) * t536;
t564 = t511 * t622 - t429;
t466 = -t516 * t669 + t595;
t563 = qJD(5) * t466 + (pkin(4) * t632 - t513 * t669) * t625 + t588 - qJD(3) * t486;
t485 = t668 * t512;
t562 = -qJD(5) * t485 - t673 + t484 * (pkin(9) * t512 - t657);
t560 = t511 * t343 + t515 * t345;
t357 = t484 * t669 + t617;
t421 = -pkin(2) * t555 - t453;
t521 = -t439 * qJ(4) + t421;
t360 = t437 * t669 + t521;
t346 = t357 * t515 - t360 * t511;
t347 = t357 * t511 + t360 * t515;
t364 = pkin(4) * t461 + pkin(9) * t638 + t385;
t448 = -pkin(2) * t659 - t547;
t526 = t448 - t646;
t368 = t460 * t669 + t526;
t558 = t364 * t515 - t368 * t511;
t557 = t364 * t511 + t368 * t515;
t554 = 0.2e1 * t592 + qJD(2);
t549 = pkin(3) * t516 - t595;
t384 = qJ(4) * t638 - t631;
t548 = -t449 * t621 - t450 * t622 + t455 * t516 - t512 * t457;
t471 = qJD(4) * t484;
t348 = -t447 + t471 + t583;
t409 = t460 * t515 + t511 * t638;
t544 = t484 * t516;
t541 = -g(1) * t416 - g(2) * t412 - g(3) * t461;
t539 = t550 * MDP(8);
t408 = -t578 + (qJD(3) * t659 + t601) * t516;
t538 = -qJ(4) * t408 - qJD(4) * t461 + t458;
t534 = -g(3) * t641 + t566;
t532 = -t421 * t484 - t664;
t367 = t437 * pkin(3) + t521;
t531 = t367 * t484 + t664;
t529 = g(1) * t415 + g(2) * t411 + g(3) * t460 - t582;
t528 = t541 - t583;
t527 = t484 * t660 + t535;
t349 = t382 * pkin(3) + t520;
t525 = -t349 + t527;
t524 = -t381 - t648;
t344 = -pkin(4) * t382 - t348;
t523 = t344 + (t431 * t669 + t658) * t431 + t541;
t522 = t367 * t439 + qJDD(4) - t529;
t410 = -t460 * t511 + t609;
t407 = qJD(3) * t461 + t512 * t601;
t397 = pkin(3) * t439 + t658;
t390 = t415 * t511 + t464 * t515;
t389 = t415 * t515 - t464 * t511;
t388 = -pkin(3) * t604 + t588;
t386 = -qJ(4) * t604 - t630;
t383 = t526 + t661;
t371 = pkin(3) * t484 - t618;
t369 = -pkin(4) * t460 - t384;
t366 = qJD(5) * t409 + t407 * t511 + t515 * t602;
t365 = -t407 * t515 - qJD(5) * t609 + (qJD(5) * t460 + t602) * t511;
t359 = pkin(3) * t407 + t538;
t358 = -pkin(3) * t602 - t548;
t355 = qJD(5) * t402 + t590;
t353 = t407 * t669 + t538;
t352 = -pkin(4) * t407 - t356;
t351 = pkin(4) * t408 - t602 * t669 - t548;
t350 = t553 - t666;
t341 = -t347 * qJD(5) + t591;
t340 = t346 * qJD(5) + t560;
t1 = [(t354 * t409 + t355 * t410 - t365 * t402 - t366 * t400) * MDP(23) + (-t354 * t410 + t366 * t402) * MDP(22) + (t349 * t383 + t367 * t359 + t348 * t384 + t373 * t356 + t350 * t385 + t371 * t358 - g(1) * (-t514 * pkin(1) - t463 * pkin(2) - pkin(3) * t412 + pkin(7) * t605 - t462 * pkin(8) - qJ(4) * t411) - g(2) * (pkin(1) * t667 + t465 * pkin(2) + t416 * pkin(3) + pkin(7) * t640 + t464 * pkin(8) + t415 * qJ(4))) * MDP(21) + (g(1) * t514 - g(2) * t667) * MDP(2) + (g(1) * t667 + g(2) * t514) * MDP(3) + t659 * t539 + (-t457 * t555 - t627 * t550 - t533 * t659 + (-t598 - t615) * t600 + t567) * MDP(10) + (-t458 * t555 + t547 * t550 - t581 * t659 + g(1) * t463 - g(2) * t465 + (-t599 + t614) * t600) * MDP(9) + ((qJDD(1) * t508 + 0.2e1 * t513 * t598) * MDP(4) + 0.2e1 * (t513 * t614 - t616 * t626) * MDP(5)) * t507 + (-t542 * t484 - t631 * t451 + t458 * t439 - t448 * t381 + t396 * t461 + t421 * t408 + (-t380 * t624 - t517 * t583) * t510 - t569) * MDP(17) + (-t548 * t484 + t589 * t451 + t458 * t437 + t448 * t382 + t396 * t460 + t421 * t407 + (-t379 * t624 + t517 * t582) * t510 - t568) * MDP(16) + (t407 * t484 - t451 * t460 + (t382 * t517 - t437 * t624) * t510) * MDP(14) + (-t349 * t461 + t356 * t484 - t359 * t439 - t367 * t408 + t381 * t383 - t384 * t451 + (t348 * t517 - t373 * t624) * t510 + t569) * MDP(20) + (-t408 * t484 + t451 * t461 + (t381 * t517 + t439 * t624) * t510) * MDP(13) + (-t349 * t460 - t358 * t484 - t359 * t437 - t367 * t407 - t382 * t383 + t385 * t451 + (-t350 * t517 + t371 * t624) * t510 + t568) * MDP(19) + (-t381 * t461 + t408 * t439) * MDP(11) + (t354 * t461 + t366 * t431 + t378 * t410 + t402 * t408) * MDP(24) + (-t378 * t461 + t408 * t431) * MDP(26) + (-t355 * t461 - t365 * t431 - t378 * t409 - t400 * t408) * MDP(25) + (t381 * t460 - t382 * t461 - t407 * t439 - t408 * t437) * MDP(12) + (-(qJD(5) * t558 + t351 * t511 + t353 * t515) * t431 + t557 * t378 - t340 * t461 - t347 * t408 + t352 * t402 + t369 * t354 - t344 * t410 + t361 * t366 + g(1) * t675 - g(2) * t389) * MDP(28) + ((-qJD(5) * t557 + t351 * t515 - t353 * t511) * t431 - t558 * t378 + t341 * t461 + t346 * t408 + t352 * t400 + t369 * t355 - t344 * t409 + t361 * t365 + g(1) * t676 - g(2) * t390) * MDP(27) + (t348 * t460 + t350 * t461 + t356 * t437 + t358 * t439 + t371 * t408 + t373 * t407 - t381 * t385 + t382 * t384 - t567) * MDP(18) + (t513 * t674 + t554 * t601) * MDP(6) + (t517 * t674 - t554 * t602) * MDP(7) + (-t451 * t517 - t484 * t624) * t510 * MDP(15) + qJDD(1) * MDP(1); (-t381 * t512 - t439 * t544) * MDP(11) + (-t378 * t512 - t431 * t544) * MDP(26) + (-t371 * t604 + t382 * t549 + t388 * t484 - t437 * t629 + t512 * t531 - t516 * t525) * MDP(19) + (t373 * t604 - t381 * t549 - t386 * t484 - t439 * t629 + t512 * t525 + t516 * t531) * MDP(20) + (t379 * t604 - pkin(2) * t382 - t456 * t437 - t441 * t484 + t532 * t512 + (-t396 + (t454 + t660) * t484 + t535) * t516) * MDP(16) + (t400 * t429 + t402 * t428 + (-t400 * t511 + t402 * t515) * t622 + (-t656 + t355 * t511 + (t400 * t515 + t402 * t511) * qJD(5)) * t516) * MDP(23) + (-t355 * t512 + t565 * t431 + (t431 * t620 + t375 + t652) * t516) * MDP(25) + (-t386 * t437 - t388 * t439 + (-t348 - t484 * t371 + (qJD(3) * t439 - t382) * pkin(8)) * t516 + (t350 - t655 + (qJD(3) * t437 - t381) * pkin(8)) * t512 + t534) * MDP(18) + ((-t381 + t648) * t516 + (-t382 + t647) * t512) * MDP(12) + (t354 * t512 + t564 * t431 + (-t431 * t619 + t637 - t650) * t516) * MDP(24) + (pkin(1) * t513 * t642 + t456 * t555 + t535 - t581) * MDP(9) + (pkin(1) * t611 + t453 * t555 + (pkin(7) * t616 + g(3)) * t641 - t566 - t606) * MDP(10) + ((t466 * t515 + t485 * t511) * t378 - t340 * t512 + t486 * t354 - g(1) * (-t464 * t635 - t465 * t511) - g(2) * (-t462 * t635 - t463 * t511) - g(3) * (-t511 * t513 + t512 * t633) * t510 + (t511 * t563 + t515 * t562) * t431 + t628 * t402 + t564 * t361 + (-t344 * t511 + t347 * t484 - t361 * t619) * t516) * MDP(28) + (-(-t466 * t511 + t485 * t515) * t378 + t341 * t512 + t486 * t355 - g(1) * (-t464 * t636 + t465 * t515) - g(2) * (-t462 * t636 + t463 * t515) - g(3) * t536 + (t511 * t562 - t515 * t563) * t431 + t628 * t400 - t565 * t361 + (t344 * t515 - t346 * t484 - t361 * t620) * t516) * MDP(27) + (t484 * t622 + t451 * t516 + (t437 * t513 - t484 * t634) * t625) * MDP(14) + (pkin(2) * t381 - t630 * t484 + t380 * t604 - t456 * t439 + t532 * t516 + (t396 - t527) * t512) * MDP(17) + (-t484 * t621 + t451 * t512 + (-t439 * t513 + t484 * t632) * t625) * MDP(13) + (-t354 * t511 * t516 + (-t516 * t619 + t564) * t402) * MDP(22) + (-t517 * t571 + t597) * MDP(6) + t484 * MDP(15) * t604 + (t513 * t571 + t500) * MDP(7) + (-t371 * t388 - t373 * t386 + t629 * t367 + (-t348 * t516 + t350 * t512 + (t371 * t516 + t373 * t512) * qJD(3) + t534) * pkin(8) + (-t349 + t535) * t549) * MDP(21) - t513 * MDP(4) * t611 + t626 * MDP(5) * t642 + t539; MDP(11) * t649 + (-t437 ^ 2 + t670) * MDP(12) + t524 * MDP(13) + (-t382 - t647) * MDP(14) + t451 * MDP(15) + (-t421 * t439 + t529 - t654) * MDP(16) + (t379 * t484 + t421 * t437 - t528) * MDP(17) + (pkin(3) * t381 - qJ(4) * t382 + (-t373 - t380) * t439 + (t371 + t618) * t437) * MDP(18) + (t397 * t437 + t522 + t654 - 0.2e1 * t666) * MDP(19) + (-t367 * t437 + t397 * t439 + t484 * t618 + 0.2e1 * t447 - t471 + t528) * MDP(20) + (-t348 * qJ(4) - t350 * pkin(3) - t367 * t397 - t371 * t380 - g(1) * (-pkin(3) * t415 + qJ(4) * t416) - g(2) * (-pkin(3) * t411 + qJ(4) * t412) - g(3) * (t646 - t661) + t618 * t373) * MDP(21) + (-t402 * t587 + t656) * MDP(22) + ((-t355 - t651) * t515 + (-t354 + t653) * t511) * MDP(23) + (t402 * t437 - t431 * t587 - t375) * MDP(24) + (-t400 * t437 - t431 * t586 + t637) * MDP(25) + t431 * t437 * MDP(26) + (qJ(4) * t355 + t346 * t437 + t617 * t400 + t523 * t511 + t515 * t672) * MDP(27) + (qJ(4) * t354 - t347 * t437 + t617 * t402 - t511 * t672 + t523 * t515) * MDP(28); t524 * MDP(18) + (t451 - t649) * MDP(19) + (-t484 ^ 2 - t670) * MDP(20) + (t522 - t655 - t666) * MDP(21) + (-t375 + t652) * MDP(27) + (t637 + t650) * MDP(28) + (-MDP(27) * t587 - MDP(28) * t586) * t431; t402 * t400 * MDP(22) + (-t400 ^ 2 + t402 ^ 2) * MDP(23) + (t608 + t653) * MDP(24) + (-t590 + t651) * MDP(25) - t378 * MDP(26) + (-g(1) * t389 - g(2) * t675 - g(3) * t409 + t347 * t431 - t361 * t402 + t591) * MDP(27) + (g(1) * t390 + g(2) * t676 - g(3) * t410 + t346 * t431 + t361 * t400 - t560) * MDP(28) + (MDP(24) * t643 - MDP(25) * t402 - MDP(27) * t347 - MDP(28) * t346) * qJD(5);];
tau = t1;
