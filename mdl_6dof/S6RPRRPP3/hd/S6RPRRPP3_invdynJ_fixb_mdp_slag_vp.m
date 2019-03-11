% Calculate vector of inverse dynamics joint torques for
% S6RPRRPP3
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
%   see S6RPRRPP3_convert_par2_MPV_fixb.m
% 
% Output:
% tau [6x1]
%   joint torques of inverse dynamics (contains inertial, gravitational Coriolis and centrifugal forces)

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 04:37
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tau = S6RPRRPP3_invdynJ_fixb_mdp_slag_vp(qJ, qJD, qJDD, g, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(3,1),zeros(9,1),zeros(26,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRRPP3_invdynJ_fixb_mdp_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RPRRPP3_invdynJ_fixb_mdp_slag_vp: qJD has to be [6x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [6 1]), ...
  'S6RPRRPP3_invdynJ_fixb_mdp_slag_vp: qJDD has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RPRRPP3_invdynJ_fixb_mdp_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S6RPRRPP3_invdynJ_fixb_mdp_slag_vp: pkin has to be [9x1] (double)');
assert(isreal(MDP) && all(size(MDP) == [26 1]), ...
  'S6RPRRPP3_invdynJ_fixb_mdp_slag_vp: MDP has to be [26x1] (double)'); 

%% Symbolic Calculation
% From invdyn_joint_fixb_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 04:37:13
% EndTime: 2019-03-09 04:37:20
% DurationCPUTime: 6.43s
% Computational Cost: add. (4217->571), mult. (8247->664), div. (0->0), fcn. (5100->10), ass. (0->233)
t515 = cos(qJ(3));
t498 = t515 * qJDD(1);
t512 = sin(qJ(3));
t601 = qJD(1) * qJD(3);
t435 = t512 * t601 + qJDD(4) - t498;
t650 = pkin(4) + qJ(6);
t580 = t650 * t435;
t508 = sin(pkin(9));
t489 = pkin(1) * t508 + pkin(7);
t457 = t489 * qJDD(1);
t681 = -qJD(2) * qJD(3) - t457;
t511 = sin(qJ(4));
t514 = cos(qJ(4));
t606 = qJD(4) * t514;
t609 = qJD(3) * t515;
t680 = t511 * t609 + t512 * t606;
t505 = qJ(1) + pkin(9);
t495 = cos(t505);
t636 = t495 * t512;
t494 = sin(t505);
t638 = t494 * t512;
t679 = g(1) * t636 + g(2) * t638;
t459 = t489 * qJD(1);
t418 = qJD(2) * t515 - t512 * t459;
t678 = qJD(3) * t418;
t419 = t512 * qJD(2) + t515 * t459;
t402 = qJD(3) * pkin(8) + t419;
t501 = t512 * pkin(8);
t503 = t515 * pkin(3);
t589 = -pkin(2) - t503;
t547 = t589 - t501;
t509 = cos(pkin(9));
t662 = pkin(1) * t509;
t536 = t547 - t662;
t407 = t536 * qJD(1);
t366 = t402 * t511 - t514 * t407;
t611 = qJD(3) * t511;
t614 = qJD(1) * t512;
t441 = t514 * t614 + t611;
t542 = pkin(5) * t441 + t366;
t603 = qJD(5) + t542;
t677 = MDP(20) - MDP(25);
t676 = MDP(21) + MDP(24);
t607 = qJD(4) * t512;
t675 = qJD(1) * t607 - qJDD(3);
t599 = qJDD(1) * t512;
t613 = qJD(1) * t515;
t383 = (qJD(3) * (qJD(4) + t613) + t599) * t511 + t675 * t514;
t428 = t435 * qJ(5);
t481 = -qJD(4) + t613;
t469 = qJD(5) * t481;
t673 = t469 - t428;
t641 = t441 * t481;
t672 = -t383 - t641;
t608 = qJD(4) * t511;
t671 = -qJD(5) * t511 - t419 + (-t511 * t613 + t608) * pkin(4);
t670 = g(1) * t495 + g(2) * t494;
t574 = MDP(17) - t677;
t669 = MDP(19) + MDP(23);
t668 = -pkin(5) * t383 + qJDD(6);
t597 = MDP(18) - MDP(21);
t666 = t511 * t574 + t514 * t597;
t605 = t514 * qJD(3);
t439 = t511 * t614 - t605;
t665 = t439 ^ 2;
t434 = t441 ^ 2;
t478 = t481 ^ 2;
t664 = 0.2e1 * t428;
t663 = pkin(5) + pkin(8);
t661 = pkin(4) * t435;
t579 = t515 * t601;
t382 = -qJD(4) * t605 + (-t579 - t599) * t514 + t675 * t511;
t660 = pkin(5) * t382;
t658 = pkin(5) * t439;
t657 = pkin(8) * t435;
t656 = g(1) * t494;
t653 = g(2) * t495;
t652 = g(3) * t512;
t651 = g(3) * t515;
t649 = pkin(8) * qJD(4);
t648 = qJ(5) * t383;
t647 = qJ(5) * t439;
t646 = qJ(5) * t514;
t367 = t514 * t402 + t511 * t407;
t360 = qJ(5) * t481 - t367;
t645 = t360 * t481;
t644 = t367 * t481;
t643 = t439 * t441;
t642 = t439 * t481;
t640 = t489 * t511;
t639 = t489 * t514;
t637 = t494 * t515;
t635 = t495 * t515;
t634 = t511 * t512;
t633 = t511 * t515;
t632 = t512 * t514;
t631 = t514 * t515;
t630 = qJDD(2) - g(3);
t587 = t515 * t605;
t629 = -t383 * t632 - t439 * t587;
t555 = qJ(6) * t511 - t646;
t541 = t555 * t515;
t628 = -qJD(1) * t541 + qJD(4) * t555 - qJD(6) * t514 + t671;
t561 = pkin(3) * t512 - pkin(8) * t515;
t445 = t561 * qJD(1);
t627 = t514 * t418 + t511 * t445;
t626 = t435 * t632 - t481 * t587;
t490 = -pkin(2) - t662;
t617 = t501 + t503;
t430 = t490 - t617;
t448 = t561 * qJD(3);
t625 = t430 * t606 + t511 * t448;
t624 = -qJ(5) * t606 + t613 * t646 + t671;
t444 = t489 * t631;
t623 = t511 * t430 + t444;
t596 = pkin(5) * t633;
t622 = -t663 * t608 - (qJ(5) * t512 - t596) * qJD(1) - t627;
t463 = t663 * t514;
t399 = t511 * t418;
t569 = -t445 * t514 + t399;
t595 = pkin(5) * t631;
t621 = qJD(4) * t463 - (-t512 * t650 + t595) * qJD(1) - t569;
t620 = t679 * t511;
t619 = t679 * t514;
t618 = pkin(4) * t634 - qJ(5) * t632;
t506 = t512 ^ 2;
t616 = -t515 ^ 2 + t506;
t615 = MDP(19) * t514;
t460 = qJD(1) * t490;
t610 = qJD(3) * t512;
t604 = -qJD(5) - t366;
t358 = t367 - t658;
t602 = -qJD(6) - t358;
t443 = t489 * t633;
t592 = -t382 * t634 + t680 * t441;
t591 = -t459 * t609 + t681 * t512;
t421 = t439 * t610;
t401 = -qJD(3) * pkin(3) - t418;
t531 = -qJ(5) * t441 + t401;
t356 = t439 * t650 + t531;
t586 = t356 * t608;
t585 = t356 * t606;
t584 = t489 * t608;
t583 = t511 * t607;
t581 = t441 * t610;
t577 = -t511 * qJ(5) - pkin(3);
t576 = -pkin(4) - t640;
t378 = -qJDD(3) * pkin(3) - qJDD(2) * t515 - t591;
t522 = qJ(5) * t382 - qJD(5) * t441 + t378;
t345 = qJD(6) * t439 + t383 * t650 + t522;
t575 = -t345 - t651;
t573 = -MDP(24) + t597;
t411 = t494 * t633 + t495 * t514;
t412 = t494 * t631 - t495 * t511;
t572 = -t411 * pkin(4) + qJ(5) * t412;
t413 = -t494 * t514 + t495 * t633;
t414 = t494 * t511 + t495 * t631;
t571 = -t413 * pkin(4) + qJ(5) * t414;
t570 = -qJ(5) + t639;
t568 = t430 * t514 - t443;
t377 = qJDD(3) * pkin(8) + qJDD(2) * t512 + t457 * t515 + t678;
t384 = qJD(1) * t448 + qJDD(1) * t536;
t565 = t511 * t377 - t514 * t384 + t402 * t606 + t407 * t608;
t564 = -t514 * t377 - t511 * t384 + t402 * t608 - t407 * t606;
t563 = t680 * pkin(4) + qJ(5) * t583 + t489 * t609;
t562 = pkin(4) * t631 + qJ(5) * t633 + t617;
t560 = g(1) * t411 - g(2) * t413;
t559 = g(1) * t412 - g(2) * t414;
t513 = sin(qJ(1));
t516 = cos(qJ(1));
t557 = g(1) * t513 - g(2) * t516;
t386 = qJ(5) * t515 - t623;
t556 = -qJD(4) * t444 - t430 * t608 + t448 * t514;
t548 = -qJDD(5) - t565;
t343 = qJD(6) * t481 - t548 - t580 - t660;
t346 = t564 + t673;
t344 = -t346 + t668;
t554 = t343 * t511 + t344 * t514;
t347 = -t548 - t661;
t553 = -t346 * t514 + t347 * t511;
t351 = t481 * t650 + t603;
t353 = qJD(6) - t360 - t658;
t552 = t351 * t514 - t353 * t511;
t551 = t351 * t511 + t353 * t514;
t359 = pkin(4) * t481 - t604;
t550 = t359 * t514 + t360 * t511;
t549 = t359 * t511 - t360 * t514;
t546 = pkin(4) * t514 - t577;
t545 = -qJ(6) * t634 - t618;
t544 = t481 * t649 - t651;
t348 = pkin(4) * t383 + t522;
t540 = -t348 + t544;
t539 = -pkin(1) * t513 - t412 * pkin(4) + t495 * pkin(7) - qJ(5) * t411;
t538 = t514 * t650 - t577;
t537 = -qJD(1) * t460 + t670;
t535 = -t401 * t481 - t657;
t365 = pkin(4) * t439 + t531;
t534 = t365 * t481 + t657;
t532 = t516 * pkin(1) + t495 * pkin(2) + pkin(3) * t635 + t414 * pkin(4) + t494 * pkin(7) + pkin(8) * t636 + qJ(5) * t413;
t517 = qJD(3) ^ 2;
t530 = 0.2e1 * qJDD(1) * t490 + t489 * t517 + t653;
t528 = 0.2e1 * qJD(3) * t460 - qJDD(3) * t489;
t527 = g(1) * t413 + g(2) * t411 + g(3) * t634 - t565;
t525 = -qJDD(5) + t527;
t362 = -t382 - t642;
t524 = g(1) * t414 + g(2) * t412 + g(3) * t632 + t564;
t523 = t365 * t441 - t525;
t521 = t356 * t441 - t525 - t660;
t520 = -t356 * t439 - t524 + t668;
t502 = t515 * pkin(4);
t468 = g(1) * t638;
t466 = pkin(8) * t635;
t464 = pkin(8) * t637;
t462 = t663 * t511;
t461 = t512 * t489;
t456 = qJDD(3) * t515 - t512 * t517;
t455 = qJDD(3) * t512 + t515 * t517;
t396 = t461 + t618;
t389 = pkin(4) * t441 + t647;
t388 = t461 - t545;
t387 = t502 - t568;
t379 = -pkin(5) * t634 - t386;
t374 = -pkin(4) * t614 + t569;
t373 = -qJ(5) * t614 - t627;
t371 = t441 * t650 + t647;
t368 = qJ(6) * t515 + t443 + t502 + (pkin(5) * t512 - t430) * t514;
t363 = (-qJ(5) * t609 - qJD(5) * t512) * t514 + t563;
t355 = t576 * t610 - t556;
t354 = (qJD(5) + t584) * t515 + t570 * t610 - t625;
t352 = qJD(3) * t541 + (qJD(6) * t511 + (qJ(6) * qJD(4) - qJD(5)) * t514) * t512 + t563;
t350 = -qJD(5) * t515 + (-pkin(5) * t632 - t443) * qJD(4) + (-t512 * t570 - t596) * qJD(3) + t625;
t349 = -pkin(5) * t583 + qJD(6) * t515 + (t595 + (-qJ(6) + t576) * t512) * qJD(3) - t556;
t1 = [(t382 * t515 + t481 * t583 + t581 + t626) * MDP(14) + (t349 * t441 - t350 * t439 - t368 * t382 - t379 * t383 + t468 + t552 * t609 + (-qJD(4) * t551 + t343 * t514 - t344 * t511 - t653) * t512) * MDP(23) + (t354 * t439 + t355 * t441 - t382 * t387 + t383 * t386 + t468 + t550 * t609 + (-qJD(4) * t549 + t346 * t511 + t347 * t514 - t653) * t512) * MDP(19) + (t383 * t515 - t435 * t634 + t680 * t481 - t421) * MDP(15) + (t345 * t388 + t356 * t352 + t343 * t368 + t351 * t349 + t344 * t379 + t353 * t350 - g(1) * (-qJ(6) * t412 + t539) - g(2) * (pkin(5) * t636 + qJ(6) * t414 + t532) - (-t512 * t663 + t589) * t656) * MDP(26) + (t439 * t583 - t592 + t629) * MDP(13) + (t512 * t530 + t515 * t528 - t468) * MDP(11) + t557 * MDP(2) + (-g(1) * t539 - g(2) * t532 + t346 * t386 + t347 * t387 + t348 * t396 + t360 * t354 + t359 * t355 + t365 * t363 - t547 * t656) * MDP(22) + (t528 * t512 + (-t530 + t656) * t515) * MDP(10) + (-t350 * t481 - t352 * t441 + t379 * t435 + t382 * t388 + (-t356 * t605 - t344) * t515 + (qJD(3) * t353 - t345 * t514 + t586) * t512 + t560) * MDP(24) + (t625 * t481 - t623 * t435 + (-t481 * t584 + (t401 * t514 + t441 * t489) * qJD(3) - t564) * t515 + (-t401 * t608 + t378 * t514 - t489 * t382 + (-t481 * t639 - t367) * qJD(3)) * t512 - t560) * MDP(18) + (-t556 * t481 + t568 * t435 + ((t401 * t511 + t439 * t489) * qJD(3) + t565) * t515 + (t401 * t606 + t378 * t511 + t489 * t383 + (-t481 * t640 - t366) * qJD(3)) * t512 + t559) * MDP(17) + (-t382 * t632 + (-t583 + t587) * t441) * MDP(12) + qJDD(1) * MDP(1) + (t354 * t481 - t363 * t441 + t382 * t396 - t386 * t435 + (-t365 * t605 + t346) * t515 + (-qJD(3) * t360 - t348 * t514 + t365 * t608) * t512 + t560) * MDP(21) + (-t435 * t515 - t481 * t610) * MDP(16) + (-t355 * t481 - t363 * t439 - t383 * t396 + t387 * t435 + (-t365 * t611 - t347) * t515 + (qJD(3) * t359 - t348 * t511 - t365 * t606) * t512 - t559) * MDP(20) + (t349 * t481 + t352 * t439 - t368 * t435 + t383 * t388 + (t356 * t611 + t343) * t515 + (-qJD(3) * t351 + t345 * t511 + t585) * t512 + t559) * MDP(25) + (qJDD(1) * t506 + 0.2e1 * t512 * t579) * MDP(5) + (t557 + (t508 ^ 2 + t509 ^ 2) * qJDD(1) * pkin(1)) * pkin(1) * MDP(4) + 0.2e1 * (t498 * t512 - t601 * t616) * MDP(6) + t455 * MDP(7) + t456 * MDP(8) + (g(1) * t516 + g(2) * t513) * MDP(3); qJDD(2) * MDP(4) + t456 * MDP(10) - t455 * MDP(11) + t592 * MDP(19) + (t592 + t629) * MDP(23) + t626 * MDP(24) + (-MDP(4) - MDP(22) - MDP(26)) * g(3) + (-t348 * MDP(22) - t345 * MDP(26) - t574 * t383 + t573 * t382 + (t549 * MDP(22) + t551 * MDP(26) - t439 * t615 + t481 * t666) * qJD(3)) * t515 + (-t383 * t615 + (qJD(3) * t365 + t553) * MDP(22) + (qJD(3) * t356 + t554) * MDP(26) - t666 * t435 + (t550 * MDP(22) + t552 * MDP(26) + t669 * t511 * t439 + (-t511 * t573 + t514 * t574) * t481) * qJD(4)) * t512 + t574 * t421 + (MDP(18) - t676) * t581; MDP(7) * t599 + MDP(8) * t498 + qJDD(3) * MDP(9) + (qJD(3) * t419 + t512 * t537 + t515 * t630 + t591) * MDP(10) + (t678 + (qJD(3) * t459 - t630) * t512 + (t537 + t681) * t515) * MDP(11) + (-t382 * t511 - t514 * t641) * MDP(12) + ((-t382 + t642) * t514 + (-t383 + t641) * t511) * MDP(13) + (-t481 * t606 + t435 * t511 + (-t441 * t512 + t481 * t631) * qJD(1)) * MDP(14) + (t481 * t608 + t435 * t514 + (t439 * t512 - t481 * t633) * qJD(1)) * MDP(15) + t481 * MDP(16) * t614 + (t366 * t614 - pkin(3) * t383 - t399 * t481 - t419 * t439 + (-t651 - t378 + (t445 + t649) * t481) * t514 + t535 * t511 + t619) * MDP(17) + (pkin(3) * t382 - t627 * t481 + t367 * t614 - t419 * t441 + t535 * t514 + (t378 - t544) * t511 - t620) * MDP(18) + (-t652 - t373 * t439 - t374 * t441 - t670 * t515 + (-t346 - t481 * t359 + (qJD(4) * t441 - t383) * pkin(8)) * t514 + (t347 - t645 + (qJD(4) * t439 - t382) * pkin(8)) * t511) * MDP(19) + (-t359 * t614 + t374 * t481 + t383 * t546 - t439 * t624 + t511 * t534 - t514 * t540 - t619) * MDP(20) + (t360 * t614 - t373 * t481 - t382 * t546 - t441 * t624 + t511 * t540 + t514 * t534 + t620) * MDP(21) + (-t360 * t373 - t359 * t374 - g(1) * t466 - g(2) * t464 - g(3) * t562 + t624 * t365 + (qJD(4) * t550 + t553) * pkin(8) + (t670 * t512 - t348) * t546) * MDP(22) + (-t652 - t382 * t462 - t383 * t463 + t621 * t441 - t622 * t439 + t552 * qJD(4) + (-qJD(1) * t552 - t670) * t515 + t554) * MDP(23) + (-t585 - t382 * t538 + t435 * t463 + t575 * t511 - t622 * t481 - t628 * t441 + (-t353 * t512 + t356 * t631) * qJD(1) + t620) * MDP(24) + (t586 - t383 * t538 - t435 * t462 + t575 * t514 + t621 * t481 + t628 * t439 + (t351 * t512 - t356 * t633) * qJD(1) + t619) * MDP(25) + (-t345 * t538 + t343 * t462 + t344 * t463 - g(1) * (pkin(5) * t635 + t466) - g(2) * (pkin(5) * t637 + t464) - g(3) * (qJ(6) * t631 + t562) + t628 * t356 + t622 * t353 + t621 * t351 + (-g(3) * pkin(5) + t538 * t670) * t512) * MDP(26) + (-MDP(5) * t512 * t515 + MDP(6) * t616) * qJD(1) ^ 2; MDP(12) * t643 + (t434 - t665) * MDP(13) + t362 * MDP(14) + t672 * MDP(15) + t435 * MDP(16) + (-t401 * t441 + t527 - t644) * MDP(17) + (t366 * t481 + t401 * t439 + t524) * MDP(18) + (pkin(4) * t382 - t648 + (-t360 - t367) * t441 + (t359 + t604) * t439) * MDP(19) + (t389 * t439 + t523 + t644 - 0.2e1 * t661) * MDP(20) + (-t365 * t439 + t389 * t441 + t481 * t604 - t469 - t524 + t664) * MDP(21) + (-t347 * pkin(4) - g(1) * t571 - g(2) * t572 + g(3) * t618 - t346 * qJ(5) - t359 * t367 + t360 * t604 - t365 * t389) * MDP(22) + (-t648 + t382 * t650 + (t353 + t602) * t441 + (t351 - t603) * t439) * MDP(23) + (t371 * t441 - t481 * t542 - 0.2e1 * t469 + t520 + t664) * MDP(24) + (-t371 * t439 + (-0.2e1 * qJD(6) - t358) * t481 + 0.2e1 * t580 - t521) * MDP(25) + (-t343 * t650 + t344 * qJ(5) - t356 * t371 - g(1) * (-qJ(6) * t413 + t571) - g(2) * (-qJ(6) * t411 + t572) - g(3) * t545 + t603 * t353 + t602 * t351) * MDP(26); (t523 - t645 - t661) * MDP(22) + ((qJD(6) + t353) * t481 - t580 + t521) * MDP(26) + t677 * (t435 - t643) + t669 * t362 + t676 * (-t434 - t478); (t435 + t643) * MDP(24) + (-t478 - t665) * MDP(25) + (-t351 * t481 + t520 - t673) * MDP(26) + t672 * MDP(23);];
tau  = t1;
