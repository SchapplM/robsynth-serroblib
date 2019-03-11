% Calculate minimal parameter regressor of Coriolis joint torque vector for
% S6RPRRPR6
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [11x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d3,d4,d6,theta2,theta5]';
% MDP [30x1]
%   Minimal dynamic parameter vector (fixed base model)
%   see S6RPRRPR6_convert_par2_MPV_fixb.m
% 
% Output:
% tauc [6x1]
%   joint torques required to compensate Coriolis and centrifugal load

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 05:18
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauc = S6RPRRPR6_coriolisvecJ_fixb_mdp_slag_vp(qJ, qJD, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(11,1),zeros(30,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRRPR6_coriolisvecJ_fixb_mdp_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RPRRPR6_coriolisvecJ_fixb_mdp_slag_vp: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6RPRRPR6_coriolisvecJ_fixb_mdp_slag_vp: pkin has to be [11x1] (double)');
assert(isreal(MDP) && all(size(MDP) == [30 1]), ...
  'S6RPRRPR6_coriolisvecJ_fixb_mdp_slag_vp: MDP has to be [30x1] (double)'); 

%% Symbolic Calculation
% From coriolisvec_joint_fixb_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 05:17:55
% EndTime: 2019-03-09 05:18:07
% DurationCPUTime: 7.70s
% Computational Cost: add. (7708->455), mult. (20062->608), div. (0->0), fcn. (15875->10), ass. (0->200)
t545 = cos(pkin(10));
t551 = cos(qJ(3));
t594 = qJD(1) * t551;
t534 = t545 * t594;
t543 = sin(pkin(10));
t548 = sin(qJ(3));
t595 = qJD(1) * t548;
t579 = t543 * t595;
t507 = t534 - t579;
t498 = qJD(4) - t507;
t549 = cos(qJ(6));
t518 = t543 * t551 + t545 * t548;
t509 = t518 * qJD(1);
t547 = sin(qJ(4));
t550 = cos(qJ(4));
t587 = t550 * qJD(3);
t477 = t509 * t547 - t587;
t479 = qJD(3) * t547 + t509 * t550;
t542 = sin(pkin(11));
t544 = cos(pkin(11));
t568 = -t477 * t544 - t479 * t542;
t607 = t549 * t568;
t424 = t477 * t542 - t479 * t544;
t546 = sin(qJ(6));
t628 = t424 * t546;
t379 = t607 + t628;
t492 = qJD(6) + t498;
t630 = t379 * t492;
t531 = qJD(3) * t534;
t493 = -qJD(3) * t579 + t531;
t590 = qJD(4) * t547;
t433 = qJD(4) * t587 + t550 * t493 - t509 * t590;
t511 = t518 * qJD(3);
t494 = qJD(1) * t511;
t537 = -pkin(2) * t545 - pkin(1);
t524 = qJD(1) * t537 + qJD(2);
t438 = -pkin(3) * t507 - pkin(8) * t509 + t524;
t635 = pkin(7) + qJ(2);
t525 = t635 * t543;
t522 = qJD(1) * t525;
t526 = t635 * t545;
t523 = qJD(1) * t526;
t470 = -t548 * t522 + t551 * t523;
t463 = qJD(3) * pkin(8) + t470;
t408 = t438 * t547 + t463 * t550;
t516 = t543 * t548 - t551 * t545;
t554 = t516 * qJD(2);
t639 = -t522 * t551 - t548 * t523;
t429 = -qJD(1) * t554 + qJD(3) * t639;
t450 = pkin(3) * t494 - pkin(8) * t493;
t444 = t550 * t450;
t553 = -qJD(4) * t408 - t429 * t547 + t444;
t348 = pkin(4) * t494 - qJ(5) * t433 - qJD(5) * t479 + t553;
t434 = t479 * qJD(4) + t493 * t547;
t589 = qJD(4) * t550;
t555 = t550 * t429 + t438 * t589 + t547 * t450 - t463 * t590;
t351 = -qJ(5) * t434 - qJD(5) * t477 + t555;
t333 = t544 * t348 - t351 * t542;
t394 = t433 * t544 - t434 * t542;
t329 = pkin(5) * t494 - pkin(9) * t394 + t333;
t334 = t542 * t348 + t544 * t351;
t393 = -t433 * t542 - t434 * t544;
t330 = pkin(9) * t393 + t334;
t407 = t550 * t438 - t463 * t547;
t387 = -qJ(5) * t479 + t407;
t374 = pkin(4) * t498 + t387;
t388 = -qJ(5) * t477 + t408;
t612 = t544 * t388;
t356 = t542 * t374 + t612;
t642 = pkin(9) * t568;
t343 = t356 + t642;
t588 = qJD(6) * t546;
t342 = t343 * t588;
t462 = -qJD(3) * pkin(3) - t639;
t418 = pkin(4) * t477 + qJD(5) + t462;
t375 = -pkin(5) * t568 + t418;
t656 = -t546 * t329 - t549 * t330 - t375 * t379 + t342;
t643 = -t549 * t424 + t546 * t568;
t655 = t494 * MDP(28) + (-t379 ^ 2 + t643 ^ 2) * MDP(25) - t379 * MDP(24) * t643;
t517 = t542 * t550 + t544 * t547;
t647 = t498 * t517;
t565 = t542 * t547 - t544 * t550;
t646 = t498 * t565;
t632 = t643 * t492;
t464 = pkin(3) * t509 - pkin(8) * t507;
t456 = t550 * t464;
t634 = -qJ(5) - pkin(8);
t577 = qJD(4) * t634;
t653 = -pkin(4) * t509 - t456 + (qJ(5) * t507 + t577) * t550 + (-qJD(5) + t639) * t547;
t600 = t547 * t464 + t550 * t639;
t620 = t507 * t547;
t652 = -qJ(5) * t620 - qJD(5) * t550 - t547 * t577 + t600;
t651 = t590 - t620;
t576 = t549 * t329 - t546 * t330;
t650 = -t375 * t643 + t576;
t649 = pkin(9) * t424;
t648 = (t543 ^ 2 + t545 ^ 2) * (MDP(7) * qJ(2) + MDP(6));
t468 = t517 * t549 - t546 * t565;
t605 = qJD(6) * t468 - t646 * t546 + t549 * t647;
t645 = t518 * qJD(2);
t510 = t516 * qJD(3);
t578 = t518 * t589;
t644 = -t510 * t547 + t578;
t602 = t542 * t652 + t544 * t653;
t601 = t542 * t653 - t544 * t652;
t473 = t525 * t551 + t548 * t526;
t638 = pkin(4) * t651 - t470;
t567 = -t517 * t546 - t549 * t565;
t606 = qJD(6) * t567 - t546 * t647 - t549 * t646;
t637 = -t468 * t494 - t492 * t606;
t575 = -t549 * t393 + t394 * t546;
t340 = qJD(6) * t643 + t575;
t636 = pkin(4) * t542;
t631 = t379 * t509;
t629 = t643 * t509;
t627 = t433 * t547;
t625 = t477 * t498;
t624 = t477 * t509;
t623 = t479 * t498;
t622 = t479 * t509;
t618 = t518 * t547;
t617 = t518 * t550;
t381 = t542 * t388;
t609 = t547 * t494;
t355 = t544 * t374 - t381;
t341 = pkin(5) * t498 + t355 + t649;
t608 = t549 * t341;
t474 = -t525 * t548 + t526 * t551;
t472 = t550 * t474;
t483 = t550 * t494;
t439 = -qJD(3) * t473 - t554;
t465 = pkin(3) * t511 + pkin(8) * t510;
t457 = t550 * t465;
t466 = pkin(3) * t516 - pkin(8) * t518 + t537;
t564 = qJ(5) * t510 - qJD(5) * t518;
t363 = pkin(4) * t511 - t439 * t547 + t457 + t564 * t550 + (-t472 + (qJ(5) * t518 - t466) * t547) * qJD(4);
t582 = t550 * t439 + t547 * t465 + t466 * t589;
t367 = -qJ(5) * t578 + (-qJD(4) * t474 + t564) * t547 + t582;
t338 = t542 * t363 + t544 * t367;
t360 = t544 * t387 - t381;
t459 = t550 * t466;
t398 = pkin(4) * t516 - qJ(5) * t617 - t474 * t547 + t459;
t599 = t547 * t466 + t472;
t409 = -qJ(5) * t618 + t599;
t370 = t542 * t398 + t544 * t409;
t598 = pkin(5) * t647 + t638;
t527 = t634 * t547;
t528 = t634 * t550;
t476 = t542 * t527 - t544 * t528;
t593 = qJD(3) * t548;
t592 = qJD(3) * t551;
t591 = qJD(4) * t518;
t585 = qJD(1) * qJD(2);
t583 = qJD(6) * t607 + t546 * t393 + t549 * t394;
t580 = -pkin(4) * t550 - pkin(3);
t337 = t544 * t363 - t367 * t542;
t359 = -t387 * t542 - t612;
t369 = t544 * t398 - t409 * t542;
t475 = t544 * t527 + t528 * t542;
t574 = t498 * t550;
t430 = qJD(1) * t645 - t522 * t593 + t523 * t592;
t440 = -t525 * t593 + t526 * t592 + t645;
t573 = -t492 * t605 + t567 * t494;
t572 = pkin(4) * t618 + t473;
t449 = -pkin(9) * t565 + t476;
t571 = pkin(5) * t509 - pkin(9) * t646 + qJD(6) * t449 - t602;
t448 = -pkin(9) * t517 + t475;
t570 = pkin(9) * t647 - qJD(6) * t448 - t601;
t332 = t546 * t341 + t549 * t343;
t452 = t517 * t518;
t453 = t565 * t518;
t569 = -t549 * t452 + t453 * t546;
t413 = -t452 * t546 - t453 * t549;
t562 = -t498 * t651 + t483;
t536 = pkin(4) * t544 + pkin(5);
t561 = t536 * t546 + t549 * t636;
t560 = t536 * t549 - t546 * t636;
t559 = pkin(4) * t644 + t440;
t557 = -t510 * t550 - t518 * t590;
t397 = pkin(4) * t434 + t430;
t556 = -pkin(8) * t494 + t462 * t498;
t339 = t424 * t588 + t583;
t485 = pkin(5) * t565 + t580;
t471 = t494 * t516;
t415 = pkin(4) * t479 - pkin(5) * t424;
t414 = pkin(5) * t452 + t572;
t411 = -t510 * t565 + t517 * t591;
t410 = t510 * t517 + t565 * t591;
t371 = -pkin(5) * t410 + t559;
t368 = -pkin(5) * t393 + t397;
t358 = -pkin(9) * t452 + t370;
t357 = pkin(5) * t516 + pkin(9) * t453 + t369;
t353 = qJD(6) * t413 - t549 * t410 - t411 * t546;
t352 = qJD(6) * t569 + t410 * t546 - t411 * t549;
t345 = t360 + t649;
t344 = t359 - t642;
t336 = pkin(9) * t410 + t338;
t335 = pkin(5) * t511 + pkin(9) * t411 + t337;
t331 = -t343 * t546 + t608;
t1 = [(t493 * t518 - t509 * t510) * MDP(8) + (-t493 * t516 - t494 * t518 - t507 * t510 - t509 * t511) * MDP(9) + (t494 * t537 + t511 * t524) * MDP(13) + (t493 * t537 - t510 * t524) * MDP(14) + (t433 * t617 + t479 * t557) * MDP(15) + (-(-t477 * t550 - t479 * t547) * t510 + (-t627 - t434 * t550 + (t477 * t547 - t479 * t550) * qJD(4)) * t518) * MDP(16) + (t433 * t516 + t479 * t511 + t483 * t518 + t498 * t557) * MDP(17) + (-t434 * t516 - t477 * t511 - t498 * t644 - t518 * t609) * MDP(18) + (t498 * t511 + t471) * MDP(19) + ((-t474 * t589 + t457) * t498 + t459 * t494 + (-t463 * t589 + t444) * t516 + t407 * t511 + t440 * t477 + t473 * t434 + t462 * t578 + ((-qJD(4) * t466 - t439) * t498 - t474 * t494 + (-qJD(4) * t438 - t429) * t516 + t430 * t518 - t462 * t510) * t547) * MDP(20) + (-(-t474 * t590 + t582) * t498 - t599 * t494 - t555 * t516 - t408 * t511 + t440 * t479 + t473 * t433 + t430 * t617 + t557 * t462) * MDP(21) + (t333 * t453 - t334 * t452 + t337 * t424 + t338 * t568 + t355 * t411 + t356 * t410 - t369 * t394 + t370 * t393) * MDP(22) + (t333 * t369 + t334 * t370 + t355 * t337 + t356 * t338 + t397 * t572 + t418 * t559) * MDP(23) + (t339 * t413 + t352 * t643) * MDP(24) + (t339 * t569 - t340 * t413 + t352 * t379 - t353 * t643) * MDP(25) + (t339 * t516 + t352 * t492 + t413 * t494 + t511 * t643) * MDP(26) + (-t340 * t516 - t353 * t492 + t379 * t511 + t494 * t569) * MDP(27) + (t492 * t511 + t471) * MDP(28) + ((t335 * t549 - t336 * t546) * t492 + (t357 * t549 - t358 * t546) * t494 + t576 * t516 + t331 * t511 - t371 * t379 + t414 * t340 - t368 * t569 + t375 * t353 + ((-t357 * t546 - t358 * t549) * t492 - t332 * t516) * qJD(6)) * MDP(29) + (-t332 * t511 + t414 * t339 + t342 * t516 + t375 * t352 + t368 * t413 + t371 * t643 + (-(-qJD(6) * t358 + t335) * t492 - t357 * t494 - t329 * t516) * t546 + (-(qJD(6) * t357 + t336) * t492 - t358 * t494 - (qJD(6) * t341 + t330) * t516) * t549) * MDP(30) + 0.2e1 * t585 * t648 + (-MDP(10) * t510 - MDP(11) * t511 - MDP(13) * t440 - MDP(14) * t439) * qJD(3); t531 * MDP(14) + (t562 - t624) * MDP(20) + (-t498 ^ 2 * t550 - t609 - t622) * MDP(21) + (t393 * t517 + t394 * t565 - t424 * t647 - t568 * t646) * MDP(22) + (-t333 * t565 + t334 * t517 - t355 * t647 - t356 * t646 - t418 * t509) * MDP(23) + (t573 + t631) * MDP(29) + (-t629 + t637) * MDP(30) + ((t543 * t594 + t545 * t595 + t509) * MDP(13) + (t507 - t579) * MDP(14)) * qJD(3) - qJD(1) ^ 2 * t648; -t507 ^ 2 * MDP(9) + (t531 + (-t507 - t579) * qJD(3)) * MDP(10) + (qJD(3) * t470 - t430) * MDP(13) + (-t507 * t524 + t516 * t585) * MDP(14) + (t479 * t574 + t627) * MDP(15) + ((t433 - t625) * t550 + (-t434 - t623) * t547) * MDP(16) + (t498 * t574 + t609 - t622) * MDP(17) + (t562 + t624) * MDP(18) + (-pkin(3) * t434 - t430 * t550 - t470 * t477 + (-pkin(8) * t589 - t456) * t498 + (t498 * t639 + t556) * t547) * MDP(20) + (-pkin(3) * t433 + t430 * t547 - t470 * t479 + (pkin(8) * t590 + t600) * t498 + t556 * t550) * MDP(21) + (-t333 * t517 - t334 * t565 + t355 * t646 - t356 * t647 + t393 * t476 - t394 * t475 + t424 * t602 + t568 * t601) * MDP(22) + (t333 * t475 + t334 * t476 + t602 * t355 + t601 * t356 + t397 * t580 + t418 * t638) * MDP(23) + (t339 * t468 + t606 * t643) * MDP(24) + (t339 * t567 - t340 * t468 + t379 * t606 - t605 * t643) * MDP(25) + (-t629 - t637) * MDP(26) + (t573 - t631) * MDP(27) + ((t448 * t549 - t449 * t546) * t494 + t485 * t340 - t368 * t567 + (t546 * t570 - t549 * t571) * t492 - t598 * t379 + t605 * t375) * MDP(29) + (-(t448 * t546 + t449 * t549) * t494 + t485 * t339 + t368 * t468 + (t546 * t571 + t549 * t570) * t492 + t598 * t643 + t606 * t375) * MDP(30) + (-t524 * MDP(13) - t498 * MDP(19) - t407 * MDP(20) + t408 * MDP(21) - t492 * MDP(28) - t331 * MDP(29) + t332 * MDP(30) - MDP(8) * t507 + t509 * MDP(9)) * t509; t479 * t477 * MDP(15) + (-t477 ^ 2 + t479 ^ 2) * MDP(16) + (t433 + t625) * MDP(17) + (-t434 + t623) * MDP(18) + t494 * MDP(19) + (t408 * t498 - t462 * t479 + t553) * MDP(20) + (t407 * t498 + t462 * t477 - t555) * MDP(21) + ((t393 * t542 - t394 * t544) * pkin(4) + (t355 - t360) * t568 + (-t356 - t359) * t424) * MDP(22) + (-t355 * t359 - t356 * t360 + (t333 * t544 + t334 * t542 - t418 * t479) * pkin(4)) * MDP(23) + (t339 - t630) * MDP(26) + (-t340 + t632) * MDP(27) + (t560 * t494 - (t344 * t549 - t345 * t546) * t492 + t415 * t379 + (-t492 * t561 - t332) * qJD(6) + t650) * MDP(29) + (-t561 * t494 + (t344 * t546 + t345 * t549) * t492 - t415 * t643 + (-t492 * t560 - t608) * qJD(6) + t656) * MDP(30) + t655; (-t424 ^ 2 - t568 ^ 2) * MDP(22) + (-t355 * t424 - t356 * t568 + t397) * MDP(23) + (t340 + t632) * MDP(29) + (t339 + t630) * MDP(30); (t583 - t630) * MDP(26) + (-t575 + t632) * MDP(27) + (t332 * t492 + t650) * MDP(29) + (t331 * t492 + t656) * MDP(30) + (MDP(26) * t628 - MDP(27) * t643 - MDP(29) * t332 - MDP(30) * t608) * qJD(6) + t655;];
tauc  = t1;
