% Calculate vector of inverse dynamics joint torques for
% S5RRRRP8
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
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d2,d3,d4]';
% MDP [26x1]
%   Minimal dynamic parameter vector (fixed base model)
%   see S5RRRRP8_convert_par2_MPV_fixb.m
% 
% Output:
% tau [5x1]
%   joint torques of inverse dynamics (contains inertial, gravitational Coriolis and centrifugal forces)

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 22:02
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tau = S5RRRRP8_invdynJ_fixb_mdp_slag_vp(qJ, qJD, qJDD, g, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(5,1),zeros(3,1),zeros(8,1),zeros(26,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRRRP8_invdynJ_fixb_mdp_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RRRRP8_invdynJ_fixb_mdp_slag_vp: qJD has to be [5x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [5 1]), ...
  'S5RRRRP8_invdynJ_fixb_mdp_slag_vp: qJDD has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RRRRP8_invdynJ_fixb_mdp_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RRRRP8_invdynJ_fixb_mdp_slag_vp: pkin has to be [8x1] (double)');
assert(isreal(MDP) && all(size(MDP) == [26 1]), ...
  'S5RRRRP8_invdynJ_fixb_mdp_slag_vp: MDP has to be [26x1] (double)'); 

%% Symbolic Calculation
% From invdyn_joint_fixb_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 22:02:21
% EndTime: 2019-12-31 22:02:30
% DurationCPUTime: 5.92s
% Computational Cost: add. (3893->461), mult. (8731->606), div. (0->0), fcn. (5964->10), ass. (0->203)
t492 = cos(qJ(3));
t490 = sin(qJ(2));
t493 = cos(qJ(2));
t579 = t492 * t493;
t516 = pkin(3) * t490 - pkin(8) * t579;
t603 = pkin(7) + pkin(8);
t545 = qJD(3) * t603;
t523 = pkin(2) * t490 - pkin(7) * t493;
t433 = t523 * qJD(1);
t489 = sin(qJ(3));
t562 = qJD(1) * t490;
t538 = t489 * t562;
t566 = pkin(6) * t538 + t492 * t433;
t624 = qJD(1) * t516 + t492 * t545 + t566;
t411 = t489 * t433;
t581 = t490 * t492;
t582 = t489 * t493;
t623 = -t411 - (-pkin(6) * t581 - pkin(8) * t582) * qJD(1) - t489 * t545;
t553 = qJD(3) * t490;
t533 = qJD(1) * t553;
t550 = qJD(1) * qJD(2);
t534 = t493 * t550;
t548 = qJDD(1) * t490;
t619 = qJD(2) * qJD(3) + t534 + t548;
t364 = (qJDD(2) - t533) * t489 + t619 * t492;
t556 = qJD(2) * t492;
t425 = -t538 + t556;
t558 = qJD(2) * t489;
t426 = t492 * t562 + t558;
t488 = sin(qJ(4));
t525 = t619 * t489 + t492 * t533;
t505 = t492 * qJDD(2) - t525;
t602 = cos(qJ(4));
t536 = t602 * qJD(4);
t551 = qJD(4) * t488;
t328 = -t602 * t364 - t425 * t536 + t426 * t551 - t488 * t505;
t513 = t488 * t425 + t426 * t602;
t329 = qJD(4) * t513 + t488 * t364 - t602 * t505;
t371 = -t602 * t425 + t426 * t488;
t369 = t371 ^ 2;
t479 = t493 * qJDD(1);
t611 = -t490 * t550 + t479;
t422 = qJDD(3) - t611;
t417 = qJDD(4) + t422;
t561 = qJD(1) * t493;
t465 = -qJD(3) + t561;
t453 = -qJD(4) + t465;
t604 = t513 ^ 2;
t622 = t417 * MDP(22) + (-t453 * t513 - t329) * MDP(21) + t371 * MDP(18) * t513 + (-t371 * t453 - t328) * MDP(20) + (-t369 + t604) * MDP(19);
t621 = t371 * qJ(5);
t584 = t488 * t489;
t512 = t492 * t602 - t584;
t609 = qJD(3) + qJD(4);
t610 = t602 * qJD(3) + t536;
t571 = -t492 * t610 + t512 * t561 + t584 * t609;
t428 = t488 * t492 + t489 * t602;
t384 = t609 * t428;
t570 = -t428 * t561 + t384;
t555 = qJD(2) * t493;
t543 = t489 * t555;
t552 = qJD(3) * t492;
t620 = t490 * t552 + t543;
t444 = -qJD(2) * pkin(2) + pkin(6) * t562;
t387 = -pkin(3) * t425 + t444;
t487 = qJ(3) + qJ(4);
t480 = sin(t487);
t481 = cos(t487);
t494 = cos(qJ(1));
t491 = sin(qJ(1));
t580 = t491 * t493;
t394 = t480 * t494 - t481 * t580;
t578 = t493 * t494;
t396 = t480 * t491 + t481 * t578;
t436 = t523 * qJD(2);
t440 = -pkin(2) * t493 - pkin(7) * t490 - pkin(1);
t385 = qJD(1) * t436 + qJDD(1) * t440;
t377 = t492 * t385;
t418 = t440 * qJD(1);
t476 = pkin(6) * t561;
t445 = qJD(2) * pkin(7) + t476;
t381 = t418 * t489 + t445 * t492;
t402 = t611 * pkin(6) + qJDD(2) * pkin(7);
t323 = pkin(3) * t422 - pkin(8) * t364 - qJD(3) * t381 - t402 * t489 + t377;
t554 = qJD(3) * t489;
t509 = t489 * t385 + t492 * t402 + t418 * t552 - t445 * t554;
t327 = pkin(8) * t505 + t509;
t380 = t492 * t418 - t445 * t489;
t354 = -pkin(8) * t426 + t380;
t349 = -pkin(3) * t465 + t354;
t355 = pkin(8) * t425 + t381;
t526 = -t488 * t323 - t602 * t327 - t349 * t536 + t355 * t551;
t595 = g(3) * t490;
t618 = g(1) * t396 - g(2) * t394 + t371 * t387 + t481 * t595 + t526;
t616 = qJ(5) * t513;
t345 = pkin(4) * t371 + qJD(5) + t387;
t615 = t345 * t513;
t424 = t492 * t440;
t600 = pkin(6) * t489;
t379 = -pkin(8) * t581 + t424 + (-pkin(3) - t600) * t493;
t467 = pkin(6) * t579;
t565 = t489 * t440 + t467;
t583 = t489 * t490;
t386 = -pkin(8) * t583 + t565;
t572 = t488 * t379 + t602 * t386;
t601 = pkin(3) * t489;
t614 = pkin(3) * t554 - t561 * t601 - t476;
t446 = t603 * t489;
t447 = t603 * t492;
t567 = -t488 * t446 + t602 * t447;
t613 = qJD(4) * t567 + t623 * t488 + t624 * t602;
t612 = -t446 * t536 - t447 * t551 - t624 * t488 + t623 * t602;
t522 = g(1) * t494 + g(2) * t491;
t474 = pkin(6) * t548;
t591 = qJDD(2) * pkin(2);
t403 = pkin(6) * t534 + t474 - t591;
t594 = g(3) * t493;
t504 = -t490 * t522 + t594;
t608 = -qJD(3) * pkin(7) * t465 + t403 + t504;
t585 = t481 * t494;
t393 = t480 * t580 + t585;
t586 = t481 * t491;
t395 = -t480 * t578 + t586;
t607 = -g(1) * t395 + g(2) * t393 + t480 * t595;
t353 = t602 * t355;
t331 = t488 * t349 + t353;
t500 = -qJD(4) * t331 + t602 * t323 - t488 * t327;
t606 = -t387 * t513 + t500 + t607;
t483 = t492 * pkin(3);
t593 = pkin(2) + t483;
t438 = pkin(4) * t480 + t601;
t592 = pkin(6) + t438;
t590 = t364 * t489;
t589 = t425 * t465;
t588 = t426 * t465;
t587 = t426 * t492;
t351 = t488 * t355;
t330 = t602 * t349 - t351;
t321 = t330 - t616;
t318 = -pkin(4) * t453 + t321;
t577 = -t321 + t318;
t576 = -t570 * qJ(5) + qJD(5) * t512 + t612;
t575 = -pkin(4) * t562 + t571 * qJ(5) - t428 * qJD(5) - t613;
t574 = t602 * t354 - t351;
t569 = t489 * t436 + t440 * t552;
t557 = qJD(2) * t490;
t568 = t492 * t436 + t557 * t600;
t439 = pkin(4) * t481 + t483;
t437 = pkin(3) * t583 + t490 * pkin(6);
t485 = t490 ^ 2;
t564 = -t493 ^ 2 + t485;
t560 = qJD(2) * t425;
t559 = qJD(2) * t426;
t388 = t620 * pkin(3) + pkin(6) * t555;
t544 = t465 * t556;
t542 = t465 * t554;
t541 = t465 * t552;
t540 = t489 * t553;
t532 = -t354 * t488 - t353;
t529 = t602 * t379 - t386 * t488;
t528 = -t602 * t446 - t447 * t488;
t527 = -qJD(3) * t418 - t402;
t524 = t602 * t555;
t521 = g(1) * t491 - g(2) * t494;
t520 = t445 * t552 - t377;
t519 = -pkin(7) * t422 + qJD(3) * t444;
t432 = pkin(2) + t439;
t484 = -qJ(5) - t603;
t518 = t432 * t493 - t484 * t490;
t515 = pkin(1) + t518;
t514 = -0.2e1 * pkin(1) * t550 - pkin(6) * qJDD(2);
t511 = t422 * t489 - t541;
t510 = t422 * t492 + t542;
t340 = t516 * qJD(2) + (-t467 + (pkin(8) * t490 - t440) * t489) * qJD(3) + t568;
t342 = -t620 * pkin(8) + (-t490 * t556 - t493 * t554) * pkin(6) + t569;
t508 = t488 * t340 + t602 * t342 + t379 * t536 - t386 * t551;
t496 = qJD(1) ^ 2;
t506 = pkin(1) * t496 + t522;
t495 = qJD(2) ^ 2;
t503 = 0.2e1 * qJDD(1) * pkin(1) - pkin(6) * t495 + t521;
t499 = -qJD(4) * t572 + t602 * t340 - t488 * t342;
t348 = -pkin(3) * t505 + t403;
t317 = t329 * pkin(4) + qJDD(5) + t348;
t472 = pkin(3) * t602 + pkin(4);
t409 = t489 * t491 + t492 * t578;
t408 = -t489 * t578 + t491 * t492;
t407 = t489 * t494 - t491 * t579;
t406 = t489 * t580 + t492 * t494;
t399 = t512 * t490;
t398 = t428 * t490;
t359 = qJ(5) * t512 + t567;
t358 = -qJ(5) * t428 + t528;
t344 = t489 * t524 - t488 * t540 - t551 * t583 + (t488 * t555 + t490 * t610) * t492;
t343 = t384 * t490 + t488 * t543 - t492 * t524;
t337 = -qJ(5) * t398 + t572;
t336 = -pkin(4) * t493 - qJ(5) * t399 + t529;
t325 = t574 - t616;
t324 = t532 + t621;
t322 = t331 - t621;
t316 = -qJ(5) * t344 - qJD(5) * t398 + t508;
t315 = pkin(4) * t557 + t343 * qJ(5) - t399 * qJD(5) + t499;
t314 = -qJ(5) * t329 - qJD(5) * t371 - t526;
t313 = t417 * pkin(4) + t328 * qJ(5) - qJD(5) * t513 + t500;
t1 = [t521 * MDP(2) + t522 * MDP(3) + (t569 * t465 - t565 * t422 - g(1) * t406 - g(2) * t408 + (t444 * t556 + (-t542 + t559) * pkin(6) + t509) * t493 + (-t444 * t554 - t381 * qJD(2) + t403 * t492 + (t364 - t544) * pkin(6)) * t490) * MDP(17) + (t364 * t581 + (t492 * t555 - t540) * t426) * MDP(11) + (t490 * t514 + t493 * t503) * MDP(9) + (-t490 * t503 + t493 * t514) * MDP(10) + (qJDD(2) * t490 + t493 * t495) * MDP(6) + (qJDD(2) * t493 - t490 * t495) * MDP(7) + (qJDD(1) * t485 + 0.2e1 * t490 * t534) * MDP(4) + ((t425 * t492 - t426 * t489) * t555 + (t492 * t505 - t590 + (-t425 * t489 - t587) * qJD(3)) * t490) * MDP(12) + (t314 * t337 + t322 * t316 + t313 * t336 + t318 * t315 + t317 * (pkin(4) * t398 + t437) + t345 * (pkin(4) * t344 + t388) + (-g(1) * t592 - g(2) * t515) * t494 + (g(1) * t515 - g(2) * t592) * t491) * MDP(26) + (-t422 * t493 - t465 * t557) * MDP(15) + (-t417 * t493 - t453 * t557) * MDP(22) + (t329 * t493 + t344 * t453 - t371 * t557 - t398 * t417) * MDP(21) + (-g(1) * t394 - g(2) * t396 + t437 * t329 + t330 * t557 + t387 * t344 + t348 * t398 + t388 * t371 + t417 * t529 - t453 * t499 - t493 * t500) * MDP(23) + ((-t364 - t544) * t493 + (t510 + t559) * t490) * MDP(13) + ((t465 * t558 - t505) * t493 + (-t511 + t560) * t490) * MDP(14) + qJDD(1) * MDP(1) + 0.2e1 * (t479 * t490 - t550 * t564) * MDP(5) + (-(-t440 * t554 + t568) * t465 + t424 * t422 - g(1) * t407 - g(2) * t409 + ((t541 - t560) * pkin(6) + (-pkin(6) * t422 + qJD(2) * t444 - t527) * t489 + t520) * t493 + (-pkin(6) * t505 + t380 * qJD(2) + t403 * t489 + t444 * t552) * t490) * MDP(16) + (-t313 * t399 - t314 * t398 - t315 * t513 - t316 * t371 + t318 * t343 - t322 * t344 + t328 * t336 - t329 * t337 + t490 * t521) * MDP(25) + (-g(1) * t393 - g(2) * t395 - t437 * t328 - t331 * t557 - t387 * t343 + t348 * t399 + t388 * t513 - t417 * t572 + t453 * t508 - t493 * t526) * MDP(24) + (t328 * t398 - t329 * t399 + t343 * t371 - t344 * t513) * MDP(19) + (-t328 * t399 - t343 * t513) * MDP(18) + (t328 * t493 + t343 * t453 + t399 * t417 + t513 * t557) * MDP(20); MDP(6) * t548 + MDP(7) * t479 + qJDD(2) * MDP(8) + (t490 * t506 - t474 - t594) * MDP(9) + (t595 + (-pkin(6) * qJDD(1) + t506) * t493) * MDP(10) + (-t465 * t587 + t590) * MDP(11) + ((t364 - t589) * t492 + (t505 + t588) * t489) * MDP(12) + ((-t426 * t490 + t465 * t579) * qJD(1) + t511) * MDP(13) + ((-t425 * t490 - t465 * t582) * qJD(1) + t510) * MDP(14) + (-pkin(2) * t525 + t566 * t465 + t519 * t489 + (-t380 * t490 + (pkin(6) * t425 - t444 * t489) * t493) * qJD(1) + (t591 - t608) * t492) * MDP(16) + (-pkin(2) * t364 - t411 * t465 + t519 * t492 + (-t444 * t579 + t381 * t490 + (-t426 * t493 + t465 * t581) * pkin(6)) * qJD(1) + t608 * t489) * MDP(17) + (-t328 * t428 - t513 * t571) * MDP(18) + (-t328 * t512 - t329 * t428 + t371 * t571 - t513 * t570) * MDP(19) + (t417 * t428 + t453 * t571) * MDP(20) + (t417 * t512 + t453 * t570) * MDP(21) + (-t593 * t329 - t348 * t512 + t417 * t528 - t481 * t594 + (g(1) * t585 + g(2) * t586) * t490 + t613 * t453 + t570 * t387 + t614 * t371) * MDP(23) + (t328 * t593 + t348 * t428 - t571 * t387 - t567 * t417 + t453 * t612 + t504 * t480 + t513 * t614) * MDP(24) + (-t313 * t428 + t314 * t512 + t318 * t571 - t322 * t570 + t328 * t358 - t329 * t359 - t371 * t576 - t493 * t522 - t513 * t575 - t595) * MDP(25) + (t314 * t359 + t313 * t358 + t317 * (-pkin(4) * t512 - t593) - g(3) * t518 + (pkin(4) * t570 + t614) * t345 + t576 * t322 + t575 * t318 + t522 * (t432 * t490 + t484 * t493)) * MDP(26) + (t465 * MDP(15) - t513 * MDP(20) + t371 * MDP(21) + t453 * MDP(22) - t330 * MDP(23) + t331 * MDP(24)) * t562 + (-t490 * t493 * MDP(4) + t564 * MDP(5)) * t496; -t426 * t425 * MDP(11) + (-t425 ^ 2 + t426 ^ 2) * MDP(12) + (t364 + t589) * MDP(13) + (t505 - t588) * MDP(14) + t422 * MDP(15) + (-g(1) * t408 + g(2) * t406 - t381 * t465 - t426 * t444 + (t527 + t595) * t489 - t520) * MDP(16) + (g(1) * t409 - g(2) * t407 + g(3) * t581 - t380 * t465 - t425 * t444 - t509) * MDP(17) + (t532 * t453 + (-t426 * t371 + t417 * t602 + t453 * t551) * pkin(3) + t606) * MDP(23) + (-t574 * t453 + (-t488 * t417 - t426 * t513 + t453 * t536) * pkin(3) + t618) * MDP(24) + (-t318 * t371 + t322 * t513 + t324 * t513 + t325 * t371 + t472 * t328 + (-t329 * t488 + (-t371 * t602 + t488 * t513) * qJD(4)) * pkin(3)) * MDP(25) + (t313 * t472 - t322 * t325 - t318 * t324 - pkin(4) * t615 - g(1) * (-t438 * t578 + t439 * t491) - g(2) * (-t438 * t580 - t439 * t494) + t438 * t595 + (t314 * t488 - t345 * t426 + (-t318 * t488 + t322 * t602) * qJD(4)) * pkin(3)) * MDP(26) + t622; (-t331 * t453 + t606) * MDP(23) + (-t330 * t453 + t618) * MDP(24) + (pkin(4) * t328 - t371 * t577) * MDP(25) + (t577 * t322 + (t313 + t607 - t615) * pkin(4)) * MDP(26) + t622; (-t369 - t604) * MDP(25) + (t318 * t513 + t322 * t371 + t317 + t504) * MDP(26);];
tau = t1;
