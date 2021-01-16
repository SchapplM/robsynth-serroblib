% Calculate vector of inverse dynamics joint torques for
% S5RRRPP6
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
%   pkin=[a2,a3,a4,a5,d1,d2,d3,theta4]';
% MDP [25x1]
%   Minimal dynamic parameter vector (fixed base model)
%   see S5RRRPP6_convert_par2_MPV_fixb.m
% 
% Output:
% tau [5x1]
%   joint torques of inverse dynamics (contains inertial, gravitational Coriolis and centrifugal forces)

% Quelle: HybrDyn-Toolbox
% Datum: 2021-01-15 22:38
% Revision: 24b2e7d74a0c1a3b64fa2f8f5ad758691ad61af3 (2021-01-15)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tau = S5RRRPP6_invdynJ_fixb_mdp_slag_vp(qJ, qJD, qJDD, g, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(5,1),zeros(3,1),zeros(8,1),zeros(25,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRRPP6_invdynJ_fixb_mdp_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RRRPP6_invdynJ_fixb_mdp_slag_vp: qJD has to be [5x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [5 1]), ...
  'S5RRRPP6_invdynJ_fixb_mdp_slag_vp: qJDD has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RRRPP6_invdynJ_fixb_mdp_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RRRPP6_invdynJ_fixb_mdp_slag_vp: pkin has to be [8x1] (double)');
assert(isreal(MDP) && all(size(MDP) == [25 1]), ...
  'S5RRRPP6_invdynJ_fixb_mdp_slag_vp: MDP has to be [25x1] (double)'); 

%% Symbolic Calculation
% From invdyn_joint_fixb_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2021-01-15 22:37:25
% EndTime: 2021-01-15 22:37:43
% DurationCPUTime: 6.70s
% Computational Cost: add. (4745->531), mult. (10574->659), div. (0->0), fcn. (6963->10), ass. (0->217)
t513 = sin(qJ(2));
t566 = qJDD(1) * t513;
t493 = pkin(6) * t566;
t516 = cos(qJ(2));
t567 = qJD(1) * qJD(2);
t552 = t516 * t567;
t421 = -qJDD(2) * pkin(2) + pkin(6) * t552 + t493;
t579 = qJD(1) * t516;
t480 = -qJD(3) + t579;
t505 = g(3) * t516;
t514 = sin(qJ(1));
t517 = cos(qJ(1));
t541 = g(1) * t517 + g(2) * t514;
t625 = t541 * t513;
t524 = -t505 + t625;
t633 = qJD(3) * pkin(7) * t480 - t421 + t524;
t511 = qJ(4) + pkin(7);
t515 = cos(qJ(3));
t463 = t511 * t515;
t509 = sin(pkin(8));
t510 = cos(pkin(8));
t512 = sin(qJ(3));
t554 = t511 * t512;
t403 = t510 * t463 - t509 * t554;
t501 = t516 * qJDD(1);
t623 = -t513 * t567 + t501;
t441 = qJDD(3) - t623;
t506 = qJ(3) + pkin(8);
t497 = sin(t506);
t632 = t403 * t441 + t497 * t524;
t580 = qJD(1) * t513;
t563 = t512 * t580;
t569 = t515 * qJD(2);
t447 = t563 - t569;
t576 = qJD(2) * t512;
t449 = t515 * t580 + t576;
t386 = t510 * t447 + t449 * t509;
t631 = t386 * t480;
t443 = t509 * t515 + t510 * t512;
t424 = t443 * qJD(3);
t587 = t443 * t579 - t424;
t573 = qJD(3) * t512;
t555 = t509 * t573;
t562 = t512 * t579;
t571 = qJD(3) * t515;
t602 = t510 * t515;
t586 = -t509 * t562 - t510 * t571 + t579 * t602 + t555;
t574 = qJD(2) * t516;
t560 = t512 * t574;
t630 = t513 * t571 + t560;
t572 = qJD(3) * t513;
t629 = -qJD(1) * t572 + qJDD(2);
t381 = (qJD(2) * (qJD(3) + t579) + t566) * t512 - t629 * t515;
t537 = -t447 * t509 + t510 * t449;
t628 = t537 ^ 2;
t455 = pkin(4) * t510 + qJ(5) * t509 + pkin(3);
t461 = -t509 * pkin(4) + qJ(5) * t510;
t400 = -t455 * t512 + t461 * t515;
t627 = -pkin(6) + t400;
t464 = -qJD(2) * pkin(2) + pkin(6) * t580;
t404 = pkin(3) * t447 + qJD(4) + t464;
t344 = pkin(4) * t386 - qJ(5) * t537 + t404;
t626 = t344 * t537;
t595 = t515 * t516;
t536 = pkin(3) * t513 - qJ(4) * t595;
t544 = pkin(2) * t513 - pkin(7) * t516;
t450 = t544 * qJD(1);
t583 = pkin(6) * t563 + t515 * t450;
t379 = qJD(1) * t536 + t583;
t431 = t512 * t450;
t599 = t513 * t515;
t600 = t512 * t516;
t390 = t431 + (-t599 * pkin(6) - qJ(4) * t600) * qJD(1);
t550 = qJD(3) * t511;
t570 = qJD(4) * t515;
t418 = -t512 * t550 + t570;
t527 = -qJD(4) * t512 - t515 * t550;
t591 = (-t379 + t527) * t510 + (t390 - t418) * t509;
t495 = pkin(6) * t579;
t545 = -t495 + (-t562 + t573) * pkin(3);
t622 = -t455 * t515 - t461 * t512;
t596 = t514 * t516;
t427 = t512 * t596 + t515 * t517;
t594 = t516 * t517;
t429 = -t512 * t594 + t514 * t515;
t621 = -g(1) * t429 + g(2) * t427;
t619 = -2 * pkin(1);
t618 = pkin(6) * t512;
t616 = g(1) * t514;
t614 = g(2) * t517;
t613 = g(3) * t513;
t489 = t511 * t513;
t612 = pkin(1) + t489;
t462 = -pkin(2) * t516 - pkin(7) * t513 - pkin(1);
t436 = t462 * qJD(1);
t465 = qJD(2) * pkin(7) + t495;
t397 = t436 * t512 + t465 * t515;
t370 = -qJ(4) * t447 + t397;
t367 = t510 * t370;
t396 = t515 * t436 - t465 * t512;
t369 = -qJ(4) * t449 + t396;
t342 = t369 * t509 + t367;
t611 = t342 * t537;
t610 = t370 * t509;
t380 = qJD(3) * t569 + (t552 + t566) * t515 + t629 * t512;
t609 = t380 * t512;
t607 = t447 * t480;
t606 = t449 * t480;
t605 = t449 * t515;
t601 = t512 * t513;
t598 = t513 * t517;
t498 = cos(t506);
t597 = t514 * t498;
t593 = t517 * t497;
t451 = t544 * qJD(2);
t398 = qJD(1) * t451 + qJDD(1) * t462;
t392 = t515 * t398;
t420 = t623 * pkin(6) + qJDD(2) * pkin(7);
t332 = pkin(3) * t441 - qJ(4) * t380 - qJD(3) * t397 - qJD(4) * t449 - t420 * t512 + t392;
t530 = t512 * t398 + t515 * t420 + t436 * t571 - t465 * t573;
t335 = -qJ(4) * t381 - qJD(4) * t447 + t530;
t323 = t510 * t332 - t509 * t335;
t324 = t509 * t332 + t510 * t335;
t483 = pkin(6) * t595;
t575 = qJD(2) * t513;
t584 = t515 * t451 + t575 * t618;
t353 = -t513 * t570 + t536 * qJD(2) + (-t483 + (qJ(4) * t513 - t462) * t512) * qJD(3) + t584;
t585 = t512 * t451 + t462 * t571;
t358 = (-pkin(6) * qJD(2) - qJ(4) * qJD(3)) * t599 + (-qJD(4) * t513 + (-pkin(6) * qJD(3) - qJ(4) * qJD(2)) * t516) * t512 + t585;
t331 = t509 * t353 + t510 * t358;
t592 = -t587 * pkin(4) + t586 * qJ(5) - qJD(5) * t443 + t545;
t365 = -pkin(3) * t480 + t369;
t341 = t509 * t365 + t367;
t590 = pkin(4) * t580 - t591;
t352 = t509 * t379 + t510 * t390;
t346 = qJ(5) * t580 + t352;
t374 = t510 * t418 + t509 * t527;
t589 = t374 - t346;
t588 = t374 - t352;
t445 = t515 * t462;
t394 = -qJ(4) * t599 + t445 + (-pkin(3) - t618) * t516;
t582 = t512 * t462 + t483;
t401 = -qJ(4) * t601 + t582;
t360 = t509 * t394 + t510 * t401;
t452 = pkin(3) * t601 + t513 * pkin(6);
t507 = t513 ^ 2;
t581 = -t516 ^ 2 + t507;
t578 = qJD(2) * t447;
t577 = qJD(2) * t449;
t343 = t369 * t510 - t610;
t568 = qJD(5) - t343;
t405 = t630 * pkin(3) + pkin(6) * t574;
t492 = pkin(3) * t515 + pkin(2);
t561 = t480 * t569;
t559 = t516 * t569;
t558 = t480 * t573;
t557 = t480 * t571;
t348 = t380 * t509 + t510 * t381;
t548 = -qJD(3) * t436 - t420;
t484 = t513 * t616;
t546 = -g(2) * t598 + t484;
t410 = t497 * t596 + t498 * t517;
t412 = t516 * t593 - t597;
t543 = g(1) * t410 - g(2) * t412;
t411 = t498 * t596 - t593;
t413 = t497 * t514 + t498 * t594;
t542 = g(1) * t411 - g(2) * t413;
t540 = t465 * t571 - t392;
t539 = -pkin(7) * t441 + qJD(3) * t464;
t330 = t353 * t510 - t358 * t509;
t340 = t365 * t510 - t610;
t349 = t380 * t510 - t381 * t509;
t359 = t394 * t510 - t401 * t509;
t322 = -t441 * pkin(4) + qJDD(5) - t323;
t534 = -pkin(6) * qJDD(2) + t567 * t619;
t532 = t441 * t512 - t557;
t531 = t441 * t515 + t558;
t519 = qJD(1) ^ 2;
t529 = pkin(1) * t519 + t541;
t518 = qJD(2) ^ 2;
t528 = pkin(6) * t518 + qJDD(1) * t619 + t614;
t402 = t463 * t509 + t510 * t554;
t525 = g(2) * t513 * t597 - t402 * t441 + (g(1) * t598 - t505) * t498;
t364 = pkin(3) * t381 + qJDD(4) + t421;
t522 = g(1) * t413 + g(2) * t411 + t498 * t613 - t324;
t521 = g(1) * t412 + g(2) * t410 - t342 * t480 + t497 * t613 + t323;
t520 = -t403 * t348 + t349 * t402 - t374 * t386 - t516 * t541 - t613;
t325 = pkin(4) * t348 - qJ(5) * t349 - qJD(5) * t537 + t364;
t491 = pkin(3) * t512 + pkin(6);
t490 = t511 * t516;
t487 = -pkin(3) * t510 - pkin(4);
t485 = pkin(3) * t509 + qJ(5);
t467 = t492 * t516;
t442 = t509 * t512 - t602;
t435 = t441 * qJ(5);
t430 = t512 * t514 + t515 * t594;
t428 = t512 * t517 - t514 * t595;
t419 = t467 + t612;
t416 = -t509 * t601 + t510 * t599;
t415 = t443 * t513;
t399 = pkin(2) - t622;
t393 = t399 * t516;
t383 = pkin(4) * t442 - qJ(5) * t443 - t492;
t376 = t393 + t612;
t372 = t424 * t513 + t509 * t560 - t510 * t559;
t371 = -t509 * t559 - t630 * t510 + t513 * t555;
t366 = pkin(4) * t415 - qJ(5) * t416 + t452;
t357 = pkin(4) * t516 - t359;
t356 = -qJ(5) * t516 + t360;
t345 = pkin(3) * t449 + pkin(4) * t537 + qJ(5) * t386;
t338 = -qJ(5) * t480 + t341;
t337 = pkin(4) * t480 + qJD(5) - t340;
t336 = -pkin(4) * t371 + qJ(5) * t372 - qJD(5) * t416 + t405;
t327 = -pkin(4) * t575 - t330;
t326 = qJ(5) * t575 - qJD(5) * t516 + t331;
t321 = -qJD(5) * t480 + t324 + t435;
t1 = [qJDD(1) * MDP(1) + (qJDD(1) * t507 + 0.2e1 * t513 * t552) * MDP(4) + 0.2e1 * (t501 * t513 - t567 * t581) * MDP(5) + (qJDD(2) * t513 + t516 * t518) * MDP(6) + (qJDD(2) * t516 - t513 * t518) * MDP(7) + (t534 * t513 + (-t528 + t616) * t516) * MDP(9) + (t513 * t528 + t516 * t534 - t484) * MDP(10) + (t380 * t599 + (-t512 * t572 + t559) * t449) * MDP(11) + ((-t447 * t515 - t449 * t512) * t574 + (-t609 - t381 * t515 + (t447 * t512 - t605) * qJD(3)) * t513) * MDP(12) + ((-t380 - t561) * t516 + (t531 + t577) * t513) * MDP(13) + ((t480 * t576 + t381) * t516 + (-t532 - t578) * t513) * MDP(14) + (-t441 * t516 - t480 * t575) * MDP(15) + (-(-t462 * t573 + t584) * t480 + t445 * t441 - g(1) * t428 - g(2) * t430 + ((t557 + t578) * pkin(6) + (-pkin(6) * t441 + qJD(2) * t464 - t548) * t512 + t540) * t516 + (pkin(6) * t381 + qJD(2) * t396 + t421 * t512 + t464 * t571) * t513) * MDP(16) + (t585 * t480 - t582 * t441 - g(1) * t427 - g(2) * t429 + (t464 * t569 + (-t558 + t577) * pkin(6) + t530) * t516 + (-t464 * t573 - t397 * qJD(2) + t421 * t515 + (t380 - t561) * pkin(6)) * t513) * MDP(17) + (-t323 * t516 - t330 * t480 + t340 * t575 + t348 * t452 + t359 * t441 + t364 * t415 - t371 * t404 + t386 * t405 + t542) * MDP(18) + (t324 * t516 + t331 * t480 - t341 * t575 + t349 * t452 - t360 * t441 + t364 * t416 - t372 * t404 + t405 * t537 - t543) * MDP(19) + (-t323 * t416 - t324 * t415 - t330 * t537 - t331 * t386 + t340 * t372 + t341 * t371 - t348 * t360 - t349 * t359 + t546) * MDP(20) + (t324 * t360 + t341 * t331 + t323 * t359 + t340 * t330 + t364 * t452 + t404 * t405 - g(1) * (-t419 * t514 + t491 * t517) - g(2) * (t419 * t517 + t491 * t514)) * MDP(21) + (t322 * t516 + t325 * t415 + t327 * t480 + t336 * t386 - t337 * t575 - t344 * t371 + t348 * t366 - t357 * t441 + t542) * MDP(22) + (-t321 * t415 + t322 * t416 - t326 * t386 + t327 * t537 - t337 * t372 + t338 * t371 - t348 * t356 + t349 * t357 + t546) * MDP(23) + (-t321 * t516 - t325 * t416 - t326 * t480 - t336 * t537 + t338 * t575 + t344 * t372 - t349 * t366 + t356 * t441 + t543) * MDP(24) + (t321 * t356 + t338 * t326 + t325 * t366 + t344 * t336 + t322 * t357 + t337 * t327 - g(1) * (-t376 * t514 - t627 * t517) - g(2) * (t376 * t517 - t627 * t514)) * MDP(25) + (-t614 + t616) * MDP(2) + t541 * MDP(3); MDP(6) * t566 + MDP(7) * t501 + qJDD(2) * MDP(8) + (t513 * t529 - t493 - t505) * MDP(9) + (t613 + (-pkin(6) * qJDD(1) + t529) * t516) * MDP(10) + (-t480 * t605 + t609) * MDP(11) + ((t380 + t607) * t515 + (-t381 + t606) * t512) * MDP(12) + ((-t449 * t513 + t480 * t595) * qJD(1) + t532) * MDP(13) + ((t447 * t513 - t480 * t600) * qJD(1) + t531) * MDP(14) + t480 * MDP(15) * t580 + (-pkin(2) * t381 + t583 * t480 + t539 * t512 + (-t396 * t513 + (-pkin(6) * t447 - t464 * t512) * t516) * qJD(1) + t633 * t515) * MDP(16) + (-pkin(2) * t380 - t431 * t480 + t539 * t515 + (-t464 * t595 + t397 * t513 + (-t449 * t516 + t480 * t599) * pkin(6)) * qJD(1) - t633 * t512) * MDP(17) + (-t340 * t580 - t348 * t492 + t364 * t442 + t386 * t545 - t404 * t587 - t480 * t591 + t525) * MDP(18) + (t341 * t580 - t349 * t492 + t364 * t443 - t404 * t586 + t480 * t588 + t537 * t545 - t632) * MDP(19) + (-t323 * t443 - t324 * t442 + t340 * t586 + t341 * t587 + t352 * t386 - t537 * t591 + t520) * MDP(20) + (t324 * t403 - t323 * t402 - t364 * t492 - g(3) * (t467 + t489) - t541 * (-t492 * t513 + t490) + t545 * t404 + t588 * t341 + t591 * t340) * MDP(21) + (t325 * t442 + t337 * t580 - t344 * t587 + t348 * t383 + t386 * t592 + t480 * t590 + t525) * MDP(22) + (-t321 * t442 + t322 * t443 - t337 * t586 + t338 * t587 + t346 * t386 + t537 * t590 + t520) * MDP(23) + (-t325 * t443 - t338 * t580 + t344 * t586 - t349 * t383 - t480 * t589 - t537 * t592 + t632) * MDP(24) + (t321 * t403 + t325 * t383 + t322 * t402 - g(3) * (t393 + t489) - t541 * (-t399 * t513 + t490) + t592 * t344 + t589 * t338 + t590 * t337) * MDP(25) + (-MDP(4) * t513 * t516 + MDP(5) * t581) * t519; t449 * t447 * MDP(11) + (-t447 ^ 2 + t449 ^ 2) * MDP(12) + (t380 - t607) * MDP(13) + (-t381 - t606) * MDP(14) + t441 * MDP(15) + (-t397 * t480 - t449 * t464 + (t548 + t613) * t512 - t540 + t621) * MDP(16) + (g(1) * t430 - g(2) * t428 + g(3) * t599 - t396 * t480 + t447 * t464 - t530) * MDP(17) + (-t537 * t404 + (-t386 * t449 + t441 * t510) * pkin(3) + t521) * MDP(18) + (-t343 * t480 + t386 * t404 + (-t441 * t509 - t449 * t537) * pkin(3) + t522) * MDP(19) + (t341 * t537 - t611 + (-t348 * t509 - t349 * t510) * pkin(3) + (-t340 + t343) * t386) * MDP(20) + (t340 * t342 - t341 * t343 + (g(3) * t601 + t323 * t510 + t324 * t509 - t404 * t449 + t621) * pkin(3)) * MDP(21) + (-t626 - t345 * t386 - qJDD(5) + (pkin(4) - t487) * t441 + t521) * MDP(22) + (t338 * t537 - t348 * t485 + t349 * t487 - t611 + (t337 - t568) * t386) * MDP(23) + (-t344 * t386 + t345 * t537 + t441 * t485 + t435 + (-0.2e1 * qJD(5) + t343) * t480 - t522) * MDP(24) + (t321 * t485 + t322 * t487 - t344 * t345 - t337 * t342 - g(1) * (t400 * t594 - t514 * t622) - g(2) * (t400 * t596 + t622 * t517) - t400 * t613 + t568 * t338) * MDP(25); (t340 * t537 + t341 * t386 + t364 + t505) * MDP(21) + (-t337 * t537 + t338 * t386 + t325 + t505) * MDP(25) + (-MDP(21) - MDP(25)) * t625 + (MDP(20) + MDP(23)) * (-t386 ^ 2 - t628) + (MDP(18) + MDP(22)) * (-t480 * t537 + t348) + (MDP(19) - MDP(24)) * (t349 + t631); (t386 * t537 - t441) * MDP(22) + (t349 - t631) * MDP(23) + (-t480 ^ 2 - t628) * MDP(24) + (t338 * t480 + t626 - g(1) * (t442 * t514 + t443 * t594) - g(2) * (-t442 * t517 + t443 * t596) - t443 * t613 + t322) * MDP(25);];
tau = t1;
