% Calculate vector of inverse dynamics joint torques for
% S5RRRPR3
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
%   pkin=[a2,a3,a4,a5,d1,d2,d3,d5,theta4]';
% MDP [24x1]
%   Minimal dynamic parameter vector (fixed base model)
%   see S5RRRPR3_convert_par2_MPV_fixb.m
% 
% Output:
% tau [5x1]
%   joint torques of inverse dynamics (contains inertial, gravitational Coriolis and centrifugal forces)

% Quelle: HybrDyn-Toolbox
% Datum: 2021-01-15 23:00
% Revision: 24b2e7d74a0c1a3b64fa2f8f5ad758691ad61af3 (2021-01-15)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tau = S5RRRPR3_invdynJ_fixb_mdp_slag_vp(qJ, qJD, qJDD, g, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(5,1),zeros(3,1),zeros(9,1),zeros(24,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRRPR3_invdynJ_fixb_mdp_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RRRPR3_invdynJ_fixb_mdp_slag_vp: qJD has to be [5x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [5 1]), ...
  'S5RRRPR3_invdynJ_fixb_mdp_slag_vp: qJDD has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RRRPR3_invdynJ_fixb_mdp_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5RRRPR3_invdynJ_fixb_mdp_slag_vp: pkin has to be [9x1] (double)');
assert(isreal(MDP) && all(size(MDP) == [24 1]), ...
  'S5RRRPR3_invdynJ_fixb_mdp_slag_vp: MDP has to be [24x1] (double)'); 

%% Symbolic Calculation
% From invdyn_joint_fixb_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2021-01-15 22:59:46
% EndTime: 2021-01-15 22:59:57
% DurationCPUTime: 4.09s
% Computational Cost: add. (3165->367), mult. (4849->455), div. (0->0), fcn. (3381->16), ass. (0->189)
t471 = sin(qJ(2));
t525 = qJD(2) * t471;
t449 = pkin(1) * t525;
t475 = cos(qJ(2));
t560 = pkin(1) * t475;
t529 = -qJD(1) * t449 + qJDD(1) * t560;
t459 = qJDD(1) + qJDD(2);
t559 = pkin(2) * t459;
t390 = -t529 - t559;
t465 = qJ(1) + qJ(2);
t454 = sin(t465);
t455 = cos(t465);
t568 = g(2) * t455 + g(3) * t454;
t570 = t390 + t568;
t569 = qJ(4) + pkin(7);
t461 = qJD(1) + qJD(2);
t467 = cos(pkin(9));
t474 = cos(qJ(3));
t537 = t467 * t474;
t516 = t461 * t537;
t466 = sin(pkin(9));
t470 = sin(qJ(3));
t541 = t466 * t470;
t385 = t461 * t541 - t516;
t473 = cos(qJ(5));
t372 = t473 * t385;
t407 = t466 * t474 + t467 * t470;
t387 = t407 * t461;
t469 = sin(qJ(5));
t545 = t387 * t469;
t333 = -t372 - t545;
t460 = qJD(3) + qJD(5);
t546 = t333 * t460;
t494 = t385 * t469 - t473 * t387;
t547 = t494 * t460;
t564 = g(2) * t454 - g(3) * t455;
t520 = qJDD(1) * t471;
t524 = qJD(2) * t475;
t391 = pkin(7) * t459 + (qJD(1) * t524 + t520) * pkin(1);
t489 = qJ(4) * t459 + qJD(4) * t461 + t391;
t561 = pkin(1) * t471;
t519 = qJD(1) * t561;
t510 = t569 * t461 + t519;
t493 = qJD(3) * t510;
t316 = qJDD(3) * pkin(3) - t470 * t489 - t474 * t493;
t319 = -t470 * t493 + t474 * t489;
t293 = t467 * t316 - t319 * t466;
t523 = qJD(3) * t470;
t513 = t461 * t523;
t420 = t466 * t513;
t522 = qJD(3) * t474;
t512 = t461 * t522;
t339 = t407 * t459 + t467 * t512 - t420;
t287 = qJDD(3) * pkin(4) - pkin(8) * t339 + t293;
t294 = t466 * t316 + t467 * t319;
t399 = t407 * qJD(3);
t425 = t459 * t537;
t338 = t399 * t461 + t459 * t541 - t425;
t288 = -pkin(8) * t338 + t294;
t374 = t510 * t470;
t368 = qJD(3) * pkin(3) - t374;
t375 = t510 * t474;
t539 = t467 * t375;
t323 = t466 * t368 + t539;
t554 = pkin(8) * t385;
t306 = t323 - t554;
t444 = pkin(3) * t474 + pkin(2);
t526 = qJD(1) * t475;
t518 = pkin(1) * t526;
t384 = -t444 * t461 + qJD(4) - t518;
t340 = pkin(4) * t385 + t384;
t462 = qJ(3) + pkin(9);
t452 = qJ(5) + t462;
t437 = sin(t452);
t438 = cos(t452);
t521 = qJD(5) * t469;
t567 = g(1) * t437 - t469 * t287 - t473 * t288 + t306 * t521 - t340 * t333 + t564 * t438;
t543 = t437 * t455;
t544 = t437 * t454;
t566 = -g(1) * t438 + g(2) * t544 - g(3) * t543 + t473 * t287 - t469 * t288 + t340 * t494;
t458 = qJDD(3) + qJDD(5);
t565 = t458 * MDP(22) + t333 * MDP(18) * t494 + (-t333 ^ 2 + t494 ^ 2) * MDP(19);
t453 = t474 * qJD(4);
t511 = qJD(3) * t569;
t396 = -t470 * t511 + t453;
t397 = -qJD(4) * t470 - t474 * t511;
t533 = -t396 * t466 + t467 * t397 + t407 * t518;
t406 = -t537 + t541;
t532 = t467 * t396 + t466 * t397 + t406 * t518;
t477 = qJD(3) ^ 2;
t563 = pkin(7) * t477 - t559;
t562 = pkin(3) * t513 + qJDD(4);
t509 = t473 * t338 + t339 * t469;
t292 = -qJD(5) * t494 + t509;
t558 = pkin(2) * t461;
t557 = pkin(3) * t466;
t556 = pkin(3) * t470;
t553 = pkin(8) * t387;
t400 = t406 * qJD(3);
t552 = pkin(8) * t400;
t551 = pkin(8) * t407;
t550 = g(1) * t474;
t542 = t461 * t470;
t362 = t466 * t375;
t536 = t470 * t474;
t322 = t467 * t368 - t362;
t304 = qJD(3) * pkin(4) + t322 - t553;
t535 = t473 * t304;
t443 = pkin(7) + t561;
t534 = -qJ(4) - t443;
t507 = qJD(3) * t534;
t517 = pkin(1) * t524;
t357 = t470 * t507 + t474 * t517 + t453;
t358 = (-qJD(4) - t517) * t470 + t474 * t507;
t318 = t467 * t357 + t466 * t358;
t325 = -t467 * t374 - t362;
t404 = t534 * t470;
t456 = t474 * qJ(4);
t405 = t443 * t474 + t456;
t349 = t466 * t404 + t467 * t405;
t428 = t569 * t470;
t429 = pkin(7) * t474 + t456;
t367 = -t466 * t428 + t467 * t429;
t531 = t454 * t444 - t455 * t569;
t463 = t470 ^ 2;
t528 = -t474 ^ 2 + t463;
t448 = pkin(3) * t523;
t515 = -qJD(5) * t372 - t469 * t338 + t473 * t339;
t514 = t461 * t525;
t371 = pkin(4) * t399 + t448;
t317 = -t357 * t466 + t467 * t358;
t324 = t374 * t466 - t539;
t348 = t467 * t404 - t405 * t466;
t366 = -t467 * t428 - t429 * t466;
t508 = t455 * t444 + t454 * t569;
t506 = t461 * t519;
t345 = -t444 * t459 - t529 + t562;
t305 = pkin(4) * t338 + t345;
t350 = t473 * t406 + t407 * t469;
t309 = -qJD(5) * t350 - t399 * t469 - t400 * t473;
t351 = -t406 * t469 + t407 * t473;
t505 = g(2) * t543 + g(3) * t544 + t305 * t351 + t340 * t309;
t450 = sin(t462);
t504 = t345 * t407 - t384 * t400 + t568 * t450;
t419 = -t518 - t558;
t503 = t419 * t522 + t570 * t470;
t502 = t371 - t519;
t403 = t406 * pkin(8);
t344 = -t403 + t367;
t499 = qJD(5) * t344 - t533 - t552;
t343 = t366 - t551;
t395 = t399 * pkin(8);
t498 = -qJD(5) * t343 + t395 - t532;
t497 = -t469 * t304 - t473 * t306;
t326 = t348 - t551;
t327 = -t403 + t349;
t496 = t326 * t473 - t327 * t469;
t495 = t326 * t469 + t327 * t473;
t379 = pkin(4) * t406 - t444;
t439 = pkin(3) * t467 + pkin(4);
t492 = t439 * t469 + t473 * t557;
t491 = t439 * t473 - t469 * t557;
t490 = -t293 * t407 - t294 * t406 + t322 * t400 - t323 * t399 - t564;
t488 = -t568 + t529;
t291 = -t387 * t521 + t515;
t487 = -t419 * t461 - t391 + t564;
t310 = qJD(5) * t351 + t473 * t399 - t400 * t469;
t486 = t305 * t350 + t340 * t310 - t438 * t568;
t451 = cos(t462);
t485 = t345 * t406 + t384 * t399 - t451 * t568;
t484 = (-t291 * t350 - t292 * t351 + t309 * t333 + t310 * t494) * MDP(19) + (t291 * t351 - t309 * t494) * MDP(18) + (t309 * t460 + t351 * t458) * MDP(20) + (-t310 * t460 - t350 * t458) * MDP(21) + 0.2e1 * (-qJD(3) * t461 * t528 + t459 * t536) * MDP(8) + (t459 * t463 + 0.2e1 * t470 * t512) * MDP(7) + (qJDD(3) * t474 - t470 * t477) * MDP(10) + (qJDD(3) * t470 + t474 * t477) * MDP(9) + t459 * MDP(4);
t483 = -t568 + t506;
t445 = -pkin(2) - t560;
t482 = pkin(1) * t514 + t443 * t477 + t445 * t459;
t480 = -pkin(7) * qJDD(3) + (t518 - t558) * qJD(3);
t479 = -qJDD(3) * t443 + (t445 * t461 - t517) * qJD(3);
t476 = cos(qJ(1));
t472 = sin(qJ(1));
t426 = -t444 - t560;
t417 = t449 + t448;
t401 = t419 * t523;
t370 = t379 - t560;
t361 = t371 + t449;
t356 = pkin(3) * t542 + pkin(4) * t387;
t308 = t325 - t553;
t307 = t324 + t554;
t302 = -t395 + t318;
t301 = t317 + t552;
t1 = [(-g(2) * t476 - g(3) * t472) * MDP(2) + (g(2) * t472 - g(3) * t476) * MDP(3) + (t294 * t349 + t323 * t318 + t293 * t348 + t322 * t317 + t345 * t426 + t384 * t417 - g(2) * (pkin(1) * t476 + t508) - g(3) * (pkin(1) * t472 + t531)) * MDP(17) + (qJD(3) * t317 + qJDD(3) * t348 + t338 * t426 + t385 * t417 + t485) * MDP(14) + (-t361 * t494 + t370 * t291 - (qJD(5) * t496 + t301 * t469 + t302 * t473) * t460 - t495 * t458 + t505) * MDP(24) + (((-qJDD(1) - t459) * t471 + (-qJD(1) - t461) * t524) * pkin(1) + t564) * MDP(6) + ((t459 * t475 - t514) * pkin(1) + t488) * MDP(5) + (t401 + t479 * t470 + (-t482 - t570) * t474) * MDP(12) + (t470 * t482 + t474 * t479 + t503) * MDP(13) + (-t361 * t333 + t370 * t292 + (-qJD(5) * t495 + t301 * t473 - t302 * t469) * t460 + t496 * t458 + t486) * MDP(23) + (-qJD(3) * t318 - qJDD(3) * t349 + t339 * t426 + t387 * t417 + t504) * MDP(15) + (-t317 * t387 - t318 * t385 - t338 * t349 - t339 * t348 + t490) * MDP(16) + qJDD(1) * MDP(1) + t484; (t379 * t291 - (t343 * t469 + t344 * t473) * t458 + (t469 * t499 + t473 * t498) * t460 - t502 * t494 + t505) * MDP(24) + (t294 * t367 + t293 * t366 - t345 * t444 - g(2) * t508 - g(3) * t531 + (t448 - t519) * t384 + t532 * t323 + t533 * t322) * MDP(17) + ((-t520 + (-qJD(2) + t461) * t526) * pkin(1) + t564) * MDP(6) + (t379 * t292 + (t343 * t473 - t344 * t469) * t458 + (t469 * t498 - t473 * t499) * t460 - t502 * t333 + t486) * MDP(23) + (t480 * t474 + (-t506 + t563) * t470 + t503) * MDP(13) + (t401 + t480 * t470 + (-t390 + t483 - t563) * t474) * MDP(12) + (t483 + t529) * MDP(5) + (-t387 * t519 - qJDD(3) * t367 - t339 * t444 + (t387 * t556 - t532) * qJD(3) + t504) * MDP(15) + (-t338 * t367 - t339 * t366 - t385 * t532 - t387 * t533 + t490) * MDP(16) + (-t385 * t519 + qJDD(3) * t366 - t338 * t444 + (t385 * t556 + t533) * qJD(3) + t485) * MDP(14) + t484; qJDD(3) * MDP(11) + (t470 * t487 - t550) * MDP(12) + (g(1) * t470 + t474 * t487) * MDP(13) + (-g(1) * t451 - qJD(3) * t324 - t384 * t387 + t564 * t450 + (qJDD(3) * t467 - t385 * t542) * pkin(3) + t293) * MDP(14) + (g(1) * t450 + qJD(3) * t325 + t384 * t385 + t564 * t451 + (-qJDD(3) * t466 - t387 * t542) * pkin(3) - t294) * MDP(15) + ((t323 + t324) * t387 + (-t322 + t325) * t385 + (-t338 * t466 - t339 * t467) * pkin(3)) * MDP(16) + (-t322 * t324 - t323 * t325 + (-t550 + t293 * t467 + t294 * t466 + (-t384 * t461 + t564) * t470) * pkin(3)) * MDP(17) + (t291 - t546) * MDP(20) + (-t292 - t547) * MDP(21) + (t491 * t458 + t356 * t333 - (t307 * t473 - t308 * t469) * t460 + (-t460 * t492 + t497) * qJD(5) + t566) * MDP(23) + (-t492 * t458 + t356 * t494 + (t307 * t469 + t308 * t473) * t460 + (-t460 * t491 - t535) * qJD(5) + t567) * MDP(24) + (MDP(10) * t474 + MDP(9) * t470) * t459 + (-MDP(7) * t536 + MDP(8) * t528) * t461 ^ 2 + t565; -t425 * MDP(14) - t420 * MDP(15) + (-t385 ^ 2 - t387 ^ 2) * MDP(16) + (t322 * t387 + t323 * t385 - t488 + t562) * MDP(17) + (t292 - t547) * MDP(23) + (t291 + t546) * MDP(24) + (MDP(14) * t541 + MDP(15) * t407 - MDP(17) * t444) * t459 + (0.2e1 * t387 * MDP(14) + (-t385 + t516) * MDP(15)) * qJD(3); (t515 - t546) * MDP(20) + (-t509 - t547) * MDP(21) + (-t460 * t497 + t566) * MDP(23) + ((-t306 * t469 + t535) * t460 + t567) * MDP(24) + (-MDP(20) * t545 + MDP(21) * t494 + MDP(23) * t497 - MDP(24) * t535) * qJD(5) + t565;];
tau = t1;
