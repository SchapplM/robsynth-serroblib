% Calculate vector of inverse dynamics joint torques for
% S5RPRPR5
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
%   pkin=[a2,a3,a4,a5,d1,d3,d5,theta2,theta4]';
% MDP [24x1]
%   Minimal dynamic parameter vector (fixed base model)
%   see S5RPRPR5_convert_par2_MPV_fixb.m
% 
% Output:
% tau [5x1]
%   joint torques of inverse dynamics (contains inertial, gravitational Coriolis and centrifugal forces)

% Quelle: HybrDyn-Toolbox
% Datum: 2022-01-23 09:26
% Revision: fd3771346c4aea32fdeb66112c511235427c26a7 (2022-01-20)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tau = S5RPRPR5_invdynJ_fixb_mdp_slag_vp(qJ, qJD, qJDD, g, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(5,1),zeros(3,1),zeros(9,1),zeros(24,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPRPR5_invdynJ_fixb_mdp_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RPRPR5_invdynJ_fixb_mdp_slag_vp: qJD has to be [5x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [5 1]), ...
  'S5RPRPR5_invdynJ_fixb_mdp_slag_vp: qJDD has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RPRPR5_invdynJ_fixb_mdp_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5RPRPR5_invdynJ_fixb_mdp_slag_vp: pkin has to be [9x1] (double)');
assert(isreal(MDP) && all(size(MDP) == [24 1]), ...
  'S5RPRPR5_invdynJ_fixb_mdp_slag_vp: MDP has to be [24x1] (double)'); 

%% Symbolic Calculation
% From invdyn_joint_fixb_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2022-01-23 09:26:06
% EndTime: 2022-01-23 09:26:12
% DurationCPUTime: 4.24s
% Computational Cost: add. (2705->395), mult. (6665->533), div. (0->0), fcn. (4900->14), ass. (0->201)
t424 = sin(pkin(8));
t418 = t424 ^ 2;
t530 = 0.2e1 * t418;
t431 = cos(qJ(3));
t490 = qJD(3) * t431;
t541 = qJ(2) * t490;
t423 = sin(pkin(9));
t425 = cos(pkin(9));
t428 = sin(qJ(3));
t379 = t423 * t431 + t425 * t428;
t443 = qJD(1) * t379;
t346 = t424 * t443;
t430 = cos(qJ(5));
t504 = t430 * t346;
t495 = qJD(1) * t424;
t473 = t428 * t495;
t457 = t423 * t473;
t493 = qJD(1) * t431;
t472 = t424 * t493;
t350 = t425 * t472 - t457;
t427 = sin(qJ(5));
t516 = t350 * t427;
t305 = t504 + t516;
t426 = cos(pkin(8));
t494 = qJD(1) * t426;
t402 = -qJD(3) + t494;
t396 = -qJD(5) + t402;
t517 = t305 * t396;
t519 = qJDD(1) * pkin(1);
t429 = sin(qJ(1));
t432 = cos(qJ(1));
t534 = g(1) * t429 - g(2) * t432;
t447 = -qJDD(2) + t519 + t534;
t480 = qJDD(1) * t428;
t468 = t424 * t480;
t482 = qJD(1) * qJD(3);
t470 = t431 * t482;
t540 = t424 * t470 + t468;
t483 = qJD(1) * qJD(2);
t484 = qJ(2) * qJDD(1);
t446 = t483 + t484;
t386 = -pkin(2) * t426 - pkin(6) * t424 - pkin(1);
t365 = t386 * qJDD(1) + qJDD(2);
t356 = t431 * t365;
t481 = qJDD(1) * t426;
t400 = -qJDD(3) + t481;
t520 = qJ(2) * t431;
t401 = t426 * t520;
t489 = qJD(4) * t424;
t492 = qJD(2) * t426;
t440 = -t428 * t492 - t431 * t489;
t512 = t424 * t431;
t476 = qJ(4) * t512;
t521 = qJ(2) * t428;
t478 = t426 * t521;
t441 = -t476 - t478;
t513 = t424 * t428;
t477 = qJ(4) * t513;
t366 = t386 * qJD(1) + qJD(2);
t491 = qJD(3) * t366;
t283 = -t428 * t491 - pkin(3) * t400 + t356 + t441 * qJDD(1) + ((-t401 + t477) * qJD(3) + t440) * qJD(1);
t435 = t441 * qJD(3) - t428 * t489;
t471 = t431 * t483;
t460 = qJDD(1) * t401 + t428 * t365 + t366 * t490 + t426 * t471;
t289 = -qJ(4) * t468 + t435 * qJD(1) + t460;
t272 = t425 * t283 - t289 * t423;
t479 = qJDD(1) * t431;
t467 = t424 * t479;
t385 = t425 * t467;
t442 = t379 * qJD(3);
t469 = t423 * t480;
t318 = -t385 + (qJD(1) * t442 + t469) * t424;
t270 = -pkin(4) * t400 + pkin(7) * t318 + t272;
t357 = t431 * t366;
t326 = t441 * qJD(1) + t357;
t316 = -pkin(3) * t402 + t326;
t327 = -qJ(4) * t473 + qJD(1) * t401 + t428 * t366;
t511 = t425 * t327;
t294 = t423 * t316 + t511;
t528 = pkin(7) * t346;
t280 = t294 - t528;
t488 = qJD(5) * t427;
t539 = t427 * t270 - t280 * t488;
t273 = t423 * t283 + t425 * t289;
t474 = t423 * t467 + t540 * t425;
t317 = qJD(3) * t457 - t474;
t271 = pkin(7) * t317 + t273;
t367 = pkin(3) * t473 + qJ(2) * t495 + qJD(4);
t325 = pkin(4) * t346 + t367;
t420 = qJ(3) + pkin(9);
t416 = qJ(5) + t420;
t407 = sin(t416);
t408 = cos(t416);
t510 = t426 * t429;
t342 = t407 * t432 - t408 * t510;
t509 = t426 * t432;
t344 = t407 * t429 + t408 * t509;
t523 = g(3) * t424;
t538 = g(1) * t344 - g(2) * t342 - t430 * t271 + t325 * t305 + t408 * t523 - t539;
t391 = -qJDD(5) + t400;
t382 = t391 * MDP(22);
t450 = -t346 * t427 + t430 * t350;
t537 = t305 * t450 * MDP(18) + (-t305 ^ 2 + t450 ^ 2) * MDP(19) - t382;
t535 = t446 * t426;
t518 = t450 * t396;
t503 = t431 * t432;
t508 = t428 * t429;
t370 = t426 * t508 + t503;
t506 = t429 * t431;
t507 = t428 * t432;
t372 = -t426 * t507 + t506;
t533 = -g(1) * t372 + g(2) * t370;
t514 = t423 * t428;
t449 = -t425 * t431 + t514;
t532 = t449 * qJD(3);
t341 = t407 * t510 + t408 * t432;
t343 = -t407 * t509 + t408 * t429;
t464 = t430 * t270 - t427 * t271;
t531 = -g(1) * t343 + g(2) * t341 - t325 * t450 + t407 * t523 + t464;
t463 = t430 * t317 + t318 * t427;
t275 = t450 * qJD(5) - t463;
t419 = t426 ^ 2;
t529 = pkin(3) * t423;
t527 = pkin(7) * t350;
t522 = MDP(9) * t424;
t433 = qJD(1) ^ 2;
t515 = t418 * t433;
t321 = t423 * t327;
t293 = t425 * t316 - t321;
t278 = -pkin(4) * t402 + t293 - t527;
t505 = t430 * t278;
t500 = t386 * t490 + t431 * t492;
t314 = t435 + t500;
t315 = (-t401 + (qJ(4) * t424 - t386) * t428) * qJD(3) + t440;
t287 = t425 * t314 + t423 * t315;
t297 = t425 * t326 - t321;
t377 = t431 * t386;
t331 = -t476 + t377 + (-pkin(3) - t521) * t426;
t499 = t428 * t386 + t401;
t336 = -t477 + t499;
t300 = t423 * t331 + t425 * t336;
t502 = -t426 * t443 + t442;
t501 = -t449 * t494 + t532;
t375 = (pkin(3) * t490 + qJD(2)) * t424;
t380 = pkin(3) * t513 + t424 * qJ(2);
t498 = t418 + t419;
t422 = t431 ^ 2;
t497 = t428 ^ 2 - t422;
t496 = MDP(10) * t424;
t486 = t400 * MDP(11);
t485 = qJD(3) + t402;
t475 = -qJD(5) * t504 + t427 * t317 - t430 * t318;
t465 = t498 * t433;
t286 = -t314 * t423 + t425 * t315;
t296 = -t326 * t423 - t511;
t299 = t425 * t331 - t336 * t423;
t462 = qJD(1) * t485;
t461 = t400 + t481;
t459 = 0.2e1 * t498;
t458 = qJD(3) * t478;
t455 = g(1) * t432 + g(2) * t429;
t453 = qJD(5) * t379 + t502;
t452 = -qJD(5) * t449 - t501;
t451 = -t427 * t278 - t430 * t280;
t363 = t379 * t424;
t364 = t449 * t424;
t319 = t430 * t363 - t364 * t427;
t320 = -t363 * t427 - t364 * t430;
t332 = t540 * pkin(3) + t446 * t424 + qJDD(4);
t448 = qJD(3) * (t402 + t494);
t410 = pkin(3) * t425 + pkin(4);
t445 = t410 * t427 + t430 * t529;
t444 = t410 * t430 - t427 * t529;
t274 = -t350 * t488 + t475;
t439 = -t402 ^ 2 - t515;
t437 = t459 * t483 - t455;
t414 = cos(t420);
t413 = sin(t420);
t409 = pkin(3) * t428 + qJ(2);
t373 = t426 * t503 + t508;
t371 = -t426 * t506 + t507;
t362 = (qJ(4) + pkin(6)) * t424 + pkin(1) + (pkin(3) * t431 + pkin(2)) * t426;
t361 = t413 * t429 + t414 * t509;
t360 = -t413 * t509 + t414 * t429;
t359 = t413 * t432 - t414 * t510;
t358 = t413 * t510 + t414 * t432;
t353 = t424 * t442;
t349 = t424 * t532;
t335 = pkin(3) * t472 + pkin(4) * t350;
t333 = pkin(4) * t363 + t380;
t328 = -pkin(4) * t349 + t375;
t298 = -pkin(4) * t317 + t332;
t295 = -pkin(7) * t363 + t300;
t292 = -pkin(4) * t426 + pkin(7) * t364 + t299;
t291 = t320 * qJD(5) - t430 * t349 - t353 * t427;
t290 = -t319 * qJD(5) + t349 * t427 - t353 * t430;
t285 = t297 - t527;
t284 = t296 + t528;
t277 = pkin(7) * t349 + t287;
t276 = pkin(7) * t353 + t286;
t1 = [qJDD(1) * MDP(1) + t534 * MDP(2) + t455 * MDP(3) + (t459 * t484 + t437) * MDP(5) + (t447 * pkin(1) + (t484 * t498 + t437) * qJ(2)) * MDP(6) + (qJDD(1) * t422 - 0.2e1 * t428 * t470) * t418 * MDP(7) + (-t428 * t479 + t497 * t482) * MDP(8) * t530 + (t428 * t448 - t461 * t431) * t522 + (t461 * t428 + t431 * t448) * t496 + (-g(1) * t371 - g(2) * t373 - t377 * t400 + (t530 + t419) * qJD(1) * t541 + (qJD(3) * t386 * t402 + t446 * t530) * t428) * MDP(12) + ((-t458 + t500) * t402 + t499 * t400 - g(1) * t370 - g(2) * t372 + (t471 + (-t428 * t482 + t479) * qJ(2)) * t530) * MDP(13) + (-g(1) * t359 - g(2) * t361 - t286 * t402 - t299 * t400 - t317 * t380 + t332 * t363 + t346 * t375 - t349 * t367) * MDP(14) + (-g(1) * t358 - g(2) * t360 + t287 * t402 + t300 * t400 - t318 * t380 - t332 * t364 + t350 * t375 - t353 * t367) * MDP(15) + (t272 * t364 - t273 * t363 - t286 * t350 - t287 * t346 + t293 * t353 + t294 * t349 + t299 * t318 + t300 * t317 + t424 * t534) * MDP(16) + (t273 * t300 + t294 * t287 + t272 * t299 + t293 * t286 + t332 * t380 + t367 * t375 - g(1) * (-t362 * t429 + t409 * t432) - g(2) * (t362 * t432 + t409 * t429)) * MDP(17) + (t274 * t320 + t290 * t450) * MDP(18) + (-t274 * t319 - t275 * t320 - t290 * t305 - t291 * t450) * MDP(19) + (-t290 * t396 - t320 * t391) * MDP(20) + (t291 * t396 + t319 * t391) * MDP(21) + (-(t292 * t430 - t295 * t427) * t391 + t328 * t305 + t333 * t275 + t298 * t319 + t325 * t291 - g(1) * t342 - g(2) * t344 + (-t276 * t430 + t277 * t427 - (-t292 * t427 - t295 * t430) * qJD(5)) * t396) * MDP(23) + (-g(1) * t341 - g(2) * t343 + t333 * t274 + t325 * t290 + t298 * t320 + t328 * t450 + ((-qJD(5) * t295 + t276) * t396 + t292 * t391) * t427 + ((qJD(5) * t292 + t277) * t396 + t295 * t391) * t430) * MDP(24) + ((t447 + t519) * MDP(4) + t486 + (-t356 + t402 * t541 + (qJ(2) * t400 + qJD(2) * t402 + t491 + t535) * t428) * MDP(12) + (-qJD(1) * t458 + t460) * MDP(13) - t272 * MDP(14) + t273 * MDP(15) - t274 * MDP(20) + t275 * MDP(21) + t382 + (-t451 * qJD(5) - t464) * MDP(23) + ((qJD(5) * t278 + t271) * t430 + t539) * MDP(24)) * t426; -MDP(4) * t481 - MDP(5) * t465 + (-qJ(2) * t465 - t447) * MDP(6) + (-t400 * t431 + t439 * t428) * MDP(12) + (t400 * t428 + t439 * t431) * MDP(13) + (-t346 * t495 + t400 * t449 + t502 * t402) * MDP(14) + (-t350 * t495 + t379 * t400 - t501 * t402) * MDP(15) + (t317 * t379 - t318 * t449 + t501 * t346 + t502 * t350) * MDP(16) + (-t272 * t449 + t273 * t379 - t502 * t293 - t501 * t294 - t367 * t495 - t534) * MDP(17) + (-(-t379 * t427 - t430 * t449) * t391 - t305 * t495 + (t452 * t427 + t453 * t430) * t396) * MDP(23) + ((t379 * t430 - t427 * t449) * t391 - t450 * t495 + (-t453 * t427 + t452 * t430) * t396) * MDP(24); t431 * t428 * MDP(7) * t515 - t497 * MDP(8) * t515 + (-t428 * t462 + t479) * t522 + (-t431 * t462 - t480) * t496 - t486 + (t356 + (-t426 * t462 - t515) * t520 + (-t485 * t366 + t523 - t535) * t428 + t533) * MDP(12) + (g(3) * t512 + g(1) * t373 - g(2) * t371 - t357 * t402 + (t485 * t494 + t515) * t521 - t460) * MDP(13) + (t413 * t523 - g(1) * t360 + g(2) * t358 + t296 * t402 - t350 * t367 + (-t346 * t472 - t400 * t425) * pkin(3) + t272) * MDP(14) + (t414 * t523 + g(1) * t361 - g(2) * t359 - t297 * t402 + t346 * t367 + (-t350 * t472 + t400 * t423) * pkin(3) - t273) * MDP(15) + ((t294 + t296) * t350 + (-t293 + t297) * t346 + (t317 * t423 + t318 * t425) * pkin(3)) * MDP(16) + (-t293 * t296 - t294 * t297 + (t273 * t423 + t272 * t425 + (g(3) * t428 - t367 * t493) * t424 + t533) * pkin(3)) * MDP(17) + (t274 - t517) * MDP(20) + (-t275 - t518) * MDP(21) + (-t444 * t391 + (t284 * t430 - t285 * t427) * t396 - t335 * t305 + (t445 * t396 + t451) * qJD(5) + t531) * MDP(23) + (t445 * t391 - (t284 * t427 + t285 * t430) * t396 - t335 * t450 + (t444 * t396 - t505) * qJD(5) + t538) * MDP(24) + t537; (-t350 * t402 + t474) * MDP(14) + (t346 * t402 + t385) * MDP(15) + (-t346 ^ 2 - t350 ^ 2) * MDP(16) + (g(3) * t426 + t293 * t350 + t294 * t346 + t332) * MDP(17) + (t275 - t518) * MDP(23) + (t274 + t517) * MDP(24) + (-MDP(15) * t469 - t455 * MDP(17) + (-MDP(14) * t514 - t379 * MDP(15)) * t482) * t424; (t475 - t517) * MDP(20) + (t463 - t518) * MDP(21) + (t451 * t396 + t531) * MDP(23) + (-(-t280 * t427 + t505) * t396 + t538) * MDP(24) + (-MDP(20) * t516 - t450 * MDP(21) + t451 * MDP(23) - MDP(24) * t505) * qJD(5) + t537;];
tau = t1;
