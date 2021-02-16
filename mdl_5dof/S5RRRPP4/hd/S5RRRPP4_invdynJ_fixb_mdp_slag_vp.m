% Calculate vector of inverse dynamics joint torques for
% S5RRRPP4
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
%   see S5RRRPP4_convert_par2_MPV_fixb.m
% 
% Output:
% tau [5x1]
%   joint torques of inverse dynamics (contains inertial, gravitational Coriolis and centrifugal forces)

% Quelle: HybrDyn-Toolbox
% Datum: 2021-01-15 22:26
% Revision: 24b2e7d74a0c1a3b64fa2f8f5ad758691ad61af3 (2021-01-15)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tau = S5RRRPP4_invdynJ_fixb_mdp_slag_vp(qJ, qJD, qJDD, g, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(5,1),zeros(3,1),zeros(8,1),zeros(25,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRRPP4_invdynJ_fixb_mdp_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RRRPP4_invdynJ_fixb_mdp_slag_vp: qJD has to be [5x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [5 1]), ...
  'S5RRRPP4_invdynJ_fixb_mdp_slag_vp: qJDD has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RRRPP4_invdynJ_fixb_mdp_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RRRPP4_invdynJ_fixb_mdp_slag_vp: pkin has to be [8x1] (double)');
assert(isreal(MDP) && all(size(MDP) == [25 1]), ...
  'S5RRRPP4_invdynJ_fixb_mdp_slag_vp: MDP has to be [25x1] (double)'); 

%% Symbolic Calculation
% From invdyn_joint_fixb_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2021-01-15 22:24:59
% EndTime: 2021-01-15 22:25:09
% DurationCPUTime: 4.30s
% Computational Cost: add. (4335->391), mult. (10173->469), div. (0->0), fcn. (7093->12), ass. (0->185)
t424 = qJ(2) + qJ(3);
t415 = sin(t424);
t416 = cos(t424);
t429 = sin(qJ(1));
t432 = cos(qJ(1));
t474 = g(1) * t432 + g(2) * t429;
t536 = -g(3) * t416 + t415 * t474;
t430 = cos(qJ(3));
t431 = cos(qJ(2));
t497 = qJD(1) * t431;
t485 = t430 * t497;
t427 = sin(qJ(3));
t428 = sin(qJ(2));
t498 = qJD(1) * t428;
t486 = t427 * t498;
t350 = -t485 + t486;
t352 = -t427 * t497 - t430 * t498;
t425 = sin(pkin(8));
t426 = cos(pkin(8));
t321 = -t426 * t350 + t352 * t425;
t421 = qJD(2) + qJD(3);
t535 = t321 * t421;
t493 = qJD(1) * qJD(2);
t483 = t431 * t493;
t492 = qJDD(1) * t428;
t533 = t483 + t492;
t463 = -t350 * t425 - t426 * t352;
t532 = t463 ^ 2;
t526 = pkin(6) + pkin(7);
t417 = t431 * pkin(2);
t517 = pkin(1) + t417;
t346 = t352 * qJ(4);
t380 = t526 * t431;
t367 = qJD(1) * t380;
t353 = t427 * t367;
t379 = t526 * t428;
t365 = qJD(1) * t379;
t502 = -t430 * t365 - t353;
t316 = t346 + t502;
t510 = t367 * t430;
t461 = t365 * t427 - t510;
t515 = qJ(4) * t350;
t452 = t461 + t515;
t490 = pkin(2) * t425 * t427;
t495 = qJD(3) * t430;
t503 = -qJD(3) * t490 - t425 * t452 + (t495 * pkin(2) - t316) * t426;
t516 = qJD(2) * pkin(2);
t358 = -t365 + t516;
t477 = t430 * t358 - t353;
t311 = t346 + t477;
t501 = -t427 * t379 + t430 * t380;
t414 = pkin(8) + t424;
t400 = sin(t414);
t401 = cos(t414);
t470 = t401 * pkin(4) + t400 * qJ(5);
t531 = g(1) * t429 - g(2) * t432;
t471 = t517 * qJDD(1);
t361 = t427 * t431 + t428 * t430;
t330 = t421 * t361;
t491 = qJDD(1) * t431;
t469 = t427 * t492 - t430 * t491;
t314 = qJD(1) * t330 + t469;
t530 = pkin(3) * t314 + qJDD(4);
t419 = qJDD(2) + qJDD(3);
t529 = -t419 * pkin(4) + qJDD(5);
t378 = t517 * qJD(1);
t332 = pkin(3) * t350 + qJD(4) - t378;
t291 = -pkin(4) * t321 - qJ(5) * t463 + t332;
t313 = qJD(3) * t485 - t421 * t486 + t427 * t491 + t533 * t430;
t331 = qJDD(2) * pkin(2) - t526 * t533;
t484 = t428 * t493;
t333 = t526 * (-t484 + t491);
t462 = -t358 * t427 - t510;
t442 = qJD(3) * t462 + t430 * t331 - t427 * t333;
t278 = pkin(3) * t419 - qJ(4) * t313 + qJD(4) * t352 + t442;
t496 = qJD(3) * t427;
t527 = (qJD(3) * t358 + t333) * t430 + t427 * t331 - t367 * t496;
t281 = -qJ(4) * t314 - qJD(4) * t350 + t527;
t506 = -t426 * t278 + t425 * t281;
t457 = -g(3) * t401 + t474 * t400 - t506;
t441 = -t291 * t463 + t457 - t529;
t360 = t427 * t428 - t430 * t431;
t487 = qJD(2) * t526;
t366 = t428 * t487;
t368 = t431 * t487;
t454 = -t430 * t366 - t427 * t368 - t379 * t495 - t380 * t496;
t292 = -qJ(4) * t330 - qJD(4) * t360 + t454;
t329 = t421 * t360;
t438 = -t501 * qJD(3) + t366 * t427 - t368 * t430;
t436 = qJ(4) * t329 - qJD(4) * t361 + t438;
t273 = t426 * t292 + t425 * t436;
t318 = -qJ(4) * t360 + t501;
t460 = -t379 * t430 - t380 * t427;
t451 = -qJ(4) * t361 + t460;
t299 = t426 * t318 + t425 * t451;
t528 = t273 * t421 + t299 * t419 + t400 * t531;
t524 = pkin(3) * t352;
t523 = pkin(3) * t415;
t522 = pkin(4) * t400;
t312 = -t462 - t515;
t307 = t426 * t312;
t286 = t311 * t425 + t307;
t514 = t286 * t463;
t513 = t286 * t421;
t511 = t312 * t425;
t287 = t311 * t426 - t511;
t512 = t287 * t421;
t509 = t401 * t429;
t508 = t401 * t432;
t507 = t426 * t427;
t269 = t425 * t278 + t426 * t281;
t306 = pkin(3) * t421 + t311;
t285 = t425 * t306 + t307;
t505 = qJD(5) + t503;
t504 = -t316 * t425 + t426 * t452 + (t425 * t430 + t507) * qJD(3) * pkin(2);
t407 = pkin(2) * t430 + pkin(3);
t345 = pkin(2) * t507 + t425 * t407;
t404 = pkin(3) * t416;
t500 = t404 + t417;
t422 = t428 ^ 2;
t499 = -t431 ^ 2 + t422;
t494 = qJD(5) - t287;
t411 = t428 * t516;
t489 = t404 + t470;
t324 = pkin(3) * t330 + t411;
t369 = -pkin(2) * t428 - t523;
t482 = t369 - t522;
t480 = t504 * t421;
t479 = t504 * t463;
t288 = t313 * t425 + t426 * t314;
t335 = pkin(3) * t360 - t517;
t475 = -t522 - t523;
t398 = pkin(2) * t484;
t472 = t398 + t530;
t468 = g(1) * t508 + g(2) * t509 + g(3) * t400 - t269;
t284 = t306 * t426 - t511;
t282 = -pkin(4) * t421 + qJD(5) - t284;
t283 = qJ(5) * t421 + t285;
t466 = -t282 * t321 + t283 * t463;
t465 = t284 * t321 + t285 * t463;
t289 = t313 * t426 - t314 * t425;
t344 = t407 * t426 - t490;
t459 = -0.2e1 * pkin(1) * t493 - pkin(6) * qJDD(2);
t347 = t398 - t471;
t296 = pkin(4) * t463 - qJ(5) * t321 - t524;
t453 = -t321 * t332 + t468;
t449 = -t352 * t350 * MDP(11) + (t350 * t421 + t313) * MDP(13) + (-t469 + (-qJD(1) * t361 - t352) * t421) * MDP(14) + (-t350 ^ 2 + t352 ^ 2) * MDP(12) + t419 * MDP(15);
t409 = t419 * qJ(5);
t447 = t291 * t321 + t409 - t468;
t446 = -t332 * t463 + t457;
t433 = qJD(2) ^ 2;
t445 = 0.2e1 * qJDD(1) * pkin(1) - pkin(6) * t433 + t531;
t434 = qJD(1) ^ 2;
t444 = pkin(1) * t434 - pkin(6) * qJDD(1) + t474;
t272 = t292 * t425 - t426 * t436;
t298 = t318 * t425 - t426 * t451;
t443 = g(1) * t509 - g(2) * t508 - t272 * t421 - t298 * t419;
t440 = t272 * t463 + t273 * t321 - t299 * t288 + t289 * t298 - t474;
t439 = pkin(4) * t288 - qJ(5) * t289 - qJD(5) * t463 + t472;
t437 = g(3) * t415 - t378 * t350 + t474 * t416 - t527;
t435 = -t378 * t352 + t442 + t536;
t420 = -qJ(4) - t526;
t412 = t421 * qJD(5);
t410 = pkin(2) * t498;
t402 = -pkin(3) * t426 - pkin(4);
t399 = pkin(3) * t425 + qJ(5);
t371 = qJ(5) * t508;
t370 = qJ(5) * t509;
t364 = pkin(1) + t500;
t357 = t432 * t364;
t340 = -pkin(4) - t344;
t339 = qJ(5) + t345;
t334 = t410 - t524;
t327 = -t360 * t425 + t361 * t426;
t326 = t426 * t360 + t361 * t425;
t304 = -t329 * t426 - t330 * t425;
t303 = -t329 * t425 + t426 * t330;
t302 = t347 + t530;
t297 = pkin(4) * t326 - qJ(5) * t327 + t335;
t295 = t296 + t410;
t274 = pkin(4) * t303 - qJ(5) * t304 - qJD(5) * t327 + t324;
t270 = -t471 + t439;
t267 = t506 + t529;
t266 = t409 + t412 + t269;
t1 = [qJDD(1) * MDP(1) + t531 * MDP(2) + t474 * MDP(3) + (qJDD(1) * t422 + 0.2e1 * t428 * t483) * MDP(4) + 0.2e1 * (t428 * t491 - t493 * t499) * MDP(5) + (qJDD(2) * t428 + t431 * t433) * MDP(6) + (qJDD(2) * t431 - t428 * t433) * MDP(7) + (t428 * t459 + t431 * t445) * MDP(9) + (-t428 * t445 + t431 * t459) * MDP(10) + (t313 * t361 + t329 * t352) * MDP(11) + (-t313 * t360 - t314 * t361 + t329 * t350 + t330 * t352) * MDP(12) + (-t329 * t421 + t361 * t419) * MDP(13) + (-t330 * t421 - t360 * t419) * MDP(14) + (-t314 * t517 - t378 * t330 + t347 * t360 + t350 * t411 + t416 * t531 + t419 * t460 + t421 * t438) * MDP(16) + (-t313 * t517 + t378 * t329 + t347 * t361 - t352 * t411 - t415 * t531 - t419 * t501 - t421 * t454) * MDP(17) + (t288 * t335 + t302 * t326 + t303 * t332 - t321 * t324 + t443) * MDP(18) + (t289 * t335 + t302 * t327 + t304 * t332 + t324 * t463 - t528) * MDP(19) + (-t269 * t326 - t284 * t304 - t285 * t303 + t327 * t506 + t440) * MDP(20) + (t269 * t299 + t285 * t273 + t506 * t298 - t284 * t272 + t302 * t335 + t332 * t324 - g(1) * (-t364 * t429 - t420 * t432) - g(2) * (-t420 * t429 + t357)) * MDP(21) + (t270 * t326 - t274 * t321 + t288 * t297 + t291 * t303 + t443) * MDP(22) + (-t266 * t326 + t267 * t327 + t282 * t304 - t283 * t303 + t440) * MDP(23) + (-t270 * t327 - t274 * t463 - t289 * t297 - t291 * t304 + t528) * MDP(24) + (-g(2) * t357 + t266 * t299 + t267 * t298 + t270 * t297 + t282 * t272 + t283 * t273 + t291 * t274 + (g(1) * t420 - g(2) * t470) * t432 + (-g(1) * (-t364 - t470) + g(2) * t420) * t429) * MDP(25); (t321 * t334 + t344 * t419 + t446 - t480) * MDP(18) + (t295 * t463 + t339 * t419 + t421 * t505 + t412 + t447) * MDP(24) + (t266 * t339 + t267 * t340 - t291 * t295 - g(1) * (t432 * t482 + t371) - g(2) * (t429 * t482 + t370) - g(3) * (t417 + t489) + t505 * t283 + t504 * t282) * MDP(25) + t449 + (-g(3) * t431 + t428 * t444) * MDP(9) + (g(3) * t428 + t431 * t444) * MDP(10) + (t295 * t321 - t340 * t419 + t441 - t480) * MDP(22) + (-t461 * t421 + (-t350 * t498 + t419 * t430 - t421 * t496) * pkin(2) + t435) * MDP(16) + (t502 * t421 + (t352 * t498 - t419 * t427 - t421 * t495) * pkin(2) + t437) * MDP(17) + (-g(3) * t500 + t269 * t345 - t284 * t504 + t285 * t503 - t332 * t334 - t344 * t506 - t369 * t474) * MDP(21) + MDP(6) * t492 + MDP(7) * t491 + qJDD(2) * MDP(8) + (-t288 * t339 + t289 * t340 + t321 * t505 + t466 + t479) * MDP(23) + (-t288 * t345 - t289 * t344 + t321 * t503 + t465 + t479) * MDP(20) + (-t334 * t463 - t345 * t419 - t421 * t503 + t453) * MDP(19) + (-MDP(4) * t428 * t431 + MDP(5) * t499) * t434; (-t421 * t462 + t435) * MDP(16) + (t421 * t477 + t437) * MDP(17) + (t513 + (-t321 * t352 + t419 * t426) * pkin(3) + t446) * MDP(18) + (t512 + (t352 * t463 - t419 * t425) * pkin(3) + t453) * MDP(19) + (-t514 - t287 * t321 + (-t288 * t425 - t289 * t426) * pkin(3) + t465) * MDP(20) + (t284 * t286 - t285 * t287 + (t269 * t425 + t332 * t352 - t426 * t506 + t536) * pkin(3)) * MDP(21) + (t296 * t321 - t402 * t419 + t441 + t513) * MDP(22) + (-t288 * t399 + t289 * t402 + t321 * t494 + t466 - t514) * MDP(23) + (t296 * t463 + t399 * t419 + 0.2e1 * t412 + t447 - t512) * MDP(24) + (t266 * t399 + t267 * t402 - t291 * t296 - t282 * t286 - g(1) * (t432 * t475 + t371) - g(2) * (t429 * t475 + t370) - g(3) * t489 + t494 * t283) * MDP(25) + t449; (t284 * t463 - t285 * t321 + t472 - t531) * MDP(21) + (-t282 * t463 - t283 * t321 + t439 - t531) * MDP(25) - (MDP(21) + MDP(25)) * t471 + (MDP(20) + MDP(23)) * (-t321 ^ 2 - t532) + (MDP(18) + MDP(22)) * (t421 * t463 + t288) + (MDP(19) - MDP(24)) * (t289 + t535); (-t321 * t463 - t419) * MDP(22) + (t289 - t535) * MDP(23) + (-t421 ^ 2 - t532) * MDP(24) + (-t283 * t421 - t441) * MDP(25);];
tau = t1;
