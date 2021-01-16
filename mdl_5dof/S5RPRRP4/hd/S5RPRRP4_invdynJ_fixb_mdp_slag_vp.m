% Calculate vector of inverse dynamics joint torques for
% S5RPRRP4
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
%   pkin=[a2,a3,a4,a5,d1,d3,d4,theta2]';
% MDP [24x1]
%   Minimal dynamic parameter vector (fixed base model)
%   see S5RPRRP4_convert_par2_MPV_fixb.m
% 
% Output:
% tau [5x1]
%   joint torques of inverse dynamics (contains inertial, gravitational Coriolis and centrifugal forces)

% Quelle: HybrDyn-Toolbox
% Datum: 2021-01-15 12:57
% Revision: d12c3222fdeb2c5f3b3c8fa5751e113be2fc3aae (2021-01-15)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tau = S5RPRRP4_invdynJ_fixb_mdp_slag_vp(qJ, qJD, qJDD, g, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(5,1),zeros(3,1),zeros(8,1),zeros(24,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPRRP4_invdynJ_fixb_mdp_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RPRRP4_invdynJ_fixb_mdp_slag_vp: qJD has to be [5x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [5 1]), ...
  'S5RPRRP4_invdynJ_fixb_mdp_slag_vp: qJDD has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RPRRP4_invdynJ_fixb_mdp_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RPRRP4_invdynJ_fixb_mdp_slag_vp: pkin has to be [8x1] (double)');
assert(isreal(MDP) && all(size(MDP) == [24 1]), ...
  'S5RPRRP4_invdynJ_fixb_mdp_slag_vp: MDP has to be [24x1] (double)'); 

%% Symbolic Calculation
% From invdyn_joint_fixb_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2021-01-15 12:56:47
% EndTime: 2021-01-15 12:56:58
% DurationCPUTime: 4.65s
% Computational Cost: add. (2862->386), mult. (6825->494), div. (0->0), fcn. (4829->10), ass. (0->201)
t412 = cos(qJ(3));
t464 = qJDD(1) * t412;
t406 = sin(pkin(8));
t411 = cos(qJ(4));
t535 = t406 * t411;
t362 = t464 * t535;
t408 = sin(qJ(4));
t409 = sin(qJ(3));
t465 = qJDD(1) * t409;
t503 = t408 * t412;
t353 = t409 * t411 + t503;
t463 = qJD(3) + qJD(4);
t537 = t463 * t353;
t288 = (qJD(1) * t537 + t408 * t465) * t406 - t362;
t401 = t406 ^ 2;
t527 = 0.2e1 * t401;
t513 = qJDD(1) * pkin(1);
t389 = qJDD(2) - t513;
t410 = sin(qJ(1));
t413 = cos(qJ(1));
t532 = -g(2) * t413 - g(3) * t410;
t536 = t389 - t532;
t407 = cos(pkin(8));
t469 = qJD(1) * qJD(2);
t470 = qJ(2) * qJDD(1);
t429 = t469 + t470;
t534 = t429 * t407;
t360 = -pkin(2) * t407 - pkin(6) * t406 - pkin(1);
t349 = t412 * t360;
t507 = t406 * t412;
t462 = pkin(7) * t507;
t518 = qJ(2) * t409;
t312 = -t462 + t349 + (-pkin(3) - t518) * t407;
t517 = qJ(2) * t412;
t381 = t407 * t517;
t489 = t409 * t360 + t381;
t508 = t406 * t409;
t319 = -pkin(7) * t508 + t489;
t494 = t408 * t312 + t411 * t319;
t481 = qJD(1) * t406;
t459 = t409 * t481;
t438 = t408 * t459;
t479 = qJD(1) * t412;
t458 = t406 * t479;
t326 = t411 * t458 - t438;
t320 = t326 * qJ(5);
t339 = t360 * qJD(1) + qJD(2);
t331 = t412 * t339;
t426 = -t407 * t518 - t462;
t308 = t426 * qJD(1) + t331;
t480 = qJD(1) * t407;
t382 = -qJD(3) + t480;
t294 = -pkin(3) * t382 + t308;
t309 = -pkin(7) * t459 + qJD(1) * t381 + t409 * t339;
t301 = t408 * t309;
t450 = t411 * t294 - t301;
t274 = -t320 + t450;
t466 = qJDD(1) * t407;
t380 = -qJDD(3) + t466;
t367 = -qJDD(4) + t380;
t474 = qJD(4) * t411;
t373 = -qJD(4) + t382;
t526 = pkin(3) * t373;
t533 = t408 * pkin(3) * t367 + t474 * t526;
t531 = -MDP(12) * t409 - MDP(13) * t412;
t361 = t367 * pkin(4);
t516 = qJ(5) * t288;
t530 = t516 - t361;
t405 = qJ(3) + qJ(4);
t392 = sin(t405);
t393 = cos(t405);
t506 = t407 * t410;
t332 = -t392 * t506 - t393 * t413;
t497 = t413 * t392;
t500 = t410 * t393;
t334 = t407 * t497 - t500;
t524 = g(1) * t406;
t529 = -g(2) * t332 - g(3) * t334 + t392 * t524;
t528 = t326 ^ 2;
t402 = t407 ^ 2;
t525 = pkin(7) * t406;
t522 = g(2) * t410;
t520 = g(3) * t413;
t519 = MDP(9) * t406;
t445 = t463 * t412;
t453 = t408 * t464;
t490 = t463 * t438;
t289 = (t453 + (qJD(1) * t445 + t465) * t411) * t406 - t490;
t515 = qJ(5) * t289;
t427 = qJD(1) * t353;
t323 = t406 * t427;
t514 = qJ(5) * t323;
t512 = t323 * t373;
t511 = t326 * t373;
t510 = (-qJ(5) - pkin(7) - pkin(6)) * t406;
t414 = qJD(1) ^ 2;
t509 = t401 * t414;
t505 = t407 * t413;
t504 = t408 * t409;
t502 = t409 * t410;
t501 = t409 * t413;
t499 = t410 * t412;
t303 = t411 * t309;
t498 = t412 * t413;
t271 = -pkin(4) * t373 + t274;
t496 = -t274 + t271;
t495 = t411 * t308 - t301;
t352 = t411 * t412 - t504;
t493 = (t463 - t480) * t352;
t492 = t407 * t427 - t537;
t476 = qJD(3) * t412;
t478 = qJD(2) * t407;
t491 = t360 * t476 + t412 * t478;
t340 = pkin(3) * t459 + qJ(2) * t481;
t359 = t412 * pkin(3) + pkin(4) * t393;
t376 = t406 * pkin(3) * t476;
t346 = t406 * qJD(2) + t376;
t383 = pkin(3) * t508;
t351 = t406 * qJ(2) + t383;
t488 = t413 * pkin(1) + t410 * qJ(2);
t486 = t401 + t402;
t404 = t412 ^ 2;
t485 = t409 ^ 2 - t404;
t484 = MDP(10) * t406;
t477 = qJD(3) * t339;
t475 = qJD(4) * t408;
t357 = t367 * MDP(18);
t473 = t380 * MDP(11);
t472 = qJD(3) + t382;
t452 = -pkin(4) * t323 - qJD(5);
t306 = -t452 + t340;
t471 = qJD(5) + t306;
t468 = qJD(1) * qJD(3);
t467 = qJDD(1) * t406;
t460 = qJ(2) * qJD(3) * t407;
t457 = qJ(2) * t466;
t456 = t412 * t469;
t455 = t412 * t468;
t338 = t360 * qJDD(1) + qJDD(2);
t330 = t412 * t338;
t437 = qJD(1) * t460;
t283 = -pkin(3) * t380 + t330 + (-pkin(7) * t467 - t437) * t412 + (-t457 - t477 + (qJD(3) * t525 - t478) * qJD(1)) * t409;
t440 = t409 * t338 + t339 * t476 + t407 * t456 + t412 * t457;
t422 = -t409 * t437 + t440;
t287 = (-t455 - t465) * t525 + t422;
t451 = t411 * t283 - t408 * t287;
t304 = t426 * qJD(3) + t491;
t305 = -t409 * t478 + (-t381 + (-t360 + t525) * t409) * qJD(3);
t449 = -t304 * t408 + t411 * t305;
t448 = -t308 * t408 - t303;
t447 = t411 * t312 - t319 * t408;
t446 = qJD(1) * t472;
t444 = pkin(3) * t458;
t443 = t380 + t466;
t442 = MDP(22) * t463;
t441 = -t408 * t283 - t411 * t287 - t294 * t474 + t309 * t475;
t315 = qJ(2) * t467 + qJD(1) * t376 + qJDD(1) * t383 + t406 * t469;
t439 = 0.2e1 * t486;
t436 = g(2) * t334 - g(3) * t332;
t333 = t407 * t500 - t497;
t335 = t392 * t410 + t393 * t505;
t435 = -g(2) * t335 - g(3) * t333;
t433 = -t520 + t522;
t432 = -t294 * t408 - t303;
t430 = qJD(3) * (t382 + t480);
t428 = t411 * t304 + t408 * t305 + t312 * t474 - t319 * t475;
t278 = pkin(4) * t289 + qJDD(5) + t315;
t321 = t323 ^ 2;
t425 = t326 * t323 * MDP(14) + (-t288 - t512) * MDP(16) + (-t511 + (-t453 + (-t463 * t479 - t465) * t411) * t406 + t490) * MDP(17) + (-t321 + t528) * MDP(15) - t357;
t424 = t439 * t469 + t520;
t423 = g(2) * t333 - g(3) * t335 + t393 * t524 + t441;
t421 = t432 * qJD(4) + t451;
t420 = t340 * t323 + t423;
t418 = t471 * t323 + t423 + t515;
t417 = t421 + t529;
t416 = -t340 * t326 + t417;
t395 = t410 * pkin(1);
t388 = pkin(3) * t411 + pkin(4);
t358 = pkin(3) * t409 + pkin(4) * t392;
t354 = pkin(2) + t359;
t344 = t407 * t498 + t502;
t343 = t407 * t501 - t499;
t342 = t407 * t499 - t501;
t341 = -t407 * t502 - t498;
t337 = t352 * t406;
t336 = t353 * t406;
t316 = pkin(4) * t326 + t444;
t313 = pkin(4) * t336 + t351;
t300 = -t463 * t406 * t504 + t445 * t535;
t299 = t537 * t406;
t291 = pkin(4) * t300 + t346;
t286 = -qJ(5) * t336 + t494;
t284 = -pkin(4) * t407 - qJ(5) * t337 + t447;
t277 = -t320 + t495;
t276 = t448 + t514;
t275 = -t432 - t514;
t270 = qJ(5) * t299 - qJD(4) * t494 - qJD(5) * t337 + t449;
t269 = -qJ(5) * t300 - qJD(5) * t336 + t428;
t268 = -qJD(5) * t323 - t441 - t515;
t267 = -qJD(5) * t326 + t421 + t530;
t1 = [qJDD(1) * MDP(1) + t532 * MDP(2) + t433 * MDP(3) + (t513 - t536) * t407 * MDP(4) + (t439 * t470 + t424 - t522) * MDP(5) + (-t389 * pkin(1) - g(2) * t488 - g(3) * t395 + (t486 * t470 + t424) * qJ(2)) * MDP(6) + (qJDD(1) * t404 - 0.2e1 * t409 * t455) * t401 * MDP(7) + (-t409 * t464 + t485 * t468) * MDP(8) * t527 + (t409 * t430 - t443 * t412) * t519 + (t443 * t409 + t412 * t430) * t484 + t407 * t473 + (-g(2) * t344 - g(3) * t342 - t330 * t407 - t349 * t380 + (t382 * t407 + (t527 + t402) * qJD(1)) * qJ(2) * t476 + (qJD(3) * t360 * t382 + t429 * t527 + (qJ(2) * t380 + qJD(2) * t382 + t477 + t534) * t407) * t409) * MDP(12) + ((-t409 * t460 + t491) * t382 + t489 * t380 + t422 * t407 + g(2) * t343 - g(3) * t341 + (t456 + (-t409 * t468 + t464) * qJ(2)) * t527) * MDP(13) + (-t288 * t337 - t299 * t326) * MDP(14) + (t288 * t336 - t289 * t337 + t299 * t323 - t300 * t326) * MDP(15) + (t288 * t407 + t299 * t373 - t337 * t367) * MDP(16) + (t289 * t407 + t300 * t373 + t336 * t367) * MDP(17) + t407 * t357 + (-t449 * t373 - t447 * t367 - t451 * t407 + t346 * t323 + t351 * t289 + t315 * t336 + t340 * t300 + (t373 * t494 - t432 * t407) * qJD(4) + t435) * MDP(19) + (-t351 * t288 - t340 * t299 + t315 * t337 + t346 * t326 + t494 * t367 + t428 * t373 - t441 * t407 + t436) * MDP(20) + (-t267 * t407 - t270 * t373 + t278 * t336 - t284 * t367 + t289 * t313 + t291 * t323 + t300 * t306 + t435) * MDP(21) + (t268 * t407 + t269 * t373 + t278 * t337 + t286 * t367 - t288 * t313 + t291 * t326 - t299 * t306 + t436) * MDP(22) + (-t267 * t337 - t268 * t336 - t269 * t323 - t270 * t326 + t271 * t299 - t275 * t300 + t284 * t288 - t286 * t289 + t406 * t532) * MDP(23) + (t268 * t286 + t275 * t269 + t267 * t284 + t271 * t270 + t278 * t313 + t306 * t291 - g(2) * (t358 * t410 + t488) - g(3) * (t354 * t506 - t410 * t510 + t395) + (-g(2) * (t354 * t407 - t510) - g(3) * (-qJ(2) - t358)) * t413) * MDP(24); -MDP(4) * t466 + t536 * MDP(6) + (t288 * t352 - t289 * t353 - t493 * t323 - t492 * t326) * MDP(23) + (t267 * t352 + t268 * t353 + t492 * t271 + t493 * t275 - t306 * t481 - t532) * MDP(24) + (-MDP(12) * t412 + MDP(13) * t409) * t380 + (MDP(19) + MDP(21)) * (-t323 * t481 - t352 * t367 - t492 * t373) + (MDP(20) + MDP(22)) * (-t326 * t481 + t353 * t367 + t493 * t373) + (-t402 * MDP(5) + (-MDP(5) + t531) * t401 - t486 * MDP(6) * qJ(2)) * t414 + t531 * t382 ^ 2; t412 * t409 * MDP(7) * t509 - t485 * MDP(8) * t509 + (-t409 * t446 + t464) * t519 + (-t412 * t446 - t465) * t484 - t473 + (-g(2) * t341 - g(3) * t343 + t330 + (-t407 * t446 - t509) * t517 + (-t472 * t339 + t524 - t534) * t409) * MDP(12) + (g(1) * t507 + g(2) * t342 - g(3) * t344 - t331 * t382 + (t472 * t480 + t509) * t518 - t440) * MDP(13) + (t448 * t373 + (-t323 * t458 - t367 * t411 + t373 * t475) * pkin(3) + t416) * MDP(19) + (-t326 * t444 - t495 * t373 + t420 + t533) * MDP(20) + (t276 * t373 - t316 * t323 - t367 * t388 - t471 * t326 + (-t303 + (-t294 + t526) * t408) * qJD(4) + t451 + t529 + t530) * MDP(21) + (-t277 * t373 - t316 * t326 + t418 + t533) * MDP(22) + (t288 * t388 + (t275 + t276) * t326 + (-t271 + t277) * t323 + (-t289 * t408 + (-t323 * t411 + t326 * t408) * qJD(4)) * pkin(3)) * MDP(23) + (t267 * t388 - t275 * t277 - t271 * t276 - t306 * t316 + t358 * t524 - g(2) * (-t358 * t506 - t359 * t413) - g(3) * (t358 * t505 - t359 * t410) + (t268 * t408 + (-t271 * t408 + t275 * t411) * qJD(4)) * pkin(3)) * MDP(24) + t425; (t432 * t373 + t416) * MDP(19) + (-t450 * t373 + t420) * MDP(20) + (t516 - t275 * t373 - 0.2e1 * t361 + (-t306 + t452) * t326 + t417) * MDP(21) + (-pkin(4) * t528 - t274 * t373 + t418) * MDP(22) + (pkin(4) * t288 - t496 * t323) * MDP(23) + (t496 * t275 + (-t306 * t326 + t267 + t529) * pkin(4)) * MDP(24) + t425; (-t490 - t511) * MDP(21) + (t512 + t362) * MDP(22) + (-t321 - t528) * MDP(23) + (g(1) * t407 + t271 * t326 + t275 * t323 + t278) * MDP(24) + (-t433 * MDP(24) + (t353 * MDP(21) - MDP(22) * t504) * qJDD(1) + (-t442 * t503 + (MDP(21) * t445 - t409 * t442) * t411) * qJD(1)) * t406;];
tau = t1;
