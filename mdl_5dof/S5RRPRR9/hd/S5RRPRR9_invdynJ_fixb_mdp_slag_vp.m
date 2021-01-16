% Calculate vector of inverse dynamics joint torques for
% S5RRPRR9
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
%   pkin=[a2,a3,a4,a5,d1,d2,d4,d5,theta3]';
% MDP [28x1]
%   Minimal dynamic parameter vector (fixed base model)
%   see S5RRPRR9_convert_par2_MPV_fixb.m
% 
% Output:
% tau [5x1]
%   joint torques of inverse dynamics (contains inertial, gravitational Coriolis and centrifugal forces)

% Quelle: HybrDyn-Toolbox
% Datum: 2021-01-15 21:48
% Revision: 24b2e7d74a0c1a3b64fa2f8f5ad758691ad61af3 (2021-01-15)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tau = S5RRPRR9_invdynJ_fixb_mdp_slag_vp(qJ, qJD, qJDD, g, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(5,1),zeros(3,1),zeros(9,1),zeros(28,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRPRR9_invdynJ_fixb_mdp_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RRPRR9_invdynJ_fixb_mdp_slag_vp: qJD has to be [5x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [5 1]), ...
  'S5RRPRR9_invdynJ_fixb_mdp_slag_vp: qJDD has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RRPRR9_invdynJ_fixb_mdp_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5RRPRR9_invdynJ_fixb_mdp_slag_vp: pkin has to be [9x1] (double)');
assert(isreal(MDP) && all(size(MDP) == [28 1]), ...
  'S5RRPRR9_invdynJ_fixb_mdp_slag_vp: MDP has to be [28x1] (double)'); 

%% Symbolic Calculation
% From invdyn_joint_fixb_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2021-01-15 21:47:50
% EndTime: 2021-01-15 21:48:07
% DurationCPUTime: 6.58s
% Computational Cost: add. (3852->461), mult. (9015->607), div. (0->0), fcn. (6781->14), ass. (0->205)
t441 = cos(qJ(2));
t529 = cos(pkin(9));
t476 = t529 * t441;
t414 = qJD(1) * t476;
t433 = sin(pkin(9));
t437 = sin(qJ(2));
t496 = qJD(1) * t437;
t388 = t433 * t496 - t414;
t537 = qJD(4) + qJD(5);
t546 = t388 + t537;
t434 = -qJ(3) - pkin(6);
t480 = qJD(2) * t434;
t451 = -qJD(3) * t437 + t441 * t480;
t482 = t434 * t437;
t348 = qJDD(2) * pkin(2) + qJD(1) * t451 + qJDD(1) * t482;
t387 = qJD(3) * t441 + t437 * t480;
t411 = t434 * t441;
t356 = qJD(1) * t387 - qJDD(1) * t411;
t301 = t529 * t348 - t433 * t356;
t299 = -qJDD(2) * pkin(3) - t301;
t380 = qJD(4) + t388;
t417 = pkin(2) * t433 + pkin(7);
t429 = qJ(2) + pkin(9);
t422 = sin(t429);
t423 = cos(t429);
t438 = sin(qJ(1));
t442 = cos(qJ(1));
t468 = g(1) * t442 + g(2) * t438;
t450 = -g(3) * t423 + t422 * t468;
t545 = -qJD(4) * t417 * t380 - t299 + t450;
t477 = t529 * t437;
t401 = t433 * t441 + t477;
t391 = t401 * qJD(1);
t436 = sin(qJ(4));
t440 = cos(qJ(4));
t491 = t440 * qJD(2);
t367 = t391 * t436 - t491;
t439 = cos(qJ(5));
t369 = qJD(2) * t436 + t391 * t440;
t435 = sin(qJ(5));
t518 = t369 * t435;
t310 = t439 * t367 + t518;
t377 = qJD(5) + t380;
t544 = t310 * t377;
t461 = t367 * t435 - t439 * t369;
t543 = t377 * t461;
t508 = t435 * t440;
t404 = t436 * t439 + t508;
t499 = t546 * t404;
t495 = qJD(4) * t436;
t517 = t388 * t436;
t542 = t495 + t517;
t421 = pkin(2) * t441 + pkin(1);
t409 = -qJD(1) * t421 + qJD(3);
t319 = pkin(3) * t388 - pkin(7) * t391 + t409;
t405 = qJD(1) * t482;
t530 = qJD(2) * pkin(2);
t397 = t405 + t530;
t406 = qJD(1) * t411;
t478 = t529 * t406;
t353 = t433 * t397 - t478;
t342 = qJD(2) * pkin(7) + t353;
t294 = t319 * t436 + t342 * t440;
t288 = -pkin(8) * t367 + t294;
t493 = qJD(5) * t435;
t285 = t288 * t493;
t394 = t433 * t406;
t352 = t397 * t529 + t394;
t341 = -qJD(2) * pkin(3) - t352;
t308 = t367 * pkin(4) + t341;
t432 = qJ(4) + qJ(5);
t428 = cos(t432);
t510 = t428 * t438;
t427 = sin(t432);
t511 = t427 * t442;
t371 = -t423 * t510 + t511;
t509 = t428 * t442;
t512 = t427 * t438;
t373 = t423 * t509 + t512;
t534 = g(3) * t422;
t541 = g(1) * t373 - g(2) * t371 + t308 * t310 + t428 * t534 + t285;
t370 = t423 * t512 + t509;
t372 = -t423 * t511 + t510;
t390 = t401 * qJD(2);
t489 = qJDD(1) * t437;
t462 = -qJDD(1) * t476 + t433 * t489;
t354 = qJD(1) * t390 + t462;
t490 = qJD(1) * qJD(2);
t481 = t437 * t490;
t445 = qJDD(1) * t401 - t433 * t481;
t355 = qJD(2) * t414 + t445;
t379 = pkin(2) * t481 - qJDD(1) * t421 + qJDD(3);
t295 = pkin(3) * t354 - pkin(7) * t355 + t379;
t292 = t440 * t295;
t302 = t433 * t348 + t529 * t356;
t300 = qJDD(2) * pkin(7) + t302;
t306 = qJD(4) * t491 + t436 * qJDD(2) + t440 * t355 - t391 * t495;
t349 = qJDD(4) + t354;
t273 = pkin(4) * t349 - pkin(8) * t306 - qJD(4) * t294 - t300 * t436 + t292;
t307 = t369 * qJD(4) - t440 * qJDD(2) + t355 * t436;
t494 = qJD(4) * t440;
t454 = t436 * t295 + t440 * t300 + t319 * t494 - t342 * t495;
t274 = -pkin(8) * t307 + t454;
t475 = t439 * t273 - t435 * t274;
t540 = -g(1) * t372 + g(2) * t370 + t308 * t461 + t427 * t534 + t475;
t344 = qJDD(5) + t349;
t539 = t344 * MDP(26) + (-t310 ^ 2 + t461 ^ 2) * MDP(23) - t310 * MDP(22) * t461;
t337 = t404 * t401;
t403 = t435 * t436 - t439 * t440;
t500 = t546 * t403;
t536 = -t344 * t404 + t377 * t500;
t474 = t306 * t435 + t439 * t307;
t279 = -qJD(5) * t461 + t474;
t535 = pkin(2) * t437;
t532 = g(3) * t441;
t531 = pkin(8) + t417;
t293 = t440 * t319 - t342 * t436;
t287 = -pkin(8) * t369 + t293;
t282 = pkin(4) * t380 + t287;
t528 = t282 * t439;
t527 = t288 * t439;
t526 = t306 * t436;
t525 = t310 * t391;
t524 = t461 * t391;
t522 = t367 * t380;
t521 = t367 * t391;
t520 = t369 * t380;
t519 = t369 * t391;
t456 = -t433 * t437 + t476;
t393 = t456 * qJD(2);
t516 = t393 * t436;
t515 = t393 * t440;
t514 = t401 * t436;
t513 = t401 * t440;
t507 = t436 * t349;
t506 = t436 * t438;
t505 = t436 * t442;
t504 = t438 * t440;
t336 = t440 * t349;
t366 = -t411 * t529 + t433 * t482;
t359 = t440 * t366;
t503 = t440 * t442;
t330 = pkin(2) * t496 + pkin(3) * t391 + pkin(7) * t388;
t358 = t405 * t529 + t394;
t502 = t436 * t330 + t440 * t358;
t351 = -pkin(3) * t456 - pkin(7) * t401 - t421;
t501 = t436 * t351 + t359;
t430 = t437 ^ 2;
t498 = -t441 ^ 2 + t430;
t492 = qJD(5) * t439;
t488 = qJDD(1) * t441;
t487 = t437 * t530;
t486 = t439 * t306 - t435 * t307 - t367 * t492;
t485 = t401 * t495;
t484 = t401 * t494;
t479 = qJD(4) * t531;
t328 = t387 * t433 - t529 * t451;
t357 = t405 * t433 - t478;
t365 = -t411 * t433 - t434 * t477;
t472 = t380 * t440;
t471 = -qJD(4) * t319 - t300;
t470 = qJD(5) * t282 + t274;
t418 = -pkin(2) * t529 - pkin(3);
t469 = pkin(4) * t542 - t357;
t467 = g(1) * t438 - g(2) * t442;
t466 = -t403 * t344 - t377 * t499;
t465 = -t342 * t494 + t292;
t323 = t440 * t330;
t399 = t531 * t440;
t464 = pkin(4) * t391 + qJD(5) * t399 - t358 * t436 + t323 + (pkin(8) * t388 + t479) * t440;
t398 = t531 * t436;
t463 = pkin(8) * t517 + qJD(5) * t398 + t436 * t479 + t502;
t276 = t282 * t435 + t527;
t460 = -t380 * t542 + t336;
t459 = -0.2e1 * pkin(1) * t490 - pkin(6) * qJDD(2);
t458 = t484 + t516;
t457 = -t485 + t515;
t329 = t387 * t529 + t433 * t451;
t331 = pkin(3) * t390 - pkin(7) * t393 + t487;
t453 = t440 * t329 + t436 * t331 + t351 * t494 - t366 * t495;
t278 = -t369 * t493 + t486;
t452 = t341 * t380 - t417 * t349;
t443 = qJD(2) ^ 2;
t447 = 0.2e1 * qJDD(1) * pkin(1) - pkin(6) * t443 + t467;
t444 = qJD(1) ^ 2;
t446 = pkin(1) * t444 - pkin(6) * qJDD(1) + t468;
t410 = -t440 * pkin(4) + t418;
t384 = t423 * t503 + t506;
t383 = -t423 * t505 + t504;
t382 = -t423 * t504 + t505;
t381 = t423 * t506 + t503;
t340 = t440 * t351;
t338 = t403 * t401;
t327 = pkin(4) * t514 + t365;
t324 = t440 * t331;
t298 = pkin(4) * t458 + t328;
t296 = -pkin(8) * t514 + t501;
t289 = -pkin(4) * t456 - pkin(8) * t513 - t366 * t436 + t340;
t284 = t393 * t508 - t435 * t485 - t493 * t514 + (t513 * t537 + t516) * t439;
t283 = -t337 * t537 - t403 * t393;
t281 = pkin(4) * t307 + t299;
t280 = -pkin(8) * t458 + t453;
t277 = -pkin(8) * t515 + pkin(4) * t390 - t329 * t436 + t324 + (-t359 + (pkin(8) * t401 - t351) * t436) * qJD(4);
t275 = -t288 * t435 + t528;
t1 = [(t437 * t459 + t441 * t447) * MDP(9) + (-t437 * t447 + t441 * t459) * MDP(10) + (-t278 * t338 - t283 * t461) * MDP(22) + (-t278 * t337 + t279 * t338 - t283 * t310 + t284 * t461) * MDP(23) + (-g(1) * t370 - g(2) * t372 - t276 * t390 + t327 * t278 - t281 * t338 + t308 * t283 - t285 * t456 - t298 * t461 + (-(-qJD(5) * t296 + t277) * t377 - t289 * t344 + t273 * t456) * t435 + (-(qJD(5) * t289 + t280) * t377 - t296 * t344 + t470 * t456) * t439) * MDP(28) + (-t278 * t456 + t283 * t377 - t338 * t344 - t390 * t461) * MDP(24) + (-t301 * t401 + t302 * t456 + t328 * t391 - t329 * t388 - t352 * t393 - t353 * t390 - t354 * t366 + t355 * t365 - t468) * MDP(13) + ((t277 * t439 - t280 * t435) * t377 + (t289 * t439 - t296 * t435) * t344 - t475 * t456 + t275 * t390 + t298 * t310 + t327 * t279 + t281 * t337 + t308 * t284 - g(1) * t371 - g(2) * t373 + ((-t289 * t435 - t296 * t439) * t377 + t276 * t456) * qJD(5)) * MDP(27) + ((-t366 * t494 + t324) * t380 + t340 * t349 - t465 * t456 + t293 * t390 + t328 * t367 + t365 * t307 + t341 * t484 - g(1) * t382 - g(2) * t384 + ((-qJD(4) * t351 - t329) * t380 - t366 * t349 - t471 * t456 + t299 * t401 + t341 * t393) * t436) * MDP(20) + (-t306 * t456 + t336 * t401 + t369 * t390 + t380 * t457) * MDP(17) + (t307 * t456 - t367 * t390 - t380 * t458 - t401 * t507) * MDP(18) + (-g(1) * t381 - g(2) * t383 - t294 * t390 + t299 * t513 + t365 * t306 + t328 * t369 + t341 * t457 - t349 * t501 - t380 * t453 + t454 * t456) * MDP(21) + (-qJDD(2) * t365 - t354 * t421 - t379 * t456 + t390 * t409 + t467 * t423 + (t388 * t535 - t328) * qJD(2)) * MDP(11) + (-t344 * t456 + t377 * t390) * MDP(26) + (t279 * t456 - t284 * t377 - t310 * t390 - t337 * t344) * MDP(25) + (-t349 * t456 + t380 * t390) * MDP(19) + (-qJDD(2) * t366 - t355 * t421 + t379 * t401 + t393 * t409 - t467 * t422 + (t391 * t535 - t329) * qJD(2)) * MDP(12) + qJDD(1) * MDP(1) + (qJDD(2) * t437 + t441 * t443) * MDP(6) + (qJDD(2) * t441 - t437 * t443) * MDP(7) + t467 * MDP(2) + t468 * MDP(3) + (qJDD(1) * t430 + 0.2e1 * t441 * t481) * MDP(4) + (t302 * t366 + t353 * t329 - t301 * t365 - t352 * t328 - t379 * t421 + t409 * t487 - g(1) * (-t421 * t438 - t434 * t442) - g(2) * (t421 * t442 - t434 * t438)) * MDP(14) + 0.2e1 * (t437 * t488 - t490 * t498) * MDP(5) + (t306 * t513 + t369 * t457) * MDP(15) + ((-t367 * t440 - t369 * t436) * t393 + (-t526 - t307 * t440 + (t367 * t436 - t369 * t440) * qJD(4)) * t401) * MDP(16); (g(3) * t437 + t441 * t446) * MDP(10) + MDP(7) * t488 + MDP(6) * t489 + qJDD(2) * MDP(8) - t380 * t391 * MDP(19) - t377 * t391 * MDP(26) + (-t293 * t391 + t418 * t307 - t323 * t380 - t357 * t367 + (t358 * t380 + t452) * t436 + t545 * t440) * MDP(20) + ((-t398 * t439 - t399 * t435) * t344 + t410 * t279 + t281 * t403 - t275 * t391 + (t435 * t463 - t439 * t464) * t377 + t469 * t310 + t499 * t308 + t450 * t428) * MDP(27) + (t278 * t404 + t461 * t500) * MDP(22) + (-(-t398 * t435 + t399 * t439) * t344 + t410 * t278 + t281 * t404 + t276 * t391 + (t435 * t464 + t439 * t463) * t377 - t469 * t461 - t500 * t308 - t450 * t427) * MDP(28) + (-t278 * t403 - t279 * t404 + t310 * t500 + t461 * t499) * MDP(23) + (t294 * t391 + t418 * t306 - t357 * t369 + t502 * t380 - t436 * t545 + t452 * t440) * MDP(21) + (t380 * t472 + t507 - t519) * MDP(17) + (t460 + t521) * MDP(18) + ((t306 - t522) * t440 + (-t307 - t520) * t436) * MDP(16) + (t524 - t536) * MDP(24) + (t466 + t525) * MDP(25) + (t369 * t472 + t526) * MDP(15) + ((t353 - t357) * t391 + (-t352 + t358) * t388 + (-t354 * t433 - t355 * t529) * pkin(2)) * MDP(13) + (t357 * qJD(2) - t409 * t391 + (qJDD(2) * t529 - t388 * t496) * pkin(2) + t450 + t301) * MDP(11) + (t437 * t446 - t532) * MDP(9) + (t352 * t357 - t353 * t358 + (t529 * t301 - t532 + t302 * t433 + (-qJD(1) * t409 + t468) * t437) * pkin(2)) * MDP(14) + (t534 + qJD(2) * t358 + t388 * t409 + t468 * t423 + (-qJDD(2) * t433 - t391 * t496) * pkin(2) - t302) * MDP(12) + (-MDP(4) * t437 * t441 + MDP(5) * t498) * t444; (0.2e1 * t391 * qJD(2) + t462) * MDP(11) + ((t414 - t388) * qJD(2) + t445) * MDP(12) + (-t388 ^ 2 - t391 ^ 2) * MDP(13) + (t352 * t391 + t353 * t388 + t379 - t467) * MDP(14) + (t460 - t521) * MDP(20) + (-t380 ^ 2 * t440 - t507 - t519) * MDP(21) + (t466 - t525) * MDP(27) + (t524 + t536) * MDP(28); t369 * t367 * MDP(15) + (-t367 ^ 2 + t369 ^ 2) * MDP(16) + (t306 + t522) * MDP(17) + (-t307 + t520) * MDP(18) + t349 * MDP(19) + (-g(1) * t383 + g(2) * t381 + t294 * t380 - t341 * t369 + (t471 + t534) * t436 + t465) * MDP(20) + (g(1) * t384 - g(2) * t382 + t293 * t380 + t341 * t367 + t440 * t534 - t454) * MDP(21) + (t278 + t544) * MDP(24) + (-t279 - t543) * MDP(25) + (-(-t287 * t435 - t527) * t377 - t276 * qJD(5) + (-t310 * t369 + t344 * t439 - t377 * t493) * pkin(4) + t540) * MDP(27) + ((-t288 * t377 - t273) * t435 + (t287 * t377 - t470) * t439 + (-t344 * t435 + t369 * t461 - t377 * t492) * pkin(4) + t541) * MDP(28) + t539; (t486 + t544) * MDP(24) + (-t474 - t543) * MDP(25) + (t276 * t377 + t540) * MDP(27) + (-t435 * t273 - t439 * t274 + t275 * t377 + t541) * MDP(28) + (-MDP(24) * t518 + MDP(25) * t461 - MDP(27) * t276 - MDP(28) * t528) * qJD(5) + t539;];
tau = t1;
