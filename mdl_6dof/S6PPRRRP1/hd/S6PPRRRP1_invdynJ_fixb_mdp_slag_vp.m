% Calculate vector of inverse dynamics joint torques for
% S6PPRRRP1
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
% pkin [12x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,alpha3,d3,d4,d5,theta1,theta2]';
% MDP [21x1]
%   Minimal dynamic parameter vector (fixed base model)
%   see S6PPRRRP1_convert_par2_MPV_fixb.m
% 
% Output:
% tau [6x1]
%   joint torques of inverse dynamics (contains inertial, gravitational Coriolis and centrifugal forces)

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-08 18:55
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tau = S6PPRRRP1_invdynJ_fixb_mdp_slag_vp(qJ, qJD, qJDD, g, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(3,1),zeros(12,1),zeros(21,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PPRRRP1_invdynJ_fixb_mdp_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6PPRRRP1_invdynJ_fixb_mdp_slag_vp: qJD has to be [6x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [6 1]), ...
  'S6PPRRRP1_invdynJ_fixb_mdp_slag_vp: qJDD has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6PPRRRP1_invdynJ_fixb_mdp_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [12 1]), ...
  'S6PPRRRP1_invdynJ_fixb_mdp_slag_vp: pkin has to be [12x1] (double)');
assert(isreal(MDP) && all(size(MDP) == [21 1]), ...
  'S6PPRRRP1_invdynJ_fixb_mdp_slag_vp: MDP has to be [21x1] (double)'); 

%% Symbolic Calculation
% From invdyn_joint_fixb_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-08 18:54:39
% EndTime: 2019-03-08 18:54:43
% DurationCPUTime: 4.62s
% Computational Cost: add. (3512->426), mult. (8867->599), div. (0->0), fcn. (8067->14), ass. (0->197)
t525 = cos(pkin(6));
t389 = t525 * qJD(1) + qJD(2);
t404 = sin(pkin(6));
t402 = sin(pkin(12));
t409 = sin(qJ(3));
t412 = cos(qJ(3));
t405 = cos(pkin(12));
t524 = cos(pkin(7));
t462 = t405 * t524;
t437 = t402 * t412 + t409 * t462;
t431 = t437 * t404;
t403 = sin(pkin(7));
t510 = t403 * t409;
t329 = qJD(1) * t431 + t389 * t510;
t408 = sin(qJ(4));
t411 = cos(qJ(4));
t454 = pkin(4) * t408 - pkin(10) * t411;
t377 = t454 * qJD(4);
t547 = t329 - t377;
t386 = t525 * qJDD(1) + qJDD(2);
t458 = t404 * t462;
t450 = qJD(1) * t458;
t493 = qJD(3) * t409;
t472 = t403 * t493;
t495 = qJD(1) * t404;
t473 = t402 * t495;
t511 = t402 * t409;
t477 = t404 * t511;
t491 = qJD(3) * t412;
t418 = -(qJDD(1) * t458 + t386 * t403) * t412 + qJDD(1) * t477 + t389 * t472 + t450 * t493 + t473 * t491;
t523 = cos(pkin(11));
t453 = t525 * t523;
t522 = sin(pkin(11));
t353 = t402 * t453 + t522 * t405;
t426 = t522 * t402 - t405 * t453;
t464 = t404 * t523;
t538 = t403 * t464 + t426 * t524;
t315 = t353 * t409 + t412 * t538;
t452 = t525 * t522;
t354 = -t402 * t452 + t523 * t405;
t427 = t523 * t402 + t405 * t452;
t463 = t404 * t522;
t537 = -t403 * t463 + t427 * t524;
t317 = t354 * t409 + t412 * t537;
t456 = t412 * t462;
t465 = t403 * t525;
t457 = t412 * t465;
t337 = -t404 * t456 - t457 + t477;
t439 = g(1) * t317 + g(2) * t315 + g(3) * t337;
t546 = t329 * qJD(3) - t418 + t439;
t480 = qJD(3) * qJD(4);
t469 = t411 * t480;
t479 = qJDD(3) * t408;
t545 = -t469 - t479;
t484 = qJD(5) * t408;
t544 = qJD(3) * t484 - qJDD(4);
t407 = sin(qJ(5));
t410 = cos(qJ(5));
t492 = qJD(3) * t411;
t331 = t407 * ((qJD(5) + t492) * qJD(4) + t479) + t544 * t410;
t393 = pkin(5) * t410 + pkin(4);
t543 = -t393 * t411 - pkin(3);
t530 = pkin(5) * t407;
t474 = pkin(9) + t530;
t328 = -t409 * t473 + (t389 * t403 + t450) * t412;
t488 = qJD(4) * t408;
t508 = t407 * t411;
t529 = pkin(9) * t407;
t542 = -t328 * t508 + t410 * t547 - t488 * t529;
t380 = -pkin(4) * t411 - pkin(10) * t408 - pkin(3);
t483 = qJD(5) * t410;
t506 = t410 * t411;
t541 = -t328 * t506 + t380 * t483 - t407 * t547;
t327 = qJD(3) * pkin(9) + t329;
t476 = t404 * t405 * t403;
t349 = -qJD(1) * t476 + t524 * t389;
t540 = -t408 * t327 + t349 * t411;
t316 = t353 * t412 - t409 * t538;
t415 = t426 * t403 - t524 * t464;
t294 = t316 * t411 + t415 * t408;
t318 = t354 * t412 - t409 * t537;
t416 = t427 * t403 + t524 * t463;
t296 = t318 * t411 + t416 * t408;
t338 = t409 * t465 + t431;
t432 = t525 * t524 - t476;
t321 = t338 * t411 + t432 * t408;
t297 = -t321 * t407 + t337 * t410;
t539 = -g(1) * (-t296 * t407 + t317 * t410) - g(2) * (-t294 * t407 + t315 * t410) - g(3) * t297;
t536 = t437 * qJDD(1);
t535 = pkin(5) * t331 + qJDD(6);
t313 = t411 * t327 + t408 * t349;
t310 = qJD(4) * pkin(10) + t313;
t390 = -qJD(5) + t492;
t533 = (pkin(9) * t390 + t310) * qJD(5) + t439;
t489 = qJD(4) * t407;
t494 = qJD(3) * t408;
t370 = t410 * t494 + t489;
t532 = t370 ^ 2;
t527 = qJ(6) + pkin(10);
t526 = qJD(3) * pkin(3);
t521 = qJ(6) * t408;
t519 = qJDD(4) * pkin(4);
t481 = t410 * qJD(4);
t330 = -qJD(5) * t481 + t407 * t544 + t410 * t545;
t518 = t330 * t407;
t348 = -qJDD(1) * t476 + t524 * t386;
t517 = t348 * t408;
t368 = t407 * t494 - t481;
t515 = t368 * t390;
t514 = t370 * t390;
t513 = t370 * t410;
t509 = t403 * t412;
t507 = t408 * t410;
t322 = t380 * qJD(3) - t328;
t284 = -t310 * t407 + t410 * t322;
t282 = -qJ(6) * t370 + t284;
t281 = -pkin(5) * t390 + t282;
t505 = -t282 + t281;
t391 = pkin(9) * t506;
t448 = pkin(5) * t408 - qJ(6) * t506;
t482 = qJD(6) * t410;
t504 = -t408 * t482 + t448 * qJD(4) + (-t391 + (-t380 + t521) * t407) * qJD(5) - t542;
t376 = t454 * qJD(3);
t503 = t407 * t376 + t410 * t540;
t502 = (-pkin(9) * qJD(4) - qJ(6) * qJD(5)) * t507 + (-qJD(6) * t408 + (-pkin(9) * qJD(5) - qJ(6) * qJD(4)) * t411) * t407 + t541;
t466 = qJD(5) * t527;
t501 = t482 - t503 + (qJ(6) * t492 - t466) * t407;
t359 = t410 * t376;
t500 = -t448 * qJD(3) - t410 * t466 - t359 + (-qJD(6) + t540) * t407;
t497 = t407 * t380 + t391;
t400 = t408 ^ 2;
t496 = -t411 ^ 2 + t400;
t490 = qJD(4) * t368;
t487 = qJD(4) * t411;
t486 = qJD(5) * t390;
t485 = qJD(5) * t407;
t478 = t411 * qJDD(3);
t436 = t456 - t511;
t303 = qJDD(3) * pkin(9) + (t386 * t409 + t389 * t491) * t403 + (t436 * qJD(3) * qJD(1) + t536) * t404;
t279 = qJDD(4) * pkin(10) + qJD(4) * t540 + t303 * t411 + t517;
t288 = qJD(3) * t377 + t380 * qJDD(3) + t418;
t475 = t410 * t279 + t407 * t288 + t322 * t483;
t471 = t403 * t491;
t470 = t390 * t481;
t285 = t310 * t410 + t322 * t407;
t298 = t321 * t410 + t337 * t407;
t449 = t408 * t303 + t327 * t487 - t348 * t411 + t349 * t488;
t309 = -qJD(4) * pkin(4) - t540;
t356 = t408 * t524 + t411 * t510;
t341 = -t356 * t407 - t410 * t509;
t446 = -t356 * t410 + t407 * t509;
t443 = -t408 * t480 + t478;
t365 = qJDD(5) - t443;
t445 = t365 * t407 - t390 * t483;
t444 = t365 * t410 + t390 * t485;
t442 = t310 * t485 - t475;
t293 = t316 * t408 - t415 * t411;
t295 = t318 * t408 - t416 * t411;
t320 = t338 * t408 - t432 * t411;
t441 = g(1) * t295 + g(2) * t293 + g(3) * t320;
t440 = g(1) * t296 + g(2) * t294 + g(3) * t321;
t438 = g(1) * t318 + g(2) * t316 + g(3) * t338;
t355 = t408 * t510 - t411 * t524;
t280 = t449 - t519;
t435 = -pkin(10) * t365 - t390 * t309;
t326 = -t328 - t526;
t430 = -pkin(9) * qJDD(4) + (t326 + t328 - t526) * qJD(4);
t429 = -g(1) * t463 + g(2) * t464 - g(3) * t525;
t287 = t410 * t288;
t425 = -t285 * qJD(5) - t407 * t279 + t287;
t424 = pkin(10) * t486 - t280 + t441;
t422 = t441 - t449;
t413 = qJD(4) ^ 2;
t417 = 0.2e1 * qJDD(3) * pkin(3) - pkin(9) * t413 + t546;
t414 = qJD(3) ^ 2;
t382 = t527 * t410;
t381 = t527 * t407;
t367 = t410 * t380;
t364 = t368 ^ 2;
t343 = -t407 * t521 + t497;
t340 = t356 * qJD(4) + t408 * t471;
t339 = -t355 * qJD(4) + t411 * t471;
t336 = -qJ(6) * t507 + t367 + (-pkin(5) - t529) * t411;
t333 = t338 * qJD(3);
t332 = (t436 * t404 + t457) * qJD(3);
t308 = t446 * qJD(5) - t339 * t407 + t410 * t472;
t307 = t341 * qJD(5) + t339 * t410 + t407 * t472;
t299 = pkin(5) * t368 + qJD(6) + t309;
t292 = -t320 * qJD(4) + t332 * t411;
t291 = t321 * qJD(4) + t332 * t408;
t283 = -qJ(6) * t368 + t285;
t277 = t280 + t535;
t276 = t297 * qJD(5) + t292 * t410 + t333 * t407;
t275 = -t298 * qJD(5) - t292 * t407 + t333 * t410;
t274 = -qJ(6) * t331 - qJD(6) * t368 - t442;
t273 = pkin(5) * t365 + qJ(6) * t330 - qJD(6) * t370 + t425;
t1 = [(qJDD(1) - g(3)) * MDP(1) + (t386 * t525 - g(3) + (t402 ^ 2 + t405 ^ 2) * t404 ^ 2 * qJDD(1)) * MDP(2) + (-qJD(3) * t333 - qJDD(3) * t337) * MDP(4) + (-qJD(3) * t332 - qJDD(3) * t338) * MDP(5) + (-t337 * t478 - qJD(4) * t291 - qJDD(4) * t320 + (-t333 * t411 + t337 * t488) * qJD(3)) * MDP(11) + (t337 * t479 - qJD(4) * t292 - qJDD(4) * t321 + (t333 * t408 + t337 * t487) * qJD(3)) * MDP(12) + (-t275 * t390 + t291 * t368 + t297 * t365 + t320 * t331) * MDP(18) + (t276 * t390 + t291 * t370 - t298 * t365 - t320 * t330) * MDP(19) + (-t275 * t370 - t276 * t368 + t297 * t330 - t298 * t331) * MDP(20) + (t273 * t297 + t274 * t298 + t275 * t281 + t276 * t283 + t277 * t320 + t291 * t299 - g(3)) * MDP(21); (t429 + t386) * MDP(2) + (-qJD(4) * t340 - qJDD(4) * t355) * MDP(11) + (-qJD(4) * t339 - qJDD(4) * t356) * MDP(12) + (-t308 * t390 + t331 * t355 + t340 * t368 + t341 * t365) * MDP(18) + (t307 * t390 - t330 * t355 + t340 * t370 + t365 * t446) * MDP(19) + (-t307 * t368 - t308 * t370 + t330 * t341 + t331 * t446) * MDP(20) + (t273 * t341 - t274 * t446 + t277 * t355 + t281 * t308 + t283 * t307 + t299 * t340 + t429) * MDP(21) + ((-qJDD(3) * MDP(5) + (-MDP(11) * t411 + MDP(12) * t408 - MDP(4)) * t414) * t409 + (t443 * MDP(11) + MDP(12) * t545 + qJDD(3) * MDP(4) - t414 * MDP(5)) * t412) * t403; qJDD(3) * MDP(3) + t546 * MDP(4) + (-t386 * t510 - t404 * t536 + (-t389 * t509 - t436 * t495 + t328) * qJD(3) + t438) * MDP(5) + (qJDD(3) * t400 + 0.2e1 * t408 * t469) * MDP(6) + 0.2e1 * (t408 * t478 - t496 * t480) * MDP(7) + (qJDD(4) * t408 + t411 * t413) * MDP(8) + (qJDD(4) * t411 - t408 * t413) * MDP(9) + (t430 * t408 + t417 * t411) * MDP(11) + (-t417 * t408 + t430 * t411) * MDP(12) + (-t330 * t507 + (-t407 * t484 + t411 * t481) * t370) * MDP(13) + ((-t368 * t410 - t370 * t407) * t487 + (t518 - t331 * t410 + (t368 * t407 - t513) * qJD(5)) * t408) * MDP(14) + ((t330 - t470) * t411 + (qJD(4) * t370 + t444) * t408) * MDP(15) + ((t390 * t489 + t331) * t411 + (-t445 - t490) * t408) * MDP(16) + (-t365 * t411 - t390 * t488) * MDP(17) + (t367 * t365 + t542 * t390 + (t380 * t486 - t438) * t407 + (pkin(9) * t490 - t287 + (-pkin(9) * t365 + qJD(4) * t309 + qJD(5) * t322 + t279) * t407 + t533 * t410) * t411 + (pkin(9) * t331 + qJD(4) * t284 + t280 * t407 + t309 * t483 - t328 * t368) * t408) * MDP(18) + (-t497 * t365 + t541 * t390 - t438 * t410 + ((pkin(9) * t370 + t309 * t410) * qJD(4) - t533 * t407 + t475) * t411 + (-t309 * t485 - t285 * qJD(4) + t280 * t410 - t328 * t370 + (-t330 - t470) * pkin(9)) * t408) * MDP(19) + (t330 * t336 - t331 * t343 - t504 * t370 - t502 * t368 + (-t281 * t410 - t283 * t407) * t487 + (-t273 * t410 - t274 * t407 + (t281 * t407 - t283 * t410) * qJD(5) + t439) * t408) * MDP(20) + (t274 * t343 + t273 * t336 - g(1) * (t317 * t543 + t474 * t318) - g(2) * (t315 * t543 + t474 * t316) - g(3) * (t337 * t543 + t474 * t338) + t502 * t283 + t504 * t281 + t299 * t474 * t487 + (t277 * t474 + (pkin(5) * t483 - t328) * t299 + t439 * t527) * t408) * MDP(21); MDP(8) * t479 + MDP(9) * t478 + qJDD(4) * MDP(10) + (qJD(4) * t313 - t326 * t494 + t422) * MDP(11) + (-t517 + (-qJD(3) * t326 - t303) * t411 + t440) * MDP(12) + (-t390 * t513 - t518) * MDP(13) + ((-t330 + t515) * t410 + (-t331 + t514) * t407) * MDP(14) + ((-t370 * t408 + t390 * t506) * qJD(3) + t445) * MDP(15) + ((t368 * t408 - t390 * t508) * qJD(3) + t444) * MDP(16) + t390 * MDP(17) * t494 + (-t284 * t494 - pkin(4) * t331 - t313 * t368 + t359 * t390 + (-t390 * t540 + t435) * t407 + t424 * t410) * MDP(18) + (pkin(4) * t330 + t285 * t494 - t313 * t370 - t503 * t390 - t424 * t407 + t435 * t410) * MDP(19) + (-t330 * t381 - t331 * t382 - t500 * t370 - t501 * t368 + (t390 * t281 + t274) * t410 + (t390 * t283 - t273) * t407 - t440) * MDP(20) + (t274 * t382 - t273 * t381 - t277 * t393 - g(1) * (-t295 * t393 + t296 * t527) - g(2) * (-t293 * t393 + t294 * t527) - g(3) * (-t320 * t393 + t321 * t527) + (-t390 * t530 - t313) * t299 + t501 * t283 + t500 * t281) * MDP(21) + (-t408 * t411 * MDP(6) + t496 * MDP(7)) * t414; t370 * t368 * MDP(13) + (-t364 + t532) * MDP(14) + (-t330 - t515) * MDP(15) + (-t331 - t514) * MDP(16) + t365 * MDP(17) + (-t285 * t390 - t309 * t370 + t425 + t539) * MDP(18) + (-t284 * t390 + t309 * t368 - g(1) * (-t296 * t410 - t317 * t407) - g(2) * (-t294 * t410 - t315 * t407) + g(3) * t298 + t442) * MDP(19) + (pkin(5) * t330 - t505 * t368) * MDP(20) + (t505 * t283 + (-t299 * t370 + t273 + t539) * pkin(5)) * MDP(21); (-t364 - t532) * MDP(20) + (t281 * t370 + t283 * t368 - t422 - t519 + t535) * MDP(21);];
tau  = t1;
