% Calculate minimal parameter regressor of Coriolis joint torque vector for
% S6RRRPPR2
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d3,d6,theta4]';
% MDP [30x1]
%   Minimal dynamic parameter vector (fixed base model)
%   see S6RRRPPR2_convert_par2_MPV_fixb.m
% 
% Output:
% tauc [6x1]
%   joint torques required to compensate Coriolis and centrifugal load

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 15:27
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauc = S6RRRPPR2_coriolisvecJ_fixb_mdp_slag_vp(qJ, qJD, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(10,1),zeros(30,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRPPR2_coriolisvecJ_fixb_mdp_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRRPPR2_coriolisvecJ_fixb_mdp_slag_vp: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RRRPPR2_coriolisvecJ_fixb_mdp_slag_vp: pkin has to be [10x1] (double)');
assert(isreal(MDP) && all(size(MDP) == [30 1]), ...
  'S6RRRPPR2_coriolisvecJ_fixb_mdp_slag_vp: MDP has to be [30x1] (double)'); 

%% Symbolic Calculation
% From coriolisvec_joint_fixb_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 15:26:50
% EndTime: 2019-03-09 15:27:00
% DurationCPUTime: 4.78s
% Computational Cost: add. (5956->415), mult. (15470->529), div. (0->0), fcn. (11421->8), ass. (0->198)
t438 = cos(qJ(3));
t439 = cos(qJ(2));
t495 = qJD(1) * t439;
t483 = t438 * t495;
t435 = sin(qJ(3));
t436 = sin(qJ(2));
t496 = qJD(1) * t436;
t484 = t435 * t496;
t389 = -t483 + t484;
t391 = -t435 * t495 - t438 * t496;
t432 = sin(pkin(10));
t433 = cos(pkin(10));
t363 = t389 * t433 - t391 * t432;
t523 = pkin(5) * t363;
t437 = cos(qJ(6));
t428 = qJD(2) + qJD(3);
t434 = sin(qJ(6));
t510 = t428 * t434;
t349 = -t437 * t363 + t510;
t460 = -t389 * t432 - t433 * t391;
t528 = qJD(6) + t460;
t538 = t349 * t528;
t351 = t363 * t434 + t428 * t437;
t537 = t351 * t528;
t536 = t363 * t428;
t535 = t528 * t437;
t526 = pkin(7) + pkin(8);
t410 = t526 * t439;
t406 = qJD(1) * t410;
t396 = t438 * t406;
t409 = t526 * t436;
t404 = qJD(1) * t409;
t472 = t404 * t435 - t396;
t521 = qJ(4) * t389;
t352 = t472 + t521;
t386 = t391 * qJ(4);
t392 = t435 * t406;
t498 = -t438 * t404 - t392;
t353 = t386 + t498;
t417 = t432 * t435 * pkin(2);
t493 = qJD(3) * t438;
t499 = -qJD(3) * t417 - t352 * t432 + (pkin(2) * t493 - t353) * t433;
t534 = t460 ^ 2;
t486 = qJD(1) * qJD(2);
t533 = -0.2e1 * t486;
t524 = pkin(5) * t460;
t532 = MDP(4) * t436;
t531 = MDP(5) * (t436 ^ 2 - t439 ^ 2);
t522 = qJD(2) * pkin(2);
t398 = -t404 + t522;
t485 = qJD(2) * t526;
t465 = qJD(1) * t485;
t399 = t436 * t465;
t530 = t438 * (qJD(3) * t398 - t399);
t470 = t528 * t434;
t503 = -qJD(5) - t499;
t508 = t433 * t435;
t500 = -t433 * t352 + t353 * t432 - (t432 * t438 + t508) * qJD(3) * pkin(2);
t475 = t438 * t398 - t392;
t346 = t386 + t475;
t403 = t435 * t439 + t436 * t438;
t529 = qJD(1) * t403;
t459 = -t398 * t435 - t396;
t347 = -t459 - t521;
t342 = t432 * t347;
t319 = t346 * t433 - t342;
t489 = qJD(5) - t319;
t373 = t428 * t403;
t368 = t373 * qJD(1);
t400 = t439 * t465;
t494 = qJD(3) * t435;
t473 = -t435 * t400 - t406 * t494;
t309 = -qJ(4) * t368 - qJD(4) * t389 + t473 + t530;
t482 = t439 * t486;
t367 = qJD(3) * t483 - t428 * t484 + t438 * t482;
t474 = t435 * t399 - t438 * t400;
t446 = qJD(3) * t459 + t474;
t310 = -qJ(4) * t367 + qJD(4) * t391 + t446;
t287 = t309 * t432 - t433 * t310;
t424 = -pkin(2) * t439 - pkin(1);
t408 = t424 * qJD(1);
t374 = pkin(3) * t389 + qJD(4) + t408;
t444 = -qJ(5) * t460 + t374;
t321 = pkin(4) * t363 + t444;
t454 = t321 * t460 + t287;
t527 = pkin(4) + pkin(9);
t525 = pkin(3) * t391;
t359 = -qJ(4) * t403 - t409 * t438 - t410 * t435;
t402 = t435 * t436 - t438 * t439;
t458 = t409 * t435 - t410 * t438;
t360 = -qJ(4) * t402 - t458;
t329 = -t433 * t359 + t360 * t432;
t520 = t287 * t329;
t492 = qJD(6) * t434;
t334 = t367 * t432 + t433 * t368;
t491 = qJD(6) * t437;
t502 = t434 * t334 + t363 * t491;
t304 = -t428 * t492 + t502;
t519 = t304 * t437;
t370 = t433 * t402 + t403 * t432;
t371 = -t402 * t432 + t403 * t433;
t457 = pkin(3) * t402 + t424;
t449 = -qJ(5) * t371 + t457;
t311 = t370 * t527 + t449;
t335 = t367 * t433 - t368 * t432;
t518 = t311 * t335;
t509 = t433 * t347;
t318 = t346 * t432 + t509;
t517 = t318 * t460;
t516 = t335 * t434;
t515 = t349 * t363;
t514 = t351 * t363;
t512 = t370 * t434;
t511 = t408 * t391;
t440 = qJD(2) ^ 2;
t507 = t436 * t440;
t333 = t437 * t335;
t506 = t439 * t440;
t441 = qJD(1) ^ 2;
t505 = t439 * t441;
t504 = -t523 + t500;
t288 = t433 * t309 + t432 * t310;
t341 = pkin(3) * t428 + t346;
t317 = t432 * t341 + t509;
t501 = t524 - t503;
t490 = t524 + t489;
t426 = t436 * t522;
t425 = pkin(2) * t496;
t286 = -t428 * qJD(5) - t288;
t421 = -pkin(3) * t433 - pkin(4);
t481 = -pkin(2) * t428 - t398;
t480 = pkin(3) * t368 + qJD(2) * t425;
t479 = pkin(3) * t373 + t426;
t478 = t500 * t460;
t477 = pkin(1) * t533;
t274 = -pkin(5) * t334 - t286;
t316 = t341 * t433 - t342;
t464 = qJD(5) - t316;
t294 = -t428 * t527 + t464 + t524;
t299 = t363 * t527 + t444;
t279 = t294 * t434 + t299 * t437;
t476 = t274 * t437 - t279 * t363;
t405 = t436 * t485;
t407 = t439 * t485;
t451 = -t438 * t405 - t435 * t407 - t409 * t493 - t410 * t494;
t322 = -qJ(4) * t373 - qJD(4) * t402 + t451;
t372 = t428 * t402;
t445 = qJD(3) * t458 + t405 * t435 - t438 * t407;
t323 = qJ(4) * t372 - qJD(4) * t403 + t445;
t290 = t322 * t432 - t433 * t323;
t423 = pkin(2) * t438 + pkin(3);
t384 = t423 * t433 - t417;
t315 = -qJ(5) * t428 - t317;
t296 = -t315 - t523;
t471 = t528 * t296;
t327 = pkin(4) * t460 + qJ(5) * t363 - t525;
t326 = t327 + t425;
t357 = t460 * pkin(9);
t379 = -pkin(4) - t384;
t376 = -pkin(9) + t379;
t467 = -qJD(6) * t376 + t326 + t357;
t418 = -pkin(9) + t421;
t466 = -qJD(6) * t418 + t327 + t357;
t278 = t294 * t437 - t299 * t434;
t314 = -pkin(4) * t428 + t464;
t463 = t314 * t363 - t315 * t460;
t462 = -t316 * t363 + t317 * t460;
t291 = t322 * t433 + t323 * t432;
t330 = t359 * t432 + t360 * t433;
t456 = t274 * t434 + t278 * t363 + (t437 * t460 + t491) * t296;
t455 = -t321 * t363 - t286;
t385 = pkin(2) * t508 + t423 * t432;
t453 = t408 * t389 - t473;
t336 = -t372 * t432 + t433 * t373;
t452 = t336 * t434 + t370 * t491;
t312 = pkin(5) * t371 + t329;
t450 = t274 * t370 + t296 * t336 - t312 * t335;
t448 = -qJ(5) * t335 - qJD(5) * t460 + t480;
t337 = -t372 * t433 - t373 * t432;
t447 = -qJ(5) * t337 - qJD(5) * t371 + t479;
t289 = pkin(4) * t334 + t448;
t443 = t287 * t371 + t290 * t460 - t291 * t363 + t329 * t335 - t330 * t334;
t332 = t437 * t334;
t305 = qJD(6) * t351 - t332;
t442 = -t391 * t389 * MDP(11) + t528 * t363 * MDP(28) + ((-t305 - t537) * t437 + (-t304 + t538) * t434) * MDP(25) + (-t470 * t528 + t333 + t514) * MDP(26) + (-t528 * t535 - t515 - t516) * MDP(27) + (-t351 * t470 + t519) * MDP(24) + t367 * MDP(13) + (-t389 ^ 2 + t391 ^ 2) * MDP(12) + (t389 * MDP(13) + (-t391 - t529) * MDP(14)) * t428;
t420 = pkin(3) * t432 + qJ(5);
t378 = qJ(5) + t385;
t328 = pkin(4) * t370 + t449;
t313 = -pkin(5) * t370 + t330;
t297 = t318 - t523;
t292 = pkin(4) * t336 + t447;
t285 = t336 * t527 + t447;
t283 = -pkin(5) * t336 + t291;
t282 = pkin(5) * t337 + t290;
t277 = t334 * t527 + t448;
t276 = pkin(5) * t335 + t287;
t275 = t437 * t276;
t1 = [(t288 * t330 - t316 * t290 + t317 * t291 + t374 * t479 + t457 * t480 + t520) * MDP(19) + (-t286 * t330 + t289 * t328 + t290 * t314 - t291 * t315 + t292 * t321 + t520) * MDP(23) + t531 * t533 + MDP(6) * t506 + (-pkin(7) * t506 + t436 * t477) * MDP(9) + (t370 * t333 - t305 * t371 - t337 * t349 + (t336 * t437 - t370 * t492) * t528) * MDP(27) + (t424 * t367 - t408 * t372 + (-t391 + t529) * t426) * MDP(17) + (t424 * t368 + t408 * t373 + (qJD(1) * t402 + t389) * t426) * MDP(16) - MDP(7) * t507 + (pkin(7) * t507 + t439 * t477) * MDP(10) + (-t288 * t370 - t316 * t337 - t317 * t336 + t443) * MDP(18) + (t286 * t370 + t314 * t337 + t315 * t336 + t443) * MDP(20) + (t304 * t371 + t335 * t512 + t337 * t351 + t452 * t528) * MDP(26) + (t304 * t512 + t351 * t452) * MDP(24) + (t275 * t371 + t278 * t337 + t283 * t349 + t313 * t305 + (-t277 * t371 - t285 * t528 - t518) * t434 + (t282 * t528 - t450) * t437 + ((-t311 * t437 - t312 * t434) * t528 - t279 * t371 + t296 * t512) * qJD(6)) * MDP(29) + (-t279 * t337 + t283 * t351 + t313 * t304 + (-(qJD(6) * t312 + t285) * t528 - t518 - (qJD(6) * t294 + t277) * t371 + t296 * qJD(6) * t370) * t437 + (-(-qJD(6) * t311 + t282) * t528 - (-qJD(6) * t299 + t276) * t371 + t450) * t434) * MDP(30) + ((-t349 * t434 + t351 * t437) * t336 + (t519 - t305 * t434 + (-t349 * t437 - t351 * t434) * qJD(6)) * t370) * MDP(25) + 0.2e1 * t482 * t532 + (t335 * t371 + t337 * t528) * MDP(28) + (t367 * t403 + t372 * t391) * MDP(11) + (-t367 * t402 - t368 * t403 + t372 * t389 + t373 * t391) * MDP(12) + (-t289 * t370 - t292 * t363 - t321 * t336 - t328 * t334) * MDP(21) + (-t289 * t371 - t292 * t460 - t321 * t337 - t328 * t335) * MDP(22) + (-t372 * MDP(13) - t373 * MDP(14) + MDP(16) * t445 - MDP(17) * t451 + t290 * MDP(21) + t291 * MDP(22)) * t428; -t505 * t532 + (-t286 * t378 + t287 * t379 - t314 * t500 + t315 * t503 - t321 * t326) * MDP(23) + (t326 * t460 - t428 * t503 + t455) * MDP(22) + (t378 * t304 + t467 * t535 + t501 * t351 + (-t376 * t335 + t504 * t528 - t471) * t434 + t476) * MDP(30) + (-t334 * t385 - t335 * t384 - t363 * t499 + t462 - t478) * MDP(18) + (t326 * t363 - t428 * t500 + t454) * MDP(21) + t442 + (t288 * t385 - t287 * t384 - t374 * (t425 - t525) + t499 * t317 + t500 * t316) * MDP(19) + t441 * t531 + (t376 * t333 + t378 * t305 + t501 * t349 + (t434 * t467 - t437 * t504) * t528 + t456) * MDP(29) + (t391 * t425 + t498 * t428 + (qJD(3) * t481 + t399) * t438 + t453) * MDP(17) + (-t389 * t425 + t511 - t472 * t428 + (t435 * t481 - t396) * qJD(3) + t474) * MDP(16) + (-t334 * t378 + t335 * t379 + t363 * t503 + t463 - t478) * MDP(20) + (MDP(9) * t436 * t441 + MDP(10) * t505) * pkin(1); (-t517 + t319 * t363 + (-t334 * t432 - t335 * t433) * pkin(3) + t462) * MDP(18) + (-t286 * t420 + t287 * t421 - t314 * t318 - t315 * t489 - t321 * t327) * MDP(23) + (t418 * t333 + t420 * t305 + (-t437 * t297 + t434 * t466) * t528 + t490 * t349 + t456) * MDP(29) + (t327 * t460 + t428 * t489 + t455) * MDP(22) + (-t334 * t420 + t335 * t421 - t363 * t489 + t463 - t517) * MDP(20) + (t420 * t304 + t466 * t535 + t490 * t351 + (t297 * t528 - t418 * t335 - t471) * t434 + t476) * MDP(30) + t442 + (-t428 * t459 + t446 + t511) * MDP(16) + (t316 * t318 - t317 * t319 + (-t287 * t433 + t288 * t432 + t374 * t391) * pkin(3)) * MDP(19) + (t428 * t475 + t453 - t530) * MDP(17) + (-t318 * t428 + t327 * t363 + t454) * MDP(21); (t316 * t460 + t317 * t363 + t480) * MDP(19) + (-t428 * t460 - t334) * MDP(21) + (-t335 + t536) * MDP(22) + (-t314 * t460 - t315 * t363 + t289) * MDP(23) + (t515 - t516) * MDP(29) + (-t333 + t514) * MDP(30) + (-MDP(29) * t535 + MDP(30) * t470) * t528 + (MDP(18) + MDP(20)) * (-t363 ^ 2 - t534); (t335 + t536) * MDP(20) - t460 * t363 * MDP(21) + (-t428 ^ 2 - t534) * MDP(22) + (t315 * t428 + t454) * MDP(23) + (-t349 * t428 + t333) * MDP(29) + (-t351 * t428 - t516) * MDP(30) - (MDP(29) * t434 + MDP(30) * t437) * t528 ^ 2; t351 * t349 * MDP(24) + (-t349 ^ 2 + t351 ^ 2) * MDP(25) + (t502 + t538) * MDP(26) + (t332 + t537) * MDP(27) + t335 * MDP(28) + (-t277 * t434 + t279 * t528 - t296 * t351 + t275) * MDP(29) + (-t276 * t434 - t277 * t437 + t278 * t528 + t296 * t349) * MDP(30) + (-MDP(26) * t510 - MDP(27) * t351 - MDP(29) * t279 - MDP(30) * t278) * qJD(6);];
tauc  = t1;
