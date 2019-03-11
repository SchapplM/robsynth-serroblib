% Calculate minimal parameter regressor of Coriolis joint torque vector for
% S6RPRRRP9
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d3,d4,d5]';
% MDP [29x1]
%   Minimal dynamic parameter vector (fixed base model)
%   see S6RPRRRP9_convert_par2_MPV_fixb.m
% 
% Output:
% tauc [6x1]
%   joint torques required to compensate Coriolis and centrifugal load

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 06:29
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauc = S6RPRRRP9_coriolisvecJ_fixb_mdp_slag_vp(qJ, qJD, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(9,1),zeros(29,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRRRP9_coriolisvecJ_fixb_mdp_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RPRRRP9_coriolisvecJ_fixb_mdp_slag_vp: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S6RPRRRP9_coriolisvecJ_fixb_mdp_slag_vp: pkin has to be [9x1] (double)');
assert(isreal(MDP) && all(size(MDP) == [29 1]), ...
  'S6RPRRRP9_coriolisvecJ_fixb_mdp_slag_vp: MDP has to be [29x1] (double)'); 

%% Symbolic Calculation
% From coriolisvec_joint_fixb_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 06:28:40
% EndTime: 2019-03-09 06:28:49
% DurationCPUTime: 5.10s
% Computational Cost: add. (4148->428), mult. (8962->581), div. (0->0), fcn. (5809->6), ass. (0->182)
t420 = sin(qJ(3));
t423 = cos(qJ(3));
t439 = pkin(3) * t423 + pkin(8) * t420;
t384 = t439 * qJD(1);
t422 = cos(qJ(4));
t366 = t422 * t384;
t518 = pkin(8) + pkin(9);
t460 = qJD(4) * t518;
t424 = -pkin(1) - pkin(7);
t403 = qJD(1) * t424 + qJD(2);
t419 = sin(qJ(4));
t505 = t419 * t423;
t463 = t403 * t505;
t503 = t420 * t422;
t465 = pkin(9) * t503;
t531 = -t463 + t366 + (pkin(4) * t423 + t465) * qJD(1) + t422 * t460;
t482 = qJD(1) * t420;
t459 = t419 * t482;
t501 = t422 * t423;
t488 = t419 * t384 + t403 * t501;
t528 = pkin(9) * t459 + t419 * t460 + t488;
t478 = qJD(3) * t422;
t481 = qJD(1) * t423;
t376 = -t419 * t481 + t478;
t421 = cos(qJ(5));
t458 = t422 * t481;
t480 = qJD(3) * t419;
t377 = t458 + t480;
t418 = sin(qJ(5));
t511 = t377 * t418;
t326 = -t421 * t376 + t511;
t324 = t326 ^ 2;
t469 = qJD(1) * qJD(3);
t452 = t423 * t469;
t436 = t376 * t418 + t421 * t377;
t519 = t436 ^ 2;
t530 = MDP(25) * t452 + (-t324 + t519) * MDP(22);
t529 = qJ(6) * t326;
t467 = qJD(4) + qJD(5);
t471 = qJD(5) * t421;
t474 = qJD(4) * t422;
t502 = t421 * t422;
t506 = t418 * t419;
t490 = t418 * t459 - t421 * t474 - t422 * t471 + t467 * t506 - t482 * t502;
t380 = t418 * t422 + t419 * t421;
t333 = t467 * t380;
t359 = t380 * qJD(1);
t489 = t420 * t359 + t333;
t479 = qJD(3) * t420;
t454 = t419 * t479;
t473 = qJD(4) * t423;
t455 = t422 * t473;
t527 = t454 - t455;
t466 = 0.2e1 * qJD(1);
t417 = t423 ^ 2;
t484 = t420 ^ 2 - t417;
t526 = MDP(8) * t484;
t525 = qJ(6) * t436;
t510 = t403 * t423;
t368 = -qJD(3) * pkin(3) - t510;
t335 = -pkin(4) * t376 + t368;
t296 = pkin(5) * t326 + qJD(6) + t335;
t524 = t296 * t436;
t523 = t335 * t436;
t379 = -t502 + t506;
t351 = t379 * t420;
t349 = t380 * t420;
t409 = qJD(4) + t482;
t401 = qJD(5) + t409;
t522 = t401 * t436;
t438 = pkin(3) * t420 - pkin(8) * t423;
t388 = qJ(2) + t438;
t373 = t422 * t388;
t504 = t419 * t424;
t451 = pkin(4) - t504;
t323 = -pkin(9) * t501 + t420 * t451 + t373;
t400 = t424 * t503;
t485 = t419 * t388 + t400;
t334 = -pkin(9) * t505 + t485;
t491 = t418 * t323 + t421 * t334;
t521 = t531 * t421;
t396 = t518 * t419;
t397 = t518 * t422;
t487 = -t418 * t396 + t421 * t397;
t472 = qJD(5) * t418;
t520 = t396 * t471 + t397 * t472 + t531 * t418 + t528 * t421;
t453 = t420 * t478;
t456 = t419 * t473;
t430 = -t453 - t456;
t468 = qJD(3) * qJD(4);
t339 = qJD(1) * t430 + t422 * t468;
t486 = qJD(1) * t455 + t419 * t468;
t431 = qJD(1) * t454 - t486;
t446 = t339 * t418 - t421 * t431;
t289 = qJD(5) * t436 + t446;
t426 = qJD(1) ^ 2;
t517 = qJ(2) * t426;
t516 = t339 * t419;
t361 = t388 * qJD(1);
t515 = t361 * t419;
t514 = t368 * t419;
t513 = t368 * t420;
t512 = t376 * t409;
t509 = t409 * t419;
t508 = t409 * t420;
t507 = t409 * t422;
t387 = t420 * t403;
t367 = qJD(3) * pkin(8) + t387;
t320 = t367 * t422 + t515;
t309 = pkin(9) * t376 + t320;
t305 = t418 * t309;
t307 = t421 * t309;
t425 = qJD(3) ^ 2;
t499 = t424 * t425;
t319 = t422 * t361 - t367 * t419;
t308 = -pkin(9) * t377 + t319;
t303 = pkin(4) * t409 + t308;
t276 = t421 * t303 - t305;
t271 = t276 - t525;
t270 = pkin(5) * t401 + t271;
t498 = t270 - t271;
t497 = -t489 * qJ(6) - qJD(6) * t379 - t520;
t496 = -pkin(5) * t481 + t490 * qJ(6) - qJD(5) * t487 - qJD(6) * t380 + t528 * t418 - t521;
t352 = t379 * t423;
t495 = -qJD(3) * t352 - t349 * t467 - t359;
t477 = qJD(3) * t423;
t494 = t379 * qJD(1) + t351 * t467 - t380 * t477;
t493 = t421 * t308 - t305;
t476 = qJD(4) * t419;
t475 = qJD(4) * t420;
t470 = t436 * MDP(21);
t462 = t409 * t505;
t461 = t421 * t339 + t376 * t471 + t418 * t431;
t415 = -pkin(4) * t422 - pkin(3);
t457 = t422 * t477;
t374 = qJD(3) * t439 + qJD(2);
t348 = t374 * qJD(1);
t342 = t422 * t348;
t428 = -qJD(4) * t320 + t342;
t282 = -pkin(9) * t339 + (pkin(4) * qJD(1) - t403 * t419) * t477 + t428;
t434 = t419 * t348 + t361 * t474 - t367 * t476 + t403 * t457;
t286 = pkin(9) * t431 + t434;
t450 = t421 * t282 - t418 * t286;
t355 = t422 * t374;
t293 = t355 + (-t400 + (pkin(9) * t423 - t388) * t419) * qJD(4) + (t423 * t451 + t465) * qJD(3);
t429 = t419 * t374 + t388 * t474 + t424 * t457 - t475 * t504;
t295 = t527 * pkin(9) + t429;
t449 = t421 * t293 - t295 * t418;
t448 = -t308 * t418 - t307;
t447 = t421 * t323 - t334 * t418;
t445 = -t421 * t396 - t397 * t418;
t375 = pkin(4) * t505 - t423 * t424;
t444 = -t376 + t478;
t443 = -t377 + t480;
t442 = qJD(1) + t475;
t441 = t418 * t282 + t421 * t286 + t303 * t471 - t309 * t472;
t398 = t420 * t452;
t440 = -t387 + (t459 + t476) * pkin(4);
t277 = t303 * t418 + t307;
t435 = t420 * t443;
t340 = -t527 * pkin(4) + t424 * t479;
t433 = t418 * t293 + t421 * t295 + t323 * t471 - t334 * t472;
t288 = t377 * t472 - t461;
t316 = -pkin(4) * t431 + t403 * t479;
t427 = -t277 * qJD(5) + t450;
t275 = t289 * pkin(5) + t316;
t414 = pkin(4) * t421 + pkin(5);
t350 = t380 * t423;
t315 = -qJ(6) * t379 + t487;
t314 = -qJ(6) * t380 + t445;
t302 = -t472 * t505 + (t467 * t501 - t454) * t421 + t430 * t418;
t300 = t333 * t423 - t418 * t454 + t421 * t453;
t290 = -qJ(6) * t350 + t491;
t287 = pkin(5) * t420 + qJ(6) * t352 + t447;
t274 = t493 - t525;
t273 = t448 + t529;
t272 = t277 - t529;
t269 = -qJ(6) * t302 - qJD(6) * t350 + t433;
t268 = pkin(5) * t477 + qJ(6) * t300 - qJD(5) * t491 + qJD(6) * t352 + t449;
t267 = -qJ(6) * t289 - qJD(6) * t326 + t441;
t266 = pkin(5) * t452 + qJ(6) * t288 - qJD(6) * t436 + t427;
t1 = [0.2e1 * t469 * t526 + (-t420 * t499 + (qJ(2) * t477 + qJD(2) * t420) * t466) * MDP(12) + (-t423 * t499 + (-qJ(2) * t479 + qJD(2) * t423) * t466) * MDP(13) + (t339 * t501 + t377 * t430) * MDP(14) + ((-t376 * t422 + t377 * t419) * t479 + (t422 * t431 - t516 + (-t419 * t376 - t377 * t422) * qJD(4)) * t423) * MDP(15) + (-t409 * t456 + t339 * t420 + (t377 * t423 + (qJD(1) * t417 - t508) * t422) * qJD(3)) * MDP(16) + (-t409 * t455 - t486 * t420 + (t376 * t423 + (qJD(1) * t484 + t508) * t419) * qJD(3)) * MDP(17) + (t409 * t477 + t398) * MDP(18) + ((-t388 * t476 + t355) * t409 + (-t424 * t486 + t368 * t474 + (t373 * qJD(1) - t409 * t504 + t319) * qJD(3)) * t423 + (t342 + (-t376 * t424 - t514) * qJD(3) + (-t515 + (-t409 * t424 - t367) * t422) * qJD(4)) * t420) * MDP(19) + (-t429 * t409 - t434 * t420 + (-t424 * t339 - t368 * t476) * t423 + ((-qJD(1) * t485 - t320) * t423 + (t424 * t377 + (-t368 + t510) * t422) * t420) * qJD(3)) * MDP(20) + (t288 * t352 - t300 * t436) * MDP(21) + (t288 * t350 + t289 * t352 + t300 * t326 - t302 * t436) * MDP(22) + (-t288 * t420 - t300 * t401 + (-qJD(1) * t352 + t436) * t477) * MDP(23) + (-t289 * t420 - t302 * t401 + (-qJD(1) * t350 - t326) * t477) * MDP(24) + (t401 * t477 + t398) * MDP(25) + (t449 * t401 + t450 * t420 + t340 * t326 + t375 * t289 + t316 * t350 + t335 * t302 + (-t277 * t420 - t401 * t491) * qJD(5) + (qJD(1) * t447 + t276) * t477) * MDP(26) + (-t433 * t401 - t441 * t420 + t340 * t436 - t375 * t288 - t316 * t352 - t335 * t300 + (-qJD(1) * t491 - t277) * t477) * MDP(27) + (t266 * t352 - t267 * t350 - t268 * t436 - t269 * t326 + t270 * t300 - t272 * t302 + t287 * t288 - t289 * t290) * MDP(28) + (t267 * t290 + t272 * t269 + t266 * t287 + t270 * t268 + t275 * (pkin(5) * t350 + t375) + t296 * (pkin(5) * t302 + t340)) * MDP(29) - 0.2e1 * MDP(7) * t398 + (MDP(6) * qJ(2) + MDP(5)) * qJD(2) * t466 + (-MDP(10) * t423 - MDP(9) * t420) * t425; -t426 * MDP(5) - MDP(6) * t517 + (-t423 * t486 - t442 * t507 + (-t420 * t376 - t462) * qJD(3)) * MDP(19) + (-t339 * t423 + t442 * t509 + (-t409 * t501 + (t377 - t458) * t420) * qJD(3)) * MDP(20) + (-t289 * t423 + t494 * t401 + (t326 * t420 - t349 * t481) * qJD(3)) * MDP(26) + (t288 * t423 - t495 * t401 + (t351 * t481 + t420 * t436) * qJD(3)) * MDP(27) + (-t288 * t349 + t289 * t351 - t326 * t495 - t436 * t494) * MDP(28) + (-t266 * t349 - t267 * t351 + t270 * t494 + t272 * t495 - t275 * t423 + t296 * t479) * MDP(29) + (t420 * MDP(12) + t423 * MDP(13)) * (-t425 - t426); t420 * MDP(13) * t517 + (t377 * t507 + t516) * MDP(14) + ((t339 + t512) * t422 + (qJD(1) * t435 - t377 * qJD(4) - t486) * t419) * MDP(15) + (t409 * t474 + (t409 * t503 + t423 * t443) * qJD(1)) * MDP(16) + (-t409 * t476 + (-t508 * t419 + t423 * t444) * qJD(1)) * MDP(17) + (-pkin(3) * t486 - t366 * t409 + (-t420 * t444 + t462) * t403 + (-pkin(8) * t507 + t514) * qJD(4) + (-t319 * t423 + (qJD(3) * t438 + t513) * t419) * qJD(1)) * MDP(19) + (-pkin(3) * t339 + t488 * t409 + t403 * t435 + (pkin(8) * t509 + t368 * t422) * qJD(4) + (t320 * t423 + (-pkin(8) * t477 + t513) * t422) * qJD(1)) * MDP(20) + (-t288 * t380 - t436 * t490) * MDP(21) + (t288 * t379 - t289 * t380 + t326 * t490 - t436 * t489) * MDP(22) + (t415 * t289 + t316 * t379 + t440 * t326 + t489 * t335) * MDP(26) + (-t415 * t288 + t316 * t380 - t490 * t335 + t436 * t440) * MDP(27) + (-t266 * t380 - t267 * t379 + t270 * t490 - t272 * t489 + t288 * t314 - t289 * t315 - t326 * t497 - t436 * t496) * MDP(28) + (t267 * t315 + t266 * t314 + t275 * (pkin(5) * t379 + t415) + (pkin(4) * t509 + pkin(5) * t489 - t387) * t296 + t497 * t272 + t496 * t270) * MDP(29) + (-t490 * MDP(23) - t489 * MDP(24) + (-t397 * t471 + (qJD(5) * t396 + t528) * t418 - t521) * MDP(26) + t520 * MDP(27)) * t401 + (-t409 * MDP(18) + (qJD(3) * t380 - t436) * MDP(23) + (-qJD(3) * t379 + t326) * MDP(24) - t401 * MDP(25) + (qJD(3) * t445 - t276) * MDP(26) + (-qJD(3) * t487 + t277) * MDP(27)) * t481 + (-t526 + (-qJ(2) * MDP(12) + t420 * MDP(7)) * t423) * t426; -t377 * t376 * MDP(14) + (-t376 ^ 2 + t377 ^ 2) * MDP(15) + (t339 - t512) * MDP(16) + (t377 * t409 + t431) * MDP(17) + MDP(18) * t452 + (-qJD(3) * t463 + t320 * t409 - t368 * t377 + t428) * MDP(19) + (t319 * t409 - t368 * t376 - t434) * MDP(20) + t326 * t470 + (t401 * t326 - t288) * MDP(23) + (-t289 + t522) * MDP(24) + (-t448 * t401 - t523 + (-t326 * t377 - t401 * t472 + t421 * t452) * pkin(4) + t427) * MDP(26) + (t493 * t401 + t335 * t326 + (-t377 * t436 - t401 * t471 - t418 * t452) * pkin(4) - t441) * MDP(27) + (-t270 * t326 + t272 * t436 + t273 * t436 + t274 * t326 + t288 * t414 + (-t289 * t418 + (-t326 * t421 + t418 * t436) * qJD(5)) * pkin(4)) * MDP(28) + (-pkin(5) * t524 + t266 * t414 - t270 * t273 - t272 * t274 + (t267 * t418 - t296 * t377 + (-t270 * t418 + t272 * t421) * qJD(5)) * pkin(4)) * MDP(29) + t530; t461 * MDP(23) + (-t446 + t522) * MDP(24) + (t277 * t401 + t450 - t523) * MDP(26) + (t276 * t401 - t441) * MDP(27) + t498 * MDP(29) * t272 + (t288 * MDP(28) + (t266 - t524) * MDP(29)) * pkin(5) + (t401 * MDP(23) + t335 * MDP(27) - MDP(28) * t498 + t470) * t326 + (-MDP(23) * t511 - t436 * MDP(24) - MDP(26) * t277) * qJD(5) + t530; (-t324 - t519) * MDP(28) + (t270 * t436 + t272 * t326 + t275) * MDP(29);];
tauc  = t1;
