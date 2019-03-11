% Calculate minimal parameter regressor of Coriolis joint torque vector for
% S6PRRRPR3
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [11x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d2,d3,d4,d6,theta1]';
% MDP [29x1]
%   Minimal dynamic parameter vector (fixed base model)
%   see S6PRRRPR3_convert_par2_MPV_fixb.m
% 
% Output:
% tauc [6x1]
%   joint torques required to compensate Coriolis and centrifugal load

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-08 23:14
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauc = S6PRRRPR3_coriolisvecJ_fixb_mdp_slag_vp(qJ, qJD, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(11,1),zeros(29,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRRRPR3_coriolisvecJ_fixb_mdp_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6PRRRPR3_coriolisvecJ_fixb_mdp_slag_vp: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6PRRRPR3_coriolisvecJ_fixb_mdp_slag_vp: pkin has to be [11x1] (double)');
assert(isreal(MDP) && all(size(MDP) == [29 1]), ...
  'S6PRRRPR3_coriolisvecJ_fixb_mdp_slag_vp: MDP has to be [29x1] (double)'); 

%% Symbolic Calculation
% From coriolisvec_joint_fixb_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-08 23:14:09
% EndTime: 2019-03-08 23:14:17
% DurationCPUTime: 3.80s
% Computational Cost: add. (3207->366), mult. (7963->496), div. (0->0), fcn. (5943->10), ass. (0->178)
t403 = cos(qJ(3));
t509 = cos(qJ(4));
t463 = t509 * t403;
t441 = qJD(2) * t463;
t399 = sin(qJ(4));
t400 = sin(qJ(3));
t476 = qJD(2) * t400;
t457 = t399 * t476;
t358 = -t441 + t457;
t402 = cos(qJ(6));
t392 = qJD(3) + qJD(4);
t398 = sin(qJ(6));
t496 = t392 * t398;
t344 = -t402 * t358 + t496;
t367 = t399 * t403 + t509 * t400;
t360 = t367 * qJD(2);
t522 = qJD(6) + t360;
t524 = t344 * t522;
t346 = t358 * t398 + t392 * t402;
t447 = t522 * t346;
t523 = t522 * t402;
t510 = -pkin(9) - pkin(8);
t401 = sin(qJ(2));
t396 = sin(pkin(6));
t479 = qJD(1) * t396;
t462 = t401 * t479;
t505 = qJD(3) * pkin(3);
t418 = t400 * t505 - t462;
t521 = MDP(17) - MDP(20);
t520 = MDP(18) - MDP(21);
t519 = MDP(5) * t403;
t518 = (t400 ^ 2 - t403 ^ 2) * MDP(6);
t453 = -qJD(2) * t510 + t462;
t397 = cos(pkin(6));
t478 = qJD(1) * t397;
t341 = t400 * t478 + t453 * t403;
t336 = t509 * t341;
t340 = -t453 * t400 + t403 * t478;
t297 = t399 * t340 + t336;
t474 = qJD(4) * t399;
t438 = pkin(3) * t474 - t297;
t335 = t399 * t341;
t298 = t509 * t340 - t335;
t456 = qJD(4) * t509;
t486 = -pkin(3) * t456 - qJD(5) + t298;
t464 = qJD(3) * t510;
t369 = t400 * t464;
t370 = t403 * t464;
t373 = t510 * t400;
t374 = t510 * t403;
t404 = cos(qJ(2));
t461 = t404 * t479;
t492 = t399 * t400;
t516 = t463 - t492;
t485 = -t509 * t369 - t399 * t370 - t373 * t456 - t374 * t474 + t516 * t461;
t427 = t399 * t373 - t509 * t374;
t484 = t427 * qJD(4) - t367 * t461 + t399 * t369 - t509 * t370;
t434 = t392 * t492;
t338 = -qJD(3) * t463 - t403 * t456 + t434;
t517 = -qJ(5) * t338 + qJD(5) * t367 - t418;
t495 = t396 * t401;
t428 = -t397 * t403 + t400 * t495;
t515 = t428 * qJD(3);
t337 = t340 + t505;
t295 = -t509 * t337 + t335;
t469 = qJD(5) + t295;
t512 = t360 ^ 2;
t511 = pkin(4) + pkin(10);
t508 = t358 * pkin(5);
t507 = t360 * pkin(5);
t506 = qJD(2) * pkin(2);
t473 = qJD(6) * t398;
t339 = t392 * t367;
t326 = t339 * qJD(2);
t472 = qJD(6) * t402;
t483 = t398 * t326 + t358 * t472;
t293 = -t392 * t473 + t483;
t502 = t293 * t402;
t296 = t399 * t337 + t336;
t501 = t296 * t392;
t482 = t392 * t441;
t325 = qJD(2) * t434 - t482;
t500 = t325 * t398;
t498 = t358 * t360;
t497 = t516 * t398;
t494 = t396 * t404;
t407 = qJD(2) ^ 2;
t493 = t396 * t407;
t406 = qJD(3) ^ 2;
t491 = t400 * t406;
t321 = t402 * t325;
t490 = t403 * t406;
t489 = t511 * t325;
t488 = -pkin(5) * t339 - t485;
t487 = -t338 * pkin(5) + t484;
t467 = qJD(2) * qJD(3);
t454 = t400 * t467;
t477 = qJD(2) * t396;
t455 = qJD(1) * t477;
t355 = pkin(3) * t454 + t401 * t455;
t481 = t507 - t486;
t475 = qJD(2) * t401;
t291 = -qJ(5) * t392 - t296;
t284 = -t291 - t508;
t471 = t284 * qJD(6);
t470 = t507 + t469;
t465 = t401 * t493;
t388 = -pkin(3) * t403 - pkin(2);
t459 = t396 * t475;
t458 = t404 * t477;
t440 = t404 * t455;
t313 = t340 * qJD(3) + t403 * t440;
t314 = -t341 * qJD(3) - t400 * t440;
t444 = t509 * t313 + t399 * t314 + t337 * t456 - t341 * t474;
t267 = -t392 * qJD(5) - t444;
t264 = -pkin(5) * t326 - t267;
t280 = -t511 * t392 + t470;
t352 = t388 * qJD(2) - t461;
t411 = -qJ(5) * t360 + t352;
t292 = t511 * t358 + t411;
t269 = t280 * t398 + t292 * t402;
t452 = t264 * t402 - t269 * t358;
t270 = t399 * t313 - t509 * t314 + t337 * t474 + t341 * t456;
t266 = -pkin(5) * t325 + t270;
t417 = qJ(5) * t325 - qJD(5) * t360 + t355;
t272 = t511 * t326 + t417;
t451 = t402 * t266 - t272 * t398;
t450 = t398 * t522;
t448 = t522 * t284;
t327 = pkin(4) * t360 + qJ(5) * t358;
t315 = pkin(3) * t476 + t327;
t353 = t360 * pkin(10);
t387 = -t509 * pkin(3) - pkin(4);
t384 = -pkin(10) + t387;
t446 = -qJD(6) * t384 + t315 + t353;
t445 = qJD(6) * t511 + t327 + t353;
t443 = t400 * t458;
t442 = t403 * t458;
t439 = t508 + t438;
t347 = -t509 * t373 - t399 * t374;
t437 = -t511 * t339 + t517;
t436 = -pkin(4) * t339 + t517;
t268 = t280 * t402 - t292 * t398;
t432 = -qJ(5) * t367 + t388;
t431 = t264 * t398 + t268 * t358 + (t284 * t360 + t471) * t402;
t357 = t397 * t400 + t403 * t495;
t318 = t399 * t357 + t428 * t509;
t430 = -t318 * t398 + t402 * t494;
t429 = t318 * t402 + t398 * t494;
t319 = t509 * t357 - t399 * t428;
t424 = t339 * t398 - t472 * t516;
t422 = t506 * qJD(2);
t305 = pkin(4) * t358 + t411;
t421 = t305 * t360 + t270;
t420 = -t352 * t360 - t270;
t419 = t352 * t358 - t444;
t416 = -t305 * t358 - t267;
t316 = t367 * pkin(5) + t347;
t415 = -t264 * t516 + t284 * t339 + t316 * t325;
t414 = t357 * qJD(3);
t413 = -0.2e1 * qJD(3) * t506;
t412 = t482 + (t358 - t457) * t392;
t410 = t414 + t443;
t409 = t442 - t515;
t322 = t402 * t326;
t294 = t346 * qJD(6) - t322;
t408 = t522 * t358 * MDP(27) + MDP(12) * t498 + ((-t294 - t447) * t402 + (-t293 + t524) * t398) * MDP(24) + (-t398 * t447 + t502) * MDP(23) + (t346 * t358 - t450 * t522 - t321) * MDP(25) + (-t344 * t358 - t522 * t523 + t500) * MDP(26) + t412 * MDP(14) + (-t358 ^ 2 + t512) * MDP(13);
t385 = pkin(3) * t399 + qJ(5);
t330 = -pkin(4) * t516 + t432;
t317 = pkin(5) * t516 + t427;
t308 = -t511 * t516 + t432;
t307 = t325 * t367;
t290 = -pkin(4) * t392 + t469;
t287 = t296 - t508;
t279 = t319 * qJD(4) + t399 * t409 + t509 * t410;
t278 = t357 * t474 + t399 * t410 - t509 * t409 + t428 * t456;
t276 = pkin(4) * t326 + t417;
t1 = [-MDP(3) * t465 - t404 * MDP(4) * t493 + (-t403 * t465 + (-t414 - 0.2e1 * t443) * qJD(3)) * MDP(10) + (t400 * t465 + (-0.2e1 * t442 + t515) * qJD(3)) * MDP(11) + (t278 * t358 + t279 * t360 - t318 * t325 - t319 * t326) * MDP(19) + (-t267 * t319 + t270 * t318 + t278 * t291 + t279 * t290 + (-t276 * t404 + t305 * t475) * t396) * MDP(22) + ((t430 * qJD(6) + t279 * t402 - t398 * t459) * t522 - t429 * t325 - t278 * t344 + t319 * t294) * MDP(28) + (-(t429 * qJD(6) + t279 * t398 + t402 * t459) * t522 - t430 * t325 - t278 * t346 + t319 * t293) * MDP(29) + t521 * ((-t326 * t404 + t358 * t475) * t396 - t279 * t392) - t520 * ((-t325 * t404 - t360 * t475) * t396 - t278 * t392); 0.2e1 * t454 * t519 - 0.2e1 * t467 * t518 + MDP(7) * t490 - MDP(8) * t491 + (-pkin(8) * t490 + t400 * t413) * MDP(10) + (pkin(8) * t491 + t403 * t413) * MDP(11) + (-t338 * t360 - t307) * MDP(12) + (-t325 * t516 - t326 * t367 + t338 * t358 - t339 * t360) * MDP(13) + (t326 * t388 + t339 * t352 - t355 * t516 + t418 * t358) * MDP(17) + (-t325 * t388 - t338 * t352 + t355 * t367 + t418 * t360) * MDP(18) + (-t267 * t516 + t270 * t367 - t290 * t338 + t291 * t339 - t325 * t347 - t326 * t427 + t485 * t358 + t484 * t360) * MDP(19) + (t276 * t516 - t305 * t339 - t326 * t330 + t436 * t358) * MDP(20) + (-t276 * t367 + t305 * t338 + t325 * t330 + t436 * t360) * MDP(21) + (-t267 * t427 + t270 * t347 + t276 * t330 + t484 * t290 + t485 * t291 - t436 * t305) * MDP(22) + (-t293 * t497 + t424 * t346) * MDP(23) + ((-t344 * t398 + t346 * t402) * t339 - (t502 - t294 * t398 + (-t344 * t402 - t346 * t398) * qJD(6)) * t516) * MDP(24) + (t293 * t367 + t325 * t497 - t338 * t346 + t424 * t522) * MDP(25) + (t516 * t321 - t294 * t367 + t338 * t344 + (t339 * t402 + t473 * t516) * t522) * MDP(26) + (-t338 * t522 - t307) * MDP(27) + (t308 * t500 + t451 * t367 - t268 * t338 + t317 * t294 - t415 * t402 + (t437 * t398 + t487 * t402) * t522 + t488 * t344 + ((-t308 * t402 - t316 * t398) * t522 - t269 * t367 - t284 * t497) * qJD(6)) * MDP(28) + (t269 * t338 + t317 * t293 + t488 * t346 + (t308 * t325 - (qJD(6) * t280 + t272) * t367 - t516 * t471 + (-qJD(6) * t316 + t437) * t522) * t402 + (-(-qJD(6) * t292 + t266) * t367 + (qJD(6) * t308 - t487) * t522 + t415) * t398) * MDP(29) + (-t338 * MDP(14) - t339 * MDP(15) - t484 * t521 + t485 * t520) * t392; (t315 * t360 - t486 * t392 + t416) * MDP(21) + (t385 * t293 + t446 * t523 + t481 * t346 + (t384 * t325 - t439 * t522 - t448) * t398 + t452) * MDP(29) + (t297 * t392 + (-t358 * t476 - t392 * t474) * pkin(3) + t420) * MDP(17) + (t298 * t392 + (-t360 * t476 - t392 * t456) * pkin(3) + t419) * MDP(18) + (-t325 * t387 - t326 * t385 + (-t291 + t438) * t360 + (t290 + t486) * t358) * MDP(19) + (t315 * t358 + t438 * t392 + t421) * MDP(20) + (-t267 * t385 + t270 * t387 + t438 * t290 + t486 * t291 - t305 * t315) * MDP(22) + t407 * t518 + (-t384 * t321 + t385 * t294 + t481 * t344 + (t446 * t398 + t439 * t402) * t522 + t431) * MDP(28) + t403 * t422 * MDP(11) + t408 + (MDP(10) * t422 - t407 * t519) * t400; (t420 + t501) * MDP(17) + (pkin(4) * t325 - qJ(5) * t326 + (-t291 - t296) * t360 + (t290 - t469) * t358) * MDP(19) + (t327 * t358 + t421 - t501) * MDP(20) + (t327 * t360 + t469 * t392 + t416) * MDP(21) + (qJ(5) * t293 + t445 * t523 + t470 * t346 + (t287 * t522 - t448 - t489) * t398 + t452) * MDP(29) + (-pkin(4) * t270 - qJ(5) * t267 - t290 * t296 - t469 * t291 - t305 * t327) * MDP(22) + t408 + (-t295 * t392 + t419) * MDP(18) + (t402 * t489 + qJ(5) * t294 + (-t402 * t287 + t445 * t398) * t522 + t470 * t344 + t431) * MDP(28); t412 * MDP(19) - MDP(20) * t498 + (-t392 ^ 2 - t512) * MDP(21) + (t291 * t392 + t421) * MDP(22) + (-t344 * t392 - t321) * MDP(28) + (-t346 * t392 + t500) * MDP(29) + (-MDP(28) * t450 - MDP(29) * t523) * t522; t346 * t344 * MDP(23) + (-t344 ^ 2 + t346 ^ 2) * MDP(24) + (t483 + t524) * MDP(25) + (t322 + t447) * MDP(26) - t325 * MDP(27) + (t269 * t522 - t284 * t346 + t451) * MDP(28) + (-t266 * t398 + t268 * t522 - t272 * t402 + t284 * t344) * MDP(29) + (-MDP(25) * t496 - t346 * MDP(26) - t269 * MDP(28) - t268 * MDP(29)) * qJD(6);];
tauc  = t1;
