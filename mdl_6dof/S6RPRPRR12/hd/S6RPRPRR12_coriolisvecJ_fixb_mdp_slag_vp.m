% Calculate minimal parameter regressor of Coriolis joint torque vector for
% S6RPRPRR12
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d3,d5,d6]';
% MDP [31x1]
%   Minimal dynamic parameter vector (fixed base model)
%   see S6RPRPRR12_convert_par2_MPV_fixb.m
% 
% Output:
% tauc [6x1]
%   joint torques required to compensate Coriolis and centrifugal load

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 04:20
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauc = S6RPRPRR12_coriolisvecJ_fixb_mdp_slag_vp(qJ, qJD, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(9,1),zeros(31,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRPRR12_coriolisvecJ_fixb_mdp_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RPRPRR12_coriolisvecJ_fixb_mdp_slag_vp: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S6RPRPRR12_coriolisvecJ_fixb_mdp_slag_vp: pkin has to be [9x1] (double)');
assert(isreal(MDP) && all(size(MDP) == [31 1]), ...
  'S6RPRPRR12_coriolisvecJ_fixb_mdp_slag_vp: MDP has to be [31x1] (double)'); 

%% Symbolic Calculation
% From coriolisvec_joint_fixb_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 04:20:04
% EndTime: 2019-03-09 04:20:11
% DurationCPUTime: 4.60s
% Computational Cost: add. (2208->399), mult. (4566->545), div. (0->0), fcn. (2780->6), ass. (0->185)
t390 = cos(qJ(3));
t457 = qJD(1) * t390;
t370 = qJD(5) + t457;
t361 = qJD(6) + t370;
t385 = sin(qJ(6));
t386 = sin(qJ(5));
t388 = cos(qJ(6));
t389 = cos(qJ(5));
t342 = t385 * t389 + t386 * t388;
t494 = qJD(5) + qJD(6);
t503 = t494 * t342;
t276 = t503 * t361;
t392 = -pkin(1) - pkin(7);
t492 = pkin(4) - t392;
t387 = sin(qJ(3));
t458 = qJD(1) * t387;
t433 = t389 * t458;
t337 = qJD(3) * t386 - t433;
t455 = qJD(3) * t389;
t339 = t386 * t458 + t455;
t477 = t339 * t385;
t280 = t388 * t337 + t477;
t505 = t280 * t361;
t411 = t337 * t385 - t388 * t339;
t504 = t361 * t411;
t362 = qJD(1) * t392 + qJD(2);
t346 = t390 * t362;
t495 = qJD(4) - t346;
t447 = pkin(4) * t457 + t495;
t391 = -pkin(3) - pkin(8);
t453 = qJD(3) * t391;
t295 = t447 + t453;
t488 = qJ(4) * t390;
t414 = pkin(8) * t387 - t488;
t462 = pkin(3) * t458 + qJD(1) * qJ(2);
t303 = qJD(1) * t414 + t462;
t271 = t295 * t386 + t303 * t389;
t266 = -pkin(9) * t337 + t271;
t449 = qJD(6) * t385;
t264 = t266 * t449;
t345 = t387 * t362;
t317 = -pkin(4) * t458 + t345;
t381 = qJD(3) * qJ(4);
t305 = t317 + t381;
t277 = pkin(5) * t337 + t305;
t502 = t277 * t280 + t264;
t446 = qJD(1) * qJD(3);
t430 = t390 * t446;
t451 = qJD(5) * t386;
t297 = -qJD(3) * t451 + qJD(5) * t433 + t386 * t430;
t413 = (qJD(3) * pkin(8) - qJD(4)) * t390;
t380 = qJD(1) * qJD(2);
t431 = t387 * t446;
t440 = pkin(3) * t430 + qJ(4) * t431 + t380;
t284 = qJD(1) * t413 + t440;
t456 = qJD(3) * t387;
t340 = t362 * t456;
t310 = -pkin(4) * t431 + t340;
t425 = -t284 * t386 + t389 * t310;
t397 = -t271 * qJD(5) + t425;
t257 = -pkin(5) * t431 - pkin(9) * t297 + t397;
t298 = t339 * qJD(5) - t389 * t430;
t450 = qJD(5) * t389;
t442 = -t389 * t284 - t295 * t450 - t386 * t310;
t402 = -t303 * t451 - t442;
t258 = -pkin(9) * t298 + t402;
t426 = t388 * t257 - t385 * t258;
t501 = t277 * t411 + t426;
t471 = t388 * t389;
t475 = t385 * t386;
t410 = -t471 + t475;
t500 = t410 * t458 - t280;
t379 = qJD(3) * qJD(4);
t454 = qJD(3) * t390;
t321 = -t362 * t454 - t379;
t490 = qJD(3) * pkin(3);
t427 = -qJD(4) + t490;
t324 = -t346 - t427;
t328 = -t345 - t381;
t481 = t328 * t390;
t499 = (t387 * (-t324 + t346) + t481) * qJD(3) + t321 * t387;
t429 = MDP(29) * t456;
t498 = (-t280 ^ 2 + t411 ^ 2) * MDP(26) - qJD(1) * t429 - t280 * MDP(25) * t411;
t383 = t387 ^ 2;
t384 = t390 ^ 2;
t497 = MDP(8) * (t383 - t384);
t403 = t361 * t410;
t496 = t387 * pkin(3) + qJ(2);
t330 = t414 + t496;
t351 = t492 * t390;
t331 = t386 * t351;
t464 = t389 * t330 + t331;
t424 = t297 * t385 + t388 * t298;
t262 = -qJD(6) * t411 + t424;
t491 = pkin(9) - t391;
t394 = qJD(1) ^ 2;
t489 = qJ(2) * t394;
t270 = t389 * t295 - t303 * t386;
t265 = -pkin(9) * t339 + t270;
t263 = pkin(5) * t370 + t265;
t487 = t263 * t388;
t486 = t266 * t388;
t485 = t297 * t389;
t300 = -pkin(4) * t430 - t321;
t484 = t300 * t386;
t483 = t300 * t389;
t480 = t337 * t370;
t479 = t337 * t387;
t478 = t339 * t370;
t476 = t370 * t391;
t474 = t386 * t387;
t473 = t386 * t390;
t472 = t387 * t389;
t470 = t389 * t390;
t469 = t390 * t394;
t393 = qJD(3) ^ 2;
t468 = t392 * t393;
t309 = t342 * t457;
t467 = -t503 - t309;
t438 = t389 * t457;
t466 = -t385 * t451 - t386 * t449 + t388 * t438 - t457 * t475 + t494 * t471;
t344 = pkin(3) * t457 + qJ(4) * t458;
t314 = pkin(8) * t457 + t344;
t465 = t389 * t314 + t386 * t317;
t439 = -pkin(5) * t389 - pkin(4);
t463 = pkin(5) * t450 - t439 * t457 + t495;
t400 = qJD(3) * t342;
t452 = qJD(4) * t390;
t448 = qJD(6) * t388;
t445 = 0.2e1 * qJD(1);
t443 = t387 * t469;
t441 = t388 * t297 - t385 * t298 - t337 * t448;
t437 = t386 * t454;
t436 = t389 * t454;
t435 = t387 * t451;
t434 = t370 * t450;
t432 = pkin(3) * t454 + qJ(4) * t456 + qJD(2);
t350 = t491 * t389;
t428 = -pkin(9) * t387 - t330;
t299 = t413 + t432;
t333 = t492 * t456;
t423 = -t299 * t386 - t389 * t333;
t422 = -t314 * t386 + t389 * t317;
t327 = -qJ(4) * t457 + t462;
t347 = -t488 + t496;
t420 = qJD(1) * t347 + t327;
t419 = t370 + t457;
t418 = qJD(6) * t263 + t258;
t334 = t492 * t454;
t348 = t491 * t386;
t405 = -pkin(5) * t387 - pkin(9) * t473;
t416 = qJD(1) * t405 - qJD(6) * t348 - t491 * t451 + t422;
t415 = pkin(9) * t438 + t350 * t494 + t465;
t255 = t263 * t385 + t486;
t332 = t389 * t351;
t274 = pkin(5) * t390 + t386 * t428 + t332;
t275 = pkin(9) * t472 + t464;
t412 = t274 * t385 + t275 * t388;
t409 = t419 * t387;
t408 = -qJD(1) * t383 + t370 * t390;
t407 = t370 * (qJD(5) * t390 + qJD(1));
t406 = t370 * t386;
t296 = -qJD(1) * t452 + t440;
t312 = t432 - t452;
t404 = -qJD(1) * t312 - t296 + t468;
t401 = t389 * t299 - t330 * t451 - t386 * t333 + t351 * t450;
t261 = -t339 * t449 + t441;
t399 = t342 * qJD(1);
t398 = -t435 + t436;
t374 = t387 * t392;
t371 = pkin(5) * t386 + qJ(4);
t356 = t386 * t431;
t349 = -pkin(4) * t387 + t374;
t319 = t387 * t439 + t374;
t316 = t342 * t387;
t315 = t385 * t474 - t387 * t471;
t311 = t327 * t457;
t287 = -pkin(5) * t398 - t334;
t273 = pkin(5) * t298 + t300;
t268 = t385 * t437 + t387 * t503 - t388 * t436;
t267 = -t494 * t387 * t410 + t390 * t400;
t260 = pkin(9) * t398 + t401;
t259 = t405 * qJD(3) + (t389 * t428 - t331) * qJD(5) + t423;
t254 = -t266 * t385 + t487;
t1 = [(-t370 * t435 - t298 * t390 + (t389 * t408 + t479) * qJD(3)) * MDP(21) + t499 * MDP(14) + (t423 * t370 - t334 * t337 + t349 * t298 + (-t305 * t455 + t425) * t390 + (-t271 * t390 - t464 * t370) * qJD(5) + (t305 * t451 - t483 + (-(-t330 * t386 + t332) * qJD(1) - t270) * qJD(3)) * t387) * MDP(23) + (-t401 * t370 - t334 * t339 + t349 * t297 + ((qJD(3) * t305 + qJD(5) * t303) * t386 + t442) * t390 + (t305 * t450 + t484 + (qJD(1) * t464 + t271) * qJD(3)) * t387) * MDP(24) + ((-t337 * t386 + t339 * t389) * t454 + (t485 - t298 * t386 + (-t337 * t389 - t339 * t386) * qJD(5)) * t387) * MDP(19) + (t297 * t474 + (t387 * t450 + t437) * t339) * MDP(18) + (-t390 * t468 + (-qJ(2) * t456 + qJD(2) * t390) * t445) * MDP(13) + (-t387 * t468 + (qJ(2) * t454 + qJD(2) * t387) * t445) * MDP(12) + ((t259 * t388 - t260 * t385) * t361 + t426 * t390 + t287 * t280 + t319 * t262 + t273 * t315 + t277 * t268 + (-t255 * t390 - t361 * t412) * qJD(6) + (-(t274 * t388 - t275 * t385) * qJD(1) - t254) * t456) * MDP(30) + (t261 * t390 + t267 * t361 + (-qJD(1) * t316 + t411) * t456) * MDP(27) + (-t262 * t390 - t268 * t361 + (qJD(1) * t315 + t280) * t456) * MDP(28) + (t319 * t261 + t264 * t390 + t277 * t267 + t273 * t316 - t287 * t411 + (-(-qJD(6) * t275 + t259) * t361 - t257 * t390) * t385 + (-(qJD(6) * t274 + t260) * t361 - t418 * t390) * t388 + (qJD(1) * t412 + t255) * t456) * MDP(31) + (t390 * t404 + t420 * t456) * MDP(16) + (-t361 - t457) * t429 + (t387 * t404 - t420 * t454) * MDP(15) + (t387 * t434 + t297 * t390 + (-t339 * t387 + t386 * t408) * qJD(3)) * MDP(20) - qJD(3) * MDP(22) * t409 + (-t261 * t315 - t262 * t316 - t267 * t280 + t268 * t411) * MDP(26) + (t261 * t316 - t267 * t411) * MDP(25) - 0.2e1 * t387 * MDP(7) * t430 + 0.2e1 * t446 * t497 + (t296 * t347 + t312 * t327 - t392 * t499) * MDP(17) + 0.2e1 * (MDP(6) * qJ(2) + MDP(5)) * t380 + (-MDP(10) * t390 - MDP(9) * t387) * t393; -t394 * MDP(5) - MDP(6) * t489 + (-qJD(1) * t327 - t499) * MDP(17) + (t298 * t387 + t386 * t407 + (t337 * t390 + t419 * t472) * qJD(3)) * MDP(23) + (t297 * t387 + t389 * t407 + (t339 * t390 - t386 * t409) * qJD(3)) * MDP(24) + (t361 * t399 + (-qJD(3) * t403 + t262) * t387 + (-qJD(3) * t500 + t276) * t390) * MDP(30) + (-qJD(1) * t403 + (-t361 * t400 + t261) * t387 + ((-t387 * t399 - t411) * qJD(3) - t494 * t403) * t390) * MDP(31) + ((-MDP(13) + MDP(16)) * t390 + (-MDP(12) + MDP(15)) * t387) * (t393 + t394); MDP(7) * t443 - t394 * t497 - qJ(2) * MDP(12) * t469 + t387 * MDP(13) * t489 + ((-t328 - t381) * t390 + (t324 + t427) * t387) * qJD(1) * MDP(14) + t311 * MDP(15) + (0.2e1 * t379 + (-t327 * t387 + t344 * t390) * qJD(1)) * MDP(16) + (-qJ(4) * t321 - qJD(4) * t328 - t327 * t344 + (t481 + (-t324 - t490) * t387) * t362) * MDP(17) + (-t339 * t406 + t485) * MDP(18) + ((-t298 - t478) * t389 + (-t297 + t480) * t386) * MDP(19) + (-t370 * t451 + (-t370 * t473 + (t339 - t455) * t387) * qJD(1)) * MDP(20) + (-t434 + t356 + (-t370 * t470 - t479) * qJD(1)) * MDP(21) + (qJ(4) * t298 + t484 - t422 * t370 + t447 * t337 + (t305 * t389 - t386 * t476) * qJD(5) + (t305 * t470 + (-t389 * t453 + t270) * t387) * qJD(1)) * MDP(23) + (qJ(4) * t297 + t483 + t465 * t370 + t447 * t339 + (-t305 * t386 - t389 * t476) * qJD(5) + (-t271 * t387 + (-t305 * t390 + t387 * t453) * t386) * qJD(1)) * MDP(24) + (-t261 * t410 - t411 * t467) * MDP(25) + (-t261 * t342 + t262 * t410 - t280 * t467 + t411 * t466) * MDP(26) + (-t309 * t361 - t276) * MDP(27) - t466 * t361 * MDP(28) + (t371 * t262 + t273 * t342 + (t385 * t415 - t388 * t416) * t361 + t463 * t280 + t466 * t277) * MDP(30) + (t371 * t261 - t273 * t410 + (t385 * t416 + t388 * t415) * t361 - t463 * t411 + t467 * t277) * MDP(31) + (t344 * MDP(15) + t370 * MDP(22) + (qJD(3) * t410 - t411) * MDP(27) + (-t280 + t400) * MDP(28) + t361 * MDP(29) + (-(t348 * t385 - t350 * t388) * qJD(3) + t254) * MDP(30) + ((-t348 * t388 - t350 * t385) * qJD(3) - t255) * MDP(31)) * t458; -MDP(15) * t443 + (-t384 * t394 - t393) * MDP(16) + (t311 + t340) * MDP(17) + t356 * MDP(24) - t276 * MDP(30) + (-t309 * MDP(30) - MDP(31) * t466) * t361 + (-MDP(24) * t370 * t389 - MDP(23) * t406) * t370 + (t328 * MDP(17) + (-t337 - t433) * MDP(23) - t339 * MDP(24) + t500 * MDP(30) + (t342 * t458 + t411) * MDP(31)) * qJD(3); t339 * t337 * MDP(18) + (-t337 ^ 2 + t339 ^ 2) * MDP(19) + (t297 + t480) * MDP(20) + (-t298 + t478) * MDP(21) - MDP(22) * t431 + (t271 * t370 - t305 * t339 + t397) * MDP(23) + (t270 * t370 + t305 * t337 - t402) * MDP(24) + (t261 + t505) * MDP(27) + (-t262 - t504) * MDP(28) + (-(-t265 * t385 - t486) * t361 - t255 * qJD(6) + (-t280 * t339 - t361 * t449 - t388 * t431) * pkin(5) + t501) * MDP(30) + ((-t266 * t361 - t257) * t385 + (t265 * t361 - t418) * t388 + (t339 * t411 - t361 * t448 + t385 * t431) * pkin(5) + t502) * MDP(31) + t498; (t441 + t505) * MDP(27) + (-t424 - t504) * MDP(28) + (t255 * t361 + t501) * MDP(30) + (t254 * t361 - t385 * t257 - t388 * t258 + t502) * MDP(31) + (-MDP(27) * t477 + MDP(28) * t411 - MDP(30) * t255 - MDP(31) * t487) * qJD(6) + t498;];
tauc  = t1;
