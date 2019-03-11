% Calculate minimal parameter regressor of Coriolis joint torque vector for
% S6RPRPRR4
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d3,d5,d6,theta2]';
% MDP [29x1]
%   Minimal dynamic parameter vector (fixed base model)
%   see S6RPRPRR4_convert_par2_MPV_fixb.m
% 
% Output:
% tauc [6x1]
%   joint torques required to compensate Coriolis and centrifugal load

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 03:46
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauc = S6RPRPRR4_coriolisvecJ_fixb_mdp_slag_vp(qJ, qJD, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(10,1),zeros(29,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRPRR4_coriolisvecJ_fixb_mdp_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RPRPRR4_coriolisvecJ_fixb_mdp_slag_vp: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RPRPRR4_coriolisvecJ_fixb_mdp_slag_vp: pkin has to be [10x1] (double)');
assert(isreal(MDP) && all(size(MDP) == [29 1]), ...
  'S6RPRPRR4_coriolisvecJ_fixb_mdp_slag_vp: MDP has to be [29x1] (double)'); 

%% Symbolic Calculation
% From coriolisvec_joint_fixb_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 03:46:34
% EndTime: 2019-03-09 03:46:42
% DurationCPUTime: 3.96s
% Computational Cost: add. (2269->378), mult. (5035->520), div. (0->0), fcn. (3164->8), ass. (0->189)
t389 = sin(pkin(10)) * pkin(1) + pkin(7);
t373 = t389 * qJD(1);
t408 = sin(qJ(3));
t411 = cos(qJ(3));
t479 = t411 * qJD(2) - t408 * t373;
t518 = qJD(4) - t479;
t511 = pkin(4) + t389;
t410 = cos(qJ(5));
t407 = sin(qJ(5));
t473 = qJD(3) * t407;
t474 = qJD(1) * t411;
t364 = t410 * t474 + t473;
t409 = cos(qJ(6));
t451 = t407 * t474;
t471 = qJD(3) * t410;
t366 = -t451 + t471;
t406 = sin(qJ(6));
t499 = t366 * t406;
t293 = t409 * t364 + t499;
t476 = qJD(1) * t408;
t388 = qJD(5) + t476;
t381 = qJD(6) + t388;
t521 = t293 * t381;
t433 = t364 * t406 - t409 * t366;
t520 = t381 * t433;
t323 = -qJD(3) * pkin(3) + t518;
t466 = qJD(5) * t411;
t453 = t407 * t466;
t519 = t408 * t471 + t453;
t513 = qJD(5) + qJD(6);
t412 = -pkin(3) - pkin(8);
t463 = pkin(4) * t476 + t518;
t302 = qJD(3) * t412 + t463;
t390 = -cos(pkin(10)) * pkin(1) - pkin(2);
t425 = -qJ(4) * t408 + t390;
t336 = t411 * t412 + t425;
t314 = t336 * qJD(1);
t278 = t302 * t407 + t314 * t410;
t275 = -pkin(9) * t364 + t278;
t465 = qJD(6) * t406;
t273 = t275 * t465;
t340 = t408 * qJD(2) + t411 * t373;
t321 = pkin(4) * t474 + t340;
t401 = qJD(3) * qJ(4);
t311 = t401 + t321;
t286 = pkin(5) * t364 + t311;
t517 = t286 * t293 + t273;
t462 = qJD(1) * qJD(3);
t450 = t408 * t462;
t318 = -qJD(5) * t364 + t407 * t450;
t387 = pkin(3) * t450;
t435 = pkin(8) * t408 - qJ(4) * t411;
t469 = qJD(4) * t408;
t417 = qJD(3) * t435 - t469;
t305 = qJD(1) * t417 + t387;
t461 = qJD(2) * qJD(3);
t470 = qJD(3) * t411;
t330 = t373 * t470 + t408 * t461;
t449 = t411 * t462;
t306 = pkin(4) * t449 + t330;
t445 = -t305 * t407 + t410 * t306;
t416 = -qJD(5) * t278 + t445;
t265 = pkin(5) * t449 - pkin(9) * t318 + t416;
t467 = qJD(5) * t410;
t319 = qJD(3) * t467 - qJD(5) * t451 - t410 * t450;
t459 = -t302 * t467 - t410 * t305 - t407 * t306;
t468 = qJD(5) * t407;
t420 = -t314 * t468 - t459;
t266 = -pkin(9) * t319 + t420;
t446 = t409 * t265 - t406 * t266;
t516 = t286 * t433 + t446;
t515 = MDP(27) * t449 + (-t293 ^ 2 + t433 ^ 2) * MDP(24) - t293 * MDP(23) * t433;
t402 = t408 ^ 2;
t403 = t411 ^ 2;
t514 = MDP(6) * (t402 - t403);
t354 = t511 * t408;
t344 = t407 * t354;
t484 = t410 * t336 + t344;
t326 = -t401 - t340;
t444 = t318 * t406 + t409 * t319;
t271 = -qJD(6) * t433 + t444;
t512 = t411 * t513;
t510 = pkin(9) - t412;
t509 = t271 * t408;
t277 = t410 * t302 - t314 * t407;
t274 = -pkin(9) * t366 + t277;
t272 = pkin(5) * t388 + t274;
t508 = t272 * t409;
t507 = t275 * t409;
t367 = t406 * t410 + t407 * t409;
t421 = t367 * t408;
t493 = t409 * t410;
t496 = t406 * t407;
t432 = -t493 + t496;
t280 = qJD(3) * t421 + t432 * t512;
t506 = t280 * t381;
t400 = qJD(3) * qJD(4);
t472 = qJD(3) * t408;
t480 = -t373 * t472 + t411 * t461;
t317 = -t400 - t480;
t298 = -pkin(4) * t450 - t317;
t505 = t298 * t407;
t504 = t298 * t410;
t503 = t318 * t410;
t502 = t319 * t408;
t501 = t364 * t388;
t500 = t366 * t388;
t498 = t388 * t410;
t497 = t388 * t412;
t495 = t407 * t408;
t413 = qJD(3) ^ 2;
t494 = t408 * t413;
t492 = t410 * t411;
t491 = t411 * t413;
t464 = qJD(6) * t409;
t458 = t409 * t318 - t406 * t319 - t364 * t464;
t270 = -t366 * t465 + t458;
t490 = t270 * t408 - t433 * t470;
t281 = t367 * t512 - t432 * t472;
t342 = t432 * t411;
t489 = t281 * t381 + t342 * t449;
t303 = t513 * t367;
t488 = -t303 * t381 - t432 * t449;
t333 = qJD(1) * t421;
t487 = -t303 - t333;
t475 = qJD(1) * t410;
t455 = t408 * t475;
t486 = -t406 * t468 - t407 * t465 + t409 * t455 - t476 * t496 + t493 * t513;
t485 = t318 * t408 + t366 * t470;
t396 = pkin(3) * t476;
t341 = qJD(1) * t435 + t396;
t483 = t407 * t321 + t410 * t341;
t482 = t519 * t388;
t355 = t511 * t411;
t457 = -pkin(5) * t410 - pkin(4);
t481 = pkin(5) * t467 - t457 * t476 + t518;
t353 = -pkin(3) * t411 + t425;
t327 = qJD(1) * t353;
t374 = qJD(1) * t390;
t414 = qJD(1) ^ 2;
t460 = t408 * t411 * t414;
t456 = t403 * t475;
t452 = t410 * t466;
t371 = t510 * t410;
t448 = pkin(9) * t411 - t336;
t443 = t410 * t321 - t341 * t407;
t395 = pkin(3) * t472;
t324 = t395 + t417;
t349 = t511 * t470;
t442 = -t324 * t407 + t410 * t349;
t441 = qJD(6) * t272 + t266;
t440 = t388 * t452;
t377 = t408 * t449;
t439 = -qJD(3) * t479 + t480;
t438 = qJD(3) * t340 - t330;
t370 = t510 * t407;
t426 = pkin(5) * t411 - pkin(9) * t495;
t437 = qJD(1) * t426 - qJD(6) * t370 - t510 * t468 + t443;
t436 = pkin(9) * t455 + t371 * t513 + t483;
t263 = t272 * t406 + t507;
t345 = t410 * t354;
t284 = pkin(5) * t408 + t407 * t448 + t345;
t285 = -pkin(9) * t492 + t484;
t434 = t284 * t406 + t285 * t409;
t431 = qJD(1) * t403 - t388 * t408;
t430 = -0.2e1 * qJD(3) * t327;
t429 = 0.2e1 * qJD(3) * t374;
t428 = t388 * t407;
t424 = t311 * t408 + t412 * t470;
t422 = -qJ(4) * t470 - t469;
t328 = qJD(1) * t422 + t387;
t350 = t395 + t422;
t423 = qJD(1) * t350 + t389 * t413 + t328;
t419 = t410 * t324 - t336 * t468 + t407 * t349 + t354 * t467;
t418 = t326 * MDP(15) - t364 * MDP(21) - t293 * MDP(28);
t415 = -t317 * t411 + t330 * t408 + (t323 * t411 + t326 * t408) * qJD(3);
t391 = pkin(5) * t407 + qJ(4);
t378 = t410 * t449;
t369 = -qJ(4) * t474 + t396;
t348 = t511 * t472;
t343 = t367 * t411;
t329 = pkin(5) * t492 + t355;
t315 = t327 * t476;
t291 = -pkin(5) * t453 + (-t389 + t457) * t472;
t282 = pkin(5) * t319 + t298;
t268 = pkin(9) * t519 + t419;
t267 = t426 * qJD(3) + (t410 * t448 - t344) * qJD(5) + t442;
t262 = -t275 * t406 + t508;
t1 = [-0.2e1 * t462 * t514 + (-t389 * t491 + t408 * t429) * MDP(10) + (t389 * t494 + t411 * t429) * MDP(11) + t415 * MDP(12) + (t408 * t430 + t411 * t423) * MDP(13) + (-t408 * t423 + t411 * t430) * MDP(14) + (t327 * t350 + t328 * t353 + t389 * t415) * MDP(15) + (-t318 * t407 * t411 + (t407 * t472 - t452) * t366) * MDP(16) + ((-t364 * t407 + t366 * t410) * t472 + (-t503 + t319 * t407 + (t364 * t410 + t366 * t407) * qJD(5)) * t411) * MDP(17) + (-t431 * t473 - t440 + t485) * MDP(18) + (-t502 + (-t364 * t411 - t456) * qJD(3) + t482) * MDP(19) + (t388 * t470 + t377) * MDP(20) + (t442 * t388 - t348 * t364 + t355 * t319 + (-t311 * t471 + t445) * t408 + (-t278 * t408 - t388 * t484) * qJD(5) + (-t311 * t468 + t504 + ((-t336 * t407 + t345) * qJD(1) + t277) * qJD(3)) * t411) * MDP(21) + (-t419 * t388 - t348 * t366 + t355 * t318 + ((qJD(3) * t311 + qJD(5) * t314) * t407 + t459) * t408 + (-t311 * t467 - t505 + (-qJD(1) * t484 - t278) * qJD(3)) * t411) * MDP(22) + (-t270 * t343 - t280 * t433) * MDP(23) + (t270 * t342 + t271 * t343 - t280 * t293 - t281 * t433) * MDP(24) + (-t343 * t449 + t490 + t506) * MDP(25) + (-t293 * t470 + t489 - t509) * MDP(26) + (t381 * t470 + t377) * MDP(27) + ((t267 * t409 - t268 * t406) * t381 + t446 * t408 + t291 * t293 + t329 * t271 - t282 * t342 - t286 * t281 + (-t263 * t408 - t381 * t434) * qJD(6) + ((t284 * t409 - t285 * t406) * qJD(1) + t262) * t470) * MDP(28) + (t329 * t270 + t273 * t408 + t286 * t280 - t282 * t343 - t291 * t433 + (-(-qJD(6) * t285 + t267) * t381 - t265 * t408) * t406 + (-(qJD(6) * t284 + t268) * t381 - t441 * t408) * t409 + (-qJD(1) * t434 - t263) * t470) * MDP(29) + MDP(7) * t491 - MDP(8) * t494 + 0.2e1 * MDP(5) * t377; (-t317 * t408 - t330 * t411) * MDP(15) + (t482 + t502) * MDP(21) + (t440 + t485) * MDP(22) + (t489 + t509) * MDP(28) + (t490 - t506) * MDP(29) + (-MDP(21) * t456 + t323 * t408 * MDP(15) + (MDP(29) * qJD(1) * t343 - t418) * t411 + t431 * MDP(22) * t407) * qJD(3) + ((-MDP(11) + MDP(14)) * t411 + (-MDP(10) + MDP(13)) * t408) * t413; -MDP(5) * t460 + t414 * t514 + (-t374 * t476 + t438) * MDP(10) - t439 * MDP(11) + (t315 - t438) * MDP(13) + (0.2e1 * t400 + (t327 * t411 + t369 * t408) * qJD(1) + t439) * MDP(14) + (-pkin(3) * t330 - qJ(4) * t317 - t323 * t340 - t326 * t518 - t327 * t369) * MDP(15) + (-t366 * t428 + t503) * MDP(16) + ((-t319 - t500) * t410 + (-t318 + t501) * t407) * MDP(17) + (-t388 * t468 + t378 + (-t366 * t411 - t388 * t495) * qJD(1)) * MDP(18) + (-t388 * t467 + (-t408 * t498 + (t364 - t473) * t411) * qJD(1)) * MDP(19) + (qJ(4) * t319 + t505 - t443 * t388 + t463 * t364 + (t311 * t410 - t407 * t497) * qJD(5) + (-t277 * t411 + t410 * t424) * qJD(1)) * MDP(21) + (qJ(4) * t318 + t504 + t483 * t388 + t463 * t366 + (-t311 * t407 - t410 * t497) * qJD(5) + (t278 * t411 - t407 * t424) * qJD(1)) * MDP(22) + (-t270 * t432 - t433 * t487) * MDP(23) + (-t270 * t367 + t271 * t432 - t293 * t487 + t433 * t486) * MDP(24) + (-t333 * t381 + t488) * MDP(25) - t486 * t381 * MDP(26) + (t391 * t271 + t282 * t367 + (t406 * t436 - t409 * t437) * t381 + t481 * t293 + t486 * t286) * MDP(28) + (t391 * t270 - t282 * t432 + (t406 * t437 + t409 * t436) * t381 - t481 * t433 + t487 * t286) * MDP(29) + (-t374 * MDP(11) - t369 * MDP(13) - t388 * MDP(20) + t433 * MDP(25) + (-qJD(3) * t367 + t293) * MDP(26) - t381 * MDP(27) + ((t370 * t406 - t371 * t409) * qJD(3) - t262) * MDP(28) + (-(-t370 * t409 - t371 * t406) * qJD(3) + t263) * MDP(29)) * t474; MDP(13) * t460 + (-t402 * t414 - t413) * MDP(14) + (t315 + t330) * MDP(15) + t378 * MDP(21) + t488 * MDP(28) + (-t333 * MDP(28) - MDP(29) * t486) * t381 + (-MDP(21) * t428 - MDP(22) * t498) * t388 + ((-t366 - t451) * MDP(22) + (-t367 * t474 + t433) * MDP(29) + t418) * qJD(3); t366 * t364 * MDP(16) + (-t364 ^ 2 + t366 ^ 2) * MDP(17) + (t318 + t501) * MDP(18) + (-t319 + t500) * MDP(19) + MDP(20) * t449 + (t278 * t388 - t311 * t366 + t416) * MDP(21) + (t277 * t388 + t311 * t364 - t420) * MDP(22) + (t270 + t521) * MDP(25) + (-t271 - t520) * MDP(26) + (-(-t274 * t406 - t507) * t381 - t263 * qJD(6) + (-t293 * t366 - t381 * t465 + t409 * t449) * pkin(5) + t516) * MDP(28) + ((-t275 * t381 - t265) * t406 + (t274 * t381 - t441) * t409 + (t366 * t433 - t381 * t464 - t406 * t449) * pkin(5) + t517) * MDP(29) + t515; (t458 + t521) * MDP(25) + (-t444 - t520) * MDP(26) + (t263 * t381 + t516) * MDP(28) + (t262 * t381 - t406 * t265 - t409 * t266 + t517) * MDP(29) + (-MDP(25) * t499 + MDP(26) * t433 - MDP(28) * t263 - MDP(29) * t508) * qJD(6) + t515;];
tauc  = t1;
