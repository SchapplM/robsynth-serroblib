% Calculate minimal parameter regressor of Coriolis joint torque vector for
% S6RPRRPP4
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d3,d4,theta2,theta5]';
% MDP [27x1]
%   Minimal dynamic parameter vector (fixed base model)
%   see S6RPRRPP4_convert_par2_MPV_fixb.m
% 
% Output:
% tauc [6x1]
%   joint torques required to compensate Coriolis and centrifugal load

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 04:41
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauc = S6RPRRPP4_coriolisvecJ_fixb_mdp_slag_vp(qJ, qJD, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(10,1),zeros(27,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRRPP4_coriolisvecJ_fixb_mdp_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RPRRPP4_coriolisvecJ_fixb_mdp_slag_vp: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RPRRPP4_coriolisvecJ_fixb_mdp_slag_vp: pkin has to be [10x1] (double)');
assert(isreal(MDP) && all(size(MDP) == [27 1]), ...
  'S6RPRRPP4_coriolisvecJ_fixb_mdp_slag_vp: MDP has to be [27x1] (double)'); 

%% Symbolic Calculation
% From coriolisvec_joint_fixb_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 04:40:55
% EndTime: 2019-03-09 04:41:04
% DurationCPUTime: 5.43s
% Computational Cost: add. (7419->435), mult. (19113->563), div. (0->0), fcn. (14372->8), ass. (0->172)
t441 = sin(pkin(9));
t443 = cos(pkin(9));
t445 = sin(qJ(3));
t447 = cos(qJ(3));
t416 = t441 * t447 + t443 * t445;
t440 = sin(pkin(10));
t442 = cos(pkin(10));
t444 = sin(qJ(4));
t446 = cos(qJ(4));
t523 = -t440 * t444 + t442 * t446;
t357 = t523 * t416;
t479 = t446 * qJD(3);
t522 = t416 * qJD(1);
t381 = t444 * t522 - t479;
t383 = qJD(3) * t444 + t446 * t522;
t331 = t442 * t381 + t383 * t440;
t485 = qJD(1) * t447;
t430 = t443 * t485;
t486 = qJD(1) * t445;
t469 = t441 * t486;
t405 = t430 - t469;
t398 = qJD(4) - t405;
t532 = t331 * t398;
t531 = (t441 ^ 2 + t443 ^ 2) * (MDP(7) * qJ(2) + MDP(6));
t415 = t440 * t446 + t442 * t444;
t404 = t415 * qJD(4);
t494 = t415 * t405 - t404;
t481 = qJD(4) * t446;
t482 = qJD(4) * t444;
t493 = t523 * t405 + t440 * t482 - t442 * t481;
t530 = t416 * qJD(2);
t414 = t441 * t445 - t447 * t443;
t408 = t414 * qJD(3);
t468 = t416 * t481;
t529 = -t408 * t444 + t468;
t460 = -t381 * t440 + t442 * t383;
t528 = t460 ^ 2;
t520 = pkin(7) + qJ(2);
t421 = t520 * t441;
t418 = qJD(1) * t421;
t422 = t520 * t443;
t419 = qJD(1) * t422;
t525 = -t418 * t447 - t445 * t419;
t366 = -qJD(3) * pkin(3) - t525;
t329 = pkin(4) * t381 + qJD(5) + t366;
t300 = pkin(5) * t331 - qJ(6) * t460 + t329;
t527 = t300 * t460;
t464 = t398 * t446;
t368 = pkin(3) * t522 - pkin(8) * t405;
t360 = t446 * t368;
t313 = -qJ(5) * t405 * t446 + pkin(4) * t522 - t444 * t525 + t360;
t490 = t444 * t368 + t446 * t525;
t509 = t405 * t444;
t321 = -qJ(5) * t509 + t490;
t519 = -qJ(5) - pkin(8);
t466 = qJD(4) * t519;
t400 = qJD(5) * t446 + t444 * t466;
t449 = -qJD(5) * t444 + t446 * t466;
t492 = (-t313 + t449) * t442 + (t321 - t400) * t440;
t377 = t421 * t447 + t445 * t422;
t373 = -t445 * t418 + t447 * t419;
t524 = -t373 + (t482 - t509) * pkin(4);
t521 = MDP(22) + MDP(25);
t435 = -pkin(2) * t443 - pkin(1);
t420 = qJD(1) * t435 + qJD(2);
t345 = -pkin(3) * t405 - pkin(8) * t522 + t420;
t367 = qJD(3) * pkin(8) + t373;
t324 = t345 * t444 + t367 * t446;
t312 = -qJ(5) * t381 + t324;
t308 = t442 * t312;
t323 = t446 * t345 - t367 * t444;
t311 = -qJ(5) * t383 + t323;
t286 = t311 * t440 + t308;
t517 = t286 * t460;
t516 = t312 * t440;
t426 = qJD(3) * t430;
t393 = -t469 * qJD(3) + t426;
t488 = qJD(4) * t479 + t446 * t393;
t341 = -t482 * t522 + t488;
t515 = t341 * t444;
t514 = t381 * t398;
t513 = t381 * t522;
t512 = t383 * t398;
t511 = t383 * t522;
t510 = t393 * t444;
t507 = t416 * t444;
t506 = t416 * t446;
t409 = t416 * qJD(3);
t394 = qJD(1) * t409;
t499 = t444 * t394;
t378 = -t421 * t445 + t422 * t447;
t374 = t446 * t378;
t387 = t446 * t394;
t451 = t414 * qJD(2);
t337 = -qJD(1) * t451 + qJD(3) * t525;
t354 = pkin(3) * t394 - pkin(8) * t393;
t351 = t446 * t354;
t465 = -t337 * t444 + t351;
t278 = pkin(4) * t394 - qJ(5) * t341 - qJD(4) * t324 - qJD(5) * t383 + t465;
t342 = t383 * qJD(4) + t510;
t473 = t446 * t337 + t345 * t481 + t444 * t354;
t453 = -t367 * t482 + t473;
t283 = -qJ(5) * t342 - qJD(5) * t381 + t453;
t498 = -t442 * t278 + t440 * t283;
t269 = t440 * t278 + t442 * t283;
t294 = t440 * t313 + t442 * t321;
t289 = qJ(6) * t522 + t294;
t353 = t442 * t400 + t440 * t449;
t497 = -t289 + t353;
t496 = -pkin(5) * t522 + t492;
t346 = -t377 * qJD(3) - t451;
t369 = pkin(3) * t409 + pkin(8) * t408;
t361 = t446 * t369;
t371 = pkin(3) * t414 - pkin(8) * t416 + t435;
t458 = qJ(5) * t408 - qJD(5) * t416;
t291 = pkin(4) * t409 - t346 * t444 + t361 + t458 * t446 + (-t374 + (qJ(5) * t416 - t371) * t444) * qJD(4);
t472 = t446 * t346 + t444 * t369 + t371 * t481;
t295 = -qJ(5) * t468 + (-qJD(4) * t378 + t458) * t444 + t472;
t273 = t440 * t291 + t442 * t295;
t495 = t494 * pkin(5) - t493 * qJ(6) + qJD(6) * t415 - t524;
t304 = pkin(4) * t398 + t311;
t285 = t440 * t304 + t308;
t363 = t446 * t371;
t318 = pkin(4) * t414 - qJ(5) * t506 - t378 * t444 + t363;
t489 = t444 * t371 + t374;
t325 = -qJ(5) * t507 + t489;
t299 = t440 * t318 + t442 * t325;
t491 = t398 * t509 + t387;
t484 = qJD(3) * t445;
t483 = qJD(3) * t447;
t480 = qJD(6) * t398;
t287 = t311 * t442 - t516;
t478 = qJD(6) - t287;
t476 = qJD(1) * qJD(2);
t474 = t394 * qJ(6) + t269;
t470 = -pkin(4) * t446 - pkin(3);
t467 = t519 * t444;
t314 = t341 * t440 + t442 * t342;
t338 = t530 * qJD(1) - t418 * t484 + t419 * t483;
t347 = -t421 * t484 + t422 * t483 + t530;
t267 = -pkin(5) * t394 + t498;
t315 = t341 * t442 - t342 * t440;
t423 = t519 * t446;
t379 = -t423 * t440 - t442 * t467;
t380 = -t442 * t423 + t440 * t467;
t463 = -t380 * t314 + t315 * t379 - t353 * t331;
t462 = pkin(4) * t507 + t377;
t272 = t291 * t442 - t295 * t440;
t284 = t304 * t442 - t516;
t298 = t318 * t442 - t325 * t440;
t456 = t529 * pkin(4) + t347;
t455 = -t408 * t446 - t416 * t482;
t317 = pkin(4) * t342 + t338;
t454 = -pkin(8) * t394 + t366 * t398;
t274 = pkin(5) * t314 - qJ(6) * t315 - qJD(6) * t460 + t317;
t434 = -pkin(4) * t442 - pkin(5);
t431 = pkin(4) * t440 + qJ(6);
t370 = -pkin(5) * t523 - qJ(6) * t415 + t470;
t356 = t415 * t416;
t327 = t404 * t416 + t523 * t408;
t326 = -qJD(4) * t357 + t408 * t415;
t306 = pkin(5) * t356 - qJ(6) * t357 + t462;
t302 = pkin(4) * t383 + pkin(5) * t460 + qJ(6) * t331;
t297 = -pkin(5) * t414 - t298;
t296 = qJ(6) * t414 + t299;
t280 = qJ(6) * t398 + t285;
t279 = -pkin(5) * t398 + qJD(6) - t284;
t275 = -pkin(5) * t326 + qJ(6) * t327 - qJD(6) * t357 + t456;
t271 = -pkin(5) * t409 - t272;
t270 = qJ(6) * t409 + qJD(6) * t414 + t273;
t266 = t474 + t480;
t1 = [(t393 * t416 - t408 * t522) * MDP(8) + (-t393 * t414 - t394 * t416 - t405 * t408 - t409 * t522) * MDP(9) + (t394 * t435 + t409 * t420) * MDP(13) + (t393 * t435 - t408 * t420) * MDP(14) + (t341 * t506 + t383 * t455) * MDP(15) + (-(-t381 * t446 - t383 * t444) * t408 + (-t515 - t342 * t446 + (t381 * t444 - t383 * t446) * qJD(4)) * t416) * MDP(16) + (t341 * t414 + t383 * t409 + t387 * t416 + t398 * t455) * MDP(17) + (-t342 * t414 - t381 * t409 - t529 * t398 - t416 * t499) * MDP(18) + (t394 * t414 + t398 * t409) * MDP(19) + ((-t378 * t481 + t361) * t398 + t363 * t394 + (-t367 * t481 + t351) * t414 + t323 * t409 + t347 * t381 + t377 * t342 + t366 * t468 + ((-qJD(4) * t371 - t346) * t398 - t378 * t394 + (-qJD(4) * t345 - t337) * t414 + t338 * t416 - t366 * t408) * t444) * MDP(20) + (-(-t378 * t482 + t472) * t398 - t489 * t394 - t453 * t414 - t324 * t409 + t347 * t383 + t377 * t341 + t338 * t506 + t455 * t366) * MDP(21) + (-t269 * t356 - t272 * t460 - t273 * t331 + t284 * t327 + t285 * t326 - t298 * t315 - t299 * t314 + t357 * t498) * MDP(22) + (t269 * t299 + t284 * t272 + t285 * t273 - t298 * t498 + t317 * t462 + t329 * t456) * MDP(23) + (-t267 * t414 - t271 * t398 + t274 * t356 + t275 * t331 - t279 * t409 - t297 * t394 - t300 * t326 + t306 * t314) * MDP(24) + (-t266 * t356 + t267 * t357 - t270 * t331 + t271 * t460 - t279 * t327 + t280 * t326 - t296 * t314 + t297 * t315) * MDP(25) + (t266 * t414 + t270 * t398 - t274 * t357 - t275 * t460 + t280 * t409 + t296 * t394 + t300 * t327 - t306 * t315) * MDP(26) + (t266 * t296 + t267 * t297 + t270 * t280 + t271 * t279 + t274 * t306 + t275 * t300) * MDP(27) + 0.2e1 * t476 * t531 + (-MDP(10) * t408 - MDP(11) * t409 - MDP(13) * t347 - MDP(14) * t346) * qJD(3); t426 * MDP(14) + (t491 - t513) * MDP(20) + (-t499 - t511) * MDP(21) + (t269 * t415 + t284 * t494 - t285 * t493 - t329 * t522 - t498 * t523) * MDP(23) + (-t331 * t522 + t394 * t523) * MDP(24) + t415 * t394 * MDP(26) + (t266 * t415 - t267 * t523 - t279 * t494 - t280 * t493 - t300 * t522) * MDP(27) + (MDP(26) * t522 - t521 * t494) * t460 + (-MDP(20) * t482 - MDP(21) * t464 + MDP(24) * t494 - MDP(26) * t493) * t398 + ((t441 * t485 + t443 * t486 + t522) * MDP(13) + (t405 - t469) * MDP(14)) * qJD(3) - qJD(1) ^ 2 * t531 + t521 * (-t415 * t314 - t315 * t523 + t493 * t331); -t405 ^ 2 * MDP(9) + (t426 + (-t405 - t469) * qJD(3)) * MDP(10) + (qJD(3) * t373 - t338) * MDP(13) + (-t405 * t420 + t414 * t476) * MDP(14) + (t383 * t464 + t515) * MDP(15) + ((t341 - t514) * t446 + (-t342 - t512) * t444) * MDP(16) + (t398 * t464 + t499 - t511) * MDP(17) + (-t398 * t482 + t491 + t513) * MDP(18) + (-pkin(3) * t342 - t338 * t446 - t373 * t381 + (-pkin(8) * t481 - t360) * t398 + (t398 * t525 + t454) * t444) * MDP(20) + (-pkin(3) * t341 + t338 * t444 - t373 * t383 + (pkin(8) * t482 + t490) * t398 + t454 * t446) * MDP(21) + (t269 * t523 + t284 * t493 + t285 * t494 + t294 * t331 + t415 * t498 - t460 * t492 + t463) * MDP(22) + (t269 * t380 + t498 * t379 + t317 * t470 + t524 * t329 + (t353 - t294) * t285 + t492 * t284) * MDP(23) + (-t274 * t523 - t300 * t494 + t314 * t370 - t331 * t495 - t379 * t394 + t398 * t496) * MDP(24) + (t266 * t523 + t267 * t415 - t279 * t493 + t280 * t494 + t289 * t331 - t460 * t496 + t463) * MDP(25) + (-t274 * t415 + t300 * t493 - t315 * t370 + t380 * t394 + t398 * t497 + t460 * t495) * MDP(26) + (t266 * t380 + t267 * t379 + t274 * t370 - t279 * t496 + t280 * t497 - t300 * t495) * MDP(27) + (-MDP(13) * t420 - MDP(19) * t398 - MDP(20) * t323 + MDP(21) * t324 + MDP(24) * t279 - MDP(26) * t280 - MDP(8) * t405 + MDP(9) * t522) * t522; t383 * t381 * MDP(15) + (-t381 ^ 2 + t383 ^ 2) * MDP(16) + (t488 + t514) * MDP(17) + (-t510 + t512) * MDP(18) + (t324 * t398 - t366 * t383 + t465) * MDP(20) + (t323 * t398 + t366 * t381 - t473) * MDP(21) + (t285 * t460 - t517) * MDP(22) + (t284 * t286 - t285 * t287) * MDP(23) + (t286 * t398 - t498 - t527) * MDP(24) + (t280 * t460 - t314 * t431 + t315 * t434 - t517) * MDP(25) + (-t287 * t398 + t302 * t460 + t474 + 0.2e1 * t480) * MDP(26) + (t266 * t431 + t267 * t434 - t279 * t286 + t280 * t478 - t300 * t302) * MDP(27) + (MDP(19) + (pkin(5) - t434) * MDP(24) + t431 * MDP(26)) * t394 + ((-MDP(18) * t522 - MDP(20) * t367) * t446 + (-MDP(17) * t522 - MDP(18) * qJD(3) - MDP(20) * t345 + MDP(21) * t367) * t444) * qJD(4) + ((-t314 * t440 - t315 * t442) * MDP(22) + (t269 * t440 - t329 * t383 - t442 * t498) * MDP(23)) * pkin(4) + ((-t284 + t287) * MDP(22) - t302 * MDP(24) + (t279 - t478) * MDP(25) - t300 * MDP(26)) * t331; (t284 * t460 + t285 * t331 + t317) * MDP(23) + (t398 * t460 + t314) * MDP(24) + (-t315 + t532) * MDP(26) + (-t279 * t460 + t280 * t331 + t274) * MDP(27) + t521 * (-t331 ^ 2 - t528); (-qJD(3) * t522 + t331 * t460) * MDP(24) + (t315 + t532) * MDP(25) + (-t398 ^ 2 - t528) * MDP(26) + (-t280 * t398 + t267 + t527) * MDP(27);];
tauc  = t1;
