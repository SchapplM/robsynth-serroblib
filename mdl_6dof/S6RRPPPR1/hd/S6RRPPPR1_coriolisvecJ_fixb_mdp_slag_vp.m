% Calculate minimal parameter regressor of Coriolis joint torque vector for
% S6RRPPPR1
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d6,theta3,theta4]';
% MDP [27x1]
%   Minimal dynamic parameter vector (fixed base model)
%   see S6RRPPPR1_convert_par2_MPV_fixb.m
% 
% Output:
% tauc [6x1]
%   joint torques required to compensate Coriolis and centrifugal load

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 08:08
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauc = S6RRPPPR1_coriolisvecJ_fixb_mdp_slag_vp(qJ, qJD, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(10,1),zeros(27,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPPPR1_coriolisvecJ_fixb_mdp_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRPPPR1_coriolisvecJ_fixb_mdp_slag_vp: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RRPPPR1_coriolisvecJ_fixb_mdp_slag_vp: pkin has to be [10x1] (double)');
assert(isreal(MDP) && all(size(MDP) == [27 1]), ...
  'S6RRPPPR1_coriolisvecJ_fixb_mdp_slag_vp: MDP has to be [27x1] (double)'); 

%% Symbolic Calculation
% From coriolisvec_joint_fixb_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 08:08:26
% EndTime: 2019-03-09 08:08:35
% DurationCPUTime: 4.87s
% Computational Cost: add. (3893->426), mult. (10162->567), div. (0->0), fcn. (7398->8), ass. (0->183)
t432 = sin(pkin(9));
t435 = sin(qJ(2));
t437 = cos(qJ(2));
t510 = cos(pkin(9));
t411 = t432 * t437 + t435 * t510;
t431 = sin(pkin(10));
t433 = cos(pkin(10));
t434 = sin(qJ(6));
t436 = cos(qJ(6));
t520 = -t431 * t436 + t433 * t434;
t337 = t520 * t411;
t395 = t411 * qJD(1);
t361 = t433 * qJD(2) - t395 * t431;
t453 = qJD(2) * t431 + t433 * t395;
t504 = t453 * t434;
t315 = t361 * t436 + t504;
t475 = t510 * t437;
t420 = qJD(1) * t475;
t490 = qJD(1) * t435;
t392 = t432 * t490 - t420;
t483 = qJD(6) - t392;
t529 = t315 * t483;
t317 = -t361 * t434 + t436 * t453;
t528 = t317 * t483;
t527 = t361 * t392;
t485 = qJD(6) * t436;
t486 = qJD(6) * t434;
t493 = t520 * t392 + t431 * t485 - t433 * t486;
t410 = t431 * t434 + t433 * t436;
t526 = t453 ^ 2;
t482 = qJD(1) * qJD(2);
t525 = -0.2e1 * t482;
t524 = MDP(5) * (t435 ^ 2 - t437 ^ 2);
t394 = t411 * qJD(2);
t380 = qJD(1) * t394;
t478 = t435 * t482;
t381 = qJD(2) * t420 - t432 * t478;
t422 = pkin(2) * t478;
t312 = pkin(3) * t380 - qJ(4) * t381 - qJD(4) * t395 + t422;
t512 = -qJ(3) - pkin(7);
t477 = qJD(2) * t512;
t389 = qJD(3) * t437 + t435 * t477;
t373 = t389 * qJD(1);
t390 = -qJD(3) * t435 + t437 * t477;
t374 = t390 * qJD(1);
t323 = t510 * t373 + t432 * t374;
t321 = qJD(2) * qJD(4) + t323;
t318 = t431 * t321;
t279 = t312 * t433 - t318;
t274 = -pkin(4) * t380 - t279;
t479 = -pkin(2) * t437 - pkin(1);
t466 = t479 * qJD(1);
t416 = qJD(3) + t466;
t325 = pkin(3) * t392 - qJ(4) * t395 + t416;
t417 = t512 * t435;
t414 = qJD(1) * t417;
t511 = qJD(2) * pkin(2);
t407 = t414 + t511;
t418 = t512 * t437;
t415 = qJD(1) * t418;
t476 = t510 * t415;
t345 = t432 * t407 - t476;
t340 = qJD(2) * qJ(4) + t345;
t300 = t431 * t325 + t433 * t340;
t290 = t392 * qJ(5) + t300;
t522 = -t290 * t392 + t274;
t333 = pkin(2) * t490 + pkin(3) * t395 + qJ(4) * t392;
t400 = t432 * t415;
t349 = t414 * t510 + t400;
t306 = t431 * t333 + t433 * t349;
t295 = t395 * qJ(5) + t306;
t488 = qJD(4) * t433;
t521 = -t295 + t488;
t519 = MDP(15) + MDP(18);
t494 = t483 * t410;
t518 = -t380 * t520 + t483 * t494;
t344 = t510 * t407 + t400;
t452 = qJD(2) * pkin(3) - qJD(4) + t344;
t441 = qJ(5) * t453 + t452;
t298 = -pkin(4) * t361 - t441;
t425 = -pkin(2) * t510 - pkin(3);
t446 = -t431 * qJ(5) + t425;
t403 = -t433 * pkin(4) + t446;
t423 = pkin(2) * t432 + qJ(4);
t502 = t380 * t423;
t517 = -t381 * t403 + (qJD(4) - t298) * t392 + t502;
t509 = qJ(5) * t433;
t515 = pkin(4) + pkin(5);
t516 = t431 * t515 - t509;
t468 = t520 * t381;
t289 = t317 * qJD(6) + t468;
t391 = t392 ^ 2;
t514 = pkin(8) * t431;
t513 = -pkin(8) + t423;
t507 = t315 * t395;
t506 = t317 * t395;
t322 = t373 * t432 - t510 * t374;
t355 = -t510 * t417 - t418 * t432;
t505 = t322 * t355;
t369 = t431 * t381;
t370 = t433 * t381;
t438 = qJD(2) ^ 2;
t497 = t435 * t438;
t496 = t437 * t438;
t439 = qJD(1) ^ 2;
t495 = t437 * t439;
t280 = t431 * t312 + t433 * t321;
t448 = -t432 * t435 + t475;
t397 = t448 * qJD(2);
t481 = t435 * t511;
t320 = pkin(3) * t394 - qJ(4) * t397 - qJD(4) * t411 + t481;
t335 = t389 * t510 + t432 * t390;
t292 = t431 * t320 + t433 * t335;
t343 = -pkin(3) * t448 - qJ(4) * t411 + t479;
t356 = t432 * t417 - t418 * t510;
t310 = t431 * t343 + t433 * t356;
t489 = qJD(4) * t453;
t487 = qJD(5) * t431;
t303 = -qJ(5) * t448 + t310;
t480 = -t361 * t485 + t381 * t410;
t474 = pkin(1) * t525;
t268 = t318 + (-pkin(8) * t381 - t312) * t433 - t515 * t380;
t270 = t380 * qJ(5) + t392 * qJD(5) + t280;
t269 = pkin(8) * t369 + t270;
t473 = t436 * t268 - t269 * t434;
t327 = t431 * t335;
t291 = t320 * t433 - t327;
t299 = t325 * t433 - t431 * t340;
t341 = t431 * t349;
t305 = t333 * t433 - t341;
t351 = t431 * t356;
t309 = t343 * t433 - t351;
t334 = t389 * t432 - t510 * t390;
t348 = t414 * t432 - t476;
t470 = t483 ^ 2;
t469 = -t392 * t516 + t348 + t487;
t276 = t394 * qJ(5) - qJD(5) * t448 + t292;
t467 = t410 * t380 - t483 * t493;
t465 = qJD(5) - t299;
t464 = pkin(4) * t431 - t509;
t463 = t268 * t434 + t269 * t436;
t273 = -pkin(8) * t453 - t392 * t515 + t465;
t278 = -pkin(8) * t361 + t290;
t265 = t273 * t436 - t278 * t434;
t266 = t273 * t434 + t278 * t436;
t286 = t351 + (-pkin(8) * t411 - t343) * t433 + t515 * t448;
t293 = t411 * t514 + t303;
t462 = t286 * t436 - t293 * t434;
t461 = t286 * t434 + t293 * t436;
t460 = -t299 * t431 + t300 * t433;
t459 = t322 * t411 + t355 * t381;
t456 = t361 * t395 + t380 * t433;
t455 = t380 * t431 + t395 * t453;
t405 = t513 * t433;
t450 = qJD(4) * t431 - qJD(6) * t405 - t341 - (pkin(8) * t392 - t333) * t433 + t515 * t395;
t404 = t513 * t431;
t449 = qJD(6) * t404 + t392 * t514 + t521;
t338 = t410 * t411;
t447 = t453 * t486 - t480;
t445 = -qJD(5) * t411 * t433 + t334;
t444 = -pkin(4) * t369 + qJD(5) * t453 - t322;
t281 = -qJ(5) * t370 - t444;
t313 = t411 * t464 + t355;
t443 = t281 * t411 + t298 * t397 + t313 * t381;
t442 = -t397 * t452 + t459;
t440 = -t502 + t381 * t425 + (-qJD(4) - t452) * t392;
t377 = t433 * t515 - t446;
t350 = t361 * t488;
t308 = -t392 * t464 + t348;
t307 = -t411 * t516 - t355;
t304 = pkin(4) * t448 - t309;
t302 = qJD(6) * t338 + t397 * t520;
t301 = -qJD(6) * t337 + t397 * t410;
t296 = -pkin(4) * t395 - t305;
t294 = t397 * t464 + t445;
t287 = -pkin(4) * t392 + t465;
t284 = t361 * t515 + t441;
t283 = -t397 * t516 - t445;
t282 = -pkin(4) * t394 - t291;
t275 = (-pkin(5) * t431 + t509) * t381 + t444;
t272 = t397 * t514 + t276;
t271 = t327 + (-pkin(8) * t397 - t320) * t433 - t515 * t394;
t1 = [0.2e1 * t437 * MDP(4) * t478 + t524 * t525 + MDP(6) * t496 - MDP(7) * t497 + (-pkin(7) * t496 + t435 * t474) * MDP(9) + (pkin(7) * t497 + t437 * t474) * MDP(10) + (t323 * t448 + t334 * t395 - t335 * t392 - t344 * t397 - t345 * t394 - t356 * t380 + t459) * MDP(11) + (t505 + t323 * t356 - t344 * t334 + t345 * t335 + (t416 + t466) * t481) * MDP(12) + (-t279 * t448 + t291 * t392 + t299 * t394 + t309 * t380 - t334 * t361 + t431 * t442) * MDP(13) + (t280 * t448 - t292 * t392 - t300 * t394 - t310 * t380 + t334 * t453 + t433 * t442) * MDP(14) + (-t291 * t453 + t292 * t361 + (-t279 * t411 - t299 * t397 - t309 * t381) * t433 + (-t280 * t411 - t300 * t397 - t310 * t381) * t431) * MDP(15) + (t279 * t309 + t280 * t310 + t291 * t299 + t292 * t300 - t334 * t452 + t505) * MDP(16) + (t274 * t448 - t282 * t392 - t287 * t394 - t294 * t361 - t304 * t380 + t431 * t443) * MDP(17) + (t276 * t361 + t282 * t453 + (t274 * t411 + t287 * t397 + t304 * t381) * t433 + (-t270 * t411 - t290 * t397 - t303 * t381) * t431) * MDP(18) + (-t270 * t448 + t276 * t392 + t290 * t394 - t294 * t453 + t303 * t380 - t443 * t433) * MDP(19) + (t270 * t303 + t274 * t304 + t276 * t290 + t281 * t313 + t282 * t287 + t294 * t298) * MDP(20) + (t301 * t317 - t338 * t447) * MDP(21) + (-t289 * t338 - t301 * t315 - t302 * t317 + t337 * t447) * MDP(22) + (t301 * t483 - t317 * t394 - t338 * t380 - t447 * t448) * MDP(23) + (-t289 * t448 - t302 * t483 + t315 * t394 + t337 * t380) * MDP(24) + (-t380 * t448 - t394 * t483) * MDP(25) + ((t271 * t436 - t272 * t434) * t483 - t462 * t380 + t473 * t448 - t265 * t394 + t283 * t315 + t307 * t289 + t275 * t337 + t284 * t302 + (-t266 * t448 - t461 * t483) * qJD(6)) * MDP(26) + (-(t271 * t434 + t272 * t436) * t483 + t461 * t380 - t463 * t448 + t266 * t394 + t283 * t317 - t307 * t447 + t275 * t338 + t284 * t301 + (-t265 * t448 - t462 * t483) * qJD(6)) * MDP(27); -t435 * MDP(4) * t495 + t439 * t524 + ((t345 - t348) * t395 + (-t344 + t349) * t392 + (-t380 * t432 - t381 * t510) * pkin(2)) * MDP(11) + (t344 * t348 - t345 * t349 + (-t322 * t510 + t323 * t432 - t416 * t490) * pkin(2)) * MDP(12) + (-t299 * t395 - t305 * t392 - t322 * t433 + t348 * t361 + t431 * t440) * MDP(13) + (t300 * t395 + t306 * t392 + t322 * t431 - t348 * t453 + t433 * t440) * MDP(14) + (t305 * t453 - t306 * t361 + t350 + (-t299 * t392 + t280) * t433 + (-t300 * t392 - t279 + t489) * t431) * MDP(15) + (-t299 * t305 - t300 * t306 + t322 * t425 + t452 * t348 + (-t279 * t431 + t280 * t433) * t423 + t460 * qJD(4)) * MDP(16) + (-t281 * t433 + t287 * t395 + t296 * t392 + t308 * t361 + (qJD(5) * t361 - t517) * t431) * MDP(17) + (-t295 * t361 - t296 * t453 + t350 + (t287 * t392 + t270) * t433 + (t489 + t522) * t431) * MDP(18) + (-t281 * t431 - t290 * t395 - t295 * t392 + (t308 + t487) * t453 + t517 * t433) * MDP(19) + (t270 * t423 * t433 + t281 * t403 - t287 * t296 - t298 * t308 + t521 * t290 + (qJD(4) * t287 - qJD(5) * t298 + t274 * t423) * t431) * MDP(20) + (-t317 * t494 + t447 * t520) * MDP(21) + (t289 * t520 + t315 * t494 - t317 * t493 + t410 * t447) * MDP(22) + (t506 - t518) * MDP(23) + (t467 - t507) * MDP(24) + t483 * t395 * MDP(25) + (-(t404 * t436 - t405 * t434) * t380 + t377 * t289 + t275 * t410 + t265 * t395 - (t434 * t449 - t436 * t450) * t483 + t469 * t315 + t493 * t284) * MDP(26) + ((t404 * t434 + t405 * t436) * t380 - t377 * t447 - t275 * t520 - t266 * t395 - (t434 * t450 + t436 * t449) * t483 + t469 * t317 - t494 * t284) * MDP(27) + (MDP(9) * t435 * t439 + MDP(10) * t495) * pkin(1); (-t395 ^ 2 - t391) * MDP(11) + (t344 * t395 + t422) * MDP(12) + t456 * MDP(13) + (-t391 * t433 - t455) * MDP(14) + (t279 * t433 + t280 * t431 + t395 * t452) * MDP(16) + (-t391 * t431 + t456) * MDP(17) + t455 * MDP(19) + (t270 * t431 - t274 * t433 - t298 * t395) * MDP(20) + (t467 + t507) * MDP(26) + (t506 + t518) * MDP(27) + t519 * t381 * (-t431 ^ 2 - t433 ^ 2) + (t345 * MDP(12) + t460 * MDP(16) + (t287 * t431 + t290 * t433) * MDP(20) + (-MDP(13) * t431 + MDP(19) * t433) * t392 + t519 * (t361 * t433 + t431 * t453)) * t392; (t299 * t453 - t300 * t361 + t322) * MDP(16) + (-t287 * t453 - t290 * t361 + t281) * MDP(20) + (-t289 - t528) * MDP(26) + (t447 + t529) * MDP(27) + (MDP(13) + MDP(17)) * (t392 * t453 + t369) + (MDP(14) - MDP(19)) * (t370 + t527) + t519 * (-t361 ^ 2 - t526); (-qJD(2) * t395 - t361 * t453) * MDP(17) + (t370 - t527) * MDP(18) + (-t391 - t526) * MDP(19) + (t298 * t453 + t522) * MDP(20) + (-t315 * t453 - t436 * t380 - t434 * t470) * MDP(26) + (-t317 * t453 + t434 * t380 - t436 * t470) * MDP(27); t317 * t315 * MDP(21) + (-t315 ^ 2 + t317 ^ 2) * MDP(22) + (t480 + t529) * MDP(23) + (-t468 + t528) * MDP(24) - t380 * MDP(25) + (t266 * t483 - t284 * t317 + t473) * MDP(26) + (t265 * t483 + t284 * t315 - t463) * MDP(27) + (-MDP(23) * t504 - MDP(24) * t317 - MDP(26) * t266 - MDP(27) * t265) * qJD(6);];
tauc  = t1;
