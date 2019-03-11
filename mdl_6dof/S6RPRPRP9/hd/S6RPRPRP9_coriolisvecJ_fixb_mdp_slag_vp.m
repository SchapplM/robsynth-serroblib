% Calculate minimal parameter regressor of Coriolis joint torque vector for
% S6RPRPRP9
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d3,d5,theta4]';
% MDP [28x1]
%   Minimal dynamic parameter vector (fixed base model)
%   see S6RPRPRP9_convert_par2_MPV_fixb.m
% 
% Output:
% tauc [6x1]
%   joint torques required to compensate Coriolis and centrifugal load

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 03:30
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauc = S6RPRPRP9_coriolisvecJ_fixb_mdp_slag_vp(qJ, qJD, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(9,1),zeros(28,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRPRP9_coriolisvecJ_fixb_mdp_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RPRPRP9_coriolisvecJ_fixb_mdp_slag_vp: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S6RPRPRP9_coriolisvecJ_fixb_mdp_slag_vp: pkin has to be [9x1] (double)');
assert(isreal(MDP) && all(size(MDP) == [28 1]), ...
  'S6RPRPRP9_coriolisvecJ_fixb_mdp_slag_vp: MDP has to be [28x1] (double)'); 

%% Symbolic Calculation
% From coriolisvec_joint_fixb_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 03:29:42
% EndTime: 2019-03-09 03:29:48
% DurationCPUTime: 4.02s
% Computational Cost: add. (3884->408), mult. (8351->543), div. (0->0), fcn. (5427->6), ass. (0->168)
t380 = sin(pkin(9));
t385 = cos(qJ(3));
t448 = qJD(1) * t385;
t422 = t380 * t448;
t381 = cos(pkin(9));
t444 = qJD(3) * t381;
t349 = -t422 + t444;
t429 = t381 * t448;
t445 = qJD(3) * t380;
t350 = t429 + t445;
t382 = sin(qJ(5));
t384 = cos(qJ(5));
t298 = -t384 * t349 + t350 * t382;
t383 = sin(qJ(3));
t449 = qJD(1) * t383;
t375 = qJD(5) + t449;
t483 = t298 * t375;
t430 = t380 * t449;
t464 = t384 * t381;
t439 = qJD(5) * t384;
t440 = qJD(5) * t382;
t476 = -t380 * t440 + t381 * t439;
t453 = -t382 * t430 + t449 * t464 + t476;
t353 = t380 * t384 + t381 * t382;
t338 = t353 * qJD(1);
t340 = t353 * qJD(5);
t452 = t383 * t338 + t340;
t469 = t349 * t382;
t398 = t350 * t384 + t469;
t482 = t398 ^ 2;
t433 = 0.2e1 * qJD(1);
t481 = MDP(8) * (t383 ^ 2 - t385 ^ 2);
t480 = MDP(6) * qJ(2) + MDP(5);
t407 = pkin(3) * t385 + qJ(4) * t383;
t355 = t407 * qJD(1);
t386 = -pkin(1) - pkin(7);
t475 = qJD(1) * t386;
t370 = qJD(2) + t475;
t467 = t380 * t385;
t308 = t381 * t355 - t370 * t467;
t432 = pkin(8) * t381 * t383;
t392 = (pkin(4) * t385 + t432) * qJD(1);
t281 = t308 + t392;
t466 = t381 * t385;
t309 = t380 * t355 + t370 * t466;
t290 = pkin(8) * t430 + t309;
t352 = t380 * t382 - t464;
t472 = pkin(8) + qJ(4);
t366 = t472 * t380;
t367 = t472 * t381;
t397 = -t366 * t384 - t367 * t382;
t479 = qJD(4) * t352 - qJD(5) * t397 + t382 * t281 + t384 * t290;
t313 = -t366 * t382 + t367 * t384;
t478 = -qJD(4) * t353 - qJD(5) * t313 - t281 * t384 + t290 * t382;
t330 = t352 * t385;
t457 = -qJD(3) * t330 - t340 * t383 - t338;
t328 = t353 * t385;
t329 = t352 * t383;
t456 = -t352 * qJD(1) + qJD(3) * t328 - qJD(5) * t329;
t360 = pkin(3) * t383 - qJ(4) * t385 + qJ(2);
t347 = t381 * t360;
t418 = -t380 * t386 + pkin(4);
t296 = -pkin(8) * t466 + t383 * t418 + t347;
t465 = t383 * t386;
t316 = t380 * t360 + t381 * t465;
t307 = -pkin(8) * t467 + t316;
t477 = t382 * t296 + t384 * t307;
t435 = MDP(23) + MDP(25);
t434 = MDP(24) - MDP(27);
t387 = qJD(3) ^ 2;
t474 = qJD(2) * t433 - t386 * t387;
t331 = qJD(3) * t407 - qJD(4) * t385 + qJD(2);
t319 = t381 * t331;
t278 = t319 + (t385 * t418 + t432) * qJD(3);
t441 = qJD(3) * t386;
t425 = t385 * t441;
t305 = t380 * t331 + t381 * t425;
t443 = qJD(3) * t383;
t428 = t380 * t443;
t288 = pkin(8) * t428 + t305;
t473 = -qJD(5) * t477 + t278 * t384 - t288 * t382;
t471 = qJD(3) * pkin(3);
t468 = t370 * t385;
t359 = t383 * t370;
t388 = qJD(1) ^ 2;
t463 = t385 * t388;
t461 = qJ(6) * t448 + t479;
t460 = -pkin(5) * t448 + t478;
t322 = -pkin(4) * t430 + t359;
t459 = -t452 * pkin(5) + t453 * qJ(6) + qJD(6) * t353 + t322;
t314 = t331 * qJD(1);
t333 = (qJD(4) + t468) * qJD(3);
t276 = t380 * t314 + t381 * t333;
t343 = t360 * qJD(1);
t344 = qJD(3) * qJ(4) + t359;
t292 = t380 * t343 + t381 * t344;
t436 = qJD(1) * qJD(3);
t421 = t383 * t436;
t411 = t382 * t421;
t454 = t349 * t439 + t380 * t411;
t410 = t384 * t421;
t451 = -t380 * t410 - t381 * t411;
t447 = qJD(3) * t397;
t446 = qJD(3) * t313;
t442 = qJD(3) * t385;
t438 = qJD(6) * t375;
t291 = t381 * t343 - t344 * t380;
t268 = pkin(4) * t449 - pkin(8) * t350 + t291;
t272 = pkin(8) * t349 + t292;
t249 = t268 * t384 - t272 * t382;
t437 = qJD(6) - t249;
t377 = -pkin(4) * t381 - pkin(3);
t427 = t382 * t443;
t426 = t384 * t443;
t374 = t383 * t441;
t420 = t385 * t436;
t419 = MDP(22) * t448;
t417 = -qJD(4) + t471;
t275 = t381 * t314 - t333 * t380;
t335 = -t417 - t468;
t416 = t335 - t468;
t348 = pkin(4) * t467 - t385 * t386;
t415 = pkin(5) * t420;
t261 = qJD(3) * t392 + t275;
t412 = t380 * t421;
t269 = pkin(8) * t412 + t276;
t414 = -t384 * t261 + t268 * t440 + t382 * t269 + t272 * t439;
t413 = qJ(6) * t420;
t406 = pkin(4) * t412;
t250 = t268 * t382 + t272 * t384;
t405 = -t275 * t380 + t276 * t381;
t402 = -t291 * t381 - t292 * t380;
t401 = -t291 * t380 + t292 * t381;
t400 = t296 * t384 - t307 * t382;
t336 = -pkin(4) * t428 + t374;
t396 = t250 * t375 - t414;
t395 = -t382 * t261 - t268 * t439 - t384 * t269 + t272 * t440;
t394 = t382 * t278 + t384 * t288 + t296 * t439 - t307 * t440;
t393 = -t335 + (t370 + t475) * t385;
t306 = -pkin(4) * t349 + t335;
t266 = t350 * t440 + t381 * t410 - t454;
t267 = qJD(5) * t398 + t451;
t354 = t370 * t443;
t391 = pkin(5) * t267 + qJ(6) * t266 - qJD(6) * t398 + t354;
t242 = t414 - t415;
t390 = t249 * t375 + t395;
t389 = -qJ(4) * t442 + (t335 + t417) * t383;
t327 = t353 * t383;
t317 = t354 - t406;
t315 = -t380 * t465 + t347;
t304 = -t380 * t425 + t319;
t295 = pkin(5) * t352 - qJ(6) * t353 + t377;
t287 = -t380 * t426 - t381 * t427 + t385 * t476;
t285 = t340 * t385 - t380 * t427 + t381 * t426;
t271 = pkin(5) * t328 + qJ(6) * t330 + t348;
t257 = pkin(5) * t398 + qJ(6) * t298;
t256 = -pkin(5) * t383 - t400;
t255 = qJ(6) * t383 + t477;
t252 = -t266 + t483;
t251 = pkin(5) * t298 - qJ(6) * t398 + t306;
t248 = pkin(5) * t287 + qJ(6) * t285 + qJD(6) * t330 + t336;
t247 = qJ(6) * t375 + t250;
t246 = -pkin(5) * t375 + t437;
t245 = -pkin(5) * t442 - t473;
t244 = qJ(6) * t442 + qJD(6) * t383 + t394;
t243 = t391 - t406;
t241 = -t395 + t413 + t438;
t1 = [0.2e1 * t436 * t481 + (-t304 * t350 + t305 * t349 + ((t315 * t381 + t316 * t380) * qJD(1) - t402) * t443) * MDP(16) + (t275 * t315 + t276 * t316 + t291 * t304 + t292 * t305 + t374 * t416) * MDP(17) + (t266 * t330 - t285 * t398) * MDP(18) + (t266 * t328 + t267 * t330 + t285 * t298 - t287 * t398) * MDP(19) + (-t285 * t375 + (-qJD(1) * t330 + t398) * t442) * MDP(20) + (-t287 * t375 + (-qJD(1) * t328 - t298) * t442) * MDP(21) + (t375 + t449) * MDP(22) * t442 + (t473 * t375 + t336 * t298 + t348 * t267 + t317 * t328 + t306 * t287 + (qJD(1) * t400 + t249) * t442) * MDP(23) + (-t394 * t375 + t336 * t398 - t348 * t266 - t317 * t330 - t306 * t285 + (-qJD(1) * t477 - t250) * t442) * MDP(24) + (t243 * t328 - t245 * t375 + t248 * t298 + t251 * t287 + t267 * t271 + (-qJD(1) * t256 - t246) * t442) * MDP(25) + (-t241 * t328 - t242 * t330 - t244 * t298 + t245 * t398 - t246 * t285 - t247 * t287 - t255 * t267 - t256 * t266) * MDP(26) + (t243 * t330 + t244 * t375 - t248 * t398 + t251 * t285 + t266 * t271 + (qJD(1) * t255 + t247) * t442) * MDP(27) + (t241 * t255 + t242 * t256 + t243 * t271 + t244 * t247 + t245 * t246 + t248 * t251) * MDP(28) + (-t387 * MDP(10) + t474 * MDP(13) + (-t275 * t381 - t276 * t380) * MDP(16) + ((qJD(1) * t315 + t291) * MDP(14) + (-qJD(1) * t316 - t292) * MDP(15)) * qJD(3)) * t385 + (-0.2e1 * MDP(7) * t420 - t387 * MDP(9) + t474 * MDP(12) + (qJD(1) * t304 + t275 + (-t349 * t386 + t380 * t393) * qJD(3)) * MDP(14) + (-qJD(1) * t305 - t276 + (t350 * t386 + t381 * t393) * qJD(3)) * MDP(15) - t266 * MDP(20) - t267 * MDP(21) - t414 * MDP(23) + t395 * MDP(24) - t242 * MDP(25) + t241 * MDP(27)) * t383 + (t480 * qJD(2) + (MDP(12) * t442 - MDP(13) * t443) * qJ(2)) * t433; (-t266 * t327 + t267 * t329 - t298 * t457 + t398 * t456) * MDP(26) + (-t241 * t329 + t242 * t327 + t246 * t456 + t247 * t457) * MDP(28) - t480 * t388 + t435 * t298 * t443 + ((-t349 * t380 + t350 * t381) * MDP(16) + t402 * MDP(17)) * qJD(1) + (-t434 * t457 - t435 * t456) * t375 + ((-t387 - t388) * MDP(13) - t243 * MDP(28) - t435 * t267 + t434 * t266 + ((t349 * t381 + t350 * t380) * MDP(16) + t401 * MDP(17) + (-t327 * t435 + t329 * t434) * qJD(1)) * qJD(3)) * t385 + (-t387 * MDP(12) + t405 * MDP(17) + (-MDP(14) * t381 + MDP(15) * t380 - MDP(12)) * t388 + ((-t349 - t422) * MDP(14) + (t350 - t429) * MDP(15) + t416 * MDP(17) + t251 * MDP(28) + t434 * t398) * qJD(3)) * t383; t383 * MDP(7) * t463 - t388 * t481 + ((t349 - t444) * t359 + (-t291 * t385 - t308 * t383 + t380 * t389) * qJD(1)) * MDP(14) + ((-t350 + t445) * t359 + (t292 * t385 + t309 * t383 + t381 * t389) * qJD(1)) * MDP(15) + (t308 * t350 - t309 * t349 + (qJD(4) * t349 - t291 * t449 + t276) * t381 + (qJD(4) * t350 - t292 * t449 - t275) * t380) * MDP(16) + (-t291 * t308 - t292 * t309 + (-t335 - t471) * t359 + t401 * qJD(4) + t405 * qJ(4)) * MDP(17) + (-t266 * t353 + t398 * t453) * MDP(18) + (t266 * t352 - t267 * t353 - t298 * t453 - t398 * t452) * MDP(19) + (t453 * t375 + (qJD(3) * t353 - t398) * t448) * MDP(20) + (-t452 * t375 + (-qJD(3) * t352 + t298) * t448) * MDP(21) - t375 * t419 + (t377 * t267 - t322 * t298 + t317 * t352 + t478 * t375 + t452 * t306 + (-t249 + t447) * t448) * MDP(23) + (-t377 * t266 - t322 * t398 + t317 * t353 + t479 * t375 + t453 * t306 + (t250 - t446) * t448) * MDP(24) + (t243 * t352 + t267 * t295 + t460 * t375 - t459 * t298 + t452 * t251 + (t246 + t447) * t448) * MDP(25) + (-t241 * t352 + t242 * t353 + t246 * t453 - t247 * t452 + t266 * t397 - t267 * t313 + t298 * t461 - t398 * t460) * MDP(26) + (-t243 * t353 + t266 * t295 - t461 * t375 + t459 * t398 - t453 * t251 + (-t247 + t446) * t448) * MDP(27) + (t241 * t313 - t242 * t397 + t243 * t295 - t246 * t460 - t247 * t461 - t251 * t459) * MDP(28) + (MDP(13) * t383 * t388 - MDP(12) * t463) * qJ(2); (-t349 ^ 2 - t350 ^ 2) * MDP(16) + (t291 * t350 - t292 * t349 + t354) * MDP(17) + (-t298 ^ 2 - t482) * MDP(26) + (-t246 * t398 + t247 * t298 + t391) * MDP(28) + (t435 * t469 + (-t382 * t434 + t384 * t435) * t350) * qJD(5) + (t350 * MDP(14) + t349 * MDP(15) + ((-MDP(28) * pkin(4) - MDP(14)) * t380 + (-t384 * t434 - MDP(15)) * t381) * qJD(3)) * t449 + t435 * (t375 * t398 + t451) + t434 * (t454 - t483); t252 * MDP(20) + (-t349 * t440 - t350 * t439 - t451) * MDP(21) + qJD(3) * t419 + t396 * MDP(23) + t390 * MDP(24) + (t396 + 0.2e1 * t415) * MDP(25) + (pkin(5) * t266 - qJ(6) * t267) * MDP(26) + (-t390 + 0.2e1 * t413 + 0.2e1 * t438) * MDP(27) + (-pkin(5) * t242 + qJ(6) * t241 - t246 * t250 + t247 * t437 - t251 * t257) * MDP(28) + (t375 * MDP(21) - t306 * MDP(23) - t251 * MDP(25) + (t247 - t250) * MDP(26) + t257 * MDP(27) + MDP(19) * t398) * t398 + (t398 * MDP(18) + t306 * MDP(24) - t257 * MDP(25) + (t246 - t437) * MDP(26) - t251 * MDP(27) - MDP(19) * t298) * t298; (t298 * t398 - t420) * MDP(25) + t252 * MDP(26) + (-t375 ^ 2 - t482) * MDP(27) + (-t247 * t375 + t251 * t398 + t242) * MDP(28);];
tauc  = t1;
