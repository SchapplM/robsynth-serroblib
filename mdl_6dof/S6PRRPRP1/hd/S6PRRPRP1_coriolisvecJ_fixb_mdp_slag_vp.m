% Calculate Coriolis joint torque vector for
% S6PRRPRP1
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [11x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d2,d3,d5,theta1,theta4]';
% MDP [26x1]
%   Minimal dynamic parameter vector (fixed base model)
%   see S6PRRPRP1_convert_par2_MPV_fixb.m
% 
% Output:
% tauc [6x1]
%   joint torques required to compensate Coriolis and centrifugal load

% Quelle: HybrDyn-Toolbox
% Datum: 2021-01-16 02:40
% Revision: 24b2e7d74a0c1a3b64fa2f8f5ad758691ad61af3 (2021-01-15)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauc = S6PRRPRP1_coriolisvecJ_fixb_mdp_slag_vp(qJ, qJD, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(11,1),zeros(26,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRRPRP1_coriolisvecJ_fixb_mdp_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6PRRPRP1_coriolisvecJ_fixb_mdp_slag_vp: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6PRRPRP1_coriolisvecJ_fixb_mdp_slag_vp: pkin has to be [11x1] (double)');
assert(isreal(MDP) && all(size(MDP) == [26 1]), ...
  'S6PRRPRP1_coriolisvecJ_fixb_mdp_slag_vp: MDP has to be [26x1] (double)'); 

%% Symbolic Calculation
% From coriolisvec_joint_fixb_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2021-01-16 02:39:01
% EndTime: 2021-01-16 02:39:15
% DurationCPUTime: 5.17s
% Computational Cost: add. (4299->417), mult. (11031->559), div. (0->0), fcn. (8353->10), ass. (0->183)
t380 = sin(qJ(2));
t376 = sin(pkin(6));
t438 = qJD(1) * t376;
t425 = t380 * t438;
t379 = sin(qJ(3));
t434 = qJD(3) * t379;
t487 = pkin(3) * t434 - t425;
t382 = cos(qJ(3));
t473 = qJ(4) + pkin(8);
t418 = qJD(3) * t473;
t343 = qJD(4) * t382 - t379 * t418;
t375 = sin(pkin(11));
t392 = -qJD(4) * t379 - t382 * t418;
t471 = cos(pkin(11));
t302 = t343 * t471 + t375 * t392;
t416 = t471 * t382;
t396 = -t375 * t379 + t416;
t383 = cos(qJ(2));
t424 = t383 * t438;
t326 = t396 * t424;
t441 = t302 - t326;
t355 = t375 * t382 + t379 * t471;
t347 = t355 * qJD(3);
t350 = t396 * qJD(3);
t486 = pkin(4) * t347 - pkin(9) * t350 + t487;
t348 = t355 * qJD(2);
t378 = sin(qJ(5));
t381 = cos(qJ(5));
t332 = qJD(3) * t378 + t348 * t381;
t463 = t332 * t378;
t485 = MDP(5) * t382;
t484 = MDP(6) * (t379 ^ 2 - t382 ^ 2);
t364 = qJD(2) * t416;
t436 = qJD(2) * t379;
t345 = t375 * t436 - t364;
t342 = qJD(5) + t345;
t483 = t342 * t463;
t377 = cos(pkin(6));
t437 = qJD(1) * t377;
t365 = t382 * t437;
t357 = qJD(2) * pkin(8) + t425;
t412 = qJ(4) * qJD(2) + t357;
t323 = -t379 * t412 + t365;
t319 = qJD(3) * pkin(3) + t323;
t423 = t379 * t437;
t324 = t382 * t412 + t423;
t417 = t471 * t324;
t276 = t375 * t319 + t417;
t273 = qJD(3) * pkin(9) + t276;
t371 = -pkin(3) * t382 - pkin(2);
t341 = qJD(2) * t371 + qJD(4) - t424;
t287 = pkin(4) * t345 - pkin(9) * t348 + t341;
t260 = t273 * t381 + t287 * t378;
t407 = qJD(4) + t424;
t291 = (-t357 * t379 + t365) * qJD(3) + (-qJ(4) * t434 + t382 * t407) * qJD(2);
t292 = (-t357 * t382 - t423) * qJD(3) + (-qJ(4) * qJD(3) * t382 - t379 * t407) * qJD(2);
t265 = t291 * t471 + t375 * t292;
t340 = qJD(2) * t347;
t430 = qJD(2) * qJD(3);
t419 = t379 * t430;
t435 = qJD(2) * t380;
t422 = t376 * t435;
t344 = pkin(3) * t419 + qJD(1) * t422;
t362 = t375 * t419;
t393 = t364 * qJD(3) - t362;
t284 = t340 * pkin(4) - pkin(9) * t393 + t344;
t282 = t381 * t284;
t388 = -qJD(5) * t260 - t265 * t378 + t282;
t431 = t381 * qJD(3);
t433 = qJD(5) * t378;
t402 = qJD(5) * t431 - t348 * t433 + t381 * t393;
t386 = -qJ(6) * t402 + t388;
t475 = pkin(5) * t340;
t245 = -qJD(6) * t332 + t386 + t475;
t330 = t348 * t378 - t431;
t254 = -qJ(6) * t330 + t260;
t470 = t254 * t342;
t482 = t245 + t470;
t481 = t326 * t378 + t486 * t381;
t442 = t343 * t375 - t355 * t424 - t471 * t392;
t314 = -pkin(4) * t396 - pkin(9) * t355 + t371;
t432 = qJD(5) * t381;
t480 = t314 * t432 + t486 * t378 + t441 * t381;
t479 = MDP(21) + MDP(23);
t478 = MDP(22) + MDP(24);
t294 = t332 * qJD(5) + t378 * t393;
t477 = t332 ^ 2;
t476 = pkin(3) * t379;
t474 = t330 * pkin(5);
t472 = qJD(2) * pkin(2);
t469 = t402 * t378;
t468 = t330 * t342;
t467 = t330 * t345;
t466 = t330 * t348;
t465 = t332 * t342;
t464 = t332 * t348;
t462 = t340 * t378;
t460 = t345 * t378;
t459 = t345 * t381;
t458 = t355 * t378;
t457 = t355 * t381;
t315 = t375 * t324;
t456 = t376 * t380;
t455 = t376 * t383;
t385 = qJD(2) ^ 2;
t454 = t376 * t385;
t384 = qJD(3) ^ 2;
t453 = t379 * t384;
t360 = t473 * t379;
t361 = t473 * t382;
t328 = -t375 * t360 + t361 * t471;
t320 = t381 * t328;
t336 = t381 * t340;
t452 = t382 * t384;
t368 = pkin(3) * t375 + pkin(9);
t451 = qJ(6) + t368;
t404 = -qJ(6) * t350 - qJD(6) * t355;
t450 = pkin(5) * t347 - t302 * t378 + t404 * t381 + (-t320 + (qJ(6) * t355 - t314) * t378) * qJD(5) + t481;
t420 = t355 * t432;
t449 = -qJ(6) * t420 + (-qJD(5) * t328 + t404) * t378 + t480;
t259 = -t273 * t378 + t381 * t287;
t253 = -qJ(6) * t332 + t259;
t251 = pkin(5) * t342 + t253;
t448 = t251 - t253;
t280 = t323 * t471 - t315;
t429 = pkin(3) * t436;
t303 = pkin(4) * t348 + pkin(9) * t345 + t429;
t298 = t381 * t303;
t414 = qJD(5) * t451;
t447 = pkin(5) * t348 + qJ(6) * t459 + t381 * t414 + t298 + (qJD(6) - t280) * t378;
t444 = t381 * t280 + t378 * t303;
t446 = qJ(6) * t460 - qJD(6) * t381 + t378 * t414 + t444;
t399 = t350 * t378 + t420;
t445 = pkin(5) * t399 + t442;
t443 = -t378 * t294 - t330 * t432;
t440 = t378 * t314 + t320;
t427 = t380 * t454;
t421 = qJD(2) * t455;
t415 = qJD(6) + t474;
t264 = t291 * t375 - t471 * t292;
t278 = t323 * t375 + t417;
t327 = t471 * t360 + t361 * t375;
t413 = t342 * t381;
t411 = t381 * t265 - t273 * t433 + t378 * t284 + t287 * t432;
t410 = t382 * t421;
t409 = t379 * t421;
t369 = -pkin(3) * t471 - pkin(4);
t275 = t319 * t471 - t315;
t406 = t264 * t355 - t328 * t340;
t395 = qJ(6) * t294 - t411;
t248 = -qJD(6) * t330 - t395;
t405 = -t251 * t342 + t248;
t403 = t336 + (-t433 - t460) * t342;
t252 = pkin(5) * t294 + t264;
t351 = t377 * t379 + t382 * t456;
t400 = t377 * t382 - t379 * t456;
t309 = t351 * t471 + t375 * t400;
t289 = -t309 * t378 - t381 * t455;
t401 = -t309 * t381 + t378 * t455;
t398 = t350 * t381 - t355 * t433;
t397 = qJD(2) * t472;
t272 = -qJD(3) * pkin(4) - t275;
t394 = t272 * t342 - t368 * t340;
t391 = t351 * qJD(3);
t390 = -0.2e1 * qJD(3) * t472;
t387 = -t391 - t409;
t359 = -t381 * pkin(5) + t369;
t353 = t451 * t381;
t352 = t451 * t378;
t329 = t330 ^ 2;
t322 = qJD(3) * t400 + t410;
t311 = t381 * t314;
t308 = t351 * t375 - t400 * t471;
t300 = pkin(5) * t458 + t327;
t279 = t322 * t471 + t375 * t387;
t277 = t322 * t375 - t387 * t471;
t269 = -pkin(5) * t460 + t278;
t268 = -qJ(6) * t458 + t440;
t267 = t272 + t415;
t266 = -pkin(5) * t396 - qJ(6) * t457 - t328 * t378 + t311;
t258 = qJD(5) * t401 - t279 * t378 + t381 * t422;
t257 = qJD(5) * t289 + t279 * t381 + t378 * t422;
t1 = [-MDP(3) * t427 - t383 * MDP(4) * t454 + (-t382 * t427 + (-t391 - 0.2e1 * t409) * qJD(3)) * MDP(10) + (t379 * t427 + (-t322 - t410) * qJD(3)) * MDP(11) + (-qJD(3) * t277 + (-t340 * t383 + t345 * t435) * t376) * MDP(12) + (-t279 * qJD(3) + (t348 * t435 - t383 * t393) * t376) * MDP(13) + (t277 * t348 - t279 * t345 + t308 * t393 - t309 * t340) * MDP(14) + (t264 * t308 + t265 * t309 - t275 * t277 + t276 * t279 + (t341 * t435 - t344 * t383) * t376) * MDP(15) + (-t257 * t330 - t258 * t332 - t289 * t402 + t294 * t401) * MDP(25) + (t245 * t289 - t248 * t401 + t251 * t258 + t252 * t308 + t254 * t257 + t267 * t277) * MDP(26) + t479 * (t258 * t342 + t277 * t330 + t289 * t340 + t294 * t308) + t478 * (-t257 * t342 + t277 * t332 + t308 * t402 + t340 * t401); 0.2e1 * t419 * t485 - 0.2e1 * t430 * t484 + MDP(7) * t452 - MDP(8) * t453 + (-pkin(8) * t452 + t379 * t390) * MDP(10) + (pkin(8) * t453 + t382 * t390) * MDP(11) + (-t345 * t425 + t340 * t371 + t341 * t347 - t344 * t396 + (t345 * t476 - t442) * qJD(3)) * MDP(12) + (-t348 * t425 + t341 * t350 + t344 * t355 - t371 * t362 + (t348 * t476 + t364 * t371 - t441) * qJD(3)) * MDP(13) + (t265 * t396 - t275 * t350 - t276 * t347 + t327 * t393 - t345 * t441 + t348 * t442 + t406) * MDP(14) + (t264 * t327 + t265 * t328 - t442 * t275 + t441 * t276 + t487 * t341 + t344 * t371) * MDP(15) + (t332 * t398 + t402 * t457) * MDP(16) + ((-t330 * t381 - t463) * t350 + (-t469 - t294 * t381 + (t330 * t378 - t332 * t381) * qJD(5)) * t355) * MDP(17) + (t332 * t347 + t336 * t355 + t342 * t398 - t396 * t402) * MDP(18) + (t294 * t396 - t330 * t347 - t340 * t458 - t342 * t399) * MDP(19) + (-t340 * t396 + t342 * t347) * MDP(20) + (t311 * t340 - (-t273 * t432 + t282) * t396 + t259 * t347 + t327 * t294 + t272 * t420 + (-t328 * t432 + t481) * t342 + t442 * t330 + ((-qJD(5) * t314 - t302) * t342 - (-qJD(5) * t287 - t265) * t396 + t272 * t350 + t406) * t378) * MDP(21) + (-t440 * t340 + t411 * t396 - t260 * t347 + t327 * t402 + t264 * t457 + (t328 * t433 - t480) * t342 + t442 * t332 + t398 * t272) * MDP(22) + (-t245 * t396 + t251 * t347 + t252 * t458 + t266 * t340 + t267 * t399 + t294 * t300 + t330 * t445 + t342 * t450) * MDP(23) + (t248 * t396 + t252 * t457 - t254 * t347 + t267 * t398 - t268 * t340 + t300 * t402 + t332 * t445 - t342 * t449) * MDP(24) + (-t266 * t402 - t268 * t294 + (-t251 * t381 - t254 * t378) * t350 - t450 * t332 - t449 * t330 + (-t245 * t381 - t248 * t378 + (t251 * t378 - t254 * t381) * qJD(5)) * t355) * MDP(25) + (t245 * t266 + t248 * t268 + t251 * t450 + t252 * t300 + t254 * t449 + t267 * t445) * MDP(26); t385 * t484 + t382 * t397 * MDP(11) + (qJD(3) * t278 - t341 * t348 - t345 * t429 - t264) * MDP(12) + (t280 * qJD(3) + t341 * t345 - t348 * t429 - t265) * MDP(13) + ((t276 - t278) * t348 + (t280 - t275) * t345 + (-t375 * t340 - t393 * t471) * pkin(3)) * MDP(14) + (t275 * t278 - t276 * t280 + (-t264 * t471 + t265 * t375 - t341 * t436) * pkin(3)) * MDP(15) + (t332 * t413 + t469) * MDP(16) + ((t402 - t467) * t381 - t483 + t443) * MDP(17) + (t342 * t413 + t462 - t464) * MDP(18) + (t403 + t466) * MDP(19) - t342 * t348 * MDP(20) + (-t259 * t348 - t264 * t381 - t278 * t330 + t369 * t294 + (-t368 * t432 - t298) * t342 + (t280 * t342 + t394) * t378) * MDP(21) + (t260 * t348 + t264 * t378 - t278 * t332 + t369 * t402 + (t368 * t433 + t444) * t342 + t394 * t381) * MDP(22) + (-t251 * t348 - t252 * t381 - t269 * t330 + t294 * t359 - t340 * t352 - t447 * t342 + (t267 * t345 + (t267 + t474) * qJD(5)) * t378) * MDP(23) + (t267 * t459 + t252 * t378 + t254 * t348 - t269 * t332 + t402 * t359 - t340 * t353 + t446 * t342 + (pkin(5) * t463 + t267 * t381) * qJD(5)) * MDP(24) + (-t294 * t353 + t446 * t330 + t447 * t332 + t352 * t402 - t482 * t378 + t405 * t381) * MDP(25) + (-t245 * t352 + t248 * t353 + t252 * t359 + (pkin(5) * t433 - t269) * t267 - t446 * t254 - t447 * t251) * MDP(26) + (MDP(10) * t397 - t385 * t485) * t379; 0.2e1 * t348 * qJD(3) * MDP(12) + (-t362 + (t364 - t345) * qJD(3)) * MDP(13) + (-t345 ^ 2 - t348 ^ 2) * MDP(14) + (t275 * t348 + t276 * t345 + t344) * MDP(15) + ((-t402 - t467) * t381 + t483 + t443) * MDP(25) + (-t267 * t348 + t405 * t378 + t482 * t381) * MDP(26) + t478 * (-t342 ^ 2 * t381 - t462 - t464) + t479 * (t403 - t466); t332 * t330 * MDP(16) + (-t329 + t477) * MDP(17) + (t402 + t468) * MDP(18) + (-t294 + t465) * MDP(19) + t340 * MDP(20) + (t260 * t342 - t272 * t332 + t388) * MDP(21) + (t259 * t342 + t272 * t330 - t411) * MDP(22) + (0.2e1 * t475 + t470 + (-t267 - t415) * t332 + t386) * MDP(23) + (-pkin(5) * t477 + t253 * t342 + (qJD(6) + t267) * t330 + t395) * MDP(24) + (-pkin(5) * t402 - t330 * t448) * MDP(25) + (t448 * t254 + (-t267 * t332 + t245) * pkin(5)) * MDP(26); (t294 + t465) * MDP(23) + (t402 - t468) * MDP(24) + (-t329 - t477) * MDP(25) + (t251 * t332 + t254 * t330 + t252) * MDP(26);];
tauc = t1;
