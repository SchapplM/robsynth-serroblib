% Calculate minimal parameter regressor of Coriolis joint torque vector for
% S6PRPRRR3
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [12x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d2,d4,d5,d6,theta1,theta3]';
% MDP [29x1]
%   Minimal dynamic parameter vector (fixed base model)
%   see S6PRPRRR3_convert_par2_MPV_fixb.m
% 
% Output:
% tauc [6x1]
%   joint torques required to compensate Coriolis and centrifugal load

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-08 20:34
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauc = S6PRPRRR3_coriolisvecJ_fixb_mdp_slag_vp(qJ, qJD, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(12,1),zeros(29,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRPRRR3_coriolisvecJ_fixb_mdp_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6PRPRRR3_coriolisvecJ_fixb_mdp_slag_vp: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [12 1]), ...
  'S6PRPRRR3_coriolisvecJ_fixb_mdp_slag_vp: pkin has to be [12x1] (double)');
assert(isreal(MDP) && all(size(MDP) == [29 1]), ...
  'S6PRPRRR3_coriolisvecJ_fixb_mdp_slag_vp: MDP has to be [29x1] (double)'); 

%% Symbolic Calculation
% From coriolisvec_joint_fixb_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-08 20:34:23
% EndTime: 2019-03-08 20:34:31
% DurationCPUTime: 4.55s
% Computational Cost: add. (4016->355), mult. (10546->490), div. (0->0), fcn. (8853->12), ass. (0->168)
t382 = cos(qJ(6));
t431 = qJD(6) * t382;
t374 = sin(pkin(12));
t380 = sin(qJ(4));
t437 = qJD(2) * t380;
t420 = t374 * t437;
t376 = cos(pkin(12));
t384 = cos(qJ(4));
t435 = qJD(2) * t384;
t422 = t376 * t435;
t346 = -t420 + t422;
t354 = t374 * t384 + t376 * t380;
t347 = qJD(2) * t354;
t379 = sin(qJ(5));
t383 = cos(qJ(5));
t317 = -t346 * t383 + t347 * t379;
t489 = t317 * t382;
t496 = t431 + t489;
t375 = sin(pkin(6));
t385 = cos(qJ(2));
t450 = t375 * t385;
t391 = t354 * t450;
t469 = pkin(8) + qJ(3);
t357 = t469 * t374;
t358 = t469 * t376;
t397 = t357 * t380 - t358 * t384;
t495 = qJD(1) * t391 - qJD(3) * t354 + qJD(4) * t397;
t449 = t384 * t376;
t353 = t374 * t380 - t449;
t390 = t353 * t450;
t452 = t357 * t384;
t494 = -qJD(1) * t390 + (qJD(3) * t374 + qJD(4) * t358) * t380 - qJD(3) * t449 + qJD(4) * t452;
t398 = t346 * t379 + t347 * t383;
t378 = sin(qJ(6));
t432 = qJD(6) * t378;
t360 = qJD(4) * t422;
t340 = -qJD(4) * t420 + t360;
t349 = t354 * qJD(4);
t341 = qJD(2) * t349;
t433 = qJD(5) * t383;
t434 = qJD(5) * t379;
t291 = t340 * t383 - t341 * t379 + t346 * t433 - t347 * t434;
t373 = qJD(4) + qJD(5);
t443 = t291 * t382 + t373 * t431;
t274 = -t398 * t432 + t443;
t273 = t274 * t382;
t309 = t373 * t378 + t382 * t398;
t464 = t291 * t378;
t275 = qJD(6) * t309 + t464;
t454 = t398 * t378;
t307 = -t373 * t382 + t454;
t493 = -t378 * t275 - t307 * t496 + t273;
t272 = t274 * t378;
t292 = qJD(5) * t398 + t340 * t379 + t341 * t383;
t284 = t378 * t292;
t485 = -qJD(6) - t317;
t444 = -t431 * t485 + t284;
t456 = t317 * t373;
t458 = t398 * t373;
t460 = t309 * t398;
t492 = (-t292 + t458) * MDP(19) - t317 ^ 2 * MDP(17) + (MDP(16) * t317 + MDP(17) * t398 + MDP(27) * t485) * t398 + (t291 + t456) * MDP(18) + (t309 * t496 + t272) * MDP(23) + (-t485 * t489 + t444 - t460) * MDP(25);
t381 = sin(qJ(2));
t439 = qJD(1) * t375;
t424 = t381 * t439;
t356 = qJD(2) * qJ(3) + t424;
t377 = cos(pkin(6));
t438 = qJD(1) * t377;
t363 = t376 * t438;
t468 = pkin(8) * qJD(2);
t326 = t363 + (-t356 - t468) * t374;
t331 = t356 * t376 + t374 * t438;
t327 = t376 * t468 + t331;
t478 = t326 * t384 - t327 * t380;
t293 = -pkin(9) * t347 + t478;
t290 = qJD(4) * pkin(4) + t293;
t400 = -t326 * t380 - t327 * t384;
t294 = pkin(9) * t346 - t400;
t463 = t294 * t379;
t266 = t290 * t383 - t463;
t264 = -pkin(5) * t373 - t266;
t491 = t264 * t317;
t368 = -pkin(3) * t376 - pkin(2);
t423 = t385 * t439;
t405 = qJD(3) - t423;
t342 = qJD(2) * t368 + t405;
t321 = -pkin(4) * t346 + t342;
t490 = t317 * t321;
t413 = t485 * t378;
t488 = -pkin(9) * t349 - t494;
t348 = t353 * qJD(4);
t487 = -pkin(9) * t348 - t495;
t296 = pkin(5) * t398 + pkin(10) * t317;
t461 = t307 * t398;
t352 = (qJD(3) + t423) * qJD(2);
t481 = t353 * t352;
t441 = t374 ^ 2 + t376 ^ 2;
t479 = t441 * MDP(7);
t286 = t382 * t292;
t477 = -t432 * t485 - t286;
t394 = pkin(4) * t349 - t424;
t462 = t294 * t383;
t267 = t290 * t379 + t462;
t277 = -pkin(9) * t341 + qJD(4) * t478 - t481;
t392 = t354 * t352;
t278 = -pkin(9) * t340 + qJD(4) * t400 - t392;
t416 = t277 * t379 - t278 * t383;
t253 = qJD(5) * t267 + t416;
t265 = pkin(10) * t373 + t267;
t279 = pkin(5) * t317 - pkin(10) * t398 + t321;
t403 = t265 * t378 - t279 * t382;
t476 = -t253 * t382 + t264 * t432 + t398 * t403;
t255 = t265 * t382 + t279 * t378;
t475 = t253 * t378 + t255 * t398 + t264 * t431;
t474 = -t321 * t398 - t416;
t415 = t278 * t379 - t294 * t434;
t252 = (qJD(5) * t290 + t277) * t383 + t415;
t311 = -pkin(9) * t354 - t358 * t380 - t452;
t312 = -pkin(9) * t353 - t397;
t282 = t311 * t379 + t312 * t383;
t322 = t353 * t383 + t354 * t379;
t323 = -t353 * t379 + t354 * t383;
t335 = pkin(4) * t353 + t368;
t283 = pkin(5) * t322 - pkin(10) * t323 + t335;
t297 = -qJD(5) * t322 - t348 * t383 - t349 * t379;
t402 = t311 * t383 - t312 * t379;
t447 = -t402 * qJD(5) + t379 * t487 - t383 * t488;
t473 = -(qJD(6) * t279 + t252) * t322 + t253 * t323 + t264 * t297 - (-qJD(6) * t283 + t447) * t485 - t282 * t292;
t471 = pkin(4) * t347;
t467 = qJD(2) * pkin(2);
t466 = t264 * t323;
t465 = t283 * t292;
t459 = t309 * t378;
t451 = t375 * t381;
t446 = t282 * qJD(5) + t379 * t488 + t383 * t487;
t436 = qJD(2) * t381;
t430 = qJD(6) * t385;
t419 = -pkin(4) * t373 - t290;
t359 = qJD(2) * t424;
t325 = pkin(4) * t341 + t359;
t418 = t441 * t352;
t369 = pkin(4) * t379 + pkin(10);
t409 = qJD(6) * t369 + t296 + t471;
t268 = t293 * t379 + t462;
t408 = pkin(4) * t434 - t268;
t269 = t293 * t383 - t463;
t407 = -pkin(4) * t433 + t269;
t298 = qJD(5) * t323 - t348 * t379 + t349 * t383;
t406 = pkin(5) * t298 - pkin(10) * t297 + t394;
t404 = -t292 * t369 + t491;
t343 = -t374 * t451 + t376 * t377;
t344 = t374 * t377 + t376 * t451;
t314 = t343 * t384 - t344 * t380;
t315 = t343 * t380 + t344 * t384;
t401 = t314 * t383 - t315 * t379;
t289 = t314 * t379 + t315 * t383;
t399 = (-t356 * t374 + t363) * t374 - t331 * t376;
t395 = t317 * t413 - t477;
t393 = t297 * t382 - t323 * t432;
t386 = qJD(2) ^ 2;
t370 = -pkin(4) * t383 - pkin(5);
t355 = t405 - t467;
t300 = -qJD(2) * t391 - qJD(4) * t315;
t299 = -qJD(2) * t390 + qJD(4) * t314;
t261 = pkin(5) * t292 - pkin(10) * t291 + t325;
t260 = t382 * t261;
t257 = qJD(5) * t289 + t299 * t379 - t300 * t383;
t256 = qJD(5) * t401 + t299 * t383 + t300 * t379;
t1 = [(-(-t256 * t378 - t289 * t431) * t485 - t289 * t284 + t257 * t307 - t401 * t275) * MDP(28) + ((t256 * t382 - t289 * t432) * t485 - t289 * t286 + t257 * t309 - t401 * t274) * MDP(29) + (-MDP(21) * t257 - MDP(22) * t256) * t373 + (-t343 * t374 + t344 * t376) * MDP(8) * t352 + (MDP(14) * t300 - MDP(15) * t299) * qJD(4) + ((-t341 * t385 - t346 * t436) * MDP(14) + (-t340 * t385 + t347 * t436) * MDP(15) + (-t292 * t385 + t317 * t436) * MDP(21) + (-t291 * t385 + t398 * t436) * MDP(22) + (-(t378 * t430 + t382 * t436) * t485 - t385 * t286) * MDP(28) + ((t378 * t436 - t382 * t430) * t485 + t385 * t284) * MDP(29) + (t355 * t381 + (-t399 - t424) * t385) * qJD(2) * MDP(8) + ((-MDP(4) + t479) * t385 + (-MDP(5) * t376 + MDP(6) * t374 - MDP(3)) * t381) * t386) * t375; (qJD(2) * t405 * t441 + t418) * MDP(7) + (-t399 * qJD(3) + qJ(3) * t418 + (t399 * t385 + (-t355 - t467) * t381) * t439) * MDP(8) + (t340 * t354 - t347 * t348) * MDP(9) + (-t340 * t353 - t341 * t354 - t346 * t348 - t347 * t349) * MDP(10) + (t368 * t341 + t342 * t349 + (qJD(2) * t353 + t346) * t424) * MDP(14) + (t368 * t340 - t342 * t348) * MDP(15) + (t291 * t323 + t297 * t398) * MDP(16) + (-t291 * t322 - t292 * t323 - t297 * t317 - t298 * t398) * MDP(17) + (t292 * t335 + t298 * t321 + t317 * t394 + t322 * t325) * MDP(21) + (t291 * t335 + t297 * t321 + t323 * t325 + t394 * t398) * MDP(22) + (t273 * t323 + t309 * t393) * MDP(23) + ((-t307 * t382 - t459) * t297 + (-t272 - t275 * t382 + (t307 * t378 - t309 * t382) * qJD(6)) * t323) * MDP(24) + (t274 * t322 + t286 * t323 + t298 * t309 - t393 * t485) * MDP(25) + (-t323 * t284 - t275 * t322 - t298 * t307 - (-t297 * t378 - t323 * t431) * t485) * MDP(26) + (t292 * t322 - t298 * t485) * MDP(27) + (-t403 * t298 + t260 * t322 - t402 * t275 + t446 * t307 + (t465 - t406 * t485 + (-t265 * t322 + t282 * t485 + t466) * qJD(6)) * t382 + t473 * t378) * MDP(28) + (-t255 * t298 - t402 * t274 + t446 * t309 + (-t465 - (-qJD(6) * t265 + t261) * t322 - qJD(6) * t466 - (qJD(6) * t282 - t406) * t485) * t378 + t473 * t382) * MDP(29) + (MDP(18) * t297 - MDP(19) * t298 - MDP(21) * t446 + MDP(22) * t447) * t373 + (-t348 * MDP(11) - t349 * MDP(12) + MDP(14) * t495 + MDP(15) * t494) * qJD(4); (qJD(2) * t399 + t359) * MDP(8) + t360 * MDP(15) + (t292 + t458) * MDP(21) + (t291 - t456) * MDP(22) + (t395 - t461) * MDP(28) + (-t382 * t485 ^ 2 - t284 - t460) * MDP(29) - t386 * t479 + ((t374 * t435 + t376 * t437 + t347) * MDP(14) + (t346 - t420) * MDP(15)) * qJD(4); -t347 * t346 * MDP(9) + (-t346 ^ 2 + t347 ^ 2) * MDP(10) + (t360 + (-t346 - t420) * qJD(4)) * MDP(11) + (-t342 * t347 - t392) * MDP(14) + (-t342 * t346 + t481) * MDP(15) + (-t317 * t471 + t268 * t373 + (t379 * t419 - t462) * qJD(5) + t474) * MDP(21) + (-t398 * t471 + t269 * t373 + t490 + (qJD(5) * t419 - t277) * t383 - t415) * MDP(22) + (t459 * t485 + t493) * MDP(24) + (t395 + t461) * MDP(26) + (t370 * t275 + t404 * t378 + t408 * t307 - (t378 * t407 - t382 * t409) * t485 + t476) * MDP(28) + (t370 * t274 + t404 * t382 + t408 * t309 - (t378 * t409 + t382 * t407) * t485 + t475) * MDP(29) + t492; ((-qJD(5) + t373) * t267 + t474) * MDP(21) + (t266 * t373 - t252 + t490) * MDP(22) + (t309 * t413 + t493) * MDP(24) + (-t413 * t485 + t286 + t461) * MDP(26) + (-pkin(5) * t275 + (-t266 * t378 + t296 * t382) * t485 - t267 * t307 + t378 * t491 - t444 * pkin(10) + t476) * MDP(28) + (-pkin(5) * t274 - (t266 * t382 + t296 * t378) * t485 - t267 * t309 + t264 * t489 + t477 * pkin(10) + t475) * MDP(29) + t492; t309 * t307 * MDP(23) + (-t307 ^ 2 + t309 ^ 2) * MDP(24) + (-t307 * t485 + t443) * MDP(25) + (-t309 * t485 - t464) * MDP(26) + t292 * MDP(27) + (-t252 * t378 - t255 * t485 - t264 * t309 + t260) * MDP(28) + (-t252 * t382 - t261 * t378 + t264 * t307 + t403 * t485) * MDP(29) + (-MDP(25) * t454 - MDP(26) * t309 - MDP(28) * t255 + MDP(29) * t403) * qJD(6);];
tauc  = t1;
