% Calculate minimal parameter regressor of Coriolis joint torque vector for
% S6RPRRPP7
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d3,d4]';
% MDP [28x1]
%   Minimal dynamic parameter vector (fixed base model)
%   see S6RPRRPP7_convert_par2_MPV_fixb.m
% 
% Output:
% tauc [6x1]
%   joint torques required to compensate Coriolis and centrifugal load

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 04:52
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauc = S6RPRRPP7_coriolisvecJ_fixb_mdp_slag_vp(qJ, qJD, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(8,1),zeros(28,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRRPP7_coriolisvecJ_fixb_mdp_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RPRRPP7_coriolisvecJ_fixb_mdp_slag_vp: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S6RPRRPP7_coriolisvecJ_fixb_mdp_slag_vp: pkin has to be [8x1] (double)');
assert(isreal(MDP) && all(size(MDP) == [28 1]), ...
  'S6RPRRPP7_coriolisvecJ_fixb_mdp_slag_vp: MDP has to be [28x1] (double)'); 

%% Symbolic Calculation
% From coriolisvec_joint_fixb_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 04:52:28
% EndTime: 2019-03-09 04:52:34
% DurationCPUTime: 4.72s
% Computational Cost: add. (2995->446), mult. (6052->573), div. (0->0), fcn. (3331->4), ass. (0->180)
t364 = sin(qJ(3));
t445 = qJD(1) * t364;
t354 = qJD(4) + t445;
t365 = cos(qJ(4));
t426 = qJD(3) * qJD(4);
t355 = t365 * t426;
t363 = sin(qJ(4));
t366 = cos(qJ(3));
t437 = qJD(4) * t366;
t410 = t363 * t437;
t431 = t365 * qJD(3);
t375 = t431 * t364 + t410;
t291 = qJD(1) * t375 - t355;
t444 = qJD(1) * t366;
t329 = t363 * t444 - t431;
t469 = t329 * t354;
t493 = -t291 - t469;
t368 = -pkin(1) - pkin(7);
t475 = qJ(5) * t365;
t481 = pkin(4) + pkin(5);
t382 = t363 * t481 - t475;
t492 = t368 - t382;
t427 = qJD(1) * qJD(3);
t407 = t366 * t427;
t491 = qJ(5) * t407 + t354 * qJD(5);
t490 = MDP(23) + MDP(26);
t423 = 0.2e1 * qJD(1);
t362 = t366 ^ 2;
t489 = MDP(8) * (t364 ^ 2 - t362);
t488 = qJ(2) * MDP(6) + MDP(5);
t415 = t363 * t445;
t346 = qJD(3) * t415;
t443 = qJD(3) * t363;
t331 = t365 * t444 + t443;
t440 = qJD(4) * t331;
t292 = -t346 + t440;
t487 = t292 * qJ(6) + t329 * qJD(6);
t486 = -t331 * t354 - t292;
t350 = qJD(1) * t368 + qJD(2);
t336 = t364 * t350;
t485 = qJD(5) * t363 + t336;
t484 = 0.2e1 * t491;
t441 = qJD(3) * t366;
t483 = qJ(5) * t441 + t364 * qJD(5);
t337 = pkin(3) * t364 - pkin(8) * t366 + qJ(2);
t312 = t337 * qJD(1);
t314 = qJD(3) * pkin(8) + t336;
t277 = t365 * t312 - t363 * t314;
t429 = qJD(5) - t277;
t424 = MDP(22) - MDP(27);
t476 = qJ(5) * t363;
t482 = -t365 * t481 - t476;
t480 = pkin(8) - qJ(6);
t478 = qJ(5) * t292;
t477 = qJ(5) * t329;
t442 = qJD(3) * t364;
t259 = t292 * pkin(4) + t291 * qJ(5) - t331 * qJD(5) + t350 * t442;
t474 = t259 * t363;
t473 = t259 * t365;
t278 = t363 * t312 + t365 * t314;
t267 = qJ(6) * t329 + t278;
t347 = t354 * qJ(5);
t263 = t267 + t347;
t472 = t263 * t354;
t271 = t347 + t278;
t471 = t271 * t354;
t470 = t291 * t363;
t468 = t329 * t363;
t467 = t329 * t365;
t465 = t331 * t363;
t464 = t331 * t365;
t463 = t350 * t366;
t462 = t354 * t365;
t461 = t363 * t364;
t460 = t363 * t366;
t459 = t364 * t365;
t458 = t364 * t368;
t457 = t365 * t366;
t370 = qJD(1) ^ 2;
t456 = t366 * t370;
t369 = qJD(3) ^ 2;
t455 = t368 * t369;
t344 = t480 * t365;
t396 = pkin(3) * t366 + pkin(8) * t364;
t333 = t396 * qJD(1);
t403 = t333 * t365 - t350 * t460;
t454 = (qJ(6) * t459 - t366 * t481) * qJD(1) - t403 - qJD(4) * t344 + qJD(6) * t363;
t449 = t363 * t333 + t350 * t457;
t280 = qJ(5) * t444 + t449;
t433 = qJD(6) * t365;
t439 = qJD(4) * t363;
t453 = qJ(6) * t415 - t439 * t480 - t280 - t433;
t452 = t354 * t382 - t485;
t394 = pkin(4) * t363 - t475;
t451 = -t354 * t394 + t485;
t448 = t363 * t337 + t365 * t458;
t446 = -t369 - t370;
t438 = qJD(4) * t365;
t436 = qJD(4) * t368;
t434 = qJD(5) * t365;
t315 = -qJD(3) * pkin(3) - t463;
t374 = qJ(5) * t331 - t315;
t274 = pkin(4) * t329 - t374;
t432 = t274 * MDP(24);
t266 = qJ(6) * t331 + t277;
t430 = qJD(5) - t266;
t264 = -t329 * t481 + qJD(6) + t374;
t428 = qJD(6) + t264;
t425 = MDP(21) + MDP(25);
t422 = pkin(8) * t354 * t363;
t421 = pkin(8) * t462;
t420 = pkin(8) * t441;
t328 = qJD(3) * t396 + qJD(2);
t412 = t368 * t441;
t418 = t363 * t328 + t337 * t438 + t365 * t412;
t409 = t364 * t436;
t417 = t337 * t439 + t363 * t412 + t365 * t409;
t411 = t354 * t438;
t414 = t363 * t441;
t416 = -qJD(1) * t462 - t354 * t414 - t364 * t411;
t286 = t364 * qJ(5) + t448;
t413 = t366 * t431;
t408 = t365 * t437;
t406 = -MDP(19) - t425;
t405 = MDP(20) - t490;
t404 = -t315 + t463;
t348 = t363 * t458;
t402 = t337 * t365 - t348;
t401 = t329 + t431;
t400 = -t331 + t443;
t399 = -t292 + t440;
t304 = t328 * qJD(1);
t398 = -t365 * t304 + t312 * t439 + t314 * t438 + t350 * t414;
t397 = t364 * t407;
t395 = pkin(4) * t365 + t476;
t260 = -t354 * t481 + t430;
t393 = t260 * t365 - t263 * t363;
t392 = t260 * t363 + t263 * t365;
t270 = -pkin(4) * t354 + t429;
t391 = t270 * t365 - t271 * t363;
t390 = t270 * t363 + t271 * t365;
t387 = qJD(1) * t362 - t354 * t364;
t386 = t328 * t365 - t417;
t385 = -t368 + t394;
t384 = -t274 * t364 + t420;
t383 = t315 * t364 - t420;
t254 = -pkin(5) * t292 - t259;
t380 = -t254 * t363 - t264 * t438;
t379 = t254 * t365 - t264 * t439;
t378 = qJ(6) * t291 + t398;
t376 = -t363 * t304 - t312 * t438 + t314 * t439 - t350 * t413;
t373 = -t363 * t409 + t418;
t255 = -t376 + t491;
t372 = t277 * t354 + t376;
t352 = t354 ^ 2;
t343 = t480 * t363;
t338 = -pkin(3) - t395;
t335 = t365 * t397;
t327 = t331 ^ 2;
t320 = pkin(3) - t482;
t305 = t329 * t442;
t296 = t385 * t366;
t287 = -pkin(4) * t364 - t402;
t285 = t492 * t366;
t283 = pkin(4) * t331 + t477;
t282 = qJ(6) * t460 + t286;
t281 = -pkin(4) * t444 - t403;
t276 = t348 + (-qJ(6) * t366 - t337) * t365 - t481 * t364;
t275 = -t331 * t481 - t477;
t272 = -t291 + t469;
t269 = (qJD(4) * t395 - t434) * t366 - t385 * t442;
t265 = -pkin(4) * t441 - t386;
t262 = t373 + t483;
t261 = (qJD(4) * t482 + t434) * t366 - t492 * t442;
t258 = -pkin(4) * t407 + t398;
t257 = qJ(6) * t408 + (qJD(6) * t366 + (-qJ(6) * qJD(3) - t436) * t364) * t363 + t418 + t483;
t256 = (qJ(6) * t442 - t328) * t365 + (qJ(6) * t439 - qJD(3) * t481 - t433) * t366 + t417;
t253 = -qJD(6) * t331 - t407 * t481 + t378;
t252 = t255 + t487;
t1 = [-0.2e1 * MDP(7) * t397 + 0.2e1 * t427 * t489 + (-t364 * t455 + (qJ(2) * t441 + qJD(2) * t364) * t423) * MDP(12) + (-t366 * t455 + (-qJ(2) * t442 + qJD(2) * t366) * t423) * MDP(13) + (-t291 * t457 - t331 * t375) * MDP(14) + ((t465 + t467) * t442 + (t470 - t292 * t365 + (-t464 + t468) * qJD(4)) * t366) * MDP(15) + (-t354 * t410 - t291 * t364 + (t331 * t366 + t365 * t387) * qJD(3)) * MDP(16) + (-t354 * t408 - t292 * t364 + (-t329 * t366 - t363 * t387) * qJD(3)) * MDP(17) + (t354 + t445) * MDP(18) * t441 + (t386 * t354 - t398 * t364 + (-t368 * t292 + t315 * t438) * t366 + ((qJD(1) * t402 + t277) * t366 + (t368 * t329 + t363 * t404) * t364) * qJD(3)) * MDP(19) + (-t373 * t354 + t376 * t364 + (t368 * t291 - t315 * t439) * t366 + ((-qJD(1) * t448 - t278) * t366 + (t368 * t331 + t365 * t404) * t364) * qJD(3)) * MDP(20) + (-t265 * t354 + t269 * t329 + t292 * t296 + (-t274 * t443 - t258) * t364 + (t274 * t438 + t474 + (-qJD(1) * t287 - t270) * qJD(3)) * t366) * MDP(21) + (-t262 * t329 + t265 * t331 - t286 * t292 - t287 * t291 - t391 * t442 + (-qJD(4) * t390 - t255 * t363 + t258 * t365) * t366) * MDP(22) + (t262 * t354 - t269 * t331 + t291 * t296 + (t274 * t431 + t255) * t364 + (t274 * t439 - t473 + (qJD(1) * t286 + t271) * qJD(3)) * t366) * MDP(23) + (t255 * t286 + t258 * t287 + t259 * t296 + t262 * t271 + t265 * t270 + t269 * t274) * MDP(24) + (-t256 * t354 - t261 * t329 - t285 * t292 + (t264 * t443 - t253) * t364 + ((-qJD(1) * t276 - t260) * qJD(3) + t380) * t366) * MDP(25) + (t257 * t354 + t261 * t331 - t285 * t291 + (-t264 * t431 + t252) * t364 + ((qJD(1) * t282 + t263) * qJD(3) + t379) * t366) * MDP(26) + (-t256 * t331 + t257 * t329 + t276 * t291 + t282 * t292 + t393 * t442 + (qJD(4) * t392 + t252 * t363 - t253 * t365) * t366) * MDP(27) + (t252 * t282 + t253 * t276 + t254 * t285 + t256 * t260 + t257 * t263 + t261 * t264) * MDP(28) + t488 * qJD(2) * t423 + (-MDP(10) * t366 - MDP(9) * t364) * t369; t305 * MDP(19) - t335 * MDP(20) + (t305 + t416) * MDP(21) + t416 * MDP(25) - t488 * t370 + (t391 * MDP(24) + t393 * MDP(28) + (-t365 * MDP(19) + t363 * t405) * t354 + t424 * (t464 + t468)) * qJD(1) + (t446 * MDP(13) - t259 * MDP(24) + t254 * MDP(28) + t406 * t292 + t405 * t291 + (t390 * MDP(24) + t392 * MDP(28) + (-t363 * MDP(19) - t365 * MDP(20)) * t354 - t424 * (-t465 + t467)) * qJD(3)) * t366 + (t446 * MDP(12) + (t329 * MDP(25) - t264 * MDP(28) + t331 * t405 + t432) * qJD(3) + (-qJD(4) * t354 * MDP(19) + (qJD(4) * t270 + t255) * MDP(24) + (qJD(4) * t260 + t252) * MDP(28) + t424 * t399) * t365 + (t258 * MDP(24) + t253 * MDP(28) - t424 * t291 + (-t271 * MDP(24) - t263 * MDP(28) + t329 * t424 + t354 * t405) * qJD(4) + t406 * t407) * t363) * t364 + t490 * (t354 * t413 + t335); t364 * MDP(7) * t456 - t370 * t489 + (t331 * t462 - t470) * MDP(14) + (t486 * t363 + t493 * t365) * MDP(15) + (t411 + (t354 * t459 + t366 * t400) * qJD(1)) * MDP(16) + (-t354 * t439 + (-t354 * t461 + t366 * t401) * qJD(1)) * MDP(17) - t354 * MDP(18) * t444 + (-pkin(3) * t292 - t403 * t354 - t401 * t336 + (t315 * t363 - t421) * qJD(4) + (-t277 * t366 + t363 * t383) * qJD(1)) * MDP(19) + (pkin(3) * t291 + t449 * t354 + t400 * t336 + (t315 * t365 + t422) * qJD(4) + (t278 * t366 + t365 * t383) * qJD(1)) * MDP(20) + (-t473 + t281 * t354 + t292 * t338 - t451 * t329 + (t274 * t363 - t421) * qJD(4) + (t270 * t366 - t363 * t384) * qJD(1)) * MDP(21) + (t280 * t329 - t281 * t331 + (pkin(8) * t399 + t270 * t354 + t255) * t365 + (t258 - t471 + (qJD(4) * t329 - t291) * pkin(8)) * t363) * MDP(22) + (-t474 - t280 * t354 + t291 * t338 + t451 * t331 + (-t274 * t365 - t422) * qJD(4) + (-t271 * t366 + t365 * t384) * qJD(1)) * MDP(23) + (t259 * t338 - t270 * t281 - t271 * t280 - t451 * t274 + (qJD(4) * t391 + t255 * t365 + t258 * t363) * pkin(8)) * MDP(24) + (-t292 * t320 + t454 * t354 + t452 * t329 + (-t264 * t461 + (-qJD(3) * t343 + t260) * t366) * qJD(1) + t379) * MDP(25) + (-t291 * t320 + t453 * t354 - t452 * t331 + (t264 * t459 + (qJD(3) * t344 - t263) * t366) * qJD(1) - t380) * MDP(26) + (t291 * t343 + t292 * t344 + t454 * t331 + t453 * t329 + (-t260 * t354 - t252) * t365 + (-t253 + t472) * t363) * MDP(27) + (t252 * t344 + t253 * t343 + t254 * t320 - t260 * t454 + t263 * t453 - t264 * t452) * MDP(28) + (MDP(13) * t364 * t370 - MDP(12) * t456) * qJ(2); t272 * MDP(16) + (-t363 * t426 + t346) * MDP(17) + t372 * MDP(20) + (pkin(4) * t291 - t478) * MDP(22) + (-t372 + t484) * MDP(23) + (-pkin(4) * t258 + qJ(5) * t255 - t270 * t278 + t271 * t429 - t274 * t283) * MDP(24) + (t267 * t354 - t378) * MDP(25) + (-t266 * t354 - t376 + t484 + t487) * MDP(26) + (-t291 * t481 + t478) * MDP(27) + (qJ(5) * t252 - t253 * t481 - t260 * t267 + t263 * t430 - t264 * t275) * MDP(28) + (-MDP(17) * t438 + (0.2e1 * pkin(4) * MDP(21) + 0.2e1 * t481 * MDP(25) + MDP(18)) * qJD(3)) * t444 + (t354 * MDP(17) - t315 * MDP(19) - t274 * MDP(21) + (t271 - t278) * MDP(22) + t283 * MDP(23) + t428 * MDP(25) - t275 * MDP(26) + (-t263 + t267) * MDP(27) + MDP(15) * t331) * t331 + (t331 * MDP(14) + t315 * MDP(20) - t283 * MDP(21) + (t270 - t429) * MDP(22) - t274 * MDP(23) + t275 * MDP(25) + t264 * MDP(26) + (-t260 + t430) * MDP(27) - MDP(15) * t329) * t329 + (MDP(19) + MDP(21)) * (t278 * t354 - t398); t272 * MDP(22) + (-t352 - t327) * MDP(23) + (t398 - t471) * MDP(24) - t352 * MDP(26) + (-t355 - t469) * MDP(27) + (t378 - t472) * MDP(28) + (-MDP(26) * t331 - MDP(28) * t428 + t329 * t425 + t432) * t331 + (MDP(27) * t410 + (MDP(27) * t459 + (-pkin(4) * MDP(24) - MDP(28) * t481 - t425) * t366) * qJD(3)) * qJD(1); t486 * MDP(25) + t493 * MDP(26) + (-t329 ^ 2 - t327) * MDP(27) + (t260 * t331 - t263 * t329 + t254) * MDP(28);];
tauc  = t1;
