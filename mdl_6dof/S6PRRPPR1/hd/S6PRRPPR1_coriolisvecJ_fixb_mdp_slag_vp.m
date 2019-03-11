% Calculate minimal parameter regressor of Coriolis joint torque vector for
% S6PRRPPR1
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [12x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d2,d3,d6,theta1,theta4,theta5]';
% MDP [24x1]
%   Minimal dynamic parameter vector (fixed base model)
%   see S6PRRPPR1_convert_par2_MPV_fixb.m
% 
% Output:
% tauc [6x1]
%   joint torques required to compensate Coriolis and centrifugal load

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-08 21:02
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauc = S6PRRPPR1_coriolisvecJ_fixb_mdp_slag_vp(qJ, qJD, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(12,1),zeros(24,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRRPPR1_coriolisvecJ_fixb_mdp_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6PRRPPR1_coriolisvecJ_fixb_mdp_slag_vp: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [12 1]), ...
  'S6PRRPPR1_coriolisvecJ_fixb_mdp_slag_vp: pkin has to be [12x1] (double)');
assert(isreal(MDP) && all(size(MDP) == [24 1]), ...
  'S6PRRPPR1_coriolisvecJ_fixb_mdp_slag_vp: MDP has to be [24x1] (double)'); 

%% Symbolic Calculation
% From coriolisvec_joint_fixb_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-08 21:02:20
% EndTime: 2019-03-08 21:02:29
% DurationCPUTime: 4.94s
% Computational Cost: add. (3645->388), mult. (9678->567), div. (0->0), fcn. (7589->12), ass. (0->177)
t409 = sin(qJ(2));
t404 = sin(pkin(6));
t457 = qJD(1) * t404;
t444 = t409 * t457;
t408 = sin(qJ(3));
t453 = qJD(3) * t408;
t497 = pkin(3) * t453 - t444;
t411 = cos(qJ(3));
t486 = cos(pkin(11));
t436 = t486 * t411;
t392 = qJD(2) * t436;
t403 = sin(pkin(11));
t454 = qJD(2) * t408;
t369 = t403 * t454 - t392;
t364 = qJD(6) + t369;
t382 = t403 * t411 + t408 * t486;
t371 = t382 * qJD(3);
t419 = -t403 * t408 + t436;
t374 = t419 * qJD(3);
t496 = pkin(4) * t371 - qJ(5) * t374 - qJD(5) * t382 + t497;
t488 = -qJ(4) - pkin(8);
t438 = qJD(3) * t488;
t365 = qJD(4) * t411 + t408 * t438;
t366 = -qJD(4) * t408 + t411 * t438;
t412 = cos(qJ(2));
t443 = t412 * t457;
t462 = t365 * t486 + t366 * t403 - t419 * t443;
t372 = t382 * qJD(2);
t402 = sin(pkin(12));
t405 = cos(pkin(12));
t356 = qJD(3) * t405 - t372 * t402;
t410 = cos(qJ(6));
t495 = t356 * t410;
t407 = sin(qJ(6));
t383 = t402 * t410 + t405 * t407;
t459 = t364 * t383;
t355 = qJD(3) * t402 + t372 * t405;
t305 = t355 * t410 + t356 * t407;
t494 = MDP(5) * t411;
t493 = MDP(6) * (t408 ^ 2 - t411 ^ 2);
t466 = -t402 * t462 + t405 * t496;
t465 = t402 * t496 + t405 * t462;
t463 = t365 * t403 - t366 * t486 - t382 * t443;
t492 = MDP(10) * t408 + MDP(11) * t411;
t361 = qJD(2) * t371;
t470 = t410 * t405;
t474 = t402 * t407;
t381 = -t470 + t474;
t460 = t364 * t381;
t491 = -t361 * t383 + t364 * t460;
t490 = pkin(9) * t405;
t395 = pkin(3) * t403 + qJ(5);
t489 = pkin(9) + t395;
t487 = qJD(2) * pkin(2);
t385 = qJD(2) * pkin(8) + t444;
t406 = cos(pkin(6));
t456 = qJD(1) * t406;
t393 = t411 * t456;
t431 = qJD(4) + t443;
t315 = (-t385 * t408 + t393) * qJD(3) + (-qJ(4) * t453 + t411 * t431) * qJD(2);
t442 = t408 * t456;
t415 = (-t385 * t411 - t442) * qJD(3) + (-qJ(4) * qJD(3) * t411 - t408 * t431) * qJD(2);
t273 = t315 * t403 - t415 * t486;
t473 = t404 * t409;
t377 = t406 * t408 + t411 * t473;
t423 = t406 * t411 - t408 * t473;
t329 = t377 * t403 - t423 * t486;
t485 = t273 * t329;
t388 = t488 * t408;
t389 = t488 * t411;
t351 = -t388 * t486 - t389 * t403;
t484 = t273 * t351;
t303 = t355 * t407 - t495;
t483 = t303 * t372;
t482 = t305 * t372;
t447 = qJD(2) * qJD(3);
t439 = t408 * t447;
t362 = qJD(3) * t392 - t403 * t439;
t480 = t362 * t402;
t479 = t362 * t405;
t478 = t369 * t402;
t477 = t374 * t402;
t476 = t382 * t402;
t475 = t382 * t405;
t435 = qJ(4) * qJD(2) + t385;
t347 = t411 * t435 + t442;
t337 = t403 * t347;
t472 = t404 * t412;
t413 = qJD(3) ^ 2;
t471 = t408 * t413;
t469 = t411 * t413;
t468 = pkin(5) * t371 - t374 * t490 + t466;
t467 = pkin(9) * t477 - t465;
t274 = t315 * t486 + t403 * t415;
t270 = qJD(3) * qJD(5) + t274;
t455 = qJD(2) * t404;
t441 = t409 * t455;
t367 = pkin(3) * t439 + qJD(1) * t441;
t287 = pkin(4) * t361 - qJ(5) * t362 - qJD(5) * t372 + t367;
t257 = t270 * t405 + t287 * t402;
t346 = -t408 * t435 + t393;
t341 = qJD(3) * pkin(3) + t346;
t437 = t486 * t347;
t295 = t341 * t403 + t437;
t290 = qJD(3) * qJ(5) + t295;
t445 = -pkin(3) * t411 - pkin(2);
t363 = qJD(2) * t445 + qJD(4) - t443;
t312 = pkin(4) * t369 - qJ(5) * t372 + t363;
t264 = t290 * t405 + t312 * t402;
t299 = t346 * t486 - t337;
t324 = pkin(3) * t454 + pkin(4) * t372 + qJ(5) * t369;
t269 = t299 * t405 + t324 * t402;
t464 = pkin(5) * t477 + t463;
t333 = -pkin(4) * t419 - qJ(5) * t382 + t445;
t352 = t388 * t403 - t389 * t486;
t292 = t333 * t402 + t352 * t405;
t450 = qJD(6) * t410;
t461 = t356 * t450 + t362 * t470;
t259 = pkin(9) * t356 + t264;
t452 = qJD(6) * t259;
t451 = qJD(6) * t382;
t440 = t412 * t455;
t256 = -t270 * t402 + t287 * t405;
t263 = -t290 * t402 + t312 * t405;
t268 = -t299 * t402 + t324 * t405;
t291 = t333 * t405 - t352 * t402;
t297 = t346 * t403 + t437;
t255 = -pkin(9) * t480 + t257;
t258 = pkin(5) * t369 - pkin(9) * t355 + t263;
t434 = -qJD(6) * t258 - t255;
t398 = -pkin(3) * t486 - pkin(4);
t433 = -t361 * t381 - t364 * t459;
t294 = t341 * t486 - t337;
t251 = t258 * t410 - t259 * t407;
t252 = t258 * t407 + t259 * t410;
t430 = -t263 * t402 + t264 * t405;
t429 = t273 * t382 + t351 * t362;
t277 = -pkin(5) * t419 - pkin(9) * t475 + t291;
t279 = -pkin(9) * t476 + t292;
t428 = t277 * t410 - t279 * t407;
t427 = t277 * t407 + t279 * t410;
t330 = t377 * t486 + t403 * t423;
t313 = -t330 * t402 - t405 * t472;
t314 = t330 * t405 - t402 * t472;
t426 = t313 * t410 - t314 * t407;
t425 = t313 * t407 + t314 * t410;
t424 = -qJD(6) * t355 - t480;
t379 = t489 * t405;
t421 = pkin(5) * t372 + qJD(5) * t402 + qJD(6) * t379 + t369 * t490 + t268;
t378 = t489 * t402;
t420 = pkin(9) * t478 - qJD(5) * t405 + qJD(6) * t378 + t269;
t289 = -qJD(3) * pkin(4) + qJD(5) - t294;
t418 = t289 * t374 + t429;
t417 = -0.2e1 * qJD(3) * t487;
t416 = -t361 * t395 + t362 * t398 + (-qJD(5) + t289) * t369;
t272 = qJD(6) * t305 + t362 * t383;
t414 = qJD(2) ^ 2;
t387 = -pkin(5) * t405 + t398;
t368 = t369 ^ 2;
t345 = -qJD(3) * t377 - t408 * t440;
t344 = qJD(3) * t423 + t411 * t440;
t328 = t381 * t382;
t327 = t383 * t382;
t321 = pkin(5) * t476 + t351;
t298 = t344 * t486 + t345 * t403;
t296 = t344 * t403 - t345 * t486;
t284 = t298 * t405 + t402 * t441;
t283 = -t298 * t402 + t405 * t441;
t282 = -pkin(5) * t478 + t297;
t281 = t374 * t383 + t450 * t475 - t451 * t474;
t280 = -t374 * t381 - t383 * t451;
t278 = -pkin(5) * t356 + t289;
t271 = t407 * t424 + t461;
t265 = pkin(5) * t480 + t273;
t254 = pkin(5) * t361 - pkin(9) * t479 + t256;
t253 = t410 * t254;
t1 = [(t296 * t372 - t298 * t369 + t329 * t362 - t330 * t361) * MDP(12) + (t274 * t330 - t294 * t296 + t295 * t298 + t485) * MDP(13) + (t283 * t369 - t296 * t356 + t313 * t361 + t329 * t480) * MDP(14) + (-t284 * t369 + t296 * t355 - t314 * t361 + t329 * t479) * MDP(15) + (-t283 * t355 + t284 * t356 + (-t313 * t405 - t314 * t402) * t362) * MDP(16) + (t256 * t313 + t257 * t314 + t263 * t283 + t264 * t284 + t289 * t296 + t485) * MDP(17) + ((-qJD(6) * t425 + t283 * t410 - t284 * t407) * t364 + t426 * t361 + t296 * t303 + t329 * t272) * MDP(23) + (-(qJD(6) * t426 + t283 * t407 + t284 * t410) * t364 - t425 * t361 + t296 * t305 + t329 * t271) * MDP(24) + (MDP(10) * t345 - MDP(11) * t344) * qJD(3) + (-t367 * t412 * MDP(13) + (t363 * t409 * MDP(13) - qJD(3) * t412 * t492) * qJD(2) + (-MDP(4) * t412 + (-MDP(10) * t411 + MDP(11) * t408 - MDP(3)) * t409) * t414) * t404; 0.2e1 * t439 * t494 - 0.2e1 * t447 * t493 + MDP(7) * t469 - MDP(8) * t471 + (-pkin(8) * t469 + t408 * t417) * MDP(10) + (pkin(8) * t471 + t411 * t417) * MDP(11) + (t274 * t419 - t294 * t374 - t295 * t371 - t352 * t361 - t369 * t462 + t372 * t463 + t429) * MDP(12) + (t274 * t352 - t463 * t294 + t462 * t295 + t363 * t497 + t367 * t445 + t484) * MDP(13) + (-t256 * t419 + t263 * t371 + t291 * t361 - t356 * t463 + t369 * t466 + t402 * t418) * MDP(14) + (t257 * t419 - t264 * t371 - t292 * t361 + t355 * t463 - t369 * t465 + t405 * t418) * MDP(15) + (-t466 * t355 + t465 * t356 + (-t256 * t382 - t263 * t374 - t291 * t362) * t405 + (-t257 * t382 - t264 * t374 - t292 * t362) * t402) * MDP(16) + (t256 * t291 + t257 * t292 + t263 * t466 + t264 * t465 + t289 * t463 + t484) * MDP(17) + (-t271 * t328 + t280 * t305) * MDP(18) + (-t271 * t327 + t272 * t328 - t280 * t303 - t281 * t305) * MDP(19) + (-t271 * t419 + t280 * t364 + t305 * t371 - t328 * t361) * MDP(20) + (t272 * t419 - t281 * t364 - t303 * t371 - t327 * t361) * MDP(21) + (-t361 * t419 + t364 * t371) * MDP(22) + (t428 * t361 - (-t255 * t407 + t253) * t419 + t251 * t371 + t321 * t272 + t265 * t327 + t278 * t281 + (t407 * t467 + t410 * t468) * t364 + t464 * t303 + (t252 * t419 - t364 * t427) * qJD(6)) * MDP(23) + (-t427 * t361 + (t254 * t407 + t255 * t410) * t419 - t252 * t371 + t321 * t271 - t265 * t328 + t278 * t280 + (-t407 * t468 + t410 * t467) * t364 + t464 * t305 + (t251 * t419 - t364 * t428) * qJD(6)) * MDP(24); ((t295 - t297) * t372 + (-t294 + t299) * t369 + (-t361 * t403 - t362 * t486) * pkin(3)) * MDP(12) + (t294 * t297 - t295 * t299 + (-t273 * t486 + t274 * t403 - t363 * t454) * pkin(3)) * MDP(13) + (-t263 * t372 - t268 * t369 - t273 * t405 + t297 * t356 + t402 * t416) * MDP(14) + (t264 * t372 + t269 * t369 + t273 * t402 - t297 * t355 + t405 * t416) * MDP(15) + (t268 * t355 - t269 * t356 + (qJD(5) * t356 - t263 * t369 + t257) * t405 + (qJD(5) * t355 - t264 * t369 - t256) * t402) * MDP(16) + (-t263 * t268 - t264 * t269 + t273 * t398 - t289 * t297 + (-t256 * t402 + t257 * t405) * t395 + t430 * qJD(5)) * MDP(17) + (t271 * t383 - t305 * t460) * MDP(18) + (-t271 * t381 - t272 * t383 + t303 * t460 - t305 * t459) * MDP(19) + (-t482 - t491) * MDP(20) + (t433 + t483) * MDP(21) - t364 * t372 * MDP(22) + ((-t378 * t410 - t379 * t407) * t361 + t387 * t272 + t265 * t381 - t251 * t372 - t282 * t303 + (t407 * t420 - t410 * t421) * t364 + t459 * t278) * MDP(23) + (-(-t378 * t407 + t379 * t410) * t361 + t387 * t271 + t265 * t383 + t252 * t372 - t282 * t305 + (t407 * t421 + t410 * t420) * t364 - t460 * t278) * MDP(24) + t492 * qJD(2) * t487 + (-t408 * t494 + t493) * t414; (-t372 ^ 2 - t368) * MDP(12) + (t294 * t372 + t367) * MDP(13) + (t356 * t372 + t361 * t405) * MDP(14) + (-t355 * t372 - t361 * t402 - t368 * t405) * MDP(15) + (t256 * t405 + t257 * t402 - t289 * t372) * MDP(17) + (t433 - t483) * MDP(23) + (-t482 + t491) * MDP(24) + (-t402 ^ 2 - t405 ^ 2) * MDP(16) * t362 + (t295 * MDP(13) + (t355 * t402 + t356 * t405) * MDP(16) + t430 * MDP(17) - MDP(14) * t478) * t369; (t355 * t369 + t480) * MDP(14) + (t356 * t369 + t479) * MDP(15) + (-t355 ^ 2 - t356 ^ 2) * MDP(16) + (t263 * t355 - t264 * t356 + t273) * MDP(17) + (t305 * t364 + t272) * MDP(23) + (t364 * t495 + (-t355 * t364 + t424) * t407 + t461) * MDP(24); -t303 ^ 2 * MDP(19) + (t303 * t364 + t461) * MDP(20) + t361 * MDP(22) + (t252 * t364 + t253) * MDP(23) + (t251 * t364 + t278 * t303) * MDP(24) + (MDP(18) * t303 + MDP(19) * t305 + MDP(21) * t364 - MDP(23) * t278) * t305 + (MDP(21) * t424 - MDP(23) * t452 + MDP(24) * t434) * t410 + (t424 * MDP(20) + (-qJD(6) * t356 - t479) * MDP(21) + t434 * MDP(23) + (-t254 + t452) * MDP(24)) * t407;];
tauc  = t1;
