% Calculate minimal parameter regressor of Coriolis joint torque vector for
% S6RRRPPR3
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d3,d6]';
% MDP [32x1]
%   Minimal dynamic parameter vector (fixed base model)
%   see S6RRRPPR3_convert_par2_MPV_fixb.m
% 
% Output:
% tauc [6x1]
%   joint torques required to compensate Coriolis and centrifugal load

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 15:30
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauc = S6RRRPPR3_coriolisvecJ_fixb_mdp_slag_vp(qJ, qJD, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(9,1),zeros(32,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRPPR3_coriolisvecJ_fixb_mdp_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRRPPR3_coriolisvecJ_fixb_mdp_slag_vp: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S6RRRPPR3_coriolisvecJ_fixb_mdp_slag_vp: pkin has to be [9x1] (double)');
assert(isreal(MDP) && all(size(MDP) == [32 1]), ...
  'S6RRRPPR3_coriolisvecJ_fixb_mdp_slag_vp: MDP has to be [32x1] (double)'); 

%% Symbolic Calculation
% From coriolisvec_joint_fixb_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 15:30:14
% EndTime: 2019-03-09 15:30:23
% DurationCPUTime: 4.13s
% Computational Cost: add. (3430->410), mult. (8258->488), div. (0->0), fcn. (5587->6), ass. (0->172)
t386 = qJD(2) + qJD(3);
t396 = cos(qJ(6));
t397 = cos(qJ(2));
t484 = cos(qJ(3));
t436 = t484 * t397;
t422 = qJD(1) * t436;
t394 = sin(qJ(3));
t395 = sin(qJ(2));
t448 = qJD(1) * t395;
t435 = t394 * t448;
t341 = -t422 + t435;
t393 = sin(qJ(6));
t464 = t341 * t393;
t318 = t396 * t386 + t464;
t358 = t394 * t397 + t395 * t484;
t449 = qJD(1) * t358;
t492 = qJD(6) + t449;
t496 = t318 * t492;
t320 = t341 * t396 - t386 * t393;
t495 = t320 * t492;
t494 = t492 * t396;
t381 = -t397 * pkin(2) - pkin(1);
t485 = -pkin(8) - pkin(7);
t366 = t485 * t397;
t362 = qJD(1) * t366;
t345 = t394 * t362;
t365 = t485 * t395;
t360 = qJD(1) * t365;
t482 = qJD(2) * pkin(2);
t349 = t360 + t482;
t453 = -t484 * t349 - t345;
t493 = qJD(4) + t453;
t441 = qJD(1) * qJD(2);
t491 = -0.2e1 * t441;
t490 = MDP(5) * (t395 ^ 2 - t397 ^ 2);
t348 = t484 * t362;
t314 = t394 * t360 - t348;
t447 = qJD(3) * t394;
t420 = pkin(2) * t447 - t314;
t315 = t360 * t484 + t345;
t434 = qJD(3) * t484;
t454 = pkin(2) * t434 + qJD(4) - t315;
t460 = t394 * t395;
t419 = t386 * t460;
t452 = t386 * t422;
t305 = qJD(1) * t419 - t452;
t458 = t396 * t305;
t488 = -t393 * t492 ^ 2 - t458;
t437 = qJD(2) * t485;
t423 = qJD(1) * t437;
t350 = t395 * t423;
t351 = t397 * t423;
t273 = t349 * t447 + t394 * t350 - t484 * t351 - t362 * t434;
t414 = qJ(5) * t305 + t273;
t262 = -qJD(5) * t449 + t414;
t322 = t394 * t365 - t366 * t484;
t361 = t395 * t437;
t363 = t397 * t437;
t281 = qJD(3) * t322 + t394 * t361 - t363 * t484;
t316 = -qJD(2) * t436 - t397 * t434 + t419;
t265 = t316 * qJ(5) - t358 * qJD(5) + t281;
t364 = t381 * qJD(1);
t295 = t341 * pkin(3) - qJ(4) * t449 + t364;
t418 = qJD(5) - t295;
t486 = -pkin(4) - pkin(9);
t267 = pkin(5) * t449 + t341 * t486 + t418;
t357 = -t436 + t460;
t428 = -t357 * pkin(3) + t358 * qJ(4) - t381;
t272 = pkin(5) * t358 + t357 * t486 + t428;
t310 = t394 * t349 - t348;
t466 = t341 * qJ(5);
t289 = t310 + t466;
t384 = t386 * qJ(4);
t285 = -t289 - t384;
t279 = pkin(5) * t386 - t285;
t317 = t386 * t358;
t383 = t386 * qJD(4);
t411 = -t349 * t434 - t484 * t350 - t394 * t351 - t362 * t447;
t271 = t383 - t411;
t306 = t317 * qJD(1);
t261 = -t306 * qJ(5) - t341 * qJD(5) - t271;
t321 = -t365 * t484 - t394 * t366;
t298 = -t358 * qJ(5) + t321;
t417 = -t261 * t357 + t298 * t305;
t487 = t279 * t317 - (qJD(6) * t272 + t265) * t492 - (qJD(6) * t267 + t262) * t358 + t417;
t340 = t449 ^ 2;
t398 = -pkin(3) - pkin(4);
t483 = pkin(4) * t449;
t481 = qJ(4) * t306;
t388 = -pkin(3) + t486;
t463 = t449 * qJ(5);
t416 = -t463 + t493;
t274 = t386 * t388 + t416;
t255 = t267 * t393 + t274 * t396;
t480 = t255 * t341;
t479 = t272 * t305;
t446 = qJD(6) * t393;
t445 = qJD(6) * t396;
t455 = t396 * t306 - t386 * t445;
t275 = -t341 * t446 + t455;
t478 = t275 * t393;
t477 = t279 * t357;
t280 = t484 * t361 + t394 * t363 + t365 * t434 + t366 * t447;
t476 = t280 * t386;
t475 = t281 * t386;
t378 = pkin(2) * t394 + qJ(4);
t474 = t306 * t378;
t473 = t306 * t393;
t472 = t453 * t386;
t471 = t310 * t386;
t470 = t318 * t341;
t469 = t320 * t341;
t465 = t341 * t386;
t461 = t393 * t305;
t399 = qJD(2) ^ 2;
t459 = t395 * t399;
t457 = t397 * t399;
t400 = qJD(1) ^ 2;
t456 = t397 * t400;
t307 = pkin(3) * t449 + t341 * qJ(4);
t451 = -t463 + t454;
t278 = -pkin(4) * t341 + t418;
t443 = -qJD(5) - t278;
t440 = pkin(2) * t448;
t439 = t395 * t482;
t433 = t395 * t441;
t432 = pkin(1) * t491;
t254 = t267 * t396 - t274 * t393;
t431 = t254 * t341 - t261 * t396;
t429 = t492 * t279;
t269 = -pkin(5) * t341 + t449 * t486 - t307;
t380 = -pkin(2) * t484 - pkin(3);
t377 = -pkin(4) + t380;
t371 = -pkin(9) + t377;
t425 = qJD(6) * t371 + t269 - t440;
t424 = qJD(6) * t388 + t269;
t421 = -t466 + t420;
t297 = t440 + t307;
t415 = t317 * t396 - t357 * t446;
t270 = t317 * pkin(3) + t316 * qJ(4) - t358 * qJD(4) + t439;
t413 = -t295 * t449 - t273;
t412 = -t364 * t449 - t273;
t266 = pkin(2) * t433 + t306 * pkin(3) + t305 * qJ(4) - qJD(4) * t449;
t409 = -t492 * t494 + t461;
t408 = -t386 * t435 + t452;
t407 = -t295 * t341 - t411;
t406 = t364 * t341 + t411;
t405 = t289 * t492 + t388 * t305 - t429;
t259 = -pkin(4) * t306 - t266;
t404 = t443 * t449 + t414;
t403 = t278 * t341 - t261;
t402 = t371 * t305 - t421 * t492 - t429;
t276 = qJD(6) * t320 + t473;
t286 = t408 + t465;
t339 = t341 ^ 2;
t401 = ((-t275 + t496) * t396 + (t276 + t495) * t393) * MDP(27) + (-t320 * t494 - t478) * MDP(26) + (-t470 - t488) * MDP(29) + (t409 + t469) * MDP(28) + t286 * MDP(13) + (-t339 + t340) * MDP(12) + (MDP(11) * t449 + MDP(30) * t492) * t341;
t392 = qJ(4) + pkin(5);
t375 = pkin(5) + t378;
t300 = t384 + t310;
t299 = t357 * qJ(5) + t322;
t296 = -pkin(3) * t386 + t493;
t292 = -pkin(4) * t357 + t428;
t291 = t305 * t358;
t287 = -t307 - t483;
t282 = -t297 - t483;
t277 = t386 * t398 + t416;
t264 = -t317 * qJ(5) - t357 * qJD(5) - t280;
t263 = -pkin(4) * t317 - t270;
t253 = -pkin(5) * t316 + t317 * t486 - t270;
t252 = -pkin(5) * t305 + t306 * t486 - t266;
t251 = t396 * t252;
t1 = [(t259 * t357 + t263 * t341 + t278 * t317 + t292 * t306) * MDP(23) + (-t266 * t428 + t270 * t295 + t271 * t322 + t273 * t321 + t280 * t300 + t281 * t296) * MDP(21) + (t266 * t357 + t270 * t341 + t295 * t317 - t306 * t428 - t475) * MDP(18) + (t259 * t358 + t263 * t449 - t278 * t316 - t292 * t305) * MDP(22) + (-t266 * t358 - t270 * t449 + t295 * t316 - t305 * t428 + t476) * MDP(20) + (-t316 * t449 - t291) * MDP(11) + (-t305 * t381 - t316 * t364 + 0.2e1 * t449 * t439 - t476) * MDP(17) + (-t262 * t358 - t264 * t341 - t265 * t449 + t277 * t316 - t285 * t317 + t299 * t306 + t417) * MDP(24) + (t305 * t357 - t306 * t358 + t316 * t341 - t317 * t449) * MDP(12) + (-t271 * t357 + t273 * t358 - t280 * t341 + t281 * t449 - t296 * t316 - t300 * t317 - t305 * t321 - t306 * t322) * MDP(19) + MDP(6) * t457 + (t259 * t292 - t261 * t299 + t262 * t298 + t263 * t278 + t264 * t285 + t265 * t277) * MDP(25) + t490 * t491 + (-t316 * t492 - t291) * MDP(30) + (t251 * t358 - t254 * t316 - t264 * t318 + t299 * t276 + (t253 * t492 - t479 + (-t274 * t358 - t298 * t492 + t477) * qJD(6)) * t396 + t487 * t393) * MDP(31) + (t255 * t316 - t264 * t320 + t299 * t275 + (-(-qJD(6) * t298 + t253) * t492 + t479 - (-qJD(6) * t274 + t252) * t358 - qJD(6) * t477) * t393 + t487 * t396) * MDP(32) + (t357 * t461 - t276 * t358 + t316 * t318 + (-t317 * t393 - t357 * t445) * t492) * MDP(29) + (t275 * t358 - t316 * t320 - t357 * t458 + t415 * t492) * MDP(28) + (-t475 + t306 * t381 + t317 * t364 + (qJD(1) * t357 + t341) * t439) * MDP(16) + ((-t318 * t396 - t320 * t393) * t317 + (-t478 - t276 * t396 + (t318 * t393 - t320 * t396) * qJD(6)) * t357) * MDP(27) - MDP(7) * t459 + (pkin(7) * t459 + t397 * t432) * MDP(10) + (-pkin(7) * t457 + t395 * t432) * MDP(9) + (-MDP(13) * t316 - MDP(14) * t317 - MDP(22) * t264 + MDP(23) * t265) * t386 + (t275 * t357 * t396 + t320 * t415) * MDP(26) + 0.2e1 * t397 * MDP(4) * t433; (t375 * t276 + t318 * t451 + t393 * t402 - t425 * t494 + t431) * MDP(31) + (-t282 * t449 + t386 * t451 + t403) * MDP(22) + (-t282 * t341 + t386 * t421 + t404) * MDP(23) + (t297 * t449 + t386 * t454 + t383 + t407) * MDP(20) + (-t305 * t380 - t474 + (t300 + t420) * t449 + (t296 - t454) * t341) * MDP(19) + (t305 * t377 + t474 + (t285 - t421) * t449 + (-t277 + t451) * t341) * MDP(24) + (t314 * t386 + (-t341 * t448 - t386 * t447) * pkin(2) + t412) * MDP(16) + (-t297 * t341 - t386 * t420 + t413) * MDP(18) + (-t480 + t375 * t275 + t451 * t320 + (t425 * t492 + t261) * t393 + t402 * t396) * MDP(32) + (t315 * t386 + (-t386 * t434 - t448 * t449) * pkin(2) + t406) * MDP(17) + (t271 * t378 + t273 * t380 - t295 * t297 + t296 * t420 + t300 * t454) * MDP(21) + (-t261 * t378 + t262 * t377 + t277 * t421 - t278 * t282 - t285 * t451) * MDP(25) + t400 * t490 + t401 - t395 * MDP(4) * t456 + (MDP(9) * t395 * t400 + MDP(10) * t456) * pkin(1); (t392 * t276 + t318 * t416 + t393 * t405 - t424 * t494 + t431) * MDP(31) + (-t480 + t392 * t275 + t416 * t320 + (t424 * t492 + t261) * t393 + t405 * t396) * MDP(32) + (t412 + t471) * MDP(16) + (-t307 * t341 + t413 + t471) * MDP(18) + (pkin(3) * t305 - t481 + (t300 - t310) * t449 + (t296 - t493) * t341) * MDP(19) + (t307 * t449 + 0.2e1 * t383 + t407 + t472) * MDP(20) + (-t287 * t449 + t386 * t416 + t403) * MDP(22) + (-t287 * t341 - t289 * t386 + t404) * MDP(23) + (t481 + t305 * t398 + (t285 + t289) * t449 + (-t277 + t416) * t341) * MDP(24) + (-qJ(4) * t261 + t262 * t398 - t277 * t289 - t278 * t287 - t285 * t416) * MDP(25) + (-pkin(3) * t273 + qJ(4) * t271 - t295 * t307 - t296 * t310 + t300 * t493) * MDP(21) + t401 + (t406 - t472) * MDP(17); t286 * MDP(19) - t340 * MDP(20) + t273 * MDP(21) + t305 * MDP(24) + t414 * MDP(25) + (-t445 * t492 + t461) * MDP(31) + (t446 * t492 + t458) * MDP(32) + (-t300 * MDP(21) - t341 * MDP(24) + t285 * MDP(25) - t318 * MDP(31) - t320 * MDP(32) + (-MDP(20) - MDP(22)) * t386) * t386 + (t295 * MDP(21) + t443 * MDP(25) - MDP(22) * t449 + (MDP(18) - MDP(23)) * t341 + (-MDP(31) * t396 + MDP(32) * t393) * t492) * t449; (t408 - t465) * MDP(22) + (t386 * t449 + t306) * MDP(23) + (-t340 - t339) * MDP(24) + (t277 * t449 + t285 * t341 + t259) * MDP(25) + (-t470 + t488) * MDP(31) + (t409 - t469) * MDP(32); t320 * t318 * MDP(26) + (-t318 ^ 2 + t320 ^ 2) * MDP(27) + (t455 + t496) * MDP(28) + (-t473 + t495) * MDP(29) - t305 * MDP(30) + (t255 * t492 - t262 * t393 - t279 * t320 + t251) * MDP(31) + (-t252 * t393 + t254 * t492 - t262 * t396 + t279 * t318) * MDP(32) + (-MDP(28) * t464 - MDP(29) * t320 - MDP(31) * t255 - MDP(32) * t254) * qJD(6);];
tauc  = t1;
