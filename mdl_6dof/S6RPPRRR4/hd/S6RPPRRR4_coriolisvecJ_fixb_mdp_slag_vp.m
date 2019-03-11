% Calculate minimal parameter regressor of Coriolis joint torque vector for
% S6RPPRRR4
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d4,d5,d6,theta3]';
% MDP [30x1]
%   Minimal dynamic parameter vector (fixed base model)
%   see S6RPPRRR4_convert_par2_MPV_fixb.m
% 
% Output:
% tauc [6x1]
%   joint torques required to compensate Coriolis and centrifugal load

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 02:27
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauc = S6RPPRRR4_coriolisvecJ_fixb_mdp_slag_vp(qJ, qJD, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(10,1),zeros(30,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPPRRR4_coriolisvecJ_fixb_mdp_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RPPRRR4_coriolisvecJ_fixb_mdp_slag_vp: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RPPRRR4_coriolisvecJ_fixb_mdp_slag_vp: pkin has to be [10x1] (double)');
assert(isreal(MDP) && all(size(MDP) == [30 1]), ...
  'S6RPPRRR4_coriolisvecJ_fixb_mdp_slag_vp: MDP has to be [30x1] (double)'); 

%% Symbolic Calculation
% From coriolisvec_joint_fixb_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 02:26:57
% EndTime: 2019-03-09 02:27:05
% DurationCPUTime: 3.94s
% Computational Cost: add. (2457->390), mult. (5065->552), div. (0->0), fcn. (3278->8), ass. (0->180)
t367 = sin(qJ(5));
t370 = cos(qJ(5));
t437 = qJD(4) * t370;
t368 = sin(qJ(4));
t442 = qJD(1) * t368;
t327 = t367 * t442 + t437;
t369 = cos(qJ(6));
t429 = t367 * qJD(4);
t328 = t370 * t442 - t429;
t366 = sin(qJ(6));
t461 = t328 * t366;
t281 = t369 * t327 + t461;
t371 = cos(qJ(4));
t441 = qJD(1) * t371;
t348 = qJD(5) + t441;
t346 = qJD(6) + t348;
t491 = t281 * t346;
t372 = -pkin(1) - pkin(2);
t347 = qJD(1) * t372 + qJD(2);
t363 = sin(pkin(10));
t364 = cos(pkin(10));
t445 = qJ(2) * qJD(1);
t315 = t363 * t347 + t364 * t445;
t304 = -qJD(1) * pkin(7) + t315;
t360 = t368 * qJD(3);
t284 = t371 * t304 + t360;
t276 = qJD(4) * pkin(8) + t284;
t314 = t347 * t364 - t363 * t445;
t303 = qJD(1) * pkin(3) - t314;
t393 = pkin(4) * t371 + pkin(8) * t368;
t277 = qJD(1) * t393 + t303;
t258 = t276 * t370 + t277 * t367;
t256 = pkin(9) * t327 + t258;
t433 = qJD(6) * t366;
t252 = t256 * t433;
t483 = qJD(3) * t371 - t368 * t304;
t275 = -qJD(4) * pkin(4) - t483;
t266 = -pkin(5) * t327 + t275;
t490 = -t266 * t281 + t252;
t385 = t327 * t366 - t369 * t328;
t438 = qJD(4) * t368;
t405 = MDP(28) * t438;
t489 = (-t281 ^ 2 + t385 ^ 2) * MDP(25) - qJD(1) * t405 - t281 * t385 * MDP(24);
t487 = t346 * t385;
t373 = qJD(4) ^ 2;
t374 = qJD(1) ^ 2;
t486 = t363 * (t373 + t374);
t361 = t368 ^ 2;
t485 = (qJD(1) * t361 - t348 * t371) * t367;
t484 = qJ(2) * MDP(6) + MDP(5);
t435 = qJD(5) * t367;
t414 = t368 * t435;
t436 = qJD(4) * t371;
t482 = t370 * t436 - t414;
t481 = MDP(10) * t368;
t480 = MDP(11) * (-t371 ^ 2 + t361);
t479 = qJD(5) + qJD(6);
t427 = qJD(1) * qJD(4);
t407 = t371 * t427;
t426 = qJD(4) * qJD(5);
t298 = qJD(1) * t414 + (-t407 + t426) * t370;
t428 = qJD(1) * qJD(2);
t409 = t364 * t428;
t268 = qJD(4) * t483 + t371 * t409;
t392 = -pkin(4) * t368 + pkin(8) * t371;
t320 = t363 * qJD(2) + qJD(4) * t392;
t305 = t320 * qJD(1);
t297 = t370 * t305;
t375 = -qJD(5) * t258 - t367 * t268 + t297;
t408 = t368 * t427;
t247 = -pkin(5) * t408 - pkin(9) * t298 + t375;
t412 = t371 * t429;
t434 = qJD(5) * t370;
t413 = t368 * t434;
t377 = t412 + t413;
t299 = qJD(1) * t377 - t367 * t426;
t422 = -t370 * t268 - t277 * t434 - t367 * t305;
t378 = -t276 * t435 - t422;
t248 = pkin(9) * t299 + t378;
t404 = t369 * t247 - t366 * t248;
t478 = -t266 * t385 + t404;
t401 = t298 * t366 - t369 * t299;
t254 = qJD(6) * t385 + t401;
t477 = pkin(8) + pkin(9);
t257 = -t276 * t367 + t370 * t277;
t255 = pkin(9) * t328 + t257;
t251 = pkin(5) * t348 + t255;
t475 = t251 * t369;
t432 = qJD(6) * t369;
t421 = t369 * t298 + t366 * t299 + t327 * t432;
t253 = t328 * t433 + t421;
t474 = t253 * t371;
t473 = t254 * t371;
t472 = t256 * t369;
t332 = t366 * t370 + t367 * t369;
t286 = t479 * t332;
t331 = t366 * t367 - t369 * t370;
t379 = t331 * t371;
t261 = qJD(4) * t379 + t286 * t368;
t471 = t261 * t346;
t459 = t367 * t368;
t423 = t366 * t459;
t457 = t370 * t368;
t262 = -qJD(6) * t423 + (t457 * t479 + t412) * t369 + t482 * t366;
t470 = t262 * t346;
t269 = qJD(4) * t360 + t304 * t436 + t368 * t409;
t469 = t269 * t367;
t468 = t269 * t370;
t467 = t275 * t370;
t466 = t298 * t367;
t465 = t298 * t371;
t464 = t299 * t371;
t463 = t327 * t348;
t462 = t328 * t348;
t460 = t348 * t370;
t458 = t367 * t371;
t456 = t370 * t371;
t333 = t392 * qJD(1);
t455 = t367 * t333 + t370 * t483;
t454 = -qJD(1) * t379 - t331 * t479;
t453 = t332 * t441 + t286;
t323 = -t363 * t458 - t364 * t370;
t415 = t368 * t437;
t452 = -t323 * qJD(5) + t363 * t415 + (t363 * t367 + t364 * t456) * qJD(1);
t324 = t363 * t456 - t364 * t367;
t416 = t368 * t429;
t451 = -t324 * qJD(5) + t363 * t416 - (t363 * t370 - t364 * t458) * qJD(1);
t448 = t364 * qJ(2) + t363 * t372;
t330 = -pkin(7) + t448;
t450 = t370 * t320 + t330 * t416;
t400 = -t363 * qJ(2) + t364 * t372;
t329 = pkin(3) - t400;
t306 = t329 + t393;
t316 = t330 * t456;
t449 = t367 * t306 + t316;
t317 = t332 * t368;
t444 = qJD(1) * t317;
t318 = t369 * t457 - t423;
t443 = qJD(1) * t318;
t440 = qJD(2) * t364;
t431 = t275 * qJD(5);
t425 = 0.2e1 * t428;
t424 = t348 * t456;
t418 = t371 * t440;
t420 = t306 * t434 + t367 * t320 + t370 * t418;
t419 = qJD(5) * t477;
t417 = t363 * t436;
t411 = t367 * t441;
t406 = MDP(21) * t438;
t403 = t330 * t348 + t276;
t402 = t370 * t333 - t367 * t483;
t399 = qJD(6) * t251 + t248;
t398 = 0.2e1 * t407;
t397 = t363 * t425;
t396 = t348 * t413;
t394 = -t284 + (t411 + t435) * pkin(5);
t341 = t477 * t370;
t382 = -pkin(5) * t368 + pkin(9) * t456;
t391 = qJD(1) * t382 + qJD(6) * t341 + t370 * t419 + t402;
t340 = t477 * t367;
t390 = pkin(9) * t411 + qJD(6) * t340 + t367 * t419 + t455;
t389 = -qJD(6) * t323 + t452;
t388 = qJD(6) * t324 - t451;
t245 = t251 * t366 + t472;
t301 = t370 * t306;
t264 = pkin(9) * t457 + t301 + (-t330 * t367 + pkin(5)) * t371;
t265 = pkin(9) * t459 + t449;
t387 = t264 * t366 + t265 * t369;
t386 = t314 * t363 - t315 * t364;
t381 = t370 * t361 * t427 + t348 * t414;
t380 = -t330 * t373 + t397;
t376 = qJD(4) * (-qJD(1) * t329 - t303 - t440);
t356 = -pkin(5) * t370 - pkin(4);
t302 = (-pkin(5) * t367 + t330) * t368;
t270 = -pkin(5) * t377 + t330 * t436 + t368 * t440;
t259 = -pkin(5) * t299 + t269;
t250 = (-t371 * t435 - t415) * t330 + t377 * pkin(9) + t420;
t249 = -t367 * t418 + t382 * qJD(4) + (-t316 + (-pkin(9) * t368 - t306) * t367) * qJD(5) + t450;
t244 = -t256 * t366 + t475;
t1 = [t398 * t481 + (-t298 * t457 + t328 * t482) * MDP(17) + (-t346 - t441) * t405 + (-t348 - t441) * t406 + ((t249 * t369 - t250 * t366) * t346 + t404 * t371 - t270 * t281 + t302 * t254 - t259 * t317 - t266 * t262 + (-t245 * t371 - t346 * t387) * qJD(6) + (-(t264 * t369 - t265 * t366) * qJD(1) - t244) * t438) * MDP(29) + (t252 * t371 + t302 * t253 - t259 * t318 + t266 * t261 + t270 * t385 + (-(-qJD(6) * t265 + t249) * t346 - t247 * t371) * t366 + (-(qJD(6) * t264 + t250) * t346 - t399 * t371) * t369 + (qJD(1) * t387 + t245) * t438) * MDP(30) + 0.2e1 * MDP(8) * t409 + (t368 * t376 + t371 * t380) * MDP(15) + (-t368 * t380 + t371 * t376) * MDP(16) - 0.2e1 * t427 * t480 + ((-t363 * t400 + t364 * t448) * qJD(1) - t386) * qJD(2) * MDP(9) + (-t253 * t318 + t261 * t385) * MDP(24) + (t253 * t317 + t254 * t318 + t261 * t281 + t262 * t385) * MDP(25) + MDP(7) * t397 + (t396 + t464 + (-t327 * t368 - t485) * qJD(4)) * MDP(20) + (t465 + (t328 * t368 - t424) * qJD(4) + t381) * MDP(19) + ((-t327 * t370 - t328 * t367) * t436 + (t466 - t299 * t370 + (t327 * t367 - t328 * t370) * qJD(5)) * t368) * MDP(18) + (-t420 * t348 + (t403 * t435 + (-t328 * t330 - t467) * qJD(4) + t422) * t371 + (-t328 * t440 + t367 * t431 - t468 + t330 * t298 + (qJD(1) * t449 + t330 * t460 + t258) * qJD(4)) * t368) * MDP(23) + ((-t306 * t435 + t450) * t348 + (-qJD(4) * t330 * t327 + t297 - t403 * t434 + (-t275 * qJD(4) - qJD(5) * t277 - t348 * t440 - t268) * t367) * t371 + (-t327 * t440 - t370 * t431 - t469 - t330 * t299 + (-(-t330 * t458 + t301) * qJD(1) - t257) * qJD(4)) * t368) * MDP(22) + (-t473 + t470 + (-t281 - t444) * t438) * MDP(27) + (t474 + t471 + (-t385 + t443) * t438) * MDP(26) + t484 * t425 + (-MDP(12) * t371 + MDP(13) * t368) * t373; t386 * qJD(1) * MDP(9) + (0.2e1 * t364 * t408 - t371 * t486) * MDP(15) + (t364 * t398 + t368 * t486) * MDP(16) + (-t327 * t417 + t451 * t348 + (-t299 * t363 + (-qJD(4) * t323 + t364 * t327) * qJD(1)) * t368) * MDP(22) + (-t328 * t417 + t452 * t348 + (t298 * t363 + (qJD(4) * t324 + t364 * t328) * qJD(1)) * t368) * MDP(23) + (-t281 * t417 + (t366 * t389 - t369 * t388) * t346 + (t363 * t254 + (-(t323 * t369 - t324 * t366) * qJD(4) + t364 * t281) * qJD(1)) * t368) * MDP(29) + (t385 * t417 + (t366 * t388 + t369 * t389) * t346 + (t363 * t253 + ((t323 * t366 + t324 * t369) * qJD(4) - t364 * t385) * qJD(1)) * t368) * MDP(30) + (-t363 * MDP(7) - t364 * MDP(8) - t484) * t374; (-t396 + t464) * MDP(22) + (t381 - t465) * MDP(23) + (-t470 - t473) * MDP(29) + (t471 - t474) * MDP(30) + (-MDP(15) * t368 - MDP(16) * t371) * t373 + (-MDP(23) * t424 + MDP(22) * t485 + (-t327 * MDP(22) - t328 * MDP(23) + (-t281 + t444) * MDP(29) + (t385 + t443) * MDP(30)) * t368) * qJD(4); (qJD(4) * t284 - t269) * MDP(15) + (t303 - t440) * t441 * MDP(16) + (-t328 * t460 + t466) * MDP(17) + ((t298 + t463) * t370 + (t299 + t462) * t367) * MDP(18) + (t348 * t434 + (t424 + (-t328 - t429) * t368) * qJD(1)) * MDP(19) + (-t348 * t435 + (-t348 * t458 + (t327 - t437) * t368) * qJD(1)) * MDP(20) + (pkin(4) * t299 - t468 - t402 * t348 + t284 * t327 + (-pkin(8) * t460 + t275 * t367) * qJD(5) + (t257 * t368 + (pkin(8) * t438 + t275 * t371) * t367) * qJD(1)) * MDP(22) + (-pkin(4) * t298 + t469 + t455 * t348 + t284 * t328 + (pkin(8) * t348 * t367 + t467) * qJD(5) + (t275 * t456 + (pkin(8) * t437 - t258) * t368) * qJD(1)) * MDP(23) + (t253 * t332 + t385 * t454) * MDP(24) + (-t253 * t331 - t254 * t332 + t281 * t454 - t385 * t453) * MDP(25) + (t356 * t254 + t259 * t331 + t453 * t266 - t281 * t394) * MDP(29) + (t356 * t253 + t259 * t332 + t454 * t266 + t385 * t394) * MDP(30) + (t454 * MDP(26) - t453 * MDP(27) + (t366 * t390 - t369 * t391) * MDP(29) + (t366 * t391 + t369 * t390) * MDP(30)) * t346 + (t303 * MDP(15) + t348 * MDP(21) + (-qJD(4) * t332 + t385) * MDP(26) + (qJD(4) * t331 + t281) * MDP(27) + t346 * MDP(28) + (-(-t340 * t369 - t341 * t366) * qJD(4) + t244) * MDP(29) + ((-t340 * t366 + t341 * t369) * qJD(4) - t245) * MDP(30)) * t442 + (-t371 * t481 + t480) * t374; t328 * t327 * MDP(17) + (-t327 ^ 2 + t328 ^ 2) * MDP(18) + (t298 - t463) * MDP(19) + (t299 - t462) * MDP(20) - qJD(1) * t406 + (t258 * t348 + t275 * t328 + t375) * MDP(22) + (t257 * t348 - t275 * t327 - t378) * MDP(23) + (t253 - t491) * MDP(26) + (-t254 + t487) * MDP(27) + (-(-t255 * t366 - t472) * t346 - t245 * qJD(6) + (-t281 * t328 - t346 * t433 - t369 * t408) * pkin(5) + t478) * MDP(29) + ((-t256 * t346 - t247) * t366 + (t255 * t346 - t399) * t369 + (t328 * t385 - t346 * t432 + t366 * t408) * pkin(5) + t490) * MDP(30) + t489; (t421 - t491) * MDP(26) + (-t401 + t487) * MDP(27) + (t245 * t346 + t478) * MDP(29) + (t244 * t346 - t366 * t247 - t369 * t248 + t490) * MDP(30) + (MDP(26) * t461 - MDP(27) * t385 - MDP(29) * t245 - MDP(30) * t475) * qJD(6) + t489;];
tauc  = t1;
