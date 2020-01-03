% Calculate minimal parameter regressor of Coriolis joint torque vector for
% S5RRRRR6
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d2,d3,d4,d5]';
% MDP [27x1]
%   Minimal dynamic parameter vector (fixed base model)
%   see S5RRRRR6_convert_par2_MPV_fixb.m
% 
% Output:
% tauc [5x1]
%   joint torques required to compensate Coriolis and centrifugal load

% Quelle: HybrDyn-Toolbox
% Datum: 2020-01-03 12:16
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauc = S5RRRRR6_coriolisvecJ_fixb_mdp_slag_vp(qJ, qJD, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(9,1),zeros(27,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRRRR6_coriolisvecJ_fixb_mdp_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RRRRR6_coriolisvecJ_fixb_mdp_slag_vp: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5RRRRR6_coriolisvecJ_fixb_mdp_slag_vp: pkin has to be [9x1] (double)');
assert(isreal(MDP) && all(size(MDP) == [27 1]), ...
  'S5RRRRR6_coriolisvecJ_fixb_mdp_slag_vp: MDP has to be [27x1] (double)'); 

%% Symbolic Calculation
% From coriolisvec_joint_fixb_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2020-01-03 12:15:24
% EndTime: 2020-01-03 12:15:30
% DurationCPUTime: 2.49s
% Computational Cost: add. (2729->263), mult. (4630->353), div. (0->0), fcn. (3270->8), ass. (0->158)
t382 = cos(qJ(3));
t373 = qJD(1) + qJD(2);
t379 = sin(qJ(2));
t460 = pkin(1) * qJD(1);
t427 = t379 * t460;
t477 = pkin(7) + pkin(8);
t416 = t477 * t373 + t427;
t319 = t416 * t382;
t381 = cos(qJ(4));
t312 = t381 * t319;
t378 = sin(qJ(3));
t318 = t416 * t378;
t313 = qJD(3) * pkin(3) - t318;
t377 = sin(qJ(4));
t401 = -t313 * t377 - t312;
t449 = t381 * t382;
t422 = t373 * t449;
t455 = t377 * t378;
t423 = t373 * t455;
t327 = -t422 + t423;
t463 = pkin(9) * t327;
t265 = -t401 - t463;
t380 = cos(qJ(5));
t320 = t380 * t327;
t344 = t377 * t382 + t378 * t381;
t329 = t344 * t373;
t376 = sin(qJ(5));
t282 = -t329 * t376 - t320;
t366 = -pkin(3) * t382 - pkin(2);
t383 = cos(qJ(2));
t426 = t383 * t460;
t330 = t366 * t373 - t426;
t287 = pkin(4) * t327 + t330;
t434 = qJD(5) * t376;
t479 = t265 * t434 - t287 * t282;
t372 = qJD(3) + qJD(4);
t438 = qJD(3) * t382;
t419 = t373 * t438;
t284 = qJD(4) * t422 - t372 * t423 + t381 * t419;
t403 = qJD(3) * t416;
t459 = pkin(1) * qJD(2);
t424 = qJD(1) * t459;
t407 = t383 * t424;
t294 = -t378 * t403 + t382 * t407;
t295 = -t378 * t407 - t382 * t403;
t412 = -t377 * t294 + t381 * t295;
t390 = t401 * qJD(4) + t412;
t246 = -pkin(9) * t284 + t390;
t429 = -qJD(4) - qJD(5);
t369 = qJD(3) - t429;
t478 = (-t265 * t369 - t246) * t376 + t479;
t301 = t372 * t344;
t285 = t301 * t373;
t436 = qJD(4) * t377;
t411 = t377 * t295 - t319 * t436;
t474 = (qJD(4) * t313 + t294) * t381;
t245 = -pkin(9) * t285 + t411 + t474;
t400 = t327 * t376 - t380 * t329;
t397 = -t376 * t245 + t380 * t246 + t287 * t400;
t247 = -qJD(5) * t320 + t380 * t284 - t376 * t285 - t329 * t434;
t387 = t400 * qJD(5) - t284 * t376 - t380 * t285;
t476 = t282 * MDP(21) * t400 + (-t282 ^ 2 + t400 ^ 2) * MDP(22) + (-t282 * t369 + t247) * MDP(23) + (-t369 * t400 + t387) * MDP(24);
t475 = MDP(7) * t378;
t473 = (t378 ^ 2 - t382 ^ 2) * MDP(8);
t421 = qJD(3) * t477;
t345 = t378 * t421;
t346 = t382 * t421;
t357 = t477 * t378;
t370 = t382 * pkin(8);
t358 = pkin(7) * t382 + t370;
t398 = t357 * t377 - t358 * t381;
t472 = t398 * qJD(4) + t344 * t426 + t377 * t345 - t381 * t346;
t343 = -t449 + t455;
t435 = qJD(4) * t381;
t471 = -t343 * t426 + t381 * t345 + t377 * t346 + t357 * t435 + t358 * t436;
t324 = t329 * pkin(9);
t310 = t377 * t319;
t410 = t381 * t313 - t310;
t264 = -t324 + t410;
t470 = t378 * MDP(12) + t382 * MDP(13);
t469 = qJD(5) - t369;
t467 = pkin(1) * t383;
t466 = pkin(3) * t369;
t465 = pkin(4) * t329;
t300 = t372 * t343;
t464 = pkin(9) * t300;
t462 = pkin(9) * t344;
t363 = pkin(1) * t379 + pkin(7);
t461 = -pkin(8) - t363;
t457 = t330 * t329;
t456 = t373 * t378;
t454 = t377 * t380;
t384 = qJD(3) ^ 2;
t451 = t378 * t384;
t450 = t380 * t265;
t448 = t382 * t384;
t297 = t343 * t380 + t344 * t376;
t256 = -t297 * qJD(5) - t300 * t380 - t301 * t376;
t360 = t379 * t424;
t439 = qJD(3) * t378;
t420 = t373 * t439;
t332 = pkin(3) * t420 + t360;
t271 = pkin(4) * t285 + t332;
t298 = -t343 * t376 + t344 * t380;
t447 = t287 * t256 + t271 * t298;
t257 = t298 * qJD(5) - t300 * t376 + t380 * t301;
t446 = t287 * t257 + t271 * t297;
t445 = -t330 * t300 + t332 * t344;
t444 = t330 * t301 + t332 * t343;
t443 = -t381 * t318 - t310;
t349 = -pkin(2) * t373 - t426;
t442 = t349 * t438 + t378 * t360;
t437 = qJD(3) * t383;
t432 = t382 * MDP(12);
t430 = -qJD(2) + t373;
t428 = pkin(3) * t456;
t425 = t383 * t459;
t367 = pkin(3) * t439;
t418 = -pkin(3) * t372 - t313;
t262 = pkin(4) * t372 + t264;
t417 = -pkin(4) * t369 - t262;
t291 = pkin(4) * t301 + t367;
t415 = qJD(3) * t461;
t409 = t318 * t377 - t312;
t406 = t291 - t427;
t299 = t301 * pkin(9);
t405 = -qJD(5) * (-t357 * t381 - t358 * t377 - t462) + t299 + t471;
t338 = t343 * pkin(9);
t404 = qJD(5) * (-t338 - t398) - t464 - t472;
t402 = -t376 * t262 - t450;
t339 = t461 * t378;
t340 = t363 * t382 + t370;
t399 = -t339 * t377 - t340 * t381;
t326 = pkin(4) * t343 + t366;
t396 = t330 * t327 - t411;
t314 = t378 * t415 + t382 * t425;
t315 = -t378 * t425 + t382 * t415;
t394 = t381 * t314 + t377 * t315 + t339 * t435 - t340 * t436;
t392 = -t367 + t427;
t389 = t399 * qJD(4) - t314 * t377 + t381 * t315;
t386 = t329 * t327 * MDP(14) + (-t327 ^ 2 + t329 ^ 2) * MDP(15) + t476 + (t327 * t372 + t284) * MDP(16);
t385 = -MDP(10) * t451 + (-t247 * t297 + t256 * t282 + t257 * t400 + t298 * t387) * MDP(22) + (t247 * t298 - t256 * t400) * MDP(21) + (-t284 * t343 - t285 * t344 + t300 * t327 - t301 * t329) * MDP(15) + (t284 * t344 - t300 * t329) * MDP(14) - 0.2e1 * t373 * qJD(3) * t473 + 0.2e1 * t419 * t475 + MDP(9) * t448 + (-t300 * MDP(16) - t301 * MDP(17)) * t372 + (t256 * MDP(23) - t257 * MDP(24)) * t369;
t368 = t379 * t459;
t365 = -pkin(2) - t467;
t364 = pkin(3) * t381 + pkin(4);
t352 = t366 - t467;
t347 = t368 + t367;
t333 = t349 * t439;
t317 = t326 - t467;
t304 = t428 + t465;
t286 = t291 + t368;
t275 = -t338 - t399;
t274 = t339 * t381 - t340 * t377 - t462;
t267 = -t324 + t443;
t266 = t409 + t463;
t252 = t389 + t464;
t251 = -t299 + t394;
t1 = [(((-qJD(1) - t373) * MDP(6) - t470 * qJD(3)) * t383 + (-qJD(1) * t432 + (t378 * MDP(13) - MDP(5) - t432) * t373) * t379) * t459 + (-t363 * t448 + t365 * t420 + t333) * MDP(12) + (t363 * t451 + t365 * t419 + t442) * MDP(13) + (-t286 * t282 - t317 * t387 + (-t251 * t376 + t252 * t380 + (-t274 * t376 - t275 * t380) * qJD(5)) * t369 + t446) * MDP(26) + t385 + (t352 * t285 + t347 * t327 + t389 * t372 + t444) * MDP(19) + (-t286 * t400 + t317 * t247 - (t251 * t380 + t252 * t376 + (t274 * t380 - t275 * t376) * qJD(5)) * t369 + t447) * MDP(27) + (t352 * t284 + t347 * t329 - t394 * t372 + t445) * MDP(20) - t360 * MDP(5); t385 + (-t326 * t387 + (t405 * t376 - t404 * t380) * t369 - t406 * t282 + t446) * MDP(26) + (-pkin(2) * t420 - pkin(7) * t448 + t333 + (t430 * t382 * t379 + t378 * t437) * t460) * MDP(12) + (-pkin(2) * t419 + pkin(7) * t451 + (-t379 * t456 + t382 * t437) * t460 + t442) * MDP(13) + (t366 * t285 - t392 * t327 + t472 * t372 + t444) * MDP(19) + (t366 * t284 - t392 * t329 + t471 * t372 + t445) * MDP(20) + (t326 * t247 + (t404 * t376 + t405 * t380) * t369 - t406 * t400 + t447) * MDP(27) + (t373 * t427 - t360) * MDP(5) + t430 * MDP(6) * t426; t386 + (-t329 * t428 + t443 * t372 + (t418 * qJD(4) - t294) * t381 + t396) * MDP(20) + (-t327 * t428 - t457 - t409 * t372 + (t418 * t377 - t312) * qJD(4) + t412) * MDP(19) + (t304 * t282 - (t266 * t380 - t267 * t376) * t369 + (-t376 * t381 - t454) * qJD(4) * t466 + ((-pkin(3) * t454 - t364 * t376) * t369 + t402) * qJD(5) + t397) * MDP(26) + (t304 * t400 + (-t429 * t377 * t466 + t266 * t369 - t246) * t376 + (-qJD(5) * t262 - t245 + (-pkin(3) * t435 - qJD(5) * t364 + t267) * t369) * t380 + t479) * MDP(27) + t470 * (-t349 * t373 - t407) + (-t382 * t475 + t473) * t373 ^ 2; (-t401 * t372 + t390 - t457) * MDP(19) + (t410 * t372 + t396 - t474) * MDP(20) + (t282 * t465 - (-t264 * t376 - t450) * t369 + (t417 * t376 - t450) * qJD(5) + t397) * MDP(26) + (t400 * t465 + (t417 * qJD(5) + t264 * t369 - t245) * t380 + t478) * MDP(27) + t386; (t469 * t402 + t397) * MDP(26) + ((-t469 * t262 - t245) * t380 + t478) * MDP(27) + t476;];
tauc = t1;
