% Calculate minimal parameter regressor of Coriolis joint torque vector for
% S6RPPRRR2
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [11x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d4,d5,d6,theta2,theta3]';
% MDP [29x1]
%   Minimal dynamic parameter vector (fixed base model)
%   see S6RPPRRR2_convert_par2_MPV_fixb.m
% 
% Output:
% tauc [6x1]
%   joint torques required to compensate Coriolis and centrifugal load

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 02:21
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauc = S6RPPRRR2_coriolisvecJ_fixb_mdp_slag_vp(qJ, qJD, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(11,1),zeros(29,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPPRRR2_coriolisvecJ_fixb_mdp_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RPPRRR2_coriolisvecJ_fixb_mdp_slag_vp: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6RPPRRR2_coriolisvecJ_fixb_mdp_slag_vp: pkin has to be [11x1] (double)');
assert(isreal(MDP) && all(size(MDP) == [29 1]), ...
  'S6RPPRRR2_coriolisvecJ_fixb_mdp_slag_vp: MDP has to be [29x1] (double)'); 

%% Symbolic Calculation
% From coriolisvec_joint_fixb_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 02:21:28
% EndTime: 2019-03-09 02:21:37
% DurationCPUTime: 4.15s
% Computational Cost: add. (3473->357), mult. (8478->471), div. (0->0), fcn. (6544->10), ass. (0->163)
t378 = cos(pkin(11));
t385 = cos(qJ(4));
t425 = qJD(1) * t385;
t365 = t378 * t425;
t376 = sin(pkin(11));
t382 = sin(qJ(4));
t426 = qJD(1) * t382;
t414 = t376 * t426;
t345 = t365 - t414;
t461 = qJD(5) + qJD(6);
t473 = t345 - t461;
t366 = sin(pkin(10)) * pkin(1) + qJ(3);
t359 = t366 * qJD(1);
t372 = t378 * qJD(2);
t457 = pkin(7) * qJD(1);
t327 = t372 + (-t359 - t457) * t376;
t338 = t376 * qJD(2) + t378 * t359;
t328 = t378 * t457 + t338;
t275 = t382 * t327 + t385 * t328;
t472 = qJD(4) * t275;
t407 = qJD(1) * (t376 ^ 2 + t378 ^ 2);
t471 = MDP(7) * t407;
t353 = t376 * t385 + t378 * t382;
t346 = t353 * qJD(1);
t381 = sin(qJ(5));
t384 = cos(qJ(5));
t420 = t384 * qJD(4);
t329 = t346 * t381 - t420;
t383 = cos(qJ(6));
t331 = qJD(4) * t381 + t346 * t384;
t380 = sin(qJ(6));
t446 = t331 * t380;
t277 = t383 * t329 + t446;
t343 = qJD(5) - t345;
t339 = qJD(6) + t343;
t470 = t277 * t339;
t396 = t329 * t380 - t383 * t331;
t469 = t339 * t396;
t436 = t380 * t384;
t355 = t381 * t383 + t436;
t429 = t473 * t355;
t468 = t376 * t425 + t378 * t426;
t424 = qJD(5) * t381;
t442 = t345 * t381;
t467 = t424 - t442;
t271 = qJD(4) * pkin(8) + t275;
t358 = -cos(pkin(10)) * pkin(1) - pkin(3) * t378 - pkin(2);
t344 = qJD(1) * t358 + qJD(3);
t289 = -pkin(4) * t345 - pkin(8) * t346 + t344;
t258 = t271 * t384 + t289 * t381;
t250 = -pkin(9) * t329 + t258;
t422 = qJD(6) * t380;
t248 = t250 * t422;
t463 = t327 * t385 - t382 * t328;
t270 = -qJD(4) * pkin(4) - t463;
t262 = pkin(5) * t329 + t270;
t466 = t262 * t277 + t248;
t364 = qJD(4) * t365;
t340 = -qJD(4) * t414 + t364;
t290 = qJD(5) * t420 + t384 * t340 - t346 * t424;
t348 = t353 * qJD(4);
t341 = qJD(1) * t348;
t352 = t376 * t382 - t385 * t378;
t388 = t352 * qJD(3);
t462 = qJD(1) * t388;
t264 = qJD(4) * t463 - t462;
t300 = pkin(4) * t341 - pkin(8) * t340;
t295 = t384 * t300;
t387 = -qJD(5) * t258 - t264 * t381 + t295;
t238 = pkin(5) * t341 - pkin(9) * t290 + t387;
t291 = qJD(5) * t331 + t340 * t381;
t423 = qJD(5) * t384;
t390 = t384 * t264 - t271 * t424 + t289 * t423 + t381 * t300;
t241 = -pkin(9) * t291 + t390;
t409 = t383 * t238 - t380 * t241;
t465 = t262 * t396 + t409;
t464 = t341 * MDP(27) + (-t277 ^ 2 + t396 ^ 2) * MDP(24) - t277 * MDP(23) * t396;
t308 = t355 * t353;
t458 = pkin(7) + t366;
t349 = t458 * t376;
t350 = t458 * t378;
t310 = t349 * t385 + t382 * t350;
t354 = t380 * t381 - t383 * t384;
t430 = t473 * t354;
t460 = -t339 * t430 - t341 * t355;
t408 = t290 * t380 + t383 * t291;
t246 = -qJD(6) * t396 + t408;
t459 = pkin(8) + pkin(9);
t257 = -t271 * t381 + t384 * t289;
t249 = -pkin(9) * t331 + t257;
t247 = pkin(5) * t343 + t249;
t456 = t247 * t383;
t455 = t250 * t383;
t454 = t277 * t346;
t453 = t396 * t346;
t452 = t290 * t381;
t450 = t329 * t343;
t449 = t329 * t346;
t448 = t331 * t343;
t447 = t331 * t346;
t443 = t341 * t381;
t347 = t352 * qJD(4);
t441 = t347 * t381;
t440 = t347 * t384;
t438 = t353 * t381;
t437 = t353 * t384;
t311 = -t349 * t382 + t350 * t385;
t302 = t384 * t311;
t333 = t384 * t341;
t421 = qJD(6) * t383;
t416 = t383 * t290 - t380 * t291 - t329 * t421;
t245 = -t331 * t422 + t416;
t435 = t245 * t352 - t348 * t396;
t411 = t353 * t424;
t255 = -t347 * t436 - t380 * t411 - t422 * t438 + (t437 * t461 - t441) * t383;
t434 = -t255 * t339 - t308 * t341;
t433 = t290 * t352 + t331 * t348;
t313 = pkin(4) * t346 - pkin(8) * t345;
t432 = t381 * t313 + t384 * t463;
t303 = pkin(4) * t352 - pkin(8) * t353 + t358;
t431 = t381 * t303 + t302;
t418 = t341 * t438;
t417 = t353 * t333;
t415 = qJD(5) * t459;
t410 = t353 * t423;
t406 = t343 * t384;
t405 = qJD(6) * t247 + t241;
t265 = qJD(3) * t468 + t472;
t404 = pkin(5) * t467 - t275;
t403 = t339 * t429 - t354 * t341;
t306 = t384 * t313;
t361 = t459 * t384;
t402 = pkin(5) * t346 + qJD(6) * t361 - t463 * t381 + t306 + (-pkin(9) * t345 + t415) * t384;
t360 = t459 * t381;
t401 = -pkin(9) * t442 + qJD(6) * t360 + t381 * t415 + t432;
t400 = -t246 * t352 - t277 * t348;
t240 = t247 * t380 + t455;
t254 = -t308 * t461 + t354 * t347;
t309 = t354 * t353;
t399 = -t254 * t339 + t309 * t341;
t398 = -t291 * t352 - t329 * t348;
t395 = (-t359 * t376 + t372) * t376 - t338 * t378;
t394 = -t343 * t467 + t333;
t393 = t410 - t441;
t392 = t411 + t440;
t391 = -pkin(8) * t341 + t270 * t343;
t282 = -qJD(4) * t310 - t388;
t314 = pkin(4) * t348 + pkin(8) * t347;
t389 = t384 * t282 + t303 * t423 - t311 * t424 + t381 * t314;
t283 = qJD(3) * t353 + qJD(4) * t311;
t370 = -pkin(5) * t384 - pkin(4);
t315 = t341 * t352;
t307 = t384 * t314;
t299 = t384 * t303;
t284 = pkin(5) * t438 + t310;
t261 = pkin(5) * t393 + t283;
t260 = -pkin(9) * t438 + t431;
t259 = pkin(5) * t352 - pkin(9) * t437 - t311 * t381 + t299;
t253 = pkin(5) * t291 + t265;
t243 = -pkin(9) * t393 + t389;
t242 = pkin(9) * t440 + pkin(5) * t348 - t282 * t381 + t307 + (-t302 + (pkin(9) * t353 - t303) * t381) * qJD(5);
t239 = -t250 * t380 + t456;
t1 = [(t340 * t353 - t346 * t347) * MDP(9) + (-t340 * t352 - t341 * t353 - t345 * t347 - t346 * t348) * MDP(10) + (t341 * t358 + t344 * t348) * MDP(14) + (t340 * t358 - t344 * t347) * MDP(15) + (t290 * t437 - t331 * t392) * MDP(16) + (-(-t329 * t384 - t331 * t381) * t347 + (-t452 - t291 * t384 + (t329 * t381 - t331 * t384) * qJD(5)) * t353) * MDP(17) + (-t343 * t392 + t417 + t433) * MDP(18) + (-t343 * t393 + t398 - t418) * MDP(19) + (t343 * t348 + t315) * MDP(20) + ((-t311 * t423 + t307) * t343 + t299 * t341 + (-t271 * t423 + t295) * t352 + t257 * t348 + t283 * t329 + t310 * t291 + t270 * t410 + ((-qJD(5) * t303 - t282) * t343 - t311 * t341 + (-qJD(5) * t289 - t264) * t352 + t265 * t353 - t270 * t347) * t381) * MDP(21) + (-t258 * t348 + t265 * t437 - t270 * t392 + t283 * t331 + t310 * t290 - t341 * t431 - t343 * t389 - t352 * t390) * MDP(22) + (-t245 * t309 - t254 * t396) * MDP(23) + (-t245 * t308 + t246 * t309 - t254 * t277 + t255 * t396) * MDP(24) + (-t399 + t435) * MDP(25) + (t400 + t434) * MDP(26) + (t339 * t348 + t315) * MDP(27) + ((t242 * t383 - t243 * t380) * t339 + (t259 * t383 - t260 * t380) * t341 + t409 * t352 + t239 * t348 + t261 * t277 + t284 * t246 + t253 * t308 + t262 * t255 + ((-t259 * t380 - t260 * t383) * t339 - t240 * t352) * qJD(6)) * MDP(28) + (-t240 * t348 + t284 * t245 + t248 * t352 - t253 * t309 + t262 * t254 - t261 * t396 + (-(-qJD(6) * t260 + t242) * t339 - t259 * t341 - t238 * t352) * t380 + (-(qJD(6) * t259 + t243) * t339 - t260 * t341 - t405 * t352) * t383) * MDP(29) + (0.2e1 * t471 + (t366 * t407 - t395) * MDP(8)) * qJD(3) + (-MDP(11) * t347 - MDP(12) * t348 - MDP(14) * t283 - MDP(15) * t282) * qJD(4); (-t398 - t418) * MDP(21) + (-t417 + t433) * MDP(22) + (-t400 + t434) * MDP(28) + (t399 + t435) * MDP(29) + (-MDP(21) * t393 + t392 * MDP(22)) * t343 + (-MDP(14) * t348 + MDP(15) * t347) * qJD(4); t364 * MDP(15) + (t394 - t449) * MDP(21) + (-t343 ^ 2 * t384 - t443 - t447) * MDP(22) + (t403 - t454) * MDP(28) + (t453 + t460) * MDP(29) + ((t346 + t468) * MDP(14) + (t345 - t414) * MDP(15)) * qJD(4) + (t395 * MDP(8) - t471) * qJD(1); -t345 ^ 2 * MDP(10) + (t364 + (-t345 - t414) * qJD(4)) * MDP(11) + (-t265 + t472) * MDP(14) + (-t344 * t345 + t462) * MDP(15) + (t331 * t406 + t452) * MDP(16) + ((t290 - t450) * t384 + (-t291 - t448) * t381) * MDP(17) + (t343 * t406 + t443 - t447) * MDP(18) + (t394 + t449) * MDP(19) + (-pkin(4) * t291 - t265 * t384 - t275 * t329 + (-pkin(8) * t423 - t306) * t343 + (t343 * t463 + t391) * t381) * MDP(21) + (-pkin(4) * t290 + t265 * t381 - t275 * t331 + (pkin(8) * t424 + t432) * t343 + t391 * t384) * MDP(22) + (t245 * t355 - t396 * t430) * MDP(23) + (-t245 * t354 - t246 * t355 - t277 * t430 - t396 * t429) * MDP(24) + (t453 - t460) * MDP(25) + (t403 + t454) * MDP(26) + ((-t360 * t383 - t361 * t380) * t341 + t370 * t246 + t253 * t354 + (t380 * t401 - t383 * t402) * t339 + t404 * t277 - t429 * t262) * MDP(28) + (-(-t360 * t380 + t361 * t383) * t341 + t370 * t245 + t253 * t355 + (t380 * t402 + t383 * t401) * t339 - t404 * t396 + t430 * t262) * MDP(29) + (MDP(10) * t346 - MDP(14) * t344 - MDP(20) * t343 - MDP(21) * t257 + MDP(22) * t258 - MDP(27) * t339 - MDP(28) * t239 + MDP(29) * t240 - MDP(9) * t345) * t346; t331 * t329 * MDP(16) + (-t329 ^ 2 + t331 ^ 2) * MDP(17) + (t290 + t450) * MDP(18) + (-t291 + t448) * MDP(19) + t341 * MDP(20) + (t258 * t343 - t270 * t331 + t387) * MDP(21) + (t257 * t343 + t270 * t329 - t390) * MDP(22) + (t245 + t470) * MDP(25) + (-t246 - t469) * MDP(26) + (-(-t249 * t380 - t455) * t339 - t240 * qJD(6) + (-t277 * t331 - t339 * t422 + t341 * t383) * pkin(5) + t465) * MDP(28) + ((-t250 * t339 - t238) * t380 + (t249 * t339 - t405) * t383 + (t331 * t396 - t339 * t421 - t341 * t380) * pkin(5) + t466) * MDP(29) + t464; (t416 + t470) * MDP(25) + (-t408 - t469) * MDP(26) + (t240 * t339 + t465) * MDP(28) + (-t380 * t238 + t239 * t339 - t383 * t241 + t466) * MDP(29) + (-MDP(25) * t446 + MDP(26) * t396 - MDP(28) * t240 - MDP(29) * t456) * qJD(6) + t464;];
tauc  = t1;
