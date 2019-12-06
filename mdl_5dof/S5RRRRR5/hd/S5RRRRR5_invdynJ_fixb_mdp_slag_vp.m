% Calculate vector of inverse dynamics joint torques for
% S5RRRRR5
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% qJDD [5x1]
%   Generalized joint accelerations
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d2,d3,d4,d5]';
% MDP [23x1]
%   Minimal dynamic parameter vector (fixed base model)
%   see S5RRRRR5_convert_par2_MPV_fixb.m
% 
% Output:
% tau [5x1]
%   joint torques of inverse dynamics (contains inertial, gravitational Coriolis and centrifugal forces)

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-05 18:59
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tau = S5RRRRR5_invdynJ_fixb_mdp_slag_vp(qJ, qJD, qJDD, g, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(5,1),zeros(3,1),zeros(9,1),zeros(23,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRRRR5_invdynJ_fixb_mdp_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RRRRR5_invdynJ_fixb_mdp_slag_vp: qJD has to be [5x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [5 1]), ...
  'S5RRRRR5_invdynJ_fixb_mdp_slag_vp: qJDD has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RRRRR5_invdynJ_fixb_mdp_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5RRRRR5_invdynJ_fixb_mdp_slag_vp: pkin has to be [9x1] (double)');
assert(isreal(MDP) && all(size(MDP) == [23 1]), ...
  'S5RRRRR5_invdynJ_fixb_mdp_slag_vp: MDP has to be [23x1] (double)'); 

%% Symbolic Calculation
% From invdyn_joint_fixb_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 18:58:37
% EndTime: 2019-12-05 18:58:42
% DurationCPUTime: 2.40s
% Computational Cost: add. (2230->288), mult. (3035->369), div. (0->0), fcn. (1907->16), ass. (0->166)
t362 = cos(qJ(2));
t455 = pkin(1) * t362;
t334 = qJDD(1) * t455;
t347 = qJDD(1) + qJDD(2);
t357 = sin(qJ(2));
t447 = pkin(1) * qJD(1);
t418 = t357 * t447;
t277 = pkin(2) * t347 - qJD(2) * t418 + t334;
t349 = qJD(1) + qJD(2);
t426 = qJD(1) * t362;
t417 = pkin(1) * t426;
t298 = pkin(2) * t349 + t417;
t356 = sin(qJ(3));
t361 = cos(qJ(3));
t408 = qJD(2) * t426;
t419 = qJDD(1) * t357;
t375 = (t408 + t419) * pkin(1);
t462 = -t356 * t277 - (qJD(3) * t298 + t375) * t361;
t353 = qJ(1) + qJ(2);
t344 = qJ(3) + t353;
t328 = sin(t344);
t329 = cos(t344);
t461 = g(2) * t329 + g(3) * t328;
t458 = -g(2) * t328 + g(3) * t329;
t359 = cos(qJ(5));
t360 = cos(qJ(4));
t431 = t359 * t360;
t354 = sin(qJ(5));
t355 = sin(qJ(4));
t438 = t354 * t355;
t285 = -t431 + t438;
t286 = t354 * t360 + t355 * t359;
t338 = qJDD(3) + t347;
t339 = qJD(3) + t349;
t348 = qJD(4) + qJD(5);
t459 = t348 * t286;
t242 = t285 * t338 + t339 * t459;
t381 = t348 * t285;
t432 = t357 * t361;
t384 = t356 * t362 + t432;
t280 = t384 * t447;
t424 = qJD(3) * t356;
t397 = pkin(2) * t424 - t280;
t330 = pkin(2) * t356 + pkin(8);
t454 = pkin(2) * t361;
t331 = -pkin(3) - t454;
t364 = qJD(4) ^ 2;
t457 = -t330 * t364 - t331 * t338 - t397 * t339;
t456 = -pkin(8) - pkin(9);
t453 = pkin(3) * t338;
t452 = pkin(3) * t339;
t451 = pkin(4) * t360;
t332 = pkin(2) + t455;
t428 = pkin(1) * t432 + t356 * t332;
t279 = pkin(8) + t428;
t449 = -pkin(9) - t279;
t448 = -pkin(9) - t330;
t271 = t298 * t356 + t361 * t418;
t264 = pkin(8) * t339 + t271;
t407 = pkin(9) * t339 + t264;
t253 = t407 * t360;
t446 = t253 * t359;
t423 = qJD(3) * t361;
t257 = t332 * t424 + (t384 * qJD(2) + t357 * t423) * pkin(1);
t445 = t257 * t339;
t444 = t271 * t339;
t352 = qJ(4) + qJ(5);
t342 = cos(t352);
t443 = t328 * t342;
t442 = t329 * t342;
t441 = t338 * t355;
t440 = t338 * t360;
t439 = t339 * t355;
t435 = t355 * t360;
t433 = t356 * t357;
t320 = pkin(1) * t433;
t400 = qJD(3) * t418;
t393 = t356 * pkin(1) * t408 + qJDD(1) * t320 + t298 * t424 + (-t277 + t400) * t361;
t240 = t393 - t453;
t315 = t356 * t418;
t270 = t298 * t361 - t315;
t263 = -t270 - t452;
t421 = qJD(4) * t360;
t430 = t240 * t355 + t263 * t421;
t422 = qJD(4) * t355;
t336 = pkin(4) * t422;
t429 = t336 + t397;
t350 = t355 ^ 2;
t427 = -t360 ^ 2 + t350;
t420 = qJD(5) * t354;
t415 = pkin(2) * t423;
t414 = t339 * t438;
t413 = t339 * t431;
t412 = t263 * t422 + t461 * t360;
t341 = sin(t353);
t343 = cos(t353);
t411 = g(2) * t343 + g(3) * t341 + t334;
t333 = -pkin(3) - t451;
t410 = qJD(4) * t456;
t409 = t339 * t421;
t406 = -g(2) * t341 + g(3) * t343;
t405 = qJD(4) * t449;
t404 = qJD(4) * t448;
t403 = t332 * t361 - t320;
t402 = qJD(1) * (-qJD(2) + t349);
t401 = qJD(2) * (-qJD(1) - t349);
t383 = t339 * t422 - t440;
t237 = t383 * pkin(4) + t240;
t255 = t333 * t339 - t270;
t399 = g(2) * t442 + g(3) * t443 + t237 * t285 + t255 * t459;
t278 = -pkin(3) - t403;
t309 = t356 * t400;
t398 = t309 + t458;
t396 = -t271 + t336;
t252 = t407 * t355;
t251 = qJD(4) * pkin(4) - t252;
t392 = -t251 * t354 - t446;
t265 = t449 * t355;
t345 = t360 * pkin(9);
t266 = t279 * t360 + t345;
t391 = t265 * t359 - t266 * t354;
t390 = t265 * t354 + t266 * t359;
t283 = t448 * t355;
t284 = t330 * t360 + t345;
t389 = t283 * t359 - t284 * t354;
t388 = t283 * t354 + t284 * t359;
t311 = t456 * t355;
t312 = pkin(8) * t360 + t345;
t387 = t311 * t359 - t312 * t354;
t386 = t311 * t354 + t312 * t359;
t379 = -t393 + t461;
t378 = -pkin(8) * t364 + t444 + t453;
t377 = -t278 * t338 - t279 * t364 - t445;
t241 = qJD(5) * t413 + t286 * t338 - t348 * t414 + t359 * t409;
t272 = -t413 + t414;
t274 = t286 * t339;
t346 = qJDD(4) + qJDD(5);
t376 = t274 * t272 * MDP(17) + (t272 * t348 + t241) * MDP(19) + (t274 * t348 - t242) * MDP(20) + (-t272 ^ 2 + t274 ^ 2) * MDP(18) + t346 * MDP(21);
t374 = -pkin(8) * qJDD(4) + (t270 - t452) * qJD(4);
t239 = pkin(8) * t338 - t309 - t462;
t373 = -t263 * t339 - t239 + t458;
t256 = t332 * t423 + (-t357 * t424 + (t361 * t362 - t433) * qJD(2)) * pkin(1);
t372 = -qJDD(4) * t279 + (t278 * t339 - t256) * qJD(4);
t340 = sin(t352);
t371 = t237 * t286 - t255 * t381 - t340 * t461;
t370 = (-t241 * t285 - t242 * t286 + t272 * t381 - t274 * t459) * MDP(18) + (t241 * t286 - t274 * t381) * MDP(17) + (t286 * t346 - t348 * t381) * MDP(19) + (-t285 * t346 - t348 * t459) * MDP(20) + 0.2e1 * (-t427 * t339 * qJD(4) + t338 * t435) * MDP(11) + (t338 * t350 + 0.2e1 * t355 * t409) * MDP(10) + (qJDD(4) * t355 + t360 * t364) * MDP(12) + (qJDD(4) * t360 - t355 * t364) * MDP(13) + t338 * MDP(7);
t369 = t347 * MDP(4) + t370;
t281 = t361 * t417 - t315;
t368 = -qJDD(4) * t330 + (t331 * t339 + t281 - t415) * qJD(4);
t231 = -t264 * t421 + qJDD(4) * pkin(4) - t239 * t355 + (-t409 - t441) * pkin(9);
t367 = t253 * t420 + g(3) * t442 + g(1) * t340 + t255 * t272 + (-t253 * t348 - t231) * t354 - g(2) * t443;
t232 = -t383 * pkin(9) + t239 * t360 - t264 * t422;
t366 = -g(1) * t342 + t392 * qJD(5) + t359 * t231 - t354 * t232 - t255 * t274 + t458 * t340;
t365 = t398 + t462;
t363 = cos(qJ(1));
t358 = sin(qJ(1));
t305 = t333 - t454;
t292 = t360 * t410;
t291 = t355 * t410;
t276 = t278 - t451;
t268 = -t355 * t415 + t360 * t404;
t267 = t355 * t404 + t360 * t415;
t254 = t336 + t257;
t247 = -t256 * t355 + t360 * t405;
t246 = t256 * t360 + t355 * t405;
t1 = [((t347 * t362 + t357 * t401) * pkin(1) + t411) * MDP(5) + (t254 * t272 + t276 * t242 + (-t390 * qJD(5) - t246 * t354 + t247 * t359) * t348 + t391 * t346 + t399) * MDP(22) + (t372 * t360 + (-t377 - t461) * t355 + t430) * MDP(16) + (t254 * t274 + t276 * t241 - (t391 * qJD(5) + t246 * t359 + t247 * t354) * t348 - t390 * t346 + t371) * MDP(23) + (((-qJDD(1) - t347) * t357 + t362 * t401) * pkin(1) + t406) * MDP(6) + qJDD(1) * MDP(1) + (t372 * t355 + (-t240 + t377) * t360 + t412) * MDP(15) + (t403 * t338 + t379 - t445) * MDP(8) + (g(2) * t363 + g(3) * t358) * MDP(2) + (-g(2) * t358 + g(3) * t363) * MDP(3) + (-t256 * t339 - t428 * t338 + t365) * MDP(9) + t369; (t281 * t339 + (-pkin(2) * t338 - t277) * t356 + ((-pkin(2) * t339 - t298) * qJD(3) - t375) * t361 + t398) * MDP(9) + (t368 * t355 + (-t240 + t457) * t360 + t412) * MDP(15) + (t368 * t360 + (-t461 - t457) * t355 + t430) * MDP(16) + (t280 * t339 + (t338 * t361 - t339 * t424) * pkin(2) + t379) * MDP(8) + (t305 * t241 - (t389 * qJD(5) + t267 * t359 + t268 * t354) * t348 - t388 * t346 - t281 * t381 + t429 * t274 + t371) * MDP(23) + (t305 * t242 + (-t388 * qJD(5) - t267 * t354 + t268 * t359) * t348 + t389 * t346 + t281 * t459 + t429 * t272 + t399) * MDP(22) + ((t362 * t402 - t419) * pkin(1) + t406) * MDP(6) + t369 + (t357 * pkin(1) * t402 + t411) * MDP(5); (t374 * t360 + (-t378 - t461) * t355 + t430) * MDP(16) + (t333 * t241 - (t387 * qJD(5) + t291 * t359 + t292 * t354) * t348 - t386 * t346 + t396 * t274 - t270 * t381 + t371) * MDP(23) + (t374 * t355 + (-t240 + t378) * t360 + t412) * MDP(15) + (t333 * t242 + (-t386 * qJD(5) - t291 * t354 + t292 * t359) * t348 + t387 * t346 + t396 * t272 + t270 * t459 + t399) * MDP(22) + (t379 + t444) * MDP(8) + (t270 * t339 + t365) * MDP(9) + t370; MDP(12) * t441 + MDP(13) * t440 + qJDD(4) * MDP(14) + (-g(1) * t360 + t373 * t355) * MDP(15) + (g(1) * t355 + t373 * t360) * MDP(16) + (-(t252 * t354 - t446) * t348 + (-t272 * t439 + t346 * t359 - t348 * t420) * pkin(4) + t366) * MDP(22) + ((-qJD(5) * t251 - t252 * t348 - t232) * t359 + (-qJD(5) * t348 * t359 - t274 * t439 - t346 * t354) * pkin(4) + t367) * MDP(23) + t376 + (-MDP(10) * t435 + t427 * MDP(11)) * t339 ^ 2; (-t392 * t348 + t366) * MDP(22) + ((-t232 + (-qJD(5) + t348) * t251) * t359 + t367) * MDP(23) + t376;];
tau = t1;
