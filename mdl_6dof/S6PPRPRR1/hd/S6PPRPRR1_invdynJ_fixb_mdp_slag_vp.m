% Calculate vector of inverse dynamics joint torques for
% S6PPRPRR1
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% qJDD [6x1]
%   Generalized joint accelerations
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [13x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,alpha3,d3,d5,d6,theta1,theta2,theta4]';
% MDP [20x1]
%   Minimal dynamic parameter vector (fixed base model)
%   see S6PPRPRR1_convert_par2_MPV_fixb.m
% 
% Output:
% tau [6x1]
%   joint torques of inverse dynamics (contains inertial, gravitational Coriolis and centrifugal forces)

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-08 18:44
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tau = S6PPRPRR1_invdynJ_fixb_mdp_slag_vp(qJ, qJD, qJDD, g, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(3,1),zeros(13,1),zeros(20,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PPRPRR1_invdynJ_fixb_mdp_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6PPRPRR1_invdynJ_fixb_mdp_slag_vp: qJD has to be [6x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [6 1]), ...
  'S6PPRPRR1_invdynJ_fixb_mdp_slag_vp: qJDD has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6PPRPRR1_invdynJ_fixb_mdp_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [13 1]), ...
  'S6PPRPRR1_invdynJ_fixb_mdp_slag_vp: pkin has to be [13x1] (double)');
assert(isreal(MDP) && all(size(MDP) == [20 1]), ...
  'S6PPRPRR1_invdynJ_fixb_mdp_slag_vp: MDP has to be [20x1] (double)'); 

%% Symbolic Calculation
% From invdyn_joint_fixb_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-08 18:43:50
% EndTime: 2019-03-08 18:43:53
% DurationCPUTime: 3.58s
% Computational Cost: add. (2396->367), mult. (6303->545), div. (0->0), fcn. (6148->16), ass. (0->181)
t366 = cos(pkin(6));
t360 = sin(pkin(7));
t369 = sin(qJ(3));
t451 = t360 * t369;
t361 = sin(pkin(6));
t358 = sin(pkin(12));
t372 = cos(qJ(3));
t363 = cos(pkin(12));
t365 = cos(pkin(7));
t447 = t363 * t365;
t398 = t358 * t372 + t369 * t447;
t473 = t398 * t361;
t301 = t366 * t451 + t473;
t371 = cos(qJ(5));
t427 = qJD(3) * qJD(5);
t418 = t371 * t427;
t368 = sin(qJ(5));
t425 = qJDD(3) * t368;
t474 = qJD(5) * qJD(6) + t418 + t425;
t343 = qJD(1) * t366 + qJD(2);
t295 = qJD(1) * t473 + t343 * t451;
t455 = t358 * t369;
t399 = t372 * t447 - t455;
t450 = t360 * t372;
t300 = t399 * t361 + t366 * t450;
t357 = sin(pkin(13));
t362 = cos(pkin(13));
t284 = t300 * t357 + t301 * t362;
t452 = t360 * t361;
t423 = t363 * t452;
t317 = t365 * t366 - t423;
t311 = t317 * t371;
t472 = -t284 * t368 + t311;
t456 = t358 * t361;
t416 = qJDD(1) * t456;
t340 = t366 * qJDD(1) + qJDD(2);
t415 = qJDD(1) * t447;
t442 = t372 * t361 * t415 + t340 * t450;
t279 = qJDD(3) * pkin(3) - t295 * qJD(3) - t369 * t416 + t442;
t428 = qJD(1) * qJD(3);
t280 = (qJD(3) * t343 * t372 + t369 * t340) * t360 + (t398 * qJDD(1) + t399 * t428) * t361;
t250 = t357 * t279 + t362 * t280;
t248 = qJDD(3) * pkin(9) + t250;
t292 = t362 * t295;
t440 = qJD(1) * t361;
t420 = t363 * t440;
t409 = t365 * t420;
t294 = t343 * t450 + t372 * t409 - t440 * t455;
t293 = qJD(3) * pkin(3) + t294;
t273 = t357 * t293 + t292;
t271 = qJD(3) * pkin(9) + t273;
t308 = t343 * t365 - t360 * t420 + qJD(4);
t258 = t271 * t371 + t308 * t368;
t307 = -qJDD(1) * t423 + t365 * t340 + qJDD(4);
t464 = t307 * t371;
t240 = -qJDD(5) * pkin(5) + qJD(5) * t258 + t248 * t368 - t464;
t438 = qJD(3) * t371;
t344 = -qJD(6) + t438;
t359 = sin(pkin(11));
t364 = cos(pkin(11));
t446 = t364 * t366;
t318 = -t358 * t359 + t363 * t446;
t448 = t361 * t365;
t302 = -t318 * t360 - t364 * t448;
t453 = t359 * t366;
t320 = -t358 * t364 - t363 * t453;
t303 = -t320 * t360 + t359 * t448;
t401 = t357 * t372 + t362 * t369;
t313 = t401 * t360;
t315 = t401 * t365;
t444 = t372 * t362;
t328 = t357 * t369 - t444;
t382 = t313 * t366 + (t315 * t363 - t328 * t358) * t361;
t321 = -t358 * t453 + t363 * t364;
t454 = t359 * t361;
t385 = t313 * t454 + t315 * t320 - t321 * t328;
t319 = t358 * t446 + t359 * t363;
t449 = t361 * t364;
t386 = -t313 * t449 + t315 * t318 - t319 * t328;
t394 = g(1) * (t303 * t371 - t368 * t385) + g(2) * (t302 * t371 - t368 * t386) + g(3) * (-t368 * t382 + t311);
t408 = pkin(5) * t368 - pkin(10) * t371;
t471 = (pkin(10) * qJD(6) + t408 * qJD(3)) * t344 - t240 - t394;
t275 = t294 * t357 + t292;
t400 = -pkin(5) * t371 - pkin(10) * t368 - pkin(4);
t468 = pkin(3) * t362;
t322 = t400 - t468;
t424 = t371 * qJDD(3);
t327 = t368 * t427 + qJDD(6) - t424;
t334 = t408 * qJD(5);
t470 = (t275 - t334) * t344 + t322 * t327;
t256 = qJD(5) * pkin(10) + t258;
t312 = t357 * t451 - t360 * t444;
t314 = t328 * t365;
t392 = g(1) * (-t312 * t454 - t314 * t320 - t321 * t401) + g(2) * (t312 * t449 - t314 * t318 - t319 * t401) + g(3) * (-t312 * t366 + (-t314 * t363 - t358 * t401) * t361);
t346 = pkin(3) * t357 + pkin(9);
t459 = t344 * t346;
t469 = (t256 + t459) * qJD(6) - t392;
t367 = sin(qJ(6));
t370 = cos(qJ(6));
t410 = t367 * qJDD(5) + t474 * t370;
t431 = qJD(6) * t368;
t417 = qJD(3) * t431;
t296 = -t367 * t417 + t410;
t466 = t296 * t367;
t465 = t307 * t368;
t463 = t317 * t368;
t429 = t370 * qJD(5);
t439 = qJD(3) * t368;
t330 = t367 * t439 - t429;
t461 = t330 * t344;
t435 = qJD(5) * t367;
t332 = t370 * t439 + t435;
t460 = t332 * t344;
t458 = t344 * t370;
t457 = t344 * t371;
t291 = t357 * t295;
t445 = t370 * t327;
t355 = t368 ^ 2;
t441 = -t371 ^ 2 + t355;
t437 = qJD(5) * t330;
t436 = qJD(5) * t346;
t434 = qJD(5) * t368;
t433 = qJD(5) * t371;
t432 = qJD(6) * t367;
t430 = qJD(6) * t370;
t419 = t344 * t435;
t249 = t279 * t362 - t357 * t280;
t272 = t293 * t362 - t291;
t270 = -qJD(3) * pkin(4) - t272;
t412 = -qJD(3) * t270 - t248;
t259 = t400 * qJD(3) - t272;
t404 = t271 * t368 - t308 * t371;
t411 = -qJDD(5) * pkin(10) + t404 * qJD(5) - qJD(6) * t259 - t248 * t371 - t465;
t407 = g(1) * t321 + g(2) * t319;
t261 = t284 * t371 + t463;
t283 = -t300 * t362 + t301 * t357;
t406 = t261 * t370 + t283 * t367;
t405 = -t261 * t367 + t283 * t370;
t305 = t313 * t371 + t365 * t368;
t403 = t305 * t370 + t312 * t367;
t402 = -t305 * t367 + t312 * t370;
t304 = t313 * t368 - t365 * t371;
t397 = -t367 * t327 + t344 * t430;
t396 = -t344 * t432 - t445;
t395 = -g(1) * t454 + g(2) * t449 - g(3) * t366;
t393 = g(1) * (t303 * t368 + t371 * t385) + g(2) * (t302 * t368 + t371 * t386) + g(3) * (t371 * t382 + t463);
t391 = -t367 * t431 + t371 * t429;
t387 = t407 - t416;
t384 = -qJD(6) * t256 + t392;
t276 = t294 * t362 - t291;
t347 = -pkin(4) - t468;
t383 = -qJDD(5) * t346 + (qJD(3) * t347 + t270 + t276) * qJD(5);
t381 = qJD(6) * t322 * t344 - g(1) * t385 - g(2) * t386 - g(3) * t382;
t380 = t393 + t411;
t255 = -qJD(5) * pkin(5) + t404;
t379 = -pkin(10) * t327 + (-t255 + t404) * t344;
t378 = qJD(5) * t255 - t276 * t344 - t327 * t346 - t411;
t376 = -g(1) * (t320 * t365 + t359 * t452) - g(2) * (t318 * t365 - t360 * t449);
t373 = qJD(5) ^ 2;
t375 = -qJD(3) * t275 + t346 * t373 - t249 + t392 + (-pkin(4) + t347) * qJDD(3);
t374 = qJD(3) ^ 2;
t351 = t370 * qJDD(5);
t337 = qJDD(5) * t371 - t368 * t373;
t336 = qJDD(5) * t368 + t371 * t373;
t316 = t332 * t434;
t310 = t328 * t360 * qJD(3);
t309 = qJD(3) * t313;
t299 = t301 * qJD(3);
t298 = t300 * qJD(3);
t297 = t370 * t417 - t351 + (t425 + (qJD(6) + t438) * qJD(5)) * t367;
t289 = t305 * qJD(5) - t310 * t368;
t288 = -t304 * qJD(5) - t310 * t371;
t282 = t298 * t362 - t299 * t357;
t281 = t298 * t357 + t299 * t362;
t246 = t261 * qJD(5) + t282 * t368;
t245 = qJD(5) * t472 + t282 * t371;
t244 = qJD(3) * t334 + t400 * qJDD(3) - t249;
t243 = t370 * t244;
t242 = t256 * t370 + t259 * t367;
t241 = -t256 * t367 + t259 * t370;
t1 = [(qJDD(1) - g(3)) * MDP(1) + (t340 * t366 - g(3) + (t358 ^ 2 + t363 ^ 2) * t361 ^ 2 * qJDD(1)) * MDP(2) + (-t249 * t283 + t250 * t284 - t272 * t281 + t273 * t282 + t307 * t317 - g(3)) * MDP(6) + (-qJD(5) * t246 + qJDD(5) * t472) * MDP(12) + (-qJD(5) * t245 - qJDD(5) * t261) * MDP(13) + (-(-t406 * qJD(6) - t245 * t367 + t281 * t370) * t344 + t405 * t327 + t246 * t330 - t472 * t297) * MDP(19) + ((t405 * qJD(6) + t245 * t370 + t281 * t367) * t344 - t406 * t327 + t246 * t332 - t472 * t296) * MDP(20) + (t300 * MDP(4) - t301 * MDP(5) + (-MDP(12) * t371 + MDP(13) * t368) * t283) * qJDD(3) + (-t299 * MDP(4) - t298 * MDP(5) + (-t281 * t371 + t283 * t434) * MDP(12) + (t281 * t368 + t283 * t433) * MDP(13)) * qJD(3); (t395 + t340) * MDP(2) + (-t249 * t312 + t250 * t313 - t272 * t309 - t273 * t310 + t307 * t365 + t395) * MDP(6) + (-qJD(5) * t289 - qJDD(5) * t304 - t312 * t424) * MDP(12) + (-qJD(5) * t288 - qJDD(5) * t305 + t312 * t425) * MDP(13) + (-(-t403 * qJD(6) - t288 * t367 + t309 * t370) * t344 + t402 * t327 + t289 * t330 + t304 * t297) * MDP(19) + ((t402 * qJD(6) + t288 * t370 + t309 * t367) * t344 - t403 * t327 + t289 * t332 + t304 * t296) * MDP(20) + ((-t309 * t371 + t312 * t434) * MDP(12) + (t309 * t368 + t312 * t433) * MDP(13)) * qJD(3) + ((qJDD(3) * t372 - t369 * t374) * MDP(4) + (-qJDD(3) * t369 - t372 * t374) * MDP(5)) * t360; qJDD(3) * MDP(3) + (-g(3) * t300 + t387 * t369 + t376 * t372 + t442) * MDP(4) + (g(3) * t301 + t294 * qJD(3) + ((-t343 * t360 - t409) * qJD(3) + t387) * t372 + (-t360 * t340 + (g(1) * t320 + g(2) * t318) * t365 + (t358 * t428 - t415 + (g(1) * t359 - g(2) * t364) * t360) * t361) * t369) * MDP(5) + (t272 * t275 - t273 * t276 + (t249 * t362 + t250 * t357 + (g(3) * t456 + t407) * t369 + (-g(3) * (t360 * t366 + t361 * t447) + t376) * t372) * pkin(3)) * MDP(6) + (qJDD(3) * t355 + 0.2e1 * t368 * t418) * MDP(7) + 0.2e1 * (t368 * t424 - t441 * t427) * MDP(8) + t336 * MDP(9) + t337 * MDP(10) + (t383 * t368 - t375 * t371) * MDP(12) + (t375 * t368 + t383 * t371) * MDP(13) + (t296 * t368 * t370 + t391 * t332) * MDP(14) + ((-t330 * t370 - t332 * t367) * t433 + (-t466 - t297 * t370 + (t330 * t367 - t332 * t370) * qJD(6)) * t368) * MDP(15) + (-t296 * t371 - t391 * t344 + t368 * t445 + t316) * MDP(16) + ((t297 + t419) * t371 + (t397 - t437) * t368) * MDP(17) + (-t327 * t371 - t344 * t434) * MDP(18) + (t470 * t370 + t381 * t367 + (t330 * t436 + t378 * t367 + t370 * t469 - t243) * t371 + (t255 * t430 + t240 * t367 - t276 * t330 + t346 * t297 + (-t367 * t459 + t241) * qJD(5)) * t368) * MDP(19) + (-t470 * t367 + t381 * t370 + (t332 * t436 + t378 * t370 + (t244 - t469) * t367) * t371 + (-t255 * t432 + t240 * t370 - t276 * t332 + t346 * t296 + (-t346 * t458 - t242) * qJD(5)) * t368) * MDP(20); (-g(1) * t303 - g(2) * t302 - g(3) * t317 + t307) * MDP(6) + t337 * MDP(12) - t336 * MDP(13) + t316 * MDP(20) + ((-t297 + t419) * MDP(19) + (t344 * t429 - t296) * MDP(20)) * t371 + ((t397 + t437) * MDP(19) + t396 * MDP(20)) * t368; MDP(9) * t425 + MDP(10) * t424 + qJDD(5) * MDP(11) + (t412 * t368 - t394 + t464) * MDP(12) + (t412 * t371 + t393 - t465) * MDP(13) + (-t332 * t458 + t466) * MDP(14) + ((t296 + t461) * t370 + (-t297 + t460) * t367) * MDP(15) + ((-t332 * t368 + t370 * t457) * qJD(3) - t397) * MDP(16) + ((t330 * t368 - t367 * t457) * qJD(3) - t396) * MDP(17) + t344 * MDP(18) * t439 + (-pkin(5) * t297 - t241 * t439 - t258 * t330 + t379 * t367 + t370 * t471) * MDP(19) + (-pkin(5) * t296 + t242 * t439 - t258 * t332 - t367 * t471 + t379 * t370) * MDP(20) + (-t368 * t371 * MDP(7) + t441 * MDP(8)) * t374; t332 * t330 * MDP(14) + (-t330 ^ 2 + t332 ^ 2) * MDP(15) + (t410 - t461) * MDP(16) + (t351 - t460) * MDP(17) + t327 * MDP(18) + (-t242 * t344 - t255 * t332 + t243) * MDP(19) + (-t241 * t344 + t255 * t330) * MDP(20) + (-MDP(17) * t417 + t384 * MDP(19) + t380 * MDP(20)) * t370 + (-MDP(16) * t417 - t474 * MDP(17) + t380 * MDP(19) + (-t244 - t384) * MDP(20)) * t367;];
tau  = t1;
