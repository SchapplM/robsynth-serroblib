% Calculate vector of inverse dynamics joint torques for
% S5RRRPR4
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
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d2,d3,d5]';
% MDP [24x1]
%   Minimal dynamic parameter vector (fixed base model)
%   see S5RRRPR4_convert_par2_MPV_fixb.m
% 
% Output:
% tau [5x1]
%   joint torques of inverse dynamics (contains inertial, gravitational Coriolis and centrifugal forces)

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 21:12
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tau = S5RRRPR4_invdynJ_fixb_mdp_slag_vp(qJ, qJD, qJDD, g, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(5,1),zeros(3,1),zeros(8,1),zeros(24,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRRPR4_invdynJ_fixb_mdp_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RRRPR4_invdynJ_fixb_mdp_slag_vp: qJD has to be [5x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [5 1]), ...
  'S5RRRPR4_invdynJ_fixb_mdp_slag_vp: qJDD has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RRRPR4_invdynJ_fixb_mdp_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RRRPR4_invdynJ_fixb_mdp_slag_vp: pkin has to be [8x1] (double)');
assert(isreal(MDP) && all(size(MDP) == [24 1]), ...
  'S5RRRPR4_invdynJ_fixb_mdp_slag_vp: MDP has to be [24x1] (double)'); 

%% Symbolic Calculation
% From invdyn_joint_fixb_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 21:11:48
% EndTime: 2019-12-31 21:11:52
% DurationCPUTime: 2.68s
% Computational Cost: add. (1730->352), mult. (2531->427), div. (0->0), fcn. (1494->10), ass. (0->186)
t359 = qJ(1) + qJ(2);
t345 = cos(t359);
t465 = g(2) * t345;
t361 = sin(qJ(3));
t346 = t361 * qJ(4);
t365 = cos(qJ(3));
t348 = t365 * pkin(3);
t443 = t348 + t346;
t429 = pkin(2) + t443;
t351 = qJDD(1) + qJDD(2);
t353 = qJD(1) + qJD(2);
t360 = sin(qJ(5));
t364 = cos(qJ(5));
t300 = t360 * t361 + t364 * t365;
t381 = t300 * qJD(5);
t435 = qJD(3) * t365;
t419 = t353 * t435;
t436 = qJD(3) * t361;
t420 = t353 * t436;
t449 = t361 * t364;
t422 = t351 * t449 + t360 * t420 + t364 * t419;
t450 = t360 * t365;
t246 = -t351 * t450 - t353 * t381 + t422;
t433 = qJD(5) * t365;
t479 = -t360 * t433 - t364 * t436;
t344 = sin(t359);
t444 = g(1) * t345 + g(2) * t344;
t357 = t361 ^ 2;
t358 = t365 ^ 2;
t441 = t357 + t358;
t477 = t353 * t441;
t362 = sin(qJ(2));
t431 = qJDD(1) * t362;
t366 = cos(qJ(2));
t438 = qJD(2) * t366;
t286 = pkin(7) * t351 + (qJD(1) * t438 + t431) * pkin(1);
t274 = t365 * t286;
t463 = pkin(1) * qJD(1);
t428 = t362 * t463;
t307 = pkin(7) * t353 + t428;
t354 = qJDD(3) * qJ(4);
t355 = qJD(3) * qJD(4);
t254 = -t307 * t436 + t274 + t354 + t355;
t273 = t361 * t286;
t291 = t307 * t435;
t415 = qJDD(4) + t273 + t291;
t460 = qJDD(3) * pkin(3);
t255 = t415 - t460;
t476 = t254 * t365 + t255 * t361;
t475 = t449 - t450;
t294 = t361 * t307;
t474 = qJD(4) + t294;
t352 = qJD(3) - qJD(5);
t473 = qJD(5) + t352;
t282 = t300 * t353;
t284 = t475 * t353;
t350 = qJDD(3) - qJDD(5);
t472 = t284 * t282 * MDP(18) - (t282 ^ 2 - t284 ^ 2) * MDP(19) - t350 * MDP(22);
t471 = pkin(3) + pkin(4);
t470 = pkin(7) - pkin(8);
t469 = pkin(1) * t366;
t468 = pkin(2) * t351;
t467 = pkin(2) * t353;
t369 = qJD(3) ^ 2;
t466 = pkin(7) * t369;
t334 = g(1) * t344;
t337 = pkin(1) * t362 + pkin(7);
t464 = -pkin(8) + t337;
t462 = pkin(7) * qJDD(3);
t461 = qJ(4) * t365;
t459 = t282 * t352;
t458 = t284 * t352;
t457 = t337 * t369;
t456 = t344 * t361;
t455 = t345 * t361;
t454 = t351 * t361;
t453 = t351 * t365;
t452 = t352 * t366;
t451 = t353 * t361;
t448 = t361 * t365;
t295 = t365 * t307;
t447 = t479 * t353;
t446 = g(1) * t456 - g(2) * t455;
t439 = qJD(2) * t362;
t426 = pkin(1) * t439;
t445 = qJD(1) * t426 - qJDD(1) * t469;
t442 = t357 - t358;
t440 = qJD(1) * t366;
t356 = qJD(3) * qJ(4);
t437 = qJD(3) * t353;
t434 = qJD(5) * t361;
t343 = t361 * qJD(4);
t432 = -pkin(8) * t451 + t474;
t430 = qJDD(3) * t337;
t427 = pkin(1) * t440;
t425 = pkin(1) * t438;
t319 = t470 * t365;
t349 = t353 ^ 2;
t424 = t349 * t448;
t325 = t365 * t334;
t410 = t353 * t428;
t423 = t365 * t410 + t427 * t436 + t325;
t289 = pkin(3) * t436 - qJ(4) * t435 - t343;
t338 = -pkin(2) - t469;
t421 = t353 * t439;
t299 = t464 * t365;
t414 = pkin(2) + t346;
t285 = t445 - t468;
t413 = -t285 - t465;
t412 = t441 * t351;
t411 = -qJD(3) * pkin(3) + qJD(4);
t308 = -t427 - t467;
t409 = t285 * t361 + t308 * t435 - t446;
t270 = -pkin(8) * t353 * t365 + t295;
t408 = t466 - t468;
t363 = sin(qJ(1));
t367 = cos(qJ(1));
t406 = g(1) * t363 - g(2) * t367;
t405 = pkin(3) * t361 - t461;
t296 = t338 - t443;
t404 = qJ(4) * t364 - t360 * t471;
t403 = -qJ(4) * t360 - t364 * t471;
t281 = t294 + t411;
t288 = t295 + t356;
t402 = t281 * t361 + t288 * t365;
t298 = t464 * t361;
t401 = t298 * t364 - t299 * t360;
t400 = t298 * t360 + t299 * t364;
t318 = t470 * t361;
t399 = t318 * t364 - t319 * t360;
t398 = t318 * t360 + t319 * t364;
t396 = t334 - t445 - t465;
t395 = g(1) * t455 + g(2) * t456 - g(3) * t365 - t273;
t394 = -t414 - t348;
t249 = (t405 * qJD(3) - t343) * t353 + t394 * t351 + t445;
t393 = t351 * t429 - t249 - t466;
t392 = t296 * t353 - t425;
t391 = -t471 * t361 + t461;
t390 = -t419 - t454;
t389 = -t353 * t434 - t453;
t245 = t390 * pkin(8) - t471 * qJDD(3) + t415;
t264 = t270 + t356;
t388 = t473 * t264 - t245;
t248 = (t420 - t453) * pkin(8) + t254;
t259 = -t471 * qJD(3) + t432;
t387 = -t473 * t259 - t248;
t271 = -pkin(4) * t436 - t289;
t380 = t471 * t365 + t414;
t242 = (t391 * qJD(3) + t343) * t353 + t380 * t351 - t445;
t256 = t380 * t353 + t427;
t379 = t360 * t435 + t364 * t434;
t260 = t379 + t479;
t277 = t300 * t344;
t279 = t300 * t345;
t385 = g(1) * t277 - g(2) * t279 + t242 * t300 + t256 * t260;
t261 = t300 * qJD(3) - t381;
t276 = t475 * t344;
t278 = t475 * t345;
t384 = g(1) * t276 - g(2) * t278 + t242 * t475 + t256 * t261;
t382 = -qJDD(4) + t395;
t272 = t289 + t426;
t378 = -t272 * t353 - t296 * t351 - t249 - t457;
t377 = t281 * t435 - t288 * t436 - t444 + t476;
t247 = t300 * t351 + t379 * t353 + t447;
t376 = (-t246 * t300 - t247 * t475 - t260 * t284 - t261 * t282) * MDP(19) + (t246 * t475 + t261 * t284) * MDP(18) + (t260 * t352 + t300 * t350) * MDP(21) + (-t261 * t352 - t350 * t475) * MDP(20) + 0.2e1 * (t351 * t448 - t442 * t437) * MDP(8) + (t351 * t357 + 0.2e1 * t361 * t419) * MDP(7) + (qJDD(3) * t365 - t361 * t369) * MDP(10) + (qJDD(3) * t361 + t365 * t369) * MDP(9) + t351 * MDP(4);
t375 = -g(1) * t278 - g(2) * t276 + g(3) * t300 - t256 * t284;
t374 = g(1) * t279 + g(2) * t277 + g(3) * t475 + t256 * t282;
t373 = pkin(1) * t421 + t338 * t351 + t457;
t372 = -t430 + (t338 * t353 - t425) * qJD(3);
t371 = (t281 * t365 - t288 * t361) * qJD(3) + t476;
t370 = -t444 * pkin(7) - t394 * t334 - t429 * t465;
t347 = t365 * pkin(4);
t341 = pkin(8) * t436;
t305 = qJD(3) * t319;
t304 = -pkin(7) * t436 + t341;
t297 = t347 + t429;
t292 = t308 * t436;
t290 = t405 * t353;
t287 = -t296 + t347;
t268 = t391 * t353;
t266 = qJD(3) * t299 + t361 * t425;
t265 = -t337 * t436 + t365 * t425 + t341;
t263 = t394 * t353 - t427;
t262 = t271 - t426;
t257 = t263 * t436;
t1 = [(t292 + t325 + t372 * t361 + (-t373 + t413) * t365) * MDP(12) + (t373 * t361 + t372 * t365 + t409) * MDP(13) + (t249 * t296 + t263 * t272 + (t402 * t438 + t406) * pkin(1) + t371 * t337 + t370) * MDP(17) + t376 + qJDD(1) * MDP(1) + t406 * MDP(2) + (g(1) * t367 + g(2) * t363) * MDP(3) + (t262 * t282 + t287 * t247 - (-t400 * qJD(5) - t265 * t360 + t266 * t364) * t352 - t401 * t350 + t385) * MDP(23) + ((t430 + (-t263 - t392) * qJD(3)) * t365 + t378 * t361 + t446) * MDP(16) + (t257 + t325 + (t392 * qJD(3) - t430) * t361 + (t378 - t465) * t365) * MDP(14) + (t337 * t412 + t425 * t477 + t377) * MDP(15) + ((t351 * t366 - t421) * pkin(1) + t396) * MDP(5) + (((-qJDD(1) - t351) * t362 + (-qJD(1) - t353) * t438) * pkin(1) + t444) * MDP(6) + (t262 * t284 + t287 * t246 + (t401 * qJD(5) + t265 * t364 + t266 * t360) * t352 + t400 * t350 + t384) * MDP(24); ((t462 + (t353 * t429 - t263 - t427) * qJD(3)) * t365 + ((-t289 + t428) * t353 + t393) * t361 + t446) * MDP(16) + (-t249 * t429 + t263 * t289 + (-t263 * t362 - t402 * t366) * t463 + t371 * pkin(7) + t370) * MDP(17) + t376 + (t396 + t410) * MDP(5) + (t271 * t282 + t297 * t247 - (-t398 * qJD(5) - t304 * t360 + t305 * t364) * t352 - t399 * t350 + (t362 * t282 + t452 * t475) * t463 + t385) * MDP(23) + (t271 * t284 + t297 * t246 + (t399 * qJD(5) + t304 * t364 + t305 * t360) * t352 + t398 * t350 + (t362 * t284 - t300 * t452) * t463 + t384) * MDP(24) + (t292 + (-pkin(2) * t437 - t462) * t361 + (-t408 + t413) * t365 + t423) * MDP(12) + (pkin(7) * t412 - t427 * t477 + t377) * MDP(15) + ((-t462 + (t427 - t467) * qJD(3)) * t365 + (t408 - t410) * t361 + t409) * MDP(13) + ((-t431 + (-qJD(2) + t353) * t440) * pkin(1) + t444) * MDP(6) + (t257 + (-t429 * t437 - t462) * t361 + (-t289 * t353 + t393 - t465) * t365 + t423) * MDP(14); -MDP(7) * t424 + t442 * MDP(8) * t349 + MDP(9) * t454 + MDP(10) * t453 + qJDD(3) * MDP(11) + (-t308 * t451 + t395) * MDP(12) + (g(3) * t361 - t274 + (-t308 * t353 + t444) * t365) * MDP(13) + (0.2e1 * t460 + (-t263 * t361 + t290 * t365) * t353 + t382) * MDP(14) + (-t405 * t351 + ((t288 - t356) * t361 + (-t281 + t411) * t365) * t353) * MDP(15) + (t274 + 0.2e1 * t354 + 0.2e1 * t355 + (t290 * t353 - g(3)) * t361 + (t263 * t353 - t444) * t365) * MDP(16) + (-t255 * pkin(3) - g(3) * t443 + t254 * qJ(4) - t263 * t290 - t281 * t295 + t474 * t288 + t444 * t405) * MDP(17) + (t459 - t246) * MDP(20) + (t247 + t458) * MDP(21) + (-t403 * t350 + t360 * t248 - t364 * t245 - t268 * t282 + (t364 * t270 + t432 * t360) * t352 + (t360 * t259 + t364 * t264 + t404 * t352) * qJD(5) - t375) * MDP(23) + (t404 * t350 + t364 * t248 + t360 * t245 - t268 * t284 + (-t360 * t270 + t432 * t364) * t352 + (t364 * t259 - t360 * t264 + t403 * t352) * qJD(5) - t374) * MDP(24) - t472; (-qJDD(3) - t424) * MDP(14) + MDP(15) * t454 + (-t349 * t357 - t369) * MDP(16) + (-qJD(3) * t288 + t263 * t451 + t291 - t382 - t460) * MDP(17) + (-t282 * t451 - t364 * t350) * MDP(23) + (-t284 * t451 + t360 * t350) * MDP(24) + (-MDP(23) * t360 - MDP(24) * t364) * t352 ^ 2; (t422 - t459) * MDP(20) + (-t447 - t458) * MDP(21) + t375 * MDP(23) + t374 * MDP(24) + (-t353 * MDP(20) * t433 + t389 * MDP(21) - t388 * MDP(23) + t387 * MDP(24)) * t364 + (t389 * MDP(20) + t390 * MDP(21) + t387 * MDP(23) + t388 * MDP(24)) * t360 + t472;];
tau = t1;
