% Calculate vector of inverse dynamics joint torques for
% S5RRPRP5
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
%   pkin=[a2,a3,a4,a5,d1,d2,d4,theta3]';
% MDP [25x1]
%   Minimal dynamic parameter vector (fixed base model)
%   see S5RRPRP5_convert_par2_MPV_fixb.m
% 
% Output:
% tau [5x1]
%   joint torques of inverse dynamics (contains inertial, gravitational Coriolis and centrifugal forces)

% Quelle: HybrDyn-Toolbox
% Datum: 2021-01-15 20:19
% Revision: 24b2e7d74a0c1a3b64fa2f8f5ad758691ad61af3 (2021-01-15)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tau = S5RRPRP5_invdynJ_fixb_mdp_slag_vp(qJ, qJD, qJDD, g, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(5,1),zeros(3,1),zeros(8,1),zeros(25,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRPRP5_invdynJ_fixb_mdp_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RRPRP5_invdynJ_fixb_mdp_slag_vp: qJD has to be [5x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [5 1]), ...
  'S5RRPRP5_invdynJ_fixb_mdp_slag_vp: qJDD has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RRPRP5_invdynJ_fixb_mdp_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RRPRP5_invdynJ_fixb_mdp_slag_vp: pkin has to be [8x1] (double)');
assert(isreal(MDP) && all(size(MDP) == [25 1]), ...
  'S5RRPRP5_invdynJ_fixb_mdp_slag_vp: MDP has to be [25x1] (double)'); 

%% Symbolic Calculation
% From invdyn_joint_fixb_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2021-01-15 20:18:28
% EndTime: 2021-01-15 20:18:41
% DurationCPUTime: 4.61s
% Computational Cost: add. (3424->371), mult. (7936->448), div. (0->0), fcn. (5721->12), ass. (0->167)
t401 = cos(pkin(8));
t400 = sin(pkin(8));
t404 = sin(qJ(2));
t406 = cos(qJ(2));
t350 = t400 * t406 + t401 * t404;
t452 = qJD(1) * qJD(2);
t445 = t404 * t452;
t419 = qJDD(1) * t350 - t400 * t445;
t444 = t406 * t452;
t309 = t401 * t444 + t419;
t462 = t401 * t406;
t433 = t400 * t404 - t462;
t341 = t433 * qJD(1);
t343 = t350 * qJD(1);
t403 = sin(qJ(4));
t451 = qJDD(1) * t404;
t447 = -t401 * t445 + (-t444 - t451) * t400;
t450 = qJDD(1) * t406;
t424 = t401 * t450 + t447;
t481 = cos(qJ(4));
t446 = qJD(4) * t481;
t455 = qJD(4) * t403;
t427 = t481 * t309 - t341 * t446 - t343 * t455 + t403 * t424;
t299 = -t481 * t341 - t343 * t403;
t396 = qJD(2) + qJD(4);
t468 = t299 * t396;
t254 = t427 - t468;
t394 = qJDD(2) + qJDD(4);
t428 = t403 * t341 - t343 * t481;
t413 = qJD(4) * t428 - t403 * t309 + t481 * t424;
t469 = t299 ^ 2;
t470 = t428 * t396;
t494 = t428 ^ 2;
t496 = t299 * t428;
t497 = t254 * MDP(17) + MDP(15) * t496 + t394 * MDP(19) + (t413 - t470) * MDP(18) + (-t469 + t494) * MDP(16);
t474 = t406 * pkin(2);
t384 = pkin(1) + t474;
t359 = -qJD(1) * t384 + qJD(3);
t315 = pkin(3) * t341 + t359;
t267 = -pkin(4) * t299 + qJ(5) * t428 + t315;
t493 = t267 * t299;
t492 = t315 * t299;
t491 = t315 * t428;
t397 = qJ(2) + pkin(8);
t391 = cos(t397);
t392 = qJ(4) + t397;
t380 = sin(t392);
t381 = cos(t392);
t459 = t381 * pkin(4) + t380 * qJ(5);
t490 = pkin(3) * t391 + t459 + t474;
t389 = t394 * pkin(4);
t484 = qJDD(5) - t389;
t489 = -t267 * t428 + t484;
t274 = -pkin(4) * t428 - qJ(5) * t299;
t402 = -qJ(3) - pkin(6);
t367 = t402 * t404;
t355 = qJD(1) * t367;
t368 = t402 * t406;
t356 = qJD(1) * t368;
t463 = t401 * t356;
t313 = -t355 * t400 + t463;
t478 = pkin(7) * t341;
t288 = t313 + t478;
t345 = t400 * t356;
t314 = t401 * t355 + t345;
t477 = pkin(7) * t343;
t289 = t314 - t477;
t382 = pkin(2) * t401 + pkin(3);
t480 = pkin(2) * t400;
t448 = t403 * t480;
t487 = qJD(4) * t448 + t403 * t288 + t289 * t481 - t382 * t446;
t458 = t403 * t382 + t481 * t480;
t385 = t394 * qJ(5);
t388 = t396 * qJD(5);
t486 = t385 + t388;
t405 = sin(qJ(1));
t407 = cos(qJ(1));
t485 = g(1) * t405 - g(2) * t407;
t440 = qJD(2) * t402;
t338 = qJD(3) * t406 + t404 * t440;
t339 = -qJD(3) * t404 + t406 * t440;
t290 = -t400 * t338 + t401 * t339;
t425 = t433 * qJD(2);
t281 = pkin(7) * t425 + t290;
t291 = t401 * t338 + t400 * t339;
t426 = t350 * qJD(2);
t282 = -pkin(7) * t426 + t291;
t316 = t401 * t367 + t368 * t400;
t292 = -pkin(7) * t350 + t316;
t317 = t400 * t367 - t401 * t368;
t293 = -pkin(7) * t433 + t317;
t429 = t292 * t481 - t403 * t293;
t251 = qJD(4) * t429 + t403 * t281 + t282 * t481;
t273 = t403 * t292 + t293 * t481;
t483 = t251 * t396 + t273 * t394 + t380 * t485;
t479 = pkin(2) * t404;
t475 = g(3) * t406;
t473 = qJD(2) * pkin(2);
t349 = t355 + t473;
t307 = t401 * t349 + t345;
t285 = qJD(2) * pkin(3) + t307 - t477;
t308 = t400 * t349 - t463;
t287 = t308 - t478;
t265 = t403 * t285 + t287 * t481;
t472 = t265 * t396;
t467 = t380 * t405;
t466 = t380 * t407;
t465 = t381 * t405;
t464 = t381 * t407;
t306 = qJDD(2) * pkin(2) + qJD(1) * t339 + qJDD(1) * t367;
t312 = qJD(1) * t338 - qJDD(1) * t368;
t276 = t400 * t306 + t401 * t312;
t461 = qJD(5) - t487;
t460 = t458 * qJD(4) + t288 * t481 - t403 * t289;
t398 = t404 ^ 2;
t457 = -t406 ^ 2 + t398;
t456 = qJD(1) * t404;
t264 = t285 * t481 - t403 * t287;
t453 = qJD(5) - t264;
t449 = pkin(2) * t445 + qJDD(3);
t387 = t404 * t473;
t318 = pkin(2) * t456 + pkin(3) * t343;
t390 = sin(t397);
t442 = -pkin(3) * t390 - pkin(4) * t380 - t479;
t275 = t401 * t306 - t312 * t400;
t263 = qJDD(2) * pkin(3) - pkin(7) * t309 + t275;
t266 = pkin(7) * t424 + t276;
t438 = t403 * t263 + t481 * t266 + t285 * t446 - t287 * t455;
t437 = -t481 * t263 + t403 * t266 + t285 * t455 + t287 * t446;
t436 = g(1) * t407 + g(2) * t405;
t431 = pkin(1) + t490;
t430 = -0.2e1 * pkin(1) * t452 - pkin(6) * qJDD(2);
t423 = t382 * t481 - t448;
t422 = -g(1) * t464 - g(2) * t465 - g(3) * t380 + t438;
t335 = -qJDD(1) * t384 + t449;
t421 = t481 * t433;
t420 = g(1) * t466 + g(2) * t467 - g(3) * t381 - t437;
t319 = pkin(3) * t426 + t387;
t408 = qJD(2) ^ 2;
t417 = 0.2e1 * qJDD(1) * pkin(1) - pkin(6) * t408 + t485;
t409 = qJD(1) ^ 2;
t416 = pkin(1) * t409 - pkin(6) * qJDD(1) + t436;
t252 = qJD(4) * t273 - t281 * t481 + t403 * t282;
t415 = g(1) * t465 - g(2) * t464 - t252 * t396 + t394 * t429;
t414 = t264 * t396 - t422;
t322 = pkin(3) * t433 - t384;
t311 = t350 * t481 - t403 * t433;
t412 = -t396 * t460 + t420;
t411 = -t420 + t489;
t286 = -t447 * pkin(3) + (-pkin(1) + (-pkin(3) * t401 - pkin(2)) * t406) * qJDD(1) + t449;
t250 = -pkin(4) * t413 - qJ(5) * t427 + qJD(5) * t428 + t286;
t395 = -pkin(7) + t402;
t361 = qJ(5) * t464;
t360 = qJ(5) * t465;
t336 = -pkin(4) - t423;
t334 = qJ(5) + t458;
t310 = t350 * t403 + t421;
t278 = qJD(4) * t311 - t403 * t425 + t426 * t481;
t277 = t350 * t455 + t396 * t421 + t403 * t426;
t271 = t310 * pkin(4) - t311 * qJ(5) + t322;
t270 = t274 + t318;
t256 = t396 * qJ(5) + t265;
t255 = -t396 * pkin(4) + t453;
t253 = t278 * pkin(4) + t277 * qJ(5) - t311 * qJD(5) + t319;
t249 = t437 + t484;
t248 = t438 + t486;
t1 = [qJDD(1) * MDP(1) + t485 * MDP(2) + t436 * MDP(3) + (qJDD(1) * t398 + 0.2e1 * t404 * t444) * MDP(4) + 0.2e1 * (t404 * t450 - t452 * t457) * MDP(5) + (qJDD(2) * t404 + t406 * t408) * MDP(6) + (qJDD(2) * t406 - t404 * t408) * MDP(7) + (t404 * t430 + t406 * t417) * MDP(9) + (-t404 * t417 + t406 * t430) * MDP(10) + (t384 * t424 + t335 * t433 + t316 * qJDD(2) + t485 * t391 + (t341 * t479 + t350 * t359 + t290) * qJD(2)) * MDP(11) + (-t317 * qJDD(2) - t384 * t309 + t335 * t350 - t485 * t390 + (t343 * t479 - t359 * t433 - t291) * qJD(2)) * MDP(12) + (-t291 * t341 + t317 * t424 - t276 * t433 - t290 * t343 - t316 * t309 - t275 * t350 + (t307 * t433 - t308 * t350) * qJD(2) - t436) * MDP(13) + (t276 * t317 + t308 * t291 + t275 * t316 + t307 * t290 - t335 * t384 + t359 * t387 - g(1) * (-t384 * t405 - t402 * t407) - g(2) * (t384 * t407 - t402 * t405)) * MDP(14) + (t277 * t428 + t311 * t427) * MDP(15) + (-t277 * t299 + t278 * t428 - t310 * t427 + t311 * t413) * MDP(16) + (-t277 * t396 + t311 * t394) * MDP(17) + (-t278 * t396 - t310 * t394) * MDP(18) + (t278 * t315 + t286 * t310 - t299 * t319 - t322 * t413 + t415) * MDP(20) + (-t277 * t315 + t286 * t311 - t319 * t428 + t322 * t427 - t483) * MDP(21) + (t250 * t310 - t253 * t299 + t267 * t278 - t271 * t413 + t415) * MDP(22) + (-t248 * t310 + t249 * t311 + t251 * t299 - t252 * t428 - t255 * t277 - t256 * t278 + t273 * t413 - t427 * t429 - t436) * MDP(23) + (-t250 * t311 + t253 * t428 + t267 * t277 - t271 * t427 + t483) * MDP(24) + (t248 * t273 - t249 * t429 + t250 * t271 + t256 * t251 + t255 * t252 + t267 * t253 + (g(1) * t395 - g(2) * t431) * t407 + (g(1) * t431 + g(2) * t395) * t405) * MDP(25); MDP(6) * t451 + MDP(7) * t450 + qJDD(2) * MDP(8) + (t404 * t416 - t475) * MDP(9) + (g(3) * t404 + t406 * t416) * MDP(10) + (-g(3) * t391 - qJD(2) * t313 - t343 * t359 + t436 * t390 + (qJDD(2) * t401 - t341 * t456) * pkin(2) + t275) * MDP(11) + (g(3) * t390 + qJD(2) * t314 + t341 * t359 + t436 * t391 + (-qJDD(2) * t400 - t343 * t456) * pkin(2) - t276) * MDP(12) + ((t308 + t313) * t343 + (t314 - t307) * t341 + (-t401 * t309 + t400 * t424) * pkin(2)) * MDP(13) + (-t307 * t313 - t308 * t314 + (-t475 + t275 * t401 + t276 * t400 + (-qJD(1) * t359 + t436) * t404) * pkin(2)) * MDP(14) + (t299 * t318 + t394 * t423 + t412 + t491) * MDP(20) + (t318 * t428 - t458 * t394 + t487 * t396 - t422 - t492) * MDP(21) + (t270 * t299 - t336 * t394 + t412 - t489) * MDP(22) + (t413 * t334 + t336 * t427 + (-t256 - t460) * t428 + (-t255 + t461) * t299) * MDP(23) + (-t270 * t428 + t334 * t394 + t396 * t461 + t422 + t486 + t493) * MDP(24) + (t248 * t334 + t249 * t336 - t267 * t270 - g(1) * (t407 * t442 + t361) - g(2) * (t405 * t442 + t360) - g(3) * t490 + t461 * t256 + t460 * t255) * MDP(25) + (-t404 * t406 * MDP(4) + MDP(5) * t457) * t409 + t497; (qJD(2) * t343 - t424) * MDP(11) + ((qJD(1) * t462 - t341) * qJD(2) + t419) * MDP(12) + (-t341 ^ 2 - t343 ^ 2) * MDP(13) + (t307 * t343 + t308 * t341 + t335 - t485) * MDP(14) + (-t469 - t494) * MDP(23) + (t255 * t428 - t256 * t299 + t250 - t485) * MDP(25) + (MDP(21) - MDP(24)) * (t427 + t468) + (MDP(20) + MDP(22)) * (-t413 - t470); (t420 + t472 + t491) * MDP(20) + (t414 - t492) * MDP(21) + (t274 * t299 + t389 - t411 + t472) * MDP(22) + (-pkin(4) * t427 + qJ(5) * t413 - (t256 - t265) * t428 - (t255 - t453) * t299) * MDP(23) + (-t274 * t428 + 0.2e1 * t385 + 0.2e1 * t388 - t414 + t493) * MDP(24) + (t248 * qJ(5) - t249 * pkin(4) - t267 * t274 - t255 * t265 - g(1) * (-pkin(4) * t466 + t361) - g(2) * (-pkin(4) * t467 + t360) - g(3) * t459 + t453 * t256) * MDP(25) + t497; (-t394 + t496) * MDP(22) + t254 * MDP(23) + (-t396 ^ 2 - t494) * MDP(24) + (-t256 * t396 + t411) * MDP(25);];
tau = t1;
