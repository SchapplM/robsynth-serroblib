% Calculate vector of inverse dynamics joint torques for
% S5RRPRR5
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
%   pkin=[a2,a3,a4,a5,d1,d2,d4,d5,theta3]';
% MDP [24x1]
%   Minimal dynamic parameter vector (fixed base model)
%   see S5RRPRR5_convert_par2_MPV_fixb.m
% 
% Output:
% tau [5x1]
%   joint torques of inverse dynamics (contains inertial, gravitational Coriolis and centrifugal forces)

% Quelle: HybrDyn-Toolbox
% Datum: 2020-01-03 12:04
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tau = S5RRPRR5_invdynJ_fixb_mdp_slag_vp(qJ, qJD, qJDD, g, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(5,1),zeros(3,1),zeros(9,1),zeros(24,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRPRR5_invdynJ_fixb_mdp_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RRPRR5_invdynJ_fixb_mdp_slag_vp: qJD has to be [5x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [5 1]), ...
  'S5RRPRR5_invdynJ_fixb_mdp_slag_vp: qJDD has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RRPRR5_invdynJ_fixb_mdp_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5RRPRR5_invdynJ_fixb_mdp_slag_vp: pkin has to be [9x1] (double)');
assert(isreal(MDP) && all(size(MDP) == [24 1]), ...
  'S5RRPRR5_invdynJ_fixb_mdp_slag_vp: MDP has to be [24x1] (double)'); 

%% Symbolic Calculation
% From invdyn_joint_fixb_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2020-01-03 12:04:04
% EndTime: 2020-01-03 12:04:09
% DurationCPUTime: 2.82s
% Computational Cost: add. (2516->308), mult. (3782->379), div. (0->0), fcn. (2849->16), ass. (0->160)
t391 = sin(qJ(2));
t451 = qJD(2) * t391;
t443 = pkin(1) * t451;
t395 = cos(qJ(2));
t481 = pkin(1) * t395;
t455 = -qJD(1) * t443 + qJDD(1) * t481;
t436 = qJDD(3) - t455;
t380 = qJDD(1) + qJDD(2);
t480 = pkin(2) * t380;
t313 = t436 - t480;
t386 = qJ(1) + qJ(2);
t375 = sin(t386);
t376 = cos(t386);
t490 = g(2) * t376 + g(3) * t375;
t491 = t313 + t490;
t452 = qJD(1) * t395;
t444 = pkin(1) * t452;
t422 = qJD(3) - t444;
t385 = qJD(1) + qJD(2);
t388 = cos(pkin(9));
t394 = cos(qJ(4));
t461 = t394 * t388;
t441 = t385 * t461;
t387 = sin(pkin(9));
t390 = sin(qJ(4));
t465 = t387 * t390;
t442 = t385 * t465;
t314 = -t441 + t442;
t393 = cos(qJ(5));
t335 = t387 * t394 + t388 * t390;
t316 = t335 * t385;
t389 = sin(qJ(5));
t469 = t316 * t389;
t273 = t393 * t314 + t469;
t384 = qJD(4) + qJD(5);
t471 = t273 * t384;
t411 = t314 * t389 - t393 * t316;
t472 = t411 * t384;
t446 = qJDD(1) * t391;
t450 = qJD(2) * t395;
t302 = qJ(3) * t380 + qJD(3) * t385 + (qJD(1) * t450 + t446) * pkin(1);
t454 = t387 ^ 2 + t388 ^ 2;
t432 = t454 * t302;
t366 = g(3) * t376;
t484 = g(2) * t375 - t366;
t482 = pkin(1) * t391;
t445 = qJD(1) * t482;
t341 = qJ(3) * t385 + t445;
t434 = pkin(7) * t385 + t341;
t303 = t434 * t387;
t304 = t434 * t388;
t412 = t303 * t390 - t304 * t394;
t264 = -pkin(8) * t314 - t412;
t363 = -pkin(3) * t388 - pkin(2);
t312 = t363 * t385 + t422;
t284 = pkin(4) * t314 + t312;
t383 = pkin(9) + qJ(4);
t374 = qJ(5) + t383;
t360 = sin(t374);
t361 = cos(t374);
t448 = qJD(5) * t389;
t489 = g(1) * t360 + t264 * t448 + t284 * t273 + t361 * t484;
t438 = qJD(4) * t441 + t335 * t380;
t280 = -t442 * qJD(4) + t438;
t435 = pkin(7) * t380 + t302;
t288 = t435 * t387;
t289 = t435 * t388;
t429 = -t394 * t288 - t390 * t289;
t248 = qJDD(4) * pkin(4) - pkin(8) * t280 + qJD(4) * t412 + t429;
t325 = t335 * qJD(4);
t350 = t380 * t461;
t419 = -t380 * t465 + t350;
t281 = t325 * t385 - t419;
t413 = -t390 * t288 + t394 * t289;
t486 = -t394 * t303 - t304 * t390;
t249 = -pkin(8) * t281 + t486 * qJD(4) + t413;
t466 = t360 * t376;
t467 = t360 * t375;
t488 = -g(1) * t361 + g(2) * t467 - g(3) * t466 + t393 * t248 - t389 * t249 + t284 * t411;
t379 = qJDD(4) + qJDD(5);
t487 = t379 * MDP(22) + (-t273 ^ 2 + t411 ^ 2) * MDP(19) - t273 * MDP(18) * t411;
t346 = (-pkin(7) - qJ(3)) * t387;
t377 = t388 * pkin(7);
t347 = qJ(3) * t388 + t377;
t458 = t390 * t346 + t394 * t347;
t485 = -t458 * qJD(4) - t422 * t335;
t362 = qJ(3) + t482;
t326 = (-pkin(7) - t362) * t387;
t327 = t362 * t388 + t377;
t460 = t390 * t326 + t394 * t327;
t430 = t280 * t389 + t393 * t281;
t251 = -t411 * qJD(5) + t430;
t334 = -t461 + t465;
t449 = qJD(4) * t394;
t483 = -t334 * t444 + (qJD(3) * t387 + qJD(4) * t347) * t390 - qJD(3) * t461 - t346 * t449;
t479 = pkin(4) * t325;
t324 = t334 * qJD(4);
t478 = pkin(8) * t324;
t477 = pkin(8) * t335;
t263 = -pkin(8) * t316 + t486;
t262 = qJD(4) * pkin(4) + t263;
t474 = t262 * t393;
t473 = t264 * t393;
t463 = t388 * MDP(7);
t457 = t376 * pkin(2) + t375 * qJ(3);
t447 = qJD(5) * t393;
t440 = t393 * t280 - t389 * t281 - t314 * t447;
t439 = t491 * t387;
t437 = t385 * t451;
t354 = pkin(1) * t450 + qJD(3);
t433 = t354 * t454;
t431 = t454 * t380;
t428 = t394 * t326 - t327 * t390;
t427 = t394 * t346 - t347 * t390;
t426 = t385 * t445;
t294 = t393 * t334 + t335 * t389;
t265 = -qJD(5) * t294 - t324 * t393 - t325 * t389;
t296 = t363 * t380 + t436;
t267 = pkin(4) * t281 + t296;
t295 = -t334 * t389 + t335 * t393;
t425 = g(2) * t466 + g(3) * t467 + t284 * t265 + t267 * t295;
t372 = sin(t383);
t424 = t296 * t335 - t312 * t324 + t372 * t490;
t423 = -t484 + t432;
t285 = t427 - t477;
t322 = t325 * pkin(8);
t418 = -qJD(5) * t285 + t322 + t483;
t329 = t334 * pkin(8);
t286 = -t329 + t458;
t417 = qJD(5) * t286 - t478 - t485;
t416 = -t262 * t389 - t473;
t270 = t428 - t477;
t271 = -t329 + t460;
t415 = t270 * t393 - t271 * t389;
t414 = t270 * t389 + t271 * t393;
t311 = pkin(4) * t334 + t363;
t409 = -t445 + t479;
t407 = -t490 + t455;
t250 = -t316 * t448 + t440;
t369 = -pkin(2) - t481;
t406 = pkin(1) * t437 + t369 * t380;
t266 = qJD(5) * t295 - t324 * t389 + t393 * t325;
t405 = t284 * t266 + t267 * t294 - t361 * t490;
t373 = cos(t383);
t404 = t296 * t334 + t312 * t325 - t373 * t490;
t403 = (-t250 * t294 - t251 * t295 - t265 * t273 + t266 * t411) * MDP(19) + (t250 * t295 - t265 * t411) * MDP(18) + (-t280 * t334 - t281 * t335 + t314 * t324 - t316 * t325) * MDP(12) + (t265 * t384 + t295 * t379) * MDP(20) + (-t266 * t384 - t294 * t379) * MDP(21) + (t280 * t335 - t316 * t324) * MDP(11) + (-qJD(4) * t324 + qJDD(4) * t335) * MDP(13) + (-qJD(4) * t325 - qJDD(4) * t334) * MDP(14) + t380 * MDP(4);
t402 = -t490 + t426;
t401 = t326 * t449 + t354 * t461 + (-qJD(4) * t327 - t354 * t387) * t390;
t399 = t422 * t454;
t398 = -t460 * qJD(4) - t335 * t354;
t396 = cos(qJ(1));
t392 = sin(qJ(1));
t364 = t375 * pkin(2);
t345 = t363 - t481;
t336 = -pkin(2) * t385 + t422;
t306 = t443 + t479;
t300 = t311 - t481;
t260 = t398 + t478;
t259 = -t322 + t401;
t1 = [t403 + qJDD(1) * MDP(1) + (((-qJDD(1) - t380) * t391 + (-qJD(1) - t385) * t450) * pkin(1) + t484) * MDP(6) + ((t380 * t395 - t437) * pkin(1) + t407) * MDP(5) + (-t306 * t411 + t300 * t250 - (qJD(5) * t415 + t259 * t393 + t260 * t389) * t384 - t414 * t379 + t425) * MDP(24) + (qJD(4) * t398 + qJDD(4) * t428 + t345 * t281 + t314 * t443 + t404) * MDP(16) + (t313 * t369 + t336 * t443 - g(2) * (pkin(1) * t396 + t457) - g(3) * (pkin(1) * t392 - qJ(3) * t376 + t364) + t341 * t433 + t362 * t432) * MDP(10) + (-t406 - t491) * t463 + (t387 * t406 + t439) * MDP(8) + (t306 * t273 + t300 * t251 + (-qJD(5) * t414 - t259 * t389 + t260 * t393) * t384 + t415 * t379 + t405) * MDP(23) + (t362 * t431 + t385 * t433 + t423) * MDP(9) + (-qJD(4) * t401 - qJDD(4) * t460 + t345 * t280 + t316 * t443 + t424) * MDP(17) + (-g(2) * t396 - g(3) * t392) * MDP(2) + (g(2) * t392 - g(3) * t396) * MDP(3); t403 + (t483 * qJD(4) - t458 * qJDD(4) + t363 * t280 - t316 * t445 + t424) * MDP(17) + ((-t446 + (-qJD(2) + t385) * t452) * pkin(1) + t484) * MDP(6) + (t311 * t251 + (t285 * t393 - t286 * t389) * t379 + (t389 * t418 - t393 * t417) * t384 + t409 * t273 + t405) * MDP(23) + (qJ(3) * t431 + t385 * t399 + t423) * MDP(9) + (-t313 * pkin(2) - t336 * t445 - g(2) * t457 - g(3) * t364 + (t432 + t366) * qJ(3) + t399 * t341) * MDP(10) + (t485 * qJD(4) + t427 * qJDD(4) + t363 * t281 - t314 * t445 + t404) * MDP(16) + ((-t426 - t480) * t387 + t439) * MDP(8) + (t402 + t455) * MDP(5) + (t311 * t250 - (t285 * t389 + t286 * t393) * t379 + (t389 * t417 + t393 * t418) * t384 - t409 * t411 + t425) * MDP(24) + (-t313 + t402 + t480) * t463; (-t454 * t341 * t385 + qJDD(3) - t407) * MDP(10) - t350 * MDP(16) + t438 * MDP(17) + (t251 - t472) * MDP(23) + (t250 - t471) * MDP(24) + (-pkin(2) * MDP(10) - t463 + (MDP(16) * t390 + MDP(8)) * t387) * t380 - t454 * MDP(9) * t385 ^ 2 + (0.2e1 * t316 * MDP(16) + (-t314 - t442) * MDP(17)) * qJD(4); t316 * t314 * MDP(11) + (-t314 ^ 2 + t316 ^ 2) * MDP(12) + ((t314 - t442) * qJD(4) + t438) * MDP(13) + t419 * MDP(14) + qJDD(4) * MDP(15) + (-g(1) * t373 - t312 * t316 + t372 * t484 + t429) * MDP(16) + (g(1) * t372 + t312 * t314 + t373 * t484 - t413) * MDP(17) + (t250 + t471) * MDP(20) + (-t251 - t472) * MDP(21) + (-(-t263 * t389 - t473) * t384 + t416 * qJD(5) + (-t273 * t316 + t393 * t379 - t384 * t448) * pkin(4) + t488) * MDP(23) + ((-t264 * t384 - t248) * t389 + (-qJD(5) * t262 + t263 * t384 - t249) * t393 + (t316 * t411 - t389 * t379 - t384 * t447) * pkin(4) + t489) * MDP(24) + t487; (t440 + t471) * MDP(20) + (-t430 - t472) * MDP(21) + (-t384 * t416 + t488) * MDP(23) + (-t393 * t249 - t389 * t248 + (-t264 * t389 + t474) * t384 + t489) * MDP(24) + (-MDP(20) * t469 + MDP(21) * t411 + MDP(23) * t416 - MDP(24) * t474) * qJD(5) + t487;];
tau = t1;
