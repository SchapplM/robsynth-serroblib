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
% MDP [23x1]
%   Minimal dynamic parameter vector (fixed base model)
%   see S5RRPRR5_convert_par2_MPV_fixb.m
% 
% Output:
% tau [5x1]
%   joint torques of inverse dynamics (contains inertial, gravitational Coriolis and centrifugal forces)

% Quelle: HybrDyn-Toolbox
% Datum: 2022-01-20 11:03
% Revision: fd3771346c4aea32fdeb66112c511235427c26a7 (2022-01-20)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tau = S5RRPRR5_invdynJ_fixb_mdp_slag_vp(qJ, qJD, qJDD, g, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(5,1),zeros(3,1),zeros(9,1),zeros(23,1)}
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
assert(isreal(MDP) && all(size(MDP) == [23 1]), ...
  'S5RRPRR5_invdynJ_fixb_mdp_slag_vp: MDP has to be [23x1] (double)'); 

%% Symbolic Calculation
% From invdyn_joint_fixb_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2022-01-20 11:03:01
% EndTime: 2022-01-20 11:03:07
% DurationCPUTime: 2.56s
% Computational Cost: add. (2487->303), mult. (3743->377), div. (0->0), fcn. (2827->16), ass. (0->159)
t396 = cos(qJ(2));
t450 = qJD(1) * t396;
t442 = pkin(1) * t450;
t421 = qJD(3) - t442;
t386 = qJD(1) + qJD(2);
t389 = cos(pkin(9));
t395 = cos(qJ(4));
t459 = t395 * t389;
t439 = t386 * t459;
t388 = sin(pkin(9));
t391 = sin(qJ(4));
t462 = t388 * t391;
t440 = t386 * t462;
t313 = -t439 + t440;
t394 = cos(qJ(5));
t334 = t388 * t395 + t389 * t391;
t315 = t334 * t386;
t390 = sin(qJ(5));
t467 = t315 * t390;
t273 = t313 * t394 + t467;
t385 = qJD(4) + qJD(5);
t469 = t273 * t385;
t411 = t313 * t390 - t315 * t394;
t470 = t411 * t385;
t381 = qJDD(1) + qJDD(2);
t392 = sin(qJ(2));
t444 = qJDD(1) * t392;
t448 = qJD(2) * t396;
t302 = qJ(3) * t381 + qJD(3) * t386 + (qJD(1) * t448 + t444) * pkin(1);
t452 = t388 ^ 2 + t389 ^ 2;
t430 = t452 * t302;
t387 = qJ(1) + qJ(2);
t376 = sin(t387);
t367 = g(1) * t376;
t377 = cos(t387);
t473 = g(2) * t377;
t488 = t367 - t473;
t333 = -t459 + t462;
t479 = pkin(1) * t392;
t443 = qJD(1) * t479;
t342 = qJ(3) * t386 + t443;
t432 = pkin(7) * t386 + t342;
t303 = t432 * t388;
t304 = t432 * t389;
t412 = t303 * t391 - t304 * t395;
t264 = -pkin(8) * t313 - t412;
t364 = -pkin(3) * t389 - pkin(2);
t311 = t364 * t386 + t421;
t284 = pkin(4) * t313 + t311;
t384 = pkin(9) + qJ(4);
t375 = qJ(5) + t384;
t361 = sin(t375);
t446 = qJD(5) * t390;
t362 = cos(t375);
t463 = t377 * t362;
t465 = t362 * t376;
t487 = g(1) * t463 + g(2) * t465 + g(3) * t361 + t264 * t446 + t284 * t273;
t437 = qJD(4) * t439 + t334 * t381;
t280 = -qJD(4) * t440 + t437;
t433 = pkin(7) * t381 + t302;
t288 = t433 * t388;
t289 = t433 * t389;
t426 = -t288 * t395 - t391 * t289;
t248 = qJDD(4) * pkin(4) - pkin(8) * t280 + qJD(4) * t412 + t426;
t324 = t334 * qJD(4);
t419 = t333 * t381;
t281 = t324 * t386 + t419;
t413 = -t391 * t288 + t395 * t289;
t484 = -t303 * t395 - t304 * t391;
t249 = -pkin(8) * t281 + qJD(4) * t484 + t413;
t464 = t376 * t361;
t466 = t361 * t377;
t486 = g(1) * t466 + g(2) * t464 - g(3) * t362 + t248 * t394 - t390 * t249 + t284 * t411;
t380 = qJDD(4) + qJDD(5);
t485 = t380 * MDP(21) + (-t273 ^ 2 + t411 ^ 2) * MDP(18) - t273 * MDP(17) * t411;
t347 = (-pkin(7) - qJ(3)) * t388;
t378 = t389 * pkin(7);
t348 = qJ(3) * t389 + t378;
t456 = t347 * t391 + t348 * t395;
t483 = -t456 * qJD(4) - t334 * t421;
t363 = qJ(3) + t479;
t325 = (-pkin(7) - t363) * t388;
t326 = t363 * t389 + t378;
t458 = t325 * t391 + t326 * t395;
t482 = g(1) * t377 + g(2) * t376;
t449 = qJD(2) * t392;
t441 = pkin(1) * t449;
t478 = pkin(1) * t396;
t453 = -qJD(1) * t441 + qJDD(1) * t478;
t481 = t367 + t453;
t427 = t280 * t390 + t281 * t394;
t251 = -qJD(5) * t411 + t427;
t447 = qJD(4) * t395;
t480 = -t333 * t442 + (qJD(3) * t388 + qJD(4) * t348) * t391 - qJD(3) * t459 - t347 * t447;
t477 = pkin(2) * t381;
t476 = pkin(4) * t324;
t323 = t333 * qJD(4);
t475 = pkin(8) * t323;
t474 = pkin(8) * t334;
t263 = -pkin(8) * t315 + t484;
t262 = qJD(4) * pkin(4) + t263;
t472 = t262 * t394;
t471 = t264 * t394;
t455 = pkin(2) * t377 + qJ(3) * t376;
t445 = qJD(5) * t394;
t438 = t280 * t394 - t281 * t390 - t313 * t445;
t436 = t386 * t449;
t435 = qJDD(3) - t453;
t434 = -pkin(2) * t376 + qJ(3) * t377;
t312 = t435 - t477;
t431 = -t312 - t473;
t354 = pkin(1) * t448 + qJD(3);
t429 = t452 * t354;
t428 = t452 * t381;
t425 = t325 * t395 - t326 * t391;
t424 = t347 * t395 - t348 * t391;
t422 = -t482 + t430;
t285 = t424 - t474;
t321 = t324 * pkin(8);
t418 = -qJD(5) * t285 + t321 + t480;
t328 = t333 * pkin(8);
t286 = -t328 + t456;
t417 = qJD(5) * t286 - t475 - t483;
t416 = -t262 * t390 - t471;
t270 = t425 - t474;
t271 = -t328 + t458;
t415 = t270 * t394 - t271 * t390;
t414 = t270 * t390 + t271 * t394;
t294 = t333 * t394 + t334 * t390;
t295 = -t333 * t390 + t334 * t394;
t310 = pkin(4) * t333 + t364;
t409 = -t443 + t476;
t265 = -qJD(5) * t294 - t323 * t394 - t324 * t390;
t296 = t364 * t381 + t435;
t267 = pkin(4) * t281 + t296;
t408 = -g(1) * t464 + g(2) * t466 + t265 * t284 + t267 * t295;
t373 = sin(t384);
t407 = t296 * t334 - t311 * t323 - t373 * t488;
t266 = qJD(5) * t295 - t323 * t390 + t324 * t394;
t406 = g(1) * t465 - g(2) * t463 + t266 * t284 + t267 * t294;
t374 = cos(t384);
t405 = t296 * t333 + t311 * t324 + t374 * t488;
t404 = t386 * t443 - t473;
t250 = -t315 * t446 + t438;
t402 = (-t250 * t294 - t251 * t295 - t265 * t273 + t266 * t411) * MDP(18) + (t250 * t295 - t265 * t411) * MDP(17) + (-t280 * t333 - t281 * t334 + t313 * t323 - t315 * t324) * MDP(11) + (t265 * t385 + t295 * t380) * MDP(19) + (-t266 * t385 - t294 * t380) * MDP(20) + (t280 * t334 - t315 * t323) * MDP(10) + (-qJD(4) * t323 + qJDD(4) * t334) * MDP(12) + (-qJD(4) * t324 - qJDD(4) * t333) * MDP(13) + t381 * MDP(4);
t401 = t325 * t447 + t354 * t459 + (-qJD(4) * t326 - t354 * t388) * t391;
t400 = t421 * t452;
t399 = -qJD(4) * t458 - t334 * t354;
t397 = cos(qJ(1));
t393 = sin(qJ(1));
t370 = -pkin(2) - t478;
t353 = t389 * t367;
t346 = t364 - t478;
t335 = -pkin(2) * t386 + t421;
t305 = t441 + t476;
t300 = t310 - t478;
t260 = t399 + t475;
t259 = -t321 + t401;
t1 = [(-t305 * t411 + t300 * t250 - (qJD(5) * t415 + t259 * t394 + t260 * t390) * t385 - t414 * t380 + t408) * MDP(23) + (((-qJDD(1) - t381) * t392 + (-qJD(1) - t386) * t448) * pkin(1) + t482) * MDP(6) + (-t473 + (t381 * t396 - t436) * pkin(1) + t481) * MDP(5) + t402 + (t363 * t428 + t386 * t429 + t422) * MDP(8) + (t312 * t370 + t335 * t441 - g(1) * (-pkin(1) * t393 + t434) - g(2) * (pkin(1) * t397 + t455) + t342 * t429 + t363 * t430) * MDP(9) + (g(1) * t393 - g(2) * t397) * MDP(2) + (g(1) * t397 + g(2) * t393) * MDP(3) + (t353 + (-pkin(1) * t436 - t370 * t381 + t431) * t389) * MDP(7) + qJDD(1) * MDP(1) + (t305 * t273 + t300 * t251 + (-qJD(5) * t414 - t259 * t390 + t260 * t394) * t385 + t415 * t380 + t406) * MDP(22) + (-qJD(4) * t401 - qJDD(4) * t458 + t346 * t280 + t315 * t441 + t407) * MDP(16) + (qJD(4) * t399 + qJDD(4) * t425 + t346 * t281 + t313 * t441 + t405) * MDP(15); t402 + (t310 * t250 - (t285 * t390 + t286 * t394) * t380 + (t390 * t417 + t394 * t418) * t385 - t409 * t411 + t408) * MDP(23) + (t404 + t481) * MDP(5) + ((-t444 + (-qJD(2) + t386) * t450) * pkin(1) + t482) * MDP(6) + (qJD(4) * t480 - t456 * qJDD(4) + t364 * t280 - t315 * t443 + t407) * MDP(16) + (t353 + (-t312 + t404 + t477) * t389) * MDP(7) + (t310 * t251 + (t285 * t394 - t286 * t390) * t380 + (t390 * t418 - t394 * t417) * t385 + t409 * t273 + t406) * MDP(22) + (qJ(3) * t428 + t386 * t400 + t422) * MDP(8) + (-t312 * pkin(2) - g(1) * t434 - g(2) * t455 + qJ(3) * t430 - t335 * t443 + t342 * t400) * MDP(9) + (qJD(4) * t483 + t424 * qJDD(4) + t364 * t281 - t313 * t443 + t405) * MDP(15); -t389 * t381 * MDP(7) + (-t342 * t386 * t452 - t367 - t431) * MDP(9) + t419 * MDP(15) + t437 * MDP(16) + (t251 - t470) * MDP(22) + (t250 - t469) * MDP(23) - t452 * MDP(8) * t386 ^ 2 + (0.2e1 * t315 * MDP(15) + (-t313 - t440) * MDP(16)) * qJD(4); t315 * t313 * MDP(10) + (-t313 ^ 2 + t315 ^ 2) * MDP(11) + ((t313 - t440) * qJD(4) + t437) * MDP(12) - t419 * MDP(13) + qJDD(4) * MDP(14) + (-g(3) * t374 - t311 * t315 + t373 * t482 + t426) * MDP(15) + (g(3) * t373 + t311 * t313 + t374 * t482 - t413) * MDP(16) + (t250 + t469) * MDP(19) + (-t251 - t470) * MDP(20) + (-(-t263 * t390 - t471) * t385 + t416 * qJD(5) + (-t273 * t315 + t380 * t394 - t385 * t446) * pkin(4) + t486) * MDP(22) + ((-t264 * t385 - t248) * t390 + (-qJD(5) * t262 + t263 * t385 - t249) * t394 + (t315 * t411 - t380 * t390 - t385 * t445) * pkin(4) + t487) * MDP(23) + t485; (t438 + t469) * MDP(19) + (-t427 - t470) * MDP(20) + (-t385 * t416 + t486) * MDP(22) + (-t394 * t249 - t390 * t248 + (-t264 * t390 + t472) * t385 + t487) * MDP(23) + (-MDP(19) * t467 + MDP(20) * t411 + t416 * MDP(22) - MDP(23) * t472) * qJD(5) + t485;];
tau = t1;
