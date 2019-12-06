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
% Datum: 2019-12-05 18:34
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
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
% StartTime: 2019-12-05 18:34:11
% EndTime: 2019-12-05 18:34:15
% DurationCPUTime: 2.54s
% Computational Cost: add. (2516->309), mult. (3782->381), div. (0->0), fcn. (2849->16), ass. (0->158)
t392 = cos(qJ(2));
t448 = qJD(1) * t392;
t440 = pkin(1) * t448;
t418 = qJD(3) - t440;
t382 = qJD(1) + qJD(2);
t385 = cos(pkin(9));
t391 = cos(qJ(4));
t457 = t391 * t385;
t437 = t382 * t457;
t384 = sin(pkin(9));
t387 = sin(qJ(4));
t460 = t384 * t387;
t438 = t382 * t460;
t311 = -t437 + t438;
t390 = cos(qJ(5));
t332 = t384 * t391 + t385 * t387;
t313 = t332 * t382;
t386 = sin(qJ(5));
t464 = t313 * t386;
t270 = t390 * t311 + t464;
t381 = qJD(4) + qJD(5);
t466 = t270 * t381;
t407 = t311 * t386 - t390 * t313;
t467 = t407 * t381;
t377 = qJDD(1) + qJDD(2);
t388 = sin(qJ(2));
t442 = qJDD(1) * t388;
t446 = qJD(2) * t392;
t299 = qJ(3) * t377 + qJD(3) * t382 + (qJD(1) * t446 + t442) * pkin(1);
t450 = t384 ^ 2 + t385 ^ 2;
t429 = t450 * t299;
t383 = qJ(1) + qJ(2);
t372 = sin(t383);
t373 = cos(t383);
t417 = g(2) * t373 + g(3) * t372;
t363 = g(2) * t372;
t477 = g(3) * t373 - t363;
t475 = pkin(1) * t388;
t441 = qJD(1) * t475;
t338 = qJ(3) * t382 + t441;
t431 = pkin(7) * t382 + t338;
t300 = t431 * t384;
t301 = t431 * t385;
t408 = t300 * t387 - t301 * t391;
t261 = -pkin(8) * t311 - t408;
t360 = -pkin(3) * t385 - pkin(2);
t309 = t360 * t382 + t418;
t281 = pkin(4) * t311 + t309;
t380 = pkin(9) + qJ(4);
t371 = qJ(5) + t380;
t357 = sin(t371);
t444 = qJD(5) * t386;
t358 = cos(t371);
t461 = t358 * t373;
t462 = t358 * t372;
t482 = g(1) * t357 - g(2) * t462 + g(3) * t461 + t261 * t444 + t281 * t270;
t435 = qJD(4) * t437 + t332 * t377;
t277 = -t438 * qJD(4) + t435;
t432 = pkin(7) * t377 + t299;
t285 = t432 * t384;
t286 = t432 * t385;
t426 = -t391 * t285 - t387 * t286;
t245 = qJDD(4) * pkin(4) - pkin(8) * t277 + qJD(4) * t408 + t426;
t322 = t332 * qJD(4);
t347 = t377 * t457;
t415 = -t377 * t460 + t347;
t278 = t322 * t382 - t415;
t409 = -t387 * t285 + t391 * t286;
t479 = -t391 * t300 - t301 * t387;
t246 = -pkin(8) * t278 + t479 * qJD(4) + t409;
t481 = -g(1) * t358 + t390 * t245 - t386 * t246 + t281 * t407 + t357 * t477;
t376 = qJDD(4) + qJDD(5);
t480 = t376 * MDP(22) + (-t270 ^ 2 + t407 ^ 2) * MDP(19) - t270 * MDP(18) * t407;
t343 = (-pkin(7) - qJ(3)) * t384;
t374 = t385 * pkin(7);
t344 = qJ(3) * t385 + t374;
t454 = t387 * t343 + t391 * t344;
t478 = -t454 * qJD(4) - t418 * t332;
t359 = qJ(3) + t475;
t323 = (-pkin(7) - t359) * t384;
t324 = t359 * t385 + t374;
t456 = t387 * t323 + t391 * t324;
t427 = t277 * t386 + t390 * t278;
t248 = -qJD(5) * t407 + t427;
t331 = -t457 + t460;
t445 = qJD(4) * t391;
t476 = -t331 * t440 + (qJD(3) * t384 + qJD(4) * t344) * t387 - qJD(3) * t457 - t343 * t445;
t474 = pkin(1) * t392;
t473 = pkin(2) * t377;
t472 = pkin(4) * t322;
t321 = t331 * qJD(4);
t471 = pkin(8) * t321;
t470 = pkin(8) * t332;
t260 = -pkin(8) * t313 + t479;
t259 = qJD(4) * pkin(4) + t260;
t469 = t259 * t390;
t468 = t261 * t390;
t453 = t417 * t385;
t447 = qJD(2) * t388;
t439 = pkin(1) * t447;
t451 = -qJD(1) * t439 + qJDD(1) * t474;
t443 = qJD(5) * t390;
t436 = t390 * t277 - t386 * t278 - t311 * t443;
t434 = t382 * t447;
t433 = qJDD(3) - t451;
t351 = pkin(1) * t446 + qJD(3);
t430 = t351 * t450;
t428 = t450 * t377;
t425 = t391 * t323 - t324 * t387;
t424 = t391 * t343 - t344 * t387;
t423 = t382 * t441;
t292 = -t331 * t386 + t332 * t390;
t263 = qJD(5) * t292 - t321 * t386 + t390 * t322;
t293 = t360 * t377 + t433;
t264 = pkin(4) * t278 + t293;
t291 = t390 * t331 + t332 * t386;
t422 = g(2) * t461 + g(3) * t462 + t281 * t263 + t264 * t291;
t370 = cos(t380);
t421 = t293 * t331 + t309 * t322 + t370 * t417;
t420 = -t477 + t429;
t419 = t451 + t417;
t282 = t424 - t470;
t319 = t322 * pkin(8);
t414 = -qJD(5) * t282 + t319 + t476;
t326 = t331 * pkin(8);
t283 = -t326 + t454;
t413 = qJD(5) * t283 - t471 - t478;
t412 = -t259 * t386 - t468;
t267 = t425 - t470;
t268 = -t326 + t456;
t411 = t267 * t390 - t268 * t386;
t410 = t267 * t386 + t268 * t390;
t308 = pkin(4) * t331 + t360;
t405 = -t441 + t472;
t403 = t423 + t473;
t247 = -t313 * t444 + t436;
t366 = -pkin(2) - t474;
t402 = -pkin(1) * t434 - t366 * t377;
t262 = -qJD(5) * t291 - t321 * t390 - t322 * t386;
t401 = t281 * t262 + t264 * t292 - t357 * t417;
t369 = sin(t380);
t400 = t293 * t332 - t309 * t321 - t369 * t417;
t399 = (-t247 * t291 - t248 * t292 - t262 * t270 + t263 * t407) * MDP(19) + (t247 * t292 - t262 * t407) * MDP(18) + (-t277 * t331 - t278 * t332 + t311 * t321 - t313 * t322) * MDP(12) + (t262 * t381 + t292 * t376) * MDP(20) + (-t263 * t381 - t291 * t376) * MDP(21) + (t277 * t332 - t313 * t321) * MDP(11) + (-qJD(4) * t321 + qJDD(4) * t332) * MDP(13) + (-qJD(4) * t322 - qJDD(4) * t331) * MDP(14) + t377 * MDP(4);
t398 = t323 * t445 + t351 * t457 + (-qJD(4) * t324 - t351 * t384) * t387;
t396 = t418 * t450;
t395 = -t456 * qJD(4) - t332 * t351;
t393 = cos(qJ(1));
t389 = sin(qJ(1));
t355 = t373 * qJ(3);
t342 = t360 - t474;
t333 = -pkin(2) * t382 + t418;
t310 = t433 - t473;
t303 = t439 + t472;
t302 = t310 * t384;
t297 = t308 - t474;
t257 = t395 + t471;
t256 = -t319 + t398;
t1 = [qJDD(1) * MDP(1) + ((-t310 + t402) * t385 + t453) * MDP(7) + t399 + (-t303 * t407 + t297 * t247 - (qJD(5) * t411 + t256 * t390 + t257 * t386) * t381 - t410 * t376 + t401) * MDP(24) + (-qJD(4) * t398 - qJDD(4) * t456 + t342 * t277 + t313 * t439 + t400) * MDP(17) + (t302 + (-t402 - t417) * t384) * MDP(8) + (t310 * t366 + t333 * t439 - g(2) * (-pkin(1) * t393 - t373 * pkin(2) - t372 * qJ(3)) - g(3) * (-pkin(1) * t389 - pkin(2) * t372 + t355) + t338 * t430 + t359 * t429) * MDP(10) + (qJD(4) * t395 + qJDD(4) * t425 + t342 * t278 + t311 * t439 + t421) * MDP(16) + (t359 * t428 + t382 * t430 + t420) * MDP(9) + (((-qJDD(1) - t377) * t388 + (-qJD(1) - t382) * t446) * pkin(1) + t477) * MDP(6) + ((t377 * t392 - t434) * pkin(1) + t419) * MDP(5) + (t303 * t270 + t297 * t248 + (-qJD(5) * t410 - t256 * t386 + t257 * t390) * t381 + t411 * t376 + t422) * MDP(23) + (g(2) * t393 + g(3) * t389) * MDP(2) + (-g(2) * t389 + g(3) * t393) * MDP(3); (t476 * qJD(4) - t454 * qJDD(4) + t360 * t277 - t313 * t441 + t400) * MDP(17) + t399 + (t308 * t248 + (t282 * t390 - t283 * t386) * t376 + (t386 * t414 - t390 * t413) * t381 + t405 * t270 + t422) * MDP(23) + (qJ(3) * t428 + t382 * t396 + t420) * MDP(9) + (-t333 * t441 - g(3) * t355 + (-t310 + t417) * pkin(2) + (t429 + t363) * qJ(3) + t396 * t338) * MDP(10) + ((-t442 + (-qJD(2) + t382) * t448) * pkin(1) + t477) * MDP(6) + (t302 + (-t403 - t417) * t384) * MDP(8) + (t478 * qJD(4) + t424 * qJDD(4) + t360 * t278 - t311 * t441 + t421) * MDP(16) + (t308 * t247 - (t282 * t386 + t283 * t390) * t376 + (t386 * t413 + t390 * t414) * t381 - t405 * t407 + t401) * MDP(24) + ((-t310 + t403) * t385 + t453) * MDP(7) + (t419 + t423) * MDP(5); (-t450 * t338 * t382 + qJDD(3) - t419) * MDP(10) - t347 * MDP(16) + t435 * MDP(17) + (t248 - t467) * MDP(23) + (t247 - t466) * MDP(24) + (-pkin(2) * MDP(10) - t385 * MDP(7) + (MDP(16) * t387 + MDP(8)) * t384) * t377 - t450 * MDP(9) * t382 ^ 2 + (0.2e1 * t313 * MDP(16) + (-t311 - t438) * MDP(17)) * qJD(4); t313 * t311 * MDP(11) + (-t311 ^ 2 + t313 ^ 2) * MDP(12) + ((t311 - t438) * qJD(4) + t435) * MDP(13) + t415 * MDP(14) + qJDD(4) * MDP(15) + (-g(1) * t370 - t309 * t313 + t369 * t477 + t426) * MDP(16) + (g(1) * t369 + t309 * t311 + t370 * t477 - t409) * MDP(17) + (t247 + t466) * MDP(20) + (-t248 - t467) * MDP(21) + (-(-t260 * t386 - t468) * t381 + t412 * qJD(5) + (-t270 * t313 + t390 * t376 - t381 * t444) * pkin(4) + t481) * MDP(23) + ((-t261 * t381 - t245) * t386 + (-qJD(5) * t259 + t260 * t381 - t246) * t390 + (t313 * t407 - t386 * t376 - t381 * t443) * pkin(4) + t482) * MDP(24) + t480; (t436 + t466) * MDP(20) + (-t427 - t467) * MDP(21) + (-t381 * t412 + t481) * MDP(23) + (-t390 * t246 - t386 * t245 + (-t261 * t386 + t469) * t381 + t482) * MDP(24) + (-MDP(20) * t464 + t407 * MDP(21) + t412 * MDP(23) - MDP(24) * t469) * qJD(5) + t480;];
tau = t1;
