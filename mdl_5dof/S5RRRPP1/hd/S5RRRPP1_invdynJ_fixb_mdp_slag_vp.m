% Calculate vector of inverse dynamics joint torques for
% S5RRRPP1
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
%   pkin=[a2,a3,a4,a5,d1,d2,d3,theta4]';
% MDP [21x1]
%   Minimal dynamic parameter vector (fixed base model)
%   see S5RRRPP1_convert_par2_MPV_fixb.m
% 
% Output:
% tau [5x1]
%   joint torques of inverse dynamics (contains inertial, gravitational Coriolis and centrifugal forces)

% Quelle: HybrDyn-Toolbox
% Datum: 2021-01-15 22:15
% Revision: 24b2e7d74a0c1a3b64fa2f8f5ad758691ad61af3 (2021-01-15)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tau = S5RRRPP1_invdynJ_fixb_mdp_slag_vp(qJ, qJD, qJDD, g, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(5,1),zeros(3,1),zeros(8,1),zeros(21,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRRPP1_invdynJ_fixb_mdp_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RRRPP1_invdynJ_fixb_mdp_slag_vp: qJD has to be [5x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [5 1]), ...
  'S5RRRPP1_invdynJ_fixb_mdp_slag_vp: qJDD has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RRRPP1_invdynJ_fixb_mdp_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RRRPP1_invdynJ_fixb_mdp_slag_vp: pkin has to be [8x1] (double)');
assert(isreal(MDP) && all(size(MDP) == [21 1]), ...
  'S5RRRPP1_invdynJ_fixb_mdp_slag_vp: MDP has to be [21x1] (double)'); 

%% Symbolic Calculation
% From invdyn_joint_fixb_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2021-01-15 22:14:56
% EndTime: 2021-01-15 22:15:04
% DurationCPUTime: 3.34s
% Computational Cost: add. (2972->355), mult. (4433->399), div. (0->0), fcn. (2695->12), ass. (0->171)
t386 = cos(pkin(8));
t391 = cos(qJ(3));
t462 = t386 * t391;
t385 = sin(pkin(8));
t388 = sin(qJ(3));
t463 = t385 * t388;
t330 = -t462 + t463;
t491 = qJ(4) + pkin(7);
t472 = pkin(3) * t391;
t364 = pkin(2) + t472;
t378 = qJDD(1) + qJDD(2);
t420 = t364 * t378;
t372 = t391 * qJD(4);
t431 = qJD(3) * t491;
t324 = -t388 * t431 + t372;
t408 = -qJD(4) * t388 - t391 * t431;
t392 = cos(qJ(2));
t452 = qJD(1) * t392;
t442 = pkin(1) * t452;
t458 = t386 * t324 + t330 * t442 + t385 * t408;
t384 = qJ(1) + qJ(2);
t374 = cos(t384);
t360 = g(2) * t374;
t389 = sin(qJ(2));
t451 = qJD(2) * t389;
t369 = pkin(1) * t451;
t476 = pkin(1) * t392;
t454 = -qJD(1) * t369 + qJDD(1) * t476;
t475 = pkin(2) * t378;
t490 = -t454 - t475 + t360;
t489 = t330 * t378;
t331 = t385 * t391 + t386 * t388;
t379 = qJD(1) + qJD(2);
t319 = t331 * t379;
t325 = t331 * qJD(3);
t287 = t379 * t325 + t489;
t448 = qJD(3) * t388;
t434 = t385 * t448;
t406 = t331 * t378 - t379 * t434;
t447 = qJD(3) * t391;
t435 = t386 * t447;
t288 = t379 * t435 + t406;
t436 = t379 * t448;
t424 = pkin(3) * t436 + qJDD(4) - t454;
t397 = pkin(4) * t287 - qJ(5) * t288 - qJD(5) * t319 + t424;
t251 = -t420 + t397;
t316 = -t364 * t379 + qJD(4) - t442;
t440 = t379 * t462;
t317 = t379 * t463 - t440;
t274 = pkin(4) * t317 - qJ(5) * t319 + t316;
t488 = t251 * t330 + t274 * t325;
t326 = -t434 + t435;
t487 = -t251 * t331 - t274 * t326;
t292 = -t420 + t424;
t486 = t292 * t330 + t316 * t325;
t485 = t292 * t331 + t316 * t326;
t373 = sin(t384);
t483 = g(1) * t374 + g(2) * t373;
t361 = g(1) * t373;
t482 = t361 - t360;
t445 = qJDD(1) * t389;
t450 = qJD(2) * t392;
t322 = pkin(7) * t378 + (qJD(1) * t450 + t445) * pkin(1);
t413 = qJ(4) * t378 + qJD(4) * t379 + t322;
t478 = pkin(1) * t389;
t443 = qJD(1) * t478;
t429 = t491 * t379 + t443;
t414 = qJD(3) * t429;
t269 = qJDD(3) * pkin(3) - t413 * t388 - t391 * t414;
t273 = -t388 * t414 + t413 * t391;
t255 = t385 * t269 + t386 * t273;
t380 = qJ(3) + pkin(8);
t370 = sin(t380);
t371 = cos(t380);
t481 = g(3) * t370 + t371 * t483 - t255;
t381 = qJDD(3) * qJ(5);
t252 = qJD(3) * qJD(5) + t255 + t381;
t254 = t386 * t269 - t385 * t273;
t470 = qJDD(3) * pkin(4);
t253 = qJDD(5) - t254 - t470;
t307 = t429 * t388;
t306 = qJD(3) * pkin(3) - t307;
t308 = t429 * t391;
t469 = t308 * t385;
t279 = t306 * t386 - t469;
t277 = -qJD(3) * pkin(4) + qJD(5) - t279;
t302 = t386 * t308;
t280 = t385 * t306 + t302;
t278 = qJD(3) * qJ(5) + t280;
t480 = -t252 * t330 + t253 * t331 + t277 * t326 - t278 * t325;
t479 = -t254 * t331 - t255 * t330 - t279 * t326 - t280 * t325;
t313 = t319 ^ 2;
t390 = sin(qJ(1));
t477 = pkin(1) * t390;
t474 = pkin(2) * t379;
t473 = pkin(3) * t388;
t471 = g(3) * t391;
t468 = t370 * t373;
t467 = t370 * t374;
t466 = t371 * t374;
t465 = t374 * t491;
t464 = t379 * t388;
t461 = t388 * t391;
t363 = pkin(7) + t478;
t460 = -qJ(4) - t363;
t459 = t324 * t385 - t331 * t442 - t386 * t408;
t337 = -t442 - t474;
t457 = t337 * t448 + t391 * t361;
t456 = g(1) * t468 - g(2) * t467;
t382 = t388 ^ 2;
t453 = -t391 ^ 2 + t382;
t282 = -t307 * t385 + t302;
t449 = qJD(3) * t282;
t283 = -t307 * t386 - t469;
t446 = qJD(5) - t283;
t441 = pkin(1) * t450;
t368 = pkin(3) * t448;
t439 = t337 * t447 + t490 * t388;
t346 = t374 * t364;
t438 = pkin(4) * t466 + qJ(5) * t467 + t346;
t437 = t379 * t451;
t433 = t491 * t388;
t430 = t460 * t388;
t428 = t373 * t491 + t346;
t427 = qJD(3) * t460;
t426 = t379 * t443;
t425 = t454 + t482;
t423 = -g(2) * t466 + t371 * t361;
t281 = pkin(4) * t325 - qJ(5) * t326 - qJD(5) * t331 + t368;
t421 = -t281 + t443;
t418 = -pkin(4) * t371 - qJ(5) * t370;
t375 = t391 * qJ(4);
t350 = pkin(7) * t391 + t375;
t305 = t386 * t350 - t385 * t433;
t417 = qJDD(3) * t305 + t456;
t394 = qJD(3) ^ 2;
t416 = 0.2e1 * (-t453 * t379 * qJD(3) + t378 * t461) * MDP(8) + (t378 * t382 + 0.2e1 * t391 * t436) * MDP(7) + (qJDD(3) * t391 - t388 * t394) * MDP(10) + (qJDD(3) * t388 + t391 * t394) * MDP(9) + t378 * MDP(4);
t415 = -t364 * t373 + t465;
t412 = g(1) * t467 + g(2) * t468 - g(3) * t371 + t254;
t304 = t350 * t385 + t386 * t433;
t410 = -qJDD(3) * t304 + t423;
t300 = t388 * t427 + t391 * t441 + t372;
t396 = (-qJD(4) - t441) * t388 + t391 * t427;
t271 = t386 * t300 + t385 * t396;
t329 = t363 * t391 + t375;
t296 = t386 * t329 + t385 * t430;
t409 = qJD(3) * t271 + qJDD(3) * t296 + t456;
t297 = pkin(4) * t330 - qJ(5) * t331 - t364;
t407 = -t337 * t379 - t322 + t483;
t405 = pkin(7) * t394 - t426 - t475;
t270 = t300 * t385 - t386 * t396;
t295 = t329 * t385 - t386 * t430;
t404 = t270 * t319 - t271 * t317 - t296 * t287 + t288 * t295 - t483;
t365 = -pkin(2) - t476;
t403 = pkin(1) * t437 + t363 * t394 + t365 * t378;
t402 = -qJD(3) * t270 - qJDD(3) * t295 + t423;
t401 = -t274 * t319 - qJDD(5) + t412;
t400 = -pkin(7) * qJDD(3) + (t442 - t474) * qJD(3);
t399 = -qJDD(3) * t363 + (t365 * t379 - t441) * qJD(3);
t398 = (-g(1) * (-t364 + t418) - g(2) * t491) * t373;
t395 = -t305 * t287 + t288 * t304 - t458 * t317 + t459 * t319 - t483;
t393 = cos(qJ(1));
t376 = t393 * pkin(1);
t358 = -pkin(3) * t386 - pkin(4);
t356 = pkin(3) * t385 + qJ(5);
t348 = -t364 - t476;
t335 = t369 + t368;
t291 = t297 - t476;
t284 = pkin(3) * t464 + pkin(4) * t319 + qJ(5) * t317;
t275 = t281 + t369;
t1 = [(g(1) * t390 - g(2) * t393) * MDP(2) + (g(1) * t393 + g(2) * t390) * MDP(3) + t416 + (t288 * t348 + t319 * t335 - t409 + t485) * MDP(15) + (-t275 * t319 - t288 * t291 + t409 + t487) * MDP(20) + (t287 * t348 + t317 * t335 + t402 + t486) * MDP(14) + (t275 * t317 + t287 * t291 + t402 + t488) * MDP(18) + (t399 * t388 + (-t403 - t490) * t391 + t457) * MDP(12) + (t399 * t391 + (t403 - t361) * t388 + t439) * MDP(13) + (t404 + t479) * MDP(16) + (t404 + t480) * MDP(19) + ((t378 * t392 - t437) * pkin(1) + t425) * MDP(5) + (((-qJDD(1) - t378) * t389 + (-qJD(1) - t379) * t450) * pkin(1) + t483) * MDP(6) + (t255 * t296 + t280 * t271 - t254 * t295 - t279 * t270 + t292 * t348 + t316 * t335 - g(1) * (t415 - t477) - g(2) * (t376 + t428)) * MDP(17) + (t252 * t296 + t278 * t271 + t251 * t291 + t274 * t275 + t253 * t295 + t277 * t270 - g(1) * (t465 - t477) - g(2) * (t376 + t438) + t398) * MDP(21) + qJDD(1) * MDP(1); (t425 + t426) * MDP(5) + ((-t445 + (-qJD(2) + t379) * t452) * pkin(1) + t483) * MDP(6) + (t400 * t388 + (-t405 - t490) * t391 + t457) * MDP(12) + (t400 * t391 + (t405 - t361) * t388 + t439) * MDP(13) + (-t317 * t443 - t287 * t364 + (t317 * t473 - t459) * qJD(3) + t410 + t486) * MDP(14) + (-t319 * t443 - t288 * t364 + (t319 * t473 - t458) * qJD(3) - t417 + t485) * MDP(15) + (t395 + t479) * MDP(16) + (t255 * t305 - t254 * t304 - t292 * t364 - g(1) * t415 - g(2) * t428 + (t368 - t443) * t316 + t458 * t280 - t459 * t279) * MDP(17) + (-t459 * qJD(3) + t287 * t297 - t421 * t317 + t410 + t488) * MDP(18) + (t395 + t480) * MDP(19) + (t458 * qJD(3) - t288 * t297 + t421 * t319 + t417 + t487) * MDP(20) + (-g(1) * t465 - g(2) * t438 + t251 * t297 + t252 * t305 + t253 * t304 - t421 * t274 + t459 * t277 + t458 * t278 + t398) * MDP(21) + t416; qJDD(3) * MDP(11) + (t407 * t388 - t471) * MDP(12) + (g(3) * t388 + t407 * t391) * MDP(13) + (t449 - t316 * t319 + (qJDD(3) * t386 - t317 * t464) * pkin(3) + t412) * MDP(14) + (qJD(3) * t283 + t316 * t317 + (-qJDD(3) * t385 - t319 * t464) * pkin(3) + t481) * MDP(15) + ((t280 - t282) * t319 + (-t279 + t283) * t317 + (-t287 * t385 - t288 * t386) * pkin(3)) * MDP(16) + (t279 * t282 - t280 * t283 + (-t471 + t254 * t386 + t255 * t385 + (-t316 * t379 + t483) * t388) * pkin(3)) * MDP(17) + (t449 - t284 * t317 + (pkin(4) - t358) * qJDD(3) + t401) * MDP(18) + (-t287 * t356 + t288 * t358 + (t278 - t282) * t319 + (t277 - t446) * t317) * MDP(19) + (qJDD(3) * t356 - t274 * t317 + t284 * t319 + t381 + (0.2e1 * qJD(5) - t283) * qJD(3) - t481) * MDP(20) + (t252 * t356 + t253 * t358 - t274 * t284 - t277 * t282 - g(3) * (-t418 + t472) + t446 * t278 + t483 * (pkin(4) * t370 - qJ(5) * t371 + t473)) * MDP(21) + (t391 * MDP(10) + t388 * MDP(9)) * t378 + (-MDP(7) * t461 + t453 * MDP(8)) * t379 ^ 2; (t279 * t319 + t280 * t317 + t424 - t482) * MDP(17) + (-t277 * t319 + t278 * t317 + t397 - t482) * MDP(21) + (-MDP(17) - MDP(21)) * t420 + (MDP(15) - MDP(20)) * ((-t317 + t440) * qJD(3) + t406) + (MDP(16) + MDP(19)) * (-t317 ^ 2 - t313) + (0.2e1 * qJD(3) * t319 + t489) * (MDP(14) + MDP(18)); (t317 * t319 - qJDD(3)) * MDP(18) + ((t317 + t440) * qJD(3) + t406) * MDP(19) + (-t313 - t394) * MDP(20) + (-qJD(3) * t278 - t401 - t470) * MDP(21);];
tau = t1;
