% Calculate vector of inverse dynamics joint torques for
% S5RRPRP7
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
%   see S5RRPRP7_convert_par2_MPV_fixb.m
% 
% Output:
% tau [5x1]
%   joint torques of inverse dynamics (contains inertial, gravitational Coriolis and centrifugal forces)

% Quelle: HybrDyn-Toolbox
% Datum: 2021-01-15 20:42
% Revision: 24b2e7d74a0c1a3b64fa2f8f5ad758691ad61af3 (2021-01-15)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tau = S5RRPRP7_invdynJ_fixb_mdp_slag_vp(qJ, qJD, qJDD, g, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(5,1),zeros(3,1),zeros(8,1),zeros(25,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRPRP7_invdynJ_fixb_mdp_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RRPRP7_invdynJ_fixb_mdp_slag_vp: qJD has to be [5x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [5 1]), ...
  'S5RRPRP7_invdynJ_fixb_mdp_slag_vp: qJDD has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RRPRP7_invdynJ_fixb_mdp_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RRPRP7_invdynJ_fixb_mdp_slag_vp: pkin has to be [8x1] (double)');
assert(isreal(MDP) && all(size(MDP) == [25 1]), ...
  'S5RRPRP7_invdynJ_fixb_mdp_slag_vp: MDP has to be [25x1] (double)'); 

%% Symbolic Calculation
% From invdyn_joint_fixb_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2021-01-15 20:41:36
% EndTime: 2021-01-15 20:41:49
% DurationCPUTime: 5.56s
% Computational Cost: add. (3733->457), mult. (8540->559), div. (0->0), fcn. (5973->10), ass. (0->177)
t367 = qJ(2) + pkin(8);
t361 = sin(t367);
t374 = sin(qJ(1));
t377 = cos(qJ(1));
t475 = -g(1) * t374 + g(2) * t377;
t477 = t361 * t475;
t376 = cos(qJ(2));
t458 = cos(pkin(8));
t415 = t458 * t376;
t351 = qJD(1) * t415;
t370 = sin(pkin(8));
t373 = sin(qJ(2));
t434 = qJD(1) * t373;
t326 = t370 * t434 - t351;
t319 = qJD(4) + t326;
t474 = -g(1) * t377 - g(2) * t374;
t476 = t361 * t474;
t416 = t458 * t373;
t338 = t370 * t376 + t416;
t460 = t376 * pkin(2);
t360 = pkin(1) + t460;
t395 = -t370 * t373 + t415;
t298 = -pkin(3) * t395 - pkin(7) * t338 - t360;
t371 = -qJ(3) - pkin(6);
t347 = t371 * t376;
t421 = t371 * t373;
t310 = -t347 * t458 + t370 * t421;
t372 = sin(qJ(4));
t375 = cos(qJ(4));
t438 = t372 * t298 + t375 * t310;
t429 = qJD(1) * qJD(2);
t420 = t373 * t429;
t385 = qJDD(1) * t338 - t370 * t420;
t381 = qJD(2) * t351 + t385;
t329 = t338 * qJD(1);
t473 = qJD(2) * t329;
t341 = qJD(1) * t347;
t332 = t370 * t341;
t340 = qJD(1) * t421;
t459 = qJD(2) * pkin(2);
t335 = t340 + t459;
t299 = t335 * t458 + t332;
t289 = -qJD(2) * pkin(3) - t299;
t311 = -t375 * qJD(2) + t329 * t372;
t313 = qJD(2) * t372 + t329 * t375;
t256 = t311 * pkin(4) - t313 * qJ(5) + t289;
t328 = t338 * qJD(2);
t427 = qJDD(1) * t373;
t405 = -qJDD(1) * t415 + t370 * t427;
t301 = qJD(1) * t328 + t405;
t296 = qJDD(4) + t301;
t468 = pkin(2) * t370;
t354 = pkin(7) + t468;
t447 = t354 * t296;
t472 = t256 * t319 - t447;
t362 = cos(t367);
t463 = g(3) * t361;
t471 = t362 * t474 - t463;
t470 = t313 ^ 2;
t469 = t319 ^ 2;
t467 = pkin(4) * t296;
t462 = g(3) * t362;
t461 = g(3) * t376;
t457 = qJ(5) * t296;
t456 = qJDD(1) * pkin(1);
t345 = -qJD(1) * t360 + qJD(3);
t274 = pkin(3) * t326 - pkin(7) * t329 + t345;
t417 = t458 * t341;
t300 = t370 * t335 - t417;
t290 = qJD(2) * pkin(7) + t300;
t255 = t274 * t372 + t290 * t375;
t249 = qJ(5) * t319 + t255;
t455 = t249 * t319;
t454 = t255 * t319;
t428 = qJD(4) * qJD(2);
t432 = qJD(4) * t372;
t268 = -t372 * qJDD(2) + t329 * t432 + (-t381 - t428) * t375;
t453 = t268 * t372;
t452 = t311 * t326;
t451 = t311 * t329;
t450 = t313 * t311;
t414 = t313 * t319;
t449 = t313 * t329;
t448 = t338 * t375;
t446 = t371 * t374;
t286 = t372 * t296;
t445 = t372 * t319;
t444 = t372 * t374;
t443 = t374 * t375;
t287 = t375 * t296;
t442 = t375 * t377;
t441 = t377 * t372;
t380 = -t375 * qJDD(2) + t372 * t381;
t269 = t313 * qJD(4) + t380;
t431 = qJD(4) * t375;
t440 = -t372 * t269 - t311 * t431;
t418 = qJD(2) * t371;
t391 = -qJD(3) * t373 + t376 * t418;
t295 = qJDD(2) * pkin(2) + qJD(1) * t391 + qJDD(1) * t421;
t324 = qJD(3) * t376 + t373 * t418;
t302 = qJD(1) * t324 - qJDD(1) * t347;
t264 = t458 * t295 - t370 * t302;
t265 = t370 * t295 + t458 * t302;
t281 = pkin(2) * t434 + pkin(3) * t329 + pkin(7) * t326;
t304 = t340 * t458 + t332;
t439 = t372 * t281 + t375 * t304;
t303 = t340 * t370 - t417;
t406 = pkin(4) * t372 - qJ(5) * t375;
t437 = -qJD(5) * t372 + t319 * t406 - t303;
t436 = (g(1) * t442 + g(2) * t443) * t361;
t368 = t373 ^ 2;
t435 = -t376 ^ 2 + t368;
t433 = qJD(4) * t354;
t254 = t274 * t375 - t290 * t372;
t430 = qJD(5) - t254;
t426 = qJDD(1) * t376;
t425 = pkin(2) * t420 + qJDD(3);
t424 = t373 * t459;
t423 = t458 * pkin(2);
t422 = t319 * t433;
t279 = t324 * t370 - t458 * t391;
t309 = -t347 * t370 - t371 * t416;
t257 = -pkin(2) * t426 + t301 * pkin(3) - pkin(7) * t381 + t425 - t456;
t263 = qJDD(2) * pkin(7) + t265;
t412 = -t375 * t257 + t372 * t263 + t274 * t432 + t290 * t431;
t320 = t362 * t444 + t442;
t322 = t362 * t441 - t443;
t411 = -g(1) * t320 + g(2) * t322;
t321 = t362 * t443 - t441;
t323 = t362 * t442 + t444;
t410 = g(1) * t321 - g(2) * t323;
t262 = -qJDD(2) * pkin(3) - t264;
t407 = t375 * pkin(4) + t372 * qJ(5);
t248 = -pkin(4) * t319 + t430;
t404 = t248 * t375 - t249 * t372;
t403 = pkin(3) + t407;
t402 = t286 + (t326 * t375 + t431) * t319;
t400 = -t319 * t432 - t326 * t445 + t287;
t399 = t422 + t462;
t398 = -0.2e1 * pkin(1) * t429 - pkin(6) * qJDD(2);
t331 = t395 * qJD(2);
t397 = t331 * t372 + t338 * t431;
t396 = -t331 * t375 + t338 * t432;
t394 = t372 * t257 + t375 * t263 + t274 * t431 - t290 * t432;
t280 = t324 * t458 + t370 * t391;
t282 = pkin(3) * t328 - pkin(7) * t331 + t424;
t393 = t375 * t280 + t372 * t282 + t298 * t431 - t310 * t432;
t392 = t289 * t319 - t447;
t318 = -qJDD(1) * t360 + t425;
t389 = -t462 - t476;
t378 = qJD(2) ^ 2;
t388 = -pkin(6) * t378 + 0.2e1 * t456 - t475;
t379 = qJD(1) ^ 2;
t387 = pkin(1) * t379 - pkin(6) * qJDD(1) - t474;
t386 = g(1) * t322 + g(2) * t320 + t372 * t463 - t412;
t384 = t256 * t313 + qJDD(5) - t386;
t383 = -g(1) * t323 - g(2) * t321 - t375 * t463 + t394;
t357 = t371 * t377;
t355 = -t423 - pkin(3);
t346 = -t370 * pkin(3) + pkin(7) * t458;
t344 = pkin(3) * t458 + t370 * pkin(7) + pkin(2);
t336 = -t423 - t403;
t307 = t344 * t376 + t346 * t373 + pkin(1);
t271 = pkin(4) * t313 + qJ(5) * t311;
t270 = t338 * t406 + t309;
t259 = pkin(4) * t395 - t298 * t375 + t310 * t372;
t258 = -qJ(5) * t395 + t438;
t251 = -pkin(4) * t329 - t281 * t375 + t304 * t372;
t250 = qJ(5) * t329 + t439;
t247 = t311 * t319 - t268;
t246 = t406 * t331 + (qJD(4) * t407 - qJD(5) * t375) * t338 + t279;
t245 = -pkin(4) * t328 + qJD(4) * t438 + t280 * t372 - t282 * t375;
t244 = qJ(5) * t328 - qJD(5) * t395 + t393;
t243 = pkin(4) * t269 + qJ(5) * t268 - qJD(5) * t313 + t262;
t242 = qJDD(5) + t412 - t467;
t241 = qJD(5) * t319 + t394 + t457;
t1 = [qJDD(1) * MDP(1) + (qJDD(1) * t368 + 0.2e1 * t376 * t420) * MDP(4) + 0.2e1 * (t373 * t426 - t429 * t435) * MDP(5) + (qJDD(2) * t373 + t376 * t378) * MDP(6) + (qJDD(2) * t376 - t373 * t378) * MDP(7) + (t373 * t398 + t376 * t388) * MDP(9) + (-t373 * t388 + t376 * t398) * MDP(10) + (-qJDD(2) * t309 - t301 * t360 - t318 * t395 + t328 * t345 - t475 * t362 + (pkin(2) * t326 * t373 - t279) * qJD(2)) * MDP(11) + (-t280 * qJD(2) - t310 * qJDD(2) + t318 * t338 + t329 * t424 + t345 * t331 - t360 * t381 + t477) * MDP(12) + (-t264 * t338 + t265 * t395 + t279 * t329 - t280 * t326 - t299 * t331 - t300 * t328 - t310 * t301 + t309 * t381 + t474) * MDP(13) + (t265 * t310 + t300 * t280 - t264 * t309 - t299 * t279 - t318 * t360 + t345 * t424 - g(1) * (-t360 * t374 - t357) - g(2) * (t360 * t377 - t446)) * MDP(14) + (-t268 * t448 - t313 * t396) * MDP(15) + ((-t311 * t375 - t313 * t372) * t331 + (t453 - t269 * t375 + (t311 * t372 - t313 * t375) * qJD(4)) * t338) * MDP(16) + (t268 * t395 + t287 * t338 + t313 * t328 - t319 * t396) * MDP(17) + (t269 * t395 - t286 * t338 - t311 * t328 - t319 * t397) * MDP(18) + (-t296 * t395 + t319 * t328) * MDP(19) + (t412 * t395 + t254 * t328 + t279 * t311 + t309 * t269 + ((-qJD(4) * t310 + t282) * t319 + t298 * t296 + t289 * qJD(4) * t338) * t375 + ((-qJD(4) * t298 - t280) * t319 - t310 * t296 + t262 * t338 + t289 * t331) * t372 + t410) * MDP(20) + (-t255 * t328 + t262 * t448 - t309 * t268 + t279 * t313 - t289 * t396 - t296 * t438 - t319 * t393 + t394 * t395 + t411) * MDP(21) + (t243 * t338 * t372 + t242 * t395 - t245 * t319 + t246 * t311 - t248 * t328 + t256 * t397 - t259 * t296 + t269 * t270 + t410) * MDP(22) + (-t244 * t311 + t245 * t313 - t258 * t269 - t259 * t268 - t477 + t404 * t331 + (-t241 * t372 + t242 * t375 + (-t248 * t372 - t249 * t375) * qJD(4)) * t338) * MDP(23) + (-t241 * t395 - t243 * t448 + t244 * t319 - t246 * t313 + t249 * t328 + t256 * t396 + t258 * t296 + t268 * t270 - t411) * MDP(24) + (t241 * t258 + t249 * t244 + t243 * t270 + t256 * t246 + t242 * t259 + t248 * t245 - g(1) * (-pkin(4) * t321 - qJ(5) * t320 - t307 * t374 - t357) - g(2) * (pkin(4) * t323 + qJ(5) * t322 + t307 * t377 - t446)) * MDP(25) - t475 * MDP(2) - t474 * MDP(3); MDP(6) * t427 + MDP(7) * t426 + qJDD(2) * MDP(8) + (t373 * t387 - t461) * MDP(9) + (g(3) * t373 + t376 * t387) * MDP(10) + (t303 * qJD(2) - t345 * t329 + (qJDD(2) * t458 - t326 * t434) * pkin(2) + t389 + t264) * MDP(11) + (qJD(2) * t304 + t326 * t345 + (-qJDD(2) * t370 - t329 * t434) * pkin(2) - t265 - t471) * MDP(12) + (-t301 * t468 - t381 * t423 - (-t300 + t303) * t329 + (t304 - t299) * t326) * MDP(13) + (t299 * t303 - t300 * t304 + (t458 * t264 - t461 + t265 * t370 + (-qJD(1) * t345 - t474) * t373) * pkin(2)) * MDP(14) + (t375 * t414 - t453) * MDP(15) + ((-t268 - t452) * t375 - t313 * t445 + t440) * MDP(16) + (t402 - t449) * MDP(17) + (t400 + t451) * MDP(18) - t319 * t329 * MDP(19) + (-t254 * t329 + t355 * t269 - t303 * t311 + (-t462 - t262 + (-t281 - t433) * t319) * t375 + (t304 * t319 + t392) * t372 + t436) * MDP(20) + (-t355 * t268 + t439 * t319 + t255 * t329 - t303 * t313 + t392 * t375 + (t262 + t399 + t476) * t372) * MDP(21) + (t248 * t329 + t251 * t319 + t269 * t336 + t437 * t311 + (-t243 - t399) * t375 + t472 * t372 + t436) * MDP(22) + (t250 * t311 - t251 * t313 + (t248 * t326 - t269 * t354 + t241 + (t313 * t354 + t248) * qJD(4)) * t375 + (-t249 * t326 - t268 * t354 + t242 + (t311 * t354 - t249) * qJD(4)) * t372 + t471) * MDP(23) + (-t249 * t329 - t250 * t319 + t268 * t336 - t437 * t313 - t472 * t375 + (-t243 + t389 - t422) * t372) * MDP(24) + (t243 * t336 - t249 * t250 - t248 * t251 - g(3) * (pkin(7) * t361 + t460) - t403 * t462 + t437 * t256 + (qJD(4) * t404 + t241 * t375 + t242 * t372) * t354 + t474 * (-t344 * t373 + t346 * t376 - t361 * t407)) * MDP(25) + (-MDP(4) * t373 * t376 + MDP(5) * t435) * t379; (t405 + 0.2e1 * t473) * MDP(11) + ((t351 - t326) * qJD(2) + t385) * MDP(12) + (-t326 ^ 2 - t329 ^ 2) * MDP(13) + (t299 * t329 + t300 * t326 + t318 + t475) * MDP(14) + (t400 - t451) * MDP(20) + (-t375 * t469 - t286 - t449) * MDP(21) + (-t319 * t445 + t287 - t451) * MDP(22) + ((t268 - t452) * t375 + t372 * t414 + t440) * MDP(23) + (t402 + t449) * MDP(24) + (-t256 * t329 + (-t242 + t455) * t375 + (t248 * t319 + t241) * t372 + t475) * MDP(25); MDP(15) * t450 + (-t311 ^ 2 + t470) * MDP(16) + t247 * MDP(17) + (-t329 * t431 - t372 * t428 - t380 + t414) * MDP(18) + t296 * MDP(19) + (-t289 * t313 + t386 + t454) * MDP(20) + (t254 * t319 + t289 * t311 - t383) * MDP(21) + (-t271 * t311 - t384 + t454 + 0.2e1 * t467) * MDP(22) + (pkin(4) * t268 - qJ(5) * t269 + (t249 - t255) * t313 + (t248 - t430) * t311) * MDP(23) + (0.2e1 * t457 - t256 * t311 + t271 * t313 + (0.2e1 * qJD(5) - t254) * t319 + t383) * MDP(24) + (t241 * qJ(5) - t242 * pkin(4) - t256 * t271 - t248 * t255 - g(1) * (-pkin(4) * t322 + qJ(5) * t323) - g(2) * (-pkin(4) * t320 + qJ(5) * t321) + t406 * t463 + t430 * t249) * MDP(25); (-qJDD(4) - t405 + t450 - t473) * MDP(22) + t247 * MDP(23) + (-t469 - t470) * MDP(24) + (t384 - t455 - t467) * MDP(25);];
tau = t1;
