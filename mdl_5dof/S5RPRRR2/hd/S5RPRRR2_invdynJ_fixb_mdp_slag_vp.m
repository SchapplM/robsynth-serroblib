% Calculate vector of inverse dynamics joint torques for
% S5RPRRR2
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
%   pkin=[a2,a3,a4,a5,d1,d3,d4,d5,theta2]';
% MDP [28x1]
%   Minimal dynamic parameter vector (fixed base model)
%   see S5RPRRR2_convert_par2_MPV_fixb.m
% 
% Output:
% tau [5x1]
%   joint torques of inverse dynamics (contains inertial, gravitational Coriolis and centrifugal forces)

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-05 18:13
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tau = S5RPRRR2_invdynJ_fixb_mdp_slag_vp(qJ, qJD, qJDD, g, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(5,1),zeros(3,1),zeros(9,1),zeros(28,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPRRR2_invdynJ_fixb_mdp_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RPRRR2_invdynJ_fixb_mdp_slag_vp: qJD has to be [5x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [5 1]), ...
  'S5RPRRR2_invdynJ_fixb_mdp_slag_vp: qJDD has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RPRRR2_invdynJ_fixb_mdp_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5RPRRR2_invdynJ_fixb_mdp_slag_vp: pkin has to be [9x1] (double)');
assert(isreal(MDP) && all(size(MDP) == [28 1]), ...
  'S5RPRRR2_invdynJ_fixb_mdp_slag_vp: MDP has to be [28x1] (double)'); 

%% Symbolic Calculation
% From invdyn_joint_fixb_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 18:12:27
% EndTime: 2019-12-05 18:12:37
% DurationCPUTime: 4.84s
% Computational Cost: add. (3565->330), mult. (8776->424), div. (0->0), fcn. (7152->16), ass. (0->160)
t388 = sin(pkin(9));
t392 = sin(qJ(3));
t389 = cos(pkin(9));
t396 = cos(qJ(3));
t456 = t396 * t389;
t344 = t388 * t392 - t456;
t336 = t344 * qJD(1);
t462 = t388 * t396;
t345 = t389 * t392 + t462;
t337 = t345 * qJD(1);
t391 = sin(qJ(4));
t395 = cos(qJ(4));
t313 = t336 * t395 + t337 * t391;
t394 = cos(qJ(5));
t390 = sin(qJ(5));
t414 = -t336 * t391 + t337 * t395;
t464 = t414 * t390;
t275 = t313 * t394 + t464;
t383 = qJDD(3) + qJDD(4);
t378 = qJDD(5) + t383;
t416 = t313 * t390 - t394 * t414;
t497 = t378 * MDP(26) + (-t275 ^ 2 + t416 ^ 2) * MDP(23) - t275 * MDP(22) * t416;
t448 = qJD(1) * t392;
t436 = t388 * t448;
t440 = qJDD(1) * t396;
t441 = qJDD(1) * t392;
t447 = qJD(3) * t396;
t437 = t388 * t440 + (qJD(1) * t447 + t441) * t389;
t321 = -qJD(3) * t436 + t437;
t339 = t345 * qJD(3);
t364 = t389 * t440;
t421 = -t388 * t441 + t364;
t322 = qJD(1) * t339 - t421;
t445 = qJD(4) * t395;
t446 = qJD(4) * t391;
t270 = t321 * t395 - t322 * t391 - t336 * t445 - t337 * t446;
t271 = qJD(4) * t414 + t321 * t391 + t322 * t395;
t443 = qJD(5) * t394;
t438 = t270 * t394 - t271 * t390 - t313 * t443;
t444 = qJD(5) * t390;
t249 = -t414 * t444 + t438;
t434 = t270 * t390 + t271 * t394;
t405 = qJD(5) * t416 - t434;
t387 = qJD(3) + qJD(4);
t465 = t313 * t387;
t466 = t414 * t387;
t439 = -qJD(4) - qJD(5);
t381 = qJD(3) - t439;
t494 = t381 * t416;
t495 = t275 * t381;
t496 = t383 * MDP(19) + t313 * MDP(15) * t414 + (-t271 + t466) * MDP(18) + (-t313 ^ 2 + t414 ^ 2) * MDP(16) + (t270 + t465) * MDP(17) + (t249 + t495) * MDP(24) + (t405 - t494) * MDP(25) + t497;
t471 = pkin(6) + qJ(2);
t358 = t471 * t388;
t346 = qJD(1) * t358;
t359 = t471 * t389;
t347 = qJD(1) * t359;
t413 = t346 * t392 - t347 * t396;
t300 = -pkin(7) * t336 - t413;
t297 = t395 * t300;
t481 = -t346 * t396 - t347 * t392;
t299 = -pkin(7) * t337 + t481;
t298 = qJD(3) * pkin(3) + t299;
t417 = -t298 * t391 - t297;
t492 = pkin(8) * t313;
t259 = -t417 - t492;
t373 = -pkin(2) * t389 - pkin(1);
t353 = qJD(1) * t373 + qJD(2);
t325 = pkin(3) * t336 + t353;
t285 = pkin(4) * t313 + t325;
t386 = pkin(9) + qJ(3);
t382 = qJ(4) + t386;
t374 = qJ(5) + t382;
t368 = sin(t374);
t369 = cos(t374);
t393 = sin(qJ(1));
t397 = cos(qJ(1));
t423 = g(1) * t397 + g(2) * t393;
t489 = g(3) * t368 + t259 * t444 + t285 * t275 + t369 * t423;
t442 = qJD(1) * qJD(2);
t476 = qJDD(1) * t471 + t442;
t329 = t476 * t388;
t330 = t476 * t389;
t429 = -t329 * t396 - t392 * t330;
t265 = qJDD(3) * pkin(3) - pkin(7) * t321 + qJD(3) * t413 + t429;
t415 = -t392 * t329 + t396 * t330;
t269 = -pkin(7) * t322 + qJD(3) * t481 + t415;
t407 = qJD(4) * t417 + t265 * t395 - t391 * t269;
t247 = pkin(4) * t383 - pkin(8) * t270 + t407;
t477 = (qJD(4) * t298 + t269) * t395 + t391 * t265 - t300 * t446;
t248 = -pkin(8) * t271 + t477;
t488 = -g(3) * t369 + t394 * t247 - t390 * t248 + t285 * t416 + t368 * t423;
t371 = sin(t382);
t372 = cos(t382);
t487 = g(3) * t371 + t325 * t313 + t372 * t423 - t477;
t483 = pkin(8) * t414;
t428 = -t358 * t396 - t359 * t392;
t310 = -pkin(7) * t345 + t428;
t453 = -t358 * t392 + t359 * t396;
t311 = -pkin(7) * t344 + t453;
t454 = t310 * t391 + t311 * t395;
t480 = qJ(2) * qJDD(1);
t422 = g(1) * t393 - g(2) * t397;
t479 = qJDD(2) - t422;
t478 = -g(3) * t372 - t325 * t414 + t371 * t423 + t407;
t475 = pkin(3) * t391;
t472 = t339 * pkin(3);
t470 = qJDD(1) * pkin(1);
t295 = t391 * t300;
t433 = t298 * t395 - t295;
t258 = t433 - t483;
t256 = pkin(4) * t387 + t258;
t469 = t256 * t394;
t461 = t389 * MDP(4);
t460 = t390 * t378;
t459 = t391 * t394;
t458 = t394 * t259;
t457 = t394 * t378;
t455 = t299 * t395 - t295;
t451 = t388 ^ 2 + t389 ^ 2;
t432 = -t299 * t391 - t297;
t431 = t310 * t395 - t311 * t391;
t426 = -qJD(5) * t256 - t248;
t424 = 0.2e1 * t451;
t420 = -t390 * t256 - t458;
t324 = -t344 * t391 + t345 * t395;
t262 = -pkin(8) * t324 + t431;
t323 = t344 * t395 + t345 * t391;
t263 = -pkin(8) * t323 + t454;
t419 = t262 * t394 - t263 * t390;
t418 = t262 * t390 + t263 * t394;
t283 = t323 * t394 + t324 * t390;
t284 = -t323 * t390 + t324 * t394;
t327 = pkin(3) * t344 + t373;
t412 = t470 - t479;
t348 = qJDD(1) * t373 + qJDD(2);
t408 = -t358 * t447 + qJD(2) * t456 + (-qJD(2) * t388 - qJD(3) * t359) * t392;
t289 = -pkin(7) * t339 + t408;
t338 = t344 * qJD(3);
t401 = -qJD(2) * t345 - qJD(3) * t453;
t290 = pkin(7) * t338 + t401;
t411 = t289 * t395 + t290 * t391 + t310 * t445 - t311 * t446;
t301 = pkin(3) * t322 + t348;
t406 = -qJD(4) * t454 - t289 * t391 + t290 * t395;
t403 = t424 * t442 - t423;
t380 = cos(t386);
t379 = sin(t386);
t376 = pkin(3) * t395 + pkin(4);
t293 = pkin(4) * t323 + t327;
t291 = pkin(3) * t337 + pkin(4) * t414;
t282 = qJD(4) * t324 - t338 * t391 + t339 * t395;
t281 = -qJD(4) * t323 - t338 * t395 - t339 * t391;
t272 = pkin(4) * t282 + t472;
t261 = t455 - t483;
t260 = t432 + t492;
t255 = pkin(4) * t271 + t301;
t254 = qJD(5) * t284 + t281 * t390 + t282 * t394;
t253 = -qJD(5) * t283 + t281 * t394 - t282 * t390;
t252 = -pkin(8) * t281 + t406;
t251 = -pkin(8) * t282 + t411;
t1 = [qJDD(1) * MDP(1) + (t249 * t284 - t253 * t416) * MDP(22) + (-t249 * t283 - t253 * t275 + t254 * t416 + t284 * t405) * MDP(23) + (t272 * t275 - t293 * t405 + t255 * t283 + t285 * t254 + (-qJD(5) * t418 - t251 * t390 + t252 * t394) * t381 + t419 * t378 + t422 * t369) * MDP(27) + t422 * MDP(2) + (-t272 * t416 + t293 * t249 + t255 * t284 + t285 * t253 - (qJD(5) * t419 + t251 * t394 + t252 * t390) * t381 - t418 * t378 - t422 * t368) * MDP(28) + t423 * MDP(3) + (t424 * t480 + t403) * MDP(6) + (qJD(3) * t401 + qJDD(3) * t428 + t373 * t322 + t353 * t339 + t348 * t344 + t380 * t422) * MDP(13) + (t270 * t324 + t281 * t414) * MDP(15) + (-t270 * t323 - t271 * t324 - t281 * t313 - t282 * t414) * MDP(16) + (-qJD(3) * t339 - qJDD(3) * t344) * MDP(11) + (-qJD(3) * t338 + qJDD(3) * t345) * MDP(10) + (t321 * t345 - t337 * t338) * MDP(8) + (-t321 * t344 - t322 * t345 + t336 * t338 - t337 * t339) * MDP(9) + (pkin(1) * t412 + (t451 * t480 + t403) * qJ(2)) * MDP(7) + (-qJD(3) * t408 - qJDD(3) * t453 + t373 * t321 - t353 * t338 + t348 * t345 - t379 * t422) * MDP(14) + (-t254 * t381 - t283 * t378) * MDP(25) + (t253 * t381 + t284 * t378) * MDP(24) + (-t282 * t387 - t323 * t383) * MDP(18) + (t281 * t387 + t324 * t383) * MDP(17) + (t327 * t270 + t325 * t281 + t301 * t324 - t371 * t422 - t383 * t454 - t387 * t411 + t414 * t472) * MDP(21) + (t327 * t271 + t325 * t282 + t301 * t323 + t313 * t472 + t372 * t422 + t383 * t431 + t387 * t406) * MDP(20) + (-MDP(5) * t388 + t461) * (t412 + t470); t479 * MDP(7) - t364 * MDP(13) + t437 * MDP(14) + (t271 + t466) * MDP(20) + (t270 - t465) * MDP(21) + (-t405 - t494) * MDP(27) + (t249 - t495) * MDP(28) + (-t461 - pkin(1) * MDP(7) + (MDP(13) * t392 + MDP(5)) * t388) * qJDD(1) + ((qJD(1) * t462 + t389 * t448 + t337) * MDP(13) + (-t336 - t436) * MDP(14)) * qJD(3) + (-MDP(7) * qJ(2) - MDP(6)) * qJD(1) ^ 2 * t451; (t455 * t387 + (-t337 * t414 - t383 * t391 - t387 * t445) * pkin(3) + t487) * MDP(21) + (t291 * t416 + (-t376 * t378 - t247 + (-t439 * t475 + t260) * t381) * t390 + (-t378 * t475 + (-pkin(3) * t445 - qJD(5) * t376 + t261) * t381 + t426) * t394 + t489) * MDP(28) + t421 * MDP(11) + qJDD(3) * MDP(12) + (t437 + (t336 - t436) * qJD(3)) * MDP(10) + (g(3) * t379 + t353 * t336 + t380 * t423 - t415) * MDP(14) + (-t432 * t387 + (-t313 * t337 + t383 * t395 - t387 * t446) * pkin(3) + t478) * MDP(20) + t337 * t336 * MDP(8) + (-g(3) * t380 - t353 * t337 + t379 * t423 + t429) * MDP(13) + (-t336 ^ 2 + t337 ^ 2) * MDP(9) + (t376 * t457 - t291 * t275 - (t260 * t394 - t261 * t390) * t381 + (-t391 * t460 + (-t390 * t395 - t459) * t381 * qJD(4)) * pkin(3) + ((-pkin(3) * t459 - t376 * t390) * t381 + t420) * qJD(5) + t488) * MDP(27) + t496; (-t387 * t417 + t478) * MDP(20) + (t387 * t433 + t487) * MDP(21) + (-(-t258 * t390 - t458) * t381 + t420 * qJD(5) + (-t275 * t414 - t381 * t444 + t457) * pkin(4) + t488) * MDP(27) + ((-t259 * t381 - t247) * t390 + (t258 * t381 + t426) * t394 + (-t381 * t443 + t414 * t416 - t460) * pkin(4) + t489) * MDP(28) + t496; (t438 + t495) * MDP(24) + (-t434 - t494) * MDP(25) + (-t381 * t420 + t488) * MDP(27) + (-t394 * t248 - t390 * t247 + (-t259 * t390 + t469) * t381 + t489) * MDP(28) + (-MDP(24) * t464 + MDP(25) * t416 + t420 * MDP(27) - MDP(28) * t469) * qJD(5) + t497;];
tau = t1;
