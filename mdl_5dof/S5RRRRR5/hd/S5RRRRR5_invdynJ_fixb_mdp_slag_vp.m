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
% Datum: 2020-01-03 12:13
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
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
% StartTime: 2020-01-03 12:13:29
% EndTime: 2020-01-03 12:13:33
% DurationCPUTime: 2.32s
% Computational Cost: add. (2230->285), mult. (3035->368), div. (0->0), fcn. (1907->16), ass. (0->167)
t355 = cos(qJ(2));
t450 = pkin(1) * t355;
t327 = qJDD(1) * t450;
t340 = qJDD(1) + qJDD(2);
t350 = sin(qJ(2));
t440 = pkin(1) * qJD(1);
t412 = t350 * t440;
t274 = pkin(2) * t340 - qJD(2) * t412 + t327;
t342 = qJD(1) + qJD(2);
t420 = qJD(1) * t355;
t411 = pkin(1) * t420;
t295 = pkin(2) * t342 + t411;
t349 = sin(qJ(3));
t354 = cos(qJ(3));
t404 = qJD(2) * t420;
t413 = qJDD(1) * t350;
t369 = (t404 + t413) * pkin(1);
t456 = -t349 * t274 - (qJD(3) * t295 + t369) * t354;
t426 = t349 * t350;
t315 = pkin(1) * t426;
t396 = qJD(3) * t412;
t418 = qJD(3) * t349;
t388 = t349 * pkin(1) * t404 + qJDD(1) * t315 + t295 * t418 + (-t274 + t396) * t354;
t331 = qJDD(3) + t340;
t448 = pkin(3) * t331;
t237 = t388 - t448;
t346 = qJ(1) + qJ(2);
t337 = qJ(3) + t346;
t321 = sin(t337);
t322 = cos(t337);
t390 = -g(2) * t322 - g(3) * t321;
t378 = -t237 + t390;
t453 = g(2) * t321 - g(3) * t322;
t352 = cos(qJ(5));
t353 = cos(qJ(4));
t424 = t352 * t353;
t347 = sin(qJ(5));
t348 = sin(qJ(4));
t431 = t347 * t348;
t282 = -t424 + t431;
t283 = t347 * t353 + t348 * t352;
t332 = qJD(3) + t342;
t341 = qJD(4) + qJD(5);
t454 = t341 * t283;
t239 = t282 * t331 + t332 * t454;
t374 = t341 * t282;
t425 = t350 * t354;
t379 = t349 * t355 + t425;
t277 = t379 * t440;
t392 = pkin(2) * t418 - t277;
t323 = pkin(2) * t349 + pkin(8);
t449 = pkin(2) * t354;
t324 = -pkin(3) - t449;
t357 = qJD(4) ^ 2;
t452 = t323 * t357 + t324 * t331 + t332 * t392;
t451 = -pkin(8) - pkin(9);
t447 = pkin(3) * t332;
t446 = pkin(4) * t353;
t325 = pkin(2) + t450;
t422 = pkin(1) * t425 + t325 * t349;
t276 = pkin(8) + t422;
t442 = -pkin(9) - t276;
t441 = -pkin(9) - t323;
t268 = t295 * t349 + t354 * t412;
t261 = pkin(8) * t332 + t268;
t403 = pkin(9) * t332 + t261;
t250 = t403 * t353;
t439 = t250 * t352;
t417 = qJD(3) * t354;
t254 = t325 * t418 + (qJD(2) * t379 + t350 * t417) * pkin(1);
t438 = t254 * t332;
t437 = t268 * t332;
t345 = qJ(4) + qJ(5);
t333 = sin(t345);
t436 = t321 * t333;
t435 = t322 * t333;
t434 = t331 * t348;
t433 = t331 * t353;
t432 = t332 * t348;
t428 = t348 * t353;
t416 = qJD(4) * t348;
t329 = pkin(4) * t416;
t423 = t329 + t392;
t343 = t348 ^ 2;
t421 = -t353 ^ 2 + t343;
t415 = qJD(4) * t353;
t414 = qJD(5) * t347;
t409 = pkin(2) * t417;
t408 = t332 * t431;
t407 = t332 * t424;
t326 = -pkin(3) - t446;
t406 = qJD(4) * t451;
t405 = t332 * t415;
t334 = sin(t346);
t336 = cos(t346);
t402 = g(2) * t334 - g(3) * t336;
t401 = qJD(4) * t442;
t400 = qJD(4) * t441;
t312 = t349 * t412;
t267 = t295 * t354 - t312;
t399 = t325 * t354 - t315;
t398 = qJD(1) * (-qJD(2) + t342);
t397 = qJD(2) * (-qJD(1) - t342);
t376 = t332 * t416 - t433;
t234 = pkin(4) * t376 + t237;
t252 = t326 * t332 - t267;
t395 = g(2) * t435 + g(3) * t436 + t234 * t283 - t252 * t374;
t260 = -t267 - t447;
t394 = t260 * t415 - t348 * t378;
t275 = -pkin(3) - t399;
t306 = t349 * t396;
t393 = t306 + t453;
t391 = -t268 + t329;
t249 = t403 * t348;
t248 = qJD(4) * pkin(4) - t249;
t387 = -t248 * t347 - t439;
t262 = t442 * t348;
t338 = t353 * pkin(9);
t263 = t276 * t353 + t338;
t386 = t262 * t352 - t263 * t347;
t385 = t262 * t347 + t263 * t352;
t280 = t441 * t348;
t281 = t323 * t353 + t338;
t384 = t280 * t352 - t281 * t347;
t383 = t280 * t347 + t281 * t352;
t308 = t451 * t348;
t309 = pkin(8) * t353 + t338;
t382 = t308 * t352 - t309 * t347;
t381 = t308 * t347 + t309 * t352;
t377 = -g(2) * t336 - g(3) * t334 + t327;
t372 = pkin(8) * t357 - t437 - t448;
t371 = t275 * t331 + t276 * t357 + t438;
t238 = qJD(5) * t407 + t283 * t331 - t341 * t408 + t352 * t405;
t269 = -t407 + t408;
t271 = t283 * t332;
t339 = qJDD(4) + qJDD(5);
t370 = t271 * t269 * MDP(17) + (t269 * t341 + t238) * MDP(19) + (t271 * t341 - t239) * MDP(20) + (-t269 ^ 2 + t271 ^ 2) * MDP(18) + t339 * MDP(21);
t368 = -pkin(8) * qJDD(4) + (t267 - t447) * qJD(4);
t236 = pkin(8) * t331 - t306 - t456;
t367 = -t260 * t332 - t236 + t453;
t253 = t325 * t417 + (-t350 * t418 + (t354 * t355 - t426) * qJD(2)) * pkin(1);
t366 = -qJDD(4) * t276 + (t275 * t332 - t253) * qJD(4);
t335 = cos(t345);
t365 = t234 * t282 + t252 * t454 + t335 * t390;
t364 = (-t238 * t282 - t239 * t283 + t269 * t374 - t271 * t454) * MDP(18) + (t238 * t283 - t271 * t374) * MDP(17) + (t283 * t339 - t341 * t374) * MDP(19) + (-t282 * t339 - t341 * t454) * MDP(20) + 0.2e1 * (-qJD(4) * t332 * t421 + t331 * t428) * MDP(11) + (t331 * t343 + 0.2e1 * t348 * t405) * MDP(10) + (qJDD(4) * t348 + t353 * t357) * MDP(12) + (qJDD(4) * t353 - t348 * t357) * MDP(13) + t331 * MDP(7);
t363 = -t388 + t390;
t362 = MDP(4) * t340 + t364;
t278 = t354 * t411 - t312;
t361 = -qJDD(4) * t323 + (t324 * t332 + t278 - t409) * qJD(4);
t228 = -t261 * t415 + qJDD(4) * pkin(4) - t236 * t348 + (-t405 - t434) * pkin(9);
t360 = t252 * t269 + t250 * t414 + g(1) * t333 + (-t250 * t341 - t228) * t347 + t453 * t335;
t229 = -pkin(9) * t376 + t236 * t353 - t261 * t416;
t359 = -g(1) * t335 + g(2) * t436 - g(3) * t435 + qJD(5) * t387 + t228 * t352 - t347 * t229 - t252 * t271;
t358 = t393 + t456;
t356 = cos(qJ(1));
t351 = sin(qJ(1));
t302 = t326 - t449;
t289 = t353 * t406;
t288 = t348 * t406;
t273 = t275 - t446;
t265 = -t348 * t409 + t353 * t400;
t264 = t348 * t400 + t353 * t409;
t256 = t260 * t416;
t251 = t329 + t254;
t244 = -t253 * t348 + t353 * t401;
t243 = t253 * t353 + t348 * t401;
t1 = [((t340 * t355 + t350 * t397) * pkin(1) + t377) * MDP(5) + (t348 * t371 + t353 * t366 + t394) * MDP(16) + (-t253 * t332 - t331 * t422 + t358) * MDP(9) + (t251 * t271 + t273 * t238 - (qJD(5) * t386 + t243 * t352 + t244 * t347) * t341 - t385 * t339 + t395) * MDP(23) + qJDD(1) * MDP(1) + t362 + (t256 + t366 * t348 + (-t371 + t378) * t353) * MDP(15) + (t251 * t269 + t273 * t239 + (-qJD(5) * t385 - t243 * t347 + t244 * t352) * t341 + t386 * t339 + t365) * MDP(22) + (t331 * t399 + t363 - t438) * MDP(8) + (-g(2) * t356 - g(3) * t351) * MDP(2) + (g(2) * t351 - g(3) * t356) * MDP(3) + (((-qJDD(1) - t340) * t350 + t355 * t397) * pkin(1) + t402) * MDP(6); (t256 + t361 * t348 + (t378 - t452) * t353) * MDP(15) + (t348 * t452 + t361 * t353 + t394) * MDP(16) + ((t355 * t398 - t413) * pkin(1) + t402) * MDP(6) + (pkin(1) * t350 * t398 + t377) * MDP(5) + t362 + (t278 * t332 + (-pkin(2) * t331 - t274) * t349 + ((-pkin(2) * t332 - t295) * qJD(3) - t369) * t354 + t393) * MDP(9) + (t277 * t332 + (t331 * t354 - t332 * t418) * pkin(2) + t363) * MDP(8) + (t302 * t239 + (-qJD(5) * t383 - t264 * t347 + t265 * t352) * t341 + t384 * t339 + t278 * t454 + t423 * t269 + t365) * MDP(22) + (t302 * t238 - (qJD(5) * t384 + t264 * t352 + t265 * t347) * t341 - t383 * t339 - t278 * t374 + t423 * t271 + t395) * MDP(23); (t326 * t239 + (-qJD(5) * t381 - t288 * t347 + t289 * t352) * t341 + t382 * t339 + t391 * t269 + t267 * t454 + t365) * MDP(22) + (t363 + t437) * MDP(8) + (t267 * t332 + t358) * MDP(9) + (t326 * t238 - (qJD(5) * t382 + t288 * t352 + t289 * t347) * t341 - t381 * t339 + t391 * t271 - t267 * t374 + t395) * MDP(23) + (t256 + t368 * t348 + (-t372 + t378) * t353) * MDP(15) + t364 + (t348 * t372 + t353 * t368 + t394) * MDP(16); MDP(12) * t434 + MDP(13) * t433 + qJDD(4) * MDP(14) + (-g(1) * t353 + t348 * t367) * MDP(15) + (g(1) * t348 + t353 * t367) * MDP(16) + (-(t249 * t347 - t439) * t341 + (-t269 * t432 + t339 * t352 - t341 * t414) * pkin(4) + t359) * MDP(22) + ((-qJD(5) * t248 - t249 * t341 - t229) * t352 + (-qJD(5) * t341 * t352 - t271 * t432 - t339 * t347) * pkin(4) + t360) * MDP(23) + t370 + (-MDP(10) * t428 + MDP(11) * t421) * t332 ^ 2; (-t341 * t387 + t359) * MDP(22) + ((-t229 + (-qJD(5) + t341) * t248) * t352 + t360) * MDP(23) + t370;];
tau = t1;
