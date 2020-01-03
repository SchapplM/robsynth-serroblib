% Calculate vector of inverse dynamics joint torques for
% S5RRPPR5
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
%   pkin=[a2,a3,a4,a5,d1,d2,d5,theta3]';
% MDP [23x1]
%   Minimal dynamic parameter vector (fixed base model)
%   see S5RRPPR5_convert_par2_MPV_fixb.m
% 
% Output:
% tau [5x1]
%   joint torques of inverse dynamics (contains inertial, gravitational Coriolis and centrifugal forces)

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 19:30
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tau = S5RRPPR5_invdynJ_fixb_mdp_slag_vp(qJ, qJD, qJDD, g, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(5,1),zeros(3,1),zeros(8,1),zeros(23,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRPPR5_invdynJ_fixb_mdp_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RRPPR5_invdynJ_fixb_mdp_slag_vp: qJD has to be [5x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [5 1]), ...
  'S5RRPPR5_invdynJ_fixb_mdp_slag_vp: qJDD has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RRPPR5_invdynJ_fixb_mdp_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RRPPR5_invdynJ_fixb_mdp_slag_vp: pkin has to be [8x1] (double)');
assert(isreal(MDP) && all(size(MDP) == [23 1]), ...
  'S5RRPPR5_invdynJ_fixb_mdp_slag_vp: MDP has to be [23x1] (double)'); 

%% Symbolic Calculation
% From invdyn_joint_fixb_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 19:30:04
% EndTime: 2019-12-31 19:30:10
% DurationCPUTime: 3.43s
% Computational Cost: add. (1761->332), mult. (4006->417), div. (0->0), fcn. (2831->10), ass. (0->148)
t348 = sin(pkin(8));
t349 = cos(pkin(8));
t352 = sin(qJ(2));
t355 = cos(qJ(2));
t315 = t348 * t355 + t349 * t352;
t304 = t315 * qJD(2);
t399 = qJDD(1) * t355;
t327 = t349 * t399;
t400 = qJDD(1) * t352;
t273 = qJD(1) * t304 + t348 * t400 - t327;
t401 = qJD(1) * qJD(2);
t393 = t352 * t401;
t325 = t348 * t393;
t366 = t315 * qJDD(1) - t325;
t392 = t355 * t401;
t274 = t349 * t392 + t366;
t351 = sin(qJ(5));
t354 = cos(qJ(5));
t406 = qJD(1) * t355;
t394 = t349 * t406;
t407 = qJD(1) * t352;
t302 = t348 * t407 - t394;
t305 = t315 * qJD(1);
t442 = t354 * t302 - t305 * t351;
t227 = t442 * qJD(5) + t351 * t273 + t354 * t274;
t340 = qJDD(2) - qJDD(5);
t341 = qJD(2) - qJD(5);
t420 = t442 * t341;
t435 = t302 * t351 + t354 * t305;
t445 = -MDP(17) * t442 * t435 + (t435 ^ 2 - t442 ^ 2) * MDP(18) + (t227 + t420) * MDP(19) - t340 * MDP(21);
t421 = t435 * t341;
t338 = t355 * pkin(2);
t438 = t338 + pkin(1);
t322 = -qJD(1) * t438 + qJD(3);
t443 = -qJ(4) * t305 + t322;
t350 = -qJ(3) - pkin(6);
t389 = qJD(2) * t350;
t298 = -qJD(3) * t352 + t355 * t389;
t323 = t350 * t352;
t269 = qJDD(2) * pkin(2) + t298 * qJD(1) + qJDD(1) * t323;
t297 = qJD(3) * t355 + t352 * t389;
t324 = t350 * t355;
t277 = t297 * qJD(1) - qJDD(1) * t324;
t232 = t349 * t269 - t348 * t277;
t391 = -qJDD(4) + t232;
t430 = -pkin(3) - pkin(4);
t225 = -pkin(7) * t274 + t430 * qJDD(2) - t391;
t344 = qJD(2) * qJD(4);
t233 = t348 * t269 + t349 * t277;
t396 = qJDD(2) * qJ(4) + t233;
t230 = t344 + t396;
t226 = pkin(7) * t273 + t230;
t237 = t430 * t302 - t443;
t353 = sin(qJ(1));
t342 = qJ(2) + pkin(8);
t336 = sin(t342);
t337 = cos(t342);
t434 = -t336 * t354 + t337 * t351;
t282 = t434 * t353;
t356 = cos(qJ(1));
t284 = t434 * t356;
t372 = t336 * t351 + t337 * t354;
t440 = -g(1) * t284 - g(2) * t282 - g(3) * t372 - t354 * t225 + t351 * t226 + t237 * t435;
t436 = g(1) * t353 - g(2) * t356;
t386 = g(1) * t356 + g(2) * t353;
t319 = qJD(1) * t323;
t320 = qJD(1) * t324;
t415 = t348 * t320;
t279 = t349 * t319 + t415;
t402 = qJD(4) - t279;
t383 = t438 * qJDD(1);
t283 = t372 * t353;
t285 = t372 * t356;
t433 = -g(1) * t285 - g(2) * t283 + g(3) * t434 + t351 * t225 + t354 * t226 + t237 * t442;
t249 = pkin(3) * t302 + t443;
t432 = -g(3) * t337 - t249 * t305 + t386 * t336 + t391;
t388 = -t354 * t273 + t274 * t351;
t228 = qJD(5) * t435 + t388;
t300 = t305 ^ 2;
t429 = pkin(7) * t302;
t428 = pkin(7) * t305;
t424 = g(3) * t355;
t422 = qJDD(2) * pkin(3);
t414 = t348 * t352;
t413 = t349 * t320;
t412 = t349 * t355;
t252 = t349 * t297 + t348 * t298;
t313 = qJD(2) * pkin(2) + t319;
t272 = t348 * t313 - t413;
t281 = t348 * t323 - t349 * t324;
t346 = t352 ^ 2;
t409 = -t355 ^ 2 + t346;
t405 = qJD(2) * t352;
t403 = -t428 + t402;
t398 = pkin(2) * t393 + qJDD(3);
t397 = pkin(2) * t405;
t263 = qJD(2) * qJ(4) + t272;
t334 = -pkin(2) * t349 - pkin(3);
t251 = t297 * t348 - t349 * t298;
t271 = t313 * t349 + t415;
t278 = t319 * t348 - t413;
t280 = -t349 * t323 - t324 * t348;
t387 = qJ(4) * t315 + t438;
t384 = qJD(4) - t271;
t382 = pkin(3) * t337 + qJ(4) * t336;
t239 = t430 * qJD(2) + t384 - t428;
t245 = t263 + t429;
t381 = t354 * t239 - t351 * t245;
t380 = -t351 * t239 - t354 * t245;
t253 = -pkin(7) * t315 + t280;
t314 = -t412 + t414;
t254 = pkin(7) * t314 + t281;
t379 = t253 * t354 - t254 * t351;
t378 = t253 * t351 + t254 * t354;
t375 = t354 * t314 - t315 * t351;
t276 = t314 * t351 + t315 * t354;
t330 = -pkin(4) + t334;
t332 = pkin(2) * t348 + qJ(4);
t374 = t330 * t354 - t332 * t351;
t373 = t330 * t351 + t332 * t354;
t371 = -pkin(2) * t407 - qJ(4) * t302;
t370 = -0.2e1 * pkin(1) * t401 - pkin(6) * qJDD(2);
t368 = qJ(4) * t274 + qJD(4) * t305 - t398;
t307 = qJD(2) * t412 - t348 * t405;
t365 = qJ(4) * t307 + qJD(4) * t315 - t397;
t357 = qJD(2) ^ 2;
t364 = 0.2e1 * qJDD(1) * pkin(1) - pkin(6) * t357 + t436;
t358 = qJD(1) ^ 2;
t363 = pkin(1) * t358 - pkin(6) * qJDD(1) + t386;
t362 = pkin(3) * t273 - t368;
t361 = t251 * t305 - t252 * t302 - t281 * t273 + t274 * t280 - t386;
t326 = t356 * t438;
t270 = pkin(3) * t314 - t387;
t255 = -qJD(2) * pkin(3) + t384;
t250 = pkin(3) * t305 - t371;
t247 = t278 + t429;
t246 = t430 * t314 + t387;
t244 = pkin(3) * t304 - t365;
t242 = pkin(7) * t304 + t252;
t241 = -pkin(7) * t307 + t251;
t240 = t430 * t305 + t371;
t236 = t430 * t304 + t365;
t235 = t276 * qJD(5) - t354 * t304 + t307 * t351;
t234 = t375 * qJD(5) + t304 * t351 + t307 * t354;
t231 = -t391 - t422;
t229 = -t383 + t362;
t224 = t430 * t273 + t368 + t383;
t1 = [qJDD(1) * MDP(1) + t436 * MDP(2) + t386 * MDP(3) + (qJDD(1) * t346 + 0.2e1 * t352 * t392) * MDP(4) + 0.2e1 * (t352 * t399 - t409 * t401) * MDP(5) + (qJDD(2) * t352 + t355 * t357) * MDP(6) + (qJDD(2) * t355 - t352 * t357) * MDP(7) + (t370 * t352 + t364 * t355) * MDP(9) + (-t364 * t352 + t370 * t355) * MDP(10) + (-t232 * t315 - t233 * t314 - t271 * t307 - t272 * t304 + t361) * MDP(11) + (t233 * t281 + t272 * t252 - t232 * t280 - t271 * t251 - (-t383 + t398) * t438 + t322 * t397 - g(1) * (-t350 * t356 - t353 * t438) - g(2) * (-t350 * t353 + t326)) * MDP(12) + (-qJD(2) * t251 - qJDD(2) * t280 + t229 * t314 + t244 * t302 + t249 * t304 + t270 * t273 + t337 * t436) * MDP(13) + (-t230 * t314 + t231 * t315 + t255 * t307 - t263 * t304 + t361) * MDP(14) + (qJD(2) * t252 + qJDD(2) * t281 - t229 * t315 - t244 * t305 - t249 * t307 - t270 * t274 + t336 * t436) * MDP(15) + (-g(2) * t326 + t229 * t270 + t230 * t281 + t231 * t280 + t249 * t244 + t255 * t251 + t263 * t252 + (g(1) * t350 - g(2) * t382) * t356 + (-g(1) * (-t438 - t382) + g(2) * t350) * t353) * MDP(16) + (t227 * t276 + t234 * t435) * MDP(17) + (t227 * t375 - t228 * t276 + t234 * t442 - t235 * t435) * MDP(18) + (-t234 * t341 - t276 * t340) * MDP(19) + (t235 * t341 - t340 * t375) * MDP(20) + (-t236 * t442 + t246 * t228 - t224 * t375 + t237 * t235 - (-t378 * qJD(5) + t241 * t354 - t242 * t351) * t341 - t379 * t340 + g(1) * t283 - g(2) * t285) * MDP(22) + (t236 * t435 + t246 * t227 + t224 * t276 + t237 * t234 + (t379 * qJD(5) + t241 * t351 + t242 * t354) * t341 + t378 * t340 - g(1) * t282 + g(2) * t284) * MDP(23); MDP(6) * t400 + MDP(7) * t399 + qJDD(2) * MDP(8) + (t363 * t352 - t424) * MDP(9) + (g(3) * t352 + t363 * t355) * MDP(10) + ((t272 - t278) * t305 + (-t271 + t279) * t302 + (-t273 * t348 - t274 * t349) * pkin(2)) * MDP(11) + (t271 * t278 - t272 * t279 + (-t424 + t232 * t349 + t233 * t348 + (-qJD(1) * t322 + t386) * t352) * pkin(2)) * MDP(12) + (qJD(2) * t278 - t250 * t302 + (pkin(3) - t334) * qJDD(2) + t432) * MDP(13) + (-t273 * t332 + t274 * t334 + (t263 - t278) * t305 + (t255 - t402) * t302) * MDP(14) + (-g(3) * t336 - qJD(2) * t279 + qJDD(2) * t332 - t249 * t302 + t250 * t305 - t386 * t337 + 0.2e1 * t344 + t396) * MDP(15) + (t230 * t332 + t231 * t334 - t249 * t250 - t255 * t278 - g(3) * (t338 + t382) + t402 * t263 + t386 * (pkin(2) * t352 + pkin(3) * t336 - qJ(4) * t337)) * MDP(16) + (t228 + t421) * MDP(20) + (-t374 * t340 + t240 * t442 + (t354 * t247 + t403 * t351) * t341 + (t373 * t341 - t380) * qJD(5) + t440) * MDP(22) + (t373 * t340 - t240 * t435 + (-t351 * t247 + t403 * t354) * t341 + (t374 * t341 + t381) * qJD(5) + t433) * MDP(23) + (-t352 * t355 * MDP(4) + t409 * MDP(5)) * t358 - t445; (t271 * t305 + t272 * t302 + t398 - t436) * MDP(12) - t327 * MDP(13) + t325 * MDP(15) + (-t255 * t305 + t263 * t302 + t362 - t436) * MDP(16) + (-t228 + t421) * MDP(22) + (-t227 + t420) * MDP(23) + (MDP(13) * t414 - t315 * MDP(15) + (-MDP(12) - MDP(16)) * t438) * qJDD(1) + ((t348 * t406 + t349 * t407 + t305) * MDP(13) + (t302 - t394) * MDP(15)) * qJD(2) + (MDP(11) + MDP(14)) * (-t302 ^ 2 - t300); (t302 * t305 - qJDD(2)) * MDP(13) + ((t302 + t394) * qJD(2) + t366) * MDP(14) + (-t300 - t357) * MDP(15) + (-qJD(2) * t263 - t422 - t432) * MDP(16) + (t305 * t442 - t354 * t340) * MDP(22) + (-t305 * t435 + t351 * t340) * MDP(23) + (-MDP(22) * t351 - MDP(23) * t354) * t341 ^ 2; (-t388 - t421) * MDP(20) + (t380 * t341 - t440) * MDP(22) + (-t381 * t341 - t433) * MDP(23) + (-MDP(20) * t435 + t380 * MDP(22) - t381 * MDP(23)) * qJD(5) + t445;];
tau = t1;
