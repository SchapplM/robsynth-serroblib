% Calculate vector of inverse dynamics joint torques for
% S6RPPRPR3
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% qJDD [6x1]
%   Generalized joint accelerations
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d4,d6,theta2,theta5]';
% MDP [23x1]
%   Minimal dynamic parameter vector (fixed base model)
%   see S6RPPRPR3_convert_par2_MPV_fixb.m
% 
% Output:
% tau [6x1]
%   joint torques of inverse dynamics (contains inertial, gravitational Coriolis and centrifugal forces)

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 01:45
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tau = S6RPPRPR3_invdynJ_fixb_mdp_slag_vp(qJ, qJD, qJDD, g, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(3,1),zeros(10,1),zeros(23,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPPRPR3_invdynJ_fixb_mdp_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RPPRPR3_invdynJ_fixb_mdp_slag_vp: qJD has to be [6x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [6 1]), ...
  'S6RPPRPR3_invdynJ_fixb_mdp_slag_vp: qJDD has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RPPRPR3_invdynJ_fixb_mdp_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RPPRPR3_invdynJ_fixb_mdp_slag_vp: pkin has to be [10x1] (double)');
assert(isreal(MDP) && all(size(MDP) == [23 1]), ...
  'S6RPPRPR3_invdynJ_fixb_mdp_slag_vp: MDP has to be [23x1] (double)'); 

%% Symbolic Calculation
% From invdyn_joint_fixb_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 01:45:04
% EndTime: 2019-03-09 01:45:09
% DurationCPUTime: 3.41s
% Computational Cost: add. (2236->332), mult. (4160->422), div. (0->0), fcn. (2765->14), ass. (0->153)
t327 = sin(qJ(4));
t330 = cos(qJ(4));
t408 = sin(pkin(10));
t409 = cos(pkin(10));
t279 = -t408 * t327 + t409 * t330;
t271 = t279 * qJD(1);
t326 = sin(qJ(6));
t329 = cos(qJ(6));
t256 = -t329 * qJD(4) + t271 * t326;
t338 = t327 * t409 + t330 * t408;
t269 = t338 * qJD(1);
t420 = qJD(6) + t269;
t424 = t256 * t420;
t258 = qJD(4) * t326 + t271 * t329;
t423 = t258 * t420;
t323 = sin(pkin(9));
t293 = pkin(1) * t323 + qJ(3);
t393 = qJDD(1) * t293;
t273 = t338 * qJD(4);
t363 = qJDD(1) * t408;
t364 = qJDD(1) * t409;
t249 = qJD(1) * t273 + t327 * t363 - t330 * t364;
t422 = -qJD(4) * qJD(6) + t249;
t317 = qJ(1) + pkin(9);
t308 = sin(t317);
t310 = cos(t317);
t370 = -g(1) * t308 + g(2) * t310;
t362 = t420 * t329;
t248 = -t271 * qJD(4) - t327 * t364 - t330 * t363;
t246 = -qJDD(6) + t248;
t403 = t246 * t326;
t421 = t362 * t420 - t403;
t324 = cos(pkin(9));
t301 = -pkin(1) * t324 - pkin(2);
t289 = -pkin(7) + t301;
t396 = qJ(5) - t289;
t274 = t396 * t327;
t282 = qJD(1) * t293;
t342 = -qJD(1) * t282 + t370;
t418 = qJDD(1) * t301;
t276 = qJD(1) * t289 + qJD(3);
t361 = qJ(5) * qJD(1) - t276;
t417 = -qJD(2) * t330 + t327 * t361;
t416 = 0.2e1 * qJD(4) * t282 + qJDD(4) * t289;
t275 = qJDD(1) * t289 + qJDD(3);
t265 = t330 * t275;
t347 = -qJ(5) * qJDD(1) - qJD(1) * qJD(5);
t226 = qJDD(4) * pkin(4) + t417 * qJD(4) - t327 * qJDD(2) + t347 * t330 + t265;
t231 = (-qJD(4) * t361 + qJDD(2)) * t330 + (-qJD(2) * qJD(4) + t275 + t347) * t327;
t215 = t226 * t409 - t231 * t408;
t213 = -qJDD(4) * pkin(5) - t215;
t267 = t330 * t276;
t389 = qJD(1) * t330;
t253 = -qJ(5) * t389 - qJD(2) * t327 + t267;
t410 = qJD(4) * pkin(4);
t252 = t253 + t410;
t367 = t408 * t417;
t229 = t252 * t409 + t367;
t224 = -qJD(4) * pkin(5) - t229;
t365 = t396 * t330;
t259 = -qJD(4) * t365 - qJD(5) * t327;
t339 = qJD(4) * t274 - qJD(5) * t330;
t235 = t259 * t409 + t339 * t408;
t314 = t327 * pkin(4);
t355 = t293 + t314;
t240 = pkin(5) * t338 - pkin(8) * t279 + t355;
t243 = -t274 * t409 - t365 * t408;
t216 = t408 * t226 + t409 * t231;
t214 = qJDD(4) * pkin(8) + t216;
t268 = qJD(1) * t314 + qJD(5) + t282;
t236 = pkin(5) * t269 - pkin(8) * t271 + t268;
t359 = qJD(6) * t236 + t214;
t415 = t213 * t279 - t224 * t273 + t243 * t246 - (qJD(6) * t240 + t235) * t420 - t338 * t359;
t299 = pkin(4) * t408 + pkin(8);
t316 = qJ(4) + pkin(10);
t307 = sin(t316);
t309 = cos(t316);
t414 = t420 * (pkin(4) * t389 + pkin(5) * t271 + pkin(8) * t269 + qJD(6) * t299) - t309 * t370 - g(3) * t307 + t213;
t411 = g(3) * t309;
t407 = t224 * t279;
t373 = t326 * qJDD(4) - t422 * t329;
t387 = qJD(6) * t326;
t227 = -t271 * t387 + t373;
t406 = t227 * t279;
t405 = t227 * t326;
t404 = t240 * t246;
t402 = t256 * t271;
t401 = t258 * t271;
t400 = t308 * t326;
t399 = t308 * t329;
t398 = t310 * t326;
t397 = t310 * t329;
t241 = t329 * t246;
t395 = qJDD(2) - g(3);
t270 = t279 * qJD(4);
t394 = t227 * t338 + t258 * t270;
t250 = t409 * t417;
t230 = t408 * t252 - t250;
t321 = t330 ^ 2;
t392 = t327 ^ 2 - t321;
t332 = qJD(4) ^ 2;
t333 = qJD(1) ^ 2;
t391 = -t332 - t333;
t390 = qJD(1) * t268;
t386 = qJD(6) * t329;
t385 = t330 * t410 + qJD(3);
t383 = qJD(1) * qJD(4);
t319 = qJD(3) * qJD(1);
t381 = qJDD(1) * t327;
t380 = qJDD(1) * t330;
t378 = qJDD(4) * t327;
t377 = qJDD(4) * t330;
t376 = t279 * t403;
t375 = t279 * t241;
t331 = cos(qJ(1));
t374 = t331 * pkin(1) + t310 * pkin(2) + t308 * qJ(3);
t277 = t319 + t393;
t372 = t330 * t383;
t328 = sin(qJ(1));
t371 = -pkin(1) * t328 + t310 * qJ(3);
t345 = qJDD(5) + t277 + (t372 + t381) * pkin(4);
t220 = -pkin(5) * t248 + pkin(8) * t249 + t345;
t225 = qJD(4) * pkin(8) + t230;
t358 = qJD(6) * t225 - t220;
t352 = -g(1) * t310 - g(2) * t308;
t350 = g(1) * t328 - g(2) * t331;
t312 = t329 * qJDD(4);
t228 = qJD(6) * t258 - t249 * t326 - t312;
t348 = -t228 * t338 - t256 * t270;
t346 = -t241 + (-t269 * t326 - t387) * t420;
t344 = t273 * t326 - t279 * t386;
t343 = t273 * t329 + t279 * t387;
t341 = t352 + t393;
t233 = t253 * t409 + t367;
t337 = t299 * t246 + (t224 + t233) * t420;
t336 = t215 * t279 + t216 * t338 - t229 * t273 + t230 * t270 + t370;
t335 = -t289 * t332 + t277 + t319 + t341;
t325 = -qJ(5) - pkin(7);
t300 = -pkin(4) * t409 - pkin(5);
t281 = -t327 * t332 + t377;
t280 = -t330 * t332 - t378;
t263 = t307 * t397 - t400;
t262 = t307 * t398 + t399;
t261 = t307 * t399 + t398;
t260 = -t307 * t400 + t397;
t242 = -t274 * t408 + t365 * t409;
t237 = pkin(5) * t270 + pkin(8) * t273 + t385;
t234 = t259 * t408 - t339 * t409;
t232 = t253 * t408 - t250;
t219 = t329 * t220;
t218 = t225 * t329 + t236 * t326;
t217 = -t225 * t326 + t236 * t329;
t1 = [qJDD(1) * MDP(1) + t350 * MDP(2) + (g(1) * t331 + g(2) * t328) * MDP(3) + (t350 + (t323 ^ 2 + t324 ^ 2) * qJDD(1) * pkin(1)) * pkin(1) * MDP(4) + (qJDD(3) + t370 + 0.2e1 * t418) * MDP(5) + (0.2e1 * t319 + t341 + t393) * MDP(6) + (t277 * t293 + t282 * qJD(3) + (qJDD(3) + t418) * t301 - g(1) * (-pkin(2) * t308 + t371) - g(2) * t374) * MDP(7) + (qJDD(1) * t321 - 0.2e1 * t327 * t372) * MDP(8) + 0.2e1 * (-t327 * t380 + t383 * t392) * MDP(9) + t281 * MDP(10) + t280 * MDP(11) + (t335 * t327 + t416 * t330) * MDP(13) + (-t416 * t327 + t335 * t330) * MDP(14) + (t234 * t271 - t235 * t269 - t242 * t249 + t243 * t248 - t336) * MDP(15) + (t216 * t243 + t230 * t235 - t215 * t242 - t229 * t234 + t345 * t355 + t268 * t385 - g(1) * (t310 * t314 + (-pkin(2) + t325) * t308 + t371) - g(2) * (t308 * t314 - t310 * t325 + t374)) * MDP(16) + (-t258 * t343 + t329 * t406) * MDP(17) + ((t256 * t329 + t258 * t326) * t273 + (-t405 - t228 * t329 + (t256 * t326 - t258 * t329) * qJD(6)) * t279) * MDP(18) + (-t343 * t420 - t375 + t394) * MDP(19) + (t344 * t420 + t348 + t376) * MDP(20) + (-t246 * t338 + t270 * t420) * MDP(21) + (-g(1) * t263 - g(2) * t261 + t217 * t270 + t219 * t338 + t242 * t228 + t234 * t256 + (t237 * t420 - t404 + (-t225 * t338 - t243 * t420 + t407) * qJD(6)) * t329 + t415 * t326) * MDP(22) + (g(1) * t262 - g(2) * t260 - t218 * t270 + t242 * t227 + t234 * t258 + (-(-qJD(6) * t243 + t237) * t420 + t404 + t358 * t338 - qJD(6) * t407) * t326 + t415 * t329) * MDP(23); t280 * MDP(13) - t281 * MDP(14) + (t248 * t279 - t249 * t338 + t269 * t273 + t270 * t271) * MDP(15) + (-t215 * t338 + t216 * t279 - t229 * t270 - t230 * t273 - g(3)) * MDP(16) + (-t348 + t376) * MDP(22) + (t375 + t394) * MDP(23) + (MDP(4) + MDP(7)) * t395 + (MDP(22) * t344 + MDP(23) * t343) * t420; -t333 * MDP(6) + (qJDD(3) + t342) * MDP(7) + (t327 * t391 + t377) * MDP(13) + (t330 * t391 - t378) * MDP(14) + (t248 * t338 + t249 * t279 - t269 * t270 + t271 * t273) * MDP(15) + (t336 - t390) * MDP(16) + (-t228 * t279 + t256 * t273 + t338 * t403) * MDP(22) + (t241 * t338 + t258 * t273 - t406) * MDP(23) + (MDP(7) * t301 + MDP(5)) * qJDD(1) + ((-qJD(1) * t329 - t270 * t326 - t338 * t386) * MDP(22) + (qJD(1) * t326 - t270 * t329 + t338 * t387) * MDP(23)) * t420; MDP(10) * t380 - MDP(11) * t381 + qJDD(4) * MDP(12) + (-t327 * t395 + t330 * t342 + t265) * MDP(13) + (t267 * qJD(4) + (-qJD(4) * t276 - t395) * t330 + (-t275 - t342) * t327) * MDP(14) + ((t230 - t232) * t271 - (t229 - t233) * t269 + (t248 * t408 + t249 * t409) * pkin(4)) * MDP(15) + (t229 * t232 - t230 * t233 + (t409 * t215 + t408 * t216 + g(3) * t327 + (t370 - t390) * t330) * pkin(4)) * MDP(16) + (t258 * t362 + t405) * MDP(17) + ((t227 - t424) * t329 + (-t228 - t423) * t326) * MDP(18) + (-t401 + t421) * MDP(19) + (t346 + t402) * MDP(20) - t420 * t271 * MDP(21) + (-t217 * t271 + t300 * t228 - t232 * t256 + t337 * t326 - t414 * t329) * MDP(22) + (t218 * t271 + t300 * t227 - t232 * t258 + t414 * t326 + t337 * t329) * MDP(23) + (MDP(8) * t327 * t330 - MDP(9) * t392) * t333; (-t269 ^ 2 - t271 ^ 2) * MDP(15) + (t229 * t271 + t230 * t269 + t345 + t352) * MDP(16) + (t346 - t402) * MDP(22) + (-t401 - t421) * MDP(23); t258 * t256 * MDP(17) + (-t256 ^ 2 + t258 ^ 2) * MDP(18) + (t373 + t424) * MDP(19) + (t312 + t423) * MDP(20) - t246 * MDP(21) + (-g(1) * t260 - g(2) * t262 + t218 * t420 - t224 * t258 + t219) * MDP(22) + (g(1) * t261 - g(2) * t263 + t217 * t420 + t224 * t256) * MDP(23) + ((-t214 + t411) * MDP(23) + (-MDP(20) * t271 - MDP(22) * t225 - MDP(23) * t236) * qJD(6)) * t329 + (-qJD(6) * t271 * MDP(19) + t422 * MDP(20) + (-t359 + t411) * MDP(22) + t358 * MDP(23)) * t326;];
tau  = t1;
