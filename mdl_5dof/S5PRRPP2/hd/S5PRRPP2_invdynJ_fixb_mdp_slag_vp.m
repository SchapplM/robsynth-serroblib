% Calculate vector of inverse dynamics joint torques for
% S5PRRPP2
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
%   pkin=[a2,a3,a4,a5,d2,d3,theta1,theta4]';
% MDP [19x1]
%   Minimal dynamic parameter vector (fixed base model)
%   see S5PRRPP2_convert_par2_MPV_fixb.m
% 
% Output:
% tau [5x1]
%   joint torques of inverse dynamics (contains inertial, gravitational Coriolis and centrifugal forces)

% Quelle: HybrDyn-Toolbox
% Datum: 2021-01-15 15:33
% Revision: 24b2e7d74a0c1a3b64fa2f8f5ad758691ad61af3 (2021-01-15)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tau = S5PRRPP2_invdynJ_fixb_mdp_slag_vp(qJ, qJD, qJDD, g, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(5,1),zeros(3,1),zeros(8,1),zeros(19,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PRRPP2_invdynJ_fixb_mdp_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5PRRPP2_invdynJ_fixb_mdp_slag_vp: qJD has to be [5x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [5 1]), ...
  'S5PRRPP2_invdynJ_fixb_mdp_slag_vp: qJDD has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5PRRPP2_invdynJ_fixb_mdp_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5PRRPP2_invdynJ_fixb_mdp_slag_vp: pkin has to be [8x1] (double)');
assert(isreal(MDP) && all(size(MDP) == [19 1]), ...
  'S5PRRPP2_invdynJ_fixb_mdp_slag_vp: MDP has to be [19x1] (double)'); 

%% Symbolic Calculation
% From invdyn_joint_fixb_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2021-01-15 15:32:40
% EndTime: 2021-01-15 15:32:48
% DurationCPUTime: 3.55s
% Computational Cost: add. (1581->317), mult. (3441->395), div. (0->0), fcn. (2375->10), ass. (0->142)
t299 = qJ(4) + pkin(6);
t302 = cos(qJ(3));
t268 = t299 * t302;
t295 = sin(pkin(8));
t297 = cos(pkin(8));
t300 = sin(qJ(3));
t340 = t299 * t300;
t231 = t297 * t268 - t295 * t340;
t291 = qJ(3) + pkin(8);
t288 = sin(t291);
t303 = cos(qJ(2));
t290 = g(3) * t303;
t301 = sin(qJ(2));
t296 = sin(pkin(7));
t298 = cos(pkin(7));
t333 = g(1) * t298 + g(2) * t296;
t325 = t333 * t301;
t393 = t290 - t325;
t394 = -qJDD(3) * t231 + t393 * t288;
t337 = qJD(3) * t299;
t248 = qJD(4) * t302 - t300 * t337;
t316 = -qJD(4) * t300 - t302 * t337;
t370 = t297 * t302;
t375 = t295 * t300;
t260 = -t370 + t375;
t321 = t260 * t303;
t361 = qJD(1) * t321 + t297 * t248 + t295 * t316;
t324 = t333 * t303;
t382 = g(3) * t301;
t312 = t324 + t382;
t261 = t295 * t302 + t297 * t300;
t254 = t261 * qJD(2);
t385 = pkin(3) * t302;
t287 = pkin(2) + t385;
t331 = t287 * qJDD(2);
t392 = MDP(14) + MDP(17);
t354 = qJD(3) * t300;
t389 = -qJD(3) * t370 + t295 * t354;
t253 = t261 * qJD(3);
t346 = MDP(12) + MDP(16);
t345 = MDP(13) - MDP(18);
t304 = qJD(3) ^ 2;
t350 = qJDD(1) * t303;
t352 = qJD(1) * qJD(2);
t334 = t301 * t352 - t350;
t387 = 0.2e1 * qJDD(2) * pkin(2) - pkin(6) * t304 + t301 * (t333 + t352) - t290 - t334;
t250 = t254 ^ 2;
t386 = pkin(3) * t300;
t381 = qJD(2) * pkin(2);
t380 = qJ(5) * t288;
t378 = qJDD(3) * pkin(4);
t359 = qJD(1) * t301;
t336 = t299 * qJD(2) + t359;
t247 = t336 * t302;
t377 = t247 * t295;
t289 = cos(t291);
t376 = t289 * t303;
t373 = t296 * t301;
t372 = t296 * t303;
t233 = t297 * t247;
t369 = t298 * t301;
t368 = t298 * t302;
t367 = t298 * t303;
t366 = t299 * t303;
t365 = t300 * t303;
t305 = qJD(2) ^ 2;
t364 = t302 * t305;
t363 = qJDD(1) - g(3);
t259 = qJDD(2) * pkin(6) + qJDD(1) * t301 + t303 * t352;
t320 = qJ(4) * qJDD(2) + qJD(2) * qJD(4) + t259;
t330 = qJD(3) * t336;
t210 = qJDD(3) * pkin(3) - t300 * t320 - t302 * t330;
t214 = -t300 * t330 + t302 * t320;
t200 = t297 * t210 - t295 * t214;
t201 = t295 * t210 + t297 * t214;
t358 = qJD(1) * t303;
t362 = t248 * t295 - t261 * t358 - t297 * t316;
t246 = t336 * t300;
t235 = qJD(3) * pkin(3) - t246;
t218 = t295 * t235 + t233;
t293 = t300 ^ 2;
t360 = -t302 ^ 2 + t293;
t357 = qJD(2) * t300;
t219 = -t246 * t295 + t233;
t355 = qJD(3) * t219;
t220 = -t246 * t297 - t377;
t353 = qJD(5) - t220;
t351 = qJD(2) * qJD(3);
t349 = qJDD(2) * t300;
t348 = qJDD(2) * t302;
t344 = pkin(3) * t354;
t341 = qJD(2) * t370;
t339 = t300 * t351;
t338 = t302 * t351;
t207 = pkin(4) * t253 + qJ(5) * t389 - qJD(5) * t261 + t344;
t335 = -t207 + t359;
t332 = g(1) * t296 - g(2) * t298;
t217 = t235 * t297 - t377;
t329 = qJDD(3) * t302 - t300 * t304;
t328 = qJDD(3) * t300 + t302 * t304;
t323 = pkin(3) * t339 + qJDD(4) + t334;
t240 = t288 * t372 + t289 * t298;
t242 = t288 * t367 - t296 * t289;
t317 = g(1) * t242 + g(2) * t240 + t288 * t382 + t200;
t230 = t268 * t295 + t297 * t340;
t314 = -g(3) * t376 - qJDD(3) * t230 + (g(1) * t369 + g(2) * t373) * t289;
t313 = qJDD(2) * t261 - t295 * t339;
t257 = -qJD(2) * t287 + qJD(4) - t358;
t241 = -t298 * t288 + t289 * t372;
t243 = t288 * t296 + t289 * t367;
t311 = g(1) * t243 + g(2) * t241 + t289 * t382 - t201;
t277 = -t358 - t381;
t310 = -pkin(6) * qJDD(3) + (t277 + t358 - t381) * qJD(3);
t251 = t295 * t357 - t341;
t215 = pkin(4) * t251 - qJ(5) * t254 + t257;
t309 = -t215 * t254 - qJDD(5) + t317;
t308 = -t277 * qJD(2) - t259 + t312;
t278 = t297 * t348;
t227 = qJD(2) * t253 + t295 * t349 - t278;
t228 = t297 * t338 + t313;
t307 = pkin(4) * t227 - qJ(5) * t228 - qJD(5) * t254 + t323;
t306 = -t231 * t227 + t228 * t230 - t361 * t251 + t254 * t362 - t312;
t292 = qJDD(3) * qJ(5);
t285 = -pkin(3) * t297 - pkin(4);
t283 = pkin(3) * t295 + qJ(5);
t281 = t296 * t385;
t273 = t303 * t287;
t271 = t298 * t366;
t270 = t296 * t366;
t245 = t260 * t301;
t244 = t261 * t301;
t229 = -t331 + t323;
t226 = pkin(4) * t260 - qJ(5) * t261 - t287;
t223 = pkin(3) * t357 + pkin(4) * t254 + qJ(5) * t251;
t222 = -qJD(2) * t321 - t301 * t253;
t221 = -t303 * t254 + t389 * t301;
t216 = qJD(3) * qJ(5) + t218;
t213 = -qJD(3) * pkin(4) + qJD(5) - t217;
t199 = qJDD(5) - t200 - t378;
t198 = -t331 + t307;
t197 = qJD(3) * qJD(5) + t201 + t292;
t1 = [t363 * MDP(1) + (-t200 * t244 - t201 * t245 + t217 * t221 + t218 * t222 - g(3)) * MDP(15) + (-t197 * t245 + t199 * t244 - t213 * t221 + t216 * t222 - g(3)) * MDP(19) + (qJDD(2) * MDP(3) - t305 * MDP(4) + (-0.2e1 * t339 + t348) * MDP(10) + (-0.2e1 * t338 - t349) * MDP(11) - t229 * MDP(15) - t198 * MDP(19) - t345 * t228 - t346 * t227) * t303 + (-t305 * MDP(3) - qJDD(2) * MDP(4) + (-t328 - t364) * MDP(10) + (t300 * t305 - t329) * MDP(11) + (t257 * MDP(15) + t215 * MDP(19) + t254 * t345) * qJD(2)) * t301 - t345 * (qJD(3) * t222 - qJDD(3) * t245) + t346 * (qJD(2) * t301 * t251 + qJD(3) * t221 - qJDD(3) * t244) + t392 * (-t221 * t254 - t222 * t251 + t245 * t227 + t228 * t244); qJDD(2) * MDP(2) + (t350 - t393) * MDP(3) + (-t301 * t363 + t324) * MDP(4) + (qJDD(2) * t293 + 0.2e1 * t300 * t338) * MDP(5) + 0.2e1 * (t300 * t348 - t351 * t360) * MDP(6) + t328 * MDP(7) + t329 * MDP(8) + (t310 * t300 + t387 * t302) * MDP(10) + (-t387 * t300 + t310 * t302) * MDP(11) + (-t251 * t359 - t227 * t287 + t229 * t260 + t253 * t257 + (t251 * t386 - t362) * qJD(3) + t314) * MDP(12) + (-t254 * t359 - t228 * t287 + t229 * t261 - t389 * t257 + (t254 * t386 - t361) * qJD(3) + t394) * MDP(13) + (-t200 * t261 - t201 * t260 + t217 * t389 - t218 * t253 + t306) * MDP(14) + (t201 * t231 - t200 * t230 - t229 * t287 - g(1) * (-t287 * t369 + t271) - g(2) * (-t287 * t373 + t270) - g(3) * (t299 * t301 + t273) + (t344 - t359) * t257 + t361 * t218 - t362 * t217) * MDP(15) + (-qJD(3) * t362 + t198 * t260 + t215 * t253 + t226 * t227 - t251 * t335 + t314) * MDP(16) + (-t197 * t260 + t199 * t261 - t213 * t389 - t216 * t253 + t306) * MDP(17) + (t361 * qJD(3) - t198 * t261 + t215 * t389 - t226 * t228 + t335 * t254 - t394) * MDP(18) + (t197 * t231 + t198 * t226 + t215 * t207 + t199 * t230 - g(1) * t271 - g(2) * t270 - g(3) * (pkin(4) * t376 + t303 * t380 + t273) + t361 * t216 + t362 * t213 + (-g(3) * t299 - t215 * qJD(1) + t333 * (pkin(4) * t289 + t287 + t380)) * t301) * MDP(19); -t300 * MDP(5) * t364 + t360 * MDP(6) * t305 + MDP(7) * t349 + MDP(8) * t348 + qJDD(3) * MDP(9) + (t300 * t308 - t302 * t332) * MDP(10) + (t300 * t332 + t302 * t308) * MDP(11) + (t355 - t254 * t257 + (qJDD(3) * t297 - t251 * t357) * pkin(3) + t317) * MDP(12) + (qJD(3) * t220 + t251 * t257 + (-qJDD(3) * t295 - t254 * t357) * pkin(3) + t311) * MDP(13) + ((t218 - t219) * t254 + (-t217 + t220) * t251 + (-t227 * t295 - t228 * t297) * pkin(3)) * MDP(14) + (-g(1) * t281 + t217 * t219 - t218 * t220 + (g(2) * t368 + t200 * t297 + t201 * t295 + (-t257 * qJD(2) + t312) * t300) * pkin(3)) * MDP(15) + (t355 - t223 * t251 + (pkin(4) - t285) * qJDD(3) + t309) * MDP(16) + (-t227 * t283 + t228 * t285 + (t216 - t219) * t254 + (t213 - t353) * t251) * MDP(17) + (qJDD(3) * t283 - t215 * t251 + t223 * t254 + t292 + (0.2e1 * qJD(5) - t220) * qJD(3) - t311) * MDP(18) + (t197 * t283 + t199 * t285 - t215 * t223 - t213 * t219 - g(1) * (-pkin(3) * t298 * t365 - pkin(4) * t242 + qJ(5) * t243 + t281) - g(2) * (-pkin(4) * t240 + qJ(5) * t241 + (-t296 * t365 - t368) * pkin(3)) + t353 * t216 - (-pkin(4) * t288 + qJ(5) * t289 - t386) * t382) * MDP(19); (t217 * t254 + t218 * t251 + t290 + t323) * MDP(15) + (-t213 * t254 + t216 * t251 + t290 + t307) * MDP(19) + t345 * ((-t251 + t341) * qJD(3) + t313) + t392 * (-t251 ^ 2 - t250) + (t325 + t331) * (-MDP(15) - MDP(19)) + (0.2e1 * qJD(3) * t254 + t375 * qJDD(2) - t278) * t346; (t251 * t254 - qJDD(3)) * MDP(16) + ((t251 + t341) * qJD(3) + t313) * MDP(17) + (-t250 - t304) * MDP(18) + (-qJD(3) * t216 - t309 - t378) * MDP(19);];
tau = t1;
