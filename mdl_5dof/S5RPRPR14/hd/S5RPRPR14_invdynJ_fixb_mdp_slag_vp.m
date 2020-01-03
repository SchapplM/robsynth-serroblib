% Calculate vector of inverse dynamics joint torques for
% S5RPRPR14
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
%   pkin=[a2,a3,a4,a5,d1,d3,d5,theta4]';
% MDP [22x1]
%   Minimal dynamic parameter vector (fixed base model)
%   see S5RPRPR14_convert_par2_MPV_fixb.m
% 
% Output:
% tau [5x1]
%   joint torques of inverse dynamics (contains inertial, gravitational Coriolis and centrifugal forces)

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 18:35
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tau = S5RPRPR14_invdynJ_fixb_mdp_slag_vp(qJ, qJD, qJDD, g, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(5,1),zeros(3,1),zeros(8,1),zeros(22,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPRPR14_invdynJ_fixb_mdp_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RPRPR14_invdynJ_fixb_mdp_slag_vp: qJD has to be [5x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [5 1]), ...
  'S5RPRPR14_invdynJ_fixb_mdp_slag_vp: qJDD has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RPRPR14_invdynJ_fixb_mdp_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RPRPR14_invdynJ_fixb_mdp_slag_vp: pkin has to be [8x1] (double)');
assert(isreal(MDP) && all(size(MDP) == [22 1]), ...
  'S5RPRPR14_invdynJ_fixb_mdp_slag_vp: MDP has to be [22x1] (double)'); 

%% Symbolic Calculation
% From invdyn_joint_fixb_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 18:35:23
% EndTime: 2019-12-31 18:35:27
% DurationCPUTime: 3.00s
% Computational Cost: add. (1606->295), mult. (3184->388), div. (0->0), fcn. (2122->10), ass. (0->140)
t297 = sin(qJ(3));
t300 = cos(qJ(3));
t373 = sin(pkin(8));
t374 = cos(pkin(8));
t251 = -t373 * t297 + t374 * t300;
t248 = t251 * qJD(1);
t296 = sin(qJ(5));
t299 = cos(qJ(5));
t229 = -t299 * qJD(3) + t248 * t296;
t309 = t374 * t297 + t373 * t300;
t381 = t309 * qJD(1);
t384 = qJD(5) + t381;
t388 = t229 * t384;
t231 = qJD(3) * t296 + t248 * t299;
t387 = t231 * t384;
t285 = t297 * pkin(3);
t359 = qJ(2) + t285;
t328 = qJDD(1) * t373;
t329 = qJDD(1) * t374;
t226 = qJD(3) * t381 + t297 * t328 - t300 * t329;
t386 = -qJD(3) * qJD(5) + t226;
t304 = qJD(1) ^ 2;
t298 = sin(qJ(1));
t301 = cos(qJ(1));
t382 = g(1) * t298 - g(2) * t301;
t312 = -qJ(2) * t304 - t382;
t327 = t384 * t299;
t225 = -qJD(3) * t248 - t297 * t329 - t300 * t328;
t223 = -qJDD(5) + t225;
t367 = t223 * t296;
t385 = t327 * t384 - t367;
t290 = qJDD(1) * qJ(2);
t318 = g(1) * t301 + g(2) * t298;
t291 = qJD(1) * qJD(2);
t339 = 0.2e1 * t291;
t380 = 0.2e1 * t290 + t339 - t318;
t302 = -pkin(1) - pkin(6);
t258 = t302 * qJDD(1) + qJDD(2);
t253 = t300 * t258;
t259 = t302 * qJD(1) + qJD(2);
t341 = qJDD(1) * t300;
t344 = qJD(1) * qJD(4);
t345 = qJD(1) * qJD(3);
t351 = qJD(3) * t297;
t209 = -t300 * t344 - t259 * t351 + qJDD(3) * pkin(3) + t253 + (t297 * t345 - t341) * qJ(4);
t326 = -qJ(4) * qJD(1) + t259;
t350 = qJD(3) * t300;
t215 = t326 * t350 + (-qJ(4) * qJDD(1) + t258 - t344) * t297;
t199 = t374 * t209 - t373 * t215;
t195 = -qJDD(3) * pkin(4) - t199;
t352 = qJD(1) * t300;
t239 = -qJ(4) * t352 + t300 * t259;
t236 = qJD(3) * pkin(3) + t239;
t238 = t326 * t297;
t334 = t373 * t238;
t210 = t374 * t236 - t334;
t206 = -qJD(3) * pkin(4) - t210;
t358 = qJ(4) - t302;
t330 = t358 * t300;
t237 = -qJD(3) * t330 - qJD(4) * t297;
t310 = -qJD(4) * t300 + t358 * t351;
t214 = t374 * t237 + t373 * t310;
t224 = pkin(4) * t309 - pkin(7) * t251 + t359;
t256 = t358 * t297;
t228 = -t374 * t256 - t373 * t330;
t331 = qJD(3) * t373;
t332 = qJD(3) * t374;
t247 = -t297 * t332 - t300 * t331;
t200 = t373 * t209 + t374 * t215;
t196 = qJDD(3) * pkin(7) + t200;
t255 = qJD(1) * t359 + qJD(4);
t212 = pkin(4) * t381 - pkin(7) * t248 + t255;
t325 = qJD(5) * t212 + t196;
t379 = t195 * t251 + t206 * t247 + t228 * t223 - (qJD(5) * t224 + t214) * t384 - t325 * t309;
t271 = t373 * pkin(3) + pkin(7);
t289 = qJ(3) + pkin(8);
t277 = sin(t289);
t278 = cos(t289);
t378 = (pkin(3) * t352 + pkin(4) * t248 + pkin(7) * t381 + qJD(5) * t271) * t384 + t382 * t278 - g(3) * t277 + t195;
t376 = g(3) * t278;
t375 = g(3) * t297;
t372 = pkin(1) * qJDD(1);
t338 = t296 * qJDD(3) - t299 * t386;
t349 = qJD(5) * t296;
t203 = -t248 * t349 + t338;
t370 = t203 * t251;
t369 = t203 * t296;
t368 = t206 * t251;
t366 = t224 * t223;
t365 = t229 * t248;
t364 = t231 * t248;
t363 = t296 * t298;
t362 = t296 * t301;
t361 = t298 * t299;
t220 = t299 * t223;
t360 = t299 * t301;
t234 = t374 * t238;
t211 = t373 * t236 + t234;
t357 = t301 * pkin(1) + t298 * qJ(2);
t294 = t300 ^ 2;
t355 = t297 ^ 2 - t294;
t303 = qJD(3) ^ 2;
t354 = -t303 - t304;
t353 = qJD(1) * t255;
t348 = qJD(5) * t299;
t347 = pkin(3) * t350 + qJD(2);
t342 = qJDD(1) * t297;
t340 = qJDD(3) * t297;
t337 = t300 * t345;
t316 = qJDD(4) + t290 + t291 + (t337 + t342) * pkin(3);
t202 = -pkin(4) * t225 + pkin(7) * t226 + t316;
t207 = qJD(3) * pkin(7) + t211;
t324 = qJD(5) * t207 - t202;
t321 = qJDD(2) - t372;
t315 = -t220 + (-t296 * t381 - t349) * t384;
t314 = t247 * t299 - t251 * t349;
t313 = 0.2e1 * qJ(2) * t345 + qJDD(3) * t302;
t217 = t374 * t239 - t334;
t308 = t271 * t223 + (t206 + t217) * t384;
t246 = t297 * t331 - t300 * t332;
t306 = t199 * t251 + t200 * t309 + t210 * t247 - t211 * t246 - t382;
t305 = -t302 * t303 + t380;
t295 = -qJ(4) - pkin(6);
t284 = t301 * qJ(2);
t281 = qJDD(3) * t300;
t280 = t299 * qJDD(3);
t272 = -t374 * pkin(3) - pkin(4);
t244 = t277 * t360 - t363;
t243 = t277 * t362 + t361;
t242 = t277 * t361 + t362;
t241 = -t277 * t363 + t360;
t227 = -t373 * t256 + t374 * t330;
t218 = -pkin(4) * t246 - pkin(7) * t247 + t347;
t216 = t373 * t239 + t234;
t213 = t373 * t237 - t374 * t310;
t204 = t231 * qJD(5) - t226 * t296 - t280;
t201 = t299 * t202;
t198 = t207 * t299 + t212 * t296;
t197 = -t207 * t296 + t212 * t299;
t1 = [qJDD(1) * MDP(1) + t382 * MDP(2) + t318 * MDP(3) + (qJDD(2) - t382 - 0.2e1 * t372) * MDP(4) + t380 * MDP(5) + (-t321 * pkin(1) - g(1) * (-pkin(1) * t298 + t284) - g(2) * t357 + (t339 + t290) * qJ(2)) * MDP(6) + (qJDD(1) * t294 - 0.2e1 * t297 * t337) * MDP(7) + 0.2e1 * (-t297 * t341 + t355 * t345) * MDP(8) + (-t297 * t303 + t281) * MDP(9) + (-t300 * t303 - t340) * MDP(10) + (t305 * t297 + t313 * t300) * MDP(12) + (-t313 * t297 + t305 * t300) * MDP(13) + (t213 * t248 - t214 * t381 + t225 * t228 - t226 * t227 - t306) * MDP(14) + (t200 * t228 + t211 * t214 - t199 * t227 - t210 * t213 + t316 * t359 + t255 * t347 - g(1) * (t301 * t285 + t284 + (-pkin(1) + t295) * t298) - g(2) * (t298 * t285 - t295 * t301 + t357)) * MDP(15) + (t314 * t231 + t299 * t370) * MDP(16) + ((-t229 * t299 - t231 * t296) * t247 + (-t369 - t204 * t299 + (t229 * t296 - t231 * t299) * qJD(5)) * t251) * MDP(17) + (t203 * t309 - t251 * t220 - t231 * t246 + t314 * t384) * MDP(18) + (t251 * t367 - t204 * t309 + t229 * t246 + (-t247 * t296 - t251 * t348) * t384) * MDP(19) + (-t223 * t309 - t246 * t384) * MDP(20) + (-g(1) * t244 - g(2) * t242 - t197 * t246 + t201 * t309 + t227 * t204 + t213 * t229 + (t218 * t384 - t366 + (-t207 * t309 - t228 * t384 + t368) * qJD(5)) * t299 + t379 * t296) * MDP(21) + (g(1) * t243 - g(2) * t241 + t198 * t246 + t227 * t203 + t213 * t231 + (-(-qJD(5) * t228 + t218) * t384 + t366 + t324 * t309 - qJD(5) * t368) * t296 + t379 * t299) * MDP(22); qJDD(1) * MDP(4) - t304 * MDP(5) + (t321 + t312) * MDP(6) + (t354 * t297 + t281) * MDP(12) + (t354 * t300 - t340) * MDP(13) + (t225 * t309 + t226 * t251 + t246 * t381 - t247 * t248) * MDP(14) + (t306 - t353) * MDP(15) + (-t204 * t251 - t229 * t247 + t309 * t367) * MDP(21) + (t220 * t309 - t231 * t247 - t370) * MDP(22) + ((-qJD(1) * t299 + t246 * t296 - t309 * t348) * MDP(21) + (qJD(1) * t296 + t246 * t299 + t309 * t349) * MDP(22)) * t384; MDP(9) * t341 - MDP(10) * t342 + qJDD(3) * MDP(11) + (t312 * t300 + t253 + t375) * MDP(12) + (g(3) * t300 + (-t258 - t312) * t297) * MDP(13) + ((t211 - t216) * t248 - (t210 - t217) * t381 + (t373 * t225 + t374 * t226) * pkin(3)) * MDP(14) + (t210 * t216 - t211 * t217 + (t374 * t199 + t373 * t200 + t375 + (-t382 - t353) * t300) * pkin(3)) * MDP(15) + (t231 * t327 + t369) * MDP(16) + ((t203 - t388) * t299 + (-t204 - t387) * t296) * MDP(17) + (-t364 + t385) * MDP(18) + (t315 + t365) * MDP(19) - t384 * t248 * MDP(20) + (-t197 * t248 + t272 * t204 - t216 * t229 + t308 * t296 - t378 * t299) * MDP(21) + (t198 * t248 + t272 * t203 - t216 * t231 + t378 * t296 + t308 * t299) * MDP(22) + (t300 * t297 * MDP(7) - t355 * MDP(8)) * t304; (-t248 ^ 2 - t381 ^ 2) * MDP(14) + (t210 * t248 + t211 * t381 + t316 - t318) * MDP(15) + (t315 - t365) * MDP(21) + (-t364 - t385) * MDP(22); t231 * t229 * MDP(16) + (-t229 ^ 2 + t231 ^ 2) * MDP(17) + (t338 + t388) * MDP(18) + (t280 + t387) * MDP(19) - t223 * MDP(20) + (-g(1) * t241 - g(2) * t243 + t198 * t384 - t206 * t231 + t201) * MDP(21) + (g(1) * t242 - g(2) * t244 + t197 * t384 + t206 * t229) * MDP(22) + ((-t196 + t376) * MDP(22) + (-MDP(19) * t248 - MDP(21) * t207 - MDP(22) * t212) * qJD(5)) * t299 + (-qJD(5) * t248 * MDP(18) + t386 * MDP(19) + (-t325 + t376) * MDP(21) + t324 * MDP(22)) * t296;];
tau = t1;
