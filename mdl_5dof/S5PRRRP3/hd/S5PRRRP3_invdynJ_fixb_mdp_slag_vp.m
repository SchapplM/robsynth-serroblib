% Calculate vector of inverse dynamics joint torques for
% S5PRRRP3
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
%   pkin=[a2,a3,a4,a5,d2,d3,d4,theta1]';
% MDP [22x1]
%   Minimal dynamic parameter vector (fixed base model)
%   see S5PRRRP3_convert_par2_MPV_fixb.m
% 
% Output:
% tau [5x1]
%   joint torques of inverse dynamics (contains inertial, gravitational Coriolis and centrifugal forces)

% Quelle: HybrDyn-Toolbox
% Datum: 2021-01-15 16:24
% Revision: 24b2e7d74a0c1a3b64fa2f8f5ad758691ad61af3 (2021-01-15)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tau = S5PRRRP3_invdynJ_fixb_mdp_slag_vp(qJ, qJD, qJDD, g, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(5,1),zeros(3,1),zeros(8,1),zeros(22,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PRRRP3_invdynJ_fixb_mdp_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5PRRRP3_invdynJ_fixb_mdp_slag_vp: qJD has to be [5x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [5 1]), ...
  'S5PRRRP3_invdynJ_fixb_mdp_slag_vp: qJDD has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5PRRRP3_invdynJ_fixb_mdp_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5PRRRP3_invdynJ_fixb_mdp_slag_vp: pkin has to be [8x1] (double)');
assert(isreal(MDP) && all(size(MDP) == [22 1]), ...
  'S5PRRRP3_invdynJ_fixb_mdp_slag_vp: MDP has to be [22x1] (double)'); 

%% Symbolic Calculation
% From invdyn_joint_fixb_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2021-01-15 16:23:32
% EndTime: 2021-01-15 16:23:38
% DurationCPUTime: 2.10s
% Computational Cost: add. (1585->267), mult. (3356->331), div. (0->0), fcn. (2250->8), ass. (0->135)
t285 = cos(qJ(3));
t359 = pkin(6) + pkin(7);
t250 = t359 * t285;
t284 = sin(qJ(3));
t334 = qJD(1) * t284;
t220 = qJD(2) * t250 + t334;
t283 = sin(qJ(4));
t214 = t283 * t220;
t324 = qJD(2) * t359;
t219 = t285 * qJD(1) - t284 * t324;
t354 = qJD(3) * pkin(3);
t217 = t219 + t354;
t358 = cos(qJ(4));
t316 = t358 * t217 - t214;
t232 = t283 * t285 + t358 * t284;
t225 = t232 * qJD(2);
t349 = t225 * qJ(5);
t191 = -t349 + t316;
t279 = qJD(3) + qJD(4);
t275 = t285 * pkin(3);
t355 = pkin(2) + t275;
t249 = t359 * t284;
t271 = t285 * qJDD(1);
t202 = qJDD(3) * pkin(3) + t271 - qJDD(2) * t249 + (-t285 * t324 - t334) * qJD(3);
t206 = t219 * qJD(3) + t284 * qJDD(1) + qJDD(2) * t250;
t365 = t358 * t202 - t283 * t206;
t337 = -t283 * t249 + t358 * t250;
t276 = qJDD(3) + qJDD(4);
t320 = t358 * qJD(4);
t357 = pkin(3) * t279;
t364 = -t283 * pkin(3) * t276 - t320 * t357;
t268 = t276 * pkin(4);
t342 = t283 * t284;
t305 = t279 * t342;
t322 = t358 * t285;
t310 = qJD(2) * t322;
t318 = qJDD(2) * t358;
t328 = qJDD(2) * t285;
t312 = t279 * t310 + t283 * t328 + t284 * t318;
t197 = qJD(2) * t305 - t312;
t353 = t197 * qJ(5);
t363 = t268 + t353;
t362 = t358 * qJD(3) + t320;
t278 = pkin(8) + qJ(2);
t270 = cos(t278);
t282 = qJ(3) + qJ(4);
t273 = sin(t282);
t345 = t270 * t273;
t269 = sin(t278);
t347 = t269 * t273;
t274 = cos(t282);
t356 = g(3) * t274;
t361 = g(1) * t345 + g(2) * t347 - t356;
t360 = t225 ^ 2;
t208 = t279 * t232;
t329 = qJDD(2) * t284;
t304 = t283 * t329 - t285 * t318;
t198 = t208 * qJD(2) + t304;
t352 = t198 * qJ(5);
t333 = qJD(2) * t284;
t223 = t283 * t333 - t310;
t351 = t223 * qJ(5);
t350 = t223 * t279;
t346 = t269 * t274;
t344 = t270 * t274;
t341 = qJDD(1) - g(3);
t189 = pkin(4) * t279 + t191;
t340 = t189 - t191;
t207 = -t362 * t285 + t305;
t339 = -t232 * t198 + t207 * t223;
t338 = t358 * t219 - t214;
t336 = pkin(4) * t274 + t275;
t280 = t284 ^ 2;
t335 = -t285 ^ 2 + t280;
t332 = qJD(4) * t283;
t248 = t355 * qJD(2);
t317 = pkin(4) * t223 + qJD(5);
t209 = -t248 + t317;
t331 = qJD(5) + t209;
t330 = qJD(2) * qJD(3);
t327 = pkin(3) * t333;
t326 = t284 * t354;
t323 = qJD(3) * t359;
t216 = t358 * t220;
t319 = t284 * t330;
t315 = -t219 * t283 - t216;
t314 = -t358 * t249 - t250 * t283;
t313 = t279 * t284;
t309 = -g(1) * t347 + g(2) * t345;
t308 = g(1) * t346 - g(2) * t344;
t307 = g(1) * t270 + g(2) * t269;
t306 = g(1) * t269 - g(2) * t270;
t231 = -t322 + t342;
t303 = -t197 * t231 + t208 * t225;
t302 = t207 * t279 - t232 * t276;
t301 = -0.2e1 * pkin(2) * t330 - pkin(6) * qJDD(3);
t300 = -t283 * t217 - t216;
t221 = pkin(3) * t319 - qJDD(2) * t355;
t237 = t284 * t323;
t238 = t285 * t323;
t299 = -t358 * t237 - t283 * t238 - t249 * t320 - t250 * t332;
t222 = t223 ^ 2;
t298 = t225 * t223 * MDP(12) + (-t283 * qJD(2) * t313 + t312 + t350) * MDP(14) - t304 * MDP(15) + (-t222 + t360) * MDP(13) + t276 * MDP(16);
t286 = qJD(3) ^ 2;
t297 = 0.2e1 * qJDD(2) * pkin(2) - pkin(6) * t286 + t306;
t287 = qJD(2) ^ 2;
t296 = pkin(2) * t287 - pkin(6) * qJDD(2) + t307;
t188 = pkin(4) * t198 + qJDD(5) + t221;
t295 = t300 * qJD(4) + t365;
t294 = -t337 * qJD(4) + t283 * t237 - t358 * t238;
t293 = t283 * t202 + t358 * t206 + t217 * t320 - t220 * t332;
t292 = g(1) * t344 + g(2) * t346 + g(3) * t273 - t293;
t291 = t295 + t361;
t290 = -t248 * t223 + t292;
t289 = t248 * t225 + t291;
t288 = t331 * t223 + t292 + t352;
t277 = -qJ(5) - t359;
t266 = t358 * pkin(3) + pkin(4);
t247 = qJDD(3) * t285 - t284 * t286;
t246 = qJDD(3) * t284 + t285 * t286;
t236 = pkin(2) + t336;
t213 = pkin(4) * t231 - t355;
t210 = pkin(4) * t225 + t327;
t205 = pkin(4) * t208 + t326;
t204 = -qJ(5) * t231 + t337;
t203 = -qJ(5) * t232 + t314;
t196 = -t208 * t279 - t231 * t276;
t194 = -t349 + t338;
t193 = t315 + t351;
t192 = -t300 - t351;
t185 = t207 * qJ(5) - t232 * qJD(5) + t294;
t184 = -qJ(5) * t208 - qJD(5) * t231 + t299;
t183 = -t223 * qJD(5) + t293 - t352;
t182 = -t225 * qJD(5) + t295 + t363;
t1 = [t341 * MDP(1) + t247 * MDP(10) - t246 * MDP(11) + (t303 + t339) * MDP(21) + (-t182 * t231 + t183 * t232 - t189 * t208 - t192 * t207 - g(3)) * MDP(22) + (MDP(17) + MDP(19)) * t196 + (MDP(18) + MDP(20)) * t302; qJDD(2) * MDP(2) + t306 * MDP(3) + t307 * MDP(4) + (qJDD(2) * t280 + 0.2e1 * t285 * t319) * MDP(5) + 0.2e1 * (t284 * t328 - t335 * t330) * MDP(6) + t246 * MDP(7) + t247 * MDP(8) + (t301 * t284 + t297 * t285) * MDP(10) + (-t297 * t284 + t301 * t285) * MDP(11) + (-t197 * t232 - t207 * t225) * MDP(12) + (-t303 + t339) * MDP(13) - t302 * MDP(14) + t196 * MDP(15) + (-t198 * t355 - t248 * t208 + t221 * t231 + t223 * t326 + t314 * t276 + t294 * t279 + t308) * MDP(17) + (t197 * t355 + t248 * t207 + t221 * t232 + t225 * t326 - t337 * t276 - t299 * t279 + t309) * MDP(18) + (t185 * t279 + t188 * t231 + t198 * t213 + t203 * t276 + t205 * t223 + t208 * t209 + t308) * MDP(19) + (-t184 * t279 + t188 * t232 - t197 * t213 - t204 * t276 + t205 * t225 - t207 * t209 + t309) * MDP(20) + (-t182 * t232 - t183 * t231 - t184 * t223 - t185 * t225 + t189 * t207 - t192 * t208 + t197 * t203 - t198 * t204 - t307) * MDP(21) + (t183 * t204 + t192 * t184 + t182 * t203 + t189 * t185 + t188 * t213 + t209 * t205 - g(1) * (-t236 * t269 - t270 * t277) - g(2) * (t236 * t270 - t269 * t277)) * MDP(22); MDP(7) * t329 + MDP(8) * t328 + qJDD(3) * MDP(9) + (-g(3) * t285 + t296 * t284 + t271) * MDP(10) + (-t341 * t284 + t296 * t285) * MDP(11) + (-t315 * t279 + (-t223 * t333 + t358 * t276 - t279 * t332) * pkin(3) + t289) * MDP(17) + (-t225 * t327 + t338 * t279 + t290 + t364) * MDP(18) + (-t193 * t279 - t210 * t223 + t266 * t276 - t331 * t225 + (-t216 + (-t217 - t357) * t283) * qJD(4) + t361 + t363 + t365) * MDP(19) + (t194 * t279 - t210 * t225 + t288 + t364) * MDP(20) + (t266 * t197 + (t192 + t193) * t225 + (-t189 + t194) * t223 + (-t198 * t283 + (-t358 * t223 + t225 * t283) * qJD(4)) * pkin(3)) * MDP(21) + (t182 * t266 - t192 * t194 - t189 * t193 - t209 * t210 - g(3) * t336 - t307 * (-pkin(3) * t284 - pkin(4) * t273) + (t183 * t283 + (-t189 * t283 + t358 * t192) * qJD(4)) * pkin(3)) * MDP(22) + t298 + (-t284 * t285 * MDP(5) + t335 * MDP(6)) * t287; (-t300 * t279 + t289) * MDP(17) + (t316 * t279 + t290) * MDP(18) + (t353 + t192 * t279 + 0.2e1 * t268 + (-t209 - t317) * t225 + t291) * MDP(19) + (-t360 * pkin(4) + t191 * t279 + t288) * MDP(20) + (pkin(4) * t197 - t340 * t223) * MDP(21) + (t340 * t192 + (-t209 * t225 + t307 * t273 + t182 - t356) * pkin(4)) * MDP(22) + t298; (t225 * t279 + t304) * MDP(19) + (t312 - t350) * MDP(20) + (-t222 - t360) * MDP(21) + (t189 * t225 + t192 * t223 + t188 - t306) * MDP(22) + (t362 * t284 * MDP(19) + (t279 * MDP(19) * t285 - MDP(20) * t313) * t283) * qJD(2);];
tau = t1;
