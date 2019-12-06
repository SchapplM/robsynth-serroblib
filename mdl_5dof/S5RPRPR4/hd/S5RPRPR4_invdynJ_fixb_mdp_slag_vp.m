% Calculate vector of inverse dynamics joint torques for
% S5RPRPR4
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
%   pkin=[a2,a3,a4,a5,d1,d3,d5,theta2,theta4]';
% MDP [20x1]
%   Minimal dynamic parameter vector (fixed base model)
%   see S5RPRPR4_convert_par2_MPV_fixb.m
% 
% Output:
% tau [5x1]
%   joint torques of inverse dynamics (contains inertial, gravitational Coriolis and centrifugal forces)

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-05 17:54
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tau = S5RPRPR4_invdynJ_fixb_mdp_slag_vp(qJ, qJD, qJDD, g, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(5,1),zeros(3,1),zeros(9,1),zeros(20,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPRPR4_invdynJ_fixb_mdp_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RPRPR4_invdynJ_fixb_mdp_slag_vp: qJD has to be [5x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [5 1]), ...
  'S5RPRPR4_invdynJ_fixb_mdp_slag_vp: qJDD has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RPRPR4_invdynJ_fixb_mdp_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5RPRPR4_invdynJ_fixb_mdp_slag_vp: pkin has to be [9x1] (double)');
assert(isreal(MDP) && all(size(MDP) == [20 1]), ...
  'S5RPRPR4_invdynJ_fixb_mdp_slag_vp: MDP has to be [20x1] (double)'); 

%% Symbolic Calculation
% From invdyn_joint_fixb_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 17:54:04
% EndTime: 2019-12-05 17:54:09
% DurationCPUTime: 1.80s
% Computational Cost: add. (1504->242), mult. (3233->327), div. (0->0), fcn. (2292->14), ass. (0->128)
t284 = sin(pkin(9));
t286 = cos(pkin(9));
t290 = sin(qJ(3));
t293 = cos(qJ(3));
t254 = -t284 * t290 + t286 * t293;
t247 = t254 * qJD(1);
t292 = cos(qJ(5));
t239 = t292 * t247;
t255 = t284 * t293 + t286 * t290;
t249 = t255 * qJD(1);
t289 = sin(qJ(5));
t335 = t249 * t289;
t209 = t239 - t335;
t280 = qJD(3) + qJD(5);
t336 = t209 * t280;
t285 = sin(pkin(8));
t269 = pkin(1) * t285 + pkin(6);
t332 = qJ(4) + t269;
t281 = qJ(1) + pkin(8);
t274 = sin(t281);
t275 = cos(t281);
t316 = -g(2) * t274 + g(3) * t275;
t276 = t293 * qJDD(2);
t260 = t269 * qJDD(1);
t301 = qJ(4) * qJDD(1) + qJD(1) * qJD(4) + qJD(2) * qJD(3) + t260;
t318 = t332 * qJD(1);
t307 = t318 * qJD(3);
t195 = qJDD(3) * pkin(3) - t301 * t290 - t293 * t307 + t276;
t199 = (qJDD(2) - t307) * t290 + t301 * t293;
t178 = t286 * t195 - t199 * t284;
t327 = qJD(1) * qJD(3);
t321 = t293 * t327;
t322 = t290 * t327;
t217 = t255 * qJDD(1) - t284 * t322 + t286 * t321;
t176 = qJDD(3) * pkin(4) - pkin(7) * t217 + t178;
t179 = t284 * t195 + t286 * t199;
t248 = t255 * qJD(3);
t216 = -qJD(1) * t248 + t254 * qJDD(1);
t177 = pkin(7) * t216 + t179;
t232 = t293 * qJD(2) - t318 * t290;
t338 = qJD(3) * pkin(3);
t227 = t232 + t338;
t233 = qJD(2) * t290 + t318 * t293;
t334 = t286 * t233;
t198 = t284 * t227 + t334;
t343 = pkin(7) * t247;
t185 = t198 + t343;
t273 = pkin(3) * t293 + pkin(2);
t287 = cos(pkin(8));
t345 = pkin(1) * t287;
t308 = -t273 - t345;
t245 = t308 * qJD(1) + qJD(4);
t218 = -pkin(4) * t247 + t245;
t277 = qJ(3) + pkin(9) + qJ(5);
t267 = sin(t277);
t268 = cos(t277);
t328 = qJD(5) * t289;
t350 = g(1) * t267 - t289 * t176 - t292 * t177 + t185 * t328 - t218 * t209 + t268 * t316;
t279 = qJDD(3) + qJDD(5);
t310 = t247 * t289 + t292 * t249;
t349 = t279 * MDP(18) + (-t209 ^ 2 + t310 ^ 2) * MDP(15) - t209 * MDP(14) * t310;
t337 = t310 * t280;
t347 = pkin(3) * t322 + t308 * qJDD(1) + qJDD(4);
t346 = -g(1) * t268 + t292 * t176 - t289 * t177 - t218 * t310 + t267 * t316;
t320 = -t292 * t216 + t217 * t289;
t182 = t310 * qJD(5) + t320;
t344 = pkin(3) * t284;
t342 = pkin(7) * t249;
t341 = g(1) * t293;
t223 = t284 * t233;
t197 = t286 * t227 - t223;
t184 = qJD(3) * pkin(4) + t197 - t342;
t333 = t292 * t184;
t331 = qJDD(2) - g(1);
t201 = t286 * t232 - t223;
t319 = qJD(3) * t332;
t236 = qJD(4) * t293 - t290 * t319;
t237 = -qJD(4) * t290 - t293 * t319;
t203 = t286 * t236 + t284 * t237;
t252 = t332 * t290;
t253 = t332 * t293;
t215 = -t284 * t252 + t286 * t253;
t282 = t290 ^ 2;
t330 = -t293 ^ 2 + t282;
t271 = -pkin(2) - t345;
t263 = qJD(1) * t271;
t326 = qJDD(1) * t293;
t324 = t290 * t338;
t323 = qJD(5) * t239 + t289 * t216 + t292 * t217;
t200 = -t232 * t284 - t334;
t202 = -t236 * t284 + t286 * t237;
t214 = -t286 * t252 - t253 * t284;
t317 = g(2) * t275 + g(3) * t274;
t291 = sin(qJ(1));
t294 = cos(qJ(1));
t315 = g(2) * t294 + g(3) * t291;
t314 = -t289 * t184 - t292 * t185;
t251 = t254 * qJD(3);
t309 = t292 * t254 - t255 * t289;
t186 = t309 * qJD(5) - t248 * t289 + t251 * t292;
t220 = t254 * t289 + t255 * t292;
t313 = t186 * t280 + t220 * t279;
t204 = -pkin(7) * t255 + t214;
t205 = pkin(7) * t254 + t215;
t312 = t204 * t292 - t205 * t289;
t311 = t204 * t289 + t205 * t292;
t270 = pkin(3) * t286 + pkin(4);
t306 = t270 * t289 + t292 * t344;
t305 = t270 * t292 - t289 * t344;
t181 = -t249 * t328 + t323;
t303 = -qJD(1) * t263 - t260 + t316;
t302 = 0.2e1 * t263 * qJD(3) - qJDD(3) * t269;
t295 = qJD(3) ^ 2;
t299 = -0.2e1 * qJDD(1) * t271 - t269 * t295 + t317;
t288 = -qJ(4) - pkin(6);
t259 = qJDD(3) * t293 - t290 * t295;
t258 = qJDD(3) * t290 + t293 * t295;
t235 = pkin(4) * t248 + t324;
t234 = pkin(3) * qJD(1) * t290 + pkin(4) * t249;
t231 = -pkin(4) * t254 + t308;
t196 = -pkin(4) * t216 + t347;
t191 = -pkin(7) * t248 + t203;
t190 = -pkin(7) * t251 + t202;
t189 = t201 - t342;
t188 = t200 - t343;
t187 = t220 * qJD(5) + t292 * t248 + t251 * t289;
t180 = -t187 * t280 + t279 * t309;
t1 = [qJDD(1) * MDP(1) + t315 * MDP(2) + (-g(2) * t291 + g(3) * t294) * MDP(3) + (t315 + (t285 ^ 2 + t287 ^ 2) * qJDD(1) * pkin(1)) * pkin(1) * MDP(4) + (qJDD(1) * t282 + 0.2e1 * t290 * t321) * MDP(5) + 0.2e1 * (t290 * t326 - t330 * t327) * MDP(6) + t258 * MDP(7) + t259 * MDP(8) + (t302 * t290 + t299 * t293) * MDP(10) + (-t299 * t290 + t302 * t293) * MDP(11) + (-t178 * t255 + t179 * t254 - t197 * t251 - t198 * t248 - t202 * t249 + t203 * t247 - t214 * t217 + t215 * t216 - t316) * MDP(12) + (t179 * t215 + t198 * t203 + t178 * t214 + t197 * t202 + t245 * t324 - g(2) * (-pkin(1) * t294 - t273 * t275 + t274 * t288) - g(3) * (-pkin(1) * t291 - t273 * t274 - t275 * t288) + t347 * t308) * MDP(13) + (t181 * t220 + t186 * t310) * MDP(14) + (t181 * t309 - t182 * t220 + t186 * t209 - t187 * t310) * MDP(15) + t313 * MDP(16) + t180 * MDP(17) + (-t235 * t209 + t231 * t182 - t196 * t309 + t218 * t187 + (-t311 * qJD(5) + t190 * t292 - t191 * t289) * t280 + t312 * t279 + t317 * t268) * MDP(19) + (t235 * t310 + t231 * t181 + t196 * t220 + t218 * t186 - (t312 * qJD(5) + t190 * t289 + t191 * t292) * t280 - t311 * t279 - t317 * t267) * MDP(20); t331 * MDP(4) + t259 * MDP(10) - t258 * MDP(11) + (t216 * t255 - t217 * t254 + t247 * t251 + t248 * t249) * MDP(12) + (t178 * t254 + t179 * t255 - t197 * t248 + t198 * t251 - g(1)) * MDP(13) + t180 * MDP(19) - t313 * MDP(20); t290 * qJDD(1) * MDP(7) + MDP(8) * t326 + qJDD(3) * MDP(9) + (t303 * t290 + t276 - t341) * MDP(10) + (-t331 * t290 + t303 * t293) * MDP(11) + ((t198 + t200) * t249 + (t197 - t201) * t247 + (t216 * t284 - t217 * t286) * pkin(3)) * MDP(12) + (-t197 * t200 - t198 * t201 + (-t341 + t178 * t286 + t179 * t284 + (-qJD(1) * t245 + t316) * t290) * pkin(3)) * MDP(13) + (t181 - t336) * MDP(16) + (-t182 + t337) * MDP(17) + (t305 * t279 + t234 * t209 - (t188 * t292 - t189 * t289) * t280 + (-t306 * t280 + t314) * qJD(5) + t346) * MDP(19) + (-t306 * t279 - t234 * t310 + (t188 * t289 + t189 * t292) * t280 + (-t305 * t280 - t333) * qJD(5) + t350) * MDP(20) + (-t290 * t293 * MDP(5) + t330 * MDP(6)) * qJD(1) ^ 2 + t349; (-t247 ^ 2 - t249 ^ 2) * MDP(12) + (t182 + t337) * MDP(19) + (t181 + t336) * MDP(20) + (t197 * t249 - t198 * t247 - t317 + t347) * MDP(13); (t323 - t336) * MDP(16) + (-t320 + t337) * MDP(17) + (-t314 * t280 + t346) * MDP(19) + ((-t185 * t289 + t333) * t280 + t350) * MDP(20) + (-MDP(16) * t335 - t310 * MDP(17) + t314 * MDP(19) - MDP(20) * t333) * qJD(5) + t349;];
tau = t1;
