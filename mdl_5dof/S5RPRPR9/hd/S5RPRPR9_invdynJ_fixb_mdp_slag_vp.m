% Calculate vector of inverse dynamics joint torques for
% S5RPRPR9
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
%   pkin=[a2,a3,a4,a5,d1,d3,d5,theta2]';
% MDP [22x1]
%   Minimal dynamic parameter vector (fixed base model)
%   see S5RPRPR9_convert_par2_MPV_fixb.m
% 
% Output:
% tau [5x1]
%   joint torques of inverse dynamics (contains inertial, gravitational Coriolis and centrifugal forces)

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 18:25
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tau = S5RPRPR9_invdynJ_fixb_mdp_slag_vp(qJ, qJD, qJDD, g, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(5,1),zeros(3,1),zeros(8,1),zeros(22,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPRPR9_invdynJ_fixb_mdp_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RPRPR9_invdynJ_fixb_mdp_slag_vp: qJD has to be [5x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [5 1]), ...
  'S5RPRPR9_invdynJ_fixb_mdp_slag_vp: qJDD has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RPRPR9_invdynJ_fixb_mdp_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RPRPR9_invdynJ_fixb_mdp_slag_vp: pkin has to be [8x1] (double)');
assert(isreal(MDP) && all(size(MDP) == [22 1]), ...
  'S5RPRPR9_invdynJ_fixb_mdp_slag_vp: MDP has to be [22x1] (double)'); 

%% Symbolic Calculation
% From invdyn_joint_fixb_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 18:24:43
% EndTime: 2019-12-31 18:24:47
% DurationCPUTime: 2.63s
% Computational Cost: add. (1059->288), mult. (2045->373), div. (0->0), fcn. (1189->10), ass. (0->146)
t263 = sin(pkin(8));
t246 = pkin(1) * t263 + pkin(6);
t236 = t246 * qJDD(1);
t364 = qJD(2) * qJD(3) + t236;
t238 = t246 * qJD(1);
t266 = sin(qJ(3));
t269 = cos(qJ(3));
t334 = t269 * qJD(2) - t266 * t238;
t361 = -qJD(4) + t334;
t363 = qJDD(1) * MDP(12);
t201 = -qJD(3) * pkin(3) - t361;
t330 = qJD(1) * t269;
t320 = qJD(1) * qJD(3);
t308 = t269 * t320;
t317 = qJDD(1) * t266;
t288 = t308 + t317;
t225 = qJDD(5) + t288;
t268 = cos(qJ(5));
t217 = t268 * t225;
t265 = sin(qJ(5));
t230 = qJD(3) * t265 + t268 * t330;
t328 = qJD(3) * t230;
t362 = (t328 - t217) * MDP(21);
t211 = t266 * qJD(2) + t269 * t238;
t203 = -qJD(3) * qJ(4) - t211;
t260 = qJ(1) + pkin(8);
t253 = sin(t260);
t254 = cos(t260);
t359 = g(1) * t254 + g(2) * t253;
t331 = qJD(1) * t266;
t321 = pkin(4) * t331 - t361;
t326 = qJD(3) * t269;
t358 = t238 * t326 + t266 * t364;
t250 = pkin(4) * t330;
t196 = -t203 + t250;
t245 = qJD(5) + t331;
t271 = -pkin(3) - pkin(7);
t357 = t271 * t225 + (t196 - t250 - t211) * t245;
t264 = cos(pkin(8));
t356 = pkin(1) * t264;
t353 = g(3) * t266;
t259 = g(3) * t269;
t352 = pkin(4) + t246;
t351 = qJ(4) * t266;
t350 = qJ(4) * t269;
t349 = qJDD(3) * pkin(3);
t327 = qJD(3) * t266;
t303 = t266 * qJDD(2) - t238 * t327 + t269 * t364;
t314 = qJDD(3) * qJ(4);
t280 = qJD(3) * qJD(4) + t303 + t314;
t309 = t266 * t320;
t316 = qJDD(1) * t269;
t187 = (-t309 + t316) * pkin(4) + t280;
t348 = t187 * t269;
t336 = t268 * qJDD(3) + t265 * t309;
t193 = -t230 * qJD(5) - t265 * t316 + t336;
t347 = t193 * t268;
t292 = qJD(3) * qJD(5) + t316;
t322 = qJD(5) * t269;
t307 = qJD(1) * t322;
t335 = -t265 * t307 - t268 * t309;
t194 = qJDD(3) * t265 + t292 * t268 + t335;
t346 = t194 * t266;
t345 = t225 * t265;
t344 = t230 * t245;
t232 = qJD(3) * t268 - t265 * t330;
t343 = t232 * t245;
t342 = t265 * t266;
t341 = t265 * t269;
t340 = t266 * t268;
t273 = qJD(1) ^ 2;
t339 = t269 * t273;
t338 = t193 * t266 + t232 * t326;
t323 = qJD(5) * t245;
t310 = t265 * t323;
t312 = t245 * t340;
t337 = qJD(3) * t312 + t269 * t310;
t261 = t266 ^ 2;
t262 = t269 ^ 2;
t333 = t261 - t262;
t300 = pkin(3) * t269 + t351;
t293 = pkin(2) + t300;
t218 = -t293 - t356;
t204 = qJD(1) * t218;
t252 = pkin(3) * t331;
t233 = -qJ(4) * t330 + t252;
t332 = qJD(1) * t233;
t247 = -pkin(2) - t356;
t239 = qJD(1) * t247;
t325 = qJD(4) * t266;
t209 = t271 * t269 + t247 - t351;
t197 = t209 * qJD(1);
t324 = qJD(5) * t197;
t318 = qJDD(1) * t218;
t315 = qJDD(2) * t269;
t313 = qJDD(3) * t246;
t311 = t265 * t327;
t220 = t352 * t269;
t285 = qJDD(4) - t315 + t358;
t190 = t285 - t349;
t305 = -qJD(3) * t203 - t190;
t244 = pkin(3) * t309;
t298 = pkin(7) * t266 - t350;
t279 = t298 * qJD(3) - t325;
t188 = t279 * qJD(1) + t209 * qJDD(1) + t244;
t195 = t271 * qJD(3) + t321;
t304 = qJD(5) * t195 + t188;
t267 = sin(qJ(1));
t270 = cos(qJ(1));
t301 = g(1) * t267 - g(2) * t270;
t299 = pkin(3) * t266 - t350;
t297 = -t324 + t259;
t184 = t195 * t265 + t197 * t268;
t295 = t245 * t265;
t291 = t301 * pkin(1);
t290 = t268 * t323 + t345;
t289 = -qJ(4) * t326 - t325;
t287 = qJD(3) * t211 - t259 - t358;
t272 = qJD(3) ^ 2;
t286 = g(1) * t253 - g(2) * t254 - t246 * t272;
t284 = -qJD(1) * t239 + t359;
t283 = -t268 * t322 + t311;
t282 = 0.2e1 * t239 * qJD(3) - t313;
t281 = -0.2e1 * t204 * qJD(3) + t313;
t278 = -0.2e1 * qJDD(1) * t247 + t286;
t277 = t280 * t269 + t190 * t266 + (t201 * t269 + t203 * t266) * qJD(3);
t191 = t289 * qJD(1) + t244 + t318;
t251 = pkin(3) * t327;
t216 = t251 + t289;
t276 = -qJD(1) * t216 - t191 + t286 - t318;
t275 = -t353 + t187 - t359 * t269 + (t298 * qJD(1) - qJD(5) * t271 + t252) * t245;
t235 = qJDD(3) * t269 - t266 * t272;
t234 = qJDD(3) * t266 + t269 * t272;
t219 = t352 * t266;
t215 = qJD(3) * t220;
t214 = t352 * t327;
t208 = -t253 * t342 + t254 * t268;
t207 = t253 * t340 + t254 * t265;
t206 = t253 * t268 + t254 * t342;
t205 = -t253 * t265 + t254 * t340;
t202 = t251 + t279;
t198 = t204 * t331;
t186 = t288 * pkin(4) + t271 * qJDD(3) + t285;
t185 = t268 * t186;
t183 = t195 * t268 - t197 * t265;
t1 = [qJDD(1) * MDP(1) + t301 * MDP(2) + (g(1) * t270 + g(2) * t267) * MDP(3) + ((t263 ^ 2 + t264 ^ 2) * pkin(1) ^ 2 * qJDD(1) + t291) * MDP(4) + (qJDD(1) * t261 + 0.2e1 * t266 * t308) * MDP(5) + 0.2e1 * (t266 * t316 - t333 * t320) * MDP(6) + t234 * MDP(7) + t235 * MDP(8) + (t282 * t266 + t278 * t269) * MDP(10) + (-t278 * t266 + t282 * t269) * MDP(11) + ((t261 + t262) * t236 + t277 - t359) * MDP(12) + (t281 * t266 - t276 * t269) * MDP(13) + (t276 * t266 + t281 * t269) * MDP(14) + (t191 * t218 + t204 * t216 + t291 + (-g(1) * pkin(6) - g(2) * t293) * t254 + (-g(2) * pkin(6) + g(1) * t293) * t253 + t277 * t246) * MDP(15) + (-t193 * t341 + t283 * t232) * MDP(16) + ((-t230 * t265 + t232 * t268) * t327 + (-t347 + t194 * t265 + (t230 * t268 + t232 * t265) * qJD(5)) * t269) * MDP(17) + (-t225 * t341 + t283 * t245 + t338) * MDP(18) + (-t346 + (-t328 - t217) * t269 + t337) * MDP(19) + (t225 * t266 + t245 * t326) * MDP(20) + ((-t202 * t265 + t215 * t268) * t245 + (-t209 * t265 + t219 * t268) * t225 + (-t188 * t265 + t185) * t266 - t214 * t230 + t220 * t194 + t268 * t348 - g(1) * t208 - g(2) * t206 + (t183 * t269 - t196 * t340) * qJD(3) + ((-t209 * t268 - t219 * t265) * t245 - t184 * t266 - t196 * t341) * qJD(5)) * MDP(21) + (-t184 * t326 + g(1) * t207 - g(2) * t205 + t220 * t193 - t214 * t232 + (-(qJD(5) * t219 + t202) * t245 - t209 * t225 - t304 * t266 - t196 * t322) * t268 + (-(-qJD(5) * t209 + t215) * t245 - t219 * t225 - t348 + (qJD(3) * t196 - t186 + t324) * t266) * t265) * MDP(22); (qJDD(2) - g(3)) * MDP(4) + (t201 * t327 + t266 * t280 - g(3)) * MDP(15) + (t337 + t346) * MDP(21) + (-t245 * t311 + t338) * MDP(22) + (MDP(10) - MDP(13)) * t235 + (-MDP(11) + MDP(14)) * t234 + (t305 * MDP(15) + t290 * MDP(22) + t362) * t269; -t266 * MDP(5) * t339 + t333 * MDP(6) * t273 + MDP(7) * t317 + MDP(8) * t316 + qJDD(3) * MDP(9) + (t284 * t266 + t287 + t315) * MDP(10) + (qJD(3) * t334 + t284 * t269 - t303 + t353) * MDP(11) - t299 * t363 + (-0.2e1 * t349 + qJDD(4) + t198 + (-qJDD(2) - t332) * t269 - t359 * t266 - t287) * MDP(13) + (0.2e1 * t314 + (-g(3) + t332) * t266 + (0.2e1 * qJD(4) - t334) * qJD(3) + (qJD(1) * t204 - t359) * t269 + t303) * MDP(14) + (-t190 * pkin(3) - g(3) * t300 + t280 * qJ(4) - t201 * t211 + t203 * t361 - t204 * t233 + t359 * t299) * MDP(15) + (-t232 * t295 + t347) * MDP(16) + ((-t194 - t343) * t268 + (-t193 + t344) * t265) * MDP(17) + (-t310 + t217 + (-t232 * t269 - t245 * t342) * qJD(1)) * MDP(18) + ((t230 * t269 - t312) * qJD(1) - t290) * MDP(19) - t245 * MDP(20) * t330 + (qJ(4) * t194 - t183 * t330 + t321 * t230 + t275 * t265 + t268 * t357) * MDP(21) + (qJ(4) * t193 + t184 * t330 + t321 * t232 - t265 * t357 + t275 * t268) * MDP(22); qJDD(3) * MDP(13) + (-t261 * t273 - t272) * MDP(14) + (t198 + t259 - t305) * MDP(15) - t362 + (-qJD(3) * t232 - t345) * MDP(22) + (MDP(13) * t339 - MDP(15) * t359 + t363) * t266 + (-t245 * MDP(22) * t268 - MDP(21) * t295) * t245; t232 * t230 * MDP(16) + (-t230 ^ 2 + t232 ^ 2) * MDP(17) + (t336 + t344) * MDP(18) + (-t335 + t343) * MDP(19) + t225 * MDP(20) + (-g(1) * t205 - g(2) * t207 + t184 * t245 - t196 * t232 + t185) * MDP(21) + (g(1) * t206 - g(2) * t208 + t183 * t245 + t196 * t230) * MDP(22) + (-MDP(18) * t307 - t292 * MDP(19) + t297 * MDP(21) - t304 * MDP(22)) * t268 + (-t292 * MDP(18) - qJDD(3) * MDP(19) - t304 * MDP(21) + (-t186 - t297) * MDP(22)) * t265;];
tau = t1;
