% Calculate vector of inverse dynamics joint torques for
% S5RRPRP1
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
% MDP [18x1]
%   Minimal dynamic parameter vector (fixed base model)
%   see S5RRPRP1_convert_par2_MPV_fixb.m
% 
% Output:
% tau [5x1]
%   joint torques of inverse dynamics (contains inertial, gravitational Coriolis and centrifugal forces)

% Quelle: HybrDyn-Toolbox
% Datum: 2021-01-15 20:09
% Revision: 24b2e7d74a0c1a3b64fa2f8f5ad758691ad61af3 (2021-01-15)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tau = S5RRPRP1_invdynJ_fixb_mdp_slag_vp(qJ, qJD, qJDD, g, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(5,1),zeros(3,1),zeros(8,1),zeros(18,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRPRP1_invdynJ_fixb_mdp_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RRPRP1_invdynJ_fixb_mdp_slag_vp: qJD has to be [5x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [5 1]), ...
  'S5RRPRP1_invdynJ_fixb_mdp_slag_vp: qJDD has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RRPRP1_invdynJ_fixb_mdp_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RRPRP1_invdynJ_fixb_mdp_slag_vp: pkin has to be [8x1] (double)');
assert(isreal(MDP) && all(size(MDP) == [18 1]), ...
  'S5RRPRP1_invdynJ_fixb_mdp_slag_vp: MDP has to be [18x1] (double)'); 

%% Symbolic Calculation
% From invdyn_joint_fixb_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2021-01-15 20:08:57
% EndTime: 2021-01-15 20:09:01
% DurationCPUTime: 1.50s
% Computational Cost: add. (1428->270), mult. (2253->323), div. (0->0), fcn. (1262->12), ass. (0->154)
t366 = qJDD(3) - g(1);
t286 = cos(qJ(2));
t355 = pkin(1) * qJD(2);
t323 = qJD(1) * t355;
t283 = sin(qJ(2));
t328 = qJDD(1) * t283;
t365 = pkin(1) * t328 + t286 * t323;
t285 = cos(qJ(4));
t275 = qJD(1) + qJD(2);
t356 = pkin(1) * qJD(1);
t325 = t286 * t356;
t232 = pkin(2) * t275 + t325;
t280 = cos(pkin(8));
t326 = t283 * t356;
t248 = t280 * t326;
t279 = sin(pkin(8));
t209 = t279 * t232 + t248;
t204 = pkin(7) * t275 + t209;
t318 = qJ(5) * t275 + t204;
t302 = t318 * t285;
t278 = qJ(1) + qJ(2);
t264 = pkin(8) + t278;
t253 = cos(t264);
t249 = g(3) * t253;
t252 = sin(t264);
t358 = g(2) * t252;
t364 = t249 - t358;
t363 = pkin(1) * t286;
t362 = pkin(2) * t280;
t282 = sin(qJ(4));
t276 = t282 ^ 2;
t361 = pkin(4) * t276;
t360 = pkin(4) * t285;
t359 = g(1) * t285;
t357 = g(2) * t253;
t354 = qJD(4) * pkin(4);
t353 = qJDD(4) * pkin(4);
t260 = pkin(2) + t363;
t344 = t279 * t283;
t316 = -pkin(1) * t344 + t260 * t280;
t219 = -pkin(3) - t316;
t214 = t219 - t360;
t352 = t214 * t275;
t221 = t279 * t325 + t248;
t351 = t221 * t275;
t259 = pkin(3) + t360;
t233 = -t259 - t362;
t274 = qJDD(1) + qJDD(2);
t350 = t233 * t274;
t349 = t252 * t282;
t348 = t274 * t282;
t347 = t274 * t285;
t346 = t275 * t282;
t345 = t275 * t285;
t343 = t280 * t283;
t342 = t282 * t285;
t336 = pkin(1) * t343 + t279 * t260;
t220 = pkin(7) + t336;
t341 = -qJ(5) - t220;
t254 = pkin(2) * t279 + pkin(7);
t340 = -qJ(5) - t254;
t266 = t285 * qJD(3);
t192 = -t282 * t318 + t266;
t189 = t192 + t354;
t339 = t189 - t192;
t247 = t279 * t326;
t223 = t280 * t325 - t247;
t332 = qJD(4) * t282;
t338 = t221 * t345 + t223 * t332;
t337 = g(3) * t349 + t282 * t357;
t277 = t285 ^ 2;
t335 = -t276 - t277;
t334 = t276 - t277;
t333 = qJD(4) * t275;
t331 = qJD(4) * t285;
t329 = qJD(4) * qJD(3);
t327 = qJDD(4) * t220;
t324 = pkin(4) * t332;
t267 = sin(t278);
t257 = pkin(2) * t267;
t281 = -qJ(5) - pkin(7);
t321 = t252 * t259 + t253 * t281 + t257;
t261 = qJDD(1) * t363;
t218 = pkin(2) * t274 - t283 * t323 + t261;
t197 = t279 * t218 + t280 * t365;
t320 = t275 * t332;
t268 = cos(t278);
t319 = g(2) * t267 - g(3) * t268;
t224 = (t280 * t286 - t344) * t355;
t317 = t219 * t275 - t224;
t208 = t232 * t280 - t247;
t315 = qJD(4) * t341;
t314 = qJD(4) * t340;
t313 = 0.2e1 * t275 * t331;
t312 = qJD(1) * (-qJD(2) + t275);
t244 = pkin(4) * t320;
t196 = t218 * t280 - t279 * t365;
t185 = -t259 * t274 + qJDD(5) - t196 + t244;
t198 = -t259 * t275 + qJD(5) - t208;
t310 = t185 * t282 + t198 * t331 + t337;
t190 = -pkin(3) * t274 - t196;
t203 = -pkin(3) * t275 - t208;
t309 = t190 * t282 + t203 * t331 + t337;
t263 = t285 * qJDD(3);
t308 = g(2) * t349 + t263 - t359;
t201 = t204 * t332;
t191 = pkin(7) * t274 + t197;
t293 = -qJ(5) * t274 - t191 - t329;
t290 = qJD(5) * t275 - t293;
t183 = -t201 + (-qJ(5) * t333 + qJDD(3)) * t282 + t290 * t285;
t307 = t183 * t285 + t364;
t306 = g(3) * t252 + t357;
t305 = -g(2) * t268 - g(3) * t267;
t258 = pkin(2) * t268;
t304 = -t252 * t281 + t253 * t259 + t258;
t288 = qJD(4) ^ 2;
t234 = qJDD(4) * t282 + t285 * t288;
t235 = qJDD(4) * t285 - t282 * t288;
t303 = 0.2e1 * (t274 * t342 - t333 * t334) * MDP(9) + (t274 * t276 + t282 * t313) * MDP(8) + t234 * MDP(10) + t235 * MDP(11) + t274 * MDP(4);
t193 = qJD(3) * t282 + t302;
t301 = t189 * t282 - t193 * t285;
t222 = (t279 * t286 + t343) * t355;
t210 = t222 + t324;
t300 = t210 * t275 + t214 * t274;
t255 = -pkin(3) - t362;
t299 = t254 * t288 + t255 * t274;
t298 = -t185 - t306;
t297 = -t190 - t306;
t296 = t261 + t305;
t295 = -t203 * t275 - t191 - t249;
t294 = -t366 * t282 + t285 * t358 + t201;
t292 = -qJDD(4) * t254 + t255 * t333;
t291 = t219 * t274 + t220 * t288 + t222 * t275;
t289 = -t249 + (-qJD(5) - t198) * t275 + t293;
t287 = cos(qJ(1));
t284 = sin(qJ(1));
t273 = t275 ^ 2;
t271 = t287 * pkin(1);
t270 = t284 * pkin(1);
t269 = t285 * qJ(5);
t265 = t285 * qJD(5);
t227 = t254 * t285 + t269;
t226 = t340 * t282;
t216 = -qJD(5) * t282 + t285 * t314;
t215 = t282 * t314 + t265;
t212 = t223 * t331;
t206 = t220 * t285 + t269;
t205 = t341 * t282;
t199 = t203 * t332;
t194 = t198 * t332;
t187 = (-qJD(5) - t224) * t282 + t285 * t315;
t186 = t224 * t285 + t282 * t315 + t265;
t182 = -qJD(4) * t302 - t282 * t290 + t263 + t353;
t1 = [(t197 * t336 + t209 * t224 + t196 * t316 - t208 * t222 - g(2) * (t258 + t271) - g(3) * (t257 + t270)) * MDP(7) + (t183 * t206 + t193 * t186 + t182 * t205 + t189 * t187 + t185 * t214 + t198 * t210 - g(2) * (t271 + t304) - g(3) * (t270 + t321)) * MDP(18) + t307 * MDP(17) + t199 * MDP(13) + (-qJD(4) * t186 - qJDD(4) * t206 + t310) * MDP(16) + t309 * MDP(14) + (qJD(4) * t187 + qJDD(4) * t205 + t194) * MDP(15) + t303 + t319 * MDP(6) + (t274 * t286 * MDP(5) + (-qJDD(1) - t274) * MDP(6) * t283 + (MDP(5) * t283 + MDP(6) * t286) * qJD(2) * (-qJD(1) - t275)) * pkin(1) + qJDD(1) * MDP(1) + t296 * MDP(5) + ((-t291 + t297) * MDP(13) - MDP(14) * t327 + (t298 - t300) * MDP(15) + (t186 * t275 + t206 * t274) * MDP(17) + (t317 * MDP(14) + MDP(16) * t352 + (-t205 * t275 - t189) * MDP(17)) * qJD(4)) * t285 + (-MDP(13) * t327 + t291 * MDP(14) + t300 * MDP(16) + (-t187 * t275 - t205 * t274 - t182) * MDP(17) + (t317 * MDP(13) + MDP(15) * t352 + (-t206 * t275 - t193) * MDP(17)) * qJD(4)) * t282 + (-g(2) * t287 - g(3) * t284) * MDP(2) + (g(2) * t284 - g(3) * t287) * MDP(3); (pkin(1) * t283 * t312 + t296) * MDP(5) + ((t286 * t312 - t328) * pkin(1) + t319) * MDP(6) + (t208 * t221 - t209 * t223 + (t196 * t280 + t197 * t279 + t305) * pkin(2)) * MDP(7) + (t199 + t292 * t282 + (t297 - t299) * t285 + t338) * MDP(13) + (t212 + t292 * t285 + (t299 - t351) * t282 + t309) * MDP(14) + (qJDD(4) * t226 + t194 + (t233 * t346 + t216) * qJD(4) + (-t244 + t298 - t350) * t285 + t338) * MDP(15) + (-qJDD(4) * t227 + t212 + (t350 - t351) * t282 + (-t215 + (t233 * t285 + t361) * t275) * qJD(4) + t310) * MDP(16) + ((-qJD(4) * t189 + t227 * t274) * t285 + (-qJD(4) * t193 - t226 * t274 - t182) * t282 + (t215 * t285 - t216 * t282 + t335 * t223 + (-t226 * t285 - t227 * t282) * qJD(4)) * t275 + t307) * MDP(17) + (t183 * t227 + t182 * t226 + t185 * t233 - g(2) * t304 - g(3) * t321 + (-t221 + t324) * t198 + (-t223 * t285 + t215) * t193 + (t223 * t282 + t216) * t189) * MDP(18) + t303; t366 * MDP(7) + (-qJD(4) * t301 + t182 * t285 + t183 * t282 - g(1)) * MDP(18) + (MDP(13) + MDP(15)) * t235 + (-MDP(14) - MDP(16)) * t234; -t273 * MDP(8) * t342 + t334 * t273 * MDP(9) + MDP(10) * t348 + MDP(11) * t347 + qJDD(4) * MDP(12) + (t282 * t295 + t308) * MDP(13) + ((-t204 * t282 + t266) * qJD(4) + (t295 - t329) * t285 + t294) * MDP(14) + (0.2e1 * t353 + (t193 - t302) * qJD(4) + (t273 * t360 + t289) * t282 + t308) * MDP(15) + (-t273 * t361 + (qJ(5) * t346 + t192) * qJD(4) + t289 * t285 + t294) * MDP(16) + (-pkin(4) * t348 + (t339 - t354) * t345) * MDP(17) + (t339 * t193 + (-t359 + t182 + (-t198 * t275 - t364) * t282) * pkin(4)) * MDP(18); (0.2e1 * t320 - t347) * MDP(15) + (t313 + t348) * MDP(16) + (t275 * t301 - t298) * MDP(18) + t335 * MDP(17) * t273;];
tau = t1;
