% Calculate vector of inverse dynamics joint torques for
% S5RRPRR7
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
%   pkin=[a2,a3,a4,a5,d1,d2,d4,d5]';
% MDP [23x1]
%   Minimal dynamic parameter vector (fixed base model)
%   see S5RRPRR7_convert_par2_MPV_fixb.m
% 
% Output:
% tau [5x1]
%   joint torques of inverse dynamics (contains inertial, gravitational Coriolis and centrifugal forces)

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 20:16
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tau = S5RRPRR7_invdynJ_fixb_mdp_slag_vp(qJ, qJD, qJDD, g, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(5,1),zeros(3,1),zeros(8,1),zeros(23,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRPRR7_invdynJ_fixb_mdp_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RRPRR7_invdynJ_fixb_mdp_slag_vp: qJD has to be [5x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [5 1]), ...
  'S5RRPRR7_invdynJ_fixb_mdp_slag_vp: qJDD has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RRPRR7_invdynJ_fixb_mdp_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RRPRR7_invdynJ_fixb_mdp_slag_vp: pkin has to be [8x1] (double)');
assert(isreal(MDP) && all(size(MDP) == [23 1]), ...
  'S5RRPRR7_invdynJ_fixb_mdp_slag_vp: MDP has to be [23x1] (double)'); 

%% Symbolic Calculation
% From invdyn_joint_fixb_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 20:15:46
% EndTime: 2019-12-31 20:15:48
% DurationCPUTime: 1.69s
% Computational Cost: add. (1292->254), mult. (1706->317), div. (0->0), fcn. (1030->12), ass. (0->130)
t280 = qJ(1) + qJ(2);
t268 = sin(t280);
t270 = cos(t280);
t358 = -g(1) * t268 + g(2) * t270;
t281 = sin(qJ(5));
t282 = sin(qJ(4));
t329 = qJD(5) * t281;
t331 = qJD(4) * t282;
t359 = -t281 * t331 - t282 * t329;
t285 = cos(qJ(5));
t286 = cos(qJ(4));
t227 = -t281 * t282 + t285 * t286;
t271 = t282 * pkin(4);
t251 = qJ(3) + t271;
t311 = -g(1) * t270 - g(2) * t268;
t275 = qJD(4) + qJD(5);
t283 = sin(qJ(2));
t354 = pkin(1) * qJD(1);
t325 = t283 * t354;
t276 = qJD(1) + qJD(2);
t351 = qJ(3) * t276;
t230 = t325 + t351;
t300 = -t230 * t276 + t358;
t272 = t276 ^ 2;
t289 = -pkin(2) - pkin(7);
t274 = qJDD(1) + qJDD(2);
t357 = pkin(2) * t274;
t287 = cos(qJ(2));
t260 = -pkin(1) * t287 - pkin(2);
t247 = -pkin(7) + t260;
t356 = -pkin(8) + t247;
t355 = -pkin(8) + t289;
t353 = pkin(1) * qJDD(1);
t352 = qJ(3) * t274;
t324 = t287 * t354;
t312 = qJD(3) - t324;
t217 = t289 * t276 + t312;
t201 = (-pkin(8) * t276 + t217) * t282;
t350 = t201 * t285;
t348 = t274 * t282;
t347 = t274 * t286;
t346 = t275 * t283;
t345 = t276 * t286;
t343 = t281 * t286;
t342 = t282 * t286;
t226 = t282 * t285 + t343;
t199 = t275 * t226;
t273 = qJDD(4) + qJDD(5);
t340 = -t199 * t275 + t227 * t273;
t322 = qJD(2) * t354;
t338 = t283 * t353 + t287 * t322;
t209 = qJD(3) * t276 + t338 + t352;
t330 = qJD(4) * t286;
t339 = t209 * t282 + t230 * t330;
t337 = t270 * pkin(2) + t268 * qJ(3);
t336 = -t283 * t322 + t287 * t353;
t290 = qJD(4) ^ 2;
t335 = -t272 - t290;
t278 = t286 ^ 2;
t334 = t282 ^ 2 - t278;
t333 = qJD(2) * t283;
t332 = qJD(2) * t287;
t328 = qJDD(4) * t247;
t327 = qJDD(4) * t282;
t326 = qJDD(4) * t289;
t323 = pkin(1) * t333;
t321 = t276 * t333;
t320 = t276 * t330;
t221 = t356 * t286;
t234 = t355 * t286;
t317 = qJDD(3) - t336;
t250 = pkin(1) * t283 + qJ(3);
t316 = -pkin(2) * t268 + t270 * qJ(3);
t315 = t275 * t286;
t314 = -t336 + t358;
t313 = t338 + t311;
t202 = -pkin(8) * t345 + t286 * t217;
t240 = pkin(1) * t332 + qJD(3);
t197 = qJD(4) * pkin(4) + t202;
t310 = -t197 * t281 - t350;
t200 = t285 * t315 + t359;
t309 = -t200 * t275 - t226 * t273;
t220 = t356 * t282;
t308 = t220 * t285 + t221 * t281;
t307 = -t220 * t281 + t221 * t285;
t233 = t355 * t282;
t306 = t233 * t285 + t234 * t281;
t305 = -t233 * t281 + t234 * t285;
t304 = t274 * t343 + t359 * t276;
t212 = t317 - t357;
t303 = -t325 + t351;
t302 = t250 * t276 + t323;
t301 = t320 + t348;
t185 = -t199 * t276 + t227 * t274;
t213 = t227 * t276;
t214 = t226 * t276;
t299 = t213 * t214 * MDP(17) + (t214 * t275 + t185) * MDP(19) + (t213 * t275 + (-t275 * t345 - t348) * t285 - t304) * MDP(20) + (t213 ^ 2 - t214 ^ 2) * MDP(18) + t273 * MDP(21);
t298 = t276 * t325 - t314;
t194 = t301 * pkin(4) + t209;
t216 = t251 * t276 + t325;
t279 = qJ(4) + qJ(5);
t267 = sin(t279);
t297 = t194 * t226 + t216 * t200 + t311 * t267;
t269 = cos(t279);
t296 = t194 * t227 - t216 * t199 + t311 * t269;
t186 = (t276 * t315 + t348) * t285 + t304;
t266 = qJDD(4) * t286;
t295 = (-t185 * t226 - t186 * t227 + t199 * t214 - t200 * t213) * MDP(18) + (t185 * t227 - t199 * t213) * MDP(17) + t340 * MDP(19) + t309 * MDP(20) + 0.2e1 * (t334 * t276 * qJD(4) - t274 * t342) * MDP(11) + (t274 * t278 - 0.2e1 * t282 * t320) * MDP(10) + (-t286 * t290 - t327) * MDP(13) + (-t282 * t290 + t266) * MDP(12) + t274 * MDP(4);
t294 = t240 * t276 - t247 * t290 + t250 * t274 + t311;
t206 = t289 * t274 + t317;
t203 = t286 * t206;
t183 = -t217 * t331 + qJDD(4) * pkin(4) + t203 + (t276 * t331 - t347) * pkin(8);
t293 = t201 * t329 + g(3) * t269 + (-t201 * t275 - t183) * t281 + t216 * t214 - t358 * t267;
t292 = t312 * t276 - t289 * t290 + t311 + t352;
t184 = -t301 * pkin(8) + t206 * t282 + t217 * t330;
t291 = g(3) * t267 + t310 * qJD(5) + t285 * t183 - t281 * t184 - t216 * t213 + t358 * t269;
t288 = cos(qJ(1));
t284 = sin(qJ(1));
t265 = pkin(4) * t330;
t264 = pkin(8) * t331;
t241 = qJD(3) + t265;
t237 = t250 + t271;
t225 = -pkin(2) * t276 + t312;
t224 = t240 + t265;
t223 = qJD(4) * t234;
t222 = -t289 * t331 + t264;
t208 = qJD(4) * t221 + t282 * t323;
t207 = -t247 * t331 + t286 * t323 + t264;
t205 = t209 * t286;
t1 = [t295 + (t224 * t214 + t237 * t186 + (-t308 * qJD(5) + t207 * t285 - t208 * t281) * t275 + t307 * t273 + t297) * MDP(22) + ((-t274 * t283 - t276 * t332) * pkin(1) - t313) * MDP(6) + ((t274 * t287 - t321) * pkin(1) - t314) * MDP(5) + qJDD(1) * MDP(1) + (t224 * t213 + t237 * t185 - (t307 * qJD(5) + t207 * t281 + t208 * t285) * t275 - t308 * t273 + t296) * MDP(23) + (g(1) * t284 - g(2) * t288) * MDP(2) + (g(1) * t288 + g(2) * t284) * MDP(3) + (t209 * t250 + t230 * t240 + t212 * t260 + t225 * t323 - g(1) * (-pkin(1) * t284 + t316) - g(2) * (pkin(1) * t288 + t337)) * MDP(9) + ((t302 * qJD(4) + t328) * t286 + t294 * t282 + t339) * MDP(15) + (t205 + (-t328 + (-t230 - t302) * qJD(4)) * t282 + t294 * t286) * MDP(16) + (pkin(1) * t321 + qJDD(3) + (-pkin(2) + t260) * t274 + t314) * MDP(7) + ((qJD(3) + t240) * t276 + (qJ(3) + t250) * t274 + t313) * MDP(8); t295 + (t241 * t214 + t251 * t186 + (-t306 * qJD(5) + t222 * t285 - t223 * t281) * t275 + t305 * t273 + (-t287 * t214 - t227 * t346) * t354 + t297) * MDP(22) + (0.2e1 * t352 + (0.2e1 * qJD(3) - t324) * t276 + t313) * MDP(8) + (t241 * t213 + t251 * t185 - (t305 * qJD(5) + t222 * t281 + t223 * t285) * t275 - t306 * t273 + (-t287 * t213 + t226 * t346) * t354 + t296) * MDP(23) + (t209 * qJ(3) + t230 * qJD(3) - t212 * pkin(2) - g(1) * t316 - g(2) * t337 + (-t225 * t283 - t230 * t287) * t354) * MDP(9) + (t276 * t324 - t313) * MDP(6) + t298 * MDP(5) + (qJDD(3) - t298 - 0.2e1 * t357) * MDP(7) + ((t303 * qJD(4) + t326) * t286 + t292 * t282 + t339) * MDP(15) + (t205 + (-t326 + (-t230 - t303) * qJD(4)) * t282 + t292 * t286) * MDP(16); t274 * MDP(7) - t272 * MDP(8) + (t212 + t300) * MDP(9) + (t335 * t282 + t266) * MDP(15) + (t335 * t286 - t327) * MDP(16) + (-t214 * t276 + t340) * MDP(22) + (-t213 * t276 + t309) * MDP(23); MDP(12) * t347 - MDP(13) * t348 + qJDD(4) * MDP(14) + (g(3) * t282 + t300 * t286 + t203) * MDP(15) + (g(3) * t286 + (-t206 - t300) * t282) * MDP(16) + (-(-t202 * t281 - t350) * t275 + (-t214 * t345 + t285 * t273 - t275 * t329) * pkin(4) + t291) * MDP(22) + ((-qJD(5) * t197 + t202 * t275 - t184) * t285 + (-qJD(5) * t285 * t275 - t213 * t345 - t281 * t273) * pkin(4) + t293) * MDP(23) + t299 + (MDP(10) * t342 - t334 * MDP(11)) * t272; (-t310 * t275 + t291) * MDP(22) + ((-t184 + (-qJD(5) + t275) * t197) * t285 + t293) * MDP(23) + t299;];
tau = t1;
