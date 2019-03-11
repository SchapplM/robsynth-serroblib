% Calculate vector of inverse dynamics joint torques for
% S6RPPPRR5
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
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d5,d6,theta4]';
% MDP [26x1]
%   Minimal dynamic parameter vector (fixed base model)
%   see S6RPPPRR5_convert_par2_MPV_fixb.m
% 
% Output:
% tau [6x1]
%   joint torques of inverse dynamics (contains inertial, gravitational Coriolis and centrifugal forces)

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 01:38
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tau = S6RPPPRR5_invdynJ_fixb_mdp_slag_vp(qJ, qJD, qJDD, g, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(3,1),zeros(9,1),zeros(26,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPPPRR5_invdynJ_fixb_mdp_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RPPPRR5_invdynJ_fixb_mdp_slag_vp: qJD has to be [6x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [6 1]), ...
  'S6RPPPRR5_invdynJ_fixb_mdp_slag_vp: qJDD has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RPPPRR5_invdynJ_fixb_mdp_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S6RPPPRR5_invdynJ_fixb_mdp_slag_vp: pkin has to be [9x1] (double)');
assert(isreal(MDP) && all(size(MDP) == [26 1]), ...
  'S6RPPPRR5_invdynJ_fixb_mdp_slag_vp: MDP has to be [26x1] (double)'); 

%% Symbolic Calculation
% From invdyn_joint_fixb_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 01:37:48
% EndTime: 2019-03-09 01:37:52
% DurationCPUTime: 2.63s
% Computational Cost: add. (1363->329), mult. (2181->427), div. (0->0), fcn. (1305->8), ass. (0->146)
t270 = cos(qJ(5));
t323 = qJD(1) * qJD(5);
t313 = t270 * t323;
t267 = sin(qJ(5));
t320 = qJDD(1) * t267;
t360 = qJD(5) * qJD(6) + t313 + t320;
t262 = sin(pkin(9));
t263 = cos(pkin(9));
t230 = qJD(2) * t263 + qJD(3) * t262;
t257 = qJD(1) * qJD(2);
t258 = qJ(2) * qJDD(1);
t310 = qJDD(3) + t257 + t258;
t222 = pkin(3) * qJDD(1) + t310;
t265 = -pkin(1) - qJ(3);
t358 = qJDD(1) * t265;
t293 = qJDD(2) + t358;
t322 = qJD(3) * qJD(1);
t223 = t293 - t322;
t340 = -t263 * t222 + t262 * t223;
t268 = sin(qJ(1));
t271 = cos(qJ(1));
t218 = t262 * t271 + t263 * t268;
t354 = g(2) * t218;
t217 = -t268 * t262 + t263 * t271;
t356 = g(1) * t217;
t359 = qJD(1) * t230 - t340 - t354 - t356;
t200 = t262 * t222 + t263 * t223;
t198 = qJDD(1) * pkin(7) + t200;
t244 = qJ(2) * qJD(1) + qJD(3);
t235 = pkin(3) * qJD(1) + t244;
t236 = qJD(1) * t265 + qJD(2);
t204 = t262 * t235 + t263 * t236;
t202 = qJD(1) * pkin(7) + t204;
t196 = qJD(4) * t267 + t202 * t270;
t332 = qJD(5) * t196;
t183 = -qJDD(5) * pkin(5) - qJDD(4) * t270 + t198 * t267 + t332;
t334 = qJD(1) * t270;
t238 = -qJD(6) + t334;
t355 = g(1) * t218;
t301 = -g(2) * t217 + t355;
t302 = pkin(5) * t267 - pkin(8) * t270;
t357 = t238 * (pkin(8) * qJD(6) + t302 * qJD(1)) + t267 * t301 - g(3) * t270 - t183;
t252 = 0.2e1 * t257;
t353 = qJDD(1) * pkin(1);
t266 = sin(qJ(6));
t269 = cos(qJ(6));
t303 = t266 * qJDD(5) + t360 * t269;
t327 = qJD(6) * t267;
t312 = qJD(1) * t327;
t191 = -t266 * t312 + t303;
t352 = t191 * t266;
t351 = t217 * t270;
t324 = t269 * qJD(5);
t335 = qJD(1) * t267;
t219 = t266 * t335 - t324;
t350 = t219 * t238;
t325 = t266 * qJD(5);
t221 = t269 * t335 + t325;
t349 = t221 * t238;
t348 = t238 * t269;
t347 = t262 * t267;
t346 = t263 * t267;
t314 = t267 * t323;
t317 = t270 * qJDD(1);
t216 = qJDD(6) + t314 - t317;
t345 = t266 * t216;
t344 = t266 * t270;
t343 = t269 * t216;
t342 = t269 * t270;
t341 = qJDD(4) - g(3);
t264 = pkin(3) + qJ(2);
t215 = t262 * t264 + t263 * t265;
t339 = t271 * pkin(1) + t268 * qJ(2);
t338 = g(1) * t268 - g(2) * t271;
t260 = t267 ^ 2;
t337 = -t270 ^ 2 + t260;
t272 = qJD(5) ^ 2;
t273 = qJD(1) ^ 2;
t336 = t272 + t273;
t195 = qJD(4) * t270 - t202 * t267;
t333 = qJD(5) * t195;
t331 = qJD(5) * t219;
t330 = qJD(5) * t221;
t329 = qJD(5) * t267;
t328 = qJD(6) * t266;
t326 = qJD(6) * t269;
t319 = qJDD(5) * t267;
t318 = qJDD(5) * t270;
t316 = t271 * qJ(3) + t339;
t315 = t238 * t325;
t311 = qJDD(2) - t338;
t189 = qJD(5) * pkin(8) + t196;
t213 = pkin(7) + t215;
t307 = t213 * t238 + t189;
t203 = t235 * t263 - t262 * t236;
t214 = -t262 * t265 + t263 * t264;
t294 = -pkin(5) * t270 - pkin(8) * t267 - pkin(4);
t190 = qJD(1) * t294 - t203;
t306 = -qJDD(5) * pkin(8) - qJD(6) * t190 - qJDD(4) * t267 - t198 * t270 - t333;
t305 = 0.2e1 * t313;
t304 = qJDD(2) - t353;
t300 = g(1) * t271 + g(2) * t268;
t299 = t238 * t324 - t191;
t246 = t269 * qJDD(5);
t192 = t269 * t312 - t246 + (t320 + (qJD(6) + t334) * qJD(5)) * t266;
t298 = -t192 + t315;
t297 = qJD(5) * t202 - t341;
t251 = t271 * qJ(2);
t296 = t265 * t268 + t251;
t295 = -qJD(6) * t189 + t354;
t292 = t302 * qJD(5);
t291 = t238 * t326 - t345;
t290 = t238 * t328 + t343;
t289 = t252 + 0.2e1 * t258 - t300;
t288 = t305 + t320;
t287 = 0.2e1 * t314 - t317;
t286 = -t270 * t336 - t319;
t285 = t267 * t336 - t318;
t283 = -t266 * t327 + t270 * t324;
t282 = t291 + t331;
t281 = -g(2) * t351 + g(3) * t267 + t306;
t201 = -qJD(1) * pkin(4) - t203;
t212 = -pkin(4) - t214;
t231 = qJD(2) * t262 - qJD(3) * t263;
t280 = -qJDD(5) * t213 + (qJD(1) * t212 + t201 - t231) * qJD(5);
t188 = -qJD(5) * pkin(5) - t195;
t279 = -pkin(8) * t216 + (-t188 - t195) * t238;
t278 = -qJD(1) * t201 - qJD(4) * qJD(5) - t198 + t301;
t277 = t188 * qJD(5) - t213 * t216 + t231 * t238 - t306;
t276 = t213 * t272 - t359 + (-pkin(4) + t212) * qJDD(1);
t275 = -t299 * t267 + (-t290 + t330) * t270;
t274 = -t267 * t298 + t270 * t282;
t233 = -t267 * t272 + t318;
t232 = t270 * t272 + t319;
t229 = qJDD(1) * t263 - t262 * t273;
t228 = -qJDD(1) * t262 - t263 * t273;
t207 = t221 * t329;
t206 = -t230 + t292;
t205 = -t214 + t294;
t194 = -t217 * t266 + t218 * t342;
t193 = -t217 * t269 - t218 * t344;
t187 = qJD(1) * t292 + qJDD(1) * t294 + t340;
t186 = t269 * t187;
t185 = t269 * t189 + t266 * t190;
t184 = -t266 * t189 + t269 * t190;
t1 = [0.2e1 * (t267 * t317 - t323 * t337) * MDP(14) + t338 * MDP(2) + (-g(1) * t296 - g(2) * t316 + qJ(2) * t310 + t244 * qJD(2) - t236 * qJD(3) + t223 * t265) * MDP(9) + (qJDD(1) * t260 + t267 * t305) * MDP(13) + (qJDD(1) * t214 + t359) * MDP(10) + (-t311 + 0.2e1 * t322 - 0.2e1 * t358) * MDP(8) + (t267 * t276 + t270 * t280) * MDP(19) + (t267 * t280 - t270 * t276) * MDP(18) + t232 * MDP(15) + t233 * MDP(16) + t300 * MDP(3) + (t200 * t215 + t204 * t231 - t340 * t214 + t203 * t230 - g(1) * (pkin(3) * t271 + t296) - g(2) * (pkin(3) * t268 + t316)) * MDP(12) + (qJDD(3) + t289) * MDP(7) + t289 * MDP(5) + (t191 * t267 * t269 + t221 * t283) * MDP(20) + (-g(2) * t194 + (t213 * t331 - t186) * t270 + (qJD(5) * t184 + t192 * t213 + t219 * t231) * t267 + (-g(1) * t351 + t205 * t216 - t206 * t238 + (t188 * t267 + t270 * t307) * qJD(6)) * t269 + (-(-qJD(6) * t205 + t213 * t329) * t238 + t183 * t267 - t355 + t277 * t270) * t266) * MDP(25) + ((t205 * t326 + t206 * t266) * t238 - t205 * t345 - t269 * t355 - g(2) * t193 + (t213 * t330 + (-qJD(6) * t307 + t187 + t356) * t266 + t277 * t269) * t270 + (-t188 * t328 + t183 * t269 + t213 * t191 + t231 * t221 + (-t213 * t348 - t185) * qJD(5)) * t267) * MDP(26) + (-t304 * pkin(1) - g(1) * (-pkin(1) * t268 + t251) - g(2) * t339 + (t252 + t258) * qJ(2)) * MDP(6) + ((-t219 * t269 - t221 * t266) * t270 * qJD(5) + (-t352 - t192 * t269 + (t219 * t266 - t221 * t269) * qJD(6)) * t267) * MDP(21) + (t311 - 0.2e1 * t353) * MDP(4) + qJDD(1) * MDP(1) + (-t216 * t270 - t238 * t329) * MDP(24) + ((t192 + t315) * t270 + (t291 - t331) * t267) * MDP(23) + (-qJD(1) * t231 - qJDD(1) * t215 - t200 + t301) * MDP(11) + (-t191 * t270 - t238 * t283 + t267 * t343 + t207) * MDP(22); (t304 - t338) * MDP(6) + ((-qJD(3) - t244) * qJD(1) + t293 - t338) * MDP(9) + t228 * MDP(10) - t229 * MDP(11) + (t340 * t262 + t200 * t263 + (-t203 * t263 - t204 * t262) * qJD(1) - t338) * MDP(12) + (t262 * t287 + t263 * t286) * MDP(18) + (t262 * t288 + t263 * t285) * MDP(19) + (t290 * t262 + t274 * t263 + ((-t262 * t344 - t263 * t269) * t238 - t219 * t347) * qJD(1)) * MDP(25) + (t291 * t262 + t275 * t263 + (-(t262 * t342 - t263 * t266) * t238 - t221 * t347) * qJD(1)) * MDP(26) + (-MDP(6) * qJ(2) - MDP(5) - MDP(7)) * t273 + (MDP(4) - MDP(8)) * qJDD(1); qJDD(1) * MDP(7) - t273 * MDP(8) + (qJD(1) * t236 - t300 + t310) * MDP(9) + t229 * MDP(10) + t228 * MDP(11) + (-t340 * t263 + t200 * t262 + (-t203 * t262 + t204 * t263) * qJD(1) - t300) * MDP(12) + (t262 * t286 - t263 * t287) * MDP(18) + (t262 * t285 - t263 * t288) * MDP(19) + (-t290 * t263 + t274 * t262 + ((-t262 * t269 + t263 * t344) * t238 + t219 * t346) * qJD(1)) * MDP(25) + (-t291 * t263 + t275 * t262 + (-(-t262 * t266 - t263 * t342) * t238 + t221 * t346) * qJD(1)) * MDP(26); t341 * MDP(12) + t233 * MDP(18) - t232 * MDP(19) + t207 * MDP(26) + (MDP(25) * t298 + MDP(26) * t299) * t270 + (MDP(25) * t282 - MDP(26) * t290) * t267; MDP(15) * t320 + MDP(16) * t317 + qJDD(5) * MDP(17) + (t267 * t278 - t270 * t297 + t332) * MDP(18) + (t267 * t297 + t270 * t278 + t333) * MDP(19) + (-t221 * t348 + t352) * MDP(20) + ((t191 + t350) * t269 + (-t192 + t349) * t266) * MDP(21) + ((-t221 * t267 + t238 * t342) * qJD(1) - t291) * MDP(22) + ((t219 * t267 - t238 * t344) * qJD(1) + t290) * MDP(23) + t238 * MDP(24) * t335 + (-pkin(5) * t192 - t184 * t335 - t196 * t219 + t279 * t266 + t357 * t269) * MDP(25) + (-pkin(5) * t191 + t185 * t335 - t196 * t221 - t357 * t266 + t279 * t269) * MDP(26) + (-MDP(13) * t267 * t270 + MDP(14) * t337) * t273; t221 * t219 * MDP(20) + (-t219 ^ 2 + t221 ^ 2) * MDP(21) + (t303 - t350) * MDP(22) + (t246 - t349) * MDP(23) + t216 * MDP(24) + (-g(1) * t193 - t185 * t238 - t188 * t221 + t186) * MDP(25) + (g(1) * t194 - t184 * t238 + t188 * t219) * MDP(26) + (-MDP(23) * t312 + MDP(25) * t295 + MDP(26) * t281) * t269 + (-MDP(22) * t312 - t360 * MDP(23) + t281 * MDP(25) + (-t187 - t295) * MDP(26)) * t266;];
tau  = t1;
