% Calculate vector of inverse dynamics joint torques for
% S5PRRRR4
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
%   pkin=[a2,a3,a4,a5,d2,d3,d4,d5,theta1]';
% MDP [21x1]
%   Minimal dynamic parameter vector (fixed base model)
%   see S5PRRRR4_convert_par2_MPV_fixb.m
% 
% Output:
% tau [5x1]
%   joint torques of inverse dynamics (contains inertial, gravitational Coriolis and centrifugal forces)

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-05 17:08
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tau = S5PRRRR4_invdynJ_fixb_mdp_slag_vp(qJ, qJD, qJDD, g, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(5,1),zeros(3,1),zeros(9,1),zeros(21,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PRRRR4_invdynJ_fixb_mdp_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5PRRRR4_invdynJ_fixb_mdp_slag_vp: qJD has to be [5x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [5 1]), ...
  'S5PRRRR4_invdynJ_fixb_mdp_slag_vp: qJDD has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5PRRRR4_invdynJ_fixb_mdp_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5PRRRR4_invdynJ_fixb_mdp_slag_vp: pkin has to be [9x1] (double)');
assert(isreal(MDP) && all(size(MDP) == [21 1]), ...
  'S5PRRRR4_invdynJ_fixb_mdp_slag_vp: MDP has to be [21x1] (double)'); 

%% Symbolic Calculation
% From invdyn_joint_fixb_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 17:08:03
% EndTime: 2019-12-05 17:08:06
% DurationCPUTime: 1.43s
% Computational Cost: add. (1079->209), mult. (1566->287), div. (0->0), fcn. (1008->12), ass. (0->124)
t350 = pkin(7) + pkin(8);
t273 = sin(qJ(5));
t274 = sin(qJ(4));
t276 = cos(qJ(5));
t277 = cos(qJ(4));
t213 = t273 * t277 + t274 * t276;
t269 = qJD(2) + qJD(3);
t204 = t213 * t269;
t275 = sin(qJ(3));
t347 = pkin(2) * t275;
t312 = qJD(3) * t347;
t278 = cos(qJ(3));
t346 = pkin(2) * t278;
t324 = qJD(2) * t312 - qJDD(2) * t346;
t266 = qJDD(2) + qJDD(3);
t345 = pkin(3) * t266;
t206 = t324 - t345;
t267 = pkin(9) + qJ(2);
t259 = qJ(3) + t267;
t247 = cos(t259);
t343 = g(2) * t247;
t349 = t206 + t343;
t328 = t276 * t277;
t332 = t273 * t274;
t212 = -t328 + t332;
t246 = sin(t259);
t325 = g(1) * t247 + g(2) * t246;
t341 = pkin(2) * qJD(2);
t314 = t275 * t341;
t303 = t350 * t269 + t314;
t196 = t277 * qJD(1) - t303 * t274;
t268 = qJD(4) + qJD(5);
t197 = qJD(1) * t274 + t303 * t277;
t344 = pkin(3) * t269;
t243 = g(1) * t246;
t251 = pkin(7) + t347;
t342 = -pkin(8) - t251;
t340 = t197 * t276;
t272 = qJ(4) + qJ(5);
t261 = sin(t272);
t339 = t246 * t261;
t262 = cos(t272);
t338 = t246 * t262;
t337 = t247 * t261;
t336 = t247 * t262;
t335 = t266 * t277;
t334 = t268 * t278;
t333 = t269 * t274;
t329 = t274 * t277;
t327 = qJDD(1) - g(3);
t321 = qJD(2) * t278;
t313 = pkin(2) * t321;
t226 = -t313 - t344;
t318 = qJD(4) * t274;
t326 = t226 * t318 + t277 * t243;
t270 = t274 ^ 2;
t323 = -t277 ^ 2 + t270;
t320 = qJD(3) * t278;
t319 = qJD(4) * t269;
t317 = qJD(4) * t277;
t316 = qJD(5) * t273;
t315 = qJDD(2) * t275;
t311 = pkin(2) * t320;
t310 = pkin(4) * t318;
t309 = t269 * t332;
t308 = t269 * t328;
t307 = t226 * t317 + t349 * t274;
t253 = -pkin(4) * t277 - pkin(3);
t306 = qJD(4) * t350;
t305 = t269 * t317;
t207 = pkin(7) * t266 + (qJD(2) * t320 + t315) * pkin(2);
t304 = pkin(8) * t266 + t207;
t301 = qJD(4) * t342;
t300 = t269 * t314;
t299 = t212 * t266;
t192 = qJD(4) * pkin(4) + t196;
t297 = -t192 * t273 - t340;
t194 = t268 * t212;
t265 = qJDD(4) + qJDD(5);
t296 = -t194 * t268 + t213 * t265;
t210 = t342 * t274;
t263 = t277 * pkin(8);
t211 = t251 * t277 + t263;
t295 = t210 * t276 - t211 * t273;
t294 = t210 * t273 + t211 * t276;
t240 = t350 * t274;
t241 = pkin(7) * t277 + t263;
t293 = -t240 * t276 - t241 * t273;
t292 = -t240 * t273 + t241 * t276;
t252 = -pkin(3) - t346;
t279 = qJD(4) ^ 2;
t291 = t251 * t279 + t252 * t266;
t290 = t243 - t324 - t343;
t190 = (t269 * t318 - t335) * pkin(4) + t206;
t205 = t253 * t269 - t313;
t289 = -g(1) * t339 + g(2) * t337 + t190 * t213 - t205 * t194;
t195 = t268 * t213;
t288 = g(1) * t338 - g(2) * t336 + t190 * t212 + t205 * t195;
t287 = -qJDD(4) * t251 + t252 * t319;
t180 = qJD(5) * t308 + t213 * t266 - t268 * t309 + t276 * t305;
t202 = -t308 + t309;
t286 = t204 * t202 * MDP(15) + (t202 * t268 + t180) * MDP(17) - t299 * MDP(18) + (-t202 ^ 2 + t204 ^ 2) * MDP(16) + t265 * MDP(19);
t285 = -t226 * t269 - t207 + t325;
t181 = t195 * t269 + t299;
t185 = -t195 * t268 - t212 * t265;
t229 = qJDD(4) * t274 + t277 * t279;
t230 = qJDD(4) * t277 - t274 * t279;
t284 = (-t180 * t212 - t181 * t213 + t194 * t202 - t195 * t204) * MDP(16) + (t180 * t213 - t194 * t204) * MDP(15) + t296 * MDP(17) + t185 * MDP(18) + 0.2e1 * (t266 * t329 - t323 * t319) * MDP(9) + (t266 * t270 + 0.2e1 * t274 * t305) * MDP(8) + t229 * MDP(10) + t230 * MDP(11) + t266 * MDP(5);
t283 = pkin(7) * t279 - t300 - t345;
t282 = -pkin(7) * qJDD(4) + (t313 - t344) * qJD(4);
t258 = t277 * qJDD(1);
t178 = qJDD(4) * pkin(4) - t197 * qJD(4) - t304 * t274 + t258;
t281 = t205 * t202 + t197 * t316 + g(2) * t338 + g(1) * t336 + g(3) * t261 + (-t197 * t268 - t178) * t273;
t179 = t196 * qJD(4) + t274 * qJDD(1) + t304 * t277;
t280 = g(1) * t337 + g(2) * t339 - g(3) * t262 + t297 * qJD(5) + t276 * t178 - t273 * t179 - t205 * t204;
t257 = cos(t267);
t256 = sin(t267);
t233 = t253 - t346;
t224 = t310 + t312;
t223 = t277 * t306;
t222 = t274 * t306;
t199 = -t274 * t311 + t277 * t301;
t198 = t274 * t301 + t277 * t311;
t1 = [t327 * MDP(1) + t230 * MDP(13) - t229 * MDP(14) + t185 * MDP(20) - t296 * MDP(21); t290 * MDP(6) + t325 * MDP(7) + ((-t291 - t349) * MDP(13) + t287 * MDP(14)) * t277 + (t287 * MDP(13) + (t291 - t243) * MDP(14)) * t274 + t284 + (t224 * t202 + t233 * t181 + (-t294 * qJD(5) - t198 * t273 + t199 * t276) * t268 + t295 * t265 + t288) * MDP(20) + (t266 * t278 * MDP(6) + (-qJDD(2) - t266) * MDP(7) * t275 + ((-MDP(13) * t277 + MDP(14) * t274 - MDP(6)) * t275 * t269 + ((-qJD(2) - t269) * MDP(7) + (-t274 * MDP(13) - t277 * MDP(14)) * qJD(4)) * t278) * qJD(3)) * pkin(2) + (g(1) * t256 - g(2) * t257) * MDP(3) + (g(1) * t257 + g(2) * t256) * MDP(4) + (t224 * t204 + t233 * t180 - (t295 * qJD(5) + t198 * t276 + t199 * t273) * t268 - t294 * t265 + t289) * MDP(21) + t326 * MDP(13) + t307 * MDP(14) + qJDD(2) * MDP(2); (t282 * t274 + (-t283 - t349) * t277 + t326) * MDP(13) + (t282 * t277 + (t283 - t243) * t274 + t307) * MDP(14) + (t202 * t310 + t253 * t181 + (-t292 * qJD(5) + t222 * t273 - t223 * t276) * t268 + t293 * t265 + (-t275 * t202 + t213 * t334) * t341 + t288) * MDP(20) + t284 + (t204 * t310 + t253 * t180 - (t293 * qJD(5) - t222 * t276 - t223 * t273) * t268 - t292 * t265 + (-t275 * t204 - t212 * t334) * t341 + t289) * MDP(21) + (t290 + t300) * MDP(6) + ((-t315 + (-qJD(3) + t269) * t321) * pkin(2) + t325) * MDP(7); t274 * t266 * MDP(10) + MDP(11) * t335 + qJDD(4) * MDP(12) + (-g(3) * t277 + t285 * t274 + t258) * MDP(13) + (-t327 * t274 + t285 * t277) * MDP(14) + (-(-t196 * t273 - t340) * t268 + (-t202 * t333 + t276 * t265 - t268 * t316) * pkin(4) + t280) * MDP(20) + ((-qJD(5) * t192 + t196 * t268 - t179) * t276 + (-qJD(5) * t276 * t268 - t204 * t333 - t273 * t265) * pkin(4) + t281) * MDP(21) + t286 + (-MDP(8) * t329 + t323 * MDP(9)) * t269 ^ 2; (-t297 * t268 + t280) * MDP(20) + ((-t179 + (-qJD(5) + t268) * t192) * t276 + t281) * MDP(21) + t286;];
tau = t1;
