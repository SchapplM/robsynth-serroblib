% Calculate vector of inverse dynamics joint torques for
% S4RRRP6
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [4x1]
%   Generalized joint coordinates (joint angles)
% qJD [4x1]
%   Generalized joint velocities
% qJDD [4x1]
%   Generalized joint accelerations
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [6x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,d1,d2,d3]';
% MDP [21x1]
%   Minimal dynamic parameter vector (fixed base model)
%   see S4RRRP6_convert_par2_MPV_fixb.m
% 
% Output:
% tau [4x1]
%   joint torques of inverse dynamics (contains inertial, gravitational Coriolis and centrifugal forces)

% Quelle: HybrDyn-Toolbox
% Datum: 2021-01-15 14:39
% Revision: 24b2e7d74a0c1a3b64fa2f8f5ad758691ad61af3 (2021-01-15)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tau = S4RRRP6_invdynJ_fixb_mdp_slag_vp(qJ, qJD, qJDD, g, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(4,1),zeros(4,1),zeros(3,1),zeros(6,1),zeros(21,1)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4RRRP6_invdynJ_fixb_mdp_slag_vp: qJ has to be [4x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [4 1]), ...
  'S4RRRP6_invdynJ_fixb_mdp_slag_vp: qJD has to be [4x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [4 1]), ...
  'S4RRRP6_invdynJ_fixb_mdp_slag_vp: qJDD has to be [4x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S4RRRP6_invdynJ_fixb_mdp_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [6 1]), ...
  'S4RRRP6_invdynJ_fixb_mdp_slag_vp: pkin has to be [6x1] (double)');
assert(isreal(MDP) && all(size(MDP) == [21 1]), ...
  'S4RRRP6_invdynJ_fixb_mdp_slag_vp: MDP has to be [21x1] (double)'); 

%% Symbolic Calculation
% From invdyn_joint_fixb_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2021-01-15 14:39:21
% EndTime: 2021-01-15 14:39:27
% DurationCPUTime: 2.50s
% Computational Cost: add. (1365->328), mult. (3088->424), div. (0->0), fcn. (1895->6), ass. (0->138)
t266 = cos(qJ(2));
t307 = qJD(1) * qJD(2);
t298 = t266 * t307;
t263 = sin(qJ(2));
t305 = qJDD(1) * t263;
t353 = qJD(2) * qJD(3) + t298 + t305;
t262 = sin(qJ(3));
t265 = cos(qJ(3));
t311 = qJD(3) * t263;
t297 = qJD(1) * t311;
t280 = (-qJDD(2) + t297) * t265;
t318 = qJD(1) * t266;
t188 = ((qJD(3) + t318) * qJD(2) + t305) * t262 + t280;
t264 = sin(qJ(1));
t267 = cos(qJ(1));
t330 = t265 * t267;
t334 = t262 * t266;
t208 = t264 * t334 + t330;
t329 = t266 * t267;
t332 = t264 * t265;
t210 = -t262 * t329 + t332;
t335 = t262 * t263;
t352 = -g(1) * t210 + g(2) * t208 + g(3) * t335;
t315 = qJD(2) * t262;
t319 = qJD(1) * t263;
t226 = t265 * t319 + t315;
t350 = t226 ^ 2;
t304 = t266 * qJDD(1);
t275 = -t263 * t307 + t304;
t221 = qJDD(3) - t275;
t349 = pkin(3) * t221;
t303 = t262 * t319;
t308 = t265 * qJD(2);
t224 = t303 - t308;
t348 = pkin(3) * t224;
t347 = pkin(3) * t262;
t346 = pkin(5) * t262;
t343 = g(3) * t263;
t342 = g(3) * t266;
t261 = qJ(4) + pkin(6);
t290 = t262 * qJDD(2) + t353 * t265;
t187 = t262 * t297 - t290;
t341 = qJ(4) * t187;
t340 = qJ(4) * t188;
t339 = t187 * t262;
t244 = -qJD(3) + t318;
t338 = t224 * t244;
t337 = t226 * t244;
t336 = t226 * t265;
t333 = t263 * t265;
t331 = t265 * t266;
t231 = -pkin(2) * t266 - pkin(6) * t263 - pkin(1);
t216 = t231 * qJD(1);
t254 = pkin(5) * t318;
t237 = qJD(2) * pkin(6) + t254;
t193 = t265 * t216 - t237 * t262;
t184 = -qJ(4) * t226 + t193;
t183 = -pkin(3) * t244 + t184;
t328 = -t184 + t183;
t288 = pkin(2) * t263 - pkin(6) * t266;
t228 = t288 * qJD(1);
t212 = t262 * t228;
t295 = qJD(3) * t261;
t309 = qJD(4) * t265;
t327 = -t262 * t295 + t309 - t212 - (-pkin(5) * t333 - qJ(4) * t334) * qJD(1);
t281 = pkin(3) * t263 - qJ(4) * t331;
t322 = pkin(5) * t303 + t265 * t228;
t326 = -t281 * qJD(1) - qJD(4) * t262 - t265 * t295 - t322;
t229 = t288 * qJD(2);
t310 = qJD(3) * t265;
t325 = t262 * t229 + t231 * t310;
t314 = qJD(2) * t263;
t324 = t265 * t229 + t314 * t346;
t323 = (g(1) * t330 + g(2) * t332) * t263;
t245 = pkin(5) * t331;
t321 = t262 * t231 + t245;
t259 = t263 ^ 2;
t320 = -t266 ^ 2 + t259;
t317 = qJD(2) * t224;
t316 = qJD(2) * t226;
t313 = qJD(2) * t266;
t312 = qJD(3) * t262;
t250 = pkin(5) + t347;
t302 = t244 * t308;
t236 = -qJD(2) * pkin(2) + pkin(5) * t319;
t294 = -qJD(4) - t348;
t197 = t236 - t294;
t301 = t197 * t310;
t300 = t244 * t312;
t299 = t244 * t310;
t252 = pkin(5) * t305;
t206 = -qJDD(2) * pkin(2) + pkin(5) * t298 + t252;
t182 = pkin(3) * t188 + qJDD(4) + t206;
t296 = -t182 - t342;
t251 = pkin(3) * t265 + pkin(2);
t293 = t251 * t266 + t261 * t263;
t195 = qJD(1) * t229 + t231 * qJDD(1);
t205 = t275 * pkin(5) + qJDD(2) * pkin(6);
t291 = t262 * t195 + t265 * t205 + t216 * t310 - t237 * t312;
t289 = -qJD(3) * pkin(6) * t244 + t206;
t287 = -g(1) * t208 - g(2) * t210;
t209 = t262 * t267 - t264 * t331;
t211 = t262 * t264 + t265 * t329;
t286 = -g(1) * t209 - g(2) * t211;
t285 = g(1) * t267 + g(2) * t264;
t284 = g(1) * t264 - g(2) * t267;
t283 = -pkin(6) * t221 + qJD(3) * t236;
t194 = t216 * t262 + t237 * t265;
t185 = -qJ(4) * t224 + t194;
t282 = t183 * t265 + t185 * t262;
t279 = t285 * t263;
t278 = -0.2e1 * pkin(1) * t307 - pkin(5) * qJDD(2);
t277 = t221 * t262 - t299;
t276 = t221 * t265 + t300;
t269 = qJD(1) ^ 2;
t274 = pkin(1) * t269 + t285;
t268 = qJD(2) ^ 2;
t273 = 0.2e1 * qJDD(1) * pkin(1) - pkin(5) * t268 + t284;
t272 = g(1) * t211 - g(2) * t209 + g(3) * t333 - t291;
t191 = t265 * t195;
t271 = -t194 * qJD(3) - t262 * t205 + t191;
t270 = t271 + t352;
t248 = g(3) * t334;
t233 = t261 * t265;
t232 = t261 * t262;
t230 = t250 * t263;
t223 = t265 * t231;
t220 = t224 ^ 2;
t218 = t318 * t347 + t254;
t204 = pkin(1) + t293;
t198 = pkin(5) * t313 + (t262 * t313 + t263 * t310) * pkin(3);
t196 = -qJ(4) * t335 + t321;
t192 = -qJ(4) * t333 + t223 + (-pkin(3) - t346) * t266;
t181 = (-pkin(5) * qJD(2) - qJ(4) * qJD(3)) * t333 + (-qJD(4) * t263 + (-pkin(5) * qJD(3) - qJ(4) * qJD(2)) * t266) * t262 + t325;
t180 = -t263 * t309 + t281 * qJD(2) + (-t245 + (qJ(4) * t263 - t231) * t262) * qJD(3) + t324;
t179 = -qJD(4) * t224 + t291 - t340;
t178 = -qJD(4) * t226 + t271 + t341 + t349;
t1 = [qJDD(1) * MDP(1) + t284 * MDP(2) + t285 * MDP(3) + (qJDD(1) * t259 + 0.2e1 * t263 * t298) * MDP(4) + 0.2e1 * (t263 * t304 - t320 * t307) * MDP(5) + (qJDD(2) * t263 + t266 * t268) * MDP(6) + (qJDD(2) * t266 - t263 * t268) * MDP(7) + (t278 * t263 + t273 * t266) * MDP(9) + (-t273 * t263 + t278 * t266) * MDP(10) + (-t187 * t333 + (-t262 * t311 + t266 * t308) * t226) * MDP(11) + ((-t224 * t265 - t226 * t262) * t313 + (t339 - t188 * t265 + (t224 * t262 - t336) * qJD(3)) * t263) * MDP(12) + ((t187 - t302) * t266 + (t276 + t316) * t263) * MDP(13) + ((t244 * t315 + t188) * t266 + (-t277 - t317) * t263) * MDP(14) + (-t221 * t266 - t244 * t314) * MDP(15) + (-(-t231 * t312 + t324) * t244 + t223 * t221 + (t237 * t310 - t191 + (t299 + t317) * pkin(5) + (-pkin(5) * t221 + qJD(2) * t236 + qJD(3) * t216 + t205) * t262) * t266 + (pkin(5) * t188 + qJD(2) * t193 + t206 * t262 + t236 * t310) * t263 + t286) * MDP(16) + (t325 * t244 - t321 * t221 + (t236 * t308 + (-t300 + t316) * pkin(5) + t291) * t266 + (-t236 * t312 - t194 * qJD(2) + t206 * t265 + (-t187 - t302) * pkin(5)) * t263 + t287) * MDP(17) + (-t180 * t244 + t188 * t230 + t192 * t221 + t198 * t224 + (t197 * t315 - t178) * t266 + (qJD(2) * t183 + t182 * t262 + t301) * t263 + t286) * MDP(18) + (t181 * t244 - t187 * t230 - t196 * t221 + t198 * t226 + (t197 * t308 + t179) * t266 + (-qJD(2) * t185 + t182 * t265 - t197 * t312) * t263 + t287) * MDP(19) + (-t180 * t226 - t181 * t224 + t187 * t192 - t188 * t196 - t282 * t313 + (-t178 * t265 - t179 * t262 + (t183 * t262 - t185 * t265) * qJD(3) + t284) * t263) * MDP(20) + (t179 * t196 + t185 * t181 + t178 * t192 + t183 * t180 + t182 * t230 + t197 * t198 - g(1) * (-t204 * t264 + t250 * t267) - g(2) * (t204 * t267 + t250 * t264)) * MDP(21); MDP(6) * t305 + MDP(7) * t304 + qJDD(2) * MDP(8) + (t274 * t263 - t252 - t342) * MDP(9) + (t343 + (-pkin(5) * qJDD(1) + t274) * t266) * MDP(10) + (-t244 * t336 - t339) * MDP(11) + ((-t187 + t338) * t265 + (-t188 + t337) * t262) * MDP(12) + ((-t226 * t263 + t244 * t331) * qJD(1) + t277) * MDP(13) + ((t224 * t263 - t244 * t334) * qJD(1) + t276) * MDP(14) + t244 * MDP(15) * t319 + (-pkin(2) * t188 + t322 * t244 + t283 * t262 + (-t289 - t342) * t265 + (-t193 * t263 + (-pkin(5) * t224 - t236 * t262) * t266) * qJD(1) + t323) * MDP(16) + (pkin(2) * t187 - t212 * t244 + t248 + t283 * t265 + (-t236 * t331 + t194 * t263 + (-t226 * t266 + t244 * t333) * pkin(5)) * qJD(1) + (-t279 + t289) * t262) * MDP(17) + (-t183 * t319 - t188 * t251 - t218 * t224 - t221 * t232 + t296 * t265 - t326 * t244 + (-t197 * t318 + (t197 + t348) * qJD(3)) * t262 + t323) * MDP(18) + (t301 + t187 * t251 - t218 * t226 - t221 * t233 + t248 + t327 * t244 + (t185 * t263 - t197 * t331) * qJD(1) + (pkin(3) * qJD(3) * t226 + t182 - t279) * t262) * MDP(19) + (-t343 - t178 * t262 + t179 * t265 - t187 * t232 - t188 * t233 - t326 * t226 - t327 * t224 - t282 * qJD(3) + (t282 * qJD(1) - t285) * t266) * MDP(20) + (t179 * t233 - t178 * t232 - t182 * t251 - g(3) * t293 - t285 * (-t251 * t263 + t261 * t266) + (pkin(3) * t312 - t218) * t197 + t327 * t185 + t326 * t183) * MDP(21) + (-t263 * t266 * MDP(4) + t320 * MDP(5)) * t269; t226 * t224 * MDP(11) + (-t220 + t350) * MDP(12) + (-t187 - t338) * MDP(13) + (-t188 - t337) * MDP(14) + t221 * MDP(15) + (-t194 * t244 - t226 * t236 + t270) * MDP(16) + (-t193 * t244 + t224 * t236 + t272) * MDP(17) + (0.2e1 * t349 + t341 - t185 * t244 + (-t197 + t294) * t226 + t270) * MDP(18) + (-pkin(3) * t350 + t340 - t184 * t244 + (qJD(4) + t197) * t224 + t272) * MDP(19) + (pkin(3) * t187 - t328 * t224) * MDP(20) + (t328 * t185 + (-t197 * t226 + t178 + t352) * pkin(3)) * MDP(21); (t280 - t337) * MDP(18) + (t290 + t338) * MDP(19) + (-t220 - t350) * MDP(20) + (t183 * t226 + t185 * t224 - t279 - t296) * MDP(21) + (t353 * MDP(18) - MDP(19) * t297) * t262;];
tau = t1;
