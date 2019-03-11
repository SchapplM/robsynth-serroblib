% Calculate minimal parameter regressor of Coriolis joint torque vector for
% S6RPPRRP3
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d4,d5,theta2]';
% MDP [25x1]
%   Minimal dynamic parameter vector (fixed base model)
%   see S6RPPRRP3_convert_par2_MPV_fixb.m
% 
% Output:
% tauc [6x1]
%   joint torques required to compensate Coriolis and centrifugal load

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 02:04
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauc = S6RPPRRP3_coriolisvecJ_fixb_mdp_slag_vp(qJ, qJD, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(9,1),zeros(25,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPPRRP3_coriolisvecJ_fixb_mdp_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RPPRRP3_coriolisvecJ_fixb_mdp_slag_vp: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S6RPPRRP3_coriolisvecJ_fixb_mdp_slag_vp: pkin has to be [9x1] (double)');
assert(isreal(MDP) && all(size(MDP) == [25 1]), ...
  'S6RPPRRP3_coriolisvecJ_fixb_mdp_slag_vp: MDP has to be [25x1] (double)'); 

%% Symbolic Calculation
% From coriolisvec_joint_fixb_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 02:03:53
% EndTime: 2019-03-09 02:03:58
% DurationCPUTime: 2.02s
% Computational Cost: add. (2066->311), mult. (4071->418), div. (0->0), fcn. (2249->6), ass. (0->134)
t255 = sin(qJ(4));
t313 = qJD(1) * t255;
t243 = qJD(5) + t313;
t346 = t255 * MDP(8);
t256 = cos(qJ(5));
t296 = MDP(21) - MDP(24);
t345 = t296 * t256;
t257 = cos(qJ(4));
t251 = t257 ^ 2;
t344 = (t255 ^ 2 - t251) * MDP(9);
t242 = -cos(pkin(9)) * pkin(1) - pkin(2) - pkin(7);
t230 = qJD(1) * t242 + qJD(3);
t249 = t257 * qJD(2);
t212 = t230 * t255 + t249;
t205 = qJD(4) * pkin(8) + t212;
t342 = -qJD(2) * t255 + t230 * t257;
t206 = t342 * qJD(4);
t244 = sin(pkin(9)) * pkin(1) + qJ(3);
t227 = pkin(4) * t255 - pkin(8) * t257 + t244;
t215 = t227 * qJD(1);
t279 = pkin(4) * t257 + pkin(8) * t255;
t231 = qJD(4) * t279 + qJD(3);
t218 = t231 * qJD(1);
t254 = sin(qJ(5));
t305 = qJD(5) * t256;
t306 = qJD(5) * t254;
t281 = t205 * t305 + t206 * t254 + t215 * t306 - t218 * t256;
t299 = qJD(1) * qJD(4);
t285 = t257 * t299;
t282 = pkin(5) * t285;
t178 = t281 - t282;
t186 = t205 * t256 + t215 * t254;
t183 = qJ(6) * t243 + t186;
t343 = -t183 * t243 + t178;
t323 = t255 * t256;
t316 = t227 * t254 + t242 * t323;
t301 = t256 * qJD(4);
t304 = qJD(5) * t257;
t264 = t254 * t304 + t255 * t301;
t298 = qJD(4) * qJD(5);
t209 = qJD(1) * t264 - t256 * t298;
t312 = qJD(1) * t257;
t233 = t254 * t312 - t301;
t311 = qJD(4) * t254;
t235 = t256 * t312 + t311;
t185 = -t205 * t254 + t215 * t256;
t300 = qJD(6) - t185;
t182 = -pkin(5) * t243 + t300;
t274 = t182 * t254 + t183 * t256;
t341 = (t233 * t256 - t235 * t254) * MDP(23) - t274 * MDP(25);
t340 = 0.2e1 * qJD(3);
t310 = qJD(4) * t255;
t207 = qJD(4) * t249 + t230 * t310;
t286 = t254 * t299;
t240 = t255 * t286;
t307 = qJD(5) * t235;
t210 = -t240 + t307;
t179 = pkin(5) * t210 + qJ(6) * t209 - qJD(6) * t235 + t207;
t339 = t179 * t254;
t338 = t179 * t256;
t204 = -qJD(4) * pkin(4) - t342;
t336 = t204 * t254;
t335 = t204 * t256;
t334 = t207 * t254;
t333 = t207 * t256;
t332 = t209 * t254;
t331 = t227 * t256;
t329 = t231 * t256;
t328 = t233 * t243;
t239 = qJD(1) * t244;
t327 = t239 * MDP(7);
t326 = t243 * t254;
t325 = t243 * t256;
t324 = t254 * t255;
t199 = t255 * t210;
t258 = qJD(4) ^ 2;
t322 = t255 * t258;
t321 = t256 * t257;
t320 = t257 * t258;
t277 = pkin(5) * t254 - qJ(6) * t256;
t319 = qJD(6) * t254 - t243 * t277 + t212;
t237 = t279 * qJD(1);
t318 = t237 * t254 + t256 * t342;
t309 = qJD(4) * t257;
t317 = -t255 * t209 + t235 * t309;
t259 = qJD(1) ^ 2;
t314 = -t258 - t259;
t308 = qJD(5) * t233;
t303 = qJD(6) * t243;
t302 = t179 * MDP(25);
t297 = MDP(20) + MDP(22);
t295 = pkin(8) * t326;
t294 = pkin(8) * t325;
t293 = pkin(8) * t309;
t292 = t243 * t324;
t291 = -t206 * t256 - t215 * t305 - t218 * t254;
t290 = t242 * t257 * t301 + t227 * t305 + t231 * t254;
t289 = t233 * t310;
t288 = t243 * t305;
t287 = t256 * t304;
t284 = MDP(19) * t309;
t283 = t242 * t254 - pkin(5);
t280 = qJ(6) * t285;
t278 = pkin(5) * t256 + qJ(6) * t254;
t266 = t205 * t306 + t291;
t177 = -t266 + t280 + t303;
t276 = t177 * t256 + t178 * t254;
t275 = t182 * t256 - t183 * t254;
t273 = t237 * t256 - t254 * t342;
t270 = -t242 + t277;
t184 = pkin(5) * t233 - qJ(6) * t235 + t204;
t269 = -t184 * t255 + t293;
t268 = t204 * t255 - t293;
t267 = t186 * t243 - t281;
t265 = t233 * t309 - t251 * t286 + t199;
t263 = -t254 * t297 - t345;
t262 = t185 * t243 + t266;
t261 = (t233 * t254 + t235 * t256) * MDP(23) + t275 * MDP(25) + (t254 * t296 - t256 * t297) * t243;
t260 = t261 * qJD(5);
t241 = t256 * t251 * t299;
t238 = -pkin(4) - t278;
t229 = qJD(4) * t292;
t208 = t270 * t257;
t196 = pkin(5) * t235 + qJ(6) * t233;
t195 = t210 * t321;
t193 = t255 * t283 - t331;
t192 = qJ(6) * t255 + t316;
t190 = -t209 + t328;
t189 = -pkin(5) * t312 - t273;
t188 = qJ(6) * t312 + t318;
t187 = (qJD(5) * t278 - qJD(6) * t256) * t257 - t270 * t310;
t181 = qJD(5) * t316 + t283 * t309 - t329;
t180 = qJ(6) * t309 + (-t242 * t306 + qJD(6)) * t255 + t290;
t1 = [qJD(1) * MDP(6) * t340 + t327 * t340 - 0.2e1 * t285 * t346 + 0.2e1 * t299 * t344 - MDP(10) * t322 - MDP(11) * t320 + (t239 * t309 - t242 * t322 + (t244 * t309 + t255 * t340) * qJD(1)) * MDP(13) + (-t239 * t310 - t242 * t320 + (-t244 * t310 + t257 * t340) * qJD(1)) * MDP(14) + (-t209 * t321 - t235 * t264) * MDP(15) + (-t195 + (-t235 * t304 + t289) * t256 + (t235 * t310 + (t209 + t308) * t257) * t254) * MDP(16) + (-t243 * t264 + t241 + t317) * MDP(17) + (-t243 * t287 - t199 + t229 + (-qJD(1) * t251 * t254 - t233 * t257) * qJD(4)) * MDP(18) + (t243 + t313) * t284 + ((-t227 * t306 + t329) * t243 + (-t242 * t288 + (t233 * t242 - t336) * qJD(4) - t281) * t255 + (t204 * t305 + t334 - t242 * t210 + (-t242 * t326 + (-t242 * t324 + t331) * qJD(1) + t185) * qJD(4)) * t257) * MDP(20) + (-t290 * t243 + ((t242 * t243 + t205) * t306 + (t235 * t242 - t335) * qJD(4) + t291) * t255 + (-t204 * t306 + t333 + t242 * t209 + (-qJD(1) * t316 - t186) * qJD(4)) * t257) * MDP(21) + (-t181 * t243 + t187 * t233 + t208 * t210 + (-t184 * t311 - t178) * t255 + (t184 * t305 + t339 + (-qJD(1) * t193 - t182) * qJD(4)) * t257) * MDP(22) + (-t180 * t233 + t181 * t235 - t192 * t210 - t193 * t209 - t275 * t310 + (-qJD(5) * t274 - t177 * t254 + t178 * t256) * t257) * MDP(23) + (t180 * t243 - t187 * t235 + t208 * t209 + (t184 * t301 + t177) * t255 + (t184 * t306 - t338 + (qJD(1) * t192 + t183) * qJD(4)) * t257) * MDP(24) + (t177 * t192 + t178 * t193 + t179 * t208 + t180 * t183 + t181 * t182 + t184 * t187) * MDP(25); (t229 + t265) * MDP(20) + (-t241 + t317) * MDP(21) + t265 * MDP(22) - t195 * MDP(23) + t241 * MDP(24) + (t258 * MDP(14) + t209 * MDP(24) + t302 + ((MDP(22) * t254 + t345) * t243 + t341) * qJD(4)) * t255 + (-t258 * MDP(13) - MDP(23) * t332 - qJD(4) * t235 * MDP(24) + (qJD(4) * t184 + t276) * MDP(25) + t260) * t257; -t259 * MDP(6) + t297 * t289 + (t261 - t327) * qJD(1) + (t314 * MDP(14) - t302 - t297 * t210 + t296 * t209 + (t243 * t263 - t341) * qJD(4)) * t257 + (t314 * MDP(13) + (-t210 * t256 - t332) * MDP(23) + t276 * MDP(25) + t260 + (t184 * MDP(25) + t235 * t296 + t263 * t312) * qJD(4)) * t255; (qJD(4) * t212 - t239 * t312 - t207) * MDP(13) + t239 * t313 * MDP(14) + (t235 * t325 - t332) * MDP(15) + ((-t209 - t328) * t256 + (-t235 * t243 - t210) * t254) * MDP(16) + (t288 + (t243 * t323 + (-t235 + t311) * t257) * qJD(1)) * MDP(17) + (-t243 * t306 + (-t292 + (t233 + t301) * t257) * qJD(1)) * MDP(18) - t243 * MDP(19) * t312 + (-pkin(4) * t210 - t333 - t273 * t243 - t212 * t233 + (-t294 + t336) * qJD(5) + (-t185 * t257 + t254 * t268) * qJD(1)) * MDP(20) + (pkin(4) * t209 + t334 + t318 * t243 - t212 * t235 + (t295 + t335) * qJD(5) + (t186 * t257 + t256 * t268) * qJD(1)) * MDP(21) + (-t338 + t189 * t243 + t210 * t238 - t319 * t233 + (t184 * t254 - t294) * qJD(5) + (t182 * t257 - t254 * t269) * qJD(1)) * MDP(22) + (t188 * t233 - t189 * t235 + (t177 + t243 * t182 + (-t210 + t307) * pkin(8)) * t256 + ((-t209 + t308) * pkin(8) + t343) * t254) * MDP(23) + (-t339 - t188 * t243 + t209 * t238 + t319 * t235 + (-t184 * t256 - t295) * qJD(5) + (-t183 * t257 + t256 * t269) * qJD(1)) * MDP(24) + (t179 * t238 - t182 * t189 - t183 * t188 - t319 * t184 + (qJD(5) * t275 + t276) * pkin(8)) * MDP(25) + (t257 * t346 - t344) * t259; t190 * MDP(17) + (-qJD(1) * t287 - t254 * t298 + t240) * MDP(18) + qJD(1) * t284 + t267 * MDP(20) + t262 * MDP(21) + (t267 + 0.2e1 * t282) * MDP(22) + (pkin(5) * t209 - qJ(6) * t210) * MDP(23) + (-t262 + 0.2e1 * t280 + 0.2e1 * t303) * MDP(24) + (-pkin(5) * t178 + qJ(6) * t177 - t182 * t186 + t183 * t300 - t184 * t196) * MDP(25) + (t243 * MDP(18) - t204 * MDP(20) - t184 * MDP(22) + (t183 - t186) * MDP(23) + t196 * MDP(24) + MDP(16) * t235) * t235 + (t235 * MDP(15) + t204 * MDP(21) - t196 * MDP(22) + (t182 - t300) * MDP(23) - t184 * MDP(24) - MDP(16) * t233) * t233; (t233 * t235 - t285) * MDP(22) + t190 * MDP(23) + (-t235 ^ 2 - t243 ^ 2) * MDP(24) + (t184 * t235 + t343) * MDP(25);];
tauc  = t1;
