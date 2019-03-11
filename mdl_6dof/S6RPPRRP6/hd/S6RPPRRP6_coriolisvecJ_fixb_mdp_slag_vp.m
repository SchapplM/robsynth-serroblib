% Calculate minimal parameter regressor of Coriolis joint torque vector for
% S6RPPRRP6
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d4,d5]';
% MDP [27x1]
%   Minimal dynamic parameter vector (fixed base model)
%   see S6RPPRRP6_convert_par2_MPV_fixb.m
% 
% Output:
% tauc [6x1]
%   joint torques required to compensate Coriolis and centrifugal load

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 02:11
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauc = S6RPPRRP6_coriolisvecJ_fixb_mdp_slag_vp(qJ, qJD, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(8,1),zeros(27,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPPRRP6_coriolisvecJ_fixb_mdp_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RPPRRP6_coriolisvecJ_fixb_mdp_slag_vp: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S6RPPRRP6_coriolisvecJ_fixb_mdp_slag_vp: pkin has to be [8x1] (double)');
assert(isreal(MDP) && all(size(MDP) == [27 1]), ...
  'S6RPPRRP6_coriolisvecJ_fixb_mdp_slag_vp: MDP has to be [27x1] (double)'); 

%% Symbolic Calculation
% From coriolisvec_joint_fixb_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 02:11:26
% EndTime: 2019-03-09 02:11:30
% DurationCPUTime: 2.43s
% Computational Cost: add. (1840->316), mult. (3498->420), div. (0->0), fcn. (1833->4), ass. (0->134)
t243 = sin(qJ(4));
t305 = qJD(1) * t243;
t233 = qJD(5) + t305;
t241 = pkin(1) + qJ(3);
t245 = cos(qJ(4));
t224 = pkin(4) * t243 - pkin(8) * t245 + t241;
t204 = t224 * qJD(1) - qJD(2);
t235 = qJ(2) * qJD(1) + qJD(3);
t230 = -pkin(7) * qJD(1) + t235;
t223 = t243 * t230;
t211 = qJD(4) * pkin(8) + t223;
t242 = sin(qJ(5));
t244 = cos(qJ(5));
t183 = t204 * t244 - t211 * t242;
t292 = qJD(6) - t183;
t179 = -pkin(5) * t233 + t292;
t184 = t204 * t242 + t211 * t244;
t180 = qJ(6) * t233 + t184;
t262 = t179 * t242 + t180 * t244;
t302 = qJD(4) * t242;
t304 = qJD(1) * t245;
t220 = t244 * t304 + t302;
t321 = t220 * t242;
t294 = t244 * qJD(4);
t218 = t242 * t304 - t294;
t322 = t218 * t244;
t287 = MDP(23) - MDP(26);
t288 = MDP(22) + MDP(24);
t335 = t288 * t242 + t287 * t244;
t249 = t335 * t233 + (-t321 + t322) * MDP(25) - t262 * MDP(27);
t240 = -pkin(7) + qJ(2);
t300 = qJD(4) * t245;
t340 = qJD(2) * t243 + t240 * t300;
t291 = qJD(1) * qJD(2);
t203 = t230 * t300 + t243 * t291;
t269 = pkin(4) * t245 + pkin(8) * t243;
t216 = t269 * qJD(4) + qJD(3);
t205 = t216 * qJD(1);
t297 = qJD(5) * t244;
t298 = qJD(5) * t242;
t272 = -t242 * t203 - t204 * t298 + t244 * t205 - t211 * t297;
t290 = qJD(1) * qJD(4);
t277 = t245 * t290;
t273 = pkin(5) * t277;
t174 = -t272 - t273;
t339 = -t180 * t233 + t174;
t338 = qJD(1) * t241;
t239 = t245 ^ 2;
t337 = (t243 ^ 2 - t239) * MDP(11);
t336 = MDP(6) * qJ(2) + MDP(5) + MDP(7);
t296 = qJD(5) * t245;
t279 = t242 * t296;
t252 = t243 * t294 + t279;
t289 = qJD(4) * qJD(5);
t192 = t252 * qJD(1) - t244 * t289;
t313 = t242 * t243;
t228 = t290 * t313;
t299 = qJD(5) * t220;
t193 = -t228 + t299;
t301 = qJD(4) * t243;
t202 = t230 * t301 - t245 * t291;
t175 = pkin(5) * t193 + qJ(6) * t192 - qJD(6) * t220 + t202;
t332 = t175 * t242;
t331 = t175 * t244;
t329 = t192 * t242;
t328 = t193 * t244;
t327 = t202 * t242;
t326 = t202 * t244;
t316 = t230 * t245;
t212 = -qJD(4) * pkin(4) - t316;
t325 = t212 * t242;
t324 = t212 * t244;
t323 = t218 * t242;
t320 = t220 * t244;
t319 = t220 * t245;
t222 = t269 * qJD(1);
t318 = t222 * t244;
t317 = t224 * t244;
t315 = t233 * t244;
t314 = t240 * t242;
t312 = t243 * t244;
t311 = t244 * t245;
t267 = pkin(5) * t242 - qJ(6) * t244;
t310 = qJD(6) * t242 - t233 * t267 + t223;
t309 = t242 * t222 + t230 * t311;
t308 = t242 * t224 + t240 * t312;
t246 = qJD(4) ^ 2;
t247 = qJD(1) ^ 2;
t306 = -t246 - t247;
t295 = qJD(6) * t233;
t231 = -qJD(2) + t338;
t293 = qJD(2) - t231;
t286 = pkin(8) * t233 * t242;
t285 = pkin(8) * t315;
t284 = pkin(8) * t300;
t282 = 0.2e1 * qJD(3) * qJD(1);
t281 = -t244 * t203 - t204 * t297 - t242 * t205;
t278 = t244 * t296;
t276 = MDP(21) * t300;
t275 = t218 * t288;
t274 = t288 * t244;
t271 = t242 * t216 + t224 * t297 + t340 * t244;
t270 = qJ(6) * t277;
t268 = pkin(5) * t244 + qJ(6) * t242;
t266 = qJD(2) + t231 + t338;
t254 = t211 * t298 + t281;
t173 = -t254 + t270 + t295;
t265 = t173 * t244 + t174 * t242;
t264 = -t173 * t242 + t174 * t244;
t263 = t179 * t244 - t180 * t242;
t260 = qJD(1) * t239 - t233 * t243;
t259 = -t240 * t246 + t282;
t258 = -t240 + t267;
t182 = pkin(5) * t218 - qJ(6) * t220 + t212;
t257 = -t182 * t243 + t284;
t256 = t212 * t243 - t284;
t255 = t184 * t233 + t272;
t253 = t182 * MDP(27) + t287 * t220;
t250 = t183 * t233 + t254;
t248 = (t320 + t323) * MDP(25) + t263 * MDP(27) + (t287 * t242 - t274) * t233;
t227 = t242 * t277;
t225 = -pkin(4) - t268;
t200 = t258 * t245;
t191 = t192 * t244;
t190 = -t317 + (-pkin(5) + t314) * t243;
t189 = qJ(6) * t243 + t308;
t187 = pkin(5) * t220 + qJ(6) * t218;
t186 = -t318 + (-pkin(5) * qJD(1) + t230 * t242) * t245;
t185 = qJ(6) * t304 + t309;
t181 = t218 * t233 - t192;
t178 = -t258 * t301 + (t268 * qJD(5) - qJD(6) * t244 - qJD(2)) * t245;
t177 = -pkin(5) * t300 + (qJD(5) * t240 * t243 - t216) * t244 + (qJD(5) * t224 + t340) * t242;
t176 = qJ(6) * t300 + (-t240 * t298 + qJD(6)) * t243 + t271;
t1 = [MDP(8) * t282 + (qJD(2) * t235 + qJD(3) * t231 + (qJ(2) * qJD(2) + qJD(3) * t241) * qJD(1)) * MDP(9) + 0.2e1 * t290 * t337 - t246 * t245 * MDP(13) + t266 * t300 * MDP(15) + (t259 * t245 - t266 * t301) * MDP(16) + (-t192 * t311 - t252 * t220) * MDP(17) + ((t321 + t322) * t301 + (t329 - t328 + (-t320 + t323) * qJD(5)) * t245) * MDP(18) + (-t233 * t279 + (t260 * t244 + t319) * qJD(4)) * MDP(19) + (-t233 * t278 + (-t218 * t245 - t260 * t242) * qJD(4)) * MDP(20) + (t233 + t305) * t276 + ((t216 * t244 - t224 * t298) * t233 + (t212 * t297 - qJD(2) * t218 - t240 * t193 + t327 + (-t233 * t314 + (-t240 * t313 + t317) * qJD(1) + t183) * qJD(4)) * t245) * MDP(22) + (-t271 * t233 + (-t212 * t298 - qJD(2) * t220 + t240 * t192 + t326 + (-t308 * qJD(1) - t184) * qJD(4)) * t245) * MDP(23) + (-t177 * t233 + t178 * t218 + t193 * t200 + (t182 * t297 + t332 + (-qJD(1) * t190 - t179) * qJD(4)) * t245) * MDP(24) + (-t176 * t218 + t177 * t220 - t189 * t193 - t190 * t192 - t263 * t301 + (-t262 * qJD(5) + t264) * t245) * MDP(25) + (t176 * t233 - t178 * t220 + t192 * t200 + (t182 * t298 - t331 + (qJD(1) * t189 + t180) * qJD(4)) * t245) * MDP(26) + (t173 * t189 + t174 * t190 + t175 * t200 + t176 * t180 + t177 * t179 + t178 * t182) * MDP(27) + (-0.2e1 * MDP(10) * t277 - t246 * MDP(12) + t259 * MDP(15) - t192 * MDP(19) - t193 * MDP(20) + ((-qJD(2) * t242 - t240 * t297) * t233 + (t218 * t240 - t325) * qJD(4) + t272) * MDP(22) + ((t233 * t240 + t211) * t298 + (t220 * t240 - t324) * qJD(4) + t281) * MDP(23) + (-t182 * t302 - t174) * MDP(24) + (t182 * t294 + t173) * MDP(26)) * t243 + 0.2e1 * t336 * t291; (t193 * t242 - t191) * MDP(25) + t264 * MDP(27) + t287 * t227 + t275 * t304 + t249 * qJD(5) - t336 * t247 + ((-qJD(3) - t235) * MDP(9) + ((-0.2e1 * MDP(15) - t274) * qJD(4) + t253) * t245 + (0.2e1 * qJD(4) * MDP(16) + t249) * t243) * qJD(1); -t247 * MDP(8) + t275 * t301 + (t293 * MDP(9) + t248) * qJD(1) + (t306 * MDP(16) - t175 * MDP(27) - t249 * qJD(4) + t287 * t192 - t288 * t193) * t245 + (t306 * MDP(15) + (-t328 - t329) * MDP(25) + t265 * MDP(27) + t248 * qJD(5) + (-t304 * t335 + t253) * qJD(4)) * t243; (t220 * t315 - t329) * MDP(17) + (-t191 - t218 * t315 + (-t220 * t233 - t193) * t242) * MDP(18) + (t233 * t297 + t227 + (t233 * t312 - t319) * qJD(1)) * MDP(19) + (-t233 * t298 + (-t233 * t313 + (t218 + t294) * t245) * qJD(1)) * MDP(20) - t233 * MDP(21) * t304 + (-pkin(4) * t193 - t326 - (-t242 * t316 + t318) * t233 - t218 * t223 + (-t285 + t325) * qJD(5) + (-t183 * t245 + t256 * t242) * qJD(1)) * MDP(22) + (pkin(4) * t192 + t327 + t309 * t233 - t220 * t223 + (t286 + t324) * qJD(5) + (t184 * t245 + t256 * t244) * qJD(1)) * MDP(23) + (-t331 + t186 * t233 + t193 * t225 - t310 * t218 + (t182 * t242 - t285) * qJD(5) + (t179 * t245 - t257 * t242) * qJD(1)) * MDP(24) + (t185 * t218 - t186 * t220 + (t173 + t233 * t179 + (-t193 + t299) * pkin(8)) * t244 + ((qJD(5) * t218 - t192) * pkin(8) + t339) * t242) * MDP(25) + (-t332 - t185 * t233 + t192 * t225 + t310 * t220 + (-t182 * t244 - t286) * qJD(5) + (-t180 * t245 + t257 * t244) * qJD(1)) * MDP(26) + (t175 * t225 - t179 * t186 - t180 * t185 - t310 * t182 + (t263 * qJD(5) + t265) * pkin(8)) * MDP(27) + (MDP(15) * t304 - MDP(16) * t305) * t293 + (t245 * t243 * MDP(10) - t337) * t247; t181 * MDP(19) + (-qJD(1) * t278 - t242 * t289 + t228) * MDP(20) + qJD(1) * t276 + t255 * MDP(22) + t250 * MDP(23) + (t255 + 0.2e1 * t273) * MDP(24) + (pkin(5) * t192 - qJ(6) * t193) * MDP(25) + (-t250 + 0.2e1 * t270 + 0.2e1 * t295) * MDP(26) + (-pkin(5) * t174 + qJ(6) * t173 - t179 * t184 + t292 * t180 - t182 * t187) * MDP(27) + (t233 * MDP(20) - t212 * MDP(22) - t182 * MDP(24) + (t180 - t184) * MDP(25) + t187 * MDP(26) + MDP(18) * t220) * t220 + (t220 * MDP(17) + t212 * MDP(23) - t187 * MDP(24) + (t179 - t292) * MDP(25) - t182 * MDP(26) - MDP(18) * t218) * t218; (t218 * t220 - t277) * MDP(24) + t181 * MDP(25) + (-t220 ^ 2 - t233 ^ 2) * MDP(26) + (t182 * t220 + t339) * MDP(27);];
tauc  = t1;
