% Calculate minimal parameter regressor of Coriolis joint torque vector for
% S5RPRPR15
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d3,d5,theta4]';
% MDP [24x1]
%   Minimal dynamic parameter vector (fixed base model)
%   see S5RPRPR15_convert_par2_MPV_fixb.m
% 
% Output:
% tauc [5x1]
%   joint torques required to compensate Coriolis and centrifugal load

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 18:37
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauc = S5RPRPR15_coriolisvecJ_fixb_mdp_slag_vp(qJ, qJD, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(8,1),zeros(24,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPRPR15_coriolisvecJ_fixb_mdp_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RPRPR15_coriolisvecJ_fixb_mdp_slag_vp: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RPRPR15_coriolisvecJ_fixb_mdp_slag_vp: pkin has to be [8x1] (double)');
assert(isreal(MDP) && all(size(MDP) == [24 1]), ...
  'S5RPRPR15_coriolisvecJ_fixb_mdp_slag_vp: MDP has to be [24x1] (double)'); 

%% Symbolic Calculation
% From coriolisvec_joint_fixb_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 18:37:29
% EndTime: 2019-12-31 18:37:34
% DurationCPUTime: 2.59s
% Computational Cost: add. (1279->267), mult. (2828->398), div. (0->0), fcn. (1799->6), ass. (0->128)
t258 = sin(pkin(8));
t259 = cos(pkin(8));
t260 = sin(qJ(5));
t262 = cos(qJ(5));
t335 = -t258 * t260 + t262 * t259;
t333 = t335 * qJD(5);
t261 = sin(qJ(3));
t314 = qJD(1) * t261;
t252 = qJD(5) + t314;
t263 = cos(qJ(3));
t313 = qJD(1) * t263;
t299 = t258 * t313;
t306 = t259 * qJD(3);
t231 = t299 - t306;
t298 = t259 * t313;
t311 = qJD(3) * t258;
t233 = t298 + t311;
t278 = t231 * t260 - t233 * t262;
t337 = t252 * t278;
t303 = 0.2e1 * qJD(1);
t336 = MDP(8) * (t261 ^ 2 - t263 ^ 2);
t264 = -pkin(1) - pkin(6);
t334 = qJD(1) * t264;
t309 = qJD(3) * t263;
t236 = t258 * t262 + t259 * t260;
t273 = t236 * qJD(5);
t332 = pkin(7) + qJ(4);
t331 = qJD(3) * pkin(3);
t266 = qJD(1) ^ 2;
t330 = qJ(2) * t266;
t218 = t262 * t231;
t328 = t233 * t260;
t189 = t218 + t328;
t329 = t189 * t252;
t251 = qJD(2) + t334;
t327 = t251 * t263;
t325 = t258 * t263;
t324 = t259 * t263;
t242 = t261 * t251;
t323 = t261 * t264;
t265 = qJD(3) ^ 2;
t320 = t264 * t265;
t285 = pkin(3) * t263 + qJ(4) * t261;
t215 = qJD(3) * t285 - qJD(4) * t263 + qJD(2);
t200 = t215 * qJD(1);
t217 = (qJD(4) + t327) * qJD(3);
t179 = t258 * t200 + t259 * t217;
t243 = pkin(3) * t261 - qJ(4) * t263 + qJ(2);
t225 = t243 * qJD(1);
t226 = qJD(3) * qJ(4) + t242;
t187 = t258 * t225 + t259 * t226;
t274 = t335 * t261;
t319 = qJD(1) * t274 + t333;
t272 = t236 * qJD(1);
t318 = t261 * t272 + t273;
t238 = t285 * qJD(1);
t197 = t258 * t238 + t251 * t324;
t304 = qJD(1) * qJD(3);
t295 = t261 * t304;
t286 = t262 * t295;
t287 = t260 * t295;
t317 = -t258 * t286 - t259 * t287;
t308 = qJD(3) * t264;
t296 = t263 * t308;
t193 = t258 * t215 + t259 * t296;
t202 = t258 * t243 + t259 * t323;
t315 = -t265 - t266;
t312 = qJD(3) * t252;
t310 = qJD(3) * t261;
t305 = t261 * MDP(14);
t302 = pkin(7) * t259 * t261;
t300 = t258 * t314;
t297 = t258 * t310;
t294 = MDP(22) * t309;
t293 = -t258 * t264 + pkin(4);
t292 = pkin(4) * t258 - t264;
t291 = -qJD(4) + t331;
t178 = t259 * t200 - t217 * t258;
t269 = (pkin(4) * t263 + t302) * qJD(1);
t172 = qJD(3) * t269 + t178;
t288 = t258 * t295;
t176 = pkin(7) * t288 + t179;
t290 = t262 * t172 - t176 * t260;
t186 = t259 * t225 - t226 * t258;
t289 = -t233 + t311;
t196 = t259 * t238 - t251 * t325;
t284 = t172 * t260 + t176 * t262;
t175 = pkin(4) * t314 - pkin(7) * t233 + t186;
t177 = -pkin(7) * t231 + t187;
t169 = t175 * t262 - t177 * t260;
t170 = t175 * t260 + t177 * t262;
t283 = -t178 * t258 + t179 * t259;
t282 = -t186 * t259 - t187 * t258;
t281 = -t186 * t258 + t187 * t259;
t229 = t259 * t243;
t188 = -pkin(7) * t324 + t261 * t293 + t229;
t195 = -pkin(7) * t325 + t202;
t280 = t188 * t262 - t195 * t260;
t279 = t188 * t260 + t195 * t262;
t277 = (-t231 - t306) * t261;
t219 = -t291 - t327;
t248 = t332 * t259;
t276 = qJD(4) * t258 + qJD(5) * t248 + t196 + t269;
t247 = t332 * t258;
t275 = pkin(7) * t300 - qJD(4) * t259 + qJD(5) * t247 + t197;
t271 = t335 * qJD(1);
t270 = -t219 + (t251 + t334) * t263;
t268 = -qJD(5) * t218 + t258 * t287 - t259 * t286;
t174 = -qJD(5) * t278 + t317;
t267 = -qJ(4) * t309 + (t219 + t291) * t261;
t173 = -qJD(5) * t328 + t268;
t254 = -pkin(4) * t259 - pkin(3);
t237 = t251 * t310;
t230 = t292 * t263;
t220 = t292 * t310;
t214 = t335 * t263;
t213 = t236 * t263;
t208 = -pkin(4) * t300 + t242;
t205 = t259 * t215;
t203 = -pkin(4) * t288 + t237;
t201 = -t258 * t323 + t229;
t194 = pkin(4) * t231 + t219;
t192 = -t258 * t296 + t205;
t184 = pkin(7) * t297 + t193;
t183 = -t260 * t261 * t306 - t262 * t297 + t263 * t333;
t182 = -qJD(3) * t274 - t263 * t273;
t180 = t205 + (t263 * t293 + t302) * qJD(3);
t1 = [-0.2e1 * t263 * MDP(7) * t295 + 0.2e1 * t304 * t336 + (-t261 * t320 + (qJ(2) * t309 + qJD(2) * t261) * t303) * MDP(12) + (-t263 * t320 + (-qJ(2) * t310 + qJD(2) * t263) * t303) * MDP(13) + ((qJD(1) * t192 + t178) * t261 + ((qJD(1) * t201 + t186) * t263 + (t231 * t264 + t258 * t270) * t261) * qJD(3)) * MDP(14) + ((-qJD(1) * t193 - t179) * t261 + ((-qJD(1) * t202 - t187) * t263 + (t233 * t264 + t259 * t270) * t261) * qJD(3)) * MDP(15) + (-t192 * t233 - t193 * t231 + (-t178 * t259 - t179 * t258) * t263 + ((t201 * t259 + t202 * t258) * qJD(1) - t282) * t310) * MDP(16) + (t178 * t201 + t179 * t202 + t186 * t192 + t187 * t193 + (t219 - t327) * t261 * t308) * MDP(17) + (t173 * t214 - t182 * t278) * MDP(18) + (-t173 * t213 - t174 * t214 - t182 * t189 + t183 * t278) * MDP(19) + (t173 * t261 + t182 * t252 + (qJD(1) * t214 - t278) * t309) * MDP(20) + (-t174 * t261 - t183 * t252 + (-qJD(1) * t213 - t189) * t309) * MDP(21) + (t252 + t314) * t294 + ((t180 * t262 - t184 * t260) * t252 + t290 * t261 - t220 * t189 + t230 * t174 + t203 * t213 + t194 * t183 + (-t170 * t261 - t252 * t279) * qJD(5) + (qJD(1) * t280 + t169) * t309) * MDP(23) + (-(t180 * t260 + t184 * t262) * t252 - t284 * t261 + t220 * t278 + t230 * t173 + t203 * t214 + t194 * t182 + (-t169 * t261 - t252 * t280) * qJD(5) + (-qJD(1) * t279 - t170) * t309) * MDP(24) + (MDP(6) * qJ(2) + MDP(5)) * qJD(2) * t303 + (-t263 * MDP(10) - t261 * MDP(9)) * t265; -t266 * MDP(5) - MDP(6) * t330 + t315 * t263 * MDP(13) + (-t259 * t266 + (t231 - t299) * qJD(3)) * t305 + ((-t231 * t259 + t233 * t258) * t309 + (t231 * t258 + t233 * t259) * qJD(1)) * MDP(16) + (t282 * qJD(1) + (t281 - t242) * t309) * MDP(17) + (-t252 * t271 + (-t236 * t312 - t174) * t263) * MDP(23) + (t252 * t272 + (-t312 * t335 - t173) * t263) * MDP(24) + (t315 * MDP(12) + (t258 * t266 + (t233 - t298) * qJD(3)) * MDP(15) + (qJD(3) * t219 + t283) * MDP(17) + (-t252 * t333 + (-t236 * t313 + t189) * qJD(3)) * MDP(23) + (t252 * t273 + (-t263 * t271 - t278) * qJD(3)) * MDP(24)) * t261; t261 * MDP(13) * t330 + (t251 * t277 + (-t186 * t263 - t196 * t261 + t258 * t267) * qJD(1)) * MDP(14) + (t289 * t242 + (t187 * t263 + t197 * t261 + t259 * t267) * qJD(1)) * MDP(15) + (t196 * t233 + t197 * t231 + (-qJD(4) * t231 - t186 * t314 + t179) * t259 + (qJD(4) * t233 - t187 * t314 - t178) * t258) * MDP(16) + (-t186 * t196 - t187 * t197 + (-t219 - t331) * t242 + t281 * qJD(4) + t283 * qJ(4)) * MDP(17) + (t173 * t236 - t278 * t319) * MDP(18) + (t173 * t335 - t174 * t236 - t189 * t319 + t278 * t318) * MDP(19) + (t319 * t252 + (qJD(3) * t236 + t278) * t313) * MDP(20) + (-t318 * t252 + (qJD(3) * t335 + t189) * t313) * MDP(21) - t252 * MDP(22) * t313 + (t254 * t174 - t208 * t189 - t203 * t335 + (t260 * t275 - t262 * t276) * t252 + t318 * t194 + ((-t247 * t262 - t248 * t260) * qJD(3) - t169) * t313) * MDP(23) + (t254 * t173 + t208 * t278 + t203 * t236 + (t260 * t276 + t262 * t275) * t252 + t319 * t194 + (-(-t247 * t260 + t248 * t262) * qJD(3) + t170) * t313) * MDP(24) + (-t336 + (-qJ(2) * MDP(12) + t261 * MDP(7)) * t263) * t266; (-t231 ^ 2 - t233 ^ 2) * MDP(16) + (t186 * t233 + t187 * t231 + t237) * MDP(17) + (t174 - t337) * MDP(23) + (t173 - t329) * MDP(24) + (MDP(15) * t277 - t289 * t305) * qJD(1); -t278 * t189 * MDP(18) + (-t189 ^ 2 + t278 ^ 2) * MDP(19) + (t268 + t329) * MDP(20) + (-t317 - t337) * MDP(21) + qJD(1) * t294 + (t170 * t252 + t194 * t278 + t290) * MDP(23) + (t169 * t252 + t189 * t194 - t284) * MDP(24) + (-MDP(20) * t328 + MDP(21) * t278 - MDP(23) * t170 - MDP(24) * t169) * qJD(5);];
tauc = t1;
