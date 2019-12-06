% Calculate minimal parameter regressor of Coriolis joint torque vector for
% S5PRRRP7
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,alpha2,d2,d3,d4,theta1]';
% MDP [20x1]
%   Minimal dynamic parameter vector (fixed base model)
%   see S5PRRRP7_convert_par2_MPV_fixb.m
% 
% Output:
% tauc [5x1]
%   joint torques required to compensate Coriolis and centrifugal load

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-05 16:57
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauc = S5PRRRP7_coriolisvecJ_fixb_mdp_slag_vp(qJ, qJD, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(9,1),zeros(20,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PRRRP7_coriolisvecJ_fixb_mdp_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5PRRRP7_coriolisvecJ_fixb_mdp_slag_vp: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5PRRRP7_coriolisvecJ_fixb_mdp_slag_vp: pkin has to be [9x1] (double)');
assert(isreal(MDP) && all(size(MDP) == [20 1]), ...
  'S5PRRRP7_coriolisvecJ_fixb_mdp_slag_vp: MDP has to be [20x1] (double)'); 

%% Symbolic Calculation
% From coriolisvec_joint_fixb_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 16:56:32
% EndTime: 2019-12-05 16:56:37
% DurationCPUTime: 1.78s
% Computational Cost: add. (1213->262), mult. (3115->385), div. (0->0), fcn. (2136->8), ass. (0->128)
t245 = sin(qJ(3));
t329 = MDP(5) * t245;
t240 = t245 ^ 2;
t248 = cos(qJ(3));
t328 = (-t248 ^ 2 + t240) * MDP(6);
t259 = pkin(3) * t245 - pkin(8) * t248;
t225 = t259 * qJD(3);
t229 = -pkin(3) * t248 - pkin(8) * t245 - pkin(2);
t244 = sin(qJ(4));
t246 = sin(qJ(2));
t247 = cos(qJ(4));
t282 = qJD(4) * t247;
t242 = sin(pkin(5));
t291 = qJD(1) * t242;
t249 = cos(qJ(2));
t304 = t248 * t249;
t327 = -(t244 * t246 + t247 * t304) * t291 + t244 * t225 + t229 * t282;
t274 = t246 * t291;
t227 = qJD(2) * pkin(7) + t274;
t243 = cos(pkin(5));
t290 = qJD(1) * t248;
t202 = -t245 * t227 + t243 * t290;
t285 = qJD(3) * t245;
t323 = pkin(7) * t244;
t326 = t247 * t225 + t285 * t323 - (-t244 * t304 + t246 * t247) * t291;
t286 = qJD(3) * t244;
t288 = qJD(2) * t245;
t222 = t247 * t288 + t286;
t325 = t222 ^ 2;
t324 = pkin(4) * t245;
t322 = -qJ(5) - pkin(8);
t321 = qJD(2) * pkin(2);
t320 = qJ(5) * t245;
t310 = t243 * t245;
t235 = qJD(1) * t310;
t289 = qJD(2) * t242;
t270 = t249 * t289;
t260 = t245 * t270;
t284 = qJD(3) * t248;
t183 = qJD(1) * t260 + qJD(3) * t235 + t227 * t284;
t319 = t183 * t244;
t318 = t183 * t247;
t194 = -qJD(3) * pkin(3) - t202;
t317 = t194 * t244;
t283 = qJD(4) * t244;
t267 = t245 * t283;
t278 = qJD(2) * qJD(3);
t265 = t248 * t278;
t276 = qJD(3) * qJD(4);
t294 = (t265 + t276) * t247;
t197 = qJD(2) * t267 - t294;
t316 = t197 * t244;
t287 = qJD(2) * t248;
t236 = -qJD(4) + t287;
t315 = t222 * t236;
t314 = t236 * t244;
t313 = t236 * t247;
t312 = t242 * t246;
t311 = t242 * t249;
t273 = t249 * t291;
t204 = t229 * qJD(2) - t273;
t309 = t244 * t204;
t308 = t244 * t248;
t307 = t245 * t247;
t250 = qJD(3) ^ 2;
t306 = t245 * t250;
t305 = t247 * t248;
t303 = t248 * t250;
t203 = t248 * t227 + t235;
t195 = qJD(3) * pkin(8) + t203;
t177 = -t195 * t244 + t247 * t204;
t173 = -qJ(5) * t222 + t177;
t168 = -pkin(4) * t236 + t173;
t302 = t168 - t173;
t237 = pkin(7) * t305;
t256 = -qJ(5) * t305 + t324;
t281 = qJD(5) * t247;
t301 = -t245 * t281 + t256 * qJD(3) + (-t237 + (-t229 + t320) * t244) * qJD(4) + t326;
t300 = (-pkin(7) * qJD(3) - qJ(5) * qJD(4)) * t307 + (-qJD(5) * t245 + (-pkin(7) * qJD(4) - qJ(5) * qJD(3)) * t248) * t244 + t327;
t224 = t259 * qJD(2);
t299 = t247 * t202 + t244 * t224;
t263 = qJD(4) * t322;
t298 = t281 - t299 + (qJ(5) * t287 + t263) * t244;
t261 = -t202 * t244 + t247 * t224;
t297 = -t256 * qJD(2) - qJD(5) * t244 + t247 * t263 - t261;
t293 = t244 * t229 + t237;
t280 = t194 * qJD(4);
t279 = t247 * qJD(3);
t277 = qJD(3) * MDP(16);
t275 = pkin(4) * t244 + pkin(7);
t271 = t246 * t289;
t269 = t244 * t284;
t268 = t236 * t283;
t266 = t245 * t282;
t264 = t244 * t276;
t182 = -t227 * t285 + (qJD(3) * t243 + t270) * t290;
t201 = (t225 + t274) * qJD(2);
t262 = t244 * t182 - t247 * t201;
t228 = -t273 - t321;
t258 = -t228 - t273;
t178 = t195 * t247 + t309;
t257 = qJD(2) * t240 - t236 * t248;
t198 = t264 + (t266 + t269) * qJD(2);
t175 = pkin(4) * t198 + t183;
t209 = t248 * t312 + t310;
t187 = -t209 * t244 - t247 * t311;
t255 = -t209 * t247 + t244 * t311;
t208 = -t243 * t248 + t245 * t312;
t254 = -t247 * t182 + t195 * t283 - t244 * t201 - t204 * t282;
t253 = qJD(3) * (-t258 - t321);
t252 = -t178 * qJD(4) - t262;
t251 = qJD(2) ^ 2;
t231 = t322 * t247;
t230 = t322 * t244;
t220 = t244 * t288 - t279;
t219 = t247 * t229;
t217 = t220 ^ 2;
t189 = -t244 * t320 + t293;
t186 = t209 * qJD(3) + t260;
t185 = -t208 * qJD(3) + t248 * t270;
t184 = -qJ(5) * t307 + t219 + (-pkin(4) - t323) * t248;
t180 = pkin(4) * t220 + qJD(5) + t194;
t174 = -qJ(5) * t220 + t178;
t171 = t187 * qJD(4) + t185 * t247 + t244 * t271;
t170 = t255 * qJD(4) - t185 * t244 + t247 * t271;
t167 = -qJ(5) * t198 - qJD(5) * t220 - t254;
t166 = qJ(5) * t197 - qJD(5) * t222 + t278 * t324 + t252;
t1 = [(-t170 * t236 + t186 * t220 + t198 * t208) * MDP(17) + (t171 * t236 + t186 * t222 - t197 * t208) * MDP(18) + (-t170 * t222 - t171 * t220 + t187 * t197 + t198 * t255) * MDP(19) + (t166 * t187 - t167 * t255 + t168 * t170 + t171 * t174 + t175 * t208 + t180 * t186) * MDP(20) + (-t186 * MDP(10) - t185 * MDP(11) + (MDP(17) * t187 + MDP(18) * t255) * t288) * qJD(3) + ((-t245 * MDP(10) - t248 * MDP(11)) * t249 * t278 + (-t249 * MDP(4) + (-MDP(10) * t248 + MDP(11) * t245 - MDP(3)) * t246) * t251) * t242; 0.2e1 * t265 * t329 - 0.2e1 * t278 * t328 + MDP(7) * t303 - MDP(8) * t306 + (-pkin(7) * t303 + t245 * t253) * MDP(10) + (pkin(7) * t306 + t248 * t253) * MDP(11) + (-t197 * t307 + (t248 * t279 - t267) * t222) * MDP(12) + ((-t220 * t247 - t222 * t244) * t284 + (t316 - t198 * t247 + (t220 * t244 - t222 * t247) * qJD(4)) * t245) * MDP(13) + (t236 * t267 + t197 * t248 + (t222 * t245 + t257 * t247) * qJD(3)) * MDP(14) + (t236 * t266 + t198 * t248 + (-t220 * t245 - t257 * t244) * qJD(3)) * MDP(15) + (-t236 - t287) * t245 * t277 + ((t229 * t283 - t326) * t236 + ((pkin(7) * t220 + t317) * qJD(3) + (t309 + (pkin(7) * t236 + t195) * t247) * qJD(4) + t262) * t248 + (-t220 * t273 + t247 * t280 + pkin(7) * t198 + t319 + ((-pkin(7) * t308 + t219) * qJD(2) + t177) * qJD(3)) * t245) * MDP(17) + (t327 * t236 + (t194 * t279 + (qJD(3) * t222 - t268) * pkin(7) - t254) * t248 + (-t222 * t273 - t244 * t280 - pkin(7) * t197 + t318 + (-pkin(7) * t313 - t293 * qJD(2) - t178) * qJD(3)) * t245) * MDP(18) + (t184 * t197 - t189 * t198 - t301 * t222 - t300 * t220 + (-t168 * t247 - t174 * t244) * t284 + (-t166 * t247 - t167 * t244 + (t168 * t244 - t174 * t247) * qJD(4)) * t245) * MDP(19) + (t166 * t184 + t167 * t189 + t300 * t174 + t301 * t168 + t180 * t275 * t284 + (t175 * t275 + (pkin(4) * t282 - t273) * t180) * t245) * MDP(20); (qJD(3) * t203 - t228 * t288 - t183) * MDP(10) + t258 * t287 * MDP(11) + (-t222 * t313 - t316) * MDP(12) + ((t220 * t236 - t197) * t247 + (-t198 + t315) * t244) * MDP(13) + (-t236 * t282 + (t236 * t305 + (-t222 + t286) * t245) * qJD(2)) * MDP(14) + (t268 + (-t236 * t308 + (t220 + t279) * t245) * qJD(2)) * MDP(15) + t236 * MDP(16) * t288 + (-pkin(3) * t198 - t318 + t261 * t236 - t203 * t220 + (pkin(8) * t313 + t317) * qJD(4) + (-t177 * t245 + (-pkin(8) * t285 - t194 * t248) * t244) * qJD(2)) * MDP(17) + (pkin(3) * t197 + t319 - t299 * t236 - t203 * t222 + (-pkin(8) * t314 + t194 * t247) * qJD(4) + (-t194 * t305 + (-pkin(8) * t279 + t178) * t245) * qJD(2)) * MDP(18) + (t197 * t230 + t198 * t231 - t297 * t222 - t298 * t220 + (t236 * t168 + t167) * t247 + (t236 * t174 - t166) * t244) * MDP(19) + (-t167 * t231 + t166 * t230 + t175 * (-pkin(4) * t247 - pkin(3)) + (-pkin(4) * t314 - t203) * t180 + t298 * t174 + t297 * t168) * MDP(20) + (-t248 * t329 + t328) * t251; (-t217 + t325) * MDP(13) + t294 * MDP(14) + (-t264 - t315) * MDP(15) + (-t178 * t236 - t194 * t222 + t252) * MDP(17) + (-t177 * t236 + t254) * MDP(18) + t302 * MDP(20) * t174 + (t197 * MDP(19) + (-t180 * t222 + t166) * MDP(20)) * pkin(4) + (t222 * MDP(12) - t236 * MDP(14) + t194 * MDP(18) - t302 * MDP(19)) * t220 + (-MDP(15) * t269 + (t277 + (-t244 * MDP(14) - t247 * MDP(15)) * qJD(4)) * t245) * qJD(2); (-t217 - t325) * MDP(19) + (t168 * t222 + t174 * t220 + t175) * MDP(20);];
tauc = t1;
