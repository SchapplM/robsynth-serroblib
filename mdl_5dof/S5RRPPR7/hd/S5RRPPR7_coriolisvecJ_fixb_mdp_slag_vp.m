% Calculate Coriolis joint torque vector for
% S5RRPPR7
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d2,d5,theta3]';
% MDP [25x1]
%   Minimal dynamic parameter vector (fixed base model)
%   see S5RRPPR7_convert_par2_MPV_fixb.m
% 
% Output:
% tauc [5x1]
%   joint torques required to compensate Coriolis and centrifugal load

% Quelle: HybrDyn-Toolbox
% Datum: 2021-01-15 20:00
% Revision: 24b2e7d74a0c1a3b64fa2f8f5ad758691ad61af3 (2021-01-15)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauc = S5RRPPR7_coriolisvecJ_fixb_mdp_slag_vp(qJ, qJD, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(8,1),zeros(25,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRPPR7_coriolisvecJ_fixb_mdp_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RRPPR7_coriolisvecJ_fixb_mdp_slag_vp: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RRPPR7_coriolisvecJ_fixb_mdp_slag_vp: pkin has to be [8x1] (double)');
assert(isreal(MDP) && all(size(MDP) == [25 1]), ...
  'S5RRPPR7_coriolisvecJ_fixb_mdp_slag_vp: MDP has to be [25x1] (double)'); 

%% Symbolic Calculation
% From coriolisvec_joint_fixb_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2021-01-15 19:59:06
% EndTime: 2021-01-15 19:59:16
% DurationCPUTime: 2.76s
% Computational Cost: add. (1498->292), mult. (3879->381), div. (0->0), fcn. (2625->6), ass. (0->134)
t252 = sin(pkin(8));
t255 = sin(qJ(2));
t290 = qJD(1) * t255;
t253 = cos(pkin(8));
t257 = cos(qJ(2));
t298 = t253 * t257;
t221 = -qJD(1) * t298 + t252 * t290;
t256 = cos(qJ(5));
t254 = sin(qJ(5));
t289 = qJD(2) * t254;
t199 = -t256 * t221 + t289;
t300 = t252 * t257;
t232 = t253 * t255 + t300;
t292 = qJD(1) * t232;
t322 = qJD(5) + t292;
t323 = t199 * t322;
t201 = qJD(2) * t256 + t221 * t254;
t273 = t322 * t201;
t282 = qJD(1) * qJD(2);
t321 = -0.2e1 * t282;
t320 = MDP(4) * t255;
t319 = t254 * t322;
t318 = (t255 ^ 2 - t257 ^ 2) * MDP(5);
t309 = -qJ(3) - pkin(6);
t238 = t309 * t257;
t235 = qJD(1) * t238;
t227 = t252 * t235;
t237 = t309 * t255;
t234 = qJD(1) * t237;
t196 = t234 * t253 + t227;
t285 = -qJD(4) + t196;
t230 = qJD(2) * pkin(2) + t234;
t299 = t253 * t235;
t194 = t252 * t230 - t299;
t191 = -qJD(2) * qJ(4) - t194;
t311 = pkin(4) * t221;
t175 = -t191 - t311;
t195 = t234 * t252 - t299;
t280 = t255 * t282;
t239 = t252 * t280;
t279 = t257 * t282;
t214 = t253 * t279 - t239;
t246 = -pkin(2) * t253 - pkin(3);
t242 = -pkin(7) + t246;
t316 = t242 * t214 + (t175 - t195 + t311) * t322;
t220 = t292 ^ 2;
t314 = pkin(3) + pkin(7);
t313 = pkin(2) * t255;
t223 = t232 * qJD(2);
t213 = qJD(1) * t223;
t312 = pkin(3) * t213;
t310 = pkin(4) * t292;
t231 = t252 * t255 - t298;
t247 = -pkin(2) * t257 - pkin(1);
t268 = -qJ(4) * t232 + t247;
t177 = t314 * t231 + t268;
t308 = t177 * t214;
t278 = qJD(2) * t309;
t218 = qJD(3) * t257 + t255 * t278;
t208 = t218 * qJD(1);
t219 = -qJD(3) * t255 + t257 * t278;
t209 = t219 * qJD(1);
t178 = t208 * t252 - t253 * t209;
t197 = -t253 * t237 - t238 * t252;
t307 = t178 * t197;
t287 = qJD(5) * t254;
t286 = qJD(5) * t256;
t294 = t254 * t213 + t221 * t286;
t182 = -qJD(2) * t287 + t294;
t306 = t182 * t256;
t291 = qJD(1) * t247;
t236 = qJD(3) + t291;
t261 = -qJ(4) * t292 + t236;
t184 = pkin(3) * t221 + t261;
t305 = t184 * t292;
t304 = t199 * t221;
t303 = t201 * t221;
t302 = t214 * t254;
t301 = t231 * t254;
t258 = qJD(2) ^ 2;
t297 = t255 * t258;
t207 = t256 * t214;
t296 = t257 * t258;
t259 = qJD(1) ^ 2;
t295 = t257 * t259;
t179 = t253 * t208 + t252 * t209;
t288 = qJD(2) * t255;
t284 = t310 - t285;
t281 = MDP(11) - MDP(16);
t249 = pkin(2) * t288;
t248 = pkin(2) * t290;
t277 = pkin(1) * t321;
t243 = pkin(2) * t280;
t276 = -qJ(4) * t214 + t243;
t275 = qJ(4) * t221 + t248;
t186 = t218 * t252 - t253 * t219;
t193 = t230 * t253 + t227;
t274 = t256 * t322;
t272 = MDP(24) * t322;
t271 = qJD(4) - t193;
t176 = -qJD(2) * qJD(4) - t179;
t169 = t314 * t221 + t261;
t170 = -t314 * qJD(2) + t271 + t310;
t162 = t169 * t256 + t170 * t254;
t270 = t169 * t254 - t170 * t256;
t187 = t218 * t253 + t219 * t252;
t198 = t237 * t252 - t238 * t253;
t267 = t223 * t254 + t231 * t286;
t266 = -qJD(4) * t292 + t276;
t226 = qJD(2) * t298 - t252 * t288;
t265 = -qJ(4) * t226 - qJD(4) * t232 + t249;
t264 = -qJD(2) * t195 + t178;
t165 = -pkin(4) * t213 - t176;
t263 = t165 + (-qJD(5) * t242 + t292 * t314 + t275) * t322;
t188 = pkin(4) * t232 + t197;
t262 = t165 * t231 + t175 * t223 - t188 * t214;
t260 = t178 * t232 + t186 * t292 - t187 * t221 + t197 * t214 - t198 * t213;
t244 = pkin(2) * t252 + qJ(4);
t216 = qJD(2) * t221;
t206 = t256 * t213;
t192 = pkin(3) * t231 + t268;
t190 = -qJD(2) * pkin(3) + t271;
t189 = -pkin(4) * t231 + t198;
t185 = pkin(3) * t292 + t275;
t183 = t201 * qJD(5) - t206;
t174 = pkin(3) * t223 + t265;
t173 = -pkin(4) * t223 + t187;
t172 = pkin(4) * t226 + t186;
t168 = t266 + t312;
t167 = pkin(4) * t214 + t178;
t166 = t256 * t167;
t164 = t314 * t223 + t265;
t163 = t314 * t213 + t266;
t1 = [0.2e1 * t279 * t320 + t318 * t321 + MDP(6) * t296 - MDP(7) * t297 + (-pkin(6) * t296 + t255 * t277) * MDP(9) + (pkin(6) * t297 + t257 * t277) * MDP(10) + (t213 * t247 + t223 * t236 + (-t186 + (qJD(1) * t231 + t221) * t313) * qJD(2)) * MDP(11) + (t214 * t247 + t226 * t236 + (0.2e1 * t292 * t313 - t187) * qJD(2)) * MDP(12) + (-t179 * t231 - t193 * t226 - t194 * t223 + t260) * MDP(13) + (t307 + t179 * t198 - t186 * t193 + t187 * t194 + (t236 + t291) * t249) * MDP(14) + (t176 * t231 + t190 * t226 + t191 * t223 + t260) * MDP(15) + (qJD(2) * t186 - t168 * t231 - t174 * t221 - t184 * t223 - t192 * t213) * MDP(16) + (qJD(2) * t187 - t168 * t232 - t174 * t292 - t184 * t226 - t192 * t214) * MDP(17) + (t168 * t192 + t174 * t184 - t176 * t198 + t186 * t190 - t187 * t191 + t307) * MDP(18) + (t182 * t301 + t267 * t201) * MDP(19) + ((-t199 * t254 + t201 * t256) * t223 + (t306 - t183 * t254 + (-t199 * t256 - t201 * t254) * qJD(5)) * t231) * MDP(20) + (t182 * t232 + t201 * t226 + t214 * t301 + t267 * t322) * MDP(21) + (t231 * t207 - t183 * t232 - t199 * t226 + (t223 * t256 - t231 * t287) * t322) * MDP(22) + (t214 * t232 + t226 * t322) * MDP(23) + (-t270 * t226 + t166 * t232 + t173 * t199 + t189 * t183 + (-t163 * t232 - t164 * t322 - t308) * t254 + (t172 * t322 - t262) * t256 + ((-t177 * t256 - t188 * t254) * t322 - t162 * t232 + t175 * t301) * qJD(5)) * MDP(24) + (-t162 * t226 + t173 * t201 + t189 * t182 + (-(qJD(5) * t188 + t164) * t322 - t308 - (qJD(5) * t170 + t163) * t232 + t175 * qJD(5) * t231) * t256 + (-(-qJD(5) * t177 + t172) * t322 - (-qJD(5) * t169 + t167) * t232 + t262) * t254) * MDP(25); -t295 * t320 + t259 * t318 + (-t221 * t248 - t236 * t292 - t264) * MDP(11) + (qJD(2) * t196 + t221 * t236 - t248 * t292 - t179) * MDP(12) + ((t194 - t195) * t292 + (-t193 + t196) * t221 + (-t213 * t252 - t214 * t253) * pkin(2)) * MDP(13) + (t193 * t195 - t194 * t196 + (-t178 * t253 + t179 * t252 - t236 * t290) * pkin(2)) * MDP(14) + (-t213 * t244 + t214 * t246 + (-t191 - t195) * t292 + (t190 + t285) * t221) * MDP(15) + (t185 * t221 + t264 + t305) * MDP(16) + (-t184 * t221 + t185 * t292 + (0.2e1 * qJD(4) - t196) * qJD(2) + t179) * MDP(17) + (-t176 * t244 + t178 * t246 - t184 * t185 - t190 * t195 + t285 * t191) * MDP(18) + (-t254 * t273 + t306) * MDP(19) + ((-t183 - t273) * t256 + (-t182 + t323) * t254) * MDP(20) + (-t319 * t322 + t207 + t303) * MDP(21) + (-t274 * t322 - t302 - t304) * MDP(22) + t322 * t221 * MDP(23) + (t244 * t183 + t284 * t199 - t221 * t270 + t263 * t254 + t256 * t316) * MDP(24) + (-t162 * t221 + t244 * t182 + t284 * t201 - t254 * t316 + t263 * t256) * MDP(25) + (t259 * t255 * MDP(9) + MDP(10) * t295) * pkin(1); -t239 * MDP(12) + (t193 * t292 + t194 * t221 + t243) * MDP(14) + (t216 + t239) * MDP(17) + (t312 - t191 * t221 + (-qJD(4) - t190) * t292 + t276) * MDP(18) + (-t302 + t304) * MDP(24) + (-t207 + t303) * MDP(25) + (MDP(25) * t319 - t256 * t272) * t322 + (-t221 * MDP(12) + t281 * t292 + (t281 * t300 + ((MDP(12) - MDP(17)) * t257 + t281 * t255) * t253) * qJD(1)) * qJD(2) + (MDP(13) + MDP(15)) * (-t221 ^ 2 - t220); (t214 + t216) * MDP(15) - t292 * t221 * MDP(16) + (-t220 - t258) * MDP(17) + (qJD(2) * t191 + t178 + t305) * MDP(18) + (-qJD(2) * t199 + t207) * MDP(24) + (-qJD(2) * t201 - t302) * MDP(25) + (-MDP(25) * t274 - t254 * t272) * t322; t201 * t199 * MDP(19) + (-t199 ^ 2 + t201 ^ 2) * MDP(20) + (t294 + t323) * MDP(21) + (t206 + t273) * MDP(22) + t214 * MDP(23) + (t162 * t322 - t163 * t254 - t175 * t201 + t166) * MDP(24) + (-t163 * t256 - t167 * t254 + t175 * t199 - t270 * t322) * MDP(25) + (-MDP(21) * t289 - t201 * MDP(22) - t162 * MDP(24) + t270 * MDP(25)) * qJD(5);];
tauc = t1;
