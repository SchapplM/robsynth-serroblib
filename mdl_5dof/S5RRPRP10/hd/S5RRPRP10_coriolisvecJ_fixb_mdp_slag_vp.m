% Calculate minimal parameter regressor of Coriolis joint torque vector for
% S5RRPRP10
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [7x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d2,d4]';
% MDP [23x1]
%   Minimal dynamic parameter vector (fixed base model)
%   see S5RRPRP10_convert_par2_MPV_fixb.m
% 
% Output:
% tauc [5x1]
%   joint torques required to compensate Coriolis and centrifugal load

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 20:11
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauc = S5RRPRP10_coriolisvecJ_fixb_mdp_slag_vp(qJ, qJD, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(7,1),zeros(23,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRPRP10_coriolisvecJ_fixb_mdp_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RRPRP10_coriolisvecJ_fixb_mdp_slag_vp: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [7 1]), ...
  'S5RRPRP10_coriolisvecJ_fixb_mdp_slag_vp: pkin has to be [7x1] (double)');
assert(isreal(MDP) && all(size(MDP) == [23 1]), ...
  'S5RRPRP10_coriolisvecJ_fixb_mdp_slag_vp: MDP has to be [23x1] (double)'); 

%% Symbolic Calculation
% From coriolisvec_joint_fixb_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 20:11:09
% EndTime: 2019-12-31 20:11:14
% DurationCPUTime: 1.91s
% Computational Cost: add. (1260->271), mult. (2861->378), div. (0->0), fcn. (1533->4), ass. (0->131)
t315 = pkin(3) + pkin(6);
t250 = cos(qJ(2));
t251 = -pkin(2) - pkin(7);
t248 = sin(qJ(2));
t271 = -qJ(3) * t248 - pkin(1);
t206 = t250 * t251 + t271;
t188 = t206 * qJD(1);
t292 = qJD(1) * t248;
t236 = pkin(6) * t292;
t260 = pkin(3) * t292 + qJD(3) + t236;
t190 = qJD(2) * t251 + t260;
t247 = sin(qJ(4));
t249 = cos(qJ(4));
t172 = t188 * t249 + t190 * t247;
t280 = qJD(1) * qJD(2);
t274 = t248 * t280;
t228 = t247 * t274;
t290 = qJD(2) * t247;
t291 = qJD(1) * t250;
t209 = t249 * t291 + t290;
t285 = qJD(4) * t209;
t181 = -t228 + t285;
t275 = t247 * t291;
t288 = qJD(2) * t249;
t211 = -t275 + t288;
t234 = pkin(2) * t274;
t265 = pkin(7) * t248 - qJ(3) * t250;
t286 = qJD(3) * t248;
t255 = qJD(2) * t265 - t286;
t179 = qJD(1) * t255 + t234;
t273 = t250 * t280;
t233 = pkin(6) * t273;
t205 = pkin(3) * t273 + t233;
t267 = -t179 * t247 + t249 * t205;
t163 = pkin(4) * t273 + qJ(5) * t181 - t172 * qJD(4) - qJD(5) * t211 + t267;
t169 = -qJ(5) * t209 + t172;
t235 = qJD(4) + t292;
t322 = t169 * t235 + t163;
t283 = qJD(4) * t249;
t294 = qJD(4) * t275 + t249 * t274;
t182 = qJD(2) * t283 - t294;
t278 = -t249 * t179 - t190 * t283 - t247 * t205;
t284 = qJD(4) * t247;
t164 = -qJ(5) * t182 - qJD(5) * t209 - t188 * t284 - t278;
t312 = t188 * t247;
t171 = t249 * t190 - t312;
t168 = -qJ(5) * t211 + t171;
t167 = pkin(4) * t235 + t168;
t321 = -t167 * t235 + t164;
t320 = -0.2e1 * t280;
t309 = t211 * t235;
t245 = t248 ^ 2;
t246 = t250 ^ 2;
t317 = (t245 - t246) * MDP(5);
t226 = t315 * t248;
t295 = t249 * t206 + t247 * t226;
t316 = t211 ^ 2;
t314 = qJD(2) * pkin(2);
t313 = t181 * t249;
t289 = qJD(2) * t248;
t217 = t315 * t289;
t243 = qJD(2) * qJD(3);
t193 = -qJD(1) * t217 + t243;
t311 = t193 * t247;
t310 = t193 * t249;
t308 = t211 * t250;
t307 = t235 * t249;
t306 = t235 * t251;
t305 = t247 * t248;
t252 = qJD(2) ^ 2;
t304 = t248 * t252;
t303 = t249 * t250;
t302 = t250 * t252;
t253 = qJD(1) ^ 2;
t301 = t250 * t253;
t300 = qJ(5) - t251;
t299 = t167 - t168;
t240 = pkin(2) * t292;
t196 = qJD(1) * t265 + t240;
t237 = pkin(6) * t291;
t218 = pkin(3) * t291 + t237;
t298 = t249 * t196 + t247 * t218;
t266 = -t196 * t247 + t249 * t218;
t297 = -qJD(5) * t249 + t284 * t300 - (pkin(4) * t250 - qJ(5) * t305) * qJD(1) - t266;
t222 = t300 * t249;
t296 = -qJ(5) * t249 * t292 - qJD(4) * t222 - qJD(5) * t247 - t298;
t227 = t315 * t250;
t224 = -pkin(2) * t250 + t271;
t200 = qJD(1) * t224;
t244 = qJD(2) * qJ(3);
t287 = qJD(2) * t250;
t282 = qJD(4) * t250;
t281 = qJD(5) * t250;
t279 = t248 * t301;
t199 = t244 + t218;
t277 = t247 * t282;
t276 = t249 * t282;
t272 = MDP(19) * t291;
t270 = pkin(1) * t320;
t269 = qJD(3) - t314;
t268 = qJ(5) * t250 - t206;
t264 = -qJD(1) * t246 + t235 * t248;
t263 = -0.2e1 * qJD(2) * t200;
t262 = t235 ^ 2;
t257 = -qJ(3) * t287 - t286;
t186 = qJD(1) * t257 + t234;
t239 = pkin(2) * t289;
t198 = t239 + t257;
t259 = pkin(6) * t252 + qJD(1) * t198 + t186;
t258 = t199 * t248 + t251 * t287;
t184 = t239 + t255;
t219 = t315 * t287;
t256 = t249 * t184 - t206 * t284 + t247 * t219 + t226 * t283;
t174 = pkin(4) * t182 + t193;
t220 = pkin(6) * t274 - t243;
t223 = t236 + t269;
t225 = -t237 - t244;
t254 = -t220 * t250 + (t223 * t250 + (t225 + t237) * t248) * qJD(2);
t229 = t249 * t273;
t221 = t300 * t247;
t215 = -qJ(3) * t291 + t240;
t214 = t249 * t226;
t208 = t209 ^ 2;
t204 = t249 * t219;
t189 = t200 * t292;
t178 = pkin(4) * t209 + qJD(5) + t199;
t176 = -qJ(5) * t303 + t295;
t175 = pkin(4) * t248 + t247 * t268 + t214;
t166 = -t249 * t281 + (t288 * t248 + t277) * qJ(5) + t256;
t165 = pkin(4) * t287 + t204 + t268 * t283 + (-qJ(5) * t289 - qJD(4) * t226 - t184 + t281) * t247;
t1 = [0.2e1 * t248 * MDP(4) * t273 + t317 * t320 + MDP(6) * t302 - MDP(7) * t304 + (-pkin(6) * t302 + t248 * t270) * MDP(9) + (pkin(6) * t304 + t250 * t270) * MDP(10) + t254 * MDP(11) + (t248 * t263 + t250 * t259) * MDP(12) + (-t248 * t259 + t250 * t263) * MDP(13) + (pkin(6) * t254 + t186 * t224 + t198 * t200) * MDP(14) + (t181 * t247 * t250 + (t247 * t289 - t276) * t211) * MDP(15) + ((-t209 * t247 + t211 * t249) * t289 + (t313 + t182 * t247 + (t209 * t249 + t211 * t247) * qJD(4)) * t250) * MDP(16) + (-t235 * t276 - t181 * t248 + (t247 * t264 + t308) * qJD(2)) * MDP(17) + (t235 * t277 - t182 * t248 + (-t209 * t250 + t249 * t264) * qJD(2)) * MDP(18) + (t235 + t292) * MDP(19) * t287 + ((-t184 * t247 + t204) * t235 - t217 * t209 + t227 * t182 + (-t199 * t288 + t267) * t248 + (-t172 * t248 - t295 * t235) * qJD(4) + (-t199 * t284 + t310 + ((-t206 * t247 + t214) * qJD(1) + t171) * qJD(2)) * t250) * MDP(20) + (-t256 * t235 - t217 * t211 - t227 * t181 + ((qJD(2) * t199 + qJD(4) * t188) * t247 + t278) * t248 + (-t199 * t283 - t311 + (-qJD(1) * t295 - t172) * qJD(2)) * t250) * MDP(21) + (-t165 * t211 - t166 * t209 + t175 * t181 - t176 * t182 + (-t167 * t247 + t169 * t249) * t289 + (t163 * t247 - t164 * t249 + (t167 * t249 + t169 * t247) * qJD(4)) * t250) * MDP(22) + (t164 * t176 + t169 * t166 + t163 * t175 + t167 * t165 + t174 * (pkin(4) * t303 + t227) + (-pkin(4) * t277 + (-pkin(4) * t249 - t315) * t289) * t178) * MDP(23); -MDP(4) * t279 + t253 * t317 + ((-t225 - t244) * t248 + (-t223 + t269) * t250) * qJD(1) * MDP(11) + (-t215 * t291 + t189) * MDP(12) + (0.2e1 * t243 + (t200 * t250 + t215 * t248) * qJD(1)) * MDP(13) + (-qJ(3) * t220 - qJD(3) * t225 - t200 * t215 + (-t225 * t248 + (-t223 - t314) * t250) * qJD(1) * pkin(6)) * MDP(14) + (-t247 * t309 - t313) * MDP(15) + ((-t182 - t309) * t249 + (t209 * t235 + t181) * t247) * MDP(16) + (-t235 * t284 + t229 + (-t235 * t305 - t308) * qJD(1)) * MDP(17) + (-t235 * t283 + (-t248 * t307 + (t209 - t290) * t250) * qJD(1)) * MDP(18) - t235 * t272 + (qJ(3) * t182 + t311 - t266 * t235 + t260 * t209 + (t199 * t249 - t247 * t306) * qJD(4) + (-t171 * t250 + t249 * t258) * qJD(1)) * MDP(20) + (-qJ(3) * t181 + t310 + t298 * t235 + t260 * t211 + (-t199 * t247 - t249 * t306) * qJD(4) + (t172 * t250 - t247 * t258) * qJD(1)) * MDP(21) + (-t181 * t222 + t182 * t221 - t296 * t209 - t297 * t211 - t321 * t247 - t322 * t249) * MDP(22) + (-t164 * t221 - t163 * t222 + t174 * (pkin(4) * t247 + qJ(3)) + (pkin(4) * t307 + t260) * t178 + t296 * t169 + t297 * t167) * MDP(23) + (MDP(9) * t248 * t253 + MDP(10) * t301) * pkin(1); MDP(12) * t279 + (-t245 * t253 - t252) * MDP(13) + (t189 + t233) * MDP(14) + t229 * MDP(20) + (t225 * MDP(14) - t209 * MDP(20) - t211 * MDP(21) - t178 * MDP(23)) * qJD(2) + ((-t209 * t292 + t181 - t285) * MDP(22) + t322 * MDP(23) - MDP(21) * t262) * t249 + (-MDP(21) * t273 + (-t182 + t309) * MDP(22) + t321 * MDP(23) - MDP(20) * t262) * t247; (-t208 + t316) * MDP(16) + t228 * MDP(17) + (t294 + t309) * MDP(18) + qJD(2) * t272 + (t172 * t235 - t199 * t211 + t267) * MDP(20) + (t171 * t235 + t278) * MDP(21) + t299 * MDP(23) * t169 + (t181 * MDP(22) + (-t178 * t211 + t163) * MDP(23)) * pkin(4) + (t211 * MDP(15) + t235 * MDP(17) + t199 * MDP(21) - MDP(22) * t299) * t209 + (-MDP(17) * t209 - MDP(18) * t288 - MDP(20) * t172 + MDP(21) * t312) * qJD(4); (-t208 - t316) * MDP(22) + (t167 * t211 + t169 * t209 + t174) * MDP(23);];
tauc = t1;
