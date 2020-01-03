% Calculate minimal parameter regressor of Coriolis joint torque vector for
% S5RPRRP13
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [7x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d3,d4]';
% MDP [24x1]
%   Minimal dynamic parameter vector (fixed base model)
%   see S5RPRRP13_convert_par2_MPV_fixb.m
% 
% Output:
% tauc [5x1]
%   joint torques required to compensate Coriolis and centrifugal load

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 19:00
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauc = S5RPRRP13_coriolisvecJ_fixb_mdp_slag_vp(qJ, qJD, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(7,1),zeros(24,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPRRP13_coriolisvecJ_fixb_mdp_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RPRRP13_coriolisvecJ_fixb_mdp_slag_vp: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [7 1]), ...
  'S5RPRRP13_coriolisvecJ_fixb_mdp_slag_vp: pkin has to be [7x1] (double)');
assert(isreal(MDP) && all(size(MDP) == [24 1]), ...
  'S5RPRRP13_coriolisvecJ_fixb_mdp_slag_vp: MDP has to be [24x1] (double)'); 

%% Symbolic Calculation
% From coriolisvec_joint_fixb_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 18:59:41
% EndTime: 2019-12-31 18:59:46
% DurationCPUTime: 2.11s
% Computational Cost: add. (1483->278), mult. (3035->389), div. (0->0), fcn. (1627->4), ass. (0->121)
t208 = sin(qJ(3));
t264 = qJD(1) * t208;
t202 = qJD(4) + t264;
t248 = 0.2e1 * qJD(1);
t210 = cos(qJ(3));
t206 = t210 ^ 2;
t297 = MDP(8) * (t208 ^ 2 - t206);
t296 = t208 * MDP(7);
t295 = MDP(6) * qJ(2) + MDP(5);
t230 = pkin(3) * t210 + pkin(7) * t208;
t189 = qJD(3) * t230 + qJD(2);
t176 = t189 * qJD(1);
t196 = pkin(3) * t208 - pkin(7) * t210 + qJ(2);
t182 = t196 * qJD(1);
t211 = -pkin(1) - pkin(6);
t292 = qJD(1) * t211;
t201 = qJD(2) + t292;
t195 = t208 * t201;
t184 = qJD(3) * pkin(7) + t195;
t209 = cos(qJ(4));
t257 = qJD(4) * t209;
t207 = sin(qJ(4));
t258 = qJD(4) * t207;
t260 = qJD(3) * t210;
t278 = t201 * t207;
t232 = -t209 * t176 + t182 * t258 + t184 * t257 + t260 * t278;
t252 = qJD(1) * qJD(3);
t238 = t210 * t252;
t233 = pkin(4) * t238;
t151 = t232 - t233;
t161 = t182 * t207 + t184 * t209;
t157 = qJ(5) * t202 + t161;
t294 = -t157 * t202 + t151;
t273 = t208 * t209;
t293 = t207 * t196 + t211 * t273;
t256 = qJD(4) * t210;
t241 = t207 * t256;
t254 = t209 * qJD(3);
t217 = t254 * t208 + t241;
t251 = qJD(3) * qJD(4);
t168 = qJD(1) * t217 - t209 * t251;
t274 = t207 * t208;
t199 = t252 * t274;
t262 = qJD(3) * t207;
t263 = qJD(1) * t210;
t192 = t209 * t263 + t262;
t259 = qJD(4) * t192;
t169 = -t199 + t259;
t261 = qJD(3) * t208;
t152 = pkin(4) * t169 + qJ(5) * t168 - qJD(5) * t192 + t201 * t261;
t290 = t152 * t207;
t289 = t152 * t209;
t287 = t168 * t207;
t286 = t169 * t209;
t277 = t201 * t210;
t185 = -qJD(3) * pkin(3) - t277;
t285 = t185 * t207;
t190 = t207 * t263 - t254;
t284 = t190 * t202;
t283 = t190 * t207;
t282 = t190 * t209;
t281 = t192 * t207;
t280 = t192 * t209;
t279 = t196 * t209;
t276 = t202 * t207;
t275 = t202 * t209;
t272 = t209 * t210;
t213 = qJD(1) ^ 2;
t271 = t210 * t213;
t212 = qJD(3) ^ 2;
t270 = t211 * t212;
t228 = pkin(4) * t207 - qJ(5) * t209;
t269 = qJD(5) * t207 - t202 * t228 + t195;
t194 = t230 * qJD(1);
t268 = t207 * t194 + t201 * t272;
t265 = -t212 - t213;
t255 = qJD(5) * t202;
t160 = t182 * t209 - t184 * t207;
t253 = qJD(5) - t160;
t250 = MDP(19) + MDP(21);
t249 = MDP(20) - MDP(23);
t247 = pkin(7) * t276;
t246 = pkin(7) * t275;
t245 = pkin(7) * t260;
t242 = t210 * t254;
t243 = t207 * t189 + t196 * t257 + t211 * t242;
t240 = t211 * t258;
t239 = t209 * t256;
t237 = MDP(18) * t263;
t236 = t207 * t211 - pkin(4);
t235 = t190 + t254;
t234 = -t192 + t262;
t231 = qJ(5) * t238;
t229 = pkin(4) * t209 + qJ(5) * t207;
t218 = -t207 * t176 - t182 * t257 + t184 * t258 - t201 * t242;
t150 = -t218 + t231 + t255;
t227 = t150 * t209 + t151 * t207;
t156 = -pkin(4) * t202 + t253;
t226 = t156 * t209 - t157 * t207;
t225 = t156 * t207 + t157 * t209;
t224 = qJD(1) * t206 - t202 * t208;
t223 = -t211 + t228;
t159 = pkin(4) * t190 - qJ(5) * t192 + t185;
t222 = -t159 * t208 + t245;
t221 = t185 * t208 - t245;
t219 = t161 * t202 - t232;
t216 = -t207 * t250 - t209 * t249;
t215 = t160 * t202 + t218;
t214 = (t280 + t283) * MDP(22) + t226 * MDP(24) + (t207 * t249 - t209 * t250) * t202;
t197 = -pkin(3) - t229;
t172 = t223 * t210;
t167 = t208 * t236 - t279;
t166 = qJ(5) * t208 + t293;
t164 = pkin(4) * t192 + qJ(5) * t190;
t163 = -t194 * t209 + (-pkin(4) * qJD(1) + t278) * t210;
t162 = qJ(5) * t263 + t268;
t158 = -t168 + t284;
t155 = (qJD(4) * t229 - qJD(5) * t209) * t210 - t223 * t261;
t154 = qJD(4) * t293 - t189 * t209 + t236 * t260;
t153 = qJ(5) * t260 + (qJD(5) - t240) * t208 + t243;
t1 = [-0.2e1 * t238 * t296 + 0.2e1 * t252 * t297 + (-t208 * t270 + (qJ(2) * t260 + qJD(2) * t208) * t248) * MDP(12) + (-t210 * t270 + (-qJ(2) * t261 + qJD(2) * t210) * t248) * MDP(13) + (-t168 * t272 - t192 * t217) * MDP(14) + ((t281 + t282) * t261 + (t287 - t286 + (-t280 + t283) * qJD(4)) * t210) * MDP(15) + (-t202 * t241 - t168 * t208 + (t192 * t210 + t209 * t224) * qJD(3)) * MDP(16) + (-t202 * t239 - t169 * t208 + (-t190 * t210 - t207 * t224) * qJD(3)) * MDP(17) + (t202 + t264) * MDP(18) * t260 + (t189 * t275 - t232 * t208 - t210 * t211 * t169 + (t185 * t272 - t202 * t293) * qJD(4) + ((t190 * t211 - t285) * t208 + (qJD(1) * t279 + t160 + (-t211 * t202 + (t201 - t292) * t208) * t207) * t210) * qJD(3)) * MDP(19) + (-(-t208 * t240 + t243) * t202 + t218 * t208 + (t211 * t168 - t185 * t258) * t210 + ((-qJD(1) * t293 - t161) * t210 + (t211 * t192 + (-t185 + t277) * t209) * t208) * qJD(3)) * MDP(20) + (-t154 * t202 + t155 * t190 + t169 * t172 + (-t159 * t262 - t151) * t208 + (t159 * t257 + t290 + (-qJD(1) * t167 - t156) * qJD(3)) * t210) * MDP(21) + (-t153 * t190 + t154 * t192 - t166 * t169 - t167 * t168 - t226 * t261 + (-qJD(4) * t225 - t150 * t207 + t151 * t209) * t210) * MDP(22) + (t153 * t202 - t155 * t192 + t168 * t172 + (t159 * t254 + t150) * t208 + (t159 * t258 - t289 + (qJD(1) * t166 + t157) * qJD(3)) * t210) * MDP(23) + (t150 * t166 + t151 * t167 + t152 * t172 + t153 * t157 + t154 * t156 + t155 * t159) * MDP(24) + t295 * qJD(2) * t248 + (-MDP(10) * t210 - MDP(9) * t208) * t212; -t295 * t213 + t250 * t190 * t261 + t214 * qJD(1) + (t265 * MDP(13) - t152 * MDP(24) - t250 * t169 + t249 * t168 + ((t281 - t282) * MDP(22) + t225 * MDP(24) + t216 * t202) * qJD(3)) * t210 + (t265 * MDP(12) + (-t286 - t287) * MDP(22) + t227 * MDP(24) + t214 * qJD(4) + (t159 * MDP(24) + t192 * t249 + t216 * t263) * qJD(3)) * t208; t271 * t296 - t213 * t297 + (t192 * t275 - t287) * MDP(14) + ((-t168 - t284) * t209 + (-t192 * t202 - t169) * t207) * MDP(15) + (t202 * t257 + (t202 * t273 + t210 * t234) * qJD(1)) * MDP(16) + (-t202 * t258 + (-t202 * t274 + t210 * t235) * qJD(1)) * MDP(17) - t202 * t237 + (-t194 * t275 - pkin(3) * t169 + (-t208 * t235 + t210 * t276) * t201 + (-t246 + t285) * qJD(4) + (-t160 * t210 + t207 * t221) * qJD(1)) * MDP(19) + (pkin(3) * t168 + t268 * t202 + t234 * t195 + (t185 * t209 + t247) * qJD(4) + (t161 * t210 + t209 * t221) * qJD(1)) * MDP(20) + (-t289 + t163 * t202 + t169 * t197 - t269 * t190 + (t159 * t207 - t246) * qJD(4) + (t156 * t210 - t207 * t222) * qJD(1)) * MDP(21) + (t162 * t190 - t163 * t192 + (t150 + t202 * t156 + (-t169 + t259) * pkin(7)) * t209 + ((qJD(4) * t190 - t168) * pkin(7) + t294) * t207) * MDP(22) + (-t290 - t162 * t202 + t168 * t197 + t269 * t192 + (-t159 * t209 - t247) * qJD(4) + (-t157 * t210 + t209 * t222) * qJD(1)) * MDP(23) + (t152 * t197 - t156 * t163 - t157 * t162 - t269 * t159 + (qJD(4) * t226 + t227) * pkin(7)) * MDP(24) + (MDP(13) * t208 * t213 - MDP(12) * t271) * qJ(2); t158 * MDP(16) + (-qJD(1) * t239 - t207 * t251 + t199) * MDP(17) + qJD(3) * t237 + t219 * MDP(19) + t215 * MDP(20) + (t219 + 0.2e1 * t233) * MDP(21) + (pkin(4) * t168 - qJ(5) * t169) * MDP(22) + (-t215 + 0.2e1 * t231 + 0.2e1 * t255) * MDP(23) + (-pkin(4) * t151 + qJ(5) * t150 - t156 * t161 + t157 * t253 - t159 * t164) * MDP(24) + (t202 * MDP(17) - t185 * MDP(19) - t159 * MDP(21) + (t157 - t161) * MDP(22) + t164 * MDP(23) + MDP(15) * t192) * t192 + (t192 * MDP(14) + t185 * MDP(20) - t164 * MDP(21) + (t156 - t253) * MDP(22) - t159 * MDP(23) - MDP(15) * t190) * t190; (t190 * t192 - t238) * MDP(21) + t158 * MDP(22) + (-t192 ^ 2 - t202 ^ 2) * MDP(23) + (t159 * t192 + t294) * MDP(24);];
tauc = t1;
