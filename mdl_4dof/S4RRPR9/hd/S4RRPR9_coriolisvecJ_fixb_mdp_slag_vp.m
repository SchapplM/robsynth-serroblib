% Calculate minimal parameter regressor of Coriolis joint torque vector for
% S4RRPR9
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [4x1]
%   Generalized joint coordinates (joint angles)
% qJD [4x1]
%   Generalized joint velocities
% pkin [7x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,d1,d2,d4,theta3]';
% MDP [21x1]
%   Minimal dynamic parameter vector (fixed base model)
%   see S4RRPR9_convert_par2_MPV_fixb.m
% 
% Output:
% tauc [4x1]
%   joint torques required to compensate Coriolis and centrifugal load

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 17:10
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauc = S4RRPR9_coriolisvecJ_fixb_mdp_slag_vp(qJ, qJD, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(4,1),zeros(7,1),zeros(21,1)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4RRPR9_coriolisvecJ_fixb_mdp_slag_vp: qJ has to be [4x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [4 1]), ...
  'S4RRPR9_coriolisvecJ_fixb_mdp_slag_vp: qJD has to be [4x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [7 1]), ...
  'S4RRPR9_coriolisvecJ_fixb_mdp_slag_vp: pkin has to be [7x1] (double)');
assert(isreal(MDP) && all(size(MDP) == [21 1]), ...
  'S4RRPR9_coriolisvecJ_fixb_mdp_slag_vp: MDP has to be [21x1] (double)'); 

%% Symbolic Calculation
% From coriolisvec_joint_fixb_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 17:10:07
% EndTime: 2019-12-31 17:10:11
% DurationCPUTime: 1.70s
% Computational Cost: add. (908->236), mult. (2415->355), div. (0->0), fcn. (1607->6), ass. (0->122)
t272 = (qJD(1) * qJD(2));
t299 = -2 * t272;
t238 = sin(qJ(2));
t298 = MDP(4) * t238;
t240 = cos(qJ(2));
t297 = MDP(5) * (t238 ^ 2 - t240 ^ 2);
t296 = pkin(6) + qJ(3);
t295 = qJD(2) * pkin(2);
t235 = sin(pkin(7));
t237 = sin(qJ(4));
t294 = t235 * t237;
t293 = t235 * t240;
t236 = cos(pkin(7));
t292 = t236 * t238;
t291 = t236 * t240;
t241 = qJD(2) ^ 2;
t290 = t238 * t241;
t289 = t240 * t241;
t242 = qJD(1) ^ 2;
t288 = t240 * t242;
t239 = cos(qJ(4));
t210 = -t239 * t236 + t294;
t246 = t210 * t240;
t287 = qJD(1) * t246 - t210 * qJD(4);
t211 = t235 * t239 + t236 * t237;
t247 = t211 * t240;
t286 = -qJD(1) * t247 + t211 * qJD(4);
t255 = pkin(2) * t238 - qJ(3) * t240;
t194 = qJD(2) * t255 - qJD(3) * t238;
t185 = t194 * qJD(1);
t283 = qJD(1) * t238;
t230 = pkin(5) * t283;
t216 = (qJD(3) - t230) * qJD(2);
t163 = t235 * t185 + t236 * t216;
t218 = -pkin(2) * t240 - qJ(3) * t238 - pkin(1);
t200 = t218 * qJD(1);
t282 = qJD(1) * t240;
t231 = pkin(5) * t282;
t222 = qJD(2) * qJ(3) + t231;
t172 = t235 * t200 + t236 * t222;
t267 = t235 * t283;
t275 = t236 * qJD(2);
t206 = t267 - t275;
t265 = t240 * t272;
t258 = t236 * t265;
t276 = qJD(4) * t239;
t285 = -t206 * t276 + t239 * t258;
t213 = t255 * qJD(1);
t177 = pkin(5) * t267 + t236 * t213;
t280 = qJD(2) * t238;
t269 = pkin(5) * t280;
t173 = t236 * t194 + t235 * t269;
t225 = pkin(5) * t291;
t182 = t235 * t218 + t225;
t281 = qJD(2) * t235;
t279 = qJD(2) * t240;
t158 = -pkin(6) * t206 + t172;
t278 = qJD(4) * t158;
t277 = qJD(4) * t238;
t274 = MDP(11) * qJD(1);
t273 = MDP(12) * qJD(1);
t271 = pkin(5) * t293;
t270 = pkin(3) * t282;
t268 = pkin(3) * t235 + pkin(5);
t266 = t236 * t283;
t264 = MDP(19) * t280;
t263 = pkin(1) * t299;
t262 = qJD(3) - t295;
t162 = t236 * t185 - t216 * t235;
t171 = t236 * t200 - t222 * t235;
t261 = t206 + t275;
t208 = t266 + t281;
t260 = -t208 + t281;
t156 = -pkin(6) * t208 + t171 - t270;
t257 = t235 * t265;
t157 = -pkin(6) * t257 + t163;
t259 = -qJD(4) * t156 - t157;
t217 = t230 + t262;
t256 = -t217 + t262;
t150 = t156 * t239 - t158 * t237;
t151 = t156 * t237 + t158 * t239;
t205 = t236 * t218;
t170 = -pkin(6) * t292 + t205 + (-pkin(5) * t235 - pkin(3)) * t240;
t175 = -pkin(6) * t235 * t238 + t182;
t254 = t170 * t239 - t175 * t237;
t253 = t170 * t237 + t175 * t239;
t252 = t206 * t237 - t208 * t239;
t251 = pkin(3) * t238 - pkin(6) * t291;
t221 = t296 * t236;
t250 = qJD(1) * t251 + qJD(3) * t235 + qJD(4) * t221 + t177;
t198 = t235 * t213;
t220 = t296 * t235;
t248 = -pkin(5) * t292 - pkin(6) * t293;
t249 = -qJD(1) * t248 + qJD(3) * t236 - qJD(4) * t220 - t198;
t245 = t251 * qJD(2);
t244 = -qJD(4) * t208 - t257;
t243 = qJD(2) * t247;
t154 = qJD(1) * t243 - qJD(4) * t252;
t229 = -pkin(3) * t236 - pkin(2);
t227 = -qJD(4) + t282;
t226 = pkin(5) * t265;
t214 = t268 * t238;
t202 = t268 * t279;
t201 = t235 * t270 + t231;
t195 = t239 * t206;
t193 = pkin(3) * t257 + t226;
t191 = t210 * t238;
t190 = t211 * t238;
t186 = t235 * t194;
t181 = t205 - t271;
t178 = -pkin(5) * t266 + t198;
t176 = pkin(3) * t206 + t217;
t174 = -t236 * t269 + t186;
t166 = t208 * t237 + t195;
t165 = qJD(2) * t248 + t186;
t161 = t276 * t292 - t277 * t294 + t243;
t160 = -qJD(2) * t246 - t211 * t277;
t159 = t245 + t173;
t155 = qJD(1) * t245 + t162;
t153 = t237 * t244 + t285;
t152 = t239 * t155;
t1 = [0.2e1 * t265 * t298 + t297 * t299 + MDP(6) * t289 - MDP(7) * t290 + (-pkin(5) * t289 + t238 * t263) * MDP(9) + (pkin(5) * t290 + t240 * t263) * MDP(10) + ((-qJD(1) * t173 - t162) * t240 + ((pkin(5) * t206 + t217 * t235) * t240 + (t171 + (t181 + 0.2e1 * t271) * qJD(1)) * t238) * qJD(2)) * MDP(11) + ((qJD(1) * t174 + t163) * t240 + ((pkin(5) * t208 + t217 * t236) * t240 + (-t172 + (-t182 + 0.2e1 * t225) * qJD(1)) * t238) * qJD(2)) * MDP(12) + (-t173 * t208 - t174 * t206 + (-t162 * t236 - t163 * t235) * t238 + (-t171 * t236 - t172 * t235 + (-t181 * t236 - t182 * t235) * qJD(1)) * t279) * MDP(13) + (t162 * t181 + t163 * t182 + t171 * t173 + t172 * t174 + (t217 + t230) * pkin(5) * t279) * MDP(14) + (-t153 * t191 - t160 * t252) * MDP(15) + (-t153 * t190 + t154 * t191 - t160 * t166 + t161 * t252) * MDP(16) + (-t153 * t240 - t160 * t227 + (-qJD(1) * t191 - t252) * t280) * MDP(17) + (t154 * t240 + t161 * t227 + (-qJD(1) * t190 - t166) * t280) * MDP(18) + (-t227 - t282) * t264 + (-(t159 * t239 - t165 * t237) * t227 - (-t157 * t237 + t152) * t240 + t202 * t166 + t214 * t154 + t193 * t190 + t176 * t161 + (t151 * t240 + t227 * t253) * qJD(4) + (qJD(1) * t254 + t150) * t280) * MDP(20) + ((t159 * t237 + t165 * t239) * t227 + (t155 * t237 + t157 * t239) * t240 - t202 * t252 + t214 * t153 - t193 * t191 + t176 * t160 + (t150 * t240 + t227 * t254) * qJD(4) + (-qJD(1) * t253 - t151) * t280) * MDP(21); -t288 * t298 + t242 * t297 + ((-qJ(3) * t281 - t171) * t238 + (-pkin(5) * t261 + t235 * t256 + t177) * t240) * t274 + ((-qJ(3) * t275 + t172) * t238 + (pkin(5) * t260 + t236 * t256 - t178) * t240) * t273 + (t177 * t208 + t178 * t206 + (-qJD(3) * t206 + t171 * t282 + t163) * t236 + (qJD(3) * t208 + t172 * t282 - t162) * t235) * MDP(13) + (-t171 * t177 - t172 * t178 + (-t171 * t235 + t172 * t236) * qJD(3) + (-t162 * t235 + t163 * t236) * qJ(3) + (-t217 - t295) * t231) * MDP(14) + (t153 * t211 - t252 * t287) * MDP(15) + (-t153 * t210 - t154 * t211 - t166 * t287 + t252 * t286) * MDP(16) + (-t287 * t227 + (qJD(2) * t211 + t252) * t283) * MDP(17) + (t286 * t227 + (-qJD(2) * t210 + t166) * t283) * MDP(18) + t227 * MDP(19) * t283 + (t229 * t154 - t201 * t166 + t193 * t210 + (t237 * t249 + t239 * t250) * t227 + t286 * t176 + ((-t220 * t239 - t221 * t237) * qJD(2) - t150) * t283) * MDP(20) + (t229 * t153 + t201 * t252 + t193 * t211 + (-t237 * t250 + t239 * t249) * t227 + t287 * t176 + (-(-t220 * t237 + t221 * t239) * qJD(2) + t151) * t283) * MDP(21) + (MDP(9) * t238 * t242 + MDP(10) * t288) * pkin(1); (-t206 ^ 2 - t208 ^ 2) * MDP(13) + (t171 * t208 + t172 * t206 + t226) * MDP(14) + (t227 * t252 + t154) * MDP(20) + (t195 * t227 + (-t257 + (-qJD(4) + t227) * t208) * t237 + t285) * MDP(21) + (t260 * t274 + t261 * t273) * t240; -t166 ^ 2 * MDP(16) + (-t166 * t227 + t285) * MDP(17) + qJD(1) * t264 + (-t151 * t227 + t152) * MDP(20) + (-t150 * t227 + t166 * t176) * MDP(21) - (MDP(15) * t166 - MDP(16) * t252 - MDP(18) * t227 - MDP(20) * t176) * t252 + (MDP(18) * t244 - MDP(20) * t278 + MDP(21) * t259) * t239 + (t244 * MDP(17) + (qJD(4) * t206 - t258) * MDP(18) + t259 * MDP(20) + (-t155 + t278) * MDP(21)) * t237;];
tauc = t1;
