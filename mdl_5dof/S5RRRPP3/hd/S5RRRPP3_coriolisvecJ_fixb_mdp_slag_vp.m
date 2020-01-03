% Calculate minimal parameter regressor of Coriolis joint torque vector for
% S5RRRPP3
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [7x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d2,d3]';
% MDP [21x1]
%   Minimal dynamic parameter vector (fixed base model)
%   see S5RRRPP3_convert_par2_MPV_fixb.m
% 
% Output:
% tauc [5x1]
%   joint torques required to compensate Coriolis and centrifugal load

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 20:54
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauc = S5RRRPP3_coriolisvecJ_fixb_mdp_slag_vp(qJ, qJD, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(7,1),zeros(21,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRRPP3_coriolisvecJ_fixb_mdp_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RRRPP3_coriolisvecJ_fixb_mdp_slag_vp: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [7 1]), ...
  'S5RRRPP3_coriolisvecJ_fixb_mdp_slag_vp: pkin has to be [7x1] (double)');
assert(isreal(MDP) && all(size(MDP) == [21 1]), ...
  'S5RRRPP3_coriolisvecJ_fixb_mdp_slag_vp: MDP has to be [21x1] (double)'); 

%% Symbolic Calculation
% From coriolisvec_joint_fixb_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 20:53:43
% EndTime: 2019-12-31 20:53:47
% DurationCPUTime: 1.15s
% Computational Cost: add. (1051->242), mult. (1751->289), div. (0->0), fcn. (767->4), ass. (0->116)
t225 = sin(qJ(3));
t222 = t225 ^ 2;
t227 = cos(qJ(3));
t223 = t227 ^ 2;
t283 = MDP(8) * (t222 - t223);
t220 = qJD(1) + qJD(2);
t226 = sin(qJ(2));
t277 = pkin(1) * qJD(1);
t249 = t226 * t277;
t193 = pkin(7) * t220 + t249;
t184 = t227 * t193;
t267 = t220 * t227;
t171 = pkin(4) * t267 + t184;
t259 = qJD(3) * qJ(4);
t239 = -qJD(5) - t259;
t161 = -t239 + t171;
t224 = -pkin(3) - qJ(5);
t282 = qJD(3) * t224;
t183 = t225 * t193;
t281 = -qJD(4) - t183;
t280 = MDP(14) + MDP(18);
t228 = cos(qJ(2));
t279 = pkin(1) * t228;
t278 = pkin(2) * t220;
t276 = pkin(1) * qJD(2);
t275 = pkin(7) * MDP(17);
t274 = qJD(3) * pkin(3);
t242 = -qJ(4) * t225 - pkin(2);
t182 = t224 * t227 + t242;
t175 = t182 - t279;
t273 = t175 * t220;
t272 = t182 * t220;
t195 = -pkin(3) * t227 + t242;
t271 = t195 * t220;
t211 = pkin(1) * t226 + pkin(7);
t229 = qJD(3) ^ 2;
t270 = t211 * t229;
t269 = t220 * MDP(5);
t268 = t220 * t225;
t266 = t225 * t227;
t246 = qJD(1) * t276;
t238 = t228 * t246;
t198 = t225 * t238;
t257 = qJD(3) * t227;
t164 = t193 * t257 + t198;
t240 = t220 * t249;
t248 = t228 * t277;
t258 = qJD(3) * t225;
t265 = t227 * t240 + t248 * t258;
t199 = t227 * t238;
t221 = qJD(3) * qJD(4);
t264 = t199 + 0.2e1 * t221;
t263 = t199 + t221;
t208 = t226 * t246;
t213 = pkin(3) * t258;
t262 = t220 * t213 + t208;
t260 = t222 + t223;
t256 = qJD(4) * t225;
t159 = t193 * t258 - t263;
t255 = t159 * MDP(14);
t254 = t211 * MDP(17);
t253 = t220 * MDP(18);
t252 = t229 * MDP(10);
t170 = -pkin(4) * t268 - t183;
t251 = qJD(4) - t170;
t250 = MDP(16) + MDP(19);
t247 = t228 * t276;
t215 = t226 * t276;
t219 = t220 ^ 2;
t245 = t219 * t266;
t243 = t220 * t257;
t244 = pkin(4) * t243 + t164;
t236 = -qJ(4) * t227 + qJ(5) * t225;
t145 = (qJD(3) * t236 - qJD(5) * t227 - t256) * t220 + t262;
t162 = qJ(5) * t258 + t227 * t239 + t213 - t256;
t157 = t215 + t162;
t241 = -t157 * t220 - t145;
t237 = (-pkin(4) * t220 - t193) * t225;
t233 = -qJ(4) * t257 - t256;
t177 = t213 + t233;
t172 = t177 + t215;
t235 = t172 * t220 + t270;
t234 = (-pkin(2) - t279) * t220 - t247;
t232 = t215 * t220 + t270;
t163 = -t248 + t271;
t185 = t195 - t279;
t231 = -t185 * t220 - t163 + t247;
t150 = qJD(3) * t237 + t263;
t151 = -qJD(3) * qJD(5) + t244;
t158 = t251 + t282;
t174 = -t274 - t281;
t176 = -t184 - t259;
t194 = -t248 - t278;
t230 = -0.2e1 * t220 * qJD(3) * t283 + 0.2e1 * t225 * MDP(7) * t243 + t229 * t227 * MDP(9) + (t194 * t257 + t225 * t208) * MDP(13) + (t150 * t227 + t151 * t225 + t158 * t257) * MDP(18) + (t164 * t225 + t174 * t257 + t176 * t258) * MDP(14) - t208 * MDP(5);
t218 = t227 * pkin(4);
t217 = t225 * pkin(4);
t214 = pkin(4) * t257;
t207 = pkin(3) * t268;
t202 = pkin(7) * t227 + t218;
t201 = pkin(7) * t225 + t217;
t190 = pkin(7) * t257 + t214;
t189 = (-pkin(4) - pkin(7)) * t258;
t187 = t211 * t227 + t218;
t186 = t211 * t225 + t217;
t180 = t194 * t258;
t178 = -qJ(4) * t267 + t207;
t167 = t220 * t236 + t207;
t166 = t211 * t257 + t225 * t247 + t214;
t165 = t227 * t247 + (-pkin(4) - t211) * t258;
t156 = t163 * t268;
t155 = t220 * t233 + t262;
t154 = -t248 + t272;
t152 = t155 * t227;
t149 = t154 * t258;
t148 = t154 * t267;
t1 = [t180 * MDP(12) + t152 * MDP(15) + (t155 * t185 + t163 * t172) * MDP(17) + t149 * MDP(20) + (t145 * t175 + t150 * t187 + t151 * t186 + t154 * t157 + t158 * t166 + t161 * t165) * MDP(21) + (-t226 * t269 + (-qJD(1) * MDP(6) + (MDP(14) * t260 - MDP(6)) * t220) * t228) * t276 + ((-t232 - t208) * MDP(12) - t255 + t235 * MDP(15) + (-t159 * t211 - t176 * t247) * MDP(17) + t165 * t253 + t241 * MDP(20)) * t227 + (-t252 + t232 * MDP(13) + (-t155 - t235) * MDP(16) + (t164 * t211 + t174 * t247) * MDP(17) + t166 * t253 + t241 * MDP(19)) * t225 + (t165 * MDP(19) - t166 * MDP(20) + (t234 * MDP(13) + t231 * MDP(16) + t174 * t254 + t186 * t253 + (-t154 - t273) * MDP(19)) * t227 + (t234 * MDP(12) + t231 * MDP(15) + t176 * t254 + (-t187 * t220 - t161) * MDP(18) + MDP(20) * t273) * t225) * qJD(3) + t230; (t145 * t182 + t150 * t202 + t151 * t201 + t154 * t162 + t158 * t190 + t161 * t189) * MDP(21) + t230 + (t155 * t195 + t163 * t177) * MDP(17) + (-MDP(13) + t250) * t225 * t240 + (t189 * MDP(19) - t190 * MDP(20) + (t176 * t275 - t163 * MDP(15) - t161 * MDP(18) + (-pkin(2) * MDP(12) - t195 * MDP(15) - t202 * MDP(18) + t182 * MDP(20)) * t220) * t225 + ((t248 - t278) * MDP(13) + (-t163 - t248 - t271) * MDP(16) + t174 * t275 + t201 * t253 + (-t154 - t248 - t272) * MDP(19)) * t227) * qJD(3) + ((-qJD(2) * t227 * MDP(12) - MDP(17) * t163 - t154 * MDP(21) + t269) * t226 + (-qJD(2) * MDP(6) + (-t174 * t225 + t176 * t227) * MDP(17) + (-t158 * t225 - t161 * t227) * MDP(21) + (-t260 * t280 + MDP(6)) * t220) * t228) * t277 + (-t252 + (-t177 * t220 - t155) * MDP(16) + t190 * t253 + (-t162 * t220 - t145) * MDP(19) + (t164 * MDP(17) + (MDP(13) - MDP(16)) * t229) * pkin(7)) * t225 + (-t255 - t145 * MDP(20) + (t177 * MDP(15) + t189 * MDP(18) - t162 * MDP(20)) * t220 + (-t159 * MDP(17) + (-MDP(12) + MDP(15)) * t229) * pkin(7)) * t227 + (t180 + t265) * MDP(12) + (t152 - t265) * MDP(15) + (t149 + t265) * MDP(20); -t198 * MDP(12) - t199 * MDP(13) + (t156 + t198) * MDP(15) + t264 * MDP(16) + (-pkin(3) * t164 - qJ(4) * t159 - t163 * t178 - t174 * t184 + t176 * t281) * MDP(17) + (t148 + t264) * MDP(19) - t244 * MDP(20) + (qJ(4) * t150 + t151 * t224 - t154 * t167 + t251 * t161 + (-qJD(5) - t171) * t158) * MDP(21) + t280 * qJD(4) * t267 + ((-t170 - t183) * MDP(19) + (0.2e1 * qJD(5) + t171) * MDP(20)) * qJD(3) + (-MDP(7) * t266 + t283) * t219 + ((-t194 * MDP(13) + (-t174 - t274) * MDP(14) - t178 * MDP(15) + t163 * MDP(16) + (-t158 - t170 + t282) * MDP(18) + t167 * MDP(20)) * t227 + (-t194 * MDP(12) + (-t176 - t259) * MDP(14) + t178 * MDP(16) + (-pkin(4) * qJD(3) + t167) * MDP(19) - t154 * MDP(20)) * t225) * t220; (qJD(3) * t176 + t156 + t164) * MDP(17) + (t154 * t268 + (-qJD(5) - t161) * qJD(3) + t244) * MDP(21) + (MDP(15) - MDP(20)) * t245 + t250 * (-t219 * t222 - t229); -MDP(19) * t245 + (-t219 * t223 - t229) * MDP(20) + (t148 + t263 + (t158 + t237) * qJD(3)) * MDP(21);];
tauc = t1;
