% Calculate Coriolis joint torque vector for
% S5RPRPR7
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d3,d5,theta2,theta4]';
% MDP [22x1]
%   Minimal dynamic parameter vector (fixed base model)
%   see S5RPRPR7_convert_par2_MPV_fixb.m
% 
% Output:
% tauc [5x1]
%   joint torques required to compensate Coriolis and centrifugal load

% Quelle: HybrDyn-Toolbox
% Datum: 2021-01-15 12:07
% Revision: d12c3222fdeb2c5f3b3c8fa5751e113be2fc3aae (2021-01-15)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauc = S5RPRPR7_coriolisvecJ_fixb_mdp_slag_vp(qJ, qJD, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(9,1),zeros(22,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPRPR7_coriolisvecJ_fixb_mdp_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RPRPR7_coriolisvecJ_fixb_mdp_slag_vp: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5RPRPR7_coriolisvecJ_fixb_mdp_slag_vp: pkin has to be [9x1] (double)');
assert(isreal(MDP) && all(size(MDP) == [22 1]), ...
  'S5RPRPR7_coriolisvecJ_fixb_mdp_slag_vp: MDP has to be [22x1] (double)'); 

%% Symbolic Calculation
% From coriolisvec_joint_fixb_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2021-01-15 12:06:32
% EndTime: 2021-01-15 12:06:37
% DurationCPUTime: 1.60s
% Computational Cost: add. (1319->228), mult. (3207->323), div. (0->0), fcn. (2170->8), ass. (0->110)
t218 = cos(qJ(3));
t283 = MDP(5) * t218;
t216 = sin(qJ(3));
t282 = (t216 ^ 2 - t218 ^ 2) * MDP(6);
t204 = sin(pkin(8)) * pkin(1) + pkin(6);
t197 = t204 * qJD(1);
t238 = qJ(4) * qJD(1) + t197;
t177 = t218 * qJD(2) - t238 * t216;
t253 = t216 * qJD(2);
t178 = t238 * t218 + t253;
t250 = qJD(1) * qJD(4);
t281 = -t178 * qJD(3) - t216 * t250;
t168 = t177 * qJD(3) + t218 * t250;
t212 = sin(pkin(9));
t276 = cos(pkin(9));
t144 = t168 * t212 - t276 * t281;
t242 = t276 * t218;
t200 = qJD(1) * t242;
t257 = qJD(1) * t216;
t187 = t212 * t257 - t200;
t185 = qJD(5) + t187;
t203 = pkin(3) * t212 + pkin(7);
t249 = pkin(3) * t257;
t195 = t212 * t218 + t276 * t216;
t259 = qJD(1) * t195;
t280 = (pkin(4) * t259 + pkin(7) * t187 + qJD(5) * t203 + t249) * t185 + t144;
t243 = t276 * t168;
t145 = t281 * t212 + t243;
t277 = qJD(3) * pkin(3);
t172 = t177 + t277;
t266 = t212 * t178;
t149 = t276 * t172 - t266;
t147 = -qJD(3) * pkin(4) - t149;
t263 = qJ(4) + t204;
t240 = qJD(3) * t263;
t179 = qJD(4) * t218 - t216 * t240;
t223 = -qJD(4) * t216 - t218 * t240;
t154 = t276 * t179 + t212 * t223;
t206 = -cos(pkin(8)) * pkin(1) - pkin(2);
t196 = -pkin(3) * t218 + t206;
t258 = qJD(1) * t196;
t186 = qJD(4) + t258;
t155 = pkin(4) * t187 - pkin(7) * t259 + t186;
t225 = -t212 * t216 + t242;
t163 = -pkin(4) * t225 - pkin(7) * t195 + t196;
t192 = t225 * qJD(3);
t193 = t263 * t218;
t241 = t263 * t216;
t166 = t276 * t193 - t212 * t241;
t189 = t195 * qJD(3);
t182 = qJD(1) * t189;
t233 = t144 * t195 - t166 * t182;
t279 = t147 * t192 - (qJD(5) * t163 + t154) * t185 + (qJD(5) * t155 + t145) * t225 + t233;
t278 = pkin(3) * t216;
t215 = sin(qJ(5));
t254 = qJD(5) * t215;
t251 = qJD(1) * qJD(3);
t245 = t216 * t251;
t199 = t212 * t245;
t183 = qJD(3) * t200 - t199;
t217 = cos(qJ(5));
t252 = t217 * qJD(3);
t261 = qJD(5) * t252 + t217 * t183;
t157 = -t254 * t259 + t261;
t275 = t157 * t215;
t274 = t163 * t182;
t267 = t259 * t215;
t174 = -t252 + t267;
t273 = t174 * t185;
t272 = t174 * t259;
t176 = qJD(3) * t215 + t217 * t259;
t271 = t176 * t185;
t270 = t176 * t259;
t269 = t182 * t215;
t268 = t183 * t215;
t219 = qJD(3) ^ 2;
t265 = t216 * t219;
t180 = t217 * t182;
t264 = t218 * t219;
t262 = -t157 * t225 + t176 * t189;
t170 = t276 * t178;
t150 = t212 * t172 + t170;
t198 = qJD(1) * t206;
t256 = qJD(1) * t218;
t255 = qJD(5) * t195;
t248 = t216 * t277;
t247 = t195 * t269;
t246 = t195 * t180;
t239 = t185 * t217;
t237 = 0.2e1 * t259;
t148 = qJD(3) * pkin(7) + t150;
t143 = t148 * t217 + t155 * t215;
t232 = t148 * t215 - t155 * t217;
t158 = t176 * qJD(5) + t268;
t231 = t158 * t225 - t174 * t189;
t229 = 0.2e1 * qJD(3) * t198;
t228 = t180 + (-t187 * t215 - t254) * t185;
t227 = -t192 * t215 - t217 * t255;
t226 = -t192 * t217 + t195 * t254;
t152 = t276 * t177 - t266;
t222 = -t203 * t182 + (t147 + t152) * t185;
t205 = -t276 * pkin(3) - pkin(4);
t201 = pkin(3) * t245;
t165 = t193 * t212 + t276 * t241;
t161 = pkin(4) * t189 - pkin(7) * t192 + t248;
t159 = pkin(4) * t182 - pkin(7) * t183 + t201;
t156 = t217 * t159;
t153 = t179 * t212 - t276 * t223;
t151 = t177 * t212 + t170;
t1 = [0.2e1 * t245 * t283 - 0.2e1 * t251 * t282 + MDP(7) * t264 - MDP(8) * t265 + (-t204 * t264 + t216 * t229) * MDP(10) + (t204 * t265 + t218 * t229) * MDP(11) + (t182 * t196 + t186 * t189 + (-t153 + (-qJD(1) * t225 + t187) * t278) * qJD(3)) * MDP(12) + (t183 * t196 + t186 * t192 + (t237 * t278 - t154) * qJD(3)) * MDP(13) + (t145 * t225 - t149 * t192 - t150 * t189 + t153 * t259 - t154 * t187 + t165 * t183 + t233) * MDP(14) + (t144 * t165 + t145 * t166 - t149 * t153 + t150 * t154 + (t186 + t258) * t248) * MDP(15) + (t157 * t195 * t217 - t226 * t176) * MDP(16) + ((-t174 * t217 - t176 * t215) * t192 + (-t275 - t158 * t217 + (t174 * t215 - t176 * t217) * qJD(5)) * t195) * MDP(17) + (-t226 * t185 + t246 + t262) * MDP(18) + (t227 * t185 + t231 - t247) * MDP(19) + (-t182 * t225 + t185 * t189) * MDP(20) + (-t232 * t189 + t153 * t174 - t156 * t225 + t165 * t158 + (t161 * t185 + t274 + (t147 * t195 + t148 * t225 - t166 * t185) * qJD(5)) * t217 + t279 * t215) * MDP(21) + (-t143 * t189 + t153 * t176 + t165 * t157 + (-(-qJD(5) * t166 + t161) * t185 - t274 + (-qJD(5) * t148 + t159) * t225 - t147 * t255) * t215 + t279 * t217) * MDP(22); (-t182 * t195 - t183 * t225 - t187 * t192 + t189 * t259) * MDP(14) + (-t144 * t225 + t145 * t195 - t149 * t189 + t150 * t192) * MDP(15) + (-t231 - t247) * MDP(21) + (-t246 + t262) * MDP(22) + (-MDP(10) * t216 - MDP(11) * t218) * t219 + (t227 * MDP(21) + t226 * MDP(22)) * t185 + (-MDP(12) * t189 - MDP(13) * t192) * qJD(3); (qJD(3) * t151 - t186 * t259 - t187 * t249 - t144) * MDP(12) + (-t243 + t186 * t187 + (-pkin(3) * t259 + qJD(4) * t212) * t257 + (-t212 * (-qJ(4) * t256 - t197 * t218 - t253) + t152) * qJD(3)) * MDP(13) + ((t150 - t151) * t259 + (-t149 + t152) * t187 + (-t182 * t212 - t276 * t183) * pkin(3)) * MDP(14) + (t149 * t151 - t150 * t152 + (-t276 * t144 + t145 * t212 - t186 * t257) * pkin(3)) * MDP(15) + (t176 * t239 + t275) * MDP(16) + ((t157 - t273) * t217 + (-t158 - t271) * t215) * MDP(17) + (t185 * t239 + t269 - t270) * MDP(18) + (t228 + t272) * MDP(19) - t185 * t259 * MDP(20) + (-t151 * t174 + t205 * t158 + t222 * t215 - t280 * t217 + t232 * t259) * MDP(21) + (t143 * t259 - t151 * t176 + t205 * t157 + t280 * t215 + t222 * t217) * MDP(22) + (-t216 * t283 + t282) * qJD(1) ^ 2 + (-MDP(10) * t257 - MDP(11) * t256) * t198; t237 * qJD(3) * MDP(12) + (-t199 + (t200 - t187) * qJD(3)) * MDP(13) + (-t187 ^ 2 - t259 ^ 2) * MDP(14) + (t149 * t259 + t150 * t187 + t201) * MDP(15) + (t228 - t272) * MDP(21) + (-t185 ^ 2 * t217 - t269 - t270) * MDP(22); t176 * t174 * MDP(16) + (-t174 ^ 2 + t176 ^ 2) * MDP(17) + (t261 + t273) * MDP(18) + (-t268 + t271) * MDP(19) + t182 * MDP(20) + (t143 * t185 - t145 * t215 - t147 * t176 + t156) * MDP(21) + (-t145 * t217 + t147 * t174 - t159 * t215 - t185 * t232) * MDP(22) + (-MDP(18) * t267 - t176 * MDP(19) - t143 * MDP(21) + t232 * MDP(22)) * qJD(5);];
tauc = t1;
