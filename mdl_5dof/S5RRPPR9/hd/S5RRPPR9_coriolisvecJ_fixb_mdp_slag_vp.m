% Calculate minimal parameter regressor of Coriolis joint torque vector for
% S5RRPPR9
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [7x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d2,d5]';
% MDP [25x1]
%   Minimal dynamic parameter vector (fixed base model)
%   see S5RRPPR9_convert_par2_MPV_fixb.m
% 
% Output:
% tauc [5x1]
%   joint torques required to compensate Coriolis and centrifugal load

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 19:42
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauc = S5RRPPR9_coriolisvecJ_fixb_mdp_slag_vp(qJ, qJD, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(7,1),zeros(25,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRPPR9_coriolisvecJ_fixb_mdp_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RRPPR9_coriolisvecJ_fixb_mdp_slag_vp: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [7 1]), ...
  'S5RRPPR9_coriolisvecJ_fixb_mdp_slag_vp: pkin has to be [7x1] (double)');
assert(isreal(MDP) && all(size(MDP) == [25 1]), ...
  'S5RRPPR9_coriolisvecJ_fixb_mdp_slag_vp: MDP has to be [25x1] (double)'); 

%% Symbolic Calculation
% From coriolisvec_joint_fixb_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 19:41:52
% EndTime: 2019-12-31 19:41:56
% DurationCPUTime: 1.92s
% Computational Cost: add. (777->263), mult. (1780->361), div. (0->0), fcn. (881->4), ass. (0->120)
t213 = sin(qJ(2));
t262 = qJD(1) * t213;
t194 = pkin(6) * t262;
t171 = -qJ(4) * t262 + t194;
t241 = qJD(3) + t171;
t216 = -pkin(2) - pkin(3);
t245 = qJD(2) * t216;
t160 = t245 + t241;
t249 = qJD(1) * qJD(2);
t283 = -0.2e1 * t249;
t280 = pkin(6) - qJ(4);
t208 = t213 ^ 2;
t215 = cos(qJ(2));
t209 = t215 ^ 2;
t282 = (t208 - t209) * MDP(5);
t261 = qJD(1) * t215;
t195 = pkin(6) * t261;
t173 = -qJ(4) * t261 + t195;
t207 = qJD(2) * qJ(3);
t166 = -t173 - t207;
t167 = -qJD(1) * pkin(1) - pkin(2) * t261 - qJ(3) * t262;
t157 = pkin(3) * t261 + qJD(4) - t167;
t230 = pkin(4) * t213 + pkin(7) * t215;
t148 = qJD(1) * t230 + t157;
t176 = -t215 * pkin(2) - t213 * qJ(3) - pkin(1);
t168 = t215 * pkin(3) - t176;
t154 = t230 + t168;
t239 = t215 * t249;
t186 = pkin(6) * t239;
t257 = qJD(4) * t213;
t258 = qJD(2) * t215;
t158 = t186 + (-qJ(4) * t258 - t257) * qJD(1);
t161 = qJD(2) * pkin(4) - t166;
t164 = t258 * t280 - t257;
t187 = qJD(5) + t262;
t281 = -(qJD(5) * t154 + t164) * t187 + (qJD(2) * t161 - qJD(5) * t148 - t158) * t213;
t279 = qJD(2) * pkin(2);
t240 = t213 * t249;
t184 = qJ(4) * t240;
t205 = qJD(2) * qJD(3);
t259 = qJD(2) * t213;
t224 = pkin(6) * t259 + qJD(4) * t215;
t151 = qJD(1) * t224 - t184 - t205;
t212 = sin(qJ(5));
t278 = t151 * t212;
t214 = cos(qJ(5));
t277 = t151 * t214;
t252 = t214 * qJD(2);
t254 = qJD(5) * t212;
t242 = t215 * t254;
t266 = qJD(1) * t242 + t214 * t240;
t152 = -qJD(5) * t252 + t266;
t276 = t152 * t212;
t169 = t212 * t261 - t252;
t275 = t169 * t187;
t274 = t169 * t215;
t260 = qJD(2) * t212;
t170 = t214 * t261 + t260;
t273 = t170 * t187;
t272 = t170 * t215;
t271 = t187 * t213;
t270 = t187 * t214;
t217 = qJD(2) ^ 2;
t269 = t213 * t217;
t268 = t215 * t217;
t218 = qJD(1) ^ 2;
t267 = t215 * t218;
t198 = t213 * qJD(3);
t265 = qJ(3) * t239 + qJD(1) * t198;
t264 = qJ(3) * t258 + t198;
t206 = -pkin(7) + t216;
t155 = qJD(2) * t206 + t241;
t256 = qJD(5) * t155;
t255 = qJD(5) * t187;
t253 = t161 * qJD(5);
t250 = -qJD(4) - t157;
t248 = t212 * t271;
t247 = t213 * t270;
t246 = t213 * t267;
t244 = t187 * t254;
t243 = t214 * t255;
t238 = MDP(23) * t261;
t237 = pkin(1) * t283;
t236 = qJD(3) - t279;
t231 = t213 * t245;
t150 = qJD(1) * t231 + t265;
t156 = t231 + t264;
t235 = qJD(1) * t156 + t150;
t234 = qJD(1) * t168 + t157;
t233 = qJD(1) * t176 + t167;
t143 = t148 * t214 - t155 * t212;
t144 = t148 * t212 + t155 * t214;
t178 = t280 * t213;
t222 = pkin(4) * t215 + t206 * t213;
t220 = t222 * qJD(2);
t229 = (-qJD(5) * t178 + t220 + t264) * t187;
t228 = qJD(1) * t209 - t271;
t227 = -t214 * MDP(24) + t212 * MDP(25);
t226 = -t212 * MDP(24) - t214 * MDP(25);
t159 = pkin(2) * t240 - t265;
t165 = pkin(2) * t259 - t264;
t225 = -pkin(6) * t217 - qJD(1) * t165 - t159;
t223 = -t161 * t213 - t206 * t258;
t174 = -pkin(6) * t240 + t205;
t175 = t194 + t236;
t177 = t195 + t207;
t219 = t174 * t215 + (t175 * t215 + (-t177 + t195) * t213) * qJD(2);
t211 = qJ(3) + pkin(4);
t202 = 0.2e1 * t205;
t191 = qJ(3) * t261;
t180 = t212 * t240;
t179 = t280 * t215;
t172 = pkin(2) * t262 - t191;
t163 = -qJ(4) * t259 + t224;
t162 = t216 * t262 + t191;
t153 = qJD(5) * t170 - t180;
t149 = qJD(1) * t222 + t191;
t146 = qJD(1) * t220 + t265;
t145 = t214 * t146;
t1 = [0.2e1 * t213 * MDP(4) * t239 + t282 * t283 + MDP(6) * t268 - MDP(7) * t269 + (-pkin(6) * t268 + t213 * t237) * MDP(9) + (pkin(6) * t269 + t215 * t237) * MDP(10) + (t215 * t225 + t233 * t259) * MDP(11) + t219 * MDP(12) + (t213 * t225 - t233 * t258) * MDP(13) + (pkin(6) * t219 + t159 * t176 + t165 * t167) * MDP(14) + (t235 * t213 + (t215 * t234 - t163) * qJD(2)) * MDP(15) + (-t235 * t215 + (t213 * t234 + t164) * qJD(2)) * MDP(16) + (t151 * t215 - t158 * t213 + (-t160 * t215 - t166 * t213) * qJD(2) + (t163 * t215 - t164 * t213 + (-t178 * t215 + t179 * t213) * qJD(2)) * qJD(1)) * MDP(17) + (t150 * t168 - t151 * t179 + t156 * t157 + t158 * t178 + t160 * t164 + t163 * t166) * MDP(18) + (-t152 * t214 * t215 + (-t213 * t252 - t242) * t170) * MDP(19) + ((t169 * t214 + t170 * t212) * t259 + (t276 - t153 * t214 + (t169 * t212 - t170 * t214) * qJD(5)) * t215) * MDP(20) + (t187 * t242 + t152 * t213 + (-t214 * t228 - t272) * qJD(2)) * MDP(21) + (t215 * t243 + t153 * t213 + (t212 * t228 + t274) * qJD(2)) * MDP(22) + (t187 + t262) * MDP(23) * t258 + (t145 * t213 - t179 * t153 + t163 * t169 + (-t213 * t256 + t229) * t214 + t281 * t212 + (-t214 * t253 + t278 + ((t154 * t214 - t178 * t212) * qJD(1) + t143) * qJD(2)) * t215) * MDP(24) + (t179 * t152 + t163 * t170 + (-t229 - (t146 - t256) * t213) * t212 + t281 * t214 + (t212 * t253 + t277 + (-(t154 * t212 + t178 * t214) * qJD(1) - t144) * qJD(2)) * t215) * MDP(25); -MDP(4) * t246 + t218 * t282 + t202 * MDP(13) + (qJ(3) * t174 + qJD(3) * t177 - t167 * t172) * MDP(14) + (qJD(2) * t171 + t184 + t202) * MDP(15) + (-qJD(2) * t173 + t186) * MDP(16) + (-qJ(3) * t151 - t157 * t162 + t158 * t216 - t160 * t173 - t166 * t241) * MDP(18) + (t170 * t270 - t276) * MDP(19) + ((-t152 - t275) * t214 + (-t153 - t273) * t212) * MDP(20) - t243 * MDP(21) + t244 * MDP(22) - t187 * t238 + (-t211 * t153 - t277 - (t149 * t214 - t173 * t212) * t187 - t241 * t169 + (-t161 * t212 - t206 * t270) * qJD(5)) * MDP(24) + (t211 * t152 + t278 + (t149 * t212 + t173 * t214) * t187 - t241 * t170 + (t187 * t206 * t212 - t161 * t214) * qJD(5)) * MDP(25) + (MDP(9) * t213 * t218 + MDP(10) * t267) * pkin(1) + ((-t167 * t213 + t172 * t215) * MDP(11) + ((t177 - t207) * t213 + (-t175 + t236) * t215) * MDP(12) + (t167 * t215 + t172 * t213) * MDP(13) + (t177 * t213 + (-t175 - t279) * t215) * pkin(6) * MDP(14) + (t250 * t215 + (-pkin(6) * qJD(2) - t162) * t213) * MDP(15) + ((-qJ(4) * qJD(2) + t162) * t215 + t250 * t213) * MDP(16) + (-t247 + (t170 - t260) * t215) * MDP(21) + (t248 + (-t169 - t252) * t215) * MDP(22) + (-t143 * t215 + t212 * t223) * MDP(24) + (t144 * t215 + t214 * t223) * MDP(25)) * qJD(1); (-qJD(2) * t177 + t186) * MDP(14) + (qJD(2) * t166 + t186) * MDP(18) + (qJD(2) * t169 - t243) * MDP(24) + (qJD(2) * t170 + t244) * MDP(25) + (-MDP(11) + MDP(16)) * t246 + (MDP(13) + MDP(15)) * (-t208 * t218 - t217) + ((-MDP(18) * qJ(4) + t226) * t258 + (t167 * MDP(14) + MDP(18) * t250 + t187 * t227) * t213) * qJD(1); t265 * MDP(18) + (-t208 - t209) * MDP(17) * t218 + t226 * t255 + ((t160 * t213 - t166 * t215) * MDP(18) + (-t248 - t274) * MDP(24) + (-t247 - t272) * MDP(25) + ((0.2e1 * MDP(15) - t227) * t215 + (MDP(18) * t216 + 0.2e1 * MDP(16)) * t213) * qJD(2)) * qJD(1); t170 * t169 * MDP(19) + (-t169 ^ 2 + t170 ^ 2) * MDP(20) + (t266 - t275) * MDP(21) + (-t180 - t273) * MDP(22) + qJD(2) * t238 + (t144 * t187 - t158 * t212 + t161 * t170 + t145) * MDP(24) + (t143 * t187 - t146 * t212 - t158 * t214 - t161 * t169) * MDP(25) + (-MDP(21) * t252 + MDP(22) * t170 - MDP(24) * t144 - MDP(25) * t143) * qJD(5);];
tauc = t1;
