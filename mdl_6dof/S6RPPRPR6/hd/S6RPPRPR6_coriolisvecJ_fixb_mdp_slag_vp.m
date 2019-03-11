% Calculate minimal parameter regressor of Coriolis joint torque vector for
% S6RPPRPR6
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d4,d6]';
% MDP [27x1]
%   Minimal dynamic parameter vector (fixed base model)
%   see S6RPPRPR6_convert_par2_MPV_fixb.m
% 
% Output:
% tauc [6x1]
%   joint torques required to compensate Coriolis and centrifugal load

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 01:51
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauc = S6RPPRPR6_coriolisvecJ_fixb_mdp_slag_vp(qJ, qJD, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(8,1),zeros(27,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPPRPR6_coriolisvecJ_fixb_mdp_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RPPRPR6_coriolisvecJ_fixb_mdp_slag_vp: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S6RPPRPR6_coriolisvecJ_fixb_mdp_slag_vp: pkin has to be [8x1] (double)');
assert(isreal(MDP) && all(size(MDP) == [27 1]), ...
  'S6RPPRPR6_coriolisvecJ_fixb_mdp_slag_vp: MDP has to be [27x1] (double)'); 

%% Symbolic Calculation
% From coriolisvec_joint_fixb_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 01:51:36
% EndTime: 2019-03-09 01:51:39
% DurationCPUTime: 1.72s
% Computational Cost: add. (804->260), mult. (1552->344), div. (0->0), fcn. (744->4), ass. (0->125)
t193 = qJ(2) * qJD(1) + qJD(3);
t181 = -pkin(7) * qJD(1) + t193;
t205 = cos(qJ(4));
t161 = (pkin(5) * qJD(1) - t181) * t205;
t247 = qJD(5) + t161;
t201 = pkin(1) + qJ(3);
t284 = qJD(1) * t201;
t206 = -pkin(4) - pkin(8);
t253 = qJD(4) * t206;
t202 = sin(qJ(6));
t204 = cos(qJ(6));
t255 = qJD(4) * t204;
t203 = sin(qJ(4));
t261 = qJD(1) * t203;
t169 = t202 * t261 + t255;
t283 = t169 * MDP(21);
t282 = qJ(2) * MDP(6) + MDP(5) + MDP(7);
t269 = t181 * t205;
t225 = -qJD(5) + t269;
t279 = qJD(4) * pkin(4);
t163 = -t225 - t279;
t171 = t203 * t181;
t258 = qJD(4) * qJ(5);
t165 = -t171 - t258;
t246 = qJD(1) * qJD(2);
t187 = t203 * t246;
t153 = -t187 + (-qJD(5) - t269) * qJD(4);
t256 = qJD(4) * t203;
t166 = t181 * t256;
t156 = -t205 * t246 + t166;
t218 = -t153 * t203 - t156 * t205;
t281 = (t163 * t203 - t165 * t205) * qJD(4) + t218;
t200 = -pkin(7) + qJ(2);
t280 = pkin(5) - t200;
t143 = t187 + (qJD(5) - t161) * qJD(4);
t277 = t143 * t202;
t276 = t143 * t204;
t251 = qJD(6) * t202;
t245 = qJD(1) * qJD(4);
t231 = t205 * t245;
t234 = t204 * t261;
t264 = qJD(6) * t234 + t202 * t231;
t150 = -qJD(4) * t251 + t264;
t275 = t150 * t203;
t274 = t150 * t204;
t257 = qJD(4) * t202;
t167 = -t234 + t257;
t260 = qJD(1) * t205;
t186 = qJD(6) + t260;
t273 = t167 * t186;
t272 = t167 * t203;
t271 = t169 * t186;
t270 = t169 * t203;
t268 = t186 * t205;
t267 = t186 * t206;
t207 = qJD(4) ^ 2;
t266 = t200 * t207;
t265 = t204 * t205;
t170 = pkin(4) * t260 + qJ(5) * t261;
t262 = MDP(27) * t204;
t259 = qJD(2) * t205;
t254 = qJD(4) * t205;
t252 = qJD(5) * t205;
t250 = qJD(6) * t204;
t249 = pkin(4) * t261 - qJD(2);
t182 = -qJD(2) + t284;
t248 = qJD(2) - t182;
t196 = qJD(3) * qJD(1);
t244 = -MDP(15) + MDP(18);
t242 = 0.2e1 * t196;
t241 = t202 * t268;
t240 = t186 * t265;
t208 = qJD(1) ^ 2;
t239 = t203 * t205 * t208;
t232 = t203 * t245;
t238 = pkin(4) * t231 + qJ(5) * t232 + t196;
t237 = t186 * t251;
t236 = qJD(6) * t268;
t235 = t203 * t250;
t233 = pkin(4) * t254 + qJ(5) * t256 + qJD(3);
t230 = MDP(25) * t261;
t229 = qJD(4) * t280;
t220 = (qJD(4) * pkin(8) - qJD(5)) * t205;
t142 = qJD(1) * t220 + t238;
t145 = t166 + (-pkin(5) * t256 - t259) * qJD(1);
t227 = -t142 * t202 + t204 * t145;
t226 = t186 + t260;
t224 = -qJ(5) * t205 + t201;
t160 = -pkin(5) * t261 + t171;
t177 = t202 * t232;
t223 = t186 * t250 - t177;
t155 = t224 * qJD(1) + t249;
t195 = t203 * pkin(4);
t172 = t195 + t224;
t222 = qJD(1) * t172 + qJD(2) + t155;
t221 = qJD(2) + t182 + t284;
t209 = pkin(8) * t203 + t224;
t146 = t209 * qJD(1) + t249;
t147 = t253 + t247;
t141 = t146 * t204 + t147 * t202;
t219 = t146 * t202 - t147 * t204;
t217 = -t163 * t205 - t165 * t203;
t164 = t195 + t209;
t174 = t280 * t205;
t215 = t164 * t204 + t174 * t202;
t198 = t203 ^ 2;
t214 = -qJD(1) * t198 + t268;
t213 = t186 * t202;
t212 = t165 * MDP(20) - t169 * MDP(27);
t211 = t242 - t266;
t148 = -qJD(1) * t252 + t238;
t159 = t233 - t252;
t210 = -qJD(1) * t159 - t148 + t266;
t199 = t205 ^ 2;
t194 = t199 * t208;
t179 = t204 * t231;
t173 = t280 * t203;
t162 = pkin(8) * t260 + t170;
t158 = qJD(2) * t203 - t205 * t229;
t157 = -t203 * t229 - t259;
t154 = t160 + t258;
t152 = t220 + t233;
t151 = t169 * qJD(6) - t179;
t149 = t155 * t260;
t1 = [MDP(8) * t242 + (qJD(2) * t193 + qJD(3) * t182 + (qJ(2) * qJD(2) + qJD(3) * t201) * qJD(1)) * MDP(9) + 0.2e1 * (t198 - t199) * MDP(11) * t245 - t207 * t205 * MDP(13) + (t211 * t205 - t221 * t256) * MDP(16) + ((-t198 - t199) * t246 - t281) * MDP(17) + (t210 * t205 + t222 * t256) * MDP(19) + (t217 * qJD(2) + t148 * t172 + t155 * t159 + t200 * t281) * MDP(20) + (t235 * t169 + t202 * t275) * MDP(21) + (t186 * t235 + t150 * t205 + (t214 * t202 - t270) * qJD(4)) * MDP(23) + (-t151 * t205 + (t214 * t204 + t272) * qJD(4)) * MDP(24) - t226 * MDP(25) * t256 + ((-t152 * t202 + t157 * t204) * t186 + t158 * t167 - t173 * t151 + (-t154 * t255 + t227) * t205 + (-t141 * t205 - t215 * t186) * qJD(6)) * MDP(26) + (-t173 * t150 + t158 * t169 + (-(qJD(6) * t174 + t152) * t186 - (qJD(6) * t147 + t142) * t205) * t204 + (-(-qJD(6) * t164 + t157) * t186 + (qJD(4) * t154 + qJD(6) * t146 - t145) * t205) * t202) * MDP(27) + (t221 * MDP(15) - t222 * MDP(18) + t202 * t283 + (-t167 * t202 + t169 * t204) * MDP(22)) * t254 + (-0.2e1 * MDP(10) * t231 - t207 * MDP(12) + t211 * MDP(15) + t210 * MDP(18) + (t274 - t151 * t202 + (-t167 * t204 - t169 * t202) * qJD(6)) * MDP(22) - t237 * MDP(24) + (t154 * t251 - t276 + (-(-t164 * t202 + t174 * t204) * qJD(1) + t219) * qJD(4)) * MDP(26) + (t154 * t250 + t277 + (t215 * qJD(1) + t141) * qJD(4)) * MDP(27)) * t203 + 0.2e1 * t282 * t246; -t196 * MDP(9) + t194 * MDP(17) - t238 * MDP(20) + t223 * MDP(26) - MDP(27) * t237 + (MDP(17) * t198 - t282) * t208 + (-t193 * MDP(9) + (-t217 + t252) * MDP(20) + (t240 - t272) * MDP(26) + (-t241 - t270) * MDP(27) + (0.2e1 * t244 * t205 + (0.2e1 * MDP(16) - 0.2e1 * MDP(19) - t262) * t203) * qJD(4)) * qJD(1); -t208 * MDP(8) + t218 * MDP(20) + (t151 * t203 + t202 * t236) * MDP(26) + (t204 * t236 + t275) * MDP(27) + (t248 * MDP(9) - t155 * MDP(20) + (MDP(26) * t202 + t262) * t186) * qJD(1) + ((MDP(26) * t167 - t212) * t205 + (t163 * MDP(20) + (MDP(26) * t204 - MDP(27) * t202) * t226) * t203) * qJD(4) + ((-MDP(16) + MDP(19)) * t205 + t244 * t203) * (t207 + t208); MDP(10) * t239 + (-t198 * t208 + t194) * MDP(11) + t248 * MDP(15) * t260 + (t182 * t261 - t187) * MDP(16) + ((-t165 - t258) * t205 + (-qJD(5) + t163 + t279) * t203) * qJD(1) * MDP(17) + (t149 + (t170 * t203 - t259) * qJD(1)) * MDP(18) + (0.2e1 * qJD(4) * qJD(5) + t187 + (-t155 * t203 + t170 * t205) * qJD(1)) * MDP(19) + (-pkin(4) * t156 - qJ(5) * t153 - t155 * t170 - t163 * t171 + t225 * t165) * MDP(20) + (-t169 * t213 + t274) * MDP(21) + ((-t151 - t271) * t204 + (-t150 + t273) * t202) * MDP(22) + (-t237 + (-t241 + (t169 - t255) * t203) * qJD(1)) * MDP(23) + ((-t240 - t272) * qJD(1) - t223) * MDP(24) + t186 * t230 + (qJ(5) * t151 + t277 - (t160 * t204 - t162 * t202) * t186 + t247 * t167 + (t154 * t204 - t202 * t267) * qJD(6) + (t154 * t265 + (-t204 * t253 - t219) * t203) * qJD(1)) * MDP(26) + (qJ(5) * t150 + t276 + (t160 * t202 + t162 * t204) * t186 + t247 * t169 + (-t154 * t202 - t204 * t267) * qJD(6) + (-t141 * t203 + (-t154 * t205 + t203 * t253) * t202) * qJD(1)) * MDP(27); -MDP(18) * t239 + (-t194 - t207) * MDP(19) + (t149 + t156) * MDP(20) + t177 * MDP(27) + ((-t167 - t234) * MDP(26) + t212) * qJD(4) + (-MDP(26) * t213 - t186 * t262) * t186; t167 * t283 + (-t167 ^ 2 + t169 ^ 2) * MDP(22) + (t264 + t273) * MDP(23) + (t179 + t271) * MDP(24) - qJD(4) * t230 + (t141 * t186 - t154 * t169 + t227) * MDP(26) + (-t142 * t204 - t145 * t202 + t154 * t167 - t186 * t219) * MDP(27) + (-MDP(23) * t257 - t169 * MDP(24) - t141 * MDP(26) + t219 * MDP(27)) * qJD(6);];
tauc  = t1;
