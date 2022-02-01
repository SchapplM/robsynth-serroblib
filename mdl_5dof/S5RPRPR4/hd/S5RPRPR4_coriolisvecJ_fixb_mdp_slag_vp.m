% Calculate Coriolis joint torque vector for
% S5RPRPR4
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
%   see S5RPRPR4_convert_par2_MPV_fixb.m
% 
% Output:
% tauc [5x1]
%   joint torques required to compensate Coriolis and centrifugal load

% Quelle: HybrDyn-Toolbox
% Datum: 2022-01-23 09:23
% Revision: fd3771346c4aea32fdeb66112c511235427c26a7 (2022-01-20)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauc = S5RPRPR4_coriolisvecJ_fixb_mdp_slag_vp(qJ, qJD, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(9,1),zeros(22,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPRPR4_coriolisvecJ_fixb_mdp_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RPRPR4_coriolisvecJ_fixb_mdp_slag_vp: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5RPRPR4_coriolisvecJ_fixb_mdp_slag_vp: pkin has to be [9x1] (double)');
assert(isreal(MDP) && all(size(MDP) == [22 1]), ...
  'S5RPRPR4_coriolisvecJ_fixb_mdp_slag_vp: MDP has to be [22x1] (double)'); 

%% Symbolic Calculation
% From coriolisvec_joint_fixb_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2022-01-23 09:23:13
% EndTime: 2022-01-23 09:23:16
% DurationCPUTime: 1.58s
% Computational Cost: add. (1158->196), mult. (2848->283), div. (0->0), fcn. (1960->8), ass. (0->95)
t243 = sin(pkin(9));
t245 = cos(pkin(9));
t250 = cos(qJ(3));
t271 = qJD(1) * t250;
t262 = t245 * t271;
t248 = sin(qJ(3));
t272 = qJD(1) * t248;
t219 = t243 * t272 - t262;
t249 = cos(qJ(5));
t210 = t249 * t219;
t228 = t243 * t250 + t245 * t248;
t221 = t228 * qJD(3);
t213 = qJD(1) * t221;
t269 = qJD(3) * t248;
t261 = qJD(1) * t269;
t232 = t243 * t261;
t214 = qJD(3) * t262 - t232;
t247 = sin(qJ(5));
t268 = qJD(5) * t247;
t275 = qJD(1) * t228;
t152 = -qJD(5) * t210 - t247 * t213 + t249 * t214 - t268 * t275;
t176 = -t247 * t275 - t210;
t255 = t219 * t247 - t249 * t275;
t253 = t255 * qJD(5) - t249 * t213 - t214 * t247;
t240 = qJD(3) + qJD(5);
t279 = t176 * t240;
t280 = t255 * t240;
t290 = t176 * t255 * MDP(16) + (-t176 ^ 2 + t255 ^ 2) * MDP(17) + (t152 - t279) * MDP(18) + (t253 - t280) * MDP(19);
t236 = sin(pkin(8)) * pkin(1) + pkin(6);
t277 = qJ(4) + t236;
t257 = t277 * qJD(1);
t200 = t250 * qJD(2) - t257 * t248;
t195 = qJD(3) * pkin(3) + t200;
t201 = qJD(2) * t248 + t257 * t250;
t278 = t245 * t201;
t166 = t243 * t195 + t278;
t282 = pkin(7) * t219;
t156 = t166 - t282;
t238 = -cos(pkin(8)) * pkin(1) - pkin(2);
t229 = -pkin(3) * t250 + t238;
t274 = qJD(1) * t229;
t218 = qJD(4) + t274;
t185 = pkin(4) * t219 + t218;
t289 = t156 * t268 - t185 * t176;
t264 = qJD(1) * qJD(4);
t183 = t200 * qJD(3) + t250 * t264;
t184 = -t201 * qJD(3) - t248 * t264;
t157 = -t183 * t243 + t245 * t184;
t150 = -pkin(7) * t214 + t157;
t158 = t245 * t183 + t243 * t184;
t151 = -pkin(7) * t213 + t158;
t288 = t249 * t150 - t247 * t151 + t185 * t255;
t286 = (t248 ^ 2 - t250 ^ 2) * MDP(6);
t265 = t250 * MDP(11);
t285 = -t248 * MDP(10) - t265;
t284 = qJD(5) - t240;
t283 = pkin(3) * t243;
t281 = pkin(7) * t275;
t191 = t243 * t201;
t168 = t245 * t200 - t191;
t258 = qJD(3) * t277;
t204 = qJD(4) * t250 - t248 * t258;
t205 = -qJD(4) * t248 - t250 * t258;
t170 = t245 * t204 + t243 * t205;
t225 = t277 * t248;
t226 = t277 * t250;
t182 = -t243 * t225 + t245 * t226;
t231 = qJD(1) * t238;
t263 = pkin(3) * t272;
t165 = t245 * t195 - t191;
t167 = -t200 * t243 - t278;
t169 = -t204 * t243 + t245 * t205;
t181 = -t245 * t225 - t226 * t243;
t155 = qJD(3) * pkin(4) + t165 - t281;
t256 = -t155 * t247 - t156 * t249;
t227 = t243 * t248 - t245 * t250;
t186 = t227 * t249 + t228 * t247;
t187 = -t227 * t247 + t228 * t249;
t251 = qJD(3) ^ 2;
t237 = pkin(3) * t245 + pkin(4);
t234 = pkin(3) * t261;
t224 = t227 * qJD(3);
t203 = pkin(3) * t269 + pkin(4) * t221;
t202 = pkin(4) * t275 + t263;
t199 = pkin(4) * t227 + t229;
t190 = pkin(4) * t213 + t234;
t172 = -pkin(7) * t227 + t182;
t171 = -pkin(7) * t228 + t181;
t164 = -pkin(7) * t221 + t170;
t163 = pkin(7) * t224 + t169;
t162 = t168 - t281;
t161 = t167 + t282;
t160 = t187 * qJD(5) + t249 * t221 - t224 * t247;
t159 = -t186 * qJD(5) - t221 * t247 - t224 * t249;
t1 = [(t213 * t229 + t218 * t221) * MDP(12) + (t214 * t229 - t218 * t224) * MDP(13) + (-t157 * t228 - t158 * t227 + t165 * t224 - t166 * t221 - t169 * t275 - t170 * t219 - t181 * t214 - t182 * t213) * MDP(14) + (t157 * t181 + t158 * t182 + t165 * t169 + t166 * t170) * MDP(15) + (t152 * t187 - t159 * t255) * MDP(16) + (-t152 * t186 + t159 * t176 + t160 * t255 + t187 * t253) * MDP(17) + (t185 * t160 - t176 * t203 + t190 * t186 - t199 * t253) * MDP(21) + (t199 * t152 + t185 * t159 + t190 * t187 - t203 * t255) * MDP(22) + (t159 * MDP(18) - t160 * MDP(19) + (t163 * t249 - t164 * t247) * MDP(21) + (-t163 * t247 - t164 * t249) * MDP(22) + ((-t171 * t247 - t172 * t249) * MDP(21) + (-t171 * t249 + t172 * t247) * MDP(22)) * qJD(5)) * t240 + (t250 * MDP(7) - t248 * MDP(8) + (-MDP(10) * t250 + MDP(11) * t248) * t236) * t251 + (t231 * t265 + t169 * MDP(12) - t170 * MDP(13) + (t238 * t265 - 0.2e1 * t286) * qJD(1) + (0.2e1 * MDP(5) * t271 + 0.2e1 * MDP(10) * t231 + ((qJD(1) * t227 + t219) * MDP(12) + 0.2e1 * t275 * MDP(13) + (t218 + t274) * MDP(15)) * pkin(3)) * t248) * qJD(3); (-t213 * t228 + t214 * t227 + t219 * t224 + t221 * t275) * MDP(14) + (-t157 * t227 + t158 * t228 - t165 * t221 - t166 * t224) * MDP(15) + t285 * t251 + (-t160 * MDP(21) - t159 * MDP(22)) * t240 + (-t221 * MDP(12) + t224 * MDP(13)) * qJD(3); (-qJD(3) * t167 - t218 * t275 - t219 * t263 + t157) * MDP(12) + (qJD(3) * t168 + t218 * t219 - t263 * t275 - t158) * MDP(13) + ((t166 + t167) * t275 + (-t165 + t168) * t219 + (-t213 * t243 - t214 * t245) * pkin(3)) * MDP(14) + (-t165 * t167 - t166 * t168 + (t157 * t245 + t158 * t243 - t218 * t272) * pkin(3)) * MDP(15) + (t202 * t176 - (t161 * t249 - t162 * t247) * t240 + ((-t237 * t247 - t249 * t283) * t240 + t256) * qJD(5) + t288) * MDP(21) + (-t249 * t151 - t247 * t150 + t202 * t255 + (t161 * t247 + t162 * t249) * t240 + (-(t237 * t249 - t247 * t283) * t240 - t155 * t249) * qJD(5) + t289) * MDP(22) + (t285 * t231 + (-t248 * t250 * MDP(5) + t286) * qJD(1)) * qJD(1) + t290; -t232 * MDP(13) + (-t219 ^ 2 - t275 ^ 2) * MDP(14) + (t165 * t275 + t166 * t219 + t234) * MDP(15) + (-t253 - t280) * MDP(21) + (t152 + t279) * MDP(22) + ((t243 * t271 + t245 * t272 + t275) * MDP(12) + (-t219 + t262) * MDP(13)) * qJD(3); (t284 * t256 + t288) * MDP(21) + ((-t156 * t240 - t150) * t247 + (-t284 * t155 - t151) * t249 + t289) * MDP(22) + t290;];
tauc = t1;
