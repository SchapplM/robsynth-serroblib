% Calculate minimal parameter regressor of Coriolis joint torque vector for
% S5RPRRR8
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d3,d4,d5]';
% MDP [23x1]
%   Minimal dynamic parameter vector (fixed base model)
%   see S5RPRRR8_convert_par2_MPV_fixb.m
% 
% Output:
% tauc [5x1]
%   joint torques required to compensate Coriolis and centrifugal load

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 19:06
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauc = S5RPRRR8_coriolisvecJ_fixb_mdp_slag_vp(qJ, qJD, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(8,1),zeros(23,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPRRR8_coriolisvecJ_fixb_mdp_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RPRRR8_coriolisvecJ_fixb_mdp_slag_vp: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RPRRR8_coriolisvecJ_fixb_mdp_slag_vp: pkin has to be [8x1] (double)');
assert(isreal(MDP) && all(size(MDP) == [23 1]), ...
  'S5RPRRR8_coriolisvecJ_fixb_mdp_slag_vp: MDP has to be [23x1] (double)'); 

%% Symbolic Calculation
% From coriolisvec_joint_fixb_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 19:06:09
% EndTime: 2019-12-31 19:06:12
% DurationCPUTime: 1.27s
% Computational Cost: add. (1021->158), mult. (1641->234), div. (0->0), fcn. (944->6), ass. (0->100)
t202 = qJD(4) + qJD(5);
t211 = cos(qJ(3));
t287 = 0.2e1 * t211;
t208 = sin(qJ(3));
t213 = qJD(4) ^ 2;
t247 = qJD(1) - qJD(3);
t234 = t247 ^ 2;
t286 = t208 * (t213 + t234);
t206 = sin(qJ(5));
t207 = sin(qJ(4));
t209 = cos(qJ(5));
t210 = cos(qJ(4));
t179 = t206 * t210 + t207 * t209;
t282 = t202 * t179;
t151 = t282 * t247;
t254 = qJD(4) * t247;
t285 = -0.2e1 * t254;
t273 = pkin(7) + pkin(8);
t261 = t209 * t210;
t283 = t202 * t261;
t276 = qJD(5) - t202;
t212 = -pkin(1) - pkin(2);
t191 = qJD(1) * t212 + qJD(2);
t256 = qJD(1) * qJ(2);
t239 = qJD(3) * t256;
t248 = qJD(1) * qJD(2);
t255 = qJD(3) * t208;
t161 = t191 * t255 + t208 * t248 + t211 * t239;
t173 = t191 * t208 + t211 * t256;
t281 = -t173 * t247 - t161;
t230 = qJ(2) * t211 + t208 * t212;
t171 = qJD(2) * t208 + qJD(3) * t230;
t280 = t171 * t247 + t161;
t265 = t206 * t207;
t158 = t202 * t265 - t283;
t279 = -t158 * MDP(19) - MDP(20) * t282;
t278 = -t208 * qJ(2) + t211 * t212;
t277 = MDP(11) * (t207 ^ 2 - t210 ^ 2);
t244 = t247 * t265;
t150 = t202 * t244 - t247 * t283;
t167 = t247 * t261 - t244;
t169 = t179 * t247;
t178 = -t261 + t265;
t275 = -(t150 * t179 + t158 * t169) * MDP(17) + (t150 * t178 - t151 * t179 - t158 * t167 - t169 * t282) * MDP(18) + t277 * t285;
t274 = 0.2e1 * qJD(2);
t272 = pkin(3) * t247;
t271 = pkin(4) * t210;
t181 = -pkin(7) + t230;
t270 = pkin(8) - t181;
t269 = qJD(4) * pkin(4);
t266 = t191 * t211;
t237 = -t247 * t273 + t173;
t156 = t237 * t210;
t262 = t209 * t156;
t259 = qJD(3) * t266 + t211 * t248;
t257 = MDP(10) * t210;
t250 = t213 * MDP(12);
t249 = t213 * MDP(13);
t246 = pkin(4) * t247 * t207;
t245 = t207 * t269;
t199 = -pkin(3) - t271;
t243 = qJD(4) * t273;
t241 = t207 * t257;
t155 = t237 * t207;
t153 = -t155 + t269;
t238 = -pkin(4) * t202 - t153;
t236 = qJD(4) * t270;
t172 = -t208 * t256 + t266;
t180 = pkin(3) - t278;
t233 = -t173 + t245;
t163 = -t172 + t272;
t170 = qJD(2) * t211 + t278 * qJD(3);
t232 = -t180 * t247 - t163 - t170;
t231 = qJD(4) * t237;
t154 = -t245 * t247 + t161;
t157 = -t199 * t247 - t172;
t227 = t154 * t179 - t157 * t158;
t226 = t154 * t178 + t157 * t282;
t225 = pkin(7) * t213 - t281;
t224 = qJD(4) * (t163 + t172 + t272);
t223 = t254 * t287;
t160 = -t208 * t239 + t259;
t145 = t160 * t210 - t207 * t231;
t146 = -t160 * t207 - t210 * t231;
t222 = -t206 * t145 + t209 * t146 + t157 * t169;
t221 = t181 * t213 - t280;
t219 = t202 * t178;
t218 = -t169 * t167 * MDP(17) + (t167 * t202 + t150) * MDP(19) + (-t169 * t202 + t151) * MDP(20) + (-t167 ^ 2 + t169 ^ 2) * MDP(18);
t217 = t157 * t167 + (t156 * t276 - t146) * t206;
t190 = t273 * t210;
t189 = t273 * t207;
t183 = t210 * t243;
t182 = t207 * t243;
t175 = t180 + t271;
t166 = t270 * t210;
t165 = t270 * t207;
t162 = t171 - t245;
t149 = -t170 * t207 + t210 * t236;
t148 = t170 * t210 + t207 * t236;
t1 = [t280 * MDP(8) + (t170 * t247 + t259) * MDP(9) + (-t175 * t151 + t162 * t167 - t226) * MDP(22) + (t175 * t150 - t162 * t169 - t227) * MDP(23) + ((-t148 * t206 + t149 * t209) * MDP(22) + (-t148 * t209 - t149 * t206) * MDP(23) + ((-t165 * t206 + t166 * t209) * MDP(22) + (-t165 * t209 - t166 * t206) * MDP(23)) * qJD(5) - t279) * t202 + (MDP(5) * t274 + (MDP(6) * t274 - MDP(9) * t255) * qJ(2)) * qJD(1) + (MDP(16) * qJD(4) * t232 - MDP(15) * t221 - t250) * t210 + (t249 + t221 * MDP(16) + (MDP(15) * t232 + 0.2e1 * t247 * t257) * qJD(4)) * t207 + t275; (t207 * t223 - t210 * t286) * MDP(15) + (t207 * t286 + t210 * t223) * MDP(16) + (t151 * t287 + (t202 ^ 2 * t178 - t247 * t167) * t208) * MDP(22) + ((-t247 * t219 - t150) * t211 + (t169 * t247 + t202 * t282) * t208) * MDP(23) + (-qJ(2) * MDP(6) - MDP(5)) * qJD(1) ^ 2 + (-MDP(8) * t208 - MDP(9) * t211) * t234; t281 * MDP(8) + (-t172 * t247 - t160) * MDP(9) + t241 * t285 + t210 * t250 - t207 * t249 + (t207 * t224 - t210 * t225) * MDP(15) + (t207 * t225 + t210 * t224) * MDP(16) + (-t199 * t151 + t233 * t167 + t172 * t282 + t226) * MDP(22) + (t199 * t150 - t233 * t169 - t172 * t219 + t227) * MDP(23) + ((t182 * t206 - t183 * t209 + (t189 * t206 - t190 * t209) * qJD(5)) * MDP(22) + (t182 * t209 + t183 * t206 - (-t189 * t209 - t190 * t206) * qJD(5)) * MDP(23) + t279) * t202 - t275; (t167 * t246 - (t155 * t206 - t262) * t202 + (t206 * t238 - t262) * qJD(5) + t222) * MDP(22) + (-t169 * t246 + (qJD(5) * t238 - t155 * t202 - t145) * t209 + t217) * MDP(23) + t218 + (t207 * MDP(15) + t210 * MDP(16)) * (t163 * t247 - t160) + (-t241 + t277) * t234; (t222 + t276 * (-t153 * t206 - t262)) * MDP(22) + ((-t276 * t153 - t145) * t209 + t217) * MDP(23) + t218;];
tauc = t1;
