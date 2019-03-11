% Calculate minimal parameter regressor of Coriolis joint torque vector for
% S6RPPRRP5
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d4,d5]';
% MDP [25x1]
%   Minimal dynamic parameter vector (fixed base model)
%   see S6RPPRRP5_convert_par2_MPV_fixb.m
% 
% Output:
% tauc [6x1]
%   joint torques required to compensate Coriolis and centrifugal load

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 02:09
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauc = S6RPPRRP5_coriolisvecJ_fixb_mdp_slag_vp(qJ, qJD, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(8,1),zeros(25,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPPRRP5_coriolisvecJ_fixb_mdp_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RPPRRP5_coriolisvecJ_fixb_mdp_slag_vp: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S6RPPRRP5_coriolisvecJ_fixb_mdp_slag_vp: pkin has to be [8x1] (double)');
assert(isreal(MDP) && all(size(MDP) == [25 1]), ...
  'S6RPPRRP5_coriolisvecJ_fixb_mdp_slag_vp: MDP has to be [25x1] (double)'); 

%% Symbolic Calculation
% From coriolisvec_joint_fixb_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 02:08:50
% EndTime: 2019-03-09 02:08:55
% DurationCPUTime: 2.00s
% Computational Cost: add. (1413->270), mult. (2759->374), div. (0->0), fcn. (1466->4), ass. (0->122)
t216 = -pkin(7) + qJ(2);
t219 = sin(qJ(4));
t221 = cos(qJ(4));
t267 = qJD(4) * t221;
t302 = qJD(2) * t219 + t216 * t267;
t217 = pkin(1) + qJ(3);
t301 = qJD(1) * t217;
t215 = t221 ^ 2;
t300 = (t219 ^ 2 - t215) * MDP(11);
t299 = qJ(2) * MDP(6) + MDP(5) + MDP(7);
t271 = qJD(1) * t219;
t209 = qJD(5) + t271;
t218 = sin(qJ(5));
t220 = cos(qJ(5));
t199 = pkin(4) * t219 - pkin(8) * t221 + t217;
t176 = t199 * qJD(1) - qJD(2);
t211 = qJ(2) * qJD(1) + qJD(3);
t206 = -pkin(7) * qJD(1) + t211;
t198 = t219 * t206;
t183 = qJD(4) * pkin(8) + t198;
t161 = t220 * t176 - t183 * t218;
t270 = qJD(1) * t221;
t251 = t220 * t270;
t195 = t218 * qJD(4) + t251;
t157 = -qJ(6) * t195 + t161;
t156 = pkin(5) * t209 + t157;
t162 = t218 * t176 + t220 * t183;
t252 = t218 * t270;
t261 = t220 * qJD(4);
t193 = t252 - t261;
t158 = -qJ(6) * t193 + t162;
t235 = t156 * t218 - t158 * t220;
t287 = t195 * t218;
t289 = t193 * t220;
t225 = (-t287 + t289) * MDP(24) + t235 * MDP(25) + (t218 * MDP(22) + t220 * MDP(23)) * t209;
t298 = t195 ^ 2;
t297 = -qJ(6) - pkin(8);
t256 = qJD(4) * qJD(5);
t210 = t220 * t256;
t264 = qJD(5) * t221;
t248 = t218 * t264;
t249 = t219 * t261;
t228 = t248 + t249;
t167 = t228 * qJD(1) - t210;
t295 = t167 * t218;
t246 = t218 * t271;
t204 = qJD(4) * t246;
t168 = t195 * qJD(5) - t204;
t294 = t168 * t220;
t259 = qJD(1) * qJD(2);
t268 = qJD(4) * t219;
t174 = t206 * t268 - t221 * t259;
t293 = t174 * t218;
t292 = t174 * t220;
t184 = -qJD(4) * pkin(4) - t206 * t221;
t291 = t184 * t220;
t290 = t193 * t218;
t288 = t195 * t209;
t286 = t195 * t220;
t285 = t195 * t221;
t284 = t209 * t220;
t283 = t218 * t209;
t282 = t218 * t219;
t281 = t218 * t221;
t280 = t219 * t220;
t279 = t220 * t221;
t278 = t156 - t157;
t243 = qJD(5) * t297;
t240 = pkin(4) * t221 + pkin(8) * t219;
t197 = t240 * qJD(1);
t275 = t218 * t197 + t206 * t279;
t277 = -qJ(6) * t246 + qJD(6) * t220 + t218 * t243 - t275;
t239 = t220 * t197 - t206 * t281;
t276 = -qJD(6) * t218 + t220 * t243 - (pkin(5) * t221 + qJ(6) * t280) * qJD(1) - t239;
t274 = t218 * t199 + t216 * t280;
t222 = qJD(4) ^ 2;
t223 = qJD(1) ^ 2;
t272 = -t222 - t223;
t266 = qJD(5) * t218;
t265 = qJD(5) * t220;
t165 = pkin(5) * t193 + qJD(6) + t184;
t263 = t165 * MDP(25);
t262 = t184 * qJD(5);
t207 = -qJD(2) + t301;
t260 = qJD(2) - t207;
t258 = qJD(1) * qJD(4);
t257 = qJD(4) * MDP(21);
t254 = 0.2e1 * qJD(3) * qJD(1);
t175 = t206 * t267 + t219 * t259;
t191 = t240 * qJD(4) + qJD(3);
t177 = t191 * qJD(1);
t253 = -t220 * t175 - t176 * t265 - t218 * t177;
t247 = t220 * t264;
t245 = t221 * t258;
t244 = pkin(5) * t218 - t216;
t242 = t209 * t216 + t183;
t241 = t218 * t191 + t199 * t265 + t220 * t302;
t238 = qJD(2) + t207 + t301;
t173 = t220 * t177;
t226 = -t162 * qJD(5) - t218 * t175 + t173;
t152 = pkin(5) * t245 + qJ(6) * t167 - qJD(6) * t195 + t226;
t229 = t183 * t266 + t253;
t153 = -qJ(6) * t168 - qJD(6) * t193 - t229;
t237 = -t152 * t220 - t153 * t218;
t236 = t156 * t220 + t158 * t218;
t233 = qJD(1) * t215 - t209 * t219;
t231 = -t216 * t222 + t254;
t230 = -pkin(8) * t267 + t184 * t219;
t159 = pkin(5) * t168 + t174;
t227 = -qJD(6) * t221 + (qJ(6) * qJD(4) - qJD(5) * t216) * t219;
t224 = (t286 + t290) * MDP(24) - t236 * MDP(25) + (-t220 * MDP(22) + t218 * MDP(23)) * t209;
t203 = t218 * t245;
t202 = t297 * t220;
t201 = t297 * t218;
t190 = t193 ^ 2;
t188 = t220 * t199;
t181 = t220 * t191;
t166 = -qJ(6) * t281 + t274;
t164 = -qJ(6) * t279 + t188 + (-t216 * t218 + pkin(5)) * t219;
t155 = -qJ(6) * t247 + t227 * t218 + t241;
t154 = pkin(5) * t267 + t181 + t227 * t220 + ((qJ(6) * t221 - t199) * qJD(5) - t302) * t218;
t1 = [MDP(8) * t254 + (qJD(2) * t211 + qJD(3) * t207 + (qJ(2) * qJD(2) + qJD(3) * t217) * qJD(1)) * MDP(9) + 0.2e1 * t258 * t300 - t222 * t221 * MDP(13) + t238 * t267 * MDP(15) + (t231 * t221 - t238 * t268) * MDP(16) + (-t167 * t279 - t228 * t195) * MDP(17) + ((t287 + t289) * t268 + (t295 - t294 + (-t286 + t290) * qJD(5)) * t221) * MDP(18) + (-t209 * t248 + (t233 * t220 + t285) * qJD(4)) * MDP(19) + (-t209 * t247 + (-t193 * t221 - t233 * t218) * qJD(4)) * MDP(20) + (t209 + t271) * t221 * t257 + ((-t199 * t266 + t181) * t209 + (t220 * t262 - qJD(2) * t193 - t216 * t168 + t293 + (-t216 * t283 + (-t216 * t282 + t188) * qJD(1) + t161) * qJD(4)) * t221) * MDP(22) + (-t241 * t209 + (-t218 * t262 - qJD(2) * t195 + t216 * t167 + t292 + (-t274 * qJD(1) - t162) * qJD(4)) * t221) * MDP(23) + (-t154 * t195 - t155 * t193 + t164 * t167 - t166 * t168 + t236 * t268 + (t235 * qJD(5) + t237) * t221) * MDP(24) + (t152 * t164 + t153 * t166 + t156 * t154 + t158 * t155 - t165 * t244 * t268 + (t159 * t244 + t165 * (pkin(5) * t265 - qJD(2))) * t221) * MDP(25) + (-0.2e1 * MDP(10) * t245 - t222 * MDP(12) + t231 * MDP(15) - t167 * MDP(19) - t168 * MDP(20) + (qJD(4) * t216 * t193 + t173 - t242 * t265 + (-qJD(2) * t209 - t184 * qJD(4) - qJD(5) * t176 - t175) * t218) * MDP(22) + (t242 * t266 + (t195 * t216 - t291) * qJD(4) + t253) * MDP(23)) * t219 + 0.2e1 * t299 * t259; t203 * MDP(23) + (-t167 * t220 + t168 * t218) * MDP(24) + t237 * MDP(25) - t299 * t223 + t225 * qJD(5) + ((-qJD(3) - t211) * MDP(9) + (-0.2e1 * qJD(4) * MDP(15) + (t193 - t261) * MDP(22) + t195 * MDP(23) + t263) * t221 + (0.2e1 * qJD(4) * MDP(16) + t225) * t219) * qJD(1); -t223 * MDP(8) + (t260 * MDP(9) + t224) * qJD(1) + (t272 * MDP(16) - t168 * MDP(22) + t167 * MDP(23) - t159 * MDP(25) - qJD(4) * t225) * t221 + (t272 * MDP(15) + (-t294 - t295) * MDP(24) + (-t152 * t218 + t153 * t220) * MDP(25) + t224 * qJD(5) + ((t193 - t252) * MDP(22) + (t195 - t251) * MDP(23) + t263) * qJD(4)) * t219; (t195 * t284 - t295) * MDP(17) + ((-t193 * t209 - t167) * t220 + (-t168 - t288) * t218) * MDP(18) + (t209 * t265 + t203 + (t209 * t280 - t285) * qJD(1)) * MDP(19) + (-t209 * t266 + (-t209 * t282 + (t193 + t261) * t221) * qJD(1)) * MDP(20) - t209 * MDP(21) * t270 + (-pkin(4) * t168 - t292 - t239 * t209 - t193 * t198 + (-pkin(8) * t284 + t184 * t218) * qJD(5) + (-t161 * t221 + t230 * t218) * qJD(1)) * MDP(22) + (pkin(4) * t167 + t293 + t275 * t209 - t195 * t198 + (pkin(8) * t283 + t291) * qJD(5) + (t162 * t221 + t230 * t220) * qJD(1)) * MDP(23) + (t167 * t201 + t168 * t202 - t276 * t195 - t277 * t193 + (-t209 * t156 + t153) * t220 + (-t209 * t158 - t152) * t218) * MDP(24) + (-t153 * t202 + t152 * t201 + t159 * (-pkin(5) * t220 - pkin(4)) + (pkin(5) * t283 - t198) * t165 + t277 * t158 + t276 * t156) * MDP(25) + (MDP(15) * t270 - MDP(16) * t271) * t260 + (t221 * t219 * MDP(10) - t300) * t223; (-t190 + t298) * MDP(18) + t210 * MDP(19) + (-t218 * t256 + t204 + t288) * MDP(20) + (t162 * t209 - t184 * t195 + t226) * MDP(22) + (t161 * t209 + t229) * MDP(23) + t278 * MDP(25) * t158 + (t167 * MDP(24) + (-t165 * t195 + t152) * MDP(25)) * pkin(5) + (t195 * MDP(17) + t209 * MDP(19) + t184 * MDP(23) - t278 * MDP(24)) * t193 + (-MDP(19) * t249 + (t257 + (-t218 * MDP(19) - t220 * MDP(20)) * qJD(5)) * t221) * qJD(1); (-t190 - t298) * MDP(24) + (t156 * t195 + t158 * t193 + t159) * MDP(25);];
tauc  = t1;
