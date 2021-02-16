% Calculate Coriolis joint torque vector for
% S5RPRRP12
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
%   see S5RPRRP12_convert_par2_MPV_fixb.m
% 
% Output:
% tauc [5x1]
%   joint torques required to compensate Coriolis and centrifugal load

% Quelle: HybrDyn-Toolbox
% Datum: 2021-01-15 19:26
% Revision: 24b2e7d74a0c1a3b64fa2f8f5ad758691ad61af3 (2021-01-15)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauc = S5RPRRP12_coriolisvecJ_fixb_mdp_slag_vp(qJ, qJD, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(7,1),zeros(24,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPRRP12_coriolisvecJ_fixb_mdp_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RPRRP12_coriolisvecJ_fixb_mdp_slag_vp: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [7 1]), ...
  'S5RPRRP12_coriolisvecJ_fixb_mdp_slag_vp: pkin has to be [7x1] (double)');
assert(isreal(MDP) && all(size(MDP) == [24 1]), ...
  'S5RPRRP12_coriolisvecJ_fixb_mdp_slag_vp: MDP has to be [24x1] (double)'); 

%% Symbolic Calculation
% From coriolisvec_joint_fixb_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2021-01-15 19:26:01
% EndTime: 2021-01-15 19:26:08
% DurationCPUTime: 2.35s
% Computational Cost: add. (1486->286), mult. (3111->405), div. (0->0), fcn. (1717->4), ass. (0->126)
t255 = 2 * qJD(1);
t222 = cos(qJ(3));
t218 = t222 ^ 2;
t220 = sin(qJ(3));
t308 = (t220 ^ 2 - t218) * MDP(8);
t307 = (MDP(6) * qJ(2) + MDP(5));
t206 = pkin(3) * t220 - pkin(7) * t222 + qJ(2);
t187 = t206 * qJD(1);
t223 = (-pkin(1) - pkin(6));
t212 = t223 * qJD(1) + qJD(2);
t205 = t220 * t212;
t190 = qJD(3) * pkin(7) + t205;
t219 = sin(qJ(4));
t221 = cos(qJ(4));
t167 = t221 * t187 - t190 * t219;
t266 = qJD(1) * t222;
t250 = t221 * t266;
t265 = qJD(3) * t219;
t201 = t250 + t265;
t162 = -qJ(5) * t201 + t167;
t267 = qJD(1) * t220;
t213 = qJD(4) + t267;
t161 = pkin(4) * t213 + t162;
t306 = t161 - t162;
t258 = t221 * qJD(3);
t214 = qJD(4) * t258;
t260 = qJD(4) * t222;
t247 = t219 * t260;
t249 = t220 * t258;
t227 = t247 + t249;
t172 = qJD(1) * t227 - t214;
t251 = t219 * t266;
t199 = t251 - t258;
t305 = -t199 * t213 - t172;
t256 = qJD(1) * qJD(3);
t244 = t220 * t256;
t210 = t219 * t244;
t173 = qJD(4) * t201 - t210;
t304 = t201 * t213 + t173;
t278 = t220 * t223;
t270 = t219 * t206 + t221 * t278;
t303 = MDP(19) + MDP(21);
t302 = pkin(4) * t199;
t301 = -qJ(5) - pkin(7);
t300 = pkin(4) * qJD(1);
t168 = t187 * t219 + t190 * t221;
t163 = -qJ(5) * t199 + t168;
t298 = t163 * t213;
t264 = qJD(3) * t220;
t165 = pkin(4) * t173 + t212 * t264;
t297 = t165 * t219;
t296 = t165 * t221;
t295 = t172 * t219;
t294 = t173 * t221;
t285 = t212 * t222;
t191 = -qJD(3) * pkin(3) - t285;
t293 = t191 * t219;
t291 = t199 * t219;
t290 = t199 * t221;
t288 = t201 * t219;
t287 = t201 * t221;
t286 = t212 * t219;
t284 = t213 * t219;
t283 = t213 * t221;
t282 = t219 * t220;
t281 = t219 * t222;
t280 = t219 * t223;
t279 = t220 * t221;
t277 = t221 * t222;
t276 = t222 * t173;
t225 = qJD(1) ^ 2;
t275 = t222 * t225;
t224 = qJD(3) ^ 2;
t274 = t223 * t224;
t234 = pkin(3) * t222 + pkin(7) * t220;
t204 = t234 * qJD(1);
t189 = t221 * t204;
t242 = qJD(4) * t301;
t273 = -t212 * t281 + t189 + (pkin(4) * t222 + qJ(5) * t279) * qJD(1) + qJD(5) * t219 - t221 * t242;
t245 = t219 * t267;
t259 = qJD(5) * t221;
t271 = t219 * t204 + t212 * t277;
t272 = qJ(5) * t245 - t219 * t242 - t259 + t271;
t268 = -t224 - t225;
t263 = qJD(3) * t222;
t262 = qJD(4) * t219;
t261 = qJD(4) * t221;
t257 = qJD(1) * MDP(18);
t253 = t219 * t278;
t197 = qJD(3) * t234 + qJD(2);
t248 = t222 * t258;
t252 = t219 * t197 + t206 * t261 + t223 * t248;
t246 = t221 * t260;
t243 = pkin(4) - t280;
t241 = -qJD(5) - t302;
t181 = t197 * qJD(1);
t176 = t221 * t181;
t240 = qJ(5) * t172 + t176;
t239 = t213 + t267;
t238 = t199 + t258;
t237 = -t201 + t265;
t236 = qJD(4) * t220 + qJD(1);
t235 = t219 * t181 + t187 * t261 - t190 * t262 + t212 * t248;
t233 = t161 * t221 + t163 * t219;
t232 = t161 * t219 - t163 * t221;
t231 = qJD(1) * t218 - t213 * t220;
t230 = -pkin(7) * t263 + t191 * t220;
t228 = -qJ(5) * t173 + t235;
t226 = (t287 + t291) * MDP(23) - t233 * MDP(24);
t215 = -pkin(4) * t221 - pkin(3);
t209 = t301 * t221;
t208 = t301 * t219;
t198 = (pkin(4) * t219 - t223) * t222;
t196 = t199 ^ 2;
t195 = t221 * t206;
t185 = t221 * t197;
t177 = -pkin(4) * t245 + t205;
t174 = t223 * t264 + (-t219 * t264 + t246) * pkin(4);
t171 = -qJ(5) * t281 + t270;
t170 = t191 - t241;
t169 = -qJ(5) * t277 + t220 * t243 + t195;
t160 = -qJ(5) * t246 + (-qJD(5) * t222 + (qJ(5) * qJD(3) - qJD(4) * t223) * t220) * t219 + t252;
t159 = qJ(5) * t249 + t185 - t270 * qJD(4) + (qJ(5) * t262 + qJD(3) * t243 - t259) * t222;
t156 = -qJD(5) * t199 + t228;
t155 = -qJD(5) * t201 - t168 * qJD(4) + (-t286 + t300) * t263 + t240;
t1 = [-0.2e1 * t222 * MDP(7) * t244 + 0.2e1 * t256 * t308 + (-t220 * t274 + (qJ(2) * t263 + qJD(2) * t220) * t255) * MDP(12) + (-t222 * t274 + (-qJ(2) * t264 + qJD(2) * t222) * t255) * MDP(13) + (-t172 * t277 - t201 * t227) * MDP(14) + ((t288 + t290) * t264 + (t295 - t294 + (-t287 + t291) * qJD(4)) * t222) * MDP(15) + (-t213 * t247 - t172 * t220 + (t201 * t222 + t221 * t231) * qJD(3)) * MDP(16) + (-t213 * t246 - t173 * t220 + (-t199 * t222 - t219 * t231) * qJD(3)) * MDP(17) + t239 * MDP(18) * t263 + (-t223 * t276 + t176 * t220 + t185 * t213 + (-t168 * t220 + t191 * t277 - t213 * t270) * qJD(4) + ((t199 * t223 - t293) * t220 + (-t213 * t280 + (t195 - t253) * qJD(1) + t167) * t222) * qJD(3)) * MDP(19) + (-(-qJD(4) * t253 + t252) * t213 - t235 * t220 + (t223 * t172 - t191 * t262) * t222 + ((-t270 * qJD(1) - t168) * t222 + (t223 * t201 + (-t191 + t285) * t221) * t220) * qJD(3)) * MDP(20) + (t159 * t213 + t173 * t198 + t174 * t199 + (-t170 * t265 + t155) * t220 + (t170 * t261 + t297 + (qJD(1) * t169 + t161) * qJD(3)) * t222) * MDP(21) + (-t160 * t213 - t172 * t198 + t174 * t201 + (-t170 * t258 - t156) * t220 + (-t170 * t262 + t296 + (-qJD(1) * t171 - t163) * qJD(3)) * t222) * MDP(22) + (-t159 * t201 - t160 * t199 + t169 * t172 - t171 * t173 + t233 * t264 + (qJD(4) * t232 - t155 * t221 - t156 * t219) * t222) * MDP(23) + (t155 * t169 + t156 * t171 + t159 * t161 + t160 * t163 + t165 * t198 + t170 * t174) * MDP(24) + (t307 * qJD(2) * t255) + (-t222 * MDP(10) - t220 * MDP(9)) * t224; -(t307 * t225) + t303 * (-t276 - t236 * t283 + (t199 * t220 - t239 * t281) * qJD(3)) + (MDP(20) + MDP(22)) * (t172 * t222 + t236 * t284 + (-t213 * t277 + (t201 - t250) * t220) * qJD(3)) + t226 * qJD(1) + (t268 * MDP(13) - t165 * MDP(24) + ((t288 - t290) * MDP(23) - t232 * MDP(24)) * qJD(3)) * t222 + (t268 * MDP(12) + (-t294 - t295) * MDP(23) + (qJD(3) * t170 - t155 * t219 + t156 * t221) * MDP(24) + t226 * qJD(4)) * t220; t220 * MDP(7) * t275 - t225 * t308 + (t201 * t283 - t295) * MDP(14) + (-t304 * t219 + t305 * t221) * MDP(15) + (t213 * t261 + (t213 * t279 + t222 * t237) * qJD(1)) * MDP(16) + (-t213 * t262 + (-t213 * t282 + t222 * t238) * qJD(1)) * MDP(17) - t213 * t222 * t257 + (-pkin(3) * t173 - t189 * t213 + (t213 * t281 - t220 * t238) * t212 + (-pkin(7) * t283 + t293) * qJD(4) + (-t167 * t222 + t219 * t230) * qJD(1)) * MDP(19) + (pkin(3) * t172 + t271 * t213 + t237 * t205 + (pkin(7) * t284 + t191 * t221) * qJD(4) + (t168 * t222 + t221 * t230) * qJD(1)) * MDP(20) + (-t296 + t173 * t215 - t177 * t199 - t273 * t213 + (t170 + t302) * t262 + (t170 * t282 + (qJD(3) * t208 - t161) * t222) * qJD(1)) * MDP(21) + (t297 - t172 * t215 - t177 * t201 + t272 * t213 + (pkin(4) * t288 + t170 * t221) * qJD(4) + (t170 * t279 + (qJD(3) * t209 + t163) * t222) * qJD(1)) * MDP(22) + (t172 * t208 + t173 * t209 + t273 * t201 + t272 * t199 + (-t161 * t213 + t156) * t221 + (-t155 - t298) * t219) * MDP(23) + (t155 * t208 - t156 * t209 + t165 * t215 + (pkin(4) * t262 - t177) * t170 - t272 * t163 - t273 * t161) * MDP(24) + (t225 * t220 * MDP(13) - MDP(12) * t275) * qJ(2); -t196 * MDP(15) + t214 * MDP(16) + t210 * MDP(17) + (t168 * t213 + t176) * MDP(19) + (t167 * t213 - t235) * MDP(20) + (t240 + t298) * MDP(21) + (t162 * t213 - t228) * MDP(22) + t172 * pkin(4) * MDP(23) + (pkin(4) * t155 + t306 * t163) * MDP(24) + (t213 * MDP(16) + t191 * MDP(20) + (qJD(5) + t170) * MDP(22) - t306 * MDP(23)) * t199 + (-t221 * MDP(16) * t267 + (t257 - MDP(19) * t286 + (-t286 + 0.2e1 * t300) * MDP(21)) * t222) * qJD(3) + (-MDP(16) * t251 - MDP(17) * t201 - t303 * t168) * qJD(4) + (t199 * MDP(14) + t213 * MDP(17) - t191 * MDP(19) + (-t170 + t241) * MDP(21) - pkin(4) * t170 * MDP(24) + (-MDP(22) * pkin(4) + MDP(15)) * t201) * t201; t304 * MDP(21) + t305 * MDP(22) + (-t201 ^ 2 - t196) * MDP(23) + (t161 * t201 + t163 * t199 + t165) * MDP(24);];
tauc = t1;
