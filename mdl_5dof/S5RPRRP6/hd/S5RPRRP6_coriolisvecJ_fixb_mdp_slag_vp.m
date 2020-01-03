% Calculate minimal parameter regressor of Coriolis joint torque vector for
% S5RPRRP6
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d3,d4,theta2]';
% MDP [20x1]
%   Minimal dynamic parameter vector (fixed base model)
%   see S5RPRRP6_convert_par2_MPV_fixb.m
% 
% Output:
% tauc [5x1]
%   joint torques required to compensate Coriolis and centrifugal load

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 18:43
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauc = S5RPRRP6_coriolisvecJ_fixb_mdp_slag_vp(qJ, qJD, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(8,1),zeros(20,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPRRP6_coriolisvecJ_fixb_mdp_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RPRRP6_coriolisvecJ_fixb_mdp_slag_vp: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RPRRP6_coriolisvecJ_fixb_mdp_slag_vp: pkin has to be [8x1] (double)');
assert(isreal(MDP) && all(size(MDP) == [20 1]), ...
  'S5RPRRP6_coriolisvecJ_fixb_mdp_slag_vp: MDP has to be [20x1] (double)'); 

%% Symbolic Calculation
% From coriolisvec_joint_fixb_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 18:43:19
% EndTime: 2019-12-31 18:43:23
% DurationCPUTime: 1.58s
% Computational Cost: add. (1127->235), mult. (2637->344), div. (0->0), fcn. (1538->6), ass. (0->113)
t226 = sin(qJ(3));
t300 = MDP(5) * t226;
t221 = t226 ^ 2;
t228 = cos(qJ(3));
t299 = (-t228 ^ 2 + t221) * MDP(6);
t215 = sin(pkin(8)) * pkin(1) + pkin(6);
t209 = t215 * qJD(1);
t298 = qJD(2) * t228 - t226 * t209;
t227 = cos(qJ(4));
t225 = sin(qJ(4));
t266 = qJD(3) * t225;
t269 = qJD(1) * t226;
t204 = t227 * t269 + t266;
t297 = t204 ^ 2;
t296 = pkin(4) * t226;
t295 = -qJ(5) - pkin(7);
t294 = qJ(5) * t226;
t175 = -qJD(3) * pkin(3) - t298;
t293 = t175 * t225;
t292 = t175 * t227;
t219 = t226 * qJD(2);
t264 = qJD(3) * t228;
t178 = qJD(3) * t219 + t209 * t264;
t291 = t178 * t225;
t290 = t178 * t227;
t268 = qJD(1) * t228;
t214 = -qJD(4) + t268;
t289 = t204 * t214;
t288 = t214 * t225;
t287 = t214 * t227;
t286 = t215 * t225;
t216 = -cos(pkin(8)) * pkin(1) - pkin(2);
t196 = -pkin(3) * t228 - pkin(7) * t226 + t216;
t179 = t196 * qJD(1);
t285 = t225 * t179;
t284 = t225 * t228;
t283 = t226 * t227;
t229 = qJD(3) ^ 2;
t282 = t226 * t229;
t281 = t227 * t228;
t280 = t228 * t229;
t183 = t228 * t209 + t219;
t176 = qJD(3) * pkin(7) + t183;
t161 = -t176 * t225 + t227 * t179;
t158 = -qJ(5) * t204 + t161;
t157 = -pkin(4) * t214 + t158;
t279 = t157 - t158;
t253 = qJD(3) * qJD(4);
t245 = t225 * t253;
t262 = qJD(4) * t227;
t248 = t226 * t262;
t251 = t225 * t264;
t172 = t245 + (t248 + t251) * qJD(1);
t256 = t227 * qJD(3);
t202 = t225 * t269 - t256;
t250 = t228 * t256;
t278 = -t172 * t283 - t202 * t250;
t239 = pkin(3) * t226 - pkin(7) * t228;
t207 = t239 * qJD(1);
t277 = t225 * t207 + t227 * t298;
t208 = t239 * qJD(3);
t276 = t196 * t262 + t225 * t208;
t243 = qJD(4) * t295;
t261 = qJD(5) * t227;
t275 = t261 - t277 + (qJ(5) * t268 + t243) * t225;
t233 = -qJ(5) * t281 + t296;
t240 = t227 * t207 - t225 * t298;
t274 = -t233 * qJD(1) - qJD(5) * t225 + t227 * t243 - t240;
t265 = qJD(3) * t226;
t273 = t227 * t208 + t265 * t286;
t206 = t215 * t281;
t272 = t225 * t196 + t206;
t255 = qJD(1) * qJD(3);
t246 = t228 * t255;
t271 = (t246 + t253) * t227;
t210 = qJD(1) * t216;
t263 = qJD(4) * t225;
t167 = pkin(4) * t202 + qJD(5) + t175;
t260 = t167 * qJD(3);
t249 = t226 * t263;
t171 = qJD(1) * t249 - t271;
t259 = t171 * MDP(19);
t258 = t175 * qJD(4);
t257 = t225 * MDP(17);
t254 = qJD(3) * MDP(16);
t177 = t298 * qJD(3);
t195 = qJD(1) * t208;
t252 = t227 * t177 + t179 * t262 + t225 * t195;
t247 = t221 * t255;
t244 = pkin(4) * t225 + t215;
t242 = t214 * t215 + t176;
t241 = t225 * t177 - t227 * t195;
t163 = pkin(4) * t172 + t178;
t162 = t176 * t227 + t285;
t159 = -qJ(5) * t202 + t162;
t238 = -t157 * t227 - t159 * t225;
t237 = t157 * t225 - t159 * t227;
t236 = qJD(1) * t221 - t214 * t228;
t234 = 0.2e1 * qJD(3) * t210;
t232 = t176 * t263 - t252;
t231 = -t162 * qJD(4) - t241;
t212 = t295 * t227;
t211 = t295 * t225;
t199 = t202 ^ 2;
t188 = t204 * t265;
t187 = t227 * t196;
t166 = -t225 * t294 + t272;
t165 = -qJ(5) * t283 + t187 + (-pkin(4) - t286) * t228;
t156 = (-qJ(5) * qJD(4) - qJD(3) * t215) * t283 + (-qJD(5) * t226 + (-qJ(5) * qJD(3) - qJD(4) * t215) * t228) * t225 + t276;
t155 = -t226 * t261 + t233 * qJD(3) + (-t206 + (-t196 + t294) * t225) * qJD(4) + t273;
t154 = -qJ(5) * t172 - qJD(5) * t202 - t232;
t153 = qJ(5) * t171 - qJD(5) * t204 + t255 * t296 + t231;
t1 = [0.2e1 * t246 * t300 - 0.2e1 * t255 * t299 + MDP(7) * t280 - MDP(8) * t282 + (-t215 * t280 + t226 * t234) * MDP(10) + (t215 * t282 + t228 * t234) * MDP(11) + (-t171 * t283 + (-t249 + t250) * t204) * MDP(12) + (-t204 * t248 + (-t204 * t264 + (qJD(4) * t202 + t171) * t226) * t225 + t278) * MDP(13) + (t171 * t228 + t214 * t249 + t236 * t256 + t188) * MDP(14) + (t214 * t248 + t172 * t228 + (-t202 * t226 - t236 * t225) * qJD(3)) * MDP(15) + (-t214 - t268) * t226 * t254 + (-(-t196 * t263 + t273) * t214 + ((t202 * t215 + t293) * qJD(3) + (t242 * t227 + t285) * qJD(4) + t241) * t228 + (t227 * t258 + t215 * t172 + t291 + ((-t215 * t284 + t187) * qJD(1) + t161) * qJD(3)) * t226) * MDP(17) + (t276 * t214 + (-t242 * t263 + (t204 * t215 + t292) * qJD(3) + t252) * t228 + (-t225 * t258 - t215 * t171 + t290 + (-t272 * qJD(1) - t215 * t287 - t162) * qJD(3)) * t226) * MDP(18) + (-t155 * t204 - t156 * t202 + t165 * t171 - t166 * t172 + t238 * t264 + (t237 * qJD(4) - t153 * t227 - t154 * t225) * t226) * MDP(19) + (t153 * t165 + t154 * t166 + t157 * t155 + t159 * t156 + t244 * t228 * t260 + (t167 * pkin(4) * t262 + t163 * t244) * t226) * MDP(20); -t247 * t257 + (-t227 * t247 + t188) * MDP(18) + t278 * MDP(19) + (-t229 * MDP(11) - t172 * MDP(17) + t171 * MDP(18) - t163 * MDP(20) + (t204 * t225 * MDP(19) - t237 * MDP(20) + (t227 * MDP(18) + t257) * t214) * qJD(3)) * t228 + (-t229 * MDP(10) + qJD(3) * t202 * MDP(17) - t225 * t259 + (-t153 * t225 + t154 * t227 + t260) * MDP(20) + ((t202 * t225 + t204 * t227) * MDP(19) + t238 * MDP(20) + (t227 * MDP(17) - t225 * MDP(18)) * t214) * qJD(4)) * t226; (qJD(3) * t183 - t210 * t269 - t178) * MDP(10) - t210 * t268 * MDP(11) + (-t171 * t225 - t204 * t287) * MDP(12) + ((t214 * t202 - t171) * t227 + (-t172 + t289) * t225) * MDP(13) + (-t214 * t262 + (t214 * t281 + (-t204 + t266) * t226) * qJD(1)) * MDP(14) + (t214 * t263 + (-t214 * t284 + (t202 + t256) * t226) * qJD(1)) * MDP(15) + t214 * MDP(16) * t269 + (-pkin(3) * t172 - t290 + t240 * t214 - t183 * t202 + (pkin(7) * t287 + t293) * qJD(4) + (-t161 * t226 + (-pkin(7) * t265 - t175 * t228) * t225) * qJD(1)) * MDP(17) + (pkin(3) * t171 + t291 - t277 * t214 - t183 * t204 + (-pkin(7) * t288 + t292) * qJD(4) + (-t175 * t281 + (-pkin(7) * t256 + t162) * t226) * qJD(1)) * MDP(18) + (t171 * t211 + t172 * t212 - t274 * t204 - t275 * t202 + (t214 * t157 + t154) * t227 + (t214 * t159 - t153) * t225) * MDP(19) + (-t154 * t212 + t153 * t211 + t163 * (-pkin(4) * t227 - pkin(3)) + (-pkin(4) * t288 - t183) * t167 + t275 * t159 + t274 * t157) * MDP(20) + (-t228 * t300 + t299) * qJD(1) ^ 2; (-t199 + t297) * MDP(13) + t271 * MDP(14) + (-t245 - t289) * MDP(15) + (-t162 * t214 - t175 * t204 + t231) * MDP(17) + (-t161 * t214 + t232) * MDP(18) + t279 * MDP(20) * t159 + (t259 + (-t167 * t204 + t153) * MDP(20)) * pkin(4) + (t204 * MDP(12) - t214 * MDP(14) + t175 * MDP(18) - t279 * MDP(19)) * t202 + (-MDP(15) * t251 + (t254 + (-t225 * MDP(14) - t227 * MDP(15)) * qJD(4)) * t226) * qJD(1); (-t199 - t297) * MDP(19) + (t157 * t204 + t159 * t202 + t163) * MDP(20);];
tauc = t1;
