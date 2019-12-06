% Calculate minimal parameter regressor of Coriolis joint torque vector for
% S5RPRRR1
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [1x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[dummy]';
% MDP [27x1]
%   Minimal dynamic parameter vector (fixed base model)
%   see S5RPRRR1_convert_par2_MPV_fixb.m
% 
% Output:
% tauc [5x1]
%   joint torques required to compensate Coriolis and centrifugal load

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-05 18:10
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauc = S5RPRRR1_coriolisvecJ_fixb_mdp_slag_vp(qJ, qJD, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(1,1),zeros(27,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPRRR1_coriolisvecJ_fixb_mdp_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RPRRR1_coriolisvecJ_fixb_mdp_slag_vp: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [1 1]), ...
  'S5RPRRR1_coriolisvecJ_fixb_mdp_slag_vp: pkin has to be [1x1] (double)');
assert(isreal(MDP) && all(size(MDP) == [27 1]), ...
  'S5RPRRR1_coriolisvecJ_fixb_mdp_slag_vp: MDP has to be [27x1] (double)'); 

%% Symbolic Calculation
% From coriolisvec_joint_fixb_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 18:09:59
% EndTime: 2019-12-05 18:10:06
% DurationCPUTime: 2.38s
% Computational Cost: add. (915->289), mult. (2525->433), div. (0->0), fcn. (1816->6), ass. (0->127)
t206 = cos(qJ(3));
t253 = qJD(1) * qJD(3);
t233 = t206 * t253;
t203 = sin(qJ(3));
t252 = qJD(2) * qJD(1);
t235 = t203 * t252;
t303 = qJ(2) * t233 + t235;
t199 = t203 ^ 2;
t302 = MDP(8) * (-t206 ^ 2 + t199);
t205 = cos(qJ(4));
t202 = sin(qJ(4));
t266 = qJD(3) * t202;
t269 = qJD(1) * t203;
t181 = t205 * t269 + t266;
t268 = qJD(1) * t206;
t195 = -qJD(4) + t268;
t201 = sin(qJ(5));
t204 = cos(qJ(5));
t164 = t181 * t201 + t204 * t195;
t197 = t205 * qJD(3);
t179 = t202 * t269 - t197;
t175 = qJD(5) + t179;
t301 = t164 * t175;
t300 = MDP(6) * qJ(2) + MDP(5);
t298 = qJ(2) * t203;
t258 = qJD(4) * t206;
t213 = -t203 * t197 - t202 * t258;
t251 = qJD(2) * qJD(4);
t232 = t205 * t251;
t254 = t205 * qJD(2);
t160 = t232 + (t213 * qJ(2) + t206 * t254) * qJD(1);
t297 = t160 * t204;
t261 = qJD(4) * t202;
t239 = t203 * t261;
t249 = qJD(3) * qJD(4);
t274 = (t233 + t249) * t205;
t169 = -qJD(1) * t239 + t274;
t296 = t169 * t202;
t259 = qJD(4) * t205;
t237 = t203 * t259;
t263 = qJD(3) * t206;
t212 = t202 * t263 + t237;
t170 = t212 * qJD(1) + t202 * t249;
t295 = t170 * t205;
t294 = t175 * t201;
t293 = t179 * t195;
t292 = t179 * t203;
t291 = t181 * t203;
t290 = t181 * t205;
t289 = t195 * t206;
t288 = t201 * t170;
t287 = t201 * t203;
t234 = t203 * t253;
t255 = qJD(5) * t204;
t245 = t204 * t169 - t195 * t255 + t201 * t234;
t256 = qJD(5) * t201;
t156 = -t181 * t256 + t245;
t286 = t202 * t156;
t285 = t202 * t203;
t284 = t202 * t206;
t283 = t203 * t170;
t282 = t203 * t205;
t207 = qJD(3) ^ 2;
t281 = t203 * t207;
t280 = t204 * t170;
t279 = t204 * t206;
t278 = t205 * t206;
t277 = t206 * t207;
t241 = t202 * t268;
t276 = -t195 * t241 + t205 * t234;
t275 = t303 * t204;
t272 = MDP(13) * t206;
t166 = t181 * t204 - t195 * t201;
t271 = MDP(21) * t166;
t270 = qJD(1) * t199;
t267 = qJD(2) * t202;
t265 = qJD(3) * t203;
t264 = qJD(3) * t204;
t262 = qJD(4) * t195;
t260 = qJD(4) * t204;
t257 = qJD(5) * t166;
t250 = qJD(3) * MDP(18);
t247 = t195 * t278;
t208 = qJD(1) ^ 2;
t246 = t203 * t206 * t208;
t244 = qJ(2) * qJD(1) * t175;
t243 = qJ(2) * t269;
t242 = t166 * t268;
t240 = t205 * t268;
t238 = t202 * t260;
t236 = t204 * t269;
t231 = MDP(17) * t266;
t230 = (qJD(4) + t195) * t203;
t229 = -t181 + t266;
t228 = -qJD(5) + t197;
t227 = qJ(2) * t246;
t225 = qJD(5) * t243;
t171 = t201 * t240 - t236;
t223 = t201 * t259 - t171;
t217 = t204 * t278 + t287;
t172 = t217 * qJD(1);
t222 = t204 * t259 - t172;
t221 = (qJ(2) * t234 - t251) * t202;
t220 = (-t195 + t268) * t203;
t219 = t270 - t289;
t218 = t270 + t289;
t174 = -t201 * t206 + t204 * t282;
t173 = t201 * t282 + t279;
t216 = -t175 * t255 - t288;
t215 = t175 * t256 - t280;
t214 = t218 * t205;
t178 = qJ(2) * t240 + t267;
t168 = t178 * t204 + t201 * t243;
t211 = t201 * t263 + t203 * t255;
t210 = -t202 * t256 + t222;
t209 = qJD(4) * t164 + t216;
t188 = t204 * t234;
t177 = qJ(2) * t241 - t254;
t167 = qJ(2) * t236 - t178 * t201;
t162 = t166 * t261;
t161 = (qJ(2) * t259 + t267) * t268 - t221;
t159 = t228 * t279 + (-t238 + (-qJD(5) * t205 + qJD(3)) * t201) * t203;
t158 = -t201 * t239 - t203 * t264 + t211 * t205 - t206 * t256;
t157 = t169 * t201 - t188 + t257;
t155 = -t168 * qJD(5) - t160 * t201 + t275;
t154 = -t178 * t256 + t297 + (t211 * qJ(2) + qJD(2) * t287) * qJD(1);
t1 = [0.2e1 * t203 * MDP(7) * t233 - 0.2e1 * t253 * t302 + MDP(9) * t277 - MDP(10) * t281 + (t169 * t282 + (t206 * t197 - t239) * t181) * MDP(14) + ((-t179 * t205 - t181 * t202) * t263 + (-t296 - t295 + (t179 * t202 - t290) * qJD(4)) * t203) * MDP(15) + (t195 * t239 - t169 * t206 + (t219 * t205 + t291) * qJD(3)) * MDP(16) + (t195 * t237 + t170 * t206 + (-t219 * t202 - t292) * qJD(3)) * MDP(17) + (-t195 - t268) * t203 * t250 + (-t177 * t265 + t161 * t206 + (t218 * t202 + t292) * qJD(2)) * MDP(19) + (-t178 * t265 + t160 * t206 + (t214 + t291) * qJD(2)) * MDP(20) + (t156 * t174 + t159 * t166) * MDP(21) + (-t156 * t173 - t157 * t174 - t158 * t166 - t159 * t164) * MDP(22) + (t156 * t285 + t159 * t175 + t212 * t166 + t170 * t174) * MDP(23) + (-t157 * t285 - t158 * t175 - t212 * t164 - t170 * t173) * MDP(24) + (t212 * t175 + t202 * t283) * MDP(25) + (t155 * t285 + t177 * t158 + t161 * t173 + t212 * t167 + ((-t201 * t278 + t203 * t204) * t175 + t164 * t284) * qJD(2)) * MDP(26) + (-t154 * t285 + t177 * t159 + t161 * t174 - t212 * t168 + (t166 * t284 - t217 * t175) * qJD(2)) * MDP(27) + 0.2e1 * t300 * t252 + (-MDP(12) * t277 + MDP(13) * t281 + (t283 + qJD(4) * t214 + (t179 * t206 + t202 * t220) * qJD(3)) * MDP(19) + (t169 * t203 - t218 * t261 + (t181 * t206 + t205 * t220) * qJD(3)) * MDP(20) + ((-t164 * t266 + t228 * t294 + t280) * t203 + ((t201 * t261 + t264) * t175 + t202 * t157 + t209 * t205) * t206) * MDP(26) + ((t228 * t175 * t204 - t166 * t266 - t288) * t203 + (-(qJD(3) * t201 - t238) * t175 + t286 + (qJD(4) * t166 + t215) * t205) * t206) * MDP(27)) * qJ(2); t276 * MDP(19) + t171 * t175 * MDP(26) + (t172 * t175 + t162) * MDP(27) - t300 * t208 + (MDP(20) * t262 + (-qJD(4) * t294 - t157) * MDP(26) + (-t175 * t260 - t156) * MDP(27)) * t205 + (-MDP(19) * t292 + (-t247 - t291) * MDP(20) + 0.2e1 * (MDP(12) * t203 + t272) * qJD(3)) * qJD(1) + (MDP(19) * t262 - MDP(20) * t234 + (-t164 * t268 + t209) * MDP(26) + (t215 - t242) * MDP(27)) * t202; -MDP(7) * t246 + t208 * t302 - 0.2e1 * MDP(12) * t235 - 0.2e1 * t252 * t272 + (-t195 * t290 + t296) * MDP(14) + ((t169 + t293) * t205 + (t181 * t195 - t170) * t202) * MDP(15) + (-t195 * t259 + (t229 * t203 + t247) * qJD(1)) * MDP(16) + (t179 * t269 + t195 * t261 + t276) * MDP(17) + t195 * MDP(18) * t269 + (-t202 * t227 + ((t177 - t254) * t203 + ((-t179 - t197) * t206 + t202 * t230) * qJ(2)) * qJD(1)) * MDP(19) + (-t205 * t227 + ((t178 + t267) * t203 + (t205 * t230 + t229 * t206) * qJ(2)) * qJD(1)) * MDP(20) + (t210 * t166 + t204 * t286) * MDP(21) + (t164 * t172 + t166 * t171 + (-t164 * t204 - t166 * t201) * t259 + (-t156 * t201 - t157 * t204 + (t164 * t201 - t166 * t204) * qJD(5)) * t202) * MDP(22) + (-t156 * t205 + t162 + (-t242 + t280) * t202 + t210 * t175) * MDP(23) + (t157 * t205 - t223 * t175 + (t195 * t164 + t216) * t202) * MDP(24) + (-t195 * t202 * t175 - t295) * MDP(25) + (-t155 * t205 + t223 * t177 - t173 * t244 + (t177 * t255 + t167 * qJD(4) + t161 * t201 + (t164 * t298 - t167 * t206) * qJD(1)) * t202) * MDP(26) + (t154 * t205 + t222 * t177 - t174 * t244 + (-t177 * t256 - t168 * qJD(4) + t161 * t204 + (t166 * t298 + t168 * t206) * qJD(1)) * t202) * MDP(27); -t179 ^ 2 * MDP(15) + (t274 - t293) * MDP(16) - qJD(4) * t231 + t221 * MDP(19) + (t177 * t195 - t232) * MDP(20) + (-MDP(19) * t195 - t164 * MDP(26) - t166 * MDP(27)) * t178 + (MDP(14) * t179 + t181 * MDP(15) - t195 * MDP(17) - t166 * MDP(23) + t164 * MDP(24) - t175 * MDP(25) - t167 * MDP(26) + t168 * MDP(27)) * t181 + (t156 * MDP(21) + (-t166 * t179 - t157 - t257) * MDP(22) + t170 * MDP(23) + t161 * MDP(27) - t175 ^ 2 * MDP(24)) * t201 + ((t156 - t301) * MDP(22) + t170 * MDP(24) - t161 * MDP(26) + (MDP(23) * t175 + t271) * t175) * t204 + ((t250 + (-t202 * MDP(16) - t205 * MDP(17)) * qJD(4)) * t203 + ((-t205 * t258 - t291) * MDP(19) + (-t213 + t292) * MDP(20)) * qJ(2) + (-t231 + (-t202 * MDP(19) - t205 * MDP(20)) * qJD(2)) * t206) * qJD(1); t164 * t271 + (-t164 ^ 2 + t166 ^ 2) * MDP(22) + (t245 + t301) * MDP(23) + (t166 * t175 - t181 * t255 + t188) * MDP(24) + t170 * MDP(25) + (-t166 * t177 + t168 * t175 - t178 * t255 + t275) * MDP(26) + (t164 * t177 + t167 * t175 - t204 * t225 - t297) * MDP(27) + (-qJD(5) * t181 * MDP(23) + (qJD(5) * t195 - t169) * MDP(24) + (-t160 - t225) * MDP(26) + (qJD(5) * t178 - t303) * MDP(27)) * t201;];
tauc = t1;
