% Calculate minimal parameter regressor of Coriolis joint torque vector for
% S4RRRR3
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [4x1]
%   Generalized joint coordinates (joint angles)
% qJD [4x1]
%   Generalized joint velocities
% pkin [7x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,d1,d2,d3,d4]';
% MDP [24x1]
%   Minimal dynamic parameter vector (fixed base model)
%   see S4RRRR3_convert_par2_MPV_fixb.m
% 
% Output:
% tauc [4x1]
%   joint torques required to compensate Coriolis and centrifugal load

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 17:24
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauc = S4RRRR3_coriolisvecJ_fixb_mdp_slag_vp(qJ, qJD, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(4,1),zeros(7,1),zeros(24,1)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4RRRR3_coriolisvecJ_fixb_mdp_slag_vp: qJ has to be [4x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [4 1]), ...
  'S4RRRR3_coriolisvecJ_fixb_mdp_slag_vp: qJD has to be [4x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [7 1]), ...
  'S4RRRR3_coriolisvecJ_fixb_mdp_slag_vp: pkin has to be [7x1] (double)');
assert(isreal(MDP) && all(size(MDP) == [24 1]), ...
  'S4RRRR3_coriolisvecJ_fixb_mdp_slag_vp: MDP has to be [24x1] (double)'); 

%% Symbolic Calculation
% From coriolisvec_joint_fixb_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 17:24:40
% EndTime: 2019-12-31 17:24:43
% DurationCPUTime: 1.45s
% Computational Cost: add. (1104->187), mult. (2948->267), div. (0->0), fcn. (2076->6), ass. (0->110)
t239 = cos(qJ(2));
t292 = pkin(5) + pkin(6);
t221 = t292 * t239;
t217 = qJD(1) * t221;
t238 = cos(qJ(3));
t207 = t238 * t217;
t236 = sin(qJ(2));
t220 = t292 * t236;
t215 = qJD(1) * t220;
t288 = qJD(2) * pkin(2);
t209 = -t215 + t288;
t235 = sin(qJ(3));
t251 = -t235 * t209 - t207;
t277 = qJD(1) * t239;
t266 = t238 * t277;
t278 = qJD(1) * t236;
t267 = t235 * t278;
t200 = -t266 + t267;
t290 = t200 * pkin(7);
t167 = -t251 - t290;
t237 = cos(qJ(4));
t194 = t237 * t200;
t202 = -t235 * t277 - t238 * t278;
t234 = sin(qJ(4));
t180 = t234 * t202 - t194;
t229 = -t239 * pkin(2) - pkin(1);
t219 = t229 * qJD(1);
t190 = t200 * pkin(3) + t219;
t274 = qJD(4) * t234;
t301 = t167 * t274 - t190 * t180;
t231 = qJD(2) + qJD(3);
t272 = qJD(1) * qJD(2);
t265 = t239 * t272;
t183 = qJD(3) * t266 - t231 * t267 + t238 * t265;
t268 = qJD(2) * t292;
t254 = qJD(1) * t268;
t210 = t236 * t254;
t211 = t239 * t254;
t258 = t235 * t210 - t238 * t211;
t245 = t251 * qJD(3) + t258;
t157 = -t183 * pkin(7) + t245;
t271 = -qJD(3) - qJD(4);
t230 = qJD(2) - t271;
t300 = (-t167 * t230 - t157) * t234 + t301;
t214 = t235 * t239 + t238 * t236;
t189 = t231 * t214;
t184 = t189 * qJD(1);
t276 = qJD(3) * t235;
t257 = -t235 * t211 - t217 * t276;
t296 = (qJD(3) * t209 - t210) * t238;
t156 = -t184 * pkin(7) + t257 + t296;
t252 = t234 * t200 + t237 * t202;
t249 = -t234 * t156 + t237 * t157 + t190 * t252;
t153 = -qJD(4) * t194 + t237 * t183 - t234 * t184 + t202 * t274;
t243 = t252 * qJD(4) - t234 * t183 - t237 * t184;
t299 = t180 * t252 * MDP(18) + (-t180 ^ 2 + t252 ^ 2) * MDP(19) + (-t180 * t230 + t153) * MDP(20) + (-t230 * t252 + t243) * MDP(21);
t298 = -0.2e1 * t272;
t297 = t236 * MDP(4);
t295 = (t236 ^ 2 - t239 ^ 2) * MDP(5);
t198 = t202 * pkin(7);
t203 = t235 * t217;
t259 = t238 * t209 - t203;
t166 = t198 + t259;
t294 = qJD(1) * t214;
t293 = qJD(4) - t230;
t291 = pkin(2) * t230;
t289 = t202 * pkin(3);
t286 = t219 * t202;
t285 = t235 * t237;
t240 = qJD(2) ^ 2;
t284 = t236 * t240;
t283 = t237 * t167;
t282 = t239 * t240;
t241 = qJD(1) ^ 2;
t281 = t239 * t241;
t280 = -t238 * t215 - t203;
t275 = qJD(3) * t238;
t270 = t236 * t288;
t269 = pkin(2) * t278;
t264 = -pkin(2) * t231 - t209;
t164 = t231 * pkin(3) + t166;
t263 = -pkin(3) * t230 - t164;
t262 = pkin(1) * t298;
t256 = t235 * t215 - t207;
t253 = -t234 * t164 - t283;
t213 = t235 * t236 - t238 * t239;
t186 = t237 * t213 + t234 * t214;
t187 = -t234 * t213 + t237 * t214;
t250 = t235 * t220 - t238 * t221;
t248 = t219 * t200 - t257;
t216 = t236 * t268;
t218 = t239 * t268;
t246 = -t238 * t216 - t235 * t218 - t220 * t275 - t221 * t276;
t244 = t250 * qJD(3) + t235 * t216 - t238 * t218;
t242 = -t202 * t200 * MDP(11) + t183 * MDP(13) + (-t200 ^ 2 + t202 ^ 2) * MDP(12) + (t200 * MDP(13) + (-t202 - t294) * MDP(14)) * t231 + t299;
t228 = t238 * pkin(2) + pkin(3);
t193 = t213 * pkin(3) + t229;
t191 = t269 - t289;
t188 = t231 * t213;
t182 = t189 * pkin(3) + t270;
t173 = -t213 * pkin(7) - t250;
t172 = -t214 * pkin(7) - t238 * t220 - t235 * t221;
t171 = t184 * pkin(3) + qJD(2) * t269;
t170 = t198 + t280;
t169 = t256 + t290;
t161 = t188 * pkin(7) + t244;
t160 = -t189 * pkin(7) + t246;
t159 = t187 * qJD(4) - t234 * t188 + t237 * t189;
t158 = -t186 * qJD(4) - t237 * t188 - t234 * t189;
t1 = [0.2e1 * t265 * t297 + t295 * t298 + MDP(6) * t282 - MDP(7) * t284 + (-pkin(5) * t282 + t236 * t262) * MDP(9) + (pkin(5) * t284 + t239 * t262) * MDP(10) + (t183 * t214 + t202 * t188) * MDP(11) + (-t183 * t213 - t214 * t184 + t188 * t200 + t202 * t189) * MDP(12) + (t229 * t184 + t219 * t189 + (qJD(1) * t213 + t200) * t270) * MDP(16) + (t229 * t183 - t219 * t188 + (-t202 + t294) * t270) * MDP(17) + (t153 * t187 - t158 * t252) * MDP(18) + (-t153 * t186 + t158 * t180 + t159 * t252 + t187 * t243) * MDP(19) + (t190 * t159 + t171 * t186 - t180 * t182 - t193 * t243) * MDP(23) + (t193 * t153 + t190 * t158 + t171 * t187 - t182 * t252) * MDP(24) + (-t188 * MDP(13) - t189 * MDP(14) + t244 * MDP(16) - t246 * MDP(17)) * t231 + (t158 * MDP(20) - t159 * MDP(21) + (-t234 * t160 + t237 * t161 + (-t172 * t234 - t173 * t237) * qJD(4)) * MDP(23) + (-t237 * t160 - t234 * t161 - (t172 * t237 - t173 * t234) * qJD(4)) * MDP(24)) * t230; (t202 * t269 + t280 * t231 + (t264 * qJD(3) + t210) * t238 + t248) * MDP(17) + (-t200 * t269 + t286 - t256 * t231 + (t264 * t235 - t207) * qJD(3) + t258) * MDP(16) - t281 * t297 + t242 + (t191 * t252 + (-t271 * t235 * t291 + t169 * t230 - t157) * t234 + (-qJD(4) * t164 - t156 + (-pkin(2) * t275 - qJD(4) * t228 + t170) * t230) * t237 + t301) * MDP(24) + (t191 * t180 - (t237 * t169 - t234 * t170) * t230 + (-t234 * t238 - t285) * qJD(3) * t291 + ((-pkin(2) * t285 - t228 * t234) * t230 + t253) * qJD(4) + t249) * MDP(23) + t241 * t295 + (t241 * t236 * MDP(9) + MDP(10) * t281) * pkin(1); (-t251 * t231 + t245 + t286) * MDP(16) + (t259 * t231 + t248 - t296) * MDP(17) + (-t180 * t289 - (-t234 * t166 - t283) * t230 + (t263 * t234 - t283) * qJD(4) + t249) * MDP(23) + (-t252 * t289 + (t263 * qJD(4) + t166 * t230 - t156) * t237 + t300) * MDP(24) + t242; (t293 * t253 + t249) * MDP(23) + ((-t293 * t164 - t156) * t237 + t300) * MDP(24) + t299;];
tauc = t1;
