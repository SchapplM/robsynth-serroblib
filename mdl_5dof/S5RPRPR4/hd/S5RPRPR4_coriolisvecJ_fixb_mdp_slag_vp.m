% Calculate minimal parameter regressor of Coriolis joint torque vector for
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
% MDP [20x1]
%   Minimal dynamic parameter vector (fixed base model)
%   see S5RPRPR4_convert_par2_MPV_fixb.m
% 
% Output:
% tauc [5x1]
%   joint torques required to compensate Coriolis and centrifugal load

% Quelle: HybrDyn-Toolbox
% Datum: 2020-01-03 11:40
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauc = S5RPRPR4_coriolisvecJ_fixb_mdp_slag_vp(qJ, qJD, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(9,1),zeros(20,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPRPR4_coriolisvecJ_fixb_mdp_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RPRPR4_coriolisvecJ_fixb_mdp_slag_vp: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5RPRPR4_coriolisvecJ_fixb_mdp_slag_vp: pkin has to be [9x1] (double)');
assert(isreal(MDP) && all(size(MDP) == [20 1]), ...
  'S5RPRPR4_coriolisvecJ_fixb_mdp_slag_vp: MDP has to be [20x1] (double)'); 

%% Symbolic Calculation
% From coriolisvec_joint_fixb_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2020-01-03 11:39:22
% EndTime: 2020-01-03 11:39:28
% DurationCPUTime: 1.38s
% Computational Cost: add. (1054->170), mult. (2594->251), div. (0->0), fcn. (1806->8), ass. (0->97)
t232 = sin(pkin(9));
t237 = sin(qJ(3));
t234 = cos(pkin(9));
t239 = cos(qJ(3));
t267 = t234 * t239;
t219 = -t232 * t237 + t267;
t212 = t219 * qJD(1);
t238 = cos(qJ(5));
t203 = t238 * t212;
t220 = t232 * t239 + t234 * t237;
t213 = t220 * qJD(3);
t206 = qJD(1) * t213;
t257 = qJD(1) * qJD(3);
t254 = t237 * t257;
t207 = -t232 * t254 + t257 * t267;
t214 = t220 * qJD(1);
t236 = sin(qJ(5));
t261 = qJD(5) * t236;
t145 = qJD(5) * t203 - t236 * t206 + t238 * t207 - t214 * t261;
t248 = t236 * t212 + t238 * t214;
t146 = t248 * qJD(5) + t238 * t206 + t207 * t236;
t169 = -t214 * t236 + t203;
t229 = qJD(3) + qJD(5);
t269 = t169 * t229;
t270 = t248 * t229;
t282 = (-t146 + t270) * MDP(17) + (-t169 ^ 2 + t248 ^ 2) * MDP(15) - t169 * t248 * MDP(14) + (t145 - t269) * MDP(16);
t225 = sin(pkin(8)) * pkin(1) + pkin(6);
t264 = qJ(4) + t225;
t250 = t264 * qJD(1);
t193 = t239 * qJD(2) - t250 * t237;
t271 = qJD(3) * pkin(3);
t188 = t193 + t271;
t194 = qJD(2) * t237 + t250 * t239;
t268 = t234 * t194;
t159 = t232 * t188 + t268;
t273 = pkin(7) * t212;
t149 = t159 + t273;
t227 = -cos(pkin(8)) * pkin(1) - pkin(2);
t245 = -pkin(3) * t239 + t227;
t243 = t245 * qJD(1);
t211 = qJD(4) + t243;
t178 = -pkin(4) * t212 + t211;
t281 = t149 * t261 - t178 * t169;
t278 = t239 * MDP(5);
t277 = (t237 ^ 2 - t239 ^ 2) * MDP(6);
t276 = qJD(5) - t229;
t256 = qJD(1) * qJD(4);
t176 = qJD(3) * t193 + t239 * t256;
t177 = -qJD(3) * t194 - t237 * t256;
t150 = -t176 * t232 + t234 * t177;
t143 = -pkin(7) * t207 + t150;
t151 = t234 * t176 + t232 * t177;
t144 = -pkin(7) * t206 + t151;
t275 = t238 * t143 - t236 * t144 - t178 * t248;
t274 = pkin(3) * t232;
t272 = pkin(7) * t214;
t184 = t232 * t194;
t240 = qJD(3) ^ 2;
t266 = t237 * t240;
t265 = t239 * t240;
t161 = t234 * t193 - t184;
t251 = qJD(3) * t264;
t197 = qJD(4) * t239 - t237 * t251;
t198 = -qJD(4) * t237 - t239 * t251;
t163 = t234 * t197 + t232 * t198;
t217 = t264 * t237;
t218 = t264 * t239;
t175 = -t232 * t217 + t234 * t218;
t222 = qJD(1) * t227;
t259 = t237 * qJD(1);
t258 = t239 * MDP(11);
t255 = t237 * t271;
t158 = t234 * t188 - t184;
t160 = -t193 * t232 - t268;
t162 = -t197 * t232 + t234 * t198;
t174 = -t234 * t217 - t218 * t232;
t148 = qJD(3) * pkin(4) + t158 - t272;
t249 = -t236 * t148 - t238 * t149;
t247 = t219 * t238 - t220 * t236;
t180 = t219 * t236 + t220 * t238;
t244 = 0.2e1 * qJD(3) * t222;
t226 = pkin(3) * t234 + pkin(4);
t224 = pkin(3) * t254;
t216 = t219 * qJD(3);
t196 = pkin(4) * t213 + t255;
t195 = pkin(3) * t259 + pkin(4) * t214;
t192 = -pkin(4) * t219 + t245;
t183 = pkin(4) * t206 + t224;
t165 = pkin(7) * t219 + t175;
t164 = -pkin(7) * t220 + t174;
t157 = -pkin(7) * t213 + t163;
t156 = -pkin(7) * t216 + t162;
t155 = t161 - t272;
t154 = t160 - t273;
t153 = t180 * qJD(5) + t238 * t213 + t216 * t236;
t152 = t247 * qJD(5) - t213 * t236 + t216 * t238;
t1 = [0.2e1 * t254 * t278 - 0.2e1 * t257 * t277 + MDP(7) * t265 - MDP(8) * t266 + (-t225 * t265 + t237 * t244) * MDP(10) + (t225 * t266 + t239 * t244) * MDP(11) + (-t150 * t220 + t151 * t219 - t158 * t216 - t159 * t213 - t162 * t214 + t163 * t212 - t174 * t207 - t175 * t206) * MDP(12) + (t150 * t174 + t151 * t175 + t158 * t162 + t159 * t163 + (t211 + t243) * t255) * MDP(13) + (t145 * t180 + t152 * t248) * MDP(14) + (t145 * t247 - t146 * t180 + t152 * t169 - t153 * t248) * MDP(15) + (t192 * t146 + t178 * t153 - t169 * t196 - t183 * t247) * MDP(19) + (t192 * t145 + t178 * t152 + t183 * t180 + t196 * t248) * MDP(20) + (t152 * MDP(16) - t153 * MDP(17) + (t156 * t238 - t157 * t236 + (-t164 * t236 - t165 * t238) * qJD(5)) * MDP(19) + (-t156 * t236 - t157 * t238 - (t164 * t238 - t165 * t236) * qJD(5)) * MDP(20)) * t229; (-t206 * t220 - t207 * t219 + t212 * t216 + t213 * t214) * MDP(12) + (t150 * t219 + t151 * t220 - t158 * t213 + t159 * t216) * MDP(13) + (-t237 * MDP(10) - t258) * t240 + (-t153 * MDP(19) - t152 * MDP(20)) * t229; ((t159 + t160) * t214 + (t158 - t161) * t212 + (-t206 * t232 - t207 * t234) * pkin(3)) * MDP(12) + (-t158 * t160 - t159 * t161 + (t150 * t234 + t151 * t232 - t211 * t259) * pkin(3)) * MDP(13) + (t195 * t169 - (t154 * t238 - t155 * t236) * t229 + ((-t226 * t236 - t238 * t274) * t229 + t249) * qJD(5) + t275) * MDP(19) + (-t238 * t144 - t236 * t143 - t195 * t248 + (t154 * t236 + t155 * t238) * t229 + (-(t226 * t238 - t236 * t274) * t229 - t238 * t148) * qJD(5) + t281) * MDP(20) + (-t237 * t278 + t277) * qJD(1) ^ 2 + (-MDP(10) * t259 - qJD(1) * t258) * t222 + t282; (-t212 ^ 2 - t214 ^ 2) * MDP(12) + (t158 * t214 - t159 * t212 + t224) * MDP(13) + (t146 + t270) * MDP(19) + (t145 + t269) * MDP(20); (t249 * t276 + t275) * MDP(19) + ((-t149 * t229 - t143) * t236 + (-t148 * t276 - t144) * t238 + t281) * MDP(20) + t282;];
tauc = t1;
