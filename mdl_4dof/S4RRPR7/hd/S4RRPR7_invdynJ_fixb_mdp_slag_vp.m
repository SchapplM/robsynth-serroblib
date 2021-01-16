% Calculate vector of inverse dynamics joint torques for
% S4RRPR7
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [4x1]
%   Generalized joint coordinates (joint angles)
% qJD [4x1]
%   Generalized joint velocities
% qJDD [4x1]
%   Generalized joint accelerations
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [7x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,d1,d2,d4,theta3]';
% MDP [21x1]
%   Minimal dynamic parameter vector (fixed base model)
%   see S4RRPR7_convert_par2_MPV_fixb.m
% 
% Output:
% tau [4x1]
%   joint torques of inverse dynamics (contains inertial, gravitational Coriolis and centrifugal forces)

% Quelle: HybrDyn-Toolbox
% Datum: 2021-01-15 10:57
% Revision: d12c3222fdeb2c5f3b3c8fa5751e113be2fc3aae (2021-01-15)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tau = S4RRPR7_invdynJ_fixb_mdp_slag_vp(qJ, qJD, qJDD, g, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(4,1),zeros(4,1),zeros(3,1),zeros(7,1),zeros(21,1)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4RRPR7_invdynJ_fixb_mdp_slag_vp: qJ has to be [4x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [4 1]), ...
  'S4RRPR7_invdynJ_fixb_mdp_slag_vp: qJD has to be [4x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [4 1]), ...
  'S4RRPR7_invdynJ_fixb_mdp_slag_vp: qJDD has to be [4x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S4RRPR7_invdynJ_fixb_mdp_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [7 1]), ...
  'S4RRPR7_invdynJ_fixb_mdp_slag_vp: pkin has to be [7x1] (double)');
assert(isreal(MDP) && all(size(MDP) == [21 1]), ...
  'S4RRPR7_invdynJ_fixb_mdp_slag_vp: MDP has to be [21x1] (double)'); 

%% Symbolic Calculation
% From invdyn_joint_fixb_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2021-01-15 10:56:43
% EndTime: 2021-01-15 10:56:49
% DurationCPUTime: 2.44s
% Computational Cost: add. (1244->270), mult. (2959->368), div. (0->0), fcn. (2051->10), ass. (0->119)
t234 = cos(qJ(2));
t284 = cos(pkin(7));
t257 = t284 * t234;
t213 = qJD(1) * t257;
t228 = sin(pkin(7));
t231 = sin(qJ(2));
t258 = t284 * t231;
t205 = t228 * t234 + t258;
t267 = qJD(1) * qJD(2);
t260 = t231 * t267;
t238 = qJDD(1) * t205 - t228 * t260;
t177 = qJD(2) * t213 + t238;
t293 = qJD(2) * qJD(4) + t177;
t229 = -qJ(3) - pkin(5);
t259 = qJD(2) * t229;
t244 = -qJD(3) * t231 + t234 * t259;
t261 = t229 * t231;
t171 = qJDD(2) * pkin(2) + qJD(1) * t244 + qJDD(1) * t261;
t194 = qJD(3) * t234 + t231 * t259;
t210 = t229 * t234;
t178 = qJD(1) * t194 - qJDD(1) * t210;
t154 = t284 * t171 - t228 * t178;
t152 = -qJDD(2) * pkin(3) - t154;
t270 = qJD(1) * t231;
t195 = t228 * t270 - t213;
t189 = qJD(4) + t195;
t198 = t205 * qJD(1);
t216 = pkin(2) * t228 + pkin(6);
t225 = qJ(2) + pkin(7);
t220 = sin(t225);
t221 = cos(t225);
t232 = sin(qJ(1));
t235 = cos(qJ(1));
t251 = g(1) * t235 + g(2) * t232;
t242 = -g(3) * t221 + t220 * t251;
t292 = t242 - t189 * (pkin(2) * t270 + pkin(3) * t198 + pkin(6) * t195 + qJD(4) * t216) - t152;
t160 = t194 * t284 + t228 * t244;
t207 = qJD(1) * t261;
t285 = qJD(2) * pkin(2);
t203 = t207 + t285;
t208 = qJD(1) * t210;
t277 = t228 * t208;
t174 = t203 * t284 + t277;
t166 = -qJD(2) * pkin(3) - t174;
t197 = t205 * qJD(2);
t265 = qJDD(1) * t231;
t249 = -qJDD(1) * t257 + t228 * t265;
t176 = qJD(1) * t197 + t249;
t172 = qJDD(4) + t176;
t219 = pkin(2) * t234 + pkin(1);
t245 = -t228 * t231 + t257;
t173 = -pkin(3) * t245 - pkin(6) * t205 - t219;
t182 = -t210 * t284 + t228 * t261;
t200 = t245 * qJD(2);
t155 = t228 * t171 + t284 * t178;
t153 = qJDD(2) * pkin(6) + t155;
t209 = -qJD(1) * t219 + qJD(3);
t158 = pkin(3) * t195 - pkin(6) * t198 + t209;
t254 = qJD(4) * t158 + t153;
t290 = t152 * t205 + t166 * t200 - t182 * t172 - (qJD(4) * t173 + t160) * t189 + t245 * t254;
t289 = pkin(2) * t231;
t288 = g(3) * t220;
t286 = g(3) * t234;
t230 = sin(qJ(4));
t233 = cos(qJ(4));
t262 = t230 * qJDD(2) + t233 * t293;
t268 = qJD(4) * t230;
t156 = -t198 * t268 + t262;
t283 = t156 * t230;
t282 = t173 * t172;
t183 = -t233 * qJD(2) + t198 * t230;
t281 = t183 * t189;
t280 = t183 * t198;
t185 = qJD(2) * t230 + t198 * t233;
t279 = t185 * t189;
t278 = t185 * t198;
t276 = t230 * t172;
t275 = t230 * t232;
t274 = t230 * t235;
t273 = t232 * t233;
t165 = t233 * t172;
t272 = t233 * t235;
t201 = t284 * t208;
t175 = t228 * t203 - t201;
t226 = t231 ^ 2;
t271 = -t234 ^ 2 + t226;
t269 = qJD(4) * t205;
t264 = qJDD(1) * t234;
t263 = t231 * t285;
t256 = t233 * t189;
t188 = pkin(2) * t260 - qJDD(1) * t219 + qJDD(3);
t151 = pkin(3) * t176 - pkin(6) * t177 + t188;
t167 = qJD(2) * pkin(6) + t175;
t255 = qJD(4) * t167 - t151;
t250 = g(1) * t232 - g(2) * t235;
t248 = t165 + (-t195 * t230 - t268) * t189;
t247 = -0.2e1 * pkin(1) * t267 - pkin(5) * qJDD(2);
t246 = t200 * t233 - t205 * t268;
t180 = t207 * t284 + t277;
t241 = -t216 * t172 + (t166 + t180) * t189;
t236 = qJD(2) ^ 2;
t240 = 0.2e1 * qJDD(1) * pkin(1) - pkin(5) * t236 + t250;
t237 = qJD(1) ^ 2;
t239 = pkin(1) * t237 - pkin(5) * qJDD(1) + t251;
t223 = t233 * qJDD(2);
t217 = -pkin(2) * t284 - pkin(3);
t193 = t221 * t272 + t275;
t192 = -t221 * t274 + t273;
t191 = -t221 * t273 + t274;
t190 = t221 * t275 + t272;
t181 = -t210 * t228 - t229 * t258;
t179 = t207 * t228 - t201;
t162 = pkin(3) * t197 - pkin(6) * t200 + t263;
t159 = t194 * t228 - t244 * t284;
t157 = t185 * qJD(4) + t177 * t230 - t223;
t150 = t158 * t230 + t167 * t233;
t149 = t158 * t233 - t167 * t230;
t148 = t233 * t151;
t1 = [qJDD(1) * MDP(1) + t250 * MDP(2) + t251 * MDP(3) + (qJDD(1) * t226 + 0.2e1 * t234 * t260) * MDP(4) + 0.2e1 * (t231 * t264 - t267 * t271) * MDP(5) + (qJDD(2) * t231 + t234 * t236) * MDP(6) + (qJDD(2) * t234 - t231 * t236) * MDP(7) + (t231 * t247 + t234 * t240) * MDP(9) + (-t231 * t240 + t234 * t247) * MDP(10) + (-qJDD(2) * t181 - t176 * t219 - t188 * t245 + t197 * t209 + t250 * t221 + (t195 * t289 - t159) * qJD(2)) * MDP(11) + (-qJDD(2) * t182 - t177 * t219 + t188 * t205 + t200 * t209 - t250 * t220 + (t198 * t289 - t160) * qJD(2)) * MDP(12) + (-t154 * t205 + t155 * t245 + t159 * t198 - t160 * t195 - t174 * t200 - t175 * t197 - t176 * t182 + t177 * t181 - t251) * MDP(13) + (t155 * t182 + t175 * t160 - t154 * t181 - t174 * t159 - t188 * t219 + t209 * t263 - g(1) * (-t219 * t232 - t229 * t235) - g(2) * (t219 * t235 - t229 * t232)) * MDP(14) + (t156 * t205 * t233 + t185 * t246) * MDP(15) + ((-t183 * t233 - t185 * t230) * t200 + (-t283 - t157 * t233 + (t183 * t230 - t185 * t233) * qJD(4)) * t205) * MDP(16) + (-t156 * t245 + t165 * t205 + t185 * t197 + t189 * t246) * MDP(17) + (-t205 * t276 + t157 * t245 - t183 * t197 + (-t200 * t230 - t233 * t269) * t189) * MDP(18) + (-t172 * t245 + t189 * t197) * MDP(19) + (-g(1) * t191 - g(2) * t193 - t148 * t245 + t149 * t197 + t181 * t157 + t159 * t183 + (t162 * t189 + t282 + (t166 * t205 + t167 * t245 - t182 * t189) * qJD(4)) * t233 + t290 * t230) * MDP(20) + (-g(1) * t190 - g(2) * t192 - t150 * t197 + t181 * t156 + t159 * t185 + (-(-qJD(4) * t182 + t162) * t189 - t282 - t255 * t245 - t166 * t269) * t230 + t290 * t233) * MDP(21); MDP(6) * t265 + MDP(7) * t264 + qJDD(2) * MDP(8) + (t231 * t239 - t286) * MDP(9) + (g(3) * t231 + t234 * t239) * MDP(10) + (t179 * qJD(2) - t209 * t198 + (qJDD(2) * t284 - t195 * t270) * pkin(2) + t242 + t154) * MDP(11) + (t288 + qJD(2) * t180 + t195 * t209 + t251 * t221 + (-qJDD(2) * t228 - t198 * t270) * pkin(2) - t155) * MDP(12) + ((t175 - t179) * t198 + (-t174 + t180) * t195 + (-t176 * t228 - t177 * t284) * pkin(2)) * MDP(13) + (t174 * t179 - t175 * t180 + (t284 * t154 - t286 + t155 * t228 + (-qJD(1) * t209 + t251) * t231) * pkin(2)) * MDP(14) + (t185 * t256 + t283) * MDP(15) + ((t156 - t281) * t233 + (-t157 - t279) * t230) * MDP(16) + (t189 * t256 + t276 - t278) * MDP(17) + (t248 + t280) * MDP(18) - t189 * t198 * MDP(19) + (-t149 * t198 + t217 * t157 - t179 * t183 + t241 * t230 + t233 * t292) * MDP(20) + (t150 * t198 + t217 * t156 - t179 * t185 - t230 * t292 + t241 * t233) * MDP(21) + (-t231 * t234 * MDP(4) + MDP(5) * t271) * t237; (0.2e1 * t198 * qJD(2) + t249) * MDP(11) + ((t213 - t195) * qJD(2) + t238) * MDP(12) + (-t195 ^ 2 - t198 ^ 2) * MDP(13) + (t174 * t198 + t175 * t195 + t188 - t250) * MDP(14) + (t248 - t280) * MDP(20) + (-t189 ^ 2 * t233 - t276 - t278) * MDP(21); t185 * t183 * MDP(15) + (-t183 ^ 2 + t185 ^ 2) * MDP(16) + (t262 + t281) * MDP(17) + (t223 + t279) * MDP(18) + t172 * MDP(19) + (-g(1) * t192 + g(2) * t190 + t150 * t189 - t166 * t185 + t148) * MDP(20) + (g(1) * t193 - g(2) * t191 + t149 * t189 + t166 * t183) * MDP(21) + ((-t153 + t288) * MDP(21) + (-MDP(18) * t198 - MDP(20) * t167 - MDP(21) * t158) * qJD(4)) * t233 + (-qJD(4) * t198 * MDP(17) - t293 * MDP(18) + (-t254 + t288) * MDP(20) + t255 * MDP(21)) * t230;];
tau = t1;
