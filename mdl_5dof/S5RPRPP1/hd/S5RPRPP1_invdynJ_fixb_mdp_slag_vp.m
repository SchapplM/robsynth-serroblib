% Calculate vector of inverse dynamics joint torques for
% S5RPRPP1
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% qJDD [5x1]
%   Generalized joint accelerations
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d3,theta2,theta4]';
% MDP [17x1]
%   Minimal dynamic parameter vector (fixed base model)
%   see S5RPRPP1_convert_par2_MPV_fixb.m
% 
% Output:
% tau [5x1]
%   joint torques of inverse dynamics (contains inertial, gravitational Coriolis and centrifugal forces)

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 18:09
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tau = S5RPRPP1_invdynJ_fixb_mdp_slag_vp(qJ, qJD, qJDD, g, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(5,1),zeros(3,1),zeros(8,1),zeros(17,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPRPP1_invdynJ_fixb_mdp_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RPRPP1_invdynJ_fixb_mdp_slag_vp: qJD has to be [5x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [5 1]), ...
  'S5RPRPP1_invdynJ_fixb_mdp_slag_vp: qJDD has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RPRPP1_invdynJ_fixb_mdp_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RPRPP1_invdynJ_fixb_mdp_slag_vp: pkin has to be [8x1] (double)');
assert(isreal(MDP) && all(size(MDP) == [17 1]), ...
  'S5RPRPP1_invdynJ_fixb_mdp_slag_vp: MDP has to be [17x1] (double)'); 

%% Symbolic Calculation
% From invdyn_joint_fixb_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 18:09:18
% EndTime: 2019-12-31 18:09:21
% DurationCPUTime: 1.66s
% Computational Cost: add. (1404->249), mult. (2849->311), div. (0->0), fcn. (1786->12), ass. (0->113)
t237 = sin(pkin(7));
t219 = pkin(1) * t237 + pkin(6);
t289 = qJ(4) + t219;
t307 = MDP(12) + MDP(15);
t232 = qJ(1) + pkin(7);
t224 = sin(t232);
t226 = cos(t232);
t265 = g(1) * t224 - g(2) * t226;
t239 = cos(pkin(7));
t221 = -pkin(1) * t239 - pkin(2);
t243 = cos(qJ(3));
t229 = t243 * pkin(3);
t305 = t221 - t229;
t231 = qJ(3) + pkin(8);
t223 = sin(t231);
t225 = cos(t231);
t304 = pkin(4) * t225 + qJ(5) * t223;
t266 = g(1) * t226 + g(2) * t224;
t236 = sin(pkin(8));
t238 = cos(pkin(8));
t241 = sin(qJ(3));
t203 = t236 * t243 + t238 * t241;
t197 = t203 * qJD(3);
t278 = qJDD(1) * t243;
t212 = t238 * t278;
t279 = qJDD(1) * t241;
t180 = qJD(1) * t197 + t236 * t279 - t212;
t280 = qJD(1) * qJD(3);
t273 = t241 * t280;
t211 = t236 * t273;
t253 = qJDD(1) * t203 - t211;
t272 = t243 * t280;
t181 = t238 * t272 + t253;
t198 = t203 * qJD(1);
t303 = pkin(4) * t180 - qJ(5) * t181 - qJD(5) * t198;
t194 = qJD(1) * t305 + qJD(4);
t284 = qJD(1) * t243;
t274 = t238 * t284;
t285 = qJD(1) * t241;
t195 = t236 * t285 - t274;
t174 = pkin(4) * t195 - qJ(5) * t198 + t194;
t227 = t243 * qJDD(2);
t207 = t219 * qJDD(1);
t251 = qJ(4) * qJDD(1) + qJD(1) * qJD(4) + qJD(2) * qJD(3) + t207;
t268 = t289 * qJD(1);
t260 = t268 * qJD(3);
t165 = qJDD(3) * pkin(3) - t241 * t251 - t243 * t260 + t227;
t168 = (qJDD(2) - t260) * t241 + t251 * t243;
t154 = t165 * t238 - t168 * t236;
t271 = -qJDD(5) + t154;
t302 = -g(3) * t225 - t174 * t198 + t266 * t223 + t271;
t193 = t198 ^ 2;
t296 = g(3) * t243;
t293 = qJDD(3) * pkin(4);
t189 = qJD(2) * t241 + t243 * t268;
t292 = t189 * t236;
t291 = t236 * t241;
t184 = t238 * t189;
t290 = t238 * t243;
t288 = qJDD(2) - g(3);
t155 = t165 * t236 + t168 * t238;
t188 = qJD(2) * t243 - t241 * t268;
t186 = qJD(3) * pkin(3) + t188;
t167 = t186 * t236 + t184;
t222 = t229 + pkin(2);
t244 = cos(qJ(1));
t287 = pkin(1) * t244 + t222 * t226;
t234 = t241 ^ 2;
t286 = -t243 ^ 2 + t234;
t210 = qJD(1) * t221;
t283 = qJD(3) * t241;
t171 = t188 * t238 - t292;
t281 = qJD(5) - t171;
t277 = pkin(3) * t273 + qJDD(4);
t276 = pkin(3) * t283;
t275 = qJDD(3) * qJ(5) + t155;
t270 = t289 * t241;
t269 = qJD(3) * t289;
t242 = sin(qJ(1));
t264 = g(1) * t242 - g(2) * t244;
t240 = -qJ(4) - pkin(6);
t263 = -pkin(1) * t242 - t226 * t240;
t166 = t186 * t238 - t292;
t259 = t277 - t265;
t256 = -qJD(4) * t241 - t243 * t269;
t254 = -qJD(1) * t210 - t207 + t266;
t252 = 0.2e1 * qJD(3) * t210 - qJDD(3) * t219;
t250 = qJDD(1) * t305 + t277;
t245 = qJD(3) ^ 2;
t249 = -0.2e1 * qJDD(1) * t221 - t219 * t245 + t265;
t190 = qJD(4) * t243 - t241 * t269;
t172 = t190 * t236 - t238 * t256;
t173 = t238 * t190 + t236 * t256;
t201 = t289 * t243;
t178 = t201 * t236 + t238 * t270;
t179 = t238 * t201 - t236 * t270;
t248 = t172 * t198 - t173 * t195 + t178 * t181 - t179 * t180 - t266;
t220 = -pkin(3) * t238 - pkin(4);
t216 = pkin(3) * t236 + qJ(5);
t206 = qJDD(3) * t243 - t241 * t245;
t205 = qJDD(3) * t241 + t243 * t245;
t202 = -t290 + t291;
t200 = qJD(3) * t290 - t236 * t283;
t177 = pkin(4) * t202 - qJ(5) * t203 + t305;
t176 = pkin(3) * t285 + pkin(4) * t198 + qJ(5) * t195;
t170 = t188 * t236 + t184;
t169 = pkin(4) * t197 - qJ(5) * t200 - qJD(5) * t203 + t276;
t162 = qJD(3) * qJ(5) + t167;
t161 = -qJD(3) * pkin(4) + qJD(5) - t166;
t156 = t250 + t303;
t153 = -t271 - t293;
t152 = qJD(3) * qJD(5) + t275;
t1 = [qJDD(1) * MDP(1) + t264 * MDP(2) + (g(1) * t244 + g(2) * t242) * MDP(3) + (t264 + (t237 ^ 2 + t239 ^ 2) * qJDD(1) * pkin(1)) * pkin(1) * MDP(4) + (qJDD(1) * t234 + 0.2e1 * t241 * t272) * MDP(5) + 0.2e1 * (t241 * t278 - t280 * t286) * MDP(6) + t205 * MDP(7) + t206 * MDP(8) + (t241 * t252 + t243 * t249) * MDP(10) + (-t241 * t249 + t243 * t252) * MDP(11) + (-t154 * t203 - t155 * t202 - t166 * t200 - t167 * t197 + t248) * MDP(12) + (t155 * t179 + t167 * t173 - t154 * t178 - t166 * t172 + t250 * t305 + t194 * t276 - g(1) * (-t222 * t224 + t263) - g(2) * (-t224 * t240 + t287)) * MDP(13) + (-qJD(3) * t172 - qJDD(3) * t178 + t156 * t202 + t169 * t195 + t174 * t197 + t177 * t180 + t225 * t265) * MDP(14) + (-t152 * t202 + t153 * t203 + t161 * t200 - t162 * t197 + t248) * MDP(15) + (qJD(3) * t173 + qJDD(3) * t179 - t156 * t203 - t169 * t198 - t174 * t200 - t177 * t181 + t223 * t265) * MDP(16) + (t152 * t179 + t162 * t173 + t156 * t177 + t174 * t169 + t153 * t178 + t161 * t172 - g(1) * t263 - g(2) * (t226 * t304 + t287) + (-g(1) * (-t222 - t304) + g(2) * t240) * t224) * MDP(17); t288 * MDP(4) + t206 * MDP(10) - t205 * MDP(11) + (-t154 * t202 + t155 * t203 - t166 * t197 + t167 * t200 - g(3)) * MDP(13) + (-qJD(3) * t197 - qJDD(3) * t202) * MDP(14) + (qJD(3) * t200 + qJDD(3) * t203) * MDP(16) + (t152 * t203 + t153 * t202 + t161 * t197 + t162 * t200 - g(3)) * MDP(17) + t307 * (-t180 * t203 + t181 * t202 - t195 * t200 + t197 * t198); MDP(7) * t279 + MDP(8) * t278 + qJDD(3) * MDP(9) + (t241 * t254 + t227 - t296) * MDP(10) + (-t241 * t288 + t243 * t254) * MDP(11) + ((t167 - t170) * t198 + (-t166 + t171) * t195 + (-t180 * t236 - t181 * t238) * pkin(3)) * MDP(12) + (t166 * t170 - t167 * t171 + (-t296 + t154 * t238 + t155 * t236 + (-qJD(1) * t194 + t266) * t241) * pkin(3)) * MDP(13) + (qJD(3) * t170 - t176 * t195 + (pkin(4) - t220) * qJDD(3) + t302) * MDP(14) + (-t180 * t216 + t181 * t220 + (t162 - t170) * t198 + (t161 - t281) * t195) * MDP(15) + (-g(3) * t223 + qJDD(3) * t216 - t174 * t195 + t176 * t198 - t266 * t225 + (0.2e1 * qJD(5) - t171) * qJD(3) + t275) * MDP(16) + (t152 * t216 + t153 * t220 - t174 * t176 - t161 * t170 - g(3) * (t229 + t304) + t281 * t162 + t266 * (pkin(3) * t241 + pkin(4) * t223 - qJ(5) * t225)) * MDP(17) + (-MDP(5) * t241 * t243 + MDP(6) * t286) * qJD(1) ^ 2; (t166 * t198 + t167 * t195 + t259) * MDP(13) - t212 * MDP(14) + t211 * MDP(16) + (-t161 * t198 + t162 * t195 + t259 + t303) * MDP(17) + ((t236 * t284 + t238 * t285 + t198) * MDP(14) + (t195 - t274) * MDP(16)) * qJD(3) + (MDP(14) * t291 - MDP(16) * t203 + (MDP(13) + MDP(17)) * t305) * qJDD(1) + t307 * (-t195 ^ 2 - t193); (t195 * t198 - qJDD(3)) * MDP(14) + ((t195 + t274) * qJD(3) + t253) * MDP(15) + (-t193 - t245) * MDP(16) + (-qJD(3) * t162 - t293 - t302) * MDP(17);];
tau = t1;
