% Calculate vector of inverse dynamics joint torques for
% S5PRRPP4
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
% pkin [7x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d2,d3,theta1]';
% MDP [19x1]
%   Minimal dynamic parameter vector (fixed base model)
%   see S5PRRPP4_convert_par2_MPV_fixb.m
% 
% Output:
% tau [5x1]
%   joint torques of inverse dynamics (contains inertial, gravitational Coriolis and centrifugal forces)

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 17:41
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tau = S5PRRPP4_invdynJ_fixb_mdp_slag_vp(qJ, qJD, qJDD, g, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(5,1),zeros(3,1),zeros(7,1),zeros(19,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PRRPP4_invdynJ_fixb_mdp_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5PRRPP4_invdynJ_fixb_mdp_slag_vp: qJD has to be [5x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [5 1]), ...
  'S5PRRPP4_invdynJ_fixb_mdp_slag_vp: qJDD has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5PRRPP4_invdynJ_fixb_mdp_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [7 1]), ...
  'S5PRRPP4_invdynJ_fixb_mdp_slag_vp: pkin has to be [7x1] (double)');
assert(isreal(MDP) && all(size(MDP) == [19 1]), ...
  'S5PRRPP4_invdynJ_fixb_mdp_slag_vp: MDP has to be [19x1] (double)'); 

%% Symbolic Calculation
% From invdyn_joint_fixb_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 17:41:17
% EndTime: 2019-12-31 17:41:19
% DurationCPUTime: 1.72s
% Computational Cost: add. (746->245), mult. (1341->287), div. (0->0), fcn. (666->4), ass. (0->119)
t232 = cos(qJ(3));
t294 = qJ(4) * t232;
t231 = sin(qJ(3));
t299 = pkin(3) + pkin(4);
t302 = t299 * t231;
t247 = t294 - t302;
t218 = t231 * qJ(4);
t262 = pkin(2) + t218;
t306 = qJDD(3) * qJ(4) + qJD(3) * qJD(4);
t272 = qJD(2) * qJD(3);
t263 = t232 * t272;
t271 = qJDD(2) * t231;
t305 = t263 + t271;
t283 = qJD(2) * t231;
t210 = pkin(6) * t283;
t176 = t232 * qJD(1) - t210;
t304 = qJD(4) - t176;
t303 = t299 * qJDD(3);
t166 = qJ(5) * t283 + t176;
t275 = qJD(4) - t166;
t161 = -t299 * qJD(3) + t275;
t282 = qJD(2) * t232;
t177 = pkin(6) * t282 + t231 * qJD(1);
t168 = -qJ(5) * t282 + t177;
t227 = qJD(3) * qJ(4);
t164 = t168 + t227;
t224 = pkin(7) + qJ(2);
t212 = sin(t224);
t213 = cos(t224);
t287 = g(1) * t213 + g(2) * t212;
t301 = 0.2e1 * t306;
t174 = t227 + t177;
t171 = -qJD(3) * pkin(3) + t304;
t222 = t232 * pkin(3);
t248 = -t262 - t222;
t172 = t248 * qJD(2);
t286 = t222 + t218;
t178 = -pkin(2) - t286;
t295 = pkin(6) * qJDD(3);
t300 = (qJD(2) * t178 + t172) * qJD(3) - t295;
t200 = g(1) * t212;
t298 = g(2) * t213;
t221 = t232 * pkin(4);
t297 = pkin(6) - qJ(5);
t296 = pkin(6) * qJD(3);
t230 = qJDD(2) * pkin(2);
t293 = qJDD(3) * pkin(3);
t292 = t212 * t231;
t291 = t212 * t232;
t290 = t213 * t231;
t289 = t213 * t232;
t235 = qJD(2) ^ 2;
t288 = t231 * t235;
t228 = t231 ^ 2;
t229 = t232 ^ 2;
t285 = t228 - t229;
t284 = t228 + t229;
t281 = qJD(3) * t164;
t183 = t297 * t232;
t280 = qJD(3) * t183;
t279 = qJD(3) * t231;
t278 = qJD(4) * t231;
t277 = qJD(5) * t231;
t276 = qJD(5) * t232;
t163 = qJD(5) + (t299 * t232 + t262) * qJD(2);
t274 = qJD(5) + t163;
t273 = qJD(1) * qJD(3);
t270 = qJDD(2) * t232;
t269 = MDP(12) + MDP(16);
t268 = MDP(14) + MDP(17);
t267 = t232 * t288;
t266 = pkin(6) * t270 + t231 * qJDD(1) + t232 * t273;
t265 = t221 + t286;
t264 = t231 * t272;
t261 = t200 - t298;
t173 = pkin(2) + t265;
t259 = qJD(2) * t173 + t163;
t257 = pkin(3) * t289 + t212 * pkin(6) + t262 * t213;
t256 = pkin(6) * t305 - t232 * qJDD(1) + t231 * t273;
t234 = qJD(3) ^ 2;
t255 = pkin(6) * t234 + t298;
t253 = pkin(3) * t231 - t294;
t252 = pkin(3) * t270 + qJ(4) * t305 + qJD(2) * t278 + t230;
t251 = t266 + t306;
t250 = -qJDD(4) - t256;
t249 = g(3) * t231 - t266;
t246 = -0.2e1 * pkin(2) * t272 - t295;
t245 = (t161 * t231 + t164 * t232) * MDP(19);
t244 = g(1) * t290 + g(2) * t292 - g(3) * t232 - t256;
t243 = -t255 + 0.2e1 * t230;
t242 = pkin(4) * t270 + qJDD(5) + t252;
t241 = -qJDD(4) + t244;
t157 = -t299 * t264 + t242;
t162 = t247 * qJD(3) + t278;
t240 = qJD(2) * t162 + qJDD(2) * t173 + t157 - t298;
t239 = qJD(3) * t177 + t244;
t238 = -t241 - t293;
t158 = pkin(3) * t264 - t252;
t170 = t253 * qJD(3) - t278;
t237 = -qJD(2) * t170 - qJDD(2) * t178 - t158 - t255;
t159 = -pkin(6) * t264 + t251;
t160 = -t250 - t293;
t236 = t159 * t232 + t160 * t231 + (t171 * t232 - t174 * t231) * qJD(3);
t197 = t213 * pkin(6);
t193 = qJ(5) * t264;
t190 = g(1) * t291;
t189 = g(1) * t292;
t186 = qJ(4) * t289;
t184 = qJ(4) * t291;
t182 = t297 * t231;
t180 = qJDD(3) * t232 - t231 * t234;
t179 = qJDD(3) * t231 + t232 * t234;
t175 = t253 * qJD(2);
t169 = -t277 + t280;
t167 = -t297 * t279 - t276;
t165 = t247 * qJD(2);
t156 = -qJ(5) * t270 + t193 + (-pkin(6) * t279 - t276) * qJD(2) + t251;
t155 = -qJ(5) * t305 - qJD(2) * t277 - t250 - t303;
t1 = [(qJDD(1) - g(3)) * MDP(1) + (t159 * t231 - t160 * t232 - g(3)) * MDP(15) + (-t155 * t232 + t156 * t231 - g(3)) * MDP(19) + ((t171 * t231 + t174 * t232) * MDP(15) + t245) * qJD(3) + (MDP(10) + t269) * t180 + (-MDP(11) + t268) * t179; qJDD(2) * MDP(2) + t261 * MDP(3) + t287 * MDP(4) + (qJDD(2) * t228 + 0.2e1 * t231 * t263) * MDP(5) + 0.2e1 * (t231 * t270 - t285 * t272) * MDP(6) + t179 * MDP(7) + t180 * MDP(8) + (t246 * t231 + t243 * t232 + t190) * MDP(10) + (-t243 * t231 + t246 * t232 - t189) * MDP(11) + (t300 * t231 + t237 * t232 + t190) * MDP(12) + (t284 * qJDD(2) * pkin(6) + t236 - t287) * MDP(13) + (t237 * t231 - t300 * t232 + t189) * MDP(14) + (t236 * pkin(6) - g(1) * t197 - g(2) * t257 + t158 * t178 + t172 * t170 - t248 * t200) * MDP(15) + (-qJDD(3) * t182 + t190 + (-t259 * t231 - t169) * qJD(3) + t240 * t232) * MDP(16) + (qJDD(3) * t183 + t189 + (t259 * t232 + t167) * qJD(3) + t240 * t231) * MDP(17) + ((-qJD(3) * t161 - qJDD(2) * t183 - t156 + (-qJD(3) * t182 - t167) * qJD(2)) * t232 + (t281 - qJDD(2) * t182 - t155 + (-t169 + t280) * qJD(2)) * t231 + t287) * MDP(18) + (t156 * t183 + t164 * t167 + t155 * t182 + t161 * t169 + t157 * t173 + t163 * t162 - g(1) * (-qJ(5) * t213 + t197) - g(2) * (pkin(4) * t289 + t257) + (-g(1) * (t248 - t221) + g(2) * qJ(5)) * t212) * MDP(19); -MDP(5) * t267 + t285 * t235 * MDP(6) + MDP(7) * t271 + MDP(8) * t270 + qJDD(3) * MDP(9) + (pkin(2) * t288 + t239) * MDP(10) + ((t176 + t210) * qJD(3) + (pkin(2) * t235 + t287) * t232 + t249) * MDP(11) + (0.2e1 * t293 - qJDD(4) + (-t172 * t231 + t175 * t232) * qJD(2) + t239) * MDP(12) - t253 * qJDD(2) * MDP(13) + (-qJD(3) * t176 - t287 * t232 + (t172 * t232 + (t175 - t296) * t231) * qJD(2) - t249 + t301) * MDP(14) + (t159 * qJ(4) - t160 * pkin(3) - t172 * t175 - t171 * t177 - g(1) * (-pkin(3) * t290 + t186) - g(2) * (-pkin(3) * t292 + t184) - g(3) * t286 + t304 * t174) * MDP(15) + (qJ(5) * t271 + qJD(3) * t168 + 0.2e1 * t303 + ((qJ(5) * qJD(3) - t165) * t232 + t274 * t231) * qJD(2) + t241) * MDP(16) + (-qJD(3) * t166 + t193 + (-g(3) + (-t165 - t296) * qJD(2)) * t231 + (-qJ(5) * qJDD(2) - t274 * qJD(2) - t287) * t232 + t266 + t301) * MDP(17) - t247 * qJDD(2) * MDP(18) + (-g(1) * t186 - g(2) * t184 - g(3) * t265 + t156 * qJ(4) - t155 * t299 - t161 * t168 - t163 * t165 + t275 * t164 + t287 * t302) * MDP(19); (-qJD(3) * t174 + t238) * MDP(15) + (-qJDD(3) * pkin(4) - qJ(5) * t263 + t238 - t281) * MDP(19) + t268 * (-t228 * t235 - t234) + t269 * (-qJDD(3) - t267) + ((-MDP(19) * qJ(5) + MDP(13) - MDP(18)) * qJDD(2) + (t172 * MDP(15) - t274 * MDP(19)) * qJD(2)) * t231; (t242 + t261) * MDP(19) - t284 * MDP(18) * t235 + (MDP(16) * t232 + t231 * MDP(17)) * qJDD(2) + (t245 + (0.2e1 * t232 * MDP(17) + (-t299 * MDP(19) - 0.2e1 * MDP(16)) * t231) * qJD(3)) * qJD(2);];
tau = t1;
