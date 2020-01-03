% Calculate vector of inverse dynamics joint torques for
% S5RPPRP2
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
%   pkin=[a2,a3,a4,a5,d1,d4,theta2,theta3]';
% MDP [19x1]
%   Minimal dynamic parameter vector (fixed base model)
%   see S5RPPRP2_convert_par2_MPV_fixb.m
% 
% Output:
% tau [5x1]
%   joint torques of inverse dynamics (contains inertial, gravitational Coriolis and centrifugal forces)

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 17:49
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tau = S5RPPRP2_invdynJ_fixb_mdp_slag_vp(qJ, qJD, qJDD, g, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(5,1),zeros(3,1),zeros(8,1),zeros(19,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPPRP2_invdynJ_fixb_mdp_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RPPRP2_invdynJ_fixb_mdp_slag_vp: qJD has to be [5x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [5 1]), ...
  'S5RPPRP2_invdynJ_fixb_mdp_slag_vp: qJDD has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RPPRP2_invdynJ_fixb_mdp_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RPPRP2_invdynJ_fixb_mdp_slag_vp: pkin has to be [8x1] (double)');
assert(isreal(MDP) && all(size(MDP) == [19 1]), ...
  'S5RPPRP2_invdynJ_fixb_mdp_slag_vp: MDP has to be [19x1] (double)'); 

%% Symbolic Calculation
% From invdyn_joint_fixb_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 17:49:28
% EndTime: 2019-12-31 17:49:31
% DurationCPUTime: 1.46s
% Computational Cost: add. (1148->224), mult. (2248->265), div. (0->0), fcn. (1536->12), ass. (0->111)
t238 = cos(pkin(7));
t296 = pkin(1) * t238;
t221 = -pkin(2) - t296;
t278 = qJDD(1) * t221;
t205 = qJDD(3) + t278;
t234 = qJ(1) + pkin(7);
t227 = sin(t234);
t229 = cos(t234);
t303 = -g(1) * t227 + g(2) * t229;
t304 = -t205 - t303;
t264 = g(1) * t229 + g(2) * t227;
t235 = sin(pkin(8));
t237 = cos(pkin(8));
t240 = sin(qJ(4));
t297 = cos(qJ(4));
t204 = t297 * t235 + t240 * t237;
t195 = t204 * qJD(1);
t302 = MDP(14) + MDP(16);
t301 = -MDP(15) + MDP(18);
t272 = t297 * t237;
t285 = t240 * t235;
t251 = t272 - t285;
t236 = sin(pkin(7));
t217 = pkin(1) * t236 + qJ(3);
t291 = pkin(6) + t217;
t199 = t291 * t235;
t200 = t291 * t237;
t252 = -t297 * t199 - t240 * t200;
t158 = t251 * qJD(3) + t252 * qJD(4);
t167 = -t240 * t199 + t297 * t200;
t233 = pkin(8) + qJ(4);
t226 = sin(t233);
t300 = -qJD(4) * t158 - qJDD(4) * t167 + t226 * t303;
t228 = cos(t233);
t201 = qJD(1) * qJD(3) + qJDD(1) * t217;
t223 = t237 * qJDD(2);
t181 = t223 + (-pkin(6) * qJDD(1) - t201) * t235;
t186 = t235 * qJDD(2) + t237 * t201;
t277 = qJDD(1) * t237;
t182 = pkin(6) * t277 + t186;
t210 = t217 * qJD(1);
t225 = t237 * qJD(2);
t290 = pkin(6) * qJD(1);
t183 = t225 + (-t210 - t290) * t235;
t269 = qJD(4) * t297;
t274 = t240 * t181 + t297 * t182 + t183 * t269;
t299 = -g(3) * t226 - t264 * t228 + t274;
t298 = t195 ^ 2;
t241 = sin(qJ(1));
t293 = g(1) * t241;
t289 = qJDD(4) * pkin(4);
t265 = qJD(1) * t272;
t271 = qJD(1) * t285;
t193 = -t265 + t271;
t288 = t193 * t195;
t287 = t235 * MDP(6);
t188 = t235 * qJD(2) + t237 * t210;
t184 = t237 * t290 + t188;
t286 = t240 * t184;
t198 = t204 * qJD(4);
t268 = qJDD(1) * t297;
t276 = qJDD(1) * t240;
t259 = t235 * t276 - t237 * t268;
t170 = qJD(1) * t198 + t259;
t281 = qJD(4) * t240;
t270 = t235 * t281;
t197 = -t237 * t269 + t270;
t284 = -t204 * t170 + t197 * t193;
t283 = t235 ^ 2 + t237 ^ 2;
t157 = t240 * t183 + t297 * t184;
t282 = qJD(4) * t157;
t156 = t297 * t183 - t286;
t280 = qJD(5) - t156;
t275 = qJDD(4) * qJ(5);
t273 = qJD(4) * t265 + t235 * t268 + t237 * t276;
t220 = pkin(3) * t237 + pkin(2);
t267 = t156 + t286;
t266 = -t297 * t181 + t240 * t182 + t183 * t281 + t184 * t269;
t242 = cos(qJ(1));
t262 = -g(2) * t242 + t293;
t261 = pkin(4) * t228 + qJ(5) * t226;
t169 = qJD(1) * t270 - t273;
t258 = t169 * t251 + t195 * t198;
t185 = -t201 * t235 + t223;
t257 = -t185 * t235 + t186 * t237;
t256 = (-t210 * t235 + t225) * t235 - t188 * t237;
t209 = -t220 - t296;
t172 = -qJD(4) * t198 + qJDD(4) * t251;
t253 = t220 + t261;
t250 = -g(3) * t228 + t264 * t226 - t266;
t191 = t209 * qJD(1) + qJD(3);
t189 = t209 * qJDD(1) + qJDD(3);
t159 = t204 * qJD(3) + t167 * qJD(4);
t248 = -qJD(4) * t159 + qJDD(4) * t252 - t303 * t228;
t162 = pkin(4) * t193 - qJ(5) * t195 + t191;
t247 = t162 * t195 + qJDD(5) - t250;
t246 = pkin(4) * t170 + qJ(5) * t169 + t189;
t239 = -pkin(6) - qJ(3);
t230 = t242 * pkin(1);
t192 = t193 ^ 2;
t171 = -qJD(4) * t197 + qJDD(4) * t204;
t168 = pkin(4) * t195 + qJ(5) * t193;
t165 = -pkin(4) * t251 - qJ(5) * t204 + t209;
t163 = pkin(4) * t198 + qJ(5) * t197 - qJD(5) * t204;
t161 = (t193 - t271) * qJD(4) + t273;
t155 = qJD(4) * qJ(5) + t157;
t154 = -qJD(4) * pkin(4) + t280;
t153 = -qJD(5) * t195 + t246;
t152 = qJDD(5) + t266 - t289;
t151 = t275 + (qJD(5) - t286) * qJD(4) + t274;
t1 = [qJDD(1) * MDP(1) + t262 * MDP(2) + (g(1) * t242 + g(2) * t241) * MDP(3) + (t262 + (t236 ^ 2 + t238 ^ 2) * qJDD(1) * pkin(1)) * pkin(1) * MDP(4) + (t201 * t283 + t257 - t264) * MDP(7) + (t205 * t221 - g(1) * (-pkin(1) * t241 - pkin(2) * t227 + qJ(3) * t229) - g(2) * (pkin(2) * t229 + qJ(3) * t227 + t230) + t257 * t217 - t256 * qJD(3)) * MDP(8) + (-t169 * t204 - t195 * t197) * MDP(9) + (-t258 + t284) * MDP(10) + t171 * MDP(11) + t172 * MDP(12) + (t170 * t209 - t189 * t251 + t191 * t198 + t248) * MDP(14) + (-t169 * t209 + t189 * t204 - t191 * t197 + t300) * MDP(15) + (-t153 * t251 + t162 * t198 + t163 * t193 + t165 * t170 + t248) * MDP(16) + (t151 * t251 + t152 * t204 - t154 * t197 - t155 * t198 - t158 * t193 + t159 * t195 - t167 * t170 + t169 * t252 - t264) * MDP(17) + (-t153 * t204 + t162 * t197 - t163 * t195 + t165 * t169 - t300) * MDP(18) + (pkin(1) * t293 - g(2) * t230 + t151 * t167 - t152 * t252 + t153 * t165 + t154 * t159 + t155 * t158 + t162 * t163 + (g(1) * t239 - g(2) * t253) * t229 + (g(1) * t253 + g(2) * t239) * t227) * MDP(19) + (t237 * MDP(5) - t287) * (-t278 + t304); (qJDD(2) - g(3)) * MDP(4) + (t185 * t237 + t186 * t235 - g(3)) * MDP(8) + (t258 + t284) * MDP(17) + (t151 * t204 - t152 * t251 + t154 * t198 - t155 * t197 - g(3)) * MDP(19) + t301 * t171 + t302 * t172; -MDP(5) * t277 + qJDD(1) * t287 - t283 * MDP(7) * qJD(1) ^ 2 + (t256 * qJD(1) - t304) * MDP(8) + (-t192 - t298) * MDP(17) + (t155 * t193 + (-qJD(5) - t154) * t195 + t246 + t303) * MDP(19) + t302 * (0.2e1 * t195 * qJD(4) + t259) + t301 * ((t193 + t271) * qJD(4) - t273); MDP(9) * t288 + (-t192 + t298) * MDP(10) + t161 * MDP(11) - t259 * MDP(12) + qJDD(4) * MDP(13) + (-t191 * t195 + t250 + t282) * MDP(14) + (t267 * qJD(4) + t191 * t193 - t299) * MDP(15) + (-t168 * t193 - t247 + t282 + 0.2e1 * t289) * MDP(16) + (pkin(4) * t169 - qJ(5) * t170 + (t155 - t157) * t195 + (t154 - t280) * t193) * MDP(17) + (0.2e1 * t275 - t162 * t193 + t168 * t195 + (0.2e1 * qJD(5) - t267) * qJD(4) + t299) * MDP(18) + (-t152 * pkin(4) - g(3) * t261 + t151 * qJ(5) - t154 * t157 + t280 * t155 - t162 * t168 + t264 * (pkin(4) * t226 - qJ(5) * t228)) * MDP(19); (-qJDD(4) + t288) * MDP(16) + t161 * MDP(17) + (-qJD(4) ^ 2 - t298) * MDP(18) + (-qJD(4) * t155 + t247 - t289) * MDP(19);];
tau = t1;
