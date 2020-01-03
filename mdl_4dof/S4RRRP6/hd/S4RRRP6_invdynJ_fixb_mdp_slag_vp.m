% Calculate vector of inverse dynamics joint torques for
% S4RRRP6
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
% pkin [6x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,d1,d2,d3]';
% MDP [19x1]
%   Minimal dynamic parameter vector (fixed base model)
%   see S4RRRP6_convert_par2_MPV_fixb.m
% 
% Output:
% tau [4x1]
%   joint torques of inverse dynamics (contains inertial, gravitational Coriolis and centrifugal forces)

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 17:19
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tau = S4RRRP6_invdynJ_fixb_mdp_slag_vp(qJ, qJD, qJDD, g, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(4,1),zeros(4,1),zeros(3,1),zeros(6,1),zeros(19,1)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4RRRP6_invdynJ_fixb_mdp_slag_vp: qJ has to be [4x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [4 1]), ...
  'S4RRRP6_invdynJ_fixb_mdp_slag_vp: qJD has to be [4x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [4 1]), ...
  'S4RRRP6_invdynJ_fixb_mdp_slag_vp: qJDD has to be [4x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S4RRRP6_invdynJ_fixb_mdp_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [6 1]), ...
  'S4RRRP6_invdynJ_fixb_mdp_slag_vp: pkin has to be [6x1] (double)');
assert(isreal(MDP) && all(size(MDP) == [19 1]), ...
  'S4RRRP6_invdynJ_fixb_mdp_slag_vp: MDP has to be [19x1] (double)'); 

%% Symbolic Calculation
% From invdyn_joint_fixb_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 17:19:15
% EndTime: 2019-12-31 17:19:18
% DurationCPUTime: 1.88s
% Computational Cost: add. (990->264), mult. (2272->362), div. (0->0), fcn. (1380->6), ass. (0->118)
t229 = sin(qJ(2));
t272 = qJD(3) * t229;
t311 = qJD(1) * t272 - qJDD(2);
t228 = sin(qJ(3));
t231 = cos(qJ(3));
t267 = qJDD(1) * t229;
t232 = cos(qJ(2));
t279 = qJD(1) * t232;
t167 = t228 * (qJD(2) * (qJD(3) + t279) + t267) + t311 * t231;
t230 = sin(qJ(1));
t233 = cos(qJ(1));
t253 = g(1) * t233 + g(2) * t230;
t292 = t228 * t232;
t184 = t230 * t292 + t231 * t233;
t289 = t232 * t233;
t186 = -t228 * t289 + t230 * t231;
t310 = -g(1) * t186 + g(2) * t184;
t218 = pkin(5) * t267;
t268 = qJD(1) * qJD(2);
t259 = t232 * t268;
t183 = -qJDD(2) * pkin(2) + pkin(5) * t259 + t218;
t214 = -qJD(3) + t279;
t299 = g(3) * t232;
t237 = -t229 * t253 + t299;
t309 = -qJD(3) * pkin(6) * t214 + t183 + t237;
t276 = qJD(2) * t228;
t280 = qJD(1) * t229;
t200 = t231 * t280 + t276;
t307 = t200 ^ 2;
t306 = pkin(3) * t228;
t305 = pkin(5) * t228;
t300 = g(3) * t229;
t298 = qJ(4) + pkin(6);
t269 = t231 * qJD(2);
t166 = -qJD(3) * t269 + (-t259 - t267) * t231 + t311 * t228;
t297 = t166 * t228;
t263 = t228 * t280;
t198 = t263 - t269;
t296 = t198 * t214;
t295 = t200 * t214;
t294 = t200 * t231;
t293 = t228 * t229;
t291 = t229 * t231;
t290 = t231 * t232;
t204 = -pkin(2) * t232 - pkin(6) * t229 - pkin(1);
t192 = t204 * qJD(1);
t220 = pkin(5) * t279;
t208 = qJD(2) * pkin(6) + t220;
t172 = t231 * t192 - t208 * t228;
t163 = -qJ(4) * t200 + t172;
t162 = -pkin(3) * t214 + t163;
t288 = -t163 + t162;
t254 = pkin(2) * t229 - pkin(6) * t232;
t202 = t254 * qJD(1);
t188 = t228 * t202;
t257 = qJD(3) * t298;
t270 = qJD(4) * t231;
t287 = -t228 * t257 + t270 - t188 - (-pkin(5) * t291 - qJ(4) * t292) * qJD(1);
t246 = pkin(3) * t229 - qJ(4) * t290;
t283 = pkin(5) * t263 + t231 * t202;
t286 = -qJD(1) * t246 - qJD(4) * t228 - t231 * t257 - t283;
t203 = t254 * qJD(2);
t271 = qJD(3) * t231;
t285 = t228 * t203 + t204 * t271;
t275 = qJD(2) * t229;
t284 = t231 * t203 + t275 * t305;
t215 = pkin(5) * t290;
t282 = t228 * t204 + t215;
t225 = t229 ^ 2;
t281 = -t232 ^ 2 + t225;
t278 = qJD(2) * t198;
t277 = qJD(2) * t200;
t274 = qJD(2) * t232;
t273 = qJD(3) * t228;
t266 = t232 * qJDD(1);
t264 = pkin(5) + t306;
t262 = t214 * t269;
t261 = t214 * t273;
t260 = t214 * t271;
t207 = -qJD(2) * pkin(2) + pkin(5) * t280;
t240 = -t229 * t268 + t266;
t182 = pkin(5) * t240 + qJDD(2) * pkin(6);
t256 = -qJD(3) * t192 - t182;
t252 = g(1) * t230 - g(2) * t233;
t174 = qJD(1) * t203 + qJDD(1) * t204;
t170 = t231 * t174;
t251 = t208 * t271 - t170;
t195 = qJDD(3) - t240;
t250 = -pkin(6) * t195 + qJD(3) * t207;
t173 = t192 * t228 + t208 * t231;
t164 = -qJ(4) * t198 + t173;
t249 = t162 * t231 + t164 * t228;
t217 = pkin(3) * t231 + pkin(2);
t248 = t217 * t232 + t229 * t298;
t244 = pkin(1) + t248;
t243 = -0.2e1 * pkin(1) * t268 - pkin(5) * qJDD(2);
t242 = t195 * t228 - t260;
t241 = t195 * t231 + t261;
t239 = t228 * t174 + t231 * t182 + t192 * t271 - t208 * t273;
t235 = qJD(1) ^ 2;
t238 = pkin(1) * t235 + t253;
t161 = pkin(3) * t167 + qJDD(4) + t183;
t234 = qJD(2) ^ 2;
t236 = 0.2e1 * qJDD(1) * pkin(1) - pkin(5) * t234 + t252;
t206 = t298 * t231;
t205 = t298 * t228;
t197 = t231 * t204;
t194 = t198 ^ 2;
t187 = t228 * t230 + t231 * t289;
t185 = t228 * t233 - t230 * t290;
t176 = pkin(3) * t198 + qJD(4) + t207;
t175 = -qJ(4) * t293 + t282;
t171 = -qJ(4) * t291 + t197 + (-pkin(3) - t305) * t232;
t160 = (-pkin(5) * qJD(2) - qJ(4) * qJD(3)) * t291 + (-qJD(4) * t229 + (-pkin(5) * qJD(3) - qJ(4) * qJD(2)) * t232) * t228 + t285;
t159 = -t229 * t270 + t246 * qJD(2) + (-t215 + (qJ(4) * t229 - t204) * t228) * qJD(3) + t284;
t158 = -qJ(4) * t167 - qJD(4) * t198 + t239;
t157 = pkin(3) * t195 + qJ(4) * t166 - t173 * qJD(3) - qJD(4) * t200 - t228 * t182 + t170;
t1 = [qJDD(1) * MDP(1) + t252 * MDP(2) + t253 * MDP(3) + (qJDD(1) * t225 + 0.2e1 * t229 * t259) * MDP(4) + 0.2e1 * (t229 * t266 - t281 * t268) * MDP(5) + (qJDD(2) * t229 + t232 * t234) * MDP(6) + (qJDD(2) * t232 - t229 * t234) * MDP(7) + (t229 * t243 + t232 * t236) * MDP(9) + (-t229 * t236 + t232 * t243) * MDP(10) + (-t166 * t291 + (-t228 * t272 + t232 * t269) * t200) * MDP(11) + ((-t198 * t231 - t200 * t228) * t274 + (t297 - t167 * t231 + (t198 * t228 - t294) * qJD(3)) * t229) * MDP(12) + ((t166 - t262) * t232 + (t241 + t277) * t229) * MDP(13) + ((t214 * t276 + t167) * t232 + (-t242 - t278) * t229) * MDP(14) + (-t195 * t232 - t214 * t275) * MDP(15) + (-(-t204 * t273 + t284) * t214 + t197 * t195 - g(1) * t185 - g(2) * t187 + ((t260 + t278) * pkin(5) + (-pkin(5) * t195 + qJD(2) * t207 - t256) * t228 + t251) * t232 + (pkin(5) * t167 + qJD(2) * t172 + t183 * t228 + t207 * t271) * t229) * MDP(16) + (t285 * t214 - t282 * t195 - g(1) * t184 - g(2) * t186 + (t207 * t269 + (-t261 + t277) * pkin(5) + t239) * t232 + (-t207 * t273 - t173 * qJD(2) + t183 * t231 + (-t166 - t262) * pkin(5)) * t229) * MDP(17) + (-t159 * t200 - t160 * t198 + t166 * t171 - t167 * t175 - t249 * t274 + (-t157 * t231 - t158 * t228 + (t162 * t228 - t164 * t231) * qJD(3) + t252) * t229) * MDP(18) + (t157 * t171 + t158 * t175 + t162 * t159 + t164 * t160 + t176 * t264 * t274 + (t176 * pkin(3) * t271 + t161 * t264) * t229 + (-g(1) * t264 - g(2) * t244) * t233 + (g(1) * t244 - g(2) * t264) * t230) * MDP(19); MDP(6) * t267 + MDP(7) * t266 + qJDD(2) * MDP(8) + (t229 * t238 - t218 - t299) * MDP(9) + (t300 + (-pkin(5) * qJDD(1) + t238) * t232) * MDP(10) + (-t214 * t294 - t297) * MDP(11) + ((-t166 + t296) * t231 + (-t167 + t295) * t228) * MDP(12) + ((-t200 * t229 + t214 * t290) * qJD(1) + t242) * MDP(13) + ((t198 * t229 - t214 * t292) * qJD(1) + t241) * MDP(14) + t214 * MDP(15) * t280 + (-pkin(2) * t167 + t283 * t214 + t250 * t228 + (-t172 * t229 + (-pkin(5) * t198 - t207 * t228) * t232) * qJD(1) - t309 * t231) * MDP(16) + (pkin(2) * t166 - t188 * t214 + t250 * t231 + (-t207 * t290 + t173 * t229 + (-t200 * t232 + t214 * t291) * pkin(5)) * qJD(1) + t309 * t228) * MDP(17) + (-t300 - t157 * t228 + t158 * t231 - t166 * t205 - t167 * t206 - t286 * t200 - t287 * t198 - t249 * qJD(3) + (qJD(1) * t249 - t253) * t232) * MDP(18) + (t158 * t206 - t157 * t205 - t161 * t217 - g(3) * t248 + (-t214 * t306 - t220) * t176 + t287 * t164 + t286 * t162 + t253 * (t217 * t229 - t232 * t298)) * MDP(19) + (-MDP(4) * t229 * t232 + MDP(5) * t281) * t235; t200 * t198 * MDP(11) + (-t194 + t307) * MDP(12) + (-t166 - t296) * MDP(13) + (-t167 - t295) * MDP(14) + t195 * MDP(15) + (-t173 * t214 - t200 * t207 + (t256 + t300) * t228 - t251 + t310) * MDP(16) + (g(1) * t187 - g(2) * t185 + g(3) * t291 - t172 * t214 + t198 * t207 - t239) * MDP(17) + (pkin(3) * t166 - t288 * t198) * MDP(18) + (t288 * t164 + (g(3) * t293 - t176 * t200 + t157 + t310) * pkin(3)) * MDP(19); (-t194 - t307) * MDP(18) + (t162 * t200 + t164 * t198 + t161 + t237) * MDP(19);];
tau = t1;
