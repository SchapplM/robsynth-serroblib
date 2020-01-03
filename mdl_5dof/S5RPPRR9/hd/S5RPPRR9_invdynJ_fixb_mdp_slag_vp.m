% Calculate vector of inverse dynamics joint torques for
% S5RPPRR9
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
%   pkin=[a2,a3,a4,a5,d1,d4,d5,theta3]';
% MDP [23x1]
%   Minimal dynamic parameter vector (fixed base model)
%   see S5RPPRR9_convert_par2_MPV_fixb.m
% 
% Output:
% tau [5x1]
%   joint torques of inverse dynamics (contains inertial, gravitational Coriolis and centrifugal forces)

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 18:03
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tau = S5RPPRR9_invdynJ_fixb_mdp_slag_vp(qJ, qJD, qJDD, g, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(5,1),zeros(3,1),zeros(8,1),zeros(23,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPPRR9_invdynJ_fixb_mdp_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RPPRR9_invdynJ_fixb_mdp_slag_vp: qJD has to be [5x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [5 1]), ...
  'S5RPPRR9_invdynJ_fixb_mdp_slag_vp: qJDD has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RPPRR9_invdynJ_fixb_mdp_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RPPRR9_invdynJ_fixb_mdp_slag_vp: pkin has to be [8x1] (double)');
assert(isreal(MDP) && all(size(MDP) == [23 1]), ...
  'S5RPPRR9_invdynJ_fixb_mdp_slag_vp: MDP has to be [23x1] (double)'); 

%% Symbolic Calculation
% From invdyn_joint_fixb_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 18:02:53
% EndTime: 2019-12-31 18:02:57
% DurationCPUTime: 1.87s
% Computational Cost: add. (1010->279), mult. (1854->374), div. (0->0), fcn. (1173->8), ass. (0->133)
t230 = sin(qJ(4));
t286 = qJD(5) * t230;
t319 = qJD(1) * t286 + qJDD(4);
t227 = cos(pkin(8));
t281 = qJD(1) * qJD(2);
t204 = t227 * t281;
t233 = -pkin(1) - pkin(2);
t198 = qJDD(1) * t233 + qJDD(2);
t226 = sin(pkin(8));
t277 = qJDD(1) * t227;
t301 = qJ(2) * t277 + t226 * t198;
t176 = t204 + t301;
t174 = -qJDD(1) * pkin(6) + t176;
t232 = cos(qJ(4));
t199 = qJD(1) * t233 + qJD(2);
t296 = qJ(2) * qJD(1);
t181 = t226 * t199 + t227 * t296;
t178 = -qJD(1) * pkin(6) + t181;
t172 = qJD(3) * t230 + t178 * t232;
t292 = qJD(4) * t172;
t159 = -qJDD(4) * pkin(4) - qJDD(3) * t232 + t174 * t230 + t292;
t200 = qJD(1) * t232 + qJD(5);
t316 = sin(qJ(1));
t317 = cos(qJ(1));
t185 = t226 * t317 - t227 * t316;
t184 = -t226 * t316 - t227 * t317;
t315 = g(1) * t184;
t254 = g(2) * t185 + t315;
t256 = -pkin(4) * t230 + pkin(7) * t232;
t318 = t200 * (pkin(7) * qJD(5) + t256 * qJD(1)) + t230 * t254 - g(3) * t232 + t159;
t314 = g(1) * t185;
t313 = g(2) * t184;
t312 = pkin(1) * qJDD(1);
t231 = cos(qJ(5));
t280 = qJD(1) * qJD(4);
t267 = t232 * t280;
t276 = qJDD(1) * t230;
t246 = t267 + t276;
t229 = sin(qJ(5));
t279 = qJD(4) * qJD(5);
t271 = t319 * t229 + t231 * t279;
t164 = -t231 * t246 + t271;
t311 = t164 * t229;
t310 = t185 * t232;
t283 = t231 * qJD(4);
t295 = qJD(1) * t230;
t186 = t229 * t295 + t283;
t309 = t186 * t200;
t284 = t229 * qJD(4);
t187 = t231 * t295 - t284;
t308 = t187 * t200;
t307 = t227 * t230;
t268 = t230 * t280;
t273 = t232 * qJDD(1);
t245 = -t268 + t273;
t183 = -qJDD(5) - t245;
t306 = t229 * t183;
t305 = t229 * t232;
t304 = t231 * t200;
t303 = t231 * t232;
t302 = qJDD(3) + g(3);
t193 = t227 * qJ(2) + t226 * t233;
t300 = t317 * pkin(1) + t316 * qJ(2);
t299 = g(1) * t316 - g(2) * t317;
t223 = t230 ^ 2;
t298 = -t232 ^ 2 + t223;
t234 = qJD(4) ^ 2;
t235 = qJD(1) ^ 2;
t297 = t234 + t235;
t294 = qJD(2) * t227;
t171 = qJD(3) * t232 - t178 * t230;
t293 = qJD(4) * t171;
t291 = qJD(4) * t186;
t290 = qJD(4) * t187;
t189 = -pkin(6) + t193;
t289 = qJD(4) * t189;
t288 = qJD(4) * t230;
t287 = qJD(5) * t229;
t285 = qJD(5) * t231;
t282 = qJ(2) * qJDD(1);
t278 = qJDD(1) * t226;
t275 = qJDD(4) * t230;
t274 = qJDD(4) * t232;
t272 = 0.2e1 * t281;
t270 = t200 * t284;
t269 = t200 * t283;
t202 = t226 * t281;
t169 = qJD(4) * pkin(7) + t172;
t263 = t189 * t200 + t169;
t262 = -qJ(2) * t278 + t198 * t227;
t180 = t199 * t227 - t226 * t296;
t192 = -t226 * qJ(2) + t227 * t233;
t177 = qJD(1) * pkin(3) - t180;
t257 = pkin(4) * t232 + pkin(7) * t230;
t170 = qJD(1) * t257 + t177;
t261 = -qJDD(4) * pkin(7) - qJD(5) * t170 - qJDD(3) * t230 - t174 * t232 - t293;
t260 = 0.2e1 * t267;
t259 = qJDD(2) - t312;
t188 = pkin(3) - t192;
t258 = -pkin(1) * t316 + t317 * qJ(2);
t255 = t313 - t314;
t175 = -t202 + t262;
t253 = t164 + t269;
t210 = t229 * t279;
t165 = t229 * t276 + qJDD(4) * t231 - t210 + (t230 * t285 + t232 * t284) * qJD(1);
t252 = -t165 + t270;
t251 = qJD(4) * t178 - t302;
t250 = -qJD(5) * t169 + t313;
t249 = t180 * t226 - t181 * t227;
t173 = qJDD(1) * pkin(3) - t175;
t248 = t285 * t200 - t306;
t247 = t231 * t183 + t200 * t287;
t244 = g(1) * t317 + g(2) * t316;
t243 = -t248 - t291;
t242 = t247 - t290;
t241 = -g(2) * t310 - g(3) * t230 + t261;
t168 = -qJD(4) * pkin(4) - t171;
t240 = pkin(7) * t183 + (t168 + t171) * t200;
t239 = qJD(1) * t177 - qJD(3) * qJD(4) - t174 - t254;
t238 = -qJDD(4) * t189 + (-qJD(1) * t188 - t177 - t294) * qJD(4);
t237 = -t168 * qJD(4) + t189 * t183 - t200 * t294 + t261;
t236 = qJDD(1) * t188 - t189 * t234 + t173 + t202 + t255;
t196 = -t230 * t234 + t274;
t195 = -t232 * t234 - t275;
t182 = qJD(2) * t226 + qJD(4) * t256;
t179 = t188 + t257;
t167 = -t184 * t303 + t185 * t229;
t166 = t184 * t305 + t185 * t231;
t163 = pkin(4) * t245 + pkin(7) * t246 + t173;
t162 = t231 * t163;
t161 = t231 * t169 + t229 * t170;
t160 = -t229 * t169 + t231 * t170;
t1 = [qJDD(1) * MDP(1) + t299 * MDP(2) + t244 * MDP(3) + (-qJDD(2) + t299 + 0.2e1 * t312) * MDP(4) + (-t244 + t272 + 0.2e1 * t282) * MDP(5) + (-t259 * pkin(1) - g(1) * t258 - g(2) * t300 + (t272 + t282) * qJ(2)) * MDP(6) + (-qJDD(1) * t192 + 0.2e1 * t202 + t255 - t262) * MDP(7) + (qJDD(1) * t193 + 0.2e1 * t204 + t254 + t301) * MDP(8) + (t176 * t193 + t175 * t192 - g(1) * (-pkin(2) * t316 + t258) - g(2) * (pkin(2) * t317 + t300) - t249 * qJD(2)) * MDP(9) + (qJDD(1) * t223 + t230 * t260) * MDP(10) + 0.2e1 * (t230 * t273 - t280 * t298) * MDP(11) + t195 * MDP(12) - t196 * MDP(13) + (t230 * t238 + t232 * t236) * MDP(15) + (-t230 * t236 + t232 * t238) * MDP(16) + (-t164 * t230 * t231 + (-t229 * t286 + t232 * t283) * t187) * MDP(17) + ((-t186 * t231 - t187 * t229) * t232 * qJD(4) + (t311 - t165 * t231 + (t186 * t229 - t187 * t231) * qJD(5)) * t230) * MDP(18) + ((t164 - t269) * t232 + (t247 + t290) * t230) * MDP(19) + ((t165 + t270) * t232 + (t248 - t291) * t230) * MDP(20) + (-t183 * t232 - t200 * t288) * MDP(21) + (-g(2) * t167 + (-t186 * t289 + t162) * t232 + (-qJD(4) * t160 - t165 * t189 - t186 * t294) * t230 + (-g(1) * t310 - t179 * t183 + t182 * t200 + (-t168 * t230 - t232 * t263) * qJD(5)) * t231 + ((-qJD(5) * t179 + t189 * t288) * t200 - t159 * t230 - t315 + t237 * t232) * t229) * MDP(22) + (-(t179 * t285 + t182 * t229) * t200 + t179 * t306 - t231 * t315 - g(2) * t166 + (-t187 * t289 + (qJD(5) * t263 - t163 + t314) * t229 + t237 * t231) * t232 + (-t187 * t294 + t168 * t287 - t159 * t231 + t189 * t164 + (t189 * t304 + t161) * qJD(4)) * t230) * MDP(23); -qJDD(1) * MDP(4) - t235 * MDP(5) + (-qJ(2) * t235 + t259 - t299) * MDP(6) + (-t226 * t235 - t277) * MDP(7) + (-t227 * t235 + t278) * MDP(8) + (qJD(1) * t249 + t175 * t227 + t176 * t226 - t299) * MDP(9) + ((0.2e1 * t268 - t273) * t227 + (-t232 * t297 - t275) * t226) * MDP(15) + ((t260 + t276) * t227 + (t230 * t297 - t274) * t226) * MDP(16) + (t247 * t227 + (t230 * t252 + t232 * t243) * t226 + (-(t226 * t231 - t227 * t305) * t200 + t186 * t307) * qJD(1)) * MDP(22) + (t248 * t227 + (t230 * t253 + t232 * t242) * t226 + ((t226 * t229 + t227 * t303) * t200 + t187 * t307) * qJD(1)) * MDP(23); t302 * MDP(9) + t196 * MDP(15) + t195 * MDP(16) + (-MDP(22) * t252 - MDP(23) * t253) * t232 + (MDP(22) * t243 + MDP(23) * t242) * t230; -MDP(12) * t276 - MDP(13) * t273 + qJDD(4) * MDP(14) + (t230 * t239 - t232 * t251 + t292) * MDP(15) + (t230 * t251 + t232 * t239 + t293) * MDP(16) + (-t187 * t304 + t311) * MDP(17) + ((t164 + t309) * t231 + (t165 + t308) * t229) * MDP(18) + ((-t187 * t230 + t200 * t303) * qJD(1) + t248) * MDP(19) + ((t186 * t230 - t200 * t305) * qJD(1) - t247) * MDP(20) + t200 * MDP(21) * t295 + (pkin(4) * t165 + t160 * t295 + t172 * t186 + t240 * t229 - t318 * t231) * MDP(22) + (-pkin(4) * t164 - t161 * t295 + t172 * t187 + t318 * t229 + t240 * t231) * MDP(23) + (-MDP(10) * t230 * t232 + MDP(11) * t298) * t235; t187 * t186 * MDP(17) + (-t186 ^ 2 + t187 ^ 2) * MDP(18) + (t271 - t309) * MDP(19) + (-t210 - t308) * MDP(20) - t183 * MDP(21) + (-g(1) * t166 + t161 * t200 + t168 * t187 + t162) * MDP(22) + (g(1) * t167 + t160 * t200 - t168 * t186) * MDP(23) + (-t246 * MDP(19) + t319 * MDP(20) + t250 * MDP(22) + t241 * MDP(23)) * t231 + (t246 * MDP(20) + t241 * MDP(22) + (-t163 - t250) * MDP(23)) * t229;];
tau = t1;
