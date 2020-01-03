% Calculate vector of inverse dynamics joint torques for
% S5RPRPP5
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
% pkin [6x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d3]';
% MDP [21x1]
%   Minimal dynamic parameter vector (fixed base model)
%   see S5RPRPP5_convert_par2_MPV_fixb.m
% 
% Output:
% tau [5x1]
%   joint torques of inverse dynamics (contains inertial, gravitational Coriolis and centrifugal forces)

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 18:16
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tau = S5RPRPP5_invdynJ_fixb_mdp_slag_vp(qJ, qJD, qJDD, g, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(5,1),zeros(3,1),zeros(6,1),zeros(21,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPRPP5_invdynJ_fixb_mdp_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RPRPP5_invdynJ_fixb_mdp_slag_vp: qJD has to be [5x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [5 1]), ...
  'S5RPRPP5_invdynJ_fixb_mdp_slag_vp: qJDD has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RPRPP5_invdynJ_fixb_mdp_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [6 1]), ...
  'S5RPRPP5_invdynJ_fixb_mdp_slag_vp: pkin has to be [6x1] (double)');
assert(isreal(MDP) && all(size(MDP) == [21 1]), ...
  'S5RPRPP5_invdynJ_fixb_mdp_slag_vp: MDP has to be [21x1] (double)'); 

%% Symbolic Calculation
% From invdyn_joint_fixb_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 18:16:47
% EndTime: 2019-12-31 18:16:49
% DurationCPUTime: 1.79s
% Computational Cost: add. (824->264), mult. (1366->302), div. (0->0), fcn. (601->4), ass. (0->129)
t241 = -pkin(3) - pkin(4);
t242 = -pkin(1) - pkin(6);
t200 = qJD(1) * t242 + qJD(2);
t239 = cos(qJ(3));
t298 = qJ(5) * qJD(1);
t180 = (t200 + t298) * t239;
t294 = qJD(4) - t180;
t172 = qJD(3) * t241 + t294;
t296 = qJD(3) * t239;
t189 = t200 * t296;
t323 = qJDD(1) * t242;
t199 = qJDD(2) + t323;
t237 = sin(qJ(3));
t191 = t237 * t199;
t232 = qJDD(3) * qJ(4);
t233 = qJD(3) * qJD(4);
t168 = t189 + t191 + t232 + t233;
t297 = qJD(3) * t237;
t188 = t200 * t297;
t192 = t239 * t199;
t276 = t192 - qJDD(4) - t188;
t315 = qJDD(3) * pkin(3);
t169 = -t276 - t315;
t270 = -t200 * t239 + qJD(4);
t319 = qJD(3) * pkin(3);
t182 = t270 - t319;
t193 = t237 * t200;
t234 = qJD(3) * qJ(4);
t184 = t193 + t234;
t245 = t168 * t237 - t169 * t239 + (t182 * t237 + t184 * t239) * qJD(3);
t240 = cos(qJ(1));
t229 = g(2) * t240;
t238 = sin(qJ(1));
t230 = g(1) * t238;
t303 = t230 - t229;
t325 = t245 - t303;
t308 = qJ(5) + t242;
t302 = g(1) * t240 + g(2) * t238;
t322 = qJDD(3) * t241;
t292 = qJD(1) * qJD(3);
t278 = t239 * t292;
t289 = qJDD(1) * t237;
t291 = qJD(1) * qJD(5);
t321 = t237 * t291 + (t278 + t289) * qJ(5);
t316 = qJ(4) * t237;
t255 = t239 * t241 - t316;
t320 = pkin(3) * t237;
t318 = pkin(1) * qJDD(1);
t244 = qJD(1) ^ 2;
t317 = qJ(2) * t244;
t314 = t237 * t238;
t313 = t237 * t240;
t312 = t238 * t239;
t221 = t239 * qJ(4);
t311 = t239 * t240;
t310 = t239 * t244;
t243 = qJD(3) ^ 2;
t309 = t242 * t243;
t307 = pkin(3) * t312 + qJ(4) * t314;
t306 = g(1) * t311 + g(2) * t312;
t288 = qJDD(1) * t239;
t295 = qJD(4) * t239;
t305 = -qJ(4) * t288 - qJD(1) * t295;
t304 = t240 * pkin(1) + t238 * qJ(2);
t235 = t237 ^ 2;
t236 = t239 ^ 2;
t301 = t235 - t236;
t300 = t235 + t236;
t299 = t243 + t244;
t281 = t237 * t241;
t265 = -qJ(2) + t281;
t174 = qJD(5) + (t265 + t221) * qJD(1);
t293 = qJD(5) + t174;
t290 = qJDD(1) * qJ(2);
t287 = qJDD(3) * t237;
t286 = qJDD(3) * t242;
t285 = MDP(14) + MDP(18);
t284 = -MDP(15) + MDP(20);
t283 = MDP(16) + MDP(19);
t282 = 0.2e1 * qJD(1) * qJD(2);
t280 = t237 * t310;
t279 = t241 * MDP(21);
t218 = t237 * t298;
t277 = qJDD(5) - t305;
t275 = qJ(2) + t320;
t274 = t221 - t320;
t273 = qJD(3) * t308;
t264 = t221 + t281;
t187 = -qJ(2) + t264;
t271 = -qJD(1) * t187 - t174;
t269 = qJDD(2) - t318;
t268 = -g(1) * t312 + g(2) * t311 + g(3) * t237 + t192;
t267 = pkin(3) * t314 + t240 * pkin(6) + t304;
t266 = -g(2) * t313 + g(3) * t239 - t191;
t262 = pkin(3) * t239 + t316;
t259 = -qJDD(4) + t268;
t183 = (t275 - t221) * qJD(1);
t196 = qJ(2) - t274;
t258 = (qJD(1) * t196 + t183) * qJD(3);
t222 = t240 * qJ(2);
t257 = pkin(3) * t313 - qJ(4) * t311 + t222;
t256 = t282 + 0.2e1 * t290;
t254 = -t188 + t259;
t247 = qJD(3) * t255 - qJD(2);
t164 = qJD(1) * t247 + qJDD(1) * t265 + t277;
t173 = t247 + t295;
t253 = qJD(1) * t173 + qJDD(1) * t187 + t164;
t202 = qJD(3) * t218;
t165 = t202 + (-qJ(5) * qJDD(1) - t291) * t239 + t322 - t276;
t166 = t168 + t321;
t175 = t218 + t184;
t252 = -t165 * t239 + t166 * t237 + t172 * t297 + t175 * t296 - t303;
t251 = qJD(3) * t262 + qJD(2);
t250 = t256 - t309;
t249 = -t254 - t315;
t167 = qJD(1) * t251 + qJDD(1) * t275 + t305;
t178 = t251 - t295;
t248 = -qJD(1) * t178 - qJDD(1) * t196 - t167 + t309;
t246 = -g(1) * t314 + 0.2e1 * t232 + 0.2e1 * t233 - t266;
t219 = qJDD(3) * t239;
t204 = t239 * t286;
t195 = t308 * t239;
t194 = t308 * t237;
t190 = t262 * qJD(1);
t181 = t255 * qJD(1);
t179 = t193 + t218;
t177 = qJD(5) * t237 + t239 * t273;
t176 = -qJD(5) * t239 + t237 * t273;
t1 = [qJDD(1) * MDP(1) + t303 * MDP(2) + t302 * MDP(3) + (qJDD(2) - t303 - 0.2e1 * t318) * MDP(4) + (t256 - t302) * MDP(5) + (-t269 * pkin(1) - g(1) * (-pkin(1) * t238 + t222) - g(2) * t304 + (t282 + t290) * qJ(2)) * MDP(6) + (qJDD(1) * t236 - 0.2e1 * t237 * t278) * MDP(7) + 0.2e1 * (-t237 * t288 + t292 * t301) * MDP(8) + (-t237 * t243 + t219) * MDP(9) + (-t239 * t243 - t287) * MDP(10) + (0.2e1 * qJ(2) * t278 + t204 + (t250 - t302) * t237) * MDP(12) + ((-0.2e1 * qJ(2) * t292 - t286) * t237 + t250 * t239 - t306) * MDP(13) + (t204 + t239 * t258 + (-t248 - t302) * t237) * MDP(14) + (-t300 * t323 - t325) * MDP(15) + ((t258 + t286) * t237 + t248 * t239 + t306) * MDP(16) + (t167 * t196 + t183 * t178 - g(1) * (t238 * t242 + t257) - g(2) * (-t221 * t238 + t267) + t245 * t242) * MDP(17) + (qJDD(3) * t195 + (t239 * t271 - t176) * qJD(3) + (-t253 - t302) * t237) * MDP(18) + (qJDD(3) * t194 + t253 * t239 + (t237 * t271 + t177) * qJD(3) + t306) * MDP(19) + ((t194 * t237 + t195 * t239) * qJDD(1) + (-t176 * t239 + t177 * t237 + (t194 * t239 - t195 * t237) * qJD(3)) * qJD(1) + t252) * MDP(20) + (t166 * t194 + t175 * t177 - t165 * t195 + t172 * t176 + t164 * t187 + t174 * t173 - g(1) * (pkin(4) * t313 + t257) - g(2) * (-qJ(5) * t240 + t267) + (-g(1) * t308 - g(2) * (pkin(4) * t237 - t221)) * t238) * MDP(21); -t244 * MDP(5) + (t269 - t303 - t317) * MDP(6) + (-qJD(1) * t183 + t325) * MDP(17) + (qJD(1) * t174 + t252) * MDP(21) + (MDP(12) + t285) * (-t237 * t299 + t219) + (-MDP(13) + t283) * (t239 * t299 + t287) + (t284 * t300 + MDP(4)) * qJDD(1); MDP(7) * t280 - t301 * t244 * MDP(8) + MDP(9) * t288 - MDP(10) * t289 + qJDD(3) * MDP(11) + (-qJ(2) * t310 + t268) * MDP(12) + ((t317 + t230) * t237 + t266) * MDP(13) + (0.2e1 * t315 + (-t183 * t239 - t190 * t237) * qJD(1) + t259) * MDP(14) + (-t262 * qJDD(1) + ((t184 - t234) * t239 + (-qJD(4) + t182 + t319) * t237) * qJD(1)) * MDP(15) + ((-t183 * t237 + t190 * t239) * qJD(1) + t246) * MDP(16) + (-t169 * pkin(3) - g(1) * t307 - g(3) * t274 + t168 * qJ(4) - t182 * t193 - t183 * t190 + t184 * t270 + t229 * t262) * MDP(17) + (qJ(5) * t288 + qJD(3) * t179 - t202 - 0.2e1 * t322 + (t181 * t237 + t239 * t293) * qJD(1) + t254) * MDP(18) + (-qJD(3) * t180 + t189 + (t174 * t237 - t181 * t239) * qJD(1) + t246 + t321) * MDP(19) + (-t255 * qJDD(1) + (-t175 + t179 + t234) * t239 * qJD(1)) * MDP(20) + (t166 * qJ(4) + t165 * t241 - t172 * t179 - t174 * t181 - g(1) * (pkin(4) * t312 + t307) - g(3) * t264 + t294 * t175 - t255 * t229) * MDP(21); (-qJD(3) * t184 + t249) * MDP(17) + (-qJDD(3) * pkin(4) - qJD(3) * t175 + t202 + t249) * MDP(21) + t283 * (-t236 * t244 - t243) + t285 * (-qJDD(3) + t280) + ((-MDP(21) * qJ(5) - t284) * qJDD(1) + (t183 * MDP(17) - MDP(21) * t293) * qJD(1)) * t239; (t277 + t302) * MDP(21) - t300 * MDP(20) * t244 + (t239 * MDP(19) - qJ(2) * MDP(21) + (-MDP(18) + t279) * t237) * qJDD(1) + ((t172 * t239 - t175 * t237 - qJD(2)) * MDP(21) + ((-MDP(21) * qJ(4) - 0.2e1 * MDP(19)) * t237 + (-0.2e1 * MDP(18) + t279) * t239) * qJD(3)) * qJD(1);];
tau = t1;
