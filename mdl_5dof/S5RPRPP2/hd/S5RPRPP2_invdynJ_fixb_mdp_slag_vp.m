% Calculate vector of inverse dynamics joint torques for
% S5RPRPP2
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
%   pkin=[a2,a3,a4,a5,d1,d3,theta2]';
% MDP [19x1]
%   Minimal dynamic parameter vector (fixed base model)
%   see S5RPRPP2_convert_par2_MPV_fixb.m
% 
% Output:
% tau [5x1]
%   joint torques of inverse dynamics (contains inertial, gravitational Coriolis and centrifugal forces)

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 18:11
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tau = S5RPRPP2_invdynJ_fixb_mdp_slag_vp(qJ, qJD, qJDD, g, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(5,1),zeros(3,1),zeros(7,1),zeros(19,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPRPP2_invdynJ_fixb_mdp_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RPRPP2_invdynJ_fixb_mdp_slag_vp: qJD has to be [5x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [5 1]), ...
  'S5RPRPP2_invdynJ_fixb_mdp_slag_vp: qJDD has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RPRPP2_invdynJ_fixb_mdp_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [7 1]), ...
  'S5RPRPP2_invdynJ_fixb_mdp_slag_vp: pkin has to be [7x1] (double)');
assert(isreal(MDP) && all(size(MDP) == [19 1]), ...
  'S5RPRPP2_invdynJ_fixb_mdp_slag_vp: MDP has to be [19x1] (double)'); 

%% Symbolic Calculation
% From invdyn_joint_fixb_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 18:11:16
% EndTime: 2019-12-31 18:11:19
% DurationCPUTime: 1.81s
% Computational Cost: add. (911->255), mult. (1636->299), div. (0->0), fcn. (835->8), ass. (0->126)
t251 = cos(qJ(3));
t316 = qJ(4) * t251;
t249 = sin(qJ(3));
t322 = pkin(3) + pkin(4);
t324 = t322 * t249;
t268 = t316 - t324;
t248 = cos(pkin(7));
t321 = pkin(1) * t248;
t287 = pkin(2) + t321;
t200 = qJDD(1) * t287;
t233 = t249 * qJ(4);
t283 = pkin(2) + t233;
t247 = sin(pkin(7));
t218 = pkin(1) * t247 + pkin(6);
t199 = t218 * qJDD(1);
t328 = qJD(2) * qJD(3) + t199;
t297 = qJD(1) * qJD(3);
t285 = t251 * t297;
t294 = qJDD(1) * t249;
t327 = t285 + t294;
t202 = t287 * qJD(1);
t201 = t218 * qJD(1);
t194 = t249 * t201;
t184 = t251 * qJD(2) - t194;
t326 = qJD(4) - t184;
t325 = t322 * qJDD(3);
t306 = qJD(1) * t249;
t175 = qJ(5) * t306 + t184;
t300 = qJD(4) - t175;
t172 = -t322 * qJD(3) + t300;
t185 = t249 * qJD(2) + t251 * t201;
t176 = -qJ(5) * qJD(1) * t251 + t185;
t243 = qJD(3) * qJ(4);
t174 = t176 + t243;
t181 = t243 + t185;
t177 = -qJD(3) * pkin(3) + t326;
t240 = qJ(1) + pkin(7);
t227 = sin(t240);
t228 = cos(t240);
t277 = g(1) * t228 + g(2) * t227;
t237 = t251 * pkin(3);
t270 = -t283 - t237;
t182 = (t270 - t321) * qJD(1);
t309 = t237 + t233;
t188 = -t287 - t309;
t292 = qJDD(3) * t218;
t323 = (qJD(1) * t188 + t182) * qJD(3) - t292;
t320 = g(1) * t227;
t317 = g(2) * t228;
t236 = t251 * pkin(4);
t315 = qJDD(3) * pkin(3);
t314 = t227 * t249;
t313 = t227 * t251;
t312 = t228 * t249;
t311 = t228 * t251;
t310 = qJ(5) - t218;
t244 = t249 ^ 2;
t245 = t251 ^ 2;
t308 = t244 - t245;
t307 = t244 + t245;
t305 = qJD(3) * t174;
t190 = t310 * t251;
t304 = qJD(3) * t190;
t303 = qJD(3) * t201;
t302 = qJD(3) * t249;
t301 = qJD(4) * t249;
t173 = qJD(5) + (t322 * t251 + t283 + t321) * qJD(1);
t299 = qJD(5) + t173;
t298 = qJ(5) * qJDD(1);
t296 = qJD(1) * qJD(5);
t293 = qJDD(1) * t251;
t291 = MDP(12) + MDP(16);
t290 = MDP(14) + MDP(17);
t255 = qJD(1) ^ 2;
t289 = t249 * t251 * t255;
t288 = t249 * qJDD(2) + t328 * t251;
t286 = t249 * t297;
t250 = sin(qJ(1));
t284 = -pkin(1) * t250 + t228 * pkin(6);
t282 = -g(3) - t303;
t183 = -t188 + t236;
t280 = qJD(1) * t183 + t173;
t278 = (-qJDD(2) + t303) * t251 + t328 * t249;
t252 = cos(qJ(1));
t276 = g(1) * t250 - g(2) * t252;
t275 = pkin(3) * t249 - t316;
t254 = qJD(3) ^ 2;
t274 = t218 * t254 + t317;
t241 = qJDD(3) * qJ(4);
t242 = qJD(3) * qJD(4);
t273 = 0.2e1 * t241 + 0.2e1 * t242 + t288;
t272 = t252 * pkin(1) + pkin(3) * t311 + t227 * pkin(6) + t283 * t228;
t271 = -qJDD(4) - t278;
t269 = pkin(3) * t293 + t327 * qJ(4) + qJD(1) * t301 + t200;
t267 = (t172 * t249 + t174 * t251) * MDP(19);
t266 = g(1) * t312 + g(2) * t314 - g(3) * t251 - t278;
t265 = -qJDD(4) + t266;
t264 = pkin(4) * t293 + qJDD(5) + t269;
t263 = -0.2e1 * t200 + t274;
t169 = -t201 * t302 + t241 + t242 + t288;
t262 = -0.2e1 * t202 * qJD(3) - t292;
t168 = -t322 * t286 + t264;
t178 = t268 * qJD(3) + t301;
t261 = qJD(1) * t178 + qJDD(1) * t183 + t168 - t317;
t260 = qJD(3) * t185 + t266;
t259 = -t265 - t315;
t171 = pkin(3) * t286 - t269;
t187 = t275 * qJD(3) - t301;
t258 = -qJD(1) * t187 - qJDD(1) * t188 - t171 - t274;
t170 = -t271 - t315;
t257 = t169 * t251 + t170 * t249 + (t177 * t251 - t181 * t249) * qJD(3);
t213 = qJ(5) * t286;
t210 = g(1) * t313;
t209 = g(1) * t314;
t206 = qJ(4) * t311;
t204 = qJ(4) * t313;
t198 = qJDD(3) * t251 - t249 * t254;
t197 = qJDD(3) * t249 + t251 * t254;
t196 = t275 * qJD(1);
t189 = t310 * t249;
t186 = t268 * qJD(1);
t180 = -qJD(5) * t249 - t304;
t179 = -qJD(5) * t251 + t310 * t302;
t167 = t213 + (-t296 - t298) * t251 + t169;
t166 = -t327 * qJ(5) - t249 * t296 - t271 - t325;
t1 = [qJDD(1) * MDP(1) + t276 * MDP(2) + (g(1) * t252 + g(2) * t250) * MDP(3) + (t276 + (t247 ^ 2 + t248 ^ 2) * qJDD(1) * pkin(1)) * pkin(1) * MDP(4) + (qJDD(1) * t244 + 0.2e1 * t249 * t285) * MDP(5) + 0.2e1 * (t249 * t293 - t308 * t297) * MDP(6) + t197 * MDP(7) + t198 * MDP(8) + (t262 * t249 - t263 * t251 + t210) * MDP(10) + (t263 * t249 + t262 * t251 - t209) * MDP(11) + (t323 * t249 + t258 * t251 + t210) * MDP(12) + (t307 * t199 + t257 - t277) * MDP(13) + (t258 * t249 - t323 * t251 + t209) * MDP(14) + (-g(1) * t284 - g(2) * t272 + t171 * t188 + t182 * t187 + t257 * t218 - t270 * t320) * MDP(15) + (qJDD(3) * t189 + t210 + (-t280 * t249 - t180) * qJD(3) + t261 * t251) * MDP(16) + (-qJDD(3) * t190 + t209 + (t280 * t251 + t179) * qJD(3) + t261 * t249) * MDP(17) + ((-qJD(3) * t172 + qJDD(1) * t190 - t167 + (qJD(3) * t189 - t179) * qJD(1)) * t251 + (t305 + qJDD(1) * t189 - t166 + (-t180 - t304) * qJD(1)) * t249 + t277) * MDP(18) + (-t167 * t190 + t174 * t179 - t166 * t189 + t172 * t180 + t168 * t183 + t173 * t178 - g(1) * (-qJ(5) * t228 + t284) - g(2) * (pkin(4) * t311 + t272) + (-g(1) * (t270 - t236) + g(2) * qJ(5)) * t227) * MDP(19); (qJDD(2) - g(3)) * MDP(4) + (t169 * t249 - t170 * t251 - g(3)) * MDP(15) + (-t166 * t251 + t167 * t249 - g(3)) * MDP(19) + ((t177 * t249 + t181 * t251) * MDP(15) + t267) * qJD(3) + (MDP(10) + t291) * t198 + (-MDP(11) + t290) * t197; -MDP(5) * t289 + t308 * MDP(6) * t255 + MDP(7) * t294 + MDP(8) * t293 + qJDD(3) * MDP(9) + (t202 * t306 + t260) * MDP(10) + (g(3) * t249 + (t184 + t194) * qJD(3) + (qJD(1) * t202 + t277) * t251 - t288) * MDP(11) + (0.2e1 * t315 - qJDD(4) + (-t182 * t249 + t196 * t251) * qJD(1) + t260) * MDP(12) + (-qJD(3) * t184 + (qJD(1) * t182 - t277) * t251 + (qJD(1) * t196 + t282) * t249 + t273) * MDP(14) + (t169 * qJ(4) - t170 * pkin(3) - t182 * t196 - t177 * t185 - g(1) * (-pkin(3) * t312 + t206) - g(2) * (-pkin(3) * t314 + t204) - g(3) * t309 + t326 * t181) * MDP(15) + (qJ(5) * t294 + qJD(3) * t176 + 0.2e1 * t325 + ((qJ(5) * qJD(3) - t186) * t251 + t299 * t249) * qJD(1) + t265) * MDP(16) + (-qJD(3) * t175 + t213 + (-qJD(1) * t186 + t282) * t249 + (-t299 * qJD(1) - t277 - t298) * t251 + t273) * MDP(17) + (t167 * qJ(4) - t166 * t322 - t172 * t176 - t173 * t186 - g(1) * t206 - g(2) * t204 - g(3) * (t236 + t309) + t300 * t174 + t277 * t324) * MDP(19) + (-t275 * MDP(13) - t268 * MDP(18)) * qJDD(1); (-qJD(3) * t181 + t259) * MDP(15) + (-qJDD(3) * pkin(4) - qJ(5) * t285 + t259 - t305) * MDP(19) + t290 * (-t244 * t255 - t254) + t291 * (-qJDD(3) - t289) + ((-MDP(19) * qJ(5) + MDP(13) - MDP(18)) * qJDD(1) + (t182 * MDP(15) - t299 * MDP(19)) * qJD(1)) * t249; (t264 - t317 + t320) * MDP(19) - t307 * MDP(18) * t255 + (MDP(16) * t251 + t249 * MDP(17)) * qJDD(1) + (t267 + (0.2e1 * t251 * MDP(17) + (-t322 * MDP(19) - 0.2e1 * MDP(16)) * t249) * qJD(3)) * qJD(1);];
tau = t1;
