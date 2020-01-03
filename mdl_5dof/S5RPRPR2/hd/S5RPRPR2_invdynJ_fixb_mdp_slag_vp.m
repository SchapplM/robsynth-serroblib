% Calculate vector of inverse dynamics joint torques for
% S5RPRPR2
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
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d3,d5,theta2,theta4]';
% MDP [18x1]
%   Minimal dynamic parameter vector (fixed base model)
%   see S5RPRPR2_convert_par2_MPV_fixb.m
% 
% Output:
% tau [5x1]
%   joint torques of inverse dynamics (contains inertial, gravitational Coriolis and centrifugal forces)

% Quelle: HybrDyn-Toolbox
% Datum: 2020-01-03 11:34
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tau = S5RPRPR2_invdynJ_fixb_mdp_slag_vp(qJ, qJD, qJDD, g, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(5,1),zeros(3,1),zeros(9,1),zeros(18,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPRPR2_invdynJ_fixb_mdp_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RPRPR2_invdynJ_fixb_mdp_slag_vp: qJD has to be [5x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [5 1]), ...
  'S5RPRPR2_invdynJ_fixb_mdp_slag_vp: qJDD has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RPRPR2_invdynJ_fixb_mdp_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5RPRPR2_invdynJ_fixb_mdp_slag_vp: pkin has to be [9x1] (double)');
assert(isreal(MDP) && all(size(MDP) == [18 1]), ...
  'S5RPRPR2_invdynJ_fixb_mdp_slag_vp: MDP has to be [18x1] (double)'); 

%% Symbolic Calculation
% From invdyn_joint_fixb_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2020-01-03 11:34:11
% EndTime: 2020-01-03 11:34:16
% DurationCPUTime: 1.22s
% Computational Cost: add. (1090->189), mult. (1784->241), div. (0->0), fcn. (1138->16), ass. (0->109)
t264 = sin(qJ(3));
t267 = cos(qJ(3));
t262 = cos(pkin(8));
t241 = pkin(1) * t262 + pkin(2);
t227 = t241 * qJD(1);
t260 = sin(pkin(8));
t326 = pkin(1) * t260;
t331 = qJD(3) * t227 + qJDD(1) * t326;
t302 = qJD(3) * t326;
t332 = -qJD(1) * t302 + t241 * qJDD(1);
t290 = -t264 * t331 + t267 * t332;
t278 = qJDD(4) - t290;
t253 = qJDD(1) + qJDD(3);
t325 = pkin(3) * t253;
t177 = t278 - t325;
t258 = qJ(1) + pkin(8);
t250 = qJ(3) + t258;
t239 = sin(t250);
t240 = cos(t250);
t289 = -g(2) * t240 - g(3) * t239;
t276 = -t177 + t289;
t259 = sin(pkin(9));
t261 = cos(pkin(9));
t263 = sin(qJ(5));
t266 = cos(qJ(5));
t210 = t259 * t266 + t261 * t263;
t257 = qJD(1) + qJD(3);
t202 = t210 * t257;
t305 = t259 ^ 2 + t261 ^ 2;
t330 = t305 * t257;
t327 = -t264 * t332 - t267 * t331;
t175 = qJ(4) * t253 + qJD(4) * t257 - t327;
t245 = t261 * qJDD(2);
t171 = -t175 * t259 + t245;
t172 = qJDD(2) * t259 + t175 * t261;
t329 = -t171 * t259 + t172 * t261;
t307 = -g(2) * t239 + g(3) * t240;
t306 = t241 * t264 + t267 * t326;
t303 = qJD(1) * t326;
t198 = t227 * t264 + t267 * t303;
t193 = qJ(4) * t257 + t198;
t184 = qJD(2) * t261 - t193 * t259;
t185 = qJD(2) * t259 + t193 * t261;
t328 = -t184 * t259 + t185 * t261;
t197 = t267 * t227 - t264 * t303;
t286 = qJD(4) - t197;
t324 = pkin(4) * t261;
t251 = t261 * pkin(7);
t318 = t198 * t257;
t203 = t306 * qJD(3);
t317 = t203 * t257;
t256 = pkin(9) + qJ(5);
t248 = sin(t256);
t316 = t239 * t248;
t315 = t240 * t248;
t314 = t241 * t267;
t313 = t259 * t263;
t311 = t261 * MDP(8);
t309 = t266 * t261;
t308 = pkin(3) * t240 + qJ(4) * t239;
t300 = t257 * t313;
t299 = t257 * t309;
t298 = t276 * t259;
t297 = qJD(5) * t299 + t210 * t253;
t242 = -pkin(3) - t324;
t296 = t305 * t253;
t295 = pkin(3) * t239 - qJ(4) * t240;
t293 = -t264 * t326 + t314;
t173 = t242 * t253 + t278;
t186 = t242 * t257 + t286;
t209 = -t309 + t313;
t206 = t209 * qJD(5);
t291 = g(2) * t315 + g(3) * t316 + t173 * t210 - t186 * t206;
t205 = -pkin(3) - t293;
t265 = sin(qJ(1));
t268 = cos(qJ(1));
t288 = -g(2) * t268 - g(3) * t265;
t222 = t253 * t309;
t287 = -t253 * t313 + t222;
t285 = -t318 - t325;
t182 = -qJD(5) * t300 + t297;
t207 = t210 * qJD(5);
t183 = t207 * t257 - t287;
t189 = -qJD(5) * t206 + qJDD(5) * t210;
t190 = -qJD(5) * t207 - qJDD(5) * t209;
t200 = -t299 + t300;
t284 = (-t182 * t209 - t183 * t210 + t200 * t206 - t202 * t207) * MDP(13) + (t182 * t210 - t202 * t206) * MDP(12) + t189 * MDP(14) + t190 * MDP(15) + t253 * MDP(5);
t204 = qJ(4) + t306;
t194 = (-pkin(7) - t204) * t259;
t195 = t204 * t261 + t251;
t283 = t194 * t266 - t195 * t263;
t282 = t194 * t263 + t195 * t266;
t281 = t205 * t253 + t317;
t218 = (-pkin(7) - qJ(4)) * t259;
t219 = qJ(4) * t261 + t251;
t280 = t218 * t266 - t219 * t263;
t279 = t218 * t263 + t219 * t266;
t277 = qJD(3) * t314 - t264 * t302;
t275 = t307 + t329;
t272 = t289 + t290;
t249 = cos(t256);
t271 = t173 * t209 + t186 * t207 + t249 * t289;
t270 = -t307 + t327;
t199 = qJD(4) + t277;
t196 = t205 - t324;
t192 = -pkin(3) * t257 + t286;
t167 = t251 * t253 + t172;
t166 = t245 + (-pkin(7) * t253 - t175) * t259;
t1 = [qJDD(1) * MDP(1) + t288 * MDP(2) + (g(2) * t265 - g(3) * t268) * MDP(3) + (t288 + (t260 ^ 2 + t262 ^ 2) * qJDD(1) * pkin(1)) * pkin(1) * MDP(4) + (t253 * t293 + t272 - t317) * MDP(6) + (-t253 * t306 - t257 * t277 + t270) * MDP(7) + (t276 - t281) * t311 + (t259 * t281 - t298) * MDP(9) + (t199 * t330 + t204 * t296 + t275) * MDP(10) + (t177 * t205 + t192 * t203 - g(2) * (pkin(2) * cos(t258) + t268 * pkin(1) + t308) - g(3) * (pkin(2) * sin(t258) + t265 * pkin(1) + t295) + t329 * t204 + t328 * t199) * MDP(11) + (t203 * t200 + t196 * t183 + t283 * qJDD(5) + (-qJD(5) * t282 - t199 * t210) * qJD(5) + t271) * MDP(17) + (t203 * t202 + t196 * t182 - t282 * qJDD(5) + (-qJD(5) * t283 + t199 * t209) * qJD(5) + t291) * MDP(18) + t284; (qJDD(2) - g(1)) * MDP(4) + (t171 * t261 + t172 * t259 - g(1)) * MDP(11) + t190 * MDP(17) - t189 * MDP(18); (t272 + t318) * MDP(6) + (t197 * t257 + t270) * MDP(7) + (t276 - t285) * t311 + (t259 * t285 - t298) * MDP(9) + (qJ(4) * t296 + t286 * t330 + t275) * MDP(10) + (-t177 * pkin(3) - t192 * t198 - g(2) * t308 - g(3) * t295 + (qJ(4) * t172 + t185 * t286) * t261 + (-qJ(4) * t171 - t184 * t286) * t259) * MDP(11) + (t242 * t183 + t280 * qJDD(5) - t198 * t200 + (-t279 * qJD(5) - t210 * t286) * qJD(5) + t271) * MDP(17) + (t242 * t182 - t279 * qJDD(5) - t198 * t202 + (-t280 * qJD(5) + t209 * t286) * qJD(5) + t291) * MDP(18) + t284; (qJDD(4) - t272) * MDP(11) - t222 * MDP(17) + t297 * MDP(18) + (-pkin(3) * MDP(11) - t311 + (MDP(17) * t263 + MDP(9)) * t259) * t253 + (0.2e1 * t202 * MDP(17) + (-t200 - t300) * MDP(18)) * qJD(5) + (-MDP(10) * t330 - MDP(11) * t328) * t257; t202 * t200 * MDP(12) + (-t200 ^ 2 + t202 ^ 2) * MDP(13) + t287 * MDP(15) + qJDD(5) * MDP(16) + (-g(1) * t249 + g(2) * t316 - g(3) * t315 + t266 * t166 - t263 * t167 - t186 * t202) * MDP(17) + (g(1) * t248 - t263 * t166 - t266 * t167 + t186 * t200 - t249 * t307) * MDP(18) + (t297 + (t200 - t300) * qJD(5)) * MDP(14);];
tau = t1;
