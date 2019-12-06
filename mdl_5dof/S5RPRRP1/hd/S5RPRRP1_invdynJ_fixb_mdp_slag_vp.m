% Calculate vector of inverse dynamics joint torques for
% S5RPRRP1
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
%   pkin=[a2,a3,a4,a5,d1,d3,d4]';
% MDP [22x1]
%   Minimal dynamic parameter vector (fixed base model)
%   see S5RPRRP1_convert_par2_MPV_fixb.m
% 
% Output:
% tau [5x1]
%   joint torques of inverse dynamics (contains inertial, gravitational Coriolis and centrifugal forces)

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-05 18:00
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tau = S5RPRRP1_invdynJ_fixb_mdp_slag_vp(qJ, qJD, qJDD, g, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(5,1),zeros(3,1),zeros(7,1),zeros(22,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPRRP1_invdynJ_fixb_mdp_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RPRRP1_invdynJ_fixb_mdp_slag_vp: qJD has to be [5x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [5 1]), ...
  'S5RPRRP1_invdynJ_fixb_mdp_slag_vp: qJDD has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RPRRP1_invdynJ_fixb_mdp_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [7 1]), ...
  'S5RPRRP1_invdynJ_fixb_mdp_slag_vp: pkin has to be [7x1] (double)');
assert(isreal(MDP) && all(size(MDP) == [22 1]), ...
  'S5RPRRP1_invdynJ_fixb_mdp_slag_vp: MDP has to be [22x1] (double)'); 

%% Symbolic Calculation
% From invdyn_joint_fixb_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 18:00:04
% EndTime: 2019-12-05 18:00:09
% DurationCPUTime: 1.77s
% Computational Cost: add. (1374->240), mult. (2592->299), div. (0->0), fcn. (1587->8), ass. (0->117)
t270 = qJD(1) ^ 2;
t264 = sin(qJ(1));
t267 = cos(qJ(1));
t312 = g(1) * t264 - g(2) * t267;
t279 = -qJ(2) * t270 - t312;
t261 = qJ(3) + qJ(4);
t245 = sin(t261);
t246 = cos(t261);
t330 = g(3) * t245 - t312 * t246;
t266 = cos(qJ(3));
t303 = qJD(1) * qJD(3);
t296 = t266 * t303;
t263 = sin(qJ(3));
t302 = qJDD(1) * t263;
t329 = t296 + t302;
t301 = qJDD(1) * t266;
t328 = t263 * t303 - t301;
t265 = cos(qJ(4));
t262 = sin(qJ(4));
t309 = qJD(1) * t263;
t298 = t262 * t309;
t308 = qJD(1) * t266;
t206 = t265 * t308 - t298;
t198 = t206 * qJ(5);
t268 = -pkin(1) - pkin(6);
t228 = t268 * qJD(1) + qJD(2);
t196 = -pkin(7) * t309 + t228 * t263;
t191 = t262 * t196;
t197 = -pkin(7) * t308 + t266 * t228;
t194 = qJD(3) * pkin(3) + t197;
t294 = t265 * t194 - t191;
t327 = t198 - t294;
t322 = pkin(7) - t268;
t220 = t322 * t263;
t221 = t322 * t266;
t314 = -t265 * t220 - t262 * t221;
t214 = t262 * t266 + t263 * t265;
t204 = t214 * qJD(1);
t255 = qJD(3) + qJD(4);
t256 = qJDD(1) * qJ(2);
t288 = g(1) * t267 + g(2) * t264;
t257 = qJD(1) * qJD(2);
t299 = 0.2e1 * t257;
t326 = 0.2e1 * t256 + t299 - t288;
t227 = t268 * qJDD(1) + qJDD(2);
t216 = t266 * t227;
t307 = qJD(3) * t263;
t183 = qJDD(3) * pkin(3) + t328 * pkin(7) - t228 * t307 + t216;
t306 = qJD(3) * t266;
t184 = -t329 * pkin(7) + t227 * t263 + t228 * t306;
t305 = qJD(4) * t262;
t325 = (qJD(4) * t194 + t184) * t265 + t262 * t183 - t196 * t305;
t324 = t206 ^ 2;
t321 = pkin(1) * qJDD(1);
t319 = qJ(5) * t204;
t222 = pkin(3) * t309 + qJD(1) * qJ(2);
t188 = pkin(4) * t204 + qJD(5) + t222;
t318 = t188 * t206;
t192 = t265 * t196;
t249 = t263 * pkin(3);
t235 = qJ(2) + t249;
t172 = pkin(4) * t255 - t327;
t317 = t172 + t327;
t186 = t255 * t214;
t215 = -t262 * t263 + t265 * t266;
t253 = qJDD(3) + qJDD(4);
t316 = -t186 * t255 + t215 * t253;
t315 = t265 * t197 - t191;
t313 = t267 * pkin(1) + t264 * qJ(2);
t260 = t266 ^ 2;
t311 = t263 ^ 2 - t260;
t269 = qJD(3) ^ 2;
t310 = -t269 - t270;
t304 = qJD(4) * t265;
t229 = pkin(3) * t306 + qJD(2);
t300 = qJDD(3) * t263;
t293 = -t197 * t262 - t192;
t292 = t220 * t262 - t265 * t221;
t290 = t255 * t266;
t289 = qJDD(2) - t321;
t195 = t329 * pkin(3) + t256 + t257;
t286 = t262 * t302 - t265 * t301;
t177 = t186 * qJD(1) + t286;
t285 = -t177 * t215 - t186 * t206;
t187 = -t262 * t307 - t263 * t305 + t265 * t290;
t284 = -t187 * t255 - t214 * t253;
t283 = -t194 * t262 - t192;
t282 = -qJD(4) * t298 - t328 * t262;
t281 = 0.2e1 * qJ(2) * t303 + qJDD(3) * t268;
t212 = t322 * t307;
t213 = qJD(3) * t221;
t280 = t262 * t212 - t265 * t213 + t220 * t305 - t221 * t304;
t178 = (qJD(1) * t290 + t302) * t265 + t282;
t278 = pkin(4) * t178 + qJDD(5) + t195;
t203 = t204 ^ 2;
t277 = t206 * t204 * MDP(14) - t286 * MDP(16) + (t206 * t255 + (-t255 * t308 - t302) * t265 - t282) * MDP(17) + (-t203 + t324) * MDP(15) + t253 * MDP(18);
t276 = t283 * qJD(4) + t265 * t183 - t262 * t184;
t275 = -qJD(4) * t314 + t265 * t212 + t213 * t262;
t166 = pkin(4) * t253 + qJ(5) * t177 - qJD(5) * t206 + t276;
t167 = -qJ(5) * t178 - qJD(5) * t204 + t325;
t174 = -t283 - t319;
t274 = t166 * t215 + t167 * t214 - t172 * t186 + t174 * t187 - t312;
t273 = -t268 * t269 + t326;
t272 = g(3) * t246 + t222 * t204 + t312 * t245 - t325;
t271 = -t222 * t206 + t276 + t330;
t254 = -qJ(5) - pkin(7) - pkin(6);
t248 = t267 * qJ(2);
t244 = qJDD(3) * t266;
t239 = pkin(3) * t265 + pkin(4);
t218 = pkin(4) * t245 + t249;
t182 = -qJ(5) * t214 + t314;
t181 = -qJ(5) * t215 + t292;
t176 = -t198 + t315;
t175 = t293 + t319;
t169 = qJ(5) * t186 - qJD(5) * t215 + t275;
t168 = -qJ(5) * t187 - qJD(5) * t214 + t280;
t1 = [qJDD(1) * MDP(1) + t312 * MDP(2) + t288 * MDP(3) + (qJDD(2) - t312 - 0.2e1 * t321) * MDP(4) + t326 * MDP(5) + (-t289 * pkin(1) - g(1) * (-pkin(1) * t264 + t248) - g(2) * t313 + (t299 + t256) * qJ(2)) * MDP(6) + (qJDD(1) * t260 - 0.2e1 * t263 * t296) * MDP(7) + 0.2e1 * (-t263 * t301 + t311 * t303) * MDP(8) + (-t263 * t269 + t244) * MDP(9) + (-t266 * t269 - t300) * MDP(10) + (t273 * t263 + t281 * t266) * MDP(12) + (-t281 * t263 + t273 * t266) * MDP(13) + t285 * MDP(14) + (t177 * t214 - t178 * t215 + t186 * t204 - t187 * t206) * MDP(15) + t316 * MDP(16) + t284 * MDP(17) + (t235 * t178 + t222 * t187 + t195 * t214 + t229 * t204 - t288 * t245 + t292 * t253 + t275 * t255) * MDP(19) + (-t235 * t177 - t222 * t186 + t195 * t215 + t229 * t206 - t288 * t246 - t314 * t253 - t280 * t255) * MDP(20) + (-t168 * t204 - t169 * t206 + t177 * t181 - t178 * t182 - t274) * MDP(21) + (t167 * t182 + t174 * t168 + t166 * t181 + t172 * t169 + t278 * (pkin(4) * t214 + t235) + t188 * (pkin(4) * t187 + t229) - g(1) * (t218 * t267 + t248 + (-pkin(1) + t254) * t264) - g(2) * (t218 * t264 - t254 * t267 + t313)) * MDP(22); qJDD(1) * MDP(4) - t270 * MDP(5) + (t289 + t279) * MDP(6) + (t310 * t263 + t244) * MDP(12) + (t310 * t266 - t300) * MDP(13) + (-qJD(1) * t204 + t316) * MDP(19) + (-qJD(1) * t206 + t284) * MDP(20) + (-t178 * t214 - t187 * t204 - t285) * MDP(21) + (-qJD(1) * t188 + t274) * MDP(22); MDP(9) * t301 - MDP(10) * t302 + qJDD(3) * MDP(11) + (g(3) * t263 + t279 * t266 + t216) * MDP(12) + (g(3) * t266 + (-t227 - t279) * t263) * MDP(13) + (-t293 * t255 + (-t204 * t308 + t253 * t265 - t255 * t305) * pkin(3) + t271) * MDP(19) + (t315 * t255 + (-t206 * t308 - t262 * t253 - t255 * t304) * pkin(3) + t272) * MDP(20) + (t177 * t239 + (t174 + t175) * t206 + (-t172 + t176) * t204 + (-t178 * t262 + (-t204 * t265 + t206 * t262) * qJD(4)) * pkin(3)) * MDP(21) + (-pkin(4) * t318 + g(3) * t218 + t166 * t239 - t172 * t175 - t174 * t176 - t312 * (pkin(3) * t266 + pkin(4) * t246) + (-t188 * t308 + t167 * t262 + (-t172 * t262 + t174 * t265) * qJD(4)) * pkin(3)) * MDP(22) + t277 + (t266 * t263 * MDP(7) - t311 * MDP(8)) * t270; (-t283 * t255 + t271) * MDP(19) + (t294 * t255 + t272) * MDP(20) + (pkin(4) * t177 - t317 * t204) * MDP(21) + (t317 * t174 + (t166 - t318 + t330) * pkin(4)) * MDP(22) + t277; (-t203 - t324) * MDP(21) + (t172 * t206 + t174 * t204 + t278 - t288) * MDP(22);];
tau = t1;
