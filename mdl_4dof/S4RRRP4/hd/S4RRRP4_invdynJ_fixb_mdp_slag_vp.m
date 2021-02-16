% Calculate vector of inverse dynamics joint torques for
% S4RRRP4
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
% MDP [21x1]
%   Minimal dynamic parameter vector (fixed base model)
%   see S4RRRP4_convert_par2_MPV_fixb.m
% 
% Output:
% tau [4x1]
%   joint torques of inverse dynamics (contains inertial, gravitational Coriolis and centrifugal forces)

% Quelle: HybrDyn-Toolbox
% Datum: 2021-01-15 14:30
% Revision: 24b2e7d74a0c1a3b64fa2f8f5ad758691ad61af3 (2021-01-15)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tau = S4RRRP4_invdynJ_fixb_mdp_slag_vp(qJ, qJD, qJDD, g, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(4,1),zeros(4,1),zeros(3,1),zeros(6,1),zeros(21,1)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4RRRP4_invdynJ_fixb_mdp_slag_vp: qJ has to be [4x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [4 1]), ...
  'S4RRRP4_invdynJ_fixb_mdp_slag_vp: qJD has to be [4x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [4 1]), ...
  'S4RRRP4_invdynJ_fixb_mdp_slag_vp: qJDD has to be [4x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S4RRRP4_invdynJ_fixb_mdp_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [6 1]), ...
  'S4RRRP4_invdynJ_fixb_mdp_slag_vp: pkin has to be [6x1] (double)');
assert(isreal(MDP) && all(size(MDP) == [21 1]), ...
  'S4RRRP4_invdynJ_fixb_mdp_slag_vp: MDP has to be [21x1] (double)'); 

%% Symbolic Calculation
% From invdyn_joint_fixb_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2021-01-15 14:30:27
% EndTime: 2021-01-15 14:30:32
% DurationCPUTime: 2.00s
% Computational Cost: add. (1225->246), mult. (2831->313), div. (0->0), fcn. (1854->8), ass. (0->125)
t263 = cos(qJ(2));
t334 = pkin(5) + pkin(6);
t226 = t334 * t263;
t221 = qJD(1) * t226;
t260 = sin(qJ(3));
t205 = t260 * t221;
t261 = sin(qJ(2));
t225 = t334 * t261;
t219 = qJD(1) * t225;
t329 = qJD(2) * pkin(2);
t211 = -t219 + t329;
t333 = cos(qJ(3));
t294 = t333 * t211 - t205;
t214 = t260 * t263 + t333 * t261;
t203 = t214 * qJD(1);
t324 = t203 * qJ(4);
t178 = -t324 + t294;
t256 = qJD(2) + qJD(3);
t253 = t263 * pkin(2);
t330 = pkin(1) + t253;
t308 = qJD(1) * qJD(2);
t297 = t263 * t308;
t307 = qJDD(1) * t261;
t191 = qJDD(2) * pkin(2) + t334 * (-t297 - t307);
t298 = t261 * t308;
t306 = qJDD(1) * t263;
t193 = t334 * (-t298 + t306);
t340 = t333 * t191 - t260 * t193;
t314 = -t260 * t225 + t333 * t226;
t254 = qJDD(2) + qJDD(3);
t299 = t333 * qJD(3);
t332 = pkin(2) * t256;
t339 = -t260 * pkin(2) * t254 - t299 * t332;
t250 = t254 * pkin(3);
t319 = t260 * t261;
t284 = t256 * t319;
t301 = t333 * t263;
t289 = qJD(1) * t301;
t296 = qJDD(1) * t333;
t290 = t256 * t289 + t260 * t306 + t261 * t296;
t180 = qJD(1) * t284 - t290;
t328 = t180 * qJ(4);
t338 = t250 + t328;
t337 = t333 * qJD(2) + t299;
t259 = qJ(2) + qJ(3);
t251 = sin(t259);
t262 = sin(qJ(1));
t318 = t262 * t251;
t264 = cos(qJ(1));
t322 = t251 * t264;
t252 = cos(t259);
t331 = g(3) * t252;
t336 = g(1) * t322 + g(2) * t318 - t331;
t335 = t203 ^ 2;
t190 = t256 * t214;
t283 = t260 * t307 - t263 * t296;
t181 = t190 * qJD(1) + t283;
t327 = t181 * qJ(4);
t311 = qJD(1) * t261;
t201 = t260 * t311 - t289;
t326 = t201 * qJ(4);
t325 = t201 * t256;
t321 = t252 * t262;
t317 = t264 * t252;
t177 = pkin(3) * t256 + t178;
t316 = t177 - t178;
t315 = -t333 * t219 - t205;
t313 = pkin(3) * t252 + t253;
t257 = t261 ^ 2;
t312 = -t263 ^ 2 + t257;
t310 = qJD(3) * t260;
t224 = t330 * qJD(1);
t295 = pkin(3) * t201 + qJD(4);
t192 = -t224 + t295;
t309 = qJD(4) + t192;
t305 = pkin(2) * t311;
t304 = t261 * t329;
t302 = qJD(2) * t334;
t209 = t333 * t221;
t293 = t219 * t260 - t209;
t292 = -t333 * t225 - t226 * t260;
t291 = t256 * t261;
t288 = -g(1) * t318 + g(2) * t322;
t287 = g(1) * t321 - g(2) * t317;
t286 = g(1) * t264 + g(2) * t262;
t285 = g(1) * t262 - g(2) * t264;
t282 = -0.2e1 * pkin(1) * t308 - pkin(5) * qJDD(2);
t281 = -t260 * t211 - t209;
t198 = pkin(2) * t298 - qJDD(1) * t330;
t220 = t261 * t302;
t222 = t263 * t302;
t278 = -t333 * t220 - t260 * t222 - t225 * t299 - t226 * t310;
t200 = t201 ^ 2;
t277 = t203 * t201 * MDP(11) + (-t260 * qJD(1) * t291 + t290 + t325) * MDP(13) - t283 * MDP(14) + (-t200 + t335) * MDP(12) + t254 * MDP(15);
t265 = qJD(2) ^ 2;
t276 = 0.2e1 * qJDD(1) * pkin(1) - pkin(5) * t265 + t285;
t266 = qJD(1) ^ 2;
t275 = pkin(1) * t266 - pkin(5) * qJDD(1) + t286;
t176 = pkin(3) * t181 + qJDD(4) + t198;
t274 = t281 * qJD(3) + t340;
t273 = -qJD(3) * t314 + t260 * t220 - t333 * t222;
t272 = t260 * t191 + t333 * t193 + t211 * t299 - t221 * t310;
t271 = g(1) * t317 + g(2) * t321 + g(3) * t251 - t272;
t270 = t274 + t336;
t269 = -t224 * t201 + t271;
t268 = t224 * t203 + t270;
t267 = t309 * t201 + t271 + t327;
t255 = -qJ(4) - t334;
t248 = t333 * pkin(2) + pkin(3);
t218 = pkin(1) + t313;
t213 = -t301 + t319;
t195 = pkin(3) * t213 - t330;
t194 = pkin(3) * t203 + t305;
t189 = -t263 * t337 + t284;
t187 = pkin(3) * t190 + t304;
t186 = -qJ(4) * t213 + t314;
t185 = -qJ(4) * t214 + t292;
t184 = -t324 + t315;
t183 = t293 + t326;
t179 = -t281 - t326;
t173 = t189 * qJ(4) - t214 * qJD(4) + t273;
t172 = -qJ(4) * t190 - qJD(4) * t213 + t278;
t171 = -t201 * qJD(4) + t272 - t327;
t170 = -t203 * qJD(4) + t274 + t338;
t1 = [qJDD(1) * MDP(1) + t285 * MDP(2) + t286 * MDP(3) + (qJDD(1) * t257 + 0.2e1 * t261 * t297) * MDP(4) + 0.2e1 * (t261 * t306 - t312 * t308) * MDP(5) + (qJDD(2) * t261 + t263 * t265) * MDP(6) + (qJDD(2) * t263 - t261 * t265) * MDP(7) + (t282 * t261 + t276 * t263) * MDP(9) + (-t276 * t261 + t282 * t263) * MDP(10) + (-t180 * t214 - t189 * t203) * MDP(11) + (t180 * t213 - t181 * t214 + t189 * t201 - t190 * t203) * MDP(12) + (-t189 * t256 + t214 * t254) * MDP(13) + (-t190 * t256 - t213 * t254) * MDP(14) + (-t181 * t330 - t224 * t190 + t198 * t213 + t201 * t304 + t292 * t254 + t273 * t256 + t287) * MDP(16) + (t180 * t330 + t224 * t189 + t198 * t214 + t203 * t304 - t314 * t254 - t278 * t256 + t288) * MDP(17) + (t173 * t256 + t176 * t213 + t181 * t195 + t185 * t254 + t187 * t201 + t190 * t192 + t287) * MDP(18) + (-t172 * t256 + t176 * t214 - t180 * t195 - t186 * t254 + t187 * t203 - t189 * t192 + t288) * MDP(19) + (-t170 * t214 - t171 * t213 - t172 * t201 - t173 * t203 + t177 * t189 - t179 * t190 + t180 * t185 - t181 * t186 - t286) * MDP(20) + (t171 * t186 + t179 * t172 + t170 * t185 + t177 * t173 + t176 * t195 + t192 * t187 - g(1) * (-t218 * t262 - t255 * t264) - g(2) * (t218 * t264 - t255 * t262)) * MDP(21); MDP(6) * t307 + MDP(7) * t306 + qJDD(2) * MDP(8) + (-g(3) * t263 + t275 * t261) * MDP(9) + (g(3) * t261 + t275 * t263) * MDP(10) + (-t293 * t256 + (-t201 * t311 + t333 * t254 - t256 * t310) * pkin(2) + t268) * MDP(16) + (-t203 * t305 + t315 * t256 + t269 + t339) * MDP(17) + (-t183 * t256 - t194 * t201 + t248 * t254 - t309 * t203 + (-t209 + (-t211 - t332) * t260) * qJD(3) + t336 + t338 + t340) * MDP(18) + (t184 * t256 - t194 * t203 + t267 + t339) * MDP(19) + (t248 * t180 + (t179 + t183) * t203 + (-t177 + t184) * t201 + (-t181 * t260 + (-t333 * t201 + t203 * t260) * qJD(3)) * pkin(2)) * MDP(20) + (t170 * t248 - t179 * t184 - t177 * t183 - t192 * t194 - g(3) * t313 - t286 * (-pkin(2) * t261 - pkin(3) * t251) + (t171 * t260 + (-t177 * t260 + t333 * t179) * qJD(3)) * pkin(2)) * MDP(21) + t277 + (-t261 * t263 * MDP(4) + t312 * MDP(5)) * t266; (-t281 * t256 + t268) * MDP(16) + (t294 * t256 + t269) * MDP(17) + (t328 + t179 * t256 + 0.2e1 * t250 + (-t192 - t295) * t203 + t270) * MDP(18) + (-t335 * pkin(3) + t178 * t256 + t267) * MDP(19) + (pkin(3) * t180 - t316 * t201) * MDP(20) + (t316 * t179 + (-t192 * t203 + t286 * t251 + t170 - t331) * pkin(3)) * MDP(21) + t277; (t203 * t256 + t283) * MDP(18) + (t290 - t325) * MDP(19) + (-t200 - t335) * MDP(20) + (t177 * t203 + t179 * t201 + t176 - t285) * MDP(21) + (t337 * t261 * MDP(18) + (t256 * MDP(18) * t263 - MDP(19) * t291) * t260) * qJD(1);];
tau = t1;
