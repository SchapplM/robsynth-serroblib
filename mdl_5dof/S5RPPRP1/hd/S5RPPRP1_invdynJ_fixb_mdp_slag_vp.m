% Calculate vector of inverse dynamics joint torques for
% S5RPPRP1
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
% MDP [18x1]
%   Minimal dynamic parameter vector (fixed base model)
%   see S5RPPRP1_convert_par2_MPV_fixb.m
% 
% Output:
% tau [5x1]
%   joint torques of inverse dynamics (contains inertial, gravitational Coriolis and centrifugal forces)

% Quelle: HybrDyn-Toolbox
% Datum: 2021-01-15 16:56
% Revision: 24b2e7d74a0c1a3b64fa2f8f5ad758691ad61af3 (2021-01-15)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tau = S5RPPRP1_invdynJ_fixb_mdp_slag_vp(qJ, qJD, qJDD, g, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(5,1),zeros(3,1),zeros(8,1),zeros(18,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPPRP1_invdynJ_fixb_mdp_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RPPRP1_invdynJ_fixb_mdp_slag_vp: qJD has to be [5x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [5 1]), ...
  'S5RPPRP1_invdynJ_fixb_mdp_slag_vp: qJDD has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RPPRP1_invdynJ_fixb_mdp_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RPPRP1_invdynJ_fixb_mdp_slag_vp: pkin has to be [8x1] (double)');
assert(isreal(MDP) && all(size(MDP) == [18 1]), ...
  'S5RPPRP1_invdynJ_fixb_mdp_slag_vp: MDP has to be [18x1] (double)'); 

%% Symbolic Calculation
% From invdyn_joint_fixb_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2021-01-15 16:55:52
% EndTime: 2021-01-15 16:55:57
% DurationCPUTime: 2.40s
% Computational Cost: add. (1142->260), mult. (2242->343), div. (0->0), fcn. (1410->10), ass. (0->129)
t242 = sin(pkin(8));
t249 = cos(qJ(4));
t291 = qJD(1) * qJD(4);
t275 = t249 * t291;
t247 = sin(qJ(4));
t286 = qJDD(1) * t247;
t333 = t275 + t286;
t334 = t242 * t333;
t239 = qJ(1) + pkin(7);
t233 = sin(t239);
t234 = cos(t239);
t308 = g(2) * t234 + g(3) * t233;
t243 = sin(pkin(7));
t224 = pkin(1) * t243 + qJ(3);
t292 = qJD(1) * qJD(3);
t206 = qJDD(1) * t224 + t292;
t244 = cos(pkin(8));
t332 = -t242 * (-qJ(5) - pkin(6)) + (pkin(4) * t249 + pkin(3)) * t244;
t314 = t244 * t247;
t193 = -t233 * t314 - t234 * t249;
t195 = -t233 * t249 + t234 * t314;
t331 = -g(2) * t193 - g(3) * t195;
t304 = qJD(1) * t244;
t220 = -qJD(4) + t304;
t330 = qJD(4) + t220;
t317 = t242 * t247;
t329 = g(1) * t317 + t331;
t287 = qJDD(1) * t244;
t219 = -qJDD(4) + t287;
t328 = pkin(4) * t219;
t326 = pkin(4) * t247;
t231 = t244 * qJDD(2);
t322 = t206 * t242;
t190 = -t231 + t322;
t323 = t190 * t242;
t320 = t233 * t247;
t237 = t242 ^ 2;
t251 = qJD(1) ^ 2;
t319 = t237 * t251;
t316 = t242 * t249;
t315 = t244 * MDP(5);
t313 = t244 * t249;
t245 = cos(pkin(7));
t229 = -pkin(1) * t245 - pkin(2);
t203 = -pkin(3) * t244 - pkin(6) * t242 + t229;
t192 = qJD(1) * t203 + qJD(3);
t211 = t224 * qJD(1);
t199 = qJD(2) * t242 + t211 * t244;
t274 = t249 * t192 - t199 * t247;
t302 = qJD(1) * t249;
t278 = t242 * t302;
t173 = -qJ(5) * t278 + t274;
t170 = -pkin(4) * t220 + t173;
t312 = -t173 + t170;
t297 = qJD(4) * t249;
t300 = qJD(3) * t249;
t311 = t203 * t297 + t244 * t300;
t207 = t224 * t313;
t310 = t247 * t203 + t207;
t248 = sin(qJ(1));
t309 = t248 * pkin(1) + t233 * pkin(2);
t307 = t244 ^ 2 + t237;
t240 = t247 ^ 2;
t241 = t249 ^ 2;
t306 = t240 - t241;
t305 = MDP(17) * t242;
t303 = qJD(1) * t247;
t301 = qJD(3) * t244;
t299 = qJD(4) * t199;
t298 = qJD(4) * t247;
t296 = qJD(5) * t242;
t295 = t219 * MDP(12);
t232 = t244 * qJD(2);
t182 = qJD(5) - t232 + (pkin(4) * t303 + t211) * t242;
t294 = qJD(5) + t182;
t293 = qJ(5) * qJDD(1);
t290 = qJD(1) * qJD(5);
t288 = qJDD(1) * t229;
t285 = qJDD(1) * t249;
t284 = MDP(13) + MDP(15);
t283 = MDP(14) + MDP(16);
t282 = qJ(5) * t316;
t281 = t247 * t319;
t250 = cos(qJ(1));
t280 = t250 * pkin(1) + t234 * pkin(2) + t233 * qJ(3);
t279 = t242 * t303;
t277 = t220 * t298;
t276 = t224 * t298;
t273 = MDP(17) * (-t240 - t241);
t202 = (t224 + t326) * t242;
t272 = qJD(1) * t202 + t182;
t191 = qJDD(2) * t242 + t206 * t244;
t271 = -qJD(4) * t192 - t191;
t270 = t219 - t287;
t269 = t219 + t287;
t189 = qJDD(1) * t203 + qJDD(3);
t268 = t247 * t189 + t249 * t191 + t192 * t297 - t199 * t298;
t267 = qJD(4) * t279;
t266 = pkin(4) * t334 + qJDD(5) - t231;
t265 = g(2) * t195 - g(3) * t193;
t194 = t233 * t313 - t234 * t247;
t196 = t234 * t313 + t320;
t264 = -g(2) * t196 - g(3) * t194;
t262 = -g(2) * t233 + g(3) * t234;
t261 = -g(2) * t250 - g(3) * t248;
t260 = t191 * t244 + t323;
t259 = -t247 * t192 - t249 * t199;
t198 = t211 * t242 - t232;
t258 = t198 * t242 + t199 * t244;
t180 = t266 + t322;
t205 = (pkin(4) * t297 + qJD(3)) * t242;
t257 = qJD(1) * t205 + qJDD(1) * t202 + t180;
t185 = t249 * t189;
t255 = qJ(5) * t267 + t247 * t271 + t185;
t254 = -t220 ^ 2 - t319;
t253 = g(1) * t316 + g(2) * t194 - g(3) * t196 - t268;
t217 = t242 * t285;
t209 = qJDD(3) + t288;
t208 = t244 * t267;
t204 = t220 * t278;
t201 = t249 * t203;
t181 = -qJ(5) * t317 + t310;
t179 = -t282 + t201 + (-t224 * t247 - pkin(4)) * t244;
t174 = -qJ(5) * t279 - t259;
t172 = -t247 * t301 - t249 * t296 + (-t207 + (qJ(5) * t242 - t203) * t247) * qJD(4);
t171 = -t247 * t296 + (-t224 * t314 - t282) * qJD(4) + t311;
t169 = (-qJ(5) * t333 - t247 * t290) * t242 + t268;
t168 = -t328 + (-t299 + (-t290 - t293) * t242) * t249 + t255;
t1 = [qJDD(1) * MDP(1) + t261 * MDP(2) + (g(2) * t248 - g(3) * t250) * MDP(3) + (t261 + (t243 ^ 2 + t245 ^ 2) * qJDD(1) * pkin(1)) * pkin(1) * MDP(4) + (-t209 - t308 - t288) * t315 + (t206 * t307 + t260 + t262) * MDP(6) + (t209 * t229 - g(2) * t280 - g(3) * (-qJ(3) * t234 + t309) + t260 * t224 + t258 * qJD(3)) * MDP(7) + (qJDD(1) * t241 - 0.2e1 * t247 * t275) * t237 * MDP(8) + 0.2e1 * (-t247 * t285 + t306 * t291) * t237 * MDP(9) + (t208 + (-t249 * t269 + t277) * t242) * MDP(10) + (t269 * t247 + (t220 + t304) * t297) * t242 * MDP(11) + t244 * t295 + (-t185 * t244 - t201 * t219 + ((qJD(1) * t237 + t220 * t244) * t224 + t258) * t297 + (-(-qJD(4) * t203 - t301) * t220 - t271 * t244 + t237 * t292 + t323 + (qJDD(1) * t237 + t219 * t244) * t224) * t247 + t264) * MDP(13) + ((-t244 * t276 + t311) * t220 + t310 * t219 + t268 * t244 + (t190 * t249 - t198 * t298) * t242 + (t224 * t285 + (-t276 + t300) * qJD(1)) * t237 + t265) * MDP(14) + (-t168 * t244 - t172 * t220 - t179 * t219 + (t257 * t247 + t272 * t297) * t242 + t264) * MDP(15) + (t169 * t244 + t171 * t220 + t181 * t219 + (t257 * t249 - t272 * t298) * t242 + t265) * MDP(16) + ((-qJD(4) * t174 - qJDD(1) * t179 - t168 + (-qJD(4) * t181 - t172) * qJD(1)) * t249 + (qJD(4) * t170 - qJDD(1) * t181 - t169 + (qJD(4) * t179 - t171) * qJD(1)) * t247 - t308) * t305 + (t169 * t181 + t174 * t171 + t168 * t179 + t170 * t172 + t180 * t202 + t182 * t205 - g(2) * (pkin(4) * t320 + t280) - g(3) * (t332 * t233 + t309) + (-g(2) * t332 - g(3) * (-qJ(3) - t326)) * t234) * MDP(18); (qJDD(2) - g(1)) * MDP(4) + (-t190 * t244 - g(1)) * MDP(7) + (-t180 * t244 - g(1)) * MDP(18) + t283 * t208 + (t284 * (t270 * t247 + (t220 - t304) * t297) + t283 * (t249 * t270 - t277) + t191 * MDP(7) + (-t168 * t247 + t169 * t249 - t170 * t297 - t174 * t298) * MDP(18)) * t242; (qJDD(3) + t308) * MDP(7) + (t168 * t249 + t169 * t247 - t170 * t298 + t174 * t297 + t308) * MDP(18) - t307 * MDP(6) * t251 + t283 * (t247 * t219 + t249 * t254) + t284 * (-t219 * t249 + t247 * t254) + (-t258 * MDP(7) + (t170 * t314 - t174 * t313 - t182 * t242) * MDP(18)) * qJD(1) + (t229 * MDP(7) + t242 * t273 - t315) * qJDD(1); t249 * MDP(8) * t281 - t306 * MDP(9) * t319 + (-t330 * t279 + t217) * MDP(10) + (-t204 - t334) * MDP(11) - t295 + (-t247 * t191 - t198 * t278 + t330 * t259 + t185 + t329) * MDP(13) + (t198 * t279 - t274 * t220 + t253) * MDP(14) + (-0.2e1 * t328 - t174 * t220 + (-pkin(4) * t281 - t299 + (-t294 * qJD(1) - t293) * t242) * t249 + t255 + t329) * MDP(15) + (-pkin(4) * t241 * t319 - t173 * t220 + (qJ(5) * t286 + (qJ(5) * t297 + t294 * t247) * qJD(1)) * t242 + t253) * MDP(16) + (-pkin(4) * t285 + (pkin(4) * qJD(4) - t312) * t303) * t305 + (t312 * t174 + (t168 + (g(1) * t247 - t182 * t302) * t242 + t331) * pkin(4)) * MDP(18); -t204 * MDP(15) + t217 * MDP(16) + (g(1) * t244 + t266) * MDP(18) + t273 * t319 + (MDP(15) * t286 + (t206 + t262) * MDP(18) + ((MDP(15) * qJD(4) + t170 * MDP(18)) * t249 + ((-qJD(4) + t220) * MDP(16) + t174 * MDP(18)) * t247) * qJD(1)) * t242;];
tau = t1;
