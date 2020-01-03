% Calculate vector of inverse dynamics joint torques for
% S5RRPRP4
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
%   pkin=[a2,a3,a4,a5,d1,d2,d4]';
% MDP [20x1]
%   Minimal dynamic parameter vector (fixed base model)
%   see S5RRPRP4_convert_par2_MPV_fixb.m
% 
% Output:
% tau [5x1]
%   joint torques of inverse dynamics (contains inertial, gravitational Coriolis and centrifugal forces)

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 19:53
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tau = S5RRPRP4_invdynJ_fixb_mdp_slag_vp(qJ, qJD, qJDD, g, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(5,1),zeros(3,1),zeros(7,1),zeros(20,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRPRP4_invdynJ_fixb_mdp_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RRPRP4_invdynJ_fixb_mdp_slag_vp: qJD has to be [5x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [5 1]), ...
  'S5RRPRP4_invdynJ_fixb_mdp_slag_vp: qJDD has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RRPRP4_invdynJ_fixb_mdp_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [7 1]), ...
  'S5RRPRP4_invdynJ_fixb_mdp_slag_vp: pkin has to be [7x1] (double)');
assert(isreal(MDP) && all(size(MDP) == [20 1]), ...
  'S5RRPRP4_invdynJ_fixb_mdp_slag_vp: MDP has to be [20x1] (double)'); 

%% Symbolic Calculation
% From invdyn_joint_fixb_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 19:53:00
% EndTime: 2019-12-31 19:53:02
% DurationCPUTime: 1.37s
% Computational Cost: add. (1188->256), mult. (1453->276), div. (0->0), fcn. (651->8), ass. (0->127)
t236 = qJD(1) + qJD(2);
t240 = sin(qJ(4));
t231 = t240 * pkin(4);
t243 = cos(qJ(4));
t277 = qJ(5) * t243 - t231;
t333 = qJ(3) - t277;
t335 = t236 * t333;
t235 = qJDD(1) + qJDD(2);
t334 = t333 * t235;
t246 = -pkin(2) - pkin(7);
t241 = sin(qJ(2));
t244 = cos(qJ(2));
t327 = pkin(1) * qJD(1);
t286 = qJD(2) * t327;
t325 = pkin(1) * qJDD(1);
t303 = -t241 * t286 + t244 * t325;
t280 = qJDD(3) - t303;
t182 = t235 * t246 + t280;
t176 = t240 * t182;
t294 = qJDD(4) * qJ(5);
t289 = t244 * t327;
t269 = qJD(3) - t289;
t190 = t246 * t236 + t269;
t317 = t190 * t243;
t171 = t294 + t176 + (qJD(5) + t317) * qJD(4);
t177 = t243 * t182;
t296 = qJD(4) * t240;
t321 = qJDD(4) * pkin(4);
t172 = t190 * t296 + qJDD(5) - t177 - t321;
t276 = qJD(5) - t317;
t326 = qJD(4) * pkin(4);
t181 = t276 - t326;
t297 = qJD(4) * qJ(5);
t318 = t190 * t240;
t183 = t297 + t318;
t249 = (t181 * t240 + t183 * t243) * qJD(4) - t172 * t243 + t171 * t240;
t239 = qJ(1) + qJ(2);
t230 = cos(t239);
t219 = g(2) * t230;
t229 = sin(t239);
t220 = g(1) * t229;
t304 = t220 - t219;
t248 = t304 - t249;
t313 = t230 * t243;
t332 = g(2) * t313 + g(3) * t240;
t268 = -g(1) * t230 - g(2) * t229;
t273 = -qJD(5) * t243 + qJD(3);
t295 = qJD(4) * t243;
t189 = pkin(4) * t295 + qJ(5) * t296 + t273;
t247 = qJD(4) ^ 2;
t311 = t246 * t247;
t331 = -t334 + (-t189 + t289) * t236 + t311;
t234 = t236 ^ 2;
t330 = pkin(1) * t241;
t329 = pkin(2) * t235;
t242 = sin(qJ(1));
t328 = g(1) * t242;
t324 = qJ(3) * t235;
t323 = qJ(3) * t236;
t290 = t241 * t327;
t180 = t290 + t335;
t319 = t180 * t236;
t199 = t290 + t323;
t316 = t199 * t236;
t223 = -pkin(1) * t244 - pkin(2);
t212 = -pkin(7) + t223;
t315 = t212 * t247;
t314 = t229 * t243;
t312 = t240 * t243;
t267 = pkin(4) * t243 + qJ(5) * t240;
t306 = t241 * t325 + t244 * t286;
t169 = t334 + (t267 * qJD(4) + t273) * t236 + t306;
t310 = t169 * t240 + t180 * t295;
t184 = qJD(3) * t236 + t306 + t324;
t309 = t184 * t240 + t199 * t295;
t228 = qJDD(4) * t243;
t299 = qJD(2) * t241;
t288 = pkin(1) * t299;
t308 = t212 * t228 + t288 * t295;
t307 = g(1) * t313 + g(2) * t314;
t305 = t230 * pkin(2) + t229 * qJ(3);
t302 = t234 + t247;
t237 = t240 ^ 2;
t238 = t243 ^ 2;
t301 = t237 - t238;
t300 = t237 + t238;
t298 = qJD(2) * t244;
t293 = qJDD(4) * t212;
t292 = qJDD(4) * t240;
t291 = qJDD(4) * t246;
t287 = pkin(1) * t298;
t285 = t180 * t296 + t307;
t284 = t177 + t332;
t283 = t184 * t243 - t307;
t282 = t236 * t299;
t281 = t236 * t295;
t214 = t230 * qJ(3);
t279 = -pkin(2) * t229 + t214;
t278 = t300 * t235;
t275 = t236 * t290;
t274 = pkin(1) * t282;
t272 = -t303 - t304;
t271 = t306 + t268;
t266 = t319 + t220;
t265 = -t316 - t220;
t264 = 0.2e1 * (t301 * t236 * qJD(4) - t235 * t312) * MDP(11) + (t235 * t238 - 0.2e1 * t240 * t281) * MDP(10) + (-t243 * t247 - t292) * MDP(13) + (-t240 * t247 + t228) * MDP(12) + t235 * MDP(4);
t263 = t181 * t243 - t183 * t240;
t260 = qJDD(3) + t272;
t258 = -t290 + t323;
t257 = t230 * pkin(7) - qJ(5) * t314 + t229 * t231 + t305;
t256 = (-t290 + t335) * qJD(4);
t185 = t189 + t287;
t195 = t333 + t330;
t255 = -t185 * t236 - t195 * t235 + t315;
t208 = qJD(3) + t287;
t215 = qJ(3) + t330;
t254 = t208 * t236 + t215 * t235 - t315;
t253 = -t272 + t275;
t251 = g(1) * (-qJ(5) * t313 + t246 * t229 + t230 * t231 + t214);
t250 = t269 * t236 - t311 + t324;
t245 = cos(qJ(1));
t232 = t245 * pkin(1);
t209 = t243 * t291;
t196 = -pkin(2) * t236 + t269;
t193 = t267 * t236;
t188 = t280 - t329;
t1 = [((qJD(3) + t208) * t236 + (qJ(3) + t215) * t235 + t271) * MDP(8) + ((-t235 * t241 - t236 * t298) * pkin(1) - t271) * MDP(6) + ((t235 * t244 - t282) * pkin(1) - t272) * MDP(5) + t264 + (-g(2) * t245 + t328) * MDP(2) + (g(1) * t245 + g(2) * t242) * MDP(3) + (t184 * t215 + t199 * t208 + t188 * t223 + t196 * t288 - g(1) * (-pkin(1) * t242 + t279) - g(2) * (t232 + t305)) * MDP(9) + ((t293 + (t195 * t236 + t288) * qJD(4)) * t240 + (-t169 + t255) * t243 + t285) * MDP(19) + (t274 + (-pkin(2) + t223) * t235 + t260) * MDP(7) + (t169 * t195 + t180 * t185 - t251 - g(2) * (t232 + t257) + (-t263 * t299 + t328) * pkin(1) + t249 * t212) * MDP(20) + (t215 * t281 + (t254 + t268) * t240 + t308 + t309) * MDP(15) + (t254 * t243 + (-t293 + (-t215 * t236 - t199 - t288) * qJD(4)) * t240 + t283) * MDP(16) + (t195 * t281 + (-t255 + t268) * t240 + t308 + t310) * MDP(17) + (-t212 * t278 - t300 * t274 + t248) * MDP(18) + qJDD(1) * MDP(1); t253 * MDP(5) + (t236 * t289 - t271) * MDP(6) + (qJDD(3) - t253 - 0.2e1 * t329) * MDP(7) + (0.2e1 * t324 + (0.2e1 * qJD(3) - t289) * t236 + t271) * MDP(8) + (t184 * qJ(3) + t199 * qJD(3) - t188 * pkin(2) - g(1) * t279 - g(2) * t305 + (-t196 * t241 - t199 * t244) * t327) * MDP(9) + (t209 + t258 * t295 + (t250 + t268) * t240 + t309) * MDP(15) + (t250 * t243 + (-t291 + (-t199 - t258) * qJD(4)) * t240 + t283) * MDP(16) + (t209 + t243 * t256 + (t268 - t331) * t240 + t310) * MDP(17) + (-t246 * t278 + t300 * t275 + t248) * MDP(18) + ((t256 + t291) * t240 + (-t169 + t331) * t243 + t285) * MDP(19) + (t169 * t333 + t180 * t189 - t251 - g(2) * t257 + (-t180 * t244 + t263 * t241) * t327 + t249 * t246) * MDP(20) + t264; -t234 * MDP(8) + (t260 - t316) * MDP(9) + (-t319 - t248) * MDP(20) + (MDP(15) + MDP(17)) * (-t302 * t240 + t228) + (-MDP(16) + MDP(19)) * (t302 * t243 + t292) + (-t300 * MDP(18) - pkin(2) * MDP(9) + MDP(7)) * t235; qJDD(4) * MDP(14) + (t265 * t243 + t284) * MDP(15) + (g(3) * t243 - t176 + (-t265 - t219) * t240) * MDP(16) + (-g(1) * t314 + 0.2e1 * t321 - qJDD(5) + (-t180 * t243 - t193 * t240) * t236 + t284) * MDP(17) + ((t183 - t297) * t243 + (-qJD(5) + t181 + t326) * t240) * t236 * MDP(18) + (0.2e1 * t294 + 0.2e1 * qJD(4) * qJD(5) + t176 + (t193 * t236 - g(3)) * t243 + (-t266 + t219) * t240) * MDP(19) + (-t172 * pkin(4) - g(3) * t277 + t171 * qJ(5) - t180 * t193 - t181 * t318 + t276 * t183 - t304 * t267) * MDP(20) + (t243 * MDP(12) - t240 * MDP(13) - t267 * MDP(18)) * t235 + (MDP(10) * t312 - t301 * MDP(11)) * t234; -qJDD(4) * MDP(17) + (-t234 * t238 - t247) * MDP(19) + (-qJD(4) * t183 + t172 - t332) * MDP(20) + (t234 * t240 * MDP(17) + t235 * MDP(18) + t266 * MDP(20)) * t243;];
tau = t1;
