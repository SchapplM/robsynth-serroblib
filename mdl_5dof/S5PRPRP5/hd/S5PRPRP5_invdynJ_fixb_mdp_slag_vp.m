% Calculate vector of inverse dynamics joint torques for
% S5PRPRP5
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
%   pkin=[a2,a3,a4,a5,d2,d4,theta1,theta3]';
% MDP [19x1]
%   Minimal dynamic parameter vector (fixed base model)
%   see S5PRPRP5_convert_par2_MPV_fixb.m
% 
% Output:
% tau [5x1]
%   joint torques of inverse dynamics (contains inertial, gravitational Coriolis and centrifugal forces)

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-05 15:39
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tau = S5PRPRP5_invdynJ_fixb_mdp_slag_vp(qJ, qJD, qJDD, g, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(5,1),zeros(3,1),zeros(8,1),zeros(19,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PRPRP5_invdynJ_fixb_mdp_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5PRPRP5_invdynJ_fixb_mdp_slag_vp: qJD has to be [5x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [5 1]), ...
  'S5PRPRP5_invdynJ_fixb_mdp_slag_vp: qJDD has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5PRPRP5_invdynJ_fixb_mdp_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5PRPRP5_invdynJ_fixb_mdp_slag_vp: pkin has to be [8x1] (double)');
assert(isreal(MDP) && all(size(MDP) == [19 1]), ...
  'S5PRPRP5_invdynJ_fixb_mdp_slag_vp: MDP has to be [19x1] (double)'); 

%% Symbolic Calculation
% From invdyn_joint_fixb_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 15:38:49
% EndTime: 2019-12-05 15:38:53
% DurationCPUTime: 2.35s
% Computational Cost: add. (1186->246), mult. (2556->300), div. (0->0), fcn. (1847->10), ass. (0->123)
t246 = sin(pkin(8));
t324 = pkin(6) + qJ(3);
t224 = t324 * t246;
t248 = cos(pkin(8));
t225 = t324 * t248;
t251 = sin(qJ(4));
t328 = cos(qJ(4));
t194 = -t251 * t224 + t225 * t328;
t245 = pkin(8) + qJ(4);
t240 = sin(t245);
t253 = cos(qJ(2));
t293 = t328 * t248;
t313 = t251 * t246;
t269 = t293 - t313;
t266 = t253 * t269;
t270 = -t224 * t328 - t251 * t225;
t311 = -qJD(1) * t266 + qJD(3) * t269 + qJD(4) * t270;
t242 = g(3) * t253;
t252 = sin(qJ(2));
t247 = sin(pkin(7));
t249 = cos(pkin(7));
t327 = g(1) * t249;
t280 = g(2) * t247 + t327;
t335 = -t280 * t252 + t242;
t342 = t311 * qJD(4) + qJDD(4) * t194 - t240 * t335;
t301 = qJDD(1) * t253;
t302 = qJD(1) * qJD(2);
t277 = t252 * t302 + qJDD(3) - t301;
t323 = qJDD(2) * pkin(2);
t212 = t277 - t323;
t308 = qJD(1) * t252;
t229 = qJD(2) * qJ(3) + t308;
t309 = t246 ^ 2 + t248 ^ 2;
t341 = qJD(2) * t229 * t309 - t212;
t340 = t248 * MDP(5) - t246 * MDP(6);
t307 = qJD(1) * t253;
t283 = qJD(3) - t307;
t339 = MDP(7) * t309;
t338 = t253 * t280;
t221 = t246 * t328 + t251 * t248;
t333 = t221 * qJD(2);
t218 = t221 * qJD(4);
t332 = MDP(14) + MDP(16);
t296 = MDP(15) - MDP(18);
t329 = t333 ^ 2;
t325 = g(3) * t252;
t322 = qJDD(4) * pkin(4);
t281 = qJD(2) * t293;
t292 = qJD(2) * t313;
t213 = -t281 + t292;
t321 = t213 * t333;
t241 = cos(t245);
t320 = t241 * t252;
t318 = t247 * t241;
t317 = t247 * t253;
t315 = t249 * t253;
t288 = pkin(6) * qJD(2) + t229;
t207 = t288 * t248;
t314 = t251 * t207;
t312 = qJDD(1) - g(3);
t310 = qJD(4) * t194 + t283 * t221;
t206 = t288 * t246;
t180 = -t251 * t206 + t207 * t328;
t305 = qJD(4) * t180;
t304 = qJD(4) * t251;
t179 = -t206 * t328 - t314;
t303 = qJD(5) - t179;
t300 = qJDD(2) * qJ(3);
t299 = qJDD(2) * t251;
t298 = qJDD(4) * qJ(5);
t208 = t300 + qJDD(1) * t252 + (qJD(3) + t307) * qJD(2);
t286 = pkin(6) * qJDD(2) + t208;
t191 = t286 * t246;
t192 = t286 * t248;
t290 = qJD(4) * t328;
t295 = -t251 * t191 + t328 * t192 - t206 * t290;
t289 = qJDD(2) * t328;
t294 = qJD(4) * t281 + t246 * t289 + t248 * t299;
t238 = pkin(3) * t248 + pkin(2);
t291 = t246 * t304;
t287 = t309 * t208;
t285 = t179 + t314;
t217 = -t248 * t290 + t291;
t173 = pkin(4) * t218 + qJ(5) * t217 - qJD(5) * t221;
t284 = -t173 + t308;
t282 = t328 * t191 + t251 * t192 - t206 * t304 + t207 * t290;
t279 = t246 * t299 - t248 * t289;
t278 = -MDP(4) + t339;
t274 = pkin(4) * t241 + qJ(5) * t240 + t238;
t271 = MDP(3) + t340;
t210 = t269 * t252;
t264 = -t325 - t338;
t219 = -qJD(2) * t238 + t283;
t262 = t283 * t309;
t202 = t240 * t317 + t241 * t249;
t204 = t240 * t315 - t318;
t261 = g(1) * t204 + g(2) * t202 + t240 * t325 - t282;
t198 = -qJDD(2) * t238 + t277;
t203 = -t249 * t240 + t241 * t317;
t205 = t240 * t247 + t241 * t315;
t260 = g(1) * t205 + g(2) * t203 + g(3) * t320 - t295;
t259 = g(2) * t252 * t318 - qJD(4) * t310 + qJDD(4) * t270 - t241 * t242 + t320 * t327;
t258 = t287 + t264;
t174 = pkin(4) * t213 - qJ(5) * t333 + t219;
t257 = t174 * t333 + qJDD(5) - t261;
t185 = qJD(2) * t291 - t294;
t186 = qJD(2) * t218 + t279;
t256 = pkin(4) * t186 + qJ(5) * t185 + t198;
t254 = qJD(2) ^ 2;
t226 = -qJD(2) * pkin(2) + t283;
t211 = t213 ^ 2;
t209 = t221 * t252;
t184 = -pkin(4) * t269 - qJ(5) * t221 - t238;
t183 = pkin(4) * t333 + qJ(5) * t213;
t182 = qJD(4) * t210 + t253 * t333;
t181 = qJD(2) * t266 - t218 * t252;
t176 = qJD(4) * qJ(5) + t180;
t175 = -qJD(4) * pkin(4) + t303;
t172 = (t213 - t292) * qJD(4) + t294;
t170 = -qJD(5) * t333 + t256;
t169 = qJDD(5) + t282 - t322;
t168 = t298 + (qJD(5) - t314) * qJD(4) + t295;
t1 = [t312 * MDP(1) - g(3) * MDP(8) + (-t181 * t213 + t182 * t333 - t185 * t209 - t186 * t210) * MDP(17) + (t168 * t210 + t169 * t209 + t175 * t182 + t176 * t181 - g(3)) * MDP(19) + (-t170 * MDP(19) + t341 * MDP(8) + t271 * qJDD(2) + t296 * t185 - t332 * t186 + t278 * t254) * t253 + (MDP(8) * t287 - t271 * t254 + t278 * qJDD(2) + (t174 * MDP(19) + t226 * MDP(8) + t296 * t333) * qJD(2)) * t252 - t296 * (qJD(4) * t181 + qJDD(4) * t210) + t332 * (qJD(2) * t252 * t213 - qJD(4) * t182 - qJDD(4) * t209); qJDD(2) * MDP(2) + (t301 - t335) * MDP(3) + (-t252 * t312 + t338) * MDP(4) + (qJD(2) * t262 + t300 * t309 + t258) * MDP(7) + (-t226 * t308 + (-t212 - t335) * pkin(2) + t258 * qJ(3) + t262 * t229) * MDP(8) + (-t185 * t221 - t217 * t333) * MDP(9) + (-t185 * t269 - t186 * t221 + t213 * t217 - t218 * t333) * MDP(10) + (-qJD(4) * t217 + qJDD(4) * t221) * MDP(11) + (-qJD(4) * t218 + qJDD(4) * t269) * MDP(12) + (-t186 * t238 - t198 * t269 - t213 * t308 + t218 * t219 + t259) * MDP(14) + (t185 * t238 + t198 * t221 - t217 * t219 - t308 * t333 - t342) * MDP(15) + (-t170 * t269 + t174 * t218 + t184 * t186 - t213 * t284 + t259) * MDP(16) + (t168 * t269 + t169 * t221 - t175 * t217 - t176 * t218 + t185 * t270 - t186 * t194 - t213 * t311 + t310 * t333 + t264) * MDP(17) + (-t170 * t221 + t174 * t217 + t184 * t185 + t284 * t333 + t342) * MDP(18) + (t168 * t194 - t169 * t270 + t170 * t184 + t174 * t173 + t311 * t176 + t310 * t175 + (-g(3) * t274 - t280 * t324) * t253 + (-g(3) * t324 - t174 * qJD(1) + t274 * t280) * t252) * MDP(19) + t340 * (t252 * (t280 + t302) - t212 + t323 - t242); -t254 * t339 + (t335 - t341) * MDP(8) + (-t211 - t329) * MDP(17) + (t176 * t213 + (-qJD(5) - t175) * t333 + t256 + t335) * MDP(19) + t332 * (0.2e1 * qJD(4) * t333 + t279) - t296 * ((t213 + t292) * qJD(4) - t294) - t340 * qJDD(2); MDP(9) * t321 + (-t211 + t329) * MDP(10) + t172 * MDP(11) - t279 * MDP(12) + qJDD(4) * MDP(13) + (-t219 * t333 + t261 + t305) * MDP(14) + (qJD(4) * t285 + t213 * t219 + t260) * MDP(15) + (-t183 * t213 - t257 + t305 + 0.2e1 * t322) * MDP(16) + (pkin(4) * t185 - qJ(5) * t186 + (t176 - t180) * t333 + (t175 - t303) * t213) * MDP(17) + (0.2e1 * t298 - t174 * t213 + t183 * t333 + (0.2e1 * qJD(5) - t285) * qJD(4) - t260) * MDP(18) + (t168 * qJ(5) - t169 * pkin(4) - t174 * t183 - t175 * t180 - g(1) * (-pkin(4) * t204 + qJ(5) * t205) - g(2) * (-pkin(4) * t202 + qJ(5) * t203) - (-pkin(4) * t240 + qJ(5) * t241) * t325 + t303 * t176) * MDP(19); (-qJDD(4) + t321) * MDP(16) + t172 * MDP(17) + (-qJD(4) ^ 2 - t329) * MDP(18) + (-qJD(4) * t176 + t257 - t322) * MDP(19);];
tau = t1;
