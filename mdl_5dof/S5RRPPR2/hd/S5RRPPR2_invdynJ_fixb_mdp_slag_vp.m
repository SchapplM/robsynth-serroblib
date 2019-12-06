% Calculate vector of inverse dynamics joint torques for
% S5RRPPR2
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
%   pkin=[a2,a3,a4,a5,d1,d2,d5,theta3,theta4]';
% MDP [18x1]
%   Minimal dynamic parameter vector (fixed base model)
%   see S5RRPPR2_convert_par2_MPV_fixb.m
% 
% Output:
% tau [5x1]
%   joint torques of inverse dynamics (contains inertial, gravitational Coriolis and centrifugal forces)

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-05 18:20
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tau = S5RRPPR2_invdynJ_fixb_mdp_slag_vp(qJ, qJD, qJDD, g, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(5,1),zeros(3,1),zeros(9,1),zeros(18,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRPPR2_invdynJ_fixb_mdp_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RRPPR2_invdynJ_fixb_mdp_slag_vp: qJD has to be [5x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [5 1]), ...
  'S5RRPPR2_invdynJ_fixb_mdp_slag_vp: qJDD has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RRPPR2_invdynJ_fixb_mdp_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5RRPPR2_invdynJ_fixb_mdp_slag_vp: pkin has to be [9x1] (double)');
assert(isreal(MDP) && all(size(MDP) == [18 1]), ...
  'S5RRPPR2_invdynJ_fixb_mdp_slag_vp: MDP has to be [18x1] (double)'); 

%% Symbolic Calculation
% From invdyn_joint_fixb_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 18:20:27
% EndTime: 2019-12-05 18:20:32
% DurationCPUTime: 1.88s
% Computational Cost: add. (1255->239), mult. (2000->330), div. (0->0), fcn. (1204->14), ass. (0->132)
t261 = cos(qJ(2));
t324 = pkin(1) * qJD(2);
t296 = qJD(1) * t324;
t258 = sin(qJ(2));
t299 = qJDD(1) * t258;
t338 = pkin(1) * t299 + t261 * t296;
t252 = qJ(1) + qJ(2);
t242 = pkin(8) + t252;
t231 = sin(t242);
t232 = cos(t242);
t337 = g(2) * t232 + g(3) * t231;
t332 = pkin(1) * t261;
t238 = qJDD(1) * t332;
t246 = qJDD(1) + qJDD(2);
t203 = pkin(2) * t246 - t258 * t296 + t238;
t254 = sin(pkin(8));
t256 = cos(pkin(8));
t187 = t254 * t203 + t338 * t256;
t249 = qJD(1) + qJD(2);
t181 = qJ(4) * t246 + qJD(4) * t249 + t187;
t253 = sin(pkin(9));
t255 = cos(pkin(9));
t180 = qJDD(3) * t253 + t181 * t255;
t179 = -t255 * qJDD(3) + t181 * t253;
t323 = t179 * t253;
t336 = t180 * t255 + t323;
t243 = sin(t252);
t244 = cos(t252);
t335 = g(2) * t244 + g(3) * t243;
t257 = sin(qJ(5));
t260 = cos(qJ(5));
t273 = MDP(17) * t257 + MDP(18) * t260;
t325 = pkin(1) * qJD(1);
t298 = t258 * t325;
t223 = t254 * t298;
t297 = t261 * t325;
t208 = t256 * t297 - t223;
t301 = qJD(4) - t208;
t334 = pkin(1) * t258;
t259 = sin(qJ(1));
t333 = pkin(1) * t259;
t262 = cos(qJ(1));
t331 = pkin(1) * t262;
t330 = pkin(2) * t243;
t329 = pkin(2) * t244;
t328 = pkin(2) * t256;
t212 = pkin(2) * t249 + t297;
t224 = t256 * t298;
t196 = t254 * t212 + t224;
t192 = qJ(4) * t249 + t196;
t189 = -t255 * qJD(3) + t192 * t253;
t322 = t189 * t253;
t274 = -pkin(4) * t255 - pkin(7) * t253 - pkin(3);
t228 = t254 * t334;
t237 = pkin(2) + t332;
t287 = t237 * t256 - t228;
t193 = t274 - t287;
t321 = t193 * t260;
t272 = t256 * t261 * t324 - qJD(2) * t228;
t202 = qJD(4) + t272;
t320 = t202 * t249;
t319 = t202 * t257;
t318 = t246 * t255;
t317 = t246 * t257;
t316 = t249 * t255;
t315 = t253 * t257;
t314 = t255 * t257;
t313 = t255 * t260;
t312 = t256 * t258;
t311 = t257 * t260;
t310 = t337 * t255;
t309 = pkin(1) * t312 + t254 * t237;
t247 = t253 ^ 2;
t248 = t255 ^ 2;
t308 = t247 + t248;
t251 = t260 ^ 2;
t307 = t257 ^ 2 - t251;
t304 = qJD(5) * t257;
t303 = qJD(5) * t260;
t217 = -qJDD(5) + t318;
t302 = t217 * MDP(16);
t220 = -qJD(5) + t316;
t300 = -qJD(5) - t220;
t294 = t238 + t335;
t293 = t249 * t303;
t292 = -g(2) * t243 + g(3) * t244;
t291 = t308 * t246;
t186 = t256 * t203 - t338 * t254;
t270 = qJDD(4) - t186;
t176 = t246 * t274 + t270;
t290 = t260 * t176 - t257 * t180;
t195 = t212 * t256 - t223;
t289 = t217 - t318;
t288 = t217 + t318;
t286 = t301 * t260;
t285 = qJD(1) * (-qJD(2) + t249);
t284 = qJD(2) * (-qJD(1) - t249);
t281 = qJD(4) - t195;
t280 = t257 * t176 + t260 * t180;
t184 = t249 * t274 + t281;
t190 = qJD(3) * t253 + t192 * t255;
t279 = t184 * t260 - t190 * t257;
t278 = -t184 * t257 - t190 * t260;
t277 = t190 * t255 + t322;
t205 = -pkin(3) - t287;
t207 = (t254 * t261 + t312) * t324;
t276 = -t205 * t246 - t207 * t249;
t206 = t254 * t297 + t224;
t233 = -pkin(3) - t328;
t275 = t206 * t249 - t233 * t246;
t271 = -pkin(3) * t231 + t232 * qJ(4) - t330;
t209 = t274 - t328;
t230 = pkin(2) * t254 + qJ(4);
t269 = t209 * t260 - t230 * t314;
t197 = t231 * t314 + t232 * t260;
t199 = -t231 * t260 + t232 * t314;
t268 = -g(2) * t199 - g(3) * t197 + (qJD(5) * t279 + t280) * t255 + t260 * t323;
t198 = t231 * t313 - t232 * t257;
t200 = -t231 * t257 - t232 * t313;
t267 = -g(2) * t200 + g(3) * t198 + t179 * t315 + t303 * t322;
t266 = g(2) * t231 - g(3) * t232 + t336;
t265 = -t232 * pkin(3) - t231 * qJ(4) - t329;
t264 = t230 * t303 + t257 * t301;
t211 = t253 * t304 * t316;
t263 = t255 * t302 + (t211 + (t304 * t220 - t260 * t288) * t253) * MDP(14) + (t288 * t257 + (t220 + t316) * t303) * t253 * MDP(15) + t246 * MDP(4) + (0.2e1 * (qJD(5) * t249 * t307 - t246 * t311) * MDP(13) + (t246 * t251 - 0.2e1 * t257 * t293) * MDP(12)) * t247;
t245 = t249 ^ 2;
t204 = qJ(4) + t309;
t191 = -pkin(3) * t249 + t281;
t183 = -pkin(3) * t246 + t270;
t182 = t183 * t253;
t171 = qJD(5) * t278 + t290;
t1 = [(t204 * t291 + t308 * t320 + t266) * MDP(10) + (t183 * t205 + t191 * t207 - g(2) * (t265 - t331) - g(3) * (t271 - t333) + t336 * t204 + t277 * t202) * MDP(11) + qJDD(1) * MDP(1) + ((t246 * t261 + t258 * t284) * pkin(1) + t294) * MDP(5) + ((-t183 + t276) * t255 + t310) * MDP(8) + (((-qJDD(1) - t246) * t258 + t261 * t284) * pkin(1) + t292) * MDP(6) + (-(-t193 * t304 + t207 * t260) * t220 - t217 * t321 + (-(-t204 * t303 - t319) * t220 + t204 * t257 * t217 - t171) * t255 + (t249 * t319 + (t293 + t317) * t204) * t247 + t267) * MDP(17) + t263 + (t182 + (-t276 - t337) * t253) * MDP(9) + (g(2) * t262 + g(3) * t259) * MDP(2) + (-g(2) * t259 + g(3) * t262) * MDP(3) + (t187 * t309 + t196 * t272 + t186 * t287 - t195 * t207 - g(2) * (-t329 - t331) - g(3) * (-t330 - t333)) * MDP(7) + ((t202 * t313 + t207 * t257) * t220 + (t193 * t257 + t204 * t313) * t217 + (t204 * t246 + t320) * t260 * t247 + (t220 * t321 + (-t322 + (-t220 * t255 - t247 * t249) * t204) * t257) * qJD(5) + t268) * MDP(18); (t285 * t334 + t294) * MDP(5) + ((t261 * t285 - t299) * pkin(1) + t292) * MDP(6) + (t195 * t206 - t196 * t208 + (t186 * t256 + t187 * t254 + t335) * pkin(2)) * MDP(7) + ((-t183 + t275) * t255 + t310) * MDP(8) + (t182 + (-t275 - t337) * t253) * MDP(9) + (t301 * t249 * t308 + t230 * t291 + t266) * MDP(10) + (t183 * t233 - t191 * t206 - g(2) * t265 - g(3) * t271 + (t180 * t230 + t190 * t301) * t255 + (t179 * t230 + t189 * t301) * t253) * MDP(11) + (-t269 * t217 - t171 * t255 + (t260 * t206 + t209 * t304 + t255 * t264) * t220 + (t230 * t317 + t249 * t264) * t247 + t267) * MDP(17) + ((t209 * t257 + t230 * t313) * t217 + (-t257 * t206 + t255 * t286) * t220 + (-t189 * t315 + t220 * t269) * qJD(5) + (t230 * t260 * t246 + (-t230 * t304 + t286) * t249) * t247 + t268) * MDP(18) + t263; (qJDD(3) - g(1)) * MDP(7) + (-t179 * t255 - g(1)) * MDP(11) + t211 * MDP(18) + (t180 * MDP(11) + (-qJD(5) * t220 * MDP(18) + MDP(17) * t289) * t257 + (t289 * MDP(18) + (t220 - t316) * MDP(17) * qJD(5)) * t260) * t253; (-t249 * t277 + t270 - t337) * MDP(11) + (-MDP(11) * pkin(3) - MDP(8) * t255 + MDP(9) * t253) * t246 + (-MDP(17) * t260 + MDP(18) * t257) * t217 + (-t248 * MDP(10) + (-MDP(10) - t273) * t247) * t245 - t273 * t220 ^ 2; -t302 + (-g(2) * t197 + g(3) * t199 + t220 * t278 + t290) * MDP(17) + (-g(2) * t198 - g(3) * t200 - t220 * t279 - t280) * MDP(18) + (MDP(17) * t278 - MDP(18) * t279) * qJD(5) + (MDP(12) * t311 - MDP(13) * t307) * t247 * t245 + ((MDP(14) * t260 - MDP(15) * t257) * t246 + t273 * g(1) + ((MDP(15) * t300 - t189 * MDP(17)) * t260 + (MDP(14) * t300 + t189 * MDP(18)) * t257) * t249) * t253;];
tau = t1;
