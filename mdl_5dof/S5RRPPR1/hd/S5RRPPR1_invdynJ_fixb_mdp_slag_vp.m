% Calculate vector of inverse dynamics joint torques for
% S5RRPPR1
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
%   see S5RRPPR1_convert_par2_MPV_fixb.m
% 
% Output:
% tau [5x1]
%   joint torques of inverse dynamics (contains inertial, gravitational Coriolis and centrifugal forces)

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-05 18:19
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tau = S5RRPPR1_invdynJ_fixb_mdp_slag_vp(qJ, qJD, qJDD, g, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(5,1),zeros(3,1),zeros(9,1),zeros(18,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRPPR1_invdynJ_fixb_mdp_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RRPPR1_invdynJ_fixb_mdp_slag_vp: qJD has to be [5x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [5 1]), ...
  'S5RRPPR1_invdynJ_fixb_mdp_slag_vp: qJDD has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RRPPR1_invdynJ_fixb_mdp_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5RRPPR1_invdynJ_fixb_mdp_slag_vp: pkin has to be [9x1] (double)');
assert(isreal(MDP) && all(size(MDP) == [18 1]), ...
  'S5RRPPR1_invdynJ_fixb_mdp_slag_vp: MDP has to be [18x1] (double)'); 

%% Symbolic Calculation
% From invdyn_joint_fixb_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 18:18:34
% EndTime: 2019-12-05 18:18:36
% DurationCPUTime: 1.46s
% Computational Cost: add. (1127->216), mult. (1797->275), div. (0->0), fcn. (1162->16), ass. (0->129)
t287 = cos(qJ(2));
t337 = pkin(1) * qJD(2);
t319 = qJD(1) * t337;
t284 = sin(qJ(2));
t322 = qJDD(1) * t284;
t356 = pkin(1) * t322 + t287 * t319;
t279 = sin(pkin(9));
t281 = cos(pkin(9));
t283 = sin(qJ(5));
t286 = cos(qJ(5));
t227 = t279 * t286 + t281 * t283;
t278 = qJ(1) + qJ(2);
t268 = pkin(8) + t278;
t253 = sin(t268);
t254 = cos(t268);
t355 = g(2) * t254 + g(3) * t253;
t277 = qJD(1) + qJD(2);
t213 = t227 * t277;
t324 = t279 ^ 2 + t281 ^ 2;
t354 = t324 * t277;
t347 = pkin(1) * t287;
t260 = qJDD(1) * t347;
t273 = qJDD(1) + qJDD(2);
t215 = pkin(2) * t273 - t284 * t319 + t260;
t280 = sin(pkin(8));
t282 = cos(pkin(8));
t194 = t280 * t215 + t356 * t282;
t186 = qJ(4) * t273 + qJD(4) * t277 + t194;
t263 = t281 * qJDD(3);
t180 = -t186 * t279 + t263;
t181 = t279 * qJDD(3) + t281 * t186;
t353 = -t180 * t279 + t181 * t281;
t352 = -g(2) * t253 + g(3) * t254;
t269 = sin(t278);
t270 = cos(t278);
t351 = g(2) * t270 + g(3) * t269;
t338 = pkin(1) * qJD(1);
t320 = t287 * t338;
t230 = pkin(2) * t277 + t320;
t321 = t284 * t338;
t245 = t282 * t321;
t208 = t280 * t230 + t245;
t204 = qJ(4) * t277 + t208;
t195 = t281 * qJD(3) - t204 * t279;
t196 = t279 * qJD(3) + t281 * t204;
t350 = -t195 * t279 + t196 * t281;
t244 = t280 * t321;
t220 = t282 * t320 - t244;
t323 = qJD(4) - t220;
t349 = pkin(1) * t284;
t285 = sin(qJ(1));
t348 = pkin(1) * t285;
t288 = cos(qJ(1));
t346 = pkin(1) * t288;
t345 = pkin(2) * t269;
t344 = pkin(2) * t270;
t343 = pkin(2) * t282;
t342 = pkin(4) * t281;
t271 = t281 * pkin(7);
t276 = pkin(9) + qJ(5);
t267 = cos(t276);
t333 = t253 * t267;
t332 = t254 * t267;
t331 = t279 * t283;
t328 = t282 * t284;
t327 = t286 * t281;
t326 = t355 * t281;
t259 = pkin(2) + t347;
t325 = pkin(1) * t328 + t280 * t259;
t249 = t280 * t349;
t317 = t277 * t331;
t316 = t277 * t327;
t315 = qJD(5) * t316 + t227 * t273;
t314 = t260 + t351;
t313 = -pkin(3) - t342;
t312 = -g(2) * t269 + g(3) * t270;
t311 = t324 * t273;
t207 = t230 * t282 - t244;
t310 = t259 * t282 - t249;
t309 = qJD(1) * (-qJD(2) + t277);
t308 = qJD(2) * (-qJD(1) - t277);
t193 = t215 * t282 - t356 * t280;
t293 = qJDD(4) - t193;
t182 = t313 * t273 + t293;
t303 = qJD(4) - t207;
t197 = t313 * t277 + t303;
t222 = t227 * qJD(5);
t226 = -t327 + t331;
t306 = g(2) * t332 + g(3) * t333 + t182 * t226 + t197 * t222;
t217 = -pkin(3) - t310;
t237 = t273 * t327;
t304 = -t273 * t331 + t237;
t191 = -qJD(5) * t317 + t315;
t192 = t222 * t277 - t304;
t221 = t226 * qJD(5);
t200 = -qJD(5) * t221 + qJDD(5) * t227;
t201 = -qJD(5) * t222 - qJDD(5) * t226;
t211 = -t316 + t317;
t302 = (-t191 * t226 - t192 * t227 + t211 * t221 - t213 * t222) * MDP(13) + (t191 * t227 - t213 * t221) * MDP(12) + t200 * MDP(14) + t201 * MDP(15) + t273 * MDP(4);
t216 = qJ(4) + t325;
t205 = (-pkin(7) - t216) * t279;
t206 = t216 * t281 + t271;
t301 = t205 * t286 - t206 * t283;
t300 = t205 * t283 + t206 * t286;
t219 = (t280 * t287 + t328) * t337;
t299 = -t217 * t273 - t219 * t277;
t218 = t280 * t320 + t245;
t255 = -pkin(3) - t343;
t298 = t218 * t277 - t255 * t273;
t251 = pkin(2) * t280 + qJ(4);
t223 = (-pkin(7) - t251) * t279;
t224 = t251 * t281 + t271;
t297 = t223 * t286 - t224 * t283;
t296 = t223 * t283 + t224 * t286;
t295 = t282 * t287 * t337 - qJD(2) * t249;
t294 = -pkin(3) * t253 + t254 * qJ(4) - t345;
t292 = -t352 + t353;
t290 = -t254 * pkin(3) - t253 * qJ(4) - t344;
t266 = sin(t276);
t289 = t182 * t227 - t197 * t221 - t266 * t355;
t232 = t313 - t343;
t214 = qJD(4) + t295;
t209 = t217 - t342;
t203 = -pkin(3) * t277 + t303;
t188 = -pkin(3) * t273 + t293;
t187 = t188 * t279;
t176 = t273 * t271 + t181;
t175 = t263 + (-pkin(7) * t273 - t186) * t279;
t1 = [qJDD(1) * MDP(1) + (g(2) * t288 + g(3) * t285) * MDP(2) + (-g(2) * t285 + g(3) * t288) * MDP(3) + ((t273 * t287 + t284 * t308) * pkin(1) + t314) * MDP(5) + (((-qJDD(1) - t273) * t284 + t287 * t308) * pkin(1) + t312) * MDP(6) + (t194 * t325 + t208 * t295 + t193 * t310 - t207 * t219 - g(2) * (-t344 - t346) - g(3) * (-t345 - t348)) * MDP(7) + ((-t188 + t299) * t281 + t326) * MDP(8) + (t187 + (-t299 - t355) * t279) * MDP(9) + (t214 * t354 + t216 * t311 + t292) * MDP(10) + (t188 * t217 + t203 * t219 - g(2) * (t290 - t346) - g(3) * (t294 - t348) + t353 * t216 + t350 * t214) * MDP(11) + (t219 * t211 + t209 * t192 + t301 * qJDD(5) + (-qJD(5) * t300 - t214 * t227) * qJD(5) + t306) * MDP(17) + (t219 * t213 + t209 * t191 - t300 * qJDD(5) + (-qJD(5) * t301 + t214 * t226) * qJD(5) + t289) * MDP(18) + t302; (t309 * t349 + t314) * MDP(5) + ((t287 * t309 - t322) * pkin(1) + t312) * MDP(6) + (t207 * t218 - t208 * t220 + (t193 * t282 + t194 * t280 + t351) * pkin(2)) * MDP(7) + ((-t188 + t298) * t281 + t326) * MDP(8) + (t187 + (-t298 - t355) * t279) * MDP(9) + (t251 * t311 + t323 * t354 + t292) * MDP(10) + (t188 * t255 - t203 * t218 - g(2) * t290 - g(3) * t294 + (t181 * t251 + t323 * t196) * t281 + (-t180 * t251 - t323 * t195) * t279) * MDP(11) + (t232 * t192 + t297 * qJDD(5) - t218 * t211 + (-qJD(5) * t296 - t323 * t227) * qJD(5) + t306) * MDP(17) + (t232 * t191 - t296 * qJDD(5) - t218 * t213 + (-qJD(5) * t297 + t323 * t226) * qJD(5) + t289) * MDP(18) + t302; (qJDD(3) - g(1)) * MDP(7) + (t180 * t281 + t181 * t279 - g(1)) * MDP(11) + t201 * MDP(17) - t200 * MDP(18); (t293 - t355) * MDP(11) - t237 * MDP(17) + t315 * MDP(18) + (-pkin(3) * MDP(11) - t281 * MDP(8) + (t283 * MDP(17) + MDP(9)) * t279) * t273 + (0.2e1 * t213 * MDP(17) + (-t211 - t317) * MDP(18)) * qJD(5) + (-MDP(10) * t354 - t350 * MDP(11)) * t277; t213 * t211 * MDP(12) + (-t211 ^ 2 + t213 ^ 2) * MDP(13) + t304 * MDP(15) + qJDD(5) * MDP(16) + (-g(1) * t267 + t286 * t175 - t283 * t176 - t197 * t213 + t352 * t266) * MDP(17) + (g(1) * t266 - g(2) * t333 + g(3) * t332 - t283 * t175 - t286 * t176 + t197 * t211) * MDP(18) + (t315 + (t211 - t317) * qJD(5)) * MDP(14);];
tau = t1;
