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
% Datum: 2020-01-03 11:56
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
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
% StartTime: 2020-01-03 11:56:03
% EndTime: 2020-01-03 11:56:06
% DurationCPUTime: 1.30s
% Computational Cost: add. (1127->213), mult. (1797->272), div. (0->0), fcn. (1162->16), ass. (0->130)
t278 = qJDD(1) + qJDD(2);
t292 = cos(qJ(2));
t352 = pkin(1) * t292;
t263 = qJDD(1) * t352;
t289 = sin(qJ(2));
t345 = pkin(1) * qJD(2);
t327 = qJD(1) * t345;
t216 = pkin(2) * t278 - t289 * t327 + t263;
t285 = sin(pkin(8));
t287 = cos(pkin(8));
t330 = qJDD(1) * t289;
t358 = pkin(1) * t330 + t292 * t327;
t194 = t216 * t287 - t358 * t285;
t297 = qJDD(4) - t194;
t189 = -pkin(3) * t278 + t297;
t283 = qJ(1) + qJ(2);
t271 = pkin(8) + t283;
t256 = sin(t271);
t257 = cos(t271);
t312 = -g(2) * t257 - g(3) * t256;
t299 = -t189 + t312;
t284 = sin(pkin(9));
t286 = cos(pkin(9));
t288 = sin(qJ(5));
t291 = cos(qJ(5));
t228 = t284 * t291 + t286 * t288;
t282 = qJD(1) + qJD(2);
t214 = t228 * t282;
t332 = t284 ^ 2 + t286 ^ 2;
t357 = t332 * t282;
t195 = t285 * t216 + t358 * t287;
t187 = qJ(4) * t278 + qJD(4) * t282 + t195;
t266 = t286 * qJDD(3);
t181 = -t187 * t284 + t266;
t182 = t284 * qJDD(3) + t286 * t187;
t356 = -t181 * t284 + t182 * t286;
t355 = g(2) * t256 - g(3) * t257;
t346 = pkin(1) * qJD(1);
t328 = t292 * t346;
t231 = pkin(2) * t282 + t328;
t329 = t289 * t346;
t246 = t287 * t329;
t209 = t285 * t231 + t246;
t205 = qJ(4) * t282 + t209;
t196 = t286 * qJD(3) - t205 * t284;
t197 = t284 * qJD(3) + t286 * t205;
t354 = -t196 * t284 + t197 * t286;
t245 = t285 * t329;
t221 = t287 * t328 - t245;
t331 = qJD(4) - t221;
t353 = pkin(1) * t289;
t351 = pkin(2) * t287;
t350 = pkin(4) * t286;
t274 = t286 * pkin(7);
t281 = pkin(9) + qJ(5);
t269 = sin(t281);
t341 = t256 * t269;
t340 = t257 * t269;
t339 = t284 * t288;
t337 = t286 * MDP(8);
t335 = t287 * t289;
t334 = t291 * t286;
t262 = pkin(2) + t352;
t333 = pkin(1) * t335 + t285 * t262;
t252 = t285 * t353;
t325 = t282 * t339;
t324 = t282 * t334;
t323 = t299 * t284;
t322 = qJD(5) * t324 + t228 * t278;
t273 = cos(t283);
t261 = pkin(2) * t273;
t321 = t257 * pkin(3) + t256 * qJ(4) + t261;
t320 = -pkin(3) - t350;
t272 = sin(t283);
t319 = g(2) * t272 - g(3) * t273;
t318 = t332 * t278;
t208 = t231 * t287 - t245;
t317 = t262 * t287 - t252;
t316 = qJD(1) * (-qJD(2) + t282);
t315 = qJD(2) * (-qJD(1) - t282);
t183 = t320 * t278 + t297;
t308 = qJD(4) - t208;
t198 = t320 * t282 + t308;
t227 = -t334 + t339;
t222 = t227 * qJD(5);
t313 = g(2) * t340 + g(3) * t341 + t183 * t228 - t198 * t222;
t218 = -pkin(3) - t317;
t311 = -g(2) * t273 - g(3) * t272;
t238 = t278 * t334;
t310 = -t278 * t339 + t238;
t260 = pkin(2) * t272;
t309 = t256 * pkin(3) - qJ(4) * t257 + t260;
t192 = -qJD(5) * t325 + t322;
t223 = t228 * qJD(5);
t193 = t282 * t223 - t310;
t201 = -qJD(5) * t222 + qJDD(5) * t228;
t202 = -qJD(5) * t223 - qJDD(5) * t227;
t212 = -t324 + t325;
t307 = (-t192 * t227 - t193 * t228 + t212 * t222 - t214 * t223) * MDP(13) + (t192 * t228 - t214 * t222) * MDP(12) + t201 * MDP(14) + t202 * MDP(15) + t278 * MDP(4);
t217 = qJ(4) + t333;
t206 = (-pkin(7) - t217) * t284;
t207 = t217 * t286 + t274;
t306 = t206 * t291 - t207 * t288;
t305 = t206 * t288 + t207 * t291;
t220 = (t285 * t292 + t335) * t345;
t304 = t218 * t278 + t220 * t282;
t219 = t285 * t328 + t246;
t258 = -pkin(3) - t351;
t303 = -t219 * t282 + t258 * t278;
t254 = pkin(2) * t285 + qJ(4);
t224 = (-pkin(7) - t254) * t284;
t225 = t254 * t286 + t274;
t302 = t224 * t291 - t225 * t288;
t301 = t224 * t288 + t225 * t291;
t300 = t287 * t292 * t345 - qJD(2) * t252;
t298 = t263 + t311;
t296 = -t355 + t356;
t270 = cos(t281);
t294 = t183 * t227 + t198 * t223 + t312 * t270;
t293 = cos(qJ(1));
t290 = sin(qJ(1));
t276 = t293 * pkin(1);
t275 = t290 * pkin(1);
t233 = t320 - t351;
t215 = qJD(4) + t300;
t210 = t218 - t350;
t204 = -pkin(3) * t282 + t308;
t177 = t278 * t274 + t182;
t176 = t266 + (-pkin(7) * t278 - t187) * t284;
t1 = [qJDD(1) * MDP(1) + (-g(2) * t293 - g(3) * t290) * MDP(2) + (g(2) * t290 - g(3) * t293) * MDP(3) + ((t278 * t292 + t289 * t315) * pkin(1) + t298) * MDP(5) + (((-qJDD(1) - t278) * t289 + t292 * t315) * pkin(1) + t319) * MDP(6) + (t195 * t333 + t209 * t300 + t194 * t317 - t208 * t220 - g(2) * (t261 + t276) - g(3) * (t260 + t275)) * MDP(7) + (t299 - t304) * t337 + (t304 * t284 - t323) * MDP(9) + (t215 * t357 + t217 * t318 + t296) * MDP(10) + (t189 * t218 + t204 * t220 - g(2) * (t276 + t321) - g(3) * (t275 + t309) + t356 * t217 + t354 * t215) * MDP(11) + (t220 * t212 + t210 * t193 + t306 * qJDD(5) + (-t305 * qJD(5) - t228 * t215) * qJD(5) + t294) * MDP(17) + (t220 * t214 + t210 * t192 - t305 * qJDD(5) + (-t306 * qJD(5) + t227 * t215) * qJD(5) + t313) * MDP(18) + t307; (t316 * t353 + t298) * MDP(5) + ((t292 * t316 - t330) * pkin(1) + t319) * MDP(6) + (t208 * t219 - t209 * t221 + (t194 * t287 + t195 * t285 + t311) * pkin(2)) * MDP(7) + (t299 - t303) * t337 + (t303 * t284 - t323) * MDP(9) + (t254 * t318 + t331 * t357 + t296) * MDP(10) + (t189 * t258 - t204 * t219 - g(2) * t321 - g(3) * t309 + (t182 * t254 + t331 * t197) * t286 + (-t181 * t254 - t331 * t196) * t284) * MDP(11) + (t233 * t193 + t302 * qJDD(5) - t219 * t212 + (-t301 * qJD(5) - t331 * t228) * qJD(5) + t294) * MDP(17) + (t233 * t192 - t301 * qJDD(5) - t219 * t214 + (-t302 * qJD(5) + t331 * t227) * qJD(5) + t313) * MDP(18) + t307; (qJDD(3) - g(1)) * MDP(7) + (t181 * t286 + t182 * t284 - g(1)) * MDP(11) + t202 * MDP(17) - t201 * MDP(18); (t297 - t312) * MDP(11) - t238 * MDP(17) + t322 * MDP(18) + (-pkin(3) * MDP(11) - t337 + (t288 * MDP(17) + MDP(9)) * t284) * t278 + (0.2e1 * t214 * MDP(17) + (-t212 - t325) * MDP(18)) * qJD(5) + (-MDP(10) * t357 - t354 * MDP(11)) * t282; t214 * t212 * MDP(12) + (-t212 ^ 2 + t214 ^ 2) * MDP(13) + t310 * MDP(15) + qJDD(5) * MDP(16) + (-g(1) * t270 + g(2) * t341 - g(3) * t340 + t291 * t176 - t288 * t177 - t198 * t214) * MDP(17) + (g(1) * t269 - t288 * t176 - t291 * t177 + t198 * t212 + t355 * t270) * MDP(18) + (t322 + (t212 - t325) * qJD(5)) * MDP(14);];
tau = t1;
