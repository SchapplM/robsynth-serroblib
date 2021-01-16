% Calculate vector of inverse dynamics joint torques for
% S5PRRPR3
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
%   pkin=[a2,a3,a4,a5,d2,d3,d5,theta1,theta4]';
% MDP [22x1]
%   Minimal dynamic parameter vector (fixed base model)
%   see S5PRRPR3_convert_par2_MPV_fixb.m
% 
% Output:
% tau [5x1]
%   joint torques of inverse dynamics (contains inertial, gravitational Coriolis and centrifugal forces)

% Quelle: HybrDyn-Toolbox
% Datum: 2021-01-15 15:42
% Revision: 24b2e7d74a0c1a3b64fa2f8f5ad758691ad61af3 (2021-01-15)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tau = S5PRRPR3_invdynJ_fixb_mdp_slag_vp(qJ, qJD, qJDD, g, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(5,1),zeros(3,1),zeros(9,1),zeros(22,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PRRPR3_invdynJ_fixb_mdp_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5PRRPR3_invdynJ_fixb_mdp_slag_vp: qJD has to be [5x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [5 1]), ...
  'S5PRRPR3_invdynJ_fixb_mdp_slag_vp: qJDD has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5PRRPR3_invdynJ_fixb_mdp_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5PRRPR3_invdynJ_fixb_mdp_slag_vp: pkin has to be [9x1] (double)');
assert(isreal(MDP) && all(size(MDP) == [22 1]), ...
  'S5PRRPR3_invdynJ_fixb_mdp_slag_vp: MDP has to be [22x1] (double)'); 

%% Symbolic Calculation
% From invdyn_joint_fixb_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2021-01-15 15:42:15
% EndTime: 2021-01-15 15:42:21
% DurationCPUTime: 2.29s
% Computational Cost: add. (1451->275), mult. (3167->364), div. (0->0), fcn. (2285->12), ass. (0->128)
t288 = sin(pkin(9));
t289 = cos(pkin(9));
t294 = cos(qJ(3));
t326 = qJD(2) * t294;
t317 = t289 * t326;
t292 = sin(qJ(3));
t327 = qJD(2) * t292;
t244 = t288 * t327 - t317;
t293 = cos(qJ(5));
t235 = t293 * t244;
t251 = t288 * t294 + t289 * t292;
t247 = t251 * qJD(2);
t291 = sin(qJ(5));
t333 = t247 * t291;
t209 = -t235 - t333;
t284 = qJD(3) + qJD(5);
t334 = t209 * t284;
t305 = t244 * t291 - t293 * t247;
t335 = t305 * t284;
t283 = pkin(8) + qJ(2);
t275 = sin(t283);
t277 = cos(t283);
t311 = g(1) * t277 + g(2) * t275;
t279 = t294 * qJDD(1);
t337 = qJ(4) + pkin(6);
t314 = t337 * qJD(3);
t312 = qJD(2) * t314;
t345 = qJD(3) * qJD(1) + qJD(2) * qJD(4) + t337 * qJDD(2);
t196 = qJDD(3) * pkin(3) - t345 * t292 - t294 * t312 + t279;
t199 = (qJDD(1) - t312) * t292 + t345 * t294;
t181 = t289 * t196 - t199 * t288;
t323 = qJD(2) * qJD(3);
t316 = t292 * t323;
t264 = t288 * t316;
t315 = t294 * t323;
t215 = t251 * qJDD(2) + t289 * t315 - t264;
t176 = qJDD(3) * pkin(4) - pkin(7) * t215 + t181;
t182 = t288 * t196 + t289 * t199;
t246 = t251 * qJD(3);
t321 = qJDD(2) * t294;
t265 = t289 * t321;
t322 = qJDD(2) * t292;
t214 = qJD(2) * t246 + t288 * t322 - t265;
t177 = -pkin(7) * t214 + t182;
t262 = t337 * t292;
t239 = t294 * qJD(1) - qJD(2) * t262;
t336 = qJD(3) * pkin(3);
t234 = t239 + t336;
t263 = t337 * t294;
t241 = qJD(1) * t292 + qJD(2) * t263;
t331 = t289 * t241;
t198 = t288 * t234 + t331;
t342 = pkin(7) * t244;
t187 = t198 - t342;
t274 = pkin(3) * t294 + pkin(2);
t259 = -t274 * qJD(2) + qJD(4);
t218 = pkin(4) * t244 + t259;
t285 = qJ(3) + pkin(9);
t280 = qJ(5) + t285;
t270 = sin(t280);
t271 = cos(t280);
t325 = qJD(5) * t291;
t348 = g(3) * t270 - t291 * t176 - t293 * t177 + t187 * t325 - t218 * t209 + t311 * t271;
t347 = -g(3) * t271 + t293 * t176 - t291 * t177 + t218 * t305 + t311 * t270;
t282 = qJDD(3) + qJDD(5);
t346 = t282 * MDP(20) + t209 * t305 * MDP(16) + (-t209 ^ 2 + t305 ^ 2) * MDP(17);
t313 = t293 * t214 + t215 * t291;
t180 = -t305 * qJD(5) + t313;
t344 = pkin(3) * t288;
t343 = pkin(3) * t292;
t341 = pkin(7) * t247;
t338 = g(3) * t294;
t227 = t288 * t241;
t332 = t288 * t292;
t197 = t289 * t234 - t227;
t186 = qJD(3) * pkin(4) + t197 - t341;
t330 = t293 * t186;
t329 = qJDD(1) - g(3);
t202 = t289 * t239 - t227;
t240 = qJD(4) * t294 - t292 * t314;
t242 = -qJD(4) * t292 - t294 * t314;
t203 = t289 * t240 + t288 * t242;
t220 = -t288 * t262 + t289 * t263;
t286 = t292 ^ 2;
t328 = -t294 ^ 2 + t286;
t320 = pkin(3) * t316 + qJDD(4);
t319 = t292 * t336;
t318 = -qJD(5) * t235 - t291 * t214 + t293 * t215;
t200 = -t239 * t288 - t331;
t201 = -t240 * t288 + t289 * t242;
t219 = -t289 * t262 - t263 * t288;
t310 = g(1) * t275 - g(2) * t277;
t250 = -t289 * t294 + t332;
t216 = t293 * t250 + t251 * t291;
t249 = t250 * qJD(3);
t183 = -t216 * qJD(5) - t246 * t291 - t249 * t293;
t217 = -t250 * t291 + t251 * t293;
t309 = t183 * t284 + t217 * t282;
t308 = -t291 * t186 - t293 * t187;
t204 = -pkin(7) * t251 + t219;
t205 = -pkin(7) * t250 + t220;
t307 = t204 * t293 - t205 * t291;
t306 = t204 * t291 + t205 * t293;
t272 = pkin(3) * t289 + pkin(4);
t303 = t272 * t291 + t293 * t344;
t302 = t272 * t293 - t291 * t344;
t301 = -0.2e1 * pkin(2) * t323 - pkin(6) * qJDD(3);
t179 = -t247 * t325 + t318;
t238 = -t274 * qJDD(2) + t320;
t295 = qJD(3) ^ 2;
t298 = 0.2e1 * qJDD(2) * pkin(2) - pkin(6) * t295 + t310;
t296 = qJD(2) ^ 2;
t297 = pkin(2) * t296 - pkin(6) * qJDD(2) + t311;
t278 = cos(t285);
t276 = sin(t285);
t261 = qJDD(3) * t294 - t292 * t295;
t260 = qJDD(3) * t292 + t294 * t295;
t226 = pkin(4) * t250 - t274;
t222 = pkin(4) * t246 + t319;
t221 = pkin(3) * t327 + pkin(4) * t247;
t192 = pkin(4) * t214 + t238;
t191 = -pkin(7) * t246 + t203;
t190 = t202 - t341;
t189 = pkin(7) * t249 + t201;
t188 = t200 + t342;
t184 = t217 * qJD(5) + t293 * t246 - t249 * t291;
t178 = -t184 * t284 - t216 * t282;
t1 = [t329 * MDP(1) + t261 * MDP(10) - t260 * MDP(11) + (-qJD(3) * t246 - qJDD(3) * t250) * MDP(12) + (qJD(3) * t249 - qJDD(3) * t251) * MDP(13) + (-t214 * t251 + t215 * t250 + t244 * t249 + t246 * t247) * MDP(14) + (-t181 * t250 + t182 * t251 - t197 * t246 - t198 * t249 - g(3)) * MDP(15) + t178 * MDP(21) - t309 * MDP(22); qJDD(2) * MDP(2) + t310 * MDP(3) + t311 * MDP(4) + (qJDD(2) * t286 + 0.2e1 * t292 * t315) * MDP(5) + 0.2e1 * (t292 * t321 - t328 * t323) * MDP(6) + t260 * MDP(7) + t261 * MDP(8) + (t301 * t292 + t298 * t294) * MDP(10) + (-t298 * t292 + t301 * t294) * MDP(11) + (qJDD(3) * t219 - t214 * t274 + t238 * t250 + t246 * t259 + t310 * t278 + (t244 * t343 + t201) * qJD(3)) * MDP(12) + (-qJDD(3) * t220 - t215 * t274 + t238 * t251 - t249 * t259 - t310 * t276 + (t247 * t343 - t203) * qJD(3)) * MDP(13) + (-t181 * t251 - t182 * t250 + t197 * t249 - t198 * t246 - t201 * t247 - t203 * t244 - t214 * t220 - t215 * t219 - t311) * MDP(14) + (t182 * t220 + t198 * t203 + t181 * t219 + t197 * t201 - t238 * t274 + t259 * t319 - g(1) * (-t274 * t275 + t277 * t337) - g(2) * (t274 * t277 + t275 * t337)) * MDP(15) + (t179 * t217 - t183 * t305) * MDP(16) + (-t179 * t216 - t180 * t217 + t183 * t209 + t184 * t305) * MDP(17) + t309 * MDP(18) + t178 * MDP(19) + (-t222 * t209 + t226 * t180 + t192 * t216 + t218 * t184 + (-t306 * qJD(5) + t189 * t293 - t191 * t291) * t284 + t307 * t282 + t310 * t271) * MDP(21) + (-t222 * t305 + t226 * t179 + t192 * t217 + t218 * t183 - (t307 * qJD(5) + t189 * t291 + t191 * t293) * t284 - t306 * t282 - t310 * t270) * MDP(22); MDP(7) * t322 + MDP(8) * t321 + qJDD(3) * MDP(9) + (t297 * t292 + t279 - t338) * MDP(10) + (-t329 * t292 + t297 * t294) * MDP(11) + (-g(3) * t278 - qJD(3) * t200 - t247 * t259 + t311 * t276 + (qJDD(3) * t289 - t244 * t327) * pkin(3) + t181) * MDP(12) + (g(3) * t276 + qJD(3) * t202 + t244 * t259 + t311 * t278 + (-qJDD(3) * t288 - t247 * t327) * pkin(3) - t182) * MDP(13) + ((t198 + t200) * t247 + (-t197 + t202) * t244 + (-t214 * t288 - t215 * t289) * pkin(3)) * MDP(14) + (-t197 * t200 - t198 * t202 + (-t338 + t181 * t289 + t182 * t288 + (-qJD(2) * t259 + t311) * t292) * pkin(3)) * MDP(15) + (t179 - t334) * MDP(18) + (-t180 - t335) * MDP(19) + (t302 * t282 + t221 * t209 - (t188 * t293 - t190 * t291) * t284 + (-t303 * t284 + t308) * qJD(5) + t347) * MDP(21) + (-t303 * t282 + t221 * t305 + (t188 * t291 + t190 * t293) * t284 + (-t302 * t284 - t330) * qJD(5) + t348) * MDP(22) + (-t292 * t294 * MDP(5) + t328 * MDP(6)) * t296 + t346; -t265 * MDP(12) - t264 * MDP(13) + (-t244 ^ 2 - t247 ^ 2) * MDP(14) + (t197 * t247 + t198 * t244 - t310 + t320) * MDP(15) + (t180 - t335) * MDP(21) + (t179 + t334) * MDP(22) + (MDP(12) * t332 + t251 * MDP(13) - t274 * MDP(15)) * qJDD(2) + ((t288 * t326 + t289 * t327 + t247) * MDP(12) + (-t244 + t317) * MDP(13)) * qJD(3); (t318 - t334) * MDP(18) + (-t313 - t335) * MDP(19) + (-t308 * t284 + t347) * MDP(21) + ((-t187 * t291 + t330) * t284 + t348) * MDP(22) + (-MDP(18) * t333 + t305 * MDP(19) + t308 * MDP(21) - MDP(22) * t330) * qJD(5) + t346;];
tau = t1;
