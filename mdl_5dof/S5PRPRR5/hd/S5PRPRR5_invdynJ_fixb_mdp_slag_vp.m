% Calculate vector of inverse dynamics joint torques for
% S5PRPRR5
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
%   pkin=[a2,a3,a4,a5,d2,d4,d5,theta1,theta3]';
% MDP [22x1]
%   Minimal dynamic parameter vector (fixed base model)
%   see S5PRPRR5_convert_par2_MPV_fixb.m
% 
% Output:
% tau [5x1]
%   joint torques of inverse dynamics (contains inertial, gravitational Coriolis and centrifugal forces)

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-05 15:55
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tau = S5PRPRR5_invdynJ_fixb_mdp_slag_vp(qJ, qJD, qJDD, g, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(5,1),zeros(3,1),zeros(9,1),zeros(22,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PRPRR5_invdynJ_fixb_mdp_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5PRPRR5_invdynJ_fixb_mdp_slag_vp: qJD has to be [5x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [5 1]), ...
  'S5PRPRR5_invdynJ_fixb_mdp_slag_vp: qJDD has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5PRPRR5_invdynJ_fixb_mdp_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5PRPRR5_invdynJ_fixb_mdp_slag_vp: pkin has to be [9x1] (double)');
assert(isreal(MDP) && all(size(MDP) == [22 1]), ...
  'S5PRPRR5_invdynJ_fixb_mdp_slag_vp: MDP has to be [22x1] (double)'); 

%% Symbolic Calculation
% From invdyn_joint_fixb_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 15:54:49
% EndTime: 2019-12-05 15:54:54
% DurationCPUTime: 2.83s
% Computational Cost: add. (1297->271), mult. (2915->358), div. (0->0), fcn. (2313->14), ass. (0->130)
t286 = cos(qJ(2));
t336 = qJD(1) * t286;
t316 = qJD(3) - t336;
t279 = cos(pkin(9));
t285 = cos(qJ(4));
t333 = qJD(2) * t285;
t323 = t279 * t333;
t277 = sin(pkin(9));
t282 = sin(qJ(4));
t334 = qJD(2) * t282;
t324 = t277 * t334;
t239 = -t323 + t324;
t284 = cos(qJ(5));
t281 = sin(qJ(5));
t250 = t277 * t285 + t279 * t282;
t362 = t250 * qJD(2);
t348 = t362 * t281;
t204 = t284 * t239 + t348;
t276 = qJD(4) + qJD(5);
t350 = t204 * t276;
t304 = t239 * t281 - t284 * t362;
t351 = t304 * t276;
t283 = sin(qJ(2));
t278 = sin(pkin(8));
t280 = cos(pkin(8));
t315 = g(1) * t280 + g(2) * t278;
t300 = t315 * t283;
t356 = g(3) * t286;
t293 = t300 - t356;
t345 = t279 * MDP(5);
t370 = -t277 * MDP(6) + t345;
t337 = qJD(1) * t283;
t257 = qJD(2) * qJ(3) + t337;
t322 = pkin(6) * qJD(2) + t257;
t232 = t322 * t277;
t233 = t322 * t279;
t307 = t232 * t282 - t233 * t285;
t197 = -pkin(7) * t239 - t307;
t266 = -pkin(3) * t279 - pkin(2);
t244 = t266 * qJD(2) + t316;
t215 = pkin(4) * t239 + t244;
t275 = pkin(9) + qJ(4);
t271 = qJ(5) + t275;
t264 = sin(t271);
t265 = cos(t271);
t332 = qJD(5) * t281;
t344 = t280 * t286;
t346 = t278 * t286;
t357 = g(3) * t283;
t369 = t215 * t204 - g(1) * (-t264 * t278 - t265 * t344) - g(2) * (t264 * t280 - t265 * t346) + t197 * t332 + t265 * t357;
t327 = qJDD(2) * t285;
t328 = qJDD(2) * t282;
t325 = qJD(4) * t323 + t277 * t327 + t279 * t328;
t211 = -qJD(4) * t324 + t325;
t329 = qJDD(2) * qJ(3);
t234 = t329 + qJDD(1) * t283 + (qJD(3) + t336) * qJD(2);
t320 = pkin(6) * qJDD(2) + t234;
t217 = t320 * t277;
t218 = t320 * t279;
t318 = -t285 * t217 - t282 * t218;
t186 = qJDD(4) * pkin(4) - pkin(7) * t211 + t307 * qJD(4) + t318;
t243 = t250 * qJD(4);
t260 = t279 * t327;
t312 = -t277 * t328 + t260;
t212 = qJD(2) * t243 - t312;
t308 = -t282 * t217 + t285 * t218;
t364 = -t285 * t232 - t233 * t282;
t187 = -pkin(7) * t212 + t364 * qJD(4) + t308;
t368 = t215 * t304 - g(1) * (-t264 * t344 + t265 * t278) - g(2) * (-t264 * t346 - t265 * t280) + t284 * t186 - t281 * t187 + t264 * t357;
t272 = qJDD(4) + qJDD(5);
t367 = t272 * MDP(20) + (-t204 ^ 2 + t304 ^ 2) * MDP(17) - t204 * MDP(16) * t304;
t366 = t315 * t286;
t339 = t277 ^ 2 + t279 ^ 2;
t365 = t339 * MDP(7);
t355 = pkin(6) + qJ(3);
t253 = t355 * t277;
t254 = t355 * t279;
t340 = -t282 * t253 + t285 * t254;
t363 = -t340 * qJD(4) - t316 * t250;
t343 = t285 * t279;
t249 = t277 * t282 - t343;
t242 = t249 * qJD(4);
t319 = t211 * t281 + t284 * t212;
t189 = -t304 * qJD(5) + t319;
t247 = t285 * t253;
t297 = t249 * t286;
t360 = -qJD(1) * t297 + (qJD(3) * t277 + qJD(4) * t254) * t282 - qJD(3) * t343 + qJD(4) * t247;
t354 = qJDD(2) * pkin(2);
t196 = -pkin(7) * t362 + t364;
t193 = qJD(4) * pkin(4) + t196;
t353 = t193 * t284;
t352 = t197 * t284;
t342 = qJDD(1) - g(3);
t331 = qJD(5) * t284;
t330 = qJD(1) * qJD(2);
t326 = t284 * t211 - t281 * t212 - t239 * t331;
t321 = t339 * t234;
t317 = -t254 * t282 - t247;
t199 = -pkin(7) * t250 + t317;
t314 = pkin(7) * t243 - qJD(5) * t199 + t360;
t200 = -pkin(7) * t249 + t340;
t313 = -pkin(7) * t242 + qJD(5) * t200 - t363;
t311 = t243 * pkin(4) - t337;
t310 = -MDP(4) + t365;
t309 = -t193 * t281 - t352;
t235 = t250 * t283;
t236 = t249 * t283;
t306 = -t235 * t284 + t236 * t281;
t305 = -t235 * t281 - t236 * t284;
t213 = t284 * t249 + t250 * t281;
t214 = -t249 * t281 + t250 * t284;
t303 = -qJDD(1) * t286 + t283 * t330 + qJDD(3);
t299 = MDP(3) + t370;
t188 = -t332 * t362 + t326;
t295 = t339 * qJD(2) * t257;
t290 = t316 * t339;
t223 = t266 * qJDD(2) + t303;
t289 = t321 - t357 - t366;
t287 = qJD(2) ^ 2;
t270 = cos(t275);
t269 = sin(t275);
t255 = -qJD(2) * pkin(2) + t316;
t238 = t303 - t354;
t227 = pkin(4) * t249 + t266;
t202 = t283 * t242 - t286 * t362;
t201 = -qJD(2) * t297 - qJD(4) * t235;
t198 = pkin(4) * t212 + t223;
t191 = t214 * qJD(5) - t242 * t281 + t284 * t243;
t190 = -t213 * qJD(5) - t242 * t284 - t243 * t281;
t1 = [t342 * MDP(1) - g(3) * MDP(8) + (qJD(4) * t202 - qJDD(4) * t235) * MDP(14) + (-qJD(4) * t201 + qJDD(4) * t236) * MDP(15) + ((-t305 * qJD(5) - t201 * t281 + t202 * t284) * t276 + t306 * t272) * MDP(21) + (-(t306 * qJD(5) + t201 * t284 + t202 * t281) * t276 - t305 * t272) * MDP(22) + ((-t238 + t295) * MDP(8) - t212 * MDP(14) - t211 * MDP(15) - t189 * MDP(21) - t188 * MDP(22) + t310 * t287 + t299 * qJDD(2)) * t286 + (MDP(8) * t321 - t299 * t287 + t310 * qJDD(2) + (t239 * MDP(14) + MDP(15) * t362 + t204 * MDP(21) - MDP(22) * t304 + t255 * MDP(8)) * qJD(2)) * t283; qJDD(2) * MDP(2) + (t342 * t286 + t300) * MDP(3) + (-t342 * t283 + t366) * MDP(4) + (t290 * qJD(2) + t339 * t329 + t289) * MDP(7) + (-t255 * t337 + (-t238 + t293) * pkin(2) + t289 * qJ(3) + t290 * t257) * MDP(8) + (t211 * t250 - t242 * t362) * MDP(9) + (-t211 * t249 - t212 * t250 + t239 * t242 - t243 * t362) * MDP(10) + (-qJD(4) * t242 + qJDD(4) * t250) * MDP(11) + (-qJD(4) * t243 - qJDD(4) * t249) * MDP(12) + (t363 * qJD(4) + t317 * qJDD(4) + t266 * t212 + t223 * t249 - t239 * t337 + t244 * t243 + t293 * t270) * MDP(14) + (t360 * qJD(4) - t340 * qJDD(4) + t266 * t211 + t223 * t250 - t244 * t242 - t269 * t293 - t337 * t362) * MDP(15) + (t188 * t214 - t190 * t304) * MDP(16) + (-t188 * t213 - t189 * t214 - t190 * t204 + t191 * t304) * MDP(17) + (t190 * t276 + t214 * t272) * MDP(18) + (-t191 * t276 - t213 * t272) * MDP(19) + ((t199 * t284 - t200 * t281) * t272 + t227 * t189 + t198 * t213 + t215 * t191 + (t314 * t281 - t313 * t284) * t276 + t311 * t204 + t293 * t265) * MDP(21) + (-(t199 * t281 + t200 * t284) * t272 + t227 * t188 + t198 * t214 + t215 * t190 + (t313 * t281 + t314 * t284) * t276 - t311 * t304 - t293 * t264) * MDP(22) + t370 * ((t315 + t330) * t283 - t238 + t354 - t356); (-t295 + t303 - t293) * MDP(8) - t260 * MDP(14) + t325 * MDP(15) + (t189 - t351) * MDP(21) + (t188 - t350) * MDP(22) - t287 * t365 + (-t345 - pkin(2) * MDP(8) + (MDP(14) * t282 + MDP(6)) * t277) * qJDD(2) + ((t277 * t333 + t279 * t334 + t362) * MDP(14) + (-t239 - t324) * MDP(15)) * qJD(4); t362 * t239 * MDP(9) + (-t239 ^ 2 + t362 ^ 2) * MDP(10) + ((t239 - t324) * qJD(4) + t325) * MDP(11) + t312 * MDP(12) + qJDD(4) * MDP(13) + (-t244 * t362 - g(1) * (-t269 * t344 + t270 * t278) - g(2) * (-t269 * t346 - t270 * t280) + t269 * t357 + t318) * MDP(14) + (t244 * t239 - g(1) * (-t269 * t278 - t270 * t344) - g(2) * (t269 * t280 - t270 * t346) + t270 * t357 - t308) * MDP(15) + (t188 + t350) * MDP(18) + (-t189 - t351) * MDP(19) + (-(-t196 * t281 - t352) * t276 + t309 * qJD(5) + (-t204 * t362 + t284 * t272 - t276 * t332) * pkin(4) + t368) * MDP(21) + ((-t197 * t276 - t186) * t281 + (-qJD(5) * t193 + t196 * t276 - t187) * t284 + (-t281 * t272 - t276 * t331 + t304 * t362) * pkin(4) + t369) * MDP(22) + t367; (t326 + t350) * MDP(18) + (-t319 - t351) * MDP(19) + (-t309 * t276 + t368) * MDP(21) + (-t284 * t187 - t281 * t186 + (-t197 * t281 + t353) * t276 + t369) * MDP(22) + (-MDP(18) * t348 + t304 * MDP(19) + t309 * MDP(21) - MDP(22) * t353) * qJD(5) + t367;];
tau = t1;
