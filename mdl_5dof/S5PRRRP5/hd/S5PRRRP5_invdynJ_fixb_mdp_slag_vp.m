% Calculate vector of inverse dynamics joint torques for
% S5PRRRP5
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
%   pkin=[a2,a3,a4,a5,d2,d3,d4,theta1]';
% MDP [20x1]
%   Minimal dynamic parameter vector (fixed base model)
%   see S5PRRRP5_convert_par2_MPV_fixb.m
% 
% Output:
% tau [5x1]
%   joint torques of inverse dynamics (contains inertial, gravitational Coriolis and centrifugal forces)

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-05 16:49
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tau = S5PRRRP5_invdynJ_fixb_mdp_slag_vp(qJ, qJD, qJDD, g, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(5,1),zeros(3,1),zeros(8,1),zeros(20,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PRRRP5_invdynJ_fixb_mdp_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5PRRRP5_invdynJ_fixb_mdp_slag_vp: qJD has to be [5x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [5 1]), ...
  'S5PRRRP5_invdynJ_fixb_mdp_slag_vp: qJDD has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5PRRRP5_invdynJ_fixb_mdp_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5PRRRP5_invdynJ_fixb_mdp_slag_vp: pkin has to be [8x1] (double)');
assert(isreal(MDP) && all(size(MDP) == [20 1]), ...
  'S5PRRRP5_invdynJ_fixb_mdp_slag_vp: MDP has to be [20x1] (double)'); 

%% Symbolic Calculation
% From invdyn_joint_fixb_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 16:49:26
% EndTime: 2019-12-05 16:49:29
% DurationCPUTime: 2.33s
% Computational Cost: add. (1360->264), mult. (3005->343), div. (0->0), fcn. (2092->10), ass. (0->134)
t299 = cos(qJ(3));
t297 = sin(qJ(2));
t353 = qJD(1) * t297;
t269 = qJD(2) * pkin(6) + t353;
t335 = pkin(7) * qJD(2) + t269;
t237 = t335 * t299;
t295 = sin(qJ(4));
t226 = t295 * t237;
t296 = sin(qJ(3));
t236 = t335 * t296;
t229 = qJD(3) * pkin(3) - t236;
t298 = cos(qJ(4));
t333 = t298 * t229 - t226;
t254 = t295 * t299 + t296 * t298;
t244 = t254 * qJD(2);
t366 = t244 * qJ(5);
t392 = t366 - t333;
t342 = qJDD(2) * t299;
t343 = qJDD(2) * t296;
t289 = qJD(3) + qJD(4);
t391 = t289 * t254;
t214 = qJD(2) * t391 + t295 * t343 - t298 * t342;
t300 = cos(qJ(2));
t293 = sin(pkin(8));
t294 = cos(pkin(8));
t328 = g(1) * t294 + g(2) * t293;
t318 = t328 * t300;
t373 = g(3) * t297;
t390 = t318 + t373;
t345 = qJD(2) * qJD(3);
t336 = t299 * t345;
t389 = -t336 - t343;
t372 = g(3) * t300;
t386 = t328 * t297;
t388 = t386 - t372;
t360 = qJDD(1) - g(3);
t387 = t360 * t300 + t386;
t286 = t299 * pkin(3);
t371 = pkin(2) + t286;
t362 = t298 * t299;
t363 = t295 * t296;
t253 = -t362 + t363;
t239 = t253 * t297;
t378 = pkin(6) + pkin(7);
t340 = qJD(3) * t378;
t259 = t296 * t340;
t260 = t299 * t340;
t352 = qJD(1) * t300;
t263 = t378 * t296;
t264 = t378 * t299;
t355 = -t295 * t263 + t298 * t264;
t384 = -t355 * qJD(4) + t254 * t352 + t259 * t295 - t298 * t260;
t316 = t300 * t253;
t347 = qJD(4) * t298;
t348 = qJD(4) * t295;
t383 = -qJD(1) * t316 + t298 * t259 + t295 * t260 + t263 * t347 + t264 * t348;
t292 = qJ(3) + qJ(4);
t284 = sin(t292);
t285 = cos(t292);
t364 = t294 * t300;
t365 = t293 * t300;
t382 = t284 * t373 - g(2) * (-t284 * t365 - t285 * t294) - g(1) * (-t284 * t364 + t285 * t293);
t346 = qJD(1) * qJD(2);
t278 = t297 * t346;
t301 = qJD(3) ^ 2;
t344 = qJDD(1) * t300;
t381 = 0.2e1 * qJDD(2) * pkin(2) - pkin(6) * t301 + t297 * (t328 + t346) - t278 + t344 - t372;
t250 = qJDD(2) * pkin(6) + qJDD(1) * t297 + t300 * t346;
t349 = qJD(3) * t299;
t216 = qJDD(3) * pkin(3) + pkin(7) * t389 - t296 * t250 - t269 * t349;
t337 = t296 * t345;
t350 = qJD(3) * t296;
t217 = -t269 * t350 + t299 * t250 + (-t337 + t342) * pkin(7);
t380 = -(qJD(4) * t229 + t217) * t298 - t295 * t216 + t237 * t348;
t379 = t244 ^ 2;
t370 = qJD(2) * pkin(2);
t338 = qJD(2) * t362;
t351 = qJD(2) * t296;
t339 = t295 * t351;
t242 = -t338 + t339;
t369 = qJ(5) * t242;
t251 = -qJD(2) * t371 - t352;
t220 = pkin(4) * t242 + qJD(5) + t251;
t367 = t220 * t244;
t228 = t298 * t237;
t302 = qJD(2) ^ 2;
t361 = t299 * t302;
t359 = -qJ(5) * t391 - qJD(5) * t253 - t383;
t326 = t289 * t363;
t221 = -t298 * t349 - t299 * t347 + t326;
t358 = qJ(5) * t221 - qJD(5) * t254 + t384;
t205 = pkin(4) * t289 - t392;
t357 = t205 + t392;
t356 = -t298 * t236 - t226;
t262 = pkin(4) * t285 + t286;
t290 = t296 ^ 2;
t354 = -t299 ^ 2 + t290;
t332 = t236 * t295 - t228;
t331 = -t298 * t263 - t264 * t295;
t329 = -qJD(4) * t338 - t295 * t342 + t298 * t389;
t327 = g(1) * t293 - g(2) * t294;
t324 = -t229 * t295 - t228;
t322 = qJDD(3) * t299 - t296 * t301;
t321 = qJDD(3) * t296 + t299 * t301;
t317 = pkin(3) * t350 - t353;
t314 = pkin(3) * t337 - qJDD(2) * t371 + t278;
t241 = t242 ^ 2;
t287 = qJDD(3) + qJDD(4);
t312 = t244 * t242 * MDP(12) + (-t329 + (t242 - t339) * t289) * MDP(14) + (t244 * t289 - t214) * MDP(15) + (-t241 + t379) * MDP(13) + t287 * MDP(16);
t270 = -t352 - t370;
t310 = -pkin(6) * qJDD(3) + (t270 + t352 - t370) * qJD(3);
t309 = t324 * qJD(4) + t298 * t216 - t295 * t217;
t307 = pkin(4) * t214 + qJDD(5) + t314;
t306 = -t270 * qJD(2) - t250 + t390;
t304 = -g(1) * (-t284 * t293 - t285 * t364) - g(2) * (t284 * t294 - t285 * t365) + t251 * t242 + t285 * t373 + t380;
t303 = -t251 * t244 + t309 + t382;
t288 = -qJ(5) - t378;
t282 = pkin(3) * t298 + pkin(4);
t261 = -pkin(3) * t296 - pkin(4) * t284;
t258 = pkin(2) + t262;
t238 = t254 * t297;
t223 = t314 - t344;
t219 = -qJ(5) * t253 + t355;
t218 = -qJ(5) * t254 + t331;
t213 = qJD(2) * t326 + t329;
t211 = t289 * t239 - t300 * t244;
t210 = -qJD(2) * t316 - t297 * t391;
t209 = -t366 + t356;
t208 = t332 + t369;
t207 = -t324 - t369;
t202 = t307 - t344;
t199 = -qJ(5) * t214 - qJD(5) * t242 - t380;
t198 = pkin(4) * t287 + qJ(5) * t213 - qJD(5) * t244 + t309;
t1 = [t360 * MDP(1) + (t211 * t289 - t238 * t287) * MDP(17) + (-t210 * t289 + t239 * t287) * MDP(18) + (-t210 * t242 - t211 * t244 - t213 * t238 + t214 * t239) * MDP(19) + (-t198 * t238 - t199 * t239 + t205 * t211 + t207 * t210 - g(3)) * MDP(20) + (qJDD(2) * MDP(3) - t302 * MDP(4) + (-0.2e1 * t337 + t342) * MDP(10) + (-0.2e1 * t336 - t343) * MDP(11) - t214 * MDP(17) + t213 * MDP(18) - t202 * MDP(20)) * t300 + (-t302 * MDP(3) - qJDD(2) * MDP(4) + (-t321 - t361) * MDP(10) + (t296 * t302 - t322) * MDP(11) + (t242 * MDP(17) + t244 * MDP(18) + t220 * MDP(20)) * qJD(2)) * t297; qJDD(2) * MDP(2) + t387 * MDP(3) + (-t360 * t297 + t318) * MDP(4) + (qJDD(2) * t290 + 0.2e1 * t296 * t336) * MDP(5) + 0.2e1 * (t296 * t342 - t354 * t345) * MDP(6) + t321 * MDP(7) + t322 * MDP(8) + (t310 * t296 + t299 * t381) * MDP(10) + (-t296 * t381 + t310 * t299) * MDP(11) + (-t213 * t254 - t221 * t244) * MDP(12) + (t213 * t253 - t214 * t254 + t221 * t242 - t244 * t391) * MDP(13) + (-t221 * t289 + t254 * t287) * MDP(14) + (-t253 * t287 - t289 * t391) * MDP(15) + (-t371 * t214 + t223 * t253 + t317 * t242 + t251 * t391 + t285 * t388 + t331 * t287 + t384 * t289) * MDP(17) + (t371 * t213 - t251 * t221 + t223 * t254 + t317 * t244 - t388 * t284 - t355 * t287 + t383 * t289) * MDP(18) + (-t198 * t254 - t199 * t253 + t205 * t221 - t207 * t391 + t213 * t218 - t214 * t219 - t359 * t242 - t358 * t244 - t390) * MDP(19) + (t199 * t219 + t198 * t218 + t202 * (pkin(4) * t253 - t371) - g(3) * (t258 * t300 - t288 * t297) + (pkin(4) * t391 + t317) * t220 + t359 * t207 + t358 * t205 + t328 * (t258 * t297 + t288 * t300)) * MDP(20); -t296 * MDP(5) * t361 + t354 * MDP(6) * t302 + MDP(7) * t343 + MDP(8) * t342 + qJDD(3) * MDP(9) + (t306 * t296 - t327 * t299) * MDP(10) + (t327 * t296 + t306 * t299) * MDP(11) + (-t332 * t289 + (-t242 * t351 + t287 * t298 - t289 * t348) * pkin(3) + t303) * MDP(17) + (t356 * t289 + (-t244 * t351 - t295 * t287 - t289 * t347) * pkin(3) + t304) * MDP(18) + (t213 * t282 + (t207 + t208) * t244 + (-t205 + t209) * t242 + (-t214 * t295 + (-t242 * t298 + t244 * t295) * qJD(4)) * pkin(3)) * MDP(19) + (t198 * t282 - t207 * t209 - t205 * t208 - pkin(4) * t367 - g(1) * (t261 * t364 + t262 * t293) - g(2) * (t261 * t365 - t262 * t294) - t261 * t373 + (-t220 * t351 + t199 * t295 + (-t205 * t295 + t207 * t298) * qJD(4)) * pkin(3)) * MDP(20) + t312; (-t324 * t289 + t303) * MDP(17) + (t333 * t289 + t304) * MDP(18) + (pkin(4) * t213 - t357 * t242) * MDP(19) + (t357 * t207 + (t198 - t367 + t382) * pkin(4)) * MDP(20) + t312; (-t241 - t379) * MDP(19) + (t205 * t244 + t207 * t242 + t307 - t387) * MDP(20);];
tau = t1;
