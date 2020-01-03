% Calculate vector of inverse dynamics joint torques for
% S5RPRRP13
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
%   pkin=[a2,a3,a4,a5,d1,d3,d4]';
% MDP [24x1]
%   Minimal dynamic parameter vector (fixed base model)
%   see S5RPRRP13_convert_par2_MPV_fixb.m
% 
% Output:
% tau [5x1]
%   joint torques of inverse dynamics (contains inertial, gravitational Coriolis and centrifugal forces)

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 19:00
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tau = S5RPRRP13_invdynJ_fixb_mdp_slag_vp(qJ, qJD, qJDD, g, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(5,1),zeros(3,1),zeros(7,1),zeros(24,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPRRP13_invdynJ_fixb_mdp_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RPRRP13_invdynJ_fixb_mdp_slag_vp: qJD has to be [5x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [5 1]), ...
  'S5RPRRP13_invdynJ_fixb_mdp_slag_vp: qJDD has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RPRRP13_invdynJ_fixb_mdp_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [7 1]), ...
  'S5RPRRP13_invdynJ_fixb_mdp_slag_vp: pkin has to be [7x1] (double)');
assert(isreal(MDP) && all(size(MDP) == [24 1]), ...
  'S5RPRRP13_invdynJ_fixb_mdp_slag_vp: MDP has to be [24x1] (double)'); 

%% Symbolic Calculation
% From invdyn_joint_fixb_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 18:59:45
% EndTime: 2019-12-31 18:59:50
% DurationCPUTime: 3.30s
% Computational Cost: add. (2003->368), mult. (3729->459), div. (0->0), fcn. (2150->6), ass. (0->149)
t277 = cos(qJ(3));
t274 = sin(qJ(3));
t374 = g(3) * t274;
t278 = cos(qJ(1));
t275 = sin(qJ(1));
t375 = g(1) * t275;
t383 = g(2) * t278 - t375;
t384 = t383 * t277 + t374;
t258 = qJD(1) * t274 + qJD(4);
t376 = pkin(7) * t277;
t311 = pkin(3) * t274 - t376;
t249 = qJ(2) + t311;
t273 = sin(qJ(4));
t279 = -pkin(1) - pkin(6);
t276 = cos(qJ(4));
t353 = t274 * t276;
t344 = t273 * t249 + t279 * t353;
t326 = qJDD(1) * t277;
t329 = qJD(1) * qJD(3);
t337 = qJD(3) * t273;
t339 = qJD(1) * t277;
t244 = t276 * t339 + t337;
t334 = qJD(4) * t244;
t356 = t273 * t274;
t208 = -t276 * qJDD(3) + t273 * t326 - t329 * t356 + t334;
t281 = qJD(1) ^ 2;
t382 = qJ(2) * t281 - t383;
t256 = t279 * qJD(1) + qJD(2);
t357 = t256 * t277;
t236 = -qJD(3) * pkin(3) - t357;
t331 = t276 * qJD(3);
t242 = t273 * t339 - t331;
t204 = pkin(4) * t242 - qJ(5) * t244 + t236;
t317 = t277 * t329;
t327 = qJDD(1) * t274;
t240 = qJDD(4) + t317 + t327;
t377 = pkin(7) * t240;
t381 = t258 * t204 - t377;
t372 = pkin(7) * qJD(4);
t322 = t258 * t372;
t380 = -t322 + t384;
t333 = qJD(4) * t273;
t292 = t274 * t331 + t277 * t333;
t207 = t292 * qJD(1) - qJD(4) * t331 - t273 * qJDD(3) - t276 * t326;
t379 = t244 ^ 2;
t378 = pkin(4) * t240;
t373 = g(3) * t277;
t371 = pkin(1) * qJDD(1);
t370 = qJ(5) * t240;
t229 = t249 * qJD(1);
t248 = t274 * t256;
t235 = qJD(3) * pkin(7) + t248;
t206 = t229 * t273 + t235 * t276;
t203 = qJ(5) * t258 + t206;
t369 = t203 * t258;
t368 = t206 * t258;
t367 = t207 * t273;
t366 = t208 * t276;
t312 = pkin(3) * t277 + pkin(7) * t274;
t241 = t312 * qJD(3) + qJD(2);
t365 = t241 * t276;
t364 = t242 * t244;
t363 = t242 * t258;
t362 = t242 * t273;
t361 = t242 * t276;
t360 = t244 * t258;
t359 = t244 * t273;
t358 = t244 * t276;
t355 = t273 * t275;
t354 = t273 * t279;
t352 = t274 * t278;
t351 = t275 * t276;
t350 = t276 * t240;
t349 = t276 * t277;
t348 = t278 * t276;
t305 = pkin(4) * t273 - qJ(5) * t276;
t347 = -qJD(5) * t273 + t258 * t305 - t248;
t247 = t312 * qJD(1);
t346 = t273 * t247 + t256 * t349;
t345 = g(2) * t277 * t348 + g(3) * t353;
t343 = t278 * pkin(1) + t275 * qJ(2);
t272 = t277 ^ 2;
t342 = t274 ^ 2 - t272;
t280 = qJD(3) ^ 2;
t341 = -t280 - t281;
t338 = qJD(3) * t242;
t336 = qJD(3) * t274;
t335 = qJD(3) * t277;
t332 = qJD(4) * t276;
t205 = t229 * t276 - t235 * t273;
t330 = qJD(5) - t205;
t328 = qJDD(1) * qJ(2);
t325 = MDP(19) + MDP(21);
t324 = MDP(20) - MDP(23);
t323 = t277 * t375;
t321 = 0.2e1 * qJD(1) * qJD(2);
t213 = t241 * qJD(1) + t249 * qJDD(1);
t254 = t279 * qJDD(1) + qJDD(2);
t219 = qJDD(3) * pkin(7) + t254 * t274 + t256 * t335;
t320 = -t273 * t213 - t276 * t219 - t229 * t332;
t319 = t277 * t279 * t331 + t273 * t241 + t249 * t332;
t316 = -pkin(4) + t354;
t314 = qJDD(2) - t371;
t313 = t276 * t213 - t273 * t219 - t229 * t333 - t235 * t332;
t230 = t274 * t355 - t348;
t232 = t273 * t352 + t351;
t310 = g(1) * t232 + g(2) * t230;
t231 = t273 * t278 + t274 * t351;
t233 = t274 * t348 - t355;
t309 = -g(1) * t233 - g(2) * t231;
t308 = g(1) * t278 + g(2) * t275;
t306 = pkin(4) * t276 + qJ(5) * t273;
t293 = -t235 * t333 - t320;
t195 = qJD(5) * t258 + t293 + t370;
t196 = qJDD(5) - t313 - t378;
t303 = t195 * t276 + t196 * t273;
t202 = -pkin(4) * t258 + t330;
t302 = t202 * t276 - t203 * t273;
t301 = t202 * t273 + t203 * t276;
t300 = pkin(3) + t306;
t298 = -t279 + t305;
t218 = -qJDD(3) * pkin(3) - t254 * t277 + t256 * t336;
t296 = t240 * t273 + t258 * t332;
t295 = -t258 * t333 + t350;
t294 = 0.2e1 * qJ(2) * t329 + qJDD(3) * t279;
t291 = -t254 + t382;
t290 = t258 * t236 - t377;
t289 = -t325 * t273 - t324 * t276;
t288 = t274 * t383 - t373;
t287 = g(1) * t230 - g(2) * t232 + t273 * t373 + t313;
t286 = -t308 + t321 + 0.2e1 * t328;
t285 = -t279 * t280 + t286;
t284 = t204 * t244 + qJDD(5) - t287;
t283 = -g(1) * t231 + g(2) * t233 - g(3) * t349 + t293;
t282 = (t358 + t362) * MDP(22) + t302 * MDP(24) + (t324 * t273 - t325 * t276) * t258;
t268 = t278 * qJ(2);
t265 = qJDD(3) * t277;
t222 = t298 * t277;
t221 = -t249 * t276 + t316 * t274;
t220 = qJ(5) * t274 + t344;
t216 = pkin(4) * t244 + qJ(5) * t242;
t212 = -t247 * t276 + (-pkin(4) * qJD(1) + t256 * t273) * t277;
t211 = qJ(5) * t339 + t346;
t201 = (t306 * qJD(4) - qJD(5) * t276) * t277 - t298 * t336;
t200 = -t207 + t363;
t199 = t344 * qJD(4) + t316 * t335 - t365;
t198 = qJ(5) * t335 + (-t279 * t333 + qJD(5)) * t274 + t319;
t197 = pkin(4) * t208 + qJ(5) * t207 - qJD(5) * t244 + t218;
t1 = [qJDD(1) * MDP(1) - t383 * MDP(2) + t308 * MDP(3) + (qJDD(2) + t383 - 0.2e1 * t371) * MDP(4) + t286 * MDP(5) + (-t314 * pkin(1) - g(1) * (-pkin(1) * t275 + t268) - g(2) * t343 + (t321 + t328) * qJ(2)) * MDP(6) + (qJDD(1) * t272 - 0.2e1 * t274 * t317) * MDP(7) + 0.2e1 * (-t274 * t326 + t342 * t329) * MDP(8) + (-t274 * t280 + t265) * MDP(9) + (-qJDD(3) * t274 - t277 * t280) * MDP(10) + (t285 * t274 + t294 * t277) * MDP(12) + (-t294 * t274 + t285 * t277) * MDP(13) + (-t207 * t349 - t292 * t244) * MDP(14) + ((t359 + t361) * t336 + (t367 - t366 + (-t358 + t362) * qJD(4)) * t277) * MDP(15) + ((-t258 * t331 - t207) * t274 + (qJD(3) * t244 + t295) * t277) * MDP(16) + ((t258 * t337 - t208) * t274 + (-t296 - t338) * t277) * MDP(17) + (t240 * t274 + t258 * t335) * MDP(18) + ((-t249 * t333 + t365) * t258 + t249 * t350 + (-t236 * t337 + (-t296 + t338) * t279 + t313) * t274 + (t236 * t332 - t279 * t208 + t218 * t273 + (-t258 * t354 + t205) * qJD(3)) * t277 + t309) * MDP(19) + (-t319 * t258 - t344 * t240 + ((t258 * t279 + t235) * t333 + (-t236 * t276 + t244 * t279) * qJD(3) + t320) * t274 + (-qJD(3) * t206 + t207 * t279 + t218 * t276 - t236 * t333) * t277 + t310) * MDP(20) + (-t199 * t258 + t201 * t242 + t208 * t222 - t221 * t240 + (-t204 * t337 - t196) * t274 + (-qJD(3) * t202 + t197 * t273 + t204 * t332) * t277 + t309) * MDP(21) + (-t198 * t242 + t199 * t244 - t207 * t221 - t208 * t220 - t302 * t336 + (-t301 * qJD(4) - t195 * t273 + t196 * t276 + t308) * t277) * MDP(22) + (t198 * t258 - t201 * t244 + t207 * t222 + t220 * t240 + (t204 * t331 + t195) * t274 + (qJD(3) * t203 - t197 * t276 + t204 * t333) * t277 - t310) * MDP(23) + (t195 * t220 + t203 * t198 + t197 * t222 + t204 * t201 + t196 * t221 + t202 * t199 - g(1) * (pkin(3) * t352 + pkin(4) * t233 + qJ(5) * t232 - t278 * t376 + t268) - g(2) * (pkin(4) * t231 + pkin(6) * t278 + qJ(5) * t230 + t343) + (-g(1) * t279 - g(2) * t311) * t275) * MDP(24); qJDD(1) * MDP(4) - t281 * MDP(5) + (t314 - t382) * MDP(6) + t265 * MDP(12) + t383 * MDP(24) + t325 * t242 * t336 + t282 * qJD(1) + (t341 * MDP(13) - t197 * MDP(24) - t325 * t208 + t324 * t207 + ((t359 - t361) * MDP(22) + t301 * MDP(24) + t289 * t258) * qJD(3)) * t277 + (t341 * MDP(12) - qJDD(3) * MDP(13) + (-t366 - t367) * MDP(22) + t303 * MDP(24) + (t204 * MDP(24) + t324 * t244) * qJD(3) + t289 * t240 + t282 * qJD(4)) * t274; MDP(9) * t326 - MDP(10) * t327 + qJDD(3) * MDP(11) + (-t291 * t277 + t374) * MDP(12) + (t291 * t274 + t373) * MDP(13) + (t258 * t358 - t367) * MDP(14) + ((-t207 - t363) * t276 + (-t208 - t360) * t273) * MDP(15) + ((-t244 * t277 + t258 * t353) * qJD(1) + t296) * MDP(16) + ((t242 * t277 - t258 * t356) * qJD(1) + t295) * MDP(17) - t258 * MDP(18) * t339 + (-t205 * t339 - t242 * t248 - pkin(3) * t208 + (-t323 - t218 + (-t247 - t372) * t258) * t276 + (t258 * t357 + t290) * t273 + t345) * MDP(19) + (pkin(3) * t207 + t346 * t258 + t206 * t339 - t244 * t248 + t290 * t276 + (t218 - t380) * t273) * MDP(20) + (t202 * t339 - t208 * t300 + t212 * t258 + t347 * t242 + (-t197 - t322 - t323) * t276 + t381 * t273 + t345) * MDP(21) + (t211 * t242 - t212 * t244 + (t195 + t258 * t202 + (-t208 + t334) * pkin(7)) * t276 + (t196 - t369 + (qJD(4) * t242 - t207) * pkin(7)) * t273 + t288) * MDP(22) + (-t203 * t339 - t207 * t300 - t211 * t258 - t347 * t244 - t381 * t276 + (-t197 + t380) * t273) * MDP(23) + (-t202 * t212 - t203 * t211 + t347 * t204 + (t302 * qJD(4) + t288 + t303) * pkin(7) + (-t197 + t384) * t300) * MDP(24) + (t277 * t274 * MDP(7) - t342 * MDP(8)) * t281; MDP(14) * t364 + (-t242 ^ 2 + t379) * MDP(15) + t200 * MDP(16) + (t360 - t208) * MDP(17) + t240 * MDP(18) + (-t236 * t244 + t287 + t368) * MDP(19) + (t205 * t258 + t236 * t242 - t283) * MDP(20) + (-t216 * t242 - t284 + t368 + 0.2e1 * t378) * MDP(21) + (pkin(4) * t207 - qJ(5) * t208 + (t203 - t206) * t244 + (t202 - t330) * t242) * MDP(22) + (0.2e1 * t370 - t204 * t242 + t216 * t244 + (0.2e1 * qJD(5) - t205) * t258 + t283) * MDP(23) + (t195 * qJ(5) - t196 * pkin(4) - t204 * t216 - t202 * t206 - g(1) * (-pkin(4) * t230 + qJ(5) * t231) - g(2) * (pkin(4) * t232 - qJ(5) * t233) + t305 * t373 + t330 * t203) * MDP(24); (-t240 + t364) * MDP(21) + t200 * MDP(22) + (-t258 ^ 2 - t379) * MDP(23) + (t284 - t369 - t378) * MDP(24);];
tau = t1;
