% Calculate vector of inverse dynamics joint torques for
% S5RPRPP3
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
%   pkin=[a2,a3,a4,a5,d1,d3,theta2]';
% MDP [22x1]
%   Minimal dynamic parameter vector (fixed base model)
%   see S5RPRPP3_convert_par2_MPV_fixb.m
% 
% Output:
% tau [5x1]
%   joint torques of inverse dynamics (contains inertial, gravitational Coriolis and centrifugal forces)

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 18:13
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tau = S5RPRPP3_invdynJ_fixb_mdp_slag_vp(qJ, qJD, qJDD, g, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(5,1),zeros(3,1),zeros(7,1),zeros(22,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPRPP3_invdynJ_fixb_mdp_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RPRPP3_invdynJ_fixb_mdp_slag_vp: qJD has to be [5x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [5 1]), ...
  'S5RPRPP3_invdynJ_fixb_mdp_slag_vp: qJDD has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RPRPP3_invdynJ_fixb_mdp_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [7 1]), ...
  'S5RPRPP3_invdynJ_fixb_mdp_slag_vp: pkin has to be [7x1] (double)');
assert(isreal(MDP) && all(size(MDP) == [22 1]), ...
  'S5RPRPP3_invdynJ_fixb_mdp_slag_vp: MDP has to be [22x1] (double)'); 

%% Symbolic Calculation
% From invdyn_joint_fixb_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 18:12:59
% EndTime: 2019-12-31 18:13:03
% DurationCPUTime: 2.74s
% Computational Cost: add. (1527->303), mult. (3407->331), div. (0->0), fcn. (2369->8), ass. (0->137)
t306 = sin(qJ(1));
t307 = cos(qJ(1));
t343 = g(1) * t306 - g(2) * t307;
t373 = qJDD(1) * pkin(1);
t325 = -qJDD(2) + t373 + t343;
t301 = sin(pkin(7));
t376 = pkin(6) + qJ(2);
t267 = t376 * t301;
t263 = qJD(1) * t267;
t302 = cos(pkin(7));
t268 = t376 * t302;
t264 = qJD(1) * t268;
t305 = sin(qJ(3));
t386 = cos(qJ(3));
t362 = -t386 * t263 - t305 * t264;
t354 = -qJD(4) + t362;
t378 = pkin(3) + qJ(5);
t341 = t378 * qJDD(3);
t262 = t301 * t386 + t305 * t302;
t259 = t262 * qJD(3);
t342 = qJDD(1) * t386;
t350 = qJDD(1) * t305;
t330 = t301 * t350 - t302 * t342;
t229 = qJD(1) * t259 + t330;
t345 = qJD(1) * t386;
t334 = t302 * t345;
t358 = qJD(1) * t305;
t347 = t301 * t358;
t254 = -t334 + t347;
t395 = t229 * qJ(5) + t254 * qJD(5);
t298 = pkin(7) + qJ(3);
t291 = sin(t298);
t292 = cos(t298);
t360 = t292 * pkin(3) + t291 * qJ(4);
t394 = t302 * MDP(4) - t301 * MDP(5);
t381 = g(1) * t307;
t331 = g(2) * t306 + t381;
t256 = t262 * qJD(1);
t392 = qJ(2) * qJDD(1);
t391 = MDP(16) - MDP(21);
t390 = -pkin(4) * t229 + qJDD(5);
t351 = qJD(1) * qJD(2);
t388 = qJDD(1) * t376 + t351;
t238 = t388 * t301;
t239 = t388 * t302;
t344 = qJD(3) * t386;
t355 = qJD(3) * t305;
t338 = -t305 * t238 + t386 * t239 - t263 * t344 - t264 * t355;
t379 = g(3) * t291;
t389 = -t292 * t331 + t338 - t379;
t348 = t386 * t302;
t219 = (qJD(2) * t301 + qJD(3) * t268) * t305 - qJD(2) * t348 + t267 * t344;
t387 = t254 ^ 2;
t249 = t256 ^ 2;
t346 = t301 * t355;
t349 = qJD(3) * t334 + t301 * t342 + t302 * t350;
t228 = qJD(1) * t346 - t349;
t385 = pkin(4) * t228;
t383 = pkin(4) * t254;
t377 = pkin(4) + t376;
t375 = qJ(4) * t229;
t374 = qJ(4) * t254;
t372 = qJDD(3) * pkin(3);
t289 = pkin(2) * t302 + pkin(1);
t266 = -qJD(1) * t289 + qJD(2);
t319 = -qJ(4) * t256 + t266;
t205 = t254 * t378 + t319;
t371 = t205 * t254;
t370 = t254 * t256;
t369 = t291 * t306;
t368 = t291 * t307;
t367 = t292 * qJ(5);
t366 = t292 * t306;
t365 = t292 * t307;
t231 = -t305 * t263 + t386 * t264;
t359 = t301 ^ 2 + t302 ^ 2;
t357 = qJD(3) * t362;
t356 = qJD(3) * t231;
t216 = -pkin(4) * t256 + t362;
t353 = qJD(4) - t216;
t217 = t231 - t383;
t352 = -qJD(5) - t217;
t340 = t359 * qJD(1) ^ 2;
t337 = t386 * t238 + t305 * t239 - t263 * t355 + t264 * t344;
t336 = 0.2e1 * t359;
t335 = g(2) * (pkin(3) * t365 + qJ(4) * t368 + t307 * t289);
t333 = -g(1) * t369 + g(2) * t368;
t332 = g(1) * t366 - g(2) * t365;
t224 = -qJD(3) * qJ(4) - t231;
t329 = -qJDD(4) - t337;
t258 = -t302 * t344 + t346;
t328 = qJ(4) * t258 - qJD(4) * t262;
t326 = -qJ(4) * t262 - t289;
t324 = -t289 - t360;
t299 = qJDD(3) * qJ(4);
t300 = qJD(3) * qJD(4);
t202 = -t299 - t300 - t338;
t232 = t267 * t386 + t305 * t268;
t233 = -t305 * t267 + t268 * t386;
t265 = -qJDD(1) * t289 + qJDD(2);
t322 = -g(1) * t368 - g(2) * t369 + g(3) * t292 + t337;
t320 = qJDD(4) + t322;
t200 = -t202 + t390;
t318 = qJD(3) * t219 - qJDD(3) * t233 + t333;
t220 = qJD(2) * t262 + qJD(3) * t233;
t317 = -qJD(3) * t220 - qJDD(3) * t232 + t332;
t218 = pkin(3) * t254 + t319;
t316 = t218 * t256 + t320;
t315 = t336 * t351 - t331;
t314 = pkin(3) * t229 + qJ(4) * t228 + t265;
t313 = t205 * t256 + t320 - t385;
t201 = -qJD(4) * t256 + t314;
t312 = t314 - t343;
t311 = 0.2e1 * t299 + 0.2e1 * t300 + t389;
t308 = qJD(3) ^ 2;
t271 = qJ(4) * t365;
t269 = qJ(4) * t366;
t261 = t301 * t305 - t348;
t242 = qJD(3) * t254;
t227 = pkin(3) * t261 + t326;
t226 = pkin(3) * t256 + t374;
t223 = -qJD(3) * pkin(3) - t354;
t222 = -t261 * pkin(4) + t233;
t221 = t262 * pkin(4) + t232;
t215 = t261 * t378 + t326;
t214 = pkin(3) * t259 + t328;
t212 = (t254 - t347) * qJD(3) + t349;
t210 = t256 * t378 + t374;
t209 = qJD(5) - t224 - t383;
t208 = -qJD(3) * t378 + t353;
t207 = -t258 * pkin(4) + t220;
t206 = -pkin(4) * t259 - t219;
t204 = qJD(5) * t261 + t259 * t378 + t328;
t203 = -t329 - t372;
t199 = -qJD(3) * qJD(5) - t329 - t341 - t385;
t198 = t201 + t395;
t1 = [qJDD(1) * MDP(1) + t343 * MDP(2) + t331 * MDP(3) + (t336 * t392 + t315) * MDP(6) + (t325 * pkin(1) + (t359 * t392 + t315) * qJ(2)) * MDP(7) + (-t228 * t262 - t256 * t258) * MDP(8) + (t228 * t261 - t229 * t262 + t254 * t258 - t256 * t259) * MDP(9) + (-qJD(3) * t258 + qJDD(3) * t262) * MDP(10) + (-qJD(3) * t259 - qJDD(3) * t261) * MDP(11) + (-t229 * t289 + t259 * t266 + t261 * t265 + t317) * MDP(13) + (t228 * t289 - t258 * t266 + t262 * t265 + t318) * MDP(14) + (t202 * t261 + t203 * t262 + t219 * t254 + t220 * t256 - t223 * t258 + t224 * t259 - t228 * t232 - t229 * t233 - t331) * MDP(15) + (-t201 * t261 - t214 * t254 - t218 * t259 - t227 * t229 - t317) * MDP(16) + (-t201 * t262 - t214 * t256 + t218 * t258 + t227 * t228 - t318) * MDP(17) + (t201 * t227 + t218 * t214 - t202 * t233 + t224 * t219 + t203 * t232 + t223 * t220 - t376 * t381 - t335 + (-g(1) * t324 - g(2) * t376) * t306) * MDP(18) + (t199 * t262 - t200 * t261 - t206 * t254 + t207 * t256 - t208 * t258 - t209 * t259 - t221 * t228 - t222 * t229 - t331) * MDP(19) + (qJD(3) * t206 + qJDD(3) * t222 - t198 * t262 - t204 * t256 + t205 * t258 + t215 * t228 - t333) * MDP(20) + (-qJD(3) * t207 - qJDD(3) * t221 + t198 * t261 + t204 * t254 + t205 * t259 + t215 * t229 + t332) * MDP(21) + (t198 * t215 + t205 * t204 + t199 * t221 + t208 * t207 + t200 * t222 + t209 * t206 - t335 + (-g(1) * t377 - g(2) * t367) * t307 + (-g(1) * (t324 - t367) - g(2) * t377) * t306) * MDP(22) + t394 * (t325 + t373); -MDP(6) * t340 + (-qJ(2) * t340 - t325) * MDP(7) + (t228 + t242) * MDP(17) + (-t224 * t254 + (-qJD(4) - t223) * t256 + t312) * MDP(18) + (t209 * t254 + (-qJD(4) - t208) * t256 + t312 + t395) * MDP(22) + (-MDP(14) + MDP(20)) * ((t254 + t347) * qJD(3) - t349) - t394 * qJDD(1) + (MDP(13) - t391) * (0.2e1 * qJD(3) * t256 + t330) + (MDP(15) + MDP(19)) * (-t387 - t249); MDP(8) * t370 + (-t387 + t249) * MDP(9) + t212 * MDP(10) - t330 * MDP(11) + qJDD(3) * MDP(12) + (-t256 * t266 - t322 + t356) * MDP(13) + (t254 * t266 + t357 - t389) * MDP(14) + (pkin(3) * t228 - t375 + (-t224 - t231) * t256 + (t223 + t354) * t254) * MDP(15) + (t226 * t254 + t316 - t356 - 0.2e1 * t372) * MDP(16) + (-t218 * t254 + t226 * t256 + t311 - t357) * MDP(17) + (-t202 * qJ(4) - t203 * pkin(3) - t218 * t226 - t223 * t231 - g(1) * (-pkin(3) * t368 + t271) - g(2) * (-pkin(3) * t369 + t269) - g(3) * t360 + t354 * t224) * MDP(18) + (-t375 + t228 * t378 + (t209 + t352) * t256 + (t208 - t353) * t254) * MDP(19) + (-qJD(3) * t216 + t210 * t256 + t311 - t371 + t390) * MDP(20) + (-t210 * t254 + (0.2e1 * qJD(5) + t217) * qJD(3) + 0.2e1 * t341 - t313) * MDP(21) + (t200 * qJ(4) - t205 * t210 - g(1) * t271 - g(2) * t269 - g(3) * (t360 + t367) + t353 * t209 + t352 * t208 + (t331 * t291 - t199) * t378) * MDP(22); (-t228 + t242) * MDP(15) + (qJD(3) * t224 + t316 - t372) * MDP(18) + t212 * MDP(19) + (-t341 + (-qJD(5) - t209) * qJD(3) + t313) * MDP(22) + (MDP(17) + MDP(20)) * (-t308 - t249) + t391 * (qJDD(3) - t370); -t330 * MDP(19) + (qJDD(3) + t370) * MDP(20) + (-t387 - t308) * MDP(21) + (-g(1) * t365 - g(2) * t366 + t200 - t371 - t379) * MDP(22) + ((-t301 * t345 - t302 * t358 + t256) * MDP(19) + t208 * MDP(22)) * qJD(3);];
tau = t1;
