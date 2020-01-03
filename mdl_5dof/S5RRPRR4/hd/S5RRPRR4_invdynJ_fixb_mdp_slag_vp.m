% Calculate vector of inverse dynamics joint torques for
% S5RRPRR4
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
%   pkin=[a2,a3,a4,a5,d1,d2,d4,d5,theta3]';
% MDP [21x1]
%   Minimal dynamic parameter vector (fixed base model)
%   see S5RRPRR4_convert_par2_MPV_fixb.m
% 
% Output:
% tau [5x1]
%   joint torques of inverse dynamics (contains inertial, gravitational Coriolis and centrifugal forces)

% Quelle: HybrDyn-Toolbox
% Datum: 2020-01-03 12:02
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tau = S5RRPRR4_invdynJ_fixb_mdp_slag_vp(qJ, qJD, qJDD, g, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(5,1),zeros(3,1),zeros(9,1),zeros(21,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRPRR4_invdynJ_fixb_mdp_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RRPRR4_invdynJ_fixb_mdp_slag_vp: qJD has to be [5x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [5 1]), ...
  'S5RRPRR4_invdynJ_fixb_mdp_slag_vp: qJDD has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RRPRR4_invdynJ_fixb_mdp_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5RRPRR4_invdynJ_fixb_mdp_slag_vp: pkin has to be [9x1] (double)');
assert(isreal(MDP) && all(size(MDP) == [21 1]), ...
  'S5RRPRR4_invdynJ_fixb_mdp_slag_vp: MDP has to be [21x1] (double)'); 

%% Symbolic Calculation
% From invdyn_joint_fixb_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2020-01-03 12:02:11
% EndTime: 2020-01-03 12:02:14
% DurationCPUTime: 1.49s
% Computational Cost: add. (1409->234), mult. (2244->317), div. (0->0), fcn. (1448->16), ass. (0->146)
t321 = sin(qJ(5));
t322 = sin(qJ(4));
t325 = cos(qJ(5));
t326 = cos(qJ(4));
t261 = t321 * t326 + t322 * t325;
t314 = qJD(1) + qJD(2);
t247 = t261 * t314;
t327 = cos(qJ(2));
t404 = pkin(1) * t327;
t300 = qJDD(1) * t404;
t312 = qJDD(1) + qJDD(2);
t323 = sin(qJ(2));
t394 = pkin(1) * qJD(2);
t369 = qJD(1) * t394;
t248 = pkin(2) * t312 - t323 * t369 + t300;
t319 = sin(pkin(9));
t320 = cos(pkin(9));
t373 = qJDD(1) * t323;
t406 = pkin(1) * t373 + t327 * t369;
t227 = t248 * t320 - t406 * t319;
t223 = -pkin(3) * t312 - t227;
t318 = qJ(1) + qJ(2);
t303 = pkin(9) + t318;
t291 = sin(t303);
t292 = cos(t303);
t351 = -g(2) * t292 - g(3) * t291;
t341 = -t223 + t351;
t381 = t325 * t326;
t385 = t321 * t322;
t260 = -t381 + t385;
t405 = g(2) * t291 - g(3) * t292;
t395 = pkin(1) * qJD(1);
t371 = t327 * t395;
t268 = pkin(2) * t314 + t371;
t372 = t323 * t395;
t288 = t320 * t372;
t240 = t319 * t268 + t288;
t362 = t240 + (pkin(7) + pkin(8)) * t314;
t225 = t326 * qJD(3) - t362 * t322;
t313 = qJD(4) + qJD(5);
t226 = qJD(3) * t322 + t326 * t362;
t403 = pkin(2) * t320;
t402 = pkin(4) * t326;
t299 = pkin(2) + t404;
t386 = t320 * t323;
t379 = pkin(1) * t386 + t319 * t299;
t250 = pkin(7) + t379;
t397 = -pkin(8) - t250;
t293 = pkin(2) * t319 + pkin(7);
t396 = -pkin(8) - t293;
t393 = t226 * t325;
t287 = t319 * t372;
t253 = t320 * t371 - t287;
t392 = t253 * t313;
t317 = qJ(4) + qJ(5);
t305 = sin(t317);
t391 = t291 * t305;
t390 = t292 * t305;
t389 = t312 * t326;
t388 = t314 * t322;
t387 = t319 * t323;
t382 = t322 * t326;
t380 = qJDD(3) - g(1);
t315 = t322 ^ 2;
t378 = -t326 ^ 2 + t315;
t376 = qJD(4) * t322;
t375 = qJD(4) * t326;
t374 = qJD(5) * t321;
t370 = pkin(4) * t376;
t367 = t314 * t385;
t366 = t314 * t381;
t228 = t319 * t248 + t406 * t320;
t365 = -pkin(3) - t402;
t364 = t314 * t375;
t224 = pkin(7) * t312 + t228;
t363 = pkin(8) * t312 + t224;
t306 = sin(t318);
t308 = cos(t318);
t361 = g(2) * t306 - g(3) * t308;
t360 = qJD(4) * t397;
t359 = qJD(4) * t396;
t239 = t268 * t320 - t287;
t358 = -pkin(1) * t387 + t299 * t320;
t357 = qJD(1) * (-qJD(2) + t314);
t356 = qJD(2) * (-qJD(1) - t314);
t211 = (t314 * t376 - t389) * pkin(4) + t223;
t229 = t314 * t365 - t239;
t233 = t313 * t260;
t354 = g(2) * t390 + g(3) * t391 + t211 * t261 - t229 * t233;
t235 = -pkin(3) * t314 - t239;
t353 = t235 * t375 - t341 * t322;
t249 = -pkin(3) - t358;
t251 = t319 * t371 + t288;
t352 = -t251 + t370;
t350 = -g(2) * t308 - g(3) * t306;
t349 = t260 * t312;
t222 = qJD(4) * pkin(4) + t225;
t347 = -t222 * t321 - t393;
t311 = qJDD(4) + qJDD(5);
t346 = -t233 * t313 + t261 * t311;
t237 = t397 * t322;
t309 = t326 * pkin(8);
t238 = t250 * t326 + t309;
t345 = t237 * t325 - t238 * t321;
t344 = t237 * t321 + t238 * t325;
t258 = t396 * t322;
t259 = t293 * t326 + t309;
t343 = t258 * t325 - t259 * t321;
t342 = t258 * t321 + t259 * t325;
t340 = t300 + t350;
t252 = (t319 * t327 + t386) * t394;
t329 = qJD(4) ^ 2;
t339 = t249 * t312 + t250 * t329 + t252 * t314;
t294 = -pkin(3) - t403;
t338 = -t251 * t314 + t293 * t329 + t294 * t312;
t212 = qJD(5) * t366 + t261 * t312 - t313 * t367 + t325 * t364;
t245 = -t366 + t367;
t337 = t247 * t245 * MDP(15) + (t245 * t313 + t212) * MDP(17) - t349 * MDP(18) + (-t245 ^ 2 + t247 ^ 2) * MDP(16) + t311 * MDP(19);
t336 = -t235 * t314 - t224 + t405;
t254 = (t320 * t327 - t387) * t394;
t335 = -qJDD(4) * t250 + (t249 * t314 - t254) * qJD(4);
t334 = -qJDD(4) * t293 + (t294 * t314 + t253) * qJD(4);
t234 = t313 * t261;
t307 = cos(t317);
t333 = t211 * t260 + t229 * t234 + t307 * t351;
t213 = t234 * t314 + t349;
t217 = -t234 * t313 - t260 * t311;
t272 = qJDD(4) * t322 + t326 * t329;
t273 = qJDD(4) * t326 - t322 * t329;
t332 = (-t212 * t260 - t213 * t261 + t233 * t245 - t234 * t247) * MDP(16) + (t212 * t261 - t233 * t247) * MDP(15) + t346 * MDP(17) + t217 * MDP(18) + 0.2e1 * (-qJD(4) * t314 * t378 + t312 * t382) * MDP(9) + (t312 * t315 + 0.2e1 * t322 * t364) * MDP(8) + t272 * MDP(10) + t273 * MDP(11) + t312 * MDP(4);
t302 = t326 * qJDD(3);
t203 = qJDD(4) * pkin(4) - t226 * qJD(4) - t363 * t322 + t302;
t331 = t229 * t245 + t226 * t374 + g(1) * t305 + (-t226 * t313 - t203) * t321 + t405 * t307;
t204 = t225 * qJD(4) + t322 * qJDD(3) + t363 * t326;
t330 = -g(1) * t307 + g(2) * t391 - g(3) * t390 + qJD(5) * t347 + t325 * t203 - t321 * t204 - t229 * t247;
t328 = cos(qJ(1));
t324 = sin(qJ(1));
t271 = t365 - t403;
t256 = t326 * t359;
t255 = t322 * t359;
t243 = t249 - t402;
t241 = t252 + t370;
t230 = t235 * t376;
t221 = -t254 * t322 + t326 * t360;
t220 = t254 * t326 + t322 * t360;
t1 = [(t241 * t247 + t243 * t212 - (qJD(5) * t345 + t220 * t325 + t221 * t321) * t313 - t344 * t311 + t354) * MDP(21) + qJDD(1) * MDP(1) + ((t312 * t327 + t323 * t356) * pkin(1) + t340) * MDP(5) + t332 + (t230 + t335 * t322 + (-t339 + t341) * t326) * MDP(13) + (t322 * t339 + t326 * t335 + t353) * MDP(14) + (((-qJDD(1) - t312) * t323 + t327 * t356) * pkin(1) + t361) * MDP(6) + (-g(2) * t328 - g(3) * t324) * MDP(2) + (g(2) * t324 - g(3) * t328) * MDP(3) + (t228 * t379 + t240 * t254 + t227 * t358 - t239 * t252 - g(2) * (pkin(1) * t328 + pkin(2) * t308) - g(3) * (pkin(1) * t324 + pkin(2) * t306)) * MDP(7) + (t241 * t245 + t243 * t213 + (-qJD(5) * t344 - t220 * t321 + t221 * t325) * t313 + t345 * t311 + t333) * MDP(20); (t323 * pkin(1) * t357 + t340) * MDP(5) + (t322 * t338 + t326 * t334 + t353) * MDP(14) + (t239 * t251 - t240 * t253 + (t227 * t320 + t228 * t319 + t350) * pkin(2)) * MDP(7) + t332 + ((t327 * t357 - t373) * pkin(1) + t361) * MDP(6) + (t230 + t334 * t322 + (-t338 + t341) * t326) * MDP(13) + (t271 * t213 + (-qJD(5) * t342 - t255 * t321 + t256 * t325) * t313 + t343 * t311 + t261 * t392 + t352 * t245 + t333) * MDP(20) + (t271 * t212 - (qJD(5) * t343 + t255 * t325 + t256 * t321) * t313 - t342 * t311 - t260 * t392 + t352 * t247 + t354) * MDP(21); t273 * MDP(13) - t272 * MDP(14) + t217 * MDP(20) - MDP(21) * t346 + MDP(7) * t380; t322 * t312 * MDP(10) + MDP(11) * t389 + qJDD(4) * MDP(12) + (-g(1) * t326 + t322 * t336 + t302) * MDP(13) + (-t322 * t380 + t326 * t336) * MDP(14) + (-(-t225 * t321 - t393) * t313 + (-t245 * t388 + t325 * t311 - t313 * t374) * pkin(4) + t330) * MDP(20) + ((-qJD(5) * t222 + t225 * t313 - t204) * t325 + (-qJD(5) * t325 * t313 - t247 * t388 - t321 * t311) * pkin(4) + t331) * MDP(21) + t337 + (-MDP(8) * t382 + t378 * MDP(9)) * t314 ^ 2; (-t313 * t347 + t330) * MDP(20) + ((-t204 + (-qJD(5) + t313) * t222) * t325 + t331) * MDP(21) + t337;];
tau = t1;
