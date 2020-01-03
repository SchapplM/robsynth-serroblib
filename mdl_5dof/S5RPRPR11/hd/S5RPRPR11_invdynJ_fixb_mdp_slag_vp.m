% Calculate vector of inverse dynamics joint torques for
% S5RPRPR11
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
%   pkin=[a2,a3,a4,a5,d1,d3,d5,theta2]';
% MDP [25x1]
%   Minimal dynamic parameter vector (fixed base model)
%   see S5RPRPR11_convert_par2_MPV_fixb.m
% 
% Output:
% tau [5x1]
%   joint torques of inverse dynamics (contains inertial, gravitational Coriolis and centrifugal forces)

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 18:28
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tau = S5RPRPR11_invdynJ_fixb_mdp_slag_vp(qJ, qJD, qJDD, g, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(5,1),zeros(3,1),zeros(8,1),zeros(25,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPRPR11_invdynJ_fixb_mdp_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RPRPR11_invdynJ_fixb_mdp_slag_vp: qJD has to be [5x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [5 1]), ...
  'S5RPRPR11_invdynJ_fixb_mdp_slag_vp: qJDD has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RPRPR11_invdynJ_fixb_mdp_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RPRPR11_invdynJ_fixb_mdp_slag_vp: pkin has to be [8x1] (double)');
assert(isreal(MDP) && all(size(MDP) == [25 1]), ...
  'S5RPRPR11_invdynJ_fixb_mdp_slag_vp: MDP has to be [25x1] (double)'); 

%% Symbolic Calculation
% From invdyn_joint_fixb_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 18:28:01
% EndTime: 2019-12-31 18:28:05
% DurationCPUTime: 2.94s
% Computational Cost: add. (1625->305), mult. (3718->376), div. (0->0), fcn. (2791->10), ass. (0->144)
t329 = qJDD(3) - qJDD(5);
t338 = cos(pkin(8));
t426 = cos(qJ(3));
t392 = t426 * t338;
t380 = qJD(1) * t392;
t337 = sin(pkin(8));
t341 = sin(qJ(3));
t408 = t337 * t341;
t391 = qJD(1) * t408;
t295 = -t380 + t391;
t305 = t426 * t337 + t341 * t338;
t297 = t305 * qJD(1);
t340 = sin(qJ(5));
t343 = cos(qJ(5));
t367 = -t343 * t295 + t297 * t340;
t401 = qJD(3) * t341;
t390 = t337 * t401;
t387 = qJDD(1) * t426;
t395 = qJDD(1) * t341;
t393 = qJD(3) * t380 + t337 * t387 + t338 * t395;
t261 = qJD(1) * t390 - t393;
t300 = t305 * qJD(3);
t374 = t337 * t395 - t338 * t387;
t262 = qJD(1) * t300 + t374;
t382 = t367 * qJD(5) + t343 * t261 - t340 * t262;
t333 = qJD(3) - qJD(5);
t415 = t367 * t333;
t435 = t295 * t340 + t343 * t297;
t443 = MDP(19) * t367 * t435 + (-t367 ^ 2 + t435 ^ 2) * MDP(20) - (t382 + t415) * MDP(21) - t329 * MDP(23);
t416 = t435 * t333;
t324 = pkin(2) * t338 + pkin(1);
t308 = -t324 * qJDD(1) + qJDD(2);
t442 = qJ(4) * t261 + t308;
t309 = -t324 * qJD(1) + qJD(2);
t441 = -qJ(4) * t297 + t309;
t418 = qJDD(1) * pkin(1);
t344 = cos(qJ(1));
t342 = sin(qJ(1));
t425 = g(1) * t342;
t436 = -g(2) * t344 + t425;
t358 = -qJDD(2) + t418 + t436;
t422 = pkin(6) + qJ(2);
t312 = t422 * t338;
t307 = qJD(1) * t312;
t292 = t341 * t307;
t311 = t422 * t337;
t306 = qJD(1) * t311;
t265 = -t426 * t306 - t292;
t397 = qJD(4) - t265;
t345 = -pkin(3) - pkin(4);
t396 = qJD(1) * qJD(2);
t429 = t422 * qJDD(1) + t396;
t279 = t429 * t337;
t280 = t429 * t338;
t389 = qJD(3) * t426;
t383 = t426 * t279 + t341 * t280 - t306 * t401 + t307 * t389;
t363 = qJDD(4) + t383;
t222 = pkin(7) * t261 + t345 * qJDD(3) + t363;
t334 = qJDD(3) * qJ(4);
t335 = qJD(3) * qJD(4);
t394 = -t341 * t279 + t426 * t280 - t306 * t389;
t227 = -t307 * t401 + t334 + t335 + t394;
t223 = pkin(7) * t262 + t227;
t232 = t345 * t295 - t441;
t332 = pkin(8) + qJ(3);
t326 = sin(t332);
t327 = cos(t332);
t288 = t326 * t343 - t327 * t340;
t273 = t288 * t342;
t410 = t327 * t344;
t411 = t326 * t344;
t275 = t340 * t410 - t343 * t411;
t364 = t326 * t340 + t327 * t343;
t438 = -g(1) * t275 + g(2) * t273 - g(3) * t364 - t343 * t222 + t340 * t223 + t232 * t435;
t268 = -t341 * t311 + t426 * t312;
t434 = t338 * MDP(4) - t337 * MDP(5);
t423 = g(2) * t342;
t378 = g(1) * t344 + t423;
t398 = -pkin(7) * t297 + t397;
t433 = qJ(2) * qJDD(1);
t274 = t364 * t342;
t276 = t364 * t344;
t432 = -g(1) * t276 - g(2) * t274 - g(3) * t288 + t340 * t222 + t343 * t223 - t232 * t367;
t245 = -t311 * t389 + qJD(2) * t392 + (-qJD(2) * t337 - qJD(3) * t312) * t341;
t431 = -qJD(3) * t245 - qJDD(3) * t268 - t326 * t436;
t385 = t261 * t340 + t343 * t262;
t225 = qJD(5) * t435 - t385;
t428 = g(3) * t326 + (t265 + t292) * qJD(3) + t378 * t327 - t394;
t427 = t297 ^ 2;
t420 = qJ(4) * t295;
t417 = qJDD(3) * pkin(3);
t414 = t295 * t297;
t266 = -t341 * t306 + t426 * t307;
t404 = t337 ^ 2 + t338 ^ 2;
t402 = qJD(3) * t266;
t400 = qJD(4) * t297;
t386 = t404 * qJD(1) ^ 2;
t381 = 0.2e1 * t404;
t243 = pkin(7) * t295 + t266;
t376 = pkin(3) * t327 + qJ(4) * t326;
t373 = qJ(4) * t343 + t340 * t345;
t372 = -qJ(4) * t340 + t343 * t345;
t235 = t345 * qJD(3) + t398;
t336 = qJD(3) * qJ(4);
t236 = t243 + t336;
t371 = t343 * t235 - t340 * t236;
t370 = -t340 * t235 - t343 * t236;
t267 = t426 * t311 + t341 * t312;
t247 = -t305 * pkin(7) + t267;
t304 = -t392 + t408;
t248 = pkin(7) * t304 + t268;
t369 = t247 * t343 - t248 * t340;
t368 = t247 * t340 + t248 * t343;
t365 = t343 * t304 - t305 * t340;
t264 = t304 * t340 + t305 * t343;
t299 = -t338 * t389 + t390;
t362 = -qJ(4) * t299 + qJD(4) * t305;
t359 = qJ(4) * t305 + t324;
t357 = t324 + t376;
t355 = g(1) * t411 - g(3) * t327 + t326 * t423 - t383;
t246 = t305 * qJD(2) + qJD(3) * t268;
t354 = -g(2) * t410 - qJD(3) * t246 - qJDD(3) * t267 + t327 * t425;
t353 = t381 * t396 - t378;
t352 = pkin(3) * t262 + t442;
t244 = pkin(3) * t295 + t441;
t349 = t244 * t297 + qJDD(4) - t355;
t289 = t295 ^ 2;
t260 = pkin(3) * t304 - t359;
t259 = pkin(3) * t297 + t420;
t258 = t336 + t266;
t257 = -qJD(3) * pkin(3) + t397;
t241 = t345 * t304 + t359;
t240 = pkin(3) * t300 - t362;
t239 = (t295 - t391) * qJD(3) + t393;
t237 = t345 * t297 - t420;
t234 = t299 * pkin(7) + t246;
t233 = pkin(7) * t300 + t245;
t231 = t345 * t300 + t362;
t230 = t264 * qJD(5) - t299 * t340 - t343 * t300;
t229 = t365 * qJD(5) - t299 * t343 + t300 * t340;
t228 = t363 - t417;
t226 = t352 - t400;
t221 = t345 * t262 + t400 - t442;
t1 = [qJDD(1) * MDP(1) + t436 * MDP(2) + t378 * MDP(3) + (t381 * t433 + t353) * MDP(6) + (t358 * pkin(1) + (t404 * t433 + t353) * qJ(2)) * MDP(7) + (-t261 * t305 - t297 * t299) * MDP(8) + (t261 * t304 - t262 * t305 + t295 * t299 - t297 * t300) * MDP(9) + (-qJD(3) * t299 + qJDD(3) * t305) * MDP(10) + (-qJD(3) * t300 - qJDD(3) * t304) * MDP(11) + (-t262 * t324 + t300 * t309 + t304 * t308 + t354) * MDP(13) + (t261 * t324 - t299 * t309 + t305 * t308 + t431) * MDP(14) + (t226 * t304 + t240 * t295 + t244 * t300 + t260 * t262 + t354) * MDP(15) + (-t227 * t304 + t228 * t305 - t245 * t295 + t246 * t297 - t257 * t299 - t258 * t300 - t261 * t267 - t262 * t268 - t378) * MDP(16) + (-t226 * t305 - t240 * t297 + t244 * t299 + t260 * t261 - t431) * MDP(17) + (t226 * t260 + t227 * t268 + t228 * t267 + t244 * t240 + t258 * t245 + t257 * t246 + (-g(1) * t422 - g(2) * t357) * t344 + (g(1) * t357 - g(2) * t422) * t342) * MDP(18) + (t229 * t435 - t264 * t382) * MDP(19) + (-t225 * t264 - t229 * t367 - t230 * t435 - t365 * t382) * MDP(20) + (-t229 * t333 - t264 * t329) * MDP(21) + (t230 * t333 - t329 * t365) * MDP(22) + (t231 * t367 + t241 * t225 - t221 * t365 + t232 * t230 - (-t368 * qJD(5) - t233 * t340 + t234 * t343) * t333 - t369 * t329 + g(1) * t274 - g(2) * t276) * MDP(24) + (t231 * t435 - t241 * t382 + t221 * t264 + t232 * t229 + (t369 * qJD(5) + t233 * t343 + t234 * t340) * t333 + t368 * t329 + g(1) * t273 + g(2) * t275) * MDP(25) + t434 * (t358 + t418); -MDP(6) * t386 + (-qJ(2) * t386 - t358) * MDP(7) + (-t289 - t427) * MDP(16) + (t258 * t295 + (-qJD(4) - t257) * t297 + t352 - t436) * MDP(18) + (-t225 + t416) * MDP(24) + (t382 - t415) * MDP(25) + (MDP(13) + MDP(15)) * (0.2e1 * qJD(3) * t297 + t374) + (-MDP(14) + MDP(17)) * ((t295 + t391) * qJD(3) - t393) - t434 * qJDD(1); MDP(8) * t414 + (-t289 + t427) * MDP(9) + t239 * MDP(10) - t374 * MDP(11) + qJDD(3) * MDP(12) + (-t297 * t309 + t355 + t402) * MDP(13) + (t295 * t309 + t428) * MDP(14) + (-t259 * t295 - t349 + t402 + 0.2e1 * t417) * MDP(15) + (pkin(3) * t261 - qJ(4) * t262 + (t258 - t266) * t297 + (t257 - t397) * t295) * MDP(16) + (-t244 * t295 + t259 * t297 + 0.2e1 * t334 + 0.2e1 * t335 - t428) * MDP(17) + (-t228 * pkin(3) - g(3) * t376 + t227 * qJ(4) - t244 * t259 - t257 * t266 + t397 * t258 + t378 * (pkin(3) * t326 - qJ(4) * t327)) * MDP(18) + (t225 + t416) * MDP(22) + (-t372 * t329 - t237 * t367 + (t343 * t243 + t398 * t340) * t333 + (t373 * t333 - t370) * qJD(5) + t438) * MDP(24) + (t373 * t329 - t237 * t435 + (-t340 * t243 + t398 * t343) * t333 + (t372 * t333 + t371) * qJD(5) + t432) * MDP(25) - t443; (-qJDD(3) + t414) * MDP(15) + t239 * MDP(16) + (-qJD(3) ^ 2 - t427) * MDP(17) + (-qJD(3) * t258 + t349 - t417) * MDP(18) + (-t297 * t367 - t343 * t329) * MDP(24) + (-t297 * t435 + t340 * t329) * MDP(25) + (-MDP(24) * t340 - MDP(25) * t343) * t333 ^ 2; (t385 - t416) * MDP(22) + (t370 * t333 - t438) * MDP(24) + (-t371 * t333 - t432) * MDP(25) + (-MDP(22) * t435 + t370 * MDP(24) - t371 * MDP(25)) * qJD(5) + t443;];
tau = t1;
