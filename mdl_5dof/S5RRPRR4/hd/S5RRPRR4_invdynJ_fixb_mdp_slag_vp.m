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
% Datum: 2019-12-05 18:32
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
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
% StartTime: 2019-12-05 18:32:13
% EndTime: 2019-12-05 18:32:16
% DurationCPUTime: 1.73s
% Computational Cost: add. (1409->236), mult. (2244->318), div. (0->0), fcn. (1448->16), ass. (0->145)
t326 = sin(qJ(5));
t327 = sin(qJ(4));
t330 = cos(qJ(5));
t331 = cos(qJ(4));
t264 = t326 * t331 + t327 * t330;
t319 = qJD(1) + qJD(2);
t250 = t264 * t319;
t332 = cos(qJ(2));
t398 = pkin(1) * qJD(2);
t372 = qJD(1) * t398;
t328 = sin(qJ(2));
t376 = qJDD(1) * t328;
t412 = pkin(1) * t376 + t332 * t372;
t385 = t330 * t331;
t389 = t326 * t327;
t263 = -t385 + t389;
t323 = qJ(1) + qJ(2);
t308 = pkin(9) + t323;
t294 = sin(t308);
t295 = cos(t308);
t411 = -g(2) * t294 + g(3) * t295;
t410 = g(2) * t295 + g(3) * t294;
t311 = sin(t323);
t313 = cos(t323);
t409 = g(2) * t313 + g(3) * t311;
t399 = pkin(1) * qJD(1);
t374 = t332 * t399;
t271 = pkin(2) * t319 + t374;
t325 = cos(pkin(9));
t375 = t328 * t399;
t291 = t325 * t375;
t324 = sin(pkin(9));
t243 = t324 * t271 + t291;
t363 = t243 + (pkin(7) + pkin(8)) * t319;
t228 = t331 * qJD(3) - t363 * t327;
t318 = qJD(4) + qJD(5);
t229 = qJD(3) * t327 + t331 * t363;
t408 = pkin(1) * t332;
t407 = pkin(2) * t325;
t406 = pkin(4) * t331;
t304 = pkin(2) + t408;
t390 = t325 * t328;
t382 = pkin(1) * t390 + t324 * t304;
t253 = pkin(7) + t382;
t401 = -pkin(8) - t253;
t296 = pkin(2) * t324 + pkin(7);
t400 = -pkin(8) - t296;
t397 = t229 * t330;
t290 = t324 * t375;
t256 = t325 * t374 - t290;
t396 = t256 * t318;
t322 = qJ(4) + qJ(5);
t312 = cos(t322);
t395 = t294 * t312;
t394 = t295 * t312;
t317 = qJDD(1) + qJDD(2);
t393 = t317 * t331;
t392 = t319 * t327;
t391 = t324 * t328;
t386 = t327 * t331;
t384 = qJDD(3) - g(1);
t305 = qJDD(1) * t408;
t251 = pkin(2) * t317 - t328 * t372 + t305;
t230 = t251 * t325 - t412 * t324;
t226 = -pkin(3) * t317 - t230;
t242 = t271 * t325 - t290;
t238 = -pkin(3) * t319 - t242;
t378 = qJD(4) * t331;
t383 = t226 * t327 + t238 * t378;
t320 = t327 ^ 2;
t381 = -t331 ^ 2 + t320;
t379 = qJD(4) * t327;
t377 = qJD(5) * t326;
t373 = pkin(4) * t379;
t370 = t319 * t389;
t369 = t319 * t385;
t368 = t238 * t379 + t410 * t331;
t231 = t324 * t251 + t412 * t325;
t367 = t305 + t409;
t366 = -pkin(3) - t406;
t365 = t319 * t378;
t227 = pkin(7) * t317 + t231;
t364 = pkin(8) * t317 + t227;
t362 = -g(2) * t311 + g(3) * t313;
t361 = qJD(4) * t401;
t360 = qJD(4) * t400;
t359 = -pkin(1) * t391 + t304 * t325;
t358 = qJD(1) * (-qJD(2) + t319);
t357 = qJD(2) * (-qJD(1) - t319);
t214 = (t319 * t379 - t393) * pkin(4) + t226;
t232 = t319 * t366 - t242;
t237 = t318 * t264;
t355 = g(2) * t394 + g(3) * t395 + t214 * t263 + t232 * t237;
t252 = -pkin(3) - t359;
t254 = t324 * t374 + t291;
t354 = -t254 + t373;
t352 = t263 * t317;
t225 = qJD(4) * pkin(4) + t228;
t350 = -t225 * t326 - t397;
t236 = t318 * t263;
t316 = qJDD(4) + qJDD(5);
t349 = -t236 * t318 + t264 * t316;
t240 = t401 * t327;
t314 = t331 * pkin(8);
t241 = t253 * t331 + t314;
t348 = t240 * t330 - t241 * t326;
t347 = t240 * t326 + t241 * t330;
t261 = t400 * t327;
t262 = t296 * t331 + t314;
t346 = t261 * t330 - t262 * t326;
t345 = t261 * t326 + t262 * t330;
t255 = (t324 * t332 + t390) * t398;
t334 = qJD(4) ^ 2;
t344 = -t252 * t317 - t253 * t334 - t255 * t319;
t297 = -pkin(3) - t407;
t343 = t254 * t319 - t296 * t334 - t297 * t317;
t215 = qJD(5) * t369 + t264 * t317 - t318 * t370 + t330 * t365;
t248 = -t369 + t370;
t342 = t250 * t248 * MDP(15) + (t248 * t318 + t215) * MDP(17) - t352 * MDP(18) + (-t248 ^ 2 + t250 ^ 2) * MDP(16) + t316 * MDP(19);
t341 = -t238 * t319 - t227 + t411;
t257 = (t325 * t332 - t391) * t398;
t340 = -qJDD(4) * t253 + (t252 * t319 - t257) * qJD(4);
t339 = -qJDD(4) * t296 + (t297 * t319 + t256) * qJD(4);
t310 = sin(t322);
t338 = t214 * t264 - t232 * t236 - t310 * t410;
t216 = t237 * t319 + t352;
t220 = -t237 * t318 - t263 * t316;
t275 = qJDD(4) * t327 + t331 * t334;
t276 = qJDD(4) * t331 - t327 * t334;
t337 = (-t215 * t263 - t216 * t264 + t236 * t248 - t237 * t250) * MDP(16) + (t215 * t264 - t236 * t250) * MDP(15) + t349 * MDP(17) + t220 * MDP(18) + 0.2e1 * (-qJD(4) * t319 * t381 + t317 * t386) * MDP(9) + (t317 * t320 + 0.2e1 * t327 * t365) * MDP(8) + t275 * MDP(10) + t276 * MDP(11) + t317 * MDP(4);
t307 = t331 * qJDD(3);
t206 = qJDD(4) * pkin(4) - t229 * qJD(4) - t364 * t327 + t307;
t336 = t229 * t377 + g(3) * t394 + g(1) * t310 + t232 * t248 + (-t229 * t318 - t206) * t326 - g(2) * t395;
t207 = t228 * qJD(4) + t327 * qJDD(3) + t364 * t331;
t335 = -g(1) * t312 + qJD(5) * t350 + t330 * t206 - t326 * t207 - t232 * t250 + t411 * t310;
t333 = cos(qJ(1));
t329 = sin(qJ(1));
t274 = t366 - t407;
t259 = t331 * t360;
t258 = t327 * t360;
t246 = t252 - t406;
t244 = t255 + t373;
t224 = -t257 * t327 + t331 * t361;
t223 = t257 * t331 + t327 * t361;
t1 = [qJDD(1) * MDP(1) + (t244 * t250 + t246 * t215 - (qJD(5) * t348 + t223 * t330 + t224 * t326) * t318 - t347 * t316 + t338) * MDP(21) + t337 + (((-qJDD(1) - t317) * t328 + t332 * t357) * pkin(1) + t362) * MDP(6) + (t340 * t327 + (-t226 + t344) * t331 + t368) * MDP(13) + ((t317 * t332 + t328 * t357) * pkin(1) + t367) * MDP(5) + (t244 * t248 + t246 * t216 + (-qJD(5) * t347 - t223 * t326 + t224 * t330) * t318 + t348 * t316 + t355) * MDP(20) + (t340 * t331 + (-t344 - t410) * t327 + t383) * MDP(14) + (g(2) * t333 + g(3) * t329) * MDP(2) + (-g(2) * t329 + g(3) * t333) * MDP(3) + (t231 * t382 + t243 * t257 + t230 * t359 - t242 * t255 - g(2) * (-pkin(1) * t333 - pkin(2) * t313) - g(3) * (-pkin(1) * t329 - pkin(2) * t311)) * MDP(7); ((t332 * t358 - t376) * pkin(1) + t362) * MDP(6) + t337 + (t339 * t327 + (-t226 + t343) * t331 + t368) * MDP(13) + (t339 * t331 + (-t343 - t410) * t327 + t383) * MDP(14) + (t274 * t216 + (-qJD(5) * t345 - t258 * t326 + t259 * t330) * t318 + t346 * t316 + t264 * t396 + t354 * t248 + t355) * MDP(20) + (t274 * t215 - (qJD(5) * t346 + t258 * t330 + t259 * t326) * t318 - t345 * t316 - t263 * t396 + t354 * t250 + t338) * MDP(21) + (t242 * t254 - t243 * t256 + (t230 * t325 + t231 * t324 + t409) * pkin(2)) * MDP(7) + (pkin(1) * t328 * t358 + t367) * MDP(5); t276 * MDP(13) - t275 * MDP(14) + t220 * MDP(20) - MDP(21) * t349 + MDP(7) * t384; t327 * t317 * MDP(10) + MDP(11) * t393 + qJDD(4) * MDP(12) + (-g(1) * t331 + t327 * t341 + t307) * MDP(13) + (-t327 * t384 + t331 * t341) * MDP(14) + (-(-t228 * t326 - t397) * t318 + (-t248 * t392 + t330 * t316 - t318 * t377) * pkin(4) + t335) * MDP(20) + ((-qJD(5) * t225 + t228 * t318 - t207) * t330 + (-qJD(5) * t330 * t318 - t250 * t392 - t326 * t316) * pkin(4) + t336) * MDP(21) + t342 + (-MDP(8) * t386 + MDP(9) * t381) * t319 ^ 2; (-t318 * t350 + t335) * MDP(20) + ((-t207 + (-qJD(5) + t318) * t225) * t330 + t336) * MDP(21) + t342;];
tau = t1;
