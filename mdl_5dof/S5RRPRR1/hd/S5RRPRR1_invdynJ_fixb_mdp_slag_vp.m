% Calculate vector of inverse dynamics joint torques for
% S5RRPRR1
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
% pkin [4x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a3,a4,d4,d5]';
% MDP [26x1]
%   Minimal dynamic parameter vector (fixed base model)
%   see S5RRPRR1_convert_par2_MPV_fixb.m
% 
% Output:
% tau [5x1]
%   joint torques of inverse dynamics (contains inertial, gravitational Coriolis and centrifugal forces)

% Quelle: HybrDyn-Toolbox
% Datum: 2019-07-18 17:22
% Revision: 08c8d617a845f5dd194efdf9aca2774760f7818f (2019-07-16)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tau = S5RRPRR1_invdynJ_fixb_mdp_slag_vp(qJ, qJD, qJDD, g, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(5,1),zeros(3,1),zeros(4,1),zeros(26,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRPRR1_invdynJ_fixb_mdp_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RRPRR1_invdynJ_fixb_mdp_slag_vp: qJD has to be [5x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [5 1]), ...
  'S5RRPRR1_invdynJ_fixb_mdp_slag_vp: qJDD has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RRPRR1_invdynJ_fixb_mdp_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [4 1]), ...
  'S5RRPRR1_invdynJ_fixb_mdp_slag_vp: pkin has to be [4x1] (double)');
assert(isreal(MDP) && all(size(MDP) == [26 1]), ...
  'S5RRPRR1_invdynJ_fixb_mdp_slag_vp: MDP has to be [26x1] (double)'); 

%% Symbolic Calculation
% From invdyn_joint_fixb_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2019-07-18 17:22:35
% EndTime: 2019-07-18 17:22:40
% DurationCPUTime: 3.32s
% Computational Cost: add. (1783->320), mult. (3988->431), div. (0->0), fcn. (2736->10), ass. (0->154)
t315 = qJ(2) + qJ(4);
t309 = sin(t315);
t321 = sin(qJ(1));
t324 = cos(qJ(1));
t349 = g(1) * t324 + g(2) * t321;
t414 = t349 * t309;
t320 = sin(qJ(2));
t319 = sin(qJ(4));
t323 = cos(qJ(2));
t387 = t319 * t323;
t407 = cos(qJ(4));
t270 = t407 * t320 + t387;
t312 = qJD(2) + qJD(4);
t245 = t312 * t270;
t360 = qJDD(1) * t407;
t373 = qJDD(1) * t320;
t345 = t319 * t373 - t323 * t360;
t229 = t245 * qJD(1) + t345;
t325 = pkin(1) + pkin(2);
t366 = qJD(1) * t407;
t296 = t323 * t366;
t381 = qJD(1) * t320;
t265 = t319 * t381 - t296;
t264 = qJD(5) + t265;
t358 = t264 ^ 2;
t413 = pkin(4) * t358;
t402 = pkin(3) + qJ(3);
t278 = t402 * t320;
t271 = qJD(1) * t278;
t317 = qJD(2) * pkin(1);
t256 = qJD(2) * pkin(2) - t271 + t317;
t279 = t402 * t323;
t272 = qJD(1) * t279;
t367 = t407 * t272;
t238 = t319 * t256 + t367;
t361 = qJD(2) * t402;
t262 = -qJD(3) * t320 - t323 * t361;
t316 = qJDD(2) * pkin(1);
t236 = qJDD(2) * pkin(2) + t262 * qJD(1) - qJDD(1) * t278 + t316;
t261 = qJD(3) * t323 - t320 * t361;
t240 = t261 * qJD(1) + qJDD(1) * t279;
t359 = -t407 * t236 + t319 * t240;
t214 = t238 * qJD(4) + t359;
t310 = cos(t315);
t404 = g(3) * t310;
t412 = t214 + t404;
t280 = t325 * t323;
t268 = -qJD(1) * t280 + qJD(3);
t302 = g(3) * t309;
t410 = t268 * t265 + t349 * t310 + t302;
t338 = -t407 * t278 - t319 * t279;
t223 = t338 * qJD(4) + t407 * t261 + t319 * t262;
t227 = qJDD(5) + t229;
t368 = t407 * t256;
t388 = t319 * t272;
t237 = -t368 + t388;
t337 = -t319 * t320 + t407 * t323;
t244 = t312 * t337;
t252 = -t319 * t278 + t407 * t279;
t254 = -pkin(4) * t270 - t280;
t311 = qJDD(2) + qJDD(4);
t378 = qJD(4) * t319;
t347 = t407 * t240 - t272 * t378;
t365 = t407 * qJD(4);
t330 = t319 * t236 + t256 * t365 + t347;
t212 = t311 * pkin(4) + t330;
t380 = qJD(1) * t323;
t267 = -t319 * t380 - t320 * t366;
t243 = pkin(4) * t267 + t268;
t356 = qJD(5) * t243 + t212;
t409 = t214 * t270 - t252 * t227 + t237 * t244 - (qJD(5) * t254 + t223) * t264 + t356 * t337;
t314 = t323 ^ 2;
t408 = 0.2e1 * t314;
t403 = g(3) * t323;
t372 = qJDD(1) * t323;
t351 = t312 * t296 + t319 * t372 + t320 * t360;
t352 = t319 * t312;
t228 = -t352 * t381 + t351;
t318 = sin(qJ(5));
t322 = cos(qJ(5));
t376 = qJD(5) * t322;
t370 = t322 * t228 + t318 * t311 + t312 * t376;
t377 = qJD(5) * t318;
t215 = t267 * t377 + t370;
t401 = t215 * t318;
t400 = t237 * t270;
t248 = -t267 * t318 - t322 * t312;
t399 = t248 * t264;
t398 = t248 * t267;
t251 = -t267 * t322 + t312 * t318;
t397 = t251 * t264;
t396 = t251 * t267;
t395 = t254 * t227;
t394 = t265 * t312;
t393 = t267 * t312;
t391 = t318 * t227;
t390 = t318 * t321;
t389 = t318 * t324;
t327 = qJD(1) ^ 2;
t386 = t320 * t327;
t385 = t321 * t322;
t225 = t322 * t227;
t384 = t322 * t324;
t273 = t325 * qJD(2) * t320;
t274 = t325 * t381;
t313 = t320 ^ 2;
t383 = t313 - t314;
t382 = t313 + t408;
t375 = qJD(1) * qJD(2);
t374 = qJD(1) * qJD(3);
t364 = t320 * t375;
t371 = pkin(1) * t364 + qJDD(3);
t369 = t325 * t407;
t363 = t323 * t375;
t357 = t264 * t322;
t247 = pkin(2) * t364 - qJDD(1) * t280 + t371;
t220 = -pkin(4) * t228 + t247;
t232 = t312 * pkin(4) + t238;
t355 = qJD(5) * t232 - t220;
t293 = t319 * t325 + pkin(4);
t353 = pkin(4) * t265 + qJD(5) * t293 + t274;
t350 = t320 * t363;
t348 = g(1) * t321 - g(2) * t324;
t241 = -t319 * t271 + t367;
t346 = t325 * t378 - t241;
t344 = -t227 * t293 + t237 * t265;
t219 = t232 * t322 + t243 * t318;
t343 = -t219 * t267 + t237 * t376 + t412 * t318;
t218 = -t232 * t318 + t243 * t322;
t342 = t218 * t267 + t237 * t377 + t322 * t414;
t341 = t225 + (-t265 * t318 - t377) * t264;
t242 = -t407 * t271 - t388;
t340 = -t325 * t365 + t242;
t336 = t244 * t322 - t270 * t377;
t335 = -pkin(4) * t227 + (-t264 + t265) * t237;
t277 = -qJ(3) * t381 + t317;
t333 = -qJD(2) * t277 * t323 - t349;
t332 = pkin(1) * t372 + t348 - t371;
t331 = t268 * t267 - t359 - t404 + t414;
t298 = t322 * t311;
t216 = t251 * qJD(5) + t228 * t318 - t298;
t329 = ((t215 - t399) * t322 + (-t216 - t397) * t318) * MDP(21) + (t251 * t357 + t401) * MDP(20) + (t341 - t398) * MDP(23) + (t264 * t357 + t391 + t396) * MDP(22) + (t228 + t394) * MDP(15) + (-t229 - t393) * MDP(16) + (-t265 ^ 2 + t267 ^ 2) * MDP(14) + t311 * MDP(17) + (-MDP(13) * t265 + t264 * MDP(24)) * t267;
t328 = qJ(3) ^ 2;
t326 = qJD(2) ^ 2;
t289 = -pkin(1) * t380 + qJD(3);
t260 = t310 * t384 + t390;
t259 = -t310 * t389 + t385;
t258 = -t310 * t385 + t389;
t257 = t310 * t390 + t384;
t253 = -t320 * t374 + t316 + (-t363 - t373) * qJ(3);
t231 = -pkin(4) * t244 + t273;
t224 = t252 * qJD(4) + t319 * t261 - t407 * t262;
t217 = t322 * t220;
t1 = [qJDD(1) * MDP(1) + (qJDD(1) * t313 + 0.2e1 * t350) * MDP(4) + 0.2e1 * (t320 * t372 - t383 * t375) * MDP(5) + (qJDD(2) * t320 + t323 * t326) * MDP(6) + (qJDD(2) * t323 - t320 * t326) * MDP(7) + (-t253 * t320 + t382 * t374 + (t382 * qJDD(1) - 0.2e1 * t350) * qJ(3) + t333) * MDP(11) + (t314 * qJDD(1) * t328 + t332 * t323 * pkin(1) + (t374 * t408 + t333) * qJ(3) + (-t253 * qJ(3) - t277 * qJD(3) + (pkin(1) * t289 - 0.2e1 * t328 * t380) * qJD(2)) * t320) * MDP(12) + (t228 * t270 - t244 * t267) * MDP(13) + (t228 * t337 - t229 * t270 - t244 * t265 + t245 * t267) * MDP(14) + (t244 * t312 + t270 * t311) * MDP(15) + (-t245 * t312 + t311 * t337) * MDP(16) + (-t224 * t312 - t229 * t280 + t245 * t268 - t247 * t337 + t265 * t273 + t311 * t338) * MDP(18) + (-t223 * t312 - t228 * t280 + t244 * t268 + t247 * t270 - t252 * t311 - t267 * t273) * MDP(19) + (t215 * t270 * t322 + t336 * t251) * MDP(20) + ((-t248 * t322 - t251 * t318) * t244 + (-t401 - t216 * t322 + (t248 * t318 - t251 * t322) * qJD(5)) * t270) * MDP(21) + (-t215 * t337 + t270 * t225 + t245 * t251 + t336 * t264) * MDP(22) + (-t270 * t391 + t216 * t337 - t245 * t248 + (-t244 * t318 - t270 * t376) * t264) * MDP(23) + (-t227 * t337 + t245 * t264) * MDP(24) + (-g(1) * t258 - g(2) * t260 - t338 * t216 - t217 * t337 + t218 * t245 + t224 * t248 + (t395 + t231 * t264 + (t232 * t337 - t252 * t264 + t400) * qJD(5)) * t322 + t409 * t318) * MDP(25) + (-g(1) * t257 - g(2) * t259 - t338 * t215 - t219 * t245 + t224 * t251 + (-(-qJD(5) * t252 + t231) * t264 - t395 - t355 * t337 - qJD(5) * t400) * t318 + t409 * t322) * MDP(26) + t349 * MDP(3) + (-t320 * MDP(10) + t310 * MDP(18) - t309 * MDP(19) + t323 * MDP(9) + MDP(2)) * t348; t329 + (g(3) * t320 + t349 * t323) * MDP(10) + (t311 * t369 + t241 * t312 - t274 * t265 + (-t367 + (-t312 * t325 - t256) * t319) * qJD(4) + t331) * MDP(18) - t323 * MDP(4) * t386 + qJDD(2) * MDP(8) + (-t216 * t369 - t412 * t322 + t344 * t318 + t346 * t248 + (t340 * t318 - t353 * t322) * t264 + t342) * MDP(25) + MDP(7) * t372 + MDP(6) * t373 + ((qJ(3) * qJD(1) * t277 + t328 * t386) * t323 + (-t403 + t253 + (-qJD(1) * t289 + t349) * t320) * pkin(1)) * MDP(12) + (t349 * t320 - t403) * MDP(9) + (-t215 * t369 + t344 * t322 - t318 * t414 + t346 * t251 + (t353 * t318 + t340 * t322) * t264 + t343) * MDP(26) + (t242 * t312 + t274 * t267 + (-t311 * t325 - t236) * t319 + (-t312 * t369 - t368) * qJD(4) - t347 + t410) * MDP(19) + t383 * MDP(5) * t327 + (-pkin(1) * t373 + (qJ(3) * t386 + (t277 - t317) * qJD(1)) * t323) * MDP(11); -t332 * MDP(12) + (t345 - t393) * MDP(18) + (t351 - t394) * MDP(19) + (t341 + t398) * MDP(25) + (-t322 * t358 - t391 + t396) * MDP(26) + ((-t313 - t314) * MDP(11) - qJ(3) * t314 * MDP(12)) * t327 + (t312 * MDP(18) * t387 + (t277 * MDP(12) + (t407 * qJD(2) + t365) * MDP(18) - MDP(19) * t352) * t320) * qJD(1); t329 + (-t237 * t312 - t330 + t410) * MDP(19) + (t331 + (-qJD(4) + t312) * t238) * MDP(18) + (-t238 * t248 + t335 * t318 + (-t412 - t413) * t322 + t342) * MDP(25) + (-t238 * t251 + t335 * t322 + (-t414 + t413) * t318 + t343) * MDP(26); t251 * t248 * MDP(20) + (-t248 ^ 2 + t251 ^ 2) * MDP(21) + (t370 + t399) * MDP(22) + (t298 + t397) * MDP(23) + t227 * MDP(24) + (-g(1) * t259 + g(2) * t257 + t219 * t264 - t237 * t251 + t217) * MDP(25) + (g(1) * t260 - g(2) * t258 + t218 * t264 + t237 * t248) * MDP(26) + ((-t212 + t302) * MDP(26) + (MDP(23) * t267 - MDP(25) * t232 - MDP(26) * t243) * qJD(5)) * t322 + (qJD(5) * t267 * MDP(22) + (-qJD(5) * t312 - t228) * MDP(23) + (-t356 + t302) * MDP(25) + t355 * MDP(26)) * t318;];
tau  = t1;
