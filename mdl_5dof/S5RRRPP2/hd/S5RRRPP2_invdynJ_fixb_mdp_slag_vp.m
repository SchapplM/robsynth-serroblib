% Calculate vector of inverse dynamics joint torques for
% S5RRRPP2
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
%   pkin=[a2,a3,a4,a5,d1,d2,d3]';
% MDP [21x1]
%   Minimal dynamic parameter vector (fixed base model)
%   see S5RRRPP2_convert_par2_MPV_fixb.m
% 
% Output:
% tau [5x1]
%   joint torques of inverse dynamics (contains inertial, gravitational Coriolis and centrifugal forces)

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 20:52
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tau = S5RRRPP2_invdynJ_fixb_mdp_slag_vp(qJ, qJD, qJDD, g, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(5,1),zeros(3,1),zeros(7,1),zeros(21,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRRPP2_invdynJ_fixb_mdp_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RRRPP2_invdynJ_fixb_mdp_slag_vp: qJD has to be [5x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [5 1]), ...
  'S5RRRPP2_invdynJ_fixb_mdp_slag_vp: qJDD has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RRRPP2_invdynJ_fixb_mdp_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [7 1]), ...
  'S5RRRPP2_invdynJ_fixb_mdp_slag_vp: pkin has to be [7x1] (double)');
assert(isreal(MDP) && all(size(MDP) == [21 1]), ...
  'S5RRRPP2_invdynJ_fixb_mdp_slag_vp: MDP has to be [21x1] (double)'); 

%% Symbolic Calculation
% From invdyn_joint_fixb_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 20:52:05
% EndTime: 2019-12-31 20:52:09
% DurationCPUTime: 2.71s
% Computational Cost: add. (1724->367), mult. (2401->409), div. (0->0), fcn. (1173->8), ass. (0->175)
t332 = cos(qJ(3));
t419 = qJ(4) * t332;
t329 = sin(qJ(3));
t428 = pkin(3) + pkin(4);
t431 = t428 * t329;
t434 = t419 - t431;
t313 = t329 * qJ(4);
t373 = pkin(2) + t313;
t392 = qJD(3) * t332;
t433 = -qJ(4) * t392 - t329 * qJD(4);
t432 = t428 * qJDD(3);
t322 = qJD(1) + qJD(2);
t330 = sin(qJ(2));
t421 = pkin(1) * qJD(1);
t386 = t330 * t421;
t272 = pkin(7) * t322 + t386;
t261 = t329 * t272;
t408 = t322 * t329;
t241 = qJ(5) * t408 - t261;
t390 = qJD(4) - t241;
t232 = -t428 * qJD(3) + t390;
t407 = t322 * t332;
t321 = qJDD(1) + qJDD(2);
t388 = qJDD(1) * t330;
t333 = cos(qJ(2));
t396 = qJD(2) * t333;
t251 = pkin(7) * t321 + (qJD(1) * t396 + t388) * pkin(1);
t246 = t332 * t251;
t323 = qJDD(3) * qJ(4);
t324 = qJD(3) * qJD(4);
t393 = qJD(3) * t329;
t226 = -t272 * t393 + t246 + t323 + t324;
t245 = t329 * t251;
t258 = t272 * t392;
t374 = qJDD(4) + t245 + t258;
t417 = qJDD(3) * pkin(3);
t227 = t374 - t417;
t430 = t226 * t332 + t227 * t329;
t262 = t332 * t272;
t242 = -qJ(5) * t407 + t262;
t325 = qJD(3) * qJ(4);
t238 = t242 + t325;
t328 = qJ(1) + qJ(2);
t311 = sin(t328);
t312 = cos(t328);
t403 = g(1) * t312 + g(2) * t311;
t300 = g(1) * t311;
t423 = g(2) * t312;
t429 = t300 - t423;
t427 = pkin(1) * t333;
t426 = pkin(2) * t322;
t336 = qJD(3) ^ 2;
t425 = pkin(7) * t336;
t331 = sin(qJ(1));
t424 = g(1) * t331;
t309 = t321 * pkin(2);
t317 = t332 * pkin(3);
t316 = t332 * pkin(4);
t422 = pkin(7) - qJ(5);
t420 = pkin(7) * qJDD(3);
t418 = qJ(5) * t321;
t401 = t317 + t313;
t378 = t316 + t401;
t264 = pkin(2) + t378;
t416 = t264 * t322;
t304 = pkin(1) * t330 + pkin(7);
t415 = t304 * t336;
t414 = t311 * t329;
t413 = t311 * t332;
t412 = t312 * t329;
t411 = t312 * t332;
t410 = t321 * t329;
t409 = t321 * t332;
t406 = t329 * t332;
t405 = -qJ(5) + t304;
t404 = -g(1) * t414 + g(2) * t412;
t397 = qJD(2) * t330;
t384 = pkin(1) * t397;
t402 = -qJD(1) * t384 + qJDD(1) * t427;
t326 = t329 ^ 2;
t327 = t332 ^ 2;
t400 = t326 - t327;
t399 = t326 + t327;
t398 = qJD(1) * t333;
t266 = t405 * t332;
t395 = qJD(3) * t266;
t394 = qJD(3) * t322;
t391 = qJD(5) * t329;
t385 = pkin(1) * t398;
t229 = t385 + qJD(5) + (t428 * t332 + t373) * t322;
t389 = qJD(5) + t229;
t387 = qJDD(3) * t304;
t383 = pkin(1) * t396;
t320 = t322 ^ 2;
t382 = t320 * t406;
t381 = t238 * t393 + t403;
t380 = t246 + 0.2e1 * t323 + 0.2e1 * t324;
t289 = g(1) * t413;
t368 = t322 * t386;
t379 = t332 * t368 + t385 * t393 + t289;
t250 = -t309 - t402;
t256 = pkin(3) * t393 + t433;
t305 = -pkin(2) - t427;
t377 = t322 * t397;
t376 = t322 * t393;
t375 = t322 * t392;
t302 = qJ(5) * t393;
t278 = t422 * t332;
t372 = -t250 - t423;
t371 = t399 * t321;
t370 = -qJD(3) * pkin(3) + qJD(4);
t297 = t312 * pkin(7);
t369 = -qJ(5) * t312 + t297;
t347 = -pkin(3) * t409 - t321 * t313 + t433 * t322 + t250;
t342 = pkin(4) * t409 + qJDD(5) - t347;
t220 = -t428 * t376 + t342;
t367 = t220 * t329 + t229 * t392 - t404;
t366 = g(1) * t412 + g(2) * t414 - g(3) * t332 - t245;
t273 = -t385 - t426;
t365 = t250 * t329 + t273 * t392 + t404;
t364 = pkin(3) * t411 + t311 * pkin(7) + t373 * t312;
t363 = t399 * t427;
t362 = -qJD(5) + t383;
t361 = -t309 + t425;
t359 = pkin(3) * t329 - t419;
t358 = 0.2e1 * (t321 * t406 - t400 * t394) * MDP(8) + (t321 * t326 + 0.2e1 * t329 * t375) * MDP(7) + (qJDD(3) * t332 - t329 * t336) * MDP(10) + (qJDD(3) * t329 + t332 * t336) * MDP(9) + t321 * MDP(4);
t334 = cos(qJ(1));
t357 = t334 * pkin(1) + t364;
t263 = t305 - t401;
t356 = t232 * t329 + t238 * t332;
t243 = -pkin(4) * t393 - t256;
t236 = t243 - t384;
t252 = -t263 + t316;
t355 = t236 * t322 + t252 * t321;
t354 = t243 * t322 + t264 * t321;
t249 = t261 + t370;
t253 = t262 + t325;
t353 = t249 * t329 + t253 * t332;
t352 = -qJDD(4) + t366;
t351 = t402 + t429;
t350 = -t373 - t317;
t349 = t263 * t322 - t383;
t348 = -t258 + t352;
t223 = pkin(3) * t376 + t347;
t274 = -pkin(2) - t401;
t346 = -t256 * t322 - t274 * t321 - t223 - t425;
t244 = t256 + t384;
t345 = -t244 * t322 - t263 * t321 - t223 - t415;
t344 = t249 * t392 - t253 * t393 - t403 + t430;
t343 = -t348 - t417;
t341 = pkin(1) * t377 + t305 * t321 + t415;
t340 = -g(1) * t297 - t350 * t300;
t339 = -t387 + (t305 * t322 - t383) * qJD(3);
t338 = (t249 * t332 - t253 * t329) * qJD(3) + t430;
t337 = (-g(1) * (t350 - t316) + g(2) * qJ(5)) * t311;
t284 = pkin(4) * t411;
t283 = qJ(4) * t411;
t281 = qJ(4) * t413;
t279 = t322 * t302;
t277 = t422 * t329;
t270 = t329 * t368;
t265 = t405 * t329;
t259 = t273 * t393;
t257 = t359 * t322;
t255 = qJD(3) * t278 - t391;
t254 = -pkin(7) * t393 - qJD(5) * t332 + t302;
t240 = t434 * t322;
t237 = t350 * t322 - t385;
t235 = t362 * t329 + t395;
t234 = -t304 * t393 + t362 * t332 + t302;
t230 = t237 * t393;
t222 = t279 + (-qJD(5) * t322 - t418) * t332 + t226;
t221 = -t322 * t391 - t432 + (-t375 - t410) * qJ(5) + t374;
t219 = t220 * t332;
t1 = [((-t265 * t321 - t221 + (-t235 + t395) * t322) * t329 + (-t234 * t322 - t266 * t321 - t222 + (-t265 * t322 - t232) * qJD(3)) * t332 + t381) * MDP(20) + (((-qJDD(1) - t321) * t330 + (-qJD(1) - t322) * t396) * pkin(1) + t403) * MDP(6) + t358 + (t222 * t266 + t238 * t234 + t221 * t265 + t232 * t235 + t220 * t252 + t229 * t236 - g(1) * (-pkin(1) * t331 + t369) - g(2) * (t284 + t357) + t337) * MDP(21) + (t223 * t263 + t237 * t244 - g(2) * t357 + (t353 * t396 + t424) * pkin(1) + t338 * t304 + t340) * MDP(17) + (-qJDD(3) * t265 + t219 + t289 + (t355 - t423) * t332 + (-t235 + (-t252 * t322 - t229) * t329) * qJD(3)) * MDP(18) + (qJDD(3) * t266 + t355 * t329 + (t252 * t407 + t234) * qJD(3) + t367) * MDP(19) + ((t321 * t333 - t377) * pkin(1) + t351) * MDP(5) + ((t387 + (-t237 - t349) * qJD(3)) * t332 + t345 * t329 - t404) * MDP(16) + (t230 + t289 + (t349 * qJD(3) - t387) * t329 + (t345 - t423) * t332) * MDP(14) + qJDD(1) * MDP(1) + (t259 + t289 + t339 * t329 + (-t341 + t372) * t332) * MDP(12) + (t341 * t329 + t339 * t332 + t365) * MDP(13) + (t322 * qJD(2) * t363 + t304 * t371 + t344) * MDP(15) + (-g(2) * t334 + t424) * MDP(2) + (g(1) * t334 + g(2) * t331) * MDP(3); (t351 + t368) * MDP(5) + ((-t388 + (-qJD(2) + t322) * t398) * pkin(1) + t403) * MDP(6) + (t259 + (-pkin(2) * t394 - t420) * t329 + (-t361 + t372) * t332 + t379) * MDP(12) + (-t270 + t361 * t329 + (-t420 + (t385 - t426) * qJD(3)) * t332 + t365) * MDP(13) + (t230 + (t274 * t394 - t420) * t329 + (t346 - t423) * t332 + t379) * MDP(14) + (-t399 * t322 * t385 + pkin(7) * t371 + t344) * MDP(15) + (t270 + (t420 + (-t274 * t322 - t237 - t385) * qJD(3)) * t332 + t346 * t329 - t404) * MDP(16) + (t223 * t274 + t237 * t256 - g(2) * t364 + (-t237 * t330 - t353 * t333) * t421 + t338 * pkin(7) + t340) * MDP(17) + (-qJDD(3) * t277 + t219 + (t354 - t423) * t332 + (-t255 + (-t229 - t416) * t329) * qJD(3) + t379) * MDP(18) + (qJDD(3) * t278 + t270 + t354 * t329 + (t254 + (-t385 + t416) * t332) * qJD(3) + t367) * MDP(19) + ((-t277 * t321 - t221) * t329 + (-qJD(3) * t232 - t278 * t321 - t222) * t332 + (-t254 * t332 - t255 * t329 + (-t277 * t332 + t278 * t329) * qJD(3) + qJD(1) * t363) * t322 + t381) * MDP(20) + (t222 * t278 + t238 * t254 + t221 * t277 + t232 * t255 + t220 * t264 + t229 * t243 - g(1) * t369 - g(2) * (t284 + t364) + t337 + (t229 * t330 - t356 * t333) * t421) * MDP(21) + t358; -MDP(7) * t382 + t400 * MDP(8) * t320 + MDP(9) * t410 + MDP(10) * t409 + qJDD(3) * MDP(11) + (-t273 * t408 + t366) * MDP(12) + (g(3) * t329 - t246 + (-t273 * t322 + t403) * t332) * MDP(13) + (0.2e1 * t417 + (-t237 * t329 + t257 * t332) * t322 + t352) * MDP(14) + (-t359 * t321 + ((t253 - t325) * t329 + (-t249 + t370) * t332) * t322) * MDP(15) + ((t257 * t322 - g(3)) * t329 + (t237 * t322 - t403) * t332 + t380) * MDP(16) + (t226 * qJ(4) - t227 * pkin(3) - t237 * t257 - t249 * t262 - g(1) * (-pkin(3) * t412 + t283) - g(2) * (-pkin(3) * t414 + t281) - g(3) * t401 + (qJD(4) + t261) * t253) * MDP(17) + (qJ(5) * t410 + qJD(3) * t242 + 0.2e1 * t432 + ((qJ(5) * qJD(3) - t240) * t332 + t389 * t329) * t322 + t348) * MDP(18) + (-qJD(3) * t241 + t279 + (-qJD(3) * t272 - t240 * t322 - g(3)) * t329 + (-t389 * t322 - t403 - t418) * t332 + t380) * MDP(19) - t434 * t321 * MDP(20) + (-g(1) * t283 - g(2) * t281 - g(3) * t378 + t222 * qJ(4) - t221 * t428 - t229 * t240 - t232 * t242 + t390 * t238 + t403 * t431) * MDP(21); (-qJD(3) * t253 + t343) * MDP(17) + (-qJDD(3) * pkin(4) - qJ(5) * t375 - qJD(3) * t238 + t343) * MDP(21) + (MDP(16) + MDP(19)) * (-t320 * t326 - t336) + (MDP(14) + MDP(18)) * (-qJDD(3) - t382) + ((t237 * MDP(17) - t389 * MDP(21)) * t322 + (-MDP(21) * qJ(5) + MDP(15) - MDP(20)) * t321) * t329; (t342 + t429) * MDP(21) + (MDP(18) * t332 + t329 * MDP(19)) * t321 - t399 * MDP(20) * t320 + (t356 * MDP(21) + (0.2e1 * t332 * MDP(19) + (-t428 * MDP(21) - 0.2e1 * MDP(18)) * t329) * qJD(3)) * t322;];
tau = t1;
