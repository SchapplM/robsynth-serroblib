% Calculate vector of inverse dynamics joint torques for
% S5RPPPR2
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
%   pkin=[a2,a3,a4,a5,d1,d5,theta2,theta3,theta4]';
% MDP [22x1]
%   Minimal dynamic parameter vector (fixed base model)
%   see S5RPPPR2_convert_par2_MPV_fixb.m
% 
% Output:
% tau [5x1]
%   joint torques of inverse dynamics (contains inertial, gravitational Coriolis and centrifugal forces)

% Quelle: HybrDyn-Toolbox
% Datum: 2022-01-23 09:00
% Revision: fd3771346c4aea32fdeb66112c511235427c26a7 (2022-01-20)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tau = S5RPPPR2_invdynJ_fixb_mdp_slag_vp(qJ, qJD, qJDD, g, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(5,1),zeros(3,1),zeros(9,1),zeros(22,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPPPR2_invdynJ_fixb_mdp_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RPPPR2_invdynJ_fixb_mdp_slag_vp: qJD has to be [5x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [5 1]), ...
  'S5RPPPR2_invdynJ_fixb_mdp_slag_vp: qJDD has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RPPPR2_invdynJ_fixb_mdp_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5RPPPR2_invdynJ_fixb_mdp_slag_vp: pkin has to be [9x1] (double)');
assert(isreal(MDP) && all(size(MDP) == [22 1]), ...
  'S5RPPPR2_invdynJ_fixb_mdp_slag_vp: MDP has to be [22x1] (double)'); 

%% Symbolic Calculation
% From invdyn_joint_fixb_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2022-01-23 08:59:51
% EndTime: 2022-01-23 08:59:55
% DurationCPUTime: 3.75s
% Computational Cost: add. (1461->343), mult. (3635->510), div. (0->0), fcn. (2868->10), ass. (0->159)
t340 = sin(qJ(1));
t342 = cos(qJ(1));
t370 = g(1) * t340 - g(2) * t342;
t333 = sin(pkin(9));
t336 = cos(pkin(9));
t338 = cos(pkin(7));
t335 = sin(pkin(7));
t337 = cos(pkin(8));
t413 = t335 * t337;
t283 = t333 * t413 + t336 * t338;
t274 = t283 * qJD(1);
t268 = qJD(5) + t274;
t339 = sin(qJ(5));
t334 = sin(pkin(8));
t341 = cos(qJ(5));
t414 = t334 * t341;
t289 = t336 * t414 - t337 * t339;
t410 = t337 * t338;
t417 = t333 * t335;
t285 = t336 * t410 + t417;
t415 = t334 * t339;
t351 = t285 * t341 + t338 * t415;
t426 = -t342 * t289 + t351 * t340;
t422 = qJDD(1) * pkin(1);
t425 = t422 + t370;
t254 = -t285 * t339 + t338 * t414;
t288 = t336 * t415 + t337 * t341;
t424 = -t254 * t342 + t340 * t288;
t392 = qJ(2) * qJDD(1);
t401 = qJD(1) * t335;
t378 = t336 * t401;
t399 = qJD(1) * t338;
t379 = t333 * t399;
t277 = t337 * t378 - t379;
t420 = t277 * t339;
t418 = t333 * t334;
t416 = t334 * t335;
t412 = t336 * t340;
t411 = t336 * t342;
t409 = t338 * t340;
t408 = t338 * t342;
t343 = qJD(1) ^ 2;
t407 = t338 * t343;
t272 = t283 * qJDD(1);
t267 = qJDD(5) + t272;
t406 = t339 * t267;
t404 = t341 * t267;
t371 = t335 * qJ(3) + pkin(1);
t298 = pkin(2) * t338 + t371;
t397 = qJD(3) * t335;
t262 = -qJD(1) * t397 - qJDD(1) * t298 + qJDD(2);
t386 = qJDD(1) * t338;
t372 = t337 * t386;
t391 = qJD(1) * qJD(2);
t377 = t338 * t391;
t236 = qJ(2) * t372 + t334 * t262 + t337 * t377;
t390 = qJD(1) * qJD(4);
t231 = (-qJ(4) * qJDD(1) - t390) * t338 + t236;
t388 = qJDD(1) * t335;
t294 = qJ(2) * t388 + t335 * t391 + qJDD(3);
t360 = pkin(3) * t334 - qJ(4) * t337;
t240 = (qJDD(1) * t360 - t337 * t390) * t335 + t294;
t217 = t336 * t231 + t333 * t240;
t282 = -qJD(1) * t298 + qJD(2);
t383 = qJ(2) * t399;
t251 = t334 * t282 + t337 * t383;
t243 = -qJ(4) * t399 + t251;
t310 = qJ(2) * t401 + qJD(3);
t260 = t360 * t401 + t310;
t224 = t336 * t243 + t333 * t260;
t265 = qJ(2) * t410 - t334 * t298;
t257 = -qJ(4) * t338 + t265;
t353 = qJ(2) + t360;
t270 = t353 * t335;
t233 = t336 * t257 + t333 * t270;
t330 = t335 ^ 2;
t403 = t338 ^ 2 + t330;
t402 = qJD(1) * t334;
t400 = qJD(1) * t337;
t398 = qJD(2) * t338;
t396 = qJD(5) * t268;
t395 = qJD(5) * t339;
t394 = t267 * MDP(20);
t389 = qJDD(1) * t334;
t387 = qJDD(1) * t337;
t311 = t333 * t386;
t373 = t336 * t387;
t273 = t335 * t373 - t311;
t380 = t341 * t402;
t364 = t335 * t380;
t375 = t334 * t388;
t384 = qJD(5) * t364 + t341 * t273 + t339 * t375;
t382 = t334 * t401;
t381 = t339 * t402;
t329 = t334 ^ 2;
t376 = t329 * t388;
t374 = t334 * t386;
t369 = t403 * t343;
t215 = pkin(6) * t375 + t217;
t235 = -qJ(2) * t374 + t262 * t337 - t334 * t377;
t234 = pkin(3) * t386 + qJDD(4) - t235;
t219 = pkin(4) * t272 - pkin(6) * t273 + t234;
t368 = -t339 * t215 + t341 * t219;
t367 = t273 * t339 - t341 * t375;
t250 = t282 * t337 - t334 * t383;
t264 = -t334 * t338 * qJ(2) - t298 * t337;
t366 = (-t337 ^ 2 - t329) * MDP(10);
t365 = 0.2e1 * t403;
t290 = t334 * t409 + t337 * t342;
t292 = t334 * t408 - t337 * t340;
t363 = -g(1) * t290 + g(2) * t292;
t362 = g(1) * t342 + g(2) * t340;
t259 = t338 * pkin(3) - t264;
t287 = -t334 * t397 + t337 * t398;
t359 = t341 * t215 + t339 * t219;
t242 = pkin(3) * t399 + qJD(4) - t250;
t220 = pkin(4) * t274 - pkin(6) * t277 + t242;
t222 = pkin(6) * t382 + t224;
t358 = t220 * t341 - t222 * t339;
t357 = -t220 * t339 - t222 * t341;
t284 = -t338 * t333 + t336 * t413;
t227 = pkin(4) * t283 - pkin(6) * t284 + t259;
t229 = pkin(6) * t416 + t233;
t356 = t227 * t341 - t229 * t339;
t355 = t227 * t339 + t229 * t341;
t216 = -t231 * t333 + t240 * t336;
t223 = -t243 * t333 + t260 * t336;
t232 = -t257 * t333 + t270 * t336;
t354 = (-qJD(4) * t337 + qJD(2)) * t335;
t352 = -t284 * t339 + t335 * t414;
t253 = t284 * t341 + t335 * t415;
t350 = t365 * t391;
t286 = t334 * t398 + t337 * t397;
t349 = qJD(1) * t286 - qJDD(1) * t264 - t235;
t348 = qJD(1) * t287 + qJDD(1) * t265 + t236;
t246 = t277 * t341 + t335 * t381;
t347 = t289 * t268;
t345 = t294 * t335 + (t391 + t392) * t330;
t327 = g(3) * t338;
t326 = t342 * qJ(2);
t324 = t340 * qJ(2);
t323 = qJDD(2) - t422;
t293 = t334 * t340 + t337 * t408;
t291 = t334 * t342 - t337 * t409;
t279 = t285 * qJD(1);
t276 = t337 * t379 - t378;
t269 = (pkin(3) * t337 + qJ(4) * t334 + pkin(2)) * t338 + t371;
t266 = -qJD(4) * t338 + t287;
t248 = t253 * qJD(5);
t247 = t352 * qJD(5);
t244 = -t364 + t420;
t239 = t336 * t266 + t333 * t354;
t238 = t266 * t333 - t336 * t354;
t228 = -pkin(4) * t416 - t232;
t226 = qJD(5) * t246 + t367;
t225 = -t277 * t395 + t384;
t221 = -pkin(4) * t382 - t223;
t214 = -pkin(4) * t375 - t216;
t1 = [qJDD(1) * MDP(1) + t370 * MDP(2) + t362 * MDP(3) + (t365 * t392 + t350 - t362) * MDP(6) + (-t323 * pkin(1) - g(1) * (-pkin(1) * t340 + t326) - g(2) * (pkin(1) * t342 + t324) + (t392 * t403 + t350) * qJ(2)) * MDP(7) + (-g(1) * t291 - g(2) * t293 + t334 * t345 + t338 * t349) * MDP(8) + (t337 * t345 + t338 * t348 + t363) * MDP(9) + (-t334 * t348 + t337 * t349 + t370) * t335 * MDP(10) + (t236 * t265 + t251 * t287 + t235 * t264 - t250 * t286 - g(1) * (-t298 * t340 + t326) - g(2) * (t298 * t342 + t324) + (qJ(2) * t294 + qJD(2) * t310) * t335) * MDP(11) + (t286 * t274 + t259 * t272 + t234 * t283 - g(1) * (t291 * t336 - t340 * t417) - g(2) * (t293 * t336 + t342 * t417) + (-qJD(1) * t238 + qJDD(1) * t232 + t216) * t416) * MDP(12) + (t286 * t277 + t259 * t273 + t234 * t284 - g(1) * (-t291 * t333 - t335 * t412) - g(2) * (-t293 * t333 + t335 * t411) + (-qJD(1) * t239 - qJDD(1) * t233 - t217) * t416) * MDP(13) + (-t216 * t284 - t217 * t283 - t232 * t273 - t233 * t272 + t238 * t277 - t239 * t274 - t363) * MDP(14) + (t217 * t233 + t224 * t239 + t216 * t232 - t223 * t238 + t234 * t259 + t242 * t286 - g(1) * (-t269 * t340 + t342 * t353) - g(2) * (t269 * t342 + t340 * t353)) * MDP(15) + (t225 * t253 + t246 * t247) * MDP(16) + (t225 * t352 - t226 * t253 - t244 * t247 - t246 * t248) * MDP(17) + (t225 * t283 + t247 * t268 + t253 * t267) * MDP(18) + (-t226 * t283 - t248 * t268 + t267 * t352) * MDP(19) + t283 * t394 + ((-t239 * t339 + t286 * t341) * t268 + t356 * t267 + t368 * t283 + t238 * t244 + t228 * t226 - t214 * t352 + t221 * t248 + g(1) * t426 - g(2) * ((t285 * t342 + t334 * t412) * t341 + t292 * t339) + (-t268 * t355 + t283 * t357) * qJD(5)) * MDP(21) + (-(t239 * t341 + t286 * t339) * t268 - t355 * t267 - t359 * t283 + t238 * t246 + t228 * t225 + t214 * t253 + t221 * t247 - g(1) * (-t254 * t340 - t288 * t342) + g(2) * t424 + (-t268 * t356 - t283 * t358) * qJD(5)) * MDP(22) + (t338 * MDP(4) - t335 * MDP(5)) * (-t323 + t425); -MDP(4) * t386 - MDP(6) * t369 + (-qJ(2) * t369 + qJDD(2) - t425) * MDP(7) + (-t334 * t369 - t372) * MDP(8) + (-t337 * t369 + t374) * MDP(9) + (t235 * t337 + t236 * t334 + (-t310 * t335 + (t250 * t334 - t251 * t337) * t338) * qJD(1) - t370) * MDP(11) + (-t333 * t376 - t272 * t337 + (-t274 * t338 + t276 * t335) * t402) * MDP(12) + (-t336 * t376 - t273 * t337 + (-t277 * t338 + t279 * t335) * t402) * MDP(13) + (t274 * t279 - t276 * t277 + (-t272 * t336 + t273 * t333) * t334) * MDP(14) + (t223 * t276 - t224 * t279 - t234 * t337 + (-t216 * t333 + t217 * t336 - t242 * t399) * t334 - t370) * MDP(15) + (-t288 * t267 + t226 * t418 - (-t279 * t339 + t338 * t380) * t268 - t276 * t244 - qJD(5) * t347) * MDP(21) + (-t289 * t267 + t225 * t418 + (t279 * t341 + t338 * t381) * t268 - t276 * t246 + t288 * t396) * MDP(22) + (MDP(5) + t366) * t388; (t327 + t294) * MDP(11) + (-t272 * t333 - t273 * t336 + (-t274 * t336 + t277 * t333) * t382) * MDP(14) + (t216 * t336 + t217 * t333 + t327) * MDP(15) + (-t336 * t226 + (-t341 * t396 - t406) * t333 + (t244 * t418 - t268 * t288) * t401) * MDP(21) + (-t336 * t225 + (t268 * t395 - t404) * t333 + (t246 * t418 - t347) * t401) * MDP(22) + ((-t337 * t407 + t389) * MDP(8) + (t334 * t407 + t387) * MDP(9) + ((t250 * t337 + t251 * t334) * qJD(1) - t362) * MDP(11) + (-t274 * t400 + t336 * t389) * MDP(12) + (-t277 * t400 - t333 * t389) * MDP(13) + ((-t242 * t337 + (-t223 * t333 + t224 * t336) * t334) * qJD(1) - t362) * MDP(15)) * t335 + (t366 + (-MDP(12) * t333 - MDP(13) * t336) * t329) * t330 * t343; t336 * MDP(12) * t386 - t311 * MDP(13) + (-t274 ^ 2 - t277 ^ 2) * MDP(14) + (-g(1) * t292 - g(2) * t290 + t223 * t277 + t224 * t274 + t234) * MDP(15) + (-t244 * t277 + t404) * MDP(21) + (-t246 * t277 - t406) * MDP(22) + ((t277 * t402 + t333 * t387) * MDP(12) + (-t274 * t402 + t373) * MDP(13) - g(3) * t334 * MDP(15)) * t335 - (MDP(21) * t339 + MDP(22) * t341) * t268 ^ 2; t246 * t244 * MDP(16) + (-t244 ^ 2 + t246 ^ 2) * MDP(17) + (t244 * t268 + t384) * MDP(18) + (t246 * t268 - t367) * MDP(19) + t394 + (-t357 * t268 - t221 * t246 + g(1) * t424 - g(2) * (-(t285 * t340 - t334 * t411) * t339 + t290 * t341) - g(3) * t352 + t368) * MDP(21) + (t358 * t268 + t221 * t244 - g(1) * (-t289 * t340 - t342 * t351) + g(2) * t426 + g(3) * t253 - t359) * MDP(22) + (-MDP(18) * t420 - MDP(19) * t246 + MDP(21) * t357 - MDP(22) * t358) * qJD(5);];
tau = t1;
