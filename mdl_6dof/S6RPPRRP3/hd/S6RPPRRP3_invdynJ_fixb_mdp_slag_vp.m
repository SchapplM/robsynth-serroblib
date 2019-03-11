% Calculate vector of inverse dynamics joint torques for
% S6RPPRRP3
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% qJDD [6x1]
%   Generalized joint accelerations
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d4,d5,theta2]';
% MDP [25x1]
%   Minimal dynamic parameter vector (fixed base model)
%   see S6RPPRRP3_convert_par2_MPV_fixb.m
% 
% Output:
% tau [6x1]
%   joint torques of inverse dynamics (contains inertial, gravitational Coriolis and centrifugal forces)

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 02:04
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tau = S6RPPRRP3_invdynJ_fixb_mdp_slag_vp(qJ, qJD, qJDD, g, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(3,1),zeros(9,1),zeros(25,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPPRRP3_invdynJ_fixb_mdp_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RPPRRP3_invdynJ_fixb_mdp_slag_vp: qJD has to be [6x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [6 1]), ...
  'S6RPPRRP3_invdynJ_fixb_mdp_slag_vp: qJDD has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RPPRRP3_invdynJ_fixb_mdp_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S6RPPRRP3_invdynJ_fixb_mdp_slag_vp: pkin has to be [9x1] (double)');
assert(isreal(MDP) && all(size(MDP) == [25 1]), ...
  'S6RPPRRP3_invdynJ_fixb_mdp_slag_vp: MDP has to be [25x1] (double)'); 

%% Symbolic Calculation
% From invdyn_joint_fixb_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 02:03:58
% EndTime: 2019-03-09 02:04:02
% DurationCPUTime: 3.76s
% Computational Cost: add. (2867->418), mult. (5017->510), div. (0->0), fcn. (2936->10), ass. (0->169)
t341 = sin(pkin(9));
t320 = pkin(1) * t341 + qJ(3);
t424 = qJDD(1) * t320;
t347 = cos(qJ(4));
t336 = qJ(1) + pkin(9);
t328 = cos(t336);
t322 = g(2) * t328;
t327 = sin(t336);
t450 = g(1) * t327;
t387 = t322 - t450;
t344 = sin(qJ(4));
t449 = g(3) * t344;
t464 = t347 * t387 + t449;
t342 = cos(pkin(9));
t323 = -pkin(1) * t342 - pkin(2);
t314 = -pkin(7) + t323;
t296 = qJD(1) * t314 + qJD(3);
t274 = -t344 * qJD(2) + t296 * t347;
t463 = qJD(4) * t274;
t315 = qJD(1) * t344 + qJD(5);
t346 = cos(qJ(5));
t400 = MDP(21) - MDP(24);
t462 = t346 * t400;
t382 = pkin(4) * t344 - pkin(8) * t347;
t292 = t320 + t382;
t343 = sin(qJ(5));
t434 = t344 * t346;
t426 = t343 * t292 + t314 * t434;
t417 = qJD(4) * t344;
t392 = t343 * t417;
t403 = qJDD(1) * t347;
t418 = qJD(4) * t343;
t420 = qJD(1) * t347;
t303 = t346 * t420 + t418;
t414 = qJD(5) * t303;
t263 = -qJD(1) * t392 - t346 * qJDD(4) + t343 * t403 + t414;
t461 = qJDD(1) * t323;
t271 = -qJD(4) * pkin(4) - t274;
t409 = t346 * qJD(4);
t301 = t343 * t420 - t409;
t247 = pkin(5) * t301 - qJ(6) * t303 + t271;
t407 = qJD(1) * qJD(4);
t390 = t347 * t407;
t404 = qJDD(1) * t344;
t298 = qJDD(5) + t390 + t404;
t452 = pkin(8) * t298;
t460 = t247 * t315 - t452;
t447 = pkin(8) * qJD(5);
t398 = t315 * t447;
t459 = -t398 + t464;
t411 = qJD(5) * t347;
t363 = t343 * t411 + t344 * t409;
t262 = qJD(1) * t363 - qJD(5) * t409 - t343 * qJDD(4) - t346 * t403;
t310 = qJD(1) * t320;
t458 = 0.2e1 * qJD(4) * t310 + qJDD(4) * t314;
t275 = t347 * qJD(2) + t344 * t296;
t272 = qJD(4) * pkin(8) + t275;
t276 = t292 * qJD(1);
t248 = -t272 * t343 + t276 * t346;
t408 = qJD(6) - t248;
t242 = -pkin(5) * t315 + t408;
t249 = t272 * t346 + t276 * t343;
t243 = qJ(6) * t315 + t249;
t371 = t242 * t343 + t243 * t346;
t457 = (t301 * t346 - t303 * t343) * MDP(23) - t371 * MDP(25);
t456 = t303 ^ 2;
t453 = pkin(5) * t298;
t448 = g(3) * t347;
t446 = qJ(6) * t298;
t445 = t243 * t315;
t444 = t249 * t315;
t443 = t262 * t343;
t383 = pkin(4) * t347 + pkin(8) * t344;
t299 = qJD(4) * t383 + qJD(3);
t441 = t299 * t346;
t440 = t301 * t303;
t439 = t301 * t315;
t438 = t303 * t315;
t437 = t303 * t346;
t436 = t314 * t343;
t435 = t343 * t344;
t258 = t344 * t263;
t433 = t346 * t298;
t432 = t346 * t347;
t431 = qJDD(2) - g(3);
t416 = qJD(4) * t347;
t430 = t301 * t416 + t258;
t429 = -t344 * t262 + t303 * t416;
t305 = t383 * qJD(1);
t428 = t346 * t274 + t343 * t305;
t374 = pkin(5) * t343 - qJ(6) * t346;
t427 = -qJD(6) * t343 + t315 * t374 - t275;
t425 = g(3) * t434 + t432 * t322;
t339 = t347 ^ 2;
t423 = t344 ^ 2 - t339;
t349 = qJD(4) ^ 2;
t350 = qJD(1) ^ 2;
t422 = -t349 - t350;
t419 = qJD(4) * t301;
t415 = qJD(5) * t301;
t413 = qJD(5) * t343;
t412 = qJD(5) * t346;
t293 = qJDD(1) * t314 + qJDD(3);
t406 = qJD(2) * qJD(4);
t395 = -t344 * qJDD(2) - t296 * t417 - t347 * t406;
t252 = -qJDD(4) * pkin(4) - t293 * t347 - t395;
t239 = pkin(5) * t263 + qJ(6) * t262 - qJD(6) * t303 + t252;
t410 = t239 * MDP(25);
t405 = qJD(3) * qJD(1);
t401 = MDP(20) + MDP(22);
t399 = t347 * t450;
t251 = qJDD(4) * pkin(8) + qJDD(2) * t347 + t293 * t344 + t463;
t261 = qJD(1) * t299 + qJDD(1) * t382 + t424;
t397 = -t346 * t251 - t343 * t261 - t276 * t412;
t396 = t347 * t314 * t409 + t292 * t412 + t343 * t299;
t348 = cos(qJ(1));
t393 = t348 * pkin(1) + t328 * pkin(2) + t327 * qJ(3);
t391 = t301 * t417;
t345 = sin(qJ(1));
t389 = -pkin(1) * t345 + t328 * qJ(3);
t388 = -pkin(5) + t436;
t385 = t401 * t343;
t384 = -t343 * t251 + t346 * t261 - t272 * t412 - t276 * t413;
t277 = t327 * t435 - t328 * t346;
t279 = t327 * t346 + t328 * t435;
t381 = g(1) * t279 + g(2) * t277;
t278 = t327 * t434 + t328 * t343;
t280 = -t327 * t343 + t328 * t434;
t380 = -g(1) * t280 - g(2) * t278;
t379 = g(1) * t328 + g(2) * t327;
t377 = g(1) * t345 - g(2) * t348;
t376 = qJDD(3) + t387;
t375 = pkin(5) * t346 + qJ(6) * t343;
t364 = -t272 * t413 - t397;
t237 = qJD(6) * t315 + t364 + t446;
t238 = qJDD(6) - t384 - t453;
t373 = t237 * t346 + t238 * t343;
t372 = t242 * t346 - t243 * t343;
t369 = pkin(4) + t375;
t367 = -t314 + t374;
t365 = t298 * t343 + t315 * t412;
t362 = -t379 + t424;
t361 = t271 * t315 - t452;
t360 = -qJD(1) * t310 + t293 + t387;
t359 = -t385 - t462;
t358 = t344 * t387 - t448;
t357 = g(1) * t277 - g(2) * t279 + t343 * t448 + t384;
t356 = t247 * t303 + qJDD(6) - t357;
t297 = t405 + t424;
t355 = -t314 * t349 + t297 + t362 + t405;
t354 = -g(1) * t278 + g(2) * t280 - g(3) * t432 + t364;
t353 = (t301 * t343 + t437) * MDP(23) + t372 * MDP(25) + (t343 * t400 - t346 * t401) * t315;
t352 = t353 * qJD(5);
t309 = qJDD(4) * t347 - t344 * t349;
t308 = -qJDD(4) * t344 - t347 * t349;
t295 = t315 * t392;
t282 = t298 * t432;
t273 = t367 * t347;
t267 = pkin(5) * t303 + qJ(6) * t301;
t265 = -t292 * t346 + t344 * t388;
t264 = qJ(6) * t344 + t426;
t255 = t263 * t432;
t254 = -pkin(5) * t420 + t274 * t343 - t305 * t346;
t253 = qJ(6) * t420 + t428;
t250 = (qJD(5) * t375 - qJD(6) * t346) * t347 - t367 * t417;
t244 = -t262 + t439;
t241 = qJD(5) * t426 + t388 * t416 - t441;
t240 = qJ(6) * t416 + (-t314 * t413 + qJD(6)) * t344 + t396;
t1 = [(-t255 + (-t303 * t411 + t391) * t346 + (t303 * t417 + (t262 + t415) * t347) * t343) * MDP(16) + (-t315 * t363 + t282 + t429) * MDP(17) + (-t258 + t295 + (-t365 - t419) * t347) * MDP(18) + (t298 * t344 + t315 * t416) * MDP(19) + ((-t292 * t413 + t441) * t315 + t292 * t433 + (-t271 * t418 + (-t365 + t419) * t314 + t384) * t344 + (t271 * t412 + t252 * t343 - t314 * t263 + (-t315 * t436 + t248) * qJD(4)) * t347 + t380) * MDP(20) + (-t396 * t315 - t426 * t298 + ((t314 * t315 + t272) * t413 + (-t271 * t346 + t303 * t314) * qJD(4) + t397) * t344 + (-qJD(4) * t249 + t252 * t346 + t262 * t314 - t271 * t413) * t347 + t381) * MDP(21) + (-t241 * t315 + t250 * t301 + t263 * t273 - t265 * t298 + (-t247 * t418 - t238) * t344 + (-qJD(4) * t242 + t239 * t343 + t247 * t412) * t347 + t380) * MDP(22) + (-t240 * t301 + t241 * t303 - t262 * t265 - t263 * t264 - t372 * t417 + (-qJD(5) * t371 - t237 * t343 + t238 * t346 + t379) * t347) * MDP(23) + (t240 * t315 - t250 * t303 + t262 * t273 + t264 * t298 + (t247 * t409 + t237) * t344 + (qJD(4) * t243 - t239 * t346 + t247 * t413) * t347 - t381) * MDP(24) + (t237 * t264 + t243 * t240 + t239 * t273 + t247 * t250 + t238 * t265 + t242 * t241 - g(1) * (pkin(5) * t280 + qJ(6) * t279 + t382 * t328 + t389) - g(2) * (pkin(5) * t278 + pkin(7) * t328 + qJ(6) * t277 + t393) + (-g(1) * (-pkin(2) - pkin(7)) - g(2) * t382) * t327) * MDP(25) + qJDD(1) * MDP(1) + (t377 + (t341 ^ 2 + t342 ^ 2) * qJDD(1) * pkin(1)) * pkin(1) * MDP(4) + (t376 + 0.2e1 * t461) * MDP(5) + (t362 + 0.2e1 * t405 + t424) * MDP(6) + (t297 * t320 + t310 * qJD(3) + (qJDD(3) + t461) * t323 - g(1) * (-pkin(2) * t327 + t389) - g(2) * t393) * MDP(7) + (qJDD(1) * t339 - 0.2e1 * t344 * t390) * MDP(8) + 0.2e1 * (-t344 * t403 + t407 * t423) * MDP(9) + t309 * MDP(10) + t308 * MDP(11) + (t355 * t344 + t347 * t458) * MDP(13) + (-t344 * t458 + t355 * t347) * MDP(14) + (-t262 * t432 - t303 * t363) * MDP(15) + t377 * MDP(2) + (g(1) * t348 + g(2) * t345) * MDP(3); t308 * MDP(13) - t309 * MDP(14) + (t295 + t430) * MDP(20) + t429 * MDP(21) + t430 * MDP(22) - t255 * MDP(23) + t282 * MDP(24) - g(3) * MDP(25) + (MDP(4) + MDP(7)) * t431 + (t262 * MDP(24) + t410 + ((t343 * MDP(22) + t462) * t315 + t457) * qJD(4)) * t344 + (-MDP(23) * t443 - qJD(4) * t303 * MDP(24) + (qJD(4) * t247 + t373) * MDP(25) + (-t346 * MDP(21) - t385) * t298 + t352) * t347; -t350 * MDP(6) + t376 * MDP(7) + t387 * MDP(25) + t401 * t391 + (MDP(7) * t323 + MDP(5)) * qJDD(1) + (-t310 * MDP(7) + t353) * qJD(1) + (qJDD(4) * MDP(13) + t422 * MDP(14) - t410 - t401 * t263 + t400 * t262 + (t315 * t359 - t457) * qJD(4)) * t347 + (t422 * MDP(13) - qJDD(4) * MDP(14) + (-t263 * t346 - t443) * MDP(23) + t373 * MDP(25) + (t247 * MDP(25) + t303 * t400) * qJD(4) + t359 * t298 + t352) * t344; MDP(10) * t403 - MDP(11) * t404 + qJDD(4) * MDP(12) + (qJD(4) * t275 + t347 * t360 + t395 + t449) * MDP(13) + (t463 + (-qJD(4) * t296 - t431) * t347 + (-t360 + t406) * t344) * MDP(14) + (t315 * t437 - t443) * MDP(15) + ((-t262 - t439) * t346 + (-t263 - t438) * t343) * MDP(16) + ((-t303 * t347 + t315 * t434) * qJD(1) + t365) * MDP(17) + (-t315 * t413 + t433 + (t301 * t347 - t315 * t435) * qJD(1)) * MDP(18) - t315 * MDP(19) * t420 + (-t248 * t420 - pkin(4) * t263 - t275 * t301 + (-t399 - t252 + (-t305 - t447) * t315) * t346 + (t274 * t315 + t361) * t343 + t425) * MDP(20) + (pkin(4) * t262 + t428 * t315 + t249 * t420 - t275 * t303 + t361 * t346 + (t252 - t459) * t343) * MDP(21) + (t242 * t420 + t254 * t315 - t263 * t369 + t427 * t301 + (-t239 - t398 - t399) * t346 + t460 * t343 + t425) * MDP(22) + (t253 * t301 - t254 * t303 + (t237 + t315 * t242 + (-t263 + t414) * pkin(8)) * t346 + (t238 - t445 + (-t262 + t415) * pkin(8)) * t343 + t358) * MDP(23) + (-t243 * t420 - t253 * t315 - t262 * t369 - t427 * t303 - t460 * t346 + (-t239 + t459) * t343) * MDP(24) + (-t242 * t254 - t243 * t253 + t427 * t247 + (qJD(5) * t372 + t358 + t373) * pkin(8) + (-t239 + t464) * t369) * MDP(25) + (MDP(8) * t344 * t347 - MDP(9) * t423) * t350; MDP(15) * t440 + (-t301 ^ 2 + t456) * MDP(16) + t244 * MDP(17) + (t438 - t263) * MDP(18) + t298 * MDP(19) + (-t271 * t303 + t357 + t444) * MDP(20) + (t248 * t315 + t271 * t301 - t354) * MDP(21) + (-t267 * t301 - t356 + t444 + 0.2e1 * t453) * MDP(22) + (pkin(5) * t262 - qJ(6) * t263 + (t243 - t249) * t303 + (t242 - t408) * t301) * MDP(23) + (0.2e1 * t446 - t247 * t301 + t267 * t303 + (0.2e1 * qJD(6) - t248) * t315 + t354) * MDP(24) + (t237 * qJ(6) - t238 * pkin(5) - t247 * t267 - t242 * t249 - g(1) * (-pkin(5) * t277 + qJ(6) * t278) - g(2) * (pkin(5) * t279 - qJ(6) * t280) + t374 * t448 + t408 * t243) * MDP(25); (-t298 + t440) * MDP(22) + t244 * MDP(23) + (-t315 ^ 2 - t456) * MDP(24) + (t356 - t445 - t453) * MDP(25);];
tau  = t1;
