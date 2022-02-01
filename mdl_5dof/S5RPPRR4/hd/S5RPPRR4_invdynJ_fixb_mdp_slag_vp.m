% Calculate vector of inverse dynamics joint torques for
% S5RPPRR4
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
%   pkin=[a2,a3,a4,a5,d1,d4,d5,theta2,theta3]';
% MDP [23x1]
%   Minimal dynamic parameter vector (fixed base model)
%   see S5RPPRR4_convert_par2_MPV_fixb.m
% 
% Output:
% tau [5x1]
%   joint torques of inverse dynamics (contains inertial, gravitational Coriolis and centrifugal forces)

% Quelle: HybrDyn-Toolbox
% Datum: 2022-01-23 09:17
% Revision: fd3771346c4aea32fdeb66112c511235427c26a7 (2022-01-20)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tau = S5RPPRR4_invdynJ_fixb_mdp_slag_vp(qJ, qJD, qJDD, g, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(5,1),zeros(3,1),zeros(9,1),zeros(23,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPPRR4_invdynJ_fixb_mdp_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RPPRR4_invdynJ_fixb_mdp_slag_vp: qJD has to be [5x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [5 1]), ...
  'S5RPPRR4_invdynJ_fixb_mdp_slag_vp: qJDD has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RPPRR4_invdynJ_fixb_mdp_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5RPPRR4_invdynJ_fixb_mdp_slag_vp: pkin has to be [9x1] (double)');
assert(isreal(MDP) && all(size(MDP) == [23 1]), ...
  'S5RPPRR4_invdynJ_fixb_mdp_slag_vp: MDP has to be [23x1] (double)'); 

%% Symbolic Calculation
% From invdyn_joint_fixb_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2022-01-23 09:17:03
% EndTime: 2022-01-23 09:17:09
% DurationCPUTime: 3.93s
% Computational Cost: add. (2055->330), mult. (5144->448), div. (0->0), fcn. (4002->14), ass. (0->165)
t362 = sin(pkin(9));
t364 = cos(pkin(9));
t363 = sin(pkin(8));
t424 = qJD(3) * t363;
t365 = cos(pkin(8));
t425 = qJD(2) * t365;
t316 = -t362 * t425 - t364 * t424;
t317 = -t362 * t424 + t364 * t425;
t367 = sin(qJ(4));
t370 = cos(qJ(4));
t328 = pkin(2) * t365 + qJ(3) * t363 + pkin(1);
t322 = t364 * t328;
t377 = -pkin(6) * t363 * t364 + (-qJ(2) * t362 - pkin(3)) * t365;
t274 = -t322 + t377;
t290 = t364 * t365 * qJ(2) - t362 * t328;
t438 = t362 * t363;
t279 = -pkin(6) * t438 + t290;
t432 = t367 * t274 + t370 * t279;
t460 = qJD(4) * t432 - t370 * t316 + t317 * t367;
t324 = t362 * t370 + t364 * t367;
t378 = qJD(1) * t324;
t296 = t363 * t378;
t369 = cos(qJ(5));
t427 = qJD(1) * t363;
t409 = t362 * t427;
t394 = t367 * t409;
t435 = t364 * t370;
t412 = t363 * t435;
t299 = qJD(1) * t412 - t394;
t366 = sin(qJ(5));
t439 = t299 * t366;
t253 = t369 * t296 + t439;
t426 = qJD(1) * t365;
t456 = qJD(4) - t426;
t340 = -qJD(5) - t456;
t442 = t253 * t340;
t319 = t324 * qJD(4);
t329 = qJDD(1) * t412;
t436 = t362 * t367;
t264 = t329 + (-qJD(1) * t319 - qJDD(1) * t436) * t363;
t413 = qJDD(1) * t365;
t343 = -qJDD(4) + t413;
t287 = -qJD(1) * t424 - qJDD(1) * t328 + qJDD(2);
t282 = t364 * t287;
t415 = qJD(1) * qJD(2);
t407 = t365 * t415;
t249 = qJDD(1) * t377 - t362 * t407 + t282;
t405 = t364 * t413;
t263 = qJ(2) * t405 + t362 * t287 + t364 * t407;
t414 = qJDD(1) * t363;
t406 = t362 * t414;
t251 = -pkin(6) * t406 + t263;
t314 = -qJD(1) * t328 + qJD(2);
t304 = t364 * t314;
t266 = qJD(1) * t377 + t304;
t410 = qJ(2) * t426;
t278 = t362 * t314 + t364 * t410;
t271 = -pkin(6) * t409 + t278;
t386 = -t266 * t367 - t271 * t370;
t375 = qJD(4) * t386 + t370 * t249 - t367 * t251;
t232 = -pkin(4) * t343 - pkin(7) * t264 + t375;
t459 = t366 * t232;
t368 = sin(qJ(1));
t371 = cos(qJ(1));
t452 = g(1) * t368 - g(2) * t371;
t399 = t370 * t266 - t271 * t367;
t240 = -pkin(7) * t299 + t399;
t238 = pkin(4) * t456 + t240;
t241 = -pkin(7) * t296 - t386;
t444 = t241 * t369;
t387 = -t238 * t366 - t444;
t457 = t387 * qJD(5);
t421 = qJD(5) * t366;
t239 = t241 * t421;
t341 = qJ(2) * t427 + qJD(3);
t315 = pkin(3) * t409 + t341;
t269 = pkin(4) * t296 + t315;
t361 = pkin(9) + qJ(4);
t355 = qJ(5) + t361;
t349 = sin(t355);
t350 = cos(t355);
t434 = t365 * t368;
t292 = t349 * t371 - t350 * t434;
t433 = t365 * t371;
t294 = t349 * t368 + t350 * t433;
t447 = g(3) * t363;
t455 = g(1) * t294 - g(2) * t292 + t269 * t253 + t350 * t447 + t239;
t326 = qJD(4) * t394;
t422 = qJD(4) * t370;
t408 = t364 * t422;
t265 = t363 * (qJD(1) * t408 + qJDD(1) * t324) - t326;
t337 = -qJDD(5) + t343;
t327 = t337 * MDP(21);
t384 = -t296 * t366 + t369 * t299;
t454 = t253 * t384 * MDP(17) + (-t253 ^ 2 + t384 ^ 2) * MDP(18) - t327;
t443 = t384 * t340;
t446 = qJDD(1) * pkin(1);
t451 = t446 + t452;
t416 = qJ(2) * qJDD(1);
t291 = t349 * t434 + t350 * t371;
t293 = -t349 * t433 + t350 * t368;
t423 = qJD(4) * t367;
t380 = -t367 * t249 - t370 * t251 - t266 * t422 + t271 * t423;
t233 = -pkin(7) * t265 - t380;
t402 = t369 * t232 - t366 * t233;
t450 = -g(1) * t293 + g(2) * t291 - t269 * t384 + t349 * t447 + t402;
t400 = t264 * t366 + t369 * t265;
t235 = qJD(5) * t384 + t400;
t445 = t238 * t369;
t441 = t296 * t456;
t440 = t299 * t456;
t437 = t362 * t365;
t431 = -t365 * t378 + t319;
t323 = t435 - t436;
t430 = t456 * t323;
t325 = pkin(3) * t438 + t363 * qJ(2);
t359 = t363 ^ 2;
t429 = t365 ^ 2 + t359;
t428 = MDP(15) * t370;
t420 = qJD(5) * t369;
t418 = t343 * MDP(14);
t417 = t363 * qJD(2);
t411 = t369 * t264 - t366 * t265 - t296 * t420;
t320 = qJ(2) * t414 + t363 * t415 + qJDD(3);
t372 = qJD(1) ^ 2;
t403 = t429 * t372;
t398 = t370 * t274 - t279 * t367;
t396 = qJD(5) * t238 + t233;
t395 = 0.2e1 * t429;
t288 = pkin(3) * t406 + t320;
t392 = g(1) * t371 + g(2) * t368;
t390 = qJD(5) * t324 + t431;
t389 = qJD(5) * t323 + t430;
t312 = t324 * t363;
t313 = t323 * t363;
t267 = t369 * t312 + t313 * t366;
t268 = -t312 * t366 + t313 * t369;
t383 = t415 + t416;
t382 = t395 * t415;
t379 = t274 * t422 - t279 * t423 + t367 * t316 + t370 * t317;
t234 = -t299 * t421 + t411;
t357 = t371 * qJ(2);
t356 = t368 * qJ(2);
t353 = cos(t361);
t352 = sin(t361);
t351 = qJDD(2) - t446;
t308 = t352 * t368 + t353 * t433;
t307 = -t352 * t433 + t353 * t368;
t306 = t352 * t371 - t353 * t434;
t305 = t352 * t434 + t353 * t371;
t302 = t363 * t408 - t423 * t438;
t301 = t363 * t319;
t289 = -qJ(2) * t437 - t322;
t280 = pkin(4) * t302 + t417;
t277 = -t362 * t410 + t304;
t276 = pkin(4) * t312 + t325;
t262 = -t383 * t437 + t282;
t246 = pkin(4) * t265 + t288;
t245 = -pkin(7) * t312 + t432;
t244 = -pkin(4) * t365 - pkin(7) * t313 + t398;
t243 = qJD(5) * t268 - t301 * t366 + t369 * t302;
t242 = -qJD(5) * t267 - t301 * t369 - t302 * t366;
t237 = pkin(7) * t301 - t460;
t236 = -pkin(7) * t302 + t379;
t1 = [qJDD(1) * MDP(1) + t452 * MDP(2) + t392 * MDP(3) + (t395 * t416 + t382 - t392) * MDP(5) + (-t351 * pkin(1) - g(1) * (-pkin(1) * t368 + t357) - g(2) * (pkin(1) * t371 + t356) + (t429 * t416 + t382) * qJ(2)) * MDP(6) + (t263 * t290 + t278 * t317 + t262 * t289 + t277 * t316 - g(1) * (-t328 * t368 + t357) - g(2) * (t328 * t371 + t356) + (t320 * qJ(2) + t341 * qJD(2)) * t363) * MDP(9) + (t264 * t313 - t299 * t301) * MDP(10) + (-t264 * t312 - t265 * t313 + t296 * t301 - t299 * t302) * MDP(11) + (-t301 * t456 - t313 * t343) * MDP(12) + (-t302 * t456 + t312 * t343) * MDP(13) + (-g(1) * t306 - g(2) * t308 + t325 * t265 + t288 * t312 + t296 * t417 + t315 * t302 - t398 * t343 - t456 * t460) * MDP(15) + (-g(1) * t305 - g(2) * t307 + t325 * t264 + t288 * t313 + t299 * t417 - t315 * t301 + t343 * t432 - t379 * t456) * MDP(16) + (t234 * t268 + t242 * t384) * MDP(17) + (-t234 * t267 - t235 * t268 - t242 * t253 - t243 * t384) * MDP(18) + (-t242 * t340 - t268 * t337) * MDP(19) + (t243 * t340 + t267 * t337) * MDP(20) + (-(t244 * t369 - t245 * t366) * t337 + t280 * t253 + t276 * t235 + t246 * t267 + t269 * t243 - g(1) * t292 - g(2) * t294 + (t236 * t366 - t237 * t369 - (-t244 * t366 - t245 * t369) * qJD(5)) * t340) * MDP(22) + (-g(1) * t291 - g(2) * t293 + t276 * t234 + t269 * t242 + t246 * t268 + t280 * t384 + ((-qJD(5) * t245 + t237) * t340 + t244 * t337) * t366 + ((qJD(5) * t244 + t236) * t340 + t245 * t337) * t369) * MDP(23) + (t362 * MDP(7) + t364 * MDP(8)) * (t320 * t363 + t359 * t383 - t392) + ((-t351 + t451) * MDP(4) + (-t316 * qJD(1) - t289 * qJDD(1) + t364 * t452 - t262) * MDP(7) + (t317 * qJD(1) + t290 * qJDD(1) - t362 * t452 + t263) * MDP(8) - t264 * MDP(12) + t265 * MDP(13) + t418 - t375 * MDP(15) - t380 * MDP(16) - t234 * MDP(19) + t235 * MDP(20) + t327 + (-t402 - t457) * MDP(22) + (t396 * t369 - t239 + t459) * MDP(23)) * t365; -MDP(4) * t413 - MDP(5) * t403 + (-qJ(2) * t403 + qJDD(2) - t451) * MDP(6) + (-t362 * t403 - t405) * MDP(7) + (t362 * t413 - t364 * t403) * MDP(8) + (t262 * t364 + t263 * t362 + (-t341 * t363 + (t277 * t362 - t278 * t364) * t365) * qJD(1) - t452) * MDP(9) + (-t296 * t427 - t323 * t343 - t431 * t456) * MDP(15) + (-t299 * t427 + t324 * t343 - t430 * t456) * MDP(16) + (-(t323 * t369 - t324 * t366) * t337 - t253 * t427 + (t366 * t389 + t369 * t390) * t340) * MDP(22) + ((t323 * t366 + t324 * t369) * t337 - t384 * t427 + (-t366 * t390 + t369 * t389) * t340) * MDP(23); (g(3) * t365 + t320) * MDP(9) + (-t326 + t440) * MDP(15) + (t329 - t441) * MDP(16) + (t235 - t443) * MDP(22) + (t234 + t442) * MDP(23) + (-t392 * MDP(9) + (-MDP(7) * t364 + MDP(8) * t362) * t372 * t365 + ((MDP(15) * t367 + MDP(8)) * t364 + (-MDP(16) * t367 + MDP(7) + t428) * t362) * qJDD(1) + ((t277 * t364 + t278 * t362) * MDP(9) + (-MDP(16) * t324 + t364 * t428) * qJD(4)) * qJD(1)) * t363; t299 * t296 * MDP(10) + (-t296 ^ 2 + t299 ^ 2) * MDP(11) + (t264 + t441) * MDP(12) + (-t265 + t440) * MDP(13) - t418 + (-g(1) * t307 + g(2) * t305 - t315 * t299 + t352 * t447 - t386 * t456 + t375) * MDP(15) + (g(1) * t308 - g(2) * t306 + t315 * t296 + t353 * t447 + t399 * t456 + t380) * MDP(16) + (t234 - t442) * MDP(19) + (-t235 - t443) * MDP(20) + ((-t240 * t366 - t444) * t340 + t457 + (-t253 * t299 - t369 * t337 + t340 * t421) * pkin(4) + t450) * MDP(22) + ((t241 * t340 - t232) * t366 + (-t240 * t340 - t396) * t369 + (-t299 * t384 + t366 * t337 + t340 * t420) * pkin(4) + t455) * MDP(23) + t454; (t411 - t442) * MDP(19) + (-t400 - t443) * MDP(20) + (t340 * t387 + t450) * MDP(22) + (-t369 * t233 - t459 - (-t241 * t366 + t445) * t340 + t455) * MDP(23) + (-MDP(19) * t439 - MDP(20) * t384 + MDP(22) * t387 - MDP(23) * t445) * qJD(5) + t454;];
tau = t1;
