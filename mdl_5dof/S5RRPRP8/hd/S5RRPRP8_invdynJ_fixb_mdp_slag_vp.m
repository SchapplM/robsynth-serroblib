% Calculate vector of inverse dynamics joint torques for
% S5RRPRP8
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
%   pkin=[a2,a3,a4,a5,d1,d2,d4]';
% MDP [23x1]
%   Minimal dynamic parameter vector (fixed base model)
%   see S5RRPRP8_convert_par2_MPV_fixb.m
% 
% Output:
% tau [5x1]
%   joint torques of inverse dynamics (contains inertial, gravitational Coriolis and centrifugal forces)

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 20:04
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tau = S5RRPRP8_invdynJ_fixb_mdp_slag_vp(qJ, qJD, qJDD, g, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(5,1),zeros(3,1),zeros(7,1),zeros(23,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRPRP8_invdynJ_fixb_mdp_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RRPRP8_invdynJ_fixb_mdp_slag_vp: qJD has to be [5x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [5 1]), ...
  'S5RRPRP8_invdynJ_fixb_mdp_slag_vp: qJDD has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RRPRP8_invdynJ_fixb_mdp_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [7 1]), ...
  'S5RRPRP8_invdynJ_fixb_mdp_slag_vp: pkin has to be [7x1] (double)');
assert(isreal(MDP) && all(size(MDP) == [23 1]), ...
  'S5RRPRP8_invdynJ_fixb_mdp_slag_vp: MDP has to be [23x1] (double)'); 

%% Symbolic Calculation
% From invdyn_joint_fixb_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 20:04:33
% EndTime: 2019-12-31 20:04:38
% DurationCPUTime: 3.11s
% Computational Cost: add. (1454->322), mult. (3108->399), div. (0->0), fcn. (1900->6), ass. (0->156)
t357 = sin(qJ(4));
t358 = sin(qJ(2));
t361 = cos(qJ(2));
t360 = cos(qJ(4));
t412 = qJD(4) * t360;
t413 = qJD(4) * t357;
t414 = qJD(2) * t361;
t465 = t357 * t414 + t358 * t412 - t361 * t413;
t363 = -pkin(2) - pkin(3);
t410 = qJDD(1) * t358;
t327 = pkin(6) * t410;
t411 = qJD(1) * qJD(2);
t401 = t361 * t411;
t400 = pkin(6) * t401 + qJDD(3) + t327;
t459 = t401 + t410;
t251 = -pkin(7) * t459 + t363 * qJDD(2) + t400;
t409 = qJDD(1) * t361;
t328 = pkin(6) * t409;
t349 = qJDD(2) * qJ(3);
t350 = qJD(2) * qJD(3);
t402 = t358 * t411;
t259 = -pkin(6) * t402 + t328 + t349 + t350;
t458 = t402 - t409;
t252 = pkin(7) * t458 + t259;
t404 = t363 * qJD(2);
t417 = qJD(1) * t358;
t332 = pkin(6) * t417;
t463 = -pkin(7) * t417 + qJD(3) + t332;
t260 = t404 + t463;
t416 = qJD(1) * t361;
t333 = pkin(6) * t416;
t294 = -pkin(7) * t416 + t333;
t351 = qJD(2) * qJ(3);
t275 = t294 + t351;
t384 = t260 * t357 + t275 * t360;
t464 = -t384 * qJD(4) + t360 * t251 - t357 * t252;
t434 = t357 * t358;
t285 = t360 * t361 + t434;
t271 = t285 * qJD(1);
t270 = t271 ^ 2;
t273 = -t357 * t416 + t360 * t417;
t347 = qJDD(2) - qJDD(4);
t448 = t273 ^ 2;
t240 = t465 * qJD(1) + qJDD(1) * t285 - t360 * t402;
t348 = qJD(2) - qJD(4);
t455 = t273 * t348 + t240;
t375 = t285 * qJD(4);
t239 = qJD(1) * t375 - t357 * t458 - t360 * t459;
t461 = t271 * t348 + t239;
t462 = -t461 * MDP(17) + t273 * t271 * MDP(15) - (-t448 + t270) * MDP(16) - t455 * MDP(18) - t347 * MDP(19);
t276 = -qJD(1) * pkin(1) - pkin(2) * t416 - qJ(3) * t417;
t258 = pkin(3) * t416 - t276;
t359 = sin(qJ(1));
t430 = t361 * t357;
t382 = -t358 * t360 + t430;
t265 = t382 * t359;
t362 = cos(qJ(1));
t429 = t361 * t362;
t407 = t357 * t429;
t432 = t358 * t362;
t267 = -t360 * t432 + t407;
t444 = g(3) * t285;
t460 = g(1) * t267 + g(2) * t265 - t258 * t273 + t444 + t464;
t420 = t361 * pkin(2) + t358 * qJ(3);
t456 = -pkin(1) - t420;
t397 = -qJ(3) * t357 + t360 * t363;
t454 = t397 * qJD(4) - t357 * t294 + t463 * t360;
t296 = qJ(3) * t360 + t357 * t363;
t453 = qJD(4) * t296 + t360 * t294 + t463 * t357;
t447 = pkin(6) - pkin(7);
t300 = t447 * t358;
t301 = t447 * t361;
t423 = t357 * t300 + t360 * t301;
t452 = g(1) * t362 + g(2) * t359;
t442 = pkin(6) * qJDD(2);
t449 = (qJD(1) * t456 + t276) * qJD(2) - t442;
t365 = qJD(1) ^ 2;
t446 = pkin(1) * t365;
t345 = g(1) * t359;
t445 = g(2) * t362;
t440 = qJ(5) * t271;
t439 = qJ(5) * t273;
t354 = qJDD(1) * pkin(1);
t438 = qJDD(2) * pkin(2);
t326 = pkin(4) * t360 + pkin(3);
t435 = t326 * t361;
t433 = t358 * t359;
t431 = t359 * t361;
t428 = t361 * t365;
t395 = t360 * t260 - t275 * t357;
t237 = t395 - t439;
t236 = -pkin(4) * t348 + t237;
t427 = -t237 + t236;
t426 = -t439 + t454;
t425 = t440 - t453;
t336 = t358 * qJD(3);
t421 = qJ(3) * t414 + t336;
t352 = t358 ^ 2;
t353 = t361 ^ 2;
t418 = t352 - t353;
t415 = qJD(2) * t358;
t408 = pkin(4) * t434;
t406 = -g(1) * t432 - g(2) * t433 + g(3) * t361;
t405 = qJD(1) * t363;
t399 = t345 - t445;
t398 = -qJD(2) * pkin(2) + qJD(3);
t393 = t360 * t300 - t301 * t357;
t283 = t361 * pkin(3) - t456;
t391 = t348 ^ 2;
t390 = t362 * pkin(1) + pkin(2) * t429 + t359 * pkin(6) + qJ(3) * t432;
t389 = -t327 - t406;
t388 = t358 * t404;
t364 = qJD(2) ^ 2;
t387 = pkin(6) * t364 + t445;
t385 = pkin(2) * t409 + qJ(3) * t459 + qJD(1) * t336 + t354;
t297 = t332 + t398;
t299 = t333 + t351;
t383 = t297 * t361 - t299 * t358;
t263 = t400 - t438;
t380 = -0.2e1 * pkin(1) * t411 - t442;
t256 = t388 + t421;
t378 = -t357 * t251 - t360 * t252 - t260 * t412 + t275 * t413;
t293 = t447 * t415;
t295 = qJD(2) * t301;
t377 = -t360 * t293 + t357 * t295 + t300 * t412 - t301 * t413;
t373 = -t387 + 0.2e1 * t354;
t371 = -t423 * qJD(4) + t293 * t357 + t360 * t295;
t250 = pkin(2) * t402 - t385;
t269 = pkin(2) * t415 - t421;
t370 = -qJD(1) * t269 - qJDD(1) * t456 - t250 - t387;
t241 = pkin(3) * t409 + qJD(1) * t388 + t385;
t369 = qJD(2) * t383 + t259 * t361 + t263 * t358;
t367 = pkin(4) * t240 + qJDD(5) + t241;
t266 = t285 * t359;
t268 = t285 * t362;
t366 = g(1) * t268 + g(2) * t266 - g(3) * t382 + t258 * t271 + t378;
t356 = -qJ(5) - pkin(7);
t341 = t362 * pkin(6);
t323 = qJ(3) * t416;
t318 = g(1) * t431;
t314 = qJ(3) * t429;
t312 = qJ(3) * t431;
t291 = -pkin(4) + t397;
t290 = pkin(2) * t417 - t323;
t264 = t358 * t405 + t323;
t254 = qJD(2) * t285 - t375;
t253 = -t360 * t415 + t465;
t246 = pkin(4) * t271 + qJD(5) + t258;
t245 = -qJ(5) * t285 + t423;
t244 = qJ(5) * t382 + t393;
t238 = t384 - t440;
t235 = -qJ(5) * t254 + qJD(5) * t382 + t371;
t234 = -qJ(5) * t253 - qJD(5) * t285 + t377;
t233 = -qJ(5) * t240 - qJD(5) * t271 - t378;
t232 = -pkin(4) * t347 + qJ(5) * t239 - qJD(5) * t273 + t464;
t1 = [qJDD(1) * MDP(1) + t399 * MDP(2) + t452 * MDP(3) + (qJDD(1) * t352 + 0.2e1 * t358 * t401) * MDP(4) + 0.2e1 * (t358 * t409 - t411 * t418) * MDP(5) + (qJDD(2) * t358 + t361 * t364) * MDP(6) + (qJDD(2) * t361 - t358 * t364) * MDP(7) + (t358 * t380 + t361 * t373 + t318) * MDP(9) + (t380 * t361 + (-t373 - t345) * t358) * MDP(10) + (t449 * t358 + t370 * t361 + t318) * MDP(11) + ((t352 + t353) * qJDD(1) * pkin(6) + t369 - t452) * MDP(12) + (-t449 * t361 + (t370 + t345) * t358) * MDP(13) + (pkin(6) * t369 - g(1) * t341 - g(2) * t390 + t276 * t269 + (t250 - t345) * t456) * MDP(14) + (t239 * t382 + t254 * t273) * MDP(15) + (t239 * t285 + t240 * t382 - t253 * t273 - t254 * t271) * MDP(16) + (-t254 * t348 + t347 * t382) * MDP(17) + (t253 * t348 + t285 * t347) * MDP(18) + (g(1) * t266 - g(2) * t268 + t283 * t240 + t241 * t285 + t258 * t253 + t256 * t271 - t347 * t393 - t348 * t371) * MDP(20) + (-g(1) * t265 + g(2) * t267 - t283 * t239 - t241 * t382 + t258 * t254 + t256 * t273 + t347 * t423 + t348 * t377) * MDP(21) + (t232 * t382 - t233 * t285 - t234 * t271 - t235 * t273 - t236 * t254 - t238 * t253 + t239 * t244 - t240 * t245 + t452) * MDP(22) + (t233 * t245 + t238 * t234 + t232 * t244 + t236 * t235 + t367 * (pkin(4) * t285 + t283) + t246 * (pkin(4) * t253 + t256) - g(1) * (t356 * t362 + t341) - g(2) * (t326 * t429 + t362 * t408 + t390) + (-g(1) * (t456 - t408 - t435) - g(2) * t356) * t359) * MDP(23); -t358 * MDP(4) * t428 + t418 * MDP(5) * t365 + MDP(6) * t410 + MDP(7) * t409 + qJDD(2) * MDP(8) + (t358 * t446 + t389) * MDP(9) + (g(3) * t358 - t328 + (t452 + t446) * t361) * MDP(10) + (0.2e1 * t438 - qJDD(3) + (-t276 * t358 + t290 * t361) * qJD(1) + t389) * MDP(11) + ((-pkin(2) * t358 + qJ(3) * t361) * qJDD(1) + ((t299 - t351) * t358 + (-t297 + t398) * t361) * qJD(1)) * MDP(12) + (t328 + 0.2e1 * t349 + 0.2e1 * t350 + (qJD(1) * t290 - g(3)) * t358 + (qJD(1) * t276 - t452) * t361) * MDP(13) + (t259 * qJ(3) + t299 * qJD(3) - t263 * pkin(2) - t276 * t290 - g(1) * (-pkin(2) * t432 + t314) - g(2) * (-pkin(2) * t433 + t312) - g(3) * t420 - t383 * qJD(1) * pkin(6)) * MDP(14) + (-t264 * t271 - t397 * t347 + t453 * t348 - t460) * MDP(20) + (-t264 * t273 + t296 * t347 + t454 * t348 - t366) * MDP(21) + (t239 * t291 - t240 * t296 + (-t238 - t425) * t273 + (t236 - t426) * t271) * MDP(22) + (t233 * t296 + t232 * t291 - t246 * (-pkin(4) * t273 + t323) - g(1) * (pkin(4) * t407 + t314) - g(2) * (pkin(4) * t359 * t430 + t312) - g(3) * (t420 + t435) + t426 * t238 + t425 * t236 + (-g(3) * pkin(4) * t357 - t246 * t405 + t452 * (pkin(2) + t326)) * t358) * MDP(23) - t462; -qJDD(2) * MDP(11) + (-t352 * t365 - t364) * MDP(13) + (-qJD(2) * t299 + t263 + t406) * MDP(14) + t406 * MDP(23) + (-MDP(11) * t428 + qJDD(1) * MDP(12) + (MDP(14) * t276 - t271 * MDP(20) - t273 * MDP(21) - MDP(23) * t246) * qJD(1)) * t358 + (-t347 * MDP(20) + t461 * MDP(22) + (-t348 * t238 + t232) * MDP(23) - MDP(21) * t391) * t360 + (t347 * MDP(21) - t455 * MDP(22) + (t348 * t236 + t233) * MDP(23) - MDP(20) * t391) * t357; (-t348 * t384 + t460) * MDP(20) + (-t348 * t395 + t366) * MDP(21) + (pkin(4) * t239 - t427 * t271) * MDP(22) + (t427 * t238 + (-t246 * t273 + t452 * t382 + t232 + t444) * pkin(4)) * MDP(23) + t462; (-t270 - t448) * MDP(22) + (t236 * t273 + t238 * t271 + t367 + t399) * MDP(23);];
tau = t1;
