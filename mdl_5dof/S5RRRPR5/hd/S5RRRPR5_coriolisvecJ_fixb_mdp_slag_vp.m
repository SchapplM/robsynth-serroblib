% Calculate minimal parameter regressor of Coriolis joint torque vector for
% S5RRRPR5
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d2,d3,d5,theta4]';
% MDP [26x1]
%   Minimal dynamic parameter vector (fixed base model)
%   see S5RRRPR5_convert_par2_MPV_fixb.m
% 
% Output:
% tauc [5x1]
%   joint torques required to compensate Coriolis and centrifugal load

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 21:15
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauc = S5RRRPR5_coriolisvecJ_fixb_mdp_slag_vp(qJ, qJD, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(9,1),zeros(26,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRRPR5_coriolisvecJ_fixb_mdp_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RRRPR5_coriolisvecJ_fixb_mdp_slag_vp: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5RRRPR5_coriolisvecJ_fixb_mdp_slag_vp: pkin has to be [9x1] (double)');
assert(isreal(MDP) && all(size(MDP) == [26 1]), ...
  'S5RRRPR5_coriolisvecJ_fixb_mdp_slag_vp: MDP has to be [26x1] (double)'); 

%% Symbolic Calculation
% From coriolisvec_joint_fixb_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 21:14:24
% EndTime: 2019-12-31 21:14:30
% DurationCPUTime: 2.67s
% Computational Cost: add. (3624->309), mult. (9497->429), div. (0->0), fcn. (6982->8), ass. (0->164)
t352 = qJD(2) + qJD(3);
t359 = cos(qJ(5));
t356 = sin(qJ(5));
t360 = cos(qJ(3));
t361 = cos(qJ(2));
t415 = qJD(1) * t361;
t407 = t360 * t415;
t357 = sin(qJ(3));
t358 = sin(qJ(2));
t416 = qJD(1) * t358;
t408 = t357 * t416;
t313 = -t407 + t408;
t315 = -t357 * t415 - t360 * t416;
t355 = sin(pkin(9));
t440 = cos(pkin(9));
t375 = -t355 * t313 - t440 * t315;
t430 = t375 * t356;
t282 = -t359 * t352 + t430;
t401 = -t440 * t313 + t315 * t355;
t446 = qJD(5) - t401;
t455 = t282 * t446;
t284 = t352 * t356 + t359 * t375;
t454 = t284 * t446;
t397 = t446 * t359;
t410 = qJD(1) * qJD(2);
t406 = t361 * t410;
t294 = qJD(3) * t407 - t352 * t408 + t360 * t406;
t327 = t357 * t361 + t358 * t360;
t447 = qJD(1) * t327;
t366 = t352 * t447;
t268 = t294 * t355 + t440 * t366;
t434 = t268 * t356;
t453 = -t446 * t397 - t434;
t452 = -0.2e1 * t410;
t450 = t358 * MDP(4);
t444 = pkin(6) + pkin(7);
t334 = t444 * t358;
t328 = qJD(1) * t334;
t441 = qJD(2) * pkin(2);
t322 = -t328 + t441;
t409 = qJD(2) * t444;
t391 = qJD(1) * t409;
t323 = t358 * t391;
t449 = (qJD(3) * t322 - t323) * t360;
t448 = (t358 ^ 2 - t361 ^ 2) * MDP(5);
t310 = t315 * qJ(4);
t335 = t444 * t361;
t330 = qJD(1) * t335;
t316 = t357 * t330;
t400 = t360 * t322 - t316;
t279 = t310 + t400;
t324 = t361 * t391;
t414 = qJD(3) * t357;
t399 = -t357 * t324 - t330 * t414;
t246 = -qJ(4) * t366 - t313 * qJD(4) + t399 + t449;
t383 = t357 * t323 - t360 * t324;
t320 = t360 * t330;
t384 = -t322 * t357 - t320;
t369 = t384 * qJD(3) + t383;
t365 = -qJ(4) * t294 + qJD(4) * t315 + t369;
t235 = t440 * t246 + t355 * t365;
t373 = t327 * qJD(3);
t299 = t327 * qJD(2) + t373;
t326 = t357 * t358 - t360 * t361;
t329 = t358 * t409;
t331 = t361 * t409;
t428 = t334 * t360;
t374 = -qJD(3) * t428 - t360 * t329 - t357 * t331 - t335 * t414;
t257 = -qJ(4) * t299 - qJD(4) * t326 + t374;
t298 = t352 * t326;
t382 = t334 * t357 - t335 * t360;
t370 = t382 * qJD(3) + t329 * t357 - t360 * t331;
t367 = qJ(4) * t298 - qJD(4) * t327 + t370;
t239 = t440 * t257 + t355 * t367;
t275 = pkin(3) * t352 + t279;
t439 = qJ(4) * t313;
t280 = -t384 - t439;
t427 = t355 * t280;
t251 = t440 * t275 - t427;
t249 = -t352 * pkin(4) - t251;
t349 = -pkin(2) * t361 - pkin(1);
t333 = t349 * qJD(1);
t300 = pkin(3) * t313 + qJD(4) + t333;
t256 = -pkin(4) * t401 - pkin(8) * t375 + t300;
t296 = t440 * t326 + t327 * t355;
t297 = -t355 * t326 + t440 * t327;
t381 = pkin(3) * t326 + t349;
t262 = pkin(4) * t296 - pkin(8) * t297 + t381;
t271 = -t440 * t298 - t355 * t299;
t234 = t246 * t355 - t440 * t365;
t288 = -qJ(4) * t326 - t382;
t372 = -qJ(4) * t327 - t335 * t357 - t428;
t264 = t440 * t288 + t355 * t372;
t389 = t234 * t297 - t264 * t268;
t445 = t249 * t271 - (qJD(5) * t262 + t239) * t446 - (qJD(5) * t256 + t235) * t296 + t389;
t443 = pkin(3) * t315;
t442 = pkin(2) * qJD(3);
t413 = qJD(5) * t356;
t269 = t440 * t294 - t355 * t366;
t412 = qJD(5) * t359;
t421 = t359 * t269 + t352 * t412;
t243 = -t375 * t413 + t421;
t438 = t243 * t356;
t437 = t249 * t401;
t436 = t249 * t297;
t435 = t262 * t268;
t433 = t269 * t356;
t432 = t282 * t375;
t431 = t284 * t375;
t429 = t333 * t315;
t426 = t355 * t357;
t362 = qJD(2) ^ 2;
t425 = t358 * t362;
t424 = t359 * t297;
t423 = t361 * t362;
t363 = qJD(1) ^ 2;
t422 = t361 * t363;
t276 = t440 * t280;
t252 = t355 * t275 + t276;
t418 = -t360 * t328 - t316;
t285 = t310 + t418;
t398 = t328 * t357 - t320;
t379 = t398 + t439;
t403 = t440 * t357;
t420 = -t285 * t355 + t440 * t379 + (t355 * t360 + t403) * t442;
t419 = -t440 * t285 - t355 * t379 + (t440 * t360 - t426) * t442;
t348 = pkin(2) * t360 + pkin(3);
t309 = pkin(2) * t403 + t355 * t348;
t351 = t358 * t441;
t350 = pkin(2) * t416;
t405 = -pkin(2) * t352 - t322;
t404 = pkin(3) * t299 + t351;
t402 = pkin(1) * t452;
t261 = pkin(4) * t375 - pkin(8) * t401 - t443;
t304 = pkin(8) + t309;
t393 = qJD(5) * t304 + t261 + t350;
t345 = pkin(3) * t355 + pkin(8);
t392 = qJD(5) * t345 + t261;
t250 = pkin(8) * t352 + t252;
t237 = t250 * t359 + t256 * t356;
t390 = t234 * t356 + t237 * t375 + t249 * t412;
t388 = -t268 * t304 - t437;
t387 = -t268 * t345 - t437;
t386 = t250 * t356 - t256 * t359;
t385 = t251 * t401 + t252 * t375;
t380 = t359 * t268 + (t356 * t401 - t413) * t446;
t378 = -t234 * t359 + t249 * t413 + t375 * t386;
t377 = t333 * t313 - t399;
t376 = t271 * t359 - t297 * t413;
t308 = -pkin(2) * t426 + t440 * t348;
t244 = t284 * qJD(5) + t433;
t368 = -t315 * t313 * MDP(11) - t446 * t375 * MDP(24) + ((t243 - t455) * t359 + (-t244 - t454) * t356) * MDP(21) + (t380 + t432) * MDP(23) + (-t431 - t453) * MDP(22) + (t284 * t397 + t438) * MDP(20) + (t313 * t352 + t294) * MDP(13) + (-t315 * t352 - t366) * MDP(14) + (-t313 ^ 2 + t315 ^ 2) * MDP(12);
t364 = pkin(3) * t366 + qJD(2) * t350;
t346 = -t440 * pkin(3) - pkin(4);
t303 = -pkin(4) - t308;
t270 = -t298 * t355 + t440 * t299;
t263 = t288 * t355 - t440 * t372;
t254 = t440 * t279 - t427;
t253 = t279 * t355 + t276;
t242 = pkin(4) * t270 - pkin(8) * t271 + t404;
t241 = t268 * pkin(4) - t269 * pkin(8) + t364;
t240 = t359 * t241;
t238 = t257 * t355 - t440 * t367;
t1 = [0.2e1 * t406 * t450 + t448 * t452 + MDP(6) * t423 - MDP(7) * t425 + (-pkin(6) * t423 + t358 * t402) * MDP(9) + (pkin(6) * t425 + t361 * t402) * MDP(10) + (t294 * t327 + t298 * t315) * MDP(11) + (-t294 * t326 + t298 * t313 + t315 * t299 - t327 * t366) * MDP(12) + (t313 * t351 + t333 * t299 + (t349 * t373 + (t358 * pkin(2) * t326 + t349 * t327) * qJD(2)) * qJD(1)) * MDP(16) + (t349 * t294 - t333 * t298 + (-t315 + t447) * t351) * MDP(17) + (-t235 * t296 + t238 * t375 + t239 * t401 - t251 * t271 - t252 * t270 + t263 * t269 + t389) * MDP(18) + (t234 * t263 + t235 * t264 - t251 * t238 + t252 * t239 + t300 * t404 + t364 * t381) * MDP(19) + (t243 * t424 + t376 * t284) * MDP(20) + ((-t282 * t359 - t284 * t356) * t271 + (-t438 - t244 * t359 + (t282 * t356 - t284 * t359) * qJD(5)) * t297) * MDP(21) + (t243 * t296 + t268 * t424 + t270 * t284 + t376 * t446) * MDP(22) + (-t297 * t434 - t244 * t296 - t270 * t282 + (-t271 * t356 - t297 * t412) * t446) * MDP(23) + (t268 * t296 + t270 * t446) * MDP(24) + (-t386 * t270 + t238 * t282 + t240 * t296 + t263 * t244 + (t242 * t446 + t435 + (-t250 * t296 - t264 * t446 + t436) * qJD(5)) * t359 + t445 * t356) * MDP(25) + (-t237 * t270 + t238 * t284 + t263 * t243 + (-(-qJD(5) * t264 + t242) * t446 - t435 - (-qJD(5) * t250 + t241) * t296 - qJD(5) * t436) * t356 + t445 * t359) * MDP(26) + (-t298 * MDP(13) - t299 * MDP(14) + t370 * MDP(16) - t374 * MDP(17)) * t352; t368 + (t303 * t244 + t388 * t356 + t420 * t282 + (-t419 * t356 - t393 * t359) * t446 + t378) * MDP(25) + (t303 * t243 + t388 * t359 + t420 * t284 + (t393 * t356 - t419 * t359) * t446 + t390) * MDP(26) + t363 * t448 + (t315 * t350 + t418 * t352 + (t405 * qJD(3) + t323) * t360 + t377) * MDP(17) + (-t313 * t350 + t429 - t398 * t352 + (t405 * t357 - t320) * qJD(3) + t383) * MDP(16) + (t235 * t309 - t234 * t308 - t300 * (t350 - t443) + t419 * t252 - t420 * t251) * MDP(19) + (-t268 * t309 - t269 * t308 + t375 * t420 + t401 * t419 + t385) * MDP(18) - t422 * t450 + (t363 * t358 * MDP(9) + MDP(10) * t422) * pkin(1); t368 + (t346 * t243 - t253 * t284 + t387 * t359 + (t359 * t254 + t392 * t356) * t446 + t390) * MDP(26) + (-t253 * t375 - t254 * t401 + (-t268 * t355 - t440 * t269) * pkin(3) + t385) * MDP(18) + (-t384 * t352 + t369 + t429) * MDP(16) + (t400 * t352 + t377 - t449) * MDP(17) + (t251 * t253 - t252 * t254 + (-t440 * t234 + t235 * t355 + t300 * t315) * pkin(3)) * MDP(19) + (t346 * t244 - t253 * t282 + t387 * t356 + (t356 * t254 - t392 * t359) * t446 + t378) * MDP(25); (-t375 ^ 2 - t401 ^ 2) * MDP(18) + (t251 * t375 - t252 * t401 + t364) * MDP(19) + (t380 - t432) * MDP(25) + (-t431 + t453) * MDP(26); t284 * t282 * MDP(20) + (-t282 ^ 2 + t284 ^ 2) * MDP(21) + (t421 + t455) * MDP(22) + (-t433 + t454) * MDP(23) + t268 * MDP(24) + (-t235 * t356 + t237 * t446 - t249 * t284 + t240) * MDP(25) + (-t235 * t359 - t241 * t356 + t249 * t282 - t386 * t446) * MDP(26) + (-MDP(22) * t430 - t284 * MDP(23) - t237 * MDP(25) + t386 * MDP(26)) * qJD(5);];
tauc = t1;
