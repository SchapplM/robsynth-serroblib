% Calculate vector of inverse dynamics joint torques for
% S5RRPRP7
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
%   pkin=[a2,a3,a4,a5,d1,d2,d4,theta3]';
% MDP [23x1]
%   Minimal dynamic parameter vector (fixed base model)
%   see S5RRPRP7_convert_par2_MPV_fixb.m
% 
% Output:
% tau [5x1]
%   joint torques of inverse dynamics (contains inertial, gravitational Coriolis and centrifugal forces)

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 20:02
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tau = S5RRPRP7_invdynJ_fixb_mdp_slag_vp(qJ, qJD, qJDD, g, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(5,1),zeros(3,1),zeros(8,1),zeros(23,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRPRP7_invdynJ_fixb_mdp_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RRPRP7_invdynJ_fixb_mdp_slag_vp: qJD has to be [5x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [5 1]), ...
  'S5RRPRP7_invdynJ_fixb_mdp_slag_vp: qJDD has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RRPRP7_invdynJ_fixb_mdp_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RRPRP7_invdynJ_fixb_mdp_slag_vp: pkin has to be [8x1] (double)');
assert(isreal(MDP) && all(size(MDP) == [23 1]), ...
  'S5RRPRP7_invdynJ_fixb_mdp_slag_vp: MDP has to be [23x1] (double)'); 

%% Symbolic Calculation
% From invdyn_joint_fixb_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 20:01:41
% EndTime: 2019-12-31 20:01:47
% DurationCPUTime: 4.27s
% Computational Cost: add. (3605->424), mult. (8260->535), div. (0->0), fcn. (5791->10), ass. (0->168)
t361 = cos(qJ(2));
t437 = cos(pkin(8));
t398 = t437 * t361;
t337 = qJD(1) * t398;
t355 = sin(pkin(8));
t358 = sin(qJ(2));
t414 = qJD(1) * t358;
t313 = -t355 * t414 + t337;
t305 = qJD(4) - t313;
t352 = qJ(2) + pkin(8);
t346 = sin(t352);
t359 = sin(qJ(1));
t362 = cos(qJ(1));
t391 = g(1) * t362 + g(2) * t359;
t455 = t391 * t346;
t324 = t355 * t361 + t437 * t358;
t440 = t361 * pkin(2);
t345 = pkin(1) + t440;
t380 = -t355 * t358 + t398;
t287 = -pkin(3) * t380 - pkin(7) * t324 - t345;
t439 = qJ(3) + pkin(6);
t331 = t439 * t358;
t332 = t439 * t361;
t297 = -t355 * t331 + t437 * t332;
t357 = sin(qJ(4));
t360 = cos(qJ(4));
t418 = t357 * t287 + t360 * t297;
t454 = g(1) * t359 - g(2) * t362;
t347 = cos(t352);
t447 = pkin(7) * t346;
t453 = pkin(3) * t347 + t447;
t327 = qJD(1) * t332;
t318 = t355 * t327;
t326 = qJD(1) * t331;
t438 = qJD(2) * pkin(2);
t321 = -t326 + t438;
t288 = t437 * t321 - t318;
t278 = -qJD(2) * pkin(3) - t288;
t315 = t324 * qJD(1);
t298 = -t360 * qJD(2) + t315 * t357;
t300 = qJD(2) * t357 + t315 * t360;
t246 = t298 * pkin(4) - t300 * qJ(5) + t278;
t314 = t324 * qJD(2);
t407 = qJDD(1) * t358;
t387 = qJDD(1) * t398 - t355 * t407;
t285 = qJD(1) * t314 + qJDD(4) - t387;
t450 = pkin(2) * t355;
t340 = pkin(7) + t450;
t429 = t340 * t285;
t452 = t305 * t246 - t429;
t409 = qJD(1) * qJD(2);
t402 = t358 * t409;
t366 = qJD(2) * t337 + t324 * qJDD(1) - t355 * t402;
t451 = t300 ^ 2;
t448 = pkin(4) * t285;
t443 = g(3) * t346;
t442 = g(3) * t347;
t441 = g(3) * t361;
t436 = qJ(5) * t285;
t435 = qJDD(1) * pkin(1);
t330 = -t345 * qJD(1) + qJD(3);
t264 = -pkin(3) * t313 - pkin(7) * t315 + t330;
t399 = t437 * t327;
t289 = t355 * t321 + t399;
t279 = qJD(2) * pkin(7) + t289;
t245 = t264 * t357 + t279 * t360;
t434 = t245 * t305;
t408 = qJD(2) * qJD(4);
t412 = qJD(4) * t357;
t258 = -t357 * qJDD(2) + t315 * t412 + (-t366 - t408) * t360;
t433 = t258 * t357;
t432 = t298 * t313;
t431 = t300 * t298;
t395 = t300 * t305;
t430 = t324 * t360;
t428 = t439 * t362;
t275 = t357 * t285;
t427 = t357 * t305;
t426 = t357 * t359;
t425 = t359 * t360;
t276 = t360 * t285;
t424 = t360 * t362;
t423 = t362 * t357;
t365 = -t360 * qJDD(2) + t357 * t366;
t259 = t300 * qJD(4) + t365;
t411 = qJD(4) * t360;
t422 = -t357 * t259 - t298 * t411;
t421 = t305 * t411 + t275;
t420 = t313 * t427 + t276;
t400 = qJD(2) * t439;
t311 = -qJD(3) * t358 - t361 * t400;
t284 = qJDD(2) * pkin(2) + t311 * qJD(1) - qJDD(1) * t331;
t310 = qJD(3) * t361 - t358 * t400;
t291 = t310 * qJD(1) + qJDD(1) * t332;
t255 = t355 * t284 + t437 * t291;
t271 = pkin(2) * t414 + pkin(3) * t315 - pkin(7) * t313;
t293 = -t437 * t326 - t318;
t419 = t357 * t271 + t360 * t293;
t292 = -t326 * t355 + t399;
t388 = pkin(4) * t357 - qJ(5) * t360;
t417 = -qJD(5) * t357 + t305 * t388 - t292;
t416 = (g(1) * t424 + g(2) * t425) * t346;
t353 = t358 ^ 2;
t415 = -t361 ^ 2 + t353;
t413 = qJD(4) * t340;
t244 = t264 * t360 - t279 * t357;
t410 = qJD(5) - t244;
t406 = qJDD(1) * t361;
t405 = pkin(2) * t402 + qJDD(3);
t404 = t358 * t438;
t403 = t437 * pkin(2);
t290 = -qJD(2) * t315 + t387;
t374 = -pkin(2) * t406 + t405 - t435;
t247 = -t290 * pkin(3) - t366 * pkin(7) + t374;
t253 = qJDD(2) * pkin(7) + t255;
t379 = t357 * t247 + t360 * t253 + t264 * t411 - t279 * t412;
t231 = qJD(5) * t305 + t379 + t436;
t238 = -pkin(4) * t305 + t410;
t397 = -t238 * t313 + t231;
t394 = -t360 * t247 + t357 * t253 + t264 * t412 + t279 * t411;
t232 = qJDD(5) + t394 - t448;
t239 = qJ(5) * t305 + t245;
t396 = t239 * t313 + t232;
t269 = t310 * t355 - t437 * t311;
t296 = t437 * t331 + t332 * t355;
t306 = t347 * t426 + t424;
t308 = t347 * t423 - t425;
t393 = -g(1) * t306 + g(2) * t308;
t307 = t347 * t425 - t423;
t309 = t347 * t424 + t426;
t392 = g(1) * t307 - g(2) * t309;
t389 = t360 * pkin(4) + t357 * qJ(5);
t254 = t437 * t284 - t355 * t291;
t386 = t238 * t360 - t239 * t357;
t385 = pkin(3) + t389;
t384 = t305 * t413 + t442;
t383 = -0.2e1 * pkin(1) * t409 - pkin(6) * qJDD(2);
t317 = t380 * qJD(2);
t382 = t317 * t357 + t324 * t411;
t381 = -t317 * t360 + t324 * t412;
t270 = t437 * t310 + t355 * t311;
t272 = pkin(3) * t314 - pkin(7) * t317 + t404;
t378 = t360 * t270 + t357 * t272 + t287 * t411 - t297 * t412;
t377 = t305 * t278 - t429;
t252 = -qJDD(2) * pkin(3) - t254;
t233 = t259 * pkin(4) + t258 * qJ(5) - t300 * qJD(5) + t252;
t376 = -t233 - t384;
t363 = qJD(2) ^ 2;
t373 = -pkin(6) * t363 + 0.2e1 * t435 + t454;
t364 = qJD(1) ^ 2;
t372 = pkin(1) * t364 - pkin(6) * qJDD(1) + t391;
t371 = g(1) * t308 + g(2) * t306 + t357 * t443 - t394;
t369 = t246 * t300 + qJDD(5) - t371;
t368 = -g(1) * t309 - g(2) * t307 - t360 * t443 + t379;
t341 = -t403 - pkin(3);
t335 = t362 * t345;
t322 = -t403 - t385;
t261 = pkin(4) * t300 + qJ(5) * t298;
t260 = t388 * t324 + t296;
t249 = pkin(4) * t380 - t287 * t360 + t297 * t357;
t248 = -qJ(5) * t380 + t418;
t241 = -pkin(4) * t315 - t271 * t360 + t293 * t357;
t240 = qJ(5) * t315 + t419;
t237 = t298 * t305 - t258;
t236 = t388 * t317 + (t389 * qJD(4) - qJD(5) * t360) * t324 + t269;
t235 = -pkin(4) * t314 + qJD(4) * t418 + t270 * t357 - t272 * t360;
t234 = qJ(5) * t314 - qJD(5) * t380 + t378;
t1 = [qJDD(1) * MDP(1) + t454 * MDP(2) + t391 * MDP(3) + (qJDD(1) * t353 + 0.2e1 * t361 * t402) * MDP(4) + 0.2e1 * (t358 * t406 - t415 * t409) * MDP(5) + (qJDD(2) * t358 + t361 * t363) * MDP(6) + (qJDD(2) * t361 - t358 * t363) * MDP(7) + (t383 * t358 + t373 * t361) * MDP(9) + (-t373 * t358 + t383 * t361) * MDP(10) + (-t254 * t324 + t255 * t380 + t269 * t315 + t270 * t313 - t288 * t317 - t289 * t314 + t297 * t290 + t296 * t366 - t391) * MDP(11) + (t255 * t297 + t289 * t270 - t254 * t296 - t288 * t269 - (-t345 * qJDD(1) + t405) * t345 + t330 * t404 - g(1) * (-t345 * t359 + t428) - g(2) * (t359 * t439 + t335)) * MDP(12) + (-t258 * t430 - t381 * t300) * MDP(13) + ((-t298 * t360 - t300 * t357) * t317 + (t433 - t259 * t360 + (t298 * t357 - t300 * t360) * qJD(4)) * t324) * MDP(14) + (t258 * t380 + t324 * t276 + t300 * t314 - t381 * t305) * MDP(15) + (t259 * t380 - t324 * t275 - t298 * t314 - t382 * t305) * MDP(16) + (-t285 * t380 + t305 * t314) * MDP(17) + (t394 * t380 + t244 * t314 + t269 * t298 + t296 * t259 + ((-qJD(4) * t297 + t272) * t305 + t287 * t285 + t278 * qJD(4) * t324) * t360 + ((-qJD(4) * t287 - t270) * t305 - t297 * t285 + t252 * t324 + t278 * t317) * t357 + t392) * MDP(18) + (-t245 * t314 + t252 * t430 - t296 * t258 + t269 * t300 - t381 * t278 - t418 * t285 - t378 * t305 + t379 * t380 + t393) * MDP(19) + (t233 * t324 * t357 + t232 * t380 - t235 * t305 + t236 * t298 - t238 * t314 + t382 * t246 - t249 * t285 + t259 * t260 + t392) * MDP(20) + (-t234 * t298 + t235 * t300 - t248 * t259 - t249 * t258 + t454 * t346 + t386 * t317 + (-t231 * t357 + t232 * t360 + (-t238 * t357 - t239 * t360) * qJD(4)) * t324) * MDP(21) + (-t231 * t380 - t233 * t430 + t234 * t305 - t236 * t300 + t239 * t314 + t381 * t246 + t248 * t285 + t258 * t260 - t393) * MDP(22) + (t231 * t248 + t239 * t234 + t233 * t260 + t246 * t236 + t232 * t249 + t238 * t235 - g(1) * (-pkin(4) * t307 - qJ(5) * t306 + t428) - g(2) * (pkin(4) * t309 + qJ(5) * t308 + t362 * t453 + t335) + (-g(1) * (-t345 - t453) - g(2) * t439) * t359) * MDP(23); MDP(6) * t407 + MDP(7) * t406 + qJDD(2) * MDP(8) + (t372 * t358 - t441) * MDP(9) + (g(3) * t358 + t372 * t361) * MDP(10) + (t290 * t450 - t366 * t403 - (-t289 + t292) * t315 + (t288 - t293) * t313) * MDP(11) + (t288 * t292 - t289 * t293 + (t437 * t254 - t441 + t255 * t355 + (-qJD(1) * t330 + t391) * t358) * pkin(2)) * MDP(12) + (t360 * t395 - t433) * MDP(13) + ((-t258 + t432) * t360 - t300 * t427 + t422) * MDP(14) + (-t305 * t313 * t360 - t300 * t315 + t421) * MDP(15) + (t298 * t315 - t305 * t412 + t420) * MDP(16) - t305 * t315 * MDP(17) + (-t244 * t315 + t341 * t259 - t292 * t298 + (-t442 - t252 + (-t271 - t413) * t305) * t360 + (t293 * t305 + t377) * t357 + t416) * MDP(18) + (-t341 * t258 + t419 * t305 + t245 * t315 - t292 * t300 + t377 * t360 + (t252 + t384 - t455) * t357) * MDP(19) + (t238 * t315 + t241 * t305 + t259 * t322 + t417 * t298 + t357 * t452 + t376 * t360 + t416) * MDP(20) + (-t443 + t240 * t298 - t241 * t300 - t391 * t347 + (-t259 * t340 + (t300 * t340 + t238) * qJD(4) + t397) * t360 + (-t258 * t340 + (t298 * t340 - t239) * qJD(4) + t396) * t357) * MDP(21) + (-t239 * t315 - t240 * t305 + t258 * t322 - t417 * t300 - t452 * t360 + (t376 + t455) * t357) * MDP(22) + (t233 * t322 - t239 * t240 - t238 * t241 - g(3) * (t440 + t447) - t385 * t442 + t417 * t246 + (t386 * qJD(4) + t231 * t360 + t232 * t357) * t340 + t391 * (pkin(2) * t358 - pkin(7) * t347 + t385 * t346)) * MDP(23) + (-MDP(4) * t358 * t361 + MDP(5) * t415) * t364; -t313 ^ 2 * MDP(11) + (-t289 * t313 + t374 - t454) * MDP(12) + t420 * MDP(18) + t422 * MDP(21) + t421 * MDP(22) - t454 * MDP(23) + (-MDP(11) * t315 + t288 * MDP(12) - t246 * MDP(23) + (-MDP(19) + MDP(22)) * t300 + (-MDP(18) - MDP(20)) * t298) * t315 + (t285 * MDP(20) + (t258 + t432) * MDP(21) + (qJD(4) * t239 - t396) * MDP(23) + (-t305 * MDP(19) - t313 * MDP(22)) * t305) * t360 + (-t285 * MDP(19) + (qJD(4) * t238 + t397) * MDP(23) + MDP(21) * t395 + (-qJD(4) * MDP(18) - t305 * MDP(20)) * t305) * t357; MDP(13) * t431 + (-t298 ^ 2 + t451) * MDP(14) + t237 * MDP(15) + (-t315 * t411 - t357 * t408 - t365 + t395) * MDP(16) + t285 * MDP(17) + (-t278 * t300 + t371 + t434) * MDP(18) + (t244 * t305 + t278 * t298 - t368) * MDP(19) + (-t261 * t298 - t369 + t434 + 0.2e1 * t448) * MDP(20) + (pkin(4) * t258 - qJ(5) * t259 + (t239 - t245) * t300 + (t238 - t410) * t298) * MDP(21) + (0.2e1 * t436 - t246 * t298 + t261 * t300 + (0.2e1 * qJD(5) - t244) * t305 + t368) * MDP(22) + (t231 * qJ(5) - t232 * pkin(4) - t246 * t261 - t238 * t245 - g(1) * (-pkin(4) * t308 + qJ(5) * t309) - g(2) * (-pkin(4) * t306 + qJ(5) * t307) + t388 * t443 + t410 * t239) * MDP(23); (-qJDD(4) + t290 + t431) * MDP(20) + t237 * MDP(21) + (-t305 ^ 2 - t451) * MDP(22) + (-t239 * t305 + t369 - t448) * MDP(23);];
tau = t1;
