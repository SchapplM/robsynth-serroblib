% Calculate Coriolis joint torque vector for
% S6PRPRRP3
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [11x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d2,d4,d5,theta1,theta3]';
% MDP [25x1]
%   Minimal dynamic parameter vector (fixed base model)
%   see S6PRPRRP3_convert_par2_MPV_fixb.m
% 
% Output:
% tauc [6x1]
%   joint torques required to compensate Coriolis and centrifugal load

% Quelle: HybrDyn-Toolbox
% Datum: 2021-01-16 01:40
% Revision: 24b2e7d74a0c1a3b64fa2f8f5ad758691ad61af3 (2021-01-15)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauc = S6PRPRRP3_coriolisvecJ_fixb_mdp_slag_vp(qJ, qJD, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(11,1),zeros(25,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRPRRP3_coriolisvecJ_fixb_mdp_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6PRPRRP3_coriolisvecJ_fixb_mdp_slag_vp: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6PRPRRP3_coriolisvecJ_fixb_mdp_slag_vp: pkin has to be [11x1] (double)');
assert(isreal(MDP) && all(size(MDP) == [25 1]), ...
  'S6PRPRRP3_coriolisvecJ_fixb_mdp_slag_vp: MDP has to be [25x1] (double)'); 

%% Symbolic Calculation
% From coriolisvec_joint_fixb_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2021-01-16 01:38:41
% EndTime: 2021-01-16 01:38:54
% DurationCPUTime: 4.74s
% Computational Cost: add. (3806->362), mult. (9758->478), div. (0->0), fcn. (7691->10), ass. (0->156)
t354 = sin(pkin(11));
t356 = cos(pkin(11));
t359 = sin(qJ(4));
t444 = cos(qJ(4));
t375 = -t359 * t354 + t444 * t356;
t441 = qJ(3) + pkin(8);
t339 = t441 * t354;
t340 = t441 * t356;
t449 = -t444 * t339 - t359 * t340;
t278 = t375 * qJD(3) + qJD(4) * t449;
t355 = sin(pkin(6));
t362 = cos(qJ(2));
t423 = t355 * t362;
t368 = t375 * t423;
t308 = qJD(1) * t368;
t461 = t278 - t308;
t330 = t375 * qJD(4);
t336 = t444 * t354 + t359 * t356;
t331 = t336 * qJD(4);
t360 = sin(qJ(2));
t407 = qJD(1) * t355;
t392 = t360 * t407;
t460 = pkin(4) * t331 - pkin(9) * t330 - t392;
t338 = qJD(2) * qJ(3) + t392;
t357 = cos(pkin(6));
t406 = qJD(1) * t357;
t345 = t356 * t406;
t439 = pkin(8) * qJD(2);
t302 = t345 + (-t338 - t439) * t354;
t315 = t356 * t338 + t354 * t406;
t303 = t356 * t439 + t315;
t265 = t359 * t302 + t444 * t303;
t459 = qJD(4) * t265;
t458 = t375 * qJD(2);
t329 = qJD(2) * t336;
t358 = sin(qJ(5));
t361 = cos(qJ(5));
t312 = qJD(4) * t358 + t329 * t361;
t430 = t312 * t358;
t367 = qJD(2) * t330;
t456 = qJD(4) * qJD(5) + t367;
t322 = qJD(5) - t458;
t455 = t322 * t430;
t391 = t362 * t407;
t334 = (qJD(3) + t391) * qJD(2);
t454 = t375 * t334;
t408 = t354 ^ 2 + t356 ^ 2;
t453 = t408 * MDP(6);
t404 = qJD(5) * t358;
t275 = t329 * t404 - t456 * t361;
t261 = qJD(4) * pkin(9) + t265;
t348 = -pkin(3) * t356 - pkin(2);
t383 = qJD(3) - t391;
t321 = t348 * qJD(2) + t383;
t268 = -pkin(4) * t458 - pkin(9) * t329 + t321;
t246 = t261 * t361 + t268 * t358;
t450 = t444 * t302 - t359 * t303;
t251 = qJD(4) * t450 + t454;
t320 = qJD(2) * t331;
t405 = qJD(2) * t360;
t390 = t355 * t405;
t343 = qJD(1) * t390;
t274 = t320 * pkin(4) - pkin(9) * t367 + t343;
t270 = t361 * t274;
t365 = -t246 * qJD(5) - t251 * t358 + t270;
t364 = qJ(6) * t275 + t365;
t443 = pkin(5) * t320;
t234 = -qJD(6) * t312 + t364 + t443;
t310 = -t361 * qJD(4) + t329 * t358;
t241 = -qJ(6) * t310 + t246;
t437 = t241 * t322;
t452 = t234 + t437;
t306 = -t359 * t339 + t444 * t340;
t369 = t336 * t423;
t410 = -qJD(1) * t369 + t336 * qJD(3) + t306 * qJD(4);
t451 = t308 * t358 + t460 * t361;
t295 = -pkin(4) * t375 - pkin(9) * t336 + t348;
t403 = qJD(5) * t361;
t448 = t295 * t403 + t460 * t358 + t461 * t361;
t447 = MDP(20) + MDP(22);
t446 = MDP(21) + MDP(23);
t445 = t312 ^ 2;
t442 = t310 * pkin(5);
t440 = -qJ(6) - pkin(9);
t438 = qJD(2) * pkin(2);
t436 = t275 * t358;
t435 = t310 * t322;
t434 = t310 * t458;
t433 = t310 * t329;
t432 = t312 * t322;
t431 = t312 * t329;
t428 = t458 * t358;
t427 = t458 * t361;
t426 = t336 * t358;
t425 = t336 * t361;
t424 = t355 * t360;
t421 = t358 * t320;
t301 = t361 * t306;
t317 = t361 * t320;
t380 = -qJ(6) * t330 - qJD(6) * t336;
t418 = pkin(5) * t331 - t278 * t358 + t380 * t361 + (-t301 + (qJ(6) * t336 - t295) * t358) * qJD(5) + t451;
t389 = t336 * t403;
t417 = -qJ(6) * t389 + (-qJD(5) * t306 + t380) * t358 + t448;
t245 = -t261 * t358 + t361 * t268;
t240 = -qJ(6) * t312 + t245;
t238 = pkin(5) * t322 + t240;
t416 = t238 - t240;
t293 = pkin(4) * t329 - pkin(9) * t458;
t288 = t361 * t293;
t388 = qJD(5) * t440;
t415 = pkin(5) * t329 - qJ(6) * t427 - t361 * t388 + t288 + (qJD(6) - t450) * t358;
t412 = t358 * t293 + t361 * t450;
t414 = -qJ(6) * t428 - qJD(6) * t361 - t358 * t388 + t412;
t374 = t330 * t358 + t389;
t413 = t374 * pkin(5) + t410;
t276 = t329 * t403 + t456 * t358;
t411 = -t358 * t276 - t310 * t403;
t409 = t358 * t295 + t301;
t387 = qJD(6) + t442;
t386 = t408 * t334;
t385 = t322 * t361;
t384 = t361 * t251 - t261 * t404 + t268 * t403 + t358 * t274;
t252 = t336 * t334 + t459;
t382 = (-t338 * t354 + t345) * t354 - t315 * t356;
t371 = qJ(6) * t276 - t384;
t235 = -qJD(6) * t310 - t371;
t381 = -t322 * t238 + t235;
t379 = t317 + (-t404 + t428) * t322;
t325 = -t354 * t424 + t356 * t357;
t326 = t354 * t357 + t356 * t424;
t284 = t359 * t325 + t444 * t326;
t271 = -t284 * t358 - t361 * t423;
t378 = -t284 * t361 + t358 * t423;
t376 = t444 * t325 - t359 * t326;
t373 = t330 * t361 - t336 * t404;
t239 = pkin(5) * t276 + t252;
t260 = -qJD(4) * pkin(4) - t450;
t372 = -pkin(9) * t320 + t322 * t260;
t366 = qJD(2) * t368;
t363 = qJD(2) ^ 2;
t350 = -pkin(5) * t361 - pkin(4);
t342 = t440 * t361;
t341 = t440 * t358;
t337 = t383 - t438;
t309 = t310 ^ 2;
t291 = t361 * t295;
t280 = pkin(5) * t426 - t449;
t263 = qJD(2) * t369 + t284 * qJD(4);
t262 = t376 * qJD(4) + t366;
t256 = pkin(5) * t428 + t265;
t255 = -qJ(6) * t426 + t409;
t254 = t260 + t387;
t253 = -pkin(5) * t375 - qJ(6) * t425 - t306 * t358 + t291;
t244 = t378 * qJD(5) - t262 * t358 + t361 * t390;
t243 = t271 * qJD(5) + t262 * t361 + t358 * t390;
t1 = [((-t325 * t354 + t326 * t356) * t334 + (t337 * t360 + (-t382 - t392) * t362) * t355 * qJD(2)) * MDP(7) + (-qJD(4) * t263 + (-t320 * t362 - t405 * t458) * t355) * MDP(13) + (t329 * t390 + (-t262 - t366) * qJD(4)) * MDP(14) + (-t243 * t310 - t244 * t312 + t271 * t275 + t276 * t378) * MDP(24) + (t234 * t271 - t235 * t378 + t238 * t244 - t239 * t376 + t241 * t243 + t254 * t263) * MDP(25) + t447 * (t244 * t322 + t263 * t310 + t271 * t320 - t276 * t376) + t446 * (-t243 * t322 + t263 * t312 + t275 * t376 + t320 * t378) + ((-t356 * MDP(5) - MDP(3)) * t360 + (-MDP(4) + t453) * t362) * t355 * t363; (qJD(2) * t383 * t408 + t386) * MDP(6) + (-t382 * qJD(3) + qJ(3) * t386 + (t382 * t362 + (-t337 - t438) * t360) * t407) * MDP(7) + (t329 * t330 + t336 * t367) * MDP(8) + (-t336 * t320 - t329 * t331 + t330 * t458 + t367 * t375) * MDP(9) + (t320 * t348 + t321 * t331) * MDP(13) + (t321 * t330 - t329 * t392 + t336 * t343 + t348 * t367) * MDP(14) + (-t275 * t425 + t373 * t312) * MDP(15) + ((-t310 * t361 - t430) * t330 + (t436 - t276 * t361 + (t310 * t358 - t312 * t361) * qJD(5)) * t336) * MDP(16) + (t275 * t375 + t312 * t331 + t336 * t317 + t373 * t322) * MDP(17) + (t276 * t375 - t310 * t331 - t374 * t322 - t336 * t421) * MDP(18) + (-t320 * t375 + t322 * t331) * MDP(19) + (t291 * t320 - (-t261 * t403 + t270) * t375 + t245 * t331 - t449 * t276 + t260 * t389 + (-t306 * t403 + t451) * t322 + t410 * t310 + ((-qJD(5) * t295 - t278) * t322 - t306 * t320 - (-qJD(5) * t268 - t251) * t375 + t252 * t336 + t260 * t330) * t358) * MDP(20) + (-t409 * t320 + t384 * t375 - t246 * t331 + t449 * t275 + t252 * t425 + (t306 * t404 - t448) * t322 + t410 * t312 + t373 * t260) * MDP(21) + (-t234 * t375 + t238 * t331 + t239 * t426 + t253 * t320 + t374 * t254 + t276 * t280 + t413 * t310 + t418 * t322) * MDP(22) + (t235 * t375 + t239 * t425 - t241 * t331 + t373 * t254 - t255 * t320 - t275 * t280 + t413 * t312 - t417 * t322) * MDP(23) + (t253 * t275 - t255 * t276 + (-t238 * t361 - t241 * t358) * t330 - t418 * t312 - t417 * t310 + (-t234 * t361 - t235 * t358 + (t238 * t358 - t241 * t361) * qJD(5)) * t336) * MDP(24) + (t234 * t253 + t235 * t255 + t418 * t238 + t239 * t280 + t417 * t241 + t413 * t254) * MDP(25) + (t330 * MDP(10) - t331 * MDP(11) - t410 * MDP(13) - t461 * MDP(14)) * qJD(4); -t363 * t453 + (t382 * qJD(2) + t343) * MDP(7) + ((t275 + t434) * t361 + t455 + t411) * MDP(24) + (-t254 * t329 + t381 * t358 + t361 * t452) * MDP(25) + t446 * (-t322 ^ 2 * t361 - t421 - t431) + t447 * (t379 - t433) + (0.2e1 * t329 * MDP(13) + 0.2e1 * t458 * MDP(14)) * qJD(4); -t458 ^ 2 * MDP(9) + (-t252 + t459) * MDP(13) + (-t321 * t458 - t454) * MDP(14) + (t312 * t385 - t436) * MDP(15) + ((-t275 + t434) * t361 - t455 + t411) * MDP(16) + (t322 * t385 + t421 - t431) * MDP(17) + (t379 + t433) * MDP(18) + (-pkin(4) * t276 - t252 * t361 - t265 * t310 + (-pkin(9) * t403 - t288) * t322 + (t322 * t450 + t372) * t358) * MDP(20) + (pkin(4) * t275 + t252 * t358 - t265 * t312 + (pkin(9) * t404 + t412) * t322 + t372 * t361) * MDP(21) + (-t239 * t361 - t256 * t310 + t276 * t350 + t320 * t341 - t415 * t322 + (-t254 * t458 + (t254 + t442) * qJD(5)) * t358) * MDP(22) + (-t254 * t427 + t239 * t358 - t256 * t312 - t275 * t350 + t320 * t342 + t414 * t322 + (pkin(5) * t430 + t254 * t361) * qJD(5)) * MDP(23) + (t275 * t341 + t276 * t342 + t414 * t310 + t415 * t312 - t358 * t452 + t381 * t361) * MDP(24) + (t234 * t341 - t235 * t342 + t239 * t350 + (pkin(5) * t404 - t256) * t254 - t414 * t241 - t415 * t238) * MDP(25) + (-t321 * MDP(13) - t322 * MDP(19) - t245 * MDP(20) + t246 * MDP(21) - t238 * MDP(22) + t241 * MDP(23) - MDP(8) * t458 + MDP(9) * t329) * t329; t312 * t310 * MDP(15) + (-t309 + t445) * MDP(16) + (-t275 + t435) * MDP(17) + (-t276 + t432) * MDP(18) + t320 * MDP(19) + (t246 * t322 - t260 * t312 + t365) * MDP(20) + (t245 * t322 + t260 * t310 - t384) * MDP(21) + (0.2e1 * t443 + t437 + (-t254 - t387) * t312 + t364) * MDP(22) + (-pkin(5) * t445 + t240 * t322 + (qJD(6) + t254) * t310 + t371) * MDP(23) + (pkin(5) * t275 - t416 * t310) * MDP(24) + (t416 * t241 + (-t254 * t312 + t234) * pkin(5)) * MDP(25); (t276 + t432) * MDP(22) + (-t275 - t435) * MDP(23) + (-t309 - t445) * MDP(24) + (t238 * t312 + t241 * t310 + t239) * MDP(25);];
tauc = t1;
