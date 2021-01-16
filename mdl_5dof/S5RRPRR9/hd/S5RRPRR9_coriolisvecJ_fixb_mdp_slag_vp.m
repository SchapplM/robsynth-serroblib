% Calculate Coriolis joint torque vector for
% S5RRPRR9
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d2,d4,d5,theta3]';
% MDP [28x1]
%   Minimal dynamic parameter vector (fixed base model)
%   see S5RRPRR9_convert_par2_MPV_fixb.m
% 
% Output:
% tauc [5x1]
%   joint torques required to compensate Coriolis and centrifugal load

% Quelle: HybrDyn-Toolbox
% Datum: 2021-01-15 21:48
% Revision: 24b2e7d74a0c1a3b64fa2f8f5ad758691ad61af3 (2021-01-15)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauc = S5RRPRR9_coriolisvecJ_fixb_mdp_slag_vp(qJ, qJD, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(9,1),zeros(28,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRPRR9_coriolisvecJ_fixb_mdp_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RRPRR9_coriolisvecJ_fixb_mdp_slag_vp: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5RRPRR9_coriolisvecJ_fixb_mdp_slag_vp: pkin has to be [9x1] (double)');
assert(isreal(MDP) && all(size(MDP) == [28 1]), ...
  'S5RRPRR9_coriolisvecJ_fixb_mdp_slag_vp: MDP has to be [28x1] (double)'); 

%% Symbolic Calculation
% From coriolisvec_joint_fixb_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2021-01-15 21:47:46
% EndTime: 2021-01-15 21:47:58
% DurationCPUTime: 4.34s
% Computational Cost: add. (3025->373), mult. (7810->507), div. (0->0), fcn. (5805->8), ass. (0->169)
t371 = cos(qJ(2));
t441 = cos(pkin(9));
t395 = t441 * t371;
t354 = qJD(1) * t395;
t365 = sin(pkin(9));
t368 = sin(qJ(2));
t412 = qJD(1) * t368;
t330 = t365 * t412 - t354;
t447 = qJD(4) + qJD(5);
t456 = t330 + t447;
t367 = sin(qJ(4));
t370 = cos(qJ(4));
t406 = t370 * qJD(2);
t343 = t365 * t371 + t441 * t368;
t414 = qJD(1) * t343;
t310 = t367 * t414 - t406;
t369 = cos(qJ(5));
t312 = qJD(2) * t367 + t370 * t414;
t366 = sin(qJ(5));
t431 = t312 * t366;
t254 = t369 * t310 + t431;
t326 = qJD(4) + t330;
t322 = qJD(5) + t326;
t455 = t254 * t322;
t383 = t310 * t366 - t369 * t312;
t454 = t322 * t383;
t346 = t366 * t370 + t367 * t369;
t416 = t456 * t346;
t411 = qJD(4) * t367;
t428 = t330 * t367;
t453 = t411 + t428;
t361 = -pkin(2) * t371 - pkin(1);
t413 = qJD(1) * t361;
t349 = qJD(3) + t413;
t273 = pkin(3) * t330 - pkin(7) * t414 + t349;
t443 = -qJ(3) - pkin(6);
t351 = t443 * t368;
t347 = qJD(1) * t351;
t442 = qJD(2) * pkin(2);
t339 = t347 + t442;
t352 = t443 * t371;
t348 = qJD(1) * t352;
t396 = t441 * t348;
t297 = t365 * t339 - t396;
t292 = qJD(2) * pkin(7) + t297;
t249 = t273 * t367 + t292 * t370;
t244 = -pkin(8) * t310 + t249;
t409 = qJD(5) * t366;
t241 = t244 * t409;
t336 = t365 * t348;
t296 = t441 * t339 + t336;
t291 = -qJD(2) * pkin(3) - t296;
t252 = t310 * pkin(4) + t291;
t452 = t252 * t254 + t241;
t405 = qJD(1) * qJD(2);
t399 = t368 * t405;
t353 = t365 * t399;
t324 = qJD(2) * t354 - t353;
t268 = qJD(4) * t406 + t370 * t324 - t411 * t414;
t332 = t343 * qJD(2);
t323 = qJD(1) * t332;
t355 = pkin(2) * t399;
t272 = pkin(3) * t323 - pkin(7) * t324 + t355;
t264 = t370 * t272;
t398 = qJD(2) * t443;
t329 = qJD(3) * t371 + t368 * t398;
t316 = t329 * qJD(1);
t375 = -qJD(3) * t368 + t371 * t398;
t317 = t375 * qJD(1);
t267 = t441 * t316 + t365 * t317;
t374 = -t249 * qJD(4) - t267 * t367 + t264;
t230 = pkin(4) * t323 - pkin(8) * t268 + t374;
t269 = t312 * qJD(4) + t324 * t367;
t410 = qJD(4) * t370;
t378 = t370 * t267 + t367 * t272 + t273 * t410 - t292 * t411;
t233 = -pkin(8) * t269 + t378;
t393 = t369 * t230 - t366 * t233;
t451 = t252 * t383 + t393;
t450 = t323 * MDP(26) + (-t254 ^ 2 + t383 ^ 2) * MDP(23) - t254 * t383 * MDP(22);
t449 = -0.2e1 * t405;
t287 = t346 * t343;
t448 = (t368 ^ 2 - t371 ^ 2) * MDP(5);
t379 = -t365 * t368 + t395;
t335 = t379 * qJD(2);
t423 = t370 * t335;
t380 = -t343 * t411 + t423;
t345 = t366 * t367 - t369 * t370;
t417 = t456 * t345;
t446 = t417 * t322 - t323 * t346;
t392 = t268 * t366 + t369 * t269;
t237 = -t383 * qJD(5) + t392;
t445 = pkin(2) * t368;
t357 = pkin(2) * t365 + pkin(7);
t444 = pkin(8) + t357;
t248 = t370 * t273 - t292 * t367;
t243 = -pkin(8) * t312 + t248;
t238 = pkin(4) * t326 + t243;
t440 = t238 * t369;
t439 = t244 * t369;
t438 = t254 * t414;
t437 = t383 * t414;
t436 = t268 * t367;
t435 = t310 * t326;
t434 = t310 * t414;
t433 = t312 * t326;
t432 = t312 * t414;
t427 = t367 * t323;
t426 = t367 * t335;
t425 = t367 * t343;
t372 = qJD(2) ^ 2;
t424 = t368 * t372;
t309 = t365 * t351 - t441 * t352;
t302 = t370 * t309;
t422 = t370 * t343;
t421 = t371 * t372;
t373 = qJD(1) ^ 2;
t420 = t371 * t373;
t404 = pkin(2) * t412;
t284 = pkin(3) * t414 + pkin(7) * t330 + t404;
t300 = t441 * t347 + t336;
t419 = t367 * t284 + t370 * t300;
t295 = -pkin(3) * t379 - pkin(7) * t343 + t361;
t418 = t367 * t295 + t302;
t408 = qJD(5) * t369;
t403 = t368 * t442;
t402 = t369 * t268 - t366 * t269 - t310 * t408;
t400 = t343 * t410;
t397 = qJD(4) * t444;
t394 = pkin(1) * t449;
t266 = t316 * t365 - t441 * t317;
t282 = t329 * t365 - t441 * t375;
t299 = t347 * t365 - t396;
t308 = -t441 * t351 - t352 * t365;
t391 = t326 * t370;
t390 = 0.2e1 * t414;
t389 = qJD(5) * t238 + t233;
t358 = -t441 * pkin(2) - pkin(3);
t388 = pkin(4) * t453 - t299;
t387 = -t322 * t416 - t345 * t323;
t277 = t370 * t284;
t341 = t444 * t370;
t386 = pkin(4) * t414 + qJD(5) * t341 - t300 * t367 + t277 + (pkin(8) * t330 + t397) * t370;
t340 = t444 * t367;
t385 = pkin(8) * t428 + qJD(5) * t340 + t367 * t397 + t419;
t232 = t238 * t366 + t439;
t384 = t266 * t343 - t309 * t323;
t382 = t370 * t323 - t326 * t453;
t381 = t400 + t426;
t283 = t441 * t329 + t365 * t375;
t285 = pkin(3) * t332 - pkin(7) * t335 + t403;
t377 = t370 * t283 + t367 * t285 + t295 * t410 - t309 * t411;
t236 = -t312 * t409 + t402;
t376 = t326 * t291 - t357 * t323;
t350 = -t370 * pkin(4) + t358;
t298 = t323 * t379;
t290 = t370 * t295;
t288 = t345 * t343;
t281 = pkin(4) * t425 + t308;
t278 = t370 * t285;
t251 = t381 * pkin(4) + t282;
t250 = -pkin(8) * t425 + t418;
t246 = pkin(4) * t269 + t266;
t245 = -pkin(4) * t379 - pkin(8) * t422 - t309 * t367 + t290;
t240 = -t409 * t425 + (t447 * t422 + t426) * t369 + t380 * t366;
t239 = -t447 * t287 - t345 * t335;
t235 = -t381 * pkin(8) + t377;
t234 = -pkin(8) * t423 + pkin(4) * t332 - t283 * t367 + t278 + (-t302 + (pkin(8) * t343 - t295) * t367) * qJD(4);
t231 = -t244 * t366 + t440;
t1 = [t448 * t449 + (-pkin(6) * t421 + t368 * t394) * MDP(9) + (pkin(6) * t424 + t371 * t394) * MDP(10) + (t323 * t361 + t332 * t349 + (-t282 + (-qJD(1) * t379 + t330) * t445) * qJD(2)) * MDP(11) + (t324 * t361 + t335 * t349 + (t390 * t445 - t283) * qJD(2)) * MDP(12) + (t267 * t379 + t282 * t414 - t283 * t330 - t296 * t335 - t297 * t332 + t308 * t324 + t384) * MDP(13) + (t266 * t308 + t267 * t309 - t282 * t296 + t283 * t297 + (t349 + t413) * t403) * MDP(14) + (t268 * t422 + t380 * t312) * MDP(15) + ((-t310 * t370 - t312 * t367) * t335 + (-t436 - t269 * t370 + (t310 * t367 - t312 * t370) * qJD(4)) * t343) * MDP(16) + (-t268 * t379 + t312 * t332 + t323 * t422 + t380 * t326) * MDP(17) + (t269 * t379 - t310 * t332 - t323 * t425 - t381 * t326) * MDP(18) + (t326 * t332 - t298) * MDP(19) + ((-t309 * t410 + t278) * t326 + t290 * t323 - (-t292 * t410 + t264) * t379 + t248 * t332 + t282 * t310 + t308 * t269 + t291 * t400 + ((-qJD(4) * t295 - t283) * t326 - (-qJD(4) * t273 - t267) * t379 + t291 * t335 + t384) * t367) * MDP(20) + (-t249 * t332 + t266 * t422 + t308 * t268 + t282 * t312 + t380 * t291 - t418 * t323 - t377 * t326 + t378 * t379) * MDP(21) + (-t236 * t288 - t239 * t383) * MDP(22) + (-t236 * t287 + t237 * t288 - t239 * t254 + t240 * t383) * MDP(23) + (-t236 * t379 + t239 * t322 - t288 * t323 - t332 * t383) * MDP(24) + (t237 * t379 - t240 * t322 - t254 * t332 - t287 * t323) * MDP(25) + (t322 * t332 - t298) * MDP(26) + ((t234 * t369 - t235 * t366) * t322 + (t245 * t369 - t250 * t366) * t323 - t393 * t379 + t231 * t332 + t251 * t254 + t281 * t237 + t246 * t287 + t252 * t240 + ((-t245 * t366 - t250 * t369) * t322 + t232 * t379) * qJD(5)) * MDP(27) + (-t232 * t332 + t281 * t236 + t252 * t239 - t241 * t379 - t246 * t288 - t251 * t383 + (-(-qJD(5) * t250 + t234) * t322 - t245 * t323 + t230 * t379) * t366 + (-(qJD(5) * t245 + t235) * t322 - t250 * t323 + t389 * t379) * t369) * MDP(28) + 0.2e1 * t371 * MDP(4) * t399 + MDP(6) * t421 - MDP(7) * t424; -t368 * MDP(4) * t420 + t373 * t448 + (qJD(2) * t299 - t330 * t404 - t266) * MDP(11) + (t300 * qJD(2) + t349 * t330 - t267) * MDP(12) + ((-t296 + t300) * t330 + (-t323 * t365 - t441 * t324) * pkin(2)) * MDP(13) + (t296 * t299 - t297 * t300 + (-t441 * t266 + t267 * t365 - t349 * t412) * pkin(2)) * MDP(14) + (t312 * t391 + t436) * MDP(15) + ((t268 - t435) * t370 + (-t269 - t433) * t367) * MDP(16) + (t326 * t391 + t427 - t432) * MDP(17) + (t382 + t434) * MDP(18) + (-t266 * t370 + t358 * t269 - t299 * t310 + (-t357 * t410 - t277) * t326 + (t300 * t326 + t376) * t367) * MDP(20) + (t266 * t367 + t358 * t268 - t299 * t312 + (t357 * t411 + t419) * t326 + t376 * t370) * MDP(21) + (t236 * t346 + t383 * t417) * MDP(22) + (-t236 * t345 - t237 * t346 + t417 * t254 + t383 * t416) * MDP(23) + (t437 - t446) * MDP(24) + (t387 + t438) * MDP(25) + ((-t340 * t369 - t341 * t366) * t323 + t350 * t237 + t246 * t345 + (t385 * t366 - t386 * t369) * t322 + t388 * t254 + t416 * t252) * MDP(27) + (-(-t340 * t366 + t341 * t369) * t323 + t350 * t236 + t246 * t346 + (t386 * t366 + t385 * t369) * t322 - t388 * t383 - t417 * t252) * MDP(28) - (t349 * MDP(11) + t404 * MDP(12) + (-t297 + t299) * MDP(13) + t326 * MDP(19) + t248 * MDP(20) - t249 * MDP(21) + t322 * MDP(26) + t231 * MDP(27) - t232 * MDP(28)) * t414 + (t373 * t368 * MDP(9) + MDP(10) * t420) * pkin(1); t390 * qJD(2) * MDP(11) + (-t353 + (t354 - t330) * qJD(2)) * MDP(12) + (-t330 ^ 2 - t414 ^ 2) * MDP(13) + (t296 * t414 + t297 * t330 + t355) * MDP(14) + (t382 - t434) * MDP(20) + (-t326 ^ 2 * t370 - t427 - t432) * MDP(21) + (t387 - t438) * MDP(27) + (t437 + t446) * MDP(28); t312 * t310 * MDP(15) + (-t310 ^ 2 + t312 ^ 2) * MDP(16) + (t268 + t435) * MDP(17) + (-t269 + t433) * MDP(18) + t323 * MDP(19) + (t249 * t326 - t291 * t312 + t374) * MDP(20) + (t248 * t326 + t291 * t310 - t378) * MDP(21) + (t236 + t455) * MDP(24) + (-t237 - t454) * MDP(25) + (-(-t243 * t366 - t439) * t322 - t232 * qJD(5) + (-t254 * t312 - t322 * t409 + t369 * t323) * pkin(4) + t451) * MDP(27) + ((-t244 * t322 - t230) * t366 + (t243 * t322 - t389) * t369 + (t312 * t383 - t322 * t408 - t366 * t323) * pkin(4) + t452) * MDP(28) + t450; (t402 + t455) * MDP(24) + (-t392 - t454) * MDP(25) + (t232 * t322 + t451) * MDP(27) + (-t366 * t230 + t231 * t322 - t369 * t233 + t452) * MDP(28) + (-MDP(24) * t431 + t383 * MDP(25) - t232 * MDP(27) - MDP(28) * t440) * qJD(5) + t450;];
tauc = t1;
