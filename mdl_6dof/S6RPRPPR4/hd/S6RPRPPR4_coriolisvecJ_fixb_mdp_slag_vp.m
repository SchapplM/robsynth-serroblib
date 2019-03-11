% Calculate minimal parameter regressor of Coriolis joint torque vector for
% S6RPRPPR4
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d3,d6,theta2,theta4]';
% MDP [29x1]
%   Minimal dynamic parameter vector (fixed base model)
%   see S6RPRPPR4_convert_par2_MPV_fixb.m
% 
% Output:
% tauc [6x1]
%   joint torques required to compensate Coriolis and centrifugal load

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 02:49
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauc = S6RPRPPR4_coriolisvecJ_fixb_mdp_slag_vp(qJ, qJD, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(10,1),zeros(29,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRPPR4_coriolisvecJ_fixb_mdp_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RPRPPR4_coriolisvecJ_fixb_mdp_slag_vp: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RPRPPR4_coriolisvecJ_fixb_mdp_slag_vp: pkin has to be [10x1] (double)');
assert(isreal(MDP) && all(size(MDP) == [29 1]), ...
  'S6RPRPPR4_coriolisvecJ_fixb_mdp_slag_vp: MDP has to be [29x1] (double)'); 

%% Symbolic Calculation
% From coriolisvec_joint_fixb_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 02:48:29
% EndTime: 2019-03-09 02:48:38
% DurationCPUTime: 4.93s
% Computational Cost: add. (3650->404), mult. (9690->528), div. (0->0), fcn. (7370->8), ass. (0->162)
t406 = sin(pkin(9));
t408 = cos(pkin(9));
t410 = sin(qJ(3));
t466 = t410 * t408;
t483 = cos(qJ(3));
t384 = t406 * t483 + t466;
t405 = sin(pkin(10));
t407 = cos(pkin(10));
t409 = sin(qJ(6));
t411 = cos(qJ(6));
t493 = -t405 * t411 + t407 * t409;
t320 = t493 * t384;
t481 = pkin(7) + qJ(2);
t390 = t481 * t406;
t385 = qJD(1) * t390;
t392 = t481 * t408;
t386 = qJD(1) * t392;
t331 = -t410 * t385 + t483 * t386;
t506 = qJD(3) * t331;
t492 = t384 * qJD(1);
t340 = -t407 * qJD(3) + t405 * t492;
t429 = qJD(3) * t405 + t407 * t492;
t473 = t429 * t409;
t302 = -t411 * t340 + t473;
t449 = qJD(1) * t483;
t396 = t408 * t449;
t467 = t410 * t406;
t451 = qJD(1) * t467;
t370 = -t396 + t451;
t456 = qJD(6) - t370;
t505 = t302 * t456;
t304 = t340 * t409 + t411 * t429;
t504 = t304 * t456;
t503 = t340 * t370;
t458 = qJD(6) * t411;
t459 = qJD(6) * t409;
t464 = t493 * t370 + t405 * t458 - t407 * t459;
t339 = -t410 * t390 + t392 * t483;
t314 = qJD(2) * t384 + qJD(3) * t339;
t502 = t407 * t384 * qJD(5) - t314;
t381 = t405 * t409 + t407 * t411;
t501 = qJD(3) * t492;
t500 = MDP(15) + MDP(19);
t499 = MDP(16) - MDP(21);
t498 = t429 ^ 2;
t424 = t408 * t483 - t467;
t416 = t424 * qJD(2);
t425 = t483 * t385 + t410 * t386;
t305 = qJD(1) * t416 + (qJD(4) - t425) * qJD(3);
t299 = t405 * t305;
t395 = qJD(3) * t396;
t356 = -qJD(3) * t451 + t395;
t376 = t384 * qJD(3);
t357 = qJD(1) * t376;
t301 = pkin(3) * t357 - qJ(4) * t356 - qJD(4) * t492;
t267 = t301 * t407 - t299;
t261 = -pkin(4) * t357 - t267;
t399 = -pkin(2) * t408 - pkin(1);
t387 = qJD(1) * t399 + qJD(2);
t312 = pkin(3) * t370 - qJ(4) * t492 + t387;
t326 = qJD(3) * qJ(4) + t331;
t286 = t405 * t312 + t407 * t326;
t276 = t370 * qJ(5) + t286;
t496 = -t276 * t370 + t261;
t327 = pkin(3) * t492 + qJ(4) * t370;
t293 = t405 * t327 - t407 * t425;
t281 = qJ(5) * t492 + t293;
t461 = qJD(4) * t407;
t495 = -t281 + t461;
t494 = -t483 * t390 - t410 * t392;
t491 = MDP(17) + MDP(20);
t490 = (t406 ^ 2 + t408 ^ 2) * (MDP(7) * qJ(2) + MDP(6));
t465 = t456 * t381;
t489 = -t357 * t493 + t456 * t465;
t427 = qJD(3) * pkin(3) - qJD(4) - t425;
t415 = qJ(5) * t429 + t427;
t284 = pkin(4) * t340 - t415;
t448 = qJ(5) * t405 + pkin(3);
t388 = -pkin(4) * t407 - t448;
t478 = qJ(4) * t357;
t486 = -t356 * t388 + (qJD(4) - t284) * t370 + t478;
t477 = qJ(5) * t407;
t484 = pkin(4) + pkin(5);
t485 = t405 * t484 - t477;
t441 = t493 * t356;
t274 = t304 * qJD(6) + t441;
t365 = t370 ^ 2;
t482 = pkin(8) * t405;
t480 = -pkin(8) + qJ(4);
t475 = t302 * t492;
t474 = t304 * t492;
t349 = t405 * t356;
t350 = t407 * t356;
t268 = t405 * t301 + t407 * t305;
t375 = t424 * qJD(3);
t308 = pkin(3) * t376 - qJ(4) * t375 - qJD(4) * t384;
t313 = qJD(3) * t494 + t416;
t279 = t405 * t308 + t407 * t313;
t328 = -pkin(3) * t424 - qJ(4) * t384 + t399;
t296 = t405 * t328 + t407 * t339;
t462 = qJD(4) * t429;
t460 = qJD(5) * t405;
t454 = qJD(1) * qJD(2);
t289 = -qJ(5) * t424 + t296;
t453 = t340 * t458 + t356 * t381;
t254 = t299 + (-pkin(8) * t356 - t301) * t407 - t484 * t357;
t256 = t357 * qJ(5) + t370 * qJD(5) + t268;
t255 = pkin(8) * t349 + t256;
t446 = t411 * t254 - t255 * t409;
t310 = t405 * t313;
t278 = t308 * t407 - t310;
t285 = t312 * t407 - t405 * t326;
t324 = t405 * t425;
t292 = t327 * t407 + t324;
t333 = t405 * t339;
t295 = t328 * t407 - t333;
t443 = t456 ^ 2;
t442 = -t370 * t485 + t331 + t460;
t262 = t376 * qJ(5) - qJD(5) * t424 + t279;
t307 = t406 * qJD(2) * t449 + t454 * t466 + t506;
t440 = t381 * t357 - t456 * t464;
t439 = qJD(5) - t285;
t438 = pkin(4) * t405 - t477;
t436 = t254 * t409 + t255 * t411;
t259 = -pkin(8) * t429 - t370 * t484 + t439;
t263 = pkin(8) * t340 + t276;
t251 = t259 * t411 - t263 * t409;
t252 = t259 * t409 + t263 * t411;
t270 = t333 + (-pkin(8) * t384 - t328) * t407 + t484 * t424;
t280 = t384 * t482 + t289;
t435 = t270 * t411 - t280 * t409;
t434 = t270 * t409 + t280 * t411;
t433 = -t285 * t405 + t286 * t407;
t391 = t480 * t407;
t423 = qJD(4) * t405 - qJD(6) * t391 + t324 - (pkin(8) * t370 - t327) * t407 + t484 * t492;
t389 = t480 * t405;
t422 = qJD(6) * t389 + t370 * t482 + t495;
t321 = t381 * t384;
t421 = t429 * t459 - t453;
t418 = -pkin(4) * t349 + qJD(5) * t429 - t307;
t265 = -qJ(5) * t350 - t418;
t298 = t384 * t438 - t494;
t420 = t265 * t384 + t284 * t375 + t298 * t356;
t419 = t307 * t384 - t356 * t494 - t375 * t427;
t414 = -pkin(3) * t356 - t478 + (-qJD(4) - t427) * t370;
t377 = t407 * t484 + t448;
t335 = t340 * t461;
t294 = -t370 * t438 + t331;
t291 = -t384 * t485 + t494;
t290 = pkin(4) * t424 - t295;
t288 = qJD(6) * t321 + t375 * t493;
t287 = -qJD(6) * t320 + t375 * t381;
t283 = -pkin(4) * t492 - t292;
t277 = t375 * t438 - t502;
t275 = -pkin(4) * t370 + t439;
t271 = -t340 * t484 + t415;
t269 = -pkin(4) * t376 - t278;
t266 = -t375 * t485 + t502;
t260 = (-pkin(5) * t405 + t477) * t356 + t418;
t258 = t375 * t482 + t262;
t257 = t310 + (-pkin(8) * t375 - t308) * t407 - t484 * t376;
t1 = [(t356 * t384 + t375 * t492) * MDP(8) + (t356 * t424 - t357 * t384 - t370 * t375 - t376 * t492) * MDP(9) + (t357 * t399 + t376 * t387) * MDP(13) + (t356 * t399 + t375 * t387) * MDP(14) + (-t267 * t424 + t278 * t370 + t285 * t376 + t295 * t357 + t314 * t340 + t405 * t419) * MDP(15) + (t268 * t424 - t279 * t370 - t286 * t376 - t296 * t357 + t314 * t429 + t407 * t419) * MDP(16) + (-t278 * t429 - t279 * t340 + (-t267 * t384 - t285 * t375 - t295 * t356) * t407 + (-t268 * t384 - t286 * t375 - t296 * t356) * t405) * MDP(17) + (t267 * t295 + t268 * t296 + t278 * t285 + t279 * t286 - t307 * t494 - t314 * t427) * MDP(18) + (t261 * t424 - t269 * t370 - t275 * t376 + t277 * t340 - t290 * t357 + t405 * t420) * MDP(19) + (-t262 * t340 + t269 * t429 + (t261 * t384 + t275 * t375 + t290 * t356) * t407 + (-t256 * t384 - t276 * t375 - t289 * t356) * t405) * MDP(20) + (-t256 * t424 + t262 * t370 + t276 * t376 - t277 * t429 + t289 * t357 - t407 * t420) * MDP(21) + (t256 * t289 + t261 * t290 + t262 * t276 + t265 * t298 + t269 * t275 + t277 * t284) * MDP(22) + (t287 * t304 - t321 * t421) * MDP(23) + (-t274 * t321 - t287 * t302 - t288 * t304 + t320 * t421) * MDP(24) + (t287 * t456 - t304 * t376 - t321 * t357 - t421 * t424) * MDP(25) + (-t274 * t424 - t288 * t456 + t302 * t376 + t320 * t357) * MDP(26) + (-t357 * t424 - t376 * t456) * MDP(27) + ((t257 * t411 - t258 * t409) * t456 - t435 * t357 + t446 * t424 - t251 * t376 + t266 * t302 + t291 * t274 + t260 * t320 + t271 * t288 + (-t252 * t424 - t434 * t456) * qJD(6)) * MDP(28) + (-(t257 * t409 + t258 * t411) * t456 + t434 * t357 - t436 * t424 + t252 * t376 + t266 * t304 - t291 * t421 + t260 * t321 + t271 * t287 + (-t251 * t424 - t435 * t456) * qJD(6)) * MDP(29) + 0.2e1 * t454 * t490 + (MDP(10) * t375 - MDP(11) * t376 - MDP(13) * t314 - MDP(14) * t313) * qJD(3); 0.2e1 * MDP(13) * t501 + (t395 + (-t370 - t451) * qJD(3)) * MDP(14) + (t267 * t407 + t268 * t405 + t370 * t433 + t427 * t492) * MDP(18) + (t256 * t405 - t261 * t407 - t284 * t492 + (t275 * t405 + t276 * t407) * t370) * MDP(22) + (t440 + t475) * MDP(28) + (t474 + t489) * MDP(29) - qJD(1) ^ 2 * t490 + t491 * ((-t340 * t407 + t405 * t429) * t370 + (-t405 ^ 2 - t407 ^ 2) * t356) + t500 * (-t340 * t492 + t357 * t407 - t365 * t405) - t499 * (t357 * t405 + t365 * t407 + t429 * t492); -t365 * MDP(9) + (t395 + (t370 - t451) * qJD(3)) * MDP(10) + (-t307 + t506) * MDP(13) + (t387 * t370 - t424 * t454) * MDP(14) + (-t292 * t370 - t307 * t407 - t331 * t340 + t405 * t414) * MDP(15) + (t293 * t370 + t307 * t405 - t331 * t429 + t407 * t414) * MDP(16) + (t292 * t429 + t293 * t340 - t335 + (-t285 * t370 + t268) * t407 + (-t286 * t370 - t267 + t462) * t405) * MDP(17) + (-pkin(3) * t307 - t285 * t292 - t286 * t293 + t427 * t331 + t433 * qJD(4) + (-t267 * t405 + t268 * t407) * qJ(4)) * MDP(18) + (-t265 * t407 + t283 * t370 - t294 * t340 + (-qJD(5) * t340 - t486) * t405) * MDP(19) + (t281 * t340 - t283 * t429 - t335 + (t275 * t370 + t256) * t407 + (t462 + t496) * t405) * MDP(20) + (-t265 * t405 - t281 * t370 + (t294 + t460) * t429 + t486 * t407) * MDP(21) + (qJ(4) * t256 * t407 + t265 * t388 - t275 * t283 - t284 * t294 + t495 * t276 + (qJ(4) * t261 + qJD(4) * t275 - qJD(5) * t284) * t405) * MDP(22) + (-t304 * t465 + t421 * t493) * MDP(23) + (t274 * t493 + t302 * t465 - t304 * t464 + t381 * t421) * MDP(24) + (t474 - t489) * MDP(25) + (t440 - t475) * MDP(26) + (-(t389 * t411 - t391 * t409) * t357 + t377 * t274 + t260 * t381 - (t409 * t422 - t411 * t423) * t456 + t442 * t302 + t464 * t271) * MDP(28) + ((t389 * t409 + t391 * t411) * t357 - t377 * t421 - t260 * t493 - (t409 * t423 + t411 * t422) * t456 + t442 * t304 - t465 * t271) * MDP(29) + (-t387 * MDP(13) - t285 * MDP(15) + t286 * MDP(16) + t275 * MDP(19) - t276 * MDP(21) + MDP(27) * t456 + t251 * MDP(28) - t252 * MDP(29) + MDP(8) * t370 + MDP(9) * t492) * t492; (t285 * t429 + t286 * t340 + t307) * MDP(18) + (-t275 * t429 + t276 * t340 + t265) * MDP(22) + (-t274 - t504) * MDP(28) + (t421 + t505) * MDP(29) + t500 * (t370 * t429 + t349) + t499 * (t350 - t503) + t491 * (-t340 ^ 2 - t498); (t340 * t429 - t501) * MDP(19) + (t350 + t503) * MDP(20) + (-t365 - t498) * MDP(21) + (t284 * t429 + t496) * MDP(22) + (-t302 * t429 - t411 * t357 - t409 * t443) * MDP(28) + (-t304 * t429 + t409 * t357 - t411 * t443) * MDP(29); t304 * t302 * MDP(23) + (-t302 ^ 2 + t304 ^ 2) * MDP(24) + (t453 + t505) * MDP(25) + (-t441 + t504) * MDP(26) - t357 * MDP(27) + (t252 * t456 - t271 * t304 + t446) * MDP(28) + (t251 * t456 + t271 * t302 - t436) * MDP(29) + (-MDP(25) * t473 - MDP(26) * t304 - MDP(28) * t252 - MDP(29) * t251) * qJD(6);];
tauc  = t1;
