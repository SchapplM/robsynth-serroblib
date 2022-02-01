% Calculate vector of inverse dynamics joint torques for
% S5RRRRR5
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
%   pkin=[a2,a3,a4,a5,d1,d2,d3,d4,d5]';
% MDP [23x1]
%   Minimal dynamic parameter vector (fixed base model)
%   see S5RRRRR5_convert_par2_MPV_fixb.m
% 
% Output:
% tau [5x1]
%   joint torques of inverse dynamics (contains inertial, gravitational Coriolis and centrifugal forces)

% Quelle: HybrDyn-Toolbox
% Datum: 2022-01-20 12:02
% Revision: fd3771346c4aea32fdeb66112c511235427c26a7 (2022-01-20)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tau = S5RRRRR5_invdynJ_fixb_mdp_slag_vp(qJ, qJD, qJDD, g, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(5,1),zeros(3,1),zeros(9,1),zeros(23,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRRRR5_invdynJ_fixb_mdp_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RRRRR5_invdynJ_fixb_mdp_slag_vp: qJD has to be [5x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [5 1]), ...
  'S5RRRRR5_invdynJ_fixb_mdp_slag_vp: qJDD has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RRRRR5_invdynJ_fixb_mdp_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5RRRRR5_invdynJ_fixb_mdp_slag_vp: pkin has to be [9x1] (double)');
assert(isreal(MDP) && all(size(MDP) == [23 1]), ...
  'S5RRRRR5_invdynJ_fixb_mdp_slag_vp: MDP has to be [23x1] (double)'); 

%% Symbolic Calculation
% From invdyn_joint_fixb_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2022-01-20 12:02:12
% EndTime: 2022-01-20 12:02:18
% DurationCPUTime: 2.52s
% Computational Cost: add. (2230->290), mult. (3035->373), div. (0->0), fcn. (1907->16), ass. (0->169)
t367 = cos(qJ(2));
t461 = pkin(1) * t367;
t339 = qJDD(1) * t461;
t352 = qJDD(1) + qJDD(2);
t362 = sin(qJ(2));
t453 = pkin(1) * qJD(1);
t422 = t362 * t453;
t280 = pkin(2) * t352 - qJD(2) * t422 + t339;
t354 = qJD(1) + qJD(2);
t430 = qJD(1) * t367;
t421 = pkin(1) * t430;
t303 = pkin(2) * t354 + t421;
t361 = sin(qJ(3));
t404 = qJD(3) * t422;
t314 = t361 * t404;
t366 = cos(qJ(3));
t413 = qJD(2) * t430;
t423 = qJDD(1) * t362;
t380 = (t413 + t423) * pkin(1);
t468 = -t361 * t280 - (qJD(3) * t303 + t380) * t366 + t314;
t438 = t361 * t362;
t325 = pkin(1) * t438;
t428 = qJD(3) * t361;
t399 = t361 * pkin(1) * t413 + qJDD(1) * t325 + t303 * t428 + (-t280 + t404) * t366;
t343 = qJDD(3) + t352;
t459 = pkin(3) * t343;
t243 = t399 - t459;
t358 = qJ(1) + qJ(2);
t349 = qJ(3) + t358;
t334 = cos(t349);
t456 = g(2) * t334;
t467 = t243 + t456;
t364 = cos(qJ(5));
t365 = cos(qJ(4));
t436 = t364 * t365;
t359 = sin(qJ(5));
t360 = sin(qJ(4));
t442 = t359 * t360;
t288 = -t436 + t442;
t289 = t359 * t365 + t360 * t364;
t344 = qJD(3) + t354;
t353 = qJD(4) + qJD(5);
t465 = t353 * t289;
t245 = t288 * t343 + t344 * t465;
t385 = t353 * t288;
t437 = t362 * t366;
t390 = t361 * t367 + t437;
t283 = t390 * t453;
t402 = pkin(2) * t428 - t283;
t333 = sin(t349);
t464 = g(1) * t334 + g(2) * t333;
t335 = pkin(2) * t361 + pkin(8);
t460 = pkin(2) * t366;
t336 = -pkin(3) - t460;
t369 = qJD(4) ^ 2;
t463 = t335 * t369 + t336 * t343 + t402 * t344;
t462 = -pkin(8) - pkin(9);
t458 = pkin(3) * t344;
t457 = pkin(4) * t365;
t323 = g(1) * t333;
t337 = pkin(2) + t461;
t433 = pkin(1) * t437 + t361 * t337;
t282 = pkin(8) + t433;
t455 = -pkin(9) - t282;
t454 = -pkin(9) - t335;
t274 = t303 * t361 + t366 * t422;
t267 = pkin(8) * t344 + t274;
t412 = pkin(9) * t344 + t267;
t256 = t412 * t365;
t452 = t256 * t364;
t427 = qJD(3) * t366;
t260 = t337 * t428 + (t390 * qJD(2) + t362 * t427) * pkin(1);
t451 = t260 * t344;
t450 = t274 * t344;
t357 = qJ(4) + qJ(5);
t345 = sin(t357);
t449 = t333 * t345;
t347 = cos(t357);
t448 = t333 * t347;
t447 = t334 * t345;
t446 = t334 * t347;
t445 = t343 * t360;
t444 = t343 * t365;
t443 = t344 * t360;
t439 = t360 * t365;
t320 = t361 * t422;
t273 = t303 * t366 - t320;
t266 = -t273 - t458;
t426 = qJD(4) * t360;
t435 = t266 * t426 + t365 * t323;
t341 = pkin(4) * t426;
t434 = t341 + t402;
t346 = sin(t358);
t348 = cos(t358);
t432 = g(1) * t348 + g(2) * t346;
t355 = t360 ^ 2;
t431 = -t365 ^ 2 + t355;
t425 = qJD(4) * t365;
t424 = qJD(5) * t359;
t419 = pkin(2) * t427;
t418 = t344 * t442;
t417 = t344 * t436;
t416 = t266 * t425 + t467 * t360;
t338 = -pkin(3) - t457;
t415 = qJD(4) * t462;
t414 = t344 * t425;
t410 = qJD(4) * t455;
t409 = qJD(4) * t454;
t407 = t337 * t366 - t325;
t406 = qJD(1) * (-qJD(2) + t354);
t405 = qJD(2) * (-qJD(1) - t354);
t281 = -pkin(3) - t407;
t403 = g(1) * t346 - g(2) * t348 + t339;
t401 = -t274 + t341;
t255 = t412 * t360;
t254 = qJD(4) * pkin(4) - t255;
t398 = -t254 * t359 - t452;
t268 = t455 * t360;
t350 = t365 * pkin(9);
t269 = t282 * t365 + t350;
t397 = t268 * t364 - t269 * t359;
t396 = t268 * t359 + t269 * t364;
t286 = t454 * t360;
t287 = t335 * t365 + t350;
t395 = t286 * t364 - t287 * t359;
t394 = t286 * t359 + t287 * t364;
t316 = t462 * t360;
t317 = pkin(8) * t365 + t350;
t393 = t316 * t364 - t317 * t359;
t392 = t316 * t359 + t317 * t364;
t389 = t344 * t426 - t444;
t240 = t389 * pkin(4) + t243;
t258 = t338 * t344 - t273;
t388 = -g(1) * t449 + g(2) * t447 + t240 * t289 - t258 * t385;
t387 = g(1) * t448 - g(2) * t446 + t240 * t288 + t258 * t465;
t383 = pkin(8) * t369 - t450 - t459;
t382 = t281 * t343 + t282 * t369 + t451;
t244 = qJD(5) * t417 + t289 * t343 - t353 * t418 + t364 * t414;
t275 = -t417 + t418;
t277 = t289 * t344;
t351 = qJDD(4) + qJDD(5);
t381 = t277 * t275 * MDP(17) + (t275 * t353 + t244) * MDP(19) + (t277 * t353 - t245) * MDP(20) + (-t275 ^ 2 + t277 ^ 2) * MDP(18) + t351 * MDP(21);
t379 = -pkin(8) * qJDD(4) + (t273 - t458) * qJD(4);
t242 = pkin(8) * t343 - t468;
t378 = -t266 * t344 - t242 + t464;
t377 = t323 - t399 - t456;
t259 = t337 * t427 + (-t362 * t428 + (t366 * t367 - t438) * qJD(2)) * pkin(1);
t376 = -qJDD(4) * t282 + (t281 * t344 - t259) * qJD(4);
t375 = (-t244 * t288 - t245 * t289 + t275 * t385 - t277 * t465) * MDP(18) + (t244 * t289 - t277 * t385) * MDP(17) + (t289 * t351 - t353 * t385) * MDP(19) + (-t288 * t351 - t353 * t465) * MDP(20) + 0.2e1 * (-t431 * t344 * qJD(4) + t343 * t439) * MDP(11) + (t343 * t355 + 0.2e1 * t360 * t414) * MDP(10) + (qJDD(4) * t360 + t365 * t369) * MDP(12) + (qJDD(4) * t365 - t360 * t369) * MDP(13) + t343 * MDP(7);
t374 = t352 * MDP(4) + t375;
t284 = t366 * t421 - t320;
t373 = -qJDD(4) * t335 + (t336 * t344 + t284 - t419) * qJD(4);
t234 = -t267 * t425 + qJDD(4) * pkin(4) - t242 * t360 + (-t414 - t445) * pkin(9);
t372 = t258 * t275 + t256 * t424 + g(2) * t448 + g(1) * t446 + g(3) * t345 + (-t256 * t353 - t234) * t359;
t235 = -t389 * pkin(9) + t242 * t365 - t267 * t426;
t371 = g(1) * t447 + g(2) * t449 - g(3) * t347 + t398 * qJD(5) + t364 * t234 - t359 * t235 - t258 * t277;
t370 = t464 + t468;
t368 = cos(qJ(1));
t363 = sin(qJ(1));
t310 = t338 - t460;
t295 = t365 * t415;
t294 = t360 * t415;
t279 = t281 - t457;
t271 = -t360 * t419 + t365 * t409;
t270 = t360 * t409 + t365 * t419;
t257 = t341 + t260;
t250 = -t259 * t360 + t365 * t410;
t249 = t259 * t365 + t360 * t410;
t1 = [(-t259 * t344 - t433 * t343 + t370) * MDP(9) + (t376 * t365 + (t382 - t323) * t360 + t416) * MDP(16) + (t407 * t343 + t377 - t451) * MDP(8) + (g(1) * t363 - g(2) * t368) * MDP(2) + (g(1) * t368 + g(2) * t363) * MDP(3) + (t376 * t360 + (-t382 - t467) * t365 + t435) * MDP(15) + (((-qJDD(1) - t352) * t362 + t367 * t405) * pkin(1) + t432) * MDP(6) + ((t352 * t367 + t362 * t405) * pkin(1) + t403) * MDP(5) + (t257 * t277 + t279 * t244 - (t397 * qJD(5) + t249 * t364 + t250 * t359) * t353 - t396 * t351 + t388) * MDP(23) + (t257 * t275 + t279 * t245 + (-t396 * qJD(5) - t249 * t359 + t250 * t364) * t353 + t397 * t351 + t387) * MDP(22) + qJDD(1) * MDP(1) + t374; (t373 * t360 + (-t467 - t463) * t365 + t435) * MDP(15) + (t373 * t365 + (-t323 + t463) * t360 + t416) * MDP(16) + (t362 * pkin(1) * t406 + t403) * MDP(5) + ((t367 * t406 - t423) * pkin(1) + t432) * MDP(6) + (t310 * t244 - (t395 * qJD(5) + t270 * t364 + t271 * t359) * t353 - t394 * t351 - t284 * t385 + t434 * t277 + t388) * MDP(23) + (t284 * t344 + t314 + (-pkin(2) * t343 - t280) * t361 + ((-pkin(2) * t344 - t303) * qJD(3) - t380) * t366 + t464) * MDP(9) + (t310 * t245 + (-t394 * qJD(5) - t270 * t359 + t271 * t364) * t353 + t395 * t351 + t284 * t465 + t434 * t275 + t387) * MDP(22) + t374 + (t283 * t344 + (t343 * t366 - t344 * t428) * pkin(2) + t377) * MDP(8); (t273 * t344 + t370) * MDP(9) + (t338 * t245 + (-t392 * qJD(5) - t294 * t359 + t295 * t364) * t353 + t393 * t351 + t401 * t275 + t273 * t465 + t387) * MDP(22) + (t379 * t365 + (t383 - t323) * t360 + t416) * MDP(16) + (t379 * t360 + (-t383 - t467) * t365 + t435) * MDP(15) + (t377 + t450) * MDP(8) + (t338 * t244 - (t393 * qJD(5) + t294 * t364 + t295 * t359) * t353 - t392 * t351 + t401 * t277 - t273 * t385 + t388) * MDP(23) + t375; MDP(12) * t445 + MDP(13) * t444 + qJDD(4) * MDP(14) + (-g(3) * t365 + t378 * t360) * MDP(15) + (g(3) * t360 + t378 * t365) * MDP(16) + (-(t255 * t359 - t452) * t353 + (-t275 * t443 + t364 * t351 - t353 * t424) * pkin(4) + t371) * MDP(22) + ((-qJD(5) * t254 - t255 * t353 - t235) * t364 + (-qJD(5) * t353 * t364 - t277 * t443 - t351 * t359) * pkin(4) + t372) * MDP(23) + t381 + (-MDP(10) * t439 + t431 * MDP(11)) * t344 ^ 2; (-t398 * t353 + t371) * MDP(22) + ((-t235 + (-qJD(5) + t353) * t254) * t364 + t372) * MDP(23) + t381;];
tau = t1;
