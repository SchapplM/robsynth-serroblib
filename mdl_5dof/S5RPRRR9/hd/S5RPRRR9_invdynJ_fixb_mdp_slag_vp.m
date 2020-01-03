% Calculate vector of inverse dynamics joint torques for
% S5RPRRR9
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
%   pkin=[a2,a3,a4,a5,d1,d3,d4,d5,theta2]';
% MDP [28x1]
%   Minimal dynamic parameter vector (fixed base model)
%   see S5RPRRR9_convert_par2_MPV_fixb.m
% 
% Output:
% tau [5x1]
%   joint torques of inverse dynamics (contains inertial, gravitational Coriolis and centrifugal forces)

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 19:08
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tau = S5RPRRR9_invdynJ_fixb_mdp_slag_vp(qJ, qJD, qJDD, g, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(5,1),zeros(3,1),zeros(9,1),zeros(28,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPRRR9_invdynJ_fixb_mdp_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RPRRR9_invdynJ_fixb_mdp_slag_vp: qJD has to be [5x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [5 1]), ...
  'S5RPRRR9_invdynJ_fixb_mdp_slag_vp: qJDD has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RPRRR9_invdynJ_fixb_mdp_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5RPRRR9_invdynJ_fixb_mdp_slag_vp: pkin has to be [9x1] (double)');
assert(isreal(MDP) && all(size(MDP) == [28 1]), ...
  'S5RPRRR9_invdynJ_fixb_mdp_slag_vp: MDP has to be [28x1] (double)'); 

%% Symbolic Calculation
% From invdyn_joint_fixb_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 19:08:10
% EndTime: 2019-12-31 19:08:17
% DurationCPUTime: 4.78s
% Computational Cost: add. (3758->362), mult. (9168->468), div. (0->0), fcn. (7375->14), ass. (0->169)
t391 = sin(pkin(9));
t395 = sin(qJ(3));
t392 = cos(pkin(9));
t399 = cos(qJ(3));
t463 = t399 * t392;
t350 = t391 * t395 - t463;
t342 = t350 * qJD(1);
t469 = t391 * t399;
t351 = t392 * t395 + t469;
t343 = t351 * qJD(1);
t394 = sin(qJ(4));
t398 = cos(qJ(4));
t316 = t398 * t342 + t343 * t394;
t397 = cos(qJ(5));
t452 = qJD(5) * t397;
t513 = t316 * t397 + t452;
t389 = pkin(9) + qJ(3);
t385 = qJ(4) + t389;
t376 = sin(t385);
t396 = sin(qJ(1));
t400 = cos(qJ(1));
t423 = g(1) * t400 + g(2) * t396;
t512 = t423 * t376;
t417 = -t342 * t394 + t398 * t343;
t457 = qJD(1) * t395;
t440 = t391 * t457;
t446 = qJDD(1) * t399;
t447 = qJDD(1) * t395;
t456 = qJD(3) * t399;
t441 = t391 * t446 + (qJD(1) * t456 + t447) * t392;
t323 = -qJD(3) * t440 + t441;
t345 = t351 * qJD(3);
t369 = t392 * t446;
t421 = -t391 * t447 + t369;
t324 = qJD(1) * t345 - t421;
t454 = qJD(4) * t398;
t455 = qJD(4) * t394;
t282 = t398 * t323 - t394 * t324 - t342 * t454 - t343 * t455;
t386 = qJDD(3) + qJDD(4);
t390 = qJD(3) + qJD(4);
t393 = sin(qJ(5));
t442 = t397 * t282 + t393 * t386 + t390 * t452;
t453 = qJD(5) * t393;
t270 = -t417 * t453 + t442;
t269 = t270 * t397;
t309 = t390 * t393 + t397 * t417;
t373 = t397 * t386;
t271 = t309 * qJD(5) + t282 * t393 - t373;
t307 = -t397 * t390 + t393 * t417;
t511 = -t393 * t271 - t513 * t307 + t269;
t268 = t270 * t393;
t283 = t417 * qJD(4) + t323 * t394 + t398 * t324;
t280 = qJDD(5) + t283;
t274 = t393 * t280;
t472 = t316 * t390;
t474 = t417 * t390;
t476 = t309 * t417;
t504 = qJD(5) + t316;
t510 = t386 * MDP(19) + (-t283 + t474) * MDP(18) - t316 ^ 2 * MDP(16) + (t316 * MDP(15) + MDP(16) * t417 - MDP(26) * t504) * t417 + (t282 + t472) * MDP(17) + (t513 * t309 + t268) * MDP(22) + (t513 * t504 + t274 - t476) * MDP(24);
t483 = pkin(6) + qJ(2);
t363 = t483 * t391;
t352 = qJD(1) * t363;
t364 = t483 * t392;
t353 = qJD(1) * t364;
t498 = -t399 * t352 - t353 * t395;
t305 = -pkin(7) * t343 + t498;
t304 = qJD(3) * pkin(3) + t305;
t416 = t352 * t395 - t353 * t399;
t306 = -pkin(7) * t342 - t416;
t479 = t306 * t394;
t286 = t304 * t398 - t479;
t284 = -pkin(4) * t390 - t286;
t509 = t284 * t316;
t448 = qJD(1) * qJD(2);
t489 = t483 * qJDD(1) + t448;
t331 = t489 * t391;
t332 = t489 * t392;
t434 = -t399 * t331 - t395 * t332;
t278 = qJDD(3) * pkin(3) - pkin(7) * t323 + t416 * qJD(3) + t434;
t418 = -t395 * t331 + t399 * t332;
t281 = -pkin(7) * t324 + t498 * qJD(3) + t418;
t478 = t306 * t398;
t287 = t304 * t394 + t478;
t490 = t287 * qJD(4) - t398 * t278 + t281 * t394;
t260 = -pkin(4) * t386 + t490;
t377 = cos(t385);
t484 = g(3) * t377;
t507 = t260 + t484;
t378 = -pkin(2) * t392 - pkin(1);
t357 = t378 * qJD(1) + qJD(2);
t327 = pkin(3) * t342 + t357;
t371 = g(3) * t376;
t491 = (qJD(4) * t304 + t281) * t398 + t278 * t394 - t306 * t455;
t503 = t316 * t327 + t377 * t423 + t371 - t491;
t500 = pkin(4) * t417;
t477 = t307 * t417;
t499 = (pkin(8) * t504 + t500) * t504;
t460 = -t395 * t363 + t399 * t364;
t497 = qJ(2) * qJDD(1);
t422 = g(1) * t396 - g(2) * t400;
t496 = qJDD(2) - t422;
t285 = pkin(8) * t390 + t287;
t288 = pkin(4) * t316 - pkin(8) * t417 + t327;
t263 = -t285 * t393 + t288 * t397;
t495 = -t263 * t417 + t284 * t453 + t397 * t512;
t264 = t285 * t397 + t288 * t393;
t494 = t264 * t417 + t284 * t452 + t393 * t507;
t493 = -t417 * t327 - t484 - t490 + t512;
t408 = -t363 * t456 + qJD(2) * t463 + (-qJD(2) * t391 - qJD(3) * t364) * t395;
t299 = -pkin(7) * t345 + t408;
t344 = t350 * qJD(3);
t404 = -t351 * qJD(2) - t460 * qJD(3);
t300 = pkin(7) * t344 + t404;
t433 = -t399 * t363 - t364 * t395;
t312 = -pkin(7) * t351 + t433;
t313 = -pkin(7) * t350 + t460;
t419 = t312 * t398 - t313 * t394;
t265 = t419 * qJD(4) + t299 * t398 + t300 * t394;
t293 = t312 * t394 + t313 * t398;
t325 = t398 * t350 + t351 * t394;
t326 = -t350 * t394 + t351 * t398;
t329 = pkin(3) * t350 + t378;
t294 = pkin(4) * t325 - pkin(8) * t326 + t329;
t297 = -t325 * qJD(4) - t344 * t398 - t345 * t394;
t259 = pkin(8) * t386 + t491;
t430 = qJD(5) * t288 + t259;
t488 = t260 * t326 - t293 * t280 + t284 * t297 - (qJD(5) * t294 + t265) * t504 - t430 * t325;
t487 = pkin(3) * t345;
t482 = qJDD(1) * pkin(1);
t481 = t284 * t326;
t480 = t294 * t280;
t475 = t309 * t393;
t468 = t392 * MDP(4);
t467 = t393 * t396;
t466 = t393 * t400;
t465 = t396 * t397;
t275 = t397 * t280;
t464 = t397 * t400;
t459 = t391 ^ 2 + t392 ^ 2;
t432 = t504 * t393;
t356 = t378 * qJDD(1) + qJDD(2);
t310 = pkin(3) * t324 + t356;
t262 = pkin(4) * t283 - pkin(8) * t282 + t310;
t429 = qJD(5) * t285 - t262;
t380 = pkin(3) * t394 + pkin(8);
t427 = pkin(3) * t343 + pkin(8) * t316 + qJD(5) * t380 + t500;
t426 = 0.2e1 * t459;
t289 = t305 * t394 + t478;
t425 = pkin(3) * t455 - t289;
t290 = t305 * t398 - t479;
t424 = -pkin(3) * t454 + t290;
t420 = -t280 * t380 + t509;
t415 = t482 - t496;
t414 = t275 - (t316 * t393 + t453) * t504;
t412 = t297 * t397 - t326 * t453;
t411 = -pkin(8) * t280 + t286 * t504 + t509;
t406 = t426 * t448 - t423;
t384 = cos(t389);
t383 = sin(t389);
t381 = -pkin(3) * t398 - pkin(4);
t339 = t377 * t464 + t467;
t338 = -t377 * t466 + t465;
t337 = -t377 * t465 + t466;
t336 = t377 * t467 + t464;
t298 = t326 * qJD(4) - t344 * t394 + t398 * t345;
t272 = pkin(4) * t298 - pkin(8) * t297 + t487;
t266 = t293 * qJD(4) + t299 * t394 - t300 * t398;
t261 = t397 * t262;
t1 = [(-t266 * t390 + t283 * t329 + t298 * t327 + t310 * t325 + t316 * t487 + t422 * t377 + t386 * t419) * MDP(20) + (-g(1) * t336 - g(2) * t338 - t264 * t298 + t266 * t309 - t419 * t270 + (-(-qJD(5) * t293 + t272) * t504 - t480 + t429 * t325 - qJD(5) * t481) * t393 + t488 * t397) * MDP(28) + (-g(1) * t337 - g(2) * t339 + t261 * t325 + t263 * t298 + t266 * t307 - t419 * t271 + (t272 * t504 + t480 + (-t285 * t325 - t293 * t504 + t481) * qJD(5)) * t397 + t488 * t393) * MDP(27) + (-t326 * t274 - t271 * t325 - t298 * t307 + (-t297 * t393 - t326 * t452) * t504) * MDP(25) + (t270 * t325 + t326 * t275 + t298 * t309 + t412 * t504) * MDP(24) + (t280 * t325 + t298 * t504) * MDP(26) + qJDD(1) * MDP(1) + (t404 * qJD(3) + t433 * qJDD(3) + t378 * t324 + t357 * t345 + t356 * t350 + t422 * t384) * MDP(13) + t422 * MDP(2) + t423 * MDP(3) + (-MDP(5) * t391 + t468) * (t415 + t482) + (-t265 * t390 + t282 * t329 - t293 * t386 + t297 * t327 + t310 * t326 - t422 * t376 + t417 * t487) * MDP(21) + (t282 * t326 + t297 * t417) * MDP(15) + (-t282 * t325 - t283 * t326 - t297 * t316 - t298 * t417) * MDP(16) + (-t408 * qJD(3) - t460 * qJDD(3) + t378 * t323 - t357 * t344 + t356 * t351 - t422 * t383) * MDP(14) + (-qJD(3) * t344 + qJDD(3) * t351) * MDP(10) + (t323 * t351 - t343 * t344) * MDP(8) + ((-t307 * t397 - t475) * t297 + (-t268 - t271 * t397 + (t307 * t393 - t309 * t397) * qJD(5)) * t326) * MDP(23) + (t326 * t269 + t412 * t309) * MDP(22) + (-t298 * t390 - t325 * t386) * MDP(18) + (t297 * t390 + t326 * t386) * MDP(17) + (t415 * pkin(1) + (t459 * t497 + t406) * qJ(2)) * MDP(7) + (t426 * t497 + t406) * MDP(6) + (-qJD(3) * t345 - qJDD(3) * t350) * MDP(11) + (-t323 * t350 - t324 * t351 + t342 * t344 - t343 * t345) * MDP(9); t496 * MDP(7) - t369 * MDP(13) + t441 * MDP(14) + (t283 + t474) * MDP(20) + (t282 - t472) * MDP(21) + (t414 - t477) * MDP(27) + (-t397 * t504 ^ 2 - t274 - t476) * MDP(28) + (-t468 - pkin(1) * MDP(7) + (MDP(13) * t395 + MDP(5)) * t391) * qJDD(1) + ((qJD(1) * t469 + t392 * t457 + t343) * MDP(13) + (-t342 - t440) * MDP(14)) * qJD(3) + (-MDP(7) * qJ(2) - MDP(6)) * qJD(1) ^ 2 * t459; (t381 * t271 - t507 * t397 + t420 * t393 + t425 * t307 + (t424 * t393 - t427 * t397) * t504 + t495) * MDP(27) + (t381 * t270 + t420 * t397 - t393 * t512 + t425 * t309 + (t427 * t393 + t424 * t397) * t504 + t494) * MDP(28) + (t290 * t390 + (-t343 * t417 - t386 * t394 - t390 * t454) * pkin(3) + t503) * MDP(21) + qJDD(3) * MDP(12) + t343 * t342 * MDP(8) + (-g(3) * t384 - t357 * t343 + t423 * t383 + t434) * MDP(13) + (t441 + (t342 - t440) * qJD(3)) * MDP(10) + (g(3) * t383 + t357 * t342 + t423 * t384 - t418) * MDP(14) + (-t475 * t504 + t511) * MDP(23) + (t289 * t390 + (-t316 * t343 + t386 * t398 - t390 * t455) * pkin(3) + t493) * MDP(20) + (t414 + t477) * MDP(25) + t421 * MDP(11) + (-t342 ^ 2 + t343 ^ 2) * MDP(9) + t510; (t287 * t390 + t493) * MDP(20) + (t286 * t390 + t503) * MDP(21) + (-t309 * t432 + t511) * MDP(23) + (-t432 * t504 + t275 + t477) * MDP(25) + (-pkin(4) * t271 - t287 * t307 + t411 * t393 + (-t507 - t499) * t397 + t495) * MDP(27) + (-pkin(4) * t270 - t287 * t309 + t411 * t397 + (-t512 + t499) * t393 + t494) * MDP(28) + t510; t309 * t307 * MDP(22) + (-t307 ^ 2 + t309 ^ 2) * MDP(23) + (t307 * t504 + t442) * MDP(24) + (t309 * t504 + t373) * MDP(25) + t280 * MDP(26) + (-g(1) * t338 + g(2) * t336 + t264 * t504 - t284 * t309 + t261) * MDP(27) + (g(1) * t339 - g(2) * t337 + t263 * t504 + t284 * t307) * MDP(28) + ((-t259 + t371) * MDP(28) + (-MDP(25) * t417 - MDP(27) * t285 - MDP(28) * t288) * qJD(5)) * t397 + (-qJD(5) * t417 * MDP(24) + (-qJD(5) * t390 - t282) * MDP(25) + (-t430 + t371) * MDP(27) + t429 * MDP(28)) * t393;];
tau = t1;
