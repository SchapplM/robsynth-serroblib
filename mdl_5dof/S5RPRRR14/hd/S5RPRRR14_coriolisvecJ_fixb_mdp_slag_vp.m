% Calculate minimal parameter regressor of Coriolis joint torque vector for
% S5RPRRR14
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [11x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,alpha2,alpha3,d1,d3,d4,d5,theta2]';
% MDP [28x1]
%   Minimal dynamic parameter vector (fixed base model)
%   see S5RPRRR14_convert_par2_MPV_fixb.m
% 
% Output:
% tauc [5x1]
%   joint torques required to compensate Coriolis and centrifugal load

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 19:22
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauc = S5RPRRR14_coriolisvecJ_fixb_mdp_slag_vp(qJ, qJD, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(11,1),zeros(28,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPRRR14_coriolisvecJ_fixb_mdp_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RPRRR14_coriolisvecJ_fixb_mdp_slag_vp: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S5RPRRR14_coriolisvecJ_fixb_mdp_slag_vp: pkin has to be [11x1] (double)');
assert(isreal(MDP) && all(size(MDP) == [28 1]), ...
  'S5RPRRR14_coriolisvecJ_fixb_mdp_slag_vp: MDP has to be [28x1] (double)'); 

%% Symbolic Calculation
% From coriolisvec_joint_fixb_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 19:19:55
% EndTime: 2019-12-31 19:20:08
% DurationCPUTime: 6.12s
% Computational Cost: add. (6431->423), mult. (21504->616), div. (0->0), fcn. (18500->12), ass. (0->184)
t395 = sin(pkin(5));
t394 = sin(pkin(6));
t501 = cos(pkin(5));
t458 = t501 * t394;
t502 = cos(qJ(3));
t427 = t502 * t458;
t397 = cos(pkin(6));
t396 = cos(pkin(11));
t463 = t502 * t396;
t448 = t397 * t463;
t518 = t395 * t448 + t427;
t393 = sin(pkin(11));
t400 = sin(qJ(3));
t487 = t397 * t400;
t419 = t393 * t487 - t463;
t414 = t395 * t419;
t365 = qJD(1) * t414;
t459 = qJD(3) * t502;
t517 = -t394 * t459 - t365;
t516 = t393 * MDP(4) + t396 * MDP(5);
t478 = qJD(1) * t395;
t462 = t393 * t478;
t515 = t518 * qJD(1) - t400 * t462;
t345 = qJD(4) - t515;
t464 = t502 * t393;
t357 = t395 * (t396 * t487 + t464) + t400 * t458;
t351 = qJD(1) * t357;
t402 = cos(qJ(4));
t455 = qJD(1) * t501;
t461 = t396 * t478;
t379 = t394 * t461;
t469 = qJD(3) - t379;
t415 = -t397 * t455 - t469;
t366 = t402 * t415;
t399 = sin(qJ(4));
t323 = t351 * t399 + t366;
t322 = qJD(5) + t323;
t512 = (t393 ^ 2 + t396 ^ 2) * MDP(6) * t395 ^ 2;
t467 = pkin(1) * t501;
t388 = t396 * t467;
t491 = t393 * t395;
t406 = t501 * pkin(2) + (-pkin(8) * t397 - qJ(2)) * t491;
t358 = t388 + t406;
t492 = t393 * t394;
t367 = (-pkin(2) * t396 - pkin(8) * t492 - pkin(1)) * t395;
t326 = -t358 * t394 + t397 * t367;
t490 = t393 * t400;
t356 = t395 * t490 - t518;
t295 = pkin(3) * t356 - pkin(9) * t357 + t326;
t457 = t501 * t397;
t488 = t395 * t396;
t371 = t394 * t488 - t457;
t412 = (t397 * t488 + t458) * pkin(8);
t480 = qJ(2) * t488 + t393 * t467;
t354 = t412 + t480;
t405 = t502 * t354 + (t358 * t397 + t367 * t394) * t400;
t299 = -t371 * pkin(9) + t405;
t511 = t399 * t295 + t402 * t299;
t444 = pkin(1) * t455;
t370 = qJ(2) * t461 + t393 * t444;
t341 = qJD(1) * t412 + t370;
t385 = t396 * t444;
t346 = qJD(1) * t406 + t385;
t361 = qJD(1) * t367 + qJD(2);
t465 = t397 * t502;
t466 = t394 * t502;
t510 = -t400 * t341 + t346 * t465 + t361 * t466;
t489 = t394 * t400;
t373 = t397 * t399 + t402 * t489;
t447 = t394 * t462;
t509 = -qJD(4) * t373 + t517 * t399 - t402 * t447;
t372 = -t397 * t402 + t399 * t489;
t508 = qJD(4) * t372 + t399 * t447 + t517 * t402;
t416 = t396 * t400 + t397 * t464;
t413 = t395 * t416;
t364 = qJD(1) * t413;
t476 = qJD(3) * t400;
t507 = -t394 * t476 + t364;
t506 = -t400 * t354 + t358 * t465 + t367 * t466;
t350 = t357 * qJD(3);
t340 = qJD(1) * t350;
t321 = -t346 * t394 + t397 * t361;
t286 = -pkin(3) * t515 - pkin(9) * t351 + t321;
t334 = t502 * t341;
t301 = t346 * t487 + t361 * t489 + t334;
t288 = -pkin(9) * t415 + t301;
t270 = t286 * t399 + t288 * t402;
t408 = qJD(2) * t414;
t283 = -qJD(1) * t408 + qJD(3) * t510;
t339 = t515 * qJD(3);
t477 = qJD(2) * t395;
t445 = t477 * t492;
t312 = pkin(3) * t340 - pkin(9) * t339 + qJD(1) * t445;
t453 = t283 * t399 - t402 * t312;
t503 = -t270 * qJD(4) - t453;
t261 = -pkin(4) * t340 - t503;
t325 = t402 * t351 - t399 * t415;
t505 = t322 * (pkin(4) * t325 + t322 * pkin(10)) + t261;
t504 = t400 * (t346 * t397 + t361 * t394) + t334;
t474 = qJD(4) * t399;
t303 = -qJD(4) * t366 + t402 * t339 - t351 * t474;
t398 = sin(qJ(5));
t401 = cos(qJ(5));
t470 = qJD(5) * t401;
t468 = t401 * t303 + t398 * t340 + t345 * t470;
t471 = qJD(5) * t398;
t274 = -t325 * t471 + t468;
t499 = t274 * t398;
t494 = t325 * t398;
t305 = -t401 * t345 + t494;
t498 = t305 * t322;
t307 = t325 * t401 + t345 * t398;
t497 = t307 * t322;
t496 = t323 * t345;
t495 = t325 * t345;
t493 = t515 * t402;
t485 = t399 * t339;
t304 = qJD(4) * t325 + t485;
t486 = t398 * t304;
t484 = t401 * t304;
t483 = t301 - t345 * (pkin(4) * t399 - pkin(10) * t402);
t320 = pkin(3) * t351 - pkin(9) * t515;
t481 = t399 * t320 + t402 * t510;
t475 = qJD(4) * t398;
t473 = qJD(4) * t401;
t472 = qJD(4) * t402;
t422 = t402 * t283 + t286 * t472 - t288 * t474 + t399 * t312;
t260 = pkin(10) * t340 + t422;
t407 = qJD(2) * t413;
t284 = qJD(1) * t407 + t504 * qJD(3);
t265 = t304 * pkin(4) - t303 * pkin(10) + t284;
t454 = -t260 * t398 + t401 * t265;
t452 = t303 * t398 - t401 * t340;
t451 = t345 * t402;
t450 = t401 * t322;
t381 = -pkin(4) * t402 - pkin(10) * t399 - pkin(3);
t449 = pkin(10) * t351 - qJD(5) * t381 + t481;
t318 = t351 * t398 + t401 * t493;
t438 = t401 * t472 - t318;
t436 = t260 * t401 + t265 * t398;
t267 = pkin(10) * t345 + t270;
t287 = pkin(3) * t415 - t510;
t273 = t323 * pkin(4) - t325 * pkin(10) + t287;
t259 = t267 * t401 + t273 * t398;
t435 = t267 * t398 - t273 * t401;
t272 = pkin(10) * t356 + t511;
t298 = t371 * pkin(3) - t506;
t327 = t357 * t399 + t371 * t402;
t328 = t357 * t402 - t371 * t399;
t278 = t327 * pkin(4) - t328 * pkin(10) + t298;
t434 = t272 * t401 + t278 * t398;
t433 = -t272 * t398 + t278 * t401;
t269 = t286 * t402 - t288 * t399;
t290 = t506 * qJD(3) - t408;
t349 = (t427 + (t448 - t490) * t395) * qJD(3);
t316 = pkin(3) * t350 - pkin(9) * t349 + t445;
t432 = -t290 * t399 + t316 * t402;
t431 = t295 * t402 - t299 * t399;
t314 = t328 * t401 + t356 * t398;
t313 = t328 * t398 - t356 * t401;
t428 = (-qJ(2) * t462 + t385) * t393 - t370 * t396;
t425 = -t322 * t470 - t486;
t424 = -t322 * t471 + t484;
t423 = -pkin(9) * t340 + t287 * t345;
t421 = t402 * t290 + t295 * t472 - t299 * t474 + t399 * t316;
t420 = -t398 * t373 - t401 * t466;
t417 = -t401 * t373 + t398 * t466;
t266 = -pkin(4) * t345 - t269;
t411 = -pkin(10) * t304 + (t266 + t269) * t322;
t291 = qJD(3) * t405 + t407;
t317 = -t401 * t351 + t398 * t493;
t311 = -qJD(4) * t327 + t349 * t402;
t310 = qJD(4) * t328 + t349 * t399;
t280 = -qJD(5) * t313 + t311 * t401 + t350 * t398;
t279 = qJD(5) * t314 + t311 * t398 - t350 * t401;
t276 = -pkin(4) * t351 - t320 * t402 + t399 * t510;
t275 = t307 * qJD(5) + t452;
t271 = -pkin(4) * t356 - t431;
t268 = t310 * pkin(4) - t311 * pkin(10) + t291;
t263 = -pkin(4) * t350 + qJD(4) * t511 - t432;
t262 = pkin(10) * t350 + t421;
t257 = -qJD(5) * t259 + t454;
t256 = -qJD(5) * t435 + t436;
t1 = [0.2e1 * qJD(2) * qJD(1) * t512 + (t339 * t357 + t349 * t351) * MDP(8) + (-t339 * t356 - t340 * t357 + t349 * t515 - t350 * t351) * MDP(9) + (-t339 * t371 - t349 * t415) * MDP(10) + (t340 * t371 + t350 * t415) * MDP(11) + (t291 * t415 + t284 * t371 + t326 * t340 + t321 * t350 + (qJD(1) * t356 - t515) * t445) * MDP(13) + (t283 * t371 + t290 * t415 + t321 * t349 + t326 * t339 + 0.2e1 * t351 * t445) * MDP(14) + (t303 * t328 + t311 * t325) * MDP(15) + (-t303 * t327 - t304 * t328 - t310 * t325 - t311 * t323) * MDP(16) + (t303 * t356 + t311 * t345 + t325 * t350 + t328 * t340) * MDP(17) + (-t304 * t356 - t310 * t345 - t323 * t350 - t327 * t340) * MDP(18) + (t340 * t356 + t345 * t350) * MDP(19) + (t432 * t345 + t431 * t340 - t453 * t356 + t269 * t350 + t291 * t323 + t298 * t304 + t284 * t327 + t287 * t310 + (-t270 * t356 - t345 * t511) * qJD(4)) * MDP(20) + (-t270 * t350 + t284 * t328 + t287 * t311 + t291 * t325 + t298 * t303 - t340 * t511 - t345 * t421 - t356 * t422) * MDP(21) + (t274 * t314 + t280 * t307) * MDP(22) + (-t274 * t313 - t275 * t314 - t279 * t307 - t280 * t305) * MDP(23) + (t274 * t327 + t280 * t322 + t304 * t314 + t307 * t310) * MDP(24) + (-t275 * t327 - t279 * t322 - t304 * t313 - t305 * t310) * MDP(25) + (t304 * t327 + t310 * t322) * MDP(26) + ((-qJD(5) * t434 - t262 * t398 + t268 * t401) * t322 + t433 * t304 + t257 * t327 - t435 * t310 + t263 * t305 + t271 * t275 + t261 * t313 + t266 * t279) * MDP(27) + (-(qJD(5) * t433 + t262 * t401 + t268 * t398) * t322 - t434 * t304 - t256 * t327 - t259 * t310 + t263 * t307 + t271 * t274 + t261 * t314 + t266 * t280) * MDP(28) + (((t396 * t480 + (qJ(2) * t491 - t388) * t393) * qJD(1) - t428) * MDP(7) - 0.2e1 * t516 * t455) * t477; t428 * MDP(7) * t478 + (t397 * t340 - t364 * t415 + (t415 * t476 + t462 * t515) * t394) * MDP(13) + (t397 * t339 + t365 * t415 + (-t351 * t462 + t415 * t459) * t394) * MDP(14) + (-t304 * t466 - t323 * t507 - t372 * t340 + t509 * t345) * MDP(20) + (-t303 * t466 - t325 * t507 - t373 * t340 + t508 * t345) * MDP(21) + (t372 * t275 + t420 * t304 + (t417 * qJD(5) + t508 * t398 - t507 * t401) * t322 - t509 * t305) * MDP(27) + (t372 * t274 + t417 * t304 + (-t420 * qJD(5) + t507 * t398 + t508 * t401) * t322 - t509 * t307) * MDP(28) + (t516 * t395 * t501 - t512) * qJD(1) ^ 2; -t515 ^ 2 * MDP(9) + (t415 * t515 + t339) * MDP(10) - t340 * MDP(11) + (-t301 * t379 + (t301 * t457 - t416 * t477) * qJD(1) + (t301 - t504) * qJD(3)) * MDP(13) + (-t510 * t379 - t321 * t515 + (t419 * t477 + t457 * t510) * qJD(1)) * MDP(14) + (t303 * t399 + t325 * t451) * MDP(15) + ((t303 - t496) * t402 + (-t304 - t495) * t399) * MDP(16) + (t399 * t340 + t345 * t451) * MDP(17) + (-t345 ^ 2 * t399 + t402 * t340) * MDP(18) + (-pkin(3) * t304 - t301 * t323 + (-t284 + (-pkin(9) * qJD(4) - t320) * t345) * t402 + (t345 * t510 + t423) * t399) * MDP(20) + (-pkin(3) * t303 + t284 * t399 - t301 * t325 + (pkin(9) * t474 + t481) * t345 + t423 * t402) * MDP(21) + (t274 * t399 * t401 + (-t399 * t471 + t438) * t307) * MDP(22) + (t305 * t318 + t307 * t317 + (-t305 * t401 - t307 * t398) * t472 + (-t499 - t275 * t401 + (t305 * t398 - t307 * t401) * qJD(5)) * t399) * MDP(23) + (-t274 * t402 + t438 * t322 + (t307 * t345 + t424) * t399) * MDP(24) + (t275 * t402 + (-t398 * t472 + t317) * t322 + (-t305 * t345 + t425) * t399) * MDP(25) + (t322 * t345 * t399 - t304 * t402) * MDP(26) + (t381 * t484 - t266 * t317 - t276 * t305 + (t398 * t449 - t401 * t483) * t322 + (t266 * t475 - t257 + (qJD(4) * t305 + t425) * pkin(9)) * t402 + (t266 * t470 + t261 * t398 - t345 * t435 + (t322 * t475 + t275) * pkin(9)) * t399) * MDP(27) + (-t381 * t486 - t266 * t318 - t276 * t307 + (t398 * t483 + t401 * t449) * t322 + (t266 * t473 + t256 + (qJD(4) * t307 - t424) * pkin(9)) * t402 + (-t266 * t471 + t261 * t401 - t345 * t259 + (t322 * t473 + t274) * pkin(9)) * t399) * MDP(28) + (-t515 * MDP(8) + (t457 * qJD(1) + t469) * MDP(11) - t321 * MDP(13) - t325 * MDP(17) + t323 * MDP(18) - t345 * MDP(19) - t269 * MDP(20) + t270 * MDP(21) + MDP(9) * t351) * t351; -t323 ^ 2 * MDP(16) + (t303 + t496) * MDP(17) + (-t485 + t495) * MDP(18) + t340 * MDP(19) + (t270 * t345 + t503) * MDP(20) + (t269 * t345 + t287 * t323 - t422) * MDP(21) + (t307 * t450 + t499) * MDP(22) + ((t274 - t498) * t401 + (-t275 - t497) * t398) * MDP(23) + (t322 * t450 + t486) * MDP(24) + (-t322 ^ 2 * t398 + t484) * MDP(25) + (-pkin(4) * t275 - t270 * t305 + t411 * t398 - t505 * t401) * MDP(27) + (-pkin(4) * t274 - t270 * t307 + t505 * t398 + t411 * t401) * MDP(28) + (MDP(15) * t323 + t325 * MDP(16) - MDP(18) * qJD(4) - t287 * MDP(20) - t307 * MDP(24) + t305 * MDP(25) - t322 * MDP(26) + MDP(27) * t435 + t259 * MDP(28)) * t325; t307 * t305 * MDP(22) + (-t305 ^ 2 + t307 ^ 2) * MDP(23) + (t468 + t498) * MDP(24) + (-t452 + t497) * MDP(25) + t304 * MDP(26) + (t259 * t322 - t266 * t307 + t454) * MDP(27) + (t266 * t305 - t322 * t435 - t436) * MDP(28) + (-MDP(24) * t494 - MDP(25) * t307 - MDP(27) * t259 + MDP(28) * t435) * qJD(5);];
tauc = t1;
