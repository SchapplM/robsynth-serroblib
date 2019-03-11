% Calculate minimal parameter regressor of Coriolis joint torque vector for
% S6RPRPPR5
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d3,d6,theta2,theta5]';
% MDP [29x1]
%   Minimal dynamic parameter vector (fixed base model)
%   see S6RPRPPR5_convert_par2_MPV_fixb.m
% 
% Output:
% tauc [6x1]
%   joint torques required to compensate Coriolis and centrifugal load

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 02:52
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauc = S6RPRPPR5_coriolisvecJ_fixb_mdp_slag_vp(qJ, qJD, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(10,1),zeros(29,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRPPR5_coriolisvecJ_fixb_mdp_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RPRPPR5_coriolisvecJ_fixb_mdp_slag_vp: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RPRPPR5_coriolisvecJ_fixb_mdp_slag_vp: pkin has to be [10x1] (double)');
assert(isreal(MDP) && all(size(MDP) == [29 1]), ...
  'S6RPRPPR5_coriolisvecJ_fixb_mdp_slag_vp: MDP has to be [29x1] (double)'); 

%% Symbolic Calculation
% From coriolisvec_joint_fixb_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 02:51:47
% EndTime: 2019-03-09 02:51:56
% DurationCPUTime: 5.66s
% Computational Cost: add. (3397->392), mult. (8794->509), div. (0->0), fcn. (6638->8), ass. (0->164)
t405 = sin(pkin(9));
t407 = cos(pkin(9));
t410 = sin(qJ(3));
t489 = cos(qJ(3));
t380 = t405 * t489 + t410 * t407;
t496 = t380 * qJD(1);
t358 = qJD(6) + t496;
t448 = t489 * t407;
t439 = qJD(1) * t448;
t392 = qJD(3) * t439;
t458 = qJD(3) * t410;
t446 = t405 * t458;
t351 = qJD(1) * t446 - t392;
t404 = sin(pkin(10));
t406 = cos(pkin(10));
t409 = sin(qJ(6));
t411 = cos(qJ(6));
t377 = t404 * t411 + t406 * t409;
t464 = t358 * t377;
t497 = -t404 * t409 + t406 * t411;
t499 = -t497 * t351 - t358 * t464;
t469 = t405 * t410;
t378 = -t448 + t469;
t323 = t497 * t378;
t456 = qJD(6) * t411;
t457 = qJD(6) * t409;
t463 = -t404 * t457 + t406 * t456 + t497 * t496;
t436 = t377 * t351 - t358 * t463;
t447 = qJD(1) * t469;
t365 = -t439 + t447;
t340 = qJD(3) * t404 - t406 * t365;
t342 = qJD(3) * t406 + t365 * t404;
t429 = t340 * t409 - t342 * t411;
t498 = t358 * t429;
t485 = pkin(7) + qJ(2);
t386 = t485 * t405;
t381 = qJD(1) * t386;
t387 = t485 * t407;
t382 = qJD(1) * t387;
t465 = -t489 * t381 - t410 * t382;
t495 = qJD(4) - t465;
t339 = t410 * t386 - t489 * t387;
t494 = (t405 ^ 2 + t407 ^ 2) * (MDP(7) * qJ(2) + MDP(6));
t332 = -t410 * t381 + t489 * t382;
t313 = -pkin(4) * t365 + t332;
t403 = qJD(3) * qJ(4);
t304 = qJD(5) + t313 + t403;
t372 = t380 * qJD(3);
t352 = qJD(1) * t372;
t481 = qJ(4) * t352;
t486 = pkin(3) + qJ(5);
t493 = -t351 * t486 + (qJD(5) - t304) * t496 + t481;
t490 = t496 ^ 2;
t492 = t351 * t404 - t406 * t490;
t444 = qJD(3) * t489;
t445 = qJD(2) * t489;
t317 = (qJD(2) * t405 + qJD(3) * t387) * t410 + t386 * t444 - t407 * t445;
t491 = t365 ^ 2;
t488 = pkin(3) * t352;
t487 = pkin(8) * t406;
t484 = -pkin(8) - t486;
t482 = qJ(4) * t351;
t480 = qJ(4) * t365;
t475 = t342 * t409;
t299 = t411 * t340 + t475;
t479 = t299 * t358;
t478 = t299 * t365;
t477 = t429 * t365;
t399 = -pkin(2) * t407 - pkin(1);
t385 = qJD(1) * t399 + qJD(2);
t414 = -qJ(4) * t496 + t385;
t314 = pkin(3) * t365 + t414;
t476 = t314 * t496;
t473 = t352 * t406;
t472 = t365 * t496;
t426 = -qJD(4) * t496 + t482;
t276 = qJD(5) * t365 + t352 * t486 + t426;
t438 = qJD(1) * t445;
t451 = qJD(1) * qJD(2);
t443 = t410 * t451;
t309 = -t381 * t458 + t382 * t444 + t405 * t438 + t407 * t443;
t286 = -pkin(4) * t351 - qJD(3) * qJD(5) + t309;
t256 = t406 * t276 + t404 * t286;
t371 = -t407 * t444 + t446;
t425 = qJ(4) * t371 - qJD(4) * t380;
t280 = qJD(5) * t378 + t372 * t486 + t425;
t318 = t380 * qJD(2) - t339 * qJD(3);
t293 = -t371 * pkin(4) + t318;
t261 = t406 * t280 + t404 * t293;
t291 = t365 * t486 + t414;
t453 = pkin(4) * t496 + t495;
t298 = -qJD(3) * t486 + t453;
t266 = t406 * t291 + t404 * t298;
t308 = t486 * t496 + t480;
t269 = t406 * t308 + t404 * t313;
t422 = -qJ(4) * t380 + t399;
t311 = t378 * t486 + t422;
t338 = t386 * t489 + t410 * t387;
t321 = t380 * pkin(4) + t338;
t274 = t406 * t311 + t404 * t321;
t460 = qJD(3) * t317;
t459 = qJD(3) * t318;
t449 = -pkin(5) * t406 - pkin(4);
t454 = -t449 * t496 + t495;
t450 = -t340 * t456 + t352 * t377;
t284 = t406 * t286;
t252 = -pkin(5) * t351 + t284 + (-pkin(8) * t352 - t276) * t404;
t253 = pkin(8) * t473 + t256;
t441 = t411 * t252 - t253 * t409;
t265 = -t291 * t404 + t406 * t298;
t440 = -t381 * t444 - t382 * t458 - t405 * t443 + t407 * t438;
t437 = t497 * t352;
t402 = qJD(3) * qJD(4);
t303 = -t402 - t440;
t434 = t252 * t409 + t253 * t411;
t255 = -t276 * t404 + t284;
t433 = t255 * t406 + t256 * t404;
t258 = pkin(5) * t496 - pkin(8) * t342 + t265;
t259 = -pkin(8) * t340 + t266;
t249 = t258 * t411 - t259 * t409;
t250 = t258 * t409 + t259 * t411;
t316 = t406 * t321;
t263 = pkin(5) * t380 + t316 + (-pkin(8) * t378 - t311) * t404;
t267 = t378 * t487 + t274;
t432 = t263 * t411 - t267 * t409;
t431 = t263 * t409 + t267 * t411;
t430 = -t265 * t406 - t266 * t404;
t428 = -t351 * t406 - t404 * t490;
t307 = t406 * t313;
t383 = t484 * t404;
t421 = qJD(5) * t406 + qJD(6) * t383 - pkin(5) * t365 + t307 + (-pkin(8) * t496 - t308) * t404;
t384 = t484 * t406;
t420 = qJD(5) * t404 - qJD(6) * t384 + t487 * t496 + t269;
t324 = t377 * t378;
t270 = -t342 * t457 + t450;
t418 = qJD(3) * t465 - t440;
t417 = qJD(3) * t332 - t309;
t285 = -pkin(4) * t352 - t303;
t322 = -pkin(4) * t378 - t339;
t416 = t285 * t378 + t304 * t372 + t322 * t352;
t271 = -qJD(6) * t429 - t437;
t397 = pkin(5) * t404 + qJ(4);
t353 = qJD(3) * t365;
t333 = t351 * t380;
t328 = pkin(3) * t378 + t422;
t327 = pkin(3) * t496 + t480;
t326 = -t403 - t332;
t325 = -qJD(3) * pkin(3) + t495;
t310 = pkin(3) * t372 + t425;
t297 = t426 + t488;
t296 = t378 * t449 - t339;
t292 = -pkin(4) * t372 - t317;
t290 = t406 * t293;
t282 = pkin(5) * t340 + t304;
t281 = t372 * t449 - t317;
t279 = qJD(6) * t324 - t497 * t372;
t278 = qJD(6) * t323 + t372 * t377;
t275 = t352 * t449 - t303;
t273 = -t311 * t404 + t316;
t268 = -t308 * t404 + t307;
t260 = -t280 * t404 + t290;
t257 = t372 * t487 + t261;
t254 = -pkin(5) * t371 + t290 + (-pkin(8) * t372 - t280) * t404;
t1 = [(-t371 * t496 - t333) * MDP(8) + (t351 * t378 - t352 * t380 + t365 * t371 - t372 * t496) * MDP(9) + (t352 * t399 + t372 * t385 - t459) * MDP(13) + (-t351 * t399 - t371 * t385 + t460) * MDP(14) + (t303 * t378 + t309 * t380 + t317 * t365 + t318 * t496 - t325 * t371 + t326 * t372 - t338 * t351 + t339 * t352) * MDP(15) + (-t297 * t378 - t310 * t365 - t314 * t372 - t328 * t352 + t459) * MDP(16) + (-t297 * t380 - t310 * t496 + t314 * t371 + t328 * t351 - t460) * MDP(17) + (t297 * t328 + t303 * t339 + t309 * t338 + t310 * t314 + t317 * t326 + t318 * t325) * MDP(18) + (t255 * t380 + t260 * t496 - t265 * t371 - t273 * t351 + t292 * t340 - t416 * t406) * MDP(19) + (-t256 * t380 - t261 * t496 + t266 * t371 + t274 * t351 + t292 * t342 + t404 * t416) * MDP(20) + (-t260 * t342 - t261 * t340 + (t256 * t378 + t266 * t372 + t274 * t352) * t406 + (-t255 * t378 - t265 * t372 - t273 * t352) * t404) * MDP(21) + (t255 * t273 + t256 * t274 + t260 * t265 + t261 * t266 + t285 * t322 + t292 * t304) * MDP(22) + (t270 * t324 - t278 * t429) * MDP(23) + (t270 * t323 - t271 * t324 - t278 * t299 + t279 * t429) * MDP(24) + (t270 * t380 + t278 * t358 - t324 * t351 + t371 * t429) * MDP(25) + (-t271 * t380 - t279 * t358 + t299 * t371 - t323 * t351) * MDP(26) + (-t358 * t371 - t333) * MDP(27) + ((t254 * t411 - t257 * t409) * t358 - t432 * t351 + t441 * t380 - t249 * t371 + t281 * t299 + t296 * t271 - t275 * t323 + t282 * t279 + (-t250 * t380 - t358 * t431) * qJD(6)) * MDP(28) + (-(t254 * t409 + t257 * t411) * t358 + t431 * t351 - t434 * t380 + t250 * t371 - t281 * t429 + t296 * t270 + t275 * t324 + t282 * t278 + (-t249 * t380 - t358 * t432) * qJD(6)) * MDP(29) + 0.2e1 * t451 * t494 + (-MDP(10) * t371 - MDP(11) * t372) * qJD(3); t392 * MDP(14) + (-t490 - t491) * MDP(15) + (t351 + t353) * MDP(17) + (t488 + t482 - t326 * t365 + (-qJD(4) - t325) * t496) * MDP(18) + (t340 * t365 + t492) * MDP(19) + (t342 * t365 - t428) * MDP(20) + ((t340 * t404 + t342 * t406) * t496 + (t404 ^ 2 + t406 ^ 2) * t352) * MDP(21) + (-t255 * t404 + t256 * t406 + t304 * t365 + t430 * t496) * MDP(22) + (t436 + t478) * MDP(28) + (-t477 - t499) * MDP(29) - qJD(1) ^ 2 * t494 + ((-t365 - t447) * MDP(14) + (0.2e1 * MDP(13) - 0.2e1 * MDP(16)) * t496) * qJD(3); MDP(8) * t472 + (t490 - t491) * MDP(9) + (t392 + (t365 - t447) * qJD(3)) * MDP(10) + (-t385 * t496 + t417) * MDP(13) + (t365 * t385 + t418) * MDP(14) + (pkin(3) * t351 - t481 + (-t326 - t332) * t496 + (t325 - t495) * t365) * MDP(15) + (t327 * t365 - t417 + t476) * MDP(16) + (-t314 * t365 + t327 * t496 + 0.2e1 * t402 - t418) * MDP(17) + (-pkin(3) * t309 - qJ(4) * t303 - t314 * t327 - t325 * t332 - t326 * t495) * MDP(18) + (t265 * t365 - t268 * t496 + t285 * t404 + t453 * t340 - t493 * t406) * MDP(19) + (-t266 * t365 + t269 * t496 + t285 * t406 + t453 * t342 + t493 * t404) * MDP(20) + (t268 * t342 + t269 * t340 + (qJD(5) * t342 - t266 * t496 - t255) * t406 + (qJD(5) * t340 + t265 * t496 - t256) * t404) * MDP(21) + (qJ(4) * t285 + qJD(5) * t430 - t265 * t268 - t266 * t269 + t304 * t453 - t433 * t486) * MDP(22) + (t270 * t497 + t429 * t464) * MDP(23) + (-t270 * t377 - t271 * t497 + t299 * t464 + t429 * t463) * MDP(24) + (-t477 + t499) * MDP(25) + (t436 - t478) * MDP(26) + t358 * t365 * MDP(27) + (-(-t383 * t409 + t384 * t411) * t351 + t397 * t271 + t275 * t377 + t249 * t365 + (t409 * t420 - t411 * t421) * t358 + t454 * t299 + t463 * t282) * MDP(28) + ((t383 * t411 + t384 * t409) * t351 + t397 * t270 + t275 * t497 - t250 * t365 + (t409 * t421 + t411 * t420) * t358 - t454 * t429 - t464 * t282) * MDP(29); (-t351 + t353) * MDP(15) - MDP(16) * t472 + (-qJD(3) ^ 2 - t490) * MDP(17) + (qJD(3) * t326 + t309 + t476) * MDP(18) + (-qJD(3) * t340 + t428) * MDP(19) + (-qJD(3) * t342 + t492) * MDP(20) + (-t340 * t406 + t342 * t404) * t496 * MDP(21) + (-qJD(3) * t304 + (-t265 * t404 + t266 * t406) * t496 + t433) * MDP(22) + (-qJD(3) * t299 + t499) * MDP(28) + (qJD(3) * t429 + t436) * MDP(29); (t342 * t496 - t473) * MDP(19) + (-t340 * t496 + t352 * t404) * MDP(20) + (-t340 ^ 2 - t342 ^ 2) * MDP(21) + (t265 * t342 + t266 * t340 + t285) * MDP(22) + (t271 - t498) * MDP(28) + (t270 - t479) * MDP(29); -t429 * t299 * MDP(23) + (-t299 ^ 2 + t429 ^ 2) * MDP(24) + (t450 + t479) * MDP(25) + (t437 - t498) * MDP(26) - t351 * MDP(27) + (t250 * t358 + t282 * t429 + t441) * MDP(28) + (t249 * t358 + t282 * t299 - t434) * MDP(29) + (-MDP(25) * t475 + MDP(26) * t429 - MDP(28) * t250 - MDP(29) * t249) * qJD(6);];
tauc  = t1;
