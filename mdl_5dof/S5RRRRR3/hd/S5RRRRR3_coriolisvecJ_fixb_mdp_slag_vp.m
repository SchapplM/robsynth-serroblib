% Calculate minimal parameter regressor of Coriolis joint torque vector for
% S5RRRRR3
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [5x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a3,a4,a5,d1,d4]';
% MDP [31x1]
%   Minimal dynamic parameter vector (fixed base model)
%   see S5RRRRR3_convert_par2_MPV_fixb.m
% 
% Output:
% tauc [5x1]
%   joint torques required to compensate Coriolis and centrifugal load

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-05 18:57
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauc = S5RRRRR3_coriolisvecJ_fixb_mdp_slag_vp(qJ, qJD, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(5,1),zeros(31,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRRRR3_coriolisvecJ_fixb_mdp_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RRRRR3_coriolisvecJ_fixb_mdp_slag_vp: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [5 1]), ...
  'S5RRRRR3_coriolisvecJ_fixb_mdp_slag_vp: pkin has to be [5x1] (double)');
assert(isreal(MDP) && all(size(MDP) == [31 1]), ...
  'S5RRRRR3_coriolisvecJ_fixb_mdp_slag_vp: MDP has to be [31x1] (double)'); 

%% Symbolic Calculation
% From coriolisvec_joint_fixb_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 18:56:38
% EndTime: 2019-12-05 18:56:48
% DurationCPUTime: 4.41s
% Computational Cost: add. (3160->335), mult. (7783->474), div. (0->0), fcn. (6150->8), ass. (0->147)
t329 = qJD(2) + qJD(3);
t334 = sin(qJ(3));
t338 = cos(qJ(2));
t428 = cos(qJ(3));
t373 = qJD(1) * t428;
t335 = sin(qJ(2));
t398 = qJD(1) * t335;
t446 = -t334 * t398 + t338 * t373;
t280 = t446 * t329;
t333 = sin(qJ(4));
t393 = qJD(4) * t333;
t445 = (-t333 * t446 + t393) * pkin(3);
t397 = qJD(1) * t338;
t306 = -t334 * t397 - t335 * t373;
t337 = cos(qJ(4));
t291 = -t306 * t333 - t337 * t329;
t336 = cos(qJ(5));
t332 = sin(qJ(5));
t354 = t306 * t337 - t329 * t333;
t415 = t354 * t332;
t257 = t336 * t291 - t415;
t301 = qJD(4) - t446;
t299 = qJD(5) + t301;
t444 = t257 * t299;
t355 = t291 * t332 + t336 * t354;
t443 = t299 * t355;
t409 = t332 * t333;
t307 = -t336 * t337 + t409;
t430 = qJD(4) + qJD(5);
t285 = t430 * t307;
t406 = t333 * t336;
t408 = t332 * t337;
t309 = t406 + t408;
t442 = t430 * t309;
t441 = -t309 * t446 + t442;
t401 = -t307 * t446 + t285;
t394 = qJD(3) * t334;
t384 = pkin(1) * t394;
t440 = t384 + t445;
t385 = pkin(1) * t397;
t273 = -pkin(2) * t446 + pkin(5) * t306 - t385;
t426 = pkin(1) * qJD(2);
t328 = t334 * t426;
t312 = pkin(5) * t329 + t328;
t262 = t273 * t333 + t312 * t337;
t391 = qJD(5) * t332;
t254 = t262 * t391;
t372 = t428 * qJD(2);
t363 = pkin(1) * t372;
t313 = -t329 * pkin(2) - t363;
t274 = t291 * pkin(3) + t313;
t439 = t257 * t274 + t254;
t310 = t334 * t338 + t428 * t335;
t288 = t329 * t310;
t281 = t288 * qJD(1);
t387 = qJD(1) * qJD(2);
t370 = t335 * t387;
t242 = pkin(1) * t370 + pkin(2) * t281 - pkin(5) * t280;
t358 = qJD(3) * t363;
t342 = -t262 * qJD(4) + t337 * t242 - t333 * t358;
t229 = t281 * pkin(3) + t342;
t228 = t336 * t229;
t435 = t337 * t273 - t312 * t333;
t230 = qJD(4) * t435 + t333 * t242 + t337 * t358;
t438 = -t332 * t230 + t274 * t355 + t228;
t437 = t281 * MDP(29) + (-t257 ^ 2 + t355 ^ 2) * MDP(26) - t257 * MDP(25) * t355;
t436 = (t335 ^ 2 - t338 ^ 2) * MDP(5);
t362 = qJD(3) * t328;
t392 = qJD(4) * t337;
t434 = t313 * t392 + t333 * t362;
t390 = qJD(5) * t333;
t433 = (-t390 - t393) * t332;
t353 = -t334 * t335 + t428 * t338;
t287 = t329 * t353;
t395 = qJD(2) * t335;
t250 = pkin(1) * t395 + pkin(2) * t288 - pkin(5) * t287;
t283 = -pkin(1) * t338 - pkin(2) * t353 - pkin(5) * t310;
t432 = t337 * t250 - t283 * t393;
t431 = -qJD(5) * t337 - t392;
t252 = -t354 * qJD(4) + t280 * t333;
t251 = t337 * t280 + t306 * t393 + t329 * t392;
t368 = t251 * t332 + t336 * t252;
t227 = -t355 * qJD(5) + t368;
t403 = t337 * t283;
t265 = -pkin(3) * t353 + t403;
t243 = pkin(3) * t301 + t435;
t364 = qJD(5) * t243 + t230;
t407 = t333 * t281;
t429 = -(qJD(5) * t265 + t250 * t333 + t283 * t392) * t299 + t364 * t353 - t283 * t407;
t427 = t337 * pkin(3);
t425 = pkin(1) * qJD(3);
t424 = t243 * t336;
t423 = t251 * t333;
t422 = t262 * t336;
t420 = t281 * t307;
t419 = t281 * t309;
t418 = t287 * t333;
t417 = t291 * t301;
t416 = t354 * t301;
t366 = t301 * t337;
t414 = t446 * t313;
t413 = t310 * t333;
t412 = t310 * t337;
t410 = t313 * t333;
t404 = t337 * t281;
t389 = qJD(5) * t336;
t381 = t336 * t251 - t332 * t252 - t291 * t389;
t380 = t333 * t428;
t379 = t337 * t428;
t378 = t428 * t329;
t371 = t428 * qJD(3);
t367 = t306 * t435 + t313 * t393;
t282 = -pkin(2) * t306 - pkin(5) * t446;
t272 = pkin(1) * t398 + t282;
t325 = pkin(1) * t334 + pkin(5);
t365 = qJD(4) * t325 + t272;
t361 = t301 * t380;
t326 = -t428 * pkin(1) - pkin(2);
t360 = -t328 + t445;
t359 = -t262 * t306 + t434;
t235 = t243 * t332 + t422;
t356 = -t281 * t325 - t414;
t352 = t310 * t392 + t418;
t351 = t287 * t337 - t310 * t393;
t226 = t354 * t391 + t381;
t234 = -t262 * t332 + t424;
t244 = pkin(3) * t252 + t362;
t348 = t234 * t306 + t244 * t307 + t441 * t274;
t347 = -t235 * t306 + t244 * t309 - t401 * t274;
t344 = -t265 * t281 - (pkin(3) * t288 - t283 * t390 + t432) * t299;
t341 = (-t226 * t307 - t227 * t309 + t401 * t257 + t355 * t441) * MDP(26) + (t226 * t309 + t355 * t401) * MDP(25) + ((t251 - t417) * t337 + (-t252 + t416) * t333) * MDP(19) + (-t401 * t299 - t306 * t355 + t419) * MDP(27) + (-t257 * t306 - t299 * t441 - t420) * MDP(28) + (-t354 * t366 + t423) * MDP(18) + (-t301 ^ 2 * t333 - t291 * t306 + t404) * MDP(21) + (t301 * t366 - t306 * t354 + t407) * MDP(20) + (-t306 * t329 - t281) * MDP(14) + (t306 ^ 2 - t446 ^ 2) * MDP(12) + (MDP(11) * t446 + t301 * MDP(22) + t299 * MDP(29)) * t306;
t339 = qJD(2) ^ 2;
t327 = -pkin(2) - t427;
t314 = t326 - t427;
t300 = t306 * pkin(3);
t296 = t306 * t385;
t295 = t446 * t385;
t278 = t307 * t310;
t277 = t309 * t310;
t276 = t337 * t282;
t268 = t333 * t282 + t337 * t363;
t267 = t281 * t353;
t255 = t272 * t337 - t300;
t253 = -t333 * t363 + t276 - t300;
t232 = t287 * t408 + (t430 * t412 + t418) * t336 + t433 * t310;
t231 = -t307 * t287 - t310 * t442;
t1 = [-0.2e1 * t387 * t436 + (t280 * t310 - t287 * t306) * MDP(11) + (t280 * t353 - t281 * t310 + t287 * t446 + t288 * t306) * MDP(12) + (t251 * t412 - t351 * t354) * MDP(18) + ((-t291 * t337 + t333 * t354) * t287 + (-t423 - t252 * t337 + (t291 * t333 + t337 * t354) * qJD(4)) * t310) * MDP(19) + (-t251 * t353 - t288 * t354 + t351 * t301 + t310 * t404) * MDP(20) + (t252 * t353 - t288 * t291 - t352 * t301 - t310 * t407) * MDP(21) + (t288 * t301 - t267) * MDP(22) + (t281 * t403 + t287 * t410 + t288 * t435 + t432 * t301 + t434 * t310 - t342 * t353) * MDP(23) + (t230 * t353 - t262 * t288 + (-qJD(4) * t283 * t301 + t313 * t287 + t310 * t362) * t337 + (-qJD(4) * t310 * t313 - t250 * t301 - t281 * t283) * t333) * MDP(24) + (-t226 * t278 - t231 * t355) * MDP(25) + (-t226 * t277 + t227 * t278 - t231 * t257 + t232 * t355) * MDP(26) + (-t226 * t353 + t231 * t299 - t278 * t281 - t288 * t355) * MDP(27) + (t227 * t353 - t232 * t299 - t257 * t288 - t277 * t281) * MDP(28) + (t288 * t299 - t267) * MDP(29) + (-t228 * t353 + t274 * t232 + t234 * t288 + t244 * t277 + (qJD(5) * t262 * t353 - t344) * t336 + t429 * t332 + (t227 * t413 + t352 * t257) * pkin(3)) * MDP(30) + (t274 * t231 - t235 * t288 - t244 * t278 - t254 * t353 + (t229 * t353 + t344) * t332 + t429 * t336 + (t226 * t413 - t352 * t355) * pkin(3)) * MDP(31) - t339 * t335 * MDP(7) + (0.2e1 * MDP(4) * t370 + t339 * MDP(6)) * t338 + (t287 * MDP(13) - t288 * MDP(14)) * t329 + ((-t446 * t395 - t281 * t338 + (-t288 * t338 - t353 * t395) * qJD(1)) * MDP(16) + (-t306 * t395 - t280 * t338 + (-t287 * t338 + t310 * t395) * qJD(1)) * MDP(17)) * pkin(1); (t314 * t226 + t325 * t420 + t347 + (-(-t332 * t380 + t336 * t379) * t425 + t442 * t325 + t255 * t332 + t272 * t406) * t299 - t440 * t355) * MDP(31) + (t314 * t227 - t325 * t419 + t348 + ((-t332 * t379 - t336 * t380) * t425 + t325 * t285 - t255 * t336 + t272 * t409) * t299 + t440 * t257) * MDP(30) + t341 + (t326 * t252 + t356 * t333 - t365 * t366 + (-t361 + (-qJD(2) * t337 + t291) * t334) * t425 + t367) * MDP(23) + (-t354 * t384 + t326 * t251 + t356 * t337 + (-pkin(1) * t337 * t371 + t365 * t333) * t301 + t359) * MDP(24) + (t295 + (t306 * t398 + (-t372 - t378) * qJD(3)) * pkin(1)) * MDP(17) + (-t296 + (t446 * t398 + (-qJD(2) - t329) * t394) * pkin(1)) * MDP(16) + (-t335 * t338 * MDP(4) + t436) * qJD(1) ^ 2; t341 + (t327 * t226 + (t253 * t332 + t268 * t336) * t299 - t360 * t355 + (-(t431 * t332 - t333 * t389 - t336 * t393) * t299 + t420) * pkin(5) + t347) * MDP(31) + (-t446 * t410 - pkin(2) * t252 - t276 * t301 + (-t301 * t392 - t407) * pkin(5) + (t361 + (-qJD(3) * t337 - t291) * t334) * t426 + t367) * MDP(23) + (t327 * t227 - (t253 * t336 - t268 * t332) * t299 + t360 * t257 + ((t431 * t336 - t433) * t299 - t419) * pkin(5) + t348) * MDP(30) + (t354 * t328 - t337 * t414 - pkin(2) * t251 + t268 * t301 + (t301 * t393 - t404) * pkin(5) + t359) * MDP(24) + (t295 + (-t371 + t378) * t426) * MDP(17) + (-t296 + (-qJD(3) + t329) * t328) * MDP(16); -t354 * t291 * MDP(18) + (-t291 ^ 2 + t354 ^ 2) * MDP(19) + (t251 + t417) * MDP(20) + (-t252 - t416) * MDP(21) + t281 * MDP(22) + (t262 * t301 + t313 * t354 + t342) * MDP(23) + (t313 * t291 + t301 * t435 - t230) * MDP(24) + (t226 + t444) * MDP(27) + (-t227 - t443) * MDP(28) + (-(-t332 * t435 - t422) * t299 - t235 * qJD(5) + (t257 * t354 + t336 * t281 - t299 * t391) * pkin(3) + t438) * MDP(30) + ((-t262 * t299 - t229) * t332 + (t299 * t435 - t364) * t336 + (-t332 * t281 - t299 * t389 - t354 * t355) * pkin(3) + t439) * MDP(31) + t437; (t381 + t444) * MDP(27) + (-t368 - t443) * MDP(28) + (t235 * t299 + t438) * MDP(30) + (-t332 * t229 - t336 * t230 + t234 * t299 + t439) * MDP(31) + (MDP(27) * t415 + t355 * MDP(28) - t235 * MDP(30) - MDP(31) * t424) * qJD(5) + t437;];
tauc = t1;
