% Calculate minimal parameter regressor of Coriolis joint torque vector for
% S5RRRRP6
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d2,d3,d4]';
% MDP [26x1]
%   Minimal dynamic parameter vector (fixed base model)
%   see S5RRRRP6_convert_par2_MPV_fixb.m
% 
% Output:
% tauc [5x1]
%   joint torques required to compensate Coriolis and centrifugal load

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 21:55
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauc = S5RRRRP6_coriolisvecJ_fixb_mdp_slag_vp(qJ, qJD, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(8,1),zeros(26,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRRRP6_coriolisvecJ_fixb_mdp_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RRRRP6_coriolisvecJ_fixb_mdp_slag_vp: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RRRRP6_coriolisvecJ_fixb_mdp_slag_vp: pkin has to be [8x1] (double)');
assert(isreal(MDP) && all(size(MDP) == [26 1]), ...
  'S5RRRRP6_coriolisvecJ_fixb_mdp_slag_vp: MDP has to be [26x1] (double)'); 

%% Symbolic Calculation
% From coriolisvec_joint_fixb_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 21:54:40
% EndTime: 2019-12-31 21:54:47
% DurationCPUTime: 3.05s
% Computational Cost: add. (3074->317), mult. (7523->431), div. (0->0), fcn. (5133->6), ass. (0->148)
t343 = cos(qJ(2));
t418 = -pkin(7) - pkin(6);
t325 = t418 * t343;
t319 = qJD(1) * t325;
t340 = sin(qJ(3));
t307 = t340 * t319;
t341 = sin(qJ(2));
t324 = t418 * t341;
t317 = qJD(1) * t324;
t417 = cos(qJ(3));
t280 = t417 * t317 + t307;
t377 = qJD(3) * t417;
t425 = -pkin(2) * t377 + t280;
t386 = qJD(1) * qJD(2);
t424 = -0.2e1 * t386;
t423 = (t341 ^ 2 - t343 ^ 2) * MDP(5);
t342 = cos(qJ(4));
t378 = qJD(1) * t417;
t390 = qJD(1) * t341;
t304 = t340 * t390 - t343 * t378;
t339 = sin(qJ(4));
t408 = t304 * t339;
t422 = -qJ(5) * t408 + t342 * qJD(5);
t404 = t340 * t343;
t306 = -qJD(1) * t404 - t341 * t378;
t274 = -pkin(3) * t306 + pkin(8) * t304;
t263 = pkin(2) * t390 + t274;
t421 = t339 * t263 + t425 * t342;
t420 = t417 * t324 + t340 * t325;
t385 = qJD(2) + qJD(3);
t358 = -t340 * t341 + t417 * t343;
t352 = t358 * qJD(3);
t283 = qJD(2) * t358 + t352;
t347 = t283 * qJD(1);
t353 = t342 * t306 - t339 * t385;
t249 = -qJD(4) * t353 + t339 * t347;
t419 = t353 ^ 2;
t416 = t342 * pkin(4);
t415 = -qJ(5) - pkin(8);
t414 = qJD(2) * pkin(2);
t370 = t342 * t385;
t388 = qJD(4) * t339;
t248 = -qJD(4) * t370 - t306 * t388 - t342 * t347;
t413 = t248 * t339;
t315 = t417 * t341 + t404;
t284 = t385 * t315;
t273 = t284 * qJD(1);
t412 = t273 * t339;
t411 = t273 * t342;
t287 = -t306 * t339 - t370;
t301 = qJD(4) + t304;
t410 = t287 * t301;
t409 = t353 * t301;
t407 = t304 * t342;
t406 = t315 * t339;
t405 = t315 * t342;
t344 = qJD(2) ^ 2;
t403 = t341 * t344;
t336 = t342 * qJ(5);
t292 = t340 * t324 - t417 * t325;
t285 = t342 * t292;
t402 = t343 * t344;
t345 = qJD(1) ^ 2;
t401 = t343 * t345;
t332 = t340 * pkin(2) + pkin(8);
t400 = -qJ(5) - t332;
t334 = -pkin(2) * t343 - pkin(1);
t323 = t334 * qJD(1);
t262 = pkin(3) * t304 + pkin(8) * t306 + t323;
t308 = t417 * t319;
t309 = t317 + t414;
t278 = t340 * t309 - t308;
t265 = pkin(8) * t385 + t278;
t234 = t342 * t262 - t265 * t339;
t227 = qJ(5) * t353 + t234;
t226 = pkin(4) * t301 + t227;
t399 = t226 - t227;
t277 = t417 * t309 + t307;
t398 = t339 * t274 + t342 * t277;
t372 = qJD(4) * t400;
t396 = t339 * t372 - t421 + t422;
t261 = t342 * t263;
t365 = -t306 * pkin(4) + t304 * t336;
t395 = t342 * t372 - t261 - t365 + (-qJD(5) + t425) * t339;
t276 = -pkin(3) * t358 - pkin(8) * t315 + t334;
t394 = t339 * t276 + t285;
t375 = qJD(4) * t415;
t393 = t339 * t375 - t398 + t422;
t373 = t342 * t274 - t277 * t339;
t392 = -t339 * qJD(5) + t342 * t375 - t365 - t373;
t389 = qJD(3) * t340;
t387 = qJD(4) * t342;
t384 = t341 * t414;
t246 = pkin(3) * t284 - pkin(8) * t283 + t384;
t381 = qJD(2) * t418;
t318 = t341 * t381;
t320 = t343 * t381;
t251 = qJD(3) * t420 + t417 * t318 + t340 * t320;
t382 = t339 * t246 + t342 * t251 + t276 * t387;
t379 = t315 * t387;
t264 = -pkin(3) * t385 - t277;
t259 = t264 * t387;
t376 = t341 * t386;
t374 = pkin(1) * t424;
t371 = t301 * t342;
t368 = qJD(1) * t381;
t310 = t341 * t368;
t311 = t343 * t368;
t243 = t309 * t389 + t340 * t310 - t417 * t311 - t319 * t377;
t333 = -t417 * pkin(2) - pkin(3);
t279 = t340 * t317 - t308;
t367 = pkin(2) * t389 - t279;
t366 = (t388 + t408) * pkin(4);
t235 = t262 * t339 + t265 * t342;
t364 = -t235 * t306 + t243 * t339 + t259;
t228 = -qJ(5) * t287 + t235;
t362 = -t226 * t342 - t228 * t339;
t361 = t264 * t304 - t273 * t332;
t360 = -qJ(5) * t283 - qJD(5) * t315;
t359 = t234 * t306 - t243 * t342 + t264 * t388;
t357 = t283 * t339 + t379;
t356 = t283 * t342 - t315 * t388;
t222 = pkin(4) * t249 + t243;
t355 = t323 * t306 - t243;
t239 = t273 * pkin(3) + (-pkin(8) * t352 + (t341 * pkin(2) - pkin(8) * t358) * qJD(2)) * qJD(1);
t242 = t309 * t377 + t417 * t310 + t340 * t311 + t319 * t389;
t354 = t339 * t239 + t342 * t242 + t262 * t387 - t265 * t388;
t237 = t342 * t239;
t351 = -qJD(4) * t235 - t242 * t339 + t237;
t214 = pkin(4) * t273 + qJ(5) * t248 + qJD(5) * t353 + t351;
t216 = -qJ(5) * t249 - qJD(5) * t287 + t354;
t350 = qJD(4) * t362 - t214 * t339 + t216 * t342 - t226 * t407 - t228 * t408;
t349 = ((-t248 - t410) * t342 + (-t249 + t409) * t339) * MDP(19) + (-t353 * t371 - t413) * MDP(18) + (-t301 ^ 2 * t339 - t287 * t306 + t411) * MDP(21) + (t301 * t371 - t306 * t353 + t412) * MDP(20) + t347 * MDP(13) + (-t304 ^ 2 + t306 ^ 2) * MDP(12) + (-MDP(11) * t304 + t301 * MDP(22)) * t306 + (t304 * MDP(13) + (-qJD(1) * t315 - t306) * MDP(14)) * t385;
t348 = t323 * t304 - t242;
t252 = t292 * qJD(3) + t340 * t318 - t417 * t320;
t322 = pkin(8) * t342 + t336;
t321 = t415 * t339;
t313 = t332 * t342 + t336;
t312 = t400 * t339;
t286 = t287 ^ 2;
t271 = t342 * t276;
t247 = t287 * pkin(4) + qJD(5) + t264;
t245 = t342 * t246;
t238 = -qJ(5) * t406 + t394;
t232 = -pkin(4) * t358 - t292 * t339 - t315 * t336 + t271;
t219 = -qJ(5) * t379 + (-qJD(4) * t292 + t360) * t339 + t382;
t217 = pkin(4) * t284 - t251 * t339 + t245 + t360 * t342 + (-t285 + (qJ(5) * t315 - t276) * t339) * qJD(4);
t1 = [0.2e1 * t343 * MDP(4) * t376 + t423 * t424 + MDP(6) * t402 - MDP(7) * t403 + (-pkin(6) * t402 + t341 * t374) * MDP(9) + (pkin(6) * t403 + t343 * t374) * MDP(10) + (-t306 * t283 + t315 * t347) * MDP(11) + (-t315 * t273 - t283 * t304 + t306 * t284 + t347 * t358) * MDP(12) + (t334 * t273 + t323 * t284 + (-qJD(1) * t358 + t304) * t384) * MDP(16) + (pkin(2) * t315 * t376 + t323 * t283 - t306 * t384 + t334 * t347) * MDP(17) + (-t248 * t405 - t353 * t356) * MDP(18) + ((-t287 * t342 + t339 * t353) * t283 + (t413 - t249 * t342 + (t287 * t339 + t342 * t353) * qJD(4)) * t315) * MDP(19) + (t248 * t358 + t273 * t405 - t284 * t353 + t301 * t356) * MDP(20) + (t249 * t358 - t273 * t406 - t284 * t287 - t301 * t357) * MDP(21) + (-t273 * t358 + t284 * t301) * MDP(22) + ((-t292 * t387 + t245) * t301 + t271 * t273 - (-t265 * t387 + t237) * t358 + t234 * t284 + t252 * t287 - t420 * t249 + t315 * t259 + ((-qJD(4) * t276 - t251) * t301 - t292 * t273 - (-qJD(4) * t262 - t242) * t358 + t243 * t315 + t264 * t283) * t339) * MDP(23) + (-(-t292 * t388 + t382) * t301 - t394 * t273 + t354 * t358 - t235 * t284 - t252 * t353 + t420 * t248 + t243 * t405 + t356 * t264) * MDP(24) + (t217 * t353 - t219 * t287 + t232 * t248 - t238 * t249 + t362 * t283 + (-t214 * t342 - t216 * t339 + (t226 * t339 - t228 * t342) * qJD(4)) * t315) * MDP(25) + (t216 * t238 + t228 * t219 + t214 * t232 + t226 * t217 + t222 * (pkin(4) * t406 - t420) + t247 * (pkin(4) * t357 + t252)) * MDP(26) + (t283 * MDP(13) - t284 * MDP(14) - t252 * MDP(16) - t251 * MDP(17)) * t385; t349 + (-t333 * t248 + t361 * t342 - t367 * t353 + (t332 * t388 + t421) * t301 + t364) * MDP(24) + (t280 * t385 + (t306 * t390 - t377 * t385) * pkin(2) + t348) * MDP(17) + (t248 * t312 - t249 * t313 - t396 * t287 + t353 * t395 + t350) * MDP(25) + (t216 * t313 + t214 * t312 + t222 * (t333 - t416) + (t308 + (pkin(2) * qJD(3) - t317) * t340 + t366) * t247 + t396 * t228 + t395 * t226) * MDP(26) + (t279 * t385 + (-t304 * t390 - t385 * t389) * pkin(2) + t355) * MDP(16) - t341 * MDP(4) * t401 + t345 * t423 + (t333 * t249 + t361 * t339 + t367 * t287 + (-t332 * t387 + t425 * t339 - t261) * t301 + t359) * MDP(23) + (t345 * t341 * MDP(9) + MDP(10) * t401) * pkin(1); t349 + (t277 * t385 + t348) * MDP(17) + (t278 * t385 + t355) * MDP(16) + (pkin(3) * t248 + t398 * t301 + t278 * t353 + t264 * t407 + (t301 * t388 - t411) * pkin(8) + t364) * MDP(24) + (-pkin(3) * t249 - t373 * t301 - t278 * t287 + t264 * t408 + (-t301 * t387 - t412) * pkin(8) + t359) * MDP(23) + (t216 * t322 + t214 * t321 + t222 * (-pkin(3) - t416) + (t366 - t278) * t247 + t393 * t228 + t392 * t226) * MDP(26) + (t248 * t321 - t249 * t322 - t393 * t287 + t353 * t392 + t350) * MDP(25); -t353 * t287 * MDP(18) + (-t286 + t419) * MDP(19) + (-t248 + t410) * MDP(20) + (-t249 - t409) * MDP(21) + t273 * MDP(22) + (t235 * t301 + t264 * t353 + t351) * MDP(23) + (t234 * t301 + t264 * t287 - t354) * MDP(24) + (pkin(4) * t248 - t399 * t287) * MDP(25) + (t399 * t228 + (t247 * t353 + t214) * pkin(4)) * MDP(26); (-t286 - t419) * MDP(25) + (-t226 * t353 + t228 * t287 + t222) * MDP(26);];
tauc = t1;
