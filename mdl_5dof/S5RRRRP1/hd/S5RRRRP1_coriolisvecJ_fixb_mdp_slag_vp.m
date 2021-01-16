% Calculate Coriolis joint torque vector for
% S5RRRRP1
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
% MDP [28x1]
%   Minimal dynamic parameter vector (fixed base model)
%   see S5RRRRP1_convert_par2_MPV_fixb.m
% 
% Output:
% tauc [5x1]
%   joint torques required to compensate Coriolis and centrifugal load

% Quelle: HybrDyn-Toolbox
% Datum: 2021-01-15 23:53
% Revision: 24b2e7d74a0c1a3b64fa2f8f5ad758691ad61af3 (2021-01-15)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauc = S5RRRRP1_coriolisvecJ_fixb_mdp_slag_vp(qJ, qJD, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(8,1),zeros(28,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRRRP1_coriolisvecJ_fixb_mdp_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RRRRP1_coriolisvecJ_fixb_mdp_slag_vp: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RRRRP1_coriolisvecJ_fixb_mdp_slag_vp: pkin has to be [8x1] (double)');
assert(isreal(MDP) && all(size(MDP) == [28 1]), ...
  'S5RRRRP1_coriolisvecJ_fixb_mdp_slag_vp: MDP has to be [28x1] (double)'); 

%% Symbolic Calculation
% From coriolisvec_joint_fixb_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2021-01-15 23:52:51
% EndTime: 2021-01-15 23:52:59
% DurationCPUTime: 3.08s
% Computational Cost: add. (3807->301), mult. (9961->385), div. (0->0), fcn. (7141->6), ass. (0->164)
t346 = sin(qJ(3));
t348 = cos(qJ(2));
t431 = cos(qJ(3));
t396 = qJD(1) * t431;
t347 = sin(qJ(2));
t405 = qJD(1) * t347;
t309 = -t346 * t405 + t348 * t396;
t342 = qJD(2) + qJD(3);
t288 = t309 * t342;
t416 = t346 * t348;
t310 = -qJD(1) * t416 - t347 * t396;
t305 = t310 * pkin(8);
t432 = pkin(6) + pkin(7);
t329 = t432 * t348;
t324 = qJD(1) * t329;
t311 = t346 * t324;
t328 = t432 * t347;
t322 = qJD(1) * t328;
t425 = qJD(2) * pkin(2);
t317 = -t322 + t425;
t385 = t431 * t317 - t311;
t263 = t305 + t385;
t254 = pkin(3) * t342 + t263;
t315 = t431 * t324;
t378 = -t346 * t317 - t315;
t426 = t309 * pkin(8);
t264 = -t378 + t426;
t430 = cos(qJ(4));
t258 = t430 * t264;
t345 = sin(qJ(4));
t376 = -t345 * t254 - t258;
t321 = t347 * t431 + t416;
t356 = t342 * t321;
t354 = t356 * qJD(1);
t399 = qJD(2) * t432;
t382 = qJD(1) * t399;
t318 = t347 * t382;
t319 = t348 * t382;
t394 = t431 * qJD(3);
t404 = qJD(3) * t346;
t363 = t317 * t394 - t318 * t431 - t346 * t319 - t324 * t404;
t239 = -pkin(8) * t354 + t363;
t384 = t346 * t318 - t431 * t319;
t367 = qJD(3) * t378 + t384;
t240 = -t288 * pkin(8) + t367;
t390 = -t345 * t239 + t430 * t240;
t369 = qJD(4) * t376 + t390;
t395 = qJD(4) * t430;
t403 = qJD(4) * t345;
t371 = t430 * t288 + t309 * t395 + t310 * t403 - t345 * t354;
t423 = t371 * qJ(5);
t361 = t369 - t423;
t374 = -t345 * t309 + t310 * t430;
t401 = t374 * qJD(5);
t219 = t361 + t401;
t338 = -t348 * pkin(2) - pkin(1);
t327 = t338 * qJD(1);
t295 = -pkin(3) * t309 + t327;
t282 = t430 * t309 + t310 * t345;
t392 = -pkin(4) * t282 + qJD(5);
t248 = t295 + t392;
t422 = t248 * t374;
t439 = t219 + t422;
t419 = t295 * t374;
t438 = t369 + t419;
t278 = t282 ^ 2;
t352 = t430 * t356;
t368 = -qJD(1) * t352 + qJD(4) * t374 - t345 * t288;
t341 = qJD(4) + t342;
t420 = t374 * t341;
t421 = t282 * t341;
t433 = t374 ^ 2;
t437 = t282 * t374 * MDP(18) + (-t278 + t433) * MDP(19) + (t368 - t420) * MDP(21) + (t371 - t421) * MDP(20);
t424 = qJ(5) * t282;
t272 = t374 * qJ(5);
t364 = t239 * t430 + t345 * t240 + t254 * t395 - t264 * t403;
t360 = -t295 * t282 - t364;
t400 = qJD(2) * qJD(1);
t435 = -0.2e1 * t400;
t434 = (t347 ^ 2 - t348 ^ 2) * MDP(5);
t393 = t347 * t400;
t335 = pkin(2) * t393;
t268 = pkin(3) * t354 + t335;
t222 = -pkin(4) * t368 + t268;
t429 = pkin(3) * t310;
t428 = pkin(3) * t341;
t418 = t327 * t310;
t256 = t345 * t264;
t417 = t345 * t346;
t349 = qJD(2) ^ 2;
t415 = t347 * t349;
t414 = t348 * t349;
t350 = qJD(1) ^ 2;
t413 = t348 * t350;
t389 = t430 * t254 - t256;
t227 = t389 + t272;
t225 = pkin(4) * t341 + t227;
t412 = t225 - t227;
t383 = t322 * t346 - t315;
t266 = t383 - t426;
t407 = -t431 * t322 - t311;
t267 = t305 + t407;
t387 = t430 * t266 - t267 * t345;
t231 = t387 - t424;
t337 = pkin(2) * t431 + pkin(3);
t398 = t430 * t346;
t287 = -t337 * t403 + (-t346 * t395 + (-t345 * t431 - t398) * qJD(3)) * pkin(2);
t411 = t231 - t287;
t408 = t345 * t266 + t430 * t267;
t232 = t272 + t408;
t286 = t337 * t395 + (-t346 * t403 + (t430 * t431 - t417) * qJD(3)) * pkin(2);
t410 = t232 - t286;
t409 = t430 * t263 - t256;
t340 = t347 * t425;
t339 = pkin(2) * t405;
t391 = pkin(1) * t435;
t388 = -t263 * t345 - t258;
t251 = -pkin(4) * t374 - t429;
t228 = -t376 + t424;
t380 = t225 * t282 - t228 * t374;
t379 = t346 * t347 - t348 * t431;
t377 = t346 * t328 - t329 * t431;
t273 = -t321 * pkin(8) - t328 * t431 - t346 * t329;
t274 = -pkin(8) * t379 - t377;
t375 = -t345 * t273 - t274 * t430;
t323 = t347 * t399;
t325 = t348 * t399;
t373 = -t431 * t323 - t346 * t325 - t328 * t394 - t329 * t404;
t246 = -pkin(8) * t356 + t373;
t294 = t342 * t379;
t365 = qJD(3) * t377 + t346 * t323 - t431 * t325;
t247 = t294 * pkin(8) + t365;
t372 = t430 * t246 + t345 * t247 + t273 * t395 - t274 * t403;
t370 = t430 * t379;
t298 = pkin(3) * t379 + t338;
t293 = t321 * t430 - t345 * t379;
t366 = qJD(4) * t375 - t345 * t246 + t430 * t247;
t362 = t310 * t309 * MDP(11) + (-t310 * t342 - t354) * MDP(14) + (-t309 ^ 2 + t310 ^ 2) * MDP(12) + t437;
t359 = -t327 * t309 - t363;
t358 = qJ(5) * t368 + t364;
t357 = (-t258 + (-t254 - t428) * t345) * qJD(4) + t390;
t218 = qJD(5) * t282 + t358;
t285 = pkin(3) * t356 + t340;
t351 = -t248 * t282 - t218;
t336 = pkin(3) * t430 + pkin(4);
t326 = t395 * t428;
t306 = pkin(2) * t398 + t345 * t337;
t303 = -pkin(2) * t417 + t337 * t430 + pkin(4);
t296 = t339 - t429;
t292 = t321 * t345 + t370;
t271 = t287 * t341;
t270 = t286 * t341;
t259 = t292 * pkin(4) + t298;
t249 = t251 + t339;
t242 = qJD(4) * t293 - t345 * t294 + t352;
t241 = qJD(4) * t370 + t294 * t430 + t321 * t403 + t345 * t356;
t235 = -t292 * qJ(5) - t375;
t234 = -t293 * qJ(5) + t273 * t430 - t345 * t274;
t233 = t242 * pkin(4) + t285;
t230 = t272 + t409;
t229 = t388 - t424;
t221 = t241 * qJ(5) - t293 * qJD(5) + t366;
t220 = -qJ(5) * t242 - qJD(5) * t292 + t372;
t1 = [0.2e1 * t348 * MDP(4) * t393 + t434 * t435 + MDP(6) * t414 - MDP(7) * t415 + (-pkin(6) * t414 + t347 * t391) * MDP(9) + (pkin(6) * t415 + t348 * t391) * MDP(10) + (t288 * t321 + t294 * t310) * MDP(11) + (-t288 * t379 - t294 * t309 + t310 * t356 - t321 * t354) * MDP(12) + (-t309 * t340 + 0.2e1 * t327 * t356 + t335 * t379) * MDP(16) + (t338 * t288 - t327 * t294 + (qJD(1) * t321 - t310) * t340) * MDP(17) + (t241 * t374 + t293 * t371) * MDP(18) + (-t241 * t282 + t242 * t374 - t292 * t371 + t293 * t368) * MDP(19) + (t295 * t242 + t268 * t292 - t282 * t285 - t298 * t368) * MDP(23) + (-t295 * t241 + t268 * t293 - t285 * t374 + t298 * t371) * MDP(24) + (t222 * t292 - t233 * t282 + t242 * t248 - t259 * t368) * MDP(25) + (t222 * t293 - t233 * t374 - t241 * t248 + t259 * t371) * MDP(26) + (-t218 * t292 - t219 * t293 + t220 * t282 + t221 * t374 + t225 * t241 - t228 * t242 - t234 * t371 + t235 * t368) * MDP(27) + (t218 * t235 + t219 * t234 + t220 * t228 + t221 * t225 + t222 * t259 + t233 * t248) * MDP(28) + (-t294 * MDP(13) - t356 * MDP(14) + t365 * MDP(16) - t373 * MDP(17)) * t342 + (-t241 * MDP(20) - t242 * MDP(21) + t366 * MDP(23) - t372 * MDP(24) + t221 * MDP(25) - t220 * MDP(26)) * t341; (t309 * t339 + t418 - t383 * t342 + (-t315 + (-pkin(2) * t342 - t317) * t346) * qJD(3) + t384) * MDP(16) + (-t282 * t410 - t303 * t371 + t306 * t368 - t374 * t411 + t380) * MDP(27) + (t218 * t306 + t219 * t303 - t225 * t411 - t410 * t228 - t248 * t249) * MDP(28) - t347 * MDP(4) * t413 + (t282 * t296 - t341 * t387 + t271 + t438) * MDP(23) + (t232 * t341 + t249 * t374 - t270 + t351) * MDP(26) + t362 + (-t231 * t341 + t249 * t282 + t271 + t439) * MDP(25) + (t407 * t342 + (t310 * t405 - t342 * t394) * pkin(2) + t359) * MDP(17) + (t296 * t374 + t341 * t408 - t270 + t360) * MDP(24) + t350 * t434 + (MDP(9) * t347 * t350 + MDP(10) * t413) * pkin(1); (t219 * t336 - t225 * t229 - t228 * t230 - t248 * t251 + (t218 * t345 + (-t225 * t345 + t228 * t430) * qJD(4)) * pkin(3)) * MDP(28) + (-t342 * t378 + t367 + t418) * MDP(16) + (t230 * t341 + t251 * t374 - t326 + t351) * MDP(26) + (-t229 * t374 - t230 * t282 - t336 * t371 + (t368 * t345 + (t282 * t430 - t345 * t374) * qJD(4)) * pkin(3) + t380) * MDP(27) + t362 + (-t229 * t341 + t251 * t282 + t357 + t401 + t422 - t423) * MDP(25) + (-t282 * t429 - t341 * t388 + t357 + t419) * MDP(23) + (t342 * t385 + t359) * MDP(17) + (t341 * t409 - t374 * t429 - t326 + t360) * MDP(24); (-t341 * t376 + t438) * MDP(23) + (t341 * t389 + t360) * MDP(24) + (t228 * t341 - (-t248 - t392) * t374 + t361) * MDP(25) + (-t433 * pkin(4) + t227 * t341 - (qJD(5) + t248) * t282 - t358) * MDP(26) + (-pkin(4) * t371 + t282 * t412) * MDP(27) + (pkin(4) * t439 + t412 * t228) * MDP(28) + t437; (-t368 - t420) * MDP(25) + (t371 + t421) * MDP(26) + (-t278 - t433) * MDP(27) + (-t225 * t374 - t228 * t282 + t222) * MDP(28);];
tauc = t1;
