% Calculate minimal parameter regressor of Coriolis joint torque vector for
% S6RPPRRR6
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d4,d5,d6]';
% MDP [30x1]
%   Minimal dynamic parameter vector (fixed base model)
%   see S6RPPRRR6_convert_par2_MPV_fixb.m
% 
% Output:
% tauc [6x1]
%   joint torques required to compensate Coriolis and centrifugal load

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 02:32
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauc = S6RPPRRR6_coriolisvecJ_fixb_mdp_slag_vp(qJ, qJD, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(9,1),zeros(30,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPPRRR6_coriolisvecJ_fixb_mdp_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RPPRRR6_coriolisvecJ_fixb_mdp_slag_vp: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S6RPPRRR6_coriolisvecJ_fixb_mdp_slag_vp: pkin has to be [9x1] (double)');
assert(isreal(MDP) && all(size(MDP) == [30 1]), ...
  'S6RPPRRR6_coriolisvecJ_fixb_mdp_slag_vp: MDP has to be [30x1] (double)'); 

%% Symbolic Calculation
% From coriolisvec_joint_fixb_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 02:31:35
% EndTime: 2019-03-09 02:31:43
% DurationCPUTime: 3.42s
% Computational Cost: add. (1850->362), mult. (3880->502), div. (0->0), fcn. (2474->6), ass. (0->155)
t321 = sin(qJ(5));
t324 = cos(qJ(5));
t374 = t324 * qJD(4);
t325 = cos(qJ(4));
t387 = qJD(1) * t325;
t280 = t321 * t387 - t374;
t323 = cos(qJ(6));
t363 = t324 * t387;
t282 = qJD(4) * t321 + t363;
t320 = sin(qJ(6));
t406 = t282 * t320;
t236 = t323 * t280 + t406;
t322 = sin(qJ(4));
t388 = qJD(1) * t322;
t309 = qJD(5) + t388;
t304 = qJD(6) + t309;
t425 = t236 * t304;
t339 = t280 * t320 - t323 * t282;
t424 = t304 * t339;
t284 = t320 * t324 + t321 * t323;
t370 = qJD(5) + qJD(6);
t328 = t370 * t284;
t318 = -pkin(7) + qJ(2);
t381 = qJD(4) * t325;
t423 = qJD(2) * t322 + t318 * t381;
t319 = pkin(1) + qJ(3);
t289 = pkin(4) * t322 - pkin(8) * t325 + t319;
t260 = qJD(1) * t289 - qJD(2);
t313 = qJD(1) * qJ(2) + qJD(3);
t302 = -qJD(1) * pkin(7) + t313;
t288 = t322 * t302;
t270 = qJD(4) * pkin(8) + t288;
t232 = t260 * t321 + t270 * t324;
t228 = -pkin(9) * t280 + t232;
t376 = qJD(6) * t320;
t224 = t228 * t376;
t271 = -qJD(4) * pkin(4) - t302 * t325;
t243 = pkin(5) * t280 + t271;
t422 = t236 * t243 + t224;
t377 = qJD(5) * t325;
t357 = t321 * t377;
t331 = -t322 * t374 - t357;
t248 = qJD(1) * t331 + qJD(5) * t374;
t346 = pkin(4) * t325 + pkin(8) * t322;
t278 = qJD(4) * t346 + qJD(3);
t261 = t278 * qJD(1);
t254 = t324 * t261;
t372 = qJD(1) * qJD(2);
t259 = t302 * t381 + t322 * t372;
t329 = -qJD(5) * t232 - t321 * t259 + t254;
t371 = qJD(1) * qJD(4);
t355 = t325 * t371;
t217 = pkin(5) * t355 - pkin(9) * t248 + t329;
t382 = qJD(4) * t322;
t361 = t321 * t382;
t249 = -qJD(1) * t361 + qJD(5) * t282;
t378 = qJD(5) * t324;
t366 = -t324 * t259 - t260 * t378 - t321 * t261;
t380 = qJD(5) * t321;
t334 = -t270 * t380 - t366;
t218 = -pkin(9) * t249 + t334;
t354 = t323 * t217 - t320 * t218;
t421 = t243 * t339 + t354;
t420 = MDP(28) * t355 + (-t236 ^ 2 + t339 ^ 2) * MDP(25) - t236 * MDP(24) * t339;
t317 = t325 ^ 2;
t419 = MDP(11) * (t322 ^ 2 - t317);
t418 = qJD(1) * t319;
t417 = MDP(6) * qJ(2) + MDP(5) + MDP(7);
t353 = t248 * t320 + t323 * t249;
t220 = -qJD(6) * t339 + t353;
t416 = pkin(8) + pkin(9);
t231 = t324 * t260 - t270 * t321;
t227 = -pkin(9) * t282 + t231;
t223 = pkin(5) * t309 + t227;
t414 = t223 * t323;
t413 = t228 * t323;
t412 = t248 * t321;
t258 = t302 * t382 - t325 * t372;
t411 = t258 * t321;
t410 = t258 * t324;
t409 = t271 * t324;
t408 = t280 * t309;
t407 = t282 * t309;
t405 = t282 * t325;
t404 = t309 * t321;
t403 = t309 * t324;
t402 = t318 * t321;
t401 = t320 * t321;
t400 = t321 * t322;
t399 = t321 * t325;
t398 = t322 * t324;
t397 = t323 * t324;
t396 = t324 * t325;
t364 = t321 * t388;
t375 = qJD(6) * t323;
t395 = t320 * t364 - t323 * t378 - t324 * t375 + t370 * t401 - t388 * t397;
t332 = t284 * qJD(1);
t394 = t322 * t332 + t328;
t285 = t346 * qJD(1);
t393 = t321 * t285 + t302 * t396;
t392 = t321 * t289 + t318 * t398;
t283 = -t397 + t401;
t385 = qJD(4) * t283;
t384 = qJD(4) * t284;
t383 = qJD(4) * t304;
t379 = qJD(5) * t322;
t303 = -qJD(2) + t418;
t373 = qJD(2) - t303;
t368 = 0.2e1 * qJD(3) * qJD(1);
t367 = t323 * t248 - t320 * t249 - t280 * t375;
t365 = qJD(5) * t416;
t359 = t309 * t380;
t358 = t318 * t379;
t356 = t324 * t377;
t352 = t309 * t318 + t270;
t351 = qJD(6) * t223 + t218;
t350 = t373 * qJD(1);
t349 = qJD(1) + t379;
t348 = t321 * t278 + t289 * t378 + t423 * t324;
t300 = t322 * t355;
t347 = -t288 + (t364 + t380) * pkin(5);
t345 = t324 * t285 - t302 * t399;
t344 = t309 * t378 + t321 * t355;
t297 = t416 * t324;
t336 = pkin(5) * t325 + pkin(9) * t398;
t343 = qJD(1) * t336 + qJD(6) * t297 + t324 * t365 + t345;
t296 = t416 * t321;
t342 = pkin(9) * t364 + qJD(6) * t296 + t321 * t365 + t393;
t341 = qJD(2) + t303 + t418;
t215 = t223 * t320 + t413;
t275 = t324 * t289;
t234 = -pkin(9) * t396 + t275 + (pkin(5) - t402) * t322;
t240 = -pkin(9) * t399 + t392;
t340 = t234 * t320 + t240 * t323;
t338 = qJD(1) * t317 - t309 * t322;
t326 = qJD(4) ^ 2;
t337 = -t318 * t326 + t368;
t335 = -pkin(8) * t381 + t271 * t322;
t219 = -t282 * t376 + t367;
t333 = qJD(1) * t283;
t330 = -t356 + t361;
t327 = qJD(1) ^ 2;
t312 = -pkin(5) * t324 - pkin(4);
t277 = (pkin(5) * t321 - t318) * t325;
t266 = t324 * t278;
t263 = t283 * t325;
t262 = t284 * t325;
t244 = -pkin(5) * t330 - qJD(2) * t325 + t318 * t382;
t229 = pkin(5) * t249 + t258;
t226 = -t376 * t399 + (t370 * t396 - t361) * t323 + t331 * t320;
t225 = t283 * t382 - t325 * t328;
t222 = pkin(9) * t330 - t321 * t358 + t348;
t221 = -t324 * t358 + t266 + t336 * qJD(4) + ((pkin(9) * t325 - t289) * qJD(5) - t423) * t321;
t214 = -t228 * t320 + t414;
t1 = [-0.2e1 * MDP(10) * t300 + 0.2e1 * t371 * t419 + (-t309 * t356 - t249 * t322 + (-t280 * t325 - t321 * t338) * qJD(4)) * MDP(20) + (t248 * t396 + t282 * t331) * MDP(17) + ((t221 * t323 - t222 * t320) * t304 + t354 * t322 + t244 * t236 + t277 * t220 + t229 * t262 + t243 * t226 + (-t215 * t322 - t304 * t340) * qJD(6) + ((t234 * t323 - t240 * t320) * qJD(1) + t214) * t381) * MDP(29) + (t219 * t322 + t225 * t304 + (-qJD(1) * t263 - t339) * t381) * MDP(26) + (-t220 * t322 - t226 * t304 + (-qJD(1) * t262 - t236) * t381) * MDP(27) + (t277 * t219 + t224 * t322 + t243 * t225 - t229 * t263 - t244 * t339 + (-(-qJD(6) * t240 + t221) * t304 - t217 * t322) * t320 + (-(qJD(6) * t234 + t222) * t304 - t351 * t322) * t323 + (-qJD(1) * t340 - t215) * t381) * MDP(30) + (t322 * t337 + t341 * t381) * MDP(15) + (t309 * t381 + t300) * MDP(21) + (t304 * t381 + t300) * MDP(28) + (t325 * t337 - t341 * t382) * MDP(16) + (qJD(2) * t313 + qJD(3) * t303 + (qJ(2) * qJD(2) + qJD(3) * t319) * qJD(1)) * MDP(9) + MDP(8) * t368 + ((t280 * t324 + t282 * t321) * t382 + (-t412 - t249 * t324 + (t280 * t321 - t282 * t324) * qJD(5)) * t325) * MDP(18) + (-t309 * t357 + t248 * t322 + (t324 * t338 + t405) * qJD(4)) * MDP(19) + (-t348 * t309 + (t352 * t380 + (t282 * t318 - t409) * qJD(4) + t366) * t322 + (-t271 * t380 - qJD(2) * t282 - t318 * t248 + t410 + (-qJD(1) * t392 - t232) * qJD(4)) * t325) * MDP(23) + ((-t289 * t380 + t266) * t309 + (qJD(4) * t318 * t280 + t254 - t352 * t378 + (-qJD(2) * t309 - qJD(4) * t271 - qJD(5) * t260 - t259) * t321) * t322 + (t271 * t378 - qJD(2) * t280 - t318 * t249 + t411 + (-t309 * t402 + (-t318 * t400 + t275) * qJD(1) + t231) * qJD(4)) * t325) * MDP(22) + (-t219 * t263 - t225 * t339) * MDP(24) + (-t219 * t262 + t220 * t263 - t225 * t236 + t226 * t339) * MDP(25) + (-MDP(12) * t322 - MDP(13) * t325) * t326 + 0.2e1 * t417 * t372; MDP(22) * t359 + t344 * MDP(23) - t417 * t327 + (MDP(29) * t394 - MDP(30) * t395) * t304 + ((-qJD(3) - t313) * MDP(9) + (0.2e1 * qJD(4) * MDP(16) + (MDP(22) * t321 + MDP(23) * t324) * t309) * t322 + (-0.2e1 * qJD(4) * MDP(15) + (t280 - t374) * MDP(22) + t282 * MDP(23) + (t236 + t385) * MDP(29) + (-t339 + t384) * MDP(30)) * t325) * qJD(1); -t327 * MDP(8) + MDP(9) * t350 + (-t249 * t325 - t349 * t403 + (t280 * t322 + (-t309 - t388) * t399) * qJD(4)) * MDP(22) + (-t248 * t325 + t349 * t404 + (-t309 * t396 + (t282 - t363) * t322) * qJD(4)) * MDP(23) + (t304 * t333 + (-t284 * t383 - t220) * t325 + ((-t284 * t387 + t236) * qJD(4) + t370 * t304 * t283) * t322) * MDP(29) + (t304 * t332 + (t283 * t383 - t219) * t325 + (t328 * t304 + (t325 * t333 - t339) * qJD(4)) * t322) * MDP(30) + (MDP(15) * t322 + MDP(16) * t325) * (-t326 - t327); t325 * MDP(15) * t350 - t373 * MDP(16) * t388 + (t282 * t403 + t412) * MDP(17) + ((t248 - t408) * t324 + (-t249 - t407) * t321) * MDP(18) + ((t309 * t398 - t405) * qJD(1) + t344) * MDP(19) + (-t359 + (-t309 * t400 + (t280 + t374) * t325) * qJD(1)) * MDP(20) + (-pkin(4) * t249 - t410 - t345 * t309 - t280 * t288 + (-pkin(8) * t403 + t271 * t321) * qJD(5) + (-t231 * t325 + t321 * t335) * qJD(1)) * MDP(22) + (-pkin(4) * t248 + t411 + t393 * t309 - t282 * t288 + (pkin(8) * t404 + t409) * qJD(5) + (t232 * t325 + t324 * t335) * qJD(1)) * MDP(23) + (t219 * t284 + t339 * t395) * MDP(24) + (-t219 * t283 - t220 * t284 + t236 * t395 + t339 * t394) * MDP(25) + (t312 * t220 + t229 * t283 + t347 * t236 + t394 * t243) * MDP(29) + (t312 * t219 + t229 * t284 - t395 * t243 - t339 * t347) * MDP(30) + (-t395 * MDP(26) - t394 * MDP(27) + (t320 * t342 - t323 * t343) * MDP(29) + (t320 * t343 + t323 * t342) * MDP(30)) * t304 + (-t309 * MDP(21) + (t339 + t384) * MDP(26) + (t236 - t385) * MDP(27) - t304 * MDP(28) + ((-t296 * t323 - t297 * t320) * qJD(4) - t214) * MDP(29) + (-(-t296 * t320 + t297 * t323) * qJD(4) + t215) * MDP(30)) * t387 + (MDP(10) * t322 * t325 - t419) * t327; t282 * t280 * MDP(17) + (-t280 ^ 2 + t282 ^ 2) * MDP(18) + (t248 + t408) * MDP(19) + (-t249 + t407) * MDP(20) + MDP(21) * t355 + (t232 * t309 - t271 * t282 + t329) * MDP(22) + (t231 * t309 + t271 * t280 - t334) * MDP(23) + (t219 + t425) * MDP(26) + (-t220 - t424) * MDP(27) + (-(-t227 * t320 - t413) * t304 - t215 * qJD(6) + (-t236 * t282 - t304 * t376 + t323 * t355) * pkin(5) + t421) * MDP(29) + ((-t228 * t304 - t217) * t320 + (t227 * t304 - t351) * t323 + (t282 * t339 - t304 * t375 - t320 * t355) * pkin(5) + t422) * MDP(30) + t420; (t367 + t425) * MDP(26) + (-t353 - t424) * MDP(27) + (t215 * t304 + t421) * MDP(29) + (t214 * t304 - t320 * t217 - t323 * t218 + t422) * MDP(30) + (-MDP(26) * t406 + MDP(27) * t339 - MDP(29) * t215 - MDP(30) * t414) * qJD(6) + t420;];
tauc  = t1;
