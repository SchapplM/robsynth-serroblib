% Calculate minimal parameter regressor of Coriolis joint torque vector for
% S5RRRPP7
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [7x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d2,d3]';
% MDP [25x1]
%   Minimal dynamic parameter vector (fixed base model)
%   see S5RRRPP7_convert_par2_MPV_fixb.m
% 
% Output:
% tauc [5x1]
%   joint torques required to compensate Coriolis and centrifugal load

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 21:06
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauc = S5RRRPP7_coriolisvecJ_fixb_mdp_slag_vp(qJ, qJD, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(7,1),zeros(25,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRRPP7_coriolisvecJ_fixb_mdp_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RRRPP7_coriolisvecJ_fixb_mdp_slag_vp: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [7 1]), ...
  'S5RRRPP7_coriolisvecJ_fixb_mdp_slag_vp: pkin has to be [7x1] (double)');
assert(isreal(MDP) && all(size(MDP) == [25 1]), ...
  'S5RRRPP7_coriolisvecJ_fixb_mdp_slag_vp: MDP has to be [25x1] (double)'); 

%% Symbolic Calculation
% From coriolisvec_joint_fixb_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 21:06:01
% EndTime: 2019-12-31 21:06:08
% DurationCPUTime: 3.43s
% Computational Cost: add. (2093->395), mult. (5013->511), div. (0->0), fcn. (2877->4), ass. (0->161)
t319 = cos(qJ(2));
t380 = qJD(1) * t319;
t303 = -qJD(3) + t380;
t317 = sin(qJ(2));
t365 = qJD(1) * qJD(2);
t354 = t317 * t365;
t428 = qJ(4) * t354 - t303 * qJD(4);
t427 = MDP(20) + MDP(23);
t426 = -0.2e1 * t365;
t425 = t317 * MDP(4);
t314 = t317 ^ 2;
t424 = (-t319 ^ 2 + t314) * MDP(5);
t289 = -pkin(2) * t319 - pkin(7) * t317 - pkin(1);
t270 = t289 * qJD(1);
t343 = pkin(2) * t317 - pkin(7) * t319;
t287 = t343 * qJD(2);
t273 = qJD(1) * t287;
t311 = pkin(6) * t380;
t295 = qJD(2) * pkin(7) + t311;
t316 = sin(qJ(3));
t318 = cos(qJ(3));
t346 = pkin(6) * t354;
t375 = qJD(3) * t318;
t376 = qJD(3) * t316;
t345 = -t270 * t376 + t318 * t273 - t295 * t375 + t316 * t346;
t223 = -pkin(3) * t354 - t345;
t243 = t316 * t270 + t318 * t295;
t298 = t303 * qJ(4);
t237 = -t298 + t243;
t423 = t237 * t303 + t223;
t364 = qJD(2) * qJD(3);
t352 = t316 * t364;
t355 = t317 * t375;
t377 = qJD(2) * t319;
t255 = (t316 * t377 + t355) * qJD(1) + t352;
t369 = t318 * qJD(2);
t381 = qJD(1) * t317;
t280 = t316 * t381 - t369;
t422 = t255 * qJ(5) + t280 * qJD(5);
t353 = t319 * t365;
t357 = t317 * t376;
t254 = qJD(1) * t357 + (-t353 - t364) * t318;
t403 = t280 * t303;
t421 = -t254 + t403;
t358 = t318 * t381;
t379 = qJD(2) * t316;
t282 = t358 + t379;
t420 = t282 * t303 - t255;
t419 = 0.2e1 * t428;
t417 = -qJD(4) * t316 - t311;
t242 = t318 * t270 - t316 * t295;
t367 = qJD(4) - t242;
t410 = qJ(4) * t316;
t415 = pkin(3) + pkin(4);
t416 = -t415 * t318 - t410;
t278 = t282 ^ 2;
t414 = pkin(6) * t316;
t413 = pkin(7) - qJ(5);
t412 = qJ(4) * t255;
t411 = qJ(4) * t280;
t409 = qJ(4) * t318;
t222 = t255 * pkin(3) + pkin(6) * t353 + t254 * qJ(4) - t282 * qJD(4);
t408 = t222 * t316;
t407 = t222 * t318;
t232 = qJ(5) * t280 + t243;
t227 = t232 - t298;
t406 = t227 * t303;
t404 = t254 * t316;
t286 = t343 * qJD(1);
t401 = t286 * t318;
t294 = -qJD(2) * pkin(2) + pkin(6) * t381;
t400 = t294 * t316;
t399 = t294 * t318;
t398 = t303 * t318;
t397 = t316 * t317;
t396 = t316 * t319;
t395 = t317 * t318;
t321 = qJD(2) ^ 2;
t394 = t317 * t321;
t393 = t318 * t319;
t392 = t319 * t321;
t322 = qJD(1) ^ 2;
t391 = t319 * t322;
t293 = t413 * t318;
t359 = -pkin(3) - t414;
t347 = -pkin(4) + t359;
t390 = -t401 + (-qJ(5) * t393 + t347 * t317) * qJD(1) - qJD(3) * t293 + qJD(5) * t316;
t372 = qJD(5) * t318;
t268 = t316 * t286;
t386 = qJ(4) * t381 + t268;
t389 = (-pkin(6) * t395 + qJ(5) * t396) * qJD(1) + t386 + t413 * t376 + t372;
t333 = -t415 * t316 + t409;
t388 = t303 * t333 + t417;
t340 = pkin(3) * t316 - t409;
t387 = t303 * t340 - t417;
t385 = t316 * t287 + t289 * t375;
t306 = pkin(6) * t393;
t384 = qJD(3) * t306 + t289 * t376;
t383 = t316 * t289 + t306;
t378 = qJD(2) * t317;
t373 = qJD(4) * t318;
t371 = t294 * qJD(3);
t370 = t317 * MDP(15);
t231 = qJ(5) * t282 + t242;
t368 = qJD(4) - t231;
t334 = qJ(4) * t282 - t294;
t228 = -t415 * t280 + qJD(5) + t334;
t366 = qJD(5) + t228;
t363 = pkin(7) * t303 * t316;
t362 = pkin(7) * t398;
t361 = pkin(7) * t378;
t360 = pkin(7) * t369;
t356 = t319 * t376;
t351 = pkin(1) * t426;
t305 = pkin(6) * t396;
t350 = t289 * t318 - t305;
t349 = t280 + t369;
t348 = -t282 + t379;
t344 = t359 * t317;
t252 = -qJ(4) * t319 + t383;
t342 = t287 * t318 - t384;
t341 = pkin(3) * t318 + t410;
t236 = pkin(3) * t303 + t367;
t339 = t236 * t318 - t237 * t316;
t338 = qJD(1) * t314 - t303 * t319;
t337 = pkin(6) + t340;
t336 = -t270 * t375 - t316 * t273 + t295 * t376;
t335 = qJ(4) * t378 - qJD(4) * t319 + t385;
t218 = -pkin(4) * t255 - t222;
t332 = -t218 * t316 - t228 * t375;
t331 = t218 * t318 - t228 * t376;
t329 = qJ(5) * t254 - t345;
t327 = -pkin(6) + t333;
t325 = -t242 * t303 + t336;
t324 = t254 + t403;
t323 = -t415 * t354 + t329;
t221 = -t318 * t346 - t336 + t428;
t313 = t319 * pkin(3);
t292 = t413 * t316;
t288 = -pkin(2) - t341;
t276 = pkin(2) - t416;
t260 = t337 * t317;
t253 = t313 - t350;
t251 = t327 * t317;
t247 = pkin(3) * t282 + t411;
t246 = qJD(1) * t344 - t401;
t245 = -pkin(6) * t358 + t386;
t241 = qJ(5) * t397 + t252;
t240 = pkin(4) * t319 + t305 + t313 + (-qJ(5) * t317 - t289) * t318;
t239 = pkin(3) * t280 - t334;
t234 = -t415 * t282 - t411;
t230 = (t341 * qJD(3) - t373) * t317 + t337 * t377;
t229 = qJD(2) * t344 - t342;
t226 = (-t317 * t369 - t356) * pkin(6) + t335;
t225 = (t416 * qJD(3) + t373) * t317 + t327 * t377;
t224 = t415 * t303 + t368;
t220 = (-pkin(6) * qJD(2) + qJ(5) * qJD(3)) * t395 + (qJD(5) * t317 + (-pkin(6) * qJD(3) + qJ(5) * qJD(2)) * t319) * t316 + t335;
t219 = (-qJ(5) * t377 - t287) * t318 + (qJ(5) * t376 + t347 * qJD(2) - t372) * t317 + t384;
t217 = -qJD(5) * t282 + t323;
t216 = t221 + t422;
t1 = [0.2e1 * t353 * t425 + t424 * t426 + MDP(6) * t392 - MDP(7) * t394 + (-pkin(6) * t392 + t317 * t351) * MDP(9) + (pkin(6) * t394 + t319 * t351) * MDP(10) + (-t254 * t395 + (t319 * t369 - t357) * t282) * MDP(11) + ((-t280 * t318 - t282 * t316) * t377 + (t404 - t255 * t318 + (t280 * t316 - t282 * t318) * qJD(3)) * t317) * MDP(12) + (t303 * t357 + t254 * t319 + (t282 * t317 + t338 * t318) * qJD(2)) * MDP(13) + (t303 * t355 + t255 * t319 + (-t280 * t317 - t338 * t316) * qJD(2)) * MDP(14) + (-t303 - t380) * qJD(2) * t370 + (-t342 * t303 - t345 * t319 + (pkin(6) * t255 + t318 * t371) * t317 + ((pkin(6) * t280 + t400) * t319 + (t350 * qJD(1) + t242 + (-t303 + t380) * t414) * t317) * qJD(2)) * MDP(16) + ((-pkin(6) * t356 + t385) * t303 - t336 * t319 + (-pkin(6) * t254 - t316 * t371) * t317 + ((pkin(6) * t282 + t399) * t319 + (-pkin(6) * t398 - t383 * qJD(1) - t243) * t317) * qJD(2)) * MDP(17) + (t229 * t303 + t230 * t280 + t255 * t260 + (t239 * t379 + t223) * t319 + (t239 * t375 + t408 + (-qJD(1) * t253 - t236) * qJD(2)) * t317) * MDP(18) + (-t226 * t280 + t229 * t282 - t252 * t255 - t253 * t254 + t339 * t377 + (-t221 * t316 + t223 * t318 + (-t236 * t316 - t237 * t318) * qJD(3)) * t317) * MDP(19) + (-t226 * t303 - t230 * t282 + t254 * t260 + (-t239 * t369 - t221) * t319 + (t239 * t376 - t407 + (qJD(1) * t252 + t237) * qJD(2)) * t317) * MDP(20) + (t221 * t252 + t222 * t260 + t223 * t253 + t226 * t237 + t229 * t236 + t230 * t239) * MDP(21) + (t219 * t303 - t225 * t280 - t251 * t255 + (-t228 * t379 + t217) * t319 + ((-qJD(1) * t240 - t224) * qJD(2) + t332) * t317) * MDP(22) + (-t220 * t303 + t225 * t282 - t251 * t254 + (t228 * t369 - t216) * t319 + ((qJD(1) * t241 + t227) * qJD(2) + t331) * t317) * MDP(23) + (-t219 * t282 + t220 * t280 + t240 * t254 + t241 * t255 + (-t224 * t318 + t227 * t316) * t377 + (t216 * t316 - t217 * t318 + (t224 * t316 + t227 * t318) * qJD(3)) * t317) * MDP(24) + (t216 * t241 + t217 * t240 + t218 * t251 + t219 * t224 + t220 * t227 + t225 * t228) * MDP(25); -t391 * t425 + t322 * t424 + (-t282 * t398 - t404) * MDP(11) + (t420 * t316 + t421 * t318) * MDP(12) + (-t303 * t375 + (t303 * t393 + t348 * t317) * qJD(1)) * MDP(13) + (t303 * t376 + (-t303 * t396 + t349 * t317) * qJD(1)) * MDP(14) + t303 * qJD(1) * t370 + (t286 * t398 - pkin(2) * t255 + (t362 + t400) * qJD(3) + (-t242 * t317 + (-t294 * t319 - t361) * t316 + (t303 * t397 - t349 * t319) * pkin(6)) * qJD(1)) * MDP(16) + (pkin(2) * t254 - t268 * t303 + (-t363 + t399) * qJD(3) + (-t294 * t393 + (t243 - t360) * t317 + (t303 * t395 + t348 * t319) * pkin(6)) * qJD(1)) * MDP(17) + (-t407 - t246 * t303 + t255 * t288 - t387 * t280 + (t239 * t316 + t362) * qJD(3) + (t236 * t317 + (-t239 * t319 - t361) * t316) * qJD(1)) * MDP(18) + (t245 * t280 - t246 * t282 + (t221 - t303 * t236 + (qJD(3) * t282 - t255) * pkin(7)) * t318 + ((qJD(3) * t280 - t254) * pkin(7) + t423) * t316) * MDP(19) + (-t408 + t245 * t303 + t254 * t288 + t387 * t282 + (-t239 * t318 + t363) * qJD(3) + (t239 * t393 + (-t237 + t360) * t317) * qJD(1)) * MDP(20) + (t222 * t288 - t236 * t246 - t237 * t245 - t387 * t239 + (t339 * qJD(3) + t221 * t318 + t223 * t316) * pkin(7)) * MDP(21) + (-t255 * t276 - t390 * t303 + t388 * t280 + (t228 * t396 + (-qJD(2) * t292 + t224) * t317) * qJD(1) + t331) * MDP(22) + (-t254 * t276 + t389 * t303 - t388 * t282 + (-t228 * t393 + (qJD(2) * t293 - t227) * t317) * qJD(1) - t332) * MDP(23) + (t254 * t292 + t255 * t293 + t390 * t282 - t389 * t280 + (t303 * t224 - t216) * t318 + (-t217 - t406) * t316) * MDP(24) + (t216 * t293 + t217 * t292 + t218 * t276 - t390 * t224 - t389 * t227 - t388 * t228) * MDP(25) + (t322 * t317 * MDP(9) + MDP(10) * t391) * pkin(1); -t324 * MDP(13) - MDP(14) * t352 + t325 * MDP(17) + (pkin(3) * t254 - t412) * MDP(19) + (-t325 + t419) * MDP(20) + (-pkin(3) * t223 + qJ(4) * t221 - t236 * t243 + t367 * t237 - t239 * t247) * MDP(21) + (-t232 * t303 - t329) * MDP(22) + (t231 * t303 - t336 + t419 + t422) * MDP(23) + (-t254 * t415 + t412) * MDP(24) + (qJ(4) * t216 - t217 * t415 - t224 * t232 + t368 * t227 - t228 * t234) * MDP(25) + (-t303 * MDP(14) - t294 * MDP(16) - t239 * MDP(18) + (t237 - t243) * MDP(19) + t247 * MDP(20) + t366 * MDP(22) - t234 * MDP(23) + (-t227 + t232) * MDP(24) + MDP(12) * t282) * t282 + (-MDP(14) * t355 + (-MDP(14) * t396 + (MDP(15) + 0.2e1 * pkin(3) * MDP(18) + 0.2e1 * t415 * MDP(22) + (MDP(17) - t427) * t318 * pkin(6)) * t317) * qJD(2)) * qJD(1) + (t282 * MDP(11) + t294 * MDP(17) - t247 * MDP(18) + (t236 - t367) * MDP(19) - t239 * MDP(20) + t234 * MDP(22) + t228 * MDP(23) + (-t224 + t368) * MDP(24) - MDP(12) * t280) * t280 + (MDP(16) + MDP(18)) * (-t243 * t303 + t345); (t239 * t282 + t423) * MDP(21) + (-t366 * t282 + t323 + t406) * MDP(25) + (MDP(18) + MDP(22)) * (t280 * t282 - t354) + t427 * (-t303 ^ 2 - t278) + (-MDP(19) + MDP(24)) * t324; t420 * MDP(22) + t421 * MDP(23) + (-t280 ^ 2 - t278) * MDP(24) + (t224 * t282 - t227 * t280 + t218) * MDP(25);];
tauc = t1;
