% Calculate Coriolis joint torque vector for
% S6PRPRRP1
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [11x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d2,d4,d5,theta1,theta3]';
% MDP [23x1]
%   Minimal dynamic parameter vector (fixed base model)
%   see S6PRPRRP1_convert_par2_MPV_fixb.m
% 
% Output:
% tauc [6x1]
%   joint torques required to compensate Coriolis and centrifugal load

% Quelle: HybrDyn-Toolbox
% Datum: 2021-01-16 01:25
% Revision: 24b2e7d74a0c1a3b64fa2f8f5ad758691ad61af3 (2021-01-15)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauc = S6PRPRRP1_coriolisvecJ_fixb_mdp_slag_vp(qJ, qJD, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(11,1),zeros(23,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRPRRP1_coriolisvecJ_fixb_mdp_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6PRPRRP1_coriolisvecJ_fixb_mdp_slag_vp: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6PRPRRP1_coriolisvecJ_fixb_mdp_slag_vp: pkin has to be [11x1] (double)');
assert(isreal(MDP) && all(size(MDP) == [23 1]), ...
  'S6PRPRRP1_coriolisvecJ_fixb_mdp_slag_vp: MDP has to be [23x1] (double)'); 

%% Symbolic Calculation
% From coriolisvec_joint_fixb_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2021-01-16 01:24:27
% EndTime: 2021-01-16 01:24:35
% DurationCPUTime: 3.30s
% Computational Cost: add. (2594->353), mult. (6436->488), div. (0->0), fcn. (4769->10), ass. (0->160)
t321 = cos(pkin(11));
t325 = sin(qJ(2));
t320 = sin(pkin(6));
t379 = qJD(1) * t320;
t362 = t325 * t379;
t305 = t321 * t362;
t319 = sin(pkin(11));
t328 = cos(qJ(2));
t361 = t328 * t379;
t273 = t319 * t361 + t305;
t324 = sin(qJ(4));
t327 = cos(qJ(4));
t344 = pkin(4) * t324 - pkin(9) * t327;
t301 = t344 * qJD(4);
t425 = -t273 + t301;
t424 = MDP(6) * t324;
t317 = t324 ^ 2;
t423 = MDP(7) * (-t327 ^ 2 + t317);
t278 = (t319 * t325 - t321 * t328) * t320;
t303 = qJD(2) * pkin(2) + t361;
t268 = t319 * t303 + t305;
t264 = qJD(2) * pkin(8) + t268;
t322 = cos(pkin(6));
t310 = qJD(1) * t322 + qJD(3);
t249 = t327 * t264 + t324 * t310;
t246 = qJD(4) * pkin(9) + t249;
t304 = t319 * t362;
t267 = t303 * t321 - t304;
t338 = -pkin(4) * t327 - pkin(9) * t324 - pkin(3);
t254 = qJD(2) * t338 - t267;
t323 = sin(qJ(5));
t326 = cos(qJ(5));
t224 = -t246 * t323 + t326 * t254;
t376 = qJD(4) * t323;
t378 = qJD(2) * t324;
t297 = t326 * t378 + t376;
t221 = -qJ(6) * t297 + t224;
t377 = qJD(2) * t327;
t311 = -qJD(5) + t377;
t220 = -pkin(5) * t311 + t221;
t422 = t220 - t221;
t275 = t321 * t361 - t304;
t415 = pkin(2) * t321;
t292 = t338 - t415;
t372 = qJD(5) * t326;
t391 = t326 * t327;
t421 = -t275 * t391 + t292 * t372 + t425 * t323;
t420 = -t324 * t264 + t310 * t327;
t373 = qJD(5) * t323;
t357 = t324 * t373;
t366 = qJD(2) * qJD(4);
t353 = t327 * t366;
t365 = qJD(4) * qJD(5);
t381 = (t353 + t365) * t326;
t271 = qJD(2) * t357 - t381;
t360 = t323 * t378;
t367 = t326 * qJD(4);
t295 = t360 - t367;
t419 = t295 * t311 - t271;
t356 = t324 * t372;
t374 = qJD(4) * t327;
t331 = t323 * t374 + t356;
t272 = qJD(2) * t331 + t323 * t365;
t418 = -t297 * t311 + t272;
t375 = qJD(4) * t324;
t393 = t323 * t327;
t312 = pkin(2) * t319 + pkin(8);
t395 = t312 * t323;
t417 = t275 * t393 + t425 * t326 + t375 * t395;
t364 = MDP(18) + MDP(20);
t339 = qJD(2) * t317 - t311 * t327;
t416 = -t311 * t357 - t339 * t367;
t414 = pkin(5) * t295;
t413 = pkin(5) * t323;
t412 = -qJ(6) - pkin(9);
t411 = qJ(6) * t324;
t394 = t323 * t254;
t225 = t246 * t326 + t394;
t222 = -qJ(6) * t295 + t225;
t410 = t222 * t311;
t276 = qJD(2) * t278;
t270 = qJD(1) * t276;
t233 = t264 * t374 - t324 * t270 + t310 * t375;
t223 = pkin(5) * t272 + t233;
t409 = t223 * t323;
t408 = t223 * t326;
t407 = t233 * t323;
t406 = t233 * t326;
t245 = -qJD(4) * pkin(4) - t420;
t405 = t245 * t323;
t404 = t245 * t326;
t403 = t272 * t327;
t402 = t275 * t295;
t401 = t275 * t297;
t399 = t295 * t324;
t396 = t311 * t326;
t392 = t324 * t326;
t299 = t312 * t391;
t337 = pkin(5) * t324 - qJ(6) * t391;
t371 = qJD(6) * t326;
t390 = -t324 * t371 + t337 * qJD(4) + (-t299 + (-t292 + t411) * t323) * qJD(5) + t417;
t389 = (-qJ(6) * qJD(5) - qJD(4) * t312) * t392 + (-qJD(6) * t324 + (-qJ(6) * qJD(4) - qJD(5) * t312) * t327) * t323 + t421;
t300 = t344 * qJD(2);
t349 = t326 * t300 - t323 * t420;
t352 = qJD(5) * t412;
t388 = qJD(2) * t337 + qJD(6) * t323 - t326 * t352 + t349;
t355 = t323 * t377;
t386 = t323 * t300 + t326 * t420;
t387 = -qJ(6) * t355 - t323 * t352 - t371 + t386;
t359 = t327 * t367;
t385 = -t272 * t392 - t295 * t359;
t382 = t323 * t292 + t299;
t370 = t245 * qJD(5);
t369 = t323 * MDP(22);
t368 = t324 * MDP(17);
t363 = MDP(19) + MDP(21);
t358 = t311 * t373;
t354 = t324 * t366;
t351 = -qJD(6) - t414;
t232 = qJD(4) * t420 - t270 * t327;
t279 = (t319 * t328 + t321 * t325) * t320;
t255 = (qJD(1) * t279 + t301) * qJD(2);
t350 = t323 * t232 - t326 * t255;
t348 = t271 * t327 + t297 * t375;
t347 = t326 * t232 - t246 * t373 + t254 * t372 + t323 * t255;
t345 = t311 * t356;
t343 = -t220 * t326 - t222 * t323;
t342 = t220 * t323 - t222 * t326;
t262 = t279 * t327 + t322 * t324;
t239 = t262 * t326 + t278 * t323;
t238 = -t262 * t323 + t278 * t326;
t261 = t279 * t324 - t322 * t327;
t336 = qJ(6) * t271 - t350;
t274 = qJD(2) * t279;
t269 = qJD(1) * t274;
t329 = qJD(4) ^ 2;
t335 = qJD(2) * t273 - t312 * t329 - t269;
t263 = -qJD(2) * pkin(3) - t267;
t334 = qJD(4) * (qJD(2) * (-pkin(3) - t415) + t263 + t275);
t333 = -qJ(6) * t272 + t347;
t332 = t339 * t323;
t330 = qJD(2) ^ 2;
t315 = -pkin(5) * t326 - pkin(4);
t308 = t412 * t326;
t307 = t412 * t323;
t294 = t295 ^ 2;
t286 = (t312 + t413) * t324;
t282 = t326 * t292;
t266 = pkin(5) * t331 + t312 * t374;
t256 = -t323 * t411 + t382;
t253 = -qJ(6) * t392 + t282 + (-pkin(5) - t395) * t327;
t242 = pkin(5) * t355 + t249;
t237 = qJD(4) * t262 - t276 * t324;
t236 = -qJD(4) * t261 - t276 * t327;
t235 = t245 - t351;
t219 = -qJD(5) * t239 - t236 * t323 + t274 * t326;
t218 = qJD(5) * t238 + t236 * t326 + t274 * t323;
t217 = -qJD(6) * t295 + t333;
t216 = pkin(5) * t354 - t225 * qJD(5) - qJD(6) * t297 + t336;
t1 = [(-t267 * t274 - t268 * t276 + t269 * t278 - t270 * t279) * MDP(5) + (-t218 * t295 - t219 * t297 + t238 * t271 - t239 * t272) * MDP(22) + (t216 * t238 + t217 * t239 + t218 * t222 + t219 * t220 + t223 * t261 + t235 * t237) * MDP(23) + (-MDP(3) * t325 - MDP(4) * t328) * t330 * t320 + t363 * (t218 * t311 + t237 * t297 - t239 * t354 - t261 * t271) + t364 * (-t219 * t311 + t237 * t295 + t238 * t354 + t261 * t272) + (-MDP(11) * t237 - MDP(12) * t236) * qJD(4) + ((-t274 * t327 + t278 * t375) * MDP(11) + (t274 * t324 + t278 * t374) * MDP(12)) * qJD(2); (t267 * t273 - t268 * t275 + (-t269 * t321 - t270 * t319) * pkin(2)) * MDP(5) + 0.2e1 * t353 * t424 - 0.2e1 * t366 * t423 + (t324 * t334 + t327 * t335) * MDP(11) + (-t324 * t335 + t327 * t334) * MDP(12) + (-t271 * t392 + (-t357 + t359) * t297) * MDP(13) + (-t297 * t356 + (-t297 * t374 + (qJD(5) * t295 + t271) * t324) * t323 + t385) * MDP(14) + (t348 - t416) * MDP(15) + (t345 + t403 + (-t332 - t399) * qJD(4)) * MDP(16) + (-t311 - t377) * qJD(4) * t368 + ((t292 * t373 - t417) * t311 + ((t295 * t312 + t405) * qJD(4) + (t394 + (t311 * t312 + t246) * t326) * qJD(5) + t350) * t327 + (t326 * t370 + t407 + t312 * t272 - t402 + ((-t312 * t393 + t282) * qJD(2) + t224) * qJD(4)) * t324) * MDP(18) + (t421 * t311 + (-t312 * t358 + (t297 * t312 + t404) * qJD(4) + t347) * t327 + (-t323 * t370 + t406 - t312 * t271 - t401 + (-qJD(2) * t382 - t312 * t396 - t225) * qJD(4)) * t324) * MDP(19) + (t266 * t295 + t272 * t286 + (t235 * t376 - t216) * t327 - t390 * t311 + (t235 * t372 + t409 - t402 + (qJD(2) * t253 + t220) * qJD(4)) * t324) * MDP(20) + (t266 * t297 - t271 * t286 + (t235 * t367 + t217) * t327 + t389 * t311 + (-t235 * t373 + t408 - t401 + (-qJD(2) * t256 - t222) * qJD(4)) * t324) * MDP(21) + (t253 * t271 - t256 * t272 - t390 * t297 - t389 * t295 + t343 * t374 + (qJD(5) * t342 - t216 * t326 - t217 * t323) * t324) * MDP(22) + (t216 * t253 + t217 * t256 + t223 * t286 + (-t275 * t324 + t266) * t235 + t389 * t222 + t390 * t220) * MDP(23) + (MDP(8) * t327 - MDP(9) * t324) * t329; t385 * MDP(22) + t363 * (t348 + t416) + t364 * (t345 - t403 + (-t332 + t399) * qJD(4)) + (-t329 * MDP(12) - t223 * MDP(23) + (-MDP(23) * t342 + t297 * t369) * qJD(4)) * t327 + (-t329 * MDP(11) - t271 * t369 + (qJD(4) * t235 - t216 * t323 + t217 * t326) * MDP(23) + ((t295 * t323 + t297 * t326) * MDP(22) + t343 * MDP(23)) * qJD(5)) * t324; (qJD(4) * t249 - t263 * t378 - t233) * MDP(11) + (-qJD(2) * t263 + t270) * t327 * MDP(12) + (-t271 * t323 - t297 * t396) * MDP(13) + (-t418 * t323 + t419 * t326) * MDP(14) + (-t311 * t372 + (t311 * t391 + (-t297 + t376) * t324) * qJD(2)) * MDP(15) + (t358 + (-t311 * t393 + (t295 + t367) * t324) * qJD(2)) * MDP(16) + t311 * qJD(2) * t368 + (-pkin(4) * t272 - t406 + t349 * t311 - t249 * t295 + (pkin(9) * t396 + t405) * qJD(5) + (-t224 * t324 + (-pkin(9) * t375 - t245 * t327) * t323) * qJD(2)) * MDP(18) + (pkin(4) * t271 + t407 - t249 * t297 - t386 * t311 + (-pkin(9) * t311 * t323 + t404) * qJD(5) + (-t245 * t391 + (-pkin(9) * t367 + t225) * t324) * qJD(2)) * MDP(19) + (-t408 - t242 * t295 + t272 * t315 + t388 * t311 + (t235 + t414) * t373 + (-t235 * t393 + (qJD(4) * t307 - t220) * t324) * qJD(2)) * MDP(20) + (t409 - t242 * t297 - t271 * t315 - t387 * t311 + (t235 * t326 + t297 * t413) * qJD(5) + (-t235 * t391 + (qJD(4) * t308 + t222) * t324) * qJD(2)) * MDP(21) + (t271 * t307 + t272 * t308 + t388 * t297 + t387 * t295 + (t220 * t311 + t217) * t326 + (-t216 + t410) * t323) * MDP(22) + (t216 * t307 - t217 * t308 + t223 * t315 + (pkin(5) * t373 - t242) * t235 - t387 * t222 - t388 * t220) * MDP(23) + (-t327 * t424 + t423) * t330; -t294 * MDP(14) + t381 * MDP(15) + (-t225 * t311 - t350) * MDP(18) + (-t224 * t311 - t347) * MDP(19) + (t336 - t410) * MDP(20) + (-t221 * t311 - t333) * MDP(21) + pkin(5) * t271 * MDP(22) + (pkin(5) * t216 + t422 * t222) * MDP(23) + (-t311 * MDP(15) + t245 * MDP(19) + (qJD(6) + t235) * MDP(21) - t422 * MDP(22)) * t295 + (-MDP(16) * t393 + (0.2e1 * MDP(20) * pkin(5) + MDP(17)) * t324) * t366 + (-MDP(15) * t360 - MDP(16) * t297 - t364 * t225) * qJD(5) + (t295 * MDP(13) - t311 * MDP(16) - t245 * MDP(18) + (-t235 + t351) * MDP(20) - pkin(5) * t235 * MDP(23) + (-pkin(5) * MDP(21) + MDP(14)) * t297) * t297; t418 * MDP(20) + t419 * MDP(21) + (-t297 ^ 2 - t294) * MDP(22) + (t220 * t297 + t222 * t295 + t223) * MDP(23);];
tauc = t1;
