% Calculate Coriolis joint torque vector for
% S5RRPRR6
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d2,d4,d5,theta3]';
% MDP [23x1]
%   Minimal dynamic parameter vector (fixed base model)
%   see S5RRPRR6_convert_par2_MPV_fixb.m
% 
% Output:
% tauc [5x1]
%   joint torques required to compensate Coriolis and centrifugal load

% Quelle: HybrDyn-Toolbox
% Datum: 2022-01-20 11:18
% Revision: fd3771346c4aea32fdeb66112c511235427c26a7 (2022-01-20)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauc = S5RRPRR6_coriolisvecJ_fixb_mdp_slag_vp(qJ, qJD, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(9,1),zeros(23,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRPRR6_coriolisvecJ_fixb_mdp_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RRPRR6_coriolisvecJ_fixb_mdp_slag_vp: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5RRPRR6_coriolisvecJ_fixb_mdp_slag_vp: pkin has to be [9x1] (double)');
assert(isreal(MDP) && all(size(MDP) == [23 1]), ...
  'S5RRPRR6_coriolisvecJ_fixb_mdp_slag_vp: MDP has to be [23x1] (double)'); 

%% Symbolic Calculation
% From coriolisvec_joint_fixb_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2022-01-20 11:17:29
% EndTime: 2022-01-20 11:17:32
% DurationCPUTime: 1.82s
% Computational Cost: add. (1613->215), mult. (2911->321), div. (0->0), fcn. (1831->8), ass. (0->135)
t294 = sin(qJ(5));
t292 = sin(pkin(9));
t295 = sin(qJ(4));
t364 = qJD(4) * t295;
t344 = t292 * t364;
t376 = t294 * t295;
t349 = t292 * t376;
t403 = -qJD(5) * t349 - t294 * t344;
t298 = cos(qJ(4));
t402 = t295 * MDP(12) + t298 * MDP(13);
t287 = t292 ^ 2;
t385 = t287 * t298;
t401 = -t295 * MDP(10) * t385 + (t295 ^ 2 - t298 ^ 2) * MDP(11) * t287;
t293 = cos(pkin(9));
t296 = sin(qJ(2));
t299 = cos(qJ(2));
t400 = t299 * MDP(6) + t296 * (t293 * MDP(7) + MDP(5));
t390 = pkin(1) * qJD(1);
t324 = -t299 * t390 + qJD(3);
t289 = qJD(1) + qJD(2);
t384 = t289 * t292;
t271 = -t293 * pkin(3) - t292 * pkin(7) - pkin(2);
t363 = qJD(4) * t298;
t366 = qJD(3) * t293;
t377 = t293 * t299;
t397 = -(t295 * t296 + t298 * t377) * t390 + t271 * t363 + t298 * t366;
t361 = qJD(4) + qJD(5);
t378 = t293 * t298;
t354 = qJ(3) * t378;
t394 = (t271 * t295 + t354) * qJD(4) + (-t295 * t377 + t296 * t298) * t390 + t295 * t366;
t392 = pkin(1) * t299;
t391 = pkin(8) * t292;
t389 = pkin(1) * qJD(2);
t269 = qJ(3) * t289 + t296 * t390;
t382 = t289 * t295;
t241 = (pkin(4) * t382 + t269) * t292;
t297 = cos(qJ(5));
t374 = t297 * t298;
t348 = t292 * t374;
t330 = t289 * t348;
t242 = t289 * t349 - t330;
t388 = t241 * t242;
t383 = t289 * t293;
t277 = -qJD(4) + t383;
t270 = -qJD(5) + t277;
t387 = t270 * t293;
t355 = qJD(1) * t389;
t265 = t289 * qJD(3) + t299 * t355;
t256 = t287 * t265;
t386 = t287 * t289;
t381 = t289 * t298;
t380 = t292 * t298;
t379 = t293 * t295;
t239 = t271 * t289 + t324;
t353 = t269 * t378;
t308 = -t239 * t295 - t353;
t360 = pkin(8) * t384;
t221 = -t295 * t360 - t308;
t375 = t297 * t221;
t342 = t269 * t364;
t331 = t296 * t355;
t346 = t239 * t363 + t265 * t378 + t295 * t331;
t373 = t265 * t385 + (-t293 * t342 + t346) * t293;
t288 = t293 ^ 2;
t372 = t288 * t265 + t256;
t370 = t403 * t289;
t369 = t287 + t288;
t362 = qJD(4) + t277;
t359 = pkin(8) * t380;
t357 = t296 * t389;
t356 = pkin(4) * t363;
t285 = pkin(1) * t296 + qJ(3);
t352 = t285 * t378;
t351 = t287 * t381;
t350 = t289 * t380;
t304 = -pkin(8) * t350 - t269 * t379;
t207 = t304 * qJD(4) + t346;
t323 = -t265 * t379 + t298 * t331;
t208 = (-t353 + (-t239 + t360) * t295) * qJD(4) + t323;
t236 = t298 * t239;
t220 = t236 + t304;
t214 = -pkin(4) * t277 + t220;
t217 = qJD(5) * t294 * t221;
t315 = t294 * t298 + t295 * t297;
t300 = t361 * t315;
t227 = t300 * t292;
t237 = (t289 * t356 + t265) * t292;
t314 = t374 - t376;
t254 = t314 * t292;
t347 = (t294 * t208 - t217 + (qJD(5) * t214 + t207) * t297) * t293 - t241 * t227 + t237 * t254;
t260 = t271 - t392;
t282 = t299 * t389 + qJD(3);
t345 = t260 * t363 + t282 * t378 + t295 * t357;
t343 = t293 * t364;
t340 = pkin(4) * t270 - t214;
t339 = t260 - t391;
t338 = t271 - t391;
t337 = t369 * t299;
t336 = t369 * t265;
t335 = t369 * t282;
t334 = -t294 * t207 + t297 * t208;
t333 = t369 * qJD(3);
t332 = pkin(4) * t350;
t322 = -t282 * t379 + t298 * t357;
t316 = -t214 * t294 - t375;
t201 = t316 * qJD(5) + t334;
t310 = t361 * t348;
t228 = t310 + t403;
t253 = t315 * t292;
t321 = -t201 * t293 + t241 * t228 + t237 * t253;
t212 = t308 * qJD(4) + t323;
t320 = t287 * t269 * t363 - t212 * t293 + t295 * t256;
t318 = qJD(5) * (t338 * t298 + (-qJ(3) * t295 - pkin(4)) * t293) + (-qJ(3) * t379 - t359) * qJD(4) + t397;
t278 = pkin(8) * t344;
t317 = qJD(5) * (t338 * t295 + t354) - t278 + t394;
t279 = t292 * t356;
t313 = -t324 * t292 - t279;
t307 = -t260 * t295 - t352;
t222 = t289 * t227;
t243 = t315 * t384;
t305 = -t242 * t243 * MDP(17) + (-t243 * t270 - t222) * MDP(19) + (t242 * t270 - t361 * t330 - t370) * MDP(20) + (t242 ^ 2 - t243 ^ 2) * MDP(18);
t223 = t289 * t310 + t370;
t303 = (t222 * t253 - t223 * t254 + t227 * t243 + t228 * t242) * MDP(18) + (-t222 * t254 + t227 * t242) * MDP(17) + (t222 * t293 + t227 * t270) * MDP(19) + (t223 * t293 + t228 * t270) * MDP(20) + (0.2e1 * t401 * t289 + t402 * t292 * (t277 + t383)) * qJD(4);
t301 = t241 * t243 + t217 + (t221 * t270 - t208) * t294;
t286 = t289 ^ 2;
t284 = t295 * t292 * pkin(4);
t267 = qJ(3) * t292 + t284;
t266 = -t289 * pkin(2) + t324;
t255 = t285 * t292 + t284;
t245 = t282 * t292 + t279;
t231 = t339 * t295 + t352;
t226 = t339 * t298 + (-t285 * t295 - pkin(4)) * t293;
t219 = t307 * qJD(4) + t278 + t322;
t218 = (-t285 * t379 - t359) * qJD(4) + t345;
t1 = [(-(-t218 * t294 + t297 * t219 + (-t226 * t294 - t231 * t297) * qJD(5)) * t270 + t245 * t243 + t255 * t223 + t321) * MDP(22) + ((t218 * t297 + t219 * t294 + (t226 * t297 - t231 * t294) * qJD(5)) * t270 - t245 * t242 - t255 * t222 + t347) * MDP(23) + (t289 * t335 + t372) * MDP(8) + t303 + ((-t285 * t343 + t345) * t277 + (t282 * t381 + (-t285 * t289 - t269) * t364) * t287 + t373) * MDP(16) + (t287 * t282 * t382 - t322 * t277 + (-t307 * t277 + t285 * t351) * qJD(4) + t320) * MDP(15) + (t269 * t335 + t285 * t336 + (t266 + (-pkin(2) - t392) * qJD(1)) * t357) * MDP(9) + t400 * (-qJD(1) - t289) * t389; (t269 * t333 + qJ(3) * t336 + ((-pkin(2) * qJD(2) - t266) * t296 - t269 * t337) * t390) * MDP(9) + (t267 * t223 + (t318 * t294 + t317 * t297) * t270 - t313 * t243 + t321) * MDP(22) + ((-t337 * t390 + t333) * t289 + t372) * MDP(8) + (t394 * t277 + (qJ(3) * t363 + t324 * t295) * t386 + t320) * MDP(15) + ((-qJ(3) * t343 + t397) * t277 + (-t342 + (-qJ(3) * t364 + t324 * t298) * t289) * t287 + t373) * MDP(16) + t303 + (-t267 * t222 + (-t317 * t294 + t318 * t297) * t270 + t313 * t242 + t347) * MDP(23) + t400 * (-qJD(2) + t289) * t390; -t369 * t286 * MDP(8) + (-t369 * t289 * t269 + t331) * MDP(9) + (t300 * t270 + (-t292 * t243 - t315 * t387) * t289) * MDP(22) + (t242 * t384 + (t361 * t270 - t387 * t289) * t314) * MDP(23) + (t295 * MDP(15) + t298 * MDP(16)) * (-t277 ^ 2 - t286 * t287); (-t269 * t351 + t308 * t277 + t212) * MDP(15) + (-t236 * t277 + (t362 * t293 + t386) * t295 * t269 - t346) * MDP(16) + ((-t220 * t294 - t375) * t270 - t243 * t332 + t388 + (t340 * t294 - t375) * qJD(5) + t334) * MDP(22) + (t242 * t332 + (t340 * qJD(5) - t220 * t270 - t207) * t297 + t301) * MDP(23) + t305 - t402 * t362 * t384 - t401 * t286; (t316 * t270 + t201 + t388) * MDP(22) + ((-t207 + (-qJD(5) - t270) * t214) * t297 + t301) * MDP(23) + t305;];
tauc = t1;
