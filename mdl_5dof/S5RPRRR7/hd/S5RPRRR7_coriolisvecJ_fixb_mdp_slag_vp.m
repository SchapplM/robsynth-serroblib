% Calculate minimal parameter regressor of Coriolis joint torque vector for
% S5RPRRR7
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d3,d4,d5,theta2]';
% MDP [25x1]
%   Minimal dynamic parameter vector (fixed base model)
%   see S5RPRRR7_convert_par2_MPV_fixb.m
% 
% Output:
% tauc [5x1]
%   joint torques required to compensate Coriolis and centrifugal load

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 19:04
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauc = S5RPRRR7_coriolisvecJ_fixb_mdp_slag_vp(qJ, qJD, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(9,1),zeros(25,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPRRR7_coriolisvecJ_fixb_mdp_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RPRRR7_coriolisvecJ_fixb_mdp_slag_vp: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5RPRRR7_coriolisvecJ_fixb_mdp_slag_vp: pkin has to be [9x1] (double)');
assert(isreal(MDP) && all(size(MDP) == [25 1]), ...
  'S5RPRRR7_coriolisvecJ_fixb_mdp_slag_vp: MDP has to be [25x1] (double)'); 

%% Symbolic Calculation
% From coriolisvec_joint_fixb_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 19:04:14
% EndTime: 2019-12-31 19:04:21
% DurationCPUTime: 2.96s
% Computational Cost: add. (1660->312), mult. (3969->440), div. (0->0), fcn. (2644->8), ass. (0->150)
t311 = sin(qJ(4));
t314 = cos(qJ(4));
t358 = t314 * qJD(3);
t312 = sin(qJ(3));
t370 = qJD(1) * t312;
t278 = t311 * t370 - t358;
t313 = cos(qJ(5));
t367 = qJD(3) * t311;
t280 = t314 * t370 + t367;
t310 = sin(qJ(5));
t388 = t280 * t310;
t229 = t313 * t278 + t388;
t315 = cos(qJ(3));
t369 = qJD(1) * t315;
t297 = -qJD(4) + t369;
t295 = -qJD(5) + t297;
t411 = t229 * t295;
t327 = t278 * t310 - t313 * t280;
t410 = t295 * t327;
t298 = sin(pkin(9)) * pkin(1) + pkin(6);
t290 = t298 * qJD(1);
t304 = t312 * qJD(2);
t256 = t315 * t290 + t304;
t248 = qJD(3) * pkin(7) + t256;
t299 = -cos(pkin(9)) * pkin(1) - pkin(2);
t273 = -pkin(3) * t315 - pkin(7) * t312 + t299;
t251 = t273 * qJD(1);
t385 = t311 * t251;
t221 = t248 * t314 + t385;
t216 = -pkin(8) * t278 + t221;
t362 = qJD(5) * t310;
t214 = t216 * t362;
t404 = qJD(2) * t315 - t312 * t290;
t247 = -qJD(3) * pkin(3) - t404;
t227 = pkin(4) * t278 + t247;
t409 = t227 * t229 + t214;
t357 = qJD(1) * qJD(3);
t346 = t315 * t357;
t364 = qJD(4) * t311;
t352 = t312 * t364;
t356 = qJD(3) * qJD(4);
t242 = -qJD(1) * t352 + (t346 + t356) * t314;
t249 = t404 * qJD(3);
t331 = pkin(3) * t312 - pkin(7) * t315;
t287 = t331 * qJD(3);
t272 = qJD(1) * t287;
t337 = t311 * t249 - t314 * t272;
t318 = -t221 * qJD(4) - t337;
t347 = t312 * t357;
t206 = pkin(4) * t347 - pkin(8) * t242 + t318;
t365 = qJD(3) * t315;
t350 = t311 * t365;
t363 = qJD(4) * t314;
t351 = t312 * t363;
t319 = t350 + t351;
t243 = t319 * qJD(1) + t311 * t356;
t354 = t314 * t249 + t251 * t363 + t311 * t272;
t321 = -t248 * t364 + t354;
t207 = -pkin(8) * t243 + t321;
t342 = t313 * t206 - t310 * t207;
t408 = t227 * t327 + t342;
t366 = qJD(3) * t312;
t344 = MDP(23) * t366;
t407 = qJD(1) * t344 + (-t229 ^ 2 + t327 ^ 2) * MDP(20) - t229 * t327 * MDP(19);
t406 = MDP(5) * t312;
t282 = t310 * t314 + t311 * t313;
t259 = t282 * t312;
t306 = t312 ^ 2;
t405 = (-t315 ^ 2 + t306) * MDP(6);
t403 = t315 * t358 - t352;
t402 = qJD(4) + qJD(5);
t340 = t242 * t310 + t313 * t243;
t209 = -t327 * qJD(5) + t340;
t401 = pkin(7) + pkin(8);
t400 = t209 * t315;
t220 = -t248 * t311 + t314 * t251;
t215 = -pkin(8) * t280 + t220;
t212 = -pkin(4) * t297 + t215;
t399 = t212 * t313;
t398 = t216 * t313;
t281 = t310 * t311 - t313 * t314;
t322 = t281 * t315;
t217 = -qJD(3) * t322 - t259 * t402;
t397 = t217 * t295;
t396 = t242 * t311;
t395 = t243 * t315;
t394 = t247 * t311;
t393 = t247 * t314;
t250 = qJD(3) * t304 + t290 * t365;
t392 = t250 * t311;
t391 = t250 * t314;
t390 = t278 * t297;
t389 = t280 * t297;
t387 = t297 * t314;
t386 = t298 * t311;
t384 = t311 * t312;
t383 = t311 * t315;
t382 = t312 * t314;
t316 = qJD(3) ^ 2;
t381 = t312 * t316;
t380 = t314 * t315;
t379 = t315 * t316;
t218 = -t362 * t384 + (t382 * t402 + t350) * t313 + t403 * t310;
t378 = t218 * t295 - t259 * t347;
t377 = qJD(1) * t322 - t281 * t402;
t376 = (-t369 + t402) * t282;
t284 = t331 * qJD(1);
t375 = t311 * t284 + t314 * t404;
t374 = t273 * t363 + t311 * t287;
t373 = t314 * t287 + t366 * t386;
t283 = t298 * t380;
t372 = t311 * t273 + t283;
t291 = qJD(1) * t299;
t361 = qJD(5) * t313;
t359 = t247 * qJD(4);
t355 = t313 * t242 - t310 * t243 - t278 * t361;
t353 = qJD(4) * t401;
t349 = t311 * t369;
t343 = MDP(16) * t366;
t208 = -t280 * t362 + t355;
t341 = -t208 * t315 - t327 * t366;
t339 = -t242 * t315 + t280 * t366;
t338 = t297 * t298 + t248;
t336 = t314 * t284 - t311 * t404;
t335 = qJD(5) * t212 + t207;
t334 = t297 * t352;
t333 = t297 * t351;
t332 = -t256 + (-t349 + t364) * pkin(4);
t293 = t401 * t314;
t323 = pkin(4) * t312 - pkin(8) * t380;
t330 = t323 * qJD(1) + qJD(5) * t293 + t314 * t353 + t336;
t292 = t401 * t311;
t329 = pkin(8) * t349 - qJD(5) * t292 - t311 * t353 - t375;
t204 = t212 * t310 + t398;
t262 = t314 * t273;
t224 = -pkin(8) * t382 + t262 + (-pkin(4) - t386) * t315;
t226 = -pkin(8) * t384 + t372;
t328 = t224 * t310 + t226 * t313;
t326 = qJD(1) * t306 - t297 * t315;
t324 = 0.2e1 * qJD(3) * t291;
t320 = t326 * t311;
t303 = -pkin(4) * t314 - pkin(3);
t265 = (pkin(4) * t311 + t298) * t312;
t260 = t281 * t312;
t235 = t319 * pkin(4) + t298 * t365;
t222 = pkin(4) * t243 + t250;
t211 = (-t312 * t358 - t315 * t364) * t298 - t319 * pkin(8) + t374;
t210 = t323 * qJD(3) + (-t283 + (pkin(8) * t312 - t273) * t311) * qJD(4) + t373;
t203 = -t216 * t310 + t399;
t1 = [0.2e1 * t346 * t406 - 0.2e1 * t357 * t405 + MDP(7) * t379 - MDP(8) * t381 + (-t298 * t379 + t312 * t324) * MDP(10) + (t298 * t381 + t315 * t324) * MDP(11) + (t242 * t382 + t280 * t403) * MDP(12) + ((-t278 * t314 - t280 * t311) * t365 + (-t396 - t243 * t314 + (t278 * t311 - t280 * t314) * qJD(4)) * t312) * MDP(13) + (t326 * t358 + t334 + t339) * MDP(14) + (t333 + t395 + (-t278 * t312 - t320) * qJD(3)) * MDP(15) + (-t297 - t369) * t343 + (-(-t273 * t364 + t373) * t297 + ((t278 * t298 + t394) * qJD(3) + (t338 * t314 + t385) * qJD(4) + t337) * t315 + (t314 * t359 + t298 * t243 + t392 + ((-t298 * t383 + t262) * qJD(1) + t220) * qJD(3)) * t312) * MDP(17) + (t374 * t297 + (-t338 * t364 + (t280 * t298 + t393) * qJD(3) + t354) * t315 + (-t311 * t359 + t298 * t242 + t391 + (-t372 * qJD(1) - t298 * t387 - t221) * qJD(3)) * t312) * MDP(18) + (-t208 * t260 - t217 * t327) * MDP(19) + (-t208 * t259 + t209 * t260 - t217 * t229 + t218 * t327) * MDP(20) + (-t260 * t347 + t341 - t397) * MDP(21) + (-t229 * t366 + t378 + t400) * MDP(22) + (-t295 - t369) * t344 + (-(t210 * t313 - t211 * t310) * t295 - t342 * t315 + t235 * t229 + t265 * t209 + t222 * t259 + t227 * t218 + (t204 * t315 + t328 * t295) * qJD(5) + ((t224 * t313 - t226 * t310) * qJD(1) + t203) * t366) * MDP(24) + (t265 * t208 - t214 * t315 + t227 * t217 - t222 * t260 - t235 * t327 + ((-qJD(5) * t226 + t210) * t295 + t206 * t315) * t310 + ((qJD(5) * t224 + t211) * t295 + t335 * t315) * t313 + (-t328 * qJD(1) - t204) * t366) * MDP(25); (t333 - t395) * MDP(17) + (-t334 + t339) * MDP(18) + (t378 - t400) * MDP(24) + (t341 + t397) * MDP(25) + (-MDP(10) * t312 - MDP(11) * t315) * t316 + (-t326 * MDP(18) * t314 - MDP(17) * t320 + (MDP(25) * qJD(1) * t260 + t278 * MDP(17) + t229 * MDP(24)) * t312) * qJD(3); (qJD(3) * t256 - t250) * MDP(10) - t291 * t369 * MDP(11) + (-t280 * t387 + t396) * MDP(12) + ((t242 + t390) * t314 + (-t243 + t389) * t311) * MDP(13) + (-t297 * t363 + (t297 * t380 + (-t280 + t367) * t312) * qJD(1)) * MDP(14) + (t297 * t364 + (-t297 * t383 + (t278 + t358) * t312) * qJD(1)) * MDP(15) + (-pkin(3) * t243 - t391 + t336 * t297 - t256 * t278 + (pkin(7) * t387 + t394) * qJD(4) + (-t220 * t312 + (-pkin(7) * t366 - t247 * t315) * t311) * qJD(1)) * MDP(17) + (-pkin(3) * t242 + t392 - t375 * t297 - t256 * t280 + (-pkin(7) * t297 * t311 + t393) * qJD(4) + (-t247 * t380 + (-pkin(7) * t358 + t221) * t312) * qJD(1)) * MDP(18) + (t208 * t282 - t327 * t377) * MDP(19) + (-t208 * t281 - t209 * t282 - t377 * t229 + t327 * t376) * MDP(20) + (t303 * t209 + t222 * t281 + t376 * t227 + t332 * t229) * MDP(24) + (t303 * t208 + t222 * t282 + t377 * t227 - t327 * t332) * MDP(25) + (-t377 * MDP(21) + t376 * MDP(22) + (t329 * t310 + t330 * t313) * MDP(24) + (-t330 * t310 + t329 * t313) * MDP(25)) * t295 + (-t291 * MDP(10) + t297 * MDP(16) + (qJD(3) * t282 + t327) * MDP(21) + (-qJD(3) * t281 + t229) * MDP(22) + t295 * MDP(23) + ((-t292 * t313 - t293 * t310) * qJD(3) - t203) * MDP(24) + (-(-t292 * t310 + t293 * t313) * qJD(3) + t204) * MDP(25)) * t370 + (-t315 * t406 + t405) * qJD(1) ^ 2; t280 * t278 * MDP(12) + (-t278 ^ 2 + t280 ^ 2) * MDP(13) + (t242 - t390) * MDP(14) + (-t243 - t389) * MDP(15) + qJD(1) * t343 + (-t221 * t297 - t247 * t280 + t318) * MDP(17) + (-t220 * t297 + t247 * t278 - t321) * MDP(18) + (t208 - t411) * MDP(21) + (-t209 + t410) * MDP(22) + ((-t215 * t310 - t398) * t295 - t204 * qJD(5) + (-t229 * t280 + t295 * t362 + t313 * t347) * pkin(4) + t408) * MDP(24) + ((t216 * t295 - t206) * t310 + (-t215 * t295 - t335) * t313 + (t280 * t327 + t295 * t361 - t310 * t347) * pkin(4) + t409) * MDP(25) + t407; (t355 - t411) * MDP(21) + (-t340 + t410) * MDP(22) + (-t204 * t295 + t408) * MDP(24) + (-t203 * t295 - t310 * t206 - t313 * t207 + t409) * MDP(25) + (-MDP(21) * t388 + t327 * MDP(22) - t204 * MDP(24) - MDP(25) * t399) * qJD(5) + t407;];
tauc = t1;
