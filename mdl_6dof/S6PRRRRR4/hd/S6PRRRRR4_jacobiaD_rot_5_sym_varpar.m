% Zeitableitung der rotatorischen Teilmatrix der analytischen Jacobi-Matrix für Segment Nr. 5 (0=Basis) von
% S6PRRRRR4
% Use Code from Maple symbolic Code Generation
%
% analytische Jacobi-Matrix: Differentieller Zusammenhang zwischen
% Endeffektorposition und verallgemeinerten Koordinaten.
% Zeitableitung der Winkeldarstellung des Endeffektors in Basis-Koordinaten
%
% Winkeldarstellung: Euler-XYZ-Winkel, rotx(alpha)*roty(beta)*rotz(gamma)
%
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [13x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,alpha3,d2,d3,d4,d5,d6,theta1]';
%
% Output:
% JaD_rot [3x6]
%   Zeitableitung der rotatorischen Teilmatrix der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-02-26 20:20
% Revision: d75aae1ac561373cd3be920984c3df29a1c2ecc8 (2019-02-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JaD_rot = S6PRRRRR4_jacobiaD_rot_5_sym_varpar(qJ, qJD, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(13,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRRRRR4_jacobiaD_rot_5_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6PRRRRR4_jacobiaD_rot_5_sym_varpar: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [13 1]), ...
  'S6PRRRRR4_jacobiaD_rot_5_sym_varpar: pkin has to be [13x1] (double)');

%% Symbolic Calculation
% From jacobiaD_rot_5_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 20:20:27
% EndTime: 2019-02-26 20:20:28
% DurationCPUTime: 1.22s
% Computational Cost: add. (6129->132), mult. (17990->254), div. (556->12), fcn. (22726->15), ass. (0->126)
t318 = sin(qJ(2));
t320 = cos(qJ(2));
t385 = cos(pkin(13));
t386 = cos(pkin(6));
t351 = t386 * t385;
t384 = sin(pkin(13));
t336 = -t384 * t318 + t320 * t351;
t299 = t336 * qJD(2);
t314 = sin(pkin(7));
t315 = sin(pkin(6));
t355 = t314 * t315 * t385;
t389 = -qJD(3) * t355 + t299;
t317 = sin(qJ(3));
t335 = -t318 * t351 - t320 * t384;
t330 = qJD(2) * t335;
t319 = cos(qJ(3));
t332 = t319 * t336;
t388 = qJD(3) * t332 + t317 * t330;
t316 = cos(pkin(7));
t331 = t336 * t317;
t327 = t316 * t331 - t319 * t335;
t370 = t316 * t319;
t253 = t327 * qJD(3) + t389 * t317 - t330 * t370;
t374 = t335 * t317;
t277 = -t316 * t332 + t319 * t355 - t374;
t275 = t277 ^ 2;
t366 = t319 * t320;
t369 = t317 * t318;
t341 = t316 * t366 - t369;
t361 = t314 * t386;
t290 = -t315 * t341 - t319 * t361;
t288 = 0.1e1 / t290 ^ 2;
t269 = t275 * t288 + 0.1e1;
t375 = t277 * t288;
t367 = t318 * t319;
t368 = t317 * t320;
t339 = t316 * t368 + t367;
t340 = t316 * t367 + t368;
t353 = qJD(3) * t361;
t273 = t317 * t353 + (qJD(2) * t340 + qJD(3) * t339) * t315;
t287 = 0.1e1 / t290;
t376 = t273 * t287 * t288;
t387 = -0.2e1 * (t253 * t375 - t275 * t376) / t269 ^ 2;
t270 = atan2(-t277, t290);
t265 = sin(t270);
t266 = cos(t270);
t247 = -t265 * t277 + t266 * t290;
t244 = 0.1e1 / t247;
t350 = t386 * t384;
t334 = t318 * t350 - t320 * t385;
t333 = t318 * t385 + t320 * t350;
t360 = t315 * t384;
t354 = t314 * t360;
t337 = -t316 * t333 + t354;
t281 = t317 * t337 - t319 * t334;
t292 = t314 * t333 + t316 * t360;
t313 = qJ(4) + qJ(5);
t310 = sin(t313);
t311 = cos(t313);
t262 = t281 * t311 + t292 * t310;
t258 = 0.1e1 / t262;
t245 = 0.1e1 / t247 ^ 2;
t259 = 0.1e1 / t262 ^ 2;
t267 = 0.1e1 / t269;
t237 = (-t253 * t287 + t273 * t375) * t267;
t349 = -t265 * t290 - t266 * t277;
t233 = t237 * t349 - t253 * t265 + t266 * t273;
t383 = t233 * t244 * t245;
t301 = t334 * qJD(2);
t312 = qJD(4) + qJD(5);
t348 = t281 * t312 + t301 * t314;
t300 = t333 * qJD(2);
t371 = t316 * t317;
t372 = t334 * t317;
t256 = t301 * t371 - t300 * t319 + (t319 * t337 + t372) * qJD(3);
t356 = t292 * t312 + t256;
t248 = t310 * t356 + t311 * t348;
t261 = t281 * t310 - t292 * t311;
t257 = t261 ^ 2;
t252 = t257 * t259 + 0.1e1;
t379 = t259 * t261;
t249 = -t310 * t348 + t311 * t356;
t380 = t249 * t258 * t259;
t382 = (t248 * t379 - t257 * t380) / t252 ^ 2;
t280 = -t319 * t354 + t333 * t370 - t372;
t381 = t245 * t280;
t378 = t265 * t280;
t377 = t266 * t280;
t373 = t334 * t314;
t276 = t280 ^ 2;
t243 = t245 * t276 + 0.1e1;
t255 = qJD(3) * t281 - t300 * t317 - t301 * t370;
t365 = 0.2e1 * (t255 * t381 - t276 * t383) / t243 ^ 2;
t364 = -0.2e1 * t382;
t363 = t261 * t380;
t362 = qJD(3) * t374;
t358 = -0.2e1 * t277 * t376;
t357 = 0.2e1 * t280 * t383;
t284 = -t317 * t333 - t334 * t370;
t352 = -qJD(3) * t284 + t300 * t371 + t301 * t319 - t312 * t373;
t285 = -t319 * t333 + t334 * t371;
t347 = t285 * t312 + t300 * t314;
t344 = -t310 * t258 + t311 * t379;
t279 = -t317 * t355 + t327;
t291 = t315 * t339 + t317 * t361;
t343 = -t279 * t287 + t291 * t375;
t283 = -t335 * t370 + t331;
t298 = t340 * t315;
t342 = -t283 * t287 + t298 * t375;
t338 = -t316 * t369 + t366;
t282 = (qJD(2) * t341 + qJD(3) * t338) * t315;
t274 = t319 * t353 + (qJD(2) * t338 + qJD(3) * t341) * t315;
t272 = t285 * t311 - t310 * t373;
t271 = t285 * t310 + t311 * t373;
t263 = t299 * t370 + t316 * t362 + t388;
t254 = t388 * t316 + t389 * t319 + t362;
t250 = 0.1e1 / t252;
t241 = 0.1e1 / t243;
t239 = t342 * t267;
t238 = t343 * t267;
t235 = t239 * t349 - t265 * t283 + t266 * t298;
t234 = t238 * t349 - t265 * t279 + t266 * t291;
t232 = t342 * t387 + (t298 * t358 - t263 * t287 + (t253 * t298 + t273 * t283 + t277 * t282) * t288) * t267;
t231 = t343 * t387 + (t291 * t358 - t254 * t287 + (t253 * t291 + t273 * t279 + t274 * t277) * t288) * t267;
t229 = t364 + 0.2e1 * (t248 * t259 * t250 + (-t250 * t380 - t259 * t382) * t261) * t261;
t1 = [0, t232, t231, 0, 0, 0; 0 (t235 * t381 - t244 * t284) * t365 + ((qJD(3) * t285 - t300 * t370 + t301 * t317) * t244 + t235 * t357 + (-t284 * t233 - t235 * t255 - (-t232 * t277 - t239 * t253 + t282 + (-t239 * t290 - t283) * t237) * t377 - (-t232 * t290 - t239 * t273 - t263 + (t239 * t277 - t298) * t237) * t378) * t245) * t241 (t234 * t381 - t244 * t281) * t365 + (t234 * t357 + t256 * t244 + (-t281 * t233 - t234 * t255 - (-t231 * t277 - t238 * t253 + t274 + (-t238 * t290 - t279) * t237) * t377 - (-t231 * t290 - t238 * t273 - t254 + (t238 * t277 - t291) * t237) * t378) * t245) * t241, 0, 0, 0; 0, 0.2e1 * (-t258 * t271 + t272 * t379) * t382 + ((t310 * t352 + t311 * t347) * t258 + 0.2e1 * t272 * t363 + (-t271 * t249 - (-t310 * t347 + t311 * t352) * t261 - t272 * t248) * t259) * t250, t344 * t280 * t364 + (t344 * t255 + ((-t258 * t312 - 0.2e1 * t363) * t311 + (t248 * t311 + (-t261 * t312 + t249) * t310) * t259) * t280) * t250, t229, t229, 0;];
JaD_rot  = t1;
