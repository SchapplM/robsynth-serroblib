% Zeitableitung der rotatorischen Teilmatrix der analytischen Jacobi-Matrix für Segment Nr. 5 (0=Basis) von
% S6PRRRPR8
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
% pkin [12x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,alpha3,d2,d3,d4,d6,theta1]';
%
% Output:
% JaD_rot [3x6]
%   Zeitableitung der rotatorischen Teilmatrix der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-02-26 20:14
% Revision: d75aae1ac561373cd3be920984c3df29a1c2ecc8 (2019-02-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JaD_rot = S6PRRRPR8_jacobiaD_rot_5_sym_varpar(qJ, qJD, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(12,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRRRPR8_jacobiaD_rot_5_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6PRRRPR8_jacobiaD_rot_5_sym_varpar: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [12 1]), ...
  'S6PRRRPR8_jacobiaD_rot_5_sym_varpar: pkin has to be [12x1] (double)');

%% Symbolic Calculation
% From jacobiaD_rot_5_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 20:14:42
% EndTime: 2019-02-26 20:14:44
% DurationCPUTime: 2.22s
% Computational Cost: add. (12413->176), mult. (37424->340), div. (781->12), fcn. (48146->15), ass. (0->138)
t318 = sin(qJ(3));
t321 = cos(qJ(3));
t314 = sin(pkin(7));
t316 = cos(pkin(7));
t319 = sin(qJ(2));
t322 = cos(qJ(2));
t385 = cos(pkin(12));
t386 = cos(pkin(6));
t347 = t386 * t385;
t384 = sin(pkin(12));
t334 = t384 * t319 - t322 * t347;
t315 = sin(pkin(6));
t355 = t315 * t385;
t329 = -t314 * t355 - t334 * t316;
t333 = -t319 * t347 - t384 * t322;
t285 = t318 * t333 + t329 * t321;
t303 = t334 * qJD(2);
t304 = t333 * qJD(2);
t369 = t316 * t318;
t261 = t285 * qJD(3) - t303 * t321 + t304 * t369;
t286 = t329 * t318 - t321 * t333;
t317 = sin(qJ(4));
t320 = cos(qJ(4));
t330 = t334 * t314 - t316 * t355;
t271 = t286 * t320 + t330 * t317;
t371 = t314 * t320;
t244 = t271 * qJD(4) + t261 * t317 + t304 * t371;
t269 = t286 * t317 - t330 * t320;
t266 = t269 ^ 2;
t365 = t319 * t321;
t366 = t318 * t322;
t337 = t316 * t366 + t365;
t356 = t314 * t386;
t298 = t337 * t315 + t318 * t356;
t370 = t314 * t322;
t307 = -t315 * t370 + t386 * t316;
t289 = t298 * t317 - t307 * t320;
t283 = 0.1e1 / t289 ^ 2;
t256 = t266 * t283 + 0.1e1;
t253 = 0.1e1 / t256;
t364 = t321 * t322;
t367 = t318 * t319;
t336 = -t316 * t367 + t364;
t339 = t316 * t364 - t367;
t348 = qJD(3) * t356;
t278 = t321 * t348 + (t336 * qJD(2) + t339 * qJD(3)) * t315;
t290 = t298 * t320 + t307 * t317;
t372 = t314 * t319;
t357 = t315 * t372;
t350 = qJD(2) * t357;
t258 = t290 * qJD(4) + t278 * t317 - t320 * t350;
t282 = 0.1e1 / t289;
t379 = t269 * t283;
t231 = (-t244 * t282 + t258 * t379) * t253;
t257 = atan2(-t269, t289);
t249 = sin(t257);
t250 = cos(t257);
t345 = -t249 * t289 - t250 * t269;
t227 = t345 * t231 - t249 * t244 + t250 * t258;
t243 = -t249 * t269 + t250 * t289;
t240 = 0.1e1 / t243;
t241 = 0.1e1 / t243 ^ 2;
t389 = t227 * t240 * t241;
t346 = t386 * t384;
t332 = t319 * t346 - t385 * t322;
t331 = t385 * t319 + t322 * t346;
t354 = t315 * t384;
t349 = t314 * t354;
t335 = -t316 * t331 + t349;
t288 = t335 * t318 - t321 * t332;
t299 = t314 * t331 + t316 * t354;
t272 = t288 * t317 - t299 * t320;
t352 = 0.2e1 * t272 * t389;
t297 = t339 * t315 + t321 * t356;
t341 = -t282 * t285 + t297 * t379;
t388 = t317 * t341;
t380 = t258 * t282 * t283;
t387 = -0.2e1 * (t244 * t379 - t266 * t380) / t256 ^ 2;
t368 = t316 * t321;
t374 = t332 * t318;
t287 = -t321 * t349 + t331 * t368 - t374;
t279 = 0.1e1 / t287;
t280 = 0.1e1 / t287 ^ 2;
t383 = t241 * t272;
t382 = t249 * t272;
t381 = t250 * t272;
t273 = t288 * t320 + t299 * t317;
t378 = t273 * t280;
t377 = t279 * t287;
t376 = t287 * t317;
t373 = t314 * t317;
t363 = qJD(4) * t317;
t362 = qJD(4) * t320;
t267 = t272 ^ 2;
t239 = t241 * t267 + 0.1e1;
t305 = t331 * qJD(2);
t306 = t332 * qJD(2);
t263 = t306 * t369 - t305 * t321 + (t335 * t321 + t374) * qJD(3);
t246 = t273 * qJD(4) + t263 * t317 + t306 * t371;
t361 = 0.2e1 * (t246 * t383 - t267 * t389) / t239 ^ 2;
t247 = -t272 * qJD(4) + t263 * t320 - t306 * t373;
t268 = t273 ^ 2;
t255 = t268 * t280 + 0.1e1;
t262 = t288 * qJD(3) - t305 * t318 - t306 * t368;
t281 = t279 * t280;
t360 = 0.2e1 * (-t262 * t268 * t281 + t247 * t378) / t255 ^ 2;
t358 = 0.2e1 * t273 * t281;
t351 = -0.2e1 * t269 * t380;
t343 = -t271 * t282 + t290 * t379;
t291 = -t334 * t321 + t333 * t369;
t274 = t291 * t317 + t333 * t371;
t302 = t336 * t315;
t294 = t302 * t317 - t320 * t357;
t342 = -t274 * t282 + t294 * t379;
t293 = -t321 * t331 + t332 * t369;
t340 = -t293 * t317 - t332 * t371;
t276 = t293 * t320 - t332 * t373;
t292 = -t318 * t331 - t332 * t368;
t338 = -t316 * t365 - t366;
t277 = -t318 * t348 + (t338 * qJD(2) - t337 * qJD(3)) * t315;
t265 = -t292 * qJD(3) + t305 * t369 + t306 * t321;
t264 = t302 * t362 + ((t338 * qJD(3) + qJD(4) * t372) * t317 + (-t337 * t317 - t320 * t370) * qJD(2)) * t315;
t260 = -t286 * qJD(3) + t303 * t318 + t304 * t368;
t259 = -t289 * qJD(4) + t278 * t320 + t317 * t350;
t251 = 0.1e1 / t255;
t248 = (t303 * t369 + t304 * t321 + (t334 * t318 + t333 * t368) * qJD(3)) * t317 + t291 * t362 + t303 * t371 - t333 * t314 * t363;
t245 = -t269 * qJD(4) + t261 * t320 - t304 * t373;
t237 = 0.1e1 / t239;
t236 = t253 * t388;
t235 = t342 * t253;
t233 = t343 * t253;
t230 = (-t249 * t285 + t250 * t297) * t317 + t345 * t236;
t229 = t345 * t235 - t249 * t274 + t250 * t294;
t228 = t345 * t233 - t249 * t271 + t250 * t290;
t226 = t342 * t387 + (t294 * t351 - t248 * t282 + (t244 * t294 + t258 * t274 + t264 * t269) * t283) * t253;
t224 = t343 * t387 + (t290 * t351 - t245 * t282 + (t244 * t290 + t258 * t271 + t259 * t269) * t283) * t253;
t223 = t387 * t388 + (t341 * t362 + (t297 * t351 - t260 * t282 + (t244 * t297 + t258 * t285 + t269 * t277) * t283) * t317) * t253;
t1 = [0, t226, t223, t224, 0, 0; 0 (t229 * t383 + t240 * t340) * t361 + ((t276 * qJD(4) + t265 * t317 + t305 * t371) * t240 + t229 * t352 + (t340 * t227 - t229 * t246 - (-t226 * t269 - t235 * t244 + t264 + (-t235 * t289 - t274) * t231) * t381 - (-t226 * t289 - t235 * t258 - t248 + (t235 * t269 - t294) * t231) * t382) * t241) * t237 (t230 * t383 + t240 * t376) * t361 + ((-t262 * t317 - t287 * t362) * t240 + t230 * t352 + (-t230 * t246 + t376 * t227 - (t297 * t362 - t223 * t269 - t236 * t244 + t277 * t317 + (-t236 * t289 - t285 * t317) * t231) * t381 - (-t285 * t362 - t223 * t289 - t236 * t258 - t260 * t317 + (t236 * t269 - t297 * t317) * t231) * t382) * t241) * t237 (t228 * t383 - t240 * t273) * t361 + (t228 * t352 + t247 * t240 + (-t273 * t227 - t228 * t246 - (-t224 * t269 - t233 * t244 + t259 + (-t233 * t289 - t271) * t231) * t381 - (-t224 * t289 - t233 * t258 - t245 + (t233 * t269 - t290) * t231) * t382) * t241) * t237, 0, 0; 0 (-t276 * t279 + t292 * t378) * t360 + ((t340 * qJD(4) + t265 * t320 - t305 * t373) * t279 + t292 * t262 * t358 + (-t276 * t262 - (t293 * qJD(3) - t305 * t368 + t306 * t318) * t273 - t292 * t247) * t280) * t251 (t288 * t378 + t320 * t377) * t360 + (t363 * t377 + (-t247 * t288 - t263 * t273) * t280 + (t288 * t358 + (t280 * t287 - t279) * t320) * t262) * t251, t272 * t279 * t360 + (t262 * t272 * t280 - t246 * t279) * t251, 0, 0;];
JaD_rot  = t1;
