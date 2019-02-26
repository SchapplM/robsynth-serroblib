% Zeitableitung der rotatorischen Teilmatrix der analytischen Jacobi-Matrix für Segment Nr. 6 (0=Basis) von
% S6RRRPPR9
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
% pkin [11x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d1,d2,d3,d6,theta4]';
%
% Output:
% JaD_rot [3x6]
%   Zeitableitung der rotatorischen Teilmatrix der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-02-26 22:08
% Revision: d75aae1ac561373cd3be920984c3df29a1c2ecc8 (2019-02-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JaD_rot = S6RRRPPR9_jacobiaD_rot_6_sym_varpar(qJ, qJD, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(11,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRPPR9_jacobiaD_rot_6_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRRPPR9_jacobiaD_rot_6_sym_varpar: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6RRRPPR9_jacobiaD_rot_6_sym_varpar: pkin has to be [11x1] (double)');

%% Symbolic Calculation
% From jacobiaD_rot_6_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 22:08:10
% EndTime: 2019-02-26 22:08:12
% DurationCPUTime: 1.77s
% Computational Cost: add. (5732->165), mult. (16796->329), div. (720->12), fcn. (21210->15), ass. (0->140)
t319 = cos(pkin(6));
t322 = sin(qJ(2));
t396 = sin(qJ(1));
t367 = t396 * t322;
t359 = t319 * t367;
t365 = qJD(2) * t396;
t324 = cos(qJ(2));
t325 = cos(qJ(1));
t377 = t325 * t324;
t286 = -qJD(1) * t359 - t322 * t365 + (qJD(2) * t319 + qJD(1)) * t377;
t321 = sin(qJ(3));
t395 = sin(pkin(6));
t314 = t325 * t395 * t321;
t366 = t396 * t324;
t378 = t325 * t322;
t340 = -t319 * t378 - t366;
t355 = t396 * t395;
t397 = cos(qJ(3));
t341 = t397 * t355;
t364 = qJD(3) * t397;
t265 = -qJD(1) * t341 - qJD(3) * t314 + t286 * t321 - t340 * t364;
t356 = t397 * t395;
t345 = t325 * t356;
t335 = -t321 * t340 + t345;
t287 = t335 ^ 2;
t363 = t322 * t395;
t334 = -t319 * t397 + t321 * t363;
t299 = 0.1e1 / t334 ^ 2;
t279 = t287 * t299 + 0.1e1;
t382 = t335 * t299;
t298 = 0.1e1 / t334;
t302 = -t319 * t321 - t322 * t356;
t362 = t324 * t395;
t357 = t321 * t362;
t289 = -qJD(2) * t357 + t302 * qJD(3);
t400 = t289 * t299;
t384 = t298 * t400;
t390 = (t265 * t382 + t287 * t384) / t279 ^ 2;
t401 = -0.2e1 * t390;
t358 = t335 * t362;
t369 = t319 * t377;
t303 = -t367 + t369;
t380 = t298 * t303;
t399 = t321 * (t299 * t358 - t380);
t306 = -t359 + t377;
t332 = -t306 * t321 + t341;
t317 = sin(pkin(11));
t318 = cos(pkin(11));
t320 = sin(qJ(6));
t323 = cos(qJ(6));
t348 = t317 * t323 - t318 * t320;
t398 = t348 * t332;
t344 = qJD(3) * t356;
t266 = (qJD(1) * t355 + qJD(3) * t340) * t321 + t286 * t397 - t325 * t344;
t280 = atan2(t335, -t334);
t275 = sin(t280);
t276 = cos(t280);
t246 = t275 * t335 - t276 * t334;
t243 = 0.1e1 / t246;
t297 = t306 * t397 + t321 * t355;
t339 = -t319 * t366 - t378;
t273 = t297 * t317 + t318 * t339;
t274 = t297 * t318 - t317 * t339;
t258 = t273 * t320 + t274 * t323;
t252 = 0.1e1 / t258;
t244 = 0.1e1 / t246 ^ 2;
t253 = 0.1e1 / t258 ^ 2;
t288 = t332 ^ 2;
t242 = t288 * t244 + 0.1e1;
t284 = t340 * qJD(1) + t339 * qJD(2);
t263 = -qJD(1) * t345 + t297 * qJD(3) + t284 * t321;
t387 = t263 * t244;
t277 = 0.1e1 / t279;
t343 = -t265 * t298 - t289 * t382;
t233 = t343 * t277;
t350 = t275 * t334 + t276 * t335;
t227 = t350 * t233 + t275 * t265 + t276 * t289;
t245 = t243 * t244;
t392 = t227 * t245;
t394 = (-t288 * t392 - t332 * t387) / t242 ^ 2;
t264 = qJD(1) * t314 + t332 * qJD(3) + t284 * t397;
t283 = -qJD(1) * t369 - qJD(2) * t377 + (t396 * qJD(1) + t319 * t365) * t322;
t247 = t264 * t317 + t283 * t318;
t248 = t264 * t318 - t283 * t317;
t230 = t258 * qJD(6) - t247 * t323 + t248 * t320;
t351 = t273 * t323 - t274 * t320;
t251 = t351 ^ 2;
t236 = t251 * t253 + 0.1e1;
t388 = t253 * t351;
t231 = t351 * qJD(6) + t247 * t320 + t248 * t323;
t391 = t231 * t252 * t253;
t393 = (-t230 * t388 - t251 * t391) / t236 ^ 2;
t389 = t244 * t332;
t386 = t275 * t332;
t385 = t276 * t332;
t383 = t335 * t298;
t381 = t335 * t302;
t379 = t339 * t321;
t376 = -0.2e1 * t394;
t375 = 0.2e1 * t394;
t374 = 0.2e1 * t393;
t373 = 0.2e1 * t390;
t372 = 0.2e1 * t245 * t332;
t371 = t244 * t386;
t370 = t244 * t385;
t368 = t339 * t397;
t361 = -0.2e1 * t351 * t391;
t360 = t298 * t401;
t354 = qJD(2) * t363;
t292 = -t340 * t397 - t314;
t271 = -t292 * t317 - t303 * t318;
t272 = -t292 * t318 + t303 * t317;
t352 = t271 * t323 - t272 * t320;
t256 = t271 * t320 + t272 * t323;
t281 = -t306 * t318 + t317 * t368;
t282 = t306 * t317 + t318 * t368;
t349 = t281 * t323 - t282 * t320;
t262 = t281 * t320 + t282 * t323;
t347 = -t317 * t320 - t318 * t323;
t346 = t324 * t356;
t342 = t292 * t298 + t299 * t381;
t269 = t347 * t332;
t338 = -qJD(3) * t379 + t397 * t283;
t337 = t275 + (-t276 * t383 - t275) * t277;
t290 = -qJD(2) * t346 + t334 * qJD(3);
t285 = t339 * qJD(1) + t340 * qJD(2);
t260 = t284 * t317 + t338 * t318;
t259 = -t284 * t318 + t338 * t317;
t250 = -t266 * t318 + t285 * t317;
t249 = -t266 * t317 - t285 * t318;
t240 = 0.1e1 / t242;
t239 = t277 * t399;
t238 = t342 * t277;
t234 = 0.1e1 / t236;
t232 = t337 * t332;
t229 = (t275 * t303 - t276 * t362) * t321 + t350 * t239;
t228 = -t350 * t238 + t275 * t292 + t276 * t302;
t225 = t342 * t373 + (-0.2e1 * t381 * t384 - t266 * t298 + (-t265 * t302 - t289 * t292 - t290 * t335) * t299) * t277;
t223 = t399 * t401 + ((t346 * t382 - t397 * t380) * qJD(3) + (0.2e1 * t358 * t384 - t285 * t298 + (t265 * t362 - t289 * t303 - t335 * t354) * t299) * t321) * t277;
t1 = [t332 * t360 + (-t263 * t298 + t332 * t400) * t277, t223, t225, 0, 0, 0; t335 * t243 * t376 + (t265 * t243 + (-t227 * t335 - t232 * t263) * t244) * t240 - (-t232 * t244 * t376 + (0.2e1 * t232 * t392 - (t233 * t277 * t383 + t373) * t371 - (-t335 * t360 + t233 + (-t233 + t343) * t277) * t370 + t337 * t387) * t240) * t332 (t229 * t389 + t243 * t379) * t375 + (t229 * t387 + (-t283 * t321 - t339 * t364) * t243 + (t229 * t372 + t244 * t379) * t227 - (-t324 * t344 + t321 * t354 + t223 * t335 + t239 * t265 + (t239 * t334 + t303 * t321) * t233) * t370 - (t303 * t364 + t223 * t334 - t239 * t289 + t285 * t321 + (-t239 * t335 + t357) * t233) * t371) * t240 (t228 * t389 + t243 * t297) * t375 + (t228 * t227 * t372 - t264 * t243 + (t297 * t227 + t228 * t263 - (t225 * t335 - t238 * t265 + t290 + (-t238 * t334 + t292) * t233) * t385 - (t225 * t334 + t238 * t289 + t266 + (t238 * t335 - t302) * t233) * t386) * t244) * t240, 0, 0, 0; (t252 * t352 - t256 * t388) * t374 + ((t256 * qJD(6) - t249 * t323 + t250 * t320) * t252 + t256 * t361 + (t352 * t231 + (t352 * qJD(6) + t249 * t320 + t250 * t323) * t351 - t256 * t230) * t253) * t234 (t252 * t349 - t262 * t388) * t374 + ((t262 * qJD(6) - t259 * t323 + t260 * t320) * t252 + t262 * t361 + (t349 * t231 + (t349 * qJD(6) + t259 * t320 + t260 * t323) * t351 - t262 * t230) * t253) * t234 (t252 * t398 + t269 * t388) * t374 + ((-qJD(6) * t269 + t348 * t263) * t252 - t269 * t361 + (t398 * t231 + (qJD(6) * t398 + t347 * t263) * t351 + t269 * t230) * t253) * t234, 0, 0, -0.2e1 * t393 - 0.2e1 * (t230 * t253 * t234 - (-t234 * t391 - t253 * t393) * t351) * t351;];
JaD_rot  = t1;
