% Zeitableitung der rotatorischen Teilmatrix der analytischen Jacobi-Matrix für Segment Nr. 6 (0=Basis) von
% S6RRPRPR6
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
%   pkin=[a2,a3,a4,a5,a6,alpha2,d1,d2,d4,d6,theta3]';
%
% Output:
% JaD_rot [3x6]
%   Zeitableitung der rotatorischen Teilmatrix der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-02-26 21:40
% Revision: d75aae1ac561373cd3be920984c3df29a1c2ecc8 (2019-02-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JaD_rot = S6RRPRPR6_jacobiaD_rot_6_sym_varpar(qJ, qJD, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(11,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPRPR6_jacobiaD_rot_6_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRPRPR6_jacobiaD_rot_6_sym_varpar: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6RRPRPR6_jacobiaD_rot_6_sym_varpar: pkin has to be [11x1] (double)');

%% Symbolic Calculation
% From jacobiaD_rot_6_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 21:40:37
% EndTime: 2019-02-26 21:40:40
% DurationCPUTime: 1.91s
% Computational Cost: add. (8421->154), mult. (24081->309), div. (726->12), fcn. (31139->15), ass. (0->130)
t310 = sin(pkin(11));
t312 = cos(pkin(11));
t316 = sin(qJ(2));
t320 = cos(qJ(2));
t301 = t316 * t310 - t320 * t312;
t313 = cos(pkin(6));
t332 = t301 * t313;
t294 = qJD(2) * t332;
t337 = t320 * t310 + t316 * t312;
t299 = t337 * t313;
t300 = t337 * qJD(2);
t317 = sin(qJ(1));
t321 = cos(qJ(1));
t359 = qJD(1) * t317;
t311 = sin(pkin(6));
t365 = t311 * t321;
t388 = -t317 * t300 - t299 * t359 + (-qJD(1) * t301 - t294) * t321 - qJD(4) * t365;
t315 = sin(qJ(4));
t338 = -t317 * t299 - t321 * t301;
t319 = cos(qJ(4));
t366 = t311 * t319;
t276 = t315 * t338 - t317 * t366;
t283 = t317 * t332 - t321 * t337;
t314 = sin(qJ(6));
t318 = cos(qJ(6));
t339 = t276 * t318 + t283 * t314;
t387 = t339 * qJD(6);
t281 = t321 * t299 - t317 * t301;
t363 = t315 * t321;
t272 = t281 * t319 - t311 * t363;
t298 = t337 * t311;
t289 = t298 * t319 + t313 * t315;
t265 = atan2(-t272, t289);
t260 = sin(t265);
t261 = cos(t265);
t236 = -t260 * t272 + t261 * t289;
t234 = 0.1e1 / t236 ^ 2;
t277 = t317 * t311 * t315 + t319 * t338;
t270 = t277 ^ 2;
t232 = t270 * t234 + 0.1e1;
t328 = -t281 * qJD(1) + t317 * t294 - t321 * t300;
t356 = qJD(4) * t319;
t357 = qJD(4) * t315;
t241 = t328 * t319 - t338 * t357 + (qJD(1) * t363 + t317 * t356) * t311;
t376 = t241 * t234;
t269 = t272 ^ 2;
t286 = 0.1e1 / t289 ^ 2;
t264 = t269 * t286 + 0.1e1;
t262 = 0.1e1 / t264;
t343 = t388 * t319;
t348 = t311 * t359;
t243 = -t281 * t357 + t315 * t348 + t343;
t288 = -t298 * t315 + t313 * t319;
t297 = t301 * t311;
t293 = qJD(2) * t297;
t267 = t288 * qJD(4) - t293 * t319;
t285 = 0.1e1 / t289;
t370 = t272 * t286;
t336 = -t243 * t285 + t267 * t370;
t224 = t336 * t262;
t341 = -t260 * t289 - t261 * t272;
t219 = t224 * t341 - t260 * t243 + t261 * t267;
t233 = 0.1e1 / t236;
t235 = t233 * t234;
t381 = t219 * t235;
t355 = 0.2e1 * (-t270 * t381 + t277 * t376) / t232 ^ 2;
t386 = t267 * t286;
t280 = -t317 * t337 - t321 * t332;
t333 = -t280 * t285 - t297 * t370;
t385 = t319 * t333;
t368 = t283 * t318;
t252 = t276 * t314 - t368;
t246 = 0.1e1 / t252;
t247 = 0.1e1 / t252 ^ 2;
t384 = -0.2e1 * t272;
t383 = 0.2e1 * t277;
t347 = qJD(1) * t366;
t240 = t277 * qJD(4) + t315 * t328 - t321 * t347;
t295 = t313 * t300;
t331 = t301 * qJD(2);
t256 = t280 * qJD(1) - t317 * t295 - t321 * t331;
t227 = t252 * qJD(6) - t240 * t318 + t256 * t314;
t245 = t339 ^ 2;
t239 = t245 * t247 + 0.1e1;
t375 = t247 * t339;
t228 = t240 * t314 + t256 * t318 + t387;
t378 = t228 * t246 * t247;
t380 = (-t227 * t375 - t245 * t378) / t239 ^ 2;
t372 = t285 * t386;
t379 = (t243 * t370 - t269 * t372) / t264 ^ 2;
t377 = t234 * t277;
t374 = t260 * t277;
t373 = t261 * t277;
t371 = t272 * t285;
t369 = t283 * t315;
t367 = t283 * t319;
t364 = t314 * t339;
t361 = t318 * t246;
t354 = 0.2e1 * t380;
t353 = -0.2e1 * t379;
t352 = t235 * t383;
t351 = t285 * t379;
t350 = t234 * t374;
t349 = t234 * t373;
t345 = -0.2e1 * t339 * t378;
t344 = t372 * t384;
t342 = qJD(6) * t369 + t328;
t271 = t281 * t315 + t319 * t365;
t340 = -t271 * t318 - t280 * t314;
t250 = -t271 * t314 + t280 * t318;
t335 = -t247 * t364 + t361;
t334 = t271 * t285 + t288 * t370;
t330 = -t260 + (t261 * t371 + t260) * t262;
t242 = t281 * t356 + t315 * t388 - t317 * t347;
t329 = -qJD(6) * t338 - t256 * t315 + t283 * t356;
t292 = t311 * t300;
t266 = -t289 * qJD(4) + t293 * t315;
t258 = t283 * qJD(1) - t321 * t295 + t317 * t331;
t254 = t314 * t369 + t318 * t338;
t253 = t314 * t338 - t315 * t368;
t237 = 0.1e1 / t239;
t230 = 0.1e1 / t232;
t229 = t262 * t385;
t226 = t334 * t262;
t223 = t330 * t277;
t221 = (-t260 * t280 - t261 * t297) * t319 + t341 * t229;
t220 = t226 * t341 + t260 * t271 + t261 * t288;
t218 = t334 * t353 + (t288 * t344 + t242 * t285 + (t243 * t288 + t266 * t272 - t267 * t271) * t286) * t262;
t216 = t353 * t385 + (-t333 * t357 + (-t297 * t344 - t258 * t285 + (-t243 * t297 + t267 * t280 - t272 * t292) * t286) * t319) * t262;
t1 = [t351 * t383 + (-t241 * t285 + t277 * t386) * t262, t216, 0, t218, 0, 0; t272 * t233 * t355 + (((qJD(4) * t281 - t348) * t315 - t343) * t233 + (t219 * t272 - t223 * t241) * t234) * t230 + (t223 * t234 * t355 + (0.2e1 * t223 * t381 - (-t224 * t262 * t371 + t353) * t350 - (t351 * t384 - t224 + (t224 - t336) * t262) * t349 - t330 * t376) * t230) * t277 (t221 * t377 - t233 * t367) * t355 + (-t221 * t376 + (-t256 * t319 - t283 * t357) * t233 + (t221 * t352 - t234 * t367) * t219 - (t297 * t357 - t216 * t272 - t229 * t243 - t292 * t319 + (-t229 * t289 - t280 * t319) * t224) * t349 - (t280 * t357 - t216 * t289 - t229 * t267 - t258 * t319 + (t229 * t272 + t297 * t319) * t224) * t350) * t230, 0 (t220 * t377 + t233 * t276) * t355 + (t220 * t219 * t352 - t240 * t233 + (t276 * t219 - t220 * t241 - (-t218 * t272 - t226 * t243 + t266 + (-t226 * t289 + t271) * t224) * t373 - (-t218 * t289 - t226 * t267 + t242 + (t226 * t272 - t288) * t224) * t374) * t234) * t230, 0, 0; (t246 * t340 - t250 * t375) * t354 + ((qJD(6) * t250 + t242 * t318 + t258 * t314) * t246 + t250 * t345 + (t340 * t228 + (qJD(6) * t340 - t242 * t314 + t258 * t318) * t339 - t250 * t227) * t247) * t237 (-t246 * t253 - t254 * t375) * t354 + (t254 * t345 + t342 * t246 * t314 - t329 * t361 + (t318 * t339 * t342 - t254 * t227 - t253 * t228 + t329 * t364) * t247) * t237, 0, t335 * t277 * t354 + (-t335 * t241 + ((qJD(6) * t246 + t345) * t314 + (-t227 * t314 + (t228 + t387) * t318) * t247) * t277) * t237, 0, -0.2e1 * t380 - 0.2e1 * (t227 * t247 * t237 - (-t237 * t378 - t247 * t380) * t339) * t339;];
JaD_rot  = t1;
