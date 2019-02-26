% Zeitableitung der rotatorischen Teilmatrix der analytischen Jacobi-Matrix für Segment Nr. 6 (0=Basis) von
% S6RRPRRP5
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
%   pkin=[a2,a3,a4,a5,a6,alpha2,d1,d2,d4,d5,theta3]';
%
% Output:
% JaD_rot [3x6]
%   Zeitableitung der rotatorischen Teilmatrix der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-02-26 21:48
% Revision: d75aae1ac561373cd3be920984c3df29a1c2ecc8 (2019-02-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JaD_rot = S6RRPRRP5_jacobiaD_rot_6_sym_varpar(qJ, qJD, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(11,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPRRP5_jacobiaD_rot_6_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRPRRP5_jacobiaD_rot_6_sym_varpar: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6RRPRRP5_jacobiaD_rot_6_sym_varpar: pkin has to be [11x1] (double)');

%% Symbolic Calculation
% From jacobiaD_rot_6_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 21:48:21
% EndTime: 2019-02-26 21:48:23
% DurationCPUTime: 1.87s
% Computational Cost: add. (8421->153), mult. (24081->307), div. (726->12), fcn. (31139->15), ass. (0->130)
t317 = sin(pkin(11));
t318 = cos(pkin(11));
t322 = sin(qJ(2));
t325 = cos(qJ(2));
t307 = t317 * t322 - t325 * t318;
t319 = cos(pkin(6));
t304 = t307 * t319;
t299 = qJD(2) * t304;
t347 = t317 * t325 + t318 * t322;
t306 = t347 * qJD(2);
t326 = cos(qJ(1));
t391 = sin(pkin(6));
t356 = t326 * t391;
t305 = t347 * t319;
t392 = sin(qJ(1));
t359 = t392 * t305;
t397 = -t392 * t306 - qJD(1) * t359 + (-qJD(1) * t307 - t299) * t326 - qJD(4) * t356;
t321 = sin(qJ(4));
t324 = cos(qJ(4));
t340 = -t326 * t305 + t392 * t307;
t277 = -t321 * t340 + t324 * t356;
t357 = t325 * t391;
t358 = t322 * t391;
t334 = -t317 * t357 - t318 * t358;
t293 = -t319 * t324 - t321 * t334;
t272 = atan2(-t277, t293);
t267 = sin(t272);
t268 = cos(t272);
t243 = -t267 * t277 + t268 * t293;
t241 = 0.1e1 / t243 ^ 2;
t341 = -t326 * t307 - t359;
t352 = t392 * t391;
t337 = -t321 * t341 + t324 * t352;
t276 = t337 ^ 2;
t239 = t241 * t276 + 0.1e1;
t283 = t321 * t352 + t324 * t341;
t333 = t340 * qJD(1) + t392 * t299 - t326 * t306;
t351 = qJD(1) * t356;
t247 = t283 * qJD(4) + t321 * t333 - t324 * t351;
t384 = t247 * t241;
t275 = t277 ^ 2;
t291 = 0.1e1 / t293 ^ 2;
t271 = t275 * t291 + 0.1e1;
t269 = 0.1e1 / t271;
t346 = qJD(1) * t352;
t370 = qJD(4) * t324;
t249 = t397 * t321 - t324 * t346 - t340 * t370;
t294 = t319 * t321 - t324 * t334;
t302 = -t317 * t358 + t318 * t357;
t298 = t302 * qJD(2);
t273 = t294 * qJD(4) + t298 * t321;
t290 = 0.1e1 / t293;
t376 = t277 * t291;
t345 = -t249 * t290 + t273 * t376;
t231 = t345 * t269;
t348 = -t267 * t293 - t268 * t277;
t226 = t348 * t231 - t267 * t249 + t268 * t273;
t240 = 0.1e1 / t243;
t242 = t240 * t241;
t389 = t226 * t242;
t368 = 0.2e1 * (-t276 * t389 - t337 * t384) / t239 ^ 2;
t396 = t273 * t291;
t285 = -t326 * t304 - t347 * t392;
t342 = -t285 * t290 + t302 * t376;
t395 = t321 * t342;
t250 = (qJD(4) * t340 + t346) * t321 + t397 * t324;
t288 = t392 * t304 - t326 * t347;
t320 = sin(qJ(5));
t323 = cos(qJ(5));
t259 = t283 * t323 - t288 * t320;
t253 = 0.1e1 / t259;
t254 = 0.1e1 / t259 ^ 2;
t394 = -0.2e1 * t277;
t393 = -0.2e1 * t337;
t248 = t337 * qJD(4) + t321 * t351 + t324 * t333;
t300 = t319 * t306;
t339 = t307 * qJD(2);
t263 = t285 * qJD(1) - t392 * t300 - t326 * t339;
t234 = t259 * qJD(5) + t248 * t320 - t263 * t323;
t258 = t283 * t320 + t288 * t323;
t252 = t258 ^ 2;
t246 = t252 * t254 + 0.1e1;
t382 = t254 * t258;
t369 = qJD(5) * t258;
t235 = t248 * t323 + t263 * t320 - t369;
t386 = t235 * t253 * t254;
t388 = (t234 * t382 - t252 * t386) / t246 ^ 2;
t378 = t290 * t396;
t387 = (t249 * t376 - t275 * t378) / t271 ^ 2;
t385 = t241 * t337;
t383 = t253 * t320;
t381 = t258 * t323;
t380 = t267 * t337;
t379 = t268 * t337;
t377 = t277 * t290;
t375 = t288 * t321;
t374 = t288 * t324;
t367 = -0.2e1 * t388;
t366 = 0.2e1 * t388;
t365 = -0.2e1 * t387;
t364 = t242 * t393;
t363 = t290 * t387;
t362 = t241 * t380;
t361 = t241 * t379;
t360 = t258 * t386;
t355 = 0.2e1 * t360;
t354 = t378 * t394;
t279 = -t321 * t356 - t324 * t340;
t349 = qJD(5) * t374 - t333;
t257 = -t279 * t323 + t285 * t320;
t256 = -t279 * t320 - t285 * t323;
t344 = t254 * t381 - t383;
t343 = -t279 * t290 + t294 * t376;
t338 = -t267 + (t268 * t377 + t267) * t269;
t335 = -qJD(4) * t375 + qJD(5) * t341 - t263 * t324;
t297 = t334 * qJD(2);
t274 = -t293 * qJD(4) + t298 * t324;
t265 = t288 * qJD(1) - t326 * t300 + t392 * t339;
t261 = t320 * t341 + t323 * t374;
t260 = t320 * t374 - t323 * t341;
t244 = 0.1e1 / t246;
t237 = 0.1e1 / t239;
t236 = t269 * t395;
t233 = t343 * t269;
t230 = t338 * t337;
t228 = (-t267 * t285 + t268 * t302) * t321 + t348 * t236;
t227 = t348 * t233 - t267 * t279 + t268 * t294;
t225 = t343 * t365 + (t294 * t354 - t250 * t290 + (t249 * t294 + t273 * t279 + t274 * t277) * t291) * t269;
t223 = t365 * t395 + (t342 * t370 + (t302 * t354 - t265 * t290 + (t249 * t302 + t273 * t285 + t277 * t297) * t291) * t321) * t269;
t1 = [t363 * t393 + (-t247 * t290 - t337 * t396) * t269, t223, 0, t225, 0, 0; t277 * t240 * t368 + (-t249 * t240 + (t226 * t277 + t230 * t247) * t241) * t237 - (-t230 * t241 * t368 + (-0.2e1 * t230 * t389 + (-t231 * t269 * t377 + t365) * t362 + (t363 * t394 - t231 + (t231 - t345) * t269) * t361 - t338 * t384) * t237) * t337 (-t228 * t385 - t240 * t375) * t368 + (-t228 * t384 + (-t263 * t321 + t288 * t370) * t240 + (t228 * t364 - t241 * t375) * t226 + (t302 * t370 - t223 * t277 - t236 * t249 + t297 * t321 + (-t236 * t293 - t285 * t321) * t231) * t361 + (-t285 * t370 - t223 * t293 - t236 * t273 - t265 * t321 + (t236 * t277 - t302 * t321) * t231) * t362) * t237, 0 (-t227 * t385 - t240 * t283) * t368 + (t227 * t226 * t364 + t248 * t240 + (-t283 * t226 - t227 * t247 + (-t225 * t277 - t233 * t249 + t274 + (-t233 * t293 - t279) * t231) * t379 + (-t225 * t293 - t233 * t273 - t250 + (t233 * t277 - t294) * t231) * t380) * t241) * t237, 0, 0; (-t253 * t256 + t257 * t382) * t366 + ((t257 * qJD(5) - t250 * t320 - t265 * t323) * t253 + t257 * t355 + (-t256 * t235 - (-t256 * qJD(5) - t250 * t323 + t265 * t320) * t258 - t257 * t234) * t254) * t244 (-t253 * t260 + t261 * t382) * t366 + (t261 * t355 + t349 * t253 * t323 + t335 * t383 + (t349 * t258 * t320 - t261 * t234 - t260 * t235 - t335 * t381) * t254) * t244, 0, -t344 * t337 * t367 + (t344 * t247 - ((-qJD(5) * t253 - 0.2e1 * t360) * t323 + (t234 * t323 + (t235 - t369) * t320) * t254) * t337) * t244, t367 + 0.2e1 * (t234 * t254 * t244 + (-t244 * t386 - t254 * t388) * t258) * t258, 0;];
JaD_rot  = t1;
