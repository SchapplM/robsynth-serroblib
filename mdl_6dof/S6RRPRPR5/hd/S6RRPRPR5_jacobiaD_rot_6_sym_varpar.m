% Zeitableitung der rotatorischen Teilmatrix der analytischen Jacobi-Matrix für Segment Nr. 6 (0=Basis) von
% S6RRPRPR5
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
%   pkin=[a2,a3,a4,a5,a6,alpha2,d1,d2,d4,d6,theta3,theta5]';
%
% Output:
% JaD_rot [3x6]
%   Zeitableitung der rotatorischen Teilmatrix der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-02-26 21:40
% Revision: d75aae1ac561373cd3be920984c3df29a1c2ecc8 (2019-02-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JaD_rot = S6RRPRPR5_jacobiaD_rot_6_sym_varpar(qJ, qJD, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(12,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPRPR5_jacobiaD_rot_6_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRPRPR5_jacobiaD_rot_6_sym_varpar: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [12 1]), ...
  'S6RRPRPR5_jacobiaD_rot_6_sym_varpar: pkin has to be [12x1] (double)');

%% Symbolic Calculation
% From jacobiaD_rot_6_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 21:40:08
% EndTime: 2019-02-26 21:40:10
% DurationCPUTime: 1.84s
% Computational Cost: add. (8842->154), mult. (24081->307), div. (726->12), fcn. (31139->15), ass. (0->131)
t320 = sin(pkin(11));
t321 = cos(pkin(11));
t324 = sin(qJ(2));
t326 = cos(qJ(2));
t307 = t324 * t320 - t326 * t321;
t322 = cos(pkin(6));
t304 = t307 * t322;
t299 = qJD(2) * t304;
t348 = t326 * t320 + t324 * t321;
t306 = t348 * qJD(2);
t327 = cos(qJ(1));
t392 = sin(pkin(6));
t357 = t327 * t392;
t305 = t348 * t322;
t393 = sin(qJ(1));
t360 = t393 * t305;
t398 = -t393 * t306 - qJD(1) * t360 + (-qJD(1) * t307 - t299) * t327 - qJD(4) * t357;
t323 = sin(qJ(4));
t325 = cos(qJ(4));
t341 = -t327 * t305 + t393 * t307;
t277 = -t323 * t341 + t325 * t357;
t358 = t326 * t392;
t359 = t324 * t392;
t335 = -t320 * t358 - t321 * t359;
t293 = -t322 * t325 - t323 * t335;
t272 = atan2(-t277, t293);
t267 = sin(t272);
t268 = cos(t272);
t243 = -t267 * t277 + t268 * t293;
t241 = 0.1e1 / t243 ^ 2;
t342 = -t327 * t307 - t360;
t353 = t393 * t392;
t338 = -t323 * t342 + t325 * t353;
t276 = t338 ^ 2;
t239 = t276 * t241 + 0.1e1;
t283 = t323 * t353 + t325 * t342;
t334 = t341 * qJD(1) + t393 * t299 - t327 * t306;
t351 = qJD(1) * t357;
t247 = t283 * qJD(4) + t323 * t334 - t325 * t351;
t385 = t247 * t241;
t275 = t277 ^ 2;
t291 = 0.1e1 / t293 ^ 2;
t271 = t275 * t291 + 0.1e1;
t269 = 0.1e1 / t271;
t347 = qJD(1) * t353;
t371 = qJD(4) * t325;
t249 = t398 * t323 - t325 * t347 - t341 * t371;
t294 = t322 * t323 - t325 * t335;
t302 = -t320 * t359 + t321 * t358;
t298 = t302 * qJD(2);
t273 = t294 * qJD(4) + t298 * t323;
t290 = 0.1e1 / t293;
t379 = t277 * t291;
t346 = -t249 * t290 + t273 * t379;
t231 = t346 * t269;
t349 = -t267 * t293 - t268 * t277;
t226 = t349 * t231 - t267 * t249 + t268 * t273;
t240 = 0.1e1 / t243;
t242 = t240 * t241;
t390 = t226 * t242;
t369 = 0.2e1 * (-t276 * t390 - t338 * t385) / t239 ^ 2;
t397 = t273 * t291;
t285 = -t327 * t304 - t348 * t393;
t343 = -t285 * t290 + t302 * t379;
t396 = t323 * t343;
t250 = (qJD(4) * t341 + t347) * t323 + t398 * t325;
t288 = t393 * t304 - t327 * t348;
t319 = pkin(12) + qJ(6);
t317 = sin(t319);
t318 = cos(t319);
t259 = t283 * t318 - t288 * t317;
t253 = 0.1e1 / t259;
t254 = 0.1e1 / t259 ^ 2;
t395 = -0.2e1 * t277;
t394 = -0.2e1 * t338;
t248 = t338 * qJD(4) + t323 * t351 + t325 * t334;
t300 = t322 * t306;
t340 = t307 * qJD(2);
t263 = t285 * qJD(1) - t393 * t300 - t327 * t340;
t234 = t259 * qJD(6) + t248 * t317 - t263 * t318;
t258 = t283 * t317 + t288 * t318;
t252 = t258 ^ 2;
t246 = t252 * t254 + 0.1e1;
t384 = t254 * t258;
t370 = qJD(6) * t258;
t235 = t248 * t318 + t263 * t317 - t370;
t387 = t235 * t253 * t254;
t389 = (t234 * t384 - t252 * t387) / t246 ^ 2;
t381 = t290 * t397;
t388 = (t249 * t379 - t275 * t381) / t271 ^ 2;
t386 = t241 * t338;
t383 = t267 * t338;
t382 = t268 * t338;
t380 = t277 * t290;
t378 = t288 * t323;
t377 = t288 * t325;
t376 = t317 * t253;
t375 = t318 * t258;
t368 = -0.2e1 * t389;
t367 = 0.2e1 * t389;
t366 = -0.2e1 * t388;
t365 = t242 * t394;
t364 = t290 * t388;
t363 = t241 * t383;
t362 = t241 * t382;
t361 = t258 * t387;
t356 = 0.2e1 * t361;
t355 = t381 * t395;
t279 = -t323 * t357 - t325 * t341;
t350 = qJD(6) * t377 - t334;
t257 = -t279 * t318 + t285 * t317;
t256 = -t279 * t317 - t285 * t318;
t345 = t254 * t375 - t376;
t344 = -t279 * t290 + t294 * t379;
t339 = -t267 + (t268 * t380 + t267) * t269;
t336 = -qJD(4) * t378 + qJD(6) * t342 - t263 * t325;
t297 = t335 * qJD(2);
t274 = -t293 * qJD(4) + t298 * t325;
t265 = t288 * qJD(1) - t327 * t300 + t393 * t340;
t261 = t317 * t342 + t318 * t377;
t260 = t317 * t377 - t318 * t342;
t244 = 0.1e1 / t246;
t237 = 0.1e1 / t239;
t236 = t269 * t396;
t233 = t344 * t269;
t230 = t339 * t338;
t228 = (-t267 * t285 + t268 * t302) * t323 + t349 * t236;
t227 = t349 * t233 - t267 * t279 + t268 * t294;
t225 = t344 * t366 + (t294 * t355 - t250 * t290 + (t249 * t294 + t273 * t279 + t274 * t277) * t291) * t269;
t223 = t366 * t396 + (t343 * t371 + (t302 * t355 - t265 * t290 + (t249 * t302 + t273 * t285 + t277 * t297) * t291) * t323) * t269;
t1 = [t364 * t394 + (-t247 * t290 - t338 * t397) * t269, t223, 0, t225, 0, 0; t277 * t240 * t369 + (-t249 * t240 + (t226 * t277 + t230 * t247) * t241) * t237 - (-t230 * t241 * t369 + (-0.2e1 * t230 * t390 + (-t231 * t269 * t380 + t366) * t363 + (t364 * t395 - t231 + (t231 - t346) * t269) * t362 - t339 * t385) * t237) * t338 (-t228 * t386 - t240 * t378) * t369 + (-t228 * t385 + (-t263 * t323 + t288 * t371) * t240 + (t228 * t365 - t241 * t378) * t226 + (t302 * t371 - t223 * t277 - t236 * t249 + t297 * t323 + (-t236 * t293 - t285 * t323) * t231) * t362 + (-t285 * t371 - t223 * t293 - t236 * t273 - t265 * t323 + (t236 * t277 - t302 * t323) * t231) * t363) * t237, 0 (-t227 * t386 - t240 * t283) * t369 + (t227 * t226 * t365 + t248 * t240 + (-t283 * t226 - t227 * t247 + (-t225 * t277 - t233 * t249 + t274 + (-t233 * t293 - t279) * t231) * t382 + (-t225 * t293 - t233 * t273 - t250 + (t233 * t277 - t294) * t231) * t383) * t241) * t237, 0, 0; (-t253 * t256 + t257 * t384) * t367 + ((t257 * qJD(6) - t250 * t317 - t265 * t318) * t253 + t257 * t356 + (-t256 * t235 - (-t256 * qJD(6) - t250 * t318 + t265 * t317) * t258 - t257 * t234) * t254) * t244 (-t253 * t260 + t261 * t384) * t367 + (t261 * t356 + t350 * t253 * t318 + t336 * t376 + (t350 * t258 * t317 - t261 * t234 - t260 * t235 - t336 * t375) * t254) * t244, 0, -t345 * t338 * t368 + (t345 * t247 - ((-qJD(6) * t253 - 0.2e1 * t361) * t318 + (t234 * t318 + (t235 - t370) * t317) * t254) * t338) * t244, 0, t368 + 0.2e1 * (t234 * t254 * t244 + (-t244 * t387 - t254 * t389) * t258) * t258;];
JaD_rot  = t1;
