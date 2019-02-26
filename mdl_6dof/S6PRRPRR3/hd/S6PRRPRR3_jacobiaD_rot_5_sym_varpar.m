% Zeitableitung der rotatorischen Teilmatrix der analytischen Jacobi-Matrix für Segment Nr. 5 (0=Basis) von
% S6PRRPRR3
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
%   pkin=[a2,a3,a4,a5,a6,alpha2,alpha3,d2,d3,d5,d6,theta1,theta4]';
%
% Output:
% JaD_rot [3x6]
%   Zeitableitung der rotatorischen Teilmatrix der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-02-26 20:05
% Revision: d75aae1ac561373cd3be920984c3df29a1c2ecc8 (2019-02-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JaD_rot = S6PRRPRR3_jacobiaD_rot_5_sym_varpar(qJ, qJD, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(13,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRRPRR3_jacobiaD_rot_5_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6PRRPRR3_jacobiaD_rot_5_sym_varpar: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [13 1]), ...
  'S6PRRPRR3_jacobiaD_rot_5_sym_varpar: pkin has to be [13x1] (double)');

%% Symbolic Calculation
% From jacobiaD_rot_5_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 20:05:24
% EndTime: 2019-02-26 20:05:26
% DurationCPUTime: 1.66s
% Computational Cost: add. (9541->145), mult. (28576->275), div. (538->12), fcn. (36806->17), ass. (0->122)
t320 = sin(pkin(13));
t327 = sin(qJ(3));
t386 = cos(pkin(13));
t388 = cos(qJ(3));
t344 = t388 * t320 + t327 * t386;
t389 = qJD(3) * t344;
t324 = cos(pkin(7));
t343 = -t327 * t320 + t388 * t386;
t300 = t343 * t324;
t325 = cos(pkin(6));
t330 = cos(qJ(2));
t387 = cos(pkin(12));
t358 = t387 * t330;
t354 = t325 * t358;
t321 = sin(pkin(12));
t328 = sin(qJ(2));
t375 = t321 * t328;
t306 = t354 - t375;
t323 = sin(pkin(6));
t322 = sin(pkin(7));
t341 = t322 * t343;
t339 = t323 * t341;
t359 = t387 * t328;
t374 = t321 * t330;
t345 = -t325 * t359 - t374;
t277 = t300 * t306 - t387 * t339 + t344 * t345;
t370 = t323 * t330;
t371 = t323 * t328;
t286 = -t300 * t370 - t325 * t341 + t344 * t371;
t260 = atan2(t277, t286);
t255 = sin(t260);
t256 = cos(t260);
t245 = t255 * t277 + t256 * t286;
t242 = 0.1e1 / t245;
t346 = t325 * t374 + t359;
t376 = t321 * t323;
t292 = t322 * t346 + t324 * t376;
t326 = sin(qJ(5));
t329 = cos(qJ(5));
t299 = t344 * t322;
t301 = t344 * t324;
t347 = t325 * t375 - t358;
t342 = t299 * t376 - t301 * t346 - t343 * t347;
t268 = t292 * t326 + t329 * t342;
t264 = 0.1e1 / t268;
t282 = 0.1e1 / t286;
t243 = 0.1e1 / t245 ^ 2;
t265 = 0.1e1 / t268 ^ 2;
t283 = 0.1e1 / t286 ^ 2;
t298 = t324 * t389;
t368 = qJD(2) * t328;
t302 = -qJD(2) * t354 + t321 * t368;
t303 = t345 * qJD(2);
t310 = t343 * qJD(3);
t338 = t322 * t389;
t337 = t323 * t338;
t252 = -t298 * t306 + t300 * t303 + t302 * t344 + t310 * t345 + t387 * t337;
t274 = t277 ^ 2;
t259 = t274 * t283 + 0.1e1;
t257 = 0.1e1 / t259;
t270 = t323 * t300 * t368 + t310 * t371 + t325 * t338 + (qJD(2) * t344 + t298) * t370;
t377 = t277 * t283;
t235 = (t252 * t282 - t270 * t377) * t257;
t352 = -t255 * t286 + t256 * t277;
t231 = t352 * t235 + t252 * t255 + t256 * t270;
t385 = t231 * t242 * t243;
t296 = t322 * t310;
t297 = qJD(3) * t300;
t304 = t346 * qJD(2);
t305 = t347 * qJD(2);
t254 = t296 * t376 - t297 * t346 + t301 * t305 - t304 * t343 + t347 * t389;
t372 = t322 * t329;
t246 = t268 * qJD(5) + t254 * t326 + t305 * t372;
t267 = -t292 * t329 + t326 * t342;
t263 = t267 ^ 2;
t250 = t263 * t265 + 0.1e1;
t379 = t265 * t267;
t366 = qJD(5) * t267;
t373 = t322 * t326;
t247 = t254 * t329 - t305 * t373 - t366;
t382 = t247 * t264 * t265;
t384 = (t246 * t379 - t263 * t382) / t250 ^ 2;
t279 = -t300 * t346 + t321 * t339 + t344 * t347;
t383 = t243 * t279;
t381 = t255 * t279;
t380 = t256 * t279;
t378 = t270 * t282 * t283;
t275 = t279 ^ 2;
t241 = t243 * t275 + 0.1e1;
t253 = t298 * t346 + t305 * t300 + t304 * t344 + t310 * t347 - t321 * t337;
t365 = 0.2e1 * (t253 * t383 - t275 * t385) / t241 ^ 2;
t364 = 0.2e1 * t384;
t363 = 0.2e1 * (t252 * t377 - t274 * t378) / t259 ^ 2;
t361 = t323 * t387;
t357 = 0.2e1 * t267 * t382;
t356 = 0.2e1 * t277 * t378;
t355 = -0.2e1 * t279 * t385;
t351 = -t264 * t326 + t329 * t379;
t276 = t299 * t361 - t301 * t306 + t343 * t345;
t285 = t325 * t299 + (t301 * t330 + t328 * t343) * t323;
t350 = -t276 * t282 + t285 * t377;
t287 = t300 * t345 - t306 * t344;
t290 = (t300 * t328 + t330 * t344) * t323;
t349 = -t282 * t287 + t290 * t377;
t289 = t301 * t347 - t343 * t346;
t348 = -t289 * t326 - t347 * t372;
t273 = t289 * t329 - t347 * t373;
t288 = -t300 * t347 - t344 * t346;
t271 = (-t298 * t328 + t310 * t330 + (t300 * t330 - t328 * t344) * qJD(2)) * t323;
t269 = t325 * t296 + (t297 * t330 - t389 * t328 + (-t301 * t328 + t330 * t343) * qJD(2)) * t323;
t262 = t297 * t347 + t301 * t304 + t305 * t343 + t346 * t389;
t261 = -t298 * t345 + t300 * t302 - t303 * t344 - t306 * t310;
t251 = t296 * t361 - t306 * t297 - t303 * t301 + t302 * t343 - t345 * t389;
t248 = 0.1e1 / t250;
t239 = 0.1e1 / t241;
t237 = t349 * t257;
t236 = t350 * t257;
t233 = -t352 * t237 + t255 * t287 + t256 * t290;
t232 = -t352 * t236 + t255 * t276 + t256 * t285;
t230 = t349 * t363 + (t290 * t356 + t261 * t282 + (-t252 * t290 - t270 * t287 - t271 * t277) * t283) * t257;
t228 = t350 * t363 + (t285 * t356 + t251 * t282 + (-t252 * t285 - t269 * t277 - t270 * t276) * t283) * t257;
t1 = [0, t230, t228, 0, 0, 0; 0 (-t233 * t383 - t242 * t288) * t365 + ((t298 * t347 - t300 * t304 + t305 * t344 - t310 * t346) * t242 + t233 * t355 + (-t288 * t231 + t233 * t253 + (t230 * t277 - t237 * t252 + t271 + (t237 * t286 + t287) * t235) * t380 + (-t230 * t286 + t237 * t270 + t261 + (t237 * t277 - t290) * t235) * t381) * t243) * t239 (-t232 * t383 - t242 * t342) * t365 + (t254 * t242 + t232 * t355 + (-t342 * t231 + t232 * t253 + (t228 * t277 - t236 * t252 + t269 + (t236 * t286 + t276) * t235) * t380 + (-t228 * t286 + t236 * t270 + t251 + (t236 * t277 - t285) * t235) * t381) * t243) * t239, 0, 0, 0; 0 (t264 * t348 + t273 * t379) * t364 + ((t273 * qJD(5) + t262 * t326 + t304 * t372) * t264 + t273 * t357 + (t348 * t247 - (t348 * qJD(5) + t262 * t329 - t304 * t373) * t267 - t273 * t246) * t265) * t248, t351 * t279 * t364 + (-t351 * t253 + ((qJD(5) * t264 + t357) * t329 + (-t246 * t329 + (-t247 + t366) * t326) * t265) * t279) * t248, 0, -0.2e1 * t384 + 0.2e1 * (t246 * t248 * t265 + (-t248 * t382 - t265 * t384) * t267) * t267, 0;];
JaD_rot  = t1;
