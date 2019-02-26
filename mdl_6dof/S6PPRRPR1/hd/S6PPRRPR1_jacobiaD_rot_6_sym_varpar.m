% Zeitableitung der rotatorischen Teilmatrix der analytischen Jacobi-Matrix für Segment Nr. 6 (0=Basis) von
% S6PPRRPR1
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
%   pkin=[a2,a3,a4,a5,a6,alpha2,alpha3,d3,d4,d6,theta1,theta2,theta5]';
%
% Output:
% JaD_rot [3x6]
%   Zeitableitung der rotatorischen Teilmatrix der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-02-26 19:40
% Revision: d75aae1ac561373cd3be920984c3df29a1c2ecc8 (2019-02-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JaD_rot = S6PPRRPR1_jacobiaD_rot_6_sym_varpar(qJ, qJD, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(13,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PPRRPR1_jacobiaD_rot_6_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6PPRRPR1_jacobiaD_rot_6_sym_varpar: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [13 1]), ...
  'S6PPRRPR1_jacobiaD_rot_6_sym_varpar: pkin has to be [13x1] (double)');

%% Symbolic Calculation
% From jacobiaD_rot_6_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 19:40:19
% EndTime: 2019-02-26 19:40:21
% DurationCPUTime: 1.67s
% Computational Cost: add. (9604->120), mult. (27507->244), div. (559->12), fcn. (36133->17), ass. (0->120)
t380 = sin(pkin(12));
t381 = sin(pkin(11));
t349 = t381 * t380;
t384 = cos(pkin(12));
t385 = cos(pkin(11));
t355 = t385 * t384;
t387 = cos(pkin(6));
t319 = -t387 * t349 + t355;
t325 = sin(qJ(3));
t327 = cos(qJ(3));
t351 = t381 * t384;
t354 = t385 * t380;
t341 = t387 * t351 + t354;
t383 = sin(pkin(6));
t350 = t381 * t383;
t382 = sin(pkin(7));
t386 = cos(pkin(7));
t389 = t341 * t386 - t382 * t350;
t303 = t319 * t325 + t389 * t327;
t318 = t387 * t354 + t351;
t340 = -t387 * t355 + t349;
t353 = t383 * t382;
t335 = -t340 * t386 - t385 * t353;
t302 = t318 * t327 + t335 * t325;
t324 = sin(qJ(4));
t326 = cos(qJ(4));
t356 = t386 * t383;
t334 = t340 * t382 - t385 * t356;
t293 = t302 * t326 + t334 * t324;
t301 = -t318 * t325 + t335 * t327;
t297 = t301 * qJD(3);
t275 = t293 * qJD(4) + t297 * t324;
t291 = t302 * t324 - t334 * t326;
t289 = t291 ^ 2;
t339 = t384 * t356 + t387 * t382;
t352 = t383 * t380;
t315 = t339 * t325 + t327 * t352;
t317 = -t384 * t353 + t387 * t386;
t308 = t315 * t324 - t317 * t326;
t306 = 0.1e1 / t308 ^ 2;
t283 = t289 * t306 + 0.1e1;
t281 = 0.1e1 / t283;
t309 = t315 * t326 + t317 * t324;
t314 = -t325 * t352 + t339 * t327;
t310 = t314 * qJD(3);
t287 = t309 * qJD(4) + t310 * t324;
t305 = 0.1e1 / t308;
t369 = t291 * t306;
t253 = (-t275 * t305 + t287 * t369) * t281;
t284 = atan2(-t291, t308);
t279 = sin(t284);
t280 = cos(t284);
t348 = -t279 * t308 - t280 * t291;
t249 = t348 * t253 - t279 * t275 + t280 * t287;
t263 = -t279 * t291 + t280 * t308;
t260 = 0.1e1 / t263;
t261 = 0.1e1 / t263 ^ 2;
t392 = t249 * t260 * t261;
t304 = t319 * t327 - t325 * t389;
t336 = t341 * t382 + t386 * t350;
t294 = t304 * t324 - t336 * t326;
t391 = 0.2e1 * t294 * t392;
t344 = -t301 * t305 + t314 * t369;
t390 = t324 * t344;
t370 = t287 * t305 * t306;
t388 = -0.2e1 * (t275 * t369 - t289 * t370) / t283 ^ 2;
t295 = t304 * t326 + t336 * t324;
t323 = pkin(13) + qJ(6);
t321 = sin(t323);
t322 = cos(t323);
t274 = t295 * t322 + t303 * t321;
t270 = 0.1e1 / t274;
t271 = 0.1e1 / t274 ^ 2;
t299 = t303 * qJD(3);
t278 = -t294 * qJD(4) - t299 * t326;
t300 = t304 * qJD(3);
t264 = t274 * qJD(6) + t278 * t321 - t300 * t322;
t273 = t295 * t321 - t303 * t322;
t269 = t273 ^ 2;
t268 = t269 * t271 + 0.1e1;
t375 = t271 * t273;
t363 = qJD(6) * t273;
t265 = t278 * t322 + t300 * t321 - t363;
t377 = t265 * t270 * t271;
t379 = (t264 * t375 - t269 * t377) / t268 ^ 2;
t378 = t261 * t294;
t376 = t270 * t321;
t374 = t273 * t322;
t277 = t295 * qJD(4) - t299 * t324;
t373 = t277 * t261;
t372 = t279 * t294;
t371 = t280 * t294;
t368 = t303 * t324;
t367 = t303 * t326;
t364 = qJD(4) * t326;
t290 = t294 ^ 2;
t259 = t290 * t261 + 0.1e1;
t362 = 0.2e1 * (-t290 * t392 + t294 * t373) / t259 ^ 2;
t361 = -0.2e1 * t379;
t359 = t273 * t377;
t358 = -0.2e1 * t291 * t370;
t357 = qJD(6) * t367 - t299;
t346 = t271 * t374 - t376;
t345 = -t293 * t305 + t309 * t369;
t342 = qJD(4) * t368 + qJD(6) * t304 - t300 * t326;
t311 = t315 * qJD(3);
t298 = t302 * qJD(3);
t288 = -t308 * qJD(4) + t310 * t326;
t286 = t304 * t321 - t322 * t367;
t285 = -t304 * t322 - t321 * t367;
t276 = -t291 * qJD(4) + t297 * t326;
t266 = 0.1e1 / t268;
t257 = 0.1e1 / t259;
t255 = t281 * t390;
t254 = t345 * t281;
t251 = (-t279 * t301 + t280 * t314) * t324 + t348 * t255;
t250 = t348 * t254 - t279 * t293 + t280 * t309;
t247 = t345 * t388 + (t309 * t358 - t276 * t305 + (t275 * t309 + t287 * t293 + t288 * t291) * t306) * t281;
t246 = t388 * t390 + (t344 * t364 + (t314 * t358 + t298 * t305 + (t275 * t314 + t287 * t301 - t291 * t311) * t306) * t324) * t281;
t1 = [0, 0, t246, t247, 0, 0; 0, 0 (t251 * t378 + t260 * t368) * t362 + ((-t300 * t324 - t303 * t364) * t260 + (-t373 + t391) * t251 + (t368 * t249 - (t314 * t364 - t246 * t291 - t255 * t275 - t311 * t324 + (-t255 * t308 - t301 * t324) * t253) * t371 - (-t301 * t364 - t246 * t308 - t255 * t287 + t298 * t324 + (t255 * t291 - t314 * t324) * t253) * t372) * t261) * t257 (t250 * t378 - t260 * t295) * t362 + (t250 * t391 + t278 * t260 + (-t295 * t249 - t250 * t277 - (-t247 * t291 - t254 * t275 + t288 + (-t254 * t308 - t293) * t253) * t371 - (-t247 * t308 - t254 * t287 - t276 + (t254 * t291 - t309) * t253) * t372) * t261) * t257, 0, 0; 0, 0, 0.2e1 * (-t270 * t285 + t286 * t375) * t379 + (0.2e1 * t286 * t359 - t357 * t270 * t322 + t342 * t376 + (-t357 * t273 * t321 - t286 * t264 - t285 * t265 - t342 * t374) * t271) * t266, t346 * t294 * t361 + (t346 * t277 + ((-qJD(6) * t270 - 0.2e1 * t359) * t322 + (t264 * t322 + (t265 - t363) * t321) * t271) * t294) * t266, 0, t361 + 0.2e1 * (t264 * t271 * t266 + (-t266 * t377 - t271 * t379) * t273) * t273;];
JaD_rot  = t1;
