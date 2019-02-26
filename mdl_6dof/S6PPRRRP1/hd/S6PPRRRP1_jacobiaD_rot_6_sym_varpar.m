% Zeitableitung der rotatorischen Teilmatrix der analytischen Jacobi-Matrix für Segment Nr. 6 (0=Basis) von
% S6PPRRRP1
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
%   pkin=[a2,a3,a4,a5,a6,alpha2,alpha3,d3,d4,d5,theta1,theta2]';
%
% Output:
% JaD_rot [3x6]
%   Zeitableitung der rotatorischen Teilmatrix der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-02-26 19:41
% Revision: d75aae1ac561373cd3be920984c3df29a1c2ecc8 (2019-02-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JaD_rot = S6PPRRRP1_jacobiaD_rot_6_sym_varpar(qJ, qJD, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(12,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PPRRRP1_jacobiaD_rot_6_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6PPRRRP1_jacobiaD_rot_6_sym_varpar: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [12 1]), ...
  'S6PPRRRP1_jacobiaD_rot_6_sym_varpar: pkin has to be [12x1] (double)');

%% Symbolic Calculation
% From jacobiaD_rot_6_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 19:41:34
% EndTime: 2019-02-26 19:41:36
% DurationCPUTime: 1.69s
% Computational Cost: add. (9293->119), mult. (27507->244), div. (559->12), fcn. (36133->17), ass. (0->118)
t378 = sin(pkin(12));
t379 = sin(pkin(11));
t348 = t379 * t378;
t382 = cos(pkin(12));
t383 = cos(pkin(11));
t354 = t383 * t382;
t385 = cos(pkin(6));
t319 = -t348 * t385 + t354;
t323 = sin(qJ(3));
t326 = cos(qJ(3));
t350 = t379 * t382;
t353 = t383 * t378;
t340 = t350 * t385 + t353;
t381 = sin(pkin(6));
t349 = t379 * t381;
t380 = sin(pkin(7));
t384 = cos(pkin(7));
t387 = t340 * t384 - t349 * t380;
t303 = t319 * t323 + t326 * t387;
t318 = t353 * t385 + t350;
t339 = -t354 * t385 + t348;
t352 = t381 * t380;
t334 = -t339 * t384 - t352 * t383;
t302 = t318 * t326 + t323 * t334;
t322 = sin(qJ(4));
t325 = cos(qJ(4));
t355 = t384 * t381;
t333 = t339 * t380 - t355 * t383;
t293 = t302 * t325 + t322 * t333;
t301 = -t318 * t323 + t326 * t334;
t296 = t301 * qJD(3);
t269 = qJD(4) * t293 + t296 * t322;
t291 = t302 * t322 - t325 * t333;
t289 = t291 ^ 2;
t338 = t355 * t382 + t380 * t385;
t351 = t381 * t378;
t315 = t323 * t338 + t326 * t351;
t317 = -t352 * t382 + t384 * t385;
t308 = t315 * t322 - t317 * t325;
t306 = 0.1e1 / t308 ^ 2;
t283 = t289 * t306 + 0.1e1;
t281 = 0.1e1 / t283;
t309 = t315 * t325 + t317 * t322;
t314 = -t323 * t351 + t326 * t338;
t310 = t314 * qJD(3);
t287 = qJD(4) * t309 + t310 * t322;
t305 = 0.1e1 / t308;
t368 = t291 * t306;
t253 = (-t269 * t305 + t287 * t368) * t281;
t284 = atan2(-t291, t308);
t279 = sin(t284);
t280 = cos(t284);
t347 = -t279 * t308 - t280 * t291;
t249 = t253 * t347 - t279 * t269 + t280 * t287;
t263 = -t279 * t291 + t280 * t308;
t260 = 0.1e1 / t263;
t261 = 0.1e1 / t263 ^ 2;
t390 = t249 * t260 * t261;
t304 = t319 * t326 - t323 * t387;
t335 = t340 * t380 + t349 * t384;
t294 = t304 * t322 - t325 * t335;
t389 = 0.2e1 * t294 * t390;
t343 = -t301 * t305 + t314 * t368;
t388 = t322 * t343;
t369 = t287 * t305 * t306;
t386 = -0.2e1 * (t269 * t368 - t289 * t369) / t283 ^ 2;
t295 = t304 * t325 + t322 * t335;
t321 = sin(qJ(5));
t324 = cos(qJ(5));
t278 = t295 * t324 + t303 * t321;
t274 = 0.1e1 / t278;
t275 = 0.1e1 / t278 ^ 2;
t377 = t261 * t294;
t298 = t303 * qJD(3);
t272 = -qJD(4) * t294 - t298 * t325;
t299 = t304 * qJD(3);
t277 = t295 * t321 - t303 * t324;
t362 = qJD(5) * t277;
t265 = t272 * t324 + t299 * t321 - t362;
t376 = t265 * t274 * t275;
t264 = qJD(5) * t278 + t272 * t321 - t299 * t324;
t273 = t277 ^ 2;
t268 = t273 * t275 + 0.1e1;
t373 = t275 * t277;
t375 = 0.1e1 / t268 ^ 2 * (t264 * t373 - t273 * t376);
t374 = t274 * t321;
t372 = t277 * t324;
t371 = t279 * t294;
t370 = t280 * t294;
t367 = t303 * t322;
t366 = t303 * t325;
t363 = qJD(4) * t325;
t290 = t294 ^ 2;
t259 = t261 * t290 + 0.1e1;
t271 = qJD(4) * t295 - t298 * t322;
t361 = 0.2e1 * (t271 * t377 - t290 * t390) / t259 ^ 2;
t359 = -0.2e1 * t375;
t358 = t277 * t376;
t357 = -0.2e1 * t291 * t369;
t356 = qJD(5) * t366 - t298;
t345 = t275 * t372 - t374;
t344 = -t293 * t305 + t309 * t368;
t341 = qJD(4) * t367 + qJD(5) * t304 - t299 * t325;
t311 = t315 * qJD(3);
t297 = t302 * qJD(3);
t288 = -qJD(4) * t308 + t310 * t325;
t286 = t304 * t321 - t324 * t366;
t285 = -t304 * t324 - t321 * t366;
t270 = -qJD(4) * t291 + t296 * t325;
t266 = 0.1e1 / t268;
t257 = 0.1e1 / t259;
t255 = t281 * t388;
t254 = t344 * t281;
t251 = (-t279 * t301 + t280 * t314) * t322 + t347 * t255;
t250 = t254 * t347 - t279 * t293 + t280 * t309;
t247 = t344 * t386 + (t309 * t357 - t270 * t305 + (t269 * t309 + t287 * t293 + t288 * t291) * t306) * t281;
t246 = t386 * t388 + (t343 * t363 + (t314 * t357 + t297 * t305 + (t269 * t314 + t287 * t301 - t291 * t311) * t306) * t322) * t281;
t1 = [0, 0, t246, t247, 0, 0; 0, 0 (t251 * t377 + t260 * t367) * t361 + ((-t299 * t322 - t303 * t363) * t260 + t251 * t389 + (-t251 * t271 + t367 * t249 - (t314 * t363 - t246 * t291 - t255 * t269 - t311 * t322 + (-t255 * t308 - t301 * t322) * t253) * t370 - (-t301 * t363 - t246 * t308 - t255 * t287 + t297 * t322 + (t255 * t291 - t314 * t322) * t253) * t371) * t261) * t257 (t250 * t377 - t260 * t295) * t361 + (t250 * t389 + t272 * t260 + (-t295 * t249 - t250 * t271 - (-t247 * t291 - t254 * t269 + t288 + (-t254 * t308 - t293) * t253) * t370 - (-t247 * t308 - t254 * t287 - t270 + (t254 * t291 - t309) * t253) * t371) * t261) * t257, 0, 0; 0, 0, 0.2e1 * (-t274 * t285 + t286 * t373) * t375 + (0.2e1 * t286 * t358 - t356 * t274 * t324 + t341 * t374 + (-t277 * t321 * t356 - t286 * t264 - t285 * t265 - t341 * t372) * t275) * t266, t345 * t294 * t359 + (t345 * t271 + ((-qJD(5) * t274 - 0.2e1 * t358) * t324 + (t264 * t324 + (t265 - t362) * t321) * t275) * t294) * t266, t359 + 0.2e1 * (t264 * t275 * t266 + (-t266 * t376 - t275 * t375) * t277) * t277, 0;];
JaD_rot  = t1;
