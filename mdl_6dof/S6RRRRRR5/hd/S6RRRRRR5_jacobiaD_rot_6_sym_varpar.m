% Zeitableitung der rotatorischen Teilmatrix der analytischen Jacobi-Matrix für Segment Nr. 6 (0=Basis) von
% S6RRRRRR5
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
%   pkin=[a2,a3,a4,a5,a6,alpha2,d1,d2,d3,d4,d5,d6]';
%
% Output:
% JaD_rot [3x6]
%   Zeitableitung der rotatorischen Teilmatrix der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-02-26 22:49
% Revision: d75aae1ac561373cd3be920984c3df29a1c2ecc8 (2019-02-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JaD_rot = S6RRRRRR5_jacobiaD_rot_6_sym_varpar(qJ, qJD, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(12,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRRRR5_jacobiaD_rot_6_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRRRRR5_jacobiaD_rot_6_sym_varpar: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [12 1]), ...
  'S6RRRRRR5_jacobiaD_rot_6_sym_varpar: pkin has to be [12x1] (double)');

%% Symbolic Calculation
% From jacobiaD_rot_6_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 22:49:26
% EndTime: 2019-02-26 22:49:28
% DurationCPUTime: 1.88s
% Computational Cost: add. (22474->153), mult. (22690->300), div. (1252->12), fcn. (28701->13), ass. (0->136)
t306 = qJ(3) + qJ(4) + qJ(5);
t303 = sin(t306);
t308 = cos(pkin(6));
t312 = cos(qJ(2));
t384 = sin(qJ(1));
t341 = t384 * t312;
t310 = sin(qJ(2));
t313 = cos(qJ(1));
t361 = t313 * t310;
t324 = -t308 * t361 - t341;
t304 = cos(t306);
t307 = sin(pkin(6));
t363 = t307 * t313;
t345 = t304 * t363;
t276 = -t303 * t324 + t345;
t365 = t307 * t310;
t347 = t303 * t365;
t286 = -t308 * t304 + t347;
t257 = atan2(-t276, t286);
t252 = sin(t257);
t253 = cos(t257);
t243 = -t252 * t276 + t253 * t286;
t241 = 0.1e1 / t243 ^ 2;
t342 = t384 * t310;
t334 = t308 * t342;
t360 = t313 * t312;
t294 = -t334 + t360;
t343 = t307 * t384;
t281 = t294 * t303 - t304 * t343;
t275 = t281 ^ 2;
t237 = t275 * t241 + 0.1e1;
t323 = -t308 * t341 - t361;
t271 = t324 * qJD(1) + t323 * qJD(2);
t305 = qJD(3) + qJD(4) + qJD(5);
t330 = t305 * t343 + t271;
t340 = qJD(1) * t363;
t366 = t304 * t305;
t247 = t294 * t366 + t330 * t303 - t304 * t340;
t377 = t241 * t281;
t274 = t276 ^ 2;
t284 = 0.1e1 / t286 ^ 2;
t256 = t274 * t284 + 0.1e1;
t254 = 0.1e1 / t256;
t339 = qJD(2) * t384;
t273 = -qJD(1) * t334 - t310 * t339 + (qJD(2) * t308 + qJD(1)) * t360;
t298 = t303 * t363;
t338 = t384 * qJD(1);
t333 = t307 * t338;
t249 = t273 * t303 - t305 * t298 - t304 * t333 - t324 * t366;
t358 = qJD(2) * t312;
t326 = t305 * t308 + t307 * t358;
t346 = t304 * t365;
t266 = t326 * t303 + t305 * t346;
t283 = 0.1e1 / t286;
t370 = t276 * t284;
t329 = -t249 * t283 + t266 * t370;
t231 = t329 * t254;
t331 = -t252 * t286 - t253 * t276;
t226 = t331 * t231 - t252 * t249 + t253 * t266;
t240 = 0.1e1 / t243;
t242 = t240 * t241;
t382 = t226 * t242;
t356 = 0.2e1 * (t247 * t377 - t275 * t382) / t237 ^ 2;
t388 = t266 * t284;
t344 = t308 * t360;
t291 = -t342 + t344;
t364 = t307 * t312;
t325 = -t283 * t291 + t364 * t370;
t387 = t303 * t325;
t250 = (t305 * t324 + t333) * t303 + t273 * t304 - t305 * t345;
t282 = t294 * t304 + t303 * t343;
t311 = cos(qJ(6));
t309 = sin(qJ(6));
t368 = t323 * t309;
t265 = t282 * t311 - t368;
t259 = 0.1e1 / t265;
t260 = 0.1e1 / t265 ^ 2;
t386 = -0.2e1 * t276;
t385 = 0.2e1 * t281;
t248 = t330 * t304 + (-t294 * t305 + t340) * t303;
t270 = -qJD(1) * t344 - t313 * t358 + (t308 * t339 + t338) * t310;
t238 = t265 * qJD(6) + t248 * t309 + t270 * t311;
t367 = t323 * t311;
t264 = t282 * t309 + t367;
t258 = t264 ^ 2;
t246 = t258 * t260 + 0.1e1;
t374 = t260 * t264;
t357 = qJD(6) * t264;
t239 = t248 * t311 - t270 * t309 - t357;
t379 = t239 * t259 * t260;
t381 = (t238 * t374 - t258 * t379) / t246 ^ 2;
t372 = t283 * t388;
t380 = (t249 * t370 - t274 * t372) / t256 ^ 2;
t378 = t241 * t247;
t376 = t252 * t281;
t375 = t253 * t281;
t373 = t264 * t311;
t371 = t276 * t283;
t369 = t323 * t303;
t362 = t309 * t259;
t359 = qJD(2) * t310;
t355 = -0.2e1 * t381;
t354 = 0.2e1 * t381;
t353 = -0.2e1 * t380;
t352 = t242 * t385;
t351 = t283 * t380;
t350 = t241 * t376;
t349 = t241 * t375;
t348 = t264 * t379;
t337 = 0.2e1 * t348;
t336 = t372 * t386;
t278 = -t304 * t324 - t298;
t332 = -qJD(6) * t304 * t323 + t271;
t263 = -t278 * t311 + t291 * t309;
t262 = -t278 * t309 - t291 * t311;
t328 = t260 * t373 - t362;
t287 = t308 * t303 + t346;
t327 = -t278 * t283 + t287 * t370;
t321 = -t252 + (t253 * t371 + t252) * t254;
t320 = qJD(6) * t294 + t270 * t304 - t305 * t369;
t272 = t323 * qJD(1) + t324 * qJD(2);
t269 = t294 * t309 + t304 * t367;
t268 = -t294 * t311 + t304 * t368;
t267 = t326 * t304 - t305 * t347;
t244 = 0.1e1 / t246;
t235 = 0.1e1 / t237;
t234 = t254 * t387;
t233 = t327 * t254;
t230 = t321 * t281;
t228 = (-t252 * t291 + t253 * t364) * t303 + t331 * t234;
t227 = t331 * t233 - t252 * t278 + t253 * t287;
t224 = t327 * t353 + (t287 * t336 - t250 * t283 + (t249 * t287 + t266 * t278 + t267 * t276) * t284) * t254;
t223 = t353 * t387 + (t325 * t366 + (t336 * t364 - t272 * t283 + (t266 * t291 + (t249 * t312 - t276 * t359) * t307) * t284) * t303) * t254;
t222 = t328 * t281 * t355 + (t328 * t247 + ((-qJD(6) * t259 - 0.2e1 * t348) * t311 + (t238 * t311 + (t239 - t357) * t309) * t260) * t281) * t244;
t221 = (t227 * t377 - t240 * t282) * t356 + (t227 * t226 * t352 + t248 * t240 + (-t282 * t226 - t227 * t247 - (-t224 * t276 - t233 * t249 + t267 + (-t233 * t286 - t278) * t231) * t375 - (-t224 * t286 - t233 * t266 - t250 + (t233 * t276 - t287) * t231) * t376) * t241) * t235;
t1 = [t351 * t385 + (-t247 * t283 + t281 * t388) * t254, t223, t224, t224, t224, 0; t276 * t240 * t356 + (-t249 * t240 + (t226 * t276 - t230 * t247) * t241) * t235 + (t230 * t241 * t356 + (0.2e1 * t230 * t382 - (-t231 * t254 * t371 + t353) * t350 - (t351 * t386 - t231 + (t231 - t329) * t254) * t349 - t321 * t378) * t235) * t281 (t228 * t377 - t240 * t369) * t356 + (-t228 * t378 + (t270 * t303 + t323 * t366) * t240 + (t228 * t352 - t241 * t369) * t226 - (-t223 * t276 - t234 * t249 + (-t303 * t359 + t312 * t366) * t307 + (-t234 * t286 - t291 * t303) * t231) * t349 - (-t291 * t366 - t223 * t286 - t234 * t266 - t272 * t303 + (t234 * t276 - t303 * t364) * t231) * t350) * t235, t221, t221, t221, 0; (-t259 * t262 + t263 * t374) * t354 + ((t263 * qJD(6) - t250 * t309 - t272 * t311) * t259 + t263 * t337 + (-t262 * t239 - (-t262 * qJD(6) - t250 * t311 + t272 * t309) * t264 - t263 * t238) * t260) * t244 (-t259 * t268 + t269 * t374) * t354 + (t269 * t337 - t332 * t259 * t311 + t320 * t362 + (-t332 * t264 * t309 - t269 * t238 - t268 * t239 - t320 * t373) * t260) * t244, t222, t222, t222, t355 + 0.2e1 * (t238 * t260 * t244 + (-t244 * t379 - t260 * t381) * t264) * t264;];
JaD_rot  = t1;
