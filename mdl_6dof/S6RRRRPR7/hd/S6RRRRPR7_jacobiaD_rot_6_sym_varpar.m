% Zeitableitung der rotatorischen Teilmatrix der analytischen Jacobi-Matrix für Segment Nr. 6 (0=Basis) von
% S6RRRRPR7
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
%   pkin=[a2,a3,a4,a5,a6,alpha2,d1,d2,d3,d4,d6,theta5]';
%
% Output:
% JaD_rot [3x6]
%   Zeitableitung der rotatorischen Teilmatrix der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-02-26 22:34
% Revision: d75aae1ac561373cd3be920984c3df29a1c2ecc8 (2019-02-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JaD_rot = S6RRRRPR7_jacobiaD_rot_6_sym_varpar(qJ, qJD, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(12,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRRPR7_jacobiaD_rot_6_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRRRPR7_jacobiaD_rot_6_sym_varpar: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [12 1]), ...
  'S6RRRRPR7_jacobiaD_rot_6_sym_varpar: pkin has to be [12x1] (double)');

%% Symbolic Calculation
% From jacobiaD_rot_6_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 22:34:07
% EndTime: 2019-02-26 22:34:09
% DurationCPUTime: 1.70s
% Computational Cost: add. (17252->153), mult. (18084->301), div. (989->12), fcn. (22873->13), ass. (0->135)
t305 = qJ(3) + qJ(4) + pkin(12);
t303 = sin(t305);
t308 = cos(pkin(6));
t312 = cos(qJ(2));
t383 = sin(qJ(1));
t341 = t383 * t312;
t310 = sin(qJ(2));
t313 = cos(qJ(1));
t360 = t313 * t310;
t324 = -t308 * t360 - t341;
t304 = cos(t305);
t307 = sin(pkin(6));
t361 = t307 * t313;
t346 = t304 * t361;
t276 = -t303 * t324 + t346;
t363 = t307 * t310;
t286 = t303 * t363 - t304 * t308;
t257 = atan2(-t276, t286);
t252 = sin(t257);
t253 = cos(t257);
t243 = -t252 * t276 + t253 * t286;
t241 = 0.1e1 / t243 ^ 2;
t342 = t383 * t310;
t334 = t308 * t342;
t359 = t313 * t312;
t294 = -t334 + t359;
t343 = t307 * t383;
t281 = t294 * t303 - t304 * t343;
t271 = t281 ^ 2;
t237 = t241 * t271 + 0.1e1;
t323 = -t308 * t341 - t360;
t273 = t324 * qJD(1) + t323 * qJD(2);
t306 = qJD(3) + qJD(4);
t330 = t306 * t343 + t273;
t340 = qJD(1) * t361;
t364 = t304 * t306;
t247 = t294 * t364 + t330 * t303 - t304 * t340;
t376 = t247 * t241;
t270 = t276 ^ 2;
t284 = 0.1e1 / t286 ^ 2;
t256 = t270 * t284 + 0.1e1;
t254 = 0.1e1 / t256;
t339 = qJD(2) * t383;
t275 = -qJD(1) * t334 - t310 * t339 + (qJD(2) * t308 + qJD(1)) * t359;
t298 = t303 * t361;
t338 = t383 * qJD(1);
t333 = t307 * t338;
t249 = t275 * t303 - t306 * t298 - t304 * t333 - t324 * t364;
t357 = qJD(2) * t312;
t326 = t306 * t308 + t307 * t357;
t345 = t306 * t363;
t268 = t326 * t303 + t304 * t345;
t283 = 0.1e1 / t286;
t368 = t276 * t284;
t329 = -t249 * t283 + t268 * t368;
t231 = t329 * t254;
t331 = -t252 * t286 - t253 * t276;
t226 = t331 * t231 - t252 * t249 + t253 * t268;
t240 = 0.1e1 / t243;
t242 = t240 * t241;
t381 = t226 * t242;
t355 = 0.2e1 * (-t271 * t381 + t281 * t376) / t237 ^ 2;
t387 = t268 * t284;
t344 = t308 * t359;
t291 = -t342 + t344;
t362 = t307 * t312;
t325 = -t283 * t291 + t362 * t368;
t386 = t303 * t325;
t250 = (t306 * t324 + t333) * t303 + t275 * t304 - t306 * t346;
t282 = t294 * t304 + t303 * t343;
t311 = cos(qJ(6));
t309 = sin(qJ(6));
t366 = t323 * t309;
t265 = t282 * t311 - t366;
t259 = 0.1e1 / t265;
t260 = 0.1e1 / t265 ^ 2;
t385 = -0.2e1 * t276;
t384 = 0.2e1 * t281;
t248 = t330 * t304 + (-t294 * t306 + t340) * t303;
t272 = -qJD(1) * t344 - t313 * t357 + (t308 * t339 + t338) * t310;
t238 = t265 * qJD(6) + t248 * t309 + t272 * t311;
t365 = t323 * t311;
t264 = t282 * t309 + t365;
t258 = t264 ^ 2;
t246 = t258 * t260 + 0.1e1;
t372 = t260 * t264;
t356 = qJD(6) * t264;
t239 = t248 * t311 - t272 * t309 - t356;
t378 = t239 * t259 * t260;
t380 = (t238 * t372 - t258 * t378) / t246 ^ 2;
t370 = t283 * t387;
t379 = (t249 * t368 - t270 * t370) / t256 ^ 2;
t377 = t241 * t281;
t375 = t252 * t281;
t374 = t253 * t281;
t373 = t259 * t309;
t371 = t264 * t311;
t369 = t276 * t283;
t367 = t323 * t303;
t358 = qJD(2) * t310;
t354 = -0.2e1 * t380;
t353 = 0.2e1 * t380;
t352 = -0.2e1 * t379;
t351 = t242 * t384;
t350 = t283 * t379;
t349 = t241 * t375;
t348 = t241 * t374;
t347 = t264 * t378;
t337 = 0.2e1 * t347;
t336 = t370 * t385;
t278 = -t304 * t324 - t298;
t332 = -qJD(6) * t304 * t323 + t273;
t263 = -t278 * t311 + t291 * t309;
t262 = -t278 * t309 - t291 * t311;
t328 = t260 * t371 - t373;
t287 = t303 * t308 + t304 * t363;
t327 = -t278 * t283 + t287 * t368;
t321 = -t252 + (t253 * t369 + t252) * t254;
t320 = qJD(6) * t294 + t272 * t304 - t306 * t367;
t274 = t323 * qJD(1) + t324 * qJD(2);
t269 = -t303 * t345 + t326 * t304;
t267 = t294 * t309 + t304 * t365;
t266 = -t294 * t311 + t304 * t366;
t244 = 0.1e1 / t246;
t235 = 0.1e1 / t237;
t234 = t254 * t386;
t232 = t327 * t254;
t230 = t321 * t281;
t228 = (-t252 * t291 + t253 * t362) * t303 + t331 * t234;
t227 = t331 * t232 - t252 * t278 + t253 * t287;
t224 = t327 * t352 + (t287 * t336 - t250 * t283 + (t249 * t287 + t268 * t278 + t269 * t276) * t284) * t254;
t223 = t352 * t386 + (t325 * t364 + (t336 * t362 - t274 * t283 + (t268 * t291 + (t249 * t312 - t276 * t358) * t307) * t284) * t303) * t254;
t222 = t328 * t281 * t354 + (t328 * t247 + ((-qJD(6) * t259 - 0.2e1 * t347) * t311 + (t238 * t311 + (t239 - t356) * t309) * t260) * t281) * t244;
t221 = (t227 * t377 - t240 * t282) * t355 + (t227 * t226 * t351 + t248 * t240 + (-t282 * t226 - t227 * t247 - (-t224 * t276 - t232 * t249 + t269 + (-t232 * t286 - t278) * t231) * t374 - (-t224 * t286 - t232 * t268 - t250 + (t232 * t276 - t287) * t231) * t375) * t241) * t235;
t1 = [t350 * t384 + (-t247 * t283 + t281 * t387) * t254, t223, t224, t224, 0, 0; t276 * t240 * t355 + (-t249 * t240 + (t226 * t276 - t230 * t247) * t241) * t235 + (t230 * t241 * t355 + (0.2e1 * t230 * t381 - (-t231 * t254 * t369 + t352) * t349 - (t350 * t385 - t231 + (t231 - t329) * t254) * t348 - t321 * t376) * t235) * t281 (t228 * t377 - t240 * t367) * t355 + (-t228 * t376 + (t272 * t303 + t323 * t364) * t240 + (t228 * t351 - t241 * t367) * t226 - (-t223 * t276 - t234 * t249 + (-t303 * t358 + t312 * t364) * t307 + (-t234 * t286 - t291 * t303) * t231) * t348 - (-t291 * t364 - t223 * t286 - t234 * t268 - t274 * t303 + (t234 * t276 - t303 * t362) * t231) * t349) * t235, t221, t221, 0, 0; (-t259 * t262 + t263 * t372) * t353 + ((t263 * qJD(6) - t250 * t309 - t274 * t311) * t259 + t263 * t337 + (-t262 * t239 - (-t262 * qJD(6) - t250 * t311 + t274 * t309) * t264 - t263 * t238) * t260) * t244 (-t259 * t266 + t267 * t372) * t353 + (t267 * t337 - t332 * t259 * t311 + t320 * t373 + (-t332 * t264 * t309 - t267 * t238 - t266 * t239 - t320 * t371) * t260) * t244, t222, t222, 0, t354 + 0.2e1 * (t238 * t244 * t260 + (-t244 * t378 - t260 * t380) * t264) * t264;];
JaD_rot  = t1;
