% Zeitableitung der rotatorischen Teilmatrix der analytischen Jacobi-Matrix für Segment Nr. 6 (0=Basis) von
% S6RRRRRP7
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
%   pkin=[a2,a3,a4,a5,a6,alpha2,d1,d2,d3,d4,d5]';
%
% Output:
% JaD_rot [3x6]
%   Zeitableitung der rotatorischen Teilmatrix der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-02-26 22:43
% Revision: d75aae1ac561373cd3be920984c3df29a1c2ecc8 (2019-02-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JaD_rot = S6RRRRRP7_jacobiaD_rot_6_sym_varpar(qJ, qJD, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(11,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRRRP7_jacobiaD_rot_6_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRRRRP7_jacobiaD_rot_6_sym_varpar: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6RRRRRP7_jacobiaD_rot_6_sym_varpar: pkin has to be [11x1] (double)');

%% Symbolic Calculation
% From jacobiaD_rot_6_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 22:43:13
% EndTime: 2019-02-26 22:43:15
% DurationCPUTime: 1.63s
% Computational Cost: add. (11947->153), mult. (18084->301), div. (989->12), fcn. (22873->13), ass. (0->135)
t306 = qJ(3) + qJ(4);
t303 = sin(t306);
t308 = cos(pkin(6));
t312 = cos(qJ(2));
t383 = sin(qJ(1));
t341 = t383 * t312;
t310 = sin(qJ(2));
t313 = cos(qJ(1));
t360 = t313 * t310;
t324 = -t308 * t360 - t341;
t304 = cos(t306);
t307 = sin(pkin(6));
t361 = t307 * t313;
t346 = t304 * t361;
t276 = -t303 * t324 + t346;
t363 = t307 * t310;
t286 = t303 * t363 - t304 * t308;
t265 = atan2(-t276, t286);
t256 = sin(t265);
t257 = cos(t265);
t243 = -t256 * t276 + t257 * t286;
t241 = 0.1e1 / t243 ^ 2;
t342 = t383 * t310;
t334 = t308 * t342;
t359 = t313 * t312;
t292 = -t334 + t359;
t343 = t307 * t383;
t281 = t292 * t303 - t304 * t343;
t275 = t281 ^ 2;
t237 = t241 * t275 + 0.1e1;
t323 = -t308 * t341 - t360;
t271 = qJD(1) * t324 + qJD(2) * t323;
t305 = qJD(3) + qJD(4);
t330 = t305 * t343 + t271;
t340 = qJD(1) * t361;
t364 = t304 * t305;
t247 = t292 * t364 + t330 * t303 - t304 * t340;
t377 = t241 * t281;
t274 = t276 ^ 2;
t284 = 0.1e1 / t286 ^ 2;
t264 = t274 * t284 + 0.1e1;
t262 = 0.1e1 / t264;
t339 = qJD(2) * t383;
t273 = -qJD(1) * t334 - t310 * t339 + (qJD(2) * t308 + qJD(1)) * t359;
t298 = t303 * t361;
t338 = t383 * qJD(1);
t333 = t307 * t338;
t249 = t273 * t303 - t305 * t298 - t304 * t333 - t324 * t364;
t357 = qJD(2) * t312;
t326 = t305 * t308 + t307 * t357;
t345 = t305 * t363;
t268 = t303 * t326 + t304 * t345;
t283 = 0.1e1 / t286;
t368 = t276 * t284;
t329 = -t249 * t283 + t268 * t368;
t231 = t329 * t262;
t331 = -t256 * t286 - t257 * t276;
t226 = t331 * t231 - t256 * t249 + t257 * t268;
t240 = 0.1e1 / t243;
t242 = t240 * t241;
t381 = t226 * t242;
t355 = 0.2e1 * (t247 * t377 - t275 * t381) / t237 ^ 2;
t387 = t268 * t284;
t344 = t308 * t359;
t289 = -t342 + t344;
t362 = t307 * t312;
t325 = -t283 * t289 + t362 * t368;
t386 = t303 * t325;
t250 = t303 * (t305 * t324 + t333) + t273 * t304 - t305 * t346;
t282 = t292 * t304 + t303 * t343;
t311 = cos(qJ(5));
t309 = sin(qJ(5));
t366 = t323 * t309;
t261 = t282 * t311 - t366;
t253 = 0.1e1 / t261;
t254 = 0.1e1 / t261 ^ 2;
t385 = -0.2e1 * t276;
t384 = 0.2e1 * t281;
t370 = t283 * t387;
t380 = (t249 * t368 - t274 * t370) / t264 ^ 2;
t248 = t330 * t304 + (-t292 * t305 + t340) * t303;
t270 = -qJD(1) * t344 - t313 * t357 + (t308 * t339 + t338) * t310;
t365 = t323 * t311;
t260 = t282 * t309 + t365;
t356 = qJD(5) * t260;
t239 = t248 * t311 - t270 * t309 - t356;
t379 = t239 * t253 * t254;
t378 = t241 * t247;
t238 = t261 * qJD(5) + t248 * t309 + t270 * t311;
t252 = t260 ^ 2;
t246 = t252 * t254 + 0.1e1;
t374 = t254 * t260;
t376 = 0.1e1 / t246 ^ 2 * (t238 * t374 - t252 * t379);
t375 = t253 * t309;
t373 = t256 * t281;
t372 = t257 * t281;
t371 = t260 * t311;
t369 = t276 * t283;
t367 = t323 * t303;
t358 = qJD(2) * t310;
t354 = -0.2e1 * t380;
t353 = t242 * t384;
t352 = -0.2e1 * t376;
t351 = 0.2e1 * t376;
t350 = t283 * t380;
t349 = t260 * t379;
t348 = t241 * t373;
t347 = t241 * t372;
t337 = 0.2e1 * t349;
t336 = t370 * t385;
t278 = -t304 * t324 - t298;
t332 = -qJD(5) * t304 * t323 + t271;
t259 = -t278 * t311 + t289 * t309;
t258 = -t278 * t309 - t289 * t311;
t328 = t254 * t371 - t375;
t287 = t303 * t308 + t304 * t363;
t327 = -t278 * t283 + t287 * t368;
t321 = -t256 + (t257 * t369 + t256) * t262;
t320 = qJD(5) * t292 + t270 * t304 - t305 * t367;
t272 = qJD(1) * t323 + qJD(2) * t324;
t269 = -t303 * t345 + t304 * t326;
t267 = t292 * t309 + t304 * t365;
t266 = -t292 * t311 + t304 * t366;
t244 = 0.1e1 / t246;
t235 = 0.1e1 / t237;
t234 = t262 * t386;
t233 = t327 * t262;
t230 = t321 * t281;
t228 = (-t256 * t289 + t257 * t362) * t303 + t331 * t234;
t227 = t331 * t233 - t256 * t278 + t257 * t287;
t224 = t327 * t354 + (t287 * t336 - t250 * t283 + (t249 * t287 + t268 * t278 + t269 * t276) * t284) * t262;
t223 = t354 * t386 + (t325 * t364 + (t336 * t362 - t272 * t283 + (t268 * t289 + (t249 * t312 - t276 * t358) * t307) * t284) * t303) * t262;
t222 = t328 * t281 * t352 + (t328 * t247 + ((-qJD(5) * t253 - 0.2e1 * t349) * t311 + (t238 * t311 + (t239 - t356) * t309) * t254) * t281) * t244;
t221 = (t227 * t377 - t240 * t282) * t355 + (t227 * t226 * t353 + t248 * t240 + (-t282 * t226 - t227 * t247 - (-t224 * t276 - t233 * t249 + t269 + (-t233 * t286 - t278) * t231) * t372 - (-t224 * t286 - t233 * t268 - t250 + (t233 * t276 - t287) * t231) * t373) * t241) * t235;
t1 = [t350 * t384 + (-t247 * t283 + t281 * t387) * t262, t223, t224, t224, 0, 0; t276 * t240 * t355 + (-t249 * t240 + (t226 * t276 - t230 * t247) * t241) * t235 + (t230 * t241 * t355 + (0.2e1 * t230 * t381 - (-t231 * t262 * t369 + t354) * t348 - (t350 * t385 - t231 + (t231 - t329) * t262) * t347 - t321 * t378) * t235) * t281 (t228 * t377 - t240 * t367) * t355 + (-t228 * t378 + (t270 * t303 + t323 * t364) * t240 + (t228 * t353 - t241 * t367) * t226 - (-t223 * t276 - t234 * t249 + (-t303 * t358 + t312 * t364) * t307 + (-t234 * t286 - t289 * t303) * t231) * t347 - (-t289 * t364 - t223 * t286 - t234 * t268 - t272 * t303 + (t234 * t276 - t303 * t362) * t231) * t348) * t235, t221, t221, 0, 0; (-t253 * t258 + t259 * t374) * t351 + ((t259 * qJD(5) - t250 * t309 - t272 * t311) * t253 + t259 * t337 + (-t258 * t239 - (-t258 * qJD(5) - t250 * t311 + t272 * t309) * t260 - t259 * t238) * t254) * t244 (-t253 * t266 + t267 * t374) * t351 + (t267 * t337 - t332 * t253 * t311 + t320 * t375 + (-t332 * t260 * t309 - t267 * t238 - t266 * t239 - t320 * t371) * t254) * t244, t222, t222, t352 + 0.2e1 * (t238 * t254 * t244 + (-t244 * t379 - t254 * t376) * t260) * t260, 0;];
JaD_rot  = t1;
