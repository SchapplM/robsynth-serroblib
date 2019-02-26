% Zeitableitung der rotatorischen Teilmatrix der analytischen Jacobi-Matrix für Segment Nr. 6 (0=Basis) von
% S6PRRRRP2
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
%   pkin=[a2,a3,a4,a5,a6,alpha2,d2,d3,d4,d5,theta1]';
%
% Output:
% JaD_rot [3x6]
%   Zeitableitung der rotatorischen Teilmatrix der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-02-26 20:15
% Revision: d75aae1ac561373cd3be920984c3df29a1c2ecc8 (2019-02-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JaD_rot = S6PRRRRP2_jacobiaD_rot_6_sym_varpar(qJ, qJD, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(11,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRRRRP2_jacobiaD_rot_6_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6PRRRRP2_jacobiaD_rot_6_sym_varpar: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6PRRRRP2_jacobiaD_rot_6_sym_varpar: pkin has to be [11x1] (double)');

%% Symbolic Calculation
% From jacobiaD_rot_6_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 20:15:43
% EndTime: 2019-02-26 20:15:45
% DurationCPUTime: 2.06s
% Computational Cost: add. (16187->161), mult. (28326->308), div. (1062->12), fcn. (36524->13), ass. (0->132)
t311 = qJ(3) + qJ(4);
t309 = cos(t311);
t314 = sin(qJ(2));
t316 = cos(qJ(2));
t376 = cos(pkin(11));
t377 = cos(pkin(6));
t338 = t377 * t376;
t375 = sin(pkin(11));
t325 = -t314 * t338 - t375 * t316;
t336 = t375 * t314 - t316 * t338;
t300 = t336 * qJD(2);
t310 = qJD(3) + qJD(4);
t312 = sin(pkin(6));
t343 = t312 * t376;
t333 = t310 * t343 + t300;
t308 = sin(t311);
t360 = t310 * t308;
t265 = -t333 * t309 + t325 * t360;
t285 = -t308 * t343 - t309 * t325;
t315 = cos(qJ(5));
t313 = sin(qJ(5));
t328 = t336 * t313;
t275 = t285 * t315 + t328;
t323 = t325 * qJD(2);
t247 = t275 * qJD(5) + t265 * t313 + t315 * t323;
t327 = t336 * t315;
t273 = t285 * t313 - t327;
t268 = t273 ^ 2;
t359 = t312 * t314;
t345 = t309 * t359;
t296 = t377 * t308 + t345;
t357 = t315 * t316;
t291 = t296 * t313 + t312 * t357;
t289 = 0.1e1 / t291 ^ 2;
t259 = t268 * t289 + 0.1e1;
t257 = 0.1e1 / t259;
t355 = qJD(2) * t312;
t326 = t377 * t310 + t316 * t355;
t346 = t308 * t359;
t282 = t326 * t309 - t310 * t346;
t358 = t313 * t316;
t292 = t296 * t315 - t312 * t358;
t344 = t314 * t355;
t261 = t292 * qJD(5) + t282 * t313 - t315 * t344;
t288 = 0.1e1 / t291;
t366 = t273 * t289;
t234 = (-t247 * t288 + t261 * t366) * t257;
t260 = atan2(-t273, t291);
t254 = sin(t260);
t255 = cos(t260);
t335 = -t254 * t291 - t255 * t273;
t230 = t335 * t234 - t254 * t247 + t255 * t261;
t246 = -t254 * t273 + t255 * t291;
t243 = 0.1e1 / t246;
t244 = 0.1e1 / t246 ^ 2;
t381 = t230 * t243 * t244;
t337 = t377 * t375;
t305 = -t314 * t337 + t376 * t316;
t342 = t312 * t375;
t287 = t305 * t309 + t308 * t342;
t324 = -t376 * t314 - t316 * t337;
t276 = t287 * t313 + t315 * t324;
t380 = -0.2e1 * t276;
t340 = 0.2e1 * t276 * t381;
t284 = t308 * t325 - t309 * t343;
t295 = t377 * t309 - t346;
t329 = -t284 * t288 + t295 * t366;
t379 = t313 * t329;
t368 = t261 * t288 * t289;
t378 = -0.2e1 * (t247 * t366 - t268 * t368) / t259 ^ 2;
t363 = t324 * t313;
t277 = t287 * t315 - t363;
t270 = 0.1e1 / t277;
t271 = 0.1e1 / t277 ^ 2;
t286 = -t305 * t308 + t309 * t342;
t283 = t286 ^ 2;
t365 = t283 * t271;
t256 = 0.1e1 + t365;
t301 = t324 * qJD(2);
t332 = t310 * t342 + t301;
t362 = t309 * t310;
t266 = -t305 * t362 - t332 * t308;
t267 = -t305 * t360 + t332 * t309;
t302 = t305 * qJD(2);
t250 = -t276 * qJD(5) + t267 * t315 + t302 * t313;
t371 = t250 * t270 * t271;
t347 = t283 * t371;
t367 = t271 * t286;
t374 = (t266 * t367 - t347) / t256 ^ 2;
t373 = t244 * t276;
t249 = t277 * qJD(5) + t267 * t313 - t302 * t315;
t372 = t249 * t244;
t370 = t254 * t276;
t369 = t255 * t276;
t364 = t286 * t313;
t361 = t309 * t315;
t356 = qJD(2) * t309;
t354 = qJD(5) * t313;
t353 = qJD(5) * t315;
t352 = t309 * qJD(5);
t269 = t276 ^ 2;
t242 = t269 * t244 + 0.1e1;
t351 = 0.2e1 * (-t269 * t381 + t276 * t372) / t242 ^ 2;
t350 = 0.2e1 * t374;
t348 = t286 * t371;
t339 = -0.2e1 * t273 * t368;
t331 = -t275 * t288 + t292 * t366;
t278 = -t309 * t328 + t315 * t325;
t293 = (t309 * t358 - t314 * t315) * t312;
t330 = -t278 * t288 + t293 * t366;
t281 = -t326 * t308 - t310 * t345;
t280 = t305 * t313 + t324 * t361;
t279 = -t305 * t315 + t309 * t363;
t264 = t333 * t308 + t325 * t362;
t263 = ((-qJD(2) + t352) * t357 + (-t316 * t360 + (qJD(5) - t356) * t314) * t313) * t312;
t262 = -t291 * qJD(5) + t282 * t315 + t313 * t344;
t252 = 0.1e1 / t256;
t251 = t328 * t360 + (-t336 * t352 + t300) * t315 + (t313 * t356 - t354) * t325;
t248 = qJD(5) * t327 + t265 * t315 - t285 * t354 - t313 * t323;
t240 = 0.1e1 / t242;
t239 = t257 * t379;
t238 = t330 * t257;
t237 = t331 * t257;
t233 = (-t254 * t284 + t255 * t295) * t313 + t335 * t239;
t232 = t335 * t238 - t254 * t278 + t255 * t293;
t231 = t335 * t237 - t254 * t275 + t255 * t292;
t229 = (t270 * t287 + t315 * t365) * t350 + (0.2e1 * t315 * t347 - t267 * t270 + (-0.2e1 * t266 * t286 * t315 + t250 * t287 + t283 * t354) * t271) * t252;
t228 = t330 * t378 + (t293 * t339 - t251 * t288 + (t247 * t293 + t261 * t278 + t263 * t273) * t289) * t257;
t226 = t331 * t378 + (t292 * t339 - t248 * t288 + (t247 * t292 + t261 * t275 + t262 * t273) * t289) * t257;
t225 = t378 * t379 + (t329 * t353 + (t295 * t339 - t264 * t288 + (t247 * t295 + t261 * t284 + t273 * t281) * t289) * t313) * t257;
t224 = (t233 * t373 - t243 * t364) * t351 + ((t266 * t313 + t286 * t353) * t243 + (-t372 + t340) * t233 + (-t364 * t230 - (t295 * t353 - t225 * t273 - t239 * t247 + t281 * t313 + (-t239 * t291 - t284 * t313) * t234) * t369 - (-t284 * t353 - t225 * t291 - t239 * t261 - t264 * t313 + (t239 * t273 - t295 * t313) * t234) * t370) * t244) * t240;
t1 = [0, t228, t225, t225, t226, 0; 0 (t232 * t373 - t243 * t279) * t351 + (t232 * t340 + (-t279 * t230 - t232 * t249 - (-t228 * t273 - t238 * t247 + t263 + (-t238 * t291 - t278) * t234) * t369 - (-t228 * t291 - t238 * t261 - t251 + (t238 * t273 - t293) * t234) * t370) * t244 + ((t324 * t352 - t301) * t315 + (qJD(5) * t305 - t302 * t309 - t324 * t360) * t313) * t243) * t240, t224, t224 (t231 * t373 - t243 * t277) * t351 + (t231 * t340 + t250 * t243 + (-t277 * t230 - t231 * t249 - (-t226 * t273 - t237 * t247 + t262 + (-t237 * t291 - t275) * t234) * t369 - (-t226 * t291 - t237 * t261 - t248 + (t237 * t273 - t292) * t234) * t370) * t244) * t240, 0; 0 (t270 * t308 * t324 + t280 * t367) * t350 + (0.2e1 * t280 * t348 + (t302 * t308 - t324 * t362) * t270 + (-(t301 * t313 - t302 * t361 + t305 * t353) * t286 - t280 * t266 - (-t308 * t250 - (t313 * t352 + t315 * t360) * t286) * t324) * t271) * t252, t229, t229, t367 * t374 * t380 + (t348 * t380 + (t249 * t286 + t266 * t276) * t271) * t252, 0;];
JaD_rot  = t1;
