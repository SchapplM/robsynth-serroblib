% Zeitableitung der rotatorischen Teilmatrix der analytischen Jacobi-Matrix für Segment Nr. 6 (0=Basis) von
% S6RRPRRR9
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
%   pkin=[a2,a3,a4,a5,a6,alpha2,d1,d2,d4,d5,d6,theta3]';
%
% Output:
% JaD_rot [3x6]
%   Zeitableitung der rotatorischen Teilmatrix der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-02-26 21:58
% Revision: d75aae1ac561373cd3be920984c3df29a1c2ecc8 (2019-02-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JaD_rot = S6RRPRRR9_jacobiaD_rot_6_sym_varpar(qJ, qJD, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(12,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPRRR9_jacobiaD_rot_6_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRPRRR9_jacobiaD_rot_6_sym_varpar: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [12 1]), ...
  'S6RRPRRR9_jacobiaD_rot_6_sym_varpar: pkin has to be [12x1] (double)');

%% Symbolic Calculation
% From jacobiaD_rot_6_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 21:58:47
% EndTime: 2019-02-26 21:58:48
% DurationCPUTime: 1.72s
% Computational Cost: add. (17252->153), mult. (18084->300), div. (989->12), fcn. (22873->13), ass. (0->136)
t302 = pkin(12) + qJ(4) + qJ(5);
t300 = sin(t302);
t305 = cos(pkin(6));
t309 = cos(qJ(2));
t381 = sin(qJ(1));
t338 = t381 * t309;
t307 = sin(qJ(2));
t310 = cos(qJ(1));
t358 = t310 * t307;
t321 = -t305 * t358 - t338;
t301 = cos(t302);
t304 = sin(pkin(6));
t361 = t304 * t310;
t342 = t301 * t361;
t273 = -t300 * t321 + t342;
t363 = t304 * t307;
t344 = t300 * t363;
t283 = -t305 * t301 + t344;
t254 = atan2(-t273, t283);
t249 = sin(t254);
t250 = cos(t254);
t240 = -t249 * t273 + t250 * t283;
t238 = 0.1e1 / t240 ^ 2;
t339 = t381 * t307;
t331 = t305 * t339;
t357 = t310 * t309;
t291 = -t331 + t357;
t340 = t304 * t381;
t278 = t291 * t300 - t301 * t340;
t268 = t278 ^ 2;
t234 = t268 * t238 + 0.1e1;
t320 = -t305 * t338 - t358;
t270 = qJD(1) * t321 + qJD(2) * t320;
t303 = qJD(4) + qJD(5);
t327 = t303 * t340 + t270;
t337 = qJD(1) * t361;
t364 = t301 * t303;
t244 = t291 * t364 + t300 * t327 - t301 * t337;
t374 = t244 * t238;
t267 = t273 ^ 2;
t281 = 0.1e1 / t283 ^ 2;
t253 = t267 * t281 + 0.1e1;
t251 = 0.1e1 / t253;
t335 = t381 * qJD(2);
t272 = -qJD(1) * t331 - t307 * t335 + (qJD(2) * t305 + qJD(1)) * t357;
t295 = t300 * t361;
t336 = t381 * qJD(1);
t330 = t304 * t336;
t246 = t272 * t300 - t303 * t295 - t301 * t330 - t321 * t364;
t355 = qJD(2) * t309;
t323 = t303 * t305 + t304 * t355;
t343 = t301 * t363;
t265 = t300 * t323 + t303 * t343;
t280 = 0.1e1 / t283;
t368 = t273 * t281;
t326 = -t246 * t280 + t265 * t368;
t228 = t326 * t251;
t328 = -t249 * t283 - t250 * t273;
t223 = t228 * t328 - t249 * t246 + t250 * t265;
t237 = 0.1e1 / t240;
t239 = t237 * t238;
t379 = t223 * t239;
t353 = 0.2e1 * (-t268 * t379 + t278 * t374) / t234 ^ 2;
t385 = t265 * t281;
t341 = t305 * t357;
t288 = -t339 + t341;
t362 = t304 * t309;
t322 = -t280 * t288 + t362 * t368;
t384 = t300 * t322;
t247 = t300 * (t303 * t321 + t330) + t272 * t301 - t303 * t342;
t279 = t291 * t301 + t300 * t340;
t308 = cos(qJ(6));
t306 = sin(qJ(6));
t366 = t320 * t306;
t262 = t279 * t308 - t366;
t256 = 0.1e1 / t262;
t257 = 0.1e1 / t262 ^ 2;
t383 = -0.2e1 * t273;
t382 = 0.2e1 * t278;
t245 = t327 * t301 + (-t291 * t303 + t337) * t300;
t269 = -qJD(1) * t341 - t310 * t355 + (t305 * t335 + t336) * t307;
t235 = qJD(6) * t262 + t245 * t306 + t269 * t308;
t365 = t320 * t308;
t261 = t279 * t306 + t365;
t255 = t261 ^ 2;
t243 = t255 * t257 + 0.1e1;
t371 = t257 * t261;
t354 = qJD(6) * t261;
t236 = t245 * t308 - t269 * t306 - t354;
t376 = t236 * t256 * t257;
t378 = (t235 * t371 - t255 * t376) / t243 ^ 2;
t370 = t280 * t385;
t377 = (t246 * t368 - t267 * t370) / t253 ^ 2;
t375 = t238 * t278;
t373 = t249 * t278;
t372 = t250 * t278;
t369 = t273 * t280;
t367 = t320 * t300;
t360 = t306 * t256;
t359 = t308 * t261;
t356 = qJD(2) * t307;
t352 = -0.2e1 * t378;
t351 = 0.2e1 * t378;
t350 = -0.2e1 * t377;
t349 = t239 * t382;
t348 = t280 * t377;
t347 = t238 * t373;
t346 = t238 * t372;
t345 = t261 * t376;
t334 = 0.2e1 * t345;
t333 = t370 * t383;
t275 = -t301 * t321 - t295;
t329 = -qJD(6) * t301 * t320 + t270;
t260 = -t275 * t308 + t288 * t306;
t259 = -t275 * t306 - t288 * t308;
t325 = t257 * t359 - t360;
t284 = t305 * t300 + t343;
t324 = -t275 * t280 + t284 * t368;
t318 = -t249 + (t250 * t369 + t249) * t251;
t317 = qJD(6) * t291 + t269 * t301 - t303 * t367;
t271 = qJD(1) * t320 + qJD(2) * t321;
t266 = t301 * t323 - t303 * t344;
t264 = t291 * t306 + t301 * t365;
t263 = -t291 * t308 + t301 * t366;
t241 = 0.1e1 / t243;
t232 = 0.1e1 / t234;
t231 = t251 * t384;
t229 = t324 * t251;
t227 = t318 * t278;
t225 = (-t249 * t288 + t250 * t362) * t300 + t328 * t231;
t224 = t229 * t328 - t249 * t275 + t250 * t284;
t221 = t324 * t350 + (t284 * t333 - t247 * t280 + (t246 * t284 + t265 * t275 + t266 * t273) * t281) * t251;
t220 = t350 * t384 + (t322 * t364 + (t333 * t362 - t271 * t280 + (t265 * t288 + (t246 * t309 - t273 * t356) * t304) * t281) * t300) * t251;
t219 = t325 * t278 * t352 + (t325 * t244 + ((-qJD(6) * t256 - 0.2e1 * t345) * t308 + (t235 * t308 + (t236 - t354) * t306) * t257) * t278) * t241;
t218 = (t224 * t375 - t237 * t279) * t353 + (t224 * t223 * t349 + t245 * t237 + (-t279 * t223 - t224 * t244 - (-t221 * t273 - t229 * t246 + t266 + (-t229 * t283 - t275) * t228) * t372 - (-t221 * t283 - t229 * t265 - t247 + (t229 * t273 - t284) * t228) * t373) * t238) * t232;
t1 = [t348 * t382 + (-t244 * t280 + t278 * t385) * t251, t220, 0, t221, t221, 0; t273 * t237 * t353 + (-t246 * t237 + (t223 * t273 - t227 * t244) * t238) * t232 + (t227 * t238 * t353 + (0.2e1 * t227 * t379 - (-t228 * t251 * t369 + t350) * t347 - (t348 * t383 - t228 + (t228 - t326) * t251) * t346 - t318 * t374) * t232) * t278 (t225 * t375 - t237 * t367) * t353 + (-t225 * t374 + (t269 * t300 + t320 * t364) * t237 + (t225 * t349 - t238 * t367) * t223 - (-t220 * t273 - t231 * t246 + (-t300 * t356 + t309 * t364) * t304 + (-t231 * t283 - t288 * t300) * t228) * t346 - (-t288 * t364 - t220 * t283 - t231 * t265 - t271 * t300 + (t231 * t273 - t300 * t362) * t228) * t347) * t232, 0, t218, t218, 0; (-t256 * t259 + t260 * t371) * t351 + ((qJD(6) * t260 - t247 * t306 - t271 * t308) * t256 + t260 * t334 + (-t259 * t236 - (-qJD(6) * t259 - t247 * t308 + t271 * t306) * t261 - t260 * t235) * t257) * t241 (-t256 * t263 + t264 * t371) * t351 + (t264 * t334 - t329 * t256 * t308 + t317 * t360 + (-t261 * t306 * t329 - t264 * t235 - t263 * t236 - t317 * t359) * t257) * t241, 0, t219, t219, t352 + 0.2e1 * (t235 * t257 * t241 + (-t241 * t376 - t257 * t378) * t261) * t261;];
JaD_rot  = t1;
