% Zeitableitung der rotatorischen Teilmatrix der analytischen Jacobi-Matrix für Segment Nr. 6 (0=Basis) von
% S6RRRPRR7
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
%   pkin=[a2,a3,a4,a5,a6,alpha2,d1,d2,d3,d5,d6,theta4]';
%
% Output:
% JaD_rot [3x6]
%   Zeitableitung der rotatorischen Teilmatrix der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-02-26 22:19
% Revision: d75aae1ac561373cd3be920984c3df29a1c2ecc8 (2019-02-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JaD_rot = S6RRRPRR7_jacobiaD_rot_6_sym_varpar(qJ, qJD, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(12,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRPRR7_jacobiaD_rot_6_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRRPRR7_jacobiaD_rot_6_sym_varpar: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [12 1]), ...
  'S6RRRPRR7_jacobiaD_rot_6_sym_varpar: pkin has to be [12x1] (double)');

%% Symbolic Calculation
% From jacobiaD_rot_6_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 22:19:19
% EndTime: 2019-02-26 22:19:21
% DurationCPUTime: 1.66s
% Computational Cost: add. (17252->153), mult. (18084->300), div. (989->12), fcn. (22873->13), ass. (0->136)
t303 = qJ(3) + pkin(12) + qJ(5);
t301 = sin(t303);
t306 = cos(pkin(6));
t310 = cos(qJ(2));
t382 = sin(qJ(1));
t339 = t382 * t310;
t308 = sin(qJ(2));
t311 = cos(qJ(1));
t359 = t311 * t308;
t322 = -t306 * t359 - t339;
t302 = cos(t303);
t305 = sin(pkin(6));
t361 = t305 * t311;
t343 = t302 * t361;
t274 = -t301 * t322 + t343;
t363 = t305 * t308;
t345 = t301 * t363;
t284 = -t306 * t302 + t345;
t255 = atan2(-t274, t284);
t250 = sin(t255);
t251 = cos(t255);
t241 = -t250 * t274 + t251 * t284;
t239 = 0.1e1 / t241 ^ 2;
t340 = t382 * t308;
t332 = t306 * t340;
t358 = t311 * t310;
t292 = -t332 + t358;
t341 = t305 * t382;
t279 = t292 * t301 - t302 * t341;
t269 = t279 ^ 2;
t235 = t269 * t239 + 0.1e1;
t321 = -t306 * t339 - t359;
t271 = t322 * qJD(1) + t321 * qJD(2);
t304 = qJD(3) + qJD(5);
t328 = t304 * t341 + t271;
t338 = qJD(1) * t361;
t364 = t302 * t304;
t245 = t292 * t364 + t328 * t301 - t302 * t338;
t375 = t239 * t279;
t268 = t274 ^ 2;
t282 = 0.1e1 / t284 ^ 2;
t254 = t268 * t282 + 0.1e1;
t252 = 0.1e1 / t254;
t337 = qJD(2) * t382;
t273 = -qJD(1) * t332 - t308 * t337 + (qJD(2) * t306 + qJD(1)) * t358;
t296 = t301 * t361;
t336 = t382 * qJD(1);
t331 = t305 * t336;
t247 = t273 * t301 - t304 * t296 - t302 * t331 - t322 * t364;
t356 = qJD(2) * t310;
t324 = t304 * t306 + t305 * t356;
t344 = t302 * t363;
t266 = t324 * t301 + t304 * t344;
t281 = 0.1e1 / t284;
t368 = t274 * t282;
t327 = -t247 * t281 + t266 * t368;
t229 = t327 * t252;
t329 = -t250 * t284 - t251 * t274;
t224 = t329 * t229 - t250 * t247 + t251 * t266;
t238 = 0.1e1 / t241;
t240 = t238 * t239;
t380 = t224 * t240;
t354 = 0.2e1 * (t245 * t375 - t269 * t380) / t235 ^ 2;
t386 = t266 * t282;
t342 = t306 * t358;
t289 = -t340 + t342;
t362 = t305 * t310;
t323 = -t281 * t289 + t362 * t368;
t385 = t301 * t323;
t248 = (t304 * t322 + t331) * t301 + t273 * t302 - t304 * t343;
t280 = t292 * t302 + t301 * t341;
t309 = cos(qJ(6));
t307 = sin(qJ(6));
t366 = t321 * t307;
t263 = t280 * t309 - t366;
t257 = 0.1e1 / t263;
t258 = 0.1e1 / t263 ^ 2;
t384 = -0.2e1 * t274;
t383 = 0.2e1 * t279;
t246 = t328 * t302 + (-t292 * t304 + t338) * t301;
t270 = -qJD(1) * t342 - t311 * t356 + (t306 * t337 + t336) * t308;
t236 = t263 * qJD(6) + t246 * t307 + t270 * t309;
t365 = t321 * t309;
t262 = t280 * t307 + t365;
t256 = t262 ^ 2;
t244 = t256 * t258 + 0.1e1;
t372 = t258 * t262;
t355 = qJD(6) * t262;
t237 = t246 * t309 - t270 * t307 - t355;
t377 = t237 * t257 * t258;
t379 = (t236 * t372 - t256 * t377) / t244 ^ 2;
t370 = t281 * t386;
t378 = (t247 * t368 - t268 * t370) / t254 ^ 2;
t376 = t239 * t245;
t374 = t250 * t279;
t373 = t251 * t279;
t371 = t262 * t309;
t369 = t274 * t281;
t367 = t321 * t301;
t360 = t307 * t257;
t357 = qJD(2) * t308;
t353 = -0.2e1 * t379;
t352 = 0.2e1 * t379;
t351 = -0.2e1 * t378;
t350 = t240 * t383;
t349 = t281 * t378;
t348 = t239 * t374;
t347 = t239 * t373;
t346 = t262 * t377;
t335 = 0.2e1 * t346;
t334 = t370 * t384;
t276 = -t302 * t322 - t296;
t330 = -qJD(6) * t302 * t321 + t271;
t261 = -t276 * t309 + t289 * t307;
t260 = -t276 * t307 - t289 * t309;
t326 = t258 * t371 - t360;
t285 = t306 * t301 + t344;
t325 = -t276 * t281 + t285 * t368;
t319 = -t250 + (t251 * t369 + t250) * t252;
t318 = qJD(6) * t292 + t270 * t302 - t304 * t367;
t272 = t321 * qJD(1) + t322 * qJD(2);
t267 = t324 * t302 - t304 * t345;
t265 = t292 * t307 + t302 * t365;
t264 = -t292 * t309 + t302 * t366;
t242 = 0.1e1 / t244;
t233 = 0.1e1 / t235;
t232 = t252 * t385;
t230 = t325 * t252;
t228 = t319 * t279;
t226 = (-t250 * t289 + t251 * t362) * t301 + t329 * t232;
t225 = t329 * t230 - t250 * t276 + t251 * t285;
t222 = t325 * t351 + (t285 * t334 - t248 * t281 + (t247 * t285 + t266 * t276 + t267 * t274) * t282) * t252;
t221 = t351 * t385 + (t323 * t364 + (t334 * t362 - t272 * t281 + (t266 * t289 + (t247 * t310 - t274 * t357) * t305) * t282) * t301) * t252;
t220 = t326 * t279 * t353 + (t326 * t245 + ((-qJD(6) * t257 - 0.2e1 * t346) * t309 + (t236 * t309 + (t237 - t355) * t307) * t258) * t279) * t242;
t219 = (t225 * t375 - t238 * t280) * t354 + (t225 * t224 * t350 + t246 * t238 + (-t280 * t224 - t225 * t245 - (-t222 * t274 - t230 * t247 + t267 + (-t230 * t284 - t276) * t229) * t373 - (-t222 * t284 - t230 * t266 - t248 + (t230 * t274 - t285) * t229) * t374) * t239) * t233;
t1 = [t349 * t383 + (-t245 * t281 + t279 * t386) * t252, t221, t222, 0, t222, 0; t274 * t238 * t354 + (-t247 * t238 + (t224 * t274 - t228 * t245) * t239) * t233 + (t228 * t239 * t354 + (0.2e1 * t228 * t380 - (-t229 * t252 * t369 + t351) * t348 - (t349 * t384 - t229 + (t229 - t327) * t252) * t347 - t319 * t376) * t233) * t279 (t226 * t375 - t238 * t367) * t354 + (-t226 * t376 + (t270 * t301 + t321 * t364) * t238 + (t226 * t350 - t239 * t367) * t224 - (-t221 * t274 - t232 * t247 + (-t301 * t357 + t310 * t364) * t305 + (-t232 * t284 - t289 * t301) * t229) * t347 - (-t289 * t364 - t221 * t284 - t232 * t266 - t272 * t301 + (t232 * t274 - t301 * t362) * t229) * t348) * t233, t219, 0, t219, 0; (-t257 * t260 + t261 * t372) * t352 + ((t261 * qJD(6) - t248 * t307 - t272 * t309) * t257 + t261 * t335 + (-t260 * t237 - (-t260 * qJD(6) - t248 * t309 + t272 * t307) * t262 - t261 * t236) * t258) * t242 (-t257 * t264 + t265 * t372) * t352 + (t265 * t335 - t330 * t257 * t309 + t318 * t360 + (-t330 * t262 * t307 - t265 * t236 - t264 * t237 - t318 * t371) * t258) * t242, t220, 0, t220, t353 + 0.2e1 * (t236 * t258 * t242 + (-t242 * t377 - t258 * t379) * t262) * t262;];
JaD_rot  = t1;
