% Zeitableitung der rotatorischen Teilmatrix der analytischen Jacobi-Matrix für Segment Nr. 5 (0=Basis) von
% S6RRRRPP9
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
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d1,d2,d3,d4]';
%
% Output:
% JaD_rot [3x6]
%   Zeitableitung der rotatorischen Teilmatrix der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-02-26 22:30
% Revision: d75aae1ac561373cd3be920984c3df29a1c2ecc8 (2019-02-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JaD_rot = S6RRRRPP9_jacobiaD_rot_5_sym_varpar(qJ, qJD, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(10,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRRPP9_jacobiaD_rot_5_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRRRPP9_jacobiaD_rot_5_sym_varpar: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RRRRPP9_jacobiaD_rot_5_sym_varpar: pkin has to be [10x1] (double)');

%% Symbolic Calculation
% From jacobiaD_rot_5_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 22:30:06
% EndTime: 2019-02-26 22:30:09
% DurationCPUTime: 2.87s
% Computational Cost: add. (9759->187), mult. (28191->360), div. (951->12), fcn. (35841->13), ass. (0->143)
t302 = cos(pkin(6));
t306 = sin(qJ(1));
t305 = sin(qJ(2));
t351 = qJD(2) * t305;
t333 = t306 * t351;
t352 = qJD(1) * t306;
t336 = t305 * t352;
t309 = cos(qJ(2));
t310 = cos(qJ(1));
t353 = t309 * t310;
t266 = -t302 * t336 - t333 + (qJD(2) * t302 + qJD(1)) * t353;
t304 = sin(qJ(3));
t308 = cos(qJ(3));
t355 = t306 * t309;
t356 = t305 * t310;
t289 = t302 * t356 + t355;
t301 = sin(pkin(6));
t359 = t301 * t310;
t322 = t289 * t304 + t308 * t359;
t338 = t301 * t352;
t245 = t322 * qJD(3) - t266 * t308 - t304 * t338;
t339 = t304 * t359;
t280 = -t289 * t308 + t339;
t303 = sin(qJ(4));
t307 = cos(qJ(4));
t357 = t305 * t306;
t330 = -t302 * t353 + t357;
t255 = t280 * t303 + t330 * t307;
t290 = t302 * t355 + t356;
t317 = t290 * qJD(1) + t289 * qJD(2);
t230 = t255 * qJD(4) - t245 * t307 + t317 * t303;
t329 = t330 * t303;
t256 = t280 * t307 - t329;
t392 = t256 * qJD(4) + t245 * t303 + t317 * t307;
t249 = t255 ^ 2;
t360 = t301 * t308;
t288 = t302 * t304 + t305 * t360;
t354 = t307 * t309;
t275 = t288 * t303 + t301 * t354;
t268 = 0.1e1 / t275 ^ 2;
t238 = t249 * t268 + 0.1e1;
t234 = 0.1e1 / t238;
t361 = t301 * t304;
t287 = t302 * t308 - t305 * t361;
t334 = qJD(2) * t301 * t309;
t274 = t287 * qJD(3) + t308 * t334;
t358 = t303 * t309;
t276 = t288 * t307 - t301 * t358;
t335 = t301 * t351;
t246 = t276 * qJD(4) + t274 * t303 - t307 * t335;
t267 = 0.1e1 / t275;
t368 = t255 * t268;
t326 = -t246 * t368 + t267 * t392;
t214 = t326 * t234;
t239 = atan2(t255, t275);
t232 = sin(t239);
t233 = cos(t239);
t328 = -t232 * t275 + t233 * t255;
t209 = t328 * t214 + t232 * t392 + t233 * t246;
t226 = t232 * t255 + t233 * t275;
t224 = 0.1e1 / t226 ^ 2;
t391 = t209 * t224;
t223 = 0.1e1 / t226;
t388 = t223 * t391;
t291 = -t302 * t357 + t353;
t282 = t291 * t308 + t306 * t361;
t257 = t282 * t303 - t290 * t307;
t381 = 0.2e1 * t257;
t332 = t381 * t388;
t265 = -t289 * qJD(1) - t290 * qJD(2);
t321 = -t291 * t304 + t306 * t360;
t337 = qJD(1) * t359;
t242 = t321 * qJD(3) + t265 * t308 + t304 * t337;
t258 = t282 * t307 + t290 * t303;
t264 = t302 * t333 + t336 + (-qJD(1) * t302 - qJD(2)) * t353;
t227 = t258 * qJD(4) + t242 * t303 + t264 * t307;
t375 = t227 * t224;
t387 = -t375 + t332;
t241 = t282 * qJD(3) + t265 * t304 - t308 * t337;
t270 = 0.1e1 / t321;
t271 = 0.1e1 / t321 ^ 2;
t272 = t270 * t271;
t371 = t241 * t272;
t386 = 0.2e1 * t258 * t371;
t385 = t246 * t268;
t323 = t267 * t322 - t287 * t368;
t384 = t303 * t323;
t382 = 0.2e1 * t255;
t250 = t257 ^ 2;
t222 = t224 * t250 + 0.1e1;
t380 = (-t250 * t388 + t257 * t375) / t222 ^ 2;
t228 = -t257 * qJD(4) + t242 * t307 - t264 * t303;
t251 = t258 ^ 2;
t240 = t251 * t271 + 0.1e1;
t367 = t258 * t271;
t378 = (t228 * t367 + t251 * t371) / t240 ^ 2;
t370 = t267 * t385;
t377 = (-t249 * t370 + t368 * t392) / t238 ^ 2;
t376 = t224 * t257;
t374 = t232 * t257;
t373 = t233 * t257;
t372 = t241 * t271;
t369 = t255 * t267;
t366 = t258 * t282;
t365 = t271 * t322;
t364 = t321 * t303;
t363 = t290 * t304;
t362 = t290 * t308;
t350 = qJD(3) * t304;
t349 = qJD(3) * t308;
t348 = qJD(4) * t307;
t347 = qJD(4) * t308;
t346 = 0.2e1 * t380;
t345 = 0.2e1 * t378;
t344 = -0.2e1 * t377;
t342 = t270 * t378;
t341 = t267 * t377;
t331 = t370 * t382;
t327 = qJD(4) * t291 + t264 * t308;
t325 = t256 * t267 - t276 * t368;
t259 = -t289 * t307 - t308 * t329;
t283 = (-t305 * t307 + t308 * t358) * t301;
t324 = -t259 * t267 - t283 * t368;
t319 = -t232 + (-t233 * t369 + t232) * t234;
t243 = qJD(3) * t339 - t266 * t304 - t289 * t349 + t308 * t338;
t273 = -t288 * qJD(3) - t304 * t334;
t261 = t291 * t303 - t307 * t362;
t260 = -t291 * t307 - t303 * t362;
t248 = ((-qJD(2) + t347) * t354 + (-t309 * t350 + (-qJD(2) * t308 + qJD(4)) * t305) * t303) * t301;
t247 = -t275 * qJD(4) + t274 * t307 + t303 * t335;
t236 = 0.1e1 / t240;
t231 = (-t330 * t347 - t266) * t307 + (t289 * qJD(4) - t308 * t317 + t330 * t350) * t303;
t220 = 0.1e1 / t222;
t219 = t234 * t384;
t218 = t324 * t234;
t217 = t325 * t234;
t212 = (t232 * t322 + t233 * t287) * t303 + t328 * t219;
t211 = t328 * t218 - t232 * t259 + t233 * t283;
t210 = t328 * t217 + t232 * t256 + t233 * t276;
t208 = t324 * t344 + (t283 * t331 - t231 * t267 + (t246 * t259 - t248 * t255 - t283 * t392) * t268) * t234;
t206 = t325 * t344 + (t276 * t331 - t230 * t267 + (-t246 * t256 - t247 * t255 - t276 * t392) * t268) * t234;
t205 = t344 * t384 + (t323 * t348 + (t287 * t331 - t243 * t267 + (-t246 * t322 - t255 * t273 - t287 * t392) * t268) * t303) * t234;
t1 = [t341 * t381 + (-t227 * t267 + t257 * t385) * t234, t208, t205, t206, 0, 0; -0.2e1 * t255 * t223 * t380 + (t392 * t223 - t255 * t391 - (t319 * t227 + ((t214 * t234 * t369 + t344) * t232 + (t341 * t382 - t214 + (t214 - t326) * t234) * t233) * t257) * t376) * t220 + (t220 * t387 + t376 * t346) * t319 * t257 (t211 * t376 - t223 * t260) * t346 + (t211 * t332 + (-t260 * t209 - t211 * t227 - (t208 * t255 + t218 * t392 + t248 + (-t218 * t275 - t259) * t214) * t373 - (-t208 * t275 - t218 * t246 - t231 + (-t218 * t255 - t283) * t214) * t374) * t224 + ((-t290 * t347 - t265) * t307 + (t290 * t350 + t327) * t303) * t223) * t220 (t212 * t376 - t223 * t364) * t346 + ((-t241 * t303 + t321 * t348) * t223 + t387 * t212 + (-t364 * t209 - (t287 * t348 + t205 * t255 + t219 * t392 + t273 * t303 + (-t219 * t275 + t303 * t322) * t214) * t373 - (t322 * t348 - t205 * t275 - t219 * t246 - t243 * t303 + (-t219 * t255 - t287 * t303) * t214) * t374) * t224) * t220 (t210 * t376 - t223 * t258) * t346 + (t210 * t332 + t228 * t223 + (-t258 * t209 - t210 * t227 - (t206 * t255 + t217 * t392 + t247 + (-t217 * t275 + t256) * t214) * t373 - (-t206 * t275 - t217 * t246 - t230 + (-t217 * t255 - t276) * t214) * t374) * t224) * t220, 0, 0; -t258 * t345 * t365 + 0.2e1 * t256 * t342 + (t228 * t365 + t230 * t270 - t243 * t367 - t256 * t372 + t322 * t386) * t236 (t261 * t270 - t363 * t367) * t345 + (-(t265 * t303 + t327 * t307) * t270 + (-(t303 * t347 + t307 * t350) * t270 + t304 * t386) * t290 + (t228 * t363 - t261 * t241 + (-t264 * t304 + t290 * t349) * t258) * t271) * t236 (t270 * t307 * t321 + t271 * t366) * t345 + (qJD(4) * t270 * t364 + (-t228 * t282 - t242 * t258) * t271 + (-0.2e1 * t272 * t366 + (-t271 * t321 + t270) * t307) * t241) * t236, -t342 * t381 + (t227 * t270 + t257 * t372) * t236, 0, 0;];
JaD_rot  = t1;
