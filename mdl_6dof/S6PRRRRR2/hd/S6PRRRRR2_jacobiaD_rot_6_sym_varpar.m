% Zeitableitung der rotatorischen Teilmatrix der analytischen Jacobi-Matrix für Segment Nr. 6 (0=Basis) von
% S6PRRRRR2
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
%   pkin=[a2,a3,a4,a5,a6,alpha2,d2,d3,d4,d5,d6,theta1]';
%
% Output:
% JaD_rot [3x6]
%   Zeitableitung der rotatorischen Teilmatrix der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-02-26 20:19
% Revision: d75aae1ac561373cd3be920984c3df29a1c2ecc8 (2019-02-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JaD_rot = S6PRRRRR2_jacobiaD_rot_6_sym_varpar(qJ, qJD, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(12,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRRRRR2_jacobiaD_rot_6_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6PRRRRR2_jacobiaD_rot_6_sym_varpar: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [12 1]), ...
  'S6PRRRRR2_jacobiaD_rot_6_sym_varpar: pkin has to be [12x1] (double)');

%% Symbolic Calculation
% From jacobiaD_rot_6_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 20:19:15
% EndTime: 2019-02-26 20:19:17
% DurationCPUTime: 1.27s
% Computational Cost: add. (9918->116), mult. (13898->233), div. (840->12), fcn. (17831->13), ass. (0->115)
t293 = sin(pkin(12));
t295 = cos(pkin(12));
t297 = sin(qJ(2));
t296 = cos(pkin(6));
t298 = cos(qJ(2));
t325 = t296 * t298;
t277 = -t293 * t297 + t295 * t325;
t273 = t277 * qJD(2);
t326 = t296 * t297;
t278 = t293 * t298 + t295 * t326;
t292 = qJ(3) + qJ(4);
t286 = sin(t292);
t290 = qJD(3) + qJD(4);
t294 = sin(pkin(6));
t329 = t294 * t295;
t316 = t286 * t329;
t288 = cos(t292);
t331 = t288 * t290;
t240 = t273 * t286 + t278 * t331 - t290 * t316;
t262 = t278 * t286 + t288 * t329;
t260 = t262 ^ 2;
t328 = t294 * t297;
t271 = t286 * t328 - t288 * t296;
t269 = 0.1e1 / t271 ^ 2;
t254 = t260 * t269 + 0.1e1;
t252 = 0.1e1 / t254;
t323 = qJD(2) * t298;
t307 = t290 * t296 + t294 * t323;
t318 = t290 * t328;
t258 = t307 * t286 + t288 * t318;
t268 = 0.1e1 / t271;
t336 = t262 * t269;
t224 = (-t240 * t268 + t258 * t336) * t252;
t255 = atan2(-t262, t271);
t250 = sin(t255);
t251 = cos(t255);
t310 = -t250 * t271 - t251 * t262;
t220 = t310 * t224 - t240 * t250 + t251 * t258;
t234 = -t250 * t262 + t251 * t271;
t231 = 0.1e1 / t234;
t232 = 0.1e1 / t234 ^ 2;
t349 = t220 * t231 * t232;
t317 = t293 * t326;
t280 = t295 * t298 - t317;
t330 = t293 * t294;
t265 = t280 * t286 - t288 * t330;
t348 = 0.2e1 * t265 * t349;
t327 = t294 * t298;
t306 = -t268 * t277 + t327 * t336;
t347 = t286 * t306;
t337 = t258 * t268 * t269;
t346 = -0.2e1 * (t240 * t336 - t260 * t337) / t254 ^ 2;
t266 = t280 * t288 + t286 * t330;
t279 = t293 * t325 + t295 * t297;
t291 = qJ(5) + qJ(6);
t285 = sin(t291);
t287 = cos(t291);
t249 = t266 * t287 + t279 * t285;
t245 = 0.1e1 / t249;
t246 = 0.1e1 / t249 ^ 2;
t276 = -qJD(2) * t317 + t295 * t323;
t289 = qJD(5) + qJD(6);
t313 = t266 * t289 - t276;
t275 = t279 * qJD(2);
t311 = t290 * t330 - t275;
t332 = t286 * t290;
t243 = -t280 * t332 + t311 * t288;
t333 = t279 * t289;
t314 = t243 + t333;
t235 = t314 * t285 + t313 * t287;
t248 = t266 * t285 - t279 * t287;
t244 = t248 ^ 2;
t239 = t244 * t246 + 0.1e1;
t341 = t246 * t248;
t236 = -t313 * t285 + t314 * t287;
t343 = t236 * t245 * t246;
t345 = (t235 * t341 - t244 * t343) / t239 ^ 2;
t344 = t232 * t265;
t342 = t245 * t285;
t340 = t248 * t287;
t339 = t250 * t265;
t338 = t251 * t265;
t335 = t279 * t286;
t334 = t279 * t288;
t324 = qJD(2) * t297;
t261 = t265 ^ 2;
t230 = t232 * t261 + 0.1e1;
t242 = t280 * t331 + t311 * t286;
t322 = 0.2e1 * (t242 * t344 - t261 * t349) / t230 ^ 2;
t321 = -0.2e1 * t345;
t319 = t248 * t343;
t315 = -0.2e1 * t262 * t337;
t312 = t288 * t333 - t275;
t309 = t246 * t340 - t342;
t264 = t278 * t288 - t316;
t272 = t286 * t296 + t288 * t328;
t308 = -t264 * t268 + t272 * t336;
t305 = -t276 * t288 + t279 * t332 + t280 * t289;
t274 = t278 * qJD(2);
t259 = -t286 * t318 + t307 * t288;
t257 = t280 * t285 - t287 * t334;
t256 = -t280 * t287 - t285 * t334;
t241 = -t278 * t332 + (-t290 * t329 + t273) * t288;
t237 = 0.1e1 / t239;
t228 = 0.1e1 / t230;
t226 = t252 * t347;
t225 = t308 * t252;
t222 = (-t250 * t277 + t251 * t327) * t286 + t310 * t226;
t221 = t310 * t225 - t250 * t264 + t251 * t272;
t219 = t308 * t346 + (t272 * t315 - t241 * t268 + (t240 * t272 + t258 * t264 + t259 * t262) * t269) * t252;
t217 = t346 * t347 + (t306 * t331 + (t315 * t327 + t268 * t274 + (t258 * t277 + (t240 * t298 - t262 * t324) * t294) * t269) * t286) * t252;
t216 = t321 + 0.2e1 * (t235 * t237 * t246 + (-t237 * t343 - t246 * t345) * t248) * t248;
t215 = t309 * t265 * t321 + (t309 * t242 + ((-t245 * t289 - 0.2e1 * t319) * t287 + (t235 * t287 + (-t248 * t289 + t236) * t285) * t246) * t265) * t237;
t214 = (t221 * t344 - t231 * t266) * t322 + (t221 * t348 + t243 * t231 + (-t266 * t220 - t221 * t242 - (-t219 * t262 - t225 * t240 + t259 + (-t225 * t271 - t264) * t224) * t338 - (-t219 * t271 - t225 * t258 - t241 + (t225 * t262 - t272) * t224) * t339) * t232) * t228;
t1 = [0, t217, t219, t219, 0, 0; 0 (t222 * t344 + t231 * t335) * t322 + ((-t276 * t286 - t279 * t331) * t231 + t222 * t348 + (-t222 * t242 + t335 * t220 - (-t217 * t262 - t226 * t240 + (-t286 * t324 + t298 * t331) * t294 + (-t226 * t271 - t277 * t286) * t224) * t338 - (-t277 * t331 - t217 * t271 - t226 * t258 + t274 * t286 + (t226 * t262 - t286 * t327) * t224) * t339) * t232) * t228, t214, t214, 0, 0; 0, 0.2e1 * (-t245 * t256 + t257 * t341) * t345 + (0.2e1 * t257 * t319 - t312 * t245 * t287 + t305 * t342 + (-t312 * t248 * t285 - t257 * t235 - t256 * t236 - t305 * t340) * t246) * t237, t215, t215, t216, t216;];
JaD_rot  = t1;
