% Zeitableitung der rotatorischen Teilmatrix der analytischen Jacobi-Matrix für Segment Nr. 5 (0=Basis) von
% S6RPRPRR13
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
%   pkin=[a2,a3,a4,a5,a6,alpha2,alpha3,d1,d3,d5,d6,theta2]';
%
% Output:
% JaD_rot [3x6]
%   Zeitableitung der rotatorischen Teilmatrix der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-02-26 20:55
% Revision: d75aae1ac561373cd3be920984c3df29a1c2ecc8 (2019-02-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JaD_rot = S6RPRPRR13_jacobiaD_rot_5_sym_varpar(qJ, qJD, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(12,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRPRR13_jacobiaD_rot_5_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RPRPRR13_jacobiaD_rot_5_sym_varpar: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [12 1]), ...
  'S6RPRPRR13_jacobiaD_rot_5_sym_varpar: pkin has to be [12x1] (double)');

%% Symbolic Calculation
% From jacobiaD_rot_5_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 20:55:39
% EndTime: 2019-02-26 20:55:40
% DurationCPUTime: 1.18s
% Computational Cost: add. (4837->121), mult. (15324->242), div. (448->12), fcn. (19446->15), ass. (0->118)
t292 = cos(pkin(6));
t287 = sin(pkin(12));
t295 = sin(qJ(1));
t338 = t295 * t287;
t323 = t292 * t338;
t290 = cos(pkin(12));
t298 = cos(qJ(1));
t335 = t298 * t290;
t279 = -t323 + t335;
t294 = sin(qJ(3));
t297 = cos(qJ(3));
t336 = t298 * t287;
t337 = t295 * t290;
t308 = t292 * t337 + t336;
t288 = sin(pkin(7));
t289 = sin(pkin(6));
t342 = t289 * t295;
t326 = t288 * t342;
t291 = cos(pkin(7));
t339 = t291 * t297;
t255 = t279 * t294 - t297 * t326 + t308 * t339;
t267 = t288 * t308 + t291 * t342;
t293 = sin(qJ(5));
t296 = cos(qJ(5));
t315 = t255 * t296 - t267 * t293;
t359 = t315 * qJD(5);
t340 = t291 * t294;
t343 = t288 * t292;
t265 = (t287 * t297 + t290 * t340) * t289 + t294 * t343;
t261 = 0.1e1 / t265;
t274 = t308 * qJD(1);
t334 = qJD(1) * t298;
t275 = -qJD(1) * t323 + t290 * t334;
t277 = t292 * t336 + t337;
t341 = t289 * t298;
t325 = t288 * t341;
t318 = qJD(3) * t325;
t281 = t297 * t318;
t322 = qJD(1) * t342;
t319 = t288 * t322;
t332 = qJD(3) * t294;
t276 = -t292 * t335 + t338;
t333 = qJD(3) * t276;
t231 = t294 * t319 + t275 * t297 - t277 * t332 - t281 + (-t274 * t294 - t297 * t333) * t291;
t311 = t276 * t340 - t277 * t297;
t254 = t294 * t325 + t311;
t248 = t254 ^ 2;
t262 = 0.1e1 / t265 ^ 2;
t245 = t248 * t262 + 0.1e1;
t264 = t297 * t343 + (-t287 * t294 + t290 * t339) * t289;
t257 = t264 * qJD(3);
t346 = t257 * t262;
t345 = t261 * t346;
t354 = (-t231 * t254 * t262 - t248 * t345) / t245 ^ 2;
t358 = 0.2e1 * t261 * t354;
t232 = t311 * qJD(3) - t274 * t339 + t297 * t319 + (-t275 + t318) * t294;
t246 = atan2(t254, t265);
t241 = sin(t246);
t242 = cos(t246);
t222 = t241 * t254 + t242 * t265;
t219 = 0.1e1 / t222;
t240 = t255 * t293 + t267 * t296;
t234 = 0.1e1 / t240;
t220 = 0.1e1 / t222 ^ 2;
t235 = 0.1e1 / t240 ^ 2;
t310 = -t291 * t308 + t326;
t256 = t279 * t297 + t310 * t294;
t249 = t256 ^ 2;
t218 = t249 * t220 + 0.1e1;
t273 = t277 * qJD(1);
t272 = t276 * qJD(1);
t321 = t289 * t334;
t307 = t272 * t291 + t288 * t321;
t229 = -t279 * t332 + t307 * t294 + (t310 * qJD(3) - t273) * t297;
t353 = t220 * t256;
t243 = 0.1e1 / t245;
t314 = -t231 * t261 - t254 * t346;
t213 = t314 * t243;
t317 = -t241 * t265 + t242 * t254;
t209 = t317 * t213 - t241 * t231 + t242 * t257;
t356 = t209 * t219 * t220;
t357 = (t229 * t353 - t249 * t356) / t218 ^ 2;
t228 = t256 * qJD(3) - t273 * t294 - t307 * t297;
t259 = -t272 * t288 + t291 * t321;
t223 = t240 * qJD(5) - t228 * t296 + t259 * t293;
t233 = t315 ^ 2;
t227 = t233 * t235 + 0.1e1;
t351 = t235 * t315;
t224 = t228 * t293 + t259 * t296 + t359;
t352 = t224 * t234 * t235;
t355 = (-t223 * t351 - t233 * t352) / t227 ^ 2;
t350 = t241 * t256;
t349 = t242 * t256;
t348 = t254 * t261;
t347 = t254 * t264;
t344 = t277 * t294;
t331 = 0.2e1 * t357;
t330 = 0.2e1 * t356;
t329 = 0.2e1 * t355;
t328 = -0.2e1 * t354;
t320 = -0.2e1 * t315 * t352;
t253 = -t344 + (-t276 * t291 - t325) * t297;
t266 = -t276 * t288 + t291 * t341;
t316 = t253 * t296 - t266 * t293;
t238 = t253 * t293 + t266 * t296;
t313 = t296 * t234 - t293 * t351;
t250 = t276 * t339 + t297 * t325 + t344;
t312 = t250 * t261 - t262 * t347;
t306 = -t241 + (-t242 * t348 + t241) * t243;
t260 = -t274 * t288 - t291 * t322;
t258 = t265 * qJD(3);
t225 = 0.1e1 / t227;
t216 = 0.1e1 / t218;
t214 = t312 * t243;
t212 = t306 * t256;
t210 = t317 * t214 + t241 * t250 + t242 * t264;
t208 = t312 * t328 + (0.2e1 * t345 * t347 - t232 * t261 + (t231 * t264 - t250 * t257 + t254 * t258) * t262) * t243;
t1 = [t256 * t358 + (-t229 * t261 + t256 * t346) * t243, 0, t208, 0, 0, 0; -0.2e1 * t254 * t219 * t357 + ((t281 + (t291 * t333 - t275) * t297 + (qJD(3) * t277 + t274 * t291 - t319) * t294) * t219 + (-t254 * t209 - t212 * t229) * t220) * t216 + (t212 * t330 * t216 + (t212 * t331 + (-(t213 * t243 * t348 + t328) * t350 - (t254 * t358 - t213 + (t213 - t314) * t243) * t349 - t306 * t229) * t216) * t220) * t256, 0 (t210 * t353 + t219 * t255) * t331 + (t210 * t256 * t330 - t228 * t219 + (t255 * t209 - t210 * t229 - (t208 * t254 - t214 * t231 - t258 + (-t214 * t265 + t250) * t213) * t349 - (-t208 * t265 - t214 * t257 - t232 + (-t214 * t254 - t264) * t213) * t350) * t220) * t216, 0, 0, 0; (t234 * t316 - t238 * t351) * t329 + ((t238 * qJD(5) - t232 * t296 + t260 * t293) * t234 + t238 * t320 + (t316 * t224 + (t316 * qJD(5) + t232 * t293 + t260 * t296) * t315 - t238 * t223) * t235) * t225, 0, t313 * t256 * t329 + (-t313 * t229 + ((qJD(5) * t234 + t320) * t293 + (-t223 * t293 + (t224 + t359) * t296) * t235) * t256) * t225, 0, -0.2e1 * t355 - 0.2e1 * (t223 * t235 * t225 - (-t225 * t352 - t235 * t355) * t315) * t315, 0;];
JaD_rot  = t1;
