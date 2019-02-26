% Zeitableitung der rotatorischen Teilmatrix der analytischen Jacobi-Matrix für Segment Nr. 4 (0=Basis) von
% S6RRRPRR15
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
%   pkin=[a2,a3,a4,a5,a6,alpha2,alpha3,d1,d2,d3,d5,d6]';
%
% Output:
% JaD_rot [3x6]
%   Zeitableitung der rotatorischen Teilmatrix der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-02-26 22:24
% Revision: d75aae1ac561373cd3be920984c3df29a1c2ecc8 (2019-02-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JaD_rot = S6RRRPRR15_jacobiaD_rot_4_sym_varpar(qJ, qJD, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(12,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRPRR15_jacobiaD_rot_4_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRRPRR15_jacobiaD_rot_4_sym_varpar: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [12 1]), ...
  'S6RRRPRR15_jacobiaD_rot_4_sym_varpar: pkin has to be [12x1] (double)');

%% Symbolic Calculation
% From jacobiaD_rot_4_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 22:24:20
% EndTime: 2019-02-26 22:24:22
% DurationCPUTime: 1.62s
% Computational Cost: add. (6571->146), mult. (20921->276), div. (667->12), fcn. (26154->13), ass. (0->138)
t263 = sin(pkin(6));
t264 = cos(pkin(7));
t267 = sin(qJ(2));
t356 = cos(qJ(3));
t317 = t356 * t267;
t266 = sin(qJ(3));
t268 = cos(qJ(2));
t335 = t266 * t268;
t288 = t264 * t317 + t335;
t289 = t264 * t335 + t317;
t265 = cos(pkin(6));
t262 = sin(pkin(7));
t339 = t262 * t266;
t323 = t265 * t339;
t213 = qJD(3) * t323 + (t288 * qJD(2) + t289 * qJD(3)) * t263;
t316 = t356 * t268;
t336 = t266 * t267;
t287 = t264 * t316 - t336;
t320 = t262 * t356;
t309 = t265 * t320;
t239 = -t287 * t263 - t309;
t237 = 0.1e1 / t239 ^ 2;
t365 = t213 * t237;
t355 = sin(qJ(1));
t261 = t355 * t267;
t307 = t265 * t261;
t269 = cos(qJ(1));
t333 = t269 * t268;
t286 = t307 - t333;
t364 = -t265 * t333 + t261;
t232 = qJD(1) * t364 + t286 * qJD(2);
t315 = t355 * t268;
t334 = t269 * t267;
t285 = t265 * t315 + t334;
t319 = t263 * t355;
t245 = t285 * t262 + t264 * t319;
t242 = 0.1e1 / t245 ^ 2;
t338 = t263 * t269;
t314 = qJD(1) * t338;
t341 = (-t232 * t262 + t264 * t314) * t242;
t308 = t262 * t319;
t224 = -t286 * t356 + (-t285 * t264 + t308) * t266;
t217 = t224 ^ 2;
t212 = t217 * t242 + 0.1e1;
t210 = 0.1e1 / t212;
t241 = 0.1e1 / t245;
t340 = t241 * t341;
t342 = t224 * t242;
t279 = t285 * t356;
t292 = t356 * t308;
t223 = t264 * t279 - t266 * t286 - t292;
t250 = t265 * t334 + t315;
t233 = t250 * qJD(1) + t285 * qJD(2);
t200 = -t233 * t356 + (t232 * t264 + t262 * t314) * t266 - t223 * qJD(3);
t351 = (t200 * t342 - t217 * t340) / t212 ^ 2;
t363 = 0.2e1 * t210 * t224 * t340 + 0.2e1 * t342 * t351;
t328 = t241 * t351;
t362 = -t210 * t341 - 0.2e1 * t328;
t299 = t320 * t338;
t281 = -t250 * t266 - t299;
t291 = t364 * t356;
t282 = t264 * t291;
t218 = t282 - t281;
t215 = t218 ^ 2;
t208 = t215 * t237 + 0.1e1;
t206 = 0.1e1 / t208;
t296 = t364 * t266;
t321 = t250 * t356;
t276 = -t264 * t296 + t321;
t234 = t285 * qJD(1) + t250 * qJD(2);
t235 = -qJD(1) * t307 - qJD(2) * t261 + (qJD(2) * t265 + qJD(1)) * t333;
t257 = t338 * t339;
t318 = t264 * t356;
t280 = -qJD(1) * t292 - qJD(3) * t257 + t234 * t318 + t235 * t266;
t201 = t276 * qJD(3) + t280;
t236 = 0.1e1 / t239;
t343 = t218 * t237;
t295 = -t201 * t236 + t213 * t343;
t187 = t295 * t206;
t209 = atan2(-t218, t239);
t204 = sin(t209);
t205 = cos(t209);
t298 = -t204 * t239 - t205 * t218;
t183 = t298 * t187 - t204 * t201 + t205 * t213;
t198 = -t204 * t218 + t205 * t239;
t195 = 0.1e1 / t198;
t196 = 0.1e1 / t198 ^ 2;
t216 = t223 ^ 2;
t194 = t196 * t216 + 0.1e1;
t192 = 0.1e1 / t194;
t349 = t192 * t196;
t199 = -qJD(1) * t299 + t224 * qJD(3) - t232 * t318 - t233 * t266;
t348 = t196 * t223;
t353 = t183 * t195 * t196;
t354 = (t199 * t348 - t216 * t353) / t194 ^ 2;
t361 = -t183 * t349 - 0.2e1 * t195 * t354;
t357 = 0.2e1 * t223;
t312 = t353 * t357;
t332 = 0.2e1 * t354;
t360 = t192 * t312 - t199 * t349 + t332 * t348;
t306 = qJD(1) * t319;
t359 = (-t250 * qJD(3) - t234 * t264 + t262 * t306) * t266 - qJD(3) * t299 + t235 * t356;
t358 = -0.2e1 * t218;
t345 = t236 * t365;
t352 = (t201 * t343 - t215 * t345) / t208 ^ 2;
t350 = t192 * t195;
t347 = t210 * t241;
t346 = t210 * t242;
t344 = t218 * t236;
t337 = t264 * t266;
t331 = -0.2e1 * t352;
t329 = t236 * t352;
t326 = t192 * t348;
t324 = t262 * t346;
t310 = t345 * t358;
t297 = t264 * t364;
t220 = -t257 + t276;
t240 = t289 * t263 + t323;
t294 = -t220 * t236 + t240 * t343;
t229 = t250 * t318 - t296;
t249 = t288 * t263;
t293 = -t229 * t236 + t249 * t343;
t290 = -t264 * t336 + t316;
t284 = -t204 + (t205 * t344 + t204) * t206;
t283 = t356 * t297;
t277 = t266 * t297 - t321;
t230 = -t285 * t266 - t286 * t318;
t231 = t286 * t337 - t279;
t226 = (t287 * qJD(2) + t290 * qJD(3)) * t263;
t214 = qJD(3) * t309 + (t290 * qJD(2) + t287 * qJD(3)) * t263;
t203 = t235 * t318 - t234 * t266 + (-t250 * t337 - t291) * qJD(3);
t202 = -qJD(3) * t282 + t359;
t191 = t293 * t206;
t189 = t294 * t206;
t184 = t298 * t189 - t204 * t220 + t205 * t240;
t182 = t293 * t331 + (t249 * t310 - t203 * t236 + (t201 * t249 + t213 * t229 + t218 * t226) * t237) * t206;
t181 = t294 * t331 + (t240 * t310 - t202 * t236 + (t201 * t240 + t213 * t220 + t214 * t218) * t237) * t206;
t1 = [t329 * t357 + (-t199 * t236 + t223 * t365) * t206, t182, t181, 0, 0, 0; (t277 * qJD(3) - t280) * t350 - (t284 * t199 + ((-t187 * t206 * t344 + t331) * t204 + (t329 * t358 - t187 + (t187 - t295) * t206) * t205) * t223) * t326 + t361 * (-t283 + t281) + t360 * t284 * t223 (t231 * qJD(3) + t232 * t266 - t233 * t318) * t350 - ((-t182 * t218 - t191 * t201 + t226 + (-t191 * t239 - t229) * t187) * t205 + (-t182 * t239 - t191 * t213 - t203 + (t191 * t218 - t249) * t187) * t204) * t326 + t361 * t230 + t360 * (t298 * t191 - t204 * t229 + t205 * t249) (t184 * t348 - t195 * t224) * t332 + (t184 * t312 + t200 * t195 + (-t224 * t183 - t184 * t199 + (-(-t181 * t218 - t189 * t201 + t214 + (-t189 * t239 - t220) * t187) * t205 - (-t181 * t239 - t189 * t213 - t202 + (t189 * t218 - t240) * t187) * t204) * t223) * t196) * t192, 0, 0, 0; (qJD(3) * t283 - t359) * t347 - (-t234 * t262 - t264 * t306) * t210 * t342 + t362 * (t257 + t277) + (-t200 * t346 + t363) * (-t262 * t364 + t264 * t338) (-t230 * qJD(3) + t232 * t356 + t233 * t337) * t347 + t233 * t224 * t324 + t362 * t231 - (-t200 * t324 + t262 * t363) * t286, t328 * t357 + (-t199 * t241 + t223 * t341) * t210, 0, 0, 0;];
JaD_rot  = t1;
