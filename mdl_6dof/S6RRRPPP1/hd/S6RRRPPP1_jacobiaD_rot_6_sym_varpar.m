% Zeitableitung der rotatorischen Teilmatrix der analytischen Jacobi-Matrix für Segment Nr. 6 (0=Basis) von
% S6RRRPPP1
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
%   pkin=[a2,a3,a4,a5,a6,alpha4,d1,d2,d3,theta4]';
%
% Output:
% JaD_rot [3x6]
%   Zeitableitung der rotatorischen Teilmatrix der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-02-26 22:02
% Revision: d75aae1ac561373cd3be920984c3df29a1c2ecc8 (2019-02-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JaD_rot = S6RRRPPP1_jacobiaD_rot_6_sym_varpar(qJ, qJD, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(10,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRPPP1_jacobiaD_rot_6_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRRPPP1_jacobiaD_rot_6_sym_varpar: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RRRPPP1_jacobiaD_rot_6_sym_varpar: pkin has to be [10x1] (double)');

%% Symbolic Calculation
% From jacobiaD_rot_6_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 22:02:48
% EndTime: 2019-02-26 22:02:50
% DurationCPUTime: 2.10s
% Computational Cost: add. (6418->167), mult. (20749->337), div. (681->12), fcn. (25776->13), ass. (0->146)
t284 = sin(qJ(3));
t286 = sin(qJ(1));
t288 = cos(qJ(2));
t287 = cos(qJ(3));
t289 = cos(qJ(1));
t343 = t287 * t289;
t267 = t286 * t284 + t288 * t343;
t335 = qJD(3) * t289;
t317 = t287 * t335;
t336 = qJD(3) * t286;
t320 = t284 * t336;
t285 = sin(qJ(2));
t348 = t285 * t286;
t324 = qJD(2) * t348;
t245 = t267 * qJD(1) - t287 * t324 - t288 * t320 - t317;
t280 = sin(pkin(10));
t282 = cos(pkin(10));
t319 = t287 * t336;
t339 = qJD(1) * t289;
t301 = t284 * t339 + t319;
t318 = t284 * t335;
t341 = qJD(1) * t286;
t302 = t287 * t341 + t318;
t244 = -t284 * t324 + t301 * t288 - t302;
t281 = sin(pkin(6));
t283 = cos(pkin(6));
t338 = qJD(2) * t288;
t322 = t286 * t338;
t304 = t285 * t339 + t322;
t372 = -t244 * t283 + t304 * t281;
t217 = t245 * t282 + t372 * t280;
t337 = qJD(2) * t289;
t323 = t285 * t337;
t340 = qJD(1) * t288;
t242 = (-qJD(3) * t288 + qJD(1)) * t343 + (t323 + (-qJD(3) + t340) * t286) * t284;
t321 = t288 * t337;
t303 = -t285 * t341 + t321;
t376 = t242 * t283 + t303 * t281;
t342 = t289 * t284;
t344 = t287 * t288;
t265 = t286 * t344 - t342;
t345 = t286 * t288;
t264 = t284 * t345 + t343;
t309 = t264 * t283 - t281 * t348;
t235 = t265 * t282 - t309 * t280;
t349 = t284 * t285;
t306 = t281 * t288 + t283 * t349;
t347 = t285 * t287;
t328 = t282 * t347;
t297 = t306 * t280 - t328;
t225 = atan2(-t235, -t297);
t220 = sin(t225);
t221 = cos(t225);
t214 = -t220 * t235 - t221 * t297;
t212 = 0.1e1 / t214 ^ 2;
t266 = t286 * t287 - t288 * t342;
t346 = t285 * t289;
t308 = t266 * t283 + t281 * t346;
t240 = t267 * t282 + t308 * t280;
t234 = t240 ^ 2;
t210 = t212 * t234 + 0.1e1;
t305 = t286 * t340 + t323;
t243 = t305 * t287 + t288 * t318 - t301;
t216 = -t243 * t282 + t376 * t280;
t366 = t216 * t212;
t232 = t235 ^ 2;
t253 = 0.1e1 / t297 ^ 2;
t224 = t232 * t253 + 0.1e1;
t222 = 0.1e1 / t224;
t350 = t283 * t288;
t256 = t282 * t344 + (t281 * t285 - t284 * t350) * t280;
t352 = t280 * t287;
t307 = -t282 * t284 - t283 * t352;
t298 = qJD(3) * t307;
t231 = t256 * qJD(2) + t285 * t298;
t252 = 0.1e1 / t297;
t361 = t235 * t253;
t312 = t217 * t252 + t231 * t361;
t203 = t312 * t222;
t313 = t220 * t297 - t221 * t235;
t199 = t313 * t203 - t217 * t220 + t221 * t231;
t211 = 0.1e1 / t214;
t369 = t199 * t211 * t212;
t334 = 0.2e1 * (-t234 * t369 + t240 * t366) / t210 ^ 2;
t375 = t231 * t253;
t374 = t267 * t281;
t261 = -t266 * t281 + t283 * t346;
t257 = 0.1e1 / t261;
t258 = 0.1e1 / t261 ^ 2;
t371 = -0.2e1 * t235;
t363 = t252 * t375;
t368 = (t217 * t361 + t232 * t363) / t224 ^ 2;
t367 = t212 * t240;
t365 = t220 * t240;
t364 = t221 * t240;
t362 = t235 * t252;
t239 = -t267 * t280 + t308 * t282;
t360 = t239 * t258;
t356 = t257 * t280;
t355 = t257 * t282;
t263 = (-t281 * t349 + t350) * t289;
t354 = t258 * t263;
t353 = t280 * t283;
t351 = t282 * t283;
t333 = 0.2e1 * t369;
t215 = t243 * t280 + t376 * t282;
t233 = t239 ^ 2;
t228 = t233 * t258 + 0.1e1;
t229 = -t242 * t281 + t303 * t283;
t259 = t257 * t258;
t332 = 0.2e1 * (-t229 * t233 * t259 + t215 * t360) / t228 ^ 2;
t331 = -0.2e1 * t368;
t330 = 0.2e1 * t239 * t259;
t329 = t252 * t368;
t327 = t283 * t348;
t325 = t284 * t337;
t316 = t229 * t330;
t315 = t363 * t371;
t314 = t240 * t333;
t246 = -t264 * t282 - t265 * t353;
t262 = t307 * t285;
t311 = t246 * t252 + t262 * t361;
t249 = t286 * t328 + (-t281 * t345 - t284 * t327) * t280;
t310 = -t249 * t252 + t256 * t361;
t300 = t284 * t341 - t317;
t299 = -t220 + (-t221 * t362 + t220) * t222;
t260 = -t264 * t281 - t327;
t251 = t297 * t289;
t250 = (t280 * t347 + t306 * t282) * t289;
t248 = t266 * t282 - t267 * t353;
t247 = -t266 * t280 - t267 * t351;
t241 = (-t282 * t287 + t284 * t353) * t285 * qJD(3) + t307 * t338;
t237 = t265 * t280 + t309 * t282;
t230 = t297 * qJD(2) + t288 * t298;
t226 = 0.1e1 / t228;
t219 = -t244 * t282 - t245 * t353;
t218 = (t287 * t322 + (t287 * t339 - t320) * t285) * t282 + ((-t288 * t339 + t324) * t281 + (-t304 * t284 - t285 * t319) * t283) * t280;
t208 = 0.1e1 / t210;
t207 = t311 * t222;
t205 = t310 * t222;
t202 = t299 * t240;
t201 = t313 * t207 - t220 * t246 + t221 * t262;
t200 = t313 * t205 + t220 * t249 + t221 * t256;
t198 = t311 * t331 + (-t262 * t315 + t219 * t252 + (t217 * t262 + t231 * t246 + t235 * t241) * t253) * t222;
t197 = t310 * t331 + (-t256 * t315 - t218 * t252 + (t217 * t256 + t230 * t235 - t231 * t249) * t253) * t222;
t1 = [-0.2e1 * t240 * t329 + (t216 * t252 + t240 * t375) * t222, t197, t198, 0, 0, 0; t235 * t211 * t334 + (-t217 * t211 + (t199 * t235 - t202 * t216) * t212) * t208 + ((t202 * t333 - t299 * t366) * t208 + (t202 * t334 + (-(t203 * t222 * t362 + t331) * t365 - (-t329 * t371 - t203 + (t203 - t312) * t222) * t364) * t208) * t212) * t240 (t200 * t367 - t211 * t251) * t334 + (t200 * t314 + (-t251 * t199 - t200 * t216 - (-t197 * t235 - t205 * t217 + t230 + (t205 * t297 + t249) * t203) * t364 - (t197 * t297 - t205 * t231 + t218 + (t205 * t235 - t256) * t203) * t365) * t212 + ((t302 * t285 - t287 * t321) * t282 + (-t305 * t281 + (t303 * t284 + t285 * t317) * t283) * t280) * t211) * t208 (t201 * t367 - t211 * t248) * t334 + ((t242 * t282 + t243 * t353) * t211 + t201 * t314 + (-t248 * t199 - t201 * t216 - (-t198 * t235 - t207 * t217 + t241 + (t207 * t297 - t246) * t203) * t364 - (t198 * t297 - t207 * t231 - t219 + (t207 * t235 - t262) * t203) * t365) * t212) * t208, 0, 0, 0; (-t237 * t257 + t260 * t360) * t332 + (t260 * t316 + t245 * t356 - t372 * t355 + (-t237 * t229 - (-t244 * t281 - t304 * t283) * t239 - t260 * t215) * t258) * t226 (t239 * t354 - t250 * t257) * t332 + (-t215 * t354 + (-t250 * t258 + t263 * t330) * t229 + ((-t281 * t282 * t341 + t325 * t351 + t337 * t352) * t257 - (-t281 * t325 - t283 * t341) * t360) * t288 + (-(t300 * t281 - t283 * t337) * t360 - t302 * t356 + (-t281 * t337 - t300 * t283) * t355) * t285) * t226 (-t247 * t257 + t360 * t374) * t332 + ((-t242 * t280 + t243 * t351) * t257 + t316 * t374 + (-t247 * t229 + (-t215 * t267 + t239 * t243) * t281) * t258) * t226, 0, 0, 0;];
JaD_rot  = t1;
