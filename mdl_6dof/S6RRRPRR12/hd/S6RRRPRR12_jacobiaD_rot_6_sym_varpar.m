% Zeitableitung der rotatorischen Teilmatrix der analytischen Jacobi-Matrix für Segment Nr. 6 (0=Basis) von
% S6RRRPRR12
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
% Datum: 2019-02-26 22:22
% Revision: d75aae1ac561373cd3be920984c3df29a1c2ecc8 (2019-02-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JaD_rot = S6RRRPRR12_jacobiaD_rot_6_sym_varpar(qJ, qJD, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(12,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRPRR12_jacobiaD_rot_6_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRRPRR12_jacobiaD_rot_6_sym_varpar: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [12 1]), ...
  'S6RRRPRR12_jacobiaD_rot_6_sym_varpar: pkin has to be [12x1] (double)');

%% Symbolic Calculation
% From jacobiaD_rot_6_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 22:22:36
% EndTime: 2019-02-26 22:22:37
% DurationCPUTime: 1.39s
% Computational Cost: add. (5847->150), mult. (14155->298), div. (744->12), fcn. (17850->13), ass. (0->132)
t285 = cos(pkin(6));
t287 = sin(qJ(2));
t361 = sin(qJ(1));
t324 = t361 * t287;
t310 = t285 * t324;
t318 = t361 * qJD(2);
t289 = cos(qJ(2));
t290 = cos(qJ(1));
t339 = t290 * t289;
t284 = sin(pkin(6));
t341 = t284 * t290;
t366 = -qJD(1) * t310 - t287 * t318 + (qJD(2) * t285 + qJD(1)) * t339 - qJD(3) * t341;
t286 = sin(qJ(3));
t288 = cos(qJ(3));
t323 = t361 * t289;
t340 = t290 * t287;
t302 = -t285 * t340 - t323;
t253 = -t286 * t302 + t288 * t341;
t343 = t284 * t287;
t264 = -t285 * t288 + t286 * t343;
t244 = atan2(-t253, t264);
t239 = sin(t244);
t240 = cos(t244);
t220 = -t239 * t253 + t240 * t264;
t218 = 0.1e1 / t220 ^ 2;
t269 = -t310 + t339;
t325 = t284 * t361;
t301 = -t269 * t286 + t288 * t325;
t250 = t301 ^ 2;
t216 = t250 * t218 + 0.1e1;
t300 = -t285 * t323 - t340;
t246 = qJD(1) * t302 + qJD(2) * t300;
t259 = t269 * t288 + t286 * t325;
t322 = qJD(1) * t341;
t224 = qJD(3) * t259 + t246 * t286 - t288 * t322;
t354 = t224 * t218;
t249 = t253 ^ 2;
t262 = 0.1e1 / t264 ^ 2;
t243 = t249 * t262 + 0.1e1;
t241 = 0.1e1 / t243;
t319 = t361 * qJD(1);
t309 = t284 * t319;
t336 = qJD(3) * t288;
t226 = t286 * t366 - t288 * t309 - t302 * t336;
t265 = t285 * t286 + t288 * t343;
t337 = qJD(2) * t289;
t321 = t284 * t337;
t251 = qJD(3) * t265 + t286 * t321;
t261 = 0.1e1 / t264;
t348 = t253 * t262;
t306 = -t226 * t261 + t251 * t348;
t208 = t306 * t241;
t307 = -t239 * t264 - t240 * t253;
t203 = t208 * t307 - t239 * t226 + t240 * t251;
t217 = 0.1e1 / t220;
t219 = t217 * t218;
t359 = t203 * t219;
t335 = 0.2e1 * (-t250 * t359 - t301 * t354) / t216 ^ 2;
t365 = t251 * t262;
t326 = t285 * t339;
t266 = -t324 + t326;
t342 = t284 * t289;
t303 = -t261 * t266 + t342 * t348;
t364 = t286 * t303;
t227 = t286 * (qJD(3) * t302 + t309) + t366 * t288;
t282 = pkin(12) + qJ(5) + qJ(6);
t280 = sin(t282);
t281 = cos(t282);
t236 = t259 * t281 - t280 * t300;
t230 = 0.1e1 / t236;
t231 = 0.1e1 / t236 ^ 2;
t363 = -0.2e1 * t253;
t362 = -0.2e1 * t301;
t245 = -qJD(1) * t326 - t290 * t337 + (t285 * t318 + t319) * t287;
t283 = qJD(5) + qJD(6);
t313 = t259 * t283 + t245;
t225 = qJD(3) * t301 + t246 * t288 + t286 * t322;
t315 = -t283 * t300 + t225;
t211 = t280 * t315 + t281 * t313;
t235 = t259 * t280 + t281 * t300;
t229 = t235 ^ 2;
t223 = t229 * t231 + 0.1e1;
t353 = t231 * t235;
t212 = -t280 * t313 + t281 * t315;
t356 = t212 * t230 * t231;
t358 = (t211 * t353 - t229 * t356) / t223 ^ 2;
t350 = t261 * t365;
t357 = (t226 * t348 - t249 * t350) / t243 ^ 2;
t355 = t218 * t301;
t352 = t239 * t301;
t351 = t240 * t301;
t349 = t253 * t261;
t347 = t300 * t286;
t346 = t300 * t288;
t345 = t280 * t230;
t344 = t281 * t235;
t338 = qJD(2) * t287;
t334 = -0.2e1 * t358;
t333 = 0.2e1 * t358;
t332 = -0.2e1 * t357;
t331 = t219 * t362;
t330 = t261 * t357;
t329 = t218 * t352;
t328 = t218 * t351;
t327 = t235 * t356;
t317 = 0.2e1 * t327;
t316 = t350 * t363;
t314 = t266 * t283 - t227;
t247 = qJD(1) * t300 + qJD(2) * t302;
t255 = -t286 * t341 - t288 * t302;
t312 = -t255 * t283 - t247;
t308 = -t283 * t346 + t246;
t305 = t344 * t231 - t345;
t304 = -t255 * t261 + t265 * t348;
t298 = -t239 + (t240 * t349 + t239) * t241;
t297 = -qJD(3) * t347 + t245 * t288 + t269 * t283;
t252 = -qJD(3) * t264 + t288 * t321;
t238 = t269 * t280 + t281 * t346;
t237 = -t269 * t281 + t280 * t346;
t234 = -t255 * t281 + t266 * t280;
t233 = -t255 * t280 - t266 * t281;
t221 = 0.1e1 / t223;
t214 = 0.1e1 / t216;
t213 = t241 * t364;
t210 = t304 * t241;
t207 = t298 * t301;
t205 = (-t239 * t266 + t240 * t342) * t286 + t307 * t213;
t204 = t210 * t307 - t239 * t255 + t240 * t265;
t202 = t304 * t332 + (t265 * t316 - t227 * t261 + (t226 * t265 + t251 * t255 + t252 * t253) * t262) * t241;
t200 = t332 * t364 + (t303 * t336 + (t316 * t342 - t247 * t261 + (t251 * t266 + (t226 * t289 - t253 * t338) * t284) * t262) * t286) * t241;
t199 = t334 + 0.2e1 * (t211 * t231 * t221 + (-t221 * t356 - t231 * t358) * t235) * t235;
t1 = [t330 * t362 + (-t224 * t261 - t301 * t365) * t241, t200, t202, 0, 0, 0; t253 * t217 * t335 + (-t226 * t217 + (t203 * t253 + t207 * t224) * t218) * t214 - (-t207 * t218 * t335 + (-0.2e1 * t207 * t359 + (-t208 * t241 * t349 + t332) * t329 + (t330 * t363 - t208 + (t208 - t306) * t241) * t328 - t298 * t354) * t214) * t301 (-t205 * t355 - t217 * t347) * t335 + (-t205 * t354 + (t245 * t286 + t300 * t336) * t217 + (t205 * t331 - t218 * t347) * t203 + (-t200 * t253 - t213 * t226 + (-t286 * t338 + t289 * t336) * t284 + (-t213 * t264 - t266 * t286) * t208) * t328 + (-t266 * t336 - t200 * t264 - t213 * t251 - t247 * t286 + (t213 * t253 - t286 * t342) * t208) * t329) * t214 (-t204 * t355 - t217 * t259) * t335 + (t204 * t203 * t331 + t225 * t217 + (-t259 * t203 - t204 * t224 + (-t202 * t253 - t210 * t226 + t252 + (-t210 * t264 - t255) * t208) * t351 + (-t202 * t264 - t210 * t251 - t227 + (t210 * t253 - t265) * t208) * t352) * t218) * t214, 0, 0, 0; (-t230 * t233 + t234 * t353) * t333 + ((t280 * t314 + t281 * t312) * t230 + t234 * t317 + (-t233 * t212 - (-t312 * t280 + t281 * t314) * t235 - t234 * t211) * t231) * t221 (-t230 * t237 + t238 * t353) * t333 + (t238 * t317 - t308 * t230 * t281 + t297 * t345 + (-t235 * t280 * t308 - t238 * t211 - t237 * t212 - t297 * t344) * t231) * t221, -t305 * t301 * t334 + (t305 * t224 - ((-t230 * t283 - 0.2e1 * t327) * t281 + (t211 * t281 + (-t235 * t283 + t212) * t280) * t231) * t301) * t221, 0, t199, t199;];
JaD_rot  = t1;
