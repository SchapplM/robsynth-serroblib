% Zeitableitung der rotatorischen Teilmatrix der analytischen Jacobi-Matrix für Segment Nr. 6 (0=Basis) von
% S6RRRRRP9
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
% Datum: 2019-02-26 22:44
% Revision: d75aae1ac561373cd3be920984c3df29a1c2ecc8 (2019-02-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JaD_rot = S6RRRRRP9_jacobiaD_rot_6_sym_varpar(qJ, qJD, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(11,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRRRP9_jacobiaD_rot_6_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRRRRP9_jacobiaD_rot_6_sym_varpar: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6RRRRRP9_jacobiaD_rot_6_sym_varpar: pkin has to be [11x1] (double)');

%% Symbolic Calculation
% From jacobiaD_rot_6_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 22:44:26
% EndTime: 2019-02-26 22:44:27
% DurationCPUTime: 1.40s
% Computational Cost: add. (5342->150), mult. (14155->298), div. (744->12), fcn. (17850->13), ass. (0->132)
t287 = cos(pkin(6));
t289 = sin(qJ(2));
t363 = sin(qJ(1));
t326 = t363 * t289;
t312 = t287 * t326;
t320 = t363 * qJD(2);
t291 = cos(qJ(2));
t292 = cos(qJ(1));
t341 = t292 * t291;
t286 = sin(pkin(6));
t343 = t286 * t292;
t368 = -qJD(1) * t312 - t289 * t320 + (qJD(2) * t287 + qJD(1)) * t341 - qJD(3) * t343;
t288 = sin(qJ(3));
t290 = cos(qJ(3));
t325 = t363 * t291;
t342 = t292 * t289;
t304 = -t287 * t342 - t325;
t255 = -t288 * t304 + t290 * t343;
t345 = t286 * t289;
t266 = -t287 * t290 + t288 * t345;
t244 = atan2(-t255, t266);
t239 = sin(t244);
t240 = cos(t244);
t222 = -t239 * t255 + t240 * t266;
t220 = 0.1e1 / t222 ^ 2;
t271 = -t312 + t341;
t327 = t286 * t363;
t303 = -t271 * t288 + t290 * t327;
t252 = t303 ^ 2;
t218 = t252 * t220 + 0.1e1;
t302 = -t287 * t325 - t342;
t248 = qJD(1) * t304 + qJD(2) * t302;
t261 = t271 * t290 + t288 * t327;
t324 = qJD(1) * t343;
t226 = qJD(3) * t261 + t248 * t288 - t290 * t324;
t356 = t226 * t220;
t251 = t255 ^ 2;
t264 = 0.1e1 / t266 ^ 2;
t243 = t251 * t264 + 0.1e1;
t241 = 0.1e1 / t243;
t321 = t363 * qJD(1);
t311 = t286 * t321;
t338 = qJD(3) * t290;
t228 = t368 * t288 - t290 * t311 - t304 * t338;
t267 = t287 * t288 + t290 * t345;
t339 = qJD(2) * t291;
t323 = t286 * t339;
t253 = qJD(3) * t267 + t288 * t323;
t263 = 0.1e1 / t266;
t350 = t255 * t264;
t308 = -t228 * t263 + t253 * t350;
t210 = t308 * t241;
t309 = -t239 * t266 - t240 * t255;
t205 = t210 * t309 - t239 * t228 + t240 * t253;
t219 = 0.1e1 / t222;
t221 = t219 * t220;
t361 = t205 * t221;
t337 = 0.2e1 * (-t252 * t361 - t303 * t356) / t218 ^ 2;
t367 = t253 * t264;
t328 = t287 * t341;
t268 = -t326 + t328;
t344 = t286 * t291;
t305 = -t263 * t268 + t344 * t350;
t366 = t288 * t305;
t229 = t288 * (qJD(3) * t304 + t311) + t368 * t290;
t285 = qJ(4) + qJ(5);
t282 = sin(t285);
t283 = cos(t285);
t238 = t261 * t283 - t282 * t302;
t232 = 0.1e1 / t238;
t233 = 0.1e1 / t238 ^ 2;
t365 = -0.2e1 * t255;
t364 = -0.2e1 * t303;
t247 = -qJD(1) * t328 - t292 * t339 + (t287 * t320 + t321) * t289;
t284 = qJD(4) + qJD(5);
t315 = t261 * t284 + t247;
t227 = qJD(3) * t303 + t248 * t290 + t288 * t324;
t317 = -t284 * t302 + t227;
t213 = t282 * t317 + t283 * t315;
t237 = t261 * t282 + t283 * t302;
t231 = t237 ^ 2;
t225 = t231 * t233 + 0.1e1;
t355 = t233 * t237;
t214 = -t282 * t315 + t283 * t317;
t358 = t214 * t232 * t233;
t360 = (t213 * t355 - t231 * t358) / t225 ^ 2;
t352 = t263 * t367;
t359 = (t228 * t350 - t251 * t352) / t243 ^ 2;
t357 = t220 * t303;
t354 = t239 * t303;
t353 = t240 * t303;
t351 = t255 * t263;
t349 = t302 * t288;
t348 = t302 * t290;
t347 = t282 * t232;
t346 = t283 * t237;
t340 = qJD(2) * t289;
t336 = -0.2e1 * t360;
t335 = 0.2e1 * t360;
t334 = -0.2e1 * t359;
t333 = t221 * t364;
t332 = t263 * t359;
t331 = t220 * t354;
t330 = t220 * t353;
t329 = t237 * t358;
t319 = 0.2e1 * t329;
t318 = t352 * t365;
t316 = t268 * t284 - t229;
t249 = qJD(1) * t302 + qJD(2) * t304;
t257 = -t288 * t343 - t290 * t304;
t314 = -t257 * t284 - t249;
t310 = -t284 * t348 + t248;
t307 = t233 * t346 - t347;
t306 = -t257 * t263 + t267 * t350;
t300 = -t239 + (t240 * t351 + t239) * t241;
t299 = -qJD(3) * t349 + t247 * t290 + t271 * t284;
t254 = -qJD(3) * t266 + t290 * t323;
t246 = t271 * t282 + t283 * t348;
t245 = -t271 * t283 + t282 * t348;
t236 = -t257 * t283 + t268 * t282;
t235 = -t257 * t282 - t268 * t283;
t223 = 0.1e1 / t225;
t216 = 0.1e1 / t218;
t215 = t241 * t366;
t212 = t306 * t241;
t209 = t300 * t303;
t207 = (-t239 * t268 + t240 * t344) * t288 + t309 * t215;
t206 = t212 * t309 - t239 * t257 + t240 * t267;
t204 = t306 * t334 + (t267 * t318 - t229 * t263 + (t228 * t267 + t253 * t257 + t254 * t255) * t264) * t241;
t202 = t334 * t366 + (t305 * t338 + (t318 * t344 - t249 * t263 + (t253 * t268 + (t228 * t291 - t255 * t340) * t286) * t264) * t288) * t241;
t201 = t336 + 0.2e1 * (t213 * t233 * t223 + (-t223 * t358 - t233 * t360) * t237) * t237;
t1 = [t332 * t364 + (-t226 * t263 - t303 * t367) * t241, t202, t204, 0, 0, 0; t255 * t219 * t337 + (-t228 * t219 + (t205 * t255 + t209 * t226) * t220) * t216 - (-t209 * t220 * t337 + (-0.2e1 * t209 * t361 + (-t210 * t241 * t351 + t334) * t331 + (t332 * t365 - t210 + (t210 - t308) * t241) * t330 - t300 * t356) * t216) * t303 (-t207 * t357 - t219 * t349) * t337 + (-t207 * t356 + (t247 * t288 + t302 * t338) * t219 + (t207 * t333 - t220 * t349) * t205 + (-t202 * t255 - t215 * t228 + (-t288 * t340 + t291 * t338) * t286 + (-t215 * t266 - t268 * t288) * t210) * t330 + (-t268 * t338 - t202 * t266 - t215 * t253 - t249 * t288 + (t215 * t255 - t288 * t344) * t210) * t331) * t216 (-t206 * t357 - t219 * t261) * t337 + (t206 * t205 * t333 + t227 * t219 + (-t261 * t205 - t206 * t226 + (-t204 * t255 - t212 * t228 + t254 + (-t212 * t266 - t257) * t210) * t353 + (-t204 * t266 - t212 * t253 - t229 + (t212 * t255 - t267) * t210) * t354) * t220) * t216, 0, 0, 0; (-t232 * t235 + t236 * t355) * t335 + ((t282 * t316 + t283 * t314) * t232 + t236 * t319 + (-t235 * t214 - (-t282 * t314 + t283 * t316) * t237 - t236 * t213) * t233) * t223 (-t232 * t245 + t246 * t355) * t335 + (t246 * t319 - t310 * t232 * t283 + t299 * t347 + (-t237 * t282 * t310 - t246 * t213 - t245 * t214 - t299 * t346) * t233) * t223, -t307 * t303 * t336 + (t307 * t226 - ((-t232 * t284 - 0.2e1 * t329) * t283 + (t213 * t283 + (-t237 * t284 + t214) * t282) * t233) * t303) * t223, t201, t201, 0;];
JaD_rot  = t1;
