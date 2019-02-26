% Zeitableitung der rotatorischen Teilmatrix der analytischen Jacobi-Matrix für Segment Nr. 6 (0=Basis) von
% S6RRRPPR6
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
%   pkin=[a2,a3,a4,a5,a6,alpha2,d1,d2,d3,d6,theta4]';
%
% Output:
% JaD_rot [3x6]
%   Zeitableitung der rotatorischen Teilmatrix der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-02-26 22:06
% Revision: d75aae1ac561373cd3be920984c3df29a1c2ecc8 (2019-02-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JaD_rot = S6RRRPPR6_jacobiaD_rot_6_sym_varpar(qJ, qJD, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(11,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRPPR6_jacobiaD_rot_6_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRRPPR6_jacobiaD_rot_6_sym_varpar: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6RRRPPR6_jacobiaD_rot_6_sym_varpar: pkin has to be [11x1] (double)');

%% Symbolic Calculation
% From jacobiaD_rot_6_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 22:06:23
% EndTime: 2019-02-26 22:06:25
% DurationCPUTime: 1.44s
% Computational Cost: add. (8428->151), mult. (13478->298), div. (726->12), fcn. (17045->13), ass. (0->129)
t271 = sin(qJ(1));
t274 = cos(qJ(1));
t268 = cos(pkin(6));
t291 = qJD(2) * t268 + qJD(1);
t270 = sin(qJ(2));
t318 = t271 * t270;
t300 = t268 * t318;
t267 = sin(pkin(6));
t308 = qJD(3) * t267;
t312 = qJD(2) * t270;
t273 = cos(qJ(2));
t314 = t274 * t273;
t345 = -qJD(1) * t300 - t271 * t312 - t274 * t308 + t291 * t314;
t254 = -t300 + t314;
t266 = qJ(3) + pkin(11);
t264 = sin(t266);
t265 = cos(t266);
t322 = t267 * t271;
t242 = t254 * t264 - t265 * t322;
t272 = cos(qJ(6));
t315 = t274 * t270;
t317 = t271 * t273;
t253 = t268 * t317 + t315;
t269 = sin(qJ(6));
t325 = t253 * t269;
t287 = t242 * t272 - t325;
t344 = qJD(6) * t287;
t252 = t268 * t315 + t317;
t320 = t267 * t274;
t238 = t252 * t265 - t264 * t320;
t323 = t267 * t270;
t250 = t268 * t264 + t265 * t323;
t225 = atan2(-t238, t250);
t212 = sin(t225);
t213 = cos(t225);
t203 = -t212 * t238 + t213 * t250;
t201 = 0.1e1 / t203 ^ 2;
t243 = t254 * t265 + t264 * t322;
t236 = t243 ^ 2;
t197 = t236 * t201 + 0.1e1;
t230 = -qJD(1) * t252 - qJD(2) * t253;
t313 = qJD(1) * t267;
t297 = t274 * t313;
t310 = qJD(3) * t264;
t208 = t264 * t297 - t254 * t310 + (t271 * t308 + t230) * t265;
t333 = t208 * t201;
t235 = t238 ^ 2;
t247 = 0.1e1 / t250 ^ 2;
t224 = t235 * t247 + 0.1e1;
t218 = 0.1e1 / t224;
t292 = t345 * t265;
t298 = t271 * t313;
t210 = -t252 * t310 + t264 * t298 + t292;
t249 = -t264 * t323 + t268 * t265;
t311 = qJD(2) * t273;
t296 = t267 * t311;
t234 = qJD(3) * t249 + t265 * t296;
t246 = 0.1e1 / t250;
t327 = t238 * t247;
t286 = -t210 * t246 + t234 * t327;
t191 = t286 * t218;
t289 = -t212 * t250 - t213 * t238;
t186 = t191 * t289 - t212 * t210 + t213 * t234;
t200 = 0.1e1 / t203;
t202 = t200 * t201;
t338 = t186 * t202;
t307 = 0.2e1 * (-t236 * t338 + t243 * t333) / t197 ^ 2;
t343 = t234 * t247;
t299 = t268 * t314;
t251 = t299 - t318;
t321 = t267 * t273;
t283 = -t246 * t251 + t321 * t327;
t342 = t265 * t283;
t324 = t253 * t272;
t223 = t242 * t269 + t324;
t215 = 0.1e1 / t223;
t216 = 0.1e1 / t223 ^ 2;
t341 = -0.2e1 * t238;
t340 = 0.2e1 * t243;
t207 = qJD(3) * t243 + t230 * t264 - t265 * t297;
t229 = -qJD(1) * t299 - t274 * t311 + t291 * t318;
t198 = qJD(6) * t223 - t207 * t272 - t229 * t269;
t214 = t287 ^ 2;
t206 = t214 * t216 + 0.1e1;
t330 = t216 * t287;
t199 = t207 * t269 - t229 * t272 + t344;
t335 = t199 * t215 * t216;
t337 = (-t198 * t330 - t214 * t335) / t206 ^ 2;
t329 = t246 * t343;
t336 = (t210 * t327 - t235 * t329) / t224 ^ 2;
t334 = t201 * t243;
t332 = t212 * t243;
t331 = t213 * t243;
t328 = t238 * t246;
t326 = t253 * t265;
t319 = t269 * t287;
t316 = t272 * t215;
t309 = qJD(3) * t265;
t306 = 0.2e1 * t337;
t305 = -0.2e1 * t336;
t304 = t202 * t340;
t303 = t246 * t336;
t302 = t201 * t332;
t301 = t201 * t331;
t294 = -0.2e1 * t287 * t335;
t293 = t329 * t341;
t290 = -qJD(6) * t253 * t264 + t230;
t237 = t252 * t264 + t265 * t320;
t288 = -t237 * t272 - t251 * t269;
t221 = -t237 * t269 + t251 * t272;
t285 = -t216 * t319 + t316;
t284 = t237 * t246 + t249 * t327;
t282 = -t212 + (t213 * t328 + t212) * t218;
t209 = t252 * t309 + t345 * t264 - t265 * t298;
t281 = qJD(6) * t254 - t229 * t264 + t253 * t309;
t233 = -qJD(3) * t250 - t264 * t296;
t231 = -qJD(1) * t253 - qJD(2) * t252;
t227 = t254 * t272 - t264 * t325;
t226 = t254 * t269 + t264 * t324;
t204 = 0.1e1 / t206;
t195 = 0.1e1 / t197;
t194 = t218 * t342;
t192 = t284 * t218;
t190 = t282 * t243;
t188 = (-t212 * t251 + t213 * t321) * t265 + t289 * t194;
t187 = t192 * t289 + t212 * t237 + t213 * t249;
t185 = t284 * t305 + (t249 * t293 + t209 * t246 + (t210 * t249 + t233 * t238 - t234 * t237) * t247) * t218;
t183 = t305 * t342 + (-t283 * t310 + (t293 * t321 - t231 * t246 + (t234 * t251 + (t210 * t273 - t238 * t312) * t267) * t247) * t265) * t218;
t1 = [t303 * t340 + (-t208 * t246 + t243 * t343) * t218, t183, t185, 0, 0, 0; t238 * t200 * t307 + (((qJD(3) * t252 - t298) * t264 - t292) * t200 + (t186 * t238 - t190 * t208) * t201) * t195 + (t190 * t201 * t307 + (0.2e1 * t190 * t338 - (-t191 * t218 * t328 + t305) * t302 - (t303 * t341 - t191 + (t191 - t286) * t218) * t301 - t282 * t333) * t195) * t243 (t188 * t334 + t200 * t326) * t307 + (-t188 * t333 + (t229 * t265 + t253 * t310) * t200 + (t188 * t304 + t201 * t326) * t186 - (-t183 * t238 - t194 * t210 + (-t265 * t312 - t273 * t310) * t267 + (-t194 * t250 - t251 * t265) * t191) * t301 - (t251 * t310 - t183 * t250 - t194 * t234 - t231 * t265 + (t194 * t238 - t265 * t321) * t191) * t302) * t195 (t187 * t334 + t200 * t242) * t307 + (t187 * t186 * t304 - t207 * t200 + (t242 * t186 - t187 * t208 - (-t185 * t238 - t192 * t210 + t233 + (-t192 * t250 + t237) * t191) * t331 - (-t185 * t250 - t192 * t234 + t209 + (t192 * t238 - t249) * t191) * t332) * t201) * t195, 0, 0, 0; (t215 * t288 - t221 * t330) * t306 + ((qJD(6) * t221 + t209 * t272 + t231 * t269) * t215 + t221 * t294 + (t288 * t199 + (qJD(6) * t288 - t209 * t269 + t231 * t272) * t287 - t221 * t198) * t216) * t204 (-t215 * t226 - t227 * t330) * t306 + (t227 * t294 + t290 * t215 * t269 + t281 * t316 + (t272 * t287 * t290 - t227 * t198 - t226 * t199 - t281 * t319) * t216) * t204, t285 * t243 * t306 + (-t285 * t208 + ((qJD(6) * t215 + t294) * t269 + (-t198 * t269 + (t199 + t344) * t272) * t216) * t243) * t204, 0, 0, -0.2e1 * t337 - 0.2e1 * (t198 * t216 * t204 - (-t204 * t335 - t216 * t337) * t287) * t287;];
JaD_rot  = t1;
