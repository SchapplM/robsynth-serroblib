% Zeitableitung der rotatorischen Teilmatrix der analytischen Jacobi-Matrix für Segment Nr. 5 (0=Basis) von
% S6RRRRPR9
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
%   pkin=[a2,a3,a4,a5,a6,alpha2,d1,d2,d3,d4,d6,theta5]';
%
% Output:
% JaD_rot [3x6]
%   Zeitableitung der rotatorischen Teilmatrix der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-02-26 22:35
% Revision: d75aae1ac561373cd3be920984c3df29a1c2ecc8 (2019-02-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JaD_rot = S6RRRRPR9_jacobiaD_rot_5_sym_varpar(qJ, qJD, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(12,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRRPR9_jacobiaD_rot_5_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRRRPR9_jacobiaD_rot_5_sym_varpar: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [12 1]), ...
  'S6RRRRPR9_jacobiaD_rot_5_sym_varpar: pkin has to be [12x1] (double)');

%% Symbolic Calculation
% From jacobiaD_rot_5_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 22:35:13
% EndTime: 2019-02-26 22:35:15
% DurationCPUTime: 1.53s
% Computational Cost: add. (11272->143), mult. (16797->286), div. (959->12), fcn. (21298->13), ass. (0->131)
t265 = qJ(3) + qJ(4);
t262 = sin(t265);
t269 = cos(pkin(6));
t271 = cos(qJ(2));
t338 = sin(qJ(1));
t298 = t338 * t271;
t270 = sin(qJ(2));
t272 = cos(qJ(1));
t316 = t272 * t270;
t282 = -t269 * t316 - t298;
t263 = cos(t265);
t267 = sin(pkin(6));
t318 = t267 * t272;
t302 = t263 * t318;
t236 = -t262 * t282 + t302;
t320 = t267 * t270;
t304 = t262 * t320;
t246 = -t269 * t263 + t304;
t225 = atan2(-t236, t246);
t220 = sin(t225);
t221 = cos(t225);
t201 = -t220 * t236 + t221 * t246;
t199 = 0.1e1 / t201 ^ 2;
t299 = t338 * t270;
t291 = t269 * t299;
t315 = t272 * t271;
t251 = -t291 + t315;
t300 = t267 * t338;
t241 = t251 * t262 - t263 * t300;
t235 = t241 ^ 2;
t197 = t235 * t199 + 0.1e1;
t281 = -t269 * t298 - t316;
t231 = t282 * qJD(1) + t281 * qJD(2);
t264 = qJD(3) + qJD(4);
t288 = t264 * t300 + t231;
t297 = qJD(1) * t318;
t322 = t263 * t264;
t207 = t251 * t322 + t288 * t262 - t263 * t297;
t332 = t207 * t199;
t234 = t236 ^ 2;
t244 = 0.1e1 / t246 ^ 2;
t224 = t234 * t244 + 0.1e1;
t222 = 0.1e1 / t224;
t295 = t338 * qJD(2);
t233 = -qJD(1) * t291 - t270 * t295 + (qJD(2) * t269 + qJD(1)) * t315;
t257 = t262 * t318;
t296 = t338 * qJD(1);
t290 = t267 * t296;
t209 = t233 * t262 - t264 * t257 - t263 * t290 - t282 * t322;
t313 = qJD(2) * t271;
t284 = t264 * t269 + t267 * t313;
t303 = t263 * t320;
t228 = t284 * t262 + t264 * t303;
t243 = 0.1e1 / t246;
t326 = t236 * t244;
t287 = -t209 * t243 + t228 * t326;
t191 = t287 * t222;
t289 = -t220 * t246 - t221 * t236;
t186 = t289 * t191 - t220 * t209 + t221 * t228;
t198 = 0.1e1 / t201;
t200 = t198 * t199;
t336 = t186 * t200;
t312 = 0.2e1 * (-t235 * t336 + t241 * t332) / t197 ^ 2;
t342 = t228 * t244;
t301 = t269 * t315;
t248 = -t299 + t301;
t319 = t267 * t271;
t283 = -t243 * t248 + t319 * t326;
t341 = t262 * t283;
t210 = (t264 * t282 + t290) * t262 + t233 * t263 - t264 * t302;
t242 = t251 * t263 + t262 * t300;
t268 = cos(pkin(12));
t266 = sin(pkin(12));
t324 = t281 * t266;
t219 = t242 * t268 - t324;
t213 = 0.1e1 / t219;
t214 = 0.1e1 / t219 ^ 2;
t340 = -0.2e1 * t236;
t339 = 0.2e1 * t241;
t328 = t243 * t342;
t335 = (t209 * t326 - t234 * t328) / t224 ^ 2;
t334 = t199 * t241;
t208 = t288 * t263 + (-t251 * t264 + t297) * t262;
t230 = -qJD(1) * t301 - t272 * t313 + (t269 * t295 + t296) * t270;
t203 = t208 * t268 - t230 * t266;
t333 = t203 * t213 * t214;
t323 = t281 * t268;
t218 = t242 * t266 + t323;
t331 = t214 * t218;
t330 = t220 * t241;
t329 = t221 * t241;
t327 = t236 * t243;
t325 = t281 * t262;
t321 = t266 * t213;
t317 = t268 * t218;
t314 = qJD(2) * t270;
t202 = t208 * t266 + t230 * t268;
t212 = t218 ^ 2;
t206 = t212 * t214 + 0.1e1;
t311 = 0.2e1 * (t202 * t331 - t212 * t333) / t206 ^ 2;
t310 = -0.2e1 * t335;
t309 = t200 * t339;
t308 = t243 * t335;
t307 = t199 * t330;
t306 = t199 * t329;
t305 = t218 * t333;
t294 = 0.2e1 * t305;
t293 = t328 * t340;
t238 = -t263 * t282 - t257;
t286 = t230 * t263 - t264 * t325;
t247 = t269 * t262 + t303;
t285 = -t238 * t243 + t247 * t326;
t279 = -t220 + (t221 * t327 + t220) * t222;
t232 = t281 * qJD(1) + t282 * qJD(2);
t229 = t284 * t263 - t264 * t304;
t227 = t251 * t266 + t263 * t323;
t226 = -t251 * t268 + t263 * t324;
t217 = -t238 * t268 + t248 * t266;
t216 = -t238 * t266 - t248 * t268;
t204 = 0.1e1 / t206;
t195 = 0.1e1 / t197;
t194 = t222 * t341;
t193 = t285 * t222;
t190 = t279 * t241;
t188 = (-t220 * t248 + t221 * t319) * t262 + t289 * t194;
t187 = t289 * t193 - t220 * t238 + t221 * t247;
t184 = t285 * t310 + (t247 * t293 - t210 * t243 + (t209 * t247 + t228 * t238 + t229 * t236) * t244) * t222;
t183 = t310 * t341 + (t283 * t322 + (t293 * t319 - t232 * t243 + (t228 * t248 + (t209 * t271 - t236 * t314) * t267) * t244) * t262) * t222;
t182 = (-t214 * t317 + t321) * t241 * t311 + (-0.2e1 * t241 * t268 * t305 - t207 * t321 + (t207 * t317 + (t202 * t268 + t203 * t266) * t241) * t214) * t204;
t181 = (t187 * t334 - t198 * t242) * t312 + (t187 * t186 * t309 + t208 * t198 + (-t242 * t186 - t187 * t207 - (-t184 * t236 - t193 * t209 + t229 + (-t193 * t246 - t238) * t191) * t329 - (-t184 * t246 - t193 * t228 - t210 + (t193 * t236 - t247) * t191) * t330) * t199) * t195;
t1 = [t308 * t339 + (-t207 * t243 + t241 * t342) * t222, t183, t184, t184, 0, 0; t236 * t198 * t312 + (-t209 * t198 + (t186 * t236 - t190 * t207) * t199) * t195 + (t190 * t199 * t312 + (0.2e1 * t190 * t336 - (-t191 * t222 * t327 + t310) * t307 - (t308 * t340 - t191 + (t191 - t287) * t222) * t306 - t279 * t332) * t195) * t241 (t188 * t334 - t198 * t325) * t312 + (-t188 * t332 + (t230 * t262 + t281 * t322) * t198 + (t188 * t309 - t199 * t325) * t186 - (-t183 * t236 - t194 * t209 + (-t262 * t314 + t271 * t322) * t267 + (-t194 * t246 - t248 * t262) * t191) * t306 - (-t248 * t322 - t183 * t246 - t194 * t228 - t232 * t262 + (t194 * t236 - t262 * t319) * t191) * t307) * t195, t181, t181, 0, 0; (-t213 * t216 + t217 * t331) * t311 + ((-t210 * t266 - t232 * t268) * t213 + t217 * t294 + (-t216 * t203 - (-t210 * t268 + t232 * t266) * t218 - t217 * t202) * t214) * t204 (-t213 * t226 + t227 * t331) * t311 + ((-t231 * t268 + t286 * t266) * t213 + t227 * t294 + (-t226 * t203 - (t231 * t266 + t286 * t268) * t218 - t227 * t202) * t214) * t204, t182, t182, 0, 0;];
JaD_rot  = t1;
