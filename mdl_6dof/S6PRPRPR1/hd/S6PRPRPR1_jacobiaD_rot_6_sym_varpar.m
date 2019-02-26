% Zeitableitung der rotatorischen Teilmatrix der analytischen Jacobi-Matrix für Segment Nr. 6 (0=Basis) von
% S6PRPRPR1
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
%   pkin=[a2,a3,a4,a5,a6,alpha2,d2,d4,d6,theta1,theta3,theta5]';
%
% Output:
% JaD_rot [3x6]
%   Zeitableitung der rotatorischen Teilmatrix der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-02-26 19:46
% Revision: d75aae1ac561373cd3be920984c3df29a1c2ecc8 (2019-02-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JaD_rot = S6PRPRPR1_jacobiaD_rot_6_sym_varpar(qJ, qJD, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(12,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRPRPR1_jacobiaD_rot_6_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6PRPRPR1_jacobiaD_rot_6_sym_varpar: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [12 1]), ...
  'S6PRPRPR1_jacobiaD_rot_6_sym_varpar: pkin has to be [12x1] (double)');

%% Symbolic Calculation
% From jacobiaD_rot_6_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 19:46:28
% EndTime: 2019-02-26 19:46:29
% DurationCPUTime: 1.25s
% Computational Cost: add. (8444->115), mult. (16325->235), div. (559->12), fcn. (21250->15), ass. (0->112)
t284 = sin(pkin(11));
t287 = cos(pkin(11));
t291 = sin(qJ(2));
t293 = cos(qJ(2));
t274 = t291 * t284 - t293 * t287;
t289 = cos(pkin(6));
t302 = t274 * t289;
t268 = qJD(2) * t302;
t307 = t293 * t284 + t291 * t287;
t273 = t307 * qJD(2);
t285 = sin(pkin(10));
t288 = cos(pkin(10));
t249 = -t288 * t268 - t285 * t273;
t271 = t307 * t289;
t255 = t288 * t271 - t285 * t274;
t283 = qJ(4) + pkin(12);
t281 = sin(t283);
t286 = sin(pkin(6));
t323 = t286 * t288;
t313 = t281 * t323;
t282 = cos(t283);
t319 = qJD(4) * t282;
t221 = -qJD(4) * t313 + t249 * t281 + t255 * t319;
t243 = t255 * t281 + t282 * t323;
t241 = t243 ^ 2;
t270 = t307 * t286;
t262 = t270 * t281 - t289 * t282;
t260 = 0.1e1 / t262 ^ 2;
t235 = t241 * t260 + 0.1e1;
t233 = 0.1e1 / t235;
t263 = t270 * t282 + t289 * t281;
t269 = t274 * t286;
t267 = qJD(2) * t269;
t239 = t263 * qJD(4) - t267 * t281;
t259 = 0.1e1 / t262;
t328 = t243 * t260;
t205 = (-t221 * t259 + t239 * t328) * t233;
t236 = atan2(-t243, t262);
t231 = sin(t236);
t232 = cos(t236);
t310 = -t231 * t262 - t232 * t243;
t201 = t310 * t205 - t231 * t221 + t232 * t239;
t215 = -t231 * t243 + t232 * t262;
t212 = 0.1e1 / t215;
t213 = 0.1e1 / t215 ^ 2;
t342 = t201 * t212 * t213;
t308 = -t285 * t271 - t288 * t274;
t324 = t285 * t286;
t303 = -t281 * t308 + t282 * t324;
t341 = -0.2e1 * t303 * t342;
t254 = -t285 * t307 - t288 * t302;
t304 = -t254 * t259 - t269 * t328;
t340 = t281 * t304;
t329 = t239 * t259 * t260;
t339 = -0.2e1 * (t221 * t328 - t241 * t329) / t235 ^ 2;
t247 = t281 * t324 + t282 * t308;
t292 = cos(qJ(6));
t257 = t285 * t302 - t288 * t307;
t290 = sin(qJ(6));
t326 = t257 * t290;
t230 = t247 * t292 - t326;
t226 = 0.1e1 / t230;
t227 = 0.1e1 / t230 ^ 2;
t309 = t285 * t268 - t288 * t273;
t224 = t303 * qJD(4) + t282 * t309;
t272 = t274 * qJD(2);
t301 = t289 * t273;
t250 = t288 * t272 + t285 * t301;
t216 = t230 * qJD(6) + t224 * t290 + t250 * t292;
t325 = t257 * t292;
t229 = t247 * t290 + t325;
t225 = t229 ^ 2;
t220 = t225 * t227 + 0.1e1;
t333 = t227 * t229;
t318 = qJD(6) * t229;
t217 = t224 * t292 - t250 * t290 - t318;
t336 = t217 * t226 * t227;
t338 = (t216 * t333 - t225 * t336) / t220 ^ 2;
t337 = t213 * t303;
t223 = t247 * qJD(4) + t281 * t309;
t335 = t223 * t213;
t334 = t226 * t290;
t332 = t229 * t292;
t331 = t231 * t303;
t330 = t232 * t303;
t327 = t257 * t281;
t242 = t303 ^ 2;
t211 = t242 * t213 + 0.1e1;
t317 = 0.2e1 * (-t242 * t342 - t303 * t335) / t211 ^ 2;
t316 = -0.2e1 * t338;
t314 = t229 * t336;
t312 = -0.2e1 * t243 * t329;
t311 = qJD(6) * t257 * t282 - t309;
t306 = t227 * t332 - t334;
t245 = t255 * t282 - t313;
t305 = -t245 * t259 + t263 * t328;
t300 = -qJD(4) * t327 + qJD(6) * t308 + t250 * t282;
t266 = t286 * t273;
t248 = t285 * t272 - t288 * t301;
t240 = -t262 * qJD(4) - t267 * t282;
t238 = t282 * t325 + t290 * t308;
t237 = t282 * t326 - t292 * t308;
t222 = -t243 * qJD(4) + t249 * t282;
t218 = 0.1e1 / t220;
t209 = 0.1e1 / t211;
t207 = t233 * t340;
t206 = t305 * t233;
t203 = (-t231 * t254 - t232 * t269) * t281 + t310 * t207;
t202 = t310 * t206 - t231 * t245 + t232 * t263;
t199 = t305 * t339 + (t263 * t312 - t222 * t259 + (t221 * t263 + t239 * t245 + t240 * t243) * t260) * t233;
t198 = t339 * t340 + (t304 * t319 + (-t269 * t312 - t248 * t259 + (-t221 * t269 + t239 * t254 - t243 * t266) * t260) * t281) * t233;
t1 = [0, t198, 0, t199, 0, 0; 0 (-t203 * t337 - t212 * t327) * t317 + ((t250 * t281 + t257 * t319) * t212 + (-t335 + t341) * t203 + (-t327 * t201 + (-t269 * t319 - t198 * t243 - t207 * t221 - t266 * t281 + (-t207 * t262 - t254 * t281) * t205) * t330 + (-t254 * t319 - t198 * t262 - t207 * t239 - t248 * t281 + (t207 * t243 + t269 * t281) * t205) * t331) * t213) * t209, 0 (-t202 * t337 - t212 * t247) * t317 + (t202 * t341 + t224 * t212 + (-t247 * t201 - t202 * t223 + (-t199 * t243 - t206 * t221 + t240 + (-t206 * t262 - t245) * t205) * t330 + (-t199 * t262 - t206 * t239 - t222 + (t206 * t243 - t263) * t205) * t331) * t213) * t209, 0, 0; 0, 0.2e1 * (-t226 * t237 + t238 * t333) * t338 + (0.2e1 * t238 * t314 + t311 * t226 * t292 + t300 * t334 + (t311 * t229 * t290 - t238 * t216 - t237 * t217 - t300 * t332) * t227) * t218, 0, -t306 * t303 * t316 + (t306 * t223 - ((-qJD(6) * t226 - 0.2e1 * t314) * t292 + (t216 * t292 + (t217 - t318) * t290) * t227) * t303) * t218, 0, t316 + 0.2e1 * (t216 * t227 * t218 + (-t218 * t336 - t227 * t338) * t229) * t229;];
JaD_rot  = t1;
