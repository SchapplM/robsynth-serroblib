% Zeitableitung der rotatorischen Teilmatrix der analytischen Jacobi-Matrix für Segment Nr. 6 (0=Basis) von
% S6PRPRPR3
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
%   pkin=[a2,a3,a4,a5,a6,alpha2,d2,d4,d6,theta1,theta3]';
%
% Output:
% JaD_rot [3x6]
%   Zeitableitung der rotatorischen Teilmatrix der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-02-26 19:47
% Revision: d75aae1ac561373cd3be920984c3df29a1c2ecc8 (2019-02-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JaD_rot = S6PRPRPR3_jacobiaD_rot_6_sym_varpar(qJ, qJD, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(11,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRPRPR3_jacobiaD_rot_6_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6PRPRPR3_jacobiaD_rot_6_sym_varpar: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6PRPRPR3_jacobiaD_rot_6_sym_varpar: pkin has to be [11x1] (double)');

%% Symbolic Calculation
% From jacobiaD_rot_6_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 19:47:39
% EndTime: 2019-02-26 19:47:40
% DurationCPUTime: 1.18s
% Computational Cost: add. (5642->115), mult. (16325->236), div. (559->12), fcn. (21250->15), ass. (0->112)
t278 = sin(qJ(4));
t281 = cos(qJ(4));
t276 = cos(pkin(6));
t273 = sin(pkin(11));
t275 = cos(pkin(11));
t279 = sin(qJ(2));
t282 = cos(qJ(2));
t302 = t282 * t273 + t279 * t275;
t263 = t302 * t276;
t266 = t279 * t273 - t282 * t275;
t274 = sin(pkin(10));
t334 = cos(pkin(10));
t297 = t334 * t263 - t274 * t266;
t333 = sin(pkin(6));
t305 = t334 * t333;
t236 = t297 * t278 + t281 * t305;
t294 = t266 * t276;
t260 = qJD(2) * t294;
t265 = t302 * qJD(2);
t242 = -t334 * t260 - t274 * t265;
t214 = -t236 * qJD(4) + t242 * t281;
t237 = -t278 * t305 + t297 * t281;
t234 = t237 ^ 2;
t309 = t282 * t333;
t310 = t279 * t333;
t289 = -t273 * t309 - t275 * t310;
t255 = t276 * t278 - t281 * t289;
t252 = 0.1e1 / t255 ^ 2;
t229 = t234 * t252 + 0.1e1;
t227 = 0.1e1 / t229;
t254 = t276 * t281 + t278 * t289;
t261 = -t273 * t310 + t275 * t309;
t259 = t261 * qJD(2);
t232 = t254 * qJD(4) + t259 * t281;
t251 = 0.1e1 / t255;
t322 = t237 * t252;
t197 = (-t214 * t251 + t232 * t322) * t227;
t230 = atan2(-t237, t255);
t225 = sin(t230);
t226 = cos(t230);
t304 = -t225 * t255 - t226 * t237;
t193 = t304 * t197 - t225 * t214 + t226 * t232;
t209 = -t225 * t237 + t226 * t255;
t206 = 0.1e1 / t209;
t207 = 0.1e1 / t209 ^ 2;
t339 = t193 * t206 * t207;
t249 = t274 * t294 - t302 * t334;
t277 = sin(qJ(6));
t280 = cos(qJ(6));
t295 = -t274 * t263 - t334 * t266;
t311 = t274 * t333;
t293 = -t278 * t295 + t281 * t311;
t303 = t249 * t277 - t280 * t293;
t338 = t303 * qJD(6);
t240 = t278 * t311 + t281 * t295;
t337 = 0.2e1 * t240 * t339;
t247 = -t274 * t302 - t334 * t294;
t298 = -t247 * t251 + t261 * t322;
t336 = t281 * t298;
t323 = t232 * t251 * t252;
t335 = -0.2e1 * (t214 * t322 - t234 * t323) / t229 ^ 2;
t320 = t249 * t280;
t222 = -t277 * t293 - t320;
t218 = 0.1e1 / t222;
t219 = 0.1e1 / t222 ^ 2;
t296 = t274 * t260 - t334 * t265;
t215 = t240 * qJD(4) + t278 * t296;
t264 = t266 * qJD(2);
t292 = t276 * t265;
t243 = t334 * t264 + t274 * t292;
t204 = t222 * qJD(6) - t215 * t280 - t243 * t277;
t217 = t303 ^ 2;
t212 = t217 * t219 + 0.1e1;
t327 = t219 * t303;
t205 = t215 * t277 - t243 * t280 + t338;
t331 = t205 * t218 * t219;
t332 = (-t204 * t327 - t217 * t331) / t212 ^ 2;
t330 = t207 * t240;
t216 = t293 * qJD(4) + t281 * t296;
t329 = t216 * t207;
t328 = t218 * t280;
t326 = t303 * t277;
t325 = t225 * t240;
t324 = t226 * t240;
t321 = t249 * t278;
t319 = t249 * t281;
t315 = qJD(4) * t278;
t235 = t240 ^ 2;
t203 = t235 * t207 + 0.1e1;
t314 = 0.2e1 * (-t235 * t339 + t240 * t329) / t203 ^ 2;
t313 = 0.2e1 * t332;
t308 = -0.2e1 * t303 * t331;
t307 = -0.2e1 * t237 * t323;
t306 = qJD(6) * t321 + t296;
t300 = -t219 * t326 + t328;
t299 = t236 * t251 + t254 * t322;
t290 = qJD(4) * t319 - qJD(6) * t295 + t243 * t278;
t258 = t289 * qJD(2);
t241 = t274 * t264 - t334 * t292;
t231 = -t255 * qJD(4) - t259 * t278;
t224 = t277 * t321 + t280 * t295;
t223 = t277 * t295 - t278 * t320;
t213 = t237 * qJD(4) + t242 * t278;
t210 = 0.1e1 / t212;
t201 = 0.1e1 / t203;
t199 = t227 * t336;
t198 = t299 * t227;
t195 = (-t225 * t247 + t226 * t261) * t281 + t304 * t199;
t194 = t304 * t198 + t225 * t236 + t226 * t254;
t192 = t299 * t335 + (t254 * t307 + t213 * t251 + (t214 * t254 + t231 * t237 - t232 * t236) * t252) * t227;
t190 = t335 * t336 + (-t298 * t315 + (t261 * t307 - t241 * t251 + (t214 * t261 + t232 * t247 + t237 * t258) * t252) * t281) * t227;
t1 = [0, t190, 0, t192, 0, 0; 0 (t195 * t330 - t206 * t319) * t314 + ((t243 * t281 - t249 * t315) * t206 + (-t329 + t337) * t195 + (-t319 * t193 - (-t261 * t315 - t190 * t237 - t199 * t214 + t258 * t281 + (-t199 * t255 - t247 * t281) * t197) * t324 - (t247 * t315 - t190 * t255 - t199 * t232 - t241 * t281 + (t199 * t237 - t261 * t281) * t197) * t325) * t207) * t201, 0 (t194 * t330 - t206 * t293) * t314 + (t194 * t337 - t215 * t206 + (-t293 * t193 - t194 * t216 - (-t192 * t237 - t198 * t214 + t231 + (-t198 * t255 + t236) * t197) * t324 - (-t192 * t255 - t198 * t232 + t213 + (t198 * t237 - t254) * t197) * t325) * t207) * t201, 0, 0; 0 (-t218 * t223 - t224 * t327) * t313 + (t224 * t308 + t306 * t218 * t277 - t290 * t328 + (t280 * t303 * t306 - t224 * t204 - t223 * t205 + t290 * t326) * t219) * t210, 0, t300 * t240 * t313 + (-t300 * t216 + ((qJD(6) * t218 + t308) * t277 + (-t204 * t277 + (t205 + t338) * t280) * t219) * t240) * t210, 0, -0.2e1 * t332 - 0.2e1 * (t204 * t219 * t210 - (-t210 * t331 - t219 * t332) * t303) * t303;];
JaD_rot  = t1;
