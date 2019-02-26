% Zeitableitung der rotatorischen Teilmatrix der analytischen Jacobi-Matrix für Segment Nr. 6 (0=Basis) von
% S6PRPPRR2
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
%   pkin=[a2,a3,a4,a5,a6,alpha2,d2,d5,d6,theta1,theta3]';
%
% Output:
% JaD_rot [3x6]
%   Zeitableitung der rotatorischen Teilmatrix der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-02-26 19:45
% Revision: d75aae1ac561373cd3be920984c3df29a1c2ecc8 (2019-02-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JaD_rot = S6PRPPRR2_jacobiaD_rot_6_sym_varpar(qJ, qJD, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(11,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRPPRR2_jacobiaD_rot_6_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6PRPPRR2_jacobiaD_rot_6_sym_varpar: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6PRPPRR2_jacobiaD_rot_6_sym_varpar: pkin has to be [11x1] (double)');

%% Symbolic Calculation
% From jacobiaD_rot_6_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 19:45:13
% EndTime: 2019-02-26 19:45:14
% DurationCPUTime: 1.19s
% Computational Cost: add. (5642->114), mult. (16325->239), div. (559->12), fcn. (21250->15), ass. (0->110)
t277 = cos(pkin(6));
t272 = sin(pkin(11));
t275 = cos(pkin(11));
t280 = sin(qJ(2));
t283 = cos(qJ(2));
t296 = t272 * t283 + t275 * t280;
t265 = t296 * t277;
t261 = qJD(2) * t265;
t269 = t272 * t280 - t275 * t283;
t267 = t269 * qJD(2);
t273 = sin(pkin(10));
t276 = cos(pkin(10));
t240 = t261 * t276 - t267 * t273;
t274 = sin(pkin(6));
t282 = cos(qJ(5));
t279 = sin(qJ(5));
t264 = t269 * t277;
t291 = t276 * t264 + t273 * t296;
t290 = t291 * t279;
t216 = -qJD(5) * t290 + (qJD(5) * t274 * t276 + t240) * t282;
t312 = t274 * t279;
t237 = t276 * t312 + t282 * t291;
t234 = t237 ^ 2;
t262 = t269 * t274;
t298 = t262 * t282 - t277 * t279;
t254 = 0.1e1 / t298 ^ 2;
t229 = t234 * t254 + 0.1e1;
t227 = 0.1e1 / t229;
t257 = t262 * t279 + t277 * t282;
t263 = t296 * t274;
t258 = qJD(2) * t263;
t232 = qJD(5) * t257 - t258 * t282;
t253 = 0.1e1 / t298;
t315 = t237 * t254;
t197 = (-t216 * t253 - t232 * t315) * t227;
t230 = atan2(t237, -t298);
t225 = sin(t230);
t226 = cos(t230);
t300 = t225 * t298 + t226 * t237;
t193 = t197 * t300 + t225 * t216 + t226 * t232;
t209 = t225 * t237 - t226 * t298;
t206 = 0.1e1 / t209;
t207 = 0.1e1 / t209 ^ 2;
t329 = t193 * t206 * t207;
t251 = t264 * t273 - t276 * t296;
t235 = t251 * t282 + t273 * t312;
t328 = 0.2e1 * t235 * t329;
t316 = t232 * t253 * t254;
t327 = (t216 * t315 + t234 * t316) / t229 ^ 2;
t249 = t265 * t276 - t269 * t273;
t293 = -t249 * t253 + t263 * t315;
t326 = t282 * t293;
t311 = t274 * t282;
t236 = -t251 * t279 + t273 * t311;
t278 = sin(qJ(6));
t281 = cos(qJ(6));
t297 = -t265 * t273 - t269 * t276;
t222 = t236 * t281 + t278 * t297;
t218 = 0.1e1 / t222;
t219 = 0.1e1 / t222 ^ 2;
t242 = t261 * t273 + t267 * t276;
t213 = -qJD(5) * t235 - t242 * t279;
t260 = t277 * t267;
t268 = t296 * qJD(2);
t299 = t260 * t273 - t268 * t276;
t204 = qJD(6) * t222 + t213 * t278 - t281 * t299;
t221 = t236 * t278 - t281 * t297;
t217 = t221 ^ 2;
t212 = t217 * t219 + 0.1e1;
t320 = t219 * t221;
t307 = qJD(6) * t221;
t205 = t213 * t281 + t278 * t299 - t307;
t324 = t205 * t218 * t219;
t325 = (t204 * t320 - t217 * t324) / t212 ^ 2;
t323 = t207 * t235;
t214 = qJD(5) * t236 + t242 * t282;
t322 = t214 * t207;
t321 = t218 * t278;
t319 = t221 * t281;
t318 = t225 * t235;
t317 = t226 * t235;
t314 = t297 * t279;
t313 = t297 * t282;
t308 = qJD(5) * t279;
t233 = t235 ^ 2;
t203 = t207 * t233 + 0.1e1;
t306 = 0.2e1 * (-t233 * t329 + t235 * t322) / t203 ^ 2;
t305 = -0.2e1 * t325;
t303 = t221 * t324;
t302 = t237 * t316;
t301 = qJD(6) * t314 - t242;
t295 = t219 * t319 - t321;
t238 = t276 * t311 - t290;
t294 = t238 * t253 + t257 * t315;
t292 = qJD(5) * t313 + qJD(6) * t251 + t279 * t299;
t259 = t274 * t267;
t241 = -t260 * t276 - t268 * t273;
t231 = qJD(5) * t298 + t258 * t279;
t224 = t251 * t278 + t281 * t314;
t223 = -t251 * t281 + t278 * t314;
t215 = qJD(5) * t237 + t240 * t279;
t210 = 0.1e1 / t212;
t201 = 0.1e1 / t203;
t199 = t227 * t326;
t198 = t294 * t227;
t195 = (t225 * t249 - t226 * t263) * t282 + t300 * t199;
t194 = -t198 * t300 + t225 * t238 + t226 * t257;
t192 = 0.2e1 * t294 * t327 + (-0.2e1 * t257 * t302 + t215 * t253 + (-t216 * t257 - t231 * t237 - t232 * t238) * t254) * t227;
t190 = -0.2e1 * t326 * t327 + (-t293 * t308 + (0.2e1 * t263 * t302 - t241 * t253 + (t216 * t263 - t232 * t249 - t237 * t259) * t254) * t282) * t227;
t1 = [0, t190, 0, 0, t192, 0; 0 (t195 * t323 + t206 * t313) * t306 + ((-t282 * t299 + t297 * t308) * t206 + (-t322 + t328) * t195 + (t313 * t193 - (t263 * t308 + t190 * t237 + t199 * t216 + t259 * t282 + (t199 * t298 + t249 * t282) * t197) * t317 - (-t249 * t308 + t190 * t298 - t199 * t232 + t241 * t282 + (-t199 * t237 + t263 * t282) * t197) * t318) * t207) * t201, 0, 0 (t194 * t323 - t206 * t236) * t306 + (t194 * t328 + t213 * t206 + (-t236 * t193 - t194 * t214 - (t192 * t237 - t198 * t216 + t231 + (-t198 * t298 + t238) * t197) * t317 - (t192 * t298 + t198 * t232 - t215 + (t198 * t237 - t257) * t197) * t318) * t207) * t201, 0; 0, 0.2e1 * (-t218 * t223 + t224 * t320) * t325 + (0.2e1 * t224 * t303 + t301 * t218 * t281 + t292 * t321 + (t221 * t278 * t301 - t224 * t204 - t223 * t205 - t292 * t319) * t219) * t210, 0, 0, t295 * t235 * t305 + (t295 * t214 + ((-qJD(6) * t218 - 0.2e1 * t303) * t281 + (t204 * t281 + (t205 - t307) * t278) * t219) * t235) * t210, t305 + 0.2e1 * (t204 * t219 * t210 + (-t210 * t324 - t219 * t325) * t221) * t221;];
JaD_rot  = t1;
