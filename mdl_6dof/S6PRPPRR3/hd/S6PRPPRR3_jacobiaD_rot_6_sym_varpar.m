% Zeitableitung der rotatorischen Teilmatrix der analytischen Jacobi-Matrix für Segment Nr. 6 (0=Basis) von
% S6PRPPRR3
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
%   pkin=[a2,a3,a4,a5,a6,alpha2,d2,d5,d6,theta1,theta4]';
%
% Output:
% JaD_rot [3x6]
%   Zeitableitung der rotatorischen Teilmatrix der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-02-26 19:46
% Revision: d75aae1ac561373cd3be920984c3df29a1c2ecc8 (2019-02-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JaD_rot = S6PRPPRR3_jacobiaD_rot_6_sym_varpar(qJ, qJD, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(11,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRPPRR3_jacobiaD_rot_6_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6PRPPRR3_jacobiaD_rot_6_sym_varpar: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6PRPPRR3_jacobiaD_rot_6_sym_varpar: pkin has to be [11x1] (double)');

%% Symbolic Calculation
% From jacobiaD_rot_6_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 19:45:50
% EndTime: 2019-02-26 19:45:51
% DurationCPUTime: 1.19s
% Computational Cost: add. (5642->117), mult. (16325->241), div. (559->12), fcn. (21250->15), ass. (0->114)
t283 = cos(pkin(6));
t286 = sin(qJ(2));
t336 = cos(pkin(10));
t310 = t336 * t286;
t280 = sin(pkin(10));
t289 = cos(qJ(2));
t320 = t280 * t289;
t273 = t283 * t310 + t320;
t279 = sin(pkin(11));
t282 = cos(pkin(11));
t309 = t336 * t289;
t321 = t280 * t286;
t297 = t283 * t309 - t321;
t252 = t273 * t282 - t279 * t297;
t285 = sin(qJ(5));
t288 = cos(qJ(5));
t281 = sin(pkin(6));
t311 = t281 * t336;
t241 = t252 * t288 + t285 * t311;
t268 = t297 * qJD(2);
t269 = t273 * qJD(2);
t244 = t268 * t282 + t269 * t279;
t217 = t241 * qJD(5) + t244 * t285;
t298 = -t252 * t285 + t288 * t311;
t237 = t298 ^ 2;
t267 = (t279 * t289 - t282 * t286) * t281;
t260 = -t267 * t285 + t283 * t288;
t258 = 0.1e1 / t260 ^ 2;
t233 = t237 * t258 + 0.1e1;
t231 = 0.1e1 / t233;
t261 = -t267 * t288 - t283 * t285;
t266 = (t279 * t286 + t282 * t289) * t281;
t264 = qJD(2) * t266;
t235 = t261 * qJD(5) + t264 * t285;
t257 = 0.1e1 / t260;
t325 = t298 * t258;
t201 = (-t217 * t257 - t235 * t325) * t231;
t234 = atan2(t298, t260);
t229 = sin(t234);
t230 = cos(t234);
t304 = -t229 * t260 + t230 * t298;
t197 = t304 * t201 - t229 * t217 + t230 * t235;
t213 = t229 * t298 + t230 * t260;
t210 = 0.1e1 / t213;
t211 = 0.1e1 / t213 ^ 2;
t340 = t197 * t210 * t211;
t275 = -t283 * t321 + t309;
t299 = -t283 * t320 - t310;
t256 = t275 * t282 - t279 * t299;
t322 = t280 * t281;
t242 = t256 * t285 + t288 * t322;
t339 = 0.2e1 * t242 * t340;
t251 = t273 * t279 + t282 * t297;
t300 = -t251 * t257 - t266 * t325;
t338 = t285 * t300;
t326 = t235 * t257 * t258;
t337 = -0.2e1 * (-t217 * t325 - t237 * t326) / t233 ^ 2;
t313 = t285 * t322;
t243 = t256 * t288 - t313;
t284 = sin(qJ(6));
t287 = cos(qJ(6));
t306 = t275 * t279 + t282 * t299;
t226 = t243 * t287 + t284 * t306;
t222 = 0.1e1 / t226;
t223 = 0.1e1 / t226 ^ 2;
t270 = t299 * qJD(2);
t271 = t275 * qJD(2);
t247 = t270 * t282 + t271 * t279;
t220 = -t242 * qJD(5) + t247 * t288;
t307 = t270 * t279 - t271 * t282;
t208 = t226 * qJD(6) + t220 * t284 - t287 * t307;
t225 = t243 * t284 - t287 * t306;
t221 = t225 ^ 2;
t216 = t221 * t223 + 0.1e1;
t330 = t223 * t225;
t318 = qJD(6) * t225;
t209 = t220 * t287 + t284 * t307 - t318;
t334 = t209 * t222 * t223;
t335 = (t208 * t330 - t221 * t334) / t216 ^ 2;
t333 = t211 * t242;
t319 = qJD(5) * t288;
t219 = -qJD(5) * t313 + t247 * t285 + t256 * t319;
t332 = t219 * t211;
t331 = t222 * t284;
t329 = t225 * t287;
t328 = t229 * t242;
t327 = t230 * t242;
t324 = t306 * t285;
t323 = t306 * t288;
t238 = t242 ^ 2;
t207 = t238 * t211 + 0.1e1;
t317 = 0.2e1 * (-t238 * t340 + t242 * t332) / t207 ^ 2;
t316 = -0.2e1 * t335;
t314 = t225 * t334;
t308 = 0.2e1 * t298 * t326;
t305 = qJD(6) * t323 + t247;
t302 = t223 * t329 - t331;
t301 = -t241 * t257 - t261 * t325;
t296 = -qJD(5) * t324 - qJD(6) * t256 + t288 * t307;
t263 = qJD(2) * t267;
t245 = t268 * t279 - t269 * t282;
t236 = -t260 * qJD(5) + t264 * t288;
t228 = -t256 * t284 + t287 * t323;
t227 = t256 * t287 + t284 * t323;
t218 = t298 * qJD(5) + t244 * t288;
t214 = 0.1e1 / t216;
t205 = 0.1e1 / t207;
t203 = t231 * t338;
t202 = t301 * t231;
t199 = (-t229 * t251 + t230 * t266) * t285 + t304 * t203;
t198 = t304 * t202 - t229 * t241 + t230 * t261;
t196 = t301 * t337 + (t261 * t308 - t218 * t257 + (t217 * t261 + t235 * t241 - t236 * t298) * t258) * t231;
t194 = t337 * t338 + (t300 * t319 + (t266 * t308 - t245 * t257 + (t217 * t266 + t235 * t251 - t263 * t298) * t258) * t285) * t231;
t1 = [0, t194, 0, 0, t196, 0; 0 (t199 * t333 - t210 * t324) * t317 + ((t285 * t307 + t306 * t319) * t210 + (-t332 + t339) * t199 + (-t324 * t197 - (t266 * t319 + t194 * t298 - t203 * t217 + t263 * t285 + (-t203 * t260 - t251 * t285) * t201) * t327 - (-t251 * t319 - t194 * t260 - t203 * t235 - t245 * t285 + (-t203 * t298 - t266 * t285) * t201) * t328) * t211) * t205, 0, 0 (t198 * t333 - t210 * t243) * t317 + (t198 * t339 + t220 * t210 + (-t243 * t197 - t198 * t219 - (t196 * t298 - t202 * t217 + t236 + (-t202 * t260 - t241) * t201) * t327 - (-t196 * t260 - t202 * t235 - t218 + (-t202 * t298 - t261) * t201) * t328) * t211) * t205, 0; 0, 0.2e1 * (-t222 * t227 + t228 * t330) * t335 + (0.2e1 * t228 * t314 + t305 * t222 * t287 + t296 * t331 + (t305 * t225 * t284 - t228 * t208 - t227 * t209 - t296 * t329) * t223) * t214, 0, 0, t302 * t242 * t316 + (t302 * t219 + ((-qJD(6) * t222 - 0.2e1 * t314) * t287 + (t208 * t287 + (t209 - t318) * t284) * t223) * t242) * t214, t316 + 0.2e1 * (t208 * t223 * t214 + (-t214 * t334 - t223 * t335) * t225) * t225;];
JaD_rot  = t1;
