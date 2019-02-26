% Zeitableitung der rotatorischen Teilmatrix der analytischen Jacobi-Matrix für Segment Nr. 6 (0=Basis) von
% S6PRPRRP1
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
%   pkin=[a2,a3,a4,a5,a6,alpha2,d2,d4,d5,theta1,theta3]';
%
% Output:
% JaD_rot [3x6]
%   Zeitableitung der rotatorischen Teilmatrix der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-02-26 19:50
% Revision: d75aae1ac561373cd3be920984c3df29a1c2ecc8 (2019-02-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JaD_rot = S6PRPRRP1_jacobiaD_rot_6_sym_varpar(qJ, qJD, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(11,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRPRRP1_jacobiaD_rot_6_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6PRPRRP1_jacobiaD_rot_6_sym_varpar: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6PRPRRP1_jacobiaD_rot_6_sym_varpar: pkin has to be [11x1] (double)');

%% Symbolic Calculation
% From jacobiaD_rot_6_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 19:50:19
% EndTime: 2019-02-26 19:50:20
% DurationCPUTime: 1.26s
% Computational Cost: add. (5642->114), mult. (16325->235), div. (559->12), fcn. (21250->15), ass. (0->110)
t281 = sin(pkin(11));
t284 = cos(pkin(11));
t289 = sin(qJ(2));
t292 = cos(qJ(2));
t274 = t289 * t281 - t292 * t284;
t286 = cos(pkin(6));
t301 = t274 * t286;
t268 = qJD(2) * t301;
t306 = t292 * t281 + t289 * t284;
t273 = t306 * qJD(2);
t282 = sin(pkin(10));
t285 = cos(pkin(10));
t249 = -t285 * t268 - t282 * t273;
t271 = t306 * t286;
t255 = t285 * t271 - t282 * t274;
t288 = sin(qJ(4));
t283 = sin(pkin(6));
t323 = t283 * t288;
t312 = t285 * t323;
t291 = cos(qJ(4));
t318 = qJD(4) * t291;
t221 = -qJD(4) * t312 + t249 * t288 + t255 * t318;
t322 = t283 * t291;
t243 = t255 * t288 + t285 * t322;
t241 = t243 ^ 2;
t270 = t306 * t283;
t262 = t270 * t288 - t286 * t291;
t260 = 0.1e1 / t262 ^ 2;
t237 = t241 * t260 + 0.1e1;
t235 = 0.1e1 / t237;
t263 = t270 * t291 + t286 * t288;
t269 = t274 * t283;
t267 = qJD(2) * t269;
t239 = qJD(4) * t263 - t267 * t288;
t259 = 0.1e1 / t262;
t326 = t243 * t260;
t205 = (-t221 * t259 + t239 * t326) * t235;
t238 = atan2(-t243, t262);
t233 = sin(t238);
t234 = cos(t238);
t309 = -t233 * t262 - t234 * t243;
t201 = t205 * t309 - t233 * t221 + t234 * t239;
t217 = -t233 * t243 + t234 * t262;
t214 = 0.1e1 / t217;
t215 = 0.1e1 / t217 ^ 2;
t340 = t201 * t214 * t215;
t307 = -t282 * t271 - t285 * t274;
t302 = t282 * t322 - t288 * t307;
t339 = -0.2e1 * t302 * t340;
t254 = -t282 * t306 - t285 * t301;
t303 = -t254 * t259 - t269 * t326;
t338 = t288 * t303;
t327 = t239 * t259 * t260;
t337 = -0.2e1 * (t221 * t326 - t241 * t327) / t237 ^ 2;
t247 = t282 * t323 + t291 * t307;
t257 = t282 * t301 - t285 * t306;
t287 = sin(qJ(5));
t290 = cos(qJ(5));
t230 = t247 * t290 - t257 * t287;
t226 = 0.1e1 / t230;
t227 = 0.1e1 / t230 ^ 2;
t308 = t282 * t268 - t285 * t273;
t224 = qJD(4) * t302 + t291 * t308;
t272 = t274 * qJD(2);
t300 = t286 * t273;
t250 = t285 * t272 + t282 * t300;
t212 = qJD(5) * t230 + t224 * t287 + t250 * t290;
t229 = t247 * t287 + t257 * t290;
t225 = t229 ^ 2;
t220 = t225 * t227 + 0.1e1;
t331 = t227 * t229;
t317 = qJD(5) * t229;
t213 = t224 * t290 - t250 * t287 - t317;
t335 = t213 * t226 * t227;
t336 = (t212 * t331 - t225 * t335) / t220 ^ 2;
t334 = t215 * t302;
t223 = qJD(4) * t247 + t288 * t308;
t333 = t223 * t215;
t332 = t226 * t287;
t330 = t229 * t290;
t329 = t233 * t302;
t328 = t234 * t302;
t325 = t257 * t288;
t324 = t257 * t291;
t242 = t302 ^ 2;
t211 = t242 * t215 + 0.1e1;
t316 = 0.2e1 * (-t242 * t340 - t302 * t333) / t211 ^ 2;
t315 = -0.2e1 * t336;
t313 = t229 * t335;
t311 = -0.2e1 * t243 * t327;
t310 = qJD(5) * t324 - t308;
t305 = t330 * t227 - t332;
t245 = t255 * t291 - t312;
t304 = -t245 * t259 + t263 * t326;
t299 = -qJD(4) * t325 + qJD(5) * t307 + t250 * t291;
t266 = t283 * t273;
t248 = t282 * t272 - t285 * t300;
t240 = -qJD(4) * t262 - t267 * t291;
t232 = t287 * t307 + t290 * t324;
t231 = t287 * t324 - t290 * t307;
t222 = -qJD(4) * t243 + t249 * t291;
t218 = 0.1e1 / t220;
t209 = 0.1e1 / t211;
t207 = t235 * t338;
t206 = t304 * t235;
t203 = (-t233 * t254 - t234 * t269) * t288 + t309 * t207;
t202 = t206 * t309 - t233 * t245 + t234 * t263;
t200 = t304 * t337 + (t263 * t311 - t222 * t259 + (t221 * t263 + t239 * t245 + t240 * t243) * t260) * t235;
t198 = t337 * t338 + (t303 * t318 + (-t269 * t311 - t248 * t259 + (-t221 * t269 + t239 * t254 - t243 * t266) * t260) * t288) * t235;
t1 = [0, t198, 0, t200, 0, 0; 0 (-t203 * t334 - t214 * t325) * t316 + ((t250 * t288 + t257 * t318) * t214 + (-t333 + t339) * t203 + (-t325 * t201 + (-t269 * t318 - t198 * t243 - t207 * t221 - t266 * t288 + (-t207 * t262 - t254 * t288) * t205) * t328 + (-t254 * t318 - t198 * t262 - t207 * t239 - t248 * t288 + (t207 * t243 + t269 * t288) * t205) * t329) * t215) * t209, 0 (-t202 * t334 - t214 * t247) * t316 + (t202 * t339 + t224 * t214 + (-t247 * t201 - t202 * t223 + (-t200 * t243 - t206 * t221 + t240 + (-t206 * t262 - t245) * t205) * t328 + (-t200 * t262 - t206 * t239 - t222 + (t206 * t243 - t263) * t205) * t329) * t215) * t209, 0, 0; 0, 0.2e1 * (-t226 * t231 + t232 * t331) * t336 + (0.2e1 * t232 * t313 + t310 * t226 * t290 + t299 * t332 + (t229 * t287 * t310 - t232 * t212 - t231 * t213 - t299 * t330) * t227) * t218, 0, -t305 * t302 * t315 + (t305 * t223 - ((-qJD(5) * t226 - 0.2e1 * t313) * t290 + (t212 * t290 + (t213 - t317) * t287) * t227) * t302) * t218, t315 + 0.2e1 * (t212 * t227 * t218 + (-t218 * t335 - t227 * t336) * t229) * t229, 0;];
JaD_rot  = t1;
