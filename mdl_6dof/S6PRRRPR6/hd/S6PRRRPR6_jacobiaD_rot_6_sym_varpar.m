% Zeitableitung der rotatorischen Teilmatrix der analytischen Jacobi-Matrix für Segment Nr. 6 (0=Basis) von
% S6PRRRPR6
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
%   pkin=[a2,a3,a4,a5,a6,alpha2,d2,d3,d4,d6,theta1]';
%
% Output:
% JaD_rot [3x6]
%   Zeitableitung der rotatorischen Teilmatrix der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-02-26 20:13
% Revision: d75aae1ac561373cd3be920984c3df29a1c2ecc8 (2019-02-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JaD_rot = S6PRRRPR6_jacobiaD_rot_6_sym_varpar(qJ, qJD, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(11,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRRRPR6_jacobiaD_rot_6_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6PRRRPR6_jacobiaD_rot_6_sym_varpar: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6PRRRPR6_jacobiaD_rot_6_sym_varpar: pkin has to be [11x1] (double)');

%% Symbolic Calculation
% From jacobiaD_rot_6_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 20:13:35
% EndTime: 2019-02-26 20:13:37
% DurationCPUTime: 1.28s
% Computational Cost: add. (4638->132), mult. (13597->274), div. (578->12), fcn. (17397->15), ass. (0->120)
t346 = -qJD(4) + qJD(6);
t277 = sin(pkin(11));
t279 = cos(pkin(11));
t284 = sin(qJ(2));
t280 = cos(pkin(6));
t288 = cos(qJ(2));
t322 = t280 * t288;
t268 = -t277 * t284 + t279 * t322;
t261 = t268 * qJD(2);
t323 = t280 * t284;
t269 = t277 * t288 + t279 * t323;
t283 = sin(qJ(3));
t278 = sin(pkin(6));
t326 = t278 * t283;
t308 = t279 * t326;
t287 = cos(qJ(3));
t319 = qJD(3) * t287;
t237 = -qJD(3) * t308 + t261 * t283 + t269 * t319;
t325 = t278 * t287;
t253 = t269 * t283 + t279 * t325;
t251 = t253 ^ 2;
t297 = -t280 * t287 + t284 * t326;
t266 = 0.1e1 / t297 ^ 2;
t247 = t251 * t266 + 0.1e1;
t245 = 0.1e1 / t247;
t273 = -t280 * t283 - t284 * t325;
t320 = qJD(2) * t288;
t307 = t278 * t320;
t257 = t273 * qJD(3) - t283 * t307;
t265 = 0.1e1 / t297;
t331 = t253 * t266;
t206 = (-t237 * t265 - t257 * t331) * t245;
t248 = atan2(t253, -t297);
t243 = sin(t248);
t244 = cos(t248);
t302 = t243 * t297 + t244 * t253;
t201 = t302 * t206 + t243 * t237 + t244 * t257;
t219 = t243 * t253 - t244 * t297;
t216 = 0.1e1 / t219;
t217 = 0.1e1 / t219 ^ 2;
t345 = t201 * t216 * t217;
t309 = t277 * t323;
t271 = t279 * t288 - t309;
t255 = -t271 * t283 + t277 * t325;
t344 = 0.2e1 * t255 * t345;
t329 = t257 * t265 * t266;
t343 = (t237 * t331 + t251 * t329) / t247 ^ 2;
t324 = t278 * t288;
t310 = t253 * t324;
t296 = -t265 * t268 + t266 * t310;
t342 = t283 * t296;
t286 = cos(qJ(4));
t341 = t346 * t286;
t282 = sin(qJ(4));
t340 = t346 * t282;
t270 = t277 * t322 + t279 * t284;
t263 = t270 * qJD(2);
t240 = t255 * qJD(3) - t263 * t287;
t256 = t271 * t287 + t277 * t326;
t242 = t256 * t286 + t270 * t282;
t264 = -qJD(2) * t309 + t279 * t320;
t220 = t242 * qJD(4) + t240 * t282 - t264 * t286;
t241 = t256 * t282 - t270 * t286;
t221 = -t241 * qJD(4) + t240 * t286 + t264 * t282;
t281 = sin(qJ(6));
t285 = cos(qJ(6));
t303 = t241 * t285 - t242 * t281;
t205 = t303 * qJD(6) + t220 * t281 + t221 * t285;
t231 = t241 * t281 + t242 * t285;
t225 = 0.1e1 / t231;
t226 = 0.1e1 / t231 ^ 2;
t204 = t231 * qJD(6) - t220 * t285 + t221 * t281;
t224 = t303 ^ 2;
t209 = t224 * t226 + 0.1e1;
t336 = t226 * t303;
t227 = t225 * t226;
t338 = t205 * t227;
t339 = (-t204 * t336 - t224 * t338) / t209 ^ 2;
t337 = t217 * t255;
t299 = t281 * t282 + t285 * t286;
t235 = t299 * t255;
t335 = t226 * t235;
t239 = -t256 * qJD(3) + t263 * t283;
t334 = t239 * t217;
t333 = t243 * t255;
t332 = t244 * t255;
t330 = t253 * t273;
t328 = t270 * t283;
t327 = t270 * t287;
t321 = qJD(2) * t284;
t314 = 0.2e1 * t339;
t252 = t255 ^ 2;
t215 = t252 * t217 + 0.1e1;
t313 = 0.2e1 * (-t252 * t345 + t255 * t334) / t215 ^ 2;
t311 = -0.2e1 * t227 * t303;
t306 = t205 * t311;
t305 = qJD(4) * t327 - t263;
t249 = -t271 * t286 - t282 * t327;
t250 = t271 * t282 - t286 * t327;
t301 = t249 * t285 - t250 * t281;
t233 = t249 * t281 + t250 * t285;
t300 = t281 * t286 - t282 * t285;
t254 = t269 * t287 - t308;
t298 = t254 * t265 + t266 * t330;
t295 = qJD(3) * t328 + qJD(4) * t271 - t264 * t287;
t262 = t269 * qJD(2);
t258 = t297 * qJD(3) - t287 * t307;
t238 = -t253 * qJD(3) + t261 * t287;
t234 = t300 * t255;
t223 = t305 * t282 + t295 * t286;
t222 = t295 * t282 - t305 * t286;
t212 = 0.1e1 / t215;
t211 = t245 * t342;
t210 = t298 * t245;
t207 = 0.1e1 / t209;
t203 = (t243 * t268 - t244 * t324) * t283 + t302 * t211;
t202 = -t302 * t210 + t243 * t254 + t244 * t273;
t200 = 0.2e1 * t298 * t343 + (-0.2e1 * t329 * t330 - t238 * t265 + (-t237 * t273 - t253 * t258 - t254 * t257) * t266) * t245;
t197 = -0.2e1 * t342 * t343 + (t296 * t319 + (0.2e1 * t310 * t329 + t262 * t265 + (-t257 * t268 + (t237 * t288 - t253 * t321) * t278) * t266) * t283) * t245;
t1 = [0, t197, t200, 0, 0, 0; 0 (t203 * t337 - t216 * t328) * t313 + ((t264 * t283 + t270 * t319) * t216 + (-t334 + t344) * t203 + (-t328 * t201 - (t197 * t253 + t211 * t237 + (t283 * t321 - t288 * t319) * t278 + (t211 * t297 + t268 * t283) * t206) * t332 - (t268 * t319 + t197 * t297 - t211 * t257 - t262 * t283 + (-t211 * t253 + t283 * t324) * t206) * t333) * t217) * t212 (t202 * t337 + t216 * t256) * t313 + (t202 * t344 - t240 * t216 + (t256 * t201 - t202 * t239 - (t200 * t253 - t210 * t237 + t258 + (-t210 * t297 + t254) * t206) * t332 - (t200 * t297 + t210 * t257 + t238 + (t210 * t253 - t273) * t206) * t333) * t217) * t212, 0, 0, 0; 0 (t225 * t301 - t233 * t336) * t314 + ((t233 * qJD(6) - t222 * t285 + t223 * t281) * t225 + t233 * t306 + (t301 * t205 + (t301 * qJD(6) + t222 * t281 + t223 * t285) * t303 - t233 * t204) * t226) * t207 (-t225 * t234 - t303 * t335) * t314 + (-t204 * t335 + (-t234 * t226 + t235 * t311) * t205 + (t300 * t225 + t299 * t336) * t239 + ((t341 * t225 + t340 * t336) * t285 + (t340 * t225 - t341 * t336) * t281) * t255) * t207 (t225 * t231 + t303 * t336) * t314 + (-t205 * t225 - t303 * t306 + (0.2e1 * t303 * t204 + t231 * t205) * t226) * t207, 0, -0.2e1 * t339 - 0.2e1 * (t204 * t226 * t207 - (-t207 * t338 - t226 * t339) * t303) * t303;];
JaD_rot  = t1;
