% Zeitableitung der rotatorischen Teilmatrix der analytischen Jacobi-Matrix für Segment Nr. 6 (0=Basis) von
% S6PRRRPR3
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
% Datum: 2019-02-26 20:11
% Revision: d75aae1ac561373cd3be920984c3df29a1c2ecc8 (2019-02-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JaD_rot = S6PRRRPR3_jacobiaD_rot_6_sym_varpar(qJ, qJD, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(11,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRRRPR3_jacobiaD_rot_6_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6PRRRPR3_jacobiaD_rot_6_sym_varpar: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6PRRRPR3_jacobiaD_rot_6_sym_varpar: pkin has to be [11x1] (double)');

%% Symbolic Calculation
% From jacobiaD_rot_6_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 20:11:44
% EndTime: 2019-02-26 20:11:45
% DurationCPUTime: 1.16s
% Computational Cost: add. (9025->112), mult. (13312->224), div. (822->12), fcn. (17117->13), ass. (0->114)
t264 = qJ(3) + qJ(4);
t262 = cos(t264);
t263 = qJD(3) + qJD(4);
t261 = sin(t264);
t267 = cos(pkin(6));
t269 = sin(qJ(2));
t325 = cos(pkin(11));
t295 = t325 * t269;
t265 = sin(pkin(11));
t271 = cos(qJ(2));
t306 = t265 * t271;
t281 = -t267 * t295 - t306;
t279 = t281 * t261;
t294 = t325 * t271;
t307 = t265 * t269;
t254 = t267 * t294 - t307;
t266 = sin(pkin(6));
t296 = t266 * t325;
t327 = t254 * qJD(2) - t263 * t296;
t217 = t327 * t262 + t263 * t279;
t278 = t281 * t262;
t240 = -t261 * t296 - t278;
t237 = t240 ^ 2;
t305 = t266 * t269;
t297 = t262 * t305;
t249 = t267 * t261 + t297;
t246 = 0.1e1 / t249 ^ 2;
t230 = t237 * t246 + 0.1e1;
t228 = 0.1e1 / t230;
t304 = t266 * t271;
t284 = qJD(2) * t304 + t263 * t267;
t298 = t261 * t305;
t235 = t284 * t262 - t263 * t298;
t245 = 0.1e1 / t249;
t314 = t240 * t246;
t200 = (-t217 * t245 + t235 * t314) * t228;
t231 = atan2(-t240, t249);
t226 = sin(t231);
t227 = cos(t231);
t288 = -t226 * t249 - t227 * t240;
t196 = t288 * t200 - t226 * t217 + t227 * t235;
t210 = -t226 * t240 + t227 * t249;
t207 = 0.1e1 / t210;
t208 = 0.1e1 / t210 ^ 2;
t331 = t196 * t207 * t208;
t256 = -t267 * t307 + t294;
t308 = t265 * t266;
t242 = t256 * t261 - t262 * t308;
t270 = cos(qJ(6));
t268 = sin(qJ(6));
t282 = -t267 * t306 - t295;
t312 = t282 * t268;
t287 = t242 * t270 + t312;
t330 = t287 * qJD(6);
t243 = t256 * t262 + t261 * t308;
t329 = 0.2e1 * t243 * t331;
t283 = -t245 * t254 + t304 * t314;
t328 = t262 * t283;
t315 = t235 * t245 * t246;
t326 = -0.2e1 * (t217 * t314 - t237 * t315) / t230 ^ 2;
t311 = t282 * t270;
t225 = t242 * t268 - t311;
t221 = 0.1e1 / t225;
t222 = 0.1e1 / t225 ^ 2;
t252 = t282 * qJD(2);
t290 = t263 * t308 + t252;
t309 = t262 * t263;
t218 = t256 * t309 + t290 * t261;
t253 = t256 * qJD(2);
t211 = t225 * qJD(6) - t218 * t270 + t253 * t268;
t220 = t287 ^ 2;
t215 = t220 * t222 + 0.1e1;
t319 = t222 * t287;
t212 = t218 * t268 + t253 * t270 + t330;
t322 = t212 * t221 * t222;
t324 = (-t211 * t319 - t220 * t322) / t215 ^ 2;
t323 = t208 * t243;
t310 = t261 * t263;
t219 = -t256 * t310 + t290 * t262;
t321 = t219 * t208;
t320 = t221 * t270;
t318 = t287 * t268;
t317 = t226 * t243;
t316 = t227 * t243;
t313 = t282 * t262;
t303 = qJD(2) * t269;
t238 = t243 ^ 2;
t206 = t238 * t208 + 0.1e1;
t302 = 0.2e1 * (-t238 * t331 + t243 * t321) / t206 ^ 2;
t301 = 0.2e1 * t324;
t293 = -0.2e1 * t287 * t322;
t292 = -0.2e1 * t240 * t315;
t289 = qJD(6) * t261 * t282 + t252;
t286 = -t222 * t318 + t320;
t239 = t262 * t296 - t279;
t248 = t267 * t262 - t298;
t285 = t239 * t245 + t248 * t314;
t280 = qJD(6) * t256 + t253 * t261 - t282 * t309;
t251 = t281 * qJD(2);
t234 = -t284 * t261 - t263 * t297;
t233 = t256 * t270 + t261 * t312;
t232 = t256 * t268 - t261 * t311;
t216 = t327 * t261 - t263 * t278;
t213 = 0.1e1 / t215;
t204 = 0.1e1 / t206;
t202 = t228 * t328;
t201 = t285 * t228;
t198 = (-t226 * t254 + t227 * t304) * t262 + t288 * t202;
t197 = t288 * t201 + t226 * t239 + t227 * t248;
t195 = t285 * t326 + (t248 * t292 + t216 * t245 + (t217 * t248 + t234 * t240 - t235 * t239) * t246) * t228;
t193 = t326 * t328 + (-t283 * t310 + (t292 * t304 - t245 * t251 + (t235 * t254 + (t217 * t271 - t240 * t303) * t266) * t246) * t262) * t228;
t192 = t286 * t243 * t301 + (-t286 * t219 + ((qJD(6) * t221 + t293) * t268 + (-t211 * t268 + (t212 + t330) * t270) * t222) * t243) * t213;
t191 = (t197 * t323 + t207 * t242) * t302 + (t197 * t329 - t218 * t207 + (t242 * t196 - t197 * t219 - (-t195 * t240 - t201 * t217 + t234 + (-t201 * t249 + t239) * t200) * t316 - (-t195 * t249 - t201 * t235 + t216 + (t201 * t240 - t248) * t200) * t317) * t208) * t204;
t1 = [0, t193, t195, t195, 0, 0; 0 (t198 * t323 - t207 * t313) * t302 + ((-t253 * t262 - t282 * t310) * t207 + (-t321 + t329) * t198 + (-t313 * t196 - (-t193 * t240 - t202 * t217 + (-t262 * t303 - t271 * t310) * t266 + (-t202 * t249 - t254 * t262) * t200) * t316 - (t254 * t310 - t193 * t249 - t202 * t235 - t251 * t262 + (t202 * t240 - t262 * t304) * t200) * t317) * t208) * t204, t191, t191, 0, 0; 0 (-t221 * t232 - t233 * t319) * t301 + (t233 * t293 + t289 * t221 * t268 + t280 * t320 + (t270 * t287 * t289 - t233 * t211 - t232 * t212 - t280 * t318) * t222) * t213, t192, t192, 0, -0.2e1 * t324 - 0.2e1 * (t211 * t222 * t213 - (-t213 * t322 - t222 * t324) * t287) * t287;];
JaD_rot  = t1;
