% Zeitableitung der rotatorischen Teilmatrix der analytischen Jacobi-Matrix für Segment Nr. 6 (0=Basis) von
% S6RRRPRP11
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
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d1,d2,d3,d5]';
%
% Output:
% JaD_rot [3x6]
%   Zeitableitung der rotatorischen Teilmatrix der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-02-26 22:15
% Revision: d75aae1ac561373cd3be920984c3df29a1c2ecc8 (2019-02-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JaD_rot = S6RRRPRP11_jacobiaD_rot_6_sym_varpar(qJ, qJD, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(10,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRPRP11_jacobiaD_rot_6_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRRPRP11_jacobiaD_rot_6_sym_varpar: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RRRPRP11_jacobiaD_rot_6_sym_varpar: pkin has to be [10x1] (double)');

%% Symbolic Calculation
% From jacobiaD_rot_6_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 22:14:59
% EndTime: 2019-02-26 22:15:00
% DurationCPUTime: 1.33s
% Computational Cost: add. (4522->150), mult. (13478->299), div. (726->12), fcn. (17045->13), ass. (0->128)
t265 = sin(qJ(1));
t261 = cos(pkin(6));
t286 = qJD(2) * t261 + qJD(1);
t264 = sin(qJ(2));
t312 = t265 * t264;
t295 = t261 * t312;
t306 = qJD(2) * t264;
t268 = cos(qJ(2));
t269 = cos(qJ(1));
t308 = t269 * t268;
t260 = sin(pkin(6));
t315 = t260 * t269;
t340 = -qJD(1) * t295 - qJD(3) * t315 - t265 * t306 + t286 * t308;
t250 = -t295 + t308;
t263 = sin(qJ(3));
t267 = cos(qJ(3));
t317 = t260 * t267;
t238 = t250 * t263 - t265 * t317;
t309 = t269 * t264;
t311 = t265 * t268;
t249 = t261 * t311 + t309;
t262 = sin(qJ(5));
t266 = cos(qJ(5));
t282 = t238 * t266 - t249 * t262;
t339 = t282 * qJD(5);
t248 = t261 * t309 + t311;
t313 = t263 * t269;
t234 = t248 * t267 - t260 * t313;
t246 = t261 * t263 + t264 * t317;
t221 = atan2(-t234, t246);
t216 = sin(t221);
t217 = cos(t221);
t199 = -t216 * t234 + t217 * t246;
t197 = 0.1e1 / t199 ^ 2;
t318 = t260 * t263;
t239 = t250 * t267 + t265 * t318;
t230 = t239 ^ 2;
t195 = t230 * t197 + 0.1e1;
t225 = -t248 * qJD(1) - t249 * qJD(2);
t303 = qJD(3) * t267;
t304 = qJD(3) * t263;
t204 = t225 * t267 - t250 * t304 + (qJD(1) * t313 + t265 * t303) * t260;
t328 = t204 * t197;
t229 = t234 ^ 2;
t243 = 0.1e1 / t246 ^ 2;
t220 = t229 * t243 + 0.1e1;
t218 = 0.1e1 / t220;
t287 = t340 * t267;
t307 = qJD(1) * t260;
t293 = t265 * t307;
t206 = -t248 * t304 + t263 * t293 + t287;
t245 = t261 * t267 - t264 * t318;
t305 = qJD(2) * t268;
t291 = t260 * t305;
t232 = t245 * qJD(3) + t267 * t291;
t242 = 0.1e1 / t246;
t322 = t234 * t243;
t281 = -t206 * t242 + t232 * t322;
t187 = t281 * t218;
t284 = -t216 * t246 - t217 * t234;
t182 = t284 * t187 - t216 * t206 + t217 * t232;
t196 = 0.1e1 / t199;
t198 = t196 * t197;
t333 = t182 * t198;
t302 = 0.2e1 * (-t230 * t333 + t239 * t328) / t195 ^ 2;
t338 = t232 * t243;
t294 = t261 * t308;
t247 = t294 - t312;
t316 = t260 * t268;
t278 = -t242 * t247 + t316 * t322;
t337 = t267 * t278;
t320 = t249 * t266;
t215 = t238 * t262 + t320;
t209 = 0.1e1 / t215;
t210 = 0.1e1 / t215 ^ 2;
t336 = -0.2e1 * t234;
t335 = 0.2e1 * t239;
t292 = t267 * t307;
t203 = t239 * qJD(3) + t225 * t263 - t269 * t292;
t224 = -qJD(1) * t294 - t269 * t305 + t286 * t312;
t191 = t215 * qJD(5) - t203 * t266 - t224 * t262;
t208 = t282 ^ 2;
t202 = t208 * t210 + 0.1e1;
t327 = t210 * t282;
t192 = t203 * t262 - t224 * t266 + t339;
t330 = t192 * t209 * t210;
t332 = (-t191 * t327 - t208 * t330) / t202 ^ 2;
t324 = t242 * t338;
t331 = (t206 * t322 - t229 * t324) / t220 ^ 2;
t329 = t197 * t239;
t326 = t216 * t239;
t325 = t217 * t239;
t323 = t234 * t242;
t321 = t249 * t263;
t319 = t249 * t267;
t314 = t262 * t282;
t310 = t266 * t209;
t301 = 0.2e1 * t332;
t300 = -0.2e1 * t331;
t299 = t198 * t335;
t298 = t242 * t331;
t297 = t197 * t326;
t296 = t197 * t325;
t289 = -0.2e1 * t282 * t330;
t288 = t324 * t336;
t285 = -qJD(5) * t321 + t225;
t233 = t248 * t263 + t267 * t315;
t283 = -t233 * t266 - t247 * t262;
t213 = -t233 * t262 + t247 * t266;
t280 = -t210 * t314 + t310;
t279 = t233 * t242 + t245 * t322;
t277 = -t216 + (t217 * t323 + t216) * t218;
t205 = t248 * t303 + t340 * t263 - t265 * t292;
t276 = qJD(5) * t250 - t224 * t263 + t249 * t303;
t231 = -t246 * qJD(3) - t263 * t291;
t226 = -t249 * qJD(1) - t248 * qJD(2);
t223 = t250 * t266 - t262 * t321;
t222 = t250 * t262 + t263 * t320;
t200 = 0.1e1 / t202;
t193 = 0.1e1 / t195;
t190 = t218 * t337;
t189 = t279 * t218;
t186 = t277 * t239;
t184 = (-t216 * t247 + t217 * t316) * t267 + t284 * t190;
t183 = t284 * t189 + t216 * t233 + t217 * t245;
t181 = t279 * t300 + (t245 * t288 + t205 * t242 + (t206 * t245 + t231 * t234 - t232 * t233) * t243) * t218;
t179 = t300 * t337 + (-t278 * t304 + (t288 * t316 - t226 * t242 + (t232 * t247 + (t206 * t268 - t234 * t306) * t260) * t243) * t267) * t218;
t1 = [t298 * t335 + (-t204 * t242 + t239 * t338) * t218, t179, t181, 0, 0, 0; t234 * t196 * t302 + (((qJD(3) * t248 - t293) * t263 - t287) * t196 + (t182 * t234 - t186 * t204) * t197) * t193 + (t186 * t197 * t302 + (0.2e1 * t186 * t333 - (-t187 * t218 * t323 + t300) * t297 - (t298 * t336 - t187 + (t187 - t281) * t218) * t296 - t277 * t328) * t193) * t239 (t184 * t329 + t196 * t319) * t302 + (-t184 * t328 + (t224 * t267 + t249 * t304) * t196 + (t184 * t299 + t197 * t319) * t182 - (-t179 * t234 - t190 * t206 + (-t267 * t306 - t268 * t304) * t260 + (-t190 * t246 - t247 * t267) * t187) * t296 - (t247 * t304 - t179 * t246 - t190 * t232 - t226 * t267 + (t190 * t234 - t267 * t316) * t187) * t297) * t193 (t183 * t329 + t196 * t238) * t302 + (t183 * t182 * t299 - t203 * t196 + (t238 * t182 - t183 * t204 - (-t181 * t234 - t189 * t206 + t231 + (-t189 * t246 + t233) * t187) * t325 - (-t181 * t246 - t189 * t232 + t205 + (t189 * t234 - t245) * t187) * t326) * t197) * t193, 0, 0, 0; (t209 * t283 - t213 * t327) * t301 + ((t213 * qJD(5) + t205 * t266 + t226 * t262) * t209 + t213 * t289 + (t283 * t192 + (t283 * qJD(5) - t205 * t262 + t226 * t266) * t282 - t213 * t191) * t210) * t200 (-t209 * t222 - t223 * t327) * t301 + (t223 * t289 + t285 * t209 * t262 + t276 * t310 + (t266 * t282 * t285 - t223 * t191 - t222 * t192 - t276 * t314) * t210) * t200, t280 * t239 * t301 + (-t280 * t204 + ((qJD(5) * t209 + t289) * t262 + (-t191 * t262 + (t192 + t339) * t266) * t210) * t239) * t200, 0, -0.2e1 * t332 - 0.2e1 * (t191 * t210 * t200 - (-t200 * t330 - t210 * t332) * t282) * t282, 0;];
JaD_rot  = t1;
