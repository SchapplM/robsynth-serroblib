% Zeitableitung der rotatorischen Teilmatrix der analytischen Jacobi-Matrix für Segment Nr. 6 (0=Basis) von
% S6RRRPPR10
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
%   pkin=[a2,a3,a4,a5,a6,alpha2,d1,d2,d3,d6,theta5]';
%
% Output:
% JaD_rot [3x6]
%   Zeitableitung der rotatorischen Teilmatrix der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-02-26 22:08
% Revision: d75aae1ac561373cd3be920984c3df29a1c2ecc8 (2019-02-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JaD_rot = S6RRRPPR10_jacobiaD_rot_6_sym_varpar(qJ, qJD, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(11,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRPPR10_jacobiaD_rot_6_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRRPPR10_jacobiaD_rot_6_sym_varpar: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6RRRPPR10_jacobiaD_rot_6_sym_varpar: pkin has to be [11x1] (double)');

%% Symbolic Calculation
% From jacobiaD_rot_6_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 22:08:52
% EndTime: 2019-02-26 22:08:53
% DurationCPUTime: 1.34s
% Computational Cost: add. (4943->151), mult. (13478->299), div. (726->12), fcn. (17045->13), ass. (0->128)
t266 = sin(qJ(1));
t263 = cos(pkin(6));
t286 = qJD(2) * t263 + qJD(1);
t265 = sin(qJ(2));
t311 = t266 * t265;
t295 = t263 * t311;
t306 = qJD(2) * t265;
t268 = cos(qJ(2));
t269 = cos(qJ(1));
t308 = t269 * t268;
t262 = sin(pkin(6));
t313 = t262 * t269;
t339 = -qJD(1) * t295 - qJD(3) * t313 - t266 * t306 + t286 * t308;
t249 = -t295 + t308;
t264 = sin(qJ(3));
t267 = cos(qJ(3));
t315 = t262 * t267;
t237 = t249 * t264 - t266 * t315;
t309 = t269 * t265;
t310 = t266 * t268;
t248 = t263 * t310 + t309;
t261 = pkin(11) + qJ(6);
t259 = sin(t261);
t260 = cos(t261);
t282 = t237 * t260 - t248 * t259;
t338 = t282 * qJD(6);
t247 = t263 * t309 + t310;
t312 = t264 * t269;
t233 = t247 * t267 - t262 * t312;
t245 = t263 * t264 + t265 * t315;
t220 = atan2(-t233, t245);
t215 = sin(t220);
t216 = cos(t220);
t198 = -t215 * t233 + t216 * t245;
t196 = 0.1e1 / t198 ^ 2;
t316 = t262 * t264;
t238 = t249 * t267 + t266 * t316;
t229 = t238 ^ 2;
t194 = t229 * t196 + 0.1e1;
t224 = -t247 * qJD(1) - t248 * qJD(2);
t303 = qJD(3) * t267;
t304 = qJD(3) * t264;
t203 = t224 * t267 - t249 * t304 + (qJD(1) * t312 + t266 * t303) * t262;
t327 = t203 * t196;
t228 = t233 ^ 2;
t242 = 0.1e1 / t245 ^ 2;
t219 = t228 * t242 + 0.1e1;
t217 = 0.1e1 / t219;
t287 = t339 * t267;
t307 = qJD(1) * t262;
t293 = t266 * t307;
t205 = -t247 * t304 + t264 * t293 + t287;
t244 = t263 * t267 - t265 * t316;
t305 = qJD(2) * t268;
t291 = t262 * t305;
t231 = t244 * qJD(3) + t267 * t291;
t241 = 0.1e1 / t245;
t321 = t233 * t242;
t281 = -t205 * t241 + t231 * t321;
t186 = t281 * t217;
t284 = -t215 * t245 - t216 * t233;
t181 = t284 * t186 - t215 * t205 + t216 * t231;
t195 = 0.1e1 / t198;
t197 = t195 * t196;
t332 = t181 * t197;
t302 = 0.2e1 * (-t229 * t332 + t238 * t327) / t194 ^ 2;
t337 = t231 * t242;
t294 = t263 * t308;
t246 = t294 - t311;
t314 = t262 * t268;
t278 = -t241 * t246 + t314 * t321;
t336 = t267 * t278;
t214 = t237 * t259 + t248 * t260;
t208 = 0.1e1 / t214;
t209 = 0.1e1 / t214 ^ 2;
t335 = -0.2e1 * t233;
t334 = 0.2e1 * t238;
t292 = t267 * t307;
t202 = t238 * qJD(3) + t224 * t264 - t269 * t292;
t223 = -qJD(1) * t294 - t269 * t305 + t286 * t311;
t189 = t214 * qJD(6) - t202 * t260 - t223 * t259;
t207 = t282 ^ 2;
t201 = t207 * t209 + 0.1e1;
t326 = t209 * t282;
t190 = t202 * t259 - t223 * t260 + t338;
t329 = t190 * t208 * t209;
t331 = (-t189 * t326 - t207 * t329) / t201 ^ 2;
t323 = t241 * t337;
t330 = (t205 * t321 - t228 * t323) / t219 ^ 2;
t328 = t196 * t238;
t325 = t215 * t238;
t324 = t216 * t238;
t322 = t233 * t241;
t320 = t248 * t264;
t319 = t248 * t267;
t318 = t259 * t282;
t317 = t260 * t208;
t301 = 0.2e1 * t331;
t300 = -0.2e1 * t330;
t299 = t197 * t334;
t298 = t241 * t330;
t297 = t196 * t325;
t296 = t196 * t324;
t289 = -0.2e1 * t282 * t329;
t288 = t323 * t335;
t285 = -qJD(6) * t320 + t224;
t232 = t247 * t264 + t267 * t313;
t283 = -t232 * t260 - t246 * t259;
t212 = -t232 * t259 + t246 * t260;
t280 = -t209 * t318 + t317;
t279 = t232 * t241 + t244 * t321;
t277 = -t215 + (t216 * t322 + t215) * t217;
t204 = t247 * t303 + t339 * t264 - t266 * t292;
t276 = qJD(6) * t249 - t223 * t264 + t248 * t303;
t230 = -t245 * qJD(3) - t264 * t291;
t225 = -t248 * qJD(1) - t247 * qJD(2);
t222 = t249 * t260 - t259 * t320;
t221 = t249 * t259 + t260 * t320;
t199 = 0.1e1 / t201;
t192 = 0.1e1 / t194;
t191 = t217 * t336;
t188 = t279 * t217;
t185 = t277 * t238;
t183 = (-t215 * t246 + t216 * t314) * t267 + t284 * t191;
t182 = t284 * t188 + t215 * t232 + t216 * t244;
t180 = t279 * t300 + (t244 * t288 + t204 * t241 + (t205 * t244 + t230 * t233 - t231 * t232) * t242) * t217;
t178 = t300 * t336 + (-t278 * t304 + (t288 * t314 - t225 * t241 + (t231 * t246 + (t205 * t268 - t233 * t306) * t262) * t242) * t267) * t217;
t1 = [t298 * t334 + (-t203 * t241 + t238 * t337) * t217, t178, t180, 0, 0, 0; t233 * t195 * t302 + (((qJD(3) * t247 - t293) * t264 - t287) * t195 + (t181 * t233 - t185 * t203) * t196) * t192 + (t185 * t196 * t302 + (0.2e1 * t185 * t332 - (-t186 * t217 * t322 + t300) * t297 - (t298 * t335 - t186 + (t186 - t281) * t217) * t296 - t277 * t327) * t192) * t238 (t183 * t328 + t195 * t319) * t302 + (-t183 * t327 + (t223 * t267 + t248 * t304) * t195 + (t183 * t299 + t196 * t319) * t181 - (-t178 * t233 - t191 * t205 + (-t267 * t306 - t268 * t304) * t262 + (-t191 * t245 - t246 * t267) * t186) * t296 - (t246 * t304 - t178 * t245 - t191 * t231 - t225 * t267 + (t191 * t233 - t267 * t314) * t186) * t297) * t192 (t182 * t328 + t195 * t237) * t302 + (t182 * t181 * t299 - t202 * t195 + (t237 * t181 - t182 * t203 - (-t180 * t233 - t188 * t205 + t230 + (-t188 * t245 + t232) * t186) * t324 - (-t180 * t245 - t188 * t231 + t204 + (t188 * t233 - t244) * t186) * t325) * t196) * t192, 0, 0, 0; (t208 * t283 - t212 * t326) * t301 + ((t212 * qJD(6) + t204 * t260 + t225 * t259) * t208 + t212 * t289 + (t283 * t190 + (t283 * qJD(6) - t204 * t259 + t225 * t260) * t282 - t212 * t189) * t209) * t199 (-t208 * t221 - t222 * t326) * t301 + (t222 * t289 + t285 * t208 * t259 + t276 * t317 + (t260 * t282 * t285 - t222 * t189 - t221 * t190 - t276 * t318) * t209) * t199, t280 * t238 * t301 + (-t280 * t203 + ((qJD(6) * t208 + t289) * t259 + (-t189 * t259 + (t190 + t338) * t260) * t209) * t238) * t199, 0, 0, -0.2e1 * t331 - 0.2e1 * (t189 * t209 * t199 - (-t199 * t329 - t209 * t331) * t282) * t282;];
JaD_rot  = t1;
