% Zeitableitung der rotatorischen Teilmatrix der analytischen Jacobi-Matrix für Segment Nr. 5 (0=Basis) von
% S6RRRRPR11
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
% Datum: 2019-02-26 22:36
% Revision: d75aae1ac561373cd3be920984c3df29a1c2ecc8 (2019-02-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JaD_rot = S6RRRRPR11_jacobiaD_rot_5_sym_varpar(qJ, qJD, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(12,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRRPR11_jacobiaD_rot_5_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRRRPR11_jacobiaD_rot_5_sym_varpar: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [12 1]), ...
  'S6RRRRPR11_jacobiaD_rot_5_sym_varpar: pkin has to be [12x1] (double)');

%% Symbolic Calculation
% From jacobiaD_rot_5_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 22:36:22
% EndTime: 2019-02-26 22:36:23
% DurationCPUTime: 1.40s
% Computational Cost: add. (4943->149), mult. (13478->297), div. (726->12), fcn. (17045->13), ass. (0->127)
t263 = cos(pkin(6));
t265 = sin(qJ(2));
t336 = sin(qJ(1));
t298 = t336 * t265;
t288 = t263 * t298;
t292 = t336 * qJD(2);
t267 = cos(qJ(2));
t268 = cos(qJ(1));
t314 = t268 * t267;
t262 = sin(pkin(6));
t316 = t262 * t268;
t341 = -qJD(1) * t288 - t265 * t292 + (qJD(2) * t263 + qJD(1)) * t314 - qJD(3) * t316;
t264 = sin(qJ(3));
t266 = cos(qJ(3));
t297 = t336 * t267;
t315 = t268 * t265;
t280 = -t263 * t315 - t297;
t232 = -t264 * t280 + t266 * t316;
t318 = t262 * t265;
t243 = -t263 * t266 + t264 * t318;
t221 = atan2(-t232, t243);
t216 = sin(t221);
t217 = cos(t221);
t199 = -t216 * t232 + t217 * t243;
t197 = 0.1e1 / t199 ^ 2;
t248 = -t288 + t314;
t299 = t262 * t336;
t279 = -t248 * t264 + t266 * t299;
t229 = t279 ^ 2;
t195 = t229 * t197 + 0.1e1;
t278 = -t263 * t297 - t315;
t225 = t280 * qJD(1) + t278 * qJD(2);
t238 = t248 * t266 + t264 * t299;
t296 = qJD(1) * t316;
t203 = t238 * qJD(3) + t225 * t264 - t266 * t296;
t329 = t203 * t197;
t228 = t232 ^ 2;
t241 = 0.1e1 / t243 ^ 2;
t220 = t228 * t241 + 0.1e1;
t218 = 0.1e1 / t220;
t293 = t336 * qJD(1);
t287 = t262 * t293;
t311 = qJD(3) * t266;
t205 = t264 * t341 - t266 * t287 - t280 * t311;
t244 = t263 * t264 + t266 * t318;
t312 = qJD(2) * t267;
t295 = t262 * t312;
t230 = t244 * qJD(3) + t264 * t295;
t240 = 0.1e1 / t243;
t323 = t232 * t241;
t284 = -t205 * t240 + t230 * t323;
t187 = t284 * t218;
t285 = -t216 * t243 - t217 * t232;
t182 = t285 * t187 - t216 * t205 + t217 * t230;
t196 = 0.1e1 / t199;
t198 = t196 * t197;
t334 = t182 * t198;
t309 = 0.2e1 * (-t229 * t334 - t279 * t329) / t195 ^ 2;
t340 = t230 * t241;
t300 = t263 * t314;
t245 = -t298 + t300;
t317 = t262 * t267;
t281 = -t240 * t245 + t317 * t323;
t339 = t264 * t281;
t206 = (qJD(3) * t280 + t287) * t264 + t341 * t266;
t261 = qJ(4) + pkin(12);
t259 = sin(t261);
t260 = cos(t261);
t215 = t238 * t260 - t259 * t278;
t209 = 0.1e1 / t215;
t210 = 0.1e1 / t215 ^ 2;
t338 = -0.2e1 * t232;
t337 = -0.2e1 * t279;
t204 = t279 * qJD(3) + t225 * t266 + t264 * t296;
t224 = -qJD(1) * t300 - t268 * t312 + (t263 * t292 + t293) * t265;
t190 = t215 * qJD(4) + t204 * t259 + t224 * t260;
t214 = t238 * t259 + t260 * t278;
t208 = t214 ^ 2;
t202 = t208 * t210 + 0.1e1;
t328 = t210 * t214;
t310 = qJD(4) * t214;
t191 = t204 * t260 - t224 * t259 - t310;
t331 = t191 * t209 * t210;
t333 = (t190 * t328 - t208 * t331) / t202 ^ 2;
t325 = t240 * t340;
t332 = (t205 * t323 - t228 * t325) / t220 ^ 2;
t330 = t197 * t279;
t327 = t216 * t279;
t326 = t217 * t279;
t324 = t232 * t240;
t322 = t278 * t264;
t321 = t278 * t266;
t320 = t259 * t209;
t319 = t260 * t214;
t313 = qJD(2) * t265;
t308 = -0.2e1 * t333;
t307 = 0.2e1 * t333;
t306 = -0.2e1 * t332;
t305 = t198 * t337;
t304 = t240 * t332;
t303 = t197 * t327;
t302 = t197 * t326;
t301 = t214 * t331;
t291 = 0.2e1 * t301;
t290 = t325 * t338;
t234 = -t264 * t316 - t266 * t280;
t286 = -qJD(4) * t321 + t225;
t213 = -t234 * t260 + t245 * t259;
t212 = -t234 * t259 - t245 * t260;
t283 = t210 * t319 - t320;
t282 = -t234 * t240 + t244 * t323;
t276 = -t216 + (t217 * t324 + t216) * t218;
t275 = -qJD(3) * t322 + qJD(4) * t248 + t224 * t266;
t231 = -t243 * qJD(3) + t266 * t295;
t226 = t278 * qJD(1) + t280 * qJD(2);
t223 = t248 * t259 + t260 * t321;
t222 = -t248 * t260 + t259 * t321;
t200 = 0.1e1 / t202;
t193 = 0.1e1 / t195;
t192 = t218 * t339;
t189 = t282 * t218;
t186 = t276 * t279;
t184 = (-t216 * t245 + t217 * t317) * t264 + t285 * t192;
t183 = t285 * t189 - t216 * t234 + t217 * t244;
t181 = t282 * t306 + (t244 * t290 - t206 * t240 + (t205 * t244 + t230 * t234 + t231 * t232) * t241) * t218;
t179 = t306 * t339 + (t281 * t311 + (t290 * t317 - t226 * t240 + (t230 * t245 + (t205 * t267 - t232 * t313) * t262) * t241) * t264) * t218;
t1 = [t304 * t337 + (-t203 * t240 - t279 * t340) * t218, t179, t181, 0, 0, 0; t232 * t196 * t309 + (-t205 * t196 + (t182 * t232 + t186 * t203) * t197) * t193 - (-t186 * t197 * t309 + (-0.2e1 * t186 * t334 + (-t187 * t218 * t324 + t306) * t303 + (t304 * t338 - t187 + (t187 - t284) * t218) * t302 - t276 * t329) * t193) * t279 (-t184 * t330 - t196 * t322) * t309 + (-t184 * t329 + (t224 * t264 + t278 * t311) * t196 + (t184 * t305 - t197 * t322) * t182 + (-t179 * t232 - t192 * t205 + (-t264 * t313 + t267 * t311) * t262 + (-t192 * t243 - t245 * t264) * t187) * t302 + (-t245 * t311 - t179 * t243 - t192 * t230 - t226 * t264 + (t192 * t232 - t264 * t317) * t187) * t303) * t193 (-t183 * t330 - t196 * t238) * t309 + (t183 * t182 * t305 + t204 * t196 + (-t238 * t182 - t183 * t203 + (-t181 * t232 - t189 * t205 + t231 + (-t189 * t243 - t234) * t187) * t326 + (-t181 * t243 - t189 * t230 - t206 + (t189 * t232 - t244) * t187) * t327) * t197) * t193, 0, 0, 0; (-t209 * t212 + t213 * t328) * t307 + ((t213 * qJD(4) - t206 * t259 - t226 * t260) * t209 + t213 * t291 + (-t212 * t191 - (-t212 * qJD(4) - t206 * t260 + t226 * t259) * t214 - t213 * t190) * t210) * t200 (-t209 * t222 + t223 * t328) * t307 + (t223 * t291 - t286 * t209 * t260 + t275 * t320 + (-t286 * t214 * t259 - t223 * t190 - t222 * t191 - t275 * t319) * t210) * t200, -t283 * t279 * t308 + (t283 * t203 - ((-qJD(4) * t209 - 0.2e1 * t301) * t260 + (t190 * t260 + (t191 - t310) * t259) * t210) * t279) * t200, t308 + 0.2e1 * (t190 * t210 * t200 + (-t200 * t331 - t210 * t333) * t214) * t214, 0, 0;];
JaD_rot  = t1;
