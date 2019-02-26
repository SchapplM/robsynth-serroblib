% Zeitableitung der rotatorischen Teilmatrix der analytischen Jacobi-Matrix für Segment Nr. 6 (0=Basis) von
% S6RRPPRR9
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
%   pkin=[a2,a3,a4,a5,a6,alpha2,d1,d2,d5,d6]';
%
% Output:
% JaD_rot [3x6]
%   Zeitableitung der rotatorischen Teilmatrix der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-02-26 21:33
% Revision: d75aae1ac561373cd3be920984c3df29a1c2ecc8 (2019-02-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JaD_rot = S6RRPPRR9_jacobiaD_rot_6_sym_varpar(qJ, qJD, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(10,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPPRR9_jacobiaD_rot_6_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRPPRR9_jacobiaD_rot_6_sym_varpar: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RRPPRR9_jacobiaD_rot_6_sym_varpar: pkin has to be [10x1] (double)');

%% Symbolic Calculation
% From jacobiaD_rot_6_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 21:33:13
% EndTime: 2019-02-26 21:33:15
% DurationCPUTime: 1.76s
% Computational Cost: add. (4522->148), mult. (13478->294), div. (726->12), fcn. (17045->13), ass. (0->125)
t260 = sin(qJ(1));
t339 = cos(pkin(6));
t340 = sin(qJ(2));
t288 = t339 * t340;
t252 = t260 * t288;
t263 = cos(qJ(2));
t341 = cos(qJ(1));
t289 = t339 * t341;
t298 = t341 * qJD(1);
t299 = qJD(2) * t340;
t257 = sin(pkin(6));
t303 = t257 * t341;
t349 = (qJD(2) * t289 + t298) * t263 - qJD(1) * t252 - t260 * t299 + qJD(5) * t303;
t282 = t341 * t288;
t318 = t260 * t263;
t246 = t282 + t318;
t259 = sin(qJ(5));
t262 = cos(qJ(5));
t280 = t246 * t262 + t259 * t303;
t228 = t280 ^ 2;
t302 = t257 * t340;
t276 = -t339 * t259 + t262 * t302;
t241 = 0.1e1 / t276 ^ 2;
t219 = t228 * t241 + 0.1e1;
t217 = 0.1e1 / t219;
t316 = qJD(1) * t260;
t301 = t257 * t316;
t347 = t349 * t262;
t202 = (qJD(5) * t246 + t301) * t259 - t347;
t244 = t259 * t302 + t339 * t262;
t319 = t257 * t263;
t300 = qJD(2) * t319;
t230 = t244 * qJD(5) - t262 * t300;
t240 = 0.1e1 / t276;
t325 = t280 * t241;
t285 = t202 * t240 - t230 * t325;
t186 = t285 * t217;
t220 = atan2(t280, -t276);
t215 = sin(t220);
t216 = cos(t220);
t286 = t215 * t276 + t216 * t280;
t181 = t286 * t186 - t215 * t202 + t216 * t230;
t198 = t215 * t280 - t216 * t276;
t196 = 0.1e1 / t198 ^ 2;
t348 = t181 * t196;
t314 = qJD(5) * t262;
t203 = t246 * t314 + t349 * t259 + t262 * t301;
t287 = t341 * t263 - t252;
t320 = t257 * t260;
t232 = t287 * t259 + t262 * t320;
t258 = sin(qJ(6));
t261 = cos(qJ(6));
t275 = t339 * t318 + t341 * t340;
t322 = t275 * t261;
t211 = t232 * t258 + t322;
t346 = 0.2e1 * t211;
t195 = 0.1e1 / t198;
t345 = t195 * t348;
t292 = t257 * t298;
t294 = -qJD(1) * t282 - t275 * qJD(2) - t263 * t316;
t205 = t232 * qJD(5) + t259 * t292 - t294 * t262;
t231 = t259 * t320 - t287 * t262;
t297 = 0.2e1 * t231 * t345;
t344 = -t196 * t205 + t297;
t343 = t230 * t241;
t245 = t260 * t340 - t263 * t289;
t304 = t280 * t319;
t281 = t240 * t245 + t241 * t304;
t342 = t262 * t281;
t212 = t232 * t261 - t258 * t275;
t208 = 0.1e1 / t212;
t209 = 0.1e1 / t212 ^ 2;
t227 = t231 ^ 2;
t194 = t196 * t227 + 0.1e1;
t333 = t196 * t231;
t338 = (t205 * t333 - t227 * t345) / t194 ^ 2;
t327 = t240 * t343;
t336 = (-t202 * t325 + t228 * t327) / t219 ^ 2;
t206 = -t231 * qJD(5) + t294 * t259 + t262 * t292;
t225 = t245 * qJD(1) - t287 * qJD(2);
t313 = qJD(6) * t211;
t191 = t206 * t261 + t225 * t258 - t313;
t335 = t191 * t208 * t209;
t207 = t211 ^ 2;
t201 = t207 * t209 + 0.1e1;
t199 = 0.1e1 / t201;
t332 = t199 * t209;
t190 = t212 * qJD(6) + t206 * t258 - t225 * t261;
t330 = t209 * t211;
t331 = 0.1e1 / t201 ^ 2 * (t190 * t330 - t207 * t335);
t329 = t215 * t231;
t328 = t216 * t231;
t326 = t280 * t240;
t324 = t280 * t244;
t323 = t275 * t259;
t321 = t275 * t262;
t315 = qJD(5) * t259;
t312 = 0.2e1 * t338;
t311 = -0.2e1 * t336;
t310 = 0.2e1 * t336;
t308 = -0.2e1 * t331;
t307 = t209 * t331;
t306 = t190 * t332;
t305 = t211 * t335;
t296 = t240 * t310;
t295 = 0.2e1 * t305;
t279 = -t246 * t259 + t262 * t303;
t214 = t245 * t258 + t261 * t279;
t213 = -t245 * t261 + t258 * t279;
t284 = -t208 * t258 + t261 * t330;
t283 = t240 * t279 + t241 * t324;
t278 = -t215 + (t216 * t326 + t215) * t217;
t277 = -qJD(6) * t323 + t294;
t273 = -t287 * qJD(6) + t225 * t259 - t275 * t314;
t229 = t276 * qJD(5) + t259 * t300;
t226 = t275 * qJD(1) + t246 * qJD(2);
t222 = -t287 * t258 - t259 * t322;
t192 = 0.1e1 / t194;
t189 = t217 * t342;
t188 = t283 * t217;
t183 = (-t215 * t245 - t216 * t319) * t262 + t286 * t189;
t182 = -t286 * t188 + t215 * t279 + t216 * t244;
t180 = t283 * t310 + (-0.2e1 * t324 * t327 + t203 * t240 + (t202 * t244 - t229 * t280 - t230 * t279) * t241) * t217;
t178 = t311 * t342 + (-t281 * t315 + (0.2e1 * t304 * t327 + t226 * t240 + (t230 * t245 + (-t202 * t263 - t280 * t299) * t257) * t241) * t262) * t217;
t1 = [-t231 * t296 + (t205 * t240 + t231 * t343) * t217, t178, 0, 0, t180, 0; -0.2e1 * t280 * t195 * t338 + ((-t246 * t315 - t259 * t301 + t347) * t195 - t280 * t348 - (t278 * t205 + ((-t186 * t217 * t326 + t311) * t215 + (-t280 * t296 - t186 + (t186 - t285) * t217) * t216) * t231) * t333) * t192 + (t344 * t192 + t333 * t312) * t278 * t231 (t183 * t333 - t195 * t321) * t312 + ((-t225 * t262 - t275 * t315) * t195 + t344 * t183 + (-t321 * t181 - (t178 * t280 - t189 * t202 + (t262 * t299 + t263 * t315) * t257 + (t189 * t276 - t245 * t262) * t186) * t328 - (t245 * t315 + t178 * t276 - t189 * t230 - t226 * t262 + (-t189 * t280 + t262 * t319) * t186) * t329) * t196) * t192, 0, 0 (t182 * t333 - t195 * t232) * t312 + (t182 * t297 + t206 * t195 + (-t232 * t181 - t182 * t205 - (t180 * t280 + t188 * t202 + t229 + (-t188 * t276 + t279) * t186) * t328 - (t180 * t276 + t188 * t230 - t203 + (t188 * t280 - t244) * t186) * t329) * t196) * t192, 0; 0.2e1 * (-t208 * t213 + t214 * t330) * t331 + ((t214 * qJD(6) - t203 * t258 - t226 * t261) * t208 + t214 * t295 + (-t213 * t191 - (-t213 * qJD(6) - t203 * t261 + t226 * t258) * t211 - t214 * t190) * t209) * t199 (t307 * t346 - t306) * t222 + (-t191 * t332 + t208 * t308) * (-t258 * t323 + t287 * t261) + ((t273 * t258 + t277 * t261) * t208 - (-t277 * t258 + t273 * t261) * t330 + t222 * t295) * t199, 0, 0, t284 * t231 * t308 + (t284 * t205 + ((-qJD(6) * t208 - 0.2e1 * t305) * t261 + (t190 * t261 + (t191 - t313) * t258) * t209) * t231) * t199, t308 + (t306 + (-t199 * t335 - t307) * t211) * t346;];
JaD_rot  = t1;
