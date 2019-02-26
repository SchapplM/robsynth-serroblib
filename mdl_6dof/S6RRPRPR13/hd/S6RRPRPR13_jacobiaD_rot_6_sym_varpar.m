% Zeitableitung der rotatorischen Teilmatrix der analytischen Jacobi-Matrix für Segment Nr. 6 (0=Basis) von
% S6RRPRPR13
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
%   pkin=[a2,a3,a4,a5,a6,alpha2,d1,d2,d4,d6,theta5]';
%
% Output:
% JaD_rot [3x6]
%   Zeitableitung der rotatorischen Teilmatrix der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-02-26 21:45
% Revision: d75aae1ac561373cd3be920984c3df29a1c2ecc8 (2019-02-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JaD_rot = S6RRPRPR13_jacobiaD_rot_6_sym_varpar(qJ, qJD, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(11,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPRPR13_jacobiaD_rot_6_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRPRPR13_jacobiaD_rot_6_sym_varpar: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6RRPRPR13_jacobiaD_rot_6_sym_varpar: pkin has to be [11x1] (double)');

%% Symbolic Calculation
% From jacobiaD_rot_6_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 21:44:50
% EndTime: 2019-02-26 21:44:52
% DurationCPUTime: 1.76s
% Computational Cost: add. (4943->150), mult. (13478->299), div. (726->12), fcn. (17045->13), ass. (0->127)
t261 = sin(qJ(2));
t262 = sin(qJ(1));
t264 = cos(qJ(2));
t265 = cos(qJ(1));
t343 = cos(pkin(6));
t295 = t265 * t343;
t247 = t262 * t261 - t264 * t295;
t260 = sin(qJ(4));
t263 = cos(qJ(4));
t259 = sin(pkin(6));
t323 = t259 * t265;
t284 = t247 * t263 + t260 * t323;
t230 = t284 ^ 2;
t324 = t259 * t264;
t281 = -t343 * t260 - t263 * t324;
t243 = 0.1e1 / t281 ^ 2;
t221 = t230 * t243 + 0.1e1;
t219 = 0.1e1 / t221;
t246 = -t260 * t324 + t343 * t263;
t317 = qJD(2) * t261;
t300 = t259 * t317;
t232 = t246 * qJD(4) - t263 * t300;
t242 = 0.1e1 / t281;
t329 = t284 * t243;
t280 = -t261 * t295 - t262 * t264;
t296 = t262 * t343;
t282 = t265 * t261 + t264 * t296;
t273 = t282 * qJD(1) - t280 * qJD(2);
t303 = t263 * t323;
t318 = qJD(1) * t262;
t344 = t247 * qJD(4) + t259 * t318;
t353 = qJD(4) * t303 - t344 * t260 + t273 * t263;
t288 = -t232 * t329 - t242 * t353;
t188 = t288 * t219;
t222 = atan2(t284, -t281);
t217 = sin(t222);
t218 = cos(t222);
t290 = t217 * t281 + t218 * t284;
t183 = t290 * t188 + t217 * t353 + t218 * t232;
t200 = t217 * t284 - t218 * t281;
t198 = 0.1e1 / t200 ^ 2;
t354 = t183 * t198;
t315 = qJD(4) * t260;
t297 = t259 * t315;
t205 = t273 * t260 + t344 * t263 + t265 * t297;
t325 = t259 * t263;
t234 = t282 * t260 + t262 * t325;
t291 = t261 * t296;
t320 = t265 * t264;
t249 = -t291 + t320;
t258 = pkin(11) + qJ(6);
t256 = sin(t258);
t257 = cos(t258);
t213 = t234 * t256 - t249 * t257;
t352 = 0.2e1 * t213;
t197 = 0.1e1 / t200;
t351 = t197 * t354;
t278 = t282 * t263;
t233 = t262 * t259 * t260 - t278;
t294 = 0.2e1 * t233 * t351;
t276 = (t343 * qJD(1) + qJD(2)) * t320 - qJD(2) * t291 - t261 * t318;
t301 = qJD(1) * t323;
t207 = t234 * qJD(4) + t260 * t301 - t276 * t263;
t335 = t207 * t198;
t350 = -t335 + t294;
t348 = t232 * t243;
t326 = t259 * t261;
t304 = t284 * t326;
t283 = t242 * t280 + t243 * t304;
t347 = t263 * t283;
t346 = -t249 * t260 * qJD(6) - t276;
t327 = t249 * t263;
t345 = qJD(4) * t327 - t282 * qJD(6);
t214 = t234 * t257 + t249 * t256;
t210 = 0.1e1 / t214;
t211 = 0.1e1 / t214 ^ 2;
t229 = t233 ^ 2;
t196 = t229 * t198 + 0.1e1;
t342 = (-t229 * t351 + t233 * t335) / t196 ^ 2;
t208 = qJD(4) * t278 + t276 * t260 - t262 * t297 + t263 * t301;
t227 = t280 * qJD(1) - t282 * qJD(2);
t191 = t214 * qJD(6) + t208 * t256 - t227 * t257;
t209 = t213 ^ 2;
t203 = t209 * t211 + 0.1e1;
t334 = t211 * t213;
t314 = qJD(6) * t213;
t192 = t208 * t257 + t227 * t256 - t314;
t338 = t192 * t210 * t211;
t341 = (t191 * t334 - t209 * t338) / t203 ^ 2;
t331 = t242 * t348;
t339 = (t230 * t331 + t329 * t353) / t221 ^ 2;
t337 = t198 * t233;
t201 = 0.1e1 / t203;
t336 = t201 * t211;
t333 = t217 * t233;
t332 = t218 * t233;
t330 = t284 * t242;
t328 = t284 * t246;
t322 = t260 * t256;
t321 = t260 * t257;
t316 = qJD(2) * t264;
t312 = 0.2e1 * t342;
t311 = -0.2e1 * t341;
t310 = -0.2e1 * t339;
t309 = 0.2e1 * t339;
t307 = t211 * t341;
t306 = t191 * t336;
t305 = t213 * t338;
t293 = t242 * t309;
t292 = 0.2e1 * t305;
t285 = -t247 * t260 + t303;
t216 = t256 * t280 + t257 * t285;
t215 = t256 * t285 - t257 * t280;
t287 = -t256 * t210 + t257 * t334;
t286 = t242 * t285 + t243 * t328;
t279 = -t217 + (t218 * t330 + t217) * t219;
t231 = t281 * qJD(4) + t260 * t300;
t228 = -qJD(1) * t291 - t262 * t317 + (qJD(2) * t343 + qJD(1)) * t320;
t224 = t249 * t321 - t282 * t256;
t194 = 0.1e1 / t196;
t193 = t219 * t347;
t190 = t286 * t219;
t185 = (-t217 * t280 - t218 * t326) * t263 + t290 * t193;
t184 = -t290 * t190 + t217 * t285 + t218 * t246;
t182 = t286 * t309 + (-0.2e1 * t328 * t331 + t205 * t242 + (-t231 * t284 - t232 * t285 - t246 * t353) * t243) * t219;
t180 = t310 * t347 + (-t283 * t315 + (0.2e1 * t304 * t331 - t228 * t242 + (t232 * t280 + (t261 * t353 + t284 * t316) * t259) * t243) * t263) * t219;
t1 = [-t233 * t293 + (t207 * t242 + t233 * t348) * t219, t180, 0, t182, 0, 0; -0.2e1 * t284 * t197 * t342 + (t353 * t197 - t284 * t354 - (t279 * t207 + ((-t188 * t219 * t330 + t310) * t217 + (-t284 * t293 - t188 + (t188 - t288) * t219) * t218) * t233) * t337) * t194 + (t350 * t194 + t337 * t312) * t279 * t233 (t185 * t337 + t197 * t327) * t312 + ((-t227 * t263 + t249 * t315) * t197 + t350 * t185 + (t327 * t183 - (t180 * t284 + t193 * t353 + (t261 * t315 - t263 * t316) * t259 + (t193 * t281 - t263 * t280) * t188) * t332 - (t280 * t315 + t180 * t281 - t193 * t232 + t228 * t263 + (-t193 * t284 + t261 * t325) * t188) * t333) * t198) * t194, 0 (t184 * t337 - t197 * t234) * t312 + (t184 * t294 + t208 * t197 + (-t234 * t183 - t184 * t207 - (t182 * t284 - t190 * t353 + t231 + (-t190 * t281 + t285) * t188) * t332 - (t182 * t281 + t190 * t232 - t205 + (t190 * t284 - t246) * t188) * t333) * t198) * t194, 0, 0; 0.2e1 * (-t210 * t215 + t216 * t334) * t341 + ((t216 * qJD(6) - t205 * t256 + t228 * t257) * t210 + t216 * t292 + (-t215 * t192 - (-t215 * qJD(6) - t205 * t257 - t228 * t256) * t213 - t216 * t191) * t211) * t201 (t307 * t352 - t306) * t224 + (-t192 * t336 + t210 * t311) * (t249 * t322 + t282 * t257) + (t224 * t292 + (t322 * t210 - t321 * t334) * t227 + (-t346 * t210 - t345 * t334) * t257 + (t345 * t210 - t346 * t334) * t256) * t201, 0, t287 * t233 * t311 + (t287 * t207 + ((-qJD(6) * t210 - 0.2e1 * t305) * t257 + (t191 * t257 + (t192 - t314) * t256) * t211) * t233) * t201, 0, t311 + (t306 + (-t201 * t338 - t307) * t213) * t352;];
JaD_rot  = t1;
