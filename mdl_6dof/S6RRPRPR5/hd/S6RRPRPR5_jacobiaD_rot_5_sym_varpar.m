% Zeitableitung der rotatorischen Teilmatrix der analytischen Jacobi-Matrix für Segment Nr. 5 (0=Basis) von
% S6RRPRPR5
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
%   pkin=[a2,a3,a4,a5,a6,alpha2,d1,d2,d4,d6,theta3,theta5]';
%
% Output:
% JaD_rot [3x6]
%   Zeitableitung der rotatorischen Teilmatrix der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-02-26 21:40
% Revision: d75aae1ac561373cd3be920984c3df29a1c2ecc8 (2019-02-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JaD_rot = S6RRPRPR5_jacobiaD_rot_5_sym_varpar(qJ, qJD, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(12,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPRPR5_jacobiaD_rot_5_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRPRPR5_jacobiaD_rot_5_sym_varpar: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [12 1]), ...
  'S6RRPRPR5_jacobiaD_rot_5_sym_varpar: pkin has to be [12x1] (double)');

%% Symbolic Calculation
% From jacobiaD_rot_5_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 21:40:09
% EndTime: 2019-02-26 21:40:10
% DurationCPUTime: 1.63s
% Computational Cost: add. (7685->143), mult. (22101->293), div. (702->12), fcn. (28611->15), ass. (0->125)
t276 = sin(pkin(11));
t278 = cos(pkin(11));
t281 = sin(qJ(2));
t284 = cos(qJ(2));
t265 = t281 * t276 - t284 * t278;
t279 = cos(pkin(6));
t262 = t265 * t279;
t257 = qJD(2) * t262;
t302 = t284 * t276 + t281 * t278;
t264 = t302 * qJD(2);
t282 = sin(qJ(1));
t285 = cos(qJ(1));
t345 = sin(pkin(6));
t312 = t285 * t345;
t263 = t302 * t279;
t327 = t282 * t263;
t350 = -t282 * t264 - qJD(1) * t327 + (-qJD(1) * t265 - t257) * t285 - qJD(4) * t312;
t280 = sin(qJ(4));
t283 = cos(qJ(4));
t303 = -t285 * t265 - t327;
t314 = t282 * t345;
t242 = t280 * t314 + t283 * t303;
t246 = t282 * t262 - t285 * t302;
t275 = sin(pkin(12));
t277 = cos(pkin(12));
t217 = t242 * t275 + t246 * t277;
t244 = t285 * t263 - t282 * t265;
t292 = -t244 * qJD(1) + t282 * t257 - t285 * t264;
t296 = -t280 * t303 + t283 * t314;
t311 = qJD(1) * t345;
t307 = t285 * t311;
t207 = t296 * qJD(4) + t280 * t307 + t283 * t292;
t243 = -t285 * t262 - t282 * t302;
t258 = t279 * t264;
t297 = t265 * qJD(2);
t222 = t243 * qJD(1) - t282 * t258 - t285 * t297;
t202 = t207 * t277 + t222 * t275;
t218 = t242 * t277 - t246 * t275;
t212 = 0.1e1 / t218;
t213 = 0.1e1 / t218 ^ 2;
t340 = t202 * t212 * t213;
t310 = 0.2e1 * t217 * t340;
t236 = t244 * t280 + t283 * t312;
t313 = t284 * t345;
t315 = t281 * t345;
t293 = -t276 * t313 - t278 * t315;
t251 = -t279 * t283 - t280 * t293;
t231 = atan2(-t236, t251);
t226 = sin(t231);
t227 = cos(t231);
t200 = -t226 * t236 + t227 * t251;
t198 = 0.1e1 / t200 ^ 2;
t235 = t296 ^ 2;
t196 = t235 * t198 + 0.1e1;
t206 = t242 * qJD(4) + t280 * t292 - t283 * t307;
t339 = t206 * t198;
t234 = t236 ^ 2;
t249 = 0.1e1 / t251 ^ 2;
t230 = t234 * t249 + 0.1e1;
t228 = 0.1e1 / t230;
t305 = t282 * t311;
t324 = qJD(4) * t283;
t208 = t244 * t324 + t350 * t280 - t283 * t305;
t252 = t279 * t280 - t283 * t293;
t260 = -t276 * t315 + t278 * t313;
t256 = t260 * qJD(2);
t232 = t252 * qJD(4) + t256 * t280;
t248 = 0.1e1 / t251;
t333 = t236 * t249;
t301 = -t208 * t248 + t232 * t333;
t190 = t301 * t228;
t304 = -t226 * t251 - t227 * t236;
t185 = t304 * t190 - t226 * t208 + t227 * t232;
t197 = 0.1e1 / t200;
t199 = t197 * t198;
t343 = t185 * t199;
t323 = 0.2e1 * (-t235 * t343 - t296 * t339) / t196 ^ 2;
t349 = t232 * t249;
t299 = -t243 * t248 + t260 * t333;
t348 = t280 * t299;
t209 = (-qJD(4) * t244 + t305) * t280 + t350 * t283;
t347 = -0.2e1 * t236;
t346 = -0.2e1 * t296;
t335 = t248 * t349;
t342 = (t208 * t333 - t234 * t335) / t230 ^ 2;
t341 = t198 * t296;
t338 = t213 * t217;
t337 = t226 * t296;
t336 = t227 * t296;
t334 = t236 * t248;
t332 = t246 * t280;
t331 = t246 * t283;
t330 = t275 * t212;
t329 = t277 * t217;
t201 = t207 * t275 - t222 * t277;
t211 = t217 ^ 2;
t205 = t211 * t213 + 0.1e1;
t322 = 0.2e1 * (t201 * t338 - t211 * t340) / t205 ^ 2;
t321 = -0.2e1 * t342;
t320 = t199 * t346;
t319 = t248 * t342;
t318 = t198 * t337;
t317 = t198 * t336;
t309 = t335 * t347;
t238 = t244 * t283 - t280 * t312;
t300 = -t238 * t248 + t252 * t333;
t298 = -qJD(4) * t332 - t222 * t283;
t295 = -t226 + (t227 * t334 + t226) * t228;
t255 = t293 * qJD(2);
t233 = -t251 * qJD(4) + t256 * t283;
t224 = t246 * qJD(1) - t285 * t258 + t282 * t297;
t220 = t275 * t303 + t277 * t331;
t219 = t275 * t331 - t277 * t303;
t216 = -t238 * t277 + t243 * t275;
t215 = -t238 * t275 - t243 * t277;
t203 = 0.1e1 / t205;
t194 = 0.1e1 / t196;
t193 = t228 * t348;
t192 = t300 * t228;
t189 = t295 * t296;
t187 = (-t226 * t243 + t227 * t260) * t280 + t304 * t193;
t186 = t304 * t192 - t226 * t238 + t227 * t252;
t184 = t300 * t321 + (t252 * t309 - t209 * t248 + (t208 * t252 + t232 * t238 + t233 * t236) * t249) * t228;
t182 = t321 * t348 + (t299 * t324 + (t260 * t309 - t224 * t248 + (t208 * t260 + t232 * t243 + t236 * t255) * t249) * t280) * t228;
t1 = [t319 * t346 + (-t206 * t248 - t296 * t349) * t228, t182, 0, t184, 0, 0; t236 * t197 * t323 + (-t208 * t197 + (t185 * t236 + t189 * t206) * t198) * t194 - (-t189 * t198 * t323 + (-0.2e1 * t189 * t343 + (-t190 * t228 * t334 + t321) * t318 + (t319 * t347 - t190 + (t190 - t301) * t228) * t317 - t295 * t339) * t194) * t296 (-t187 * t341 - t197 * t332) * t323 + (-t187 * t339 + (-t222 * t280 + t246 * t324) * t197 + (t187 * t320 - t198 * t332) * t185 + (t260 * t324 - t182 * t236 - t193 * t208 + t255 * t280 + (-t193 * t251 - t243 * t280) * t190) * t317 + (-t243 * t324 - t182 * t251 - t193 * t232 - t224 * t280 + (t193 * t236 - t260 * t280) * t190) * t318) * t194, 0 (-t186 * t341 - t197 * t242) * t323 + (t186 * t185 * t320 + t207 * t197 + (-t242 * t185 - t186 * t206 + (-t184 * t236 - t192 * t208 + t233 + (-t192 * t251 - t238) * t190) * t336 + (-t184 * t251 - t192 * t232 - t209 + (t192 * t236 - t252) * t190) * t337) * t198) * t194, 0, 0; (-t212 * t215 + t216 * t338) * t322 + ((-t209 * t275 - t224 * t277) * t212 + t216 * t310 + (-t215 * t202 - (-t209 * t277 + t224 * t275) * t217 - t216 * t201) * t213) * t203 (-t212 * t219 + t220 * t338) * t322 + ((t298 * t275 - t277 * t292) * t212 + t220 * t310 + (-t219 * t202 - (t275 * t292 + t298 * t277) * t217 - t220 * t201) * t213) * t203, 0 -(-t213 * t329 + t330) * t296 * t322 + (t296 * t277 * t310 - t206 * t330 + (t206 * t329 - (t201 * t277 + t202 * t275) * t296) * t213) * t203, 0, 0;];
JaD_rot  = t1;
