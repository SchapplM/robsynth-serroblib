% Zeitableitung der rotatorischen Teilmatrix der analytischen Jacobi-Matrix für Segment Nr. 6 (0=Basis) von
% S6RRPRRP13
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
%   pkin=[a2,a3,a4,a5,a6,alpha2,d1,d2,d4,d5]';
%
% Output:
% JaD_rot [3x6]
%   Zeitableitung der rotatorischen Teilmatrix der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-02-26 21:53
% Revision: d75aae1ac561373cd3be920984c3df29a1c2ecc8 (2019-02-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JaD_rot = S6RRPRRP13_jacobiaD_rot_6_sym_varpar(qJ, qJD, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(10,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPRRP13_jacobiaD_rot_6_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRPRRP13_jacobiaD_rot_6_sym_varpar: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RRPRRP13_jacobiaD_rot_6_sym_varpar: pkin has to be [10x1] (double)');

%% Symbolic Calculation
% From jacobiaD_rot_6_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 21:52:49
% EndTime: 2019-02-26 21:52:51
% DurationCPUTime: 1.83s
% Computational Cost: add. (4522->147), mult. (13478->296), div. (726->12), fcn. (17045->13), ass. (0->123)
t261 = sin(qJ(2));
t262 = sin(qJ(1));
t265 = cos(qJ(2));
t266 = cos(qJ(1));
t344 = cos(pkin(6));
t296 = t266 * t344;
t281 = -t261 * t296 - t262 * t265;
t297 = t262 * t344;
t283 = t266 * t261 + t265 * t297;
t258 = sin(pkin(6));
t324 = t258 * t266;
t357 = t283 * qJD(1) - t281 * qJD(2) + qJD(4) * t324;
t249 = t262 * t261 - t265 * t296;
t260 = sin(qJ(4));
t264 = cos(qJ(4));
t285 = t249 * t264 + t260 * t324;
t232 = t285 ^ 2;
t325 = t258 * t265;
t282 = -t344 * t260 - t264 * t325;
t245 = 0.1e1 / t282 ^ 2;
t223 = t232 * t245 + 0.1e1;
t221 = 0.1e1 / t223;
t248 = -t260 * t325 + t344 * t264;
t318 = qJD(2) * t261;
t301 = t258 * t318;
t234 = t248 * qJD(4) - t264 * t301;
t244 = 0.1e1 / t282;
t330 = t285 * t245;
t319 = qJD(1) * t262;
t345 = t249 * qJD(4) + t258 * t319;
t355 = -t345 * t260 + t357 * t264;
t289 = -t234 * t330 - t244 * t355;
t190 = t289 * t221;
t224 = atan2(t285, -t282);
t219 = sin(t224);
t220 = cos(t224);
t291 = t219 * t282 + t220 * t285;
t185 = t291 * t190 + t219 * t355 + t220 * t234;
t202 = t219 * t285 - t220 * t282;
t200 = 0.1e1 / t202 ^ 2;
t356 = t185 * t200;
t207 = t357 * t260 + t345 * t264;
t326 = t258 * t264;
t236 = t283 * t260 + t262 * t326;
t292 = t261 * t297;
t321 = t266 * t265;
t251 = -t292 + t321;
t259 = sin(qJ(5));
t263 = cos(qJ(5));
t215 = t236 * t259 - t251 * t263;
t354 = 0.2e1 * t215;
t199 = 0.1e1 / t202;
t353 = t199 * t356;
t346 = -t262 * t258 * t260 + t283 * t264;
t293 = -0.2e1 * t346 * t353;
t277 = (t344 * qJD(1) + qJD(2)) * t321 - qJD(2) * t292 - t261 * t319;
t302 = qJD(1) * t324;
t209 = t236 * qJD(4) + t260 * t302 - t277 * t264;
t336 = t209 * t200;
t352 = -t336 + t293;
t350 = t234 * t245;
t327 = t258 * t261;
t305 = t285 * t327;
t284 = t244 * t281 + t245 * t305;
t349 = t264 * t284;
t348 = -t251 * t260 * qJD(5) - t277;
t328 = t251 * t264;
t347 = qJD(4) * t328 - t283 * qJD(5);
t216 = t236 * t263 + t251 * t259;
t212 = 0.1e1 / t216;
t213 = 0.1e1 / t216 ^ 2;
t231 = t346 ^ 2;
t198 = t231 * t200 + 0.1e1;
t343 = (-t231 * t353 - t336 * t346) / t198 ^ 2;
t210 = t346 * qJD(4) + t277 * t260 + t264 * t302;
t229 = t281 * qJD(1) - t283 * qJD(2);
t194 = t216 * qJD(5) + t210 * t259 - t229 * t263;
t211 = t215 ^ 2;
t205 = t211 * t213 + 0.1e1;
t335 = t213 * t215;
t315 = qJD(5) * t215;
t195 = t210 * t263 + t229 * t259 - t315;
t339 = t195 * t212 * t213;
t342 = (t194 * t335 - t211 * t339) / t205 ^ 2;
t332 = t244 * t350;
t340 = (t232 * t332 + t330 * t355) / t223 ^ 2;
t338 = t200 * t346;
t203 = 0.1e1 / t205;
t337 = t203 * t213;
t334 = t219 * t346;
t333 = t220 * t346;
t331 = t285 * t244;
t329 = t285 * t248;
t323 = t260 * t259;
t322 = t260 * t263;
t317 = qJD(2) * t265;
t316 = qJD(4) * t260;
t313 = 0.2e1 * t343;
t312 = -0.2e1 * t342;
t311 = -0.2e1 * t340;
t310 = 0.2e1 * t340;
t308 = t213 * t342;
t307 = t194 * t337;
t306 = t215 * t339;
t295 = t244 * t310;
t294 = 0.2e1 * t306;
t286 = -t249 * t260 + t264 * t324;
t218 = t259 * t281 + t263 * t286;
t217 = t259 * t286 - t263 * t281;
t288 = -t259 * t212 + t263 * t335;
t287 = t244 * t286 + t245 * t329;
t280 = -t219 + (t220 * t331 + t219) * t221;
t233 = t282 * qJD(4) + t260 * t301;
t230 = -qJD(1) * t292 - t262 * t318 + (qJD(2) * t344 + qJD(1)) * t321;
t226 = t251 * t322 - t283 * t259;
t196 = 0.1e1 / t198;
t193 = t221 * t349;
t192 = t287 * t221;
t187 = (-t219 * t281 - t220 * t327) * t264 + t291 * t193;
t186 = -t291 * t192 + t219 * t286 + t220 * t248;
t184 = t287 * t310 + (-0.2e1 * t329 * t332 + t207 * t244 + (-t233 * t285 - t234 * t286 - t248 * t355) * t245) * t221;
t182 = t311 * t349 + (-t284 * t316 + (0.2e1 * t305 * t332 - t230 * t244 + (t234 * t281 + (t261 * t355 + t285 * t317) * t258) * t245) * t264) * t221;
t1 = [t346 * t295 + (t209 * t244 - t346 * t350) * t221, t182, 0, t184, 0, 0; -0.2e1 * t285 * t199 * t343 + (t355 * t199 - t285 * t356 + (t280 * t209 - ((-t190 * t221 * t331 + t311) * t219 + (-t285 * t295 - t190 + (t190 - t289) * t221) * t220) * t346) * t338) * t196 - (t196 * t352 - t338 * t313) * t280 * t346 (-t187 * t338 + t199 * t328) * t313 + ((-t229 * t264 + t251 * t316) * t199 + t352 * t187 + (t328 * t185 + (t182 * t285 + t193 * t355 + (t261 * t316 - t264 * t317) * t258 + (t193 * t282 - t264 * t281) * t190) * t333 + (t281 * t316 + t182 * t282 - t193 * t234 + t230 * t264 + (-t193 * t285 + t261 * t326) * t190) * t334) * t200) * t196, 0 (-t186 * t338 - t199 * t236) * t313 + (t186 * t293 + t210 * t199 + (-t236 * t185 - t186 * t209 + (t184 * t285 - t192 * t355 + t233 + (-t192 * t282 + t286) * t190) * t333 + (t184 * t282 + t192 * t234 - t207 + (t192 * t285 - t248) * t190) * t334) * t200) * t196, 0, 0; 0.2e1 * (-t212 * t217 + t218 * t335) * t342 + ((t218 * qJD(5) - t207 * t259 + t230 * t263) * t212 + t218 * t294 + (-t217 * t195 - (-t217 * qJD(5) - t207 * t263 - t230 * t259) * t215 - t218 * t194) * t213) * t203 (t308 * t354 - t307) * t226 + (-t195 * t337 + t212 * t312) * (t251 * t323 + t283 * t263) + (t226 * t294 + (t323 * t212 - t322 * t335) * t229 + (-t348 * t212 - t347 * t335) * t263 + (t347 * t212 - t348 * t335) * t259) * t203, 0, -t288 * t346 * t312 + (t288 * t209 - ((-qJD(5) * t212 - 0.2e1 * t306) * t263 + (t194 * t263 + (t195 - t315) * t259) * t213) * t346) * t203, t312 + (t307 + (-t203 * t339 - t308) * t215) * t354, 0;];
JaD_rot  = t1;
