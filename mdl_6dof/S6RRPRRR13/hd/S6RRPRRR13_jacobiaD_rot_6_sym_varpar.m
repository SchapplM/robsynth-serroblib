% Zeitableitung der rotatorischen Teilmatrix der analytischen Jacobi-Matrix für Segment Nr. 6 (0=Basis) von
% S6RRPRRR13
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
%   pkin=[a2,a3,a4,a5,a6,alpha2,d1,d2,d4,d5,d6]';
%
% Output:
% JaD_rot [3x6]
%   Zeitableitung der rotatorischen Teilmatrix der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-02-26 22:01
% Revision: d75aae1ac561373cd3be920984c3df29a1c2ecc8 (2019-02-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JaD_rot = S6RRPRRR13_jacobiaD_rot_6_sym_varpar(qJ, qJD, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(11,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPRRR13_jacobiaD_rot_6_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRPRRR13_jacobiaD_rot_6_sym_varpar: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6RRPRRR13_jacobiaD_rot_6_sym_varpar: pkin has to be [11x1] (double)');

%% Symbolic Calculation
% From jacobiaD_rot_6_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 22:01:09
% EndTime: 2019-02-26 22:01:11
% DurationCPUTime: 1.84s
% Computational Cost: add. (5342->149), mult. (14155->296), div. (744->12), fcn. (17850->13), ass. (0->130)
t283 = sin(qJ(2));
t284 = sin(qJ(1));
t286 = cos(qJ(2));
t287 = cos(qJ(1));
t369 = cos(pkin(6));
t321 = t287 * t369;
t302 = -t283 * t321 - t284 * t286;
t322 = t284 * t369;
t304 = t287 * t283 + t286 * t322;
t281 = sin(pkin(6));
t348 = t281 * t287;
t382 = t304 * qJD(1) - t302 * qJD(2) + qJD(4) * t348;
t268 = t284 * t283 - t286 * t321;
t282 = sin(qJ(4));
t285 = cos(qJ(4));
t306 = t268 * t285 + t282 * t348;
t251 = t306 ^ 2;
t349 = t281 * t286;
t303 = -t369 * t282 - t285 * t349;
t264 = 0.1e1 / t303 ^ 2;
t242 = t251 * t264 + 0.1e1;
t240 = 0.1e1 / t242;
t267 = -t282 * t349 + t369 * t285;
t342 = qJD(2) * t283;
t325 = t281 * t342;
t253 = t267 * qJD(4) - t285 * t325;
t263 = 0.1e1 / t303;
t355 = t306 * t264;
t343 = qJD(1) * t284;
t370 = t268 * qJD(4) + t281 * t343;
t380 = -t370 * t282 + t382 * t285;
t310 = -t253 * t355 - t263 * t380;
t209 = t310 * t240;
t243 = atan2(t306, -t303);
t238 = sin(t243);
t239 = cos(t243);
t312 = t238 * t303 + t239 * t306;
t204 = t312 * t209 + t238 * t380 + t239 * t253;
t221 = t238 * t306 - t239 * t303;
t219 = 0.1e1 / t221 ^ 2;
t381 = t204 * t219;
t226 = t382 * t282 + t370 * t285;
t350 = t281 * t285;
t255 = t304 * t282 + t284 * t350;
t313 = t283 * t322;
t345 = t287 * t286;
t270 = -t313 + t345;
t280 = qJ(5) + qJ(6);
t277 = sin(t280);
t278 = cos(t280);
t234 = t255 * t277 - t270 * t278;
t379 = 0.2e1 * t234;
t218 = 0.1e1 / t221;
t378 = t218 * t381;
t372 = -t284 * t281 * t282 + t304 * t285;
t318 = -0.2e1 * t372 * t378;
t298 = (t369 * qJD(1) + qJD(2)) * t345 - qJD(2) * t313 - t283 * t343;
t326 = qJD(1) * t348;
t228 = t255 * qJD(4) + t282 * t326 - t298 * t285;
t361 = t228 * t219;
t377 = -t361 + t318;
t375 = t253 * t264;
t351 = t281 * t283;
t330 = t306 * t351;
t305 = t263 * t302 + t264 * t330;
t374 = t285 * t305;
t279 = qJD(5) + qJD(6);
t353 = t270 * t279;
t373 = -t282 * t353 - t298;
t352 = t270 * t285;
t371 = qJD(4) * t352 - t304 * t279;
t235 = t255 * t278 + t270 * t277;
t231 = 0.1e1 / t235;
t232 = 0.1e1 / t235 ^ 2;
t250 = t372 ^ 2;
t217 = t250 * t219 + 0.1e1;
t368 = (-t250 * t378 - t361 * t372) / t217 ^ 2;
t248 = t302 * qJD(1) - t304 * qJD(2);
t315 = t255 * t279 - t248;
t229 = t372 * qJD(4) + t298 * t282 + t285 * t326;
t316 = t229 + t353;
t212 = t316 * t277 + t315 * t278;
t230 = t234 ^ 2;
t224 = t230 * t232 + 0.1e1;
t360 = t232 * t234;
t213 = -t315 * t277 + t316 * t278;
t364 = t213 * t231 * t232;
t367 = (t212 * t360 - t230 * t364) / t224 ^ 2;
t357 = t263 * t375;
t365 = (t251 * t357 + t355 * t380) / t242 ^ 2;
t363 = t219 * t372;
t222 = 0.1e1 / t224;
t362 = t222 * t232;
t359 = t238 * t372;
t358 = t239 * t372;
t356 = t306 * t263;
t354 = t306 * t267;
t347 = t282 * t277;
t346 = t282 * t278;
t341 = qJD(2) * t286;
t340 = qJD(4) * t282;
t338 = 0.2e1 * t368;
t337 = -0.2e1 * t367;
t336 = -0.2e1 * t365;
t335 = 0.2e1 * t365;
t333 = t232 * t367;
t332 = t212 * t362;
t331 = t234 * t364;
t320 = t263 * t335;
t319 = 0.2e1 * t331;
t317 = t279 * t302 - t226;
t249 = -qJD(1) * t313 - t284 * t342 + (qJD(2) * t369 + qJD(1)) * t345;
t307 = -t268 * t282 + t285 * t348;
t314 = t279 * t307 + t249;
t309 = -t277 * t231 + t278 * t360;
t308 = t263 * t307 + t264 * t354;
t301 = -t238 + (t239 * t356 + t238) * t240;
t252 = t303 * qJD(4) + t282 * t325;
t245 = t270 * t346 - t304 * t277;
t237 = t277 * t302 + t278 * t307;
t236 = t277 * t307 - t278 * t302;
t215 = 0.1e1 / t217;
t214 = t240 * t374;
t211 = t308 * t240;
t206 = (-t238 * t302 - t239 * t351) * t285 + t312 * t214;
t205 = -t312 * t211 + t238 * t307 + t239 * t267;
t203 = t308 * t335 + (-0.2e1 * t354 * t357 + t226 * t263 + (-t252 * t306 - t253 * t307 - t267 * t380) * t264) * t240;
t201 = t336 * t374 + (-t305 * t340 + (0.2e1 * t330 * t357 - t249 * t263 + (t253 * t302 + (t283 * t380 + t306 * t341) * t281) * t264) * t285) * t240;
t200 = t337 + (t332 + (-t222 * t364 - t333) * t234) * t379;
t1 = [t372 * t320 + (t228 * t263 - t372 * t375) * t240, t201, 0, t203, 0, 0; -0.2e1 * t306 * t218 * t368 + (t380 * t218 - t306 * t381 + (t301 * t228 - ((-t209 * t240 * t356 + t336) * t238 + (-t306 * t320 - t209 + (t209 - t310) * t240) * t239) * t372) * t363) * t215 - (t377 * t215 - t363 * t338) * t301 * t372 (-t206 * t363 + t218 * t352) * t338 + ((-t248 * t285 + t270 * t340) * t218 + t377 * t206 + (t352 * t204 + (t201 * t306 + t214 * t380 + (t283 * t340 - t285 * t341) * t281 + (t214 * t303 - t285 * t302) * t209) * t358 + (t302 * t340 + t201 * t303 - t214 * t253 + t249 * t285 + (-t214 * t306 + t283 * t350) * t209) * t359) * t219) * t215, 0 (-t205 * t363 - t218 * t255) * t338 + (t205 * t318 + t229 * t218 + (-t255 * t204 - t205 * t228 + (t203 * t306 - t211 * t380 + t252 + (-t211 * t303 + t307) * t209) * t358 + (t203 * t303 + t211 * t253 - t226 + (t211 * t306 - t267) * t209) * t359) * t219) * t215, 0, 0; 0.2e1 * (-t231 * t236 + t237 * t360) * t367 + ((t317 * t277 + t314 * t278) * t231 + t237 * t319 + (-t236 * t213 - (-t314 * t277 + t317 * t278) * t234 - t237 * t212) * t232) * t222 (t333 * t379 - t332) * t245 + (-t213 * t362 + t231 * t337) * (t270 * t347 + t304 * t278) + (t245 * t319 + (t347 * t231 - t346 * t360) * t248 + (-t373 * t231 - t371 * t360) * t278 + (t371 * t231 - t373 * t360) * t277) * t222, 0, -t309 * t372 * t337 + (t309 * t228 - ((-t231 * t279 - 0.2e1 * t331) * t278 + (t212 * t278 + (-t234 * t279 + t213) * t277) * t232) * t372) * t222, t200, t200;];
JaD_rot  = t1;
