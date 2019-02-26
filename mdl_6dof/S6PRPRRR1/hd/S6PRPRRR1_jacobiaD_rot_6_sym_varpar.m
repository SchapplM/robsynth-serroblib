% Zeitableitung der rotatorischen Teilmatrix der analytischen Jacobi-Matrix für Segment Nr. 6 (0=Basis) von
% S6PRPRRR1
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
%   pkin=[a2,a3,a4,a5,a6,alpha2,d2,d4,d5,d6,theta1,theta3]';
%
% Output:
% JaD_rot [3x6]
%   Zeitableitung der rotatorischen Teilmatrix der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-02-26 19:53
% Revision: d75aae1ac561373cd3be920984c3df29a1c2ecc8 (2019-02-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JaD_rot = S6PRPRRR1_jacobiaD_rot_6_sym_varpar(qJ, qJD, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(12,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRPRRR1_jacobiaD_rot_6_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6PRPRRR1_jacobiaD_rot_6_sym_varpar: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [12 1]), ...
  'S6PRPRRR1_jacobiaD_rot_6_sym_varpar: pkin has to be [12x1] (double)');

%% Symbolic Calculation
% From jacobiaD_rot_6_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 19:53:44
% EndTime: 2019-02-26 19:53:45
% DurationCPUTime: 1.58s
% Computational Cost: add. (12871->119), mult. (23866->239), div. (822->12), fcn. (31086->15), ass. (0->117)
t319 = sin(pkin(12));
t322 = cos(pkin(12));
t326 = sin(qJ(2));
t328 = cos(qJ(2));
t309 = t326 * t319 - t328 * t322;
t324 = cos(pkin(6));
t337 = t309 * t324;
t302 = qJD(2) * t337;
t341 = t328 * t319 + t326 * t322;
t307 = t341 * qJD(2);
t320 = sin(pkin(11));
t323 = cos(pkin(11));
t283 = -t323 * t302 - t320 * t307;
t305 = t341 * t324;
t289 = t323 * t305 - t320 * t309;
t318 = qJ(4) + qJ(5);
t315 = sin(t318);
t317 = qJD(4) + qJD(5);
t321 = sin(pkin(6));
t358 = t321 * t323;
t349 = t315 * t358;
t316 = cos(t318);
t360 = t316 * t317;
t255 = t283 * t315 + t289 * t360 - t317 * t349;
t277 = t289 * t315 + t316 * t358;
t275 = t277 ^ 2;
t304 = t341 * t321;
t296 = t304 * t315 - t324 * t316;
t294 = 0.1e1 / t296 ^ 2;
t269 = t275 * t294 + 0.1e1;
t267 = 0.1e1 / t269;
t303 = t309 * t321;
t347 = -qJD(2) * t303 + t317 * t324;
t273 = t304 * t360 + t347 * t315;
t293 = 0.1e1 / t296;
t365 = t277 * t294;
t239 = (-t255 * t293 + t273 * t365) * t267;
t272 = atan2(-t277, t296);
t265 = sin(t272);
t266 = cos(t272);
t344 = -t265 * t296 - t266 * t277;
t235 = t344 * t239 - t265 * t255 + t266 * t273;
t249 = -t265 * t277 + t266 * t296;
t246 = 0.1e1 / t249;
t247 = 0.1e1 / t249 ^ 2;
t379 = t235 * t246 * t247;
t342 = -t320 * t305 - t323 * t309;
t359 = t320 * t321;
t280 = t315 * t342 - t316 * t359;
t378 = 0.2e1 * t280 * t379;
t288 = -t320 * t341 - t323 * t337;
t338 = -t288 * t293 - t303 * t365;
t377 = t315 * t338;
t366 = t273 * t293 * t294;
t376 = -0.2e1 * (t255 * t365 - t275 * t366) / t269 ^ 2;
t281 = t315 * t359 + t316 * t342;
t327 = cos(qJ(6));
t291 = t320 * t337 - t323 * t341;
t325 = sin(qJ(6));
t363 = t291 * t325;
t264 = t281 * t327 - t363;
t260 = 0.1e1 / t264;
t261 = 0.1e1 / t264 ^ 2;
t343 = t320 * t302 - t323 * t307;
t346 = t317 * t359 + t343;
t361 = t315 * t317;
t258 = t346 * t316 - t342 * t361;
t306 = t309 * qJD(2);
t336 = t324 * t307;
t284 = t323 * t306 + t320 * t336;
t250 = t264 * qJD(6) + t258 * t325 + t284 * t327;
t362 = t291 * t327;
t263 = t281 * t325 + t362;
t259 = t263 ^ 2;
t254 = t259 * t261 + 0.1e1;
t370 = t261 * t263;
t354 = qJD(6) * t263;
t251 = t258 * t327 - t284 * t325 - t354;
t373 = t251 * t260 * t261;
t375 = (t250 * t370 - t259 * t373) / t254 ^ 2;
t374 = t247 * t280;
t257 = t346 * t315 + t342 * t360;
t372 = t257 * t247;
t371 = t260 * t325;
t369 = t263 * t327;
t368 = t265 * t280;
t367 = t266 * t280;
t364 = t291 * t315;
t276 = t280 ^ 2;
t245 = t276 * t247 + 0.1e1;
t353 = 0.2e1 * (-t276 * t379 + t280 * t372) / t245 ^ 2;
t352 = -0.2e1 * t375;
t350 = t263 * t373;
t348 = -0.2e1 * t277 * t366;
t345 = qJD(6) * t291 * t316 - t343;
t340 = t261 * t369 - t371;
t279 = t289 * t316 - t349;
t297 = t304 * t316 + t324 * t315;
t339 = -t279 * t293 + t297 * t365;
t335 = qJD(6) * t342 + t284 * t316 - t291 * t361;
t300 = t321 * t307;
t282 = t320 * t306 - t323 * t336;
t274 = -t304 * t361 + t347 * t316;
t271 = t316 * t362 + t325 * t342;
t270 = t316 * t363 - t327 * t342;
t256 = -t289 * t361 + (-t317 * t358 + t283) * t316;
t252 = 0.1e1 / t254;
t243 = 0.1e1 / t245;
t241 = t267 * t377;
t240 = t339 * t267;
t237 = (-t265 * t288 - t266 * t303) * t315 + t344 * t241;
t236 = t344 * t240 - t265 * t279 + t266 * t297;
t233 = t339 * t376 + (t297 * t348 - t256 * t293 + (t255 * t297 + t273 * t279 + t274 * t277) * t294) * t267;
t232 = t376 * t377 + (t338 * t360 + (-t303 * t348 - t282 * t293 + (-t255 * t303 + t273 * t288 - t277 * t300) * t294) * t315) * t267;
t231 = t340 * t280 * t352 + (t340 * t257 + ((-qJD(6) * t260 - 0.2e1 * t350) * t327 + (t250 * t327 + (t251 - t354) * t325) * t261) * t280) * t252;
t230 = (t236 * t374 - t246 * t281) * t353 + (t236 * t378 + t258 * t246 + (-t281 * t235 - t236 * t257 - (-t233 * t277 - t240 * t255 + t274 + (-t240 * t296 - t279) * t239) * t367 - (-t233 * t296 - t240 * t273 - t256 + (t240 * t277 - t297) * t239) * t368) * t247) * t243;
t1 = [0, t232, 0, t233, t233, 0; 0 (t237 * t374 - t246 * t364) * t353 + ((t284 * t315 + t291 * t360) * t246 + (-t372 + t378) * t237 + (-t364 * t235 - (-t303 * t360 - t232 * t277 - t241 * t255 - t300 * t315 + (-t241 * t296 - t288 * t315) * t239) * t367 - (-t288 * t360 - t232 * t296 - t241 * t273 - t282 * t315 + (t241 * t277 + t303 * t315) * t239) * t368) * t247) * t243, 0, t230, t230, 0; 0, 0.2e1 * (-t260 * t270 + t271 * t370) * t375 + (0.2e1 * t271 * t350 + t345 * t260 * t327 + t335 * t371 + (t345 * t263 * t325 - t271 * t250 - t270 * t251 - t335 * t369) * t261) * t252, 0, t231, t231, t352 + 0.2e1 * (t250 * t261 * t252 + (-t252 * t373 - t261 * t375) * t263) * t263;];
JaD_rot  = t1;
