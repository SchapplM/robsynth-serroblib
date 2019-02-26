% Zeitableitung der rotatorischen Teilmatrix der analytischen Jacobi-Matrix für Segment Nr. 5 (0=Basis) von
% S6RPRPRR11
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
% pkin [13x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,alpha3,d1,d3,d5,d6,theta2,theta4]';
%
% Output:
% JaD_rot [3x6]
%   Zeitableitung der rotatorischen Teilmatrix der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-02-26 20:54
% Revision: d75aae1ac561373cd3be920984c3df29a1c2ecc8 (2019-02-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JaD_rot = S6RPRPRR11_jacobiaD_rot_5_sym_varpar(qJ, qJD, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(13,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRPRR11_jacobiaD_rot_5_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RPRPRR11_jacobiaD_rot_5_sym_varpar: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [13 1]), ...
  'S6RPRPRR11_jacobiaD_rot_5_sym_varpar: pkin has to be [13x1] (double)');

%% Symbolic Calculation
% From jacobiaD_rot_5_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 20:54:32
% EndTime: 2019-02-26 20:54:33
% DurationCPUTime: 1.39s
% Computational Cost: add. (5146->112), mult. (15324->225), div. (448->12), fcn. (19446->15), ass. (0->109)
t355 = sin(pkin(12));
t360 = cos(pkin(6));
t329 = t360 * t355;
t358 = cos(pkin(12));
t361 = sin(qJ(1));
t363 = cos(qJ(1));
t284 = t363 * t329 + t361 * t358;
t297 = sin(qJ(3));
t356 = sin(pkin(7));
t357 = sin(pkin(6));
t327 = t357 * t356;
t320 = t363 * t327;
t291 = t297 * t320;
t331 = t360 * t358;
t283 = -t363 * t331 + t361 * t355;
t359 = cos(pkin(7));
t335 = t283 * t359;
t362 = cos(qJ(3));
t369 = -t284 * t362 + t297 * t335 + t291;
t313 = t362 * t320;
t333 = t359 * t362;
t368 = -t283 * t333 - t313;
t326 = t357 * t355;
t328 = t359 * t357;
t367 = t358 * t328 + t360 * t356;
t273 = t367 * t297 + t362 * t326;
t265 = t273 * qJD(3);
t304 = -t297 * t326 + t367 * t362;
t270 = 0.1e1 / t304 ^ 2;
t344 = t265 * t270;
t285 = -t361 * t329 + t363 * t358;
t309 = t361 * t331 + t363 * t355;
t307 = t309 * t359;
t317 = t361 * t327;
t312 = t362 * t317;
t366 = -t285 * t297 - t362 * t307 + t312;
t281 = t309 * qJD(1);
t282 = t285 * qJD(1);
t239 = (qJD(1) * t317 - qJD(3) * t284 - t359 * t281) * t297 + t282 * t362 + t368 * qJD(3);
t365 = qJD(1) * t312 + t369 * qJD(3) - t281 * t333 - t282 * t297;
t257 = t284 * t297 - t368;
t254 = atan2(-t257, -t304);
t249 = sin(t254);
t250 = cos(t254);
t230 = -t249 * t257 - t250 * t304;
t227 = 0.1e1 / t230;
t263 = t285 * t362 + (-t307 + t317) * t297;
t318 = t361 * t328;
t275 = t309 * t356 + t318;
t296 = pkin(13) + qJ(5);
t294 = sin(t296);
t295 = cos(t296);
t248 = t263 * t295 + t275 * t294;
t242 = 0.1e1 / t248;
t269 = 0.1e1 / t304;
t228 = 0.1e1 / t230 ^ 2;
t243 = 0.1e1 / t248 ^ 2;
t364 = -0.2e1 * t366;
t256 = t366 ^ 2;
t226 = t228 * t256 + 0.1e1;
t280 = t284 * qJD(1);
t306 = qJD(1) * t335;
t236 = -qJD(1) * t313 + t263 * qJD(3) - t280 * t297 - t362 * t306;
t350 = t228 * t366;
t255 = t257 ^ 2;
t253 = t255 * t270 + 0.1e1;
t251 = 0.1e1 / t253;
t323 = t257 * t344 - t269 * t365;
t221 = t323 * t251;
t325 = t249 * t304 - t250 * t257;
t217 = t325 * t221 + t249 * t365 + t250 * t265;
t353 = t217 * t227 * t228;
t354 = (-t236 * t350 - t256 * t353) / t226 ^ 2;
t343 = t269 * t344;
t352 = (-t257 * t270 * t365 + t255 * t343) / t253 ^ 2;
t224 = 0.1e1 / t226;
t351 = t224 * t228;
t237 = qJD(1) * t291 + t366 * qJD(3) - t280 * t362 + t297 * t306;
t274 = -t283 * t356 + t363 * t328;
t267 = t274 * qJD(1);
t247 = t263 * t294 - t275 * t295;
t341 = qJD(5) * t247;
t232 = t237 * t295 + t267 * t294 - t341;
t349 = t232 * t242 * t243;
t231 = t248 * qJD(5) + t237 * t294 - t267 * t295;
t241 = t247 ^ 2;
t235 = t241 * t243 + 0.1e1;
t347 = t243 * t247;
t348 = 0.1e1 / t235 ^ 2 * (t231 * t347 - t241 * t349);
t346 = t257 * t269;
t345 = t257 * t273;
t340 = 0.2e1 * t354;
t339 = -0.2e1 * t352;
t338 = -0.2e1 * t348;
t337 = t269 * t352;
t336 = t247 * t349;
t334 = t353 * t364;
t246 = t274 * t294 + t295 * t369;
t245 = -t274 * t295 + t294 * t369;
t322 = -t242 * t294 + t295 * t347;
t321 = -t269 * t369 + t270 * t345;
t315 = -t249 + (-t250 * t346 + t249) * t251;
t268 = -qJD(1) * t318 - t281 * t356;
t264 = t304 * qJD(3);
t233 = 0.1e1 / t235;
t222 = t321 * t251;
t218 = t325 * t222 + t249 * t369 + t250 * t273;
t216 = t321 * t339 + (0.2e1 * t343 * t345 + t239 * t269 + (t257 * t264 - t265 * t369 - t273 * t365) * t270) * t251;
t1 = [-t337 * t364 + (t236 * t269 - t344 * t366) * t251, 0, t216, 0, 0, 0; -(-t217 * t351 - 0.2e1 * t227 * t354) * t257 + (t365 * t227 + (t315 * t236 - ((t221 * t251 * t346 + t339) * t249 + (0.2e1 * t257 * t337 - t221 + (t221 - t323) * t251) * t250) * t366) * t350) * t224 - (t224 * t334 - t236 * t351 - t350 * t340) * t315 * t366, 0 (-t218 * t350 - t227 * t263) * t340 + (t218 * t334 + t237 * t227 + (-t263 * t217 - t218 * t236 - (-(-t216 * t257 + t222 * t365 + t264 + (t222 * t304 + t369) * t221) * t250 - (t216 * t304 - t222 * t265 - t239 + (t222 * t257 - t273) * t221) * t249) * t366) * t228) * t224, 0, 0, 0; 0.2e1 * (-t242 * t245 + t246 * t347) * t348 + ((t246 * qJD(5) - t239 * t294 - t268 * t295) * t242 + 0.2e1 * t246 * t336 + (-t245 * t232 - (-t245 * qJD(5) - t239 * t295 + t268 * t294) * t247 - t246 * t231) * t243) * t233, 0, -t322 * t366 * t338 + (t322 * t236 - ((-qJD(5) * t242 - 0.2e1 * t336) * t295 + (t231 * t295 + (t232 - t341) * t294) * t243) * t366) * t233, 0, t338 + 0.2e1 * (t231 * t243 * t233 + (-t233 * t349 - t243 * t348) * t247) * t247, 0;];
JaD_rot  = t1;
