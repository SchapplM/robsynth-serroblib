% Zeitableitung der rotatorischen Teilmatrix der analytischen Jacobi-Matrix für Segment Nr. 5 (0=Basis) von
% S6PRRRPR5
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
%   pkin=[a2,a3,a4,a5,a6,alpha2,alpha3,d2,d3,d4,d6,theta1,theta5]';
%
% Output:
% JaD_rot [3x6]
%   Zeitableitung der rotatorischen Teilmatrix der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-02-26 20:13
% Revision: d75aae1ac561373cd3be920984c3df29a1c2ecc8 (2019-02-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JaD_rot = S6PRRRPR5_jacobiaD_rot_5_sym_varpar(qJ, qJD, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(13,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRRRPR5_jacobiaD_rot_5_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6PRRRPR5_jacobiaD_rot_5_sym_varpar: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [13 1]), ...
  'S6PRRRPR5_jacobiaD_rot_5_sym_varpar: pkin has to be [13x1] (double)');

%% Symbolic Calculation
% From jacobiaD_rot_5_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 20:12:56
% EndTime: 2019-02-26 20:12:57
% DurationCPUTime: 1.22s
% Computational Cost: add. (5650->131), mult. (16984->252), div. (538->12), fcn. (21480->15), ass. (0->123)
t300 = sin(qJ(2));
t302 = cos(qJ(2));
t366 = cos(pkin(12));
t367 = cos(pkin(6));
t332 = t367 * t366;
t365 = sin(pkin(12));
t318 = -t365 * t300 + t302 * t332;
t282 = t318 * qJD(2);
t296 = sin(pkin(7));
t297 = sin(pkin(6));
t335 = t296 * t297 * t366;
t370 = -qJD(3) * t335 + t282;
t299 = sin(qJ(3));
t317 = -t300 * t332 - t365 * t302;
t312 = qJD(2) * t317;
t301 = cos(qJ(3));
t314 = t301 * t318;
t369 = qJD(3) * t314 + t299 * t312;
t298 = cos(pkin(7));
t313 = t318 * t299;
t309 = t298 * t313 - t301 * t317;
t350 = t298 * t301;
t236 = t309 * qJD(3) + t370 * t299 - t312 * t350;
t355 = t317 * t299;
t260 = -t298 * t314 + t301 * t335 - t355;
t258 = t260 ^ 2;
t346 = t301 * t302;
t349 = t299 * t300;
t323 = t298 * t346 - t349;
t340 = t296 * t367;
t273 = -t323 * t297 - t301 * t340;
t271 = 0.1e1 / t273 ^ 2;
t252 = t258 * t271 + 0.1e1;
t356 = t260 * t271;
t347 = t300 * t301;
t348 = t299 * t302;
t321 = t298 * t348 + t347;
t322 = t298 * t347 + t348;
t333 = qJD(3) * t340;
t256 = t299 * t333 + (t322 * qJD(2) + t321 * qJD(3)) * t297;
t270 = 0.1e1 / t273;
t357 = t256 * t270 * t271;
t368 = -0.2e1 * (t236 * t356 - t258 * t357) / t252 ^ 2;
t253 = atan2(-t260, t273);
t248 = sin(t253);
t249 = cos(t253);
t230 = -t248 * t260 + t249 * t273;
t227 = 0.1e1 / t230;
t331 = t367 * t365;
t316 = t300 * t331 - t366 * t302;
t315 = t366 * t300 + t302 * t331;
t339 = t297 * t365;
t334 = t296 * t339;
t319 = -t298 * t315 + t334;
t264 = t319 * t299 - t301 * t316;
t275 = t296 * t315 + t298 * t339;
t295 = qJ(4) + pkin(13);
t293 = sin(t295);
t294 = cos(t295);
t245 = t264 * t294 + t275 * t293;
t241 = 0.1e1 / t245;
t228 = 0.1e1 / t230 ^ 2;
t242 = 0.1e1 / t245 ^ 2;
t250 = 0.1e1 / t252;
t220 = (-t236 * t270 + t256 * t356) * t250;
t330 = -t248 * t273 - t249 * t260;
t216 = t330 * t220 - t248 * t236 + t249 * t256;
t364 = t216 * t227 * t228;
t283 = t315 * qJD(2);
t284 = t316 * qJD(2);
t351 = t298 * t299;
t354 = t316 * t299;
t239 = t284 * t351 - t283 * t301 + (t319 * t301 + t354) * qJD(3);
t352 = t294 * t296;
t231 = t245 * qJD(4) + t239 * t293 + t284 * t352;
t244 = t264 * t293 - t275 * t294;
t240 = t244 ^ 2;
t235 = t240 * t242 + 0.1e1;
t360 = t242 * t244;
t345 = qJD(4) * t244;
t353 = t293 * t296;
t232 = t239 * t294 - t284 * t353 - t345;
t361 = t232 * t241 * t242;
t363 = (t231 * t360 - t240 * t361) / t235 ^ 2;
t263 = -t301 * t334 + t315 * t350 - t354;
t362 = t228 * t263;
t359 = t248 * t263;
t358 = t249 * t263;
t259 = t263 ^ 2;
t226 = t259 * t228 + 0.1e1;
t238 = t264 * qJD(3) - t283 * t299 - t284 * t350;
t344 = 0.2e1 * (t238 * t362 - t259 * t364) / t226 ^ 2;
t343 = -0.2e1 * t363;
t342 = t244 * t361;
t341 = qJD(3) * t355;
t337 = 0.2e1 * t263 * t364;
t336 = -0.2e1 * t260 * t357;
t327 = -t293 * t241 + t294 * t360;
t262 = -t299 * t335 + t309;
t274 = t321 * t297 + t299 * t340;
t326 = -t262 * t270 + t274 * t356;
t266 = -t317 * t350 + t313;
t281 = t322 * t297;
t325 = -t266 * t270 + t281 * t356;
t268 = -t301 * t315 + t316 * t351;
t324 = -t268 * t293 - t316 * t352;
t255 = t268 * t294 - t316 * t353;
t267 = -t299 * t315 - t316 * t350;
t320 = -t298 * t349 + t346;
t265 = (t323 * qJD(2) + t320 * qJD(3)) * t297;
t257 = t301 * t333 + (t320 * qJD(2) + t323 * qJD(3)) * t297;
t247 = -t267 * qJD(3) + t283 * t351 + t284 * t301;
t246 = t282 * t350 + t298 * t341 + t369;
t237 = t369 * t298 + t370 * t301 + t341;
t233 = 0.1e1 / t235;
t224 = 0.1e1 / t226;
t222 = t325 * t250;
t221 = t326 * t250;
t218 = t330 * t222 - t248 * t266 + t249 * t281;
t217 = t330 * t221 - t248 * t262 + t249 * t274;
t215 = t325 * t368 + (t281 * t336 - t246 * t270 + (t236 * t281 + t256 * t266 + t260 * t265) * t271) * t250;
t214 = t326 * t368 + (t274 * t336 - t237 * t270 + (t236 * t274 + t256 * t262 + t257 * t260) * t271) * t250;
t1 = [0, t215, t214, 0, 0, 0; 0 (t218 * t362 - t227 * t267) * t344 + ((t268 * qJD(3) - t283 * t350 + t284 * t299) * t227 + t218 * t337 + (-t267 * t216 - t218 * t238 - (-t215 * t260 - t222 * t236 + t265 + (-t222 * t273 - t266) * t220) * t358 - (-t215 * t273 - t222 * t256 - t246 + (t222 * t260 - t281) * t220) * t359) * t228) * t224 (t217 * t362 - t227 * t264) * t344 + (t217 * t337 + t239 * t227 + (-t264 * t216 - t217 * t238 - (-t214 * t260 - t221 * t236 + t257 + (-t221 * t273 - t262) * t220) * t358 - (-t214 * t273 - t221 * t256 - t237 + (t221 * t260 - t274) * t220) * t359) * t228) * t224, 0, 0, 0; 0, 0.2e1 * (t241 * t324 + t255 * t360) * t363 + ((t255 * qJD(4) + t247 * t293 + t283 * t352) * t241 + 0.2e1 * t255 * t342 + (t324 * t232 - (t324 * qJD(4) + t247 * t294 - t283 * t353) * t244 - t255 * t231) * t242) * t233, t327 * t263 * t343 + (t327 * t238 + ((-qJD(4) * t241 - 0.2e1 * t342) * t294 + (t231 * t294 + (t232 - t345) * t293) * t242) * t263) * t233, t343 + 0.2e1 * (t231 * t242 * t233 + (-t233 * t361 - t242 * t363) * t244) * t244, 0, 0;];
JaD_rot  = t1;
