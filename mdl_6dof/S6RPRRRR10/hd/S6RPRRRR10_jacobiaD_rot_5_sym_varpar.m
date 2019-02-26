% Zeitableitung der rotatorischen Teilmatrix der analytischen Jacobi-Matrix für Segment Nr. 5 (0=Basis) von
% S6RPRRRR10
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
%   pkin=[a2,a3,a4,a5,a6,alpha2,alpha3,d1,d3,d4,d5,d6,theta2]';
%
% Output:
% JaD_rot [3x6]
%   Zeitableitung der rotatorischen Teilmatrix der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-02-26 21:20
% Revision: d75aae1ac561373cd3be920984c3df29a1c2ecc8 (2019-02-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JaD_rot = S6RPRRRR10_jacobiaD_rot_5_sym_varpar(qJ, qJD, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(13,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRRRR10_jacobiaD_rot_5_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RPRRRR10_jacobiaD_rot_5_sym_varpar: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [13 1]), ...
  'S6RPRRRR10_jacobiaD_rot_5_sym_varpar: pkin has to be [13x1] (double)');

%% Symbolic Calculation
% From jacobiaD_rot_5_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 21:19:50
% EndTime: 2019-02-26 21:19:52
% DurationCPUTime: 1.46s
% Computational Cost: add. (5639->113), mult. (16379->226), div. (466->12), fcn. (20741->15), ass. (0->112)
t377 = sin(pkin(13));
t382 = cos(pkin(6));
t348 = t382 * t377;
t380 = cos(pkin(13));
t383 = sin(qJ(1));
t385 = cos(qJ(1));
t302 = t385 * t348 + t383 * t380;
t316 = sin(qJ(3));
t378 = sin(pkin(7));
t379 = sin(pkin(6));
t346 = t379 * t378;
t339 = t385 * t346;
t309 = t316 * t339;
t350 = t382 * t380;
t301 = -t385 * t350 + t383 * t377;
t381 = cos(pkin(7));
t358 = t301 * t381;
t384 = cos(qJ(3));
t391 = -t302 * t384 + t316 * t358 + t309;
t332 = t384 * t339;
t352 = t381 * t384;
t390 = -t301 * t352 - t332;
t345 = t379 * t377;
t347 = t381 * t379;
t389 = t380 * t347 + t382 * t378;
t291 = t389 * t316 + t384 * t345;
t283 = t291 * qJD(3);
t323 = -t316 * t345 + t389 * t384;
t288 = 0.1e1 / t323 ^ 2;
t366 = t283 * t288;
t303 = -t383 * t348 + t385 * t380;
t328 = t383 * t350 + t385 * t377;
t326 = t328 * t381;
t336 = t383 * t346;
t331 = t384 * t336;
t388 = -t303 * t316 - t384 * t326 + t331;
t299 = t328 * qJD(1);
t300 = t303 * qJD(1);
t257 = (qJD(1) * t336 - qJD(3) * t302 - t381 * t299) * t316 + t300 * t384 + t390 * qJD(3);
t387 = qJD(1) * t331 + t391 * qJD(3) - t299 * t352 - t300 * t316;
t275 = t302 * t316 - t390;
t272 = atan2(-t275, -t323);
t267 = sin(t272);
t268 = cos(t272);
t250 = -t267 * t275 - t268 * t323;
t245 = 0.1e1 / t250;
t281 = t303 * t384 + (-t326 + t336) * t316;
t337 = t383 * t347;
t293 = t328 * t378 + t337;
t315 = qJ(4) + qJ(5);
t312 = sin(t315);
t313 = cos(t315);
t266 = t281 * t313 + t293 * t312;
t260 = 0.1e1 / t266;
t287 = 0.1e1 / t323;
t246 = 0.1e1 / t250 ^ 2;
t261 = 0.1e1 / t266 ^ 2;
t386 = -0.2e1 * t388;
t274 = t388 ^ 2;
t244 = t246 * t274 + 0.1e1;
t298 = t302 * qJD(1);
t325 = qJD(1) * t358;
t254 = -qJD(1) * t332 + t281 * qJD(3) - t298 * t316 - t384 * t325;
t372 = t246 * t388;
t273 = t275 ^ 2;
t271 = t273 * t288 + 0.1e1;
t269 = 0.1e1 / t271;
t342 = t275 * t366 - t287 * t387;
t239 = t342 * t269;
t344 = t267 * t323 - t268 * t275;
t235 = t344 * t239 + t267 * t387 + t268 * t283;
t375 = t235 * t245 * t246;
t376 = (-t254 * t372 - t274 * t375) / t244 ^ 2;
t365 = t287 * t366;
t374 = (-t275 * t288 * t387 + t273 * t365) / t271 ^ 2;
t242 = 0.1e1 / t244;
t373 = t242 * t246;
t292 = -t301 * t378 + t385 * t347;
t314 = qJD(4) + qJD(5);
t353 = -t292 * qJD(1) + t281 * t314;
t255 = qJD(1) * t309 + t388 * qJD(3) - t298 * t384 + t316 * t325;
t356 = t293 * t314 + t255;
t249 = -t353 * t312 + t356 * t313;
t371 = t249 * t260 * t261;
t248 = t356 * t312 + t353 * t313;
t265 = t281 * t312 - t293 * t313;
t259 = t265 ^ 2;
t253 = t259 * t261 + 0.1e1;
t369 = t261 * t265;
t370 = 0.1e1 / t253 ^ 2 * (t248 * t369 - t259 * t371);
t368 = t275 * t287;
t367 = t275 * t291;
t363 = 0.2e1 * t376;
t362 = -0.2e1 * t374;
t361 = -0.2e1 * t370;
t360 = t287 * t374;
t359 = t265 * t371;
t357 = t375 * t386;
t355 = t292 * t314 - t257;
t354 = qJD(1) * t337 + t299 * t378 + t314 * t391;
t341 = -t260 * t312 + t313 * t369;
t340 = -t287 * t391 + t288 * t367;
t334 = -t267 + (-t268 * t368 + t267) * t269;
t282 = t323 * qJD(3);
t264 = t292 * t312 + t313 * t391;
t263 = -t292 * t313 + t312 * t391;
t251 = 0.1e1 / t253;
t240 = t340 * t269;
t236 = t344 * t240 + t267 * t391 + t268 * t291;
t234 = t340 * t362 + (0.2e1 * t365 * t367 + t257 * t287 + (t275 * t282 - t283 * t391 - t291 * t387) * t288) * t269;
t232 = t361 + 0.2e1 * (t248 * t261 * t251 + (-t251 * t371 - t261 * t370) * t265) * t265;
t1 = [-t360 * t386 + (t254 * t287 - t366 * t388) * t269, 0, t234, 0, 0, 0; -(-t235 * t373 - 0.2e1 * t245 * t376) * t275 + (t387 * t245 + (t334 * t254 - ((t239 * t269 * t368 + t362) * t267 + (0.2e1 * t275 * t360 - t239 + (t239 - t342) * t269) * t268) * t388) * t372) * t242 - (t242 * t357 - t254 * t373 - t372 * t363) * t334 * t388, 0 (-t236 * t372 - t245 * t281) * t363 + (t236 * t357 + t255 * t245 + (-t281 * t235 - t236 * t254 - (-(-t234 * t275 + t240 * t387 + t282 + (t240 * t323 + t391) * t239) * t268 - (t234 * t323 - t240 * t283 - t257 + (t240 * t275 - t291) * t239) * t267) * t388) * t246) * t242, 0, 0, 0; 0.2e1 * (-t260 * t263 + t264 * t369) * t370 + ((t355 * t312 + t354 * t313) * t260 + 0.2e1 * t264 * t359 + (-t263 * t249 - (-t354 * t312 + t355 * t313) * t265 - t264 * t248) * t261) * t251, 0, -t341 * t388 * t361 + (t341 * t254 - ((-t260 * t314 - 0.2e1 * t359) * t313 + (t248 * t313 + (-t265 * t314 + t249) * t312) * t261) * t388) * t251, t232, t232, 0;];
JaD_rot  = t1;
