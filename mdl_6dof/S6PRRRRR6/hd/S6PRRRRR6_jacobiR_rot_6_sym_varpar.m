% Rotatorische Teilmatrix der Rotationsmatrix-Jacobi-Matrix für Segment Nr. 6 (0=Basis) von
% S6PRRRRR6
% Use Code from Maple symbolic Code Generation
%
% Rotationsmatrix-Jacobi-Matrix: Differentieller Zusammenhang zwischen
% gestapelter Endeffektor-Rotationsmatrix und verallgemeinerten Koordinaten.
%
%
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% pkin [14x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,alpha3,alpha4,d2,d3,d4,d5,d6,theta1]';
%
% Output:
% JR_rot [9x6]
%   Jacobi-Matrix der Endeffektor-Rotationsmatrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-02-22 10:05
% Revision: 2b76964ad985d937eecd005a1a368749e6b3dc4d (2019-02-18)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JR_rot = S6PRRRRR6_jacobiR_rot_6_sym_varpar(qJ, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(14,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRRRRR6_jacobiR_rot_6_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [14 1]), ...
  'S6PRRRRR6_jacobiR_rot_6_sym_varpar: pkin has to be [14x1] (double)');

%% Symbolic Calculation
% From jacobiR_rot_6_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-22 10:05:19
% EndTime: 2019-02-22 10:05:20
% DurationCPUTime: 0.40s
% Computational Cost: add. (813->105), mult. (2422->221), div. (0->0), fcn. (3274->18), ass. (0->102)
t344 = sin(pkin(8));
t345 = sin(pkin(7));
t391 = t344 * t345;
t352 = sin(qJ(5));
t390 = t344 * t352;
t357 = cos(qJ(5));
t389 = t344 * t357;
t346 = sin(pkin(6));
t388 = t345 * t346;
t348 = cos(pkin(8));
t387 = t345 * t348;
t350 = cos(pkin(6));
t386 = t345 * t350;
t349 = cos(pkin(7));
t385 = t346 * t349;
t353 = sin(qJ(4));
t384 = t348 * t353;
t358 = cos(qJ(4));
t383 = t348 * t358;
t354 = sin(qJ(3));
t382 = t349 * t354;
t359 = cos(qJ(3));
t381 = t349 * t359;
t355 = sin(qJ(2));
t380 = t350 * t355;
t360 = cos(qJ(2));
t379 = t350 * t360;
t351 = sin(qJ(6));
t378 = t351 * t357;
t377 = t354 * t355;
t376 = t354 * t360;
t375 = t355 * t359;
t356 = cos(qJ(6));
t374 = t356 * t357;
t373 = t359 * t360;
t347 = cos(pkin(14));
t372 = t347 * t388;
t371 = t355 * t388;
t343 = sin(pkin(14));
t337 = -t343 * t355 + t347 * t379;
t338 = t343 * t360 + t347 * t380;
t324 = -t337 * t354 - t338 * t381;
t370 = t324 * t348 + t338 * t391;
t339 = -t343 * t379 - t347 * t355;
t340 = -t343 * t380 + t347 * t360;
t326 = -t339 * t354 - t340 * t381;
t369 = t326 * t348 + t340 * t391;
t368 = -t337 * t345 - t347 * t385;
t367 = -t339 * t345 + t343 * t385;
t366 = t339 * t349 + t343 * t388;
t365 = t350 * t349 - t360 * t388;
t335 = (-t349 * t375 - t376) * t346;
t364 = t335 * t348 + t344 * t371;
t363 = t368 * t344;
t362 = t367 * t344;
t361 = t365 * t344;
t336 = (-t349 * t377 + t373) * t346;
t333 = t354 * t386 + (t349 * t376 + t375) * t346;
t332 = t359 * t386 + (t349 * t373 - t377) * t346;
t328 = -t335 * t344 + t348 * t371;
t327 = t339 * t359 - t340 * t382;
t325 = t337 * t359 - t338 * t382;
t323 = t340 * t359 + t366 * t354;
t322 = -t340 * t354 + t366 * t359;
t321 = t337 * t382 + t338 * t359 - t354 * t372;
t320 = -t338 * t354 + (t337 * t349 - t372) * t359;
t319 = -t332 * t344 + t365 * t348;
t316 = t336 * t358 + t364 * t353;
t315 = t336 * t353 - t364 * t358;
t314 = -t326 * t344 + t340 * t387;
t313 = -t324 * t344 + t338 * t387;
t312 = t332 * t358 - t333 * t384;
t311 = t332 * t353 + t333 * t383;
t310 = -t322 * t344 + t367 * t348;
t309 = -t320 * t344 + t368 * t348;
t308 = t333 * t358 + (t332 * t348 + t361) * t353;
t307 = -t332 * t383 + t333 * t353 - t358 * t361;
t306 = t312 * t357 + t333 * t390;
t305 = t316 * t357 + t328 * t352;
t304 = t322 * t358 - t323 * t384;
t303 = t322 * t353 + t323 * t383;
t302 = t320 * t358 - t321 * t384;
t301 = t320 * t353 + t321 * t383;
t300 = t327 * t358 + t369 * t353;
t299 = t327 * t353 - t369 * t358;
t298 = t325 * t358 + t370 * t353;
t297 = t325 * t353 - t370 * t358;
t296 = t323 * t358 + (t322 * t348 + t362) * t353;
t295 = -t322 * t383 + t323 * t353 - t358 * t362;
t294 = t321 * t358 + (t320 * t348 + t363) * t353;
t293 = -t320 * t383 + t321 * t353 - t358 * t363;
t292 = t308 * t357 + t319 * t352;
t291 = -t308 * t352 + t319 * t357;
t290 = t304 * t357 + t323 * t390;
t289 = t302 * t357 + t321 * t390;
t288 = t300 * t357 + t314 * t352;
t287 = t298 * t357 + t313 * t352;
t286 = t296 * t357 + t310 * t352;
t285 = -t296 * t352 + t310 * t357;
t284 = t294 * t357 + t309 * t352;
t283 = -t294 * t352 + t309 * t357;
t1 = [0, t288 * t356 + t299 * t351, t290 * t356 + t303 * t351, -t295 * t374 + t296 * t351, t285 * t356, -t286 * t351 + t295 * t356; 0, t287 * t356 + t297 * t351, t289 * t356 + t301 * t351, -t293 * t374 + t294 * t351, t283 * t356, -t284 * t351 + t293 * t356; 0, t305 * t356 + t315 * t351, t306 * t356 + t311 * t351, -t307 * t374 + t308 * t351, t291 * t356, -t292 * t351 + t307 * t356; 0, -t288 * t351 + t299 * t356, -t290 * t351 + t303 * t356, t295 * t378 + t296 * t356, -t285 * t351, -t286 * t356 - t295 * t351; 0, -t287 * t351 + t297 * t356, -t289 * t351 + t301 * t356, t293 * t378 + t294 * t356, -t283 * t351, -t284 * t356 - t293 * t351; 0, -t305 * t351 + t315 * t356, -t306 * t351 + t311 * t356, t307 * t378 + t308 * t356, -t291 * t351, -t292 * t356 - t307 * t351; 0, t300 * t352 - t314 * t357, t304 * t352 - t323 * t389, -t295 * t352, t286, 0; 0, t298 * t352 - t313 * t357, t302 * t352 - t321 * t389, -t293 * t352, t284, 0; 0, t316 * t352 - t328 * t357, t312 * t352 - t333 * t389, -t307 * t352, t292, 0;];
JR_rot  = t1;
