% Rotatorische Teilmatrix der Rotationsmatrix-Jacobi-Matrix für Segment Nr. 7 (0=Basis) von
% S7RRRRRRR1
% Use Code from Maple symbolic Code Generation
%
% Rotationsmatrix-Jacobi-Matrix: Differentieller Zusammenhang zwischen
% gestapelter Endeffektor-Rotationsmatrix und verallgemeinerten Koordinaten.
%
%
% Input:
% qJ [7x1]
%   Generalized joint coordinates (joint angles)
% pkin [4x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[d1,d3,d5,d7]';
%
% Output:
% JR_rot [9x7]
%   Jacobi-Matrix der Endeffektor-Rotationsmatrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-02-26 22:54
% Revision: d75aae1ac561373cd3be920984c3df29a1c2ecc8 (2019-02-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JR_rot = S7RRRRRRR1_jacobiR_rot_7_sym_varpar(qJ, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(7,1),zeros(4,1)}
assert(isreal(qJ) && all(size(qJ) == [7 1]), ...
  'S7RRRRRRR1_jacobiR_rot_7_sym_varpar: qJ has to be [7x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [4 1]), ...
  'S7RRRRRRR1_jacobiR_rot_7_sym_varpar: pkin has to be [4x1] (double)');

%% Symbolic Calculation
% From jacobiR_rot_7_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 22:54:31
% EndTime: 2019-02-26 22:54:32
% DurationCPUTime: 0.78s
% Computational Cost: add. (405->92), mult. (1175->199), div. (0->0), fcn. (1631->14), ass. (0->92)
t343 = cos(qJ(3));
t336 = sin(qJ(3));
t345 = cos(qJ(1));
t350 = t345 * t336;
t338 = sin(qJ(1));
t344 = cos(qJ(2));
t356 = t338 * t344;
t323 = t343 * t356 + t350;
t342 = cos(qJ(4));
t335 = sin(qJ(4));
t337 = sin(qJ(2));
t360 = t337 * t335;
t305 = t323 * t342 + t338 * t360;
t349 = t345 * t343;
t322 = t336 * t356 - t349;
t334 = sin(qJ(5));
t341 = cos(qJ(5));
t286 = t305 * t341 - t322 * t334;
t359 = t337 * t342;
t304 = t323 * t335 - t338 * t359;
t333 = sin(qJ(6));
t340 = cos(qJ(6));
t273 = t286 * t340 + t304 * t333;
t285 = t305 * t334 + t322 * t341;
t332 = sin(qJ(7));
t339 = cos(qJ(7));
t378 = t273 * t332 + t285 * t339;
t377 = -t273 * t339 + t285 * t332;
t374 = t286 * t333 - t304 * t340;
t369 = t332 * t334;
t368 = t332 * t340;
t367 = t333 * t335;
t366 = t333 * t341;
t365 = t334 * t339;
t364 = t334 * t342;
t363 = t335 * t340;
t362 = t336 * t337;
t361 = t336 * t344;
t358 = t337 * t343;
t357 = t337 * t345;
t355 = t339 * t340;
t354 = t340 * t341;
t353 = t341 * t342;
t352 = t344 * t335;
t351 = t344 * t342;
t348 = t336 * t360;
t347 = t334 * t362;
t346 = t341 * t362;
t321 = t342 * t358 - t352;
t320 = t335 * t358 + t351;
t327 = -t338 * t336 + t344 * t349;
t326 = -t338 * t343 - t344 * t350;
t325 = t343 * t351 + t360;
t324 = t343 * t352 - t359;
t317 = t321 * t345;
t316 = t320 * t345;
t315 = t321 * t338;
t314 = t320 * t338;
t313 = (-t334 * t343 - t336 * t353) * t337;
t312 = (-t336 * t364 + t341 * t343) * t337;
t311 = t327 * t342 + t335 * t357;
t310 = t327 * t335 - t342 * t357;
t309 = t325 * t341 - t334 * t361;
t308 = t325 * t334 + t341 * t361;
t303 = t321 * t341 - t347;
t302 = t321 * t334 + t346;
t301 = -t317 * t341 + t345 * t347;
t300 = -t317 * t334 - t345 * t346;
t299 = -t315 * t341 + t338 * t347;
t298 = -t315 * t334 - t338 * t346;
t297 = t313 * t340 - t333 * t348;
t296 = t326 * t353 - t327 * t334;
t295 = t326 * t364 + t327 * t341;
t294 = -t322 * t353 - t323 * t334;
t293 = -t322 * t364 + t323 * t341;
t292 = -t320 * t354 + t321 * t333;
t291 = t311 * t341 + t326 * t334;
t290 = t311 * t334 - t326 * t341;
t289 = t309 * t340 + t324 * t333;
t284 = t303 * t340 + t320 * t333;
t283 = -t303 * t333 + t320 * t340;
t282 = t301 * t340 - t316 * t333;
t281 = t299 * t340 - t314 * t333;
t280 = t296 * t340 + t326 * t367;
t279 = t294 * t340 - t322 * t367;
t278 = -t310 * t354 + t311 * t333;
t277 = -t304 * t354 + t305 * t333;
t276 = t291 * t340 + t310 * t333;
t275 = -t291 * t333 + t310 * t340;
t271 = t276 * t339 - t290 * t332;
t270 = -t276 * t332 - t290 * t339;
t1 = [t377, t282 * t339 - t300 * t332, t280 * t339 - t295 * t332, t278 * t339 + t310 * t369, -t290 * t355 - t291 * t332, t275 * t339, t270; t271, t281 * t339 - t298 * t332, t279 * t339 - t293 * t332, t277 * t339 + t304 * t369, -t285 * t355 - t286 * t332, -t374 * t339, -t378; 0, t289 * t339 - t308 * t332, t297 * t339 - t312 * t332, t292 * t339 + t320 * t369, -t302 * t355 - t303 * t332, t283 * t339, -t284 * t332 - t302 * t339; t378, -t282 * t332 - t300 * t339, -t280 * t332 - t295 * t339, -t278 * t332 + t310 * t365, t290 * t368 - t291 * t339, -t275 * t332, -t271; t270, -t281 * t332 - t298 * t339, -t279 * t332 - t293 * t339, -t277 * t332 + t304 * t365, t285 * t368 - t286 * t339, t374 * t332, t377; 0, -t289 * t332 - t308 * t339, -t297 * t332 - t312 * t339, -t292 * t332 + t320 * t365, t302 * t368 - t303 * t339, -t283 * t332, -t284 * t339 + t302 * t332; t374, -t301 * t333 - t316 * t340, -t296 * t333 + t326 * t363, t310 * t366 + t311 * t340, t290 * t333, -t276, 0; t275, -t299 * t333 - t314 * t340, -t294 * t333 - t322 * t363, t304 * t366 + t305 * t340, t285 * t333, -t273, 0; 0, -t309 * t333 + t324 * t340, -t313 * t333 - t340 * t348, t320 * t366 + t321 * t340, t302 * t333, -t284, 0;];
JR_rot  = t1;
