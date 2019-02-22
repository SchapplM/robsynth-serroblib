% Rotatorische Teilmatrix der Rotationsmatrix-Jacobi-Matrix für Segment Nr. 6 (0=Basis) von
% S6RPRRRR12
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
%   pkin=[a2,a3,a4,a5,a6,alpha2,alpha3,alpha4,d1,d3,d4,d5,d6,theta2]';
%
% Output:
% JR_rot [9x6]
%   Jacobi-Matrix der Endeffektor-Rotationsmatrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-02-22 11:05
% Revision: 2b76964ad985d937eecd005a1a368749e6b3dc4d (2019-02-18)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JR_rot = S6RPRRRR12_jacobiR_rot_6_sym_varpar(qJ, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(14,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRRRR12_jacobiR_rot_6_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [14 1]), ...
  'S6RPRRRR12_jacobiR_rot_6_sym_varpar: pkin has to be [14x1] (double)');

%% Symbolic Calculation
% From jacobiR_rot_6_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-22 11:05:35
% EndTime: 2019-02-22 11:05:36
% DurationCPUTime: 0.75s
% Computational Cost: add. (824->80), mult. (2410->158), div. (0->0), fcn. (3270->18), ass. (0->82)
t344 = cos(pkin(6));
t341 = cos(pkin(14));
t354 = cos(qJ(1));
t365 = t354 * t341;
t337 = sin(pkin(14));
t349 = sin(qJ(1));
t369 = t349 * t337;
t331 = -t344 * t365 + t369;
t366 = t354 * t337;
t368 = t349 * t341;
t332 = t344 * t366 + t368;
t348 = sin(qJ(3));
t353 = cos(qJ(3));
t339 = sin(pkin(7));
t340 = sin(pkin(6));
t374 = t340 * t354;
t364 = t339 * t374;
t343 = cos(pkin(7));
t371 = t343 * t348;
t321 = t331 * t371 - t332 * t353 + t348 * t364;
t347 = sin(qJ(4));
t352 = cos(qJ(4));
t320 = (t331 * t343 + t364) * t353 + t332 * t348;
t342 = cos(pkin(8));
t329 = -t331 * t339 + t343 * t374;
t338 = sin(pkin(8));
t381 = t329 * t338;
t363 = t320 * t342 + t381;
t300 = t321 * t352 + t363 * t347;
t311 = t320 * t338 - t329 * t342;
t346 = sin(qJ(5));
t351 = cos(qJ(5));
t290 = t300 * t351 - t311 * t346;
t345 = sin(qJ(6));
t388 = t290 * t345;
t350 = cos(qJ(6));
t387 = t290 * t350;
t288 = t300 * t346 + t311 * t351;
t383 = t321 * t347;
t378 = t338 * t346;
t377 = t338 * t351;
t376 = t339 * t344;
t375 = t340 * t349;
t373 = t342 * t347;
t372 = t342 * t352;
t370 = t345 * t351;
t367 = t350 * t351;
t333 = -t344 * t368 - t366;
t361 = -t333 * t339 + t343 * t375;
t360 = t333 * t343 + t339 * t375;
t359 = -t340 * t341 * t339 + t344 * t343;
t357 = t361 * t338;
t356 = t359 * t338;
t334 = -t344 * t369 + t365;
t322 = -t334 * t348 + t360 * t353;
t355 = -t322 * t338 + t361 * t342;
t328 = t348 * t376 + (t337 * t353 + t341 * t371) * t340;
t327 = t353 * t376 + (t341 * t343 * t353 - t337 * t348) * t340;
t323 = t334 * t353 + t360 * t348;
t317 = -t327 * t338 + t359 * t342;
t314 = t327 * t352 - t328 * t373;
t313 = t327 * t347 + t328 * t372;
t309 = t328 * t352 + (t327 * t342 + t356) * t347;
t308 = -t327 * t372 + t328 * t347 - t352 * t356;
t307 = t314 * t351 + t328 * t378;
t306 = t322 * t352 - t323 * t373;
t305 = t322 * t347 + t323 * t372;
t304 = -t320 * t352 + t321 * t373;
t303 = -t320 * t347 - t321 * t372;
t302 = t323 * t352 + (t322 * t342 + t357) * t347;
t301 = -t322 * t372 + t323 * t347 - t352 * t357;
t299 = -t363 * t352 + t383;
t297 = t320 * t372 + t352 * t381 - t383;
t296 = t309 * t351 + t317 * t346;
t295 = -t309 * t346 + t317 * t351;
t294 = t306 * t351 + t323 * t378;
t293 = t304 * t351 - t321 * t378;
t292 = t302 * t351 + t355 * t346;
t291 = t302 * t346 - t355 * t351;
t287 = t292 * t350 + t301 * t345;
t286 = -t292 * t345 + t301 * t350;
t1 = [t299 * t345 + t387, 0, t294 * t350 + t305 * t345, -t301 * t367 + t302 * t345, -t291 * t350, t286; t287, 0, t293 * t350 + t303 * t345, -t297 * t367 - t300 * t345, t288 * t350, t297 * t350 + t388; 0, 0, t307 * t350 + t313 * t345, -t308 * t367 + t309 * t345, t295 * t350, -t296 * t345 + t308 * t350; t299 * t350 - t388, 0, -t294 * t345 + t305 * t350, t301 * t370 + t302 * t350, t291 * t345, -t287; t286, 0, -t293 * t345 + t303 * t350, t297 * t370 - t300 * t350, -t288 * t345, -t297 * t345 + t387; 0, 0, -t307 * t345 + t313 * t350, t308 * t370 + t309 * t350, -t295 * t345, -t296 * t350 - t308 * t345; t288, 0, t306 * t346 - t323 * t377, -t301 * t346, t292, 0; t291, 0, t304 * t346 + t321 * t377, -t297 * t346, -t290, 0; 0, 0, t314 * t346 - t328 * t377, -t308 * t346, t296, 0;];
JR_rot  = t1;
