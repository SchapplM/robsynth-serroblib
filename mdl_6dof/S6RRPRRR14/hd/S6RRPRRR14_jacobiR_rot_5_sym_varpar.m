% Rotatorische Teilmatrix der Rotationsmatrix-Jacobi-Matrix für Segment Nr. 5 (0=Basis) von
% S6RRPRRR14
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
%   pkin=[a2,a3,a4,a5,a6,alpha2,alpha3,alpha4,d1,d2,d4,d5,d6,theta3]';
%
% Output:
% JR_rot [9x6]
%   Jacobi-Matrix der Endeffektor-Rotationsmatrix

% Quelle: HybrDyn-Toolbox (ehem. IRT-Maple-Toolbox)
% Datum: 2018-12-10 18:38
% Revision: bb42a8b95257d9bc83910d26e849f5825122f662 (2018-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für mechatronische Systeme, Universität Hannover

function JR_rot = S6RRPRRR14_jacobiR_rot_5_sym_varpar(qJ, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(14,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPRRR14_jacobiR_rot_5_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [14 1]), ...
  'S6RRPRRR14_jacobiR_rot_5_sym_varpar: pkin has to be [14x1] (double)');

%% Symbolic Calculation
% From jacobiR_rot_5_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2018-12-10 18:38:22
% EndTime: 2018-12-10 18:38:23
% DurationCPUTime: 0.50s
% Computational Cost: add. (1558->91), mult. (1569->146), div. (0->0), fcn. (1610->28), ass. (0->91)
t374 = sin(qJ(2));
t375 = sin(qJ(1));
t379 = cos(qJ(1));
t391 = pkin(6) + qJ(2);
t384 = cos(t391) / 0.2e1;
t392 = pkin(6) - qJ(2);
t388 = cos(t392);
t380 = t388 / 0.2e1 + t384;
t333 = t374 * t375 - t379 * t380;
t382 = sin(t391) / 0.2e1;
t386 = sin(t392);
t347 = t382 - t386 / 0.2e1;
t378 = cos(qJ(2));
t334 = t347 * t379 + t375 * t378;
t362 = pkin(7) + pkin(14);
t351 = sin(t362) / 0.2e1;
t363 = pkin(7) - pkin(14);
t360 = sin(t363);
t340 = t351 + t360 / 0.2e1;
t352 = cos(t363) / 0.2e1;
t361 = cos(t362);
t342 = t352 + t361 / 0.2e1;
t364 = sin(pkin(14));
t367 = sin(pkin(6));
t393 = t367 * t379;
t314 = t333 * t342 + t334 * t364 + t340 * t393;
t341 = t351 - t360 / 0.2e1;
t343 = t352 - t361 / 0.2e1;
t368 = cos(pkin(14));
t315 = t333 * t341 - t334 * t368 + t343 * t393;
t366 = sin(pkin(7));
t370 = cos(pkin(7));
t329 = -t333 * t366 + t370 * t393;
t389 = pkin(8) + qJ(4);
t381 = sin(t389) / 0.2e1;
t390 = pkin(8) - qJ(4);
t385 = sin(t390);
t345 = t381 - t385 / 0.2e1;
t383 = cos(t389) / 0.2e1;
t387 = cos(t390);
t348 = t383 - t387 / 0.2e1;
t377 = cos(qJ(4));
t295 = t314 * t345 + t315 * t377 - t329 * t348;
t365 = sin(pkin(8));
t369 = cos(pkin(8));
t303 = t314 * t365 - t329 * t369;
t372 = sin(qJ(5));
t376 = cos(qJ(5));
t400 = t295 * t376 - t303 * t372;
t399 = t295 * t372 + t303 * t376;
t344 = t381 + t385 / 0.2e1;
t349 = t387 / 0.2e1 + t383;
t373 = sin(qJ(4));
t293 = -t314 * t349 + t315 * t373 - t329 * t344;
t398 = t344 * t366;
t397 = t348 * t366;
t350 = t384 - t388 / 0.2e1;
t396 = t350 * t366;
t395 = t366 * t369;
t394 = t367 * t375;
t338 = t347 * t375 - t378 * t379;
t336 = -t374 * t379 - t375 * t380;
t316 = t336 * t342 + t338 * t364 + t340 * t394;
t317 = t336 * t341 - t338 * t368 + t343 * t394;
t331 = -t336 * t366 + t370 * t394;
t297 = t316 * t345 + t317 * t377 - t331 * t348;
t346 = t382 + t386 / 0.2e1;
t371 = cos(pkin(6));
t324 = t340 * t371 + t342 * t346 + t350 * t364;
t325 = t341 * t346 + t343 * t371 - t350 * t368;
t332 = -t346 * t366 + t370 * t371;
t301 = t324 * t345 + t325 * t377 - t332 * t348;
t328 = t341 * t350 + t346 * t368;
t327 = t342 * t350 - t346 * t364;
t322 = t336 * t368 + t338 * t341;
t321 = -t336 * t364 + t338 * t342;
t320 = -t333 * t368 - t334 * t341;
t319 = t333 * t364 - t334 * t342;
t318 = -t327 * t365 - t350 * t395;
t309 = -t324 * t365 + t332 * t369;
t307 = -t321 * t365 - t338 * t395;
t306 = -t319 * t365 + t334 * t395;
t305 = -t316 * t365 + t331 * t369;
t302 = t327 * t345 + t328 * t377 + t348 * t396;
t300 = t324 * t349 - t325 * t373 + t332 * t344;
t299 = t321 * t345 + t322 * t377 + t338 * t397;
t298 = t319 * t345 + t320 * t377 - t334 * t397;
t296 = -t316 * t349 + t317 * t373 - t331 * t344;
t292 = t297 * t376 + t305 * t372;
t291 = -t297 * t372 + t305 * t376;
t1 = [t400, t299 * t376 + t307 * t372, 0, -t296 * t376, t291, 0; t292, t298 * t376 + t306 * t372, 0, t293 * t376, t399, 0; 0, t302 * t376 + t318 * t372, 0, t300 * t376, -t301 * t372 + t309 * t376, 0; -t399, -t299 * t372 + t307 * t376, 0, t296 * t372, -t292, 0; t291, -t298 * t372 + t306 * t376, 0, -t293 * t372, t400, 0; 0, -t302 * t372 + t318 * t376, 0, -t300 * t372, -t301 * t376 - t309 * t372, 0; t293, -t321 * t349 + t322 * t373 + t338 * t398, 0, t297, 0, 0; t296, -t319 * t349 + t320 * t373 - t334 * t398, 0, -t295, 0, 0; 0, -t327 * t349 + t328 * t373 + t344 * t396, 0, t301, 0, 0;];
JR_rot  = t1;
