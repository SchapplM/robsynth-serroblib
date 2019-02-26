% Rotatorische Teilmatrix der Rotationsmatrix-Jacobi-Matrix für Segment Nr. 6 (0=Basis) von
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

% Quelle: HybrDyn-Toolbox
% Datum: 2019-02-26 22:56
% Revision: d75aae1ac561373cd3be920984c3df29a1c2ecc8 (2019-02-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JR_rot = S6RRPRRR14_jacobiR_rot_6_sym_varpar(qJ, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(14,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPRRR14_jacobiR_rot_6_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [14 1]), ...
  'S6RRPRRR14_jacobiR_rot_6_sym_varpar: pkin has to be [14x1] (double)');

%% Symbolic Calculation
% From jacobiR_rot_6_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 22:55:47
% EndTime: 2019-02-26 22:55:48
% DurationCPUTime: 0.79s
% Computational Cost: add. (825->85), mult. (2425->180), div. (0->0), fcn. (3286->18), ass. (0->92)
t368 = sin(qJ(2));
t369 = sin(qJ(1));
t373 = cos(qJ(2));
t374 = cos(qJ(1));
t412 = cos(pkin(6));
t394 = t374 * t412;
t353 = t368 * t394 + t369 * t373;
t358 = sin(pkin(14));
t362 = cos(pkin(14));
t352 = t369 * t368 - t373 * t394;
t360 = sin(pkin(7));
t364 = cos(pkin(7));
t361 = sin(pkin(6));
t402 = t361 * t374;
t389 = t352 * t364 + t360 * t402;
t337 = -t353 * t362 + t389 * t358;
t367 = sin(qJ(4));
t372 = cos(qJ(4));
t348 = -t352 * t360 + t364 * t402;
t359 = sin(pkin(8));
t363 = cos(pkin(8));
t413 = t353 * t358 + t389 * t362;
t417 = t348 * t359 + t363 * t413;
t315 = t337 * t372 + t367 * t417;
t326 = -t348 * t363 + t359 * t413;
t366 = sin(qJ(5));
t371 = cos(qJ(5));
t305 = t315 * t371 - t326 * t366;
t314 = t337 * t367 - t372 * t417;
t365 = sin(qJ(6));
t370 = cos(qJ(6));
t423 = t305 * t365 - t314 * t370;
t422 = t305 * t370 + t314 * t365;
t303 = t315 * t366 + t326 * t371;
t393 = t412 * t360;
t399 = t364 * t373;
t377 = t362 * t393 + (-t358 * t368 + t362 * t399) * t361;
t385 = -t361 * t373 * t360 + t412 * t364;
t416 = t385 * t359 + t377 * t363;
t395 = t369 * t412;
t355 = -t368 * t395 + t374 * t373;
t354 = -t374 * t368 - t373 * t395;
t403 = t361 * t369;
t387 = t354 * t364 + t360 * t403;
t380 = t355 * t358 - t387 * t362;
t388 = -t354 * t360 + t364 * t403;
t415 = -t388 * t359 + t380 * t363;
t407 = t358 * t364;
t406 = t360 * t359;
t405 = t360 * t363;
t404 = t361 * t368;
t401 = t362 * t364;
t400 = t364 * t368;
t398 = t365 * t371;
t397 = t370 * t371;
t396 = t360 * t404;
t339 = t352 * t358 - t353 * t401;
t391 = t339 * t363 + t353 * t406;
t341 = -t354 * t358 - t355 * t401;
t390 = t341 * t363 + t355 * t406;
t350 = (-t358 * t373 - t362 * t400) * t361;
t386 = t350 * t363 + t359 * t396;
t375 = t380 * t359 + t388 * t363;
t351 = (-t358 * t400 + t362 * t373) * t361;
t347 = t362 * t404 + (t361 * t399 + t393) * t358;
t343 = -t350 * t359 + t363 * t396;
t342 = t354 * t362 - t355 * t407;
t340 = -t352 * t362 - t353 * t407;
t338 = t355 * t362 + t387 * t358;
t334 = -t377 * t359 + t385 * t363;
t331 = -t341 * t359 + t355 * t405;
t330 = -t339 * t359 + t353 * t405;
t329 = t351 * t372 + t386 * t367;
t328 = t351 * t367 - t386 * t372;
t324 = t347 * t372 + t416 * t367;
t323 = t347 * t367 - t416 * t372;
t322 = t329 * t371 + t343 * t366;
t321 = t342 * t372 + t390 * t367;
t320 = t342 * t367 - t390 * t372;
t319 = t340 * t372 + t391 * t367;
t318 = t340 * t367 - t391 * t372;
t317 = t338 * t372 - t415 * t367;
t316 = t338 * t367 + t415 * t372;
t311 = t324 * t371 + t334 * t366;
t310 = -t324 * t366 + t334 * t371;
t309 = t321 * t371 + t331 * t366;
t308 = t319 * t371 + t330 * t366;
t307 = t317 * t371 + t375 * t366;
t306 = t317 * t366 - t375 * t371;
t302 = t307 * t370 + t316 * t365;
t301 = -t307 * t365 + t316 * t370;
t1 = [t422, t309 * t370 + t320 * t365, 0, -t316 * t397 + t317 * t365, -t306 * t370, t301; t302, t308 * t370 + t318 * t365, 0, t314 * t397 - t315 * t365, t303 * t370, t423; 0, t322 * t370 + t328 * t365, 0, -t323 * t397 + t324 * t365, t310 * t370, -t311 * t365 + t323 * t370; -t423, -t309 * t365 + t320 * t370, 0, t316 * t398 + t317 * t370, t306 * t365, -t302; t301, -t308 * t365 + t318 * t370, 0, -t314 * t398 - t315 * t370, -t303 * t365, t422; 0, -t322 * t365 + t328 * t370, 0, t323 * t398 + t324 * t370, -t310 * t365, -t311 * t370 - t323 * t365; t303, t321 * t366 - t331 * t371, 0, -t316 * t366, t307, 0; t306, t319 * t366 - t330 * t371, 0, t314 * t366, -t305, 0; 0, t329 * t366 - t343 * t371, 0, -t323 * t366, t311, 0;];
JR_rot  = t1;
