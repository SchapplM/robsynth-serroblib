% Rotatorische Teilmatrix der Rotationsmatrix-Jacobi-Matrix für Segment Nr. 6 (0=Basis) von
% S6RRRRRR10
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
%   pkin=[a2,a3,a4,a5,a6,alpha2,alpha3,alpha4,d1,d2,d3,d4,d5,d6]';
%
% Output:
% JR_rot [9x6]
%   Jacobi-Matrix der Endeffektor-Rotationsmatrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-02-22 12:39
% Revision: 2b76964ad985d937eecd005a1a368749e6b3dc4d (2019-02-18)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JR_rot = S6RRRRRR10_jacobiR_rot_6_sym_varpar(qJ, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(14,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRRRR10_jacobiR_rot_6_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [14 1]), ...
  'S6RRRRRR10_jacobiR_rot_6_sym_varpar: pkin has to be [14x1] (double)');

%% Symbolic Calculation
% From jacobiR_rot_6_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-22 12:39:02
% EndTime: 2019-02-22 12:39:03
% DurationCPUTime: 1.16s
% Computational Cost: add. (985->110), mult. (2912->225), div. (0->0), fcn. (3942->18), ass. (0->110)
t415 = sin(qJ(2));
t416 = sin(qJ(1));
t421 = cos(qJ(2));
t422 = cos(qJ(1));
t460 = cos(pkin(6));
t436 = t422 * t460;
t399 = t416 * t415 - t421 * t436;
t400 = t415 * t436 + t416 * t421;
t414 = sin(qJ(3));
t420 = cos(qJ(3));
t407 = sin(pkin(7));
t408 = sin(pkin(6));
t450 = t408 * t422;
t438 = t407 * t450;
t410 = cos(pkin(7));
t447 = t410 * t414;
t382 = t399 * t447 - t400 * t420 + t414 * t438;
t413 = sin(qJ(4));
t419 = cos(qJ(4));
t381 = (t399 * t410 + t438) * t420 + t400 * t414;
t409 = cos(pkin(8));
t395 = -t399 * t407 + t410 * t450;
t406 = sin(pkin(8));
t459 = t395 * t406;
t434 = t381 * t409 + t459;
t352 = t382 * t419 + t434 * t413;
t368 = t381 * t406 - t395 * t409;
t412 = sin(qJ(5));
t418 = cos(qJ(5));
t340 = t352 * t418 - t368 * t412;
t411 = sin(qJ(6));
t467 = t340 * t411;
t417 = cos(qJ(6));
t466 = t340 * t417;
t338 = t352 * t412 + t368 * t418;
t462 = t382 * t413;
t456 = t406 * t407;
t455 = t406 * t412;
t454 = t406 * t418;
t453 = t407 * t408;
t452 = t407 * t409;
t451 = t408 * t416;
t449 = t409 * t413;
t448 = t409 * t419;
t446 = t410 * t420;
t445 = t411 * t418;
t444 = t414 * t415;
t443 = t414 * t421;
t442 = t415 * t420;
t441 = t417 * t418;
t440 = t420 * t421;
t439 = t415 * t453;
t437 = t416 * t460;
t435 = t460 * t407;
t385 = t399 * t414 - t400 * t446;
t433 = t385 * t409 + t400 * t456;
t401 = -t422 * t415 - t421 * t437;
t402 = -t415 * t437 + t422 * t421;
t387 = -t401 * t414 - t402 * t446;
t432 = t387 * t409 + t402 * t456;
t430 = -t401 * t407 + t410 * t451;
t429 = t401 * t410 + t407 * t451;
t397 = (-t410 * t442 - t443) * t408;
t428 = t397 * t409 + t406 * t439;
t427 = t410 * t460 - t421 * t453;
t425 = t430 * t406;
t424 = t427 * t406;
t383 = -t402 * t414 + t429 * t420;
t423 = -t383 * t406 + t430 * t409;
t398 = (-t410 * t444 + t440) * t408;
t394 = t414 * t435 + (t410 * t443 + t442) * t408;
t393 = t420 * t435 + (t410 * t440 - t444) * t408;
t389 = -t397 * t406 + t409 * t439;
t388 = t401 * t420 - t402 * t447;
t386 = -t399 * t420 - t400 * t447;
t384 = t402 * t420 + t429 * t414;
t378 = -t393 * t406 + t427 * t409;
t375 = -t387 * t406 + t402 * t452;
t374 = -t385 * t406 + t400 * t452;
t373 = t398 * t419 + t413 * t428;
t372 = t398 * t413 - t419 * t428;
t371 = t393 * t419 - t394 * t449;
t370 = t393 * t413 + t394 * t448;
t366 = t394 * t419 + (t393 * t409 + t424) * t413;
t365 = -t393 * t448 + t394 * t413 - t419 * t424;
t364 = t371 * t418 + t394 * t455;
t363 = t373 * t418 + t389 * t412;
t362 = t383 * t419 - t384 * t449;
t361 = t383 * t413 + t384 * t448;
t360 = -t381 * t419 + t382 * t449;
t359 = -t381 * t413 - t382 * t448;
t358 = t388 * t419 + t413 * t432;
t357 = t388 * t413 - t419 * t432;
t356 = t386 * t419 + t413 * t433;
t355 = t386 * t413 - t419 * t433;
t354 = t384 * t419 + (t383 * t409 + t425) * t413;
t353 = -t383 * t448 + t384 * t413 - t419 * t425;
t351 = -t419 * t434 + t462;
t349 = t381 * t448 + t419 * t459 - t462;
t348 = t366 * t418 + t378 * t412;
t347 = -t366 * t412 + t378 * t418;
t346 = t362 * t418 + t384 * t455;
t345 = t360 * t418 - t382 * t455;
t344 = t358 * t418 + t375 * t412;
t343 = t356 * t418 + t374 * t412;
t342 = t354 * t418 + t412 * t423;
t341 = t354 * t412 - t418 * t423;
t337 = t342 * t417 + t353 * t411;
t336 = -t342 * t411 + t353 * t417;
t1 = [t351 * t411 + t466, t344 * t417 + t357 * t411, t346 * t417 + t361 * t411, -t353 * t441 + t354 * t411, -t341 * t417, t336; t337, t343 * t417 + t355 * t411, t345 * t417 + t359 * t411, -t349 * t441 - t352 * t411, t338 * t417, t349 * t417 + t467; 0, t363 * t417 + t372 * t411, t364 * t417 + t370 * t411, -t365 * t441 + t366 * t411, t347 * t417, -t348 * t411 + t365 * t417; t351 * t417 - t467, -t344 * t411 + t357 * t417, -t346 * t411 + t361 * t417, t353 * t445 + t354 * t417, t341 * t411, -t337; t336, -t343 * t411 + t355 * t417, -t345 * t411 + t359 * t417, t349 * t445 - t352 * t417, -t338 * t411, -t349 * t411 + t466; 0, -t363 * t411 + t372 * t417, -t364 * t411 + t370 * t417, t365 * t445 + t366 * t417, -t347 * t411, -t348 * t417 - t365 * t411; t338, t358 * t412 - t375 * t418, t362 * t412 - t384 * t454, -t353 * t412, t342, 0; t341, t356 * t412 - t374 * t418, t360 * t412 + t382 * t454, -t349 * t412, -t340, 0; 0, t373 * t412 - t389 * t418, t371 * t412 - t394 * t454, -t365 * t412, t348, 0;];
JR_rot  = t1;
