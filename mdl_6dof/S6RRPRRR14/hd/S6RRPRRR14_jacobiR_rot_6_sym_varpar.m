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

% Quelle: HybrDyn-Toolbox (ehem. IRT-Maple-Toolbox)
% Datum: 2018-12-10 18:38
% Revision: bb42a8b95257d9bc83910d26e849f5825122f662 (2018-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für mechatronische Systeme, Universität Hannover

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
% StartTime: 2018-12-10 18:38:24
% EndTime: 2018-12-10 18:38:25
% DurationCPUTime: 0.91s
% Computational Cost: add. (3171->112), mult. (3207->186), div. (0->0), fcn. (3286->30), ass. (0->106)
t456 = sin(qJ(2));
t457 = sin(qJ(1));
t462 = cos(qJ(1));
t481 = pkin(6) + qJ(2);
t474 = cos(t481) / 0.2e1;
t482 = pkin(6) - qJ(2);
t478 = cos(t482);
t463 = t478 / 0.2e1 + t474;
t416 = t457 * t456 - t462 * t463;
t472 = sin(t481) / 0.2e1;
t476 = sin(t482);
t429 = t472 - t476 / 0.2e1;
t461 = cos(qJ(2));
t417 = t462 * t429 + t457 * t461;
t443 = pkin(7) + pkin(14);
t432 = sin(t443) / 0.2e1;
t444 = pkin(7) - pkin(14);
t441 = sin(t444);
t423 = t432 + t441 / 0.2e1;
t433 = cos(t444) / 0.2e1;
t442 = cos(t443);
t425 = t433 + t442 / 0.2e1;
t445 = sin(pkin(14));
t448 = sin(pkin(6));
t485 = t448 * t462;
t394 = t416 * t425 + t417 * t445 + t423 * t485;
t424 = t432 - t441 / 0.2e1;
t426 = t433 - t442 / 0.2e1;
t449 = cos(pkin(14));
t395 = t416 * t424 - t417 * t449 + t426 * t485;
t447 = sin(pkin(7));
t451 = cos(pkin(7));
t411 = -t416 * t447 + t451 * t485;
t479 = pkin(8) + qJ(4);
t471 = sin(t479) / 0.2e1;
t480 = pkin(8) - qJ(4);
t475 = sin(t480);
t427 = t471 - t475 / 0.2e1;
t473 = cos(t479) / 0.2e1;
t477 = cos(t480);
t430 = t473 - t477 / 0.2e1;
t460 = cos(qJ(4));
t368 = t394 * t427 + t395 * t460 - t411 * t430;
t446 = sin(pkin(8));
t450 = cos(pkin(8));
t383 = t394 * t446 - t411 * t450;
t454 = sin(qJ(5));
t459 = cos(qJ(5));
t357 = t368 * t459 - t383 * t454;
t455 = sin(qJ(4));
t468 = t471 + t475 / 0.2e1;
t469 = t477 / 0.2e1 + t473;
t364 = t394 * t469 - t395 * t455 + t411 * t468;
t453 = sin(qJ(6));
t458 = cos(qJ(6));
t492 = t357 * t453 + t364 * t458;
t491 = t357 * t458 - t364 * t453;
t355 = t368 * t454 + t383 * t459;
t488 = t430 * t447;
t487 = t447 * t450;
t486 = t448 * t457;
t484 = t453 * t459;
t483 = t458 * t459;
t419 = -t462 * t456 - t457 * t463;
t421 = t457 * t429 - t462 * t461;
t396 = t419 * t425 + t421 * t445 + t423 * t486;
t413 = -t419 * t447 + t451 * t486;
t470 = -t396 * t446 + t413 * t450;
t397 = t419 * t424 - t421 * t449 + t426 * t486;
t466 = t396 * t427 + t397 * t460 - t413 * t430;
t428 = t472 + t476 / 0.2e1;
t431 = t474 - t478 / 0.2e1;
t452 = cos(pkin(6));
t405 = t452 * t423 + t428 * t425 + t431 * t445;
t406 = t428 * t424 + t452 * t426 - t431 * t449;
t415 = -t428 * t447 + t452 * t451;
t465 = t405 * t427 + t406 * t460 - t415 * t430;
t464 = t447 * t468;
t410 = t431 * t424 + t428 * t449;
t409 = t431 * t425 - t428 * t445;
t403 = t419 * t449 + t421 * t424;
t402 = -t419 * t445 + t421 * t425;
t401 = -t416 * t449 - t417 * t424;
t400 = t416 * t445 - t417 * t425;
t399 = -t409 * t446 - t431 * t487;
t389 = -t405 * t446 + t415 * t450;
t386 = -t402 * t446 - t421 * t487;
t385 = -t400 * t446 + t417 * t487;
t382 = t409 * t427 + t410 * t460 + t431 * t488;
t381 = -t409 * t469 + t410 * t455 + t431 * t464;
t377 = -t405 * t469 + t406 * t455 - t415 * t468;
t376 = t402 * t427 + t403 * t460 + t421 * t488;
t375 = -t402 * t469 + t403 * t455 + t421 * t464;
t374 = t400 * t427 + t401 * t460 - t417 * t488;
t373 = -t400 * t469 + t401 * t455 - t417 * t464;
t372 = t382 * t459 + t399 * t454;
t369 = -t396 * t469 + t397 * t455 - t413 * t468;
t363 = t389 * t454 + t459 * t465;
t362 = t389 * t459 - t454 * t465;
t361 = t376 * t459 + t386 * t454;
t360 = t374 * t459 + t385 * t454;
t359 = t470 * t454 + t459 * t466;
t358 = t454 * t466 - t470 * t459;
t354 = t359 * t458 + t369 * t453;
t353 = -t359 * t453 + t369 * t458;
t1 = [t491, t361 * t458 + t375 * t453, 0, -t369 * t483 + t453 * t466, -t358 * t458, t353; t354, t360 * t458 + t373 * t453, 0, -t364 * t483 - t368 * t453, t355 * t458, t492; 0, t372 * t458 + t381 * t453, 0, -t377 * t483 + t453 * t465, t362 * t458, -t363 * t453 + t377 * t458; -t492, -t361 * t453 + t375 * t458, 0, t369 * t484 + t458 * t466, t358 * t453, -t354; t353, -t360 * t453 + t373 * t458, 0, t364 * t484 - t368 * t458, -t355 * t453, t491; 0, -t372 * t453 + t381 * t458, 0, t377 * t484 + t458 * t465, -t362 * t453, -t363 * t458 - t377 * t453; t355, t376 * t454 - t386 * t459, 0, -t369 * t454, t359, 0; t358, t374 * t454 - t385 * t459, 0, -t364 * t454, -t357, 0; 0, t382 * t454 - t399 * t459, 0, -t377 * t454, t363, 0;];
JR_rot  = t1;
