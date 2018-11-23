% Rotatorische Teilmatrix der Rotationsmatrix-Jacobi-Matrix für Segment Nr. 5 (0=Basis) von
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

% Quelle: HybrDyn-Toolbox (ehem. IRT-Maple-Toolbox)
% Datum: 2018-11-23 11:28
% Revision: 76f9d5e39f14dc242b53c0d9d3d9db48bd8f37c0
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für mechatronische Systeme, Universität Hannover

function JR_rot = S6RRRRRR10_jacobiR_rot_5_sym_varpar(qJ, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(14,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRRRR10_jacobiR_rot_5_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [14 1]), ...
  'S6RRRRRR10_jacobiR_rot_5_sym_varpar: pkin has to be [14x1] (double)');

%% Symbolic Calculation
% From jacobiR_rot_5_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2018-11-23 11:27:19
% EndTime: 2018-11-23 11:27:20
% DurationCPUTime: 0.60s
% Computational Cost: add. (1968->102), mult. (1967->172), div. (0->0), fcn. (2007->28), ass. (0->96)
t403 = sin(qJ(2));
t404 = sin(qJ(1));
t409 = cos(qJ(1));
t427 = pkin(6) + qJ(2);
t416 = cos(t427) / 0.2e1;
t428 = pkin(6) - qJ(2);
t422 = cos(t428);
t410 = t422 / 0.2e1 + t416;
t365 = t404 * t403 - t409 * t410;
t413 = sin(t427) / 0.2e1;
t419 = sin(t428);
t377 = t413 - t419 / 0.2e1;
t408 = cos(qJ(2));
t366 = t409 * t377 + t404 * t408;
t425 = pkin(7) + qJ(3);
t412 = sin(t425) / 0.2e1;
t426 = pkin(7) - qJ(3);
t418 = sin(t426);
t375 = t412 - t418 / 0.2e1;
t415 = cos(t425) / 0.2e1;
t421 = cos(t426);
t380 = t415 - t421 / 0.2e1;
t407 = cos(qJ(3));
t396 = sin(pkin(6));
t429 = t396 * t409;
t342 = -t365 * t375 + t366 * t407 + t380 * t429;
t374 = t412 + t418 / 0.2e1;
t381 = t421 / 0.2e1 + t415;
t402 = sin(qJ(3));
t344 = t365 * t381 + t366 * t402 + t374 * t429;
t395 = sin(pkin(7));
t398 = cos(pkin(7));
t361 = -t365 * t395 + t398 * t429;
t423 = pkin(8) + qJ(4);
t411 = sin(t423) / 0.2e1;
t424 = pkin(8) - qJ(4);
t417 = sin(t424);
t373 = t411 - t417 / 0.2e1;
t414 = cos(t423) / 0.2e1;
t420 = cos(t424);
t378 = t414 - t420 / 0.2e1;
t406 = cos(qJ(4));
t320 = t342 * t406 - t344 * t373 + t361 * t378;
t394 = sin(pkin(8));
t397 = cos(pkin(8));
t332 = t344 * t394 - t361 * t397;
t400 = sin(qJ(5));
t405 = cos(qJ(5));
t438 = -t320 * t405 - t332 * t400;
t437 = t320 * t400 - t332 * t405;
t372 = t411 + t417 / 0.2e1;
t379 = t420 / 0.2e1 + t414;
t401 = sin(qJ(4));
t319 = -t342 * t401 - t344 * t379 - t361 * t372;
t436 = t372 * t395;
t435 = t378 * t395;
t382 = t416 - t422 / 0.2e1;
t434 = t382 * t395;
t433 = t394 * t400;
t432 = t394 * t405;
t431 = t395 * t397;
t430 = t396 * t404;
t370 = t404 * t377 - t409 * t408;
t368 = -t409 * t403 - t404 * t410;
t346 = t368 * t381 + t370 * t402 + t374 * t430;
t348 = -t368 * t375 + t370 * t407 + t380 * t430;
t363 = -t368 * t395 + t398 * t430;
t323 = t346 * t373 - t348 * t406 - t363 * t378;
t376 = t413 + t419 / 0.2e1;
t399 = cos(pkin(6));
t355 = t399 * t374 + t376 * t381 + t382 * t402;
t356 = t376 * t375 - t399 * t380 - t382 * t407;
t364 = -t376 * t395 + t399 * t398;
t329 = t355 * t373 + t356 * t406 - t364 * t378;
t360 = t382 * t375 + t376 * t407;
t359 = -t376 * t402 + t382 * t381;
t353 = t368 * t407 + t370 * t375;
t352 = -t368 * t402 + t370 * t381;
t351 = -t365 * t407 - t366 * t375;
t350 = t365 * t402 - t366 * t381;
t349 = -t359 * t394 - t382 * t431;
t338 = -t355 * t394 + t364 * t397;
t336 = -t352 * t394 - t370 * t431;
t335 = -t350 * t394 + t366 * t431;
t334 = -t346 * t394 + t363 * t397;
t331 = t359 * t373 + t360 * t406 + t378 * t434;
t330 = t355 * t406 - t356 * t373;
t328 = t355 * t379 - t356 * t401 + t364 * t372;
t327 = t346 * t406 + t348 * t373;
t326 = -t342 * t373 - t344 * t406;
t325 = t352 * t373 + t353 * t406 + t370 * t435;
t324 = t350 * t373 + t351 * t406 - t366 * t435;
t322 = -t346 * t379 - t348 * t401 - t363 * t372;
t318 = t323 * t405 + t334 * t400;
t317 = -t323 * t400 + t334 * t405;
t1 = [t438, t325 * t405 + t336 * t400, t327 * t405 - t348 * t433, -t322 * t405, t317, 0; t318, t324 * t405 + t335 * t400, t326 * t405 + t342 * t433, t319 * t405, -t437, 0; 0, t331 * t405 + t349 * t400, t330 * t405 + t356 * t433, t328 * t405, -t329 * t400 + t338 * t405, 0; t437, -t325 * t400 + t336 * t405, -t327 * t400 - t348 * t432, t322 * t400, -t318, 0; t317, -t324 * t400 + t335 * t405, -t326 * t400 + t342 * t432, -t319 * t400, t438, 0; 0, -t331 * t400 + t349 * t405, -t330 * t400 + t356 * t432, -t328 * t400, -t329 * t405 - t338 * t400, 0; t319, -t352 * t379 + t353 * t401 + t370 * t436, t346 * t401 - t348 * t379, t323, 0, 0; t322, -t350 * t379 + t351 * t401 - t366 * t436, t342 * t379 - t344 * t401, t320, 0, 0; 0, -t359 * t379 + t360 * t401 + t372 * t434, t355 * t401 + t356 * t379, t329, 0, 0;];
JR_rot  = t1;
