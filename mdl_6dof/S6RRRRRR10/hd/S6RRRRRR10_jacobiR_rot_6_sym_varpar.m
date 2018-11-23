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

% Quelle: HybrDyn-Toolbox (ehem. IRT-Maple-Toolbox)
% Datum: 2018-11-23 11:28
% Revision: 76f9d5e39f14dc242b53c0d9d3d9db48bd8f37c0
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für mechatronische Systeme, Universität Hannover

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
% StartTime: 2018-11-23 11:27:21
% EndTime: 2018-11-23 11:27:22
% DurationCPUTime: 1.09s
% Computational Cost: add. (3844->129), mult. (3865->224), div. (0->0), fcn. (3942->30), ass. (0->117)
t502 = sin(qJ(2));
t503 = sin(qJ(1));
t509 = cos(qJ(1));
t533 = pkin(6) + qJ(2);
t522 = cos(t533) / 0.2e1;
t534 = pkin(6) - qJ(2);
t528 = cos(t534);
t510 = t528 / 0.2e1 + t522;
t464 = t503 * t502 - t509 * t510;
t519 = sin(t533) / 0.2e1;
t525 = sin(t534);
t475 = t519 - t525 / 0.2e1;
t508 = cos(qJ(2));
t465 = t509 * t475 + t503 * t508;
t531 = pkin(7) + qJ(3);
t518 = sin(t531) / 0.2e1;
t532 = pkin(7) - qJ(3);
t524 = sin(t532);
t472 = t518 + t524 / 0.2e1;
t521 = cos(t531) / 0.2e1;
t527 = cos(t532);
t479 = t527 / 0.2e1 + t521;
t501 = sin(qJ(3));
t494 = sin(pkin(6));
t537 = t494 * t509;
t440 = t464 * t479 + t465 * t501 + t472 * t537;
t493 = sin(pkin(7));
t496 = cos(pkin(7));
t459 = -t464 * t493 + t496 * t537;
t492 = sin(pkin(8));
t495 = cos(pkin(8));
t428 = t440 * t492 - t459 * t495;
t499 = sin(qJ(5));
t505 = cos(qJ(5));
t473 = t518 - t524 / 0.2e1;
t478 = t521 - t527 / 0.2e1;
t507 = cos(qJ(3));
t438 = -t464 * t473 + t465 * t507 + t478 * t537;
t529 = pkin(8) + qJ(4);
t517 = sin(t529) / 0.2e1;
t530 = pkin(8) - qJ(4);
t523 = sin(t530);
t471 = t517 - t523 / 0.2e1;
t520 = cos(t529) / 0.2e1;
t526 = cos(t530);
t476 = t520 - t526 / 0.2e1;
t506 = cos(qJ(4));
t514 = t438 * t506 - t440 * t471 + t459 * t476;
t392 = t428 * t499 + t505 * t514;
t477 = t526 / 0.2e1 + t520;
t500 = sin(qJ(4));
t515 = t517 + t523 / 0.2e1;
t402 = t438 * t500 + t440 * t477 + t459 * t515;
t498 = sin(qJ(6));
t504 = cos(qJ(6));
t546 = t392 * t498 - t402 * t504;
t545 = -t392 * t504 - t402 * t498;
t391 = t428 * t505 - t499 * t514;
t542 = t476 * t493;
t541 = t492 * t499;
t540 = t492 * t505;
t539 = t493 * t495;
t538 = t494 * t503;
t536 = t498 * t505;
t535 = t504 * t505;
t467 = -t509 * t502 - t503 * t510;
t469 = t503 * t475 - t509 * t508;
t442 = t467 * t479 + t469 * t501 + t472 * t538;
t461 = -t467 * t493 + t496 * t538;
t516 = -t442 * t492 + t461 * t495;
t444 = -t467 * t473 + t469 * t507 + t478 * t538;
t513 = t442 * t471 - t444 * t506 - t461 * t476;
t474 = t519 + t525 / 0.2e1;
t480 = t522 - t528 / 0.2e1;
t497 = cos(pkin(6));
t452 = t497 * t472 + t474 * t479 + t480 * t501;
t453 = t474 * t473 - t497 * t478 - t480 * t507;
t463 = -t474 * t493 + t497 * t496;
t512 = t452 * t471 + t453 * t506 - t463 * t476;
t511 = t493 * t515;
t458 = t480 * t473 + t474 * t507;
t457 = -t474 * t501 + t480 * t479;
t450 = t467 * t507 + t469 * t473;
t449 = -t467 * t501 + t469 * t479;
t448 = -t464 * t507 - t465 * t473;
t447 = t464 * t501 - t465 * t479;
t446 = -t457 * t492 - t480 * t539;
t434 = -t452 * t492 + t463 * t495;
t431 = -t449 * t492 - t469 * t539;
t430 = -t447 * t492 + t465 * t539;
t427 = t457 * t471 + t458 * t506 + t480 * t542;
t426 = -t457 * t477 + t458 * t500 + t480 * t511;
t424 = t452 * t506 - t453 * t471;
t423 = t452 * t500 + t453 * t477;
t420 = -t452 * t477 + t453 * t500 - t463 * t515;
t419 = t442 * t506 + t444 * t471;
t418 = t442 * t500 - t444 * t477;
t417 = -t438 * t471 - t440 * t506;
t416 = t438 * t477 - t440 * t500;
t415 = t449 * t471 + t450 * t506 + t469 * t542;
t414 = -t449 * t477 + t450 * t500 + t469 * t511;
t413 = t447 * t471 + t448 * t506 - t465 * t542;
t412 = -t447 * t477 + t448 * t500 - t465 * t511;
t411 = t424 * t505 + t453 * t541;
t410 = t427 * t505 + t446 * t499;
t407 = -t442 * t477 - t444 * t500 - t461 * t515;
t401 = t434 * t499 + t505 * t512;
t400 = t434 * t505 - t499 * t512;
t399 = t419 * t505 - t444 * t541;
t398 = t417 * t505 + t438 * t541;
t397 = t415 * t505 + t431 * t499;
t396 = t413 * t505 + t430 * t499;
t395 = t499 * t516 + t505 * t513;
t394 = t499 * t513 - t505 * t516;
t390 = t395 * t504 + t407 * t498;
t389 = -t395 * t498 + t407 * t504;
t1 = [t545, t397 * t504 + t414 * t498, t399 * t504 + t418 * t498, -t407 * t535 + t498 * t513, -t394 * t504, t389; t390, t396 * t504 + t412 * t498, t398 * t504 + t416 * t498, -t402 * t535 + t498 * t514, t391 * t504, -t546; 0, t410 * t504 + t426 * t498, t411 * t504 + t423 * t498, -t420 * t535 + t498 * t512, t400 * t504, -t401 * t498 + t420 * t504; t546, -t397 * t498 + t414 * t504, -t399 * t498 + t418 * t504, t407 * t536 + t504 * t513, t394 * t498, -t390; t389, -t396 * t498 + t412 * t504, -t398 * t498 + t416 * t504, t402 * t536 + t504 * t514, -t391 * t498, t545; 0, -t410 * t498 + t426 * t504, -t411 * t498 + t423 * t504, t420 * t536 + t504 * t512, -t400 * t498, -t401 * t504 - t420 * t498; t391, t415 * t499 - t431 * t505, t419 * t499 + t444 * t540, -t407 * t499, t395, 0; t394, t413 * t499 - t430 * t505, t417 * t499 - t438 * t540, -t402 * t499, t392, 0; 0, t427 * t499 - t446 * t505, t424 * t499 - t453 * t540, -t420 * t499, t401, 0;];
JR_rot  = t1;
