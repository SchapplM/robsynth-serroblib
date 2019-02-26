% Zeitableitung der rotatorischen Teilmatrix der geometrischen Jacobi-Matrix für Segment Nr. 5 (0=Basis) von
% S6RRRRRR10
% Use Code from Maple symbolic Code Generation
%
% geometrische Jacobi-Matrix: Differentieller Zusammenhang zwischen
% Endeffektorposition und verallgemeinerten Koordinaten.
%
%
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [14x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,alpha3,alpha4,d1,d2,d3,d4,d5,d6]';
%
% Output:
% JgD_rot [3x6]
%   Zeitableitung der rotatorischen Teilmatrix der geometrischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox (ehem. IRT-Maple-Toolbox)
% Datum: 2018-11-23 11:28
% Revision: 76f9d5e39f14dc242b53c0d9d3d9db48bd8f37c0
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für mechatronische Systeme, Universität Hannover

function JgD_rot = S6RRRRRR10_jacobigD_rot_5_floatb_twist_sym_varpar(qJ, qJD, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(14,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRRRR10_jacobigD_rot_5_floatb_twist_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRRRRR10_jacobigD_rot_5_floatb_twist_sym_varpar: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [14 1]), ...
  'S6RRRRRR10_jacobigD_rot_5_floatb_twist_sym_varpar: pkin has to be [14x1] (double)');

%% Symbolic Calculation
% From jacobigD_rot_5_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2018-11-23 11:27:20
% EndTime: 2018-11-23 11:27:20
% DurationCPUTime: 0.29s
% Computational Cost: add. (563->102), mult. (708->160), div. (0->0), fcn. (559->26), ass. (0->91)
t493 = pkin(8) - qJ(4);
t531 = sin(t493) / 0.2e1;
t492 = pkin(8) + qJ(4);
t530 = cos(t492) / 0.2e1;
t529 = qJD(2) / 0.2e1;
t528 = qJD(3) / 0.2e1;
t496 = pkin(6) + qJ(2);
t490 = cos(t496);
t475 = t490 * t529;
t497 = pkin(6) - qJ(2);
t491 = cos(t497);
t521 = qJD(2) * t491;
t461 = t475 - t521 / 0.2e1;
t499 = sin(pkin(7));
t527 = t461 * t499;
t500 = sin(pkin(6));
t507 = sin(qJ(1));
t526 = t500 * t507;
t511 = cos(qJ(1));
t525 = t500 * t511;
t524 = qJD(1) * t507;
t523 = qJD(1) * t511;
t484 = sin(t496);
t522 = qJD(2) * t484;
t520 = qJD(2) * t507;
t519 = qJD(2) * t511;
t494 = pkin(7) + qJ(3);
t482 = sin(t494);
t518 = qJD(3) * t482;
t495 = pkin(7) - qJ(3);
t489 = cos(t495);
t517 = qJD(3) * t489;
t505 = sin(qJ(3));
t516 = qJD(3) * t505;
t509 = cos(qJ(3));
t515 = qJD(3) * t509;
t514 = qJD(4) * cos(qJ(4));
t513 = t500 * t524;
t512 = t500 * t523;
t477 = t484 / 0.2e1;
t485 = sin(t497);
t466 = t477 - t485 / 0.2e1;
t510 = cos(qJ(2));
t449 = t466 * t511 + t507 * t510;
t451 = -t466 * t507 + t510 * t511;
t479 = t491 / 0.2e1;
t470 = t479 + t490 / 0.2e1;
t506 = sin(qJ(2));
t448 = t470 * t511 - t506 * t507;
t450 = -t470 * t507 - t506 * t511;
t504 = sin(qJ(4));
t503 = cos(pkin(6));
t502 = cos(pkin(7));
t501 = cos(pkin(8));
t498 = sin(pkin(8));
t488 = cos(t494);
t487 = cos(t493);
t483 = sin(t495);
t480 = sin(t492);
t478 = t489 / 0.2e1;
t476 = t482 / 0.2e1;
t474 = t485 * t529;
t473 = t488 * t528;
t472 = t483 * t528;
t471 = t479 - t490 / 0.2e1;
t469 = t478 - t488 / 0.2e1;
t468 = t478 + t488 / 0.2e1;
t467 = t487 / 0.2e1 + t530;
t465 = t477 + t485 / 0.2e1;
t464 = t476 - t483 / 0.2e1;
t463 = t476 + t483 / 0.2e1;
t462 = t480 / 0.2e1 + t531;
t460 = t475 + t521 / 0.2e1;
t459 = t474 - t522 / 0.2e1;
t458 = t474 + t522 / 0.2e1;
t457 = t473 - t517 / 0.2e1;
t456 = t473 + t517 / 0.2e1;
t455 = t472 - t518 / 0.2e1;
t454 = t472 + t518 / 0.2e1;
t453 = (t530 - t487 / 0.2e1) * qJD(4);
t452 = (t531 - t480 / 0.2e1) * qJD(4);
t447 = qJD(1) * t451 + t511 * t460 - t506 * t520;
t446 = qJD(1) * t450 + t511 * t459 - t510 * t520;
t445 = -qJD(1) * t449 - t507 * t460 - t506 * t519;
t444 = -qJD(1) * t448 - t507 * t459 - t510 * t519;
t443 = -t446 * t499 + t502 * t513;
t442 = -t444 * t499 + t502 * t512;
t441 = t455 * t465 + t457 * t503 - t458 * t505 + t461 * t468 - t471 * t515;
t440 = -t449 * t515 + t446 * t468 - t447 * t505 + t448 * t455 + (-t457 * t511 + t463 * t524) * t500;
t439 = -t451 * t515 + t444 * t468 - t445 * t505 + t450 * t455 + (t457 * t507 + t463 * t523) * t500;
t1 = [0, t512, t442, -t439 * t498 + t442 * t501 (-t451 * t516 + t444 * t464 + t445 * t509 + t450 * t456 + (t454 * t507 + t469 * t523) * t500) * t504 + (t450 * t464 + t451 * t509 + t469 * t526) * t514 - t439 * t467 - (t450 * t468 - t451 * t505 + t463 * t526) * t452 - t442 * t462 - (-t450 * t499 + t502 * t526) * t453, 0; 0, t513, t443, -t440 * t498 + t443 * t501 (-t449 * t516 + t446 * t464 + t447 * t509 + t448 * t456 + (-t454 * t511 + t469 * t524) * t500) * t504 + (t448 * t464 + t449 * t509 - t469 * t525) * t514 - t440 * t467 - (t448 * t468 - t449 * t505 - t463 * t525) * t452 - t443 * t462 - (-t448 * t499 - t502 * t525) * t453, 0; 0, 0, -t527, -t441 * t498 - t501 * t527 (t454 * t503 + t456 * t465 + t458 * t509 + t461 * t464 - t471 * t516) * t504 + (t464 * t465 + t469 * t503 + t471 * t509) * t514 - t441 * t467 - (t463 * t503 + t465 * t468 - t471 * t505) * t452 + t462 * t527 - (-t465 * t499 + t502 * t503) * t453, 0;];
JgD_rot  = t1;
