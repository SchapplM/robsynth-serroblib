% Zeitableitung der rotatorischen Teilmatrix der geometrischen Jacobi-Matrix für Segment Nr. 6 (0=Basis) von
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

% Quelle: HybrDyn-Toolbox
% Datum: 2019-02-26 22:53
% Revision: d75aae1ac561373cd3be920984c3df29a1c2ecc8 (2019-02-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JgD_rot = S6RRRRRR10_jacobigD_rot_6_sym_varpar(qJ, qJD, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(14,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRRRR10_jacobigD_rot_6_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRRRRR10_jacobigD_rot_6_sym_varpar: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [14 1]), ...
  'S6RRRRRR10_jacobigD_rot_6_sym_varpar: pkin has to be [14x1] (double)');

%% Symbolic Calculation
% From jacobigD_rot_6_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 22:53:00
% EndTime: 2019-02-26 22:53:01
% DurationCPUTime: 0.39s
% Computational Cost: add. (364->81), mult. (1171->158), div. (0->0), fcn. (1294->16), ass. (0->78)
t504 = sin(qJ(3));
t509 = cos(qJ(3));
t501 = cos(pkin(6));
t510 = cos(qJ(2));
t511 = cos(qJ(1));
t532 = t511 * t510;
t505 = sin(qJ(2));
t506 = sin(qJ(1));
t536 = t506 * t505;
t517 = t501 * t536 - t532;
t533 = t511 * t505;
t535 = t506 * t510;
t494 = -t501 * t535 - t533;
t497 = sin(pkin(7));
t500 = cos(pkin(7));
t498 = sin(pkin(6));
t542 = t498 * t506;
t520 = t494 * t500 + t497 * t542;
t480 = t504 * t520 - t509 * t517;
t493 = t501 * t533 + t535;
t492 = t501 * t532 - t536;
t541 = t498 * t511;
t521 = -t492 * t500 + t497 * t541;
t548 = -t493 * t509 + t504 * t521;
t496 = sin(pkin(8));
t503 = sin(qJ(4));
t545 = t496 * t503;
t544 = t497 * t498;
t543 = t497 * t501;
t499 = cos(pkin(8));
t540 = t499 * t503;
t539 = t504 * t505;
t538 = t504 * t510;
t537 = t505 * t509;
t534 = t509 * t510;
t531 = qJD(1) * t498;
t502 = sin(qJ(5));
t530 = qJD(4) * t502;
t529 = t506 * t531;
t528 = t511 * t531;
t527 = qJD(3) * t543;
t526 = qJD(2) * t505 * t544;
t525 = t496 * t526;
t477 = -t493 * t504 - t509 * t521;
t489 = -t492 * t497 - t500 * t541;
t524 = t477 * t499 + t489 * t496;
t479 = t504 * t517 + t509 * t520;
t490 = -t494 * t497 + t500 * t542;
t523 = t479 * t499 + t490 * t496;
t519 = t500 * t534 - t539;
t487 = t498 * t519 + t509 * t543;
t491 = t501 * t500 - t510 * t544;
t522 = t487 * t499 + t491 * t496;
t518 = t500 * t538 + t537;
t483 = -qJD(1) * t492 + qJD(2) * t517;
t516 = t483 * t500 + t497 * t528;
t485 = qJD(1) * t494 - qJD(2) * t493;
t515 = t485 * t500 + t497 * t529;
t508 = cos(qJ(4));
t514 = t503 * t524 - t508 * t548;
t513 = t480 * t508 + t503 * t523;
t488 = t498 * t518 + t504 * t543;
t512 = t488 * t508 + t503 * t522;
t507 = cos(qJ(5));
t486 = -qJD(1) * t517 + qJD(2) * t492;
t484 = -qJD(1) * t493 + qJD(2) * t494;
t482 = -t485 * t497 + t500 * t529;
t481 = -t483 * t497 + t500 * t528;
t476 = t509 * t527 + (t519 * qJD(3) + (-t500 * t539 + t534) * qJD(2)) * t498;
t475 = -t504 * t527 + (-t518 * qJD(3) + (-t500 * t537 - t538) * qJD(2)) * t498;
t474 = -t475 * t496 + t499 * t526;
t473 = qJD(3) * t477 + t486 * t509 + t504 * t515;
t472 = t548 * qJD(3) - t486 * t504 + t515 * t509;
t471 = qJD(3) * t479 + t484 * t509 + t504 * t516;
t470 = -t480 * qJD(3) - t484 * t504 + t516 * t509;
t469 = -t472 * t496 + t482 * t499;
t468 = -t470 * t496 + t481 * t499;
t1 = [0, t528, t481, t468, t471 * t503 + (-t470 * t499 - t481 * t496) * t508 + t513 * qJD(4) (t470 * t540 + t471 * t508 + t481 * t545) * t502 - t468 * t507 + (t513 * t507 + (-t479 * t496 + t490 * t499) * t502) * qJD(5) + (-t480 * t503 + t508 * t523) * t530; 0, t529, t482, t469, t473 * t503 + (-t472 * t499 - t482 * t496) * t508 + t514 * qJD(4) (t472 * t540 + t473 * t508 + t482 * t545) * t502 - t469 * t507 + (t514 * t507 + (-t477 * t496 + t489 * t499) * t502) * qJD(5) + (t503 * t548 + t508 * t524) * t530; 0, 0, t526, t474, t476 * t503 + (-t475 * t499 - t525) * t508 + t512 * qJD(4) (t475 * t540 + t476 * t508 + t503 * t525) * t502 - t474 * t507 + (t512 * t507 + (-t487 * t496 + t491 * t499) * t502) * qJD(5) + (-t488 * t503 + t508 * t522) * t530;];
JgD_rot  = t1;
