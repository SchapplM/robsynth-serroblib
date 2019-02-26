% Zeitableitung der translatorischen Teilmatrix der analytischen Jacobi-Matrix für Segment Nr. 6 (0=Basis) von
% S6RRRPRP9
% Use Code from Maple symbolic Code Generation
%
% analytische Jacobi-Matrix: Differentieller Zusammenhang zwischen
% Endeffektorposition und verallgemeinerten Koordinaten.
%
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% r_i_i_C [3x1]
%   Ortsvektor vom KörperKS-Ursprung zum gesuchten Punkt
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d3,d5]';
%
% Output:
% JaD_transl [3x6]
%   Zeitableitung der translatorischen Teilmatrix der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-02-26 22:13
% Revision: d75aae1ac561373cd3be920984c3df29a1c2ecc8 (2019-02-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JaD_transl = S6RRRPRP9_jacobiaD_transl_6_sym_varpar(qJ, qJD, r_i_i_C, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(3,1),zeros(9,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRPRP9_jacobiaD_transl_6_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRRPRP9_jacobiaD_transl_6_sym_varpar: qJD has to be [6x1] (double)');
assert(isa(r_i_i_C,'double') && isreal(r_i_i_C) && all(size(r_i_i_C) == [3 1]), ...
	'S6RRRPRP9_jacobiaD_transl_6_sym_varpar: Position vector r_i_i_C has to be [3x1] double');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S6RRRPRP9_jacobiaD_transl_6_sym_varpar: pkin has to be [9x1] (double)');

%% Symbolic Calculation
% From jacobiaD_transl_6_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 22:13:46
% EndTime: 2019-02-26 22:13:46
% DurationCPUTime: 0.64s
% Computational Cost: add. (677->110), mult. (2049->174), div. (0->0), fcn. (2002->8), ass. (0->70)
t470 = sin(qJ(3));
t474 = cos(qJ(3));
t475 = cos(qJ(2));
t472 = sin(qJ(1));
t512 = qJD(3) * t474;
t476 = cos(qJ(1));
t517 = qJD(1) * t476;
t484 = t470 * t517 + t472 * t512;
t511 = qJD(3) * t476;
t501 = t470 * t511;
t471 = sin(qJ(2));
t516 = qJD(2) * t472;
t503 = t471 * t516;
t519 = qJD(1) * t472;
t445 = -t470 * t503 - t474 * t519 + t475 * t484 - t501;
t520 = t476 * t474;
t451 = t472 * t470 + t475 * t520;
t523 = t472 * t475;
t507 = t470 * t523;
t446 = t451 * qJD(1) - qJD(3) * t507 + (-t503 - t511) * t474;
t469 = sin(qJ(5));
t473 = cos(qJ(5));
t448 = t507 + t520;
t521 = t476 * t470;
t449 = t474 * t523 - t521;
t497 = t448 * t469 + t449 * t473;
t428 = qJD(5) * t497 - t445 * t473 + t446 * t469;
t528 = r_i_i_C(3) + qJ(6);
t530 = r_i_i_C(1) + pkin(5);
t498 = t448 * t473 - t449 * t469;
t535 = qJD(5) * t498 + t445 * t469 + t446 * t473;
t546 = t497 * qJD(6) - t530 * t428 + t528 * t535;
t494 = t469 * t474 - t470 * t473;
t508 = -r_i_i_C(2) - pkin(9) + pkin(8);
t531 = pkin(4) + pkin(3);
t537 = -t474 * qJ(4) + t470 * t531;
t477 = -t508 * qJD(2) + t537 * qJD(3) - t470 * qJD(4) - t494 * qJD(6);
t514 = qJD(2) * t476;
t502 = t471 * t514;
t518 = qJD(1) * t475;
t443 = (-qJD(3) * t475 + qJD(1)) * t520 + (t502 + (-qJD(3) + t518) * t472) * t470;
t444 = t475 * t501 + (t472 * t518 + t502) * t474 - t484;
t450 = -t472 * t474 + t475 * t521;
t495 = t450 * t469 + t451 * t473;
t424 = qJD(5) * t495 + t443 * t473 - t444 * t469;
t496 = t450 * t473 - t451 * t469;
t425 = qJD(5) * t496 - t443 * t469 - t444 * t473;
t545 = t495 * qJD(6) - t530 * t424 + t425 * t528;
t492 = -qJ(4) * t470 - t474 * t531;
t486 = -pkin(2) + t492;
t542 = t475 * t494;
t515 = qJD(2) * t475;
t540 = t470 * t515 + t471 * t512;
t538 = qJD(3) - qJD(5);
t536 = -t471 * pkin(2) + t508 * t475;
t493 = t469 * t470 + t473 * t474;
t487 = t493 * qJD(5);
t506 = t474 * t515;
t513 = qJD(3) * t470;
t439 = t471 * t487 - t540 * t473 + (-t471 * t513 + t506) * t469;
t533 = t528 * (-t473 * t506 + (t494 * qJD(5) + t473 * t513) * t471 - t540 * t469) + t530 * t439 - t493 * t471 * qJD(6);
t490 = qJD(1) * t494;
t489 = qJD(2) * t494;
t488 = qJD(2) * t493;
t485 = -pkin(2) * t475 - t471 * t508 - pkin(1);
t483 = t475 * t488;
t482 = qJD(2) * t486;
t479 = t538 * t494;
t478 = t538 * t493;
t1 = [t498 * qJD(6) - t445 * qJ(4) - t448 * qJD(4) - t531 * t446 - t530 * t535 - t528 * t428 - t536 * t516 + (-t472 * pkin(7) + t476 * t485) * qJD(1), t530 * (-t476 * t483 + (-t476 * t479 + t493 * t519) * t471) + t528 * (-t514 * t542 + (t472 * t490 + t476 * t478) * t471) + (t476 * t482 - t508 * t519) * t475 + (t477 * t476 - t486 * t519) * t471, -t444 * qJ(4) + t451 * qJD(4) + t531 * t443 - t545, -t443, t545, t424; -t496 * qJD(6) - t443 * qJ(4) + t450 * qJD(4) - t531 * t444 + t530 * t425 + t528 * t424 + t536 * t514 + (t476 * pkin(7) + t472 * t485) * qJD(1), t530 * (-t472 * t483 + (-t472 * t479 - t493 * t517) * t471) - t528 * (t489 * t523 + (t476 * t490 + (-qJD(3) * t493 + t487) * t472) * t471) + (t472 * t482 + t508 * t517) * t475 + (t477 * t472 + t486 * t517) * t471, t446 * qJ(4) + t449 * qJD(4) - t531 * t445 - t546, t445, t546, t428; 0, t530 * (-t471 * t488 + t538 * t542) - t528 * (t471 * t489 + t475 * t478) + t471 * t482 - t477 * t475, -t537 * t515 + (qJD(3) * t492 + qJD(4) * t474) * t471 + t533, t540, -t533, t439;];
JaD_transl  = t1;
