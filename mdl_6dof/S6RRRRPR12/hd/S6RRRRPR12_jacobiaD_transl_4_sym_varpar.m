% Zeitableitung der translatorischen Teilmatrix der analytischen Jacobi-Matrix für Segment Nr. 4 (0=Basis) von
% S6RRRRPR12
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
% pkin [13x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,alpha3,d1,d2,d3,d4,d6,theta5]';
%
% Output:
% JaD_transl [3x6]
%   Zeitableitung der translatorischen Teilmatrix der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-02-26 22:37
% Revision: d75aae1ac561373cd3be920984c3df29a1c2ecc8 (2019-02-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JaD_transl = S6RRRRPR12_jacobiaD_transl_4_sym_varpar(qJ, qJD, r_i_i_C, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(3,1),zeros(13,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRRPR12_jacobiaD_transl_4_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRRRPR12_jacobiaD_transl_4_sym_varpar: qJD has to be [6x1] (double)');
assert(isa(r_i_i_C,'double') && isreal(r_i_i_C) && all(size(r_i_i_C) == [3 1]), ...
	'S6RRRRPR12_jacobiaD_transl_4_sym_varpar: Position vector r_i_i_C has to be [3x1] double');
assert(isreal(pkin) && all(size(pkin) == [13 1]), ...
  'S6RRRRPR12_jacobiaD_transl_4_sym_varpar: pkin has to be [13x1] (double)');

%% Symbolic Calculation
% From jacobiaD_transl_4_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 22:37:01
% EndTime: 2019-02-26 22:37:03
% DurationCPUTime: 0.95s
% Computational Cost: add. (730->138), mult. (2280->247), div. (0->0), fcn. (2372->12), ass. (0->84)
t501 = sin(qJ(2));
t502 = sin(qJ(1));
t505 = cos(qJ(2));
t506 = cos(qJ(1));
t552 = cos(pkin(6));
t527 = t506 * t552;
t486 = t502 * t501 - t505 * t527;
t496 = sin(pkin(7));
t498 = cos(pkin(7));
t497 = sin(pkin(6));
t545 = t497 * t506;
t519 = t486 * t498 + t496 * t545;
t528 = t502 * t552;
t523 = t501 * t528;
t536 = qJD(2) * t501;
t538 = t506 * t505;
t478 = -qJD(1) * t523 - t502 * t536 + (qJD(2) * t552 + qJD(1)) * t538;
t487 = t501 * t527 + t502 * t505;
t500 = sin(qJ(3));
t504 = cos(qJ(3));
t510 = t506 * t501 + t505 * t528;
t477 = t510 * qJD(1) + t487 * qJD(2);
t537 = qJD(1) * t497;
t530 = t502 * t537;
t512 = -t477 * t498 + t496 * t530;
t453 = (-qJD(3) * t487 + t512) * t500 + (-t519 * qJD(3) + t478) * t504;
t550 = t477 * t496;
t469 = t498 * t530 + t550;
t499 = sin(qJ(4));
t503 = cos(qJ(4));
t561 = t453 * t499 - t469 * t503;
t560 = -t453 * t503 - t469 * t499;
t464 = -t487 * t504 + t519 * t500;
t481 = -t486 * t496 + t498 * t545;
t559 = -t464 * t499 + t481 * t503;
t558 = t464 * t503 + t481 * t499;
t553 = r_i_i_C(3) + pkin(11);
t511 = t523 - t538;
t475 = t486 * qJD(1) + t511 * qJD(2);
t551 = t475 * t496;
t547 = t496 * t497;
t546 = t497 * t502;
t544 = t498 * t500;
t543 = t498 * t504;
t542 = t500 * t501;
t541 = t500 * t505;
t540 = t501 * t504;
t539 = t504 * t505;
t535 = qJD(2) * t505;
t534 = qJD(4) * t499;
t533 = qJD(4) * t503;
t531 = pkin(10) * t498 + pkin(9);
t529 = t506 * t537;
t526 = t552 * t496;
t525 = t496 * t529;
t524 = t536 * t547;
t522 = qJD(3) * t526;
t520 = t503 * r_i_i_C(1) - t499 * r_i_i_C(2) + pkin(3);
t473 = -t486 * t504 - t487 * t544;
t474 = -t504 * t510 + t511 * t544;
t518 = t496 * t546 - t498 * t510;
t517 = t498 * t539 - t542;
t516 = -t498 * t540 - t541;
t515 = t498 * t541 + t540;
t514 = t498 * t542 - t539;
t513 = qJD(4) * (-t499 * r_i_i_C(1) - t503 * r_i_i_C(2));
t508 = t500 * t511 + t518 * t504;
t466 = t518 * t500 - t504 * t511;
t507 = t464 * qJD(3) - t478 * t500 + t512 * t504;
t485 = t552 * t498 - t505 * t547;
t484 = t514 * t497;
t483 = t496 * t510 + t498 * t546;
t480 = t515 * t497 + t500 * t526;
t476 = t487 * qJD(1) + t510 * qJD(2);
t471 = (-t515 * qJD(2) + t516 * qJD(3)) * t497;
t467 = t498 * t529 - t551;
t461 = t504 * t522 + (-t514 * qJD(2) + t517 * qJD(3)) * t497;
t459 = -t478 * t544 - t477 * t504 + (t486 * t500 - t487 * t543) * qJD(3);
t457 = t476 * t544 + t475 * t504 + (t500 * t510 + t511 * t543) * qJD(3);
t451 = -t476 * t504 + (t475 * t498 + t525) * t500 + t508 * qJD(3);
t450 = t466 * qJD(3) - t475 * t543 - t476 * t500 - t504 * t525;
t449 = t451 * t503 + t467 * t499 + (-t466 * t499 + t483 * t503) * qJD(4);
t448 = -t451 * t499 + t467 * t503 + (-t466 * t503 - t483 * t499) * qJD(4);
t1 = [t560 * r_i_i_C(1) + t561 * r_i_i_C(2) - t453 * pkin(3) - t478 * pkin(2) - pkin(10) * t550 + t553 * t507 + (t559 * r_i_i_C(1) - t558 * r_i_i_C(2)) * qJD(4) + (-t506 * pkin(1) - t531 * t546) * qJD(1) (t457 * t503 - t474 * t534) * r_i_i_C(1) + (-t457 * t499 - t474 * t533) * r_i_i_C(2) + t457 * pkin(3) + t475 * pkin(2) + t553 * (t474 * qJD(3) + t475 * t500 - t476 * t543) + ((-t476 * t499 - t511 * t533) * r_i_i_C(1) + (-t476 * t503 + t511 * t534) * r_i_i_C(2) - t476 * pkin(10)) * t496, -t520 * t450 + t553 * t451 + t508 * t513, t448 * r_i_i_C(1) - t449 * r_i_i_C(2), 0, 0; -pkin(10) * t551 - t476 * pkin(2) + t451 * pkin(3) + t449 * r_i_i_C(1) + t448 * r_i_i_C(2) + t553 * t450 + (-pkin(1) * t502 + t531 * t545) * qJD(1) (t459 * t503 - t473 * t534) * r_i_i_C(1) + (-t459 * t499 - t473 * t533) * r_i_i_C(2) + t459 * pkin(3) - t477 * pkin(2) + t553 * (t473 * qJD(3) - t477 * t500 + t478 * t543) + ((t478 * t499 + t487 * t533) * r_i_i_C(1) + (t478 * t503 - t487 * t534) * r_i_i_C(2) + t478 * pkin(10)) * t496, t553 * t453 + (-t487 * t500 - t519 * t504) * t513 + t520 * t507, -t561 * r_i_i_C(1) + t560 * r_i_i_C(2) + (t558 * r_i_i_C(1) + t559 * r_i_i_C(2)) * qJD(4), 0, 0; 0 (t471 * t503 + t484 * t534) * r_i_i_C(1) + (-t471 * t499 + t484 * t533) * r_i_i_C(2) + t471 * pkin(3) + (-t553 * (-t517 * qJD(2) + t514 * qJD(3)) - pkin(2) * t536 + ((t499 * t535 + t501 * t533) * r_i_i_C(1) + (-t501 * t534 + t503 * t535) * r_i_i_C(2) + pkin(10) * t535) * t496) * t497, t553 * t461 + (t517 * t497 + t504 * t526) * t513 + t520 * (-t500 * t522 + (t516 * qJD(2) - t515 * qJD(3)) * t497) (-t461 * t499 + t503 * t524) * r_i_i_C(1) + (-t461 * t503 - t499 * t524) * r_i_i_C(2) + ((-t480 * t503 - t485 * t499) * r_i_i_C(1) + (t480 * t499 - t485 * t503) * r_i_i_C(2)) * qJD(4), 0, 0;];
JaD_transl  = t1;
