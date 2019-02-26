% Zeitableitung der translatorischen Teilmatrix der analytischen Jacobi-Matrix für Segment Nr. 6 (0=Basis) von
% S6RRRRPR13
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
% pkin [11x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d1,d2,d3,d4,d6]';
%
% Output:
% JaD_transl [3x6]
%   Zeitableitung der translatorischen Teilmatrix der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-02-26 22:37
% Revision: d75aae1ac561373cd3be920984c3df29a1c2ecc8 (2019-02-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JaD_transl = S6RRRRPR13_jacobiaD_transl_6_sym_varpar(qJ, qJD, r_i_i_C, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(3,1),zeros(11,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRRPR13_jacobiaD_transl_6_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRRRPR13_jacobiaD_transl_6_sym_varpar: qJD has to be [6x1] (double)');
assert(isa(r_i_i_C,'double') && isreal(r_i_i_C) && all(size(r_i_i_C) == [3 1]), ...
	'S6RRRRPR13_jacobiaD_transl_6_sym_varpar: Position vector r_i_i_C has to be [3x1] double');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6RRRRPR13_jacobiaD_transl_6_sym_varpar: pkin has to be [11x1] (double)');

%% Symbolic Calculation
% From jacobiaD_transl_6_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 22:37:39
% EndTime: 2019-02-26 22:37:40
% DurationCPUTime: 1.36s
% Computational Cost: add. (1624->179), mult. (4756->297), div. (0->0), fcn. (5097->12), ass. (0->100)
t503 = sin(qJ(2));
t504 = sin(qJ(1));
t508 = cos(qJ(2));
t561 = cos(pkin(6));
t563 = cos(qJ(1));
t526 = t561 * t563;
t484 = t503 * t526 + t504 * t508;
t502 = sin(qJ(3));
t507 = cos(qJ(3));
t499 = sin(pkin(6));
t537 = t499 * t563;
t471 = t484 * t507 - t502 * t537;
t483 = t504 * t503 - t508 * t526;
t501 = sin(qJ(4));
t506 = cos(qJ(4));
t450 = t471 * t501 - t483 * t506;
t530 = t504 * t561;
t485 = t563 * t503 + t508 * t530;
t464 = qJD(1) * t485 + qJD(2) * t484;
t574 = -qJD(4) * t450 + t464 * t501;
t451 = t471 * t506 + t483 * t501;
t573 = qJD(4) * t451 - t464 * t506;
t500 = sin(qJ(6));
t505 = cos(qJ(6));
t572 = t450 * t505 - t451 * t500;
t571 = t450 * t500 + t451 * t505;
t570 = (t571 * r_i_i_C(1) + t572 * r_i_i_C(2)) * qJD(6);
t551 = t503 * t507;
t482 = t499 * t551 + t561 * t502;
t548 = t506 * t508;
t468 = t482 * t501 + t499 * t548;
t552 = t501 * t508;
t558 = t482 * t506;
t469 = -t499 * t552 + t558;
t569 = ((t468 * t500 + t469 * t505) * r_i_i_C(1) + (t468 * t505 - t469 * t500) * r_i_i_C(2)) * qJD(6);
t540 = -r_i_i_C(3) - pkin(11) + pkin(10);
t568 = t507 * pkin(3) + t502 * t540 + pkin(2);
t567 = -pkin(3) * t502 + t507 * t540;
t564 = pkin(4) + pkin(5);
t516 = t505 * r_i_i_C(1) - t500 * r_i_i_C(2) + t564;
t517 = t500 * r_i_i_C(1) + t505 * r_i_i_C(2) + qJ(5);
t511 = t501 * t517 + t506 * t516 + pkin(3);
t555 = t499 * t504;
t554 = t501 * t503;
t553 = t501 * t507;
t550 = t504 * t507;
t549 = t506 * t507;
t547 = t507 * t508;
t546 = qJD(1) * t504;
t545 = qJD(2) * t503;
t544 = qJD(3) * t502;
t543 = qJD(3) * t508;
t542 = qJD(4) * t499;
t541 = qJD(4) * t507;
t539 = t499 * t503 * t502;
t536 = t499 * t546;
t535 = t499 * t545;
t534 = qJD(2) * t499 * t508;
t533 = t502 * t543;
t532 = t508 * t542;
t531 = t563 * qJD(1);
t527 = t503 * t530;
t465 = -qJD(1) * t527 - t504 * t545 + (qJD(2) * t526 + t531) * t508;
t528 = t507 * t537;
t529 = qJD(3) * t528 - t465 * t507;
t463 = qJD(1) * t484 + qJD(2) * t485;
t525 = t485 * t541 - t463;
t524 = t483 * t541 + t465;
t486 = t563 * t508 - t527;
t474 = t486 * t507 + t502 * t555;
t454 = t474 * t501 - t485 * t506;
t455 = t474 * t506 + t485 * t501;
t521 = t454 * t505 - t455 * t500;
t520 = t454 * t500 + t455 * t505;
t514 = qJD(3) * t567;
t462 = qJD(2) * t527 + t503 * t546 + (-qJD(1) * t526 - t563 * qJD(2)) * t508;
t513 = qJD(4) * t486 + t462 * t507 + t485 * t544;
t512 = qJD(4) * t484 - t464 * t507 + t483 * t544;
t510 = -t471 * qJD(3) - t465 * t502 + t507 * t536;
t509 = t501 * qJD(5) + ((-t500 * t506 + t501 * t505) * r_i_i_C(1) + (-t500 * t501 - t505 * t506) * r_i_i_C(2)) * qJD(6) + (-t501 * t516 + t506 * t517) * qJD(4);
t476 = (t506 * t547 + t554) * t499;
t475 = (t501 * t547 - t503 * t506) * t499;
t467 = -qJD(3) * t539 + (t561 * qJD(3) + t534) * t507;
t459 = -t485 * t549 + t486 * t501;
t458 = -t485 * t553 - t486 * t506;
t457 = -t483 * t549 + t484 * t501;
t456 = -t483 * t553 - t484 * t506;
t447 = -qJD(4) * t468 + t467 * t506 + t501 * t535;
t446 = qJD(4) * t558 - t506 * t535 + (t467 - t532) * t501;
t445 = (qJD(3) * t484 - t536) * t502 + t529;
t443 = -t484 * t544 + t502 * t536 - t529;
t441 = -t463 * t507 - t486 * t544 + (qJD(3) * t550 + t502 * t531) * t499;
t440 = -qJD(1) * t528 + t474 * qJD(3) - t463 * t502;
t433 = t443 * t506 + t574;
t432 = t443 * t501 + t573;
t431 = -qJD(4) * t454 + t441 * t506 - t462 * t501;
t430 = qJD(4) * t455 + t441 * t501 + t462 * t506;
t429 = qJD(6) * t521 + t430 * t500 + t431 * t505;
t428 = -qJD(6) * t520 + t430 * t505 - t431 * t500;
t1 = [-t465 * pkin(2) + t445 * pkin(3) - t464 * pkin(9) - t450 * qJD(5) + t517 * (t445 * t501 - t573) + t516 * (t445 * t506 - t574) + (-t572 * r_i_i_C(1) + t571 * r_i_i_C(2)) * qJD(6) + (-t563 * pkin(1) - pkin(8) * t555) * qJD(1) + t540 * t510, -t463 * pkin(9) + t458 * qJD(5) + t517 * (t501 * t513 - t506 * t525) + t516 * (t501 * t525 + t506 * t513) + ((t458 * t505 - t459 * t500) * r_i_i_C(1) + (-t458 * t500 - t459 * t505) * r_i_i_C(2)) * qJD(6) - t485 * t514 + t568 * t462, t540 * t441 - t511 * t440 + t509 * (-t486 * t502 + t499 * t550) t455 * qJD(5) + t517 * t431 - t516 * t430 + (r_i_i_C(1) * t520 + r_i_i_C(2) * t521) * qJD(6), t430, t428 * r_i_i_C(1) - t429 * r_i_i_C(2); -t463 * pkin(2) + t441 * pkin(3) - t462 * pkin(9) + t429 * r_i_i_C(1) + t428 * r_i_i_C(2) + t430 * qJ(5) + t454 * qJD(5) + t564 * t431 + (-pkin(1) * t504 + pkin(8) * t537) * qJD(1) + t540 * t440, t465 * pkin(9) + t456 * qJD(5) + t517 * (t501 * t512 - t506 * t524) + t516 * (t501 * t524 + t506 * t512) + ((t456 * t505 - t457 * t500) * r_i_i_C(1) + (-t456 * t500 - t457 * t505) * r_i_i_C(2)) * qJD(6) - t483 * t514 - t568 * t464, t540 * t443 + t511 * t510 + t509 * (-t484 * t502 - t528) t451 * qJD(5) - t516 * t432 + t517 * t433 + t570, t432 (t432 * t505 - t433 * t500) * r_i_i_C(1) + (-t432 * t500 - t433 * t505) * r_i_i_C(2) - t570; 0, t475 * qJD(5) - t517 * (-t532 * t549 - t542 * t554) + ((t475 * t505 - t476 * t500) * r_i_i_C(1) + (-t475 * t500 - t476 * t505) * r_i_i_C(2)) * qJD(6) + (-t517 * (t501 * t533 + (t501 * t551 + t548) * qJD(2)) + t516 * ((qJD(2) - t541) * t552 + (-t533 + (-qJD(2) * t507 + qJD(4)) * t503) * t506) + t567 * t543 + (pkin(9) * t508 - t503 * t568) * qJD(2)) * t499, t540 * t467 + t511 * (-t482 * qJD(3) - t502 * t534) + t509 * (t561 * t507 - t539) t469 * qJD(5) - t516 * t446 + t517 * t447 + t569, t446 (t446 * t505 - t447 * t500) * r_i_i_C(1) + (-t446 * t500 - t447 * t505) * r_i_i_C(2) - t569;];
JaD_transl  = t1;
