% Zeitableitung der translatorischen Teilmatrix der analytischen Jacobi-Matrix für Segment Nr. 6 (0=Basis) von
% S6PRRPRR4
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
%   pkin=[a2,a3,a4,a5,a6,alpha2,d2,d3,d5,d6,theta1]';
%
% Output:
% JaD_transl [3x6]
%   Zeitableitung der translatorischen Teilmatrix der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-02-26 20:06
% Revision: d75aae1ac561373cd3be920984c3df29a1c2ecc8 (2019-02-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JaD_transl = S6PRRPRR4_jacobiaD_transl_6_sym_varpar(qJ, qJD, r_i_i_C, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(3,1),zeros(11,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRRPRR4_jacobiaD_transl_6_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6PRRPRR4_jacobiaD_transl_6_sym_varpar: qJD has to be [6x1] (double)');
assert(isa(r_i_i_C,'double') && isreal(r_i_i_C) && all(size(r_i_i_C) == [3 1]), ...
	'S6PRRPRR4_jacobiaD_transl_6_sym_varpar: Position vector r_i_i_C has to be [3x1] double');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6PRRPRR4_jacobiaD_transl_6_sym_varpar: pkin has to be [11x1] (double)');

%% Symbolic Calculation
% From jacobiaD_transl_6_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 20:05:56
% EndTime: 2019-02-26 20:05:56
% DurationCPUTime: 0.71s
% Computational Cost: add. (905->115), mult. (2786->204), div. (0->0), fcn. (3016->12), ass. (0->72)
t517 = sin(qJ(5));
t518 = sin(qJ(3));
t521 = cos(qJ(5));
t522 = cos(qJ(3));
t534 = t517 * t518 + t521 * t522;
t565 = qJD(5) - qJD(3);
t524 = t565 * t534;
t515 = cos(pkin(6));
t519 = sin(qJ(2));
t513 = sin(pkin(6));
t557 = t513 * t522;
t506 = t515 * t518 + t519 * t557;
t523 = cos(qJ(2));
t556 = t513 * t523;
t548 = qJD(2) * t556;
t494 = qJD(3) * t506 + t518 * t548;
t558 = t513 * t518;
t505 = -t515 * t522 + t519 * t558;
t495 = -qJD(3) * t505 + t522 * t548;
t536 = t505 * t521 - t506 * t517;
t472 = qJD(5) * t536 + t494 * t517 + t495 * t521;
t489 = t505 * t517 + t506 * t521;
t516 = sin(qJ(6));
t520 = cos(qJ(6));
t545 = -t516 * r_i_i_C(1) - t520 * r_i_i_C(2);
t528 = qJD(6) * t545;
t546 = -t520 * r_i_i_C(1) + t516 * r_i_i_C(2);
t533 = pkin(5) - t546;
t560 = r_i_i_C(3) + pkin(10);
t574 = -t533 * (qJD(5) * t489 - t494 * t521 + t495 * t517) + t560 * t472 + t536 * t528;
t512 = sin(pkin(11));
t514 = cos(pkin(11));
t555 = t515 * t519;
t529 = t512 * t555 - t514 * t523;
t493 = t512 * t558 - t522 * t529;
t554 = t515 * t523;
t530 = t512 * t554 + t514 * t519;
t499 = t530 * qJD(2);
t483 = qJD(3) * t493 - t499 * t518;
t531 = t512 * t557 + t518 * t529;
t484 = qJD(3) * t531 - t499 * t522;
t539 = -t493 * t517 - t521 * t531;
t464 = qJD(5) * t539 + t483 * t517 + t484 * t521;
t480 = t493 * t521 - t517 * t531;
t573 = -t533 * (qJD(5) * t480 - t483 * t521 + t484 * t517) + t560 * t464 + t539 * t528;
t502 = t512 * t523 + t514 * t555;
t491 = t502 * t522 - t514 * t558;
t550 = t514 * t554;
t553 = qJD(2) * t519;
t497 = -qJD(2) * t550 + t512 * t553;
t481 = t491 * qJD(3) - t497 * t518;
t490 = t502 * t518 + t514 * t557;
t482 = -qJD(3) * t490 - t497 * t522;
t540 = t490 * t521 - t491 * t517;
t460 = qJD(5) * t540 + t481 * t517 + t482 * t521;
t477 = t490 * t517 + t491 * t521;
t572 = -t533 * (qJD(5) * t477 - t481 * t521 + t482 * t517) + t560 * t460 + t540 * t528;
t561 = pkin(3) + pkin(4);
t527 = t518 * qJ(4) + t522 * t561 + pkin(2);
t552 = qJD(6) * t534 * t556;
t549 = t513 * t553;
t535 = t517 * t522 - t518 * t521;
t532 = -pkin(8) + pkin(9) - t545;
t526 = t518 * qJD(4) + (qJ(4) * t522 - t518 * t561) * qJD(3);
t525 = t565 * t535;
t501 = -t512 * t519 + t550;
t500 = t529 * qJD(2);
t498 = t502 * qJD(2);
t486 = t534 * t530;
t485 = t534 * t501;
t474 = (-t523 * t525 - t534 * t553) * t513;
t1 = [0, t533 * (t500 * t534 + t525 * t530) + t532 * t499 + t560 * (t500 * t535 - t524 * t530) + ((t486 * t516 + t520 * t529) * r_i_i_C(1) + (t486 * t520 - t516 * t529) * r_i_i_C(2)) * qJD(6) - t526 * t530 + t527 * t500, t484 * qJ(4) + t493 * qJD(4) - t561 * t483 - t573, t483, t573 (-t464 * t516 + t500 * t520) * r_i_i_C(1) + (-t464 * t520 - t500 * t516) * r_i_i_C(2) + ((-t480 * t520 + t516 * t530) * r_i_i_C(1) + (t480 * t516 + t520 * t530) * r_i_i_C(2)) * qJD(6); 0, t533 * (-t498 * t534 - t501 * t525) + t532 * t497 + t560 * (-t498 * t535 + t501 * t524) + ((-t485 * t516 - t502 * t520) * r_i_i_C(1) + (-t485 * t520 + t502 * t516) * r_i_i_C(2)) * qJD(6) + t526 * t501 - t527 * t498, t482 * qJ(4) + t491 * qJD(4) - t561 * t481 - t572, t481, t572 (-t460 * t516 - t498 * t520) * r_i_i_C(1) + (-t460 * t520 + t498 * t516) * r_i_i_C(2) + ((-t477 * t520 - t501 * t516) * r_i_i_C(1) + (t477 * t516 - t501 * t520) * r_i_i_C(2)) * qJD(6); 0 (t474 * t520 - t516 * t552) * r_i_i_C(1) + (-t474 * t516 - t520 * t552) * r_i_i_C(2) + t474 * pkin(5) + (-t560 * (-t523 * t524 + t535 * t553) + t546 * t519 * qJD(6) + t526 * t523 + (-t527 * t519 - t532 * t523) * qJD(2)) * t513, t495 * qJ(4) + t506 * qJD(4) - t561 * t494 - t574, t494, t574 (-t472 * t516 - t520 * t549) * r_i_i_C(1) + (-t472 * t520 + t516 * t549) * r_i_i_C(2) + ((-t489 * t520 - t516 * t556) * r_i_i_C(1) + (t489 * t516 - t520 * t556) * r_i_i_C(2)) * qJD(6);];
JaD_transl  = t1;
