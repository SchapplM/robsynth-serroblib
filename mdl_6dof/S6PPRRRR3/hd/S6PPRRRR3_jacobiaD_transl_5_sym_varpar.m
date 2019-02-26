% Zeitableitung der translatorischen Teilmatrix der analytischen Jacobi-Matrix für Segment Nr. 5 (0=Basis) von
% S6PPRRRR3
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
% pkin [14x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,alpha3,alpha4,d3,d4,d5,d6,theta1,theta2]';
%
% Output:
% JaD_transl [3x6]
%   Zeitableitung der translatorischen Teilmatrix der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-02-26 19:44
% Revision: d75aae1ac561373cd3be920984c3df29a1c2ecc8 (2019-02-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JaD_transl = S6PPRRRR3_jacobiaD_transl_5_sym_varpar(qJ, qJD, r_i_i_C, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(3,1),zeros(14,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PPRRRR3_jacobiaD_transl_5_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6PPRRRR3_jacobiaD_transl_5_sym_varpar: qJD has to be [6x1] (double)');
assert(isa(r_i_i_C,'double') && isreal(r_i_i_C) && all(size(r_i_i_C) == [3 1]), ...
	'S6PPRRRR3_jacobiaD_transl_5_sym_varpar: Position vector r_i_i_C has to be [3x1] double');
assert(isreal(pkin) && all(size(pkin) == [14 1]), ...
  'S6PPRRRR3_jacobiaD_transl_5_sym_varpar: pkin has to be [14x1] (double)');

%% Symbolic Calculation
% From jacobiaD_transl_5_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 19:43:59
% EndTime: 2019-02-26 19:43:59
% DurationCPUTime: 0.66s
% Computational Cost: add. (1005->131), mult. (3235->244), div. (0->0), fcn. (3808->16), ass. (0->82)
t510 = sin(pkin(14));
t514 = sin(pkin(6));
t515 = cos(pkin(14));
t522 = sin(qJ(3));
t525 = cos(qJ(3));
t518 = cos(pkin(7));
t545 = t518 * t522;
t513 = sin(pkin(7));
t519 = cos(pkin(6));
t551 = t513 * t519;
t496 = (t510 * t525 + t515 * t545) * t514 + t522 * t551;
t521 = sin(qJ(4));
t524 = cos(qJ(4));
t495 = t525 * t551 + (t515 * t518 * t525 - t510 * t522) * t514;
t550 = t514 * t515;
t501 = -t513 * t550 + t518 * t519;
t512 = sin(pkin(8));
t517 = cos(pkin(8));
t534 = t495 * t517 + t501 * t512;
t480 = t496 * t524 + t534 * t521;
t516 = cos(pkin(13));
t511 = sin(pkin(13));
t555 = t511 * t519;
t505 = -t510 * t555 + t515 * t516;
t504 = -t510 * t516 - t515 * t555;
t552 = t513 * t514;
t531 = t504 * t518 + t511 * t552;
t491 = t505 * t525 + t531 * t522;
t490 = -t505 * t522 + t531 * t525;
t549 = t514 * t518;
t498 = -t504 * t513 + t511 * t549;
t535 = t490 * t517 + t498 * t512;
t474 = t491 * t524 + t535 * t521;
t548 = t516 * t519;
t502 = -t511 * t510 + t515 * t548;
t540 = t516 * t552;
t503 = t510 * t548 + t511 * t515;
t557 = t503 * t525;
t489 = t502 * t545 - t522 * t540 + t557;
t532 = -t502 * t518 + t540;
t558 = t503 * t522;
t488 = -t532 * t525 - t558;
t497 = -t502 * t513 - t516 * t549;
t536 = t488 * t517 + t497 * t512;
t472 = t489 * t524 + t536 * t521;
t562 = r_i_i_C(3) + pkin(11);
t520 = sin(qJ(5));
t554 = t512 * t520;
t523 = cos(qJ(5));
t553 = t512 * t523;
t547 = t517 * t521;
t546 = t517 * t524;
t544 = qJD(3) * t522;
t543 = qJD(3) * t525;
t542 = qJD(5) * t520;
t541 = qJD(5) * t523;
t538 = t513 * t543;
t537 = t518 * t543;
t533 = r_i_i_C(1) * t523 - r_i_i_C(2) * t520 + pkin(4);
t477 = t488 * t524 - t489 * t547;
t478 = t490 * t524 - t491 * t547;
t483 = t495 * t524 - t496 * t547;
t529 = qJD(5) * (-r_i_i_C(1) * t520 - r_i_i_C(2) * t523);
t528 = -t489 * t521 + t536 * t524;
t527 = -t491 * t521 + t535 * t524;
t526 = -t496 * t521 + t534 * t524;
t494 = t496 * qJD(3);
t493 = t510 * t514 * t544 - t519 * t538 - t537 * t550;
t492 = -t495 * t512 + t501 * t517;
t487 = t491 * qJD(3);
t486 = -t511 * t514 * t538 - t504 * t537 + t505 * t544;
t485 = (t532 * t522 - t557) * qJD(3);
t484 = -t502 * t537 + (t525 * t540 + t558) * qJD(3);
t482 = -t490 * t512 + t498 * t517;
t481 = -t488 * t512 + t497 * t517;
t476 = t493 * t547 - t494 * t524 + (-t495 * t521 - t496 * t546) * qJD(4);
t470 = t526 * qJD(4) - t493 * t524 - t494 * t547;
t468 = t486 * t547 - t487 * t524 + (-t490 * t521 - t491 * t546) * qJD(4);
t466 = t484 * t547 + t485 * t524 + (-t488 * t521 - t489 * t546) * qJD(4);
t464 = t527 * qJD(4) - t486 * t524 - t487 * t547;
t462 = t528 * qJD(4) - t484 * t524 + t485 * t547;
t1 = [0, 0 (t468 * t523 - t478 * t542) * r_i_i_C(1) + (-t468 * t520 - t478 * t541) * r_i_i_C(2) + t468 * pkin(4) - t487 * pkin(3) + t562 * (t478 * qJD(4) - t486 * t546 - t487 * t521) + ((-t486 * t520 + t491 * t541) * r_i_i_C(1) + (-t486 * t523 - t491 * t542) * r_i_i_C(2) - t486 * pkin(10)) * t512, t562 * t464 + t527 * t529 + t533 * (-t474 * qJD(4) + t486 * t521 - t487 * t546) (-t464 * t520 + t487 * t553) * r_i_i_C(1) + (-t464 * t523 - t487 * t554) * r_i_i_C(2) + ((-t474 * t523 - t482 * t520) * r_i_i_C(1) + (t474 * t520 - t482 * t523) * r_i_i_C(2)) * qJD(5), 0; 0, 0 (t466 * t523 - t477 * t542) * r_i_i_C(1) + (-t466 * t520 - t477 * t541) * r_i_i_C(2) + t466 * pkin(4) + t485 * pkin(3) + t562 * (t477 * qJD(4) - t484 * t546 + t485 * t521) + ((-t484 * t520 + t489 * t541) * r_i_i_C(1) + (-t484 * t523 - t489 * t542) * r_i_i_C(2) - t484 * pkin(10)) * t512, t562 * t462 + t528 * t529 + t533 * (-t472 * qJD(4) + t484 * t521 + t485 * t546) (-t462 * t520 - t485 * t553) * r_i_i_C(1) + (-t462 * t523 + t485 * t554) * r_i_i_C(2) + ((-t472 * t523 - t481 * t520) * r_i_i_C(1) + (t472 * t520 - t481 * t523) * r_i_i_C(2)) * qJD(5), 0; 0, 0 (t476 * t523 - t483 * t542) * r_i_i_C(1) + (-t476 * t520 - t483 * t541) * r_i_i_C(2) + t476 * pkin(4) - t494 * pkin(3) + t562 * (t483 * qJD(4) - t493 * t546 - t494 * t521) + ((-t493 * t520 + t496 * t541) * r_i_i_C(1) + (-t493 * t523 - t496 * t542) * r_i_i_C(2) - t493 * pkin(10)) * t512, t562 * t470 + t526 * t529 + t533 * (-t480 * qJD(4) + t493 * t521 - t494 * t546) (-t470 * t520 + t494 * t553) * r_i_i_C(1) + (-t470 * t523 - t494 * t554) * r_i_i_C(2) + ((-t480 * t523 - t492 * t520) * r_i_i_C(1) + (t480 * t520 - t492 * t523) * r_i_i_C(2)) * qJD(5), 0;];
JaD_transl  = t1;
