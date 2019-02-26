% Zeitableitung der translatorischen Teilmatrix der analytischen Jacobi-Matrix für Segment Nr. 6 (0=Basis) von
% S6RRPRRR5
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
% pkin [12x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d1,d2,d4,d5,d6,theta3]';
%
% Output:
% JaD_transl [3x6]
%   Zeitableitung der translatorischen Teilmatrix der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-02-26 21:56
% Revision: d75aae1ac561373cd3be920984c3df29a1c2ecc8 (2019-02-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JaD_transl = S6RRPRRR5_jacobiaD_transl_6_sym_varpar(qJ, qJD, r_i_i_C, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(3,1),zeros(12,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPRRR5_jacobiaD_transl_6_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRPRRR5_jacobiaD_transl_6_sym_varpar: qJD has to be [6x1] (double)');
assert(isa(r_i_i_C,'double') && isreal(r_i_i_C) && all(size(r_i_i_C) == [3 1]), ...
	'S6RRPRRR5_jacobiaD_transl_6_sym_varpar: Position vector r_i_i_C has to be [3x1] double');
assert(isreal(pkin) && all(size(pkin) == [12 1]), ...
  'S6RRPRRR5_jacobiaD_transl_6_sym_varpar: pkin has to be [12x1] (double)');

%% Symbolic Calculation
% From jacobiaD_transl_6_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 21:56:22
% EndTime: 2019-02-26 21:56:23
% DurationCPUTime: 1.08s
% Computational Cost: add. (1327->143), mult. (3394->238), div. (0->0), fcn. (3744->14), ass. (0->88)
t491 = cos(pkin(6));
t489 = sin(pkin(12));
t555 = cos(pkin(12));
t558 = cos(qJ(2));
t519 = t558 * t555;
t494 = sin(qJ(2));
t541 = qJD(2) * t494;
t562 = -qJD(2) * t519 + t489 * t541;
t464 = t562 * t491;
t530 = t494 * t555;
t509 = t558 * t489 + t530;
t469 = t509 * t491;
t531 = qJD(2) * t558;
t471 = -qJD(2) * t530 - t489 * t531;
t495 = sin(qJ(1));
t498 = cos(qJ(1));
t508 = -t494 * t489 + t519;
t543 = qJD(1) * t495;
t440 = t469 * t543 - t495 * t471 + (-qJD(1) * t508 + t464) * t498;
t493 = sin(qJ(4));
t497 = cos(qJ(4));
t516 = t469 * t498 + t495 * t508;
t490 = sin(pkin(6));
t534 = t490 * t543;
t550 = t490 * t498;
t537 = t497 * t550;
t434 = (-qJD(4) * t516 + t534) * t493 - qJD(4) * t537 - t440 * t497;
t507 = t491 * t508;
t451 = -t495 * t509 + t498 * t507;
t487 = qJD(5) + qJD(6);
t567 = t451 * t487 - t434;
t496 = cos(qJ(5));
t483 = pkin(5) * t496 + pkin(4);
t488 = qJ(5) + qJ(6);
t485 = sin(t488);
t486 = cos(t488);
t514 = r_i_i_C(1) * t486 - r_i_i_C(2) * t485 + t483;
t556 = r_i_i_C(3) + pkin(11) + pkin(10);
t492 = sin(qJ(5));
t561 = qJD(5) * t492 * pkin(5) + (r_i_i_C(1) * t485 + r_i_i_C(2) * t486) * t487;
t566 = t497 * t561 + (t514 * t493 - t556 * t497) * qJD(4);
t465 = qJD(2) * t469;
t542 = qJD(1) * t498;
t563 = qJD(1) * t507 + qJD(2) * t508;
t441 = -t498 * t465 - t495 * t563 - t509 * t542;
t447 = t493 * t550 - t497 * t516;
t564 = -t447 * t487 + t441;
t559 = t556 * t493 + t514 * t497 + pkin(3);
t557 = pkin(2) * t491;
t553 = t485 * t487;
t552 = t486 * t487;
t551 = t490 * t495;
t548 = t494 * t495;
t547 = t494 * t498;
t438 = -t495 * t465 + t498 * t563 - t509 * t543;
t515 = -t469 * t495 + t498 * t508;
t449 = t493 * t551 + t497 * t515;
t526 = -t449 * t487 + t438;
t502 = -t516 * qJD(1) + t495 * t464 + t498 * t471;
t512 = -t493 * t515 + t497 * t551;
t533 = t490 * t542;
t432 = t512 * qJD(4) + t493 * t533 + t497 * t502;
t454 = -t495 * t507 - t498 * t509;
t529 = t454 * t487 - t432;
t427 = t529 * t485 + t526 * t486;
t428 = t526 * t485 - t529 * t486;
t546 = t427 * r_i_i_C(1) - t428 * r_i_i_C(2);
t545 = (t485 * t567 - t486 * t564) * r_i_i_C(1) + (t485 * t564 + t486 * t567) * r_i_i_C(2);
t468 = t509 * t490;
t457 = t468 * t497 + t491 * t493;
t463 = qJD(2) * t468;
t521 = t457 * t487 - t463;
t462 = t562 * t490;
t517 = -t468 * t493 + t491 * t497;
t444 = t517 * qJD(4) - t462 * t497;
t467 = t508 * t490;
t522 = t467 * t487 - t444;
t544 = (t522 * t485 - t521 * t486) * r_i_i_C(1) + (t521 * t485 + t522 * t486) * r_i_i_C(2);
t540 = qJD(5) * t496;
t539 = pkin(2) * t541;
t536 = t558 * t495;
t535 = t558 * t498;
t501 = t447 * qJD(4) + t440 * t493 + t497 * t534;
t484 = t558 * pkin(2) + pkin(1);
t472 = -t490 * qJD(3) + t531 * t557;
t470 = t494 * t557 + (-pkin(8) - qJ(3)) * t490;
t431 = t449 * qJD(4) + t493 * t502 - t497 * t533;
t1 = [t495 * t539 + t440 * pkin(3) + t441 * pkin(9) - t434 * t483 - t498 * t472 + (r_i_i_C(1) * t567 + r_i_i_C(2) * t564) * t486 + (r_i_i_C(1) * t564 - r_i_i_C(2) * t567) * t485 + t556 * t501 + (t470 * t495 - t484 * t498) * qJD(1) + (t441 * t492 + (-t447 * t492 + t451 * t496) * qJD(5)) * pkin(5) (t485 * t502 + t515 * t552) * r_i_i_C(1) + (t486 * t502 - t515 * t553) * r_i_i_C(2) + t502 * pkin(9) + (t492 * t502 + t515 * t540) * pkin(5) - t559 * t438 + ((t491 * t548 - t535) * qJD(2) + (-t491 * t535 + t548) * qJD(1)) * pkin(2) - t566 * t454, t533, -t514 * t431 + t556 * t432 - t512 * t561 (-t432 * t492 + t438 * t496 + (-t449 * t496 + t454 * t492) * qJD(5)) * pkin(5) + t546, t546; -t498 * t539 + t502 * pkin(3) + t438 * pkin(9) + t428 * r_i_i_C(1) + t427 * r_i_i_C(2) + t432 * t483 - t495 * t472 + t556 * t431 + (-t470 * t498 - t484 * t495) * qJD(1) + (t438 * t492 + (-t449 * t492 - t454 * t496) * qJD(5)) * pkin(5) (-t440 * t485 + t516 * t552) * r_i_i_C(1) + (-t440 * t486 - t516 * t553) * r_i_i_C(2) - t440 * pkin(9) + (-t440 * t492 + t516 * t540) * pkin(5) + t559 * t441 + ((-t491 * t547 - t536) * qJD(2) + (-t491 * t536 - t547) * qJD(1)) * pkin(2) - t566 * t451, t534, t556 * t434 - t561 * (-t493 * t516 - t537) + t514 * t501 (-t434 * t492 - t441 * t496 + (t447 * t496 + t451 * t492) * qJD(5)) * pkin(5) + t545, t545; 0 (-t462 * t485 + t468 * t552) * r_i_i_C(1) + (-t462 * t486 - t468 * t553) * r_i_i_C(2) - t462 * pkin(9) - t490 * t539 + (-t462 * t492 + t468 * t540) * pkin(5) - t559 * t463 - t566 * t467, 0, t556 * t444 - t561 * t517 + t514 * (-t457 * qJD(4) + t462 * t493) (-t444 * t492 + t463 * t496 + (-t457 * t496 + t467 * t492) * qJD(5)) * pkin(5) + t544, t544;];
JaD_transl  = t1;
