% Zeitableitung der translatorischen Teilmatrix der analytischen Jacobi-Matrix für Segment Nr. 5 (0=Basis) von
% S6RRRPRR15
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
%   pkin=[a2,a3,a4,a5,a6,alpha2,alpha3,d1,d2,d3,d5,d6]';
%
% Output:
% JaD_transl [3x6]
%   Zeitableitung der translatorischen Teilmatrix der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-02-26 22:24
% Revision: d75aae1ac561373cd3be920984c3df29a1c2ecc8 (2019-02-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JaD_transl = S6RRRPRR15_jacobiaD_transl_5_sym_varpar(qJ, qJD, r_i_i_C, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(3,1),zeros(12,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRPRR15_jacobiaD_transl_5_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRRPRR15_jacobiaD_transl_5_sym_varpar: qJD has to be [6x1] (double)');
assert(isa(r_i_i_C,'double') && isreal(r_i_i_C) && all(size(r_i_i_C) == [3 1]), ...
	'S6RRRPRR15_jacobiaD_transl_5_sym_varpar: Position vector r_i_i_C has to be [3x1] double');
assert(isreal(pkin) && all(size(pkin) == [12 1]), ...
  'S6RRRPRR15_jacobiaD_transl_5_sym_varpar: pkin has to be [12x1] (double)');

%% Symbolic Calculation
% From jacobiaD_transl_5_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 22:24:27
% EndTime: 2019-02-26 22:24:28
% DurationCPUTime: 0.86s
% Computational Cost: add. (901->146), mult. (2800->248), div. (0->0), fcn. (2920->12), ass. (0->91)
t483 = sin(qJ(2));
t484 = sin(qJ(1));
t487 = cos(qJ(2));
t488 = cos(qJ(1));
t540 = cos(pkin(6));
t512 = t488 * t540;
t465 = t483 * t512 + t484 * t487;
t513 = t484 * t540;
t492 = t488 * t483 + t487 * t513;
t454 = t492 * qJD(1) + t465 * qJD(2);
t507 = t483 * t513;
t523 = qJD(2) * t483;
t525 = t488 * t487;
t455 = -qJD(1) * t507 - t484 * t523 + (qJD(2) * t540 + qJD(1)) * t525;
t464 = t484 * t483 - t487 * t512;
t482 = sin(qJ(3));
t486 = cos(qJ(3));
t478 = sin(pkin(7));
t479 = sin(pkin(6));
t524 = qJD(1) * t479;
t515 = t484 * t524;
t510 = t478 * t515;
t532 = t479 * t488;
t518 = t478 * t532;
t522 = qJD(3) * t482;
t480 = cos(pkin(7));
t530 = t480 * t486;
t531 = t480 * t482;
t535 = t465 * t486;
t430 = (t464 * t531 - t535) * qJD(3) - t454 * t530 - t455 * t482 + t486 * t510 + t518 * t522;
t538 = t454 * t478;
t445 = t480 * t515 + t538;
t481 = sin(qJ(5));
t485 = cos(qJ(5));
t549 = t430 * t481 - t445 * t485;
t548 = t430 * t485 + t445 * t481;
t500 = t464 * t480 + t518;
t536 = t465 * t482;
t545 = (t500 * t486 + t536) * qJD(3) - (-t454 * t480 + t510) * t482 - t455 * t486;
t438 = t464 * t530 + t486 * t518 + t536;
t458 = -t464 * t478 + t480 * t532;
t544 = -t438 * t485 - t458 * t481;
t543 = t438 * t481 - t458 * t485;
t493 = t507 - t525;
t452 = t464 * qJD(1) + t493 * qJD(2);
t539 = t452 * t478;
t534 = t478 * t479;
t533 = t479 * t484;
t529 = t482 * t483;
t528 = t482 * t487;
t527 = t483 * t486;
t526 = t486 * t487;
t521 = qJD(5) * t481;
t520 = qJD(5) * t485;
t519 = r_i_i_C(3) + pkin(11) + pkin(3);
t517 = t480 * t526;
t516 = pkin(10) * t480 + pkin(9);
t514 = t488 * t524;
t511 = t540 * t478;
t509 = t478 * t514;
t508 = t523 * t534;
t506 = qJD(3) * t511;
t505 = t485 * r_i_i_C(1) - t481 * r_i_i_C(2);
t504 = -t481 * r_i_i_C(1) - t485 * r_i_i_C(2);
t503 = qJ(4) - t504;
t502 = pkin(4) + pkin(10) + t505;
t450 = -t464 * t482 + t465 * t530;
t451 = -t482 * t492 - t493 * t530;
t499 = t478 * t533 - t480 * t492;
t498 = t517 - t529;
t497 = t480 * t527 + t528;
t496 = t480 * t528 + t527;
t495 = qJD(5) * t504;
t491 = t505 * qJD(5) + qJD(4);
t489 = t499 * t482 - t486 * t493;
t463 = t540 * t480 - t487 * t534;
t461 = t497 * t479;
t460 = t478 * t492 + t480 * t533;
t456 = -t498 * t479 - t486 * t511;
t453 = t465 * qJD(1) + t492 * qJD(2);
t446 = (-qJD(2) * t517 - qJD(3) * t526 + (qJD(3) * t480 + qJD(2)) * t529) * t479;
t443 = t480 * t514 - t539;
t441 = -t482 * t493 - t499 * t486;
t436 = t482 * t506 + (t497 * qJD(2) + t496 * qJD(3)) * t479;
t434 = t455 * t530 - t454 * t482 + (-t464 * t486 - t465 * t531) * qJD(3);
t432 = -t453 * t530 + t452 * t482 + (-t486 * t492 + t493 * t531) * qJD(3);
t427 = t493 * t522 + (t452 * t480 + t509) * t482 + (t499 * qJD(3) - t453) * t486;
t426 = t489 * qJD(3) - t452 * t530 - t453 * t482 - t486 * t509;
t425 = t426 * t481 + t443 * t485 + (t441 * t485 - t460 * t481) * qJD(5);
t424 = t426 * t485 - t443 * t481 + (-t441 * t481 - t460 * t485) * qJD(5);
t1 = [t549 * r_i_i_C(1) + t548 * r_i_i_C(2) - t445 * pkin(4) + t430 * qJ(4) - t438 * qJD(4) - t455 * pkin(2) - pkin(10) * t538 + (t544 * r_i_i_C(1) + t543 * r_i_i_C(2)) * qJD(5) + t519 * t545 + (-t488 * pkin(1) - t516 * t533) * qJD(1) (t432 * t481 + t451 * t520) * r_i_i_C(1) + (t432 * t485 - t451 * t521) * r_i_i_C(2) + t432 * qJ(4) + t451 * qJD(4) + t452 * pkin(2) + t519 * (-t451 * qJD(3) + t452 * t486 + t453 * t531) + (-t502 * t453 - t493 * t495) * t478, -t519 * t426 + t503 * t427 + t491 * t489, t426, t424 * r_i_i_C(1) - t425 * r_i_i_C(2), 0; -pkin(10) * t539 - t453 * pkin(2) + t443 * pkin(4) + t425 * r_i_i_C(1) + t424 * r_i_i_C(2) + t426 * qJ(4) + t441 * qJD(4) + t519 * t427 + (-pkin(1) * t484 + t516 * t532) * qJD(1) (t434 * t481 + t450 * t520) * r_i_i_C(1) + (t434 * t485 - t450 * t521) * r_i_i_C(2) + t434 * qJ(4) + t450 * qJD(4) - t454 * pkin(2) + t519 * (-t450 * qJD(3) - t454 * t486 - t455 * t531) + (t502 * t455 + t465 * t495) * t478, t491 * (-t500 * t482 + t535) - t503 * t545 + t519 * t430, -t430, -t548 * r_i_i_C(1) + t549 * r_i_i_C(2) + (-t543 * r_i_i_C(1) + t544 * r_i_i_C(2)) * qJD(5), 0; 0 (-t446 * t481 + t461 * t520) * r_i_i_C(1) + (-t446 * t485 - t461 * t521) * r_i_i_C(2) - t446 * qJ(4) + t461 * qJD(4) + (-t519 * (t496 * qJD(2) + t497 * qJD(3)) - pkin(2) * t523 + (t502 * t487 * qJD(2) + t483 * t495) * t478) * t479, t491 * (t496 * t479 + t482 * t511) + t503 * (t486 * t506 + (t498 * qJD(3) + (-t480 * t529 + t526) * qJD(2)) * t479) - t519 * t436, t436 (t436 * t485 - t481 * t508) * r_i_i_C(1) + (-t436 * t481 - t485 * t508) * r_i_i_C(2) + ((-t456 * t481 - t463 * t485) * r_i_i_C(1) + (-t456 * t485 + t463 * t481) * r_i_i_C(2)) * qJD(5), 0;];
JaD_transl  = t1;
