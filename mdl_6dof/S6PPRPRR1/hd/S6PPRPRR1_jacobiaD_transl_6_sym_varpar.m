% Zeitableitung der translatorischen Teilmatrix der analytischen Jacobi-Matrix für Segment Nr. 6 (0=Basis) von
% S6PPRPRR1
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
%   pkin=[a2,a3,a4,a5,a6,alpha2,alpha3,d3,d5,d6,theta1,theta2,theta4]';
%
% Output:
% JaD_transl [3x6]
%   Zeitableitung der translatorischen Teilmatrix der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-02-26 19:39
% Revision: d75aae1ac561373cd3be920984c3df29a1c2ecc8 (2019-02-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JaD_transl = S6PPRPRR1_jacobiaD_transl_6_sym_varpar(qJ, qJD, r_i_i_C, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(3,1),zeros(13,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PPRPRR1_jacobiaD_transl_6_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6PPRPRR1_jacobiaD_transl_6_sym_varpar: qJD has to be [6x1] (double)');
assert(isa(r_i_i_C,'double') && isreal(r_i_i_C) && all(size(r_i_i_C) == [3 1]), ...
	'S6PPRPRR1_jacobiaD_transl_6_sym_varpar: Position vector r_i_i_C has to be [3x1] double');
assert(isreal(pkin) && all(size(pkin) == [13 1]), ...
  'S6PPRPRR1_jacobiaD_transl_6_sym_varpar: pkin has to be [13x1] (double)');

%% Symbolic Calculation
% From jacobiaD_transl_6_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 19:39:51
% EndTime: 2019-02-26 19:39:51
% DurationCPUTime: 0.59s
% Computational Cost: add. (972->109), mult. (3071->210), div. (0->0), fcn. (3566->16), ass. (0->75)
t484 = sin(pkin(13));
t489 = cos(pkin(13));
t499 = cos(qJ(3));
t517 = qJD(3) * t499;
t496 = sin(qJ(3));
t518 = qJD(3) * t496;
t528 = t484 * t518 - t489 * t517;
t487 = sin(pkin(7));
t462 = t528 * t487;
t492 = cos(pkin(7));
t464 = t528 * t492;
t476 = -t484 * t517 - t489 * t518;
t485 = sin(pkin(12));
t488 = sin(pkin(6));
t490 = cos(pkin(12));
t493 = cos(pkin(6));
t453 = t493 * t462 + (t464 * t490 - t476 * t485) * t488;
t527 = pkin(10) + r_i_i_C(3);
t526 = pkin(3) * qJD(3);
t486 = sin(pkin(11));
t525 = t486 * t488;
t524 = t486 * t493;
t523 = t487 * t488;
t491 = cos(pkin(11));
t522 = t488 * t491;
t521 = t488 * t492;
t520 = t491 * t493;
t494 = sin(qJ(6));
t516 = qJD(6) * t494;
t497 = cos(qJ(6));
t515 = qJD(6) * t497;
t471 = -t486 * t485 + t490 * t520;
t457 = -t471 * t487 - t491 * t521;
t495 = sin(qJ(5));
t498 = cos(qJ(5));
t508 = t499 * t484 + t496 * t489;
t467 = t508 * t487;
t469 = t508 * t492;
t472 = t485 * t520 + t486 * t490;
t477 = t496 * t484 - t499 * t489;
t504 = -t467 * t522 + t471 * t469 - t472 * t477;
t434 = t457 * t495 + t498 * t504;
t512 = t457 * t498 - t495 * t504;
t473 = -t491 * t485 - t490 * t524;
t458 = -t473 * t487 + t486 * t521;
t474 = -t485 * t524 + t491 * t490;
t503 = t467 * t525 + t473 * t469 - t474 * t477;
t436 = t458 * t495 + t498 * t503;
t511 = t458 * t498 - t495 * t503;
t470 = -t490 * t523 + t493 * t492;
t502 = t493 * t467 + (t469 * t490 - t477 * t485) * t488;
t450 = t470 * t495 + t498 * t502;
t510 = t470 * t498 - t495 * t502;
t507 = t497 * r_i_i_C(1) - t494 * r_i_i_C(2) + pkin(5);
t506 = qJD(6) * (-t494 * r_i_i_C(1) - t497 * r_i_i_C(2));
t505 = qJD(3) * t508;
t438 = t462 * t522 - t471 * t464 + t472 * t476;
t442 = t462 * t525 + t473 * t464 - t474 * t476;
t501 = t495 * t527 + t507 * t498 + pkin(4);
t500 = t498 * t506 + (-t507 * t495 + t498 * t527) * qJD(5);
t475 = t477 * qJD(3);
t468 = t477 * t492;
t466 = t477 * t487;
t465 = t492 * t505;
t463 = t487 * t505;
t455 = -t493 * t466 + (-t468 * t490 - t485 * t508) * t488;
t451 = -t493 * t463 + (-t465 * t490 + t475 * t485) * t488;
t447 = -t466 * t525 - t473 * t468 - t474 * t508;
t444 = t466 * t522 - t471 * t468 - t472 * t508;
t440 = -t463 * t525 - t473 * t465 + t474 * t475;
t437 = t463 * t522 - t471 * t465 + t472 * t475;
t432 = qJD(5) * t510 - t453 * t498;
t430 = qJD(5) * t511 - t442 * t498;
t428 = qJD(5) * t512 + t438 * t498;
t1 = [0, 0 (-t442 * t494 + t503 * t515) * r_i_i_C(1) + (-t442 * t497 - t503 * t516) * r_i_i_C(2) - t442 * pkin(9) + (-t474 * t499 + (-t473 * t492 - t486 * t523) * t496) * t526 + t501 * t440 + t500 * t447, 0, t527 * t430 + t511 * t506 + t507 * (-qJD(5) * t436 + t442 * t495) (-t430 * t494 - t440 * t497) * r_i_i_C(1) + (-t430 * t497 + t440 * t494) * r_i_i_C(2) + ((-t436 * t497 + t447 * t494) * r_i_i_C(1) + (t436 * t494 + t447 * t497) * r_i_i_C(2)) * qJD(6); 0, 0 (t438 * t494 + t504 * t515) * r_i_i_C(1) + (t438 * t497 - t504 * t516) * r_i_i_C(2) + t438 * pkin(9) + (-t472 * t499 + (-t471 * t492 + t487 * t522) * t496) * t526 + t501 * t437 + t500 * t444, 0, t527 * t428 + t512 * t506 + t507 * (-qJD(5) * t434 - t438 * t495) (-t428 * t494 - t437 * t497) * r_i_i_C(1) + (-t428 * t497 + t437 * t494) * r_i_i_C(2) + ((-t434 * t497 + t444 * t494) * r_i_i_C(1) + (t434 * t494 + t444 * t497) * r_i_i_C(2)) * qJD(6); 0, 0 (-t453 * t494 + t502 * t515) * r_i_i_C(1) + (-t453 * t497 - t502 * t516) * r_i_i_C(2) - t453 * pkin(9) + (-t487 * t493 * t496 + (-t490 * t492 * t496 - t485 * t499) * t488) * t526 + t501 * t451 + t500 * t455, 0, t527 * t432 + t510 * t506 + t507 * (-qJD(5) * t450 + t453 * t495) (-t432 * t494 - t451 * t497) * r_i_i_C(1) + (-t432 * t497 + t451 * t494) * r_i_i_C(2) + ((-t450 * t497 + t455 * t494) * r_i_i_C(1) + (t450 * t494 + t455 * t497) * r_i_i_C(2)) * qJD(6);];
JaD_transl  = t1;
