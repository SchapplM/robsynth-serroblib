% Zeitableitung der translatorischen Teilmatrix der analytischen Jacobi-Matrix für Segment Nr. 4 (0=Basis) von
% S6PRRRRR6
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
%   pkin=[a2,a3,a4,a5,a6,alpha2,alpha3,alpha4,d2,d3,d4,d5,d6,theta1]';
%
% Output:
% JaD_transl [3x6]
%   Zeitableitung der translatorischen Teilmatrix der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-02-26 20:22
% Revision: d75aae1ac561373cd3be920984c3df29a1c2ecc8 (2019-02-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JaD_transl = S6PRRRRR6_jacobiaD_transl_4_sym_varpar(qJ, qJD, r_i_i_C, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(3,1),zeros(14,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRRRRR6_jacobiaD_transl_4_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6PRRRRR6_jacobiaD_transl_4_sym_varpar: qJD has to be [6x1] (double)');
assert(isa(r_i_i_C,'double') && isreal(r_i_i_C) && all(size(r_i_i_C) == [3 1]), ...
	'S6PRRRRR6_jacobiaD_transl_4_sym_varpar: Position vector r_i_i_C has to be [3x1] double');
assert(isreal(pkin) && all(size(pkin) == [14 1]), ...
  'S6PRRRRR6_jacobiaD_transl_4_sym_varpar: pkin has to be [14x1] (double)');

%% Symbolic Calculation
% From jacobiaD_transl_4_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 20:21:54
% EndTime: 2019-02-26 20:21:54
% DurationCPUTime: 0.82s
% Computational Cost: add. (594->150), mult. (2011->271), div. (0->0), fcn. (2160->14), ass. (0->89)
t509 = pkin(11) + r_i_i_C(3);
t456 = sin(qJ(3));
t459 = cos(qJ(3));
t447 = sin(pkin(14));
t451 = cos(pkin(14));
t460 = cos(qJ(2));
t454 = cos(pkin(6));
t457 = sin(qJ(2));
t493 = t454 * t457;
t466 = t447 * t493 - t451 * t460;
t453 = cos(pkin(7));
t492 = t454 * t460;
t467 = t447 * t492 + t451 * t457;
t449 = sin(pkin(7));
t450 = sin(pkin(6));
t502 = t449 * t450;
t468 = t447 * t502 - t453 * t467;
t417 = t468 * t456 - t459 * t466;
t448 = sin(pkin(8));
t452 = cos(pkin(8));
t455 = sin(qJ(4));
t458 = cos(qJ(4));
t461 = t509 * t448 - (t455 * r_i_i_C(1) + t458 * r_i_i_C(2)) * t452;
t508 = t449 * pkin(10);
t438 = t447 * t460 + t451 * t493;
t433 = t438 * qJD(2);
t507 = t433 * t456;
t506 = t438 * t459;
t504 = t448 * t455;
t503 = t448 * t458;
t501 = t449 * t452;
t500 = t449 * t454;
t499 = t450 * t453;
t498 = t450 * t460;
t497 = t452 * t455;
t496 = t452 * t458;
t495 = t453 * t456;
t494 = t453 * t459;
t491 = t456 * t457;
t490 = t456 * t460;
t489 = t457 * t459;
t488 = t459 * t460;
t487 = qJD(2) * t457;
t486 = qJD(2) * t460;
t485 = qJD(3) * t456;
t484 = qJD(3) * t459;
t483 = qJD(4) * t457;
t482 = t449 * t504;
t481 = t449 * t503;
t480 = t451 * t502;
t479 = t451 * t492;
t478 = t456 * t500;
t477 = t450 * t487;
t475 = t448 * t449 * t477;
t474 = t458 * r_i_i_C(1) - t455 * r_i_i_C(2) + pkin(3);
t437 = -t447 * t457 + t479;
t420 = -t437 * t456 - t438 * t494;
t471 = -t437 * t459 + t438 * t495;
t470 = -t437 * t453 + t480;
t422 = t456 * t467 + t466 * t494;
t469 = t459 * t467 - t466 * t495;
t465 = t453 * t488 - t491;
t464 = -t453 * t489 - t490;
t463 = t453 * t490 + t489;
t462 = t453 * t491 - t488;
t436 = -t449 * t498 + t454 * t453;
t435 = t466 * qJD(2);
t434 = t467 * qJD(2);
t432 = -qJD(2) * t479 + t447 * t487;
t429 = t462 * t450;
t428 = t464 * t450;
t427 = t447 * t499 + t449 * t467;
t426 = -t437 * t449 - t451 * t499;
t425 = t463 * t450 + t478;
t424 = t465 * t450 + t459 * t500;
t416 = t456 * t466 + t468 * t459;
t415 = t437 * t495 - t456 * t480 + t506;
t414 = -t438 * t456 - t470 * t459;
t413 = -t477 * t495 - t450 * t457 * t485 + (t450 * t486 + (t453 * t498 + t500) * qJD(3)) * t459;
t412 = -qJD(3) * t478 + (t464 * qJD(2) - t463 * qJD(3)) * t450;
t411 = t422 * qJD(3) + t434 * t495 + t435 * t459;
t410 = t469 * qJD(3) + t434 * t494 - t435 * t456;
t409 = t420 * qJD(3) + t432 * t495 - t433 * t459;
t408 = t471 * qJD(3) + t432 * t494 + t507;
t407 = t435 * t495 + t466 * t485 + (t468 * qJD(3) - t434) * t459;
t406 = -t417 * qJD(3) + t434 * t456 + t435 * t494;
t405 = -t432 * t459 - t438 * t485 - t480 * t484 + (t437 * t484 - t507) * t453;
t404 = -t433 * t494 + t432 * t456 + (t470 * t456 - t506) * qJD(3);
t1 = [0 (t410 * t497 + t411 * t458 - t434 * t482) * r_i_i_C(1) + (t410 * t496 - t411 * t455 - t434 * t481) * r_i_i_C(2) + t411 * pkin(3) + t435 * pkin(2) - t434 * t508 + ((t422 * t496 + t455 * t469 - t466 * t481) * r_i_i_C(1) + (-t422 * t497 + t458 * t469 + t466 * t482) * r_i_i_C(2)) * qJD(4) + t509 * (-t410 * t448 - t434 * t501) t474 * t406 + t461 * t407 + ((-t416 * t455 - t417 * t496) * r_i_i_C(1) + (-t416 * t458 + t417 * t497) * r_i_i_C(2)) * qJD(4) (t406 * t496 - t407 * t455 - t435 * t481) * r_i_i_C(1) + (-t406 * t497 - t407 * t458 + t435 * t482) * r_i_i_C(2) + ((-t416 * t497 - t417 * t458 - t427 * t504) * r_i_i_C(1) + (-t416 * t496 + t417 * t455 - t427 * t503) * r_i_i_C(2)) * qJD(4), 0, 0; 0 (t408 * t497 + t409 * t458 - t432 * t482) * r_i_i_C(1) + (t408 * t496 - t409 * t455 - t432 * t481) * r_i_i_C(2) + t409 * pkin(3) - t433 * pkin(2) - t432 * t508 + ((t420 * t496 + t438 * t481 + t455 * t471) * r_i_i_C(1) + (-t420 * t497 - t438 * t482 + t458 * t471) * r_i_i_C(2)) * qJD(4) + t509 * (-t408 * t448 - t432 * t501) t474 * t404 + t461 * t405 + ((-t414 * t455 - t415 * t496) * r_i_i_C(1) + (-t414 * t458 + t415 * t497) * r_i_i_C(2)) * qJD(4) (t404 * t496 - t405 * t455 + t433 * t481) * r_i_i_C(1) + (-t404 * t497 - t405 * t458 - t433 * t482) * r_i_i_C(2) + ((-t414 * t497 - t415 * t458 - t426 * t504) * r_i_i_C(1) + (-t414 * t496 + t415 * t455 - t426 * t503) * r_i_i_C(2)) * qJD(4), 0, 0; 0 ((t428 * t496 + t429 * t455) * r_i_i_C(1) + (-t428 * t497 + t429 * t458) * r_i_i_C(2)) * qJD(4) + (t474 * (-t463 * qJD(2) + t464 * qJD(3)) - t461 * (-t465 * qJD(2) + t462 * qJD(3)) - pkin(2) * t487 + ((t509 * t452 + pkin(10)) * t486 + ((t455 * t486 + t458 * t483) * r_i_i_C(1) + (-t455 * t483 + t458 * t486) * r_i_i_C(2)) * t448) * t449) * t450, t474 * t412 + t461 * t413 + ((-t424 * t455 - t425 * t496) * r_i_i_C(1) + (-t424 * t458 + t425 * t497) * r_i_i_C(2)) * qJD(4) (t412 * t496 - t413 * t455 + t458 * t475) * r_i_i_C(1) + (-t412 * t497 - t413 * t458 - t455 * t475) * r_i_i_C(2) + ((-t424 * t497 - t425 * t458 - t436 * t504) * r_i_i_C(1) + (-t424 * t496 + t425 * t455 - t436 * t503) * r_i_i_C(2)) * qJD(4), 0, 0;];
JaD_transl  = t1;
