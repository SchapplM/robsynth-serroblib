% Zeitableitung der translatorischen Teilmatrix der analytischen Jacobi-Matrix für Segment Nr. 5 (0=Basis) von
% S6PRRPRR3
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
%   pkin=[a2,a3,a4,a5,a6,alpha2,alpha3,d2,d3,d5,d6,theta1,theta4]';
%
% Output:
% JaD_transl [3x6]
%   Zeitableitung der translatorischen Teilmatrix der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-02-26 20:05
% Revision: d75aae1ac561373cd3be920984c3df29a1c2ecc8 (2019-02-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JaD_transl = S6PRRPRR3_jacobiaD_transl_5_sym_varpar(qJ, qJD, r_i_i_C, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(3,1),zeros(13,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRRPRR3_jacobiaD_transl_5_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6PRRPRR3_jacobiaD_transl_5_sym_varpar: qJD has to be [6x1] (double)');
assert(isa(r_i_i_C,'double') && isreal(r_i_i_C) && all(size(r_i_i_C) == [3 1]), ...
	'S6PRRPRR3_jacobiaD_transl_5_sym_varpar: Position vector r_i_i_C has to be [3x1] double');
assert(isreal(pkin) && all(size(pkin) == [13 1]), ...
  'S6PRRPRR3_jacobiaD_transl_5_sym_varpar: pkin has to be [13x1] (double)');

%% Symbolic Calculation
% From jacobiaD_transl_5_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 20:05:26
% EndTime: 2019-02-26 20:05:27
% DurationCPUTime: 0.98s
% Computational Cost: add. (642->156), mult. (2061->282), div. (0->0), fcn. (2236->14), ass. (0->84)
t449 = sin(pkin(7));
t447 = sin(pkin(13));
t496 = cos(pkin(13));
t498 = cos(qJ(3));
t472 = t498 * t496;
t455 = sin(qJ(3));
t483 = qJD(3) * t455;
t501 = -qJD(3) * t472 + t447 * t483;
t416 = t501 * t449;
t452 = cos(pkin(7));
t418 = t501 * t452;
t474 = t455 * t496;
t475 = t498 * qJD(3);
t435 = -qJD(3) * t474 - t447 * t475;
t450 = sin(pkin(6));
t453 = cos(pkin(6));
t456 = sin(qJ(2));
t458 = cos(qJ(2));
t462 = t498 * t447 + t474;
t423 = t462 * t452;
t461 = -t455 * t447 + t472;
t468 = t423 * t456 - t458 * t461;
t502 = t453 * t416 + (t468 * qJD(2) + t418 * t458 - t435 * t456) * t450;
t499 = r_i_i_C(3) + pkin(10);
t497 = pkin(3) * t452;
t448 = sin(pkin(12));
t495 = t448 * t450;
t494 = t449 * t450;
t454 = sin(qJ(5));
t493 = t449 * t454;
t457 = cos(qJ(5));
t492 = t449 * t457;
t451 = cos(pkin(12));
t491 = t450 * t451;
t490 = t450 * t452;
t488 = t453 * t456;
t487 = t453 * t458;
t486 = t455 * t458;
t485 = qJD(2) * t456;
t484 = qJD(2) * t458;
t482 = qJD(5) * t454;
t481 = qJD(5) * t457;
t480 = pkin(3) * t483;
t479 = t451 * t487;
t478 = t452 * t498;
t477 = t498 * t456;
t473 = t485 * t494;
t422 = t461 * t452;
t470 = t422 * t458 - t456 * t462;
t469 = t423 * t458 + t456 * t461;
t467 = t457 * r_i_i_C(1) - t454 * r_i_i_C(2) + pkin(4);
t431 = t448 * t458 + t451 * t488;
t465 = t448 * t487 + t451 * t456;
t464 = t448 * t488 - t451 * t458;
t463 = qJD(5) * (-t454 * r_i_i_C(1) - t457 * r_i_i_C(2));
t460 = qJD(3) * t462;
t424 = -qJD(2) * t479 + t448 * t485;
t425 = t431 * qJD(2);
t430 = -t448 * t456 + t479;
t392 = t416 * t491 - t430 * t418 - t425 * t423 - t424 * t461 + t431 * t435;
t426 = t465 * qJD(2);
t427 = t464 * qJD(2);
t459 = t416 * t495 - t418 * t465 - t427 * t423 + t426 * t461 + t435 * t464;
t446 = t498 * pkin(3) + pkin(2);
t436 = -t449 * qJD(4) + t475 * t497;
t434 = t461 * qJD(3);
t429 = t453 * t452 - t458 * t494;
t428 = t455 * t497 + (-pkin(9) - qJ(4)) * t449;
t421 = t462 * t449;
t420 = t461 * t449;
t419 = t452 * t460;
t417 = t449 * t460;
t415 = t448 * t490 + t449 * t465;
t414 = -t430 * t449 - t451 * t490;
t413 = t468 * t450;
t412 = t423 * t464 - t461 * t465;
t411 = -t431 * t423 + t430 * t461;
t410 = t453 * t421 + t469 * t450;
t408 = t421 * t495 - t423 * t465 - t461 * t464;
t406 = -t421 * t491 + t430 * t423 + t431 * t461;
t404 = (-t469 * qJD(2) + t418 * t456 + t435 * t458) * t450;
t399 = -t418 * t464 + t426 * t423 + t427 * t461 - t435 * t465;
t397 = t431 * t418 + t424 * t423 - t425 * t461 + t430 * t435;
t1 = [0 (t399 * t457 - t426 * t493) * r_i_i_C(1) + (-t399 * t454 - t426 * t492) * r_i_i_C(2) + t399 * pkin(4) + t427 * t446 + t465 * t480 + t426 * t428 + t464 * t436 - t499 * (-t419 * t464 + t426 * t422 - t427 * t462 + t434 * t465) + ((-t412 * t454 - t464 * t492) * r_i_i_C(1) + (-t412 * t457 + t464 * t493) * r_i_i_C(2)) * qJD(5), -t499 * t459 + (t420 * t495 - t422 * t465 + t462 * t464) * t463 + t467 * (-t417 * t495 + t419 * t465 + t427 * t422 + t426 * t462 + t434 * t464) + (t427 * t478 + t426 * t455 + (t498 * t464 + (-t448 * t494 + t452 * t465) * t455) * qJD(3)) * pkin(3), -t427 * t449 (-t427 * t492 + t454 * t459) * r_i_i_C(1) + (t427 * t493 + t457 * t459) * r_i_i_C(2) + ((-t408 * t457 - t415 * t454) * r_i_i_C(1) + (t408 * t454 - t415 * t457) * r_i_i_C(2)) * qJD(5), 0; 0 (t397 * t457 - t424 * t493) * r_i_i_C(1) + (-t397 * t454 - t424 * t492) * r_i_i_C(2) + t397 * pkin(4) - t425 * t446 - t430 * t480 + t424 * t428 - t431 * t436 - t499 * (t431 * t419 + t424 * t422 + t425 * t462 - t430 * t434) + ((-t411 * t454 + t431 * t492) * r_i_i_C(1) + (-t411 * t457 - t431 * t493) * r_i_i_C(2)) * qJD(5), t499 * t392 + (-t420 * t491 + t430 * t422 - t431 * t462) * t463 + t467 * (t417 * t491 - t430 * t419 - t425 * t422 + t424 * t462 - t431 * t434) + (-t425 * t478 + t424 * t455 + (-t498 * t431 + (-t430 * t452 + t449 * t491) * t455) * qJD(3)) * pkin(3), t425 * t449 (-t392 * t454 + t425 * t492) * r_i_i_C(1) + (-t392 * t457 - t425 * t493) * r_i_i_C(2) + ((-t406 * t457 - t414 * t454) * r_i_i_C(1) + (t406 * t454 - t414 * t457) * r_i_i_C(2)) * qJD(5), 0; 0 (t404 * t457 + t413 * t482) * r_i_i_C(1) + (-t404 * t454 + t413 * t481) * r_i_i_C(2) + t404 * pkin(4) + (t499 * (t470 * qJD(2) - t419 * t456 + t434 * t458) - t458 * t480 - t456 * t436 + (-t458 * t428 - t456 * t446) * qJD(2) + ((t454 * t484 + t456 * t481) * r_i_i_C(1) + (-t456 * t482 + t457 * t484) * r_i_i_C(2)) * t449) * t450, -t499 * t502 + (t453 * t420 + t470 * t450) * t463 + t467 * (-t453 * t417 + (-t419 * t458 - t434 * t456 + (-t422 * t456 - t458 * t462) * qJD(2)) * t450) + (-t449 * t453 * t483 + ((-t452 * t486 - t477) * qJD(3) + (-t452 * t477 - t486) * qJD(2)) * t450) * pkin(3), t473 (t454 * t502 + t457 * t473) * r_i_i_C(1) + (-t454 * t473 + t457 * t502) * r_i_i_C(2) + ((-t410 * t457 - t429 * t454) * r_i_i_C(1) + (t410 * t454 - t429 * t457) * r_i_i_C(2)) * qJD(5), 0;];
JaD_transl  = t1;
