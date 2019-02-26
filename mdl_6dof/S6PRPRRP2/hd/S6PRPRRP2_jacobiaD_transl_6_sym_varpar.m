% Zeitableitung der translatorischen Teilmatrix der analytischen Jacobi-Matrix für Segment Nr. 6 (0=Basis) von
% S6PRPRRP2
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
%   pkin=[a2,a3,a4,a5,a6,alpha2,d2,d4,d5,theta1,theta3]';
%
% Output:
% JaD_transl [3x6]
%   Zeitableitung der translatorischen Teilmatrix der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-02-26 19:51
% Revision: d75aae1ac561373cd3be920984c3df29a1c2ecc8 (2019-02-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JaD_transl = S6PRPRRP2_jacobiaD_transl_6_sym_varpar(qJ, qJD, r_i_i_C, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(3,1),zeros(11,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRPRRP2_jacobiaD_transl_6_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6PRPRRP2_jacobiaD_transl_6_sym_varpar: qJD has to be [6x1] (double)');
assert(isa(r_i_i_C,'double') && isreal(r_i_i_C) && all(size(r_i_i_C) == [3 1]), ...
	'S6PRPRRP2_jacobiaD_transl_6_sym_varpar: Position vector r_i_i_C has to be [3x1] double');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6PRPRRP2_jacobiaD_transl_6_sym_varpar: pkin has to be [11x1] (double)');

%% Symbolic Calculation
% From jacobiaD_transl_6_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 19:51:07
% EndTime: 2019-02-26 19:51:07
% DurationCPUTime: 0.53s
% Computational Cost: add. (822->101), mult. (2522->179), div. (0->0), fcn. (2797->12), ass. (0->71)
t451 = sin(qJ(4));
t454 = cos(qJ(4));
t497 = pkin(9) + r_i_i_C(2);
t463 = qJD(4) * (-pkin(4) * t451 + t454 * t497);
t449 = cos(pkin(6));
t445 = sin(pkin(11));
t452 = sin(qJ(2));
t491 = cos(pkin(11));
t495 = cos(qJ(2));
t462 = t445 * t495 + t452 * t491;
t434 = t462 * t449;
t478 = t495 * t491;
t485 = qJD(2) * t452;
t499 = -qJD(2) * t478 + t445 * t485;
t461 = -t445 * t452 + t478;
t464 = t454 * pkin(4) + t451 * t497 + pkin(3);
t496 = -r_i_i_C(1) - pkin(5);
t493 = r_i_i_C(3) + qJ(6);
t492 = pkin(2) * qJD(2);
t447 = sin(pkin(6));
t490 = t447 * t451;
t489 = t447 * t454;
t488 = t449 * t452;
t450 = sin(qJ(5));
t487 = t450 * t454;
t484 = qJD(4) * t451;
t483 = qJD(5) * t454;
t430 = t499 * t449;
t436 = t462 * qJD(2);
t446 = sin(pkin(10));
t448 = cos(pkin(10));
t413 = t430 * t448 + t436 * t446;
t456 = t461 * t449;
t418 = -t446 * t462 + t448 * t456;
t477 = t418 * t483 + t413;
t415 = t430 * t446 - t436 * t448;
t421 = -t446 * t456 - t448 * t462;
t476 = t421 * t483 - t415;
t428 = t499 * t447;
t432 = t461 * t447;
t475 = -t432 * t483 - t428;
t453 = cos(qJ(5));
t470 = t434 * t448 + t446 * t461;
t466 = t448 * t490 - t454 * t470;
t474 = -t418 * t450 - t453 * t466;
t469 = -t434 * t446 + t448 * t461;
t408 = t446 * t490 + t454 * t469;
t473 = t408 * t453 - t421 * t450;
t433 = t462 * t447;
t424 = t433 * t454 + t449 * t451;
t472 = t424 * t453 - t432 * t450;
t471 = -t433 * t451 + t449 * t454;
t467 = -t448 * t489 - t451 * t470;
t465 = t446 * t489 - t451 * t469;
t460 = t450 * t493 - t453 * t496 + pkin(4);
t431 = qJD(2) * t434;
t435 = t461 * qJD(2);
t411 = -t431 * t448 - t435 * t446;
t459 = qJD(5) * t470 + t411 * t454 - t418 * t484;
t414 = t431 * t446 - t435 * t448;
t458 = qJD(5) * t469 + t414 * t454 - t421 * t484;
t429 = qJD(2) * t433;
t457 = qJD(5) * t433 - t429 * t454 - t432 * t484;
t455 = qJD(6) * t450 + (t450 * t496 + t453 * t493) * qJD(5);
t404 = qJD(4) * t471 - t428 * t454;
t402 = qJD(4) * t465 + t415 * t454;
t400 = qJD(4) * t467 - t413 * t454;
t395 = qJD(5) * t472 + t404 * t450 - t429 * t453;
t389 = qJD(5) * t473 + t402 * t450 + t414 * t453;
t387 = qJD(5) * t474 + t400 * t450 + t411 * t453;
t1 = [0 -(-t421 * t487 + t453 * t469) * qJD(6) + t415 * pkin(8) - t496 * (-t450 * t476 + t453 * t458) + t493 * (t450 * t458 + t453 * t476) + (t446 * t488 - t448 * t495) * t492 + t421 * t463 + t464 * t414, 0, t497 * t402 + t455 * t465 + t460 * (-qJD(4) * t408 - t415 * t451) t473 * qJD(6) + t493 * (t402 * t453 - t414 * t450 + (-t408 * t450 - t421 * t453) * qJD(5)) + t496 * t389, t389; 0 -(-t418 * t487 + t453 * t470) * qJD(6) - t413 * pkin(8) - t496 * (-t450 * t477 + t453 * t459) + t493 * (t450 * t459 + t453 * t477) + (-t446 * t495 - t448 * t488) * t492 + t418 * t463 + t464 * t411, 0, t497 * t400 + t455 * t467 + t460 * (qJD(4) * t466 + t413 * t451) t474 * qJD(6) + t493 * (t400 * t453 - t411 * t450 + (-t418 * t453 + t450 * t466) * qJD(5)) + t496 * t387, t387; 0 -(-t432 * t487 + t433 * t453) * qJD(6) - t428 * pkin(8) - t447 * pkin(2) * t485 - t496 * (t450 * t475 + t453 * t457) + t493 * (t450 * t457 - t453 * t475) + t432 * t463 - t464 * t429, 0, t497 * t404 + t455 * t471 + t460 * (-t424 * qJD(4) + t428 * t451) t472 * qJD(6) + t493 * (t404 * t453 + t429 * t450 + (-t424 * t450 - t432 * t453) * qJD(5)) + t496 * t395, t395;];
JaD_transl  = t1;
