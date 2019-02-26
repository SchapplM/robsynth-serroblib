% Zeitableitung der translatorischen Teilmatrix der analytischen Jacobi-Matrix für Segment Nr. 5 (0=Basis) von
% S6RRRRPP9
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
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d1,d2,d3,d4]';
%
% Output:
% JaD_transl [3x6]
%   Zeitableitung der translatorischen Teilmatrix der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-02-26 22:30
% Revision: d75aae1ac561373cd3be920984c3df29a1c2ecc8 (2019-02-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JaD_transl = S6RRRRPP9_jacobiaD_transl_5_sym_varpar(qJ, qJD, r_i_i_C, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(3,1),zeros(10,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRRPP9_jacobiaD_transl_5_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRRRPP9_jacobiaD_transl_5_sym_varpar: qJD has to be [6x1] (double)');
assert(isa(r_i_i_C,'double') && isreal(r_i_i_C) && all(size(r_i_i_C) == [3 1]), ...
	'S6RRRRPP9_jacobiaD_transl_5_sym_varpar: Position vector r_i_i_C has to be [3x1] double');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RRRRPP9_jacobiaD_transl_5_sym_varpar: pkin has to be [10x1] (double)');

%% Symbolic Calculation
% From jacobiaD_transl_5_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 22:30:09
% EndTime: 2019-02-26 22:30:10
% DurationCPUTime: 0.76s
% Computational Cost: add. (793->118), mult. (2339->192), div. (0->0), fcn. (2413->10), ass. (0->74)
t438 = cos(pkin(6));
t442 = sin(qJ(1));
t441 = sin(qJ(2));
t483 = t441 * t442;
t472 = t438 * t483;
t477 = qJD(2) * t441;
t445 = cos(qJ(2));
t446 = cos(qJ(1));
t479 = t445 * t446;
t417 = -qJD(1) * t472 - t442 * t477 + (qJD(2) * t438 + qJD(1)) * t479;
t481 = t442 * t445;
t482 = t441 * t446;
t429 = t438 * t482 + t481;
t440 = sin(qJ(3));
t444 = cos(qJ(3));
t437 = sin(pkin(6));
t478 = qJD(1) * t437;
t470 = t442 * t478;
t486 = t437 * t446;
t473 = t444 * t486;
t406 = (-qJD(3) * t429 + t470) * t440 - qJD(3) * t473 + t417 * t444;
t430 = t438 * t481 + t482;
t416 = t430 * qJD(1) + t429 * qJD(2);
t439 = sin(qJ(4));
t443 = cos(qJ(4));
t423 = -t429 * t444 + t440 * t486;
t428 = -t438 * t479 + t483;
t500 = t423 * t443 - t428 * t439;
t505 = qJD(4) * t500 - t406 * t439 + t416 * t443;
t501 = t423 * t439 + t428 * t443;
t504 = qJD(4) * t501 + t406 * t443 + t416 * t439;
t497 = r_i_i_C(1) + pkin(10);
t499 = t444 * pkin(3) + t497 * t440 + pkin(2);
t498 = -pkin(3) * t440 + t497 * t444;
t494 = r_i_i_C(3) + qJ(5);
t496 = r_i_i_C(2) - pkin(4);
t452 = t494 * t439 - t496 * t443 + pkin(3);
t456 = t472 - t479;
t488 = t437 * t442;
t425 = t440 * t488 - t444 * t456;
t491 = t425 * t439;
t487 = t437 * t444;
t485 = t439 * t444;
t484 = t439 * t445;
t480 = t443 * t445;
t476 = qJD(3) * t440;
t475 = qJD(3) * t445;
t474 = qJD(4) * t444;
t469 = t446 * t478;
t468 = t437 * t477;
t467 = qJD(2) * t437 * t445;
t465 = -qJD(2) + t474;
t415 = t429 * qJD(1) + t430 * qJD(2);
t464 = t430 * t474 - t415;
t463 = t428 * t474 + t417;
t460 = t425 * t443 + t430 * t439;
t427 = t438 * t440 + t441 * t487;
t459 = -t427 * t443 + t437 * t484;
t458 = t440 * t456 + t442 * t487;
t457 = -t437 * t440 * t441 + t438 * t444;
t453 = qJD(3) * t498;
t414 = t428 * qJD(1) + t456 * qJD(2);
t451 = -qJD(4) * t456 + t414 * t444 + t430 * t476;
t450 = qJD(4) * t429 - t416 * t444 + t428 * t476;
t449 = t440 * t475 + (qJD(2) * t444 - qJD(4)) * t441;
t448 = qJD(5) * t439 + (t496 * t439 + t494 * t443) * qJD(4);
t447 = t423 * qJD(3) - t417 * t440 + t444 * t470;
t420 = t457 * qJD(3) + t444 * t467;
t409 = -t459 * qJD(4) + t420 * t439 - t443 * t468;
t404 = t458 * qJD(3) - t415 * t444 + t440 * t469;
t403 = t425 * qJD(3) - t415 * t440 - t444 * t469;
t394 = -t414 * t439 - qJD(4) * t491 + (qJD(4) * t430 + t404) * t443;
t393 = t460 * qJD(4) + t404 * t439 + t414 * t443;
t1 = [t501 * qJD(5) - t406 * pkin(3) - t417 * pkin(2) - t416 * pkin(9) + t497 * t447 + t496 * t504 + t494 * t505 + (-pkin(1) * t446 - pkin(8) * t488) * qJD(1) -(t430 * t485 - t443 * t456) * qJD(5) - t415 * pkin(9) - t496 * (t464 * t439 + t451 * t443) + t494 * (t451 * t439 - t464 * t443) - t430 * t453 + t499 * t414, -t452 * t403 + t497 * t404 + t448 * t458, t460 * qJD(5) + t496 * t393 + t494 * t394, t393, 0; -(t430 * t443 - t491) * qJD(5) + t404 * pkin(3) - t415 * pkin(2) - t414 * pkin(9) + t497 * t403 - t496 * t394 + t494 * t393 + (-pkin(1) * t442 + pkin(8) * t486) * qJD(1) -(t428 * t485 + t429 * t443) * qJD(5) + t417 * pkin(9) - t496 * (t463 * t439 + t450 * t443) + t494 * (t450 * t439 - t463 * t443) - t428 * t453 - t499 * t416, t497 * t406 + t448 * (-t429 * t440 - t473) + t452 * t447, -t500 * qJD(5) + t494 * t504 - t496 * t505, -t505, 0; 0 (t496 * (t449 * t443 + t465 * t484) - t494 * (t449 * t439 - t465 * t480) - (t441 * t443 - t444 * t484) * qJD(5) + t498 * t475 + (t445 * pkin(9) - t441 * t499) * qJD(2)) * t437, t497 * t420 + t448 * t457 + t452 * (-t427 * qJD(3) - t440 * t467) -t459 * qJD(5) + t494 * (t439 * t468 + t420 * t443 + (-t427 * t439 - t437 * t480) * qJD(4)) + t496 * t409, t409, 0;];
JaD_transl  = t1;
