% Zeitableitung der translatorischen Teilmatrix der analytischen Jacobi-Matrix für Segment Nr. 4 (0=Basis) von
% S6RRRPRR13
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
%   pkin=[a2,a3,a4,a5,a6,alpha2,alpha3,d1,d2,d3,d5,d6,theta4]';
%
% Output:
% JaD_transl [3x6]
%   Zeitableitung der translatorischen Teilmatrix der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-02-26 22:23
% Revision: d75aae1ac561373cd3be920984c3df29a1c2ecc8 (2019-02-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JaD_transl = S6RRRPRR13_jacobiaD_transl_4_sym_varpar(qJ, qJD, r_i_i_C, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(3,1),zeros(13,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRPRR13_jacobiaD_transl_4_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRRPRR13_jacobiaD_transl_4_sym_varpar: qJD has to be [6x1] (double)');
assert(isa(r_i_i_C,'double') && isreal(r_i_i_C) && all(size(r_i_i_C) == [3 1]), ...
	'S6RRRPRR13_jacobiaD_transl_4_sym_varpar: Position vector r_i_i_C has to be [3x1] double');
assert(isreal(pkin) && all(size(pkin) == [13 1]), ...
  'S6RRRPRR13_jacobiaD_transl_4_sym_varpar: pkin has to be [13x1] (double)');

%% Symbolic Calculation
% From jacobiaD_transl_4_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 22:23:07
% EndTime: 2019-02-26 22:23:08
% DurationCPUTime: 0.63s
% Computational Cost: add. (578->109), mult. (1811->183), div. (0->0), fcn. (1844->12), ass. (0->72)
t446 = sin(qJ(2));
t447 = sin(qJ(1));
t497 = cos(pkin(6));
t474 = t447 * t497;
t468 = t446 * t474;
t449 = cos(qJ(2));
t450 = cos(qJ(1));
t480 = t450 * t449;
t482 = t447 * t446;
t427 = -qJD(1) * t468 - qJD(2) * t482 + (qJD(2) * t497 + qJD(1)) * t480;
t441 = sin(pkin(7));
t442 = sin(pkin(6));
t488 = t442 * t450;
t478 = t441 * t488;
t506 = -qJD(3) * t478 + t427;
t498 = r_i_i_C(3) + qJ(4);
t440 = sin(pkin(13));
t443 = cos(pkin(13));
t505 = -t443 * r_i_i_C(1) + t440 * r_i_i_C(2) - pkin(3);
t445 = sin(qJ(3));
t448 = cos(qJ(3));
t455 = t468 - t480;
t444 = cos(pkin(7));
t454 = t450 * t446 + t449 * t474;
t489 = t442 * t447;
t460 = t441 * t489 - t444 * t454;
t504 = t460 * t445 - t448 * t455;
t473 = t450 * t497;
t429 = t446 * t473 + t447 * t449;
t426 = t454 * qJD(1) + t429 * qJD(2);
t428 = -t449 * t473 + t482;
t479 = qJD(1) * t442;
t476 = t447 * t479;
t471 = t441 * t476;
t486 = t444 * t448;
t487 = t444 * t445;
t493 = t429 * t448;
t503 = (t428 * t487 - t493) * qJD(3) - t426 * t486 + t448 * t471 - t506 * t445;
t494 = t428 * t444;
t502 = (-qJD(3) * t429 - t426 * t444 + t471) * t445 + (-qJD(3) * t494 + t506) * t448;
t500 = t441 * pkin(10);
t424 = t428 * qJD(1) + t455 * qJD(2);
t496 = t424 * t441;
t495 = t426 * t441;
t491 = t440 * t441;
t490 = t441 * t443;
t485 = t445 * t446;
t484 = t445 * t449;
t483 = t446 * t448;
t481 = t448 * t449;
t477 = pkin(10) * t444 + pkin(9);
t475 = t450 * t479;
t472 = t497 * t441;
t470 = t441 * t475;
t467 = qJD(3) * t472;
t464 = t428 * t445 - t429 * t486;
t462 = t478 + t494;
t461 = t445 * t454 + t455 * t486;
t459 = t444 * t481 - t485;
t458 = t444 * t483 + t484;
t457 = t444 * t484 + t483;
t456 = t444 * t485 - t481;
t451 = t445 * t455 + t460 * t448;
t425 = t429 * qJD(1) + t454 * qJD(2);
t419 = -t444 * t476 - t495;
t418 = t444 * t475 - t496;
t416 = t445 * t467 + (t458 * qJD(2) + t457 * qJD(3)) * t442;
t415 = t464 * qJD(3) - t426 * t448 - t427 * t487;
t413 = t461 * qJD(3) + t424 * t448 + t425 * t487;
t407 = -t425 * t448 + (t424 * t444 + t470) * t445 + t451 * qJD(3);
t406 = t504 * qJD(3) - t424 * t486 - t425 * t445 - t448 * t470;
t1 = [(t419 * t440 - t443 * t502) * r_i_i_C(1) + (t419 * t443 + t440 * t502) * r_i_i_C(2) - t502 * pkin(3) - (t429 * t445 + t462 * t448) * qJD(4) - t427 * pkin(2) - pkin(10) * t495 + t498 * t503 + (-t450 * pkin(1) - t477 * t489) * qJD(1) (t413 * t443 - t425 * t491) * r_i_i_C(1) + (-t413 * t440 - t425 * t490) * r_i_i_C(2) + t413 * pkin(3) - t461 * qJD(4) + t424 * pkin(2) - t425 * t500 + t498 * (-t425 * t486 + t424 * t445 + (-t448 * t454 + t455 * t487) * qJD(3)) t504 * qJD(4) + t406 * t505 + t498 * t407, t406, 0, 0; (t407 * t443 + t418 * t440) * r_i_i_C(1) + (-t407 * t440 + t418 * t443) * r_i_i_C(2) + t407 * pkin(3) - t451 * qJD(4) - t425 * pkin(2) - pkin(10) * t496 + t498 * t406 + (-t447 * pkin(1) + t477 * t488) * qJD(1) (t415 * t443 + t427 * t491) * r_i_i_C(1) + (-t415 * t440 + t427 * t490) * r_i_i_C(2) + t415 * pkin(3) - t464 * qJD(4) - t426 * pkin(2) + t427 * t500 + t498 * (t427 * t486 - t426 * t445 + (-t428 * t448 - t429 * t487) * qJD(3)) -(t462 * t445 - t493) * qJD(4) + t498 * t502 - t505 * t503, -t503, 0, 0; 0 (-t498 * (-t459 * qJD(2) + t456 * qJD(3)) - t505 * (-t457 * qJD(2) - t458 * qJD(3)) + t458 * qJD(4) + (-t446 * pkin(2) + (r_i_i_C(1) * t440 + r_i_i_C(2) * t443 + pkin(10)) * t449 * t441) * qJD(2)) * t442 -(-t457 * t442 - t445 * t472) * qJD(4) + t498 * (t448 * t467 + (-t456 * qJD(2) + t459 * qJD(3)) * t442) + t505 * t416, t416, 0, 0;];
JaD_transl  = t1;
