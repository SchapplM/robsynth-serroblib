% Zeitableitung der translatorischen Teilmatrix der analytischen Jacobi-Matrix für Segment Nr. 6 (0=Basis) von
% S6RRRRPR7
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
%   pkin=[a2,a3,a4,a5,a6,alpha2,d1,d2,d3,d4,d6,theta5]';
%
% Output:
% JaD_transl [3x6]
%   Zeitableitung der translatorischen Teilmatrix der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-02-26 22:34
% Revision: d75aae1ac561373cd3be920984c3df29a1c2ecc8 (2019-02-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JaD_transl = S6RRRRPR7_jacobiaD_transl_6_sym_varpar(qJ, qJD, r_i_i_C, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(3,1),zeros(12,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRRPR7_jacobiaD_transl_6_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRRRPR7_jacobiaD_transl_6_sym_varpar: qJD has to be [6x1] (double)');
assert(isa(r_i_i_C,'double') && isreal(r_i_i_C) && all(size(r_i_i_C) == [3 1]), ...
	'S6RRRRPR7_jacobiaD_transl_6_sym_varpar: Position vector r_i_i_C has to be [3x1] double');
assert(isreal(pkin) && all(size(pkin) == [12 1]), ...
  'S6RRRRPR7_jacobiaD_transl_6_sym_varpar: pkin has to be [12x1] (double)');

%% Symbolic Calculation
% From jacobiaD_transl_6_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 22:34:09
% EndTime: 2019-02-26 22:34:10
% DurationCPUTime: 0.79s
% Computational Cost: add. (1324->141), mult. (1913->216), div. (0->0), fcn. (1883->14), ass. (0->93)
t508 = pkin(11) + r_i_i_C(3);
t449 = sin(qJ(6));
t453 = cos(qJ(6));
t467 = t453 * r_i_i_C(1) - t449 * r_i_i_C(2);
t464 = pkin(5) + t467;
t480 = qJD(6) * t453;
t481 = qJD(6) * t449;
t517 = -r_i_i_C(1) * t481 - t480 * r_i_i_C(2);
t448 = cos(pkin(6));
t451 = sin(qJ(2));
t456 = cos(qJ(1));
t488 = t456 * t451;
t452 = sin(qJ(1));
t455 = cos(qJ(2));
t489 = t452 * t455;
t425 = t448 * t488 + t489;
t446 = qJ(3) + qJ(4);
t440 = pkin(12) + t446;
t437 = sin(t440);
t438 = cos(t440);
t447 = sin(pkin(6));
t491 = t447 * t456;
t413 = -t425 * t438 + t437 * t491;
t487 = t456 * t455;
t474 = t448 * t487;
t490 = t452 * t451;
t424 = -t474 + t490;
t516 = -t413 * t449 - t424 * t453;
t515 = t413 * t453 - t424 * t449;
t445 = qJD(3) + qJD(4);
t450 = sin(qJ(3));
t502 = pkin(3) * qJD(3);
t441 = sin(t446);
t507 = pkin(4) * t441;
t420 = -t445 * t507 - t450 * t502;
t514 = -(t464 * t437 - t508 * t438) * t445 + t420;
t476 = t448 * t490;
t427 = -t476 + t487;
t485 = qJD(1) * t456;
t472 = t447 * t485;
t512 = -t427 * t445 + t472;
t511 = t449 * r_i_i_C(1) + t453 * r_i_i_C(2);
t442 = cos(t446);
t454 = cos(qJ(3));
t431 = t454 * pkin(3) + pkin(4) * t442;
t429 = pkin(2) + t431;
t509 = t508 * t437 + t464 * t438 + t429;
t430 = t450 * pkin(3) + t507;
t503 = pkin(8) + t430;
t426 = t448 * t489 + t488;
t409 = t426 * qJD(1) + t425 * qJD(2);
t501 = t409 * t449;
t500 = t409 * t453;
t496 = t437 * t445;
t495 = t442 * t445;
t494 = t447 * t451;
t493 = t447 * t452;
t492 = t447 * t455;
t486 = qJD(1) * t452;
t484 = qJD(2) * t451;
t483 = qJD(2) * t455;
t482 = qJD(6) * t438;
t478 = t437 * t494;
t477 = t438 * t494;
t475 = t438 * t491;
t473 = t447 * t486;
t471 = t447 * t484;
t469 = qJD(2) * t448 + qJD(1);
t410 = -qJD(1) * t476 - t452 * t484 + t469 * t487;
t470 = -t410 * t438 + t445 * t475;
t408 = t425 * qJD(1) + t426 * qJD(2);
t466 = t445 * t493 - t408;
t465 = t445 * t491 - t410;
t463 = -t425 * t445 + t473;
t462 = t445 * t448 + t447 * t483;
t398 = t465 * t437 + t463 * t438;
t396 = t466 * t437 - t512 * t438;
t397 = -t427 * t496 + t437 * t472 + t466 * t438;
t460 = t517 * (-t427 * t437 + t438 * t493) + t508 * t397 - t464 * t396;
t399 = -t425 * t496 + t437 * t473 - t470;
t459 = t517 * (-t425 * t437 - t475) + t508 * t399 + t464 * t398;
t406 = t462 * t438 - t445 * t478;
t458 = t517 * (t448 * t438 - t478) + t508 * t406 + t464 * (-t462 * t437 - t445 * t477);
t457 = t511 * t482 - t514;
t444 = -qJ(5) - pkin(10) - pkin(9);
t421 = pkin(4) * t495 + t454 * t502;
t419 = t448 * t437 + t477;
t415 = t427 * t438 + t437 * t493;
t407 = -qJD(1) * t474 - t456 * t483 + t469 * t490;
t401 = -t463 * t437 + t470;
t389 = t397 * t453 - t407 * t449 + (-t415 * t449 + t426 * t453) * qJD(6);
t388 = -t397 * t449 - t407 * t453 + (-t415 * t453 - t426 * t449) * qJD(6);
t1 = [(t401 * t453 - t501) * r_i_i_C(1) + (-t401 * t449 - t500) * r_i_i_C(2) + t401 * pkin(5) - t410 * t429 - t425 * t420 + t409 * t444 - t424 * qJD(5) + t421 * t491 + t508 * t398 + (t516 * r_i_i_C(1) - t515 * r_i_i_C(2)) * qJD(6) + (-t456 * pkin(1) - t503 * t493) * qJD(1) (-t408 * t449 + t427 * t480) * r_i_i_C(1) + (-t408 * t453 - t427 * t481) * r_i_i_C(2) + t408 * t444 + t427 * qJD(5) + t509 * t407 + t457 * t426, t408 * t430 - t427 * t421 + (t420 * t452 + t431 * t485) * t447 + t460 (-t466 * t441 + t512 * t442) * pkin(4) + t460, -t407, t388 * r_i_i_C(1) - t389 * r_i_i_C(2); t421 * t493 + t397 * pkin(5) + t389 * r_i_i_C(1) + t388 * r_i_i_C(2) + t426 * qJD(5) + t407 * t444 - t408 * t429 + t427 * t420 + t508 * t396 + (-pkin(1) * t452 + t503 * t491) * qJD(1) (t410 * t449 + t425 * t480) * r_i_i_C(1) + (t410 * t453 - t425 * t481) * r_i_i_C(2) - t410 * t444 + t425 * qJD(5) - t509 * t409 + t457 * t424, -t410 * t430 - t425 * t421 + (-t420 * t456 + t431 * t486) * t447 + t459 (t465 * t441 + t463 * t442) * pkin(4) + t459, t409 (-t399 * t449 + t500) * r_i_i_C(1) + (-t399 * t453 - t501) * r_i_i_C(2) + (t515 * r_i_i_C(1) + t516 * r_i_i_C(2)) * qJD(6); 0 ((-qJD(2) * t509 + t467 * qJD(6) + qJD(5)) * t451 + (-qJD(2) * t444 + t511 * (qJD(2) - t482) + t514) * t455) * t447, t448 * t420 + (-t421 * t451 - t430 * t483) * t447 + t458 (-t462 * t441 - t494 * t495) * pkin(4) + t458, t471 (-t406 * t449 + t453 * t471) * r_i_i_C(1) + (-t406 * t453 - t449 * t471) * r_i_i_C(2) + ((-t419 * t453 + t449 * t492) * r_i_i_C(1) + (t419 * t449 + t453 * t492) * r_i_i_C(2)) * qJD(6);];
JaD_transl  = t1;
