% Zeitableitung der translatorischen Teilmatrix der analytischen Jacobi-Matrix für Segment Nr. 6 (0=Basis) von
% S6RRRRPR10
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
%   pkin=[a2,a3,a4,a5,a6,alpha2,d1,d2,d3,d4,d6]';
%
% Output:
% JaD_transl [3x6]
%   Zeitableitung der translatorischen Teilmatrix der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-02-26 22:35
% Revision: d75aae1ac561373cd3be920984c3df29a1c2ecc8 (2019-02-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JaD_transl = S6RRRRPR10_jacobiaD_transl_6_sym_varpar(qJ, qJD, r_i_i_C, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(3,1),zeros(11,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRRPR10_jacobiaD_transl_6_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRRRPR10_jacobiaD_transl_6_sym_varpar: qJD has to be [6x1] (double)');
assert(isa(r_i_i_C,'double') && isreal(r_i_i_C) && all(size(r_i_i_C) == [3 1]), ...
	'S6RRRRPR10_jacobiaD_transl_6_sym_varpar: Position vector r_i_i_C has to be [3x1] double');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6RRRRPR10_jacobiaD_transl_6_sym_varpar: pkin has to be [11x1] (double)');

%% Symbolic Calculation
% From jacobiaD_transl_6_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 22:35:55
% EndTime: 2019-02-26 22:35:56
% DurationCPUTime: 0.81s
% Computational Cost: add. (1172->117), mult. (2134->191), div. (0->0), fcn. (2129->12), ass. (0->76)
t441 = sin(qJ(6));
t445 = cos(qJ(6));
t463 = t445 * r_i_i_C(1) - t441 * r_i_i_C(2);
t505 = t463 * qJD(6) + qJD(5);
t462 = -t441 * r_i_i_C(1) - t445 * r_i_i_C(2);
t506 = qJ(5) - t462;
t478 = pkin(4) + pkin(11) + r_i_i_C(3);
t443 = sin(qJ(2));
t444 = sin(qJ(1));
t447 = cos(qJ(2));
t448 = cos(qJ(1));
t493 = cos(pkin(6));
t467 = t448 * t493;
t420 = t443 * t467 + t444 * t447;
t439 = qJ(3) + qJ(4);
t437 = cos(t439);
t440 = sin(pkin(6));
t484 = t440 * t448;
t429 = t437 * t484;
t436 = sin(t439);
t409 = t420 * t436 + t429;
t464 = t447 * t467;
t483 = t444 * t443;
t419 = -t464 + t483;
t503 = t409 * t441 + t419 * t445;
t502 = -t409 * t445 + t419 * t441;
t438 = qJD(3) + qJD(4);
t442 = sin(qJ(3));
t489 = t437 * t438;
t450 = -(-t478 * t438 + t505) * t436 + qJD(3) * t442 * pkin(3) - t506 * t489;
t446 = cos(qJ(3));
t435 = t446 * pkin(3) + pkin(2);
t500 = t436 * t506 + t478 * t437 + t435;
t494 = -pkin(5) - pkin(10) - pkin(9);
t490 = t436 * t438;
t488 = t440 * t443;
t487 = t440 * t444;
t486 = t440 * t446;
t485 = t440 * t447;
t482 = t448 * t447;
t481 = qJD(1) * t440;
t480 = qJD(2) * t443;
t479 = qJD(2) * t447;
t475 = t436 * t488;
t474 = t437 * t488;
t473 = t436 * t484;
t472 = t444 * t481;
t471 = t448 * t481;
t470 = t440 * t480;
t469 = t440 * t479;
t468 = t444 * t493;
t459 = qJD(2) * t493 + qJD(1);
t465 = t443 * t468;
t408 = -qJD(1) * t465 - t444 * t480 + t459 * t482;
t466 = -t408 * t437 + t438 * t429;
t421 = t448 * t443 + t447 * t468;
t406 = t420 * qJD(1) + t421 * qJD(2);
t460 = t438 * t487 - t406;
t457 = t462 * qJD(6);
t456 = t463 - t494;
t455 = t493 * t438 + t469;
t393 = t408 * t436 + t420 * t489 - t437 * t472 - t438 * t473;
t422 = -t465 + t482;
t391 = t422 * t489 + t460 * t436 - t437 * t471;
t392 = -t422 * t490 + t436 * t471 + t460 * t437;
t453 = t505 * (t422 * t437 + t436 * t487) + t506 * t392 - t478 * t391;
t452 = t505 * (t420 * t437 - t473) + t506 * (-t420 * t490 + t436 * t472 - t466) - t478 * t393;
t401 = t455 * t436 + t438 * t474;
t451 = t505 * (t493 * t436 + t474) + t506 * (t455 * t437 - t438 * t475) - t478 * t401;
t417 = -t493 * t437 + t475;
t412 = t422 * t436 - t437 * t487;
t407 = t421 * qJD(1) + t420 * qJD(2);
t405 = -qJD(1) * t464 - t448 * t479 + t459 * t483;
t382 = t391 * t441 - t405 * t445 + (t412 * t445 - t421 * t441) * qJD(6);
t381 = t391 * t445 + t405 * t441 + (-t412 * t441 - t421 * t445) * qJD(6);
t1 = [-t409 * qJD(5) - t408 * t435 - t506 * t393 - t456 * t407 + (t502 * r_i_i_C(1) + t503 * r_i_i_C(2)) * qJD(6) + (-t448 * pkin(1) - pkin(8) * t487) * qJD(1) + t478 * ((t420 * t438 - t472) * t436 + t466) + (-t442 * t472 + (t420 * t442 + t446 * t484) * qJD(3)) * pkin(3), t500 * t405 - t456 * t406 + t450 * t421 + t422 * t457 (t446 * t471 + t406 * t442 + (-t422 * t446 - t442 * t487) * qJD(3)) * pkin(3) + t453, t453, t391, t381 * r_i_i_C(1) - t382 * r_i_i_C(2); t382 * r_i_i_C(1) + t381 * r_i_i_C(2) + t391 * qJ(5) + t412 * qJD(5) - t406 * t435 + t494 * t405 + (-pkin(1) * t444 + pkin(8) * t484) * qJD(1) + t478 * t392 + (t442 * t471 + (-t422 * t442 + t444 * t486) * qJD(3)) * pkin(3), -t407 * t500 + t456 * t408 + t450 * t419 + t420 * t457 (t446 * t472 - t408 * t442 + (-t420 * t446 + t442 * t484) * qJD(3)) * pkin(3) + t452, t452, t393 (t393 * t445 - t407 * t441) * r_i_i_C(1) + (-t393 * t441 - t407 * t445) * r_i_i_C(2) + (-t503 * r_i_i_C(1) + t502 * r_i_i_C(2)) * qJD(6); 0 ((-qJD(2) * t500 + t457) * t443 + (t456 * qJD(2) - t450) * t447) * t440 (-t442 * t469 + (-t493 * t442 - t443 * t486) * qJD(3)) * pkin(3) + t451, t451, t401 (t401 * t445 - t441 * t470) * r_i_i_C(1) + (-t401 * t441 - t445 * t470) * r_i_i_C(2) + ((-t417 * t441 + t445 * t485) * r_i_i_C(1) + (-t417 * t445 - t441 * t485) * r_i_i_C(2)) * qJD(6);];
JaD_transl  = t1;
