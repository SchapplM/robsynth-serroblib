% Zeitableitung der translatorischen Teilmatrix der analytischen Jacobi-Matrix für Segment Nr. 5 (0=Basis) von
% S6PRRRPR8
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
%   pkin=[a2,a3,a4,a5,a6,alpha2,alpha3,d2,d3,d4,d6,theta1]';
%
% Output:
% JaD_transl [3x6]
%   Zeitableitung der translatorischen Teilmatrix der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-02-26 20:14
% Revision: d75aae1ac561373cd3be920984c3df29a1c2ecc8 (2019-02-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JaD_transl = S6PRRRPR8_jacobiaD_transl_5_sym_varpar(qJ, qJD, r_i_i_C, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(3,1),zeros(12,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRRRPR8_jacobiaD_transl_5_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6PRRRPR8_jacobiaD_transl_5_sym_varpar: qJD has to be [6x1] (double)');
assert(isa(r_i_i_C,'double') && isreal(r_i_i_C) && all(size(r_i_i_C) == [3 1]), ...
	'S6PRRRPR8_jacobiaD_transl_5_sym_varpar: Position vector r_i_i_C has to be [3x1] double');
assert(isreal(pkin) && all(size(pkin) == [12 1]), ...
  'S6PRRRPR8_jacobiaD_transl_5_sym_varpar: pkin has to be [12x1] (double)');

%% Symbolic Calculation
% From jacobiaD_transl_5_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 20:14:44
% EndTime: 2019-02-26 20:14:45
% DurationCPUTime: 0.58s
% Computational Cost: add. (725->118), mult. (2362->212), div. (0->0), fcn. (2557->12), ass. (0->79)
t436 = sin(pkin(12));
t439 = cos(pkin(12));
t444 = sin(qJ(2));
t441 = cos(pkin(6));
t447 = cos(qJ(2));
t477 = t441 * t447;
t430 = -t436 * t444 + t439 * t477;
t443 = sin(qJ(3));
t446 = cos(qJ(3));
t478 = t441 * t444;
t457 = t436 * t478 - t439 * t447;
t440 = cos(pkin(7));
t458 = t436 * t477 + t439 * t444;
t437 = sin(pkin(7));
t438 = sin(pkin(6));
t486 = t437 * t438;
t459 = t436 * t486 - t440 * t458;
t415 = t459 * t443 - t446 * t457;
t431 = t436 * t447 + t439 * t478;
t460 = -t430 * t440 + t439 * t486;
t494 = -t431 * t446 + t460 * t443;
t493 = -r_i_i_C(1) - pkin(10);
t492 = r_i_i_C(2) - pkin(4);
t491 = pkin(9) * t437;
t490 = r_i_i_C(3) + qJ(5);
t485 = t437 * t441;
t442 = sin(qJ(4));
t484 = t437 * t442;
t445 = cos(qJ(4));
t483 = t437 * t445;
t482 = t437 * t447;
t481 = t438 * t440;
t480 = t440 * t443;
t479 = t440 * t446;
t476 = t443 * t444;
t475 = t443 * t447;
t474 = t444 * t446;
t473 = t446 * t447;
t472 = qJD(2) * t438;
t471 = t444 * t486;
t469 = t437 * t472;
t468 = qJD(3) * t485;
t467 = t444 * t469;
t466 = t447 * t469;
t422 = -t430 * t437 - t439 * t481;
t465 = t422 * t442 - t445 * t494;
t423 = t436 * t481 + t437 * t458;
t464 = t415 * t445 + t423 * t442;
t454 = t440 * t475 + t474;
t421 = t454 * t438 + t443 * t485;
t429 = -t438 * t482 + t441 * t440;
t463 = t421 * t445 + t429 * t442;
t418 = t430 * t446 - t431 * t480;
t462 = -t418 * t442 + t431 * t483;
t419 = -t446 * t458 + t457 * t480;
t461 = -t419 * t442 - t457 * t483;
t456 = t440 * t473 - t476;
t455 = -t440 * t474 - t475;
t453 = t440 * t476 - t473;
t424 = t453 * t438;
t452 = t424 * t442 + t445 * t471;
t451 = t490 * t442 - t492 * t445 + pkin(3);
t450 = -t431 * t443 - t460 * t446;
t449 = t443 * t457 + t459 * t446;
t448 = qJD(5) * t442 + (t492 * t442 + t490 * t445) * qJD(4);
t428 = t457 * qJD(2);
t427 = t458 * qJD(2);
t426 = t431 * qJD(2);
t425 = t430 * qJD(2);
t417 = (-t454 * qJD(2) + t455 * qJD(3)) * t438;
t411 = t446 * t468 + (-t453 * qJD(2) + t456 * qJD(3)) * t438;
t409 = t427 * t480 + t428 * t446 + (t443 * t458 + t457 * t479) * qJD(3);
t407 = -t425 * t480 - t426 * t446 + (-t430 * t443 - t431 * t479) * qJD(3);
t403 = t449 * qJD(3) - t427 * t446 + t428 * t480;
t401 = t450 * qJD(3) + t425 * t446 - t426 * t480;
t398 = t463 * qJD(4) + t411 * t442 - t445 * t467;
t392 = t464 * qJD(4) + t403 * t442 + t428 * t483;
t390 = t465 * qJD(4) + t401 * t442 - t426 * t483;
t1 = [0, -t461 * qJD(5) + t409 * pkin(3) + t428 * pkin(2) - t427 * t491 - t493 * (t419 * qJD(3) - t427 * t479 + t428 * t443) - t492 * (t461 * qJD(4) + t409 * t445 - t427 * t484) + t490 * (t427 * t483 + t409 * t442 + (t419 * t445 - t457 * t484) * qJD(4)) -t493 * t403 + t448 * t449 + t451 * (-qJD(3) * t415 + t427 * t443 + t428 * t479) t464 * qJD(5) + t490 * (-t428 * t484 + t403 * t445 + (-t415 * t442 + t423 * t445) * qJD(4)) + t492 * t392, t392, 0; 0, -t462 * qJD(5) + t407 * pkin(3) - t426 * pkin(2) + t425 * t491 - t493 * (t418 * qJD(3) + t425 * t479 - t426 * t443) - t492 * (t462 * qJD(4) + t407 * t445 + t425 * t484) + t490 * (-t425 * t483 + t407 * t442 + (t418 * t445 + t431 * t484) * qJD(4)) -t493 * t401 + t448 * t450 + t451 * (qJD(3) * t494 - t425 * t443 - t426 * t479) t465 * qJD(5) + t490 * (t426 * t484 + t401 * t445 + (t422 * t445 + t442 * t494) * qJD(4)) + t492 * t390, t390, 0; 0, -t452 * qJD(5) + t417 * pkin(3) + t493 * (-t456 * qJD(2) + t453 * qJD(3)) * t438 - t492 * (t452 * qJD(4) + t417 * t445 + t442 * t466) + t490 * (-t445 * t466 + t417 * t442 + (-t424 * t445 + t442 * t471) * qJD(4)) + (-pkin(2) * t444 + pkin(9) * t482) * t472, -t493 * t411 + t448 * (t456 * t438 + t446 * t485) + t451 * (-t443 * t468 + (t455 * qJD(2) - t454 * qJD(3)) * t438) t463 * qJD(5) + t490 * (t442 * t467 + t411 * t445 + (-t421 * t442 + t429 * t445) * qJD(4)) + t492 * t398, t398, 0;];
JaD_transl  = t1;
