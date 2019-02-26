% Zeitableitung der translatorischen Teilmatrix der analytischen Jacobi-Matrix für Segment Nr. 6 (0=Basis) von
% S6RRPRPR5
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
%   pkin=[a2,a3,a4,a5,a6,alpha2,d1,d2,d4,d6,theta3,theta5]';
%
% Output:
% JaD_transl [3x6]
%   Zeitableitung der translatorischen Teilmatrix der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-02-26 21:40
% Revision: d75aae1ac561373cd3be920984c3df29a1c2ecc8 (2019-02-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JaD_transl = S6RRPRPR5_jacobiaD_transl_6_sym_varpar(qJ, qJD, r_i_i_C, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(3,1),zeros(12,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPRPR5_jacobiaD_transl_6_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRPRPR5_jacobiaD_transl_6_sym_varpar: qJD has to be [6x1] (double)');
assert(isa(r_i_i_C,'double') && isreal(r_i_i_C) && all(size(r_i_i_C) == [3 1]), ...
	'S6RRPRPR5_jacobiaD_transl_6_sym_varpar: Position vector r_i_i_C has to be [3x1] double');
assert(isreal(pkin) && all(size(pkin) == [12 1]), ...
  'S6RRPRPR5_jacobiaD_transl_6_sym_varpar: pkin has to be [12x1] (double)');

%% Symbolic Calculation
% From jacobiaD_transl_6_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 21:40:10
% EndTime: 2019-02-26 21:40:11
% DurationCPUTime: 0.74s
% Computational Cost: add. (1046->119), mult. (2775->199), div. (0->0), fcn. (3058->14), ass. (0->83)
t463 = sin(pkin(6));
t470 = cos(qJ(1));
t510 = t463 * t470;
t464 = cos(pkin(6));
t462 = sin(pkin(11));
t514 = cos(pkin(11));
t517 = cos(qJ(2));
t489 = t517 * t514;
t467 = sin(qJ(2));
t504 = qJD(2) * t467;
t521 = -qJD(2) * t489 + t462 * t504;
t434 = t521 * t464;
t492 = t467 * t514;
t477 = t517 * t462 + t492;
t439 = t477 * t464;
t493 = t517 * qJD(2);
t441 = -qJD(2) * t492 - t462 * t493;
t468 = sin(qJ(1));
t476 = -t467 * t462 + t489;
t506 = qJD(1) * t468;
t520 = t439 * t506 - t468 * t441 + (-qJD(1) * t476 + t434) * t470;
t528 = qJD(4) * t510 + t520;
t422 = t470 * t439 + t468 * t476;
t497 = t463 * t506;
t527 = -qJD(4) * t422 + t497;
t466 = sin(qJ(4));
t469 = cos(qJ(4));
t416 = t422 * t469 - t466 * t510;
t475 = t476 * t464;
t421 = -t468 * t477 + t470 * t475;
t460 = pkin(12) + qJ(6);
t458 = sin(t460);
t459 = cos(t460);
t526 = t416 * t458 + t421 * t459;
t525 = t416 * t459 - t421 * t458;
t487 = t458 * r_i_i_C(1) + t459 * r_i_i_C(2);
t481 = qJD(6) * t487;
t456 = cos(pkin(12)) * pkin(5) + pkin(4);
t488 = t459 * r_i_i_C(1) - t458 * r_i_i_C(2);
t485 = t456 + t488;
t515 = r_i_i_C(3) + pkin(10) + qJ(5);
t524 = (t485 * t466 - t515 * t469) * qJD(4) - t466 * qJD(5) + t469 * t481;
t522 = qJD(1) * t475 + t476 * qJD(2);
t404 = t527 * t466 - t528 * t469;
t518 = t515 * t466 + t485 * t469 + pkin(3);
t516 = pkin(2) * t464;
t511 = t463 * t468;
t508 = t467 * t468;
t507 = t467 * t470;
t505 = qJD(1) * t470;
t501 = pkin(2) * t504;
t500 = sin(pkin(12)) * pkin(5) + pkin(9);
t499 = t517 * t468;
t498 = t517 * t470;
t496 = t463 * t505;
t438 = t477 * t463;
t427 = t438 * t469 + t464 * t466;
t486 = -t438 * t466 + t464 * t469;
t425 = -t468 * t439 + t470 * t476;
t483 = t422 * t466 + t469 * t510;
t418 = -t425 * t466 + t469 * t511;
t419 = t425 * t469 + t466 * t511;
t482 = qJD(6) * t488;
t478 = t487 + t500;
t403 = -t528 * t466 - t527 * t469;
t409 = -t422 * qJD(1) + t468 * t434 + t470 * t441;
t457 = t517 * pkin(2) + pkin(1);
t442 = -t463 * qJD(3) + t493 * t516;
t440 = t467 * t516 + (-pkin(8) - qJ(3)) * t463;
t437 = t476 * t463;
t435 = qJD(2) * t439;
t433 = qJD(2) * t438;
t432 = t521 * t463;
t424 = -t468 * t475 - t470 * t477;
t414 = t486 * qJD(4) - t432 * t469;
t413 = t427 * qJD(4) - t432 * t466;
t411 = -t470 * t435 - t468 * t522 - t477 * t505;
t408 = -t468 * t435 + t470 * t522 - t477 * t506;
t402 = t418 * qJD(4) + t409 * t469 + t466 * t496;
t401 = t419 * qJD(4) + t409 * t466 - t469 * t496;
t400 = t402 * t459 + t408 * t458 + (-t419 * t458 - t424 * t459) * qJD(6);
t399 = -t402 * t458 + t408 * t459 + (-t419 * t459 + t424 * t458) * qJD(6);
t1 = [-t483 * qJD(5) + t520 * pkin(3) + t468 * t501 - t470 * t442 - t485 * t404 + t478 * t411 - t515 * t403 + (t526 * r_i_i_C(1) + t525 * r_i_i_C(2)) * qJD(6) + (t468 * t440 - t470 * t457) * qJD(1), t425 * t482 + t478 * t409 - t518 * t408 + ((t464 * t508 - t498) * qJD(2) + (-t464 * t498 + t508) * qJD(1)) * pkin(2) - t524 * t424, t496, t419 * qJD(5) - t485 * t401 + t515 * t402 - t418 * t481, t401, t399 * r_i_i_C(1) - t400 * r_i_i_C(2); -t470 * t501 + t409 * pkin(3) + t400 * r_i_i_C(1) + t399 * r_i_i_C(2) - t418 * qJD(5) + t402 * t456 - t468 * t442 + t500 * t408 + t515 * t401 + (-t440 * t470 - t457 * t468) * qJD(1), t422 * t482 - t478 * t520 + t518 * t411 + ((-t464 * t507 - t499) * qJD(2) + (-t464 * t499 - t507) * qJD(1)) * pkin(2) - t524 * t421, t497, t416 * qJD(5) - t485 * t403 + t515 * t404 + t483 * t481, t403 (-t404 * t458 - t411 * t459) * r_i_i_C(1) + (-t404 * t459 + t411 * t458) * r_i_i_C(2) + (-t525 * r_i_i_C(1) + t526 * r_i_i_C(2)) * qJD(6); 0, -t478 * t432 - t433 * t518 - t524 * t437 + t438 * t482 - t463 * t501, 0, t427 * qJD(5) - t485 * t413 + t515 * t414 - t486 * t481, t413 (-t414 * t458 + t433 * t459) * r_i_i_C(1) + (-t414 * t459 - t433 * t458) * r_i_i_C(2) + ((-t427 * t459 + t437 * t458) * r_i_i_C(1) + (t427 * t458 + t437 * t459) * r_i_i_C(2)) * qJD(6);];
JaD_transl  = t1;
