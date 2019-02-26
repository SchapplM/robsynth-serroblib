% Zeitableitung der translatorischen Teilmatrix der analytischen Jacobi-Matrix für Segment Nr. 5 (0=Basis) von
% S6PRRRPR5
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
%   pkin=[a2,a3,a4,a5,a6,alpha2,alpha3,d2,d3,d4,d6,theta1,theta5]';
%
% Output:
% JaD_transl [3x6]
%   Zeitableitung der translatorischen Teilmatrix der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-02-26 20:13
% Revision: d75aae1ac561373cd3be920984c3df29a1c2ecc8 (2019-02-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JaD_transl = S6PRRRPR5_jacobiaD_transl_5_sym_varpar(qJ, qJD, r_i_i_C, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(3,1),zeros(13,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRRRPR5_jacobiaD_transl_5_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6PRRRPR5_jacobiaD_transl_5_sym_varpar: qJD has to be [6x1] (double)');
assert(isa(r_i_i_C,'double') && isreal(r_i_i_C) && all(size(r_i_i_C) == [3 1]), ...
	'S6PRRRPR5_jacobiaD_transl_5_sym_varpar: Position vector r_i_i_C has to be [3x1] double');
assert(isreal(pkin) && all(size(pkin) == [13 1]), ...
  'S6PRRRPR5_jacobiaD_transl_5_sym_varpar: pkin has to be [13x1] (double)');

%% Symbolic Calculation
% From jacobiaD_transl_5_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 20:12:57
% EndTime: 2019-02-26 20:12:58
% DurationCPUTime: 0.65s
% Computational Cost: add. (577->142), mult. (1681->254), div. (0->0), fcn. (1784->14), ass. (0->84)
t419 = sin(qJ(4));
t467 = t419 * pkin(4);
t466 = r_i_i_C(3) + qJ(5) + pkin(10);
t412 = sin(pkin(12));
t415 = cos(pkin(12));
t424 = cos(qJ(2));
t417 = cos(pkin(6));
t421 = sin(qJ(2));
t455 = t417 * t421;
t402 = t412 * t424 + t415 * t455;
t423 = cos(qJ(3));
t465 = t402 * t423;
t411 = qJ(4) + pkin(13);
t409 = sin(t411);
t413 = sin(pkin(7));
t464 = t409 * t413;
t410 = cos(t411);
t463 = t410 * t413;
t414 = sin(pkin(6));
t462 = t413 * t414;
t420 = sin(qJ(3));
t461 = t413 * t420;
t422 = cos(qJ(4));
t460 = t413 * t422;
t459 = t414 * t415;
t416 = cos(pkin(7));
t458 = t414 * t416;
t457 = t416 * t420;
t456 = t416 * t423;
t454 = t417 * t424;
t453 = t420 * t421;
t452 = t420 * t424;
t451 = t421 * t423;
t450 = t423 * t424;
t449 = qJD(2) * t421;
t448 = qJD(2) * t424;
t447 = qJD(4) * t409;
t446 = qJD(4) * t410;
t445 = qJD(4) * t421;
t444 = qJD(4) * t422;
t443 = qJD(4) * t467;
t442 = t415 * t454;
t441 = t417 * t413 * t423;
t440 = qJD(3) * t461;
t439 = t449 * t462;
t408 = t422 * pkin(4) + pkin(3);
t438 = -t410 * r_i_i_C(1) + t409 * r_i_i_C(2) - t408;
t401 = -t412 * t421 + t442;
t437 = -t401 * t420 - t402 * t456;
t387 = t401 * t423 - t402 * t457;
t436 = t401 * t416 - t413 * t459;
t432 = t412 * t455 - t415 * t424;
t433 = t412 * t454 + t415 * t421;
t435 = t420 * t433 + t432 * t456;
t388 = -t423 * t433 + t432 * t457;
t434 = t412 * t462 - t416 * t433;
t431 = t416 * t450 - t453;
t430 = t416 * t451 + t452;
t429 = t416 * t452 + t451;
t428 = t416 * t453 - t450;
t427 = qJD(4) * (-t409 * r_i_i_C(1) - t410 * r_i_i_C(2) - t467);
t426 = -t402 * t420 + t423 * t436;
t425 = t420 * t432 + t423 * t434;
t384 = t420 * t434 - t423 * t432;
t400 = t417 * t416 - t424 * t462;
t399 = t432 * qJD(2);
t398 = t433 * qJD(2);
t397 = t402 * qJD(2);
t396 = -qJD(2) * t442 + t412 * t449;
t395 = t428 * t414;
t392 = t412 * t458 + t413 * t433;
t391 = -t401 * t413 - t415 * t458;
t390 = t414 * t429 + t417 * t461;
t386 = (-qJD(2) * t429 - qJD(3) * t430) * t414;
t382 = t420 * t436 + t465;
t380 = qJD(3) * t441 + (-t428 * qJD(2) + qJD(3) * t431) * t414;
t379 = t417 * t440 + (qJD(2) * t430 + qJD(3) * t429) * t414;
t378 = qJD(3) * t435 + t398 * t457 + t399 * t423;
t376 = qJD(3) * t437 + t396 * t457 - t397 * t423;
t374 = qJD(3) * t425 - t398 * t423 + t399 * t457;
t373 = qJD(3) * t384 - t398 * t420 - t399 * t456;
t372 = qJD(3) * t426 - t396 * t423 - t397 * t457;
t371 = -t396 * t420 + t397 * t456 - t440 * t459 + (t401 * t457 + t465) * qJD(3);
t1 = [0 (t378 * t410 - t388 * t447) * r_i_i_C(1) + (-t378 * t409 - t388 * t446) * r_i_i_C(2) + t378 * t408 - t388 * t443 - t435 * qJD(5) + t399 * pkin(2) + t466 * (qJD(3) * t388 - t398 * t456 + t399 * t420) + ((-t398 * t409 - t432 * t446) * r_i_i_C(1) + (-t398 * t410 + t432 * t447) * r_i_i_C(2) - t398 * pkin(9) + (-t398 * t419 - t432 * t444) * pkin(4)) * t413, t384 * qJD(5) + t373 * t438 + t374 * t466 + t425 * t427 (-t374 * t409 - t399 * t463) * r_i_i_C(1) + (-t374 * t410 + t399 * t464) * r_i_i_C(2) + ((-t384 * t410 - t392 * t409) * r_i_i_C(1) + (t384 * t409 - t392 * t410) * r_i_i_C(2)) * qJD(4) + (-t399 * t460 - t374 * t419 + (-t384 * t422 - t392 * t419) * qJD(4)) * pkin(4), t373, 0; 0 (t376 * t410 - t387 * t447) * r_i_i_C(1) + (-t376 * t409 - t387 * t446) * r_i_i_C(2) + t376 * t408 - t387 * t443 - t437 * qJD(5) - t397 * pkin(2) + t466 * (qJD(3) * t387 - t396 * t456 - t397 * t420) + ((-t396 * t409 + t402 * t446) * r_i_i_C(1) + (-t396 * t410 - t402 * t447) * r_i_i_C(2) - t396 * pkin(9) + (-t396 * t419 + t402 * t444) * pkin(4)) * t413, t382 * qJD(5) + t371 * t438 + t372 * t466 + t426 * t427 (-t372 * t409 + t397 * t463) * r_i_i_C(1) + (-t372 * t410 - t397 * t464) * r_i_i_C(2) + ((-t382 * t410 - t391 * t409) * r_i_i_C(1) + (t382 * t409 - t391 * t410) * r_i_i_C(2)) * qJD(4) + (t397 * t460 - t372 * t419 + (-t382 * t422 - t391 * t419) * qJD(4)) * pkin(4), t371, 0; 0 (t386 * t410 + t395 * t447) * r_i_i_C(1) + (-t386 * t409 + t395 * t446) * r_i_i_C(2) + t386 * t408 + t395 * t443 + (-t466 * (-qJD(2) * t431 + t428 * qJD(3)) + t430 * qJD(5) - pkin(2) * t449 + ((t409 * t448 + t410 * t445) * r_i_i_C(1) + (-t409 * t445 + t410 * t448) * r_i_i_C(2) + pkin(9) * t448 + (t419 * t448 + t421 * t444) * pkin(4)) * t413) * t414, t390 * qJD(5) + t466 * t380 + t438 * t379 + (t414 * t431 + t441) * t427 (-t380 * t409 + t410 * t439) * r_i_i_C(1) + (-t380 * t410 - t409 * t439) * r_i_i_C(2) + ((-t390 * t410 - t400 * t409) * r_i_i_C(1) + (t390 * t409 - t400 * t410) * r_i_i_C(2)) * qJD(4) + (t422 * t439 - t380 * t419 + (-t390 * t422 - t400 * t419) * qJD(4)) * pkin(4), t379, 0;];
JaD_transl  = t1;
