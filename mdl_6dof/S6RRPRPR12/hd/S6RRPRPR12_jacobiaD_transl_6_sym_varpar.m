% Zeitableitung der translatorischen Teilmatrix der analytischen Jacobi-Matrix für Segment Nr. 6 (0=Basis) von
% S6RRPRPR12
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
%   pkin=[a2,a3,a4,a5,a6,alpha2,d1,d2,d4,d6,theta5]';
%
% Output:
% JaD_transl [3x6]
%   Zeitableitung der translatorischen Teilmatrix der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-02-26 21:44
% Revision: d75aae1ac561373cd3be920984c3df29a1c2ecc8 (2019-02-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JaD_transl = S6RRPRPR12_jacobiaD_transl_6_sym_varpar(qJ, qJD, r_i_i_C, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(3,1),zeros(11,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPRPR12_jacobiaD_transl_6_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRPRPR12_jacobiaD_transl_6_sym_varpar: qJD has to be [6x1] (double)');
assert(isa(r_i_i_C,'double') && isreal(r_i_i_C) && all(size(r_i_i_C) == [3 1]), ...
	'S6RRPRPR12_jacobiaD_transl_6_sym_varpar: Position vector r_i_i_C has to be [3x1] double');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6RRPRPR12_jacobiaD_transl_6_sym_varpar: pkin has to be [11x1] (double)');

%% Symbolic Calculation
% From jacobiaD_transl_6_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 21:44:24
% EndTime: 2019-02-26 21:44:24
% DurationCPUTime: 0.58s
% Computational Cost: add. (732->117), mult. (1617->183), div. (0->0), fcn. (1604->12), ass. (0->68)
t398 = cos(pkin(6));
t406 = cos(qJ(2));
t407 = cos(qJ(1));
t435 = t406 * t407;
t429 = t398 * t435;
t402 = sin(qJ(2));
t403 = sin(qJ(1));
t438 = t402 * t403;
t381 = -t429 + t438;
t396 = qJ(4) + pkin(11);
t394 = sin(t396);
t395 = cos(t396);
t397 = sin(pkin(6));
t439 = t397 * t407;
t431 = t395 * t439;
t377 = -t381 * t394 + t431;
t436 = t403 * t406;
t437 = t402 * t407;
t382 = t398 * t437 + t436;
t400 = sin(qJ(6));
t404 = cos(qJ(6));
t451 = -t377 * t400 - t382 * t404;
t450 = t377 * t404 - t382 * t400;
t383 = t398 * t436 + t437;
t370 = t383 * qJD(1) + t382 * qJD(2);
t434 = qJD(1) * t397;
t428 = t403 * t434;
t449 = (qJD(4) * t381 + t428) * t394 - qJD(4) * t431 - t370 * t395;
t418 = t381 * t395 + t394 * t439;
t363 = t418 * qJD(4) + t370 * t394 + t395 * t428;
t421 = r_i_i_C(1) * t404 - r_i_i_C(2) * t400;
t410 = -t421 * qJD(6) - qJD(5);
t401 = sin(qJ(4));
t419 = pkin(5) + t421;
t448 = pkin(10) + r_i_i_C(3);
t409 = t401 * pkin(4) + t419 * t394 - t448 * t395 + qJ(3);
t405 = cos(qJ(4));
t447 = t405 * pkin(4);
t446 = -pkin(2) - qJ(5) - pkin(9);
t445 = pkin(8) + pkin(3) + t447;
t442 = t397 * t402;
t441 = t397 * t403;
t440 = t397 * t406;
t433 = qJD(2) * t402;
t432 = qJD(2) * t406;
t430 = t398 * t438;
t427 = t407 * t434;
t426 = t397 * t433;
t425 = t397 * t432;
t422 = qJD(2) * t398 + qJD(1);
t420 = -r_i_i_C(1) * t400 - r_i_i_C(2) * t404;
t375 = t383 * t394 + t395 * t441;
t417 = t383 * t395 - t394 * t441;
t416 = t394 * t398 + t395 * t440;
t415 = t394 * t440 - t395 * t398;
t414 = qJD(6) * t420;
t413 = -t420 - t446;
t408 = qJD(3) + t394 * t414 + (t448 * t394 + t419 * t395 + t447) * qJD(4);
t384 = -t430 + t435;
t372 = t416 * qJD(4) - t394 * t426;
t371 = -qJD(1) * t430 - t403 * t433 + t422 * t435;
t369 = t382 * qJD(1) + t383 * qJD(2);
t368 = -qJD(1) * t429 - t407 * t432 + t422 * t438;
t366 = t417 * qJD(4) - t368 * t394 + t395 * t427;
t365 = t375 * qJD(4) + t368 * t395 + t394 * t427;
t360 = t366 * t404 - t369 * t400 + (-t375 * t400 + t384 * t404) * qJD(6);
t359 = -t366 * t400 - t369 * t404 + (-t375 * t404 - t384 * t400) * qJD(6);
t1 = [-t370 * qJ(3) - t381 * qJD(3) - t382 * qJD(5) - t419 * t363 - t413 * t371 - t448 * t449 + (r_i_i_C(1) * t451 - t450 * r_i_i_C(2)) * qJD(6) + (-t407 * pkin(1) - t445 * t441) * qJD(1) + (-t370 * t401 + (-t381 * t405 - t401 * t439) * qJD(4)) * pkin(4), t413 * t368 - t409 * t369 + t410 * t383 + t408 * t384, -t368, t448 * t366 + t417 * t414 - t419 * t365 + (-t401 * t427 - t368 * t405 + (-t383 * t401 - t405 * t441) * qJD(4)) * pkin(4), -t369, r_i_i_C(1) * t359 - t360 * r_i_i_C(2); t366 * pkin(5) + t360 * r_i_i_C(1) + t359 * r_i_i_C(2) - t368 * qJ(3) + t383 * qJD(3) + t384 * qJD(5) + t446 * t369 + t448 * t365 + (-pkin(1) * t403 + t445 * t439) * qJD(1) + (-t368 * t401 + (t383 * t405 - t401 * t441) * qJD(4)) * pkin(4), -t413 * t370 + t409 * t371 + t410 * t381 + t408 * t382, t370, t448 * t363 + t418 * t414 - t419 * t449 + (-t401 * t428 + t370 * t405 + (-t381 * t401 + t405 * t439) * qJD(4)) * pkin(4), t371 (-t363 * t400 + t371 * t404) * r_i_i_C(1) + (-t363 * t404 - t371 * t400) * r_i_i_C(2) + (t450 * r_i_i_C(1) + r_i_i_C(2) * t451) * qJD(6); 0 ((t409 * qJD(2) - t410) * t406 + (-t413 * qJD(2) + t408) * t402) * t397, t426, -t448 * t372 - t416 * t414 + t419 * (t415 * qJD(4) + t395 * t426) + (t405 * t426 + (-t398 * t405 + t401 * t440) * qJD(4)) * pkin(4), t425 (t372 * t400 + t404 * t425) * r_i_i_C(1) + (t372 * t404 - t400 * t425) * r_i_i_C(2) + ((-t400 * t442 + t404 * t415) * r_i_i_C(1) + (-t400 * t415 - t404 * t442) * r_i_i_C(2)) * qJD(6);];
JaD_transl  = t1;
