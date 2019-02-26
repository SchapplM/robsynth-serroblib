% Zeitableitung der translatorischen Teilmatrix der analytischen Jacobi-Matrix für Segment Nr. 6 (0=Basis) von
% S6RRPRRP13
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
%   pkin=[a2,a3,a4,a5,a6,alpha2,d1,d2,d4,d5]';
%
% Output:
% JaD_transl [3x6]
%   Zeitableitung der translatorischen Teilmatrix der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-02-26 21:53
% Revision: d75aae1ac561373cd3be920984c3df29a1c2ecc8 (2019-02-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JaD_transl = S6RRPRRP13_jacobiaD_transl_6_sym_varpar(qJ, qJD, r_i_i_C, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(3,1),zeros(10,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPRRP13_jacobiaD_transl_6_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRPRRP13_jacobiaD_transl_6_sym_varpar: qJD has to be [6x1] (double)');
assert(isa(r_i_i_C,'double') && isreal(r_i_i_C) && all(size(r_i_i_C) == [3 1]), ...
	'S6RRPRRP13_jacobiaD_transl_6_sym_varpar: Position vector r_i_i_C has to be [3x1] double');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RRPRRP13_jacobiaD_transl_6_sym_varpar: pkin has to be [10x1] (double)');

%% Symbolic Calculation
% From jacobiaD_transl_6_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 21:52:51
% EndTime: 2019-02-26 21:52:51
% DurationCPUTime: 0.51s
% Computational Cost: add. (621->106), mult. (1799->163), div. (0->0), fcn. (1806->10), ass. (0->68)
t401 = cos(qJ(5));
t397 = sin(qJ(5));
t445 = t397 * r_i_i_C(2);
t446 = pkin(5) + r_i_i_C(1);
t409 = (t446 * t401 - t445) * qJD(5);
t418 = -t401 * r_i_i_C(2) - t446 * t397;
t448 = -pkin(2) - pkin(9);
t412 = t418 + t448;
t395 = cos(pkin(6));
t398 = sin(qJ(4));
t402 = cos(qJ(4));
t394 = sin(pkin(6));
t403 = cos(qJ(2));
t438 = t394 * t403;
t381 = t395 * t402 - t398 * t438;
t399 = sin(qJ(2));
t404 = cos(qJ(1));
t433 = t404 * t399;
t400 = sin(qJ(1));
t434 = t400 * t403;
t383 = t395 * t433 + t434;
t432 = t404 * t403;
t435 = t400 * t399;
t382 = -t395 * t432 + t435;
t437 = t394 * t404;
t427 = t402 * t437;
t416 = -t382 * t398 + t427;
t450 = -t383 * t401 - t397 * t416;
t421 = -t383 * t397 + t401 * t416;
t384 = t395 * t434 + t433;
t370 = t384 * qJD(1) + t383 * qJD(2);
t415 = t382 * t402 + t398 * t437;
t431 = qJD(1) * t394;
t426 = t400 * t431;
t362 = t415 * qJD(4) + t370 * t398 + t402 * t426;
t393 = t401 * pkin(5) + pkin(4);
t419 = t401 * r_i_i_C(1) + t393 - t445;
t443 = r_i_i_C(3) + qJ(6) + pkin(10);
t406 = t419 * t398 - t443 * t402 + qJ(3);
t447 = pkin(3) + pkin(8);
t440 = t394 * t399;
t439 = t394 * t400;
t430 = qJD(2) * t399;
t428 = t395 * t435;
t425 = t404 * t431;
t424 = qJD(2) * t438;
t423 = t394 * t430;
t371 = -qJD(1) * t428 - t400 * t430 + (qJD(2) * t395 + qJD(1)) * t432;
t422 = -t362 * t397 + t371 * t401;
t417 = -t381 * t401 - t397 * t440;
t375 = t384 * t398 + t402 * t439;
t374 = t384 * t402 - t398 * t439;
t414 = t395 * t398 + t402 * t438;
t413 = t428 - t432;
t372 = t414 * qJD(4) - t398 * t423;
t410 = t372 * t397 + t401 * t424;
t408 = qJD(5) * t418;
t369 = t383 * qJD(1) + t384 * qJD(2);
t407 = -t369 * t397 + (-t375 * t397 - t401 * t413) * qJD(5);
t360 = -t370 * t402 - qJD(4) * t427 + (qJD(4) * t382 + t426) * t398;
t368 = t382 * qJD(1) + t413 * qJD(2);
t365 = t374 * qJD(4) - t368 * t398 + t402 * t425;
t358 = -t365 * t397 - t369 * t401 + (-t375 * t401 + t397 * t413) * qJD(5);
t405 = -t402 * qJD(6) + qJD(3) + t398 * t408 + (t443 * t398 + t419 * t402) * qJD(4);
t373 = t381 * qJD(4) - t402 * t423;
t364 = t375 * qJD(4) + t368 * t402 + t398 * t425;
t359 = t365 * t401 + t407;
t1 = [t415 * qJD(6) - t370 * qJ(3) - t382 * qJD(3) - t443 * t360 - t419 * t362 + t412 * t371 + (-t404 * pkin(1) - t447 * t439) * qJD(1) + (-t421 * r_i_i_C(2) + t446 * t450) * qJD(5), -t412 * t368 - t406 * t369 - t384 * t409 - t405 * t413, -t368, t375 * qJD(6) - t419 * t364 + t443 * t365 + t374 * t408, -t359 * r_i_i_C(2) + t446 * t358, t364; t359 * r_i_i_C(1) + t358 * r_i_i_C(2) - t368 * qJ(3) + t384 * qJD(3) - t374 * qJD(6) + t365 * t393 + t448 * t369 + t443 * t364 + (-t400 * pkin(1) + t447 * t437) * qJD(1) + t407 * pkin(5), t412 * t370 + t406 * t371 - t382 * t409 + t405 * t383, t370, -qJD(6) * t416 - t419 * t360 + t443 * t362 + t415 * t408, t422 * r_i_i_C(1) + (-t362 * t401 - t371 * t397) * r_i_i_C(2) + (t421 * r_i_i_C(1) + t450 * r_i_i_C(2)) * qJD(5) + (t421 * qJD(5) + t422) * pkin(5), t360; 0 ((t406 * qJD(2) + t409) * t403 + (t412 * qJD(2) + t405) * t399) * t394, t423, t381 * qJD(6) - t443 * t372 - t419 * t373 - t414 * t408, t410 * r_i_i_C(1) + (t372 * t401 - t397 * t424) * r_i_i_C(2) + (t417 * r_i_i_C(1) + (t381 * t397 - t401 * t440) * r_i_i_C(2)) * qJD(5) + (t417 * qJD(5) + t410) * pkin(5), t373;];
JaD_transl  = t1;
