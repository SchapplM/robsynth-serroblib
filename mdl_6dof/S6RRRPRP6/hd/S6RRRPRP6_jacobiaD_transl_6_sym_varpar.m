% Zeitableitung der translatorischen Teilmatrix der analytischen Jacobi-Matrix für Segment Nr. 6 (0=Basis) von
% S6RRRPRP6
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
%   pkin=[a2,a3,a4,a5,a6,alpha2,d1,d2,d3,d5,theta4]';
%
% Output:
% JaD_transl [3x6]
%   Zeitableitung der translatorischen Teilmatrix der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-02-26 22:12
% Revision: d75aae1ac561373cd3be920984c3df29a1c2ecc8 (2019-02-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JaD_transl = S6RRRPRP6_jacobiaD_transl_6_sym_varpar(qJ, qJD, r_i_i_C, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(3,1),zeros(11,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRPRP6_jacobiaD_transl_6_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRRPRP6_jacobiaD_transl_6_sym_varpar: qJD has to be [6x1] (double)');
assert(isa(r_i_i_C,'double') && isreal(r_i_i_C) && all(size(r_i_i_C) == [3 1]), ...
	'S6RRRPRP6_jacobiaD_transl_6_sym_varpar: Position vector r_i_i_C has to be [3x1] double');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6RRRPRP6_jacobiaD_transl_6_sym_varpar: pkin has to be [11x1] (double)');

%% Symbolic Calculation
% From jacobiaD_transl_6_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 22:12:06
% EndTime: 2019-02-26 22:12:06
% DurationCPUTime: 0.77s
% Computational Cost: add. (883->136), mult. (1858->212), div. (0->0), fcn. (1858->12), ass. (0->76)
t399 = sin(qJ(1));
t402 = cos(qJ(2));
t446 = cos(pkin(6));
t451 = cos(qJ(1));
t420 = t446 * t451;
t398 = sin(qJ(2));
t426 = t399 * t446;
t421 = t398 * t426;
t427 = t451 * qJD(1);
t438 = qJD(2) * t398;
t364 = -qJD(1) * t421 - t399 * t438 + (qJD(2) * t420 + t427) * t402;
t393 = sin(pkin(6));
t432 = t393 * t451;
t460 = -qJD(3) * t432 + t364;
t375 = t398 * t420 + t399 * t402;
t442 = t393 * t399;
t430 = qJD(1) * t442;
t459 = -qJD(3) * t375 + t430;
t392 = qJ(3) + pkin(11);
t390 = sin(t392);
t391 = cos(t392);
t368 = t375 * t391 - t390 * t432;
t416 = t402 * t420;
t439 = t398 * t399;
t374 = -t416 + t439;
t396 = sin(qJ(5));
t400 = cos(qJ(5));
t458 = t368 * t396 - t374 * t400;
t457 = t368 * t400 + t374 * t396;
t452 = pkin(5) + r_i_i_C(1);
t455 = t400 * r_i_i_C(2) + t452 * t396;
t397 = sin(qJ(3));
t388 = pkin(5) * t400 + pkin(4);
t450 = t396 * r_i_i_C(2);
t415 = r_i_i_C(1) * t400 + t388 - t450;
t447 = r_i_i_C(3) + qJ(6) + pkin(10);
t456 = -(t397 * pkin(3) + t415 * t390 - t447 * t391) * qJD(3) + t390 * qJD(6);
t358 = t459 * t390 + t460 * t391;
t401 = cos(qJ(3));
t389 = pkin(3) * t401 + pkin(2);
t454 = t447 * t390 + t415 * t391 + t389;
t443 = t393 * t398;
t441 = t393 * t401;
t440 = t393 * t402;
t436 = qJD(5) * t391;
t435 = qJD(5) * t396;
t434 = qJD(5) * t400;
t431 = t451 * t402;
t429 = qJD(2) * t440;
t428 = t393 * t438;
t423 = t393 * t427;
t376 = t451 * t398 + t402 * t426;
t363 = t376 * qJD(1) + t375 * qJD(2);
t419 = -t358 * t396 + t363 * t400;
t373 = t446 * t390 + t391 * t443;
t413 = -t373 * t400 + t396 * t440;
t377 = t431 - t421;
t370 = -t377 * t390 + t391 * t442;
t371 = t377 * t391 + t390 * t442;
t409 = t375 * t390 + t391 * t432;
t406 = -t390 * t443 + t446 * t391;
t366 = t406 * qJD(3) + t391 * t429;
t408 = -t366 * t396 + t400 * t428;
t407 = qJD(5) * t455;
t357 = t460 * t390 - t459 * t391;
t361 = -qJD(1) * t416 - qJD(2) * t431 + (qJD(2) * t446 + qJD(1)) * t439;
t405 = -t361 * t396 + (-t371 * t396 + t376 * t400) * qJD(5);
t362 = t375 * qJD(1) + t376 * qJD(2);
t356 = t370 * qJD(3) - t362 * t391 + t390 * t423;
t353 = -t356 * t396 - t361 * t400 + (-t371 * t400 - t376 * t396) * qJD(5);
t403 = t455 * t436 - t456;
t395 = -qJ(4) - pkin(9);
t365 = t373 * qJD(3) + t390 * t429;
t355 = t371 * qJD(3) - t362 * t390 - t391 * t423;
t354 = t356 * t400 + t405;
t1 = [-t409 * qJD(6) - t364 * t389 - t374 * qJD(4) - t415 * t358 + (t395 - t455) * t363 - t447 * t357 + (-t451 * pkin(1) - pkin(8) * t442) * qJD(1) + (-t397 * t430 + (t375 * t397 + t401 * t432) * qJD(3)) * pkin(3) + (t457 * r_i_i_C(2) + t452 * t458) * qJD(5) (-t362 * t400 - t377 * t435) * r_i_i_C(2) + t362 * t395 + t377 * qJD(4) + t454 * t361 + t403 * t376 + t452 * (-t362 * t396 + t377 * t434) t371 * qJD(6) + t447 * t356 - t415 * t355 - t370 * t407 + (t401 * t423 + t362 * t397 + (-t377 * t401 - t397 * t442) * qJD(3)) * pkin(3), -t361, -r_i_i_C(2) * t354 + t452 * t353, t355; t354 * r_i_i_C(1) + t353 * r_i_i_C(2) + t376 * qJD(4) - t370 * qJD(6) + t356 * t388 + t361 * t395 - t362 * t389 + t447 * t355 + (-t399 * pkin(1) + pkin(8) * t432) * qJD(1) + t405 * pkin(5) + (t397 * t423 + (-t377 * t397 + t399 * t441) * qJD(3)) * pkin(3) (t364 * t400 - t375 * t435) * r_i_i_C(2) - t364 * t395 + t375 * qJD(4) - t454 * t363 + t403 * t374 + t452 * (t364 * t396 + t375 * t434) t368 * qJD(6) + t447 * t358 - t415 * t357 + t409 * t407 + (t401 * t430 - t364 * t397 + (-t375 * t401 + t397 * t432) * qJD(3)) * pkin(3), t363, t419 * r_i_i_C(1) + (-t358 * t400 - t363 * t396) * r_i_i_C(2) + (-r_i_i_C(1) * t457 + t458 * r_i_i_C(2)) * qJD(5) + (-qJD(5) * t457 + t419) * pkin(5), t357; 0 ((qJD(4) + (t452 * t400 - t450) * qJD(5) - t454 * qJD(2)) * t398 + (-qJD(2) * t395 + t455 * (qJD(2) - t436) + t456) * t402) * t393, t373 * qJD(6) + t447 * t366 - t415 * t365 - t406 * t407 + (-t397 * t429 + (-t446 * t397 - t398 * t441) * qJD(3)) * pkin(3), t428, t408 * r_i_i_C(1) + (-t366 * t400 - t396 * t428) * r_i_i_C(2) + (t413 * r_i_i_C(1) + (t373 * t396 + t400 * t440) * r_i_i_C(2)) * qJD(5) + (t413 * qJD(5) + t408) * pkin(5), t365;];
JaD_transl  = t1;
