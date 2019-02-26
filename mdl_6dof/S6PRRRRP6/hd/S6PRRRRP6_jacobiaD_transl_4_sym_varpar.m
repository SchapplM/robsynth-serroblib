% Zeitableitung der translatorischen Teilmatrix der analytischen Jacobi-Matrix für Segment Nr. 4 (0=Basis) von
% S6PRRRRP6
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
%   pkin=[a2,a3,a4,a5,a6,alpha2,alpha3,d2,d3,d4,d5,theta1]';
%
% Output:
% JaD_transl [3x6]
%   Zeitableitung der translatorischen Teilmatrix der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-02-26 20:18
% Revision: d75aae1ac561373cd3be920984c3df29a1c2ecc8 (2019-02-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JaD_transl = S6PRRRRP6_jacobiaD_transl_4_sym_varpar(qJ, qJD, r_i_i_C, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(3,1),zeros(12,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRRRRP6_jacobiaD_transl_4_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6PRRRRP6_jacobiaD_transl_4_sym_varpar: qJD has to be [6x1] (double)');
assert(isa(r_i_i_C,'double') && isreal(r_i_i_C) && all(size(r_i_i_C) == [3 1]), ...
	'S6PRRRRP6_jacobiaD_transl_4_sym_varpar: Position vector r_i_i_C has to be [3x1] double');
assert(isreal(pkin) && all(size(pkin) == [12 1]), ...
  'S6PRRRRP6_jacobiaD_transl_4_sym_varpar: pkin has to be [12x1] (double)');

%% Symbolic Calculation
% From jacobiaD_transl_4_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 20:18:08
% EndTime: 2019-02-26 20:18:09
% DurationCPUTime: 0.59s
% Computational Cost: add. (380->110), mult. (1266->209), div. (0->0), fcn. (1334->12), ass. (0->67)
t400 = sin(qJ(3));
t403 = cos(qJ(3));
t393 = sin(pkin(12));
t396 = cos(pkin(12));
t404 = cos(qJ(2));
t398 = cos(pkin(6));
t401 = sin(qJ(2));
t429 = t398 * t401;
t412 = t393 * t429 - t396 * t404;
t397 = cos(pkin(7));
t428 = t398 * t404;
t413 = t393 * t428 + t396 * t401;
t394 = sin(pkin(7));
t395 = sin(pkin(6));
t436 = t394 * t395;
t414 = t393 * t436 - t397 * t413;
t372 = t414 * t400 - t403 * t412;
t388 = t393 * t404 + t396 * t429;
t419 = t396 * t428;
t387 = -t393 * t401 + t419;
t415 = -t387 * t397 + t396 * t436;
t440 = -t388 * t403 + t415 * t400;
t439 = r_i_i_C(3) + pkin(10);
t435 = t394 * t398;
t399 = sin(qJ(4));
t434 = t394 * t399;
t402 = cos(qJ(4));
t433 = t394 * t402;
t432 = t395 * t397;
t431 = t397 * t400;
t430 = t397 * t403;
t427 = t400 * t401;
t426 = t400 * t404;
t425 = t401 * t403;
t424 = t403 * t404;
t423 = qJD(2) * t401;
t422 = qJD(2) * t404;
t421 = qJD(4) * t399;
t420 = qJD(4) * t402;
t418 = qJD(3) * t435;
t417 = t423 * t436;
t416 = t402 * r_i_i_C(1) - t399 * r_i_i_C(2) + pkin(3);
t375 = t387 * t403 - t388 * t431;
t376 = -t403 * t413 + t412 * t431;
t411 = t397 * t424 - t427;
t410 = -t397 * t425 - t426;
t409 = t397 * t426 + t425;
t408 = t397 * t427 - t424;
t407 = qJD(4) * (-t399 * r_i_i_C(1) - t402 * r_i_i_C(2));
t406 = -t388 * t400 - t415 * t403;
t405 = t400 * t412 + t414 * t403;
t386 = t398 * t397 - t404 * t436;
t385 = t412 * qJD(2);
t384 = t413 * qJD(2);
t383 = t388 * qJD(2);
t382 = -qJD(2) * t419 + t393 * t423;
t381 = t408 * t395;
t380 = t393 * t432 + t394 * t413;
t379 = -t387 * t394 - t396 * t432;
t378 = t409 * t395 + t400 * t435;
t374 = (-t409 * qJD(2) + t410 * qJD(3)) * t395;
t368 = t403 * t418 + (-t408 * qJD(2) + t411 * qJD(3)) * t395;
t366 = t384 * t431 + t385 * t403 + (t400 * t413 + t412 * t430) * qJD(3);
t364 = t382 * t431 - t383 * t403 + (-t387 * t400 - t388 * t430) * qJD(3);
t362 = t405 * qJD(3) - t384 * t403 + t385 * t431;
t360 = t406 * qJD(3) - t382 * t403 - t383 * t431;
t1 = [0 (t366 * t402 - t376 * t421) * r_i_i_C(1) + (-t366 * t399 - t376 * t420) * r_i_i_C(2) + t366 * pkin(3) + t385 * pkin(2) + t439 * (t376 * qJD(3) - t384 * t430 + t385 * t400) + ((-t384 * t399 - t412 * t420) * r_i_i_C(1) + (-t384 * t402 + t412 * t421) * r_i_i_C(2) - t384 * pkin(9)) * t394, t439 * t362 + t405 * t407 + t416 * (-t372 * qJD(3) + t384 * t400 + t385 * t430) (-t362 * t399 - t385 * t433) * r_i_i_C(1) + (-t362 * t402 + t385 * t434) * r_i_i_C(2) + ((-t372 * t402 - t380 * t399) * r_i_i_C(1) + (t372 * t399 - t380 * t402) * r_i_i_C(2)) * qJD(4), 0, 0; 0 (t364 * t402 - t375 * t421) * r_i_i_C(1) + (-t364 * t399 - t375 * t420) * r_i_i_C(2) + t364 * pkin(3) - t383 * pkin(2) + t439 * (t375 * qJD(3) - t382 * t430 - t383 * t400) + ((-t382 * t399 + t388 * t420) * r_i_i_C(1) + (-t382 * t402 - t388 * t421) * r_i_i_C(2) - t382 * pkin(9)) * t394, t439 * t360 + t406 * t407 + t416 * (t440 * qJD(3) + t382 * t400 - t383 * t430) (-t360 * t399 + t383 * t433) * r_i_i_C(1) + (-t360 * t402 - t383 * t434) * r_i_i_C(2) + ((-t379 * t399 + t402 * t440) * r_i_i_C(1) + (-t379 * t402 - t399 * t440) * r_i_i_C(2)) * qJD(4), 0, 0; 0 (t374 * t402 + t381 * t421) * r_i_i_C(1) + (-t374 * t399 + t381 * t420) * r_i_i_C(2) + t374 * pkin(3) + (-t439 * (-t411 * qJD(2) + t408 * qJD(3)) - pkin(2) * t423 + ((t399 * t422 + t401 * t420) * r_i_i_C(1) + (-t401 * t421 + t402 * t422) * r_i_i_C(2) + pkin(9) * t422) * t394) * t395, t439 * t368 + (t411 * t395 + t403 * t435) * t407 + t416 * (-t400 * t418 + (t410 * qJD(2) - t409 * qJD(3)) * t395) (-t368 * t399 + t402 * t417) * r_i_i_C(1) + (-t368 * t402 - t399 * t417) * r_i_i_C(2) + ((-t378 * t402 - t386 * t399) * r_i_i_C(1) + (t378 * t399 - t386 * t402) * r_i_i_C(2)) * qJD(4), 0, 0;];
JaD_transl  = t1;
