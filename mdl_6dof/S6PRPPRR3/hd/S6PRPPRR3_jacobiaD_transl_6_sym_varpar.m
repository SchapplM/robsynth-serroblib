% Zeitableitung der translatorischen Teilmatrix der analytischen Jacobi-Matrix für Segment Nr. 6 (0=Basis) von
% S6PRPPRR3
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
%   pkin=[a2,a3,a4,a5,a6,alpha2,d2,d5,d6,theta1,theta4]';
%
% Output:
% JaD_transl [3x6]
%   Zeitableitung der translatorischen Teilmatrix der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-02-26 19:46
% Revision: d75aae1ac561373cd3be920984c3df29a1c2ecc8 (2019-02-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JaD_transl = S6PRPPRR3_jacobiaD_transl_6_sym_varpar(qJ, qJD, r_i_i_C, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(3,1),zeros(11,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRPPRR3_jacobiaD_transl_6_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6PRPPRR3_jacobiaD_transl_6_sym_varpar: qJD has to be [6x1] (double)');
assert(isa(r_i_i_C,'double') && isreal(r_i_i_C) && all(size(r_i_i_C) == [3 1]), ...
	'S6PRPPRR3_jacobiaD_transl_6_sym_varpar: Position vector r_i_i_C has to be [3x1] double');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6PRPPRR3_jacobiaD_transl_6_sym_varpar: pkin has to be [11x1] (double)');

%% Symbolic Calculation
% From jacobiaD_transl_6_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 19:45:51
% EndTime: 2019-02-26 19:45:52
% DurationCPUTime: 0.46s
% Computational Cost: add. (446->92), mult. (1402->170), div. (0->0), fcn. (1504->12), ass. (0->59)
t393 = sin(pkin(11));
t395 = sin(pkin(6));
t396 = cos(pkin(11));
t401 = sin(qJ(2));
t404 = cos(qJ(2));
t378 = (-t393 * t404 + t396 * t401) * t395;
t427 = -pkin(3) - pkin(2);
t426 = pkin(9) + r_i_i_C(3);
t400 = sin(qJ(5));
t424 = t395 * t400;
t403 = cos(qJ(5));
t423 = t395 * t403;
t398 = cos(pkin(6));
t421 = t398 * t401;
t420 = t398 * t404;
t419 = qJD(2) * t401;
t418 = qJD(2) * t404;
t399 = sin(qJ(6));
t417 = qJD(6) * t399;
t402 = cos(qJ(6));
t416 = qJD(6) * t402;
t394 = sin(pkin(10));
t415 = t394 * t419;
t414 = t395 * t419;
t397 = cos(pkin(10));
t413 = t397 * t418;
t379 = -t398 * t413 + t415;
t384 = t394 * t404 + t397 * t421;
t380 = t384 * qJD(2);
t351 = -t379 * t393 - t380 * t396;
t352 = -t379 * t396 + t380 * t393;
t385 = t394 * t420 + t397 * t401;
t381 = t385 * qJD(2);
t382 = -t398 * t415 + t413;
t355 = -t381 * t393 - t382 * t396;
t356 = -t381 * t396 + t382 * t393;
t383 = t394 * t401 - t397 * t420;
t361 = -t383 * t396 + t384 * t393;
t386 = -t394 * t421 + t397 * t404;
t365 = -t385 * t396 + t386 * t393;
t368 = t378 * t403 - t398 * t400;
t412 = -t378 * t400 - t398 * t403;
t362 = t383 * t393 + t384 * t396;
t366 = t385 * t393 + t386 * t396;
t411 = t402 * r_i_i_C(1) - t399 * r_i_i_C(2) + pkin(5);
t410 = -t362 * t400 + t397 * t423;
t348 = t362 * t403 + t397 * t424;
t409 = -t366 * t400 - t394 * t423;
t408 = -t366 * t403 + t394 * t424;
t407 = qJD(6) * (-t399 * r_i_i_C(1) - t402 * r_i_i_C(2));
t406 = t426 * t400 + t411 * t403 + pkin(4);
t405 = t403 * t407 + (-t411 * t400 + t426 * t403) * qJD(5);
t377 = (t393 * t401 + t396 * t404) * t395;
t374 = qJD(2) * t378;
t373 = -t395 * t396 * t418 - t393 * t414;
t346 = t412 * qJD(5) - t373 * t403;
t344 = t409 * qJD(5) + t356 * t403;
t342 = t410 * qJD(5) + t352 * t403;
t1 = [0 (-t356 * t399 - t366 * t416) * r_i_i_C(1) + (-t356 * t402 + t366 * t417) * r_i_i_C(2) - t356 * pkin(8) - t381 * qJ(3) + t386 * qJD(3) + t427 * t382 + t406 * t355 + t405 * t365, t382, 0, t426 * t344 + t409 * t407 + t411 * (t408 * qJD(5) - t356 * t400) (-t344 * t399 + t355 * t402) * r_i_i_C(1) + (-t344 * t402 - t355 * t399) * r_i_i_C(2) + ((-t365 * t399 + t402 * t408) * r_i_i_C(1) + (-t365 * t402 - t399 * t408) * r_i_i_C(2)) * qJD(6); 0 (-t352 * t399 - t362 * t416) * r_i_i_C(1) + (-t352 * t402 + t362 * t417) * r_i_i_C(2) - t352 * pkin(8) - t379 * qJ(3) + t384 * qJD(3) + t427 * t380 + t406 * t351 + t405 * t361, t380, 0, t426 * t342 + t410 * t407 + t411 * (-t348 * qJD(5) - t352 * t400) (-t342 * t399 + t351 * t402) * r_i_i_C(1) + (-t342 * t402 - t351 * t399) * r_i_i_C(2) + ((-t348 * t402 - t361 * t399) * r_i_i_C(1) + (t348 * t399 - t361 * t402) * r_i_i_C(2)) * qJD(6); 0 (t373 * t399 - t378 * t416) * r_i_i_C(1) + (t373 * t402 + t378 * t417) * r_i_i_C(2) + t373 * pkin(8) + (t401 * qJD(3) + (qJ(3) * t404 + t427 * t401) * qJD(2)) * t395 - t406 * t374 + t405 * t377, t414, 0, t426 * t346 + t412 * t407 + t411 * (-t368 * qJD(5) + t373 * t400) (-t346 * t399 - t374 * t402) * r_i_i_C(1) + (-t346 * t402 + t374 * t399) * r_i_i_C(2) + ((-t368 * t402 - t377 * t399) * r_i_i_C(1) + (t368 * t399 - t377 * t402) * r_i_i_C(2)) * qJD(6);];
JaD_transl  = t1;
