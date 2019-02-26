% Zeitableitung der translatorischen Teilmatrix der analytischen Jacobi-Matrix für Segment Nr. 5 (0=Basis) von
% S6RRRRPR9
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
% Datum: 2019-02-26 22:35
% Revision: d75aae1ac561373cd3be920984c3df29a1c2ecc8 (2019-02-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JaD_transl = S6RRRRPR9_jacobiaD_transl_5_sym_varpar(qJ, qJD, r_i_i_C, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(3,1),zeros(12,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRRPR9_jacobiaD_transl_5_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRRRPR9_jacobiaD_transl_5_sym_varpar: qJD has to be [6x1] (double)');
assert(isa(r_i_i_C,'double') && isreal(r_i_i_C) && all(size(r_i_i_C) == [3 1]), ...
	'S6RRRRPR9_jacobiaD_transl_5_sym_varpar: Position vector r_i_i_C has to be [3x1] double');
assert(isreal(pkin) && all(size(pkin) == [12 1]), ...
  'S6RRRRPR9_jacobiaD_transl_5_sym_varpar: pkin has to be [12x1] (double)');

%% Symbolic Calculation
% From jacobiaD_transl_5_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 22:35:15
% EndTime: 2019-02-26 22:35:15
% DurationCPUTime: 0.43s
% Computational Cost: add. (797->95), mult. (1447->150), div. (0->0), fcn. (1413->12), ass. (0->62)
t431 = r_i_i_C(3) + qJ(5);
t387 = sin(pkin(12));
t389 = cos(pkin(12));
t406 = r_i_i_C(1) * t389 - r_i_i_C(2) * t387 + pkin(4);
t394 = cos(qJ(3));
t382 = t394 * pkin(3) + pkin(2);
t386 = qJ(3) + qJ(4);
t383 = sin(t386);
t384 = cos(t386);
t434 = t431 * t383 + t406 * t384 + t382;
t385 = qJD(3) + qJD(4);
t430 = t383 * t385;
t429 = t384 * t385;
t388 = sin(pkin(6));
t392 = sin(qJ(2));
t428 = t388 * t392;
t393 = sin(qJ(1));
t427 = t388 * t393;
t426 = t388 * t394;
t396 = cos(qJ(1));
t425 = t388 * t396;
t424 = t393 * t392;
t395 = cos(qJ(2));
t423 = t393 * t395;
t422 = t396 * t392;
t421 = t396 * t395;
t420 = qJD(1) * t388;
t419 = qJD(2) * t392;
t418 = qJD(2) * t395;
t417 = t384 * t428;
t390 = cos(pkin(6));
t416 = t390 * t424;
t415 = t383 * t425;
t414 = t384 * t425;
t413 = t390 * t421;
t412 = t393 * t420;
t411 = t396 * t420;
t410 = t388 * t418;
t408 = qJD(2) * t390 + qJD(1);
t363 = -qJD(1) * t416 - t393 * t419 + t408 * t421;
t409 = -t363 * t384 + t385 * t414;
t370 = t390 * t422 + t423;
t404 = t390 * t423 + t422;
t361 = t370 * qJD(1) + t404 * qJD(2);
t407 = t385 * t427 - t361;
t397 = -pkin(10) - pkin(9);
t405 = t387 * r_i_i_C(1) + t389 * r_i_i_C(2) - t397;
t403 = t385 * t390 + t410;
t351 = t363 * t383 + t370 * t429 - t384 * t412 - t385 * t415;
t372 = -t416 + t421;
t349 = t372 * t429 + t407 * t383 - t384 * t411;
t350 = -t372 * t430 + t383 * t411 + t407 * t384;
t402 = -(-t372 * t384 - t383 * t427) * qJD(5) + t431 * t350 - t406 * t349;
t401 = -(-t370 * t384 + t415) * qJD(5) + t431 * (-t370 * t430 + t383 * t412 - t409) - t406 * t351;
t358 = t403 * t383 + t385 * t417;
t400 = -(-t390 * t383 - t417) * qJD(5) + t431 * (t403 * t384 - t428 * t430) - t406 * t358;
t391 = sin(qJ(3));
t398 = -qJD(3) * t391 * pkin(3) + t383 * qJD(5) + (-t406 * t383 + t431 * t384) * t385;
t362 = t404 * qJD(1) + t370 * qJD(2);
t360 = -qJD(1) * t413 - t396 * t418 + t408 * t424;
t354 = (t370 * t385 - t412) * t383 + t409;
t1 = [(t354 * t389 - t362 * t387) * r_i_i_C(1) + (-t354 * t387 - t362 * t389) * r_i_i_C(2) + t354 * pkin(4) - (t370 * t383 + t414) * qJD(5) - t363 * t382 + t362 * t397 - t431 * t351 + (-t396 * pkin(1) - pkin(8) * t427) * qJD(1) + (-t391 * t412 + (t370 * t391 + t394 * t425) * qJD(3)) * pkin(3), t434 * t360 - t405 * t361 - t398 * t404 (t394 * t411 + t361 * t391 + (-t372 * t394 - t391 * t427) * qJD(3)) * pkin(3) + t402, t402, t349, 0; (t350 * t389 - t360 * t387) * r_i_i_C(1) + (-t350 * t387 - t360 * t389) * r_i_i_C(2) + t350 * pkin(4) - (-t372 * t383 + t384 * t427) * qJD(5) - t361 * t382 + t360 * t397 + t431 * t349 + (-t393 * pkin(1) + pkin(8) * t425) * qJD(1) + (t391 * t411 + (-t372 * t391 + t393 * t426) * qJD(3)) * pkin(3), t405 * t363 - t434 * t362 + t398 * (t413 - t424) (t394 * t412 - t363 * t391 + (-t370 * t394 + t391 * t425) * qJD(3)) * pkin(3) + t401, t401, t351, 0; 0 (-t434 * t419 + (t405 * qJD(2) + t398) * t395) * t388 (-t391 * t410 + (-t390 * t391 - t392 * t426) * qJD(3)) * pkin(3) + t400, t400, t358, 0;];
JaD_transl  = t1;
