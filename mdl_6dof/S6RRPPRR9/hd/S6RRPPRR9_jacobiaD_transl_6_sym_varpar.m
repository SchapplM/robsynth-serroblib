% Zeitableitung der translatorischen Teilmatrix der analytischen Jacobi-Matrix für Segment Nr. 6 (0=Basis) von
% S6RRPPRR9
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
%   pkin=[a2,a3,a4,a5,a6,alpha2,d1,d2,d5,d6]';
%
% Output:
% JaD_transl [3x6]
%   Zeitableitung der translatorischen Teilmatrix der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-02-26 21:33
% Revision: d75aae1ac561373cd3be920984c3df29a1c2ecc8 (2019-02-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JaD_transl = S6RRPPRR9_jacobiaD_transl_6_sym_varpar(qJ, qJD, r_i_i_C, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(3,1),zeros(10,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPPRR9_jacobiaD_transl_6_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRPPRR9_jacobiaD_transl_6_sym_varpar: qJD has to be [6x1] (double)');
assert(isa(r_i_i_C,'double') && isreal(r_i_i_C) && all(size(r_i_i_C) == [3 1]), ...
	'S6RRPPRR9_jacobiaD_transl_6_sym_varpar: Position vector r_i_i_C has to be [3x1] double');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RRPPRR9_jacobiaD_transl_6_sym_varpar: pkin has to be [10x1] (double)');

%% Symbolic Calculation
% From jacobiaD_transl_6_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 21:33:15
% EndTime: 2019-02-26 21:33:15
% DurationCPUTime: 0.49s
% Computational Cost: add. (512->99), mult. (1517->159), div. (0->0), fcn. (1514->10), ass. (0->65)
t388 = sin(qJ(2));
t389 = sin(qJ(1));
t392 = cos(qJ(2));
t393 = cos(qJ(1));
t431 = cos(pkin(6));
t412 = t393 * t431;
t373 = t388 * t412 + t389 * t392;
t387 = sin(qJ(5));
t391 = cos(qJ(5));
t385 = sin(pkin(6));
t425 = t385 * t393;
t368 = -t373 * t387 + t391 * t425;
t408 = t392 * t412;
t424 = t389 * t388;
t372 = -t408 + t424;
t386 = sin(qJ(6));
t390 = cos(qJ(6));
t440 = -t368 * t386 + t372 * t390;
t439 = t368 * t390 + t372 * t386;
t406 = t386 * r_i_i_C(1) + t390 * r_i_i_C(2);
t401 = qJD(6) * t406;
t407 = -t390 * r_i_i_C(1) + t386 * r_i_i_C(2);
t404 = pkin(5) - t407;
t434 = pkin(10) + r_i_i_C(3);
t394 = -(t434 * t387 + t404 * t391) * qJD(5) - qJD(4) + t387 * t401;
t405 = qJD(2) * t431 + qJD(1);
t413 = t389 * t431;
t409 = t388 * t413;
t421 = qJD(2) * t388;
t423 = t393 * t392;
t362 = -qJD(1) * t409 - t389 * t421 + t405 * t423;
t427 = t385 * t391;
t416 = qJD(5) * t427;
t422 = qJD(1) * t385;
t418 = t389 * t422;
t437 = (qJD(5) * t373 + t418) * t387 - t362 * t391 - t393 * t416;
t403 = t373 * t391 + t387 * t425;
t354 = t403 * qJD(5) + t362 * t387 + t391 * t418;
t433 = -pkin(2) - qJ(4);
t435 = t404 * t387 - t434 * t391 - t433;
t432 = pkin(9) - qJ(3);
t428 = t385 * t389;
t426 = t385 * t392;
t420 = qJD(2) * t392;
t419 = pkin(3) + pkin(4) + pkin(8);
t417 = t393 * t422;
t415 = t385 * t421;
t414 = t385 * t420;
t411 = t431 * t387;
t375 = -t409 + t423;
t366 = t375 * t387 + t389 * t427;
t402 = t375 * t391 - t387 * t428;
t400 = t406 + t432;
t374 = t393 * t388 + t392 * t413;
t397 = -t385 * t388 * t387 - t431 * t391;
t396 = t407 * qJD(6) + qJD(3);
t363 = qJD(5) * t411 - t387 * t414 - t388 * t416;
t361 = t374 * qJD(1) + t373 * qJD(2);
t360 = t373 * qJD(1) + t374 * qJD(2);
t359 = -qJD(1) * t408 - t393 * t420 + t405 * t424;
t357 = t402 * qJD(5) - t360 * t387 + t391 * t417;
t356 = t366 * qJD(5) + t360 * t391 + t387 * t417;
t351 = t357 * t390 + t359 * t386 + (-t366 * t386 - t374 * t390) * qJD(6);
t350 = -t357 * t386 + t359 * t390 + (-t366 * t390 + t374 * t386) * qJD(6);
t1 = [-t372 * qJD(3) - t373 * qJD(4) + t433 * t362 - t404 * t354 + t400 * t361 - t434 * t437 + (r_i_i_C(1) * t440 - t439 * r_i_i_C(2)) * qJD(6) + (-t393 * pkin(1) - t419 * t428) * qJD(1), t435 * t359 + t400 * t360 + t394 * t374 + t396 * t375, -t359, -t360, -t404 * t356 + t434 * t357 - t402 * t401, t350 * r_i_i_C(1) - t351 * r_i_i_C(2); t357 * pkin(5) + t351 * r_i_i_C(1) + t350 * r_i_i_C(2) + t374 * qJD(3) + t375 * qJD(4) + t433 * t360 + t432 * t359 + t434 * t356 + (-pkin(1) * t389 + t419 * t425) * qJD(1), -t361 * t435 - t400 * t362 + t394 * t372 + t396 * t373, t361, t362, t434 * t354 - t403 * t401 - t404 * t437 (-t354 * t386 - t361 * t390) * r_i_i_C(1) + (-t354 * t390 + t361 * t386) * r_i_i_C(2) + (t439 * r_i_i_C(1) + r_i_i_C(2) * t440) * qJD(6); 0 ((-qJD(2) * t435 + t396) * t388 + (-t400 * qJD(2) - t394) * t392) * t385, t415, t414, -t434 * t363 - (t388 * t427 - t411) * t401 + t404 * (t397 * qJD(5) + t391 * t414) (t363 * t386 - t390 * t415) * r_i_i_C(1) + (t363 * t390 + t386 * t415) * r_i_i_C(2) + ((-t386 * t426 + t390 * t397) * r_i_i_C(1) + (-t386 * t397 - t390 * t426) * r_i_i_C(2)) * qJD(6);];
JaD_transl  = t1;
