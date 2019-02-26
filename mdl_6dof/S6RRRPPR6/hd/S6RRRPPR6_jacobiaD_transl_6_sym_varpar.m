% Zeitableitung der translatorischen Teilmatrix der analytischen Jacobi-Matrix für Segment Nr. 6 (0=Basis) von
% S6RRRPPR6
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
%   pkin=[a2,a3,a4,a5,a6,alpha2,d1,d2,d3,d6,theta4]';
%
% Output:
% JaD_transl [3x6]
%   Zeitableitung der translatorischen Teilmatrix der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-02-26 22:06
% Revision: d75aae1ac561373cd3be920984c3df29a1c2ecc8 (2019-02-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JaD_transl = S6RRRPPR6_jacobiaD_transl_6_sym_varpar(qJ, qJD, r_i_i_C, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(3,1),zeros(11,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRPPR6_jacobiaD_transl_6_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRRPPR6_jacobiaD_transl_6_sym_varpar: qJD has to be [6x1] (double)');
assert(isa(r_i_i_C,'double') && isreal(r_i_i_C) && all(size(r_i_i_C) == [3 1]), ...
	'S6RRRPPR6_jacobiaD_transl_6_sym_varpar: Position vector r_i_i_C has to be [3x1] double');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6RRRPPR6_jacobiaD_transl_6_sym_varpar: pkin has to be [11x1] (double)');

%% Symbolic Calculation
% From jacobiaD_transl_6_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 22:06:25
% EndTime: 2019-02-26 22:06:25
% DurationCPUTime: 0.57s
% Computational Cost: add. (836->115), mult. (1772->186), div. (0->0), fcn. (1768->12), ass. (0->68)
t389 = sin(qJ(2));
t390 = sin(qJ(1));
t393 = cos(qJ(2));
t394 = cos(qJ(1));
t433 = cos(pkin(6));
t410 = t394 * t433;
t369 = t389 * t410 + t390 * t393;
t384 = qJ(3) + pkin(11);
t382 = sin(t384);
t383 = cos(t384);
t385 = sin(pkin(6));
t423 = t385 * t394;
t443 = t369 * t383 - t382 * t423;
t408 = t393 * t410;
t422 = t390 * t389;
t368 = -t408 + t422;
t387 = sin(qJ(6));
t391 = cos(qJ(6));
t439 = t369 * t382 + t383 * t423;
t442 = t368 * t391 + t387 * t439;
t441 = t368 * t387 - t391 * t439;
t388 = sin(qJ(3));
t407 = t391 * r_i_i_C(1) - t387 * r_i_i_C(2);
t397 = t407 * qJD(6) + qJD(5);
t406 = -t387 * r_i_i_C(1) - t391 * r_i_i_C(2);
t404 = qJ(5) - t406;
t417 = r_i_i_C(3) + pkin(10) + pkin(4);
t395 = (t388 * pkin(3) + t417 * t382 - t404 * t383) * qJD(3) - t397 * t382;
t405 = qJD(2) * t433 + qJD(1);
t411 = t390 * t433;
t409 = t389 * t411;
t419 = qJD(2) * t389;
t421 = t394 * t393;
t357 = -qJD(1) * t409 - t390 * t419 + t405 * t421;
t420 = qJD(1) * t385;
t415 = t390 * t420;
t438 = qJD(3) * t439 - t357 * t383 - t382 * t415;
t392 = cos(qJ(3));
t381 = t392 * pkin(3) + pkin(2);
t437 = t404 * t382 + t417 * t383 + t381;
t434 = -pkin(5) - qJ(4) - pkin(9);
t371 = -t409 + t421;
t428 = t371 * t382;
t427 = t385 * t389;
t426 = t385 * t390;
t425 = t385 * t392;
t424 = t385 * t393;
t418 = qJD(2) * t393;
t414 = t394 * t420;
t413 = t385 * t418;
t412 = t385 * t419;
t402 = t371 * t383 + t382 * t426;
t401 = t407 - t434;
t370 = t394 * t389 + t393 * t411;
t366 = t382 * t427 - t433 * t383;
t399 = t433 * t382 + t383 * t427;
t398 = t406 * qJD(6) + qJD(4);
t350 = t443 * qJD(3) + t357 * t382 - t383 * t415;
t363 = -t383 * t426 + t428;
t358 = t399 * qJD(3) + t382 * t413;
t356 = t370 * qJD(1) + t369 * qJD(2);
t355 = t369 * qJD(1) + t370 * qJD(2);
t354 = -qJD(1) * t408 - t394 * t418 + t405 * t422;
t349 = t382 * t414 - qJD(3) * t428 + (qJD(3) * t426 - t355) * t383;
t348 = t402 * qJD(3) - t355 * t382 - t383 * t414;
t347 = t348 * t387 - t354 * t391 + (t363 * t391 - t370 * t387) * qJD(6);
t346 = t348 * t391 + t354 * t387 + (-t363 * t387 - t370 * t391) * qJD(6);
t1 = [-t368 * qJD(4) - t439 * qJD(5) - t357 * t381 - t404 * t350 - t401 * t356 + (t441 * r_i_i_C(1) + t442 * r_i_i_C(2)) * qJD(6) + (-t394 * pkin(1) - pkin(8) * t426) * qJD(1) + t417 * t438 + (-t388 * t415 + (t369 * t388 + t392 * t423) * qJD(3)) * pkin(3), t437 * t354 - t401 * t355 + t395 * t370 + t398 * t371, t397 * t402 + t404 * t349 - t417 * t348 + (t392 * t414 + t355 * t388 + (-t371 * t392 - t388 * t426) * qJD(3)) * pkin(3), -t354, t348, t346 * r_i_i_C(1) - t347 * r_i_i_C(2); t347 * r_i_i_C(1) + t346 * r_i_i_C(2) + t348 * qJ(5) + t370 * qJD(4) + t363 * qJD(5) - t355 * t381 + t434 * t354 + (-pkin(1) * t390 + pkin(8) * t423) * qJD(1) + t417 * t349 + (t388 * t414 + (-t371 * t388 + t390 * t425) * qJD(3)) * pkin(3), -t356 * t437 + t401 * t357 + t395 * t368 + t398 * t369, t397 * t443 - t404 * t438 - t417 * t350 + (t392 * t415 - t357 * t388 + (-t369 * t392 + t388 * t423) * qJD(3)) * pkin(3), t356, t350 (t350 * t391 - t356 * t387) * r_i_i_C(1) + (-t350 * t387 - t356 * t391) * r_i_i_C(2) + (-t442 * r_i_i_C(1) + t441 * r_i_i_C(2)) * qJD(6); 0 ((-qJD(2) * t437 + t398) * t389 + (t401 * qJD(2) - t395) * t393) * t385, t397 * t399 + t404 * (-t366 * qJD(3) + t383 * t413) - t417 * t358 + (-t388 * t413 + (-t433 * t388 - t389 * t425) * qJD(3)) * pkin(3), t412, t358 (t358 * t391 - t387 * t412) * r_i_i_C(1) + (-t358 * t387 - t391 * t412) * r_i_i_C(2) + ((-t366 * t387 + t391 * t424) * r_i_i_C(1) + (-t366 * t391 - t387 * t424) * r_i_i_C(2)) * qJD(6);];
JaD_transl  = t1;
