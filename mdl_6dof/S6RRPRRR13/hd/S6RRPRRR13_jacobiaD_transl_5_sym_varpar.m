% Zeitableitung der translatorischen Teilmatrix der analytischen Jacobi-Matrix für Segment Nr. 5 (0=Basis) von
% S6RRPRRR13
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
%   pkin=[a2,a3,a4,a5,a6,alpha2,d1,d2,d4,d5,d6]';
%
% Output:
% JaD_transl [3x6]
%   Zeitableitung der translatorischen Teilmatrix der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-02-26 22:01
% Revision: d75aae1ac561373cd3be920984c3df29a1c2ecc8 (2019-02-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JaD_transl = S6RRPRRR13_jacobiaD_transl_5_sym_varpar(qJ, qJD, r_i_i_C, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(3,1),zeros(11,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPRRR13_jacobiaD_transl_5_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRPRRR13_jacobiaD_transl_5_sym_varpar: qJD has to be [6x1] (double)');
assert(isa(r_i_i_C,'double') && isreal(r_i_i_C) && all(size(r_i_i_C) == [3 1]), ...
	'S6RRPRRR13_jacobiaD_transl_5_sym_varpar: Position vector r_i_i_C has to be [3x1] double');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6RRPRRR13_jacobiaD_transl_5_sym_varpar: pkin has to be [11x1] (double)');

%% Symbolic Calculation
% From jacobiaD_transl_5_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 22:01:01
% EndTime: 2019-02-26 22:01:02
% DurationCPUTime: 0.45s
% Computational Cost: add. (478->92), mult. (1424->150), div. (0->0), fcn. (1424->10), ass. (0->61)
t380 = cos(pkin(6));
t387 = cos(qJ(2));
t388 = cos(qJ(1));
t413 = t388 * t387;
t383 = sin(qJ(2));
t384 = sin(qJ(1));
t416 = t384 * t383;
t370 = -t380 * t413 + t416;
t382 = sin(qJ(4));
t386 = cos(qJ(4));
t379 = sin(pkin(6));
t417 = t379 * t388;
t409 = t386 * t417;
t366 = -t370 * t382 + t409;
t414 = t388 * t383;
t415 = t384 * t387;
t371 = t380 * t414 + t415;
t381 = sin(qJ(5));
t385 = cos(qJ(5));
t428 = -t366 * t381 - t371 * t385;
t427 = t366 * t385 - t371 * t381;
t403 = t385 * r_i_i_C(1) - t381 * r_i_i_C(2);
t394 = t403 * qJD(5);
t372 = t380 * t415 + t414;
t359 = t372 * qJD(1) + t371 * qJD(2);
t412 = qJD(1) * t379;
t408 = t384 * t412;
t426 = (qJD(4) * t370 + t408) * t382 - qJD(4) * t409 - t359 * t386;
t399 = t370 * t386 + t382 * t417;
t352 = t399 * qJD(4) + t359 * t382 + t386 * t408;
t401 = pkin(4) + t403;
t423 = -r_i_i_C(3) - pkin(10);
t390 = t401 * t382 + t423 * t386 + qJ(3);
t425 = -pkin(2) - pkin(9);
t424 = pkin(3) + pkin(8);
t420 = t379 * t383;
t419 = t379 * t384;
t418 = t379 * t387;
t411 = qJD(2) * t383;
t410 = t380 * t416;
t407 = t388 * t412;
t406 = qJD(2) * t418;
t405 = t379 * t411;
t402 = -t381 * r_i_i_C(1) - t385 * r_i_i_C(2);
t400 = -t402 - t425;
t364 = t372 * t382 + t386 * t419;
t398 = t372 * t386 - t382 * t419;
t397 = t380 * t382 + t386 * t418;
t396 = -t380 * t386 + t382 * t418;
t395 = t410 - t413;
t393 = qJD(5) * t402;
t389 = qJD(3) + t382 * t393 + (-t423 * t382 + t401 * t386) * qJD(4);
t361 = t397 * qJD(4) - t382 * t405;
t360 = -qJD(1) * t410 - t384 * t411 + (qJD(2) * t380 + qJD(1)) * t413;
t358 = t371 * qJD(1) + t372 * qJD(2);
t357 = t370 * qJD(1) + t395 * qJD(2);
t355 = t398 * qJD(4) - t357 * t382 + t386 * t407;
t354 = t364 * qJD(4) + t357 * t386 + t382 * t407;
t349 = t355 * t385 - t358 * t381 + (-t364 * t381 - t385 * t395) * qJD(5);
t348 = -t355 * t381 - t358 * t385 + (-t364 * t385 + t381 * t395) * qJD(5);
t1 = [-t359 * qJ(3) - t370 * qJD(3) - t401 * t352 - t400 * t360 + t423 * t426 + (r_i_i_C(1) * t428 - t427 * r_i_i_C(2)) * qJD(5) + (-t388 * pkin(1) - t424 * t419) * qJD(1), t400 * t357 - t390 * t358 - t372 * t394 - t389 * t395, -t357, -t401 * t354 - t423 * t355 + t398 * t393, t348 * r_i_i_C(1) - t349 * r_i_i_C(2), 0; t355 * pkin(4) + t349 * r_i_i_C(1) + t348 * r_i_i_C(2) - t357 * qJ(3) + t372 * qJD(3) + t425 * t358 - t423 * t354 + (-pkin(1) * t384 + t424 * t417) * qJD(1), -t400 * t359 + t390 * t360 - t370 * t394 + t389 * t371, t359, -t423 * t352 + t399 * t393 - t401 * t426 (-t352 * t381 + t360 * t385) * r_i_i_C(1) + (-t352 * t385 - t360 * t381) * r_i_i_C(2) + (t427 * r_i_i_C(1) + r_i_i_C(2) * t428) * qJD(5), 0; 0 ((t390 * qJD(2) + t394) * t387 + (-t400 * qJD(2) + t389) * t383) * t379, t405, t423 * t361 - t397 * t393 + t401 * (t396 * qJD(4) + t386 * t405) (t361 * t381 + t385 * t406) * r_i_i_C(1) + (t361 * t385 - t381 * t406) * r_i_i_C(2) + ((-t381 * t420 + t385 * t396) * r_i_i_C(1) + (-t381 * t396 - t385 * t420) * r_i_i_C(2)) * qJD(5), 0;];
JaD_transl  = t1;
