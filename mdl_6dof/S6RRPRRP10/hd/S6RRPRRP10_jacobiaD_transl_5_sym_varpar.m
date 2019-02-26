% Zeitableitung der translatorischen Teilmatrix der analytischen Jacobi-Matrix für Segment Nr. 5 (0=Basis) von
% S6RRPRRP10
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
%   pkin=[a2,a3,a4,a5,a6,alpha2,d1,d2,d4,d5,theta3]';
%
% Output:
% JaD_transl [3x6]
%   Zeitableitung der translatorischen Teilmatrix der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-02-26 21:51
% Revision: d75aae1ac561373cd3be920984c3df29a1c2ecc8 (2019-02-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JaD_transl = S6RRPRRP10_jacobiaD_transl_5_sym_varpar(qJ, qJD, r_i_i_C, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(3,1),zeros(11,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPRRP10_jacobiaD_transl_5_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRPRRP10_jacobiaD_transl_5_sym_varpar: qJD has to be [6x1] (double)');
assert(isa(r_i_i_C,'double') && isreal(r_i_i_C) && all(size(r_i_i_C) == [3 1]), ...
	'S6RRPRRP10_jacobiaD_transl_5_sym_varpar: Position vector r_i_i_C has to be [3x1] double');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6RRPRRP10_jacobiaD_transl_5_sym_varpar: pkin has to be [11x1] (double)');

%% Symbolic Calculation
% From jacobiaD_transl_5_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 21:51:06
% EndTime: 2019-02-26 21:51:06
% DurationCPUTime: 0.58s
% Computational Cost: add. (656->106), mult. (1384->172), div. (0->0), fcn. (1389->12), ass. (0->69)
t392 = sin(qJ(1));
t388 = cos(pkin(6));
t407 = qJD(2) * t388 + qJD(1);
t391 = sin(qJ(2));
t426 = t391 * t392;
t415 = t388 * t426;
t421 = qJD(2) * t391;
t394 = cos(qJ(2));
t395 = cos(qJ(1));
t423 = t394 * t395;
t361 = -qJD(1) * t415 - t392 * t421 + t407 * t423;
t424 = t392 * t394;
t425 = t391 * t395;
t372 = t388 * t425 + t424;
t385 = pkin(11) + qJ(4);
t383 = sin(t385);
t384 = cos(t385);
t387 = sin(pkin(6));
t422 = qJD(1) * t387;
t412 = t392 * t422;
t427 = t387 * t395;
t416 = t384 * t427;
t355 = (-qJD(4) * t372 + t412) * t383 - qJD(4) * t416 + t361 * t384;
t373 = t388 * t424 + t425;
t360 = t373 * qJD(1) + t372 * qJD(2);
t390 = sin(qJ(5));
t393 = cos(qJ(5));
t444 = t355 * t390 - t360 * t393;
t443 = -t355 * t393 - t360 * t390;
t405 = r_i_i_C(1) * t393 - r_i_i_C(2) * t390;
t403 = pkin(4) + t405;
t437 = r_i_i_C(3) + pkin(10);
t442 = (t403 * t383 - t437 * t384) * qJD(4);
t366 = -t372 * t384 + t383 * t427;
t414 = t388 * t423;
t371 = -t414 + t426;
t441 = -t366 * t390 - t371 * t393;
t440 = t366 * t393 - t371 * t390;
t404 = r_i_i_C(1) * t390 + r_i_i_C(2) * t393;
t382 = cos(pkin(11)) * pkin(3) + pkin(2);
t438 = t437 * t383 + t403 * t384 + t382;
t430 = t387 * t391;
t429 = t387 * t392;
t428 = t387 * t394;
t420 = qJD(2) * t394;
t419 = qJD(5) * t384;
t418 = qJD(5) * t390;
t417 = qJD(5) * t393;
t413 = pkin(3) * sin(pkin(11)) + pkin(8);
t411 = t395 * t422;
t410 = t387 * t420;
t409 = t387 * t421;
t374 = -t415 + t423;
t402 = -t374 * t383 + t384 * t429;
t368 = t374 * t384 + t383 * t429;
t370 = t383 * t388 + t384 * t430;
t401 = -t383 * t430 + t384 * t388;
t400 = qJD(5) * t404;
t397 = t366 * qJD(4) - t361 * t383 + t384 * t412;
t396 = t404 * t419 + t442;
t389 = -pkin(9) - qJ(3);
t363 = t401 * qJD(4) + t384 * t410;
t359 = t372 * qJD(1) + t373 * qJD(2);
t358 = -qJD(1) * t414 - t395 * t420 + t407 * t426;
t353 = t402 * qJD(4) - t359 * t384 + t383 * t411;
t352 = t368 * qJD(4) - t359 * t383 - t384 * t411;
t351 = t353 * t393 - t358 * t390 + (-t368 * t390 + t373 * t393) * qJD(5);
t350 = -t353 * t390 - t358 * t393 + (-t368 * t393 - t373 * t390) * qJD(5);
t1 = [t443 * r_i_i_C(1) + t444 * r_i_i_C(2) - t355 * pkin(4) - t361 * t382 + t360 * t389 - t371 * qJD(3) + t437 * t397 + (t441 * r_i_i_C(1) - t440 * r_i_i_C(2)) * qJD(5) + (-t395 * pkin(1) - t413 * t429) * qJD(1) (-t359 * t390 + t374 * t417) * r_i_i_C(1) + (-t359 * t393 - t374 * t418) * r_i_i_C(2) + t359 * t389 + t374 * qJD(3) + t438 * t358 + t396 * t373, -t358, -t403 * t352 + t437 * t353 - t402 * t400, r_i_i_C(1) * t350 - r_i_i_C(2) * t351, 0; t353 * pkin(4) + t351 * r_i_i_C(1) + t350 * r_i_i_C(2) + t373 * qJD(3) + t358 * t389 - t359 * t382 + t437 * t352 + (-pkin(1) * t392 + t413 * t427) * qJD(1) (t361 * t390 + t372 * t417) * r_i_i_C(1) + (t361 * t393 - t372 * t418) * r_i_i_C(2) - t361 * t389 + t372 * qJD(3) - t438 * t360 + t396 * t371, t360, t437 * t355 - (-t372 * t383 - t416) * t400 + t403 * t397, -t444 * r_i_i_C(1) + t443 * r_i_i_C(2) + (t440 * r_i_i_C(1) + t441 * r_i_i_C(2)) * qJD(5), 0; 0 ((-qJD(2) * t438 + t405 * qJD(5) + qJD(3)) * t391 + (-qJD(2) * t389 - t442 + t404 * (qJD(2) - t419)) * t394) * t387, t409, t437 * t363 - t401 * t400 + t403 * (-t370 * qJD(4) - t383 * t410) (-t363 * t390 + t393 * t409) * r_i_i_C(1) + (-t363 * t393 - t390 * t409) * r_i_i_C(2) + ((-t370 * t393 + t390 * t428) * r_i_i_C(1) + (t370 * t390 + t393 * t428) * r_i_i_C(2)) * qJD(5), 0;];
JaD_transl  = t1;
