% Zeitableitung der translatorischen Teilmatrix der analytischen Jacobi-Matrix für Segment Nr. 6 (0=Basis) von
% S6RRPRPR10
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
%   pkin=[a2,a3,a4,a5,a6,alpha2,d1,d2,d4,d6,theta3]';
%
% Output:
% JaD_transl [3x6]
%   Zeitableitung der translatorischen Teilmatrix der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-02-26 21:43
% Revision: d75aae1ac561373cd3be920984c3df29a1c2ecc8 (2019-02-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JaD_transl = S6RRPRPR10_jacobiaD_transl_6_sym_varpar(qJ, qJD, r_i_i_C, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(3,1),zeros(11,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPRPR10_jacobiaD_transl_6_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRPRPR10_jacobiaD_transl_6_sym_varpar: qJD has to be [6x1] (double)');
assert(isa(r_i_i_C,'double') && isreal(r_i_i_C) && all(size(r_i_i_C) == [3 1]), ...
	'S6RRPRPR10_jacobiaD_transl_6_sym_varpar: Position vector r_i_i_C has to be [3x1] double');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6RRPRPR10_jacobiaD_transl_6_sym_varpar: pkin has to be [11x1] (double)');

%% Symbolic Calculation
% From jacobiaD_transl_6_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 21:43:02
% EndTime: 2019-02-26 21:43:03
% DurationCPUTime: 0.48s
% Computational Cost: add. (806->98), mult. (1673->158), div. (0->0), fcn. (1681->12), ass. (0->66)
t388 = sin(qJ(2));
t389 = sin(qJ(1));
t391 = cos(qJ(2));
t392 = cos(qJ(1));
t431 = cos(pkin(6));
t408 = t392 * t431;
t368 = t388 * t408 + t389 * t391;
t383 = pkin(11) + qJ(4);
t381 = sin(t383);
t382 = cos(t383);
t385 = sin(pkin(6));
t422 = t385 * t392;
t440 = t368 * t382 - t381 * t422;
t406 = t391 * t408;
t421 = t389 * t388;
t367 = -t406 + t421;
t387 = sin(qJ(6));
t390 = cos(qJ(6));
t436 = t368 * t381 + t382 * t422;
t439 = t367 * t390 + t387 * t436;
t438 = t367 * t387 - t390 * t436;
t405 = t390 * r_i_i_C(1) - t387 * r_i_i_C(2);
t395 = t405 * qJD(6) + qJD(5);
t404 = -t387 * r_i_i_C(1) - t390 * r_i_i_C(2);
t402 = qJ(5) - t404;
t416 = r_i_i_C(3) + pkin(10) + pkin(4);
t393 = (t416 * t381 - t402 * t382) * qJD(4) - t395 * t381;
t403 = qJD(2) * t431 + qJD(1);
t409 = t389 * t431;
t407 = t388 * t409;
t418 = qJD(2) * t388;
t420 = t392 * t391;
t356 = -qJD(1) * t407 - t389 * t418 + t403 * t420;
t419 = qJD(1) * t385;
t413 = t389 * t419;
t435 = qJD(4) * t436 - t356 * t382 - t381 * t413;
t380 = cos(pkin(11)) * pkin(3) + pkin(2);
t433 = t402 * t381 + t416 * t382 + t380;
t432 = -pkin(5) - pkin(9) - qJ(3);
t370 = -t407 + t420;
t426 = t370 * t381;
t425 = t385 * t388;
t424 = t385 * t389;
t423 = t385 * t391;
t417 = qJD(2) * t391;
t414 = pkin(3) * sin(pkin(11)) + pkin(8);
t412 = t392 * t419;
t411 = t385 * t417;
t410 = t385 * t418;
t400 = t370 * t382 + t381 * t424;
t399 = t405 - t432;
t369 = t392 * t388 + t391 * t409;
t365 = t381 * t425 - t431 * t382;
t397 = t431 * t381 + t382 * t425;
t396 = t404 * qJD(6) + qJD(3);
t349 = t440 * qJD(4) + t356 * t381 - t382 * t413;
t362 = -t382 * t424 + t426;
t357 = qJD(4) * t397 + t381 * t411;
t355 = qJD(1) * t369 + qJD(2) * t368;
t354 = qJD(1) * t368 + qJD(2) * t369;
t353 = -qJD(1) * t406 - t392 * t417 + t403 * t421;
t348 = t381 * t412 - qJD(4) * t426 + (qJD(4) * t424 - t354) * t382;
t347 = qJD(4) * t400 - t354 * t381 - t382 * t412;
t346 = t347 * t387 - t353 * t390 + (t362 * t390 - t369 * t387) * qJD(6);
t345 = t347 * t390 + t353 * t387 + (-t362 * t387 - t369 * t390) * qJD(6);
t1 = [-t367 * qJD(3) - t436 * qJD(5) - t356 * t380 - t402 * t349 - t399 * t355 + (t438 * r_i_i_C(1) + t439 * r_i_i_C(2)) * qJD(6) + t416 * t435 + (-t392 * pkin(1) - t414 * t424) * qJD(1), t433 * t353 - t399 * t354 + t393 * t369 + t396 * t370, -t353, -t416 * t347 + t402 * t348 + t395 * t400, t347, t345 * r_i_i_C(1) - t346 * r_i_i_C(2); t346 * r_i_i_C(1) + t345 * r_i_i_C(2) + t347 * qJ(5) + t369 * qJD(3) + t362 * qJD(5) - t354 * t380 + t432 * t353 + t416 * t348 + (-pkin(1) * t389 + t414 * t422) * qJD(1), -t355 * t433 + t356 * t399 + t367 * t393 + t368 * t396, t355, -t416 * t349 + t395 * t440 - t402 * t435, t349 (t349 * t390 - t355 * t387) * r_i_i_C(1) + (-t349 * t387 - t355 * t390) * r_i_i_C(2) + (-t439 * r_i_i_C(1) + t438 * r_i_i_C(2)) * qJD(6); 0 ((-qJD(2) * t433 + t396) * t388 + (t399 * qJD(2) - t393) * t391) * t385, t410, t395 * t397 + t402 * (-qJD(4) * t365 + t382 * t411) - t416 * t357, t357 (t357 * t390 - t387 * t410) * r_i_i_C(1) + (-t357 * t387 - t390 * t410) * r_i_i_C(2) + ((-t365 * t387 + t390 * t423) * r_i_i_C(1) + (-t365 * t390 - t387 * t423) * r_i_i_C(2)) * qJD(6);];
JaD_transl  = t1;
