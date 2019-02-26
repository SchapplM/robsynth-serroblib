% Zeitableitung der translatorischen Teilmatrix der analytischen Jacobi-Matrix für Segment Nr. 6 (0=Basis) von
% S6PRPRPR3
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
%   pkin=[a2,a3,a4,a5,a6,alpha2,d2,d4,d6,theta1,theta3]';
%
% Output:
% JaD_transl [3x6]
%   Zeitableitung der translatorischen Teilmatrix der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-02-26 19:47
% Revision: d75aae1ac561373cd3be920984c3df29a1c2ecc8 (2019-02-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JaD_transl = S6PRPRPR3_jacobiaD_transl_6_sym_varpar(qJ, qJD, r_i_i_C, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(3,1),zeros(11,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRPRPR3_jacobiaD_transl_6_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6PRPRPR3_jacobiaD_transl_6_sym_varpar: qJD has to be [6x1] (double)');
assert(isa(r_i_i_C,'double') && isreal(r_i_i_C) && all(size(r_i_i_C) == [3 1]), ...
	'S6PRPRPR3_jacobiaD_transl_6_sym_varpar: Position vector r_i_i_C has to be [3x1] double');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6PRPRPR3_jacobiaD_transl_6_sym_varpar: pkin has to be [11x1] (double)');

%% Symbolic Calculation
% From jacobiaD_transl_6_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 19:47:40
% EndTime: 2019-02-26 19:47:41
% DurationCPUTime: 0.50s
% Computational Cost: add. (525->74), mult. (1644->141), div. (0->0), fcn. (1774->12), ass. (0->58)
t390 = sin(qJ(4));
t393 = cos(qJ(4));
t389 = sin(qJ(6));
t392 = cos(qJ(6));
t410 = t392 * r_i_i_C(1) - t389 * r_i_i_C(2);
t398 = t410 * qJD(6) + qJD(5);
t409 = -t389 * r_i_i_C(1) - t392 * r_i_i_C(2);
t406 = qJ(5) - t409;
t416 = pkin(4) + pkin(9) + r_i_i_C(3);
t430 = (t416 * t390 - t406 * t393) * qJD(4) - t398 * t390;
t388 = cos(pkin(6));
t384 = sin(pkin(11));
t391 = sin(qJ(2));
t423 = cos(pkin(11));
t425 = cos(qJ(2));
t400 = t425 * t384 + t391 * t423;
t372 = t400 * t388;
t411 = t425 * t423;
t417 = qJD(2) * t391;
t428 = -qJD(2) * t411 + t384 * t417;
t399 = -t391 * t384 + t411;
t385 = sin(pkin(10));
t387 = cos(pkin(10));
t358 = t372 * t387 + t385 * t399;
t386 = sin(pkin(6));
t421 = t386 * t390;
t427 = t358 * t393 - t387 * t421;
t395 = t406 * t390 + t416 * t393 + pkin(3);
t424 = pkin(2) * qJD(2);
t420 = t386 * t393;
t419 = t388 * t391;
t369 = t428 * t388;
t374 = t400 * qJD(2);
t408 = t369 * t387 + t374 * t385;
t354 = t385 * t369 - t374 * t387;
t371 = t400 * t386;
t407 = t371 * t393 + t388 * t390;
t362 = t371 * t390 - t388 * t393;
t361 = -t372 * t385 + t387 * t399;
t404 = -pkin(5) - pkin(8) - t410;
t346 = t358 * t390 + t387 * t420;
t403 = -t361 * t390 + t385 * t420;
t402 = t361 * t393 + t385 * t421;
t401 = qJD(6) * t409;
t397 = t399 * t388;
t396 = qJD(2) * t372;
t373 = t399 * qJD(2);
t370 = t399 * t386;
t368 = qJD(2) * t371;
t367 = t428 * t386;
t360 = -t385 * t397 - t387 * t400;
t357 = -t385 * t400 + t387 * t397;
t353 = -t373 * t387 + t385 * t396;
t350 = -t385 * t373 - t387 * t396;
t344 = t407 * qJD(4) - t367 * t390;
t342 = t402 * qJD(4) + t354 * t390;
t340 = t427 * qJD(4) - t390 * t408;
t1 = [0, t361 * t401 - t404 * t354 + (t385 * t419 - t425 * t387) * t424 + t395 * t353 - t430 * t360, 0, t398 * t402 + t406 * (t403 * qJD(4) + t354 * t393) - t416 * t342, t342 (t342 * t392 + t353 * t389) * r_i_i_C(1) + (-t342 * t389 + t353 * t392) * r_i_i_C(2) + ((t360 * t392 + t389 * t403) * r_i_i_C(1) + (-t360 * t389 + t392 * t403) * r_i_i_C(2)) * qJD(6); 0, t358 * t401 + t404 * t408 + (-t425 * t385 - t387 * t419) * t424 + t395 * t350 - t430 * t357, 0, t398 * t427 + t406 * (-t346 * qJD(4) - t393 * t408) - t416 * t340, t340 (t340 * t392 + t350 * t389) * r_i_i_C(1) + (-t340 * t389 + t350 * t392) * r_i_i_C(2) + ((-t346 * t389 + t357 * t392) * r_i_i_C(1) + (-t346 * t392 - t357 * t389) * r_i_i_C(2)) * qJD(6); 0, -t386 * pkin(2) * t417 + t404 * t367 - t395 * t368 - t430 * t370 + t371 * t401, 0, t398 * t407 + t406 * (-t362 * qJD(4) - t367 * t393) - t416 * t344, t344 (t344 * t392 - t368 * t389) * r_i_i_C(1) + (-t344 * t389 - t368 * t392) * r_i_i_C(2) + ((-t362 * t389 + t370 * t392) * r_i_i_C(1) + (-t362 * t392 - t370 * t389) * r_i_i_C(2)) * qJD(6);];
JaD_transl  = t1;
