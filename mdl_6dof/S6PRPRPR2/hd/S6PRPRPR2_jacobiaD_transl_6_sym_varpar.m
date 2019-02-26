% Zeitableitung der translatorischen Teilmatrix der analytischen Jacobi-Matrix für Segment Nr. 6 (0=Basis) von
% S6PRPRPR2
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
%   pkin=[a2,a3,a4,a5,a6,alpha2,d2,d4,d6,theta1,theta3,theta5]';
%
% Output:
% JaD_transl [3x6]
%   Zeitableitung der translatorischen Teilmatrix der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-02-26 19:47
% Revision: d75aae1ac561373cd3be920984c3df29a1c2ecc8 (2019-02-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JaD_transl = S6PRPRPR2_jacobiaD_transl_6_sym_varpar(qJ, qJD, r_i_i_C, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(3,1),zeros(12,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRPRPR2_jacobiaD_transl_6_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6PRPRPR2_jacobiaD_transl_6_sym_varpar: qJD has to be [6x1] (double)');
assert(isa(r_i_i_C,'double') && isreal(r_i_i_C) && all(size(r_i_i_C) == [3 1]), ...
	'S6PRPRPR2_jacobiaD_transl_6_sym_varpar: Position vector r_i_i_C has to be [3x1] double');
assert(isreal(pkin) && all(size(pkin) == [12 1]), ...
  'S6PRPRPR2_jacobiaD_transl_6_sym_varpar: pkin has to be [12x1] (double)');

%% Symbolic Calculation
% From jacobiaD_transl_6_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 19:47:05
% EndTime: 2019-02-26 19:47:05
% DurationCPUTime: 0.49s
% Computational Cost: add. (566->80), mult. (1521->147), div. (0->0), fcn. (1646->14), ass. (0->62)
t395 = sin(qJ(4));
t397 = cos(qJ(4));
t387 = pkin(12) + qJ(6);
t385 = sin(t387);
t386 = cos(t387);
t413 = t385 * r_i_i_C(1) + t386 * r_i_i_C(2);
t405 = qJD(6) * t413;
t414 = t386 * r_i_i_C(1) - t385 * r_i_i_C(2);
t410 = cos(pkin(12)) * pkin(5) + pkin(4) + t414;
t429 = r_i_i_C(3) + pkin(9) + qJ(5);
t434 = (t410 * t395 - t429 * t397) * qJD(4) - t395 * qJD(5) + t397 * t405;
t393 = cos(pkin(6));
t389 = sin(pkin(11));
t396 = sin(qJ(2));
t427 = cos(pkin(11));
t430 = cos(qJ(2));
t403 = t430 * t389 + t396 * t427;
t372 = t403 * t393;
t415 = t430 * t427;
t421 = qJD(2) * t396;
t432 = -qJD(2) * t415 + t389 * t421;
t402 = -t396 * t389 + t415;
t390 = sin(pkin(10));
t392 = cos(pkin(10));
t358 = t392 * t372 + t390 * t402;
t391 = sin(pkin(6));
t425 = t391 * t395;
t347 = t358 * t397 - t392 * t425;
t399 = t429 * t395 + t410 * t397 + pkin(3);
t428 = pkin(2) * qJD(2);
t424 = t391 * t397;
t423 = t393 * t396;
t369 = t432 * t393;
t374 = t403 * qJD(2);
t412 = t392 * t369 + t390 * t374;
t354 = t390 * t369 - t392 * t374;
t371 = t403 * t391;
t363 = t371 * t397 + t393 * t395;
t411 = -t371 * t395 + t393 * t397;
t361 = -t390 * t372 + t392 * t402;
t408 = -t358 * t395 - t392 * t424;
t407 = -t361 * t395 + t390 * t424;
t349 = t361 * t397 + t390 * t425;
t406 = qJD(6) * t414;
t404 = -sin(pkin(12)) * pkin(5) - pkin(8) - t413;
t401 = t402 * t393;
t400 = qJD(2) * t372;
t373 = t402 * qJD(2);
t370 = t402 * t391;
t368 = qJD(2) * t371;
t367 = t432 * t391;
t360 = -t390 * t401 - t392 * t403;
t357 = -t390 * t403 + t392 * t401;
t353 = -t392 * t373 + t390 * t400;
t350 = -t390 * t373 - t392 * t400;
t345 = t411 * qJD(4) - t367 * t397;
t344 = t363 * qJD(4) - t367 * t395;
t343 = t407 * qJD(4) + t354 * t397;
t342 = t349 * qJD(4) + t354 * t395;
t341 = t408 * qJD(4) - t397 * t412;
t340 = t347 * qJD(4) - t395 * t412;
t1 = [0, t361 * t406 - t404 * t354 + (t390 * t423 - t430 * t392) * t428 + t399 * t353 - t434 * t360, 0, t349 * qJD(5) - t410 * t342 + t429 * t343 - t407 * t405, t342 (-t343 * t385 - t353 * t386) * r_i_i_C(1) + (-t343 * t386 + t353 * t385) * r_i_i_C(2) + ((-t349 * t386 + t360 * t385) * r_i_i_C(1) + (t349 * t385 + t360 * t386) * r_i_i_C(2)) * qJD(6); 0, t358 * t406 + t404 * t412 + (-t430 * t390 - t392 * t423) * t428 + t399 * t350 - t434 * t357, 0, t347 * qJD(5) - t410 * t340 + t429 * t341 - t408 * t405, t340 (-t341 * t385 - t350 * t386) * r_i_i_C(1) + (-t341 * t386 + t350 * t385) * r_i_i_C(2) + ((-t347 * t386 + t357 * t385) * r_i_i_C(1) + (t347 * t385 + t357 * t386) * r_i_i_C(2)) * qJD(6); 0, -t391 * pkin(2) * t421 + t404 * t367 - t399 * t368 - t434 * t370 + t371 * t406, 0, t363 * qJD(5) - t410 * t344 + t429 * t345 - t411 * t405, t344 (-t345 * t385 + t368 * t386) * r_i_i_C(1) + (-t345 * t386 - t368 * t385) * r_i_i_C(2) + ((-t363 * t386 + t370 * t385) * r_i_i_C(1) + (t363 * t385 + t370 * t386) * r_i_i_C(2)) * qJD(6);];
JaD_transl  = t1;
