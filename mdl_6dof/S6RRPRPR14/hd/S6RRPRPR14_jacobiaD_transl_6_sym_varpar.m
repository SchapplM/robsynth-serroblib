% Zeitableitung der translatorischen Teilmatrix der analytischen Jacobi-Matrix für Segment Nr. 6 (0=Basis) von
% S6RRPRPR14
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
%   pkin=[a2,a3,a4,a5,a6,alpha2,d1,d2,d4,d6]';
%
% Output:
% JaD_transl [3x6]
%   Zeitableitung der translatorischen Teilmatrix der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-02-26 21:45
% Revision: d75aae1ac561373cd3be920984c3df29a1c2ecc8 (2019-02-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JaD_transl = S6RRPRPR14_jacobiaD_transl_6_sym_varpar(qJ, qJD, r_i_i_C, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(3,1),zeros(10,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPRPR14_jacobiaD_transl_6_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRPRPR14_jacobiaD_transl_6_sym_varpar: qJD has to be [6x1] (double)');
assert(isa(r_i_i_C,'double') && isreal(r_i_i_C) && all(size(r_i_i_C) == [3 1]), ...
	'S6RRPRPR14_jacobiaD_transl_6_sym_varpar: Position vector r_i_i_C has to be [3x1] double');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RRPRPR14_jacobiaD_transl_6_sym_varpar: pkin has to be [10x1] (double)');

%% Symbolic Calculation
% From jacobiaD_transl_6_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 21:45:28
% EndTime: 2019-02-26 21:45:28
% DurationCPUTime: 0.46s
% Computational Cost: add. (579->103), mult. (1713->165), div. (0->0), fcn. (1716->10), ass. (0->67)
t380 = cos(pkin(6));
t387 = cos(qJ(2));
t388 = cos(qJ(1));
t417 = t388 * t387;
t407 = t380 * t417;
t383 = sin(qJ(2));
t384 = sin(qJ(1));
t420 = t384 * t383;
t365 = -t407 + t420;
t382 = sin(qJ(4));
t386 = cos(qJ(4));
t379 = sin(pkin(6));
t421 = t379 * t388;
t357 = t365 * t386 + t382 * t421;
t418 = t388 * t383;
t419 = t384 * t387;
t366 = t380 * t418 + t419;
t381 = sin(qJ(6));
t385 = cos(qJ(6));
t431 = t357 * t381 - t366 * t385;
t430 = t357 * t385 + t366 * t381;
t397 = t381 * r_i_i_C(1) + t385 * r_i_i_C(2);
t394 = qJD(6) * t397;
t398 = t385 * r_i_i_C(1) - t381 * r_i_i_C(2);
t391 = t398 * qJD(6) + qJD(5);
t396 = qJ(5) + t397;
t410 = -r_i_i_C(3) - pkin(10) - pkin(4);
t429 = t410 * t382 + t396 * t386 - qJ(3);
t428 = pkin(3) + pkin(8);
t367 = t380 * t419 + t418;
t351 = t367 * qJD(1) + t366 * qJD(2);
t427 = t351 * t382;
t424 = t379 * t383;
t423 = t379 * t384;
t422 = t379 * t387;
t416 = qJD(1) * t384;
t415 = qJD(2) * t383;
t414 = qJD(2) * t387;
t413 = qJD(4) * t382;
t412 = qJD(4) * t386;
t411 = -pkin(2) - pkin(5) - pkin(9);
t409 = t380 * t420;
t408 = t386 * t421;
t406 = t379 * t416;
t405 = qJD(1) * t421;
t404 = t386 * t416;
t403 = t379 * t414;
t402 = t379 * t413;
t401 = t379 * t415;
t400 = -qJD(4) * t408 - t351 * t386;
t399 = qJD(2) * t380 + qJD(1);
t395 = t367 * t382 + t386 * t423;
t363 = t380 * t382 + t386 * t422;
t393 = t398 - t411;
t389 = qJD(3) - t391 * t386 + (t396 * t382 - t410 * t386) * qJD(4);
t368 = -t409 + t417;
t355 = -t367 * t386 + t382 * t423;
t354 = t380 * t412 - t386 * t401 - t387 * t402;
t352 = -qJD(1) * t409 - t384 * t415 + t399 * t417;
t350 = t366 * qJD(1) + t367 * qJD(2);
t349 = -qJD(1) * t407 - t388 * t414 + t399 * t420;
t345 = -t349 * t382 - t384 * t402 + (qJD(4) * t367 + t405) * t386;
t344 = t395 * qJD(4) + t349 * t386 + t382 * t405;
t340 = (qJD(4) * t365 + t406) * t382 + t400;
t339 = t344 * t381 - t350 * t385 + (t355 * t385 - t368 * t381) * qJD(6);
t338 = t344 * t385 + t350 * t381 + (-t355 * t381 - t368 * t385) * qJD(6);
t1 = [-t351 * qJ(3) - t365 * qJD(3) + t357 * qJD(5) - t396 * (t365 * t413 + t382 * t406 + t400) + (t430 * r_i_i_C(1) - t431 * r_i_i_C(2)) * qJD(6) - t393 * t352 + t410 * (t365 * t412 + t427 + (t388 * t413 + t404) * t379) + (-t388 * pkin(1) - t428 * t423) * qJD(1), t393 * t349 + t429 * t350 + t367 * t394 + t389 * t368, -t349, t410 * t344 + t396 * t345 + t391 * t395, t344, t338 * r_i_i_C(1) - t339 * r_i_i_C(2); t339 * r_i_i_C(1) + t338 * r_i_i_C(2) - t349 * qJ(3) + t344 * qJ(5) + t367 * qJD(3) + t355 * qJD(5) + t411 * t350 - t410 * t345 + (-pkin(1) * t384 + t428 * t421) * qJD(1), -t393 * t351 - t352 * t429 + t365 * t394 + t389 * t366, t351, -t391 * (-t365 * t382 + t408) + t396 * (t357 * qJD(4) + t379 * t404 + t427) + t410 * t340, t340 (t340 * t385 - t352 * t381) * r_i_i_C(1) + (-t340 * t381 - t352 * t385) * r_i_i_C(2) + (t431 * r_i_i_C(1) + t430 * r_i_i_C(2)) * qJD(6); 0 ((-qJD(2) * t429 - t394) * t387 + (-t393 * qJD(2) + t389) * t383) * t379, t401, t391 * (t380 * t386 - t382 * t422) + t410 * t354 - t396 * (t363 * qJD(4) - t382 * t401) t354 (t354 * t385 - t381 * t403) * r_i_i_C(1) + (-t354 * t381 - t385 * t403) * r_i_i_C(2) + ((-t363 * t381 - t385 * t424) * r_i_i_C(1) + (-t363 * t385 + t381 * t424) * r_i_i_C(2)) * qJD(6);];
JaD_transl  = t1;
