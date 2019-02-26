% Zeitableitung der translatorischen Teilmatrix der analytischen Jacobi-Matrix für Segment Nr. 6 (0=Basis) von
% S6RRRPPR10
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
%   pkin=[a2,a3,a4,a5,a6,alpha2,d1,d2,d3,d6,theta5]';
%
% Output:
% JaD_transl [3x6]
%   Zeitableitung der translatorischen Teilmatrix der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-02-26 22:08
% Revision: d75aae1ac561373cd3be920984c3df29a1c2ecc8 (2019-02-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JaD_transl = S6RRRPPR10_jacobiaD_transl_6_sym_varpar(qJ, qJD, r_i_i_C, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(3,1),zeros(11,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRPPR10_jacobiaD_transl_6_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRRPPR10_jacobiaD_transl_6_sym_varpar: qJD has to be [6x1] (double)');
assert(isa(r_i_i_C,'double') && isreal(r_i_i_C) && all(size(r_i_i_C) == [3 1]), ...
	'S6RRRPPR10_jacobiaD_transl_6_sym_varpar: Position vector r_i_i_C has to be [3x1] double');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6RRRPPR10_jacobiaD_transl_6_sym_varpar: pkin has to be [11x1] (double)');

%% Symbolic Calculation
% From jacobiaD_transl_6_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 22:08:53
% EndTime: 2019-02-26 22:08:54
% DurationCPUTime: 0.50s
% Computational Cost: add. (727->102), mult. (1878->165), div. (0->0), fcn. (1892->12), ass. (0->68)
t386 = cos(qJ(2));
t387 = cos(qJ(1));
t427 = cos(pkin(6));
t403 = t387 * t427;
t401 = t386 * t403;
t383 = sin(qJ(2));
t384 = sin(qJ(1));
t417 = t384 * t383;
t362 = -t401 + t417;
t378 = pkin(11) + qJ(6);
t376 = sin(t378);
t377 = cos(t378);
t363 = t383 * t403 + t384 * t386;
t382 = sin(qJ(3));
t385 = cos(qJ(3));
t380 = sin(pkin(6));
t419 = t380 * t387;
t432 = t363 * t382 + t385 * t419;
t435 = t362 * t377 + t376 * t432;
t434 = t362 * t376 - t377 * t432;
t400 = t377 * r_i_i_C(1) - t376 * r_i_i_C(2);
t390 = t400 * qJD(6) + qJD(4);
t399 = -t376 * r_i_i_C(1) - t377 * r_i_i_C(2);
t405 = sin(pkin(11)) * pkin(5) + qJ(4);
t391 = -t399 + t405;
t411 = r_i_i_C(3) + pkin(10) + qJ(5) + pkin(3);
t388 = (t411 * t382 - t391 * t385) * qJD(3) - t385 * qJD(5) - t390 * t382;
t398 = qJD(2) * t427 + qJD(1);
t404 = t384 * t427;
t402 = t383 * t404;
t415 = qJD(2) * t383;
t416 = t387 * t386;
t351 = -qJD(1) * t402 - t384 * t415 + t398 * t416;
t422 = t380 * t384;
t410 = t382 * t422;
t431 = -qJD(1) * t410 + qJD(3) * t432 - t351 * t385;
t429 = t391 * t382 + t411 * t385 + pkin(2);
t428 = -pkin(9) - cos(pkin(11)) * pkin(5) - pkin(4);
t365 = -t402 + t416;
t423 = t365 * t382;
t421 = t380 * t385;
t420 = t380 * t386;
t418 = t382 * t387;
t414 = qJD(2) * t386;
t413 = qJD(3) * t385;
t409 = t380 * t418;
t408 = qJD(1) * t421;
t407 = t380 * t415;
t406 = t380 * t414;
t396 = t363 * t385 - t409;
t358 = t365 * t385 + t410;
t395 = t399 * qJD(6);
t394 = t400 - t428;
t364 = t387 * t383 + t386 * t404;
t360 = t380 * t383 * t382 - t427 * t385;
t392 = t427 * t382 + t383 * t421;
t344 = -qJD(3) * t409 + t351 * t382 + t363 * t413 - t384 * t408;
t357 = -t384 * t421 + t423;
t353 = -t360 * qJD(3) + t385 * t406;
t352 = t392 * qJD(3) + t382 * t406;
t350 = t364 * qJD(1) + t363 * qJD(2);
t349 = t363 * qJD(1) + t364 * qJD(2);
t348 = -qJD(1) * t401 - t387 * t414 + t398 * t417;
t343 = -t349 * t385 - qJD(3) * t423 + (qJD(1) * t418 + t384 * t413) * t380;
t342 = t358 * qJD(3) - t349 * t382 - t387 * t408;
t341 = t342 * t376 - t348 * t377 + (t357 * t377 - t364 * t376) * qJD(6);
t340 = t342 * t377 + t348 * t376 + (-t357 * t376 - t364 * t377) * qJD(6);
t1 = [-t396 * qJD(5) - t432 * qJD(4) - t351 * pkin(2) - t394 * t350 - t391 * t344 + (t434 * r_i_i_C(1) + t435 * r_i_i_C(2)) * qJD(6) + (-t387 * pkin(1) - pkin(8) * t422) * qJD(1) + t411 * t431, t429 * t348 - t394 * t349 + t388 * t364 + t365 * t395, -t357 * qJD(5) - t411 * t342 + t391 * t343 + t390 * t358, t342, t343, t340 * r_i_i_C(1) - t341 * r_i_i_C(2); -t349 * pkin(2) + t341 * r_i_i_C(1) + t340 * r_i_i_C(2) + t357 * qJD(4) + t358 * qJD(5) + t428 * t348 + t405 * t342 + (-pkin(1) * t384 + pkin(8) * t419) * qJD(1) + t411 * t343, -t350 * t429 + t394 * t351 + t388 * t362 + t363 * t395, -qJD(5) * t432 - t411 * t344 + t390 * t396 - t391 * t431, t344, -t431 (t344 * t377 - t350 * t376) * r_i_i_C(1) + (-t344 * t376 - t350 * t377) * r_i_i_C(2) + (-t435 * r_i_i_C(1) + t434 * r_i_i_C(2)) * qJD(6); 0 ((-qJD(2) * t429 + t395) * t383 + (t394 * qJD(2) - t388) * t386) * t380, -t360 * qJD(5) - t411 * t352 + t391 * t353 + t390 * t392, t352, t353 (t352 * t377 - t376 * t407) * r_i_i_C(1) + (-t352 * t376 - t377 * t407) * r_i_i_C(2) + ((-t360 * t376 + t377 * t420) * r_i_i_C(1) + (-t360 * t377 - t376 * t420) * r_i_i_C(2)) * qJD(6);];
JaD_transl  = t1;
