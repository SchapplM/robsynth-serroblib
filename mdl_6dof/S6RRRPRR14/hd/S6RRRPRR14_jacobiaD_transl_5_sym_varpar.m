% Zeitableitung der translatorischen Teilmatrix der analytischen Jacobi-Matrix für Segment Nr. 5 (0=Basis) von
% S6RRRPRR14
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
%   pkin=[a2,a3,a4,a5,a6,alpha2,d1,d2,d3,d5,d6]';
%
% Output:
% JaD_transl [3x6]
%   Zeitableitung der translatorischen Teilmatrix der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-02-26 22:23
% Revision: d75aae1ac561373cd3be920984c3df29a1c2ecc8 (2019-02-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JaD_transl = S6RRRPRR14_jacobiaD_transl_5_sym_varpar(qJ, qJD, r_i_i_C, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(3,1),zeros(11,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRPRR14_jacobiaD_transl_5_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRRPRR14_jacobiaD_transl_5_sym_varpar: qJD has to be [6x1] (double)');
assert(isa(r_i_i_C,'double') && isreal(r_i_i_C) && all(size(r_i_i_C) == [3 1]), ...
	'S6RRRPRR14_jacobiaD_transl_5_sym_varpar: Position vector r_i_i_C has to be [3x1] double');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6RRRPRR14_jacobiaD_transl_5_sym_varpar: pkin has to be [11x1] (double)');

%% Symbolic Calculation
% From jacobiaD_transl_5_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 22:23:46
% EndTime: 2019-02-26 22:23:47
% DurationCPUTime: 0.48s
% Computational Cost: add. (545->92), mult. (1620->157), div. (0->0), fcn. (1626->10), ass. (0->64)
t376 = cos(qJ(2));
t377 = cos(qJ(1));
t415 = cos(pkin(6));
t393 = t377 * t415;
t391 = t376 * t393;
t372 = sin(qJ(2));
t373 = sin(qJ(1));
t405 = t372 * t373;
t356 = -t391 + t405;
t370 = sin(qJ(5));
t374 = cos(qJ(5));
t357 = t372 * t393 + t373 * t376;
t371 = sin(qJ(3));
t375 = cos(qJ(3));
t369 = sin(pkin(6));
t407 = t369 * t377;
t420 = t357 * t371 + t375 * t407;
t423 = t356 * t374 + t370 * t420;
t422 = t356 * t370 - t374 * t420;
t390 = r_i_i_C(1) * t374 - r_i_i_C(2) * t370;
t380 = t390 * qJD(5) + qJD(4);
t389 = -r_i_i_C(1) * t370 - r_i_i_C(2) * t374;
t387 = qJ(4) - t389;
t400 = pkin(3) + pkin(10) + r_i_i_C(3);
t378 = (t400 * t371 - t387 * t375) * qJD(3) - t380 * t371;
t388 = qJD(2) * t415 + qJD(1);
t394 = t373 * t415;
t392 = t372 * t394;
t403 = qJD(2) * t372;
t404 = t377 * t376;
t345 = -qJD(1) * t392 - t373 * t403 + t388 * t404;
t410 = t369 * t373;
t399 = t371 * t410;
t419 = -qJD(1) * t399 + qJD(3) * t420 - t345 * t375;
t417 = t387 * t371 + t400 * t375 + pkin(2);
t416 = -pkin(4) - pkin(9);
t359 = -t392 + t404;
t411 = t359 * t371;
t409 = t369 * t375;
t408 = t369 * t376;
t406 = t371 * t377;
t402 = qJD(2) * t376;
t401 = qJD(3) * t375;
t398 = t369 * t406;
t397 = qJD(1) * t409;
t396 = t369 * t403;
t395 = t369 * t402;
t386 = t390 - t416;
t384 = t359 * t375 + t399;
t383 = t389 * qJD(5);
t358 = t377 * t372 + t376 * t394;
t354 = t369 * t372 * t371 - t415 * t375;
t381 = t415 * t371 + t372 * t409;
t338 = -qJD(3) * t398 + t345 * t371 + t357 * t401 - t373 * t397;
t351 = -t373 * t409 + t411;
t346 = t381 * qJD(3) + t371 * t395;
t344 = t358 * qJD(1) + t357 * qJD(2);
t343 = t357 * qJD(1) + t358 * qJD(2);
t342 = -qJD(1) * t391 - t377 * t402 + t388 * t405;
t337 = -t343 * t375 - qJD(3) * t411 + (qJD(1) * t406 + t373 * t401) * t369;
t336 = t384 * qJD(3) - t343 * t371 - t377 * t397;
t335 = t336 * t370 - t342 * t374 + (t351 * t374 - t358 * t370) * qJD(5);
t334 = t336 * t374 + t342 * t370 + (-t351 * t370 - t358 * t374) * qJD(5);
t1 = [-t345 * pkin(2) - t420 * qJD(4) - t387 * t338 - t386 * t344 + (t422 * r_i_i_C(1) + t423 * r_i_i_C(2)) * qJD(5) + (-pkin(1) * t377 - pkin(8) * t410) * qJD(1) + t400 * t419, t417 * t342 - t386 * t343 + t378 * t358 + t359 * t383, -t400 * t336 + t387 * t337 + t380 * t384, t336, r_i_i_C(1) * t334 - t335 * r_i_i_C(2), 0; -t343 * pkin(2) + t335 * r_i_i_C(1) + t334 * r_i_i_C(2) + t336 * qJ(4) + t351 * qJD(4) + t416 * t342 + (-pkin(1) * t373 + pkin(8) * t407) * qJD(1) + t400 * t337, -t344 * t417 + t386 * t345 + t378 * t356 + t357 * t383, t380 * (t357 * t375 - t398) - t387 * t419 - t400 * t338, t338 (t338 * t374 - t344 * t370) * r_i_i_C(1) + (-t338 * t370 - t344 * t374) * r_i_i_C(2) + (-t423 * r_i_i_C(1) + t422 * r_i_i_C(2)) * qJD(5), 0; 0 ((-qJD(2) * t417 + t383) * t372 + (t386 * qJD(2) - t378) * t376) * t369, t380 * t381 + t387 * (-t354 * qJD(3) + t375 * t395) - t400 * t346, t346 (t346 * t374 - t370 * t396) * r_i_i_C(1) + (-t346 * t370 - t374 * t396) * r_i_i_C(2) + ((-t354 * t370 + t374 * t408) * r_i_i_C(1) + (-t354 * t374 - t370 * t408) * r_i_i_C(2)) * qJD(5), 0;];
JaD_transl  = t1;
