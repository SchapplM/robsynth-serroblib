% Zeitableitung der translatorischen Teilmatrix der analytischen Jacobi-Matrix für Segment Nr. 6 (0=Basis) von
% S6RRPRRR7
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
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d4,d5,d6]';
%
% Output:
% JaD_transl [3x6]
%   Zeitableitung der translatorischen Teilmatrix der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-02-26 21:57
% Revision: d75aae1ac561373cd3be920984c3df29a1c2ecc8 (2019-02-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JaD_transl = S6RRPRRR7_jacobiaD_transl_6_sym_varpar(qJ, qJD, r_i_i_C, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(3,1),zeros(10,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPRRR7_jacobiaD_transl_6_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRPRRR7_jacobiaD_transl_6_sym_varpar: qJD has to be [6x1] (double)');
assert(isa(r_i_i_C,'double') && isreal(r_i_i_C) && all(size(r_i_i_C) == [3 1]), ...
	'S6RRPRRR7_jacobiaD_transl_6_sym_varpar: Position vector r_i_i_C has to be [3x1] double');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RRPRRR7_jacobiaD_transl_6_sym_varpar: pkin has to be [10x1] (double)');

%% Symbolic Calculation
% From jacobiaD_transl_6_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 21:57:37
% EndTime: 2019-02-26 21:57:37
% DurationCPUTime: 0.55s
% Computational Cost: add. (572->90), mult. (1353->141), div. (0->0), fcn. (1269->10), ass. (0->61)
t382 = sin(qJ(2));
t385 = cos(qJ(2));
t432 = sin(qJ(4));
t433 = cos(qJ(4));
t358 = t382 * t433 - t385 * t432;
t379 = qJD(5) + qJD(6);
t381 = sin(qJ(5));
t380 = qJ(5) + qJ(6);
t377 = sin(t380);
t378 = cos(t380);
t430 = t378 * r_i_i_C(2);
t400 = t377 * r_i_i_C(1) + t430;
t428 = pkin(5) * qJD(5);
t389 = t400 * t379 + t381 * t428;
t444 = t358 * t389;
t357 = t382 * t432 + t385 * t433;
t414 = qJD(2) * t432;
t415 = qJD(2) * t433;
t350 = t357 * qJD(4) - t382 * t414 - t385 * t415;
t383 = sin(qJ(1));
t386 = cos(qJ(1));
t439 = t358 * qJD(1);
t348 = t350 * t383 - t386 * t439;
t412 = qJD(4) * t432;
t413 = qJD(4) * t433;
t351 = (-t412 + t414) * t385 + (t413 - t415) * t382;
t394 = qJD(1) * t357;
t349 = t351 * t383 + t386 * t394;
t384 = cos(qJ(5));
t376 = t384 * pkin(5) + pkin(4);
t398 = t378 * r_i_i_C(1) - t377 * r_i_i_C(2) + t376;
t429 = r_i_i_C(3) + pkin(10) + pkin(9);
t443 = -t398 * t348 + t429 * t349 - t383 * t444;
t346 = -t351 * t386 + t383 * t394;
t356 = t357 * t386;
t347 = qJD(2) * t356 + (-t382 * t412 - t385 * t413) * t386 - t439 * t383;
t442 = -t429 * t346 + t398 * t347 - t386 * t444;
t441 = -t429 * t350 - t398 * t351 + t389 * t357;
t434 = pkin(2) + pkin(3);
t438 = -qJ(3) * t385 + t434 * t382;
t440 = t438 * qJD(2) - t382 * qJD(3);
t431 = pkin(5) * t381;
t426 = t377 * t379;
t425 = t378 * t379;
t421 = qJD(1) * t386;
t399 = -t356 * t379 - t421;
t411 = t379 * t383 + t346;
t344 = t411 * t377 + t399 * t378;
t345 = t399 * t377 - t411 * t378;
t424 = t344 * r_i_i_C(1) - t345 * r_i_i_C(2);
t354 = t357 * t383;
t422 = qJD(1) * t383;
t369 = t377 * t422;
t410 = -t379 * t386 - t349;
t423 = ((-t354 * t379 - t422) * t378 + t410 * t377) * r_i_i_C(1) + (t354 * t426 + t410 * t378 + t369) * r_i_i_C(2);
t420 = qJD(2) * t386;
t409 = pkin(7) - pkin(8) - t431;
t397 = -qJ(3) * t382 - t434 * t385;
t395 = -pkin(1) + t397;
t352 = t358 * r_i_i_C(2) * t426;
t1 = [t369 * r_i_i_C(1) - t398 * t349 + ((t354 * t377 - t378 * t386) * r_i_i_C(1) + (t354 * t378 + t377 * t386) * r_i_i_C(2)) * t379 - t429 * t348 + (t354 * t381 - t384 * t386) * t428 + t440 * t383 + (t395 * t386 + (-t409 + t430) * t383) * qJD(1) (-qJ(3) * t420 + t434 * t422) * t382 + (-qJ(3) * t422 + (-t434 * qJD(2) + qJD(3)) * t386) * t385 - t442, -t382 * t422 + t385 * t420, t442 (-t384 * t421 + t346 * t381 + (-t356 * t384 + t381 * t383) * qJD(5)) * pkin(5) + t424, t424; t345 * r_i_i_C(1) + t344 * r_i_i_C(2) - t346 * t376 - t429 * t347 + (-t356 * t381 - t383 * t384) * t428 - t440 * t386 + (t395 * t383 + t409 * t386) * qJD(1), -t438 * t421 + (t397 * qJD(2) + qJD(3) * t385) * t383 - t443, t383 * qJD(2) * t385 + t382 * t421, t443 (-t384 * t422 - t349 * t381 + (-t354 * t384 - t381 * t386) * qJD(5)) * pkin(5) + t423, t423; 0, -t440 - t441, qJD(2) * t382, t441, t352 + (-r_i_i_C(1) * t425 - t384 * t428) * t358 + (t400 + t431) * t350, t350 * t430 + t352 + (t350 * t377 - t358 * t425) * r_i_i_C(1);];
JaD_transl  = t1;
