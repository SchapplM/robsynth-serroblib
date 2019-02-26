% Zeitableitung der translatorischen Teilmatrix der analytischen Jacobi-Matrix für Segment Nr. 6 (0=Basis) von
% S6RRRRPR5
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
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d3,d4,d6]';
%
% Output:
% JaD_transl [3x6]
%   Zeitableitung der translatorischen Teilmatrix der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-02-26 22:33
% Revision: d75aae1ac561373cd3be920984c3df29a1c2ecc8 (2019-02-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JaD_transl = S6RRRRPR5_jacobiaD_transl_6_sym_varpar(qJ, qJD, r_i_i_C, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(3,1),zeros(10,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRRPR5_jacobiaD_transl_6_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRRRPR5_jacobiaD_transl_6_sym_varpar: qJD has to be [6x1] (double)');
assert(isa(r_i_i_C,'double') && isreal(r_i_i_C) && all(size(r_i_i_C) == [3 1]), ...
	'S6RRRRPR5_jacobiaD_transl_6_sym_varpar: Position vector r_i_i_C has to be [3x1] double');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RRRRPR5_jacobiaD_transl_6_sym_varpar: pkin has to be [10x1] (double)');

%% Symbolic Calculation
% From jacobiaD_transl_6_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 22:33:04
% EndTime: 2019-02-26 22:33:05
% DurationCPUTime: 0.80s
% Computational Cost: add. (972->126), mult. (1504->185), div. (0->0), fcn. (1371->10), ass. (0->81)
t368 = sin(qJ(4));
t372 = cos(qJ(4));
t374 = cos(qJ(1));
t422 = qJD(4) * t374;
t370 = sin(qJ(1));
t425 = qJD(1) * t370;
t388 = t368 * t422 + t372 * t425;
t448 = pkin(4) + pkin(5);
t447 = pkin(10) + r_i_i_C(3);
t367 = sin(qJ(6));
t371 = cos(qJ(6));
t397 = t367 * t368 + t371 * t372;
t441 = qJD(4) - qJD(6);
t379 = t441 * t397;
t398 = t367 * t372 - t368 * t371;
t446 = t441 * t398;
t423 = qJD(4) * t372;
t393 = -qJ(5) * t423 - qJD(5) * t368;
t366 = qJ(2) + qJ(3);
t364 = cos(t366);
t427 = t372 * t374;
t429 = t368 * t370;
t328 = t364 * t429 + t427;
t426 = t374 * t368;
t428 = t370 * t372;
t329 = t364 * t428 - t426;
t401 = t328 * t367 + t329 * t371;
t402 = t328 * t371 - t329 * t367;
t443 = (t401 * r_i_i_C(1) + t402 * r_i_i_C(2)) * qJD(6);
t363 = sin(t366);
t365 = qJD(2) + qJD(3);
t369 = sin(qJ(2));
t433 = pkin(2) * qJD(2);
t419 = t369 * t433;
t420 = pkin(9) - t447;
t442 = -t419 + (-pkin(3) * t363 + t420 * t364) * t365;
t436 = pkin(2) * t369;
t434 = pkin(9) * t364;
t432 = qJ(5) * t368;
t424 = qJD(1) * t374;
t387 = t368 * t424 + t370 * t423;
t430 = t365 * t370;
t417 = t363 * t430;
t326 = t387 * t364 - t368 * t417 - t388;
t323 = t326 * t371;
t431 = t364 * t365;
t418 = t448 * t368;
t416 = t374 * t365 * t363;
t413 = t364 * t425;
t411 = qJD(4) * t429;
t409 = t372 * t422;
t408 = -t367 * r_i_i_C(1) - qJ(5);
t407 = r_i_i_C(2) * t367 - t448;
t392 = t398 * t365;
t384 = t364 * t392;
t403 = -r_i_i_C(1) * (-t379 * t363 + t384) - (t363 * t446 + t397 * t431) * r_i_i_C(2);
t330 = t364 * t426 - t428;
t331 = t364 * t427 + t429;
t400 = t330 * t371 - t331 * t367;
t399 = t330 * t367 + t331 * t371;
t396 = t371 * r_i_i_C(2) - t408;
t395 = r_i_i_C(1) * t371 - t407;
t394 = -t372 * t448 - t432;
t391 = t397 * t365;
t390 = -pkin(3) + t394;
t385 = t364 * t391;
t389 = (t370 * t384 + (-t379 * t370 + t398 * t424) * t363) * r_i_i_C(2) + (-t370 * t385 + (-t370 * t446 - t397 * t424) * t363) * r_i_i_C(1) + t424 * t434 + t447 * t417 + t448 * t363 * t411;
t373 = cos(qJ(2));
t386 = -pkin(2) * t373 - pkin(3) * t364 - t420 * t363 - pkin(1);
t383 = t447 * (t416 + t413) + (-t385 * r_i_i_C(1) + t384 * r_i_i_C(2)) * t374 + ((-t379 * t374 - t398 * t425) * r_i_i_C(2) + (-t374 * t446 + t397 * t425) * r_i_i_C(1) + (t432 + pkin(3)) * t425 + t448 * t388) * t363;
t382 = t390 * t363 - t364 * t447;
t378 = t393 * t363 + (-pkin(9) * t363 + t390 * t364) * t365;
t377 = -t373 * t433 + t378;
t376 = pkin(9) * t431 + t382 * t365 + (-t391 * r_i_i_C(1) + t392 * r_i_i_C(2)) * t363 + (r_i_i_C(1) * t446 + r_i_i_C(2) * t379 - qJD(4) * t418 - t393) * t364;
t375 = -pkin(8) - pkin(7);
t327 = t331 * qJD(1) - t364 * t411 - t372 * t417 - t409;
t325 = t388 * t364 + t372 * t416 - t387;
t324 = t328 * qJD(1) - t364 * t409 + t368 * t416 - t411;
t314 = t400 * qJD(6) - t324 * t367 - t325 * t371;
t313 = -t399 * qJD(6) - t324 * t371 + t325 * t367;
t1 = [-t323 * r_i_i_C(2) - t328 * qJD(5) + t408 * t326 - t395 * t327 + (-t402 * r_i_i_C(1) + t401 * r_i_i_C(2)) * qJD(6) - t442 * t370 + (t370 * t375 + t386 * t374) * qJD(1), t383 + t377 * t374 + (-t434 + t436) * t425, -pkin(9) * t413 + t378 * t374 + t383, t331 * qJD(5) - t396 * t325 + t395 * t324 + (t399 * r_i_i_C(1) + t400 * r_i_i_C(2)) * qJD(6), -t324, r_i_i_C(1) * t313 - r_i_i_C(2) * t314; t314 * r_i_i_C(1) + t313 * r_i_i_C(2) - t324 * qJ(5) + t330 * qJD(5) - t448 * t325 + t442 * t374 + (t386 * t370 - t374 * t375) * qJD(1) (t382 - t436) * t424 + t377 * t370 + t389 (t390 * t430 - t424 * t447) * t364 + ((-pkin(9) * t365 + t393) * t370 + t390 * t424) * t363 + t389, -t323 * r_i_i_C(1) + t329 * qJD(5) + t407 * t326 + t396 * t327 + t443, t326 (-t327 * t367 + t323) * r_i_i_C(1) + (-t326 * t367 - t327 * t371) * r_i_i_C(2) - t443; 0, t376 - t419, t376 (qJ(5) * t372 - t418) * t431 + (t394 * qJD(4) + qJD(5) * t372) * t363 - t403, t363 * t423 + t368 * t431, t403;];
JaD_transl  = t1;
