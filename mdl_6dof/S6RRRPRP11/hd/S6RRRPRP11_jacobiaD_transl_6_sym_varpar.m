% Zeitableitung der translatorischen Teilmatrix der analytischen Jacobi-Matrix für Segment Nr. 6 (0=Basis) von
% S6RRRPRP11
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
%   pkin=[a2,a3,a4,a5,a6,alpha2,d1,d2,d3,d5]';
%
% Output:
% JaD_transl [3x6]
%   Zeitableitung der translatorischen Teilmatrix der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-02-26 22:15
% Revision: d75aae1ac561373cd3be920984c3df29a1c2ecc8 (2019-02-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JaD_transl = S6RRRPRP11_jacobiaD_transl_6_sym_varpar(qJ, qJD, r_i_i_C, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(3,1),zeros(10,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRPRP11_jacobiaD_transl_6_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRRPRP11_jacobiaD_transl_6_sym_varpar: qJD has to be [6x1] (double)');
assert(isa(r_i_i_C,'double') && isreal(r_i_i_C) && all(size(r_i_i_C) == [3 1]), ...
	'S6RRRPRP11_jacobiaD_transl_6_sym_varpar: Position vector r_i_i_C has to be [3x1] double');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RRRPRP11_jacobiaD_transl_6_sym_varpar: pkin has to be [10x1] (double)');

%% Symbolic Calculation
% From jacobiaD_transl_6_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 22:15:00
% EndTime: 2019-02-26 22:15:01
% DurationCPUTime: 0.74s
% Computational Cost: add. (714->113), mult. (2093->176), div. (0->0), fcn. (2106->10), ass. (0->72)
t380 = sin(qJ(3));
t384 = cos(qJ(3));
t379 = sin(qJ(5));
t383 = cos(qJ(5));
t434 = pkin(5) + r_i_i_C(1);
t400 = -t383 * r_i_i_C(2) - t434 * t379;
t397 = -qJ(4) + t400;
t417 = pkin(3) + r_i_i_C(3) + qJ(6) + pkin(10);
t442 = -(t417 * t380 + t397 * t384) * qJD(3) + t384 * qJD(6);
t387 = t397 * t380 - t417 * t384 - pkin(2);
t385 = cos(qJ(2));
t429 = cos(pkin(6));
t433 = cos(qJ(1));
t405 = t429 * t433;
t401 = t385 * t405;
t381 = sin(qJ(2));
t382 = sin(qJ(1));
t422 = t382 * t381;
t363 = -t401 + t422;
t364 = t381 * t405 + t382 * t385;
t377 = sin(pkin(6));
t416 = t377 * t433;
t372 = t384 * t416;
t409 = t364 * t380 + t372;
t440 = t363 * t383 + t379 * t409;
t439 = t363 * t379 - t383 * t409;
t432 = t379 * r_i_i_C(2);
t389 = qJD(4) + (t434 * t383 - t432) * qJD(5);
t410 = t382 * t429;
t406 = t381 * t410;
t411 = t433 * qJD(1);
t421 = qJD(2) * t381;
t352 = -qJD(1) * t406 - t382 * t421 + (qJD(2) * t405 + t411) * t385;
t425 = t377 * t382;
t414 = qJD(1) * t425;
t435 = qJD(3) * t409 - t352 * t384 - t380 * t414;
t430 = -t383 * pkin(5) - pkin(4) - pkin(9);
t415 = t433 * t385;
t366 = t415 - t406;
t426 = t366 * t380;
t424 = t377 * t384;
t423 = t377 * t385;
t420 = qJD(3) * t384;
t419 = qJD(5) * t380;
t413 = t377 * t421;
t412 = qJD(2) * t423;
t408 = -t430 - t432;
t407 = t380 * t416;
t345 = -qJD(3) * t407 + t352 * t380 + t364 * t420 - t384 * t414;
t365 = t433 * t381 + t385 * t410;
t351 = qJD(1) * t365 + qJD(2) * t364;
t404 = t345 * t383 - t351 * t379;
t361 = t377 * t381 * t380 - t429 * t384;
t399 = -t361 * t379 + t383 * t423;
t359 = t366 * t384 + t380 * t425;
t398 = t383 * r_i_i_C(1) + t408;
t395 = t364 * t384 - t407;
t390 = t429 * t380 + t381 * t424;
t353 = qJD(3) * t390 + t380 * t412;
t392 = t353 * t383 - t379 * t413;
t391 = t400 * qJD(5);
t350 = qJD(1) * t364 + qJD(2) * t365;
t343 = -qJD(1) * t372 + t359 * qJD(3) - t350 * t380;
t358 = -t382 * t424 + t426;
t388 = t343 * t379 + (t358 * t383 - t365 * t379) * qJD(5);
t349 = -qJD(1) * t401 - qJD(2) * t415 + (qJD(2) * t429 + qJD(1)) * t422;
t341 = t343 * t383 + t349 * t379 + (-t358 * t379 - t365 * t383) * qJD(5);
t386 = -t389 * t380 - t442;
t354 = -qJD(3) * t361 + t384 * t412;
t344 = -t350 * t384 - qJD(3) * t426 + (t380 * t411 + t382 * t420) * t377;
t342 = -t349 * t383 + t388;
t1 = [-t395 * qJD(6) - t409 * qJD(4) - t352 * pkin(2) - t398 * t351 + t397 * t345 + (-t433 * pkin(1) - pkin(8) * t425) * qJD(1) + t417 * t435 + (t440 * r_i_i_C(2) + t434 * t439) * qJD(5), -t387 * t349 - t398 * t350 + t386 * t365 + t366 * t391, -t358 * qJD(6) - t417 * t343 - t344 * t397 + t389 * t359, t343, -t342 * r_i_i_C(2) + t434 * t341, t344; -t350 * pkin(2) + t342 * r_i_i_C(1) + t341 * r_i_i_C(2) + t343 * qJ(4) + t358 * qJD(4) + t359 * qJD(6) + t430 * t349 + (-pkin(1) * t382 + pkin(8) * t416) * qJD(1) + t417 * t344 + t388 * pkin(5), t351 * t387 + t352 * t398 + t363 * t386 + t364 * t391, -qJD(6) * t409 - t417 * t345 + t389 * t395 + t397 * t435, t345, t404 * r_i_i_C(1) + (-t345 * t379 - t351 * t383) * r_i_i_C(2) + (-r_i_i_C(1) * t440 + t439 * r_i_i_C(2)) * qJD(5) + (-qJD(5) * t440 + t404) * pkin(5), -t435; 0 ((qJD(2) * t387 + t391) * t381 + ((-qJD(5) * t432 + qJD(4)) * t380 + t408 * qJD(2) + ((qJD(2) + t419) * r_i_i_C(1) + pkin(5) * t419) * t383 + t442) * t385) * t377, -t361 * qJD(6) - t417 * t353 - t354 * t397 + t389 * t390, t353, t392 * r_i_i_C(1) + (-t353 * t379 - t383 * t413) * r_i_i_C(2) + (t399 * r_i_i_C(1) + (-t361 * t383 - t379 * t423) * r_i_i_C(2)) * qJD(5) + (t399 * qJD(5) + t392) * pkin(5), t354;];
JaD_transl  = t1;
