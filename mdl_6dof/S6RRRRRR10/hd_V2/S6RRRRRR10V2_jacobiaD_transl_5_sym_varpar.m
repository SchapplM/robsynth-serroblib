% Zeitableitung der translatorischen Teilmatrix der analytischen Jacobi-Matrix für Segment Nr. 5 (0=Basis) von
% S6RRRRRR10V2
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
% pkin [6x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,d1,d4,d6]';
%
% Output:
% JaD_transl [3x6]
%   Zeitableitung der translatorischen Teilmatrix der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-04-11 14:56
% Revision: b693519ea345eb34ae9622239e7f1167217e9d53 (2019-04-09)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JaD_transl = S6RRRRRR10V2_jacobiaD_transl_5_sym_varpar(qJ, qJD, r_i_i_C, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(3,1),zeros(6,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRRRR10V2_jacobiaD_transl_5_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRRRRR10V2_jacobiaD_transl_5_sym_varpar: qJD has to be [6x1] (double)');
assert(isa(r_i_i_C,'double') && isreal(r_i_i_C) && all(size(r_i_i_C) == [3 1]), ...
	'S6RRRRRR10V2_jacobiaD_transl_5_sym_varpar: Position vector r_i_i_C has to be [3x1] double');
assert(isreal(pkin) && all(size(pkin) == [6 1]), ...
  'S6RRRRRR10V2_jacobiaD_transl_5_sym_varpar: pkin has to be [6x1] (double)');

%% Symbolic Calculation
% From jacobiaD_transl_5_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-04-11 14:56:34
% EndTime: 2019-04-11 14:56:34
% DurationCPUTime: 0.67s
% Computational Cost: add. (562->105), mult. (860->177), div. (0->0), fcn. (777->10), ass. (0->81)
t365 = sin(qJ(5));
t369 = cos(qJ(5));
t366 = sin(qJ(4));
t413 = qJD(4) * t366;
t363 = qJD(2) + qJD(3);
t370 = cos(qJ(4));
t435 = -qJD(5) * t370 + t363;
t434 = -t365 * t435 + t369 * t413;
t437 = t365 * t413 + t369 * t435;
t429 = r_i_i_C(3) * t366;
t436 = pkin(3) + t429;
t364 = qJ(2) + qJ(3);
t361 = sin(t364);
t368 = sin(qJ(1));
t362 = cos(t364);
t416 = qJD(1) * t362;
t392 = -qJD(4) + t416;
t372 = cos(qJ(1));
t421 = t363 * t372;
t433 = t361 * t421 + t392 * t368;
t367 = sin(qJ(2));
t431 = pkin(2) * t367;
t430 = pkin(5) * t361;
t428 = t361 * pkin(3);
t427 = t362 * pkin(5);
t426 = t365 * r_i_i_C(1);
t425 = pkin(2) * qJD(2);
t424 = t361 * t370;
t423 = t362 * t363;
t422 = t363 * t368;
t420 = t368 * t366;
t419 = t368 * t370;
t418 = t372 * t366;
t417 = t372 * t370;
t415 = qJD(1) * t368;
t414 = qJD(1) * t372;
t412 = qJD(4) * t370;
t411 = qJD(5) * t361;
t410 = qJD(5) * t365;
t409 = qJD(5) * t369;
t407 = pkin(5) * t416;
t406 = r_i_i_C(3) * t412;
t405 = t367 * t425;
t403 = t362 * t420;
t402 = t361 * t422;
t395 = t363 * t370 - qJD(5);
t384 = t395 * t365;
t373 = -t361 * t437 + t362 * t384;
t385 = t395 * t369;
t374 = t361 * t434 - t362 * t385;
t380 = t362 * t369 + t365 * t424;
t381 = -t362 * t365 + t369 * t424;
t401 = (t373 * t368 + t380 * t414) * r_i_i_C(2) + (t374 * t368 - t381 * t414) * r_i_i_C(1) + t372 * t407;
t399 = t361 * t415;
t371 = cos(qJ(2));
t396 = -t371 * pkin(2) - pkin(3) * t362 - pkin(1);
t393 = -qJD(4) * t362 + qJD(1);
t391 = (t373 * t372 - t380 * t415) * r_i_i_C(2) + (t374 * t372 + t381 * t415) * r_i_i_C(1) + t436 * t399;
t390 = t436 * t363;
t389 = -t369 * r_i_i_C(1) + t365 * r_i_i_C(2);
t388 = t369 * r_i_i_C(2) + t426;
t382 = t393 * t372;
t344 = t366 * t382 - t370 * t433;
t387 = t372 * t411 + t344;
t350 = t362 * t417 + t420;
t346 = t350 * qJD(1) - qJD(4) * t403 - t370 * t402 - t372 * t412;
t386 = -t368 * t411 - t346;
t383 = -pkin(5) - t388;
t348 = t362 * t419 - t418;
t379 = -qJD(5) * t348 + t361 * t414 + t362 * t422;
t378 = -qJD(5) * t350 + t362 * t421 - t399;
t377 = -t361 * t390 + (t361 * t384 + t362 * t437) * r_i_i_C(2) + (-t361 * t385 - t362 * t434) * r_i_i_C(1) + t362 * t406 + pkin(5) * t423;
t376 = -t361 * t406 + (-t362 * t436 - t430) * t363;
t375 = -t371 * t425 + t376;
t349 = -t362 * t418 + t419;
t347 = -t403 - t417;
t345 = t393 * t419 + (-t392 * t372 + t402) * t366;
t343 = t366 * t433 + t370 * t382;
t336 = t378 * t365 + t387 * t369;
t335 = -t387 * t365 + t378 * t369;
t1 = [(-t346 * t369 + t348 * t410) * r_i_i_C(1) + (t346 * t365 + t348 * t409) * r_i_i_C(2) + t345 * r_i_i_C(3) + (t383 * t361 + t396) * t414 + (t405 + t389 * t411 + (t383 * t362 + t428) * t363) * t368 (-t427 + t431) * t415 + t375 * t372 + t391, -t368 * t407 + t376 * t372 + t391, t344 * r_i_i_C(3) + (-t343 * t365 - t349 * t409) * r_i_i_C(2) + (t343 * t369 - t349 * t410) * r_i_i_C(1), t335 * r_i_i_C(1) - t336 * r_i_i_C(2), 0; t336 * r_i_i_C(1) + t335 * r_i_i_C(2) - t343 * r_i_i_C(3) + (-t405 + (t427 - t428) * t363) * t372 + (t396 - t430) * t415 (-t361 * t436 - t431) * t414 + t375 * t368 + t401, -t368 * t362 * t390 + (-pkin(3) * t414 - pkin(5) * t422 + (-t366 * t414 - t368 * t412) * r_i_i_C(3)) * t361 + t401, t346 * r_i_i_C(3) + (-t345 * t365 - t347 * t409) * r_i_i_C(2) + (t345 * t369 - t347 * t410) * r_i_i_C(1) (t379 * r_i_i_C(1) + t386 * r_i_i_C(2)) * t369 + (t386 * r_i_i_C(1) - t379 * r_i_i_C(2)) * t365, 0; 0, t377 - t405, t377 (r_i_i_C(3) * t370 + t389 * t366) * t423 + (t388 * t366 * qJD(5) + (t389 * t370 - t429) * qJD(4)) * t361 (-r_i_i_C(2) * t385 - t395 * t426) * t362 + (t437 * r_i_i_C(1) + t434 * r_i_i_C(2)) * t361, 0;];
JaD_transl  = t1;
