% Zeitableitung der translatorischen Teilmatrix der analytischen Jacobi-Matrix für Segment Nr. 5 (0=Basis) von
% S6RRRPPR9
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
%   pkin=[a2,a3,a4,a5,a6,alpha2,d1,d2,d3,d6,theta4]';
%
% Output:
% JaD_transl [3x6]
%   Zeitableitung der translatorischen Teilmatrix der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-02-26 22:08
% Revision: d75aae1ac561373cd3be920984c3df29a1c2ecc8 (2019-02-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JaD_transl = S6RRRPPR9_jacobiaD_transl_5_sym_varpar(qJ, qJD, r_i_i_C, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(3,1),zeros(11,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRPPR9_jacobiaD_transl_5_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRRPPR9_jacobiaD_transl_5_sym_varpar: qJD has to be [6x1] (double)');
assert(isa(r_i_i_C,'double') && isreal(r_i_i_C) && all(size(r_i_i_C) == [3 1]), ...
	'S6RRRPPR9_jacobiaD_transl_5_sym_varpar: Position vector r_i_i_C has to be [3x1] double');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6RRRPPR9_jacobiaD_transl_5_sym_varpar: pkin has to be [11x1] (double)');

%% Symbolic Calculation
% From jacobiaD_transl_5_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 22:08:12
% EndTime: 2019-02-26 22:08:12
% DurationCPUTime: 0.53s
% Computational Cost: add. (529->96), mult. (1587->159), div. (0->0), fcn. (1585->10), ass. (0->62)
t370 = sin(qJ(3));
t373 = cos(qJ(3));
t369 = cos(pkin(6));
t371 = sin(qJ(2));
t375 = cos(qJ(1));
t405 = t375 * t371;
t372 = sin(qJ(1));
t374 = cos(qJ(2));
t406 = t372 * t374;
t354 = t369 * t405 + t406;
t367 = sin(pkin(6));
t403 = qJD(1) * t367;
t418 = -qJD(3) * t354 + t372 * t403;
t386 = qJD(2) * t369 + qJD(1);
t407 = t372 * t371;
t396 = t369 * t407;
t402 = qJD(2) * t371;
t404 = t375 * t374;
t348 = -qJD(1) * t396 - t372 * t402 + t386 * t404;
t409 = t367 * t375;
t419 = -qJD(3) * t409 + t348;
t340 = t418 * t370 + t419 * t373;
t355 = t369 * t406 + t405;
t347 = t355 * qJD(1) + t354 * qJD(2);
t366 = sin(pkin(11));
t368 = cos(pkin(11));
t420 = t340 * t366 - t347 * t368;
t398 = qJD(5) * t366;
t413 = r_i_i_C(2) + qJ(4);
t376 = (-pkin(3) * t370 + t413 * t373) * qJD(3) + t370 * qJD(4) + t373 * t398;
t417 = t373 * pkin(3) + t413 * t370 + pkin(2);
t415 = pkin(4) + r_i_i_C(1);
t412 = r_i_i_C(3) + qJ(5);
t410 = t367 * t372;
t408 = t371 * t373;
t401 = qJD(2) * t374;
t399 = qJD(3) * t370;
t397 = t368 * qJD(5);
t395 = t369 * t404;
t392 = t375 * t403;
t391 = t367 * t401;
t389 = t374 * t399;
t387 = t354 * t373 - t370 * t409;
t385 = t354 * t370 + t373 * t409;
t356 = -t396 + t404;
t351 = -t356 * t370 + t373 * t410;
t352 = t356 * t373 + t370 * t410;
t384 = t367 * t408 + t369 * t370;
t383 = -t367 * t371 * t370 + t369 * t373;
t345 = -qJD(1) * t395 - t375 * t401 + t386 * t407;
t382 = t345 * t373 + t355 * t399;
t353 = t395 - t407;
t381 = -t347 * t373 - t353 * t399;
t377 = -t412 * t366 - t415 * t368 - pkin(3);
t339 = t419 * t370 - t418 * t373;
t350 = t383 * qJD(3) + t373 * t391;
t349 = t384 * qJD(3) + t370 * t391;
t346 = t354 * qJD(1) + t355 * qJD(2);
t338 = t351 * qJD(3) - t346 * t373 + t370 * t392;
t337 = t352 * qJD(3) - t346 * t370 - t373 * t392;
t329 = t338 * t366 + t345 * t368;
t1 = [-(t353 * t368 + t387 * t366) * qJD(5) - t340 * pkin(3) - t385 * qJD(4) - t348 * pkin(2) - t347 * pkin(9) - t413 * t339 + t415 * (-t340 * t368 - t347 * t366) - t412 * t420 + (-t375 * pkin(1) - pkin(8) * t410) * qJD(1), -t356 * t397 - t346 * pkin(9) + t415 * (-t346 * t366 + t382 * t368) + t412 * (t346 * t368 + t382 * t366) - t376 * t355 + t417 * t345, t352 * qJD(4) + t377 * t337 + t413 * t338 + t351 * t398, t337, t329, 0; -(-t352 * t366 + t355 * t368) * qJD(5) + t338 * pkin(3) - t351 * qJD(4) - t346 * pkin(2) - t345 * pkin(9) + t413 * t337 + t415 * (t338 * t368 - t345 * t366) + t412 * t329 + (-t372 * pkin(1) + pkin(8) * t409) * qJD(1), -t354 * t397 + t348 * pkin(9) + t415 * (t348 * t366 + t381 * t368) + t412 * (-t348 * t368 + t381 * t366) + t376 * t353 - t417 * t347, t387 * qJD(4) + t377 * t339 + t413 * t340 - t385 * t398, t339, t420, 0; 0 (t415 * (-t368 * t389 + (t366 * t374 - t368 * t408) * qJD(2)) - t412 * (t366 * t389 + (t366 * t408 + t368 * t374) * qJD(2)) - t371 * t397 + t376 * t374 + (t374 * pkin(9) - t371 * t417) * qJD(2)) * t367, t384 * qJD(4) + t377 * t349 + t413 * t350 + t383 * t398, t349, -t367 * t368 * t402 + t350 * t366, 0;];
JaD_transl  = t1;
