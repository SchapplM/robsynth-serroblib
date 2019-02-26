% Zeitableitung der translatorischen Teilmatrix der analytischen Jacobi-Matrix für Segment Nr. 6 (0=Basis) von
% S6RRPPRR5
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
%   pkin=[a2,a3,a4,a5,a6,alpha2,d1,d2,d5,d6]';
%
% Output:
% JaD_transl [3x6]
%   Zeitableitung der translatorischen Teilmatrix der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-02-26 21:30
% Revision: d75aae1ac561373cd3be920984c3df29a1c2ecc8 (2019-02-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JaD_transl = S6RRPPRR5_jacobiaD_transl_6_sym_varpar(qJ, qJD, r_i_i_C, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(3,1),zeros(10,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPPRR5_jacobiaD_transl_6_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRPPRR5_jacobiaD_transl_6_sym_varpar: qJD has to be [6x1] (double)');
assert(isa(r_i_i_C,'double') && isreal(r_i_i_C) && all(size(r_i_i_C) == [3 1]), ...
	'S6RRPPRR5_jacobiaD_transl_6_sym_varpar: Position vector r_i_i_C has to be [3x1] double');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RRPPRR5_jacobiaD_transl_6_sym_varpar: pkin has to be [10x1] (double)');

%% Symbolic Calculation
% From jacobiaD_transl_6_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 21:30:50
% EndTime: 2019-02-26 21:30:50
% DurationCPUTime: 0.54s
% Computational Cost: add. (516->98), mult. (1526->159), div. (0->0), fcn. (1516->10), ass. (0->66)
t365 = sin(qJ(2));
t366 = sin(qJ(1));
t369 = cos(qJ(2));
t370 = cos(qJ(1));
t408 = cos(pkin(6));
t388 = t370 * t408;
t351 = t365 * t388 + t366 * t369;
t364 = sin(qJ(5));
t368 = cos(qJ(5));
t362 = sin(pkin(6));
t402 = t362 * t370;
t344 = t351 * t368 + t364 * t402;
t385 = t369 * t388;
t401 = t366 * t365;
t350 = -t385 + t401;
t363 = sin(qJ(6));
t367 = cos(qJ(6));
t417 = t344 * t363 + t350 * t367;
t416 = t344 * t367 - t350 * t363;
t383 = t363 * r_i_i_C(1) + t367 * r_i_i_C(2);
t378 = qJD(6) * t383;
t384 = -t367 * r_i_i_C(1) + t363 * r_i_i_C(2);
t381 = pkin(5) - t384;
t411 = pkin(10) + r_i_i_C(3);
t415 = (t381 * t364 - t411 * t368) * qJD(5) + t368 * t378;
t382 = qJD(2) * t408 + qJD(1);
t389 = t366 * t408;
t386 = t365 * t389;
t398 = qJD(2) * t365;
t400 = t370 * t369;
t340 = -qJD(1) * t386 - t366 * t398 + t382 * t400;
t380 = -t351 * t364 + t368 * t402;
t399 = qJD(1) * t362;
t393 = t366 * t399;
t334 = t380 * qJD(5) + t340 * t368 - t364 * t393;
t395 = -pkin(2) - pkin(3) - pkin(4);
t412 = t411 * t364 + t381 * t368 - t395;
t410 = pkin(8) - qJ(4);
t409 = pkin(9) - qJ(3);
t405 = t362 * t366;
t404 = t362 * t368;
t403 = t362 * t369;
t397 = qJD(2) * t369;
t396 = qJD(4) * t362;
t394 = t364 * t405;
t392 = t370 * t399;
t391 = t362 * t397;
t390 = t362 * t398;
t353 = -t386 + t400;
t379 = -t353 * t364 - t366 * t404;
t377 = t383 + t409;
t376 = -t362 * t365 * t364 - t408 * t368;
t349 = -t408 * t364 + t365 * t404;
t375 = t370 * t365 + t369 * t389;
t374 = t384 * qJD(6) + qJD(3);
t372 = -t344 * qJD(5) - t340 * t364 - t368 * t393;
t347 = t353 * t368 - t394;
t342 = t376 * qJD(5) + t368 * t391;
t339 = t375 * qJD(1) + t351 * qJD(2);
t338 = t351 * qJD(1) + t375 * qJD(2);
t337 = -qJD(1) * t385 - t370 * t397 + t382 * t401;
t332 = t379 * qJD(5) - t338 * t368 - t364 * t392;
t331 = -t338 * t364 - qJD(5) * t394 + (qJD(5) * t353 + t392) * t368;
t330 = t332 * t367 + t337 * t363 + (-t347 * t363 - t367 * t375) * qJD(6);
t329 = -t332 * t363 + t337 * t367 + (-t347 * t367 + t363 * t375) * qJD(6);
t1 = [-t370 * t396 - t350 * qJD(3) - t381 * t334 + t377 * t339 + t411 * t372 + (t417 * r_i_i_C(1) + t416 * r_i_i_C(2)) * qJD(6) + t395 * t340 + (-t370 * pkin(1) - t410 * t405) * qJD(1), t412 * t337 + t377 * t338 + t374 * t353 + t375 * t415, -t337, -t392, -t381 * t331 + t411 * t332 - t379 * t378, t329 * r_i_i_C(1) - t330 * r_i_i_C(2); -t366 * t396 + t332 * pkin(5) + t330 * r_i_i_C(1) + t329 * r_i_i_C(2) + t375 * qJD(3) + t409 * t337 + t411 * t331 + t395 * t338 + (-pkin(1) * t366 + t410 * t402) * qJD(1), -t339 * t412 - t377 * t340 + t415 * t350 + t374 * t351, t339, -t393, t411 * t334 + t381 * t372 - t380 * t378 (-t334 * t363 - t339 * t367) * r_i_i_C(1) + (-t334 * t367 + t339 * t363) * r_i_i_C(2) + (-t416 * r_i_i_C(1) + t417 * r_i_i_C(2)) * qJD(6); 0 (t374 * t365 - t415 * t369 + (-t365 * t412 - t377 * t369) * qJD(2)) * t362, t390, 0, t411 * t342 - t376 * t378 + t381 * (-t349 * qJD(5) - t364 * t391) (-t342 * t363 - t367 * t390) * r_i_i_C(1) + (-t342 * t367 + t363 * t390) * r_i_i_C(2) + ((-t349 * t367 - t363 * t403) * r_i_i_C(1) + (t349 * t363 - t367 * t403) * r_i_i_C(2)) * qJD(6);];
JaD_transl  = t1;
