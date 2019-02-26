% Zeitableitung der translatorischen Teilmatrix der analytischen Jacobi-Matrix für Segment Nr. 6 (0=Basis) von
% S6RRPPRR7
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
% Datum: 2019-02-26 21:32
% Revision: d75aae1ac561373cd3be920984c3df29a1c2ecc8 (2019-02-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JaD_transl = S6RRPPRR7_jacobiaD_transl_6_sym_varpar(qJ, qJD, r_i_i_C, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(3,1),zeros(10,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPPRR7_jacobiaD_transl_6_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRPPRR7_jacobiaD_transl_6_sym_varpar: qJD has to be [6x1] (double)');
assert(isa(r_i_i_C,'double') && isreal(r_i_i_C) && all(size(r_i_i_C) == [3 1]), ...
	'S6RRPPRR7_jacobiaD_transl_6_sym_varpar: Position vector r_i_i_C has to be [3x1] double');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RRPPRR7_jacobiaD_transl_6_sym_varpar: pkin has to be [10x1] (double)');

%% Symbolic Calculation
% From jacobiaD_transl_6_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 21:31:58
% EndTime: 2019-02-26 21:31:59
% DurationCPUTime: 0.49s
% Computational Cost: add. (516->98), mult. (1526->155), div. (0->0), fcn. (1516->10), ass. (0->66)
t363 = cos(pkin(6));
t370 = cos(qJ(2));
t371 = cos(qJ(1));
t399 = t370 * t371;
t391 = t363 * t399;
t366 = sin(qJ(2));
t367 = sin(qJ(1));
t402 = t366 * t367;
t351 = -t391 + t402;
t365 = sin(qJ(5));
t369 = cos(qJ(5));
t362 = sin(pkin(6));
t403 = t362 * t371;
t345 = t351 * t369 + t365 * t403;
t400 = t367 * t370;
t401 = t366 * t371;
t352 = t363 * t401 + t400;
t364 = sin(qJ(6));
t368 = cos(qJ(6));
t413 = t345 * t364 - t352 * t368;
t412 = t345 * t368 + t352 * t364;
t384 = r_i_i_C(1) * t368 - r_i_i_C(2) * t364;
t377 = t384 * qJD(6);
t353 = t363 * t400 + t401;
t340 = t353 * qJD(1) + t352 * qJD(2);
t381 = -t351 * t365 + t369 * t403;
t398 = qJD(1) * t362;
t390 = t367 * t398;
t335 = t381 * qJD(5) + t340 * t369 - t365 * t390;
t382 = pkin(5) + t384;
t410 = -pkin(4) - qJ(3);
t411 = pkin(10) + r_i_i_C(3);
t374 = t411 * t365 + t382 * t369 - t410;
t409 = pkin(8) - qJ(4);
t406 = t362 * t366;
t405 = t362 * t367;
t404 = t362 * t370;
t397 = qJD(2) * t366;
t396 = qJD(2) * t370;
t395 = qJD(4) * t362;
t394 = -pkin(2) - pkin(3) - pkin(9);
t393 = t365 * t405;
t392 = t363 * t402;
t389 = t371 * t398;
t388 = t362 * t396;
t387 = t362 * t397;
t385 = qJD(2) * t363 + qJD(1);
t383 = -r_i_i_C(1) * t364 - r_i_i_C(2) * t368;
t380 = -t353 * t365 - t369 * t405;
t379 = t363 * t365 + t369 * t404;
t378 = -t363 * t369 + t365 * t404;
t376 = qJD(6) * t383;
t375 = -t383 - t394;
t373 = -t345 * qJD(5) - t340 * t365 - t369 * t390;
t372 = qJD(3) + t369 * t376 + (-t382 * t365 + t411 * t369) * qJD(5);
t354 = -t392 + t399;
t348 = t353 * t369 - t393;
t343 = t378 * qJD(5) + t369 * t387;
t341 = -qJD(1) * t392 - t367 * t397 + t385 * t399;
t339 = t352 * qJD(1) + t353 * qJD(2);
t338 = -qJD(1) * t391 - t371 * t396 + t385 * t402;
t333 = t380 * qJD(5) - t338 * t369 - t365 * t389;
t332 = -t338 * t365 - qJD(5) * t393 + (qJD(5) * t353 + t389) * t369;
t331 = t333 * t368 - t339 * t364 + (-t348 * t364 + t354 * t368) * qJD(6);
t330 = -t333 * t364 - t339 * t368 + (-t348 * t368 - t354 * t364) * qJD(6);
t1 = [-t371 * t395 - t351 * qJD(3) + t410 * t340 + t411 * t373 - t382 * t335 + (t413 * r_i_i_C(1) + t412 * r_i_i_C(2)) * qJD(6) - t375 * t341 + (-t371 * pkin(1) - t409 * t405) * qJD(1), t375 * t338 - t374 * t339 - t353 * t377 + t372 * t354, -t338, -t389, -t382 * t332 + t411 * t333 + t380 * t376, r_i_i_C(1) * t330 - t331 * r_i_i_C(2); -t367 * t395 + t333 * pkin(5) + t331 * r_i_i_C(1) + t330 * r_i_i_C(2) + t353 * qJD(3) + t410 * t338 + t411 * t332 + t394 * t339 + (-pkin(1) * t367 + t409 * t403) * qJD(1), -t375 * t340 + t374 * t341 - t351 * t377 + t372 * t352, t340, -t390, t411 * t335 + t382 * t373 + t381 * t376 (-t335 * t364 + t341 * t368) * r_i_i_C(1) + (-t335 * t368 - t341 * t364) * r_i_i_C(2) + (-t412 * r_i_i_C(1) + t413 * r_i_i_C(2)) * qJD(6); 0 ((t374 * qJD(2) + t377) * t370 + (-t375 * qJD(2) + t372) * t366) * t362, t387, 0, t411 * t343 + t378 * t376 + t382 * (t379 * qJD(5) - t365 * t387) (-t343 * t364 + t368 * t388) * r_i_i_C(1) + (-t343 * t368 - t364 * t388) * r_i_i_C(2) + ((-t364 * t406 + t368 * t379) * r_i_i_C(1) + (-t364 * t379 - t368 * t406) * r_i_i_C(2)) * qJD(6);];
JaD_transl  = t1;
