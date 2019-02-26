% Zeitableitung der translatorischen Teilmatrix der analytischen Jacobi-Matrix für Segment Nr. 6 (0=Basis) von
% S6PRRPPR4
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
%   pkin=[a2,a3,a4,a5,a6,alpha2,d2,d3,d6,theta1,theta4]';
%
% Output:
% JaD_transl [3x6]
%   Zeitableitung der translatorischen Teilmatrix der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-02-26 20:00
% Revision: d75aae1ac561373cd3be920984c3df29a1c2ecc8 (2019-02-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JaD_transl = S6PRRPPR4_jacobiaD_transl_6_sym_varpar(qJ, qJD, r_i_i_C, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(3,1),zeros(11,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRRPPR4_jacobiaD_transl_6_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6PRRPPR4_jacobiaD_transl_6_sym_varpar: qJD has to be [6x1] (double)');
assert(isa(r_i_i_C,'double') && isreal(r_i_i_C) && all(size(r_i_i_C) == [3 1]), ...
	'S6PRRPPR4_jacobiaD_transl_6_sym_varpar: Position vector r_i_i_C has to be [3x1] double');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6PRRPPR4_jacobiaD_transl_6_sym_varpar: pkin has to be [11x1] (double)');

%% Symbolic Calculation
% From jacobiaD_transl_6_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 20:00:01
% EndTime: 2019-02-26 20:00:02
% DurationCPUTime: 0.72s
% Computational Cost: add. (541->121), mult. (1742->225), div. (0->0), fcn. (1827->12), ass. (0->77)
t373 = sin(qJ(3));
t376 = cos(qJ(3));
t397 = -r_i_i_C(3) - pkin(9) + qJ(4);
t415 = t373 * qJD(4) + (-pkin(3) * t373 + t397 * t376) * qJD(3);
t369 = sin(pkin(10));
t377 = cos(qJ(2));
t410 = cos(pkin(10));
t411 = cos(pkin(6));
t387 = t411 * t410;
t385 = t377 * t387;
t374 = sin(qJ(2));
t400 = qJD(2) * t374;
t352 = -qJD(2) * t385 + t369 * t400;
t370 = sin(pkin(6));
t391 = t370 * t410;
t414 = -qJD(3) * t391 - t352;
t392 = t369 * t411;
t359 = -t374 * t392 + t410 * t377;
t357 = t369 * t377 + t374 * t387;
t409 = t357 * t376;
t368 = sin(pkin(11));
t408 = t368 * t376;
t407 = t369 * t370;
t406 = t370 * t373;
t405 = t370 * t377;
t371 = cos(pkin(11));
t404 = t371 * t376;
t403 = t371 * t377;
t402 = t374 * t376;
t401 = t376 * t377;
t399 = qJD(3) * t373;
t396 = t374 * t406;
t395 = t370 * t400;
t394 = qJD(2) * t405;
t393 = t377 * t399;
t372 = sin(qJ(6));
t375 = cos(qJ(6));
t386 = t372 * r_i_i_C(1) + t375 * r_i_i_C(2) + qJ(5);
t384 = t375 * r_i_i_C(1) - t372 * r_i_i_C(2) + pkin(4) + pkin(5);
t345 = t359 * t376 + t369 * t406;
t353 = t357 * qJD(2);
t356 = t369 * t374 - t385;
t383 = -t353 * t376 + t356 * t399;
t355 = t359 * qJD(2);
t358 = t410 * t374 + t377 * t392;
t382 = -t355 * t376 + t358 * t399;
t361 = t370 * t402 + t411 * t373;
t381 = -t376 * pkin(3) - t397 * t373 - pkin(2);
t379 = -t386 * t368 - t384 * t371 - pkin(3);
t378 = t368 * qJD(5) + ((t368 * t375 - t371 * t372) * r_i_i_C(1) + (-t368 * t372 - t371 * t375) * r_i_i_C(2)) * qJD(6);
t354 = t358 * qJD(2);
t349 = (t368 * t374 + t371 * t401) * t370;
t348 = (t368 * t401 - t371 * t374) * t370;
t347 = -qJD(3) * t396 + (t411 * qJD(3) + t394) * t376;
t346 = t361 * qJD(3) + t373 * t394;
t343 = -t373 * t391 + t409;
t341 = t361 * t371 - t368 * t405;
t340 = t361 * t368 + t370 * t403;
t337 = -t358 * t404 + t359 * t368;
t336 = -t358 * t408 - t359 * t371;
t335 = -t356 * t404 + t357 * t368;
t334 = -t356 * t408 - t357 * t371;
t333 = t347 * t371 + t368 * t395;
t332 = t347 * t368 - t371 * t395;
t331 = t345 * t371 + t358 * t368;
t330 = t345 * t368 - t358 * t371;
t329 = t343 * t371 + t356 * t368;
t328 = t343 * t368 - t356 * t371;
t327 = -t359 * t399 + (qJD(3) * t407 - t354) * t376;
t326 = t345 * qJD(3) - t354 * t373;
t325 = -t357 * t399 + t414 * t376;
t324 = qJD(3) * t409 + t414 * t373;
t319 = t327 * t371 + t355 * t368;
t318 = t327 * t368 - t355 * t371;
t317 = t325 * t371 + t353 * t368;
t316 = t325 * t368 - t353 * t371;
t1 = [0, -t354 * pkin(8) + t336 * qJD(5) + t386 * (t354 * t371 + t382 * t368) + t384 * (-t354 * t368 + t382 * t371) + ((t336 * t375 - t337 * t372) * r_i_i_C(1) + (-t336 * t372 - t337 * t375) * r_i_i_C(2)) * qJD(6) - t415 * t358 + t381 * t355, t345 * qJD(4) + t397 * t327 + t378 * (-t359 * t373 + t376 * t407) + t379 * t326, t326, t318 (t318 * t375 - t319 * t372) * r_i_i_C(1) + (-t318 * t372 - t319 * t375) * r_i_i_C(2) + ((-t330 * t372 - t331 * t375) * r_i_i_C(1) + (-t330 * t375 + t331 * t372) * r_i_i_C(2)) * qJD(6); 0, -t352 * pkin(8) + t334 * qJD(5) + t386 * (t352 * t371 + t383 * t368) + t384 * (-t352 * t368 + t383 * t371) + ((t334 * t375 - t335 * t372) * r_i_i_C(1) + (-t334 * t372 - t335 * t375) * r_i_i_C(2)) * qJD(6) - t415 * t356 + t381 * t353, t343 * qJD(4) + t397 * t325 + t378 * (-t357 * t373 - t376 * t391) + t379 * t324, t324, t316 (t316 * t375 - t317 * t372) * r_i_i_C(1) + (-t316 * t372 - t317 * t375) * r_i_i_C(2) + ((-t328 * t372 - t329 * t375) * r_i_i_C(1) + (-t328 * t375 + t329 * t372) * r_i_i_C(2)) * qJD(6); 0, t348 * qJD(5) + ((t348 * t375 - t349 * t372) * r_i_i_C(1) + (-t348 * t372 - t349 * t375) * r_i_i_C(2)) * qJD(6) + (-t386 * (t368 * t393 + (t368 * t402 + t403) * qJD(2)) + t384 * (-t371 * t393 + (t368 * t377 - t371 * t402) * qJD(2)) + t415 * t377 + (pkin(8) * t377 + t381 * t374) * qJD(2)) * t370, t361 * qJD(4) + t397 * t347 + t378 * (t411 * t376 - t396) + t379 * t346, t346, t332 (t332 * t375 - t333 * t372) * r_i_i_C(1) + (-t332 * t372 - t333 * t375) * r_i_i_C(2) + ((-t340 * t372 - t341 * t375) * r_i_i_C(1) + (-t340 * t375 + t341 * t372) * r_i_i_C(2)) * qJD(6);];
JaD_transl  = t1;
