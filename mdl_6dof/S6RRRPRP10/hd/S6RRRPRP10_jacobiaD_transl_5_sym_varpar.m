% Zeitableitung der translatorischen Teilmatrix der analytischen Jacobi-Matrix für Segment Nr. 5 (0=Basis) von
% S6RRRPRP10
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
%   pkin=[a2,a3,a4,a5,a6,alpha2,d1,d2,d3,d5,theta4]';
%
% Output:
% JaD_transl [3x6]
%   Zeitableitung der translatorischen Teilmatrix der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-02-26 22:14
% Revision: d75aae1ac561373cd3be920984c3df29a1c2ecc8 (2019-02-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JaD_transl = S6RRRPRP10_jacobiaD_transl_5_sym_varpar(qJ, qJD, r_i_i_C, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(3,1),zeros(11,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRPRP10_jacobiaD_transl_5_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRRPRP10_jacobiaD_transl_5_sym_varpar: qJD has to be [6x1] (double)');
assert(isa(r_i_i_C,'double') && isreal(r_i_i_C) && all(size(r_i_i_C) == [3 1]), ...
	'S6RRRPRP10_jacobiaD_transl_5_sym_varpar: Position vector r_i_i_C has to be [3x1] double');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6RRRPRP10_jacobiaD_transl_5_sym_varpar: pkin has to be [11x1] (double)');

%% Symbolic Calculation
% From jacobiaD_transl_5_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 22:14:22
% EndTime: 2019-02-26 22:14:22
% DurationCPUTime: 0.47s
% Computational Cost: add. (600->95), mult. (1491->158), div. (0->0), fcn. (1502->12), ass. (0->66)
t386 = sin(qJ(1));
t382 = cos(pkin(6));
t401 = qJD(2) * t382 + qJD(1);
t385 = sin(qJ(2));
t419 = t386 * t385;
t410 = t382 * t419;
t414 = qJD(2) * t385;
t388 = cos(qJ(2));
t389 = cos(qJ(1));
t416 = t389 * t388;
t352 = -qJD(1) * t410 - t386 * t414 + t401 * t416;
t381 = sin(pkin(6));
t420 = t381 * t389;
t433 = -qJD(3) * t420 + t352;
t417 = t389 * t385;
t418 = t386 * t388;
t363 = t382 * t417 + t418;
t415 = qJD(1) * t381;
t432 = -qJD(3) * t363 + t386 * t415;
t384 = sin(qJ(3));
t387 = cos(qJ(3));
t356 = t363 * t387 - t384 * t420;
t409 = t382 * t416;
t362 = -t409 + t419;
t379 = pkin(11) + qJ(5);
t377 = sin(t379);
t378 = cos(t379);
t431 = t356 * t377 - t362 * t378;
t430 = t356 * t378 + t362 * t377;
t399 = t377 * r_i_i_C(1) + t378 * r_i_i_C(2);
t395 = qJD(5) * t399;
t376 = cos(pkin(11)) * pkin(4) + pkin(3);
t400 = t378 * r_i_i_C(1) - t377 * r_i_i_C(2);
t398 = t376 + t400;
t426 = r_i_i_C(3) + pkin(10) + qJ(4);
t390 = (t398 * t384 - t426 * t387) * qJD(3) - t384 * qJD(4) + t387 * t395;
t346 = t432 * t384 + t433 * t387;
t427 = t426 * t384 + t398 * t387 + pkin(2);
t423 = t381 * t386;
t422 = t381 * t387;
t421 = t381 * t388;
t413 = qJD(2) * t388;
t408 = -sin(pkin(11)) * pkin(4) - pkin(9);
t406 = t389 * t415;
t405 = t381 * t414;
t404 = t381 * t413;
t397 = t363 * t384 + t387 * t420;
t365 = -t410 + t416;
t358 = -t365 * t384 + t386 * t422;
t359 = t365 * t387 + t384 * t423;
t361 = t382 * t384 + t385 * t422;
t396 = -t381 * t385 * t384 + t382 * t387;
t364 = t382 * t418 + t417;
t394 = t400 * qJD(5);
t392 = t399 - t408;
t345 = t433 * t384 - t432 * t387;
t354 = t396 * qJD(3) + t387 * t404;
t353 = t361 * qJD(3) + t384 * t404;
t351 = t364 * qJD(1) + t363 * qJD(2);
t350 = t363 * qJD(1) + t364 * qJD(2);
t349 = -qJD(1) * t409 - t389 * t413 + t401 * t419;
t344 = t358 * qJD(3) - t350 * t387 + t384 * t406;
t343 = t359 * qJD(3) - t350 * t384 - t387 * t406;
t342 = t344 * t378 - t349 * t377 + (-t359 * t377 + t364 * t378) * qJD(5);
t341 = -t344 * t377 - t349 * t378 + (-t359 * t378 - t364 * t377) * qJD(5);
t1 = [-t397 * qJD(4) - t352 * pkin(2) - t398 * t346 - t392 * t351 - t426 * t345 + (t431 * r_i_i_C(1) + t430 * r_i_i_C(2)) * qJD(5) + (-t389 * pkin(1) - pkin(8) * t423) * qJD(1), t427 * t349 - t392 * t350 + t390 * t364 + t365 * t394, t359 * qJD(4) - t398 * t343 + t426 * t344 - t358 * t395, t343, t341 * r_i_i_C(1) - t342 * r_i_i_C(2), 0; -t350 * pkin(2) + t342 * r_i_i_C(1) + t341 * r_i_i_C(2) - t358 * qJD(4) + t344 * t376 + t408 * t349 + t426 * t343 + (-pkin(1) * t386 + pkin(8) * t420) * qJD(1), -t351 * t427 + t392 * t352 + t390 * t362 + t363 * t394, t356 * qJD(4) - t398 * t345 + t426 * t346 + t397 * t395, t345 (-t346 * t377 + t351 * t378) * r_i_i_C(1) + (-t346 * t378 - t351 * t377) * r_i_i_C(2) + (-t430 * r_i_i_C(1) + t431 * r_i_i_C(2)) * qJD(5), 0; 0 ((-qJD(2) * t427 + t394) * t385 + (t392 * qJD(2) - t390) * t388) * t381, t361 * qJD(4) - t398 * t353 + t426 * t354 - t396 * t395, t353 (-t354 * t377 + t378 * t405) * r_i_i_C(1) + (-t354 * t378 - t377 * t405) * r_i_i_C(2) + ((-t361 * t378 + t377 * t421) * r_i_i_C(1) + (t361 * t377 + t378 * t421) * r_i_i_C(2)) * qJD(5), 0;];
JaD_transl  = t1;
