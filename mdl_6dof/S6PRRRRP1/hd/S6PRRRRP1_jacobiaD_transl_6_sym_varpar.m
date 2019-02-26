% Zeitableitung der translatorischen Teilmatrix der analytischen Jacobi-Matrix für Segment Nr. 6 (0=Basis) von
% S6PRRRRP1
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
%   pkin=[a2,a3,a4,a5,a6,alpha2,d2,d3,d4,d5,theta1]';
%
% Output:
% JaD_transl [3x6]
%   Zeitableitung der translatorischen Teilmatrix der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-02-26 20:15
% Revision: d75aae1ac561373cd3be920984c3df29a1c2ecc8 (2019-02-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JaD_transl = S6PRRRRP1_jacobiaD_transl_6_sym_varpar(qJ, qJD, r_i_i_C, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(3,1),zeros(11,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRRRRP1_jacobiaD_transl_6_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6PRRRRP1_jacobiaD_transl_6_sym_varpar: qJD has to be [6x1] (double)');
assert(isa(r_i_i_C,'double') && isreal(r_i_i_C) && all(size(r_i_i_C) == [3 1]), ...
	'S6PRRRRP1_jacobiaD_transl_6_sym_varpar: Position vector r_i_i_C has to be [3x1] double');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6PRRRRP1_jacobiaD_transl_6_sym_varpar: pkin has to be [11x1] (double)');

%% Symbolic Calculation
% From jacobiaD_transl_6_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 20:15:14
% EndTime: 2019-02-26 20:15:15
% DurationCPUTime: 0.66s
% Computational Cost: add. (776->97), mult. (1336->164), div. (0->0), fcn. (1336->12), ass. (0->69)
t374 = sin(qJ(5));
t377 = cos(qJ(5));
t425 = pkin(5) + r_i_i_C(1);
t431 = -t374 * r_i_i_C(2) + t425 * t377;
t370 = qJ(3) + qJ(4);
t367 = sin(t370);
t368 = cos(t370);
t369 = qJD(3) + qJD(4);
t375 = sin(qJ(3));
t394 = pkin(4) + t431;
t421 = r_i_i_C(3) + qJ(6) + pkin(10);
t392 = t377 * r_i_i_C(2) + t425 * t374;
t427 = t392 * qJD(5);
t381 = (-t421 * t369 + t427) * t368 + qJD(3) * t375 * pkin(3) + (t394 * t369 - qJD(6)) * t367;
t376 = sin(qJ(2));
t379 = cos(qJ(2));
t371 = sin(pkin(11));
t420 = cos(pkin(6));
t405 = t371 * t420;
t419 = cos(pkin(11));
t359 = -t376 * t405 + t419 * t379;
t418 = t367 * t369;
t417 = t368 * t369;
t372 = sin(pkin(6));
t416 = t371 * t372;
t415 = t372 * t376;
t414 = t372 * t379;
t413 = qJD(2) * t376;
t412 = qJD(5) * t374;
t411 = qJD(5) * t377;
t409 = t367 * t415;
t408 = t368 * t415;
t407 = t372 * t413;
t406 = qJD(2) * t414;
t404 = t372 * t419;
t401 = t367 * t404;
t358 = t419 * t376 + t379 * t405;
t354 = t358 * qJD(2);
t400 = t369 * t416 - t354;
t399 = t420 * t419;
t393 = t379 * t399;
t352 = -qJD(2) * t393 + t371 * t413;
t357 = t371 * t379 + t376 * t399;
t337 = -t357 * t418 + (-t369 * t404 - t352) * t368;
t353 = t357 * qJD(2);
t398 = -t337 * t374 + t353 * t377;
t339 = -t359 * t418 + t400 * t368;
t355 = t359 * qJD(2);
t397 = -t339 * t374 + t355 * t377;
t347 = t357 * t368 - t401;
t356 = t371 * t376 - t393;
t396 = -t347 * t377 - t356 * t374;
t349 = t359 * t368 + t367 * t416;
t395 = -t349 * t377 - t358 * t374;
t351 = t420 * t367 + t408;
t391 = -t351 * t377 + t374 * t414;
t386 = t420 * t369 + t406;
t345 = t386 * t368 - t369 * t409;
t388 = -t345 * t374 + t377 * t407;
t378 = cos(qJ(3));
t385 = -t378 * pkin(3) - t421 * t367 - t394 * t368 - pkin(2);
t336 = -t352 * t367 + t357 * t417 - t369 * t401;
t384 = t347 * qJD(6) - (-t357 * t367 - t368 * t404) * t427 + t421 * t337 - t394 * t336;
t338 = t359 * t417 + t400 * t367;
t383 = t349 * qJD(6) - (-t359 * t367 + t368 * t416) * t427 + t421 * t339 - t394 * t338;
t344 = t386 * t367 + t369 * t408;
t382 = t351 * qJD(6) - (t420 * t368 - t409) * t427 + t421 * t345 - t394 * t344;
t380 = -pkin(9) - pkin(8);
t1 = [0 (-t354 * t377 - t359 * t412) * r_i_i_C(2) + t354 * t380 + t385 * t355 + t381 * t358 + t425 * (-t354 * t374 + t359 * t411) (t354 * t375 + (-t359 * t378 - t375 * t416) * qJD(3)) * pkin(3) + t383, t383, t397 * r_i_i_C(1) + (-t339 * t377 - t355 * t374) * r_i_i_C(2) + (t395 * r_i_i_C(1) + (t349 * t374 - t358 * t377) * r_i_i_C(2)) * qJD(5) + (t395 * qJD(5) + t397) * pkin(5), t338; 0 (-t352 * t377 - t357 * t412) * r_i_i_C(2) + t352 * t380 + t385 * t353 + t381 * t356 + t425 * (-t352 * t374 + t357 * t411) (t352 * t375 + (-t357 * t378 + t375 * t404) * qJD(3)) * pkin(3) + t384, t384, t398 * r_i_i_C(1) + (-t337 * t377 - t353 * t374) * r_i_i_C(2) + (t396 * r_i_i_C(1) + (t347 * t374 - t356 * t377) * r_i_i_C(2)) * qJD(5) + (t396 * qJD(5) + t398) * pkin(5), t336; 0 ((t385 * qJD(2) + t431 * qJD(5)) * t376 + ((-t380 + t392) * qJD(2) - t381) * t379) * t372 (-t375 * t406 + (-t420 * t375 - t378 * t415) * qJD(3)) * pkin(3) + t382, t382, t388 * r_i_i_C(1) + (-t345 * t377 - t374 * t407) * r_i_i_C(2) + (t391 * r_i_i_C(1) + (t351 * t374 + t377 * t414) * r_i_i_C(2)) * qJD(5) + (t391 * qJD(5) + t388) * pkin(5), t344;];
JaD_transl  = t1;
