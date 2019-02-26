% Zeitableitung der translatorischen Teilmatrix der analytischen Jacobi-Matrix für Segment Nr. 6 (0=Basis) von
% S6PRPPRR2
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
%   pkin=[a2,a3,a4,a5,a6,alpha2,d2,d5,d6,theta1,theta3]';
%
% Output:
% JaD_transl [3x6]
%   Zeitableitung der translatorischen Teilmatrix der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-02-26 19:45
% Revision: d75aae1ac561373cd3be920984c3df29a1c2ecc8 (2019-02-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JaD_transl = S6PRPPRR2_jacobiaD_transl_6_sym_varpar(qJ, qJD, r_i_i_C, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(3,1),zeros(11,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRPPRR2_jacobiaD_transl_6_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6PRPPRR2_jacobiaD_transl_6_sym_varpar: qJD has to be [6x1] (double)');
assert(isa(r_i_i_C,'double') && isreal(r_i_i_C) && all(size(r_i_i_C) == [3 1]), ...
	'S6PRPPRR2_jacobiaD_transl_6_sym_varpar: Position vector r_i_i_C has to be [3x1] double');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6PRPPRR2_jacobiaD_transl_6_sym_varpar: pkin has to be [11x1] (double)');

%% Symbolic Calculation
% From jacobiaD_transl_6_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 19:45:14
% EndTime: 2019-02-26 19:45:14
% DurationCPUTime: 0.46s
% Computational Cost: add. (456->79), mult. (1438->146), div. (0->0), fcn. (1546->12), ass. (0->58)
t382 = cos(pkin(6));
t377 = sin(pkin(11));
t380 = cos(pkin(11));
t385 = sin(qJ(2));
t388 = cos(qJ(2));
t400 = t388 * t377 + t385 * t380;
t407 = t400 * t382;
t384 = sin(qJ(5));
t387 = cos(qJ(5));
t383 = sin(qJ(6));
t386 = cos(qJ(6));
t401 = t383 * r_i_i_C(1) + t386 * r_i_i_C(2);
t393 = qJD(6) * t401;
t402 = t386 * r_i_i_C(1) - t383 * r_i_i_C(2);
t399 = pkin(5) + t402;
t421 = pkin(9) + r_i_i_C(3);
t427 = qJD(4) - t384 * t393 + (t421 * t384 + t399 * t387) * qJD(5);
t368 = t385 * t377 - t388 * t380;
t394 = qJD(6) * t402;
t366 = t368 * qJD(2);
t405 = qJD(2) * t388;
t406 = qJD(2) * t385;
t367 = -t377 * t405 - t380 * t406;
t378 = sin(pkin(10));
t381 = cos(pkin(10));
t398 = (t377 * t406 - t380 * t405) * t382;
t345 = t381 * t367 + t378 * t398;
t353 = -t381 * t368 - t378 * t407;
t425 = -t378 * t367 + t381 * t398;
t424 = t378 * t368 - t381 * t407;
t422 = t399 * t384 - t421 * t387 + qJ(4);
t420 = pkin(2) * qJD(2);
t379 = sin(pkin(6));
t417 = t379 * t384;
t416 = t379 * t387;
t413 = t382 * t384;
t412 = t382 * t385;
t404 = qJD(5) * t387;
t403 = t378 * t417;
t364 = t368 * t379;
t355 = t364 * t384 + t382 * t387;
t397 = pkin(3) + pkin(8) + t401;
t392 = t368 * t382;
t349 = -t378 * t400 - t381 * t392;
t396 = t349 * t384 + t381 * t416;
t395 = -t349 * t387 + t381 * t417;
t352 = t378 * t392 - t381 * t400;
t337 = -t352 * t384 + t378 * t416;
t365 = t400 * t379;
t391 = qJD(2) * t407;
t363 = t379 * t366;
t362 = qJD(2) * t365;
t344 = t381 * t366 + t378 * t391;
t341 = t378 * t366 - t381 * t391;
t334 = qJD(5) * t413 - t362 * t384 - t364 * t404;
t332 = t395 * qJD(5) - t341 * t384;
t330 = qJD(5) * t403 + t344 * t384 + t352 * t404;
t1 = [0, t352 * t394 + t397 * t344 + (t378 * t412 - t381 * t388) * t420 + t422 * t345 + t427 * t353, 0, -t344, -t421 * t330 - (-t352 * t387 - t403) * t393 + t399 * (-t337 * qJD(5) - t344 * t387) (t330 * t383 + t345 * t386) * r_i_i_C(1) + (t330 * t386 - t345 * t383) * r_i_i_C(2) + ((-t337 * t386 - t353 * t383) * r_i_i_C(1) + (t337 * t383 - t353 * t386) * r_i_i_C(2)) * qJD(6); 0, t349 * t394 + t397 * t341 + (-t378 * t388 - t381 * t412) * t420 - t422 * t425 - t427 * t424, 0, -t341, t421 * t332 - t395 * t393 + t399 * (t396 * qJD(5) - t341 * t387) (-t332 * t383 - t386 * t425) * r_i_i_C(1) + (-t332 * t386 + t383 * t425) * r_i_i_C(2) + ((t383 * t424 + t386 * t396) * r_i_i_C(1) + (-t383 * t396 + t386 * t424) * r_i_i_C(2)) * qJD(6); 0, -t379 * pkin(2) * t406 - t397 * t362 - t422 * t363 - t364 * t394 + t427 * t365, 0, t362, -t421 * t334 - (t364 * t387 - t413) * t393 + t399 * (-t355 * qJD(5) + t362 * t387) (t334 * t383 - t363 * t386) * r_i_i_C(1) + (t334 * t386 + t363 * t383) * r_i_i_C(2) + ((-t355 * t386 - t365 * t383) * r_i_i_C(1) + (t355 * t383 - t365 * t386) * r_i_i_C(2)) * qJD(6);];
JaD_transl  = t1;
