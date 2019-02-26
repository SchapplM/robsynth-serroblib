% Zeitableitung der translatorischen Teilmatrix der analytischen Jacobi-Matrix für Segment Nr. 6 (0=Basis) von
% S6PRRPRP2
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
%   pkin=[a2,a3,a4,a5,a6,alpha2,d2,d3,d5,theta1,theta4]';
%
% Output:
% JaD_transl [3x6]
%   Zeitableitung der translatorischen Teilmatrix der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-02-26 20:01
% Revision: d75aae1ac561373cd3be920984c3df29a1c2ecc8 (2019-02-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JaD_transl = S6PRRPRP2_jacobiaD_transl_6_sym_varpar(qJ, qJD, r_i_i_C, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(3,1),zeros(11,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRRPRP2_jacobiaD_transl_6_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6PRRPRP2_jacobiaD_transl_6_sym_varpar: qJD has to be [6x1] (double)');
assert(isa(r_i_i_C,'double') && isreal(r_i_i_C) && all(size(r_i_i_C) == [3 1]), ...
	'S6PRRPRP2_jacobiaD_transl_6_sym_varpar: Position vector r_i_i_C has to be [3x1] double');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6PRRPRP2_jacobiaD_transl_6_sym_varpar: pkin has to be [11x1] (double)');

%% Symbolic Calculation
% From jacobiaD_transl_6_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 20:01:50
% EndTime: 2019-02-26 20:01:50
% DurationCPUTime: 0.45s
% Computational Cost: add. (699->106), mult. (1453->180), div. (0->0), fcn. (1497->12), ass. (0->67)
t380 = sin(qJ(2));
t383 = cos(qJ(2));
t375 = sin(pkin(10));
t422 = cos(pkin(6));
t406 = t375 * t422;
t421 = cos(pkin(10));
t365 = -t380 * t406 + t421 * t383;
t374 = qJ(3) + pkin(11);
t372 = sin(t374);
t373 = cos(t374);
t376 = sin(pkin(6));
t417 = t376 * t380;
t357 = t422 * t372 + t373 * t417;
t381 = cos(qJ(5));
t378 = sin(qJ(5));
t415 = t378 * t383;
t428 = -t357 * t381 + t376 * t415;
t379 = sin(qJ(3));
t425 = pkin(9) + r_i_i_C(2);
t427 = -pkin(3) * t379 - pkin(4) * t372 + t425 * t373;
t412 = qJD(3) * t383;
t426 = (qJD(2) * t373 - qJD(5)) * t380 + t372 * t412;
t424 = r_i_i_C(1) + pkin(5);
t423 = r_i_i_C(3) + qJ(6);
t419 = t373 * t378;
t418 = t375 * t376;
t416 = t376 * t383;
t414 = qJD(2) * t380;
t413 = qJD(3) * t372;
t411 = qJD(5) * t373;
t409 = qJD(2) * t416;
t408 = t376 * t414;
t405 = t376 * t421;
t398 = t422 * t421;
t394 = t383 * t398;
t358 = -qJD(2) * t394 + t375 * t414;
t362 = t375 * t380 - t394;
t400 = t362 * t411 - t358;
t364 = t421 * t380 + t383 * t406;
t360 = t364 * qJD(2);
t399 = t364 * t411 - t360;
t363 = t375 * t383 + t380 * t398;
t389 = -t363 * t373 + t372 * t405;
t397 = t362 * t378 - t381 * t389;
t353 = t365 * t373 + t372 * t418;
t396 = t353 * t381 + t364 * t378;
t395 = (qJD(2) - t411) * t383;
t393 = -t365 * t372 + t373 * t418;
t382 = cos(qJ(3));
t392 = -t382 * pkin(3) - t373 * pkin(4) - t425 * t372 - pkin(2);
t391 = -t363 * t372 - t373 * t405;
t390 = -t372 * t417 + t422 * t373;
t388 = t423 * t378 + t424 * t381 + pkin(4);
t359 = t363 * qJD(2);
t387 = qJD(5) * t363 - t359 * t373 + t362 * t413;
t361 = t365 * qJD(2);
t386 = qJD(5) * t365 - t361 * t373 + t364 * t413;
t385 = qJD(3) * t427;
t384 = qJD(6) * t378 + (-t424 * t378 + t423 * t381) * qJD(5);
t377 = -qJ(4) - pkin(8);
t349 = t390 * qJD(3) + t373 * t409;
t347 = t393 * qJD(3) - t360 * t373;
t345 = t391 * qJD(3) - t358 * t373;
t340 = -qJD(5) * t428 + t349 * t378 - t381 * t408;
t334 = t396 * qJD(5) + t347 * t378 - t361 * t381;
t332 = t397 * qJD(5) + t345 * t378 - t359 * t381;
t1 = [0 -(t364 * t419 + t365 * t381) * qJD(6) + t360 * t377 + t365 * qJD(4) + t424 * (t399 * t378 + t386 * t381) + t423 * (t386 * t378 - t399 * t381) + t392 * t361 - t364 * t385, t425 * t347 + (t360 * t379 + (-t365 * t382 - t379 * t418) * qJD(3)) * pkin(3) + t384 * t393 + t388 * (-t353 * qJD(3) + t360 * t372) t361, t396 * qJD(6) + t423 * (t347 * t381 + t361 * t378 + (-t353 * t378 + t364 * t381) * qJD(5)) - t424 * t334, t334; 0 -(t362 * t419 + t363 * t381) * qJD(6) + t358 * t377 + t363 * qJD(4) + t424 * (t400 * t378 + t387 * t381) + t423 * (t387 * t378 - t400 * t381) + t392 * t359 - t362 * t385, t425 * t345 + (t358 * t379 + (-t363 * t382 + t379 * t405) * qJD(3)) * pkin(3) + t384 * t391 + t388 * (t389 * qJD(3) + t358 * t372) t359, t397 * qJD(6) + t423 * (t345 * t381 + t359 * t378 + (t362 * t381 + t378 * t389) * qJD(5)) - t424 * t332, t332; 0 (t424 * (t378 * t395 - t381 * t426) - t423 * (t378 * t426 + t381 * t395) - (-t373 * t415 + t380 * t381) * qJD(6) + t380 * qJD(4) + t427 * t412 + (-t383 * t377 + t392 * t380) * qJD(2)) * t376, t425 * t349 + (-t379 * t409 + (-t422 * t379 - t382 * t417) * qJD(3)) * pkin(3) + t384 * t390 + t388 * (-t357 * qJD(3) - t372 * t409) t408, -t428 * qJD(6) + t423 * (t378 * t408 + t349 * t381 + (-t357 * t378 - t381 * t416) * qJD(5)) - t424 * t340, t340;];
JaD_transl  = t1;
