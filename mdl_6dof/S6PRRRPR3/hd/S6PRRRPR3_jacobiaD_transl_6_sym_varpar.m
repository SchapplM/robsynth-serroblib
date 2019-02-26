% Zeitableitung der translatorischen Teilmatrix der analytischen Jacobi-Matrix für Segment Nr. 6 (0=Basis) von
% S6PRRRPR3
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
%   pkin=[a2,a3,a4,a5,a6,alpha2,d2,d3,d4,d6,theta1]';
%
% Output:
% JaD_transl [3x6]
%   Zeitableitung der translatorischen Teilmatrix der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-02-26 20:11
% Revision: d75aae1ac561373cd3be920984c3df29a1c2ecc8 (2019-02-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JaD_transl = S6PRRRPR3_jacobiaD_transl_6_sym_varpar(qJ, qJD, r_i_i_C, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(3,1),zeros(11,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRRRPR3_jacobiaD_transl_6_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6PRRRPR3_jacobiaD_transl_6_sym_varpar: qJD has to be [6x1] (double)');
assert(isa(r_i_i_C,'double') && isreal(r_i_i_C) && all(size(r_i_i_C) == [3 1]), ...
	'S6PRRRPR3_jacobiaD_transl_6_sym_varpar: Position vector r_i_i_C has to be [3x1] double');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6PRRRPR3_jacobiaD_transl_6_sym_varpar: pkin has to be [11x1] (double)');

%% Symbolic Calculation
% From jacobiaD_transl_6_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 20:11:45
% EndTime: 2019-02-26 20:11:45
% DurationCPUTime: 0.54s
% Computational Cost: add. (734->80), mult. (1250->142), div. (0->0), fcn. (1245->12), ass. (0->59)
t383 = sin(qJ(6));
t386 = cos(qJ(6));
t404 = t386 * r_i_i_C(1) - t383 * r_i_i_C(2);
t433 = t404 * qJD(6) + qJD(5);
t403 = -t383 * r_i_i_C(1) - t386 * r_i_i_C(2);
t434 = qJ(5) - t403;
t417 = pkin(4) + pkin(10) + r_i_i_C(3);
t380 = qJ(3) + qJ(4);
t377 = sin(t380);
t379 = qJD(3) + qJD(4);
t384 = sin(qJ(3));
t378 = cos(t380);
t422 = t378 * t379;
t390 = -(-t417 * t379 + t433) * t377 + qJD(3) * t384 * pkin(3) - t434 * t422;
t385 = sin(qJ(2));
t388 = cos(qJ(2));
t381 = sin(pkin(11));
t425 = cos(pkin(6));
t410 = t381 * t425;
t424 = cos(pkin(11));
t368 = -t385 * t410 + t424 * t388;
t423 = t377 * t379;
t382 = sin(pkin(6));
t421 = t381 * t382;
t420 = t382 * t385;
t419 = t382 * t388;
t418 = qJD(2) * t385;
t414 = t377 * t420;
t413 = t378 * t420;
t412 = t382 * t418;
t411 = qJD(2) * t419;
t409 = t382 * t424;
t406 = t377 * t409;
t405 = t378 * t409;
t367 = t424 * t385 + t388 * t410;
t363 = t367 * qJD(2);
t401 = t379 * t421 - t363;
t400 = t425 * t424;
t398 = t388 * t400;
t397 = t403 * qJD(6);
t396 = pkin(5) + pkin(9) + pkin(8) + t404;
t395 = t425 * t379 + t411;
t366 = t381 * t388 + t385 * t400;
t387 = cos(qJ(3));
t394 = -pkin(3) * t387 - t377 * t434 - t417 * t378 - pkin(2);
t342 = t368 * t422 + t401 * t377;
t393 = t433 * (t368 * t378 + t377 * t421) + t434 * (-t368 * t423 + t401 * t378) - t417 * t342;
t361 = -qJD(2) * t398 + t381 * t418;
t340 = -t361 * t377 + t366 * t422 - t379 * t406;
t392 = t433 * (t366 * t378 - t406) + t434 * (-t361 * t378 - t366 * t423 - t379 * t405) - t417 * t340;
t348 = t395 * t377 + t379 * t413;
t391 = t433 * (t425 * t377 + t413) + t434 * (t395 * t378 - t379 * t414) - t417 * t348;
t365 = t381 * t385 - t398;
t364 = t368 * qJD(2);
t362 = t366 * qJD(2);
t359 = -t425 * t378 + t414;
t354 = t368 * t377 - t378 * t421;
t352 = t366 * t377 + t405;
t1 = [0, -t396 * t363 + t394 * t364 + t390 * t367 + t368 * t397 (t363 * t384 + (-t368 * t387 - t384 * t421) * qJD(3)) * pkin(3) + t393, t393, t342 (t342 * t386 - t364 * t383) * r_i_i_C(1) + (-t342 * t383 - t364 * t386) * r_i_i_C(2) + ((-t354 * t383 - t367 * t386) * r_i_i_C(1) + (-t354 * t386 + t367 * t383) * r_i_i_C(2)) * qJD(6); 0, -t396 * t361 + t394 * t362 + t390 * t365 + t366 * t397 (t361 * t384 + (-t366 * t387 + t384 * t409) * qJD(3)) * pkin(3) + t392, t392, t340 (t340 * t386 - t362 * t383) * r_i_i_C(1) + (-t340 * t383 - t362 * t386) * r_i_i_C(2) + ((-t352 * t383 - t365 * t386) * r_i_i_C(1) + (-t352 * t386 + t365 * t383) * r_i_i_C(2)) * qJD(6); 0 ((t394 * qJD(2) + t397) * t385 + (t396 * qJD(2) - t390) * t388) * t382 (-t384 * t411 + (-t425 * t384 - t387 * t420) * qJD(3)) * pkin(3) + t391, t391, t348 (t348 * t386 - t383 * t412) * r_i_i_C(1) + (-t348 * t383 - t386 * t412) * r_i_i_C(2) + ((-t359 * t383 + t386 * t419) * r_i_i_C(1) + (-t359 * t386 - t383 * t419) * r_i_i_C(2)) * qJD(6);];
JaD_transl  = t1;
