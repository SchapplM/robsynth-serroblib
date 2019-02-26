% Zeitableitung der translatorischen Teilmatrix der analytischen Jacobi-Matrix für Segment Nr. 5 (0=Basis) von
% S6PPPRRR1
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
% pkin [14x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,alpha3,alpha4,d4,d5,d6,theta1,theta2,theta3]';
%
% Output:
% JaD_transl [3x6]
%   Zeitableitung der translatorischen Teilmatrix der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-02-26 19:38
% Revision: d75aae1ac561373cd3be920984c3df29a1c2ecc8 (2019-02-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JaD_transl = S6PPPRRR1_jacobiaD_transl_5_sym_varpar(qJ, qJD, r_i_i_C, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(3,1),zeros(14,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PPPRRR1_jacobiaD_transl_5_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6PPPRRR1_jacobiaD_transl_5_sym_varpar: qJD has to be [6x1] (double)');
assert(isa(r_i_i_C,'double') && isreal(r_i_i_C) && all(size(r_i_i_C) == [3 1]), ...
	'S6PPPRRR1_jacobiaD_transl_5_sym_varpar: Position vector r_i_i_C has to be [3x1] double');
assert(isreal(pkin) && all(size(pkin) == [14 1]), ...
  'S6PPPRRR1_jacobiaD_transl_5_sym_varpar: pkin has to be [14x1] (double)');

%% Symbolic Calculation
% From jacobiaD_transl_5_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 19:38:45
% EndTime: 2019-02-26 19:38:45
% DurationCPUTime: 0.22s
% Computational Cost: add. (410->57), mult. (1254->116), div. (0->0), fcn. (1562->16), ass. (0->58)
t391 = sin(pkin(14));
t392 = sin(pkin(13));
t396 = sin(pkin(6));
t397 = cos(pkin(14));
t398 = cos(pkin(13));
t401 = cos(pkin(7));
t420 = t398 * t401;
t395 = sin(pkin(7));
t402 = cos(pkin(6));
t422 = t395 * t402;
t382 = t396 * t392 * t397 + (t396 * t420 + t422) * t391;
t404 = sin(qJ(4));
t406 = cos(qJ(4));
t381 = t397 * t422 + (-t391 * t392 + t397 * t420) * t396;
t423 = t395 * t396;
t386 = -t398 * t423 + t402 * t401;
t394 = sin(pkin(8));
t400 = cos(pkin(8));
t411 = t381 * t400 + t386 * t394;
t368 = t382 * t406 + t411 * t404;
t399 = cos(pkin(12));
t393 = sin(pkin(12));
t424 = t393 * t402;
t390 = -t392 * t424 + t399 * t398;
t389 = -t399 * t392 - t398 * t424;
t409 = t389 * t401 + t393 * t423;
t376 = t390 * t397 + t409 * t391;
t375 = -t390 * t391 + t409 * t397;
t421 = t396 * t401;
t384 = -t389 * t395 + t393 * t421;
t412 = t375 * t400 + t384 * t394;
t364 = t376 * t406 + t412 * t404;
t419 = t399 * t402;
t388 = t392 * t419 + t393 * t398;
t387 = -t393 * t392 + t398 * t419;
t410 = t387 * t401 - t399 * t423;
t374 = t388 * t397 + t410 * t391;
t373 = -t388 * t391 + t410 * t397;
t383 = -t387 * t395 - t399 * t421;
t413 = t373 * t400 + t383 * t394;
t362 = t374 * t406 + t413 * t404;
t428 = -pkin(10) - r_i_i_C(3);
t418 = qJD(4) * t404;
t417 = qJD(4) * t406;
t416 = t394 * t417;
t415 = t400 * t417;
t403 = sin(qJ(5));
t405 = cos(qJ(5));
t414 = t403 * r_i_i_C(1) + t405 * r_i_i_C(2);
t408 = qJD(5) * t414;
t407 = qJD(4) * (t405 * r_i_i_C(1) - t403 * r_i_i_C(2) + pkin(4));
t377 = -t381 * t394 + t386 * t400;
t370 = -t375 * t394 + t384 * t400;
t369 = -t373 * t394 + t383 * t400;
t365 = -t381 * t415 + t382 * t418 - t386 * t416;
t359 = -t375 * t415 + t376 * t418 - t384 * t416;
t357 = -t373 * t415 + t374 * t418 - t383 * t416;
t1 = [0, 0, 0, t428 * t359 - (-t376 * t404 + t412 * t406) * t408 - t364 * t407, t414 * t359 + ((-t364 * t405 - t370 * t403) * r_i_i_C(1) + (t364 * t403 - t370 * t405) * r_i_i_C(2)) * qJD(5), 0; 0, 0, 0, t428 * t357 - (-t374 * t404 + t413 * t406) * t408 - t362 * t407, t414 * t357 + ((-t362 * t405 - t369 * t403) * r_i_i_C(1) + (t362 * t403 - t369 * t405) * r_i_i_C(2)) * qJD(5), 0; 0, 0, 0, t428 * t365 - (-t382 * t404 + t411 * t406) * t408 - t368 * t407, t414 * t365 + ((-t368 * t405 - t377 * t403) * r_i_i_C(1) + (t368 * t403 - t377 * t405) * r_i_i_C(2)) * qJD(5), 0;];
JaD_transl  = t1;
