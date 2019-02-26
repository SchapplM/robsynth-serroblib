% Zeitableitung der translatorischen Teilmatrix der analytischen Jacobi-Matrix für Segment Nr. 5 (0=Basis) von
% S6RRRPRR8
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
% pkin [12x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d1,d2,d3,d5,d6,theta4]';
%
% Output:
% JaD_transl [3x6]
%   Zeitableitung der translatorischen Teilmatrix der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-02-26 22:20
% Revision: d75aae1ac561373cd3be920984c3df29a1c2ecc8 (2019-02-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JaD_transl = S6RRRPRR8_jacobiaD_transl_5_sym_varpar(qJ, qJD, r_i_i_C, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(3,1),zeros(12,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRPRR8_jacobiaD_transl_5_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRRPRR8_jacobiaD_transl_5_sym_varpar: qJD has to be [6x1] (double)');
assert(isa(r_i_i_C,'double') && isreal(r_i_i_C) && all(size(r_i_i_C) == [3 1]), ...
	'S6RRRPRR8_jacobiaD_transl_5_sym_varpar: Position vector r_i_i_C has to be [3x1] double');
assert(isreal(pkin) && all(size(pkin) == [12 1]), ...
  'S6RRRPRR8_jacobiaD_transl_5_sym_varpar: pkin has to be [12x1] (double)');

%% Symbolic Calculation
% From jacobiaD_transl_5_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 22:19:49
% EndTime: 2019-02-26 22:19:50
% DurationCPUTime: 0.68s
% Computational Cost: add. (686->123), mult. (1483->200), div. (0->0), fcn. (1476->12), ass. (0->71)
t393 = sin(qJ(1));
t388 = cos(pkin(6));
t409 = qJD(2) * t388 + qJD(1);
t392 = sin(qJ(2));
t427 = t393 * t392;
t417 = t388 * t427;
t422 = qJD(2) * t392;
t396 = cos(qJ(2));
t397 = cos(qJ(1));
t424 = t397 * t396;
t362 = -qJD(1) * t417 - t393 * t422 + t409 * t424;
t425 = t397 * t392;
t426 = t393 * t396;
t373 = t388 * t425 + t426;
t386 = qJ(3) + pkin(12);
t384 = sin(t386);
t385 = cos(t386);
t387 = sin(pkin(6));
t423 = qJD(1) * t387;
t414 = t393 * t423;
t428 = t387 * t397;
t416 = t385 * t428;
t356 = (-qJD(3) * t373 + t414) * t384 - qJD(3) * t416 + t362 * t385;
t374 = t388 * t426 + t425;
t361 = t374 * qJD(1) + t373 * qJD(2);
t390 = sin(qJ(5));
t394 = cos(qJ(5));
t447 = t356 * t390 - t361 * t394;
t446 = -t356 * t394 - t361 * t390;
t391 = sin(qJ(3));
t407 = t394 * r_i_i_C(1) - t390 * r_i_i_C(2);
t405 = pkin(4) + t407;
t440 = pkin(10) + r_i_i_C(3);
t445 = (t391 * pkin(3) + t405 * t384 - t440 * t385) * qJD(3);
t367 = -t373 * t385 + t384 * t428;
t415 = t388 * t424;
t372 = -t415 + t427;
t444 = -t367 * t390 - t372 * t394;
t443 = t367 * t394 - t372 * t390;
t406 = t390 * r_i_i_C(1) + t394 * r_i_i_C(2);
t395 = cos(qJ(3));
t383 = t395 * pkin(3) + pkin(2);
t442 = t440 * t384 + t405 * t385 + t383;
t432 = t387 * t392;
t431 = t387 * t393;
t430 = t387 * t395;
t429 = t387 * t396;
t421 = qJD(2) * t396;
t420 = qJD(5) * t385;
t419 = qJD(5) * t390;
t418 = qJD(5) * t394;
t413 = t397 * t423;
t412 = t387 * t421;
t411 = t387 * t422;
t375 = -t417 + t424;
t404 = -t375 * t384 + t385 * t431;
t369 = t375 * t385 + t384 * t431;
t371 = t388 * t384 + t385 * t432;
t403 = -t384 * t432 + t388 * t385;
t402 = qJD(5) * t406;
t399 = t367 * qJD(3) - t362 * t384 + t385 * t414;
t398 = t406 * t420 + t445;
t389 = -qJ(4) - pkin(9);
t364 = t403 * qJD(3) + t385 * t412;
t360 = t373 * qJD(1) + t374 * qJD(2);
t359 = -qJD(1) * t415 - t397 * t421 + t409 * t427;
t354 = t404 * qJD(3) - t360 * t385 + t384 * t413;
t353 = t369 * qJD(3) - t360 * t384 - t385 * t413;
t352 = t354 * t394 - t359 * t390 + (-t369 * t390 + t374 * t394) * qJD(5);
t351 = -t354 * t390 - t359 * t394 + (-t369 * t394 - t374 * t390) * qJD(5);
t1 = [t446 * r_i_i_C(1) + t447 * r_i_i_C(2) - t356 * pkin(4) - t362 * t383 + t361 * t389 - t372 * qJD(4) + t440 * t399 + (t444 * r_i_i_C(1) - t443 * r_i_i_C(2)) * qJD(5) + (-t397 * pkin(1) - pkin(8) * t431) * qJD(1) + (-t391 * t414 + (t373 * t391 + t395 * t428) * qJD(3)) * pkin(3) (-t360 * t390 + t375 * t418) * r_i_i_C(1) + (-t360 * t394 - t375 * t419) * r_i_i_C(2) + t360 * t389 + t375 * qJD(4) + t442 * t359 + t398 * t374, t440 * t354 - t404 * t402 - t405 * t353 + (t395 * t413 + t360 * t391 + (-t375 * t395 - t391 * t431) * qJD(3)) * pkin(3), -t359, t351 * r_i_i_C(1) - t352 * r_i_i_C(2), 0; t354 * pkin(4) + t352 * r_i_i_C(1) + t351 * r_i_i_C(2) + t374 * qJD(4) + t359 * t389 - t360 * t383 + t440 * t353 + (-pkin(1) * t393 + pkin(8) * t428) * qJD(1) + (t391 * t413 + (-t375 * t391 + t393 * t430) * qJD(3)) * pkin(3) (t362 * t390 + t373 * t418) * r_i_i_C(1) + (t362 * t394 - t373 * t419) * r_i_i_C(2) - t362 * t389 + t373 * qJD(4) - t442 * t361 + t398 * t372, t440 * t356 - (-t373 * t384 - t416) * t402 + t405 * t399 + (t395 * t414 - t362 * t391 + (-t373 * t395 + t391 * t428) * qJD(3)) * pkin(3), t361, -t447 * r_i_i_C(1) + t446 * r_i_i_C(2) + (t443 * r_i_i_C(1) + t444 * r_i_i_C(2)) * qJD(5), 0; 0 ((-qJD(2) * t442 + t407 * qJD(5) + qJD(4)) * t392 + (-qJD(2) * t389 - t445 + t406 * (qJD(2) - t420)) * t396) * t387, t440 * t364 - t403 * t402 + t405 * (-t371 * qJD(3) - t384 * t412) + (-t391 * t412 + (-t388 * t391 - t392 * t430) * qJD(3)) * pkin(3), t411 (-t364 * t390 + t394 * t411) * r_i_i_C(1) + (-t364 * t394 - t390 * t411) * r_i_i_C(2) + ((-t371 * t394 + t390 * t429) * r_i_i_C(1) + (t371 * t390 + t394 * t429) * r_i_i_C(2)) * qJD(5), 0;];
JaD_transl  = t1;
