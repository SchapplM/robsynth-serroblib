% Zeitableitung der translatorischen Teilmatrix der analytischen Jacobi-Matrix für Segment Nr. 5 (0=Basis) von
% S6RRRRPR11
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
%   pkin=[a2,a3,a4,a5,a6,alpha2,d1,d2,d3,d4,d6,theta5]';
%
% Output:
% JaD_transl [3x6]
%   Zeitableitung der translatorischen Teilmatrix der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-02-26 22:36
% Revision: d75aae1ac561373cd3be920984c3df29a1c2ecc8 (2019-02-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JaD_transl = S6RRRRPR11_jacobiaD_transl_5_sym_varpar(qJ, qJD, r_i_i_C, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(3,1),zeros(12,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRRPR11_jacobiaD_transl_5_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRRRPR11_jacobiaD_transl_5_sym_varpar: qJD has to be [6x1] (double)');
assert(isa(r_i_i_C,'double') && isreal(r_i_i_C) && all(size(r_i_i_C) == [3 1]), ...
	'S6RRRRPR11_jacobiaD_transl_5_sym_varpar: Position vector r_i_i_C has to be [3x1] double');
assert(isreal(pkin) && all(size(pkin) == [12 1]), ...
  'S6RRRRPR11_jacobiaD_transl_5_sym_varpar: pkin has to be [12x1] (double)');

%% Symbolic Calculation
% From jacobiaD_transl_5_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 22:36:23
% EndTime: 2019-02-26 22:36:24
% DurationCPUTime: 0.62s
% Computational Cost: add. (669->126), mult. (1706->207), div. (0->0), fcn. (1716->12), ass. (0->69)
t386 = sin(qJ(1));
t389 = cos(qJ(2));
t426 = cos(pkin(6));
t432 = cos(qJ(1));
t400 = t426 * t432;
t385 = sin(qJ(2));
t406 = t386 * t426;
t401 = t385 * t406;
t407 = t432 * qJD(1);
t419 = qJD(2) * t385;
t353 = -qJD(1) * t401 - t386 * t419 + (qJD(2) * t400 + t407) * t389;
t381 = sin(pkin(6));
t412 = t381 * t432;
t439 = -qJD(3) * t412 + t353;
t364 = t385 * t400 + t386 * t389;
t423 = t381 * t386;
t438 = qJD(1) * t423 - qJD(3) * t364;
t384 = sin(qJ(3));
t388 = cos(qJ(3));
t357 = t364 * t388 - t384 * t412;
t398 = t389 * t400;
t420 = t386 * t385;
t363 = -t398 + t420;
t380 = qJ(4) + pkin(12);
t378 = sin(t380);
t379 = cos(t380);
t437 = t357 * t378 - t363 * t379;
t436 = t357 * t379 + t363 * t378;
t387 = cos(qJ(4));
t428 = t387 * pkin(4);
t377 = pkin(3) + t428;
t399 = t379 * r_i_i_C(1) - t378 * r_i_i_C(2);
t397 = t377 + t399;
t427 = r_i_i_C(3) + qJ(5) + pkin(10);
t435 = -(t397 * t384 - t427 * t388) * qJD(3) + t384 * qJD(5);
t383 = sin(qJ(4));
t395 = t383 * pkin(4) + t378 * r_i_i_C(1) + t379 * r_i_i_C(2);
t347 = t438 * t384 + t439 * t388;
t433 = t427 * t384 + t397 * t388 + pkin(2);
t422 = t381 * t388;
t421 = t381 * t389;
t417 = qJD(4) * t378;
t416 = qJD(4) * t379;
t415 = qJD(4) * t387;
t414 = qJD(4) * t388;
t411 = t432 * t389;
t409 = t381 * t419;
t408 = qJD(2) * t421;
t402 = t381 * t407;
t366 = t411 - t401;
t359 = -t366 * t384 + t386 * t422;
t360 = t366 * t388 + t384 * t423;
t394 = t364 * t384 + t388 * t412;
t393 = -t381 * t385 * t384 + t426 * t388;
t362 = t426 * t384 + t385 * t422;
t346 = t439 * t384 - t438 * t388;
t392 = qJD(4) * t395;
t365 = t432 * t385 + t389 * t406;
t390 = t395 * t414 - t435;
t355 = t393 * qJD(3) + t388 * t408;
t354 = t362 * qJD(3) + t384 * t408;
t352 = t365 * qJD(1) + t364 * qJD(2);
t351 = t364 * qJD(1) + t365 * qJD(2);
t350 = -qJD(1) * t398 - qJD(2) * t411 + (qJD(2) * t426 + qJD(1)) * t420;
t345 = t359 * qJD(3) - t351 * t388 + t384 * t402;
t344 = t360 * qJD(3) - t351 * t384 - t388 * t402;
t343 = t345 * t379 - t350 * t378 + (-t360 * t378 + t365 * t379) * qJD(4);
t342 = -t345 * t378 - t350 * t379 + (-t360 * t379 - t365 * t378) * qJD(4);
t1 = [-t394 * qJD(5) - t353 * pkin(2) - t397 * t347 + (-pkin(9) - t395) * t352 - t427 * t346 + (-t432 * pkin(1) - pkin(8) * t423) * qJD(1) + (t437 * r_i_i_C(1) + t436 * r_i_i_C(2) + (t357 * t383 - t363 * t387) * pkin(4)) * qJD(4) (-t351 * t378 + t366 * t416) * r_i_i_C(1) + (-t351 * t379 - t366 * t417) * r_i_i_C(2) - t351 * pkin(9) + (-t351 * t383 + t366 * t415) * pkin(4) + t433 * t350 + t390 * t365, t360 * qJD(5) - t397 * t344 + t427 * t345 - t359 * t392, t342 * r_i_i_C(1) - t343 * r_i_i_C(2) + (-t345 * t383 - t350 * t387 + (-t360 * t387 - t365 * t383) * qJD(4)) * pkin(4), t344, 0; -t351 * pkin(2) - t350 * pkin(9) + t343 * r_i_i_C(1) + t342 * r_i_i_C(2) - t359 * qJD(5) + t345 * t377 + t427 * t344 + (-pkin(1) * t386 + pkin(8) * t412) * qJD(1) + (-t350 * t383 + (-t360 * t383 + t365 * t387) * qJD(4)) * pkin(4) (t353 * t378 + t364 * t416) * r_i_i_C(1) + (t353 * t379 - t364 * t417) * r_i_i_C(2) + t353 * pkin(9) + (t353 * t383 + t364 * t415) * pkin(4) - t433 * t352 + t390 * t363, t357 * qJD(5) - t397 * t346 + t427 * t347 + t394 * t392 (-t347 * t378 + t352 * t379) * r_i_i_C(1) + (-t347 * t379 - t352 * t378) * r_i_i_C(2) + (-t436 * r_i_i_C(1) + t437 * r_i_i_C(2)) * qJD(4) + (-t347 * t383 + t352 * t387 + (-t357 * t387 - t363 * t383) * qJD(4)) * pkin(4), t346, 0; 0 (((t399 + t428) * qJD(4) - t433 * qJD(2)) * t385 + (qJD(2) * pkin(9) + t395 * (qJD(2) - t414) + t435) * t389) * t381, t362 * qJD(5) - t397 * t354 + t427 * t355 - t393 * t392 (-t355 * t378 + t379 * t409) * r_i_i_C(1) + (-t355 * t379 - t378 * t409) * r_i_i_C(2) + ((-t362 * t379 + t378 * t421) * r_i_i_C(1) + (t362 * t378 + t379 * t421) * r_i_i_C(2)) * qJD(4) + (t387 * t409 - t355 * t383 + (-t362 * t387 + t383 * t421) * qJD(4)) * pkin(4), t354, 0;];
JaD_transl  = t1;
