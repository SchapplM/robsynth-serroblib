% Zeitableitung der translatorischen Teilmatrix der analytischen Jacobi-Matrix für Segment Nr. 6 (0=Basis) von
% S6PRPRPR1
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
%   pkin=[a2,a3,a4,a5,a6,alpha2,d2,d4,d6,theta1,theta3,theta5]';
%
% Output:
% JaD_transl [3x6]
%   Zeitableitung der translatorischen Teilmatrix der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-02-26 19:46
% Revision: d75aae1ac561373cd3be920984c3df29a1c2ecc8 (2019-02-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JaD_transl = S6PRPRPR1_jacobiaD_transl_6_sym_varpar(qJ, qJD, r_i_i_C, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(3,1),zeros(12,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRPRPR1_jacobiaD_transl_6_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6PRPRPR1_jacobiaD_transl_6_sym_varpar: qJD has to be [6x1] (double)');
assert(isa(r_i_i_C,'double') && isreal(r_i_i_C) && all(size(r_i_i_C) == [3 1]), ...
	'S6PRPRPR1_jacobiaD_transl_6_sym_varpar: Position vector r_i_i_C has to be [3x1] double');
assert(isreal(pkin) && all(size(pkin) == [12 1]), ...
  'S6PRPRPR1_jacobiaD_transl_6_sym_varpar: pkin has to be [12x1] (double)');

%% Symbolic Calculation
% From jacobiaD_transl_6_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 19:46:29
% EndTime: 2019-02-26 19:46:30
% DurationCPUTime: 0.54s
% Computational Cost: add. (608->99), mult. (1513->178), div. (0->0), fcn. (1626->14), ass. (0->61)
t390 = qJ(4) + pkin(12);
t388 = sin(t390);
t389 = cos(t390);
t398 = sin(qJ(4));
t397 = sin(qJ(6));
t400 = cos(qJ(6));
t408 = (t397 * r_i_i_C(1) + t400 * r_i_i_C(2)) * qJD(6);
t413 = t400 * r_i_i_C(1) - t397 * r_i_i_C(2) + pkin(5);
t432 = pkin(9) + r_i_i_C(3);
t436 = (t398 * pkin(4) + t413 * t388 - t432 * t389) * qJD(4) + t389 * t408;
t395 = cos(pkin(6));
t391 = sin(pkin(11));
t399 = sin(qJ(2));
t428 = cos(pkin(11));
t431 = cos(qJ(2));
t407 = t431 * t391 + t399 * t428;
t376 = t407 * t395;
t416 = t431 * t428;
t422 = qJD(2) * t399;
t434 = -qJD(2) * t416 + t391 * t422;
t406 = -t399 * t391 + t416;
t401 = cos(qJ(4));
t403 = t401 * pkin(4) + t432 * t388 + t413 * t389 + pkin(3);
t429 = pkin(2) * qJD(2);
t392 = sin(pkin(10));
t393 = sin(pkin(6));
t427 = t392 * t393;
t394 = cos(pkin(10));
t426 = t393 * t394;
t425 = t393 * t398;
t424 = t395 * t399;
t421 = qJD(6) * t397;
t420 = qJD(6) * t400;
t373 = t434 * t395;
t378 = t407 * qJD(2);
t355 = t394 * t373 + t392 * t378;
t357 = t392 * t373 - t394 * t378;
t375 = t407 * t393;
t366 = t375 * t389 + t395 * t388;
t414 = -t375 * t388 + t395 * t389;
t361 = t394 * t376 + t392 * t406;
t362 = t392 * t376 - t394 * t406;
t411 = -t361 * t388 - t389 * t426;
t410 = -t361 * t389 + t388 * t426;
t409 = t362 * t388 + t389 * t427;
t352 = -t362 * t389 + t388 * t427;
t405 = t406 * t395;
t404 = qJD(2) * t376;
t396 = -qJ(5) - pkin(8);
t377 = t406 * qJD(2);
t374 = t406 * t393;
t372 = qJD(2) * t375;
t371 = t434 * t393;
t363 = -t392 * t405 - t394 * t407;
t360 = -t392 * t407 + t394 * t405;
t356 = -t394 * t377 + t392 * t404;
t353 = -t392 * t377 - t394 * t404;
t348 = t414 * qJD(4) - t371 * t389;
t346 = t409 * qJD(4) + t357 * t389;
t344 = t411 * qJD(4) - t355 * t389;
t1 = [0 (t357 * t397 - t362 * t420) * r_i_i_C(1) + (t357 * t400 + t362 * t421) * r_i_i_C(2) - t357 * t396 - t362 * qJD(5) + (t392 * t424 - t431 * t394) * t429 + t403 * t356 - t436 * t363, 0, t432 * t346 - t409 * t408 + t413 * (-t352 * qJD(4) - t357 * t388) + (-t357 * t398 + (t362 * t401 - t392 * t425) * qJD(4)) * pkin(4), -t356 (-t346 * t397 - t356 * t400) * r_i_i_C(1) + (-t346 * t400 + t356 * t397) * r_i_i_C(2) + ((-t352 * t400 + t363 * t397) * r_i_i_C(1) + (t352 * t397 + t363 * t400) * r_i_i_C(2)) * qJD(6); 0 (-t355 * t397 + t361 * t420) * r_i_i_C(1) + (-t355 * t400 - t361 * t421) * r_i_i_C(2) + t355 * t396 + t361 * qJD(5) + (-t431 * t392 - t394 * t424) * t429 + t403 * t353 - t436 * t360, 0, t432 * t344 - t411 * t408 + t413 * (t410 * qJD(4) + t355 * t388) + (t355 * t398 + (-t361 * t401 + t394 * t425) * qJD(4)) * pkin(4), -t353 (-t344 * t397 - t353 * t400) * r_i_i_C(1) + (-t344 * t400 + t353 * t397) * r_i_i_C(2) + ((t360 * t397 + t400 * t410) * r_i_i_C(1) + (t360 * t400 - t397 * t410) * r_i_i_C(2)) * qJD(6); 0 (-t371 * t397 + t375 * t420) * r_i_i_C(1) + (-t371 * t400 - t375 * t421) * r_i_i_C(2) + t371 * t396 + t375 * qJD(5) - t393 * pkin(2) * t422 - t403 * t372 - t436 * t374, 0, t432 * t348 - t414 * t408 + t413 * (-t366 * qJD(4) + t371 * t388) + (t371 * t398 + (-t375 * t401 - t395 * t398) * qJD(4)) * pkin(4), t372 (-t348 * t397 + t372 * t400) * r_i_i_C(1) + (-t348 * t400 - t372 * t397) * r_i_i_C(2) + ((-t366 * t400 + t374 * t397) * r_i_i_C(1) + (t366 * t397 + t374 * t400) * r_i_i_C(2)) * qJD(6);];
JaD_transl  = t1;
