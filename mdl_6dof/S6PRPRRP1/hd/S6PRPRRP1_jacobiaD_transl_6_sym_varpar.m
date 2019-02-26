% Zeitableitung der translatorischen Teilmatrix der analytischen Jacobi-Matrix für Segment Nr. 6 (0=Basis) von
% S6PRPRRP1
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
%   pkin=[a2,a3,a4,a5,a6,alpha2,d2,d4,d5,theta1,theta3]';
%
% Output:
% JaD_transl [3x6]
%   Zeitableitung der translatorischen Teilmatrix der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-02-26 19:50
% Revision: d75aae1ac561373cd3be920984c3df29a1c2ecc8 (2019-02-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JaD_transl = S6PRPRRP1_jacobiaD_transl_6_sym_varpar(qJ, qJD, r_i_i_C, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(3,1),zeros(11,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRPRRP1_jacobiaD_transl_6_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6PRPRRP1_jacobiaD_transl_6_sym_varpar: qJD has to be [6x1] (double)');
assert(isa(r_i_i_C,'double') && isreal(r_i_i_C) && all(size(r_i_i_C) == [3 1]), ...
	'S6PRPRRP1_jacobiaD_transl_6_sym_varpar: Position vector r_i_i_C has to be [3x1] double');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6PRPRRP1_jacobiaD_transl_6_sym_varpar: pkin has to be [11x1] (double)');

%% Symbolic Calculation
% From jacobiaD_transl_6_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 19:50:20
% EndTime: 2019-02-26 19:50:20
% DurationCPUTime: 0.59s
% Computational Cost: add. (589->92), mult. (1794->167), div. (0->0), fcn. (1946->12), ass. (0->66)
t393 = sin(qJ(4));
t396 = cos(qJ(4));
t392 = sin(qJ(5));
t395 = cos(qJ(5));
t439 = pkin(5) + r_i_i_C(1);
t403 = qJD(5) * (t395 * r_i_i_C(2) + t439 * t392);
t410 = -t392 * r_i_i_C(2) + t439 * t395 + pkin(4);
t436 = r_i_i_C(3) + qJ(6) + pkin(9);
t397 = -(t410 * t393 - t436 * t396) * qJD(4) + t393 * qJD(6) - t396 * t403;
t390 = cos(pkin(6));
t386 = sin(pkin(11));
t394 = sin(qJ(2));
t434 = cos(pkin(11));
t438 = cos(qJ(2));
t402 = t438 * t386 + t394 * t434;
t373 = t402 * t390;
t420 = t438 * t434;
t428 = qJD(2) * t394;
t442 = -qJD(2) * t420 + t386 * t428;
t401 = -t394 * t386 + t420;
t389 = cos(pkin(10));
t387 = sin(pkin(10));
t412 = t389 * t373 + t387 * t401;
t388 = sin(pkin(6));
t432 = t388 * t393;
t348 = -t389 * t432 + t412 * t396;
t398 = t436 * t393 + t410 * t396 + pkin(3);
t435 = pkin(2) * qJD(2);
t431 = t388 * t396;
t430 = t390 * t394;
t427 = qJD(5) * t392;
t426 = qJD(5) * t395;
t370 = t442 * t390;
t375 = t402 * qJD(2);
t353 = t389 * t370 + t387 * t375;
t408 = -t389 * t431 - t393 * t412;
t342 = t408 * qJD(4) - t353 * t396;
t374 = t401 * qJD(2);
t399 = qJD(2) * t373;
t351 = -t387 * t374 - t389 * t399;
t419 = -t342 * t392 - t351 * t395;
t355 = t387 * t370 - t389 * t375;
t411 = -t387 * t373 + t389 * t401;
t407 = t387 * t431 - t393 * t411;
t344 = t407 * qJD(4) + t355 * t396;
t354 = -t389 * t374 + t387 * t399;
t418 = -t344 * t392 - t354 * t395;
t368 = t442 * t388;
t372 = t402 * t388;
t413 = -t372 * t393 + t390 * t396;
t346 = t413 * qJD(4) - t368 * t396;
t369 = qJD(2) * t372;
t417 = -t346 * t392 + t369 * t395;
t400 = t401 * t390;
t358 = -t387 * t402 + t389 * t400;
t416 = -t348 * t395 + t358 * t392;
t350 = t387 * t432 + t396 * t411;
t361 = -t387 * t400 - t389 * t402;
t415 = -t350 * t395 + t361 * t392;
t364 = t372 * t396 + t390 * t393;
t371 = t401 * t388;
t414 = -t364 * t395 + t371 * t392;
t345 = t364 * qJD(4) - t368 * t393;
t343 = t350 * qJD(4) + t355 * t393;
t341 = t348 * qJD(4) - t353 * t393;
t1 = [0 (t355 * t395 - t411 * t427) * r_i_i_C(2) + t355 * pkin(8) + (t387 * t430 - t438 * t389) * t435 + t398 * t354 + t397 * t361 + t439 * (t355 * t392 + t411 * t426) 0, t350 * qJD(6) - t410 * t343 + t436 * t344 - t407 * t403, t418 * r_i_i_C(1) + (-t344 * t395 + t354 * t392) * r_i_i_C(2) + (t415 * r_i_i_C(1) + (t350 * t392 + t361 * t395) * r_i_i_C(2)) * qJD(5) + (t415 * qJD(5) + t418) * pkin(5), t343; 0 (-t353 * t395 - t412 * t427) * r_i_i_C(2) - t353 * pkin(8) + (-t438 * t387 - t389 * t430) * t435 + t398 * t351 + t397 * t358 + t439 * (-t353 * t392 + t412 * t426) 0, t348 * qJD(6) - t410 * t341 + t436 * t342 - t408 * t403, t419 * r_i_i_C(1) + (-t342 * t395 + t351 * t392) * r_i_i_C(2) + (t416 * r_i_i_C(1) + (t348 * t392 + t358 * t395) * r_i_i_C(2)) * qJD(5) + (t416 * qJD(5) + t419) * pkin(5), t341; 0 (-t368 * t395 - t372 * t427) * r_i_i_C(2) - t368 * pkin(8) - t388 * pkin(2) * t428 - t398 * t369 + t397 * t371 + t439 * (-t368 * t392 + t372 * t426) 0, t364 * qJD(6) - t410 * t345 + t436 * t346 - t413 * t403, t417 * r_i_i_C(1) + (-t346 * t395 - t369 * t392) * r_i_i_C(2) + (t414 * r_i_i_C(1) + (t364 * t392 + t371 * t395) * r_i_i_C(2)) * qJD(5) + (t414 * qJD(5) + t417) * pkin(5), t345;];
JaD_transl  = t1;
