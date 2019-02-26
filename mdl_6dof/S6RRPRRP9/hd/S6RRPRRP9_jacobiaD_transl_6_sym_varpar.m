% Zeitableitung der translatorischen Teilmatrix der analytischen Jacobi-Matrix für Segment Nr. 6 (0=Basis) von
% S6RRPRRP9
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
%   pkin=[a2,a3,a4,a5,a6,alpha2,d1,d2,d4,d5,theta3]';
%
% Output:
% JaD_transl [3x6]
%   Zeitableitung der translatorischen Teilmatrix der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-02-26 21:50
% Revision: d75aae1ac561373cd3be920984c3df29a1c2ecc8 (2019-02-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JaD_transl = S6RRPRRP9_jacobiaD_transl_6_sym_varpar(qJ, qJD, r_i_i_C, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(3,1),zeros(11,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPRRP9_jacobiaD_transl_6_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRPRRP9_jacobiaD_transl_6_sym_varpar: qJD has to be [6x1] (double)');
assert(isa(r_i_i_C,'double') && isreal(r_i_i_C) && all(size(r_i_i_C) == [3 1]), ...
	'S6RRPRRP9_jacobiaD_transl_6_sym_varpar: Position vector r_i_i_C has to be [3x1] double');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6RRPRRP9_jacobiaD_transl_6_sym_varpar: pkin has to be [11x1] (double)');

%% Symbolic Calculation
% From jacobiaD_transl_6_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 21:50:30
% EndTime: 2019-02-26 21:50:30
% DurationCPUTime: 0.66s
% Computational Cost: add. (853->119), mult. (1759->184), div. (0->0), fcn. (1771->12), ass. (0->73)
t398 = sin(qJ(1));
t400 = cos(qJ(2));
t443 = cos(pkin(6));
t448 = cos(qJ(1));
t418 = t443 * t448;
t397 = sin(qJ(2));
t424 = t398 * t443;
t419 = t397 * t424;
t425 = t448 * qJD(1);
t436 = qJD(2) * t397;
t363 = -qJD(1) * t419 - t398 * t436 + (qJD(2) * t418 + t425) * t400;
t393 = sin(pkin(6));
t430 = t393 * t448;
t458 = -qJD(4) * t430 + t363;
t374 = t397 * t418 + t398 * t400;
t439 = t393 * t398;
t457 = qJD(1) * t439 - qJD(4) * t374;
t391 = pkin(11) + qJ(4);
t389 = sin(t391);
t390 = cos(t391);
t367 = t374 * t390 - t389 * t430;
t414 = t400 * t418;
t437 = t398 * t397;
t373 = -t414 + t437;
t396 = sin(qJ(5));
t399 = cos(qJ(5));
t456 = t367 * t396 - t373 * t399;
t455 = t367 * t399 + t373 * t396;
t449 = pkin(5) + r_i_i_C(1);
t452 = t399 * r_i_i_C(2) + t449 * t396;
t388 = t399 * pkin(5) + pkin(4);
t446 = t396 * r_i_i_C(2);
t413 = t399 * r_i_i_C(1) + t388 - t446;
t444 = r_i_i_C(3) + qJ(6) + pkin(10);
t454 = -(t413 * t389 - t444 * t390) * qJD(4) + t389 * qJD(6);
t453 = pkin(8) + pkin(3) * sin(pkin(11));
t357 = t457 * t389 + t458 * t390;
t387 = cos(pkin(11)) * pkin(3) + pkin(2);
t450 = t444 * t389 + t413 * t390 + t387;
t440 = t393 * t397;
t438 = t393 * t400;
t434 = qJD(5) * t390;
t433 = qJD(5) * t396;
t432 = qJD(5) * t399;
t429 = t448 * t400;
t427 = qJD(2) * t438;
t426 = t393 * t436;
t420 = t393 * t425;
t375 = t448 * t397 + t400 * t424;
t362 = t375 * qJD(1) + t374 * qJD(2);
t417 = -t357 * t396 + t362 * t399;
t372 = t443 * t389 + t390 * t440;
t411 = -t372 * t399 + t396 * t438;
t376 = t429 - t419;
t369 = -t376 * t389 + t390 * t439;
t370 = t376 * t390 + t389 * t439;
t407 = t374 * t389 + t390 * t430;
t404 = -t389 * t440 + t443 * t390;
t365 = t404 * qJD(4) + t390 * t427;
t406 = -t365 * t396 + t399 * t426;
t405 = qJD(5) * t452;
t356 = t458 * t389 - t457 * t390;
t360 = -qJD(1) * t414 - qJD(2) * t429 + (qJD(2) * t443 + qJD(1)) * t437;
t403 = -t360 * t396 + (-t370 * t396 + t375 * t399) * qJD(5);
t361 = t374 * qJD(1) + t375 * qJD(2);
t355 = t369 * qJD(4) - t361 * t390 + t389 * t420;
t352 = -t355 * t396 - t360 * t399 + (-t370 * t399 - t375 * t396) * qJD(5);
t401 = t452 * t434 - t454;
t395 = -pkin(9) - qJ(3);
t364 = t372 * qJD(4) + t389 * t427;
t354 = t370 * qJD(4) - t361 * t389 - t390 * t420;
t353 = t355 * t399 + t403;
t1 = [-t407 * qJD(6) - t363 * t387 - t373 * qJD(3) - t413 * t357 + (t395 - t452) * t362 - t444 * t356 + (-t448 * pkin(1) - t453 * t439) * qJD(1) + (t455 * r_i_i_C(2) + t449 * t456) * qJD(5) (-t361 * t399 - t376 * t433) * r_i_i_C(2) + t361 * t395 + t376 * qJD(3) + t450 * t360 + t401 * t375 + t449 * (-t361 * t396 + t376 * t432) -t360, t370 * qJD(6) - t413 * t354 + t444 * t355 - t369 * t405, -t353 * r_i_i_C(2) + t449 * t352, t354; t353 * r_i_i_C(1) + t352 * r_i_i_C(2) + t375 * qJD(3) - t369 * qJD(6) + t355 * t388 + t360 * t395 - t361 * t387 + t444 * t354 + (-pkin(1) * t398 + t453 * t430) * qJD(1) + t403 * pkin(5) (t363 * t399 - t374 * t433) * r_i_i_C(2) - t363 * t395 + t374 * qJD(3) - t450 * t362 + t401 * t373 + t449 * (t363 * t396 + t374 * t432) t362, t367 * qJD(6) - t413 * t356 + t444 * t357 + t407 * t405, t417 * r_i_i_C(1) + (-t357 * t399 - t362 * t396) * r_i_i_C(2) + (-r_i_i_C(1) * t455 + t456 * r_i_i_C(2)) * qJD(5) + (-qJD(5) * t455 + t417) * pkin(5), t356; 0 ((qJD(3) + (t449 * t399 - t446) * qJD(5) - t450 * qJD(2)) * t397 + (-qJD(2) * t395 + t452 * (qJD(2) - t434) + t454) * t400) * t393, t426, t372 * qJD(6) - t413 * t364 + t444 * t365 - t404 * t405, t406 * r_i_i_C(1) + (-t365 * t399 - t396 * t426) * r_i_i_C(2) + (t411 * r_i_i_C(1) + (t372 * t396 + t399 * t438) * r_i_i_C(2)) * qJD(5) + (qJD(5) * t411 + t406) * pkin(5), t364;];
JaD_transl  = t1;
