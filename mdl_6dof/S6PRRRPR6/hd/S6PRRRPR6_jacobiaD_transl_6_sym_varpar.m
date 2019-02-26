% Zeitableitung der translatorischen Teilmatrix der analytischen Jacobi-Matrix für Segment Nr. 6 (0=Basis) von
% S6PRRRPR6
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
% Datum: 2019-02-26 20:13
% Revision: d75aae1ac561373cd3be920984c3df29a1c2ecc8 (2019-02-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JaD_transl = S6PRRRPR6_jacobiaD_transl_6_sym_varpar(qJ, qJD, r_i_i_C, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(3,1),zeros(11,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRRRPR6_jacobiaD_transl_6_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6PRRRPR6_jacobiaD_transl_6_sym_varpar: qJD has to be [6x1] (double)');
assert(isa(r_i_i_C,'double') && isreal(r_i_i_C) && all(size(r_i_i_C) == [3 1]), ...
	'S6PRRRPR6_jacobiaD_transl_6_sym_varpar: Position vector r_i_i_C has to be [3x1] double');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6PRRRPR6_jacobiaD_transl_6_sym_varpar: pkin has to be [11x1] (double)');

%% Symbolic Calculation
% From jacobiaD_transl_6_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 20:13:37
% EndTime: 2019-02-26 20:13:37
% DurationCPUTime: 0.86s
% Computational Cost: add. (988->143), mult. (3026->250), div. (0->0), fcn. (3289->12), ass. (0->85)
t406 = sin(pkin(11));
t411 = sin(qJ(2));
t415 = cos(qJ(2));
t461 = cos(pkin(11));
t462 = cos(pkin(6));
t432 = t462 * t461;
t424 = t415 * t432;
t391 = t406 * t411 - t424;
t409 = sin(qJ(4));
t413 = cos(qJ(4));
t392 = t406 * t415 + t411 * t432;
t410 = sin(qJ(3));
t414 = cos(qJ(3));
t407 = sin(pkin(6));
t438 = t407 * t461;
t421 = -t392 * t414 + t410 * t438;
t363 = -t391 * t413 - t409 * t421;
t364 = t391 * t409 - t413 * t421;
t408 = sin(qJ(6));
t412 = cos(qJ(6));
t470 = ((t363 * t408 + t364 * t412) * r_i_i_C(1) + (t363 * t412 - t364 * t408) * r_i_i_C(2)) * qJD(6);
t439 = t406 * t462;
t394 = -t411 * t439 + t461 * t415;
t458 = t407 * t410;
t374 = t394 * t414 + t406 * t458;
t393 = t461 * t411 + t415 * t439;
t365 = t374 * t409 - t393 * t413;
t366 = t374 * t413 + t393 * t409;
t469 = ((t365 * t408 + t366 * t412) * r_i_i_C(1) + (t365 * t412 - t366 * t408) * r_i_i_C(2)) * qJD(6);
t454 = t411 * t414;
t396 = t407 * t454 + t462 * t410;
t452 = t413 * t415;
t377 = t396 * t409 + t407 * t452;
t455 = t409 * t415;
t460 = t396 * t413;
t378 = -t407 * t455 + t460;
t468 = ((t377 * t408 + t378 * t412) * r_i_i_C(1) + (t377 * t412 - t378 * t408) * r_i_i_C(2)) * qJD(6);
t445 = -r_i_i_C(3) - pkin(10) + pkin(9);
t467 = -pkin(3) * t410 + t445 * t414;
t459 = t406 * t407;
t457 = t409 * t411;
t456 = t409 * t414;
t453 = t413 * t414;
t451 = t414 * t415;
t450 = qJD(2) * t411;
t449 = qJD(3) * t410;
t448 = qJD(3) * t415;
t447 = qJD(4) * t407;
t446 = qJD(4) * t414;
t444 = t411 * t458;
t443 = t407 * t450;
t442 = qJD(2) * t407 * t415;
t441 = t410 * t448;
t440 = t415 * t447;
t435 = t414 * t438;
t387 = -qJD(2) * t424 + t406 * t450;
t434 = t391 * t446 - t387;
t389 = t393 * qJD(2);
t433 = t393 * t446 - t389;
t425 = t408 * r_i_i_C(1) + t412 * r_i_i_C(2) + qJ(5);
t423 = t412 * r_i_i_C(1) - t408 * r_i_i_C(2) + pkin(4) + pkin(5);
t422 = -t414 * pkin(3) - t445 * t410 - pkin(2);
t420 = qJD(3) * t467;
t388 = t392 * qJD(2);
t419 = qJD(4) * t392 - t388 * t414 + t391 * t449;
t390 = t394 * qJD(2);
t418 = qJD(4) * t394 - t390 * t414 + t393 * t449;
t417 = t425 * t409 + t423 * t413 + pkin(3);
t416 = t409 * qJD(5) + ((-t408 * t413 + t409 * t412) * r_i_i_C(1) + (-t408 * t409 - t412 * t413) * r_i_i_C(2)) * qJD(6) + (-t423 * t409 + t425 * t413) * qJD(4);
t382 = (t413 * t451 + t457) * t407;
t381 = (t409 * t451 - t411 * t413) * t407;
t376 = -qJD(3) * t444 + (t462 * qJD(3) + t442) * t414;
t370 = -t393 * t453 + t394 * t409;
t369 = -t393 * t456 - t394 * t413;
t368 = -t391 * t453 + t392 * t409;
t367 = -t391 * t456 - t392 * t413;
t362 = -t394 * t449 + (qJD(3) * t459 - t389) * t414;
t360 = -qJD(3) * t435 - t387 * t414 - t392 * t449;
t356 = -t377 * qJD(4) + t376 * t413 + t409 * t443;
t355 = qJD(4) * t460 - t413 * t443 + (t376 - t440) * t409;
t350 = -t365 * qJD(4) + t362 * t413 + t390 * t409;
t349 = t366 * qJD(4) + t362 * t409 - t390 * t413;
t348 = -t363 * qJD(4) + t360 * t413 + t388 * t409;
t347 = t364 * qJD(4) + t360 * t409 - t388 * t413;
t1 = [0, -t389 * pkin(8) + t369 * qJD(5) + t425 * (t418 * t409 - t433 * t413) + t423 * (t433 * t409 + t418 * t413) + ((t369 * t412 - t370 * t408) * r_i_i_C(1) + (-t369 * t408 - t370 * t412) * r_i_i_C(2)) * qJD(6) - t393 * t420 + t422 * t390, t445 * t362 + t417 * (-t374 * qJD(3) + t389 * t410) + t416 * (-t394 * t410 + t414 * t459) t366 * qJD(5) - t423 * t349 + t425 * t350 + t469, t349 (t349 * t412 - t350 * t408) * r_i_i_C(1) + (-t349 * t408 - t350 * t412) * r_i_i_C(2) - t469; 0, -t387 * pkin(8) + t367 * qJD(5) + t425 * (t419 * t409 - t434 * t413) + t423 * (t434 * t409 + t419 * t413) + ((t367 * t412 - t368 * t408) * r_i_i_C(1) + (-t367 * t408 - t368 * t412) * r_i_i_C(2)) * qJD(6) - t391 * t420 + t422 * t388, t445 * t360 + t417 * (t421 * qJD(3) + t387 * t410) + t416 * (-t392 * t410 - t435) t364 * qJD(5) - t423 * t347 + t425 * t348 + t470, t347 (t347 * t412 - t348 * t408) * r_i_i_C(1) + (-t347 * t408 - t348 * t412) * r_i_i_C(2) - t470; 0, t381 * qJD(5) - t425 * (-t440 * t453 - t447 * t457) + ((t381 * t412 - t382 * t408) * r_i_i_C(1) + (-t381 * t408 - t382 * t412) * r_i_i_C(2)) * qJD(6) + (-t425 * (t409 * t441 + (t409 * t454 + t452) * qJD(2)) + t423 * ((qJD(2) - t446) * t455 + (-t441 + (-qJD(2) * t414 + qJD(4)) * t411) * t413) + t467 * t448 + (pkin(8) * t415 + t422 * t411) * qJD(2)) * t407, t445 * t376 + t417 * (-t396 * qJD(3) - t410 * t442) + t416 * (t462 * t414 - t444) t378 * qJD(5) - t423 * t355 + t425 * t356 + t468, t355 (t355 * t412 - t356 * t408) * r_i_i_C(1) + (-t355 * t408 - t356 * t412) * r_i_i_C(2) - t468;];
JaD_transl  = t1;
