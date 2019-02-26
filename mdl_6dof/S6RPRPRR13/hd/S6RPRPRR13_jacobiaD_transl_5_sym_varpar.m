% Zeitableitung der translatorischen Teilmatrix der analytischen Jacobi-Matrix für Segment Nr. 5 (0=Basis) von
% S6RPRPRR13
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
%   pkin=[a2,a3,a4,a5,a6,alpha2,alpha3,d1,d3,d5,d6,theta2]';
%
% Output:
% JaD_transl [3x6]
%   Zeitableitung der translatorischen Teilmatrix der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-02-26 20:55
% Revision: d75aae1ac561373cd3be920984c3df29a1c2ecc8 (2019-02-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JaD_transl = S6RPRPRR13_jacobiaD_transl_5_sym_varpar(qJ, qJD, r_i_i_C, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(3,1),zeros(12,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRPRR13_jacobiaD_transl_5_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RPRPRR13_jacobiaD_transl_5_sym_varpar: qJD has to be [6x1] (double)');
assert(isa(r_i_i_C,'double') && isreal(r_i_i_C) && all(size(r_i_i_C) == [3 1]), ...
	'S6RPRPRR13_jacobiaD_transl_5_sym_varpar: Position vector r_i_i_C has to be [3x1] double');
assert(isreal(pkin) && all(size(pkin) == [12 1]), ...
  'S6RPRPRR13_jacobiaD_transl_5_sym_varpar: pkin has to be [12x1] (double)');

%% Symbolic Calculation
% From jacobiaD_transl_5_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 20:55:40
% EndTime: 2019-02-26 20:55:40
% DurationCPUTime: 0.47s
% Computational Cost: add. (485->92), mult. (1598->158), div. (0->0), fcn. (1696->12), ass. (0->74)
t408 = cos(pkin(6));
t403 = sin(pkin(12));
t414 = cos(qJ(1));
t441 = t414 * t403;
t406 = cos(pkin(12));
t411 = sin(qJ(1));
t442 = t411 * t406;
t420 = t408 * t442 + t441;
t388 = t420 * qJD(1);
t443 = t411 * t403;
t432 = t408 * t443;
t439 = qJD(1) * t414;
t389 = -qJD(1) * t432 + t406 * t439;
t440 = t414 * t406;
t391 = -t408 * t440 + t443;
t404 = sin(pkin(7));
t410 = sin(qJ(3));
t413 = cos(qJ(3));
t405 = sin(pkin(6));
t447 = t405 * t411;
t431 = qJD(1) * t447;
t427 = t404 * t431;
t438 = qJD(3) * t410;
t429 = t405 * t438;
t407 = cos(pkin(7));
t444 = t407 * t413;
t445 = t407 * t410;
t392 = t408 * t441 + t442;
t449 = t392 * t413;
t366 = t414 * t404 * t429 + qJD(3) * (t391 * t445 - t449) - t388 * t444 - t389 * t410 + t413 * t427;
t452 = t388 * t404;
t377 = t407 * t431 + t452;
t409 = sin(qJ(5));
t412 = cos(qJ(5));
t462 = t366 * t409 - t377 * t412;
t461 = t366 * t412 + t377 * t409;
t446 = t405 * t414;
t435 = t404 * t446;
t422 = t391 * t407 + t435;
t450 = t392 * t410;
t458 = (t413 * t422 + t450) * qJD(3) - (-t388 * t407 + t427) * t410 - t389 * t413;
t368 = t391 * t444 + t413 * t435 + t450;
t380 = -t391 * t404 + t407 * t446;
t457 = -t368 * t412 - t380 * t409;
t456 = t368 * t409 - t380 * t412;
t386 = t391 * qJD(1);
t453 = t386 * t404;
t448 = t404 * t408;
t437 = t405 * qJD(2);
t436 = r_i_i_C(3) + pkin(10) + pkin(3);
t434 = t406 * t444;
t433 = t413 * t448;
t430 = t405 * t439;
t428 = pkin(9) * t407 + qJ(2);
t426 = t404 * t430;
t425 = t412 * r_i_i_C(1) - t409 * r_i_i_C(2);
t424 = t409 * r_i_i_C(1) + t412 * r_i_i_C(2) + qJ(4);
t421 = t404 * t447 - t407 * t420;
t418 = qJD(5) * t425 + qJD(4);
t394 = -t432 + t440;
t416 = t394 * t413 + t410 * t421;
t415 = t410 * t448 + (t403 * t413 + t406 * t445) * t405;
t390 = -t405 * t406 * t404 + t408 * t407;
t387 = t392 * qJD(1);
t382 = t404 * t420 + t407 * t447;
t378 = -t433 + (t403 * t410 - t434) * t405;
t375 = t407 * t430 - t453;
t374 = t415 * qJD(3);
t371 = t394 * t410 - t413 * t421;
t363 = -t394 * t438 + (t386 * t407 + t426) * t410 + (qJD(3) * t421 - t387) * t413;
t362 = qJD(3) * t416 - t386 * t444 - t387 * t410 - t413 * t426;
t361 = t362 * t409 + t375 * t412 + (t371 * t412 - t382 * t409) * qJD(5);
t360 = t362 * t412 - t375 * t409 + (-t371 * t409 - t382 * t412) * qJD(5);
t1 = [t462 * r_i_i_C(1) + t461 * r_i_i_C(2) - t377 * pkin(4) + t366 * qJ(4) - t368 * qJD(4) - t389 * pkin(2) - pkin(9) * t452 + t414 * t437 + (t457 * r_i_i_C(1) + t456 * r_i_i_C(2)) * qJD(5) + t436 * t458 + (-t414 * pkin(1) - t428 * t447) * qJD(1), t430, -t362 * t436 + t363 * t424 + t416 * t418, t362, t360 * r_i_i_C(1) - t361 * r_i_i_C(2), 0; -pkin(9) * t453 + t411 * t437 - t387 * pkin(2) + t375 * pkin(4) + t361 * r_i_i_C(1) + t360 * r_i_i_C(2) + t362 * qJ(4) + t371 * qJD(4) + t436 * t363 + (-pkin(1) * t411 + t428 * t446) * qJD(1), t431, t418 * (-t410 * t422 + t449) - t424 * t458 + t436 * t366, -t366, -t461 * r_i_i_C(1) + t462 * r_i_i_C(2) + (-t456 * r_i_i_C(1) + t457 * r_i_i_C(2)) * qJD(5), 0; 0, 0, t418 * t415 - t436 * t374 - t424 * (t403 * t429 + (-t405 * t434 - t433) * qJD(3)) t374, t425 * t374 + ((-t378 * t409 - t390 * t412) * r_i_i_C(1) + (-t378 * t412 + t390 * t409) * r_i_i_C(2)) * qJD(5), 0;];
JaD_transl  = t1;
