% Zeitableitung der translatorischen Teilmatrix der analytischen Jacobi-Matrix für Segment Nr. 4 (0=Basis) von
% S6RPRRRR11
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
% pkin [13x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,alpha3,d1,d3,d4,d5,d6,theta2]';
%
% Output:
% JaD_transl [3x6]
%   Zeitableitung der translatorischen Teilmatrix der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-02-26 21:20
% Revision: d75aae1ac561373cd3be920984c3df29a1c2ecc8 (2019-02-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JaD_transl = S6RPRRRR11_jacobiaD_transl_4_sym_varpar(qJ, qJD, r_i_i_C, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(3,1),zeros(13,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRRRR11_jacobiaD_transl_4_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RPRRRR11_jacobiaD_transl_4_sym_varpar: qJD has to be [6x1] (double)');
assert(isa(r_i_i_C,'double') && isreal(r_i_i_C) && all(size(r_i_i_C) == [3 1]), ...
	'S6RPRRRR11_jacobiaD_transl_4_sym_varpar: Position vector r_i_i_C has to be [3x1] double');
assert(isreal(pkin) && all(size(pkin) == [13 1]), ...
  'S6RPRRRR11_jacobiaD_transl_4_sym_varpar: pkin has to be [13x1] (double)');

%% Symbolic Calculation
% From jacobiaD_transl_4_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 21:20:27
% EndTime: 2019-02-26 21:20:28
% DurationCPUTime: 0.44s
% Computational Cost: add. (393->78), mult. (1301->142), div. (0->0), fcn. (1376->12), ass. (0->64)
t406 = sin(qJ(3));
t400 = sin(pkin(7));
t401 = sin(pkin(6));
t410 = cos(qJ(1));
t440 = t401 * t410;
t429 = t400 * t440;
t403 = cos(pkin(7));
t402 = cos(pkin(13));
t404 = cos(pkin(6));
t438 = t404 * t410;
t399 = sin(pkin(13));
t407 = sin(qJ(1));
t443 = t399 * t407;
t419 = t402 * t438 - t443;
t452 = t419 * t403;
t422 = -t452 + t429;
t461 = t422 * t406;
t430 = t404 * t443;
t436 = qJD(1) * t410;
t387 = -qJD(1) * t430 + t402 * t436;
t437 = t407 * t402;
t390 = t399 * t438 + t437;
t409 = cos(qJ(3));
t418 = t399 * t410 + t404 * t437;
t386 = t418 * qJD(1);
t441 = t401 * t407;
t427 = qJD(1) * t441;
t416 = -t386 * t403 + t400 * t427;
t426 = t409 * t429;
t366 = -qJD(3) * t426 + (-qJD(3) * t390 + t416) * t406 - (-qJD(3) * t452 - t387) * t409;
t447 = t386 * t400;
t378 = t403 * t427 + t447;
t405 = sin(qJ(4));
t408 = cos(qJ(4));
t460 = t366 * t405 - t378 * t408;
t459 = -t366 * t408 - t378 * t405;
t371 = -t390 * t409 + t461;
t451 = t419 * t400 + t403 * t440;
t458 = -t371 * t405 + t451 * t408;
t457 = t371 * t408 + t451 * t405;
t439 = t402 * t403;
t442 = t400 * t404;
t450 = (-t399 * t406 + t409 * t439) * t401 + t409 * t442;
t380 = (t399 * t409 + t406 * t439) * t401 + t406 * t442;
t448 = r_i_i_C(3) + pkin(10);
t434 = t401 * qJD(2);
t425 = t405 * r_i_i_C(1) + t408 * r_i_i_C(2);
t423 = t408 * r_i_i_C(1) - t405 * r_i_i_C(2) + pkin(3);
t421 = t400 * t441 - t403 * t418;
t417 = qJD(4) * t425;
t392 = t402 * t410 - t430;
t413 = -t392 * t406 + t421 * t409;
t373 = t392 * t409 + t421 * t406;
t411 = t371 * qJD(3) - t387 * t406 + t416 * t409;
t388 = -t400 * t401 * t402 + t403 * t404;
t385 = t390 * qJD(1);
t383 = t400 * t418 + t403 * t441;
t376 = t451 * qJD(1);
t374 = t450 * qJD(3);
t364 = qJD(1) * t461 + t413 * qJD(3) - t385 * t409;
t363 = t373 * qJD(3) - t385 * t406 + (t409 * t452 - t426) * qJD(1);
t362 = t364 * t408 + t376 * t405 + (-t373 * t405 + t383 * t408) * qJD(4);
t361 = -t364 * t405 + t376 * t408 + (-t373 * t408 - t383 * t405) * qJD(4);
t1 = [t459 * r_i_i_C(1) + t460 * r_i_i_C(2) - t366 * pkin(3) - t387 * pkin(2) - pkin(9) * t447 + t410 * t434 + t448 * t411 + (t458 * r_i_i_C(1) - t457 * r_i_i_C(2)) * qJD(4) + (-t410 * pkin(1) + (-pkin(9) * t403 - qJ(2)) * t441) * qJD(1), t401 * t436, -t423 * t363 + t448 * t364 - t413 * t417, r_i_i_C(1) * t361 - t362 * r_i_i_C(2), 0, 0; t407 * t434 - t385 * pkin(2) + t364 * pkin(3) + t362 * r_i_i_C(1) + t361 * r_i_i_C(2) + t448 * t363 + (-t407 * pkin(1) + pkin(9) * t451 + qJ(2) * t440) * qJD(1), t427, t448 * t366 - (-t390 * t406 - t422 * t409) * t417 + t423 * t411, -t460 * r_i_i_C(1) + t459 * r_i_i_C(2) + (t457 * r_i_i_C(1) + t458 * r_i_i_C(2)) * qJD(4), 0, 0; 0, 0, -t423 * t380 * qJD(3) + t448 * t374 - t450 * t417, -t425 * t374 + ((-t380 * t408 - t388 * t405) * r_i_i_C(1) + (t380 * t405 - t388 * t408) * r_i_i_C(2)) * qJD(4), 0, 0;];
JaD_transl  = t1;
