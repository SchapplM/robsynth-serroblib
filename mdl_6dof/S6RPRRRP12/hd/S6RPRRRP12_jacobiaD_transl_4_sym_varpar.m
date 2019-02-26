% Zeitableitung der translatorischen Teilmatrix der analytischen Jacobi-Matrix für Segment Nr. 4 (0=Basis) von
% S6RPRRRP12
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
%   pkin=[a2,a3,a4,a5,a6,alpha2,alpha3,d1,d3,d4,d5,theta2]';
%
% Output:
% JaD_transl [3x6]
%   Zeitableitung der translatorischen Teilmatrix der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-02-26 21:14
% Revision: d75aae1ac561373cd3be920984c3df29a1c2ecc8 (2019-02-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JaD_transl = S6RPRRRP12_jacobiaD_transl_4_sym_varpar(qJ, qJD, r_i_i_C, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(3,1),zeros(12,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRRRP12_jacobiaD_transl_4_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RPRRRP12_jacobiaD_transl_4_sym_varpar: qJD has to be [6x1] (double)');
assert(isa(r_i_i_C,'double') && isreal(r_i_i_C) && all(size(r_i_i_C) == [3 1]), ...
	'S6RPRRRP12_jacobiaD_transl_4_sym_varpar: Position vector r_i_i_C has to be [3x1] double');
assert(isreal(pkin) && all(size(pkin) == [12 1]), ...
  'S6RPRRRP12_jacobiaD_transl_4_sym_varpar: pkin has to be [12x1] (double)');

%% Symbolic Calculation
% From jacobiaD_transl_4_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 21:14:24
% EndTime: 2019-02-26 21:14:25
% DurationCPUTime: 0.44s
% Computational Cost: add. (393->78), mult. (1301->141), div. (0->0), fcn. (1376->12), ass. (0->65)
t408 = sin(qJ(3));
t402 = sin(pkin(7));
t403 = sin(pkin(6));
t412 = cos(qJ(1));
t444 = t403 * t412;
t431 = t402 * t444;
t405 = cos(pkin(7));
t406 = cos(pkin(6));
t404 = cos(pkin(12));
t439 = t412 * t404;
t401 = sin(pkin(12));
t409 = sin(qJ(1));
t442 = t409 * t401;
t421 = t406 * t439 - t442;
t455 = t421 * t405;
t424 = -t455 + t431;
t464 = t424 * t408;
t432 = t406 * t442;
t438 = qJD(1) * t412;
t389 = -qJD(1) * t432 + t404 * t438;
t440 = t412 * t401;
t441 = t409 * t404;
t392 = t406 * t440 + t441;
t411 = cos(qJ(3));
t420 = t406 * t441 + t440;
t388 = t420 * qJD(1);
t445 = t403 * t409;
t429 = qJD(1) * t445;
t418 = -t388 * t405 + t402 * t429;
t428 = t411 * t431;
t368 = -qJD(3) * t428 + (-qJD(3) * t392 + t418) * t408 - (-qJD(3) * t455 - t389) * t411;
t450 = t388 * t402;
t380 = t405 * t429 + t450;
t407 = sin(qJ(4));
t410 = cos(qJ(4));
t463 = t368 * t407 - t380 * t410;
t462 = -t368 * t410 - t380 * t407;
t373 = -t392 * t411 + t464;
t454 = t421 * t402 + t405 * t444;
t461 = -t373 * t407 + t454 * t410;
t460 = t373 * t410 + t454 * t407;
t443 = t404 * t405;
t446 = t402 * t406;
t453 = t403 * (-t401 * t408 + t411 * t443) + t411 * t446;
t382 = (t401 * t411 + t408 * t443) * t403 + t408 * t446;
t451 = r_i_i_C(3) + pkin(10);
t436 = t403 * qJD(2);
t427 = t407 * r_i_i_C(1) + t410 * r_i_i_C(2);
t425 = t410 * r_i_i_C(1) - t407 * r_i_i_C(2) + pkin(3);
t423 = t402 * t445 - t405 * t420;
t419 = qJD(4) * t427;
t394 = -t432 + t439;
t415 = -t394 * t408 + t423 * t411;
t375 = t394 * t411 + t423 * t408;
t413 = t373 * qJD(3) - t389 * t408 + t418 * t411;
t390 = -t403 * t404 * t402 + t406 * t405;
t387 = t392 * qJD(1);
t385 = t402 * t420 + t405 * t445;
t378 = t454 * qJD(1);
t376 = t453 * qJD(3);
t366 = qJD(1) * t464 + t415 * qJD(3) - t387 * t411;
t365 = t375 * qJD(3) - t387 * t408 + (t411 * t455 - t428) * qJD(1);
t364 = t366 * t410 + t378 * t407 + (-t375 * t407 + t385 * t410) * qJD(4);
t363 = -t366 * t407 + t378 * t410 + (-t375 * t410 - t385 * t407) * qJD(4);
t1 = [t462 * r_i_i_C(1) + t463 * r_i_i_C(2) - t368 * pkin(3) - t389 * pkin(2) - pkin(9) * t450 + t412 * t436 + t451 * t413 + (t461 * r_i_i_C(1) - t460 * r_i_i_C(2)) * qJD(4) + (-t412 * pkin(1) + (-pkin(9) * t405 - qJ(2)) * t445) * qJD(1), t403 * t438, -t425 * t365 + t451 * t366 - t415 * t419, t363 * r_i_i_C(1) - t364 * r_i_i_C(2), 0, 0; t409 * t436 - t387 * pkin(2) + t366 * pkin(3) + t364 * r_i_i_C(1) + t363 * r_i_i_C(2) + t451 * t365 + (-t409 * pkin(1) + pkin(9) * t454 + qJ(2) * t444) * qJD(1), t429, t451 * t368 - (-t392 * t408 - t424 * t411) * t419 + t425 * t413, -t463 * r_i_i_C(1) + t462 * r_i_i_C(2) + (t460 * r_i_i_C(1) + t461 * r_i_i_C(2)) * qJD(4), 0, 0; 0, 0, -t425 * t382 * qJD(3) + t451 * t376 - t453 * t419, -t427 * t376 + ((-t382 * t410 - t390 * t407) * r_i_i_C(1) + (t382 * t407 - t390 * t410) * r_i_i_C(2)) * qJD(4), 0, 0;];
JaD_transl  = t1;
