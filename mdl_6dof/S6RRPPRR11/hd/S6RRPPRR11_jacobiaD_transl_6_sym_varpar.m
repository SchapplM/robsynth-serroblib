% Zeitableitung der translatorischen Teilmatrix der analytischen Jacobi-Matrix für Segment Nr. 6 (0=Basis) von
% S6RRPPRR11
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
%   pkin=[a2,a3,a4,a5,a6,alpha2,d1,d2,d5,d6,theta4]';
%
% Output:
% JaD_transl [3x6]
%   Zeitableitung der translatorischen Teilmatrix der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-02-26 21:34
% Revision: d75aae1ac561373cd3be920984c3df29a1c2ecc8 (2019-02-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JaD_transl = S6RRPPRR11_jacobiaD_transl_6_sym_varpar(qJ, qJD, r_i_i_C, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(3,1),zeros(11,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPPRR11_jacobiaD_transl_6_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRPPRR11_jacobiaD_transl_6_sym_varpar: qJD has to be [6x1] (double)');
assert(isa(r_i_i_C,'double') && isreal(r_i_i_C) && all(size(r_i_i_C) == [3 1]), ...
	'S6RRPPRR11_jacobiaD_transl_6_sym_varpar: Position vector r_i_i_C has to be [3x1] double');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6RRPPRR11_jacobiaD_transl_6_sym_varpar: pkin has to be [11x1] (double)');

%% Symbolic Calculation
% From jacobiaD_transl_6_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 21:34:21
% EndTime: 2019-02-26 21:34:22
% DurationCPUTime: 0.47s
% Computational Cost: add. (702->100), mult. (1518->156), div. (0->0), fcn. (1517->12), ass. (0->66)
t402 = cos(pkin(6));
t408 = cos(qJ(2));
t409 = cos(qJ(1));
t437 = t408 * t409;
t431 = t402 * t437;
t405 = sin(qJ(2));
t406 = sin(qJ(1));
t440 = t405 * t406;
t384 = -t431 + t440;
t399 = pkin(11) + qJ(5);
t397 = sin(t399);
t398 = cos(t399);
t401 = sin(pkin(6));
t441 = t401 * t409;
t433 = t398 * t441;
t380 = -t384 * t397 + t433;
t438 = t406 * t408;
t439 = t405 * t409;
t385 = t402 * t439 + t438;
t404 = sin(qJ(6));
t407 = cos(qJ(6));
t452 = -t380 * t404 - t385 * t407;
t451 = t380 * t407 - t385 * t404;
t386 = t402 * t438 + t439;
t373 = t386 * qJD(1) + t385 * qJD(2);
t436 = qJD(1) * t401;
t430 = t406 * t436;
t450 = (qJD(5) * t384 + t430) * t397 - qJD(5) * t433 - t373 * t398;
t420 = t384 * t398 + t397 * t441;
t366 = t420 * qJD(5) + t373 * t397 + t398 * t430;
t423 = r_i_i_C(1) * t407 - r_i_i_C(2) * t404;
t412 = -t423 * qJD(6) - qJD(4);
t421 = pkin(5) + t423;
t426 = -sin(pkin(11)) * pkin(4) - qJ(3);
t449 = -r_i_i_C(3) - pkin(10);
t411 = t421 * t397 + t449 * t398 - t426;
t448 = -pkin(2) - pkin(9) - qJ(4);
t447 = pkin(8) + cos(pkin(11)) * pkin(4) + pkin(3);
t444 = t401 * t405;
t443 = t401 * t406;
t442 = t401 * t408;
t435 = qJD(2) * t405;
t434 = qJD(2) * t408;
t432 = t402 * t440;
t429 = t409 * t436;
t428 = t401 * t435;
t427 = t401 * t434;
t424 = qJD(2) * t402 + qJD(1);
t422 = -t404 * r_i_i_C(1) - t407 * r_i_i_C(2);
t378 = t386 * t397 + t398 * t443;
t419 = t386 * t398 - t397 * t443;
t418 = t397 * t402 + t398 * t442;
t417 = t397 * t442 - t398 * t402;
t416 = qJD(6) * t422;
t415 = -t422 - t448;
t410 = qJD(3) + t397 * t416 + (-t449 * t397 + t421 * t398) * qJD(5);
t387 = -t432 + t437;
t375 = t418 * qJD(5) - t397 * t428;
t374 = -qJD(1) * t432 - t406 * t435 + t424 * t437;
t372 = t385 * qJD(1) + t386 * qJD(2);
t371 = -qJD(1) * t431 - t409 * t434 + t424 * t440;
t369 = t419 * qJD(5) - t371 * t397 + t398 * t429;
t368 = t378 * qJD(5) + t371 * t398 + t397 * t429;
t363 = t369 * t407 - t372 * t404 + (-t378 * t404 + t387 * t407) * qJD(6);
t362 = -t369 * t404 - t372 * t407 + (-t378 * t407 - t387 * t404) * qJD(6);
t1 = [-t384 * qJD(3) - t385 * qJD(4) - t421 * t366 - t415 * t374 + t426 * t373 + t449 * t450 + (t452 * r_i_i_C(1) - t451 * r_i_i_C(2)) * qJD(6) + (-t409 * pkin(1) - t447 * t443) * qJD(1), t415 * t371 - t411 * t372 + t412 * t386 + t410 * t387, -t371, -t372, -t421 * t368 - t449 * t369 + t419 * t416, r_i_i_C(1) * t362 - t363 * r_i_i_C(2); t369 * pkin(5) + t363 * r_i_i_C(1) + t362 * r_i_i_C(2) + t386 * qJD(3) + t387 * qJD(4) + t448 * t372 + t426 * t371 - t449 * t368 + (-pkin(1) * t406 + t447 * t441) * qJD(1), -t415 * t373 + t411 * t374 + t412 * t384 + t410 * t385, t373, t374, -t449 * t366 + t420 * t416 - t421 * t450 (-t366 * t404 + t374 * t407) * r_i_i_C(1) + (-t366 * t407 - t374 * t404) * r_i_i_C(2) + (t451 * r_i_i_C(1) + t452 * r_i_i_C(2)) * qJD(6); 0 ((t411 * qJD(2) - t412) * t408 + (-t415 * qJD(2) + t410) * t405) * t401, t428, t427, t449 * t375 - t418 * t416 + t421 * (t417 * qJD(5) + t398 * t428) (t375 * t404 + t407 * t427) * r_i_i_C(1) + (t375 * t407 - t404 * t427) * r_i_i_C(2) + ((-t404 * t444 + t407 * t417) * r_i_i_C(1) + (-t404 * t417 - t407 * t444) * r_i_i_C(2)) * qJD(6);];
JaD_transl  = t1;
