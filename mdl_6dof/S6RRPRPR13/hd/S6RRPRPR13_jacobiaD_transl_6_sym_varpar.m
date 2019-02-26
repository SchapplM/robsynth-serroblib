% Zeitableitung der translatorischen Teilmatrix der analytischen Jacobi-Matrix für Segment Nr. 6 (0=Basis) von
% S6RRPRPR13
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
%   pkin=[a2,a3,a4,a5,a6,alpha2,d1,d2,d4,d6,theta5]';
%
% Output:
% JaD_transl [3x6]
%   Zeitableitung der translatorischen Teilmatrix der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-02-26 21:45
% Revision: d75aae1ac561373cd3be920984c3df29a1c2ecc8 (2019-02-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JaD_transl = S6RRPRPR13_jacobiaD_transl_6_sym_varpar(qJ, qJD, r_i_i_C, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(3,1),zeros(11,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPRPR13_jacobiaD_transl_6_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRPRPR13_jacobiaD_transl_6_sym_varpar: qJD has to be [6x1] (double)');
assert(isa(r_i_i_C,'double') && isreal(r_i_i_C) && all(size(r_i_i_C) == [3 1]), ...
	'S6RRPRPR13_jacobiaD_transl_6_sym_varpar: Position vector r_i_i_C has to be [3x1] double');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6RRPRPR13_jacobiaD_transl_6_sym_varpar: pkin has to be [11x1] (double)');

%% Symbolic Calculation
% From jacobiaD_transl_6_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 21:44:52
% EndTime: 2019-02-26 21:44:52
% DurationCPUTime: 0.44s
% Computational Cost: add. (634->101), mult. (1584->158), div. (0->0), fcn. (1592->12), ass. (0->64)
t396 = cos(pkin(6));
t398 = sin(qJ(4));
t401 = cos(qJ(4));
t395 = sin(pkin(6));
t402 = cos(qJ(2));
t433 = t395 * t402;
t378 = t396 * t401 - t398 * t433;
t393 = pkin(11) + qJ(6);
t391 = sin(t393);
t392 = cos(t393);
t416 = t392 * r_i_i_C(1) - t391 * r_i_i_C(2);
t409 = t416 * qJD(6);
t399 = sin(qJ(2));
t403 = cos(qJ(1));
t428 = t403 * t399;
t400 = sin(qJ(1));
t429 = t400 * t402;
t380 = t396 * t428 + t429;
t427 = t403 * t402;
t430 = t400 * t399;
t379 = -t396 * t427 + t430;
t432 = t395 * t403;
t422 = t401 * t432;
t413 = -t379 * t398 + t422;
t441 = -t380 * t392 - t391 * t413;
t440 = -t380 * t391 + t392 * t413;
t381 = t396 * t429 + t428;
t367 = t381 * qJD(1) + t380 * qJD(2);
t412 = t379 * t401 + t398 * t432;
t426 = qJD(1) * t395;
t421 = t400 * t426;
t359 = t412 * qJD(4) + t367 * t398 + t401 * t421;
t390 = cos(pkin(11)) * pkin(5) + pkin(4);
t414 = t390 + t416;
t438 = r_i_i_C(3) + pkin(10) + qJ(5);
t405 = t414 * t398 - t438 * t401 + qJ(3);
t439 = pkin(3) + pkin(8);
t435 = t395 * t399;
t434 = t395 * t400;
t425 = qJD(2) * t399;
t423 = t396 * t430;
t420 = t403 * t426;
t419 = qJD(2) * t433;
t418 = t395 * t425;
t417 = -sin(pkin(11)) * pkin(5) - pkin(2) - pkin(9);
t415 = -t391 * r_i_i_C(1) - t392 * r_i_i_C(2);
t372 = t381 * t398 + t401 * t434;
t371 = t381 * t401 - t398 * t434;
t411 = t396 * t398 + t401 * t433;
t410 = t423 - t427;
t408 = qJD(6) * t415;
t406 = -t415 - t417;
t357 = -t367 * t401 - qJD(4) * t422 + (qJD(4) * t379 + t421) * t398;
t404 = -t401 * qJD(5) + qJD(3) + t398 * t408 + (t438 * t398 + t414 * t401) * qJD(4);
t370 = t378 * qJD(4) - t401 * t418;
t369 = t411 * qJD(4) - t398 * t418;
t368 = -qJD(1) * t423 - t400 * t425 + (qJD(2) * t396 + qJD(1)) * t427;
t366 = t380 * qJD(1) + t381 * qJD(2);
t365 = t379 * qJD(1) + t410 * qJD(2);
t362 = t371 * qJD(4) - t365 * t398 + t401 * t420;
t361 = t372 * qJD(4) + t365 * t401 + t398 * t420;
t356 = t362 * t392 - t366 * t391 + (-t372 * t391 - t392 * t410) * qJD(6);
t355 = -t362 * t391 - t366 * t392 + (-t372 * t392 + t391 * t410) * qJD(6);
t1 = [t412 * qJD(5) - t367 * qJ(3) - t379 * qJD(3) - t438 * t357 - t414 * t359 + (t441 * r_i_i_C(1) - t440 * r_i_i_C(2)) * qJD(6) - t406 * t368 + (-t403 * pkin(1) - t439 * t434) * qJD(1), t406 * t365 - t405 * t366 - t381 * t409 - t404 * t410, -t365, t372 * qJD(5) - t414 * t361 + t438 * t362 + t371 * t408, t361, t355 * r_i_i_C(1) - t356 * r_i_i_C(2); t356 * r_i_i_C(1) + t355 * r_i_i_C(2) - t365 * qJ(3) + t381 * qJD(3) - t371 * qJD(5) + t362 * t390 + t438 * t361 + t417 * t366 + (-pkin(1) * t400 + t439 * t432) * qJD(1), -t367 * t406 + t368 * t405 - t379 * t409 + t380 * t404, t367, -qJD(5) * t413 - t414 * t357 + t438 * t359 + t412 * t408, t357 (-t359 * t391 + t368 * t392) * r_i_i_C(1) + (-t359 * t392 - t368 * t391) * r_i_i_C(2) + (t440 * r_i_i_C(1) + t441 * r_i_i_C(2)) * qJD(6); 0 ((qJD(2) * t405 + t409) * t402 + (-qJD(2) * t406 + t404) * t399) * t395, t418, t378 * qJD(5) - t438 * t369 - t414 * t370 - t411 * t408, t370 (t369 * t391 + t392 * t419) * r_i_i_C(1) + (t369 * t392 - t391 * t419) * r_i_i_C(2) + ((-t378 * t392 - t391 * t435) * r_i_i_C(1) + (t378 * t391 - t392 * t435) * r_i_i_C(2)) * qJD(6);];
JaD_transl  = t1;
