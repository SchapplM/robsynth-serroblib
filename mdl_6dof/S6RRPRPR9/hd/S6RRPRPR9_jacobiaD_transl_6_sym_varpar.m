% Zeitableitung der translatorischen Teilmatrix der analytischen Jacobi-Matrix für Segment Nr. 6 (0=Basis) von
% S6RRPRPR9
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
%   pkin=[a2,a3,a4,a5,a6,alpha2,d1,d2,d4,d6,theta3,theta5]';
%
% Output:
% JaD_transl [3x6]
%   Zeitableitung der translatorischen Teilmatrix der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-02-26 21:42
% Revision: d75aae1ac561373cd3be920984c3df29a1c2ecc8 (2019-02-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JaD_transl = S6RRPRPR9_jacobiaD_transl_6_sym_varpar(qJ, qJD, r_i_i_C, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(3,1),zeros(12,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPRPR9_jacobiaD_transl_6_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRPRPR9_jacobiaD_transl_6_sym_varpar: qJD has to be [6x1] (double)');
assert(isa(r_i_i_C,'double') && isreal(r_i_i_C) && all(size(r_i_i_C) == [3 1]), ...
	'S6RRPRPR9_jacobiaD_transl_6_sym_varpar: Position vector r_i_i_C has to be [3x1] double');
assert(isreal(pkin) && all(size(pkin) == [12 1]), ...
  'S6RRPRPR9_jacobiaD_transl_6_sym_varpar: pkin has to be [12x1] (double)');

%% Symbolic Calculation
% From jacobiaD_transl_6_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 21:42:25
% EndTime: 2019-02-26 21:42:25
% DurationCPUTime: 0.49s
% Computational Cost: add. (836->103), mult. (1544->163), div. (0->0), fcn. (1557->14), ass. (0->69)
t401 = sin(qJ(1));
t441 = cos(pkin(6));
t413 = qJD(2) * t441 + qJD(1);
t400 = sin(qJ(2));
t420 = t401 * t441;
t417 = t400 * t420;
t431 = qJD(2) * t400;
t402 = cos(qJ(2));
t403 = cos(qJ(1));
t433 = t403 * t402;
t363 = -qJD(1) * t417 - t401 * t431 + t413 * t433;
t397 = sin(pkin(6));
t435 = t397 * t403;
t449 = -qJD(4) * t435 + t363;
t419 = t403 * t441;
t374 = t400 * t419 + t401 * t402;
t432 = qJD(1) * t397;
t448 = -qJD(4) * t374 + t401 * t432;
t394 = pkin(11) + qJ(4);
t390 = sin(t394);
t392 = cos(t394);
t367 = t374 * t392 - t390 * t435;
t416 = t402 * t419;
t434 = t401 * t400;
t373 = -t416 + t434;
t393 = pkin(12) + qJ(6);
t389 = sin(t393);
t391 = cos(t393);
t447 = t367 * t389 - t373 * t391;
t446 = t367 * t391 + t373 * t389;
t414 = t389 * r_i_i_C(1) + t391 * r_i_i_C(2);
t410 = qJD(6) * t414;
t387 = cos(pkin(12)) * pkin(5) + pkin(4);
t415 = t391 * r_i_i_C(1) - t389 * r_i_i_C(2);
t412 = t387 + t415;
t442 = r_i_i_C(3) + pkin(10) + qJ(5);
t404 = (t412 * t390 - t442 * t392) * qJD(4) - t390 * qJD(5) + t392 * t410;
t357 = t448 * t390 + t449 * t392;
t388 = cos(pkin(11)) * pkin(3) + pkin(2);
t443 = t442 * t390 + t412 * t392 + t388;
t438 = t397 * t400;
t437 = t397 * t401;
t436 = t397 * t402;
t430 = qJD(2) * t402;
t427 = pkin(3) * sin(pkin(11)) + pkin(8);
t425 = t403 * t432;
t424 = t397 * t430;
t422 = t397 * t431;
t421 = -sin(pkin(12)) * pkin(5) - pkin(9) - qJ(3);
t411 = t374 * t390 + t392 * t435;
t376 = -t417 + t433;
t369 = -t376 * t390 + t392 * t437;
t370 = t376 * t392 + t390 * t437;
t375 = t403 * t400 + t402 * t420;
t408 = -t390 * t438 + t441 * t392;
t372 = t441 * t390 + t392 * t438;
t407 = t414 - t421;
t406 = t415 * qJD(6) + qJD(3);
t356 = t449 * t390 - t448 * t392;
t365 = t408 * qJD(4) + t392 * t424;
t364 = t372 * qJD(4) + t390 * t424;
t362 = t375 * qJD(1) + t374 * qJD(2);
t361 = t374 * qJD(1) + t375 * qJD(2);
t360 = -qJD(1) * t416 - t403 * t430 + t413 * t434;
t355 = t369 * qJD(4) - t361 * t392 + t390 * t425;
t354 = t370 * qJD(4) - t361 * t390 - t392 * t425;
t353 = t355 * t391 - t360 * t389 + (-t370 * t389 + t375 * t391) * qJD(6);
t352 = -t355 * t389 - t360 * t391 + (-t370 * t391 - t375 * t389) * qJD(6);
t1 = [-t411 * qJD(5) - t363 * t388 - t373 * qJD(3) - t412 * t357 - t407 * t362 - t442 * t356 + (t447 * r_i_i_C(1) + t446 * r_i_i_C(2)) * qJD(6) + (-t403 * pkin(1) - t427 * t437) * qJD(1), t443 * t360 - t407 * t361 + t404 * t375 + t406 * t376, -t360, t370 * qJD(5) - t412 * t354 + t442 * t355 - t369 * t410, t354, t352 * r_i_i_C(1) - t353 * r_i_i_C(2); t353 * r_i_i_C(1) + t352 * r_i_i_C(2) + t375 * qJD(3) - t369 * qJD(5) + t355 * t387 - t361 * t388 + t421 * t360 + t442 * t354 + (-pkin(1) * t401 + t427 * t435) * qJD(1), -t362 * t443 + t407 * t363 + t404 * t373 + t406 * t374, t362, t367 * qJD(5) - t412 * t356 + t442 * t357 + t411 * t410, t356 (-t357 * t389 + t362 * t391) * r_i_i_C(1) + (-t357 * t391 - t362 * t389) * r_i_i_C(2) + (-t446 * r_i_i_C(1) + t447 * r_i_i_C(2)) * qJD(6); 0 ((-qJD(2) * t443 + t406) * t400 + (t407 * qJD(2) - t404) * t402) * t397, t422, t372 * qJD(5) - t412 * t364 + t442 * t365 - t408 * t410, t364 (-t365 * t389 + t391 * t422) * r_i_i_C(1) + (-t365 * t391 - t389 * t422) * r_i_i_C(2) + ((-t372 * t391 + t389 * t436) * r_i_i_C(1) + (t372 * t389 + t391 * t436) * r_i_i_C(2)) * qJD(6);];
JaD_transl  = t1;
