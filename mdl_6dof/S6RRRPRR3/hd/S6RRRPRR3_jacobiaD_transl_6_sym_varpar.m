% Zeitableitung der translatorischen Teilmatrix der analytischen Jacobi-Matrix für Segment Nr. 6 (0=Basis) von
% S6RRRPRR3
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
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d3,d5,d6]';
%
% Output:
% JaD_transl [3x6]
%   Zeitableitung der translatorischen Teilmatrix der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-02-26 22:17
% Revision: d75aae1ac561373cd3be920984c3df29a1c2ecc8 (2019-02-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JaD_transl = S6RRRPRR3_jacobiaD_transl_6_sym_varpar(qJ, qJD, r_i_i_C, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(3,1),zeros(10,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRPRR3_jacobiaD_transl_6_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRRPRR3_jacobiaD_transl_6_sym_varpar: qJD has to be [6x1] (double)');
assert(isa(r_i_i_C,'double') && isreal(r_i_i_C) && all(size(r_i_i_C) == [3 1]), ...
	'S6RRRPRR3_jacobiaD_transl_6_sym_varpar: Position vector r_i_i_C has to be [3x1] double');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RRRPRR3_jacobiaD_transl_6_sym_varpar: pkin has to be [10x1] (double)');

%% Symbolic Calculation
% From jacobiaD_transl_6_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 22:17:09
% EndTime: 2019-02-26 22:17:10
% DurationCPUTime: 0.60s
% Computational Cost: add. (1051->98), mult. (1398->150), div. (0->0), fcn. (1290->10), ass. (0->78)
t406 = qJ(2) + qJ(3);
t403 = sin(t406);
t405 = qJD(2) + qJD(3);
t404 = cos(t406);
t458 = qJ(4) * t404;
t468 = pkin(3) + pkin(4);
t480 = -t405 * t458 - (-t468 * t405 + qJD(4)) * t403;
t408 = sin(qJ(2));
t460 = pkin(2) * qJD(2);
t450 = t408 * t460;
t479 = t450 + t480;
t409 = sin(qJ(1));
t412 = cos(qJ(1));
t465 = sin(qJ(5));
t444 = t404 * t465;
t429 = qJD(5) * t444;
t466 = cos(qJ(5));
t443 = t405 * t466;
t418 = t403 * t443 + t429;
t376 = t403 * t465 + t404 * t466;
t421 = qJD(1) * t376;
t437 = qJD(5) * t466;
t430 = t403 * t437;
t439 = t412 * t465;
t432 = t404 * t439;
t363 = -t405 * t432 + t409 * t421 + (t418 - t430) * t412;
t445 = t403 * t466;
t431 = qJD(1) * t445;
t440 = t412 * t466;
t441 = t409 * t465;
t364 = -t409 * t431 + (t405 - qJD(5)) * t403 * t439 + (qJD(1) * t441 + t405 * t440 - t412 * t437) * t404;
t467 = -r_i_i_C(3) - pkin(10);
t407 = sin(qJ(6));
t410 = cos(qJ(6));
t469 = -t410 * r_i_i_C(1) + t407 * r_i_i_C(2) - pkin(5);
t478 = t467 * t363 - t469 * t364;
t442 = t405 * t465;
t370 = t376 * qJD(5) - t403 * t442 - t404 * t443;
t365 = qJD(1) * t432 + t370 * t409 - t412 * t431;
t419 = t404 * t442 + t430;
t434 = t409 * t445;
t366 = -t405 * t434 + t412 * t421 + (t419 - t429) * t409;
t477 = t469 * t365 - t467 * t366;
t476 = t469 * (-t418 + t419) + t467 * t370;
t451 = qJD(6) * t410;
t452 = qJD(6) * t407;
t471 = -r_i_i_C(1) * t452 - t451 * r_i_i_C(2);
t464 = pkin(2) * t408;
t462 = t410 * r_i_i_C(2);
t461 = -pkin(9) + pkin(8) + pkin(7);
t459 = qJ(4) * t403;
t457 = t405 * t409;
t456 = t405 * t412;
t455 = qJD(1) * t409;
t454 = qJD(1) * t412;
t453 = qJD(4) * t404;
t447 = t468 * t412;
t446 = t404 * t457;
t438 = t403 * t455;
t373 = t376 * t409;
t428 = t373 * t410 + t407 * t412;
t427 = t373 * t407 - t410 * t412;
t425 = -t468 * t404 - t459;
t424 = qJD(6) * (t407 * r_i_i_C(1) + t462);
t411 = cos(qJ(2));
t422 = qJD(1) * (-t411 * pkin(2) - pkin(1) + t425);
t417 = t405 * t425 - t411 * t460;
t372 = t404 * t441 - t434;
t416 = t471 * t372 + t409 * t453 + t454 * t458 - t477;
t374 = -t403 * t440 + t432;
t415 = t471 * t374 + t412 * t453 + t468 * t438 - t478;
t414 = t471 * t376 - t476 - t480;
t400 = t407 * t455;
t377 = t445 - t444;
t375 = t376 * t412;
t356 = -t407 * t454 - t363 * t410 + (-t375 * t407 - t409 * t410) * qJD(6);
t355 = -t410 * t454 + t363 * t407 + (-t375 * t410 + t407 * t409) * qJD(6);
t1 = [t400 * r_i_i_C(1) + t469 * t366 + t467 * t365 + (r_i_i_C(1) * t427 + r_i_i_C(2) * t428) * qJD(6) + t412 * t422 + ((-t461 + t462) * qJD(1) + t479) * t409 (-t458 + t464) * t455 + t417 * t412 + t415, -t456 * t459 + (-qJ(4) * t455 - t405 * t447) * t404 + t415, t404 * t456 - t438, t374 * t424 + t478, t355 * r_i_i_C(1) - t356 * r_i_i_C(2); -t363 * pkin(5) + t356 * r_i_i_C(1) + t355 * r_i_i_C(2) + t467 * t364 + t409 * t422 + (t461 * qJD(1) - t479) * t412 (-t468 * t403 - t464) * t454 + t417 * t409 + t416, -t468 * t446 + (-qJ(4) * t457 - qJD(1) * t447) * t403 + t416, t403 * t454 + t446, t372 * t424 + t477 (-t366 * t407 - t410 * t455) * r_i_i_C(1) + (-t366 * t410 + t400) * r_i_i_C(2) + (-r_i_i_C(1) * t428 + r_i_i_C(2) * t427) * qJD(6); 0, t414 - t450, t414, t405 * t403, t376 * t424 + t476 (t370 * t410 + t377 * t452) * r_i_i_C(2) + (t370 * t407 - t377 * t451) * r_i_i_C(1);];
JaD_transl  = t1;
