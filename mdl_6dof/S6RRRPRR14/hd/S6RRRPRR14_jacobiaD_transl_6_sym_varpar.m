% Zeitableitung der translatorischen Teilmatrix der analytischen Jacobi-Matrix für Segment Nr. 6 (0=Basis) von
% S6RRRPRR14
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
%   pkin=[a2,a3,a4,a5,a6,alpha2,d1,d2,d3,d5,d6]';
%
% Output:
% JaD_transl [3x6]
%   Zeitableitung der translatorischen Teilmatrix der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-02-26 22:23
% Revision: d75aae1ac561373cd3be920984c3df29a1c2ecc8 (2019-02-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JaD_transl = S6RRRPRR14_jacobiaD_transl_6_sym_varpar(qJ, qJD, r_i_i_C, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(3,1),zeros(11,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRPRR14_jacobiaD_transl_6_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRRPRR14_jacobiaD_transl_6_sym_varpar: qJD has to be [6x1] (double)');
assert(isa(r_i_i_C,'double') && isreal(r_i_i_C) && all(size(r_i_i_C) == [3 1]), ...
	'S6RRRPRR14_jacobiaD_transl_6_sym_varpar: Position vector r_i_i_C has to be [3x1] double');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6RRRPRR14_jacobiaD_transl_6_sym_varpar: pkin has to be [11x1] (double)');

%% Symbolic Calculation
% From jacobiaD_transl_6_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 22:23:42
% EndTime: 2019-02-26 22:23:42
% DurationCPUTime: 0.57s
% Computational Cost: add. (914->113), mult. (2230->180), div. (0->0), fcn. (2242->12), ass. (0->75)
t411 = sin(qJ(3));
t415 = cos(qJ(3));
t407 = qJD(5) + qJD(6);
t414 = cos(qJ(5));
t408 = qJ(5) + qJ(6);
t405 = sin(t408);
t406 = cos(t408);
t433 = t406 * r_i_i_C(1) - t405 * r_i_i_C(2);
t462 = qJD(5) * pkin(5);
t420 = t433 * t407 + t414 * t462 + qJD(4);
t410 = sin(qJ(5));
t432 = -t405 * r_i_i_C(1) - t406 * r_i_i_C(2);
t429 = qJ(4) - t432;
t422 = pkin(5) * t410 + t429;
t450 = pkin(3) + r_i_i_C(3) + pkin(11) + pkin(10);
t418 = (t450 * t411 - t422 * t415) * qJD(3) - t420 * t411;
t413 = sin(qJ(1));
t416 = cos(qJ(2));
t461 = cos(pkin(6));
t464 = cos(qJ(1));
t435 = t461 * t464;
t412 = sin(qJ(2));
t443 = t413 * t461;
t436 = t412 * t443;
t444 = t464 * qJD(1);
t452 = qJD(2) * t412;
t380 = -qJD(1) * t436 - t413 * t452 + (qJD(2) * t435 + t444) * t416;
t392 = t412 * t435 + t413 * t416;
t409 = sin(pkin(6));
t449 = t409 * t464;
t400 = t415 * t449;
t438 = t392 * t411 + t400;
t459 = t409 * t413;
t447 = qJD(1) * t459;
t467 = qJD(3) * t438 - t380 * t415 - t411 * t447;
t465 = t422 * t411 + t450 * t415 + pkin(2);
t463 = -t414 * pkin(5) - pkin(4) - pkin(9);
t448 = t464 * t416;
t394 = t448 - t436;
t460 = t394 * t411;
t458 = t409 * t415;
t457 = t409 * t416;
t456 = t413 * t412;
t430 = t416 * t435;
t377 = -qJD(1) * t430 - qJD(2) * t448 + (qJD(2) * t461 + qJD(1)) * t456;
t386 = -t413 * t458 + t460;
t440 = t386 * t407 - t377;
t393 = t464 * t412 + t416 * t443;
t378 = t392 * qJD(1) + t393 * qJD(2);
t428 = t394 * t415 + t411 * t459;
t371 = -qJD(1) * t400 + t428 * qJD(3) - t378 * t411;
t442 = -t393 * t407 + t371;
t367 = -t440 * t405 + t442 * t406;
t368 = t442 * t405 + t440 * t406;
t455 = t367 * r_i_i_C(1) - t368 * r_i_i_C(2);
t379 = t393 * qJD(1) + t392 * qJD(2);
t439 = -t407 * t438 - t379;
t437 = t411 * t449;
t451 = qJD(3) * t415;
t373 = -qJD(3) * t437 + t380 * t411 + t392 * t451 - t415 * t447;
t391 = -t430 + t456;
t441 = t391 * t407 - t373;
t454 = (t439 * t405 - t441 * t406) * r_i_i_C(1) + (t441 * t405 + t439 * t406) * r_i_i_C(2);
t389 = t409 * t412 * t411 - t461 * t415;
t446 = t409 * t452;
t427 = -t389 * t407 - t446;
t423 = t461 * t411 + t412 * t458;
t445 = qJD(2) * t457;
t381 = t423 * qJD(3) + t411 * t445;
t431 = t407 * t457 + t381;
t453 = (t427 * t405 + t431 * t406) * r_i_i_C(1) + (-t431 * t405 + t427 * t406) * r_i_i_C(2);
t426 = t433 - t463;
t421 = t432 * t407 - t410 * t462;
t372 = -t378 * t415 - qJD(3) * t460 + (t411 * t444 + t413 * t451) * t409;
t1 = [-t380 * pkin(2) - t438 * qJD(4) - t429 * t373 + ((t391 * t405 - t406 * t438) * r_i_i_C(1) + (t391 * t406 + t405 * t438) * r_i_i_C(2)) * t407 - t426 * t379 + (-t464 * pkin(1) - pkin(8) * t459) * qJD(1) + t450 * t467 + (-t373 * t410 + (t391 * t410 - t414 * t438) * qJD(5)) * pkin(5), t465 * t377 - t426 * t378 + t418 * t393 + t421 * t394, -t450 * t371 + t422 * t372 + t420 * t428, t371 (t371 * t414 + t377 * t410 + (-t386 * t410 - t393 * t414) * qJD(5)) * pkin(5) + t455, t455; -t378 * pkin(2) + t368 * r_i_i_C(1) + t367 * r_i_i_C(2) + t371 * qJ(4) + t386 * qJD(4) + t463 * t377 + (-pkin(1) * t413 + pkin(8) * t449) * qJD(1) + t450 * t372 + (t371 * t410 + (t386 * t414 - t393 * t410) * qJD(5)) * pkin(5), -t379 * t465 + t426 * t380 + t418 * t391 + t421 * t392, -t450 * t373 + t420 * (t392 * t415 - t437) - t422 * t467, t373 (t373 * t414 - t379 * t410 + (-t391 * t414 - t410 * t438) * qJD(5)) * pkin(5) + t454, t454; 0 ((-qJD(2) * t465 + t421) * t412 + (t426 * qJD(2) - t418) * t416) * t409, -t450 * t381 + t420 * t423 + t422 * (-t389 * qJD(3) + t415 * t445) t381 (-t410 * t446 + t381 * t414 + (-t389 * t410 + t414 * t457) * qJD(5)) * pkin(5) + t453, t453;];
JaD_transl  = t1;
