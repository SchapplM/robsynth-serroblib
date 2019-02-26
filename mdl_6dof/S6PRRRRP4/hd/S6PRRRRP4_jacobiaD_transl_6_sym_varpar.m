% Zeitableitung der translatorischen Teilmatrix der analytischen Jacobi-Matrix für Segment Nr. 6 (0=Basis) von
% S6PRRRRP4
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
%   pkin=[a2,a3,a4,a5,a6,alpha2,d2,d3,d4,d5,theta1]';
%
% Output:
% JaD_transl [3x6]
%   Zeitableitung der translatorischen Teilmatrix der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-02-26 20:17
% Revision: d75aae1ac561373cd3be920984c3df29a1c2ecc8 (2019-02-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JaD_transl = S6PRRRRP4_jacobiaD_transl_6_sym_varpar(qJ, qJD, r_i_i_C, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(3,1),zeros(11,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRRRRP4_jacobiaD_transl_6_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6PRRRRP4_jacobiaD_transl_6_sym_varpar: qJD has to be [6x1] (double)');
assert(isa(r_i_i_C,'double') && isreal(r_i_i_C) && all(size(r_i_i_C) == [3 1]), ...
	'S6PRRRRP4_jacobiaD_transl_6_sym_varpar: Position vector r_i_i_C has to be [3x1] double');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6PRRRRP4_jacobiaD_transl_6_sym_varpar: pkin has to be [11x1] (double)');

%% Symbolic Calculation
% From jacobiaD_transl_6_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 20:16:56
% EndTime: 2019-02-26 20:16:57
% DurationCPUTime: 0.58s
% Computational Cost: add. (946->110), mult. (1948->180), div. (0->0), fcn. (2026->12), ass. (0->73)
t477 = pkin(5) + r_i_i_C(1);
t474 = r_i_i_C(3) + qJ(6);
t427 = cos(qJ(4));
t417 = pkin(4) * t427 + pkin(3);
t428 = cos(qJ(3));
t421 = qJ(4) + qJ(5);
t418 = sin(t421);
t424 = sin(qJ(4));
t476 = pkin(4) * t424;
t444 = -qJD(4) * t476 + qJD(6) * t418;
t425 = sin(qJ(3));
t464 = qJD(3) * t425;
t475 = r_i_i_C(2) + pkin(10) + pkin(9);
t481 = -t417 * t464 + (t475 * qJD(3) + t444) * t428;
t426 = sin(qJ(2));
t423 = sin(pkin(6));
t468 = t423 * t425;
t473 = cos(pkin(6));
t440 = -t426 * t468 + t473 * t428;
t429 = cos(qJ(2));
t466 = t423 * t429;
t459 = qJD(2) * t466;
t397 = t440 * qJD(3) + t428 * t459;
t420 = qJD(4) + qJD(5);
t480 = -t420 * t466 + t397;
t422 = sin(pkin(11));
t456 = t422 * t473;
t472 = cos(pkin(11));
t408 = -t426 * t456 + t472 * t429;
t479 = (qJD(2) * t428 - t420) * t426 + t429 * t464;
t471 = t418 * t420;
t419 = cos(t421);
t470 = t419 * t420;
t469 = t420 * t428;
t467 = t423 * t428;
t465 = qJD(2) * t426;
t463 = qJD(4) * t427;
t462 = qJD(6) * t419;
t460 = t423 * t465;
t455 = t423 * t472;
t447 = t473 * t472;
t445 = t429 * t447;
t401 = -qJD(2) * t445 + t422 * t465;
t406 = t422 * t429 + t426 * t447;
t441 = -t406 * t425 - t428 * t455;
t386 = t441 * qJD(3) - t401 * t428;
t405 = t422 * t426 - t445;
t453 = t405 * t420 + t386;
t407 = t472 * t426 + t429 * t456;
t403 = t407 * qJD(2);
t443 = -t408 * t425 + t422 * t467;
t388 = t443 * qJD(3) - t403 * t428;
t452 = t407 * t420 + t388;
t449 = t405 * t469 - t401;
t448 = t407 * t469 - t403;
t446 = (qJD(2) - t469) * t429;
t395 = t408 * t428 + t422 * t468;
t442 = -t428 * t417 - t475 * t425 - pkin(2);
t410 = t473 * t425 + t426 * t467;
t439 = -t406 * t428 + t425 * t455;
t402 = t406 * qJD(2);
t368 = -t402 * t419 + t453 * t418 - t439 * t470;
t438 = -(-t405 * t418 + t419 * t439) * qJD(6) + t474 * (t402 * t418 + t453 * t419 + t439 * t471) - t477 * t368;
t404 = t408 * qJD(2);
t370 = t395 * t470 - t404 * t419 + t452 * t418;
t437 = -(-t395 * t419 - t407 * t418) * qJD(6) + t474 * (-t395 * t471 + t404 * t418 + t452 * t419) - t477 * t370;
t379 = t410 * t470 + t480 * t418 - t419 * t460;
t436 = -(-t410 * t419 + t418 * t466) * qJD(6) + t474 * (-t410 * t471 + t418 * t460 + t480 * t419) - t477 * t379;
t435 = t474 * t418 + t477 * t419 + t417;
t434 = -t402 * t428 + t405 * t464 + t406 * t420;
t433 = -t404 * t428 + t407 * t464 + t408 * t420;
t432 = (-t477 * t418 + t474 * t419) * t420 + t444;
t1 = [0, -t408 * t462 - t403 * pkin(8) + t477 * (t448 * t418 + t433 * t419) + t474 * (t433 * t418 - t448 * t419) + (-t403 * t424 + t408 * t463) * pkin(4) - t481 * t407 + t442 * t404, t475 * t388 + t435 * (-t395 * qJD(3) + t403 * t425) + t432 * t443 (-t388 * t424 + t404 * t427 + (-t395 * t427 - t407 * t424) * qJD(4)) * pkin(4) + t437, t437, t370; 0, -t406 * t462 - t401 * pkin(8) + t477 * (t449 * t418 + t434 * t419) + t474 * (t434 * t418 - t449 * t419) + (-t401 * t424 + t406 * t463) * pkin(4) - t481 * t405 + t442 * t402, t475 * t386 + t435 * (t439 * qJD(3) + t401 * t425) + t432 * t441 (-t386 * t424 + t402 * t427 + (-t405 * t424 + t427 * t439) * qJD(4)) * pkin(4) + t438, t438, t368; 0 (t477 * (t418 * t446 - t479 * t419) - t474 * (t479 * t418 + t419 * t446) + (pkin(4) * t463 - t462) * t426 + t481 * t429 + ((pkin(8) + t476) * t429 + t442 * t426) * qJD(2)) * t423, t475 * t397 + t435 * (-t410 * qJD(3) - t425 * t459) + t432 * t440 (t427 * t460 - t397 * t424 + (-t410 * t427 + t424 * t466) * qJD(4)) * pkin(4) + t436, t436, t379;];
JaD_transl  = t1;
