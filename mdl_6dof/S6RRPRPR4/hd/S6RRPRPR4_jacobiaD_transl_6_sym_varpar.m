% Zeitableitung der translatorischen Teilmatrix der analytischen Jacobi-Matrix für Segment Nr. 6 (0=Basis) von
% S6RRPRPR4
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
% Datum: 2019-02-26 21:39
% Revision: d75aae1ac561373cd3be920984c3df29a1c2ecc8 (2019-02-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JaD_transl = S6RRPRPR4_jacobiaD_transl_6_sym_varpar(qJ, qJD, r_i_i_C, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(3,1),zeros(12,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPRPR4_jacobiaD_transl_6_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRPRPR4_jacobiaD_transl_6_sym_varpar: qJD has to be [6x1] (double)');
assert(isa(r_i_i_C,'double') && isreal(r_i_i_C) && all(size(r_i_i_C) == [3 1]), ...
	'S6RRPRPR4_jacobiaD_transl_6_sym_varpar: Position vector r_i_i_C has to be [3x1] double');
assert(isreal(pkin) && all(size(pkin) == [12 1]), ...
  'S6RRPRPR4_jacobiaD_transl_6_sym_varpar: pkin has to be [12x1] (double)');

%% Symbolic Calculation
% From jacobiaD_transl_6_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 21:39:35
% EndTime: 2019-02-26 21:39:36
% DurationCPUTime: 0.89s
% Computational Cost: add. (1128->148), mult. (2753->244), div. (0->0), fcn. (3018->14), ass. (0->82)
t469 = cos(pkin(6));
t467 = sin(pkin(11));
t521 = cos(pkin(11));
t524 = cos(qJ(2));
t495 = t524 * t521;
t473 = sin(qJ(2));
t509 = qJD(2) * t473;
t528 = -qJD(2) * t495 + t467 * t509;
t443 = t528 * t469;
t498 = t473 * t521;
t486 = t467 * t524 + t498;
t448 = t486 * t469;
t499 = t524 * qJD(2);
t450 = -qJD(2) * t498 - t467 * t499;
t474 = sin(qJ(1));
t477 = cos(qJ(1));
t485 = -t473 * t467 + t495;
t511 = qJD(1) * t474;
t418 = t448 * t511 - t474 * t450 + t477 * (-qJD(1) * t485 + t443);
t430 = t477 * t448 + t474 * t485;
t466 = qJ(4) + pkin(12);
t464 = sin(t466);
t465 = cos(t466);
t468 = sin(pkin(6));
t502 = t468 * t511;
t515 = t468 * t477;
t505 = t465 * t515;
t412 = t464 * (-qJD(4) * t430 + t502) - qJD(4) * t505 - t418 * t465;
t444 = qJD(2) * t448;
t510 = qJD(1) * t477;
t484 = t485 * t469;
t529 = qJD(1) * t484 + t485 * qJD(2);
t419 = -t477 * t444 - t529 * t474 - t486 * t510;
t471 = sin(qJ(6));
t475 = cos(qJ(6));
t535 = t412 * t471 + t419 * t475;
t534 = -t412 * t475 + t419 * t471;
t425 = -t430 * t465 + t464 * t515;
t429 = -t474 * t486 + t477 * t484;
t533 = -t425 * t471 + t429 * t475;
t532 = t425 * t475 + t429 * t471;
t472 = sin(qJ(4));
t489 = qJD(6) * (t471 * r_i_i_C(1) + t475 * r_i_i_C(2));
t492 = t475 * r_i_i_C(1) - t471 * r_i_i_C(2) + pkin(5);
t525 = r_i_i_C(3) + pkin(10);
t531 = (t472 * pkin(4) + t492 * t464 - t465 * t525) * qJD(4) + t465 * t489;
t476 = cos(qJ(4));
t462 = t476 * pkin(4) + pkin(3);
t527 = t464 * t525 + t492 * t465 + t462;
t523 = pkin(2) * t469;
t516 = t468 * t474;
t513 = t473 * t474;
t512 = t473 * t477;
t508 = qJD(6) * t471;
t507 = qJD(6) * t475;
t506 = pkin(2) * t509;
t504 = t524 * t474;
t503 = t524 * t477;
t501 = t468 * t510;
t447 = t486 * t468;
t435 = t447 * t465 + t469 * t464;
t493 = -t447 * t464 + t469 * t465;
t431 = t474 * t448 - t477 * t485;
t490 = t431 * t464 + t465 * t516;
t427 = -t431 * t465 + t464 * t516;
t480 = -qJD(1) * t430 + t474 * t443 + t477 * t450;
t479 = qJD(4) * t425 + t418 * t464 + t465 * t502;
t470 = -qJ(5) - pkin(9);
t463 = pkin(2) * t524 + pkin(1);
t451 = -t468 * qJD(3) + t499 * t523;
t449 = t473 * t523 + (-pkin(8) - qJ(3)) * t468;
t446 = t485 * t468;
t442 = qJD(2) * t447;
t441 = t528 * t468;
t432 = -t474 * t484 - t477 * t486;
t422 = qJD(4) * t493 - t441 * t465;
t416 = -t474 * t444 + t529 * t477 - t486 * t511;
t410 = qJD(4) * t490 + t464 * t501 + t465 * t480;
t409 = qJD(4) * t427 + t464 * t480 - t465 * t501;
t408 = t410 * t475 + t416 * t471 + (-t427 * t471 - t432 * t475) * qJD(6);
t407 = -t410 * t471 + t416 * t475 + (-t427 * t475 + t432 * t471) * qJD(6);
t1 = [t534 * r_i_i_C(1) + t535 * r_i_i_C(2) - t412 * pkin(5) + t418 * t462 - t419 * t470 + t429 * qJD(5) + t474 * t506 - t477 * t451 + t525 * t479 + (t533 * r_i_i_C(1) - t532 * r_i_i_C(2)) * qJD(6) + (t474 * t449 - t477 * t463) * qJD(1) + (-t472 * t502 + (t430 * t472 + t476 * t515) * qJD(4)) * pkin(4) (-t431 * t507 + t471 * t480) * r_i_i_C(1) + (t431 * t508 + t475 * t480) * r_i_i_C(2) - t480 * t470 - t431 * qJD(5) - t527 * t416 + ((t469 * t513 - t503) * qJD(2) + (-t469 * t503 + t513) * qJD(1)) * pkin(2) - t531 * t432, t501, t525 * t410 - t490 * t489 - t492 * t409 + (t476 * t501 - t480 * t472 + (t431 * t476 - t472 * t516) * qJD(4)) * pkin(4), t416, t407 * r_i_i_C(1) - t408 * r_i_i_C(2); -t477 * t506 + t410 * pkin(5) + t408 * r_i_i_C(1) + t407 * r_i_i_C(2) - t432 * qJD(5) - t416 * t470 + t480 * t462 - t474 * t451 + t525 * t409 + (-t449 * t477 - t463 * t474) * qJD(1) + (t472 * t501 + (t431 * t472 + t476 * t516) * qJD(4)) * pkin(4) (-t418 * t471 + t430 * t507) * r_i_i_C(1) + (-t418 * t475 - t430 * t508) * r_i_i_C(2) + t418 * t470 + t430 * qJD(5) + t527 * t419 + ((-t469 * t512 - t504) * qJD(2) + (-t469 * t504 - t512) * qJD(1)) * pkin(2) - t531 * t429, t502, t525 * t412 - (-t430 * t464 - t505) * t489 + t492 * t479 + (t476 * t502 + t418 * t472 + (-t430 * t476 + t472 * t515) * qJD(4)) * pkin(4), -t419, -t535 * r_i_i_C(1) + t534 * r_i_i_C(2) + (t532 * r_i_i_C(1) + t533 * r_i_i_C(2)) * qJD(6); 0 (-t441 * t471 + t447 * t507) * r_i_i_C(1) + (-t441 * t475 - t447 * t508) * r_i_i_C(2) + t441 * t470 + t447 * qJD(5) - t468 * t506 - t527 * t442 - t531 * t446, 0, t525 * t422 - t493 * t489 + t492 * (-qJD(4) * t435 + t441 * t464) + (t441 * t472 + (-t447 * t476 - t469 * t472) * qJD(4)) * pkin(4), t442 (-t422 * t471 + t442 * t475) * r_i_i_C(1) + (-t422 * t475 - t442 * t471) * r_i_i_C(2) + ((-t435 * t475 + t446 * t471) * r_i_i_C(1) + (t435 * t471 + t446 * t475) * r_i_i_C(2)) * qJD(6);];
JaD_transl  = t1;
