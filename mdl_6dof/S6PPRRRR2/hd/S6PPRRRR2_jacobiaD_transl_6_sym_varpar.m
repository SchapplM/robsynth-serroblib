% Zeitableitung der translatorischen Teilmatrix der analytischen Jacobi-Matrix für Segment Nr. 6 (0=Basis) von
% S6PPRRRR2
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
%   pkin=[a2,a3,a4,a5,a6,alpha2,alpha3,d3,d4,d5,d6,theta1,theta2]';
%
% Output:
% JaD_transl [3x6]
%   Zeitableitung der translatorischen Teilmatrix der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-02-26 19:43
% Revision: d75aae1ac561373cd3be920984c3df29a1c2ecc8 (2019-02-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JaD_transl = S6PPRRRR2_jacobiaD_transl_6_sym_varpar(qJ, qJD, r_i_i_C, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(3,1),zeros(13,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PPRRRR2_jacobiaD_transl_6_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6PPRRRR2_jacobiaD_transl_6_sym_varpar: qJD has to be [6x1] (double)');
assert(isa(r_i_i_C,'double') && isreal(r_i_i_C) && all(size(r_i_i_C) == [3 1]), ...
	'S6PPRRRR2_jacobiaD_transl_6_sym_varpar: Position vector r_i_i_C has to be [3x1] double');
assert(isreal(pkin) && all(size(pkin) == [13 1]), ...
  'S6PPRRRR2_jacobiaD_transl_6_sym_varpar: pkin has to be [13x1] (double)');

%% Symbolic Calculation
% From jacobiaD_transl_6_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 19:43:18
% EndTime: 2019-02-26 19:43:19
% DurationCPUTime: 0.59s
% Computational Cost: add. (1049->100), mult. (2843->174), div. (0->0), fcn. (3310->16), ass. (0->80)
t517 = sin(pkin(13));
t518 = sin(pkin(12));
t492 = t518 * t517;
t521 = cos(pkin(13));
t522 = cos(pkin(12));
t498 = t522 * t521;
t524 = cos(pkin(6));
t476 = -t524 * t498 + t492;
t523 = cos(pkin(7));
t474 = t476 * t523;
t519 = sin(pkin(7));
t520 = sin(pkin(6));
t496 = t520 * t519;
t485 = t522 * t496;
t531 = t474 + t485;
t494 = t518 * t521;
t497 = t522 * t517;
t477 = t524 * t494 + t497;
t493 = t518 * t520;
t530 = t477 * t523 - t519 * t493;
t499 = t523 * t520;
t529 = t521 * t499 + t524 * t519;
t451 = t524 * t497 + t494;
t467 = sin(qJ(3));
t526 = cos(qJ(3));
t434 = t451 * t467 + t531 * t526;
t528 = t530 * t526;
t495 = t520 * t517;
t442 = t467 * t495 - t529 * t526;
t464 = qJ(5) + qJ(6);
t461 = sin(t464);
t462 = cos(t464);
t463 = qJD(5) + qJD(6);
t465 = sin(qJ(5));
t527 = qJD(5) * t465 * pkin(5) + (t461 * r_i_i_C(1) + t462 * r_i_i_C(2)) * t463;
t525 = r_i_i_C(3) + pkin(11) + pkin(10);
t516 = t461 * t463;
t515 = t462 * t463;
t508 = t451 * t526;
t435 = -t531 * t467 + t508;
t444 = t476 * t519 - t522 * t499;
t466 = sin(qJ(4));
t469 = cos(qJ(4));
t427 = t435 * t469 + t444 * t466;
t511 = qJD(3) * t467;
t431 = -t485 * t511 + (-t467 * t474 + t508) * qJD(3);
t504 = t427 * t463 - t431;
t430 = t434 * qJD(3);
t491 = -t435 * t466 + t444 * t469;
t421 = t491 * qJD(4) - t430 * t469;
t507 = -t434 * t463 - t421;
t514 = (t507 * t461 - t504 * t462) * r_i_i_C(1) + (t504 * t461 + t507 * t462) * r_i_i_C(2);
t452 = -t524 * t492 + t498;
t437 = t452 * t526 - t530 * t467;
t445 = t477 * t519 + t523 * t493;
t429 = t437 * t469 + t445 * t466;
t433 = t437 * qJD(3);
t503 = t429 * t463 - t433;
t432 = t528 * qJD(3) + t452 * t511;
t490 = -t437 * t466 + t445 * t469;
t423 = t490 * qJD(4) - t432 * t469;
t436 = t452 * t467 + t528;
t506 = -t436 * t463 - t423;
t513 = (t506 * t461 - t503 * t462) * r_i_i_C(1) + (t503 * t461 + t506 * t462) * r_i_i_C(2);
t443 = t529 * t467 + t526 * t495;
t450 = -t521 * t496 + t524 * t523;
t439 = t443 * t469 + t450 * t466;
t441 = t443 * qJD(3);
t502 = t439 * t463 - t441;
t440 = t442 * qJD(3);
t489 = -t443 * t466 + t450 * t469;
t425 = t489 * qJD(4) - t440 * t469;
t505 = -t442 * t463 - t425;
t512 = (t505 * t461 - t502 * t462) * r_i_i_C(1) + (t502 * t461 + t505 * t462) * r_i_i_C(2);
t468 = cos(qJ(5));
t510 = qJD(5) * t468;
t488 = t468 * pkin(5) + r_i_i_C(1) * t462 - r_i_i_C(2) * t461 + pkin(4);
t479 = -t525 * t466 - t488 * t469 - pkin(3);
t471 = t527 * t469 + (t488 * t466 - t525 * t469) * qJD(4);
t1 = [0, 0 (-t432 * t461 + t437 * t515) * r_i_i_C(1) + (-t432 * t462 - t437 * t516) * r_i_i_C(2) - t432 * pkin(9) + (-t432 * t465 + t437 * t510) * pkin(5) + t479 * t433 + t471 * t436, t525 * t423 - t527 * t490 + t488 * (-t429 * qJD(4) + t432 * t466) (-t423 * t465 + t433 * t468 + (-t429 * t468 - t436 * t465) * qJD(5)) * pkin(5) + t513, t513; 0, 0 (-t430 * t461 + t435 * t515) * r_i_i_C(1) + (-t430 * t462 - t435 * t516) * r_i_i_C(2) - t430 * pkin(9) + (-t430 * t465 + t435 * t510) * pkin(5) + t479 * t431 + t471 * t434, t525 * t421 - t527 * t491 + t488 * (-t427 * qJD(4) + t430 * t466) (-t421 * t465 + t431 * t468 + (-t427 * t468 - t434 * t465) * qJD(5)) * pkin(5) + t514, t514; 0, 0 (-t440 * t461 + t443 * t515) * r_i_i_C(1) + (-t440 * t462 - t443 * t516) * r_i_i_C(2) - t440 * pkin(9) + (-t440 * t465 + t443 * t510) * pkin(5) + t479 * t441 + t471 * t442, t525 * t425 - t527 * t489 + t488 * (-t439 * qJD(4) + t440 * t466) (-t425 * t465 + t441 * t468 + (-t439 * t468 - t442 * t465) * qJD(5)) * pkin(5) + t512, t512;];
JaD_transl  = t1;
