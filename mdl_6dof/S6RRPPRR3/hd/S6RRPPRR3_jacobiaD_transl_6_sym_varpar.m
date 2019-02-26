% Zeitableitung der translatorischen Teilmatrix der analytischen Jacobi-Matrix für Segment Nr. 6 (0=Basis) von
% S6RRPPRR3
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
%   pkin=[a2,a3,a4,a5,a6,alpha2,d1,d2,d5,d6,theta3,theta4]';
%
% Output:
% JaD_transl [3x6]
%   Zeitableitung der translatorischen Teilmatrix der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-02-26 21:29
% Revision: d75aae1ac561373cd3be920984c3df29a1c2ecc8 (2019-02-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JaD_transl = S6RRPPRR3_jacobiaD_transl_6_sym_varpar(qJ, qJD, r_i_i_C, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(3,1),zeros(12,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPPRR3_jacobiaD_transl_6_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRPPRR3_jacobiaD_transl_6_sym_varpar: qJD has to be [6x1] (double)');
assert(isa(r_i_i_C,'double') && isreal(r_i_i_C) && all(size(r_i_i_C) == [3 1]), ...
	'S6RRPPRR3_jacobiaD_transl_6_sym_varpar: Position vector r_i_i_C has to be [3x1] double');
assert(isreal(pkin) && all(size(pkin) == [12 1]), ...
  'S6RRPPRR3_jacobiaD_transl_6_sym_varpar: pkin has to be [12x1] (double)');

%% Symbolic Calculation
% From jacobiaD_transl_6_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 21:29:36
% EndTime: 2019-02-26 21:29:37
% DurationCPUTime: 0.82s
% Computational Cost: add. (1075->131), mult. (2590->217), div. (0->0), fcn. (2849->14), ass. (0->80)
t469 = cos(pkin(6));
t467 = sin(pkin(11));
t520 = cos(pkin(11));
t522 = cos(qJ(2));
t493 = t522 * t520;
t472 = sin(qJ(2));
t508 = qJD(2) * t472;
t526 = -qJD(2) * t493 + t467 * t508;
t442 = t526 * t469;
t497 = t472 * t520;
t484 = t522 * t467 + t497;
t447 = t484 * t469;
t498 = t522 * qJD(2);
t449 = -qJD(2) * t497 - t467 * t498;
t473 = sin(qJ(1));
t475 = cos(qJ(1));
t483 = -t472 * t467 + t493;
t510 = qJD(1) * t473;
t417 = t447 * t510 - t473 * t449 + (-qJD(1) * t483 + t442) * t475;
t429 = t447 * t475 + t473 * t483;
t465 = pkin(12) + qJ(5);
t463 = sin(t465);
t464 = cos(t465);
t468 = sin(pkin(6));
t501 = t468 * t510;
t514 = t468 * t475;
t504 = t464 * t514;
t411 = (-qJD(5) * t429 + t501) * t463 - qJD(5) * t504 - t417 * t464;
t443 = qJD(2) * t447;
t509 = qJD(1) * t475;
t482 = t469 * t483;
t527 = qJD(1) * t482 + qJD(2) * t483;
t418 = -t475 * t443 - t527 * t473 - t484 * t509;
t471 = sin(qJ(6));
t474 = cos(qJ(6));
t533 = t411 * t471 + t418 * t474;
t532 = -t411 * t474 + t418 * t471;
t424 = -t429 * t464 + t463 * t514;
t428 = -t473 * t484 + t475 * t482;
t531 = -t424 * t471 + t428 * t474;
t530 = t424 * t474 + t428 * t471;
t487 = qJD(6) * (r_i_i_C(1) * t471 + r_i_i_C(2) * t474);
t490 = r_i_i_C(1) * t474 - r_i_i_C(2) * t471 + pkin(5);
t523 = pkin(10) + r_i_i_C(3);
t529 = (t490 * t463 - t523 * t464) * qJD(5) + t464 * t487;
t461 = cos(pkin(12)) * pkin(4) + pkin(3);
t524 = t523 * t463 + t490 * t464 + t461;
t521 = pkin(2) * t469;
t515 = t468 * t473;
t512 = t472 * t473;
t511 = t472 * t475;
t507 = qJD(6) * t471;
t506 = qJD(6) * t474;
t505 = pkin(2) * t508;
t503 = t522 * t473;
t502 = t522 * t475;
t500 = t468 * t509;
t494 = -t472 * t521 + (pkin(4) * sin(pkin(12)) + pkin(8) + qJ(3)) * t468;
t446 = t484 * t468;
t434 = t446 * t464 + t463 * t469;
t491 = -t446 * t463 + t464 * t469;
t430 = t473 * t447 - t475 * t483;
t488 = t430 * t463 + t464 * t515;
t426 = -t430 * t464 + t463 * t515;
t478 = -t429 * qJD(1) + t473 * t442 + t449 * t475;
t477 = t424 * qJD(5) + t417 * t463 + t464 * t501;
t470 = -pkin(9) - qJ(4);
t462 = t522 * pkin(2) + pkin(1);
t450 = -t468 * qJD(3) + t498 * t521;
t445 = t483 * t468;
t441 = qJD(2) * t446;
t440 = t526 * t468;
t431 = -t473 * t482 - t475 * t484;
t421 = t491 * qJD(5) - t440 * t464;
t415 = -t473 * t443 + t527 * t475 - t484 * t510;
t409 = t488 * qJD(5) + t463 * t500 + t464 * t478;
t408 = t426 * qJD(5) + t463 * t478 - t464 * t500;
t407 = t409 * t474 + t415 * t471 + (-t426 * t471 - t431 * t474) * qJD(6);
t406 = -t409 * t471 + t415 * t474 + (-t426 * t474 + t431 * t471) * qJD(6);
t1 = [t532 * r_i_i_C(1) + t533 * r_i_i_C(2) - t411 * pkin(5) + t417 * t461 - t418 * t470 + t428 * qJD(4) + t473 * t505 - t475 * t450 + t523 * t477 + (t531 * r_i_i_C(1) - t530 * r_i_i_C(2)) * qJD(6) + (-t475 * t462 - t494 * t473) * qJD(1) (-t430 * t506 + t471 * t478) * r_i_i_C(1) + (t430 * t507 + t474 * t478) * r_i_i_C(2) - t478 * t470 - t430 * qJD(4) - t524 * t415 + ((t469 * t512 - t502) * qJD(2) + (-t469 * t502 + t512) * qJD(1)) * pkin(2) - t529 * t431, t500, t415, -t490 * t408 + t523 * t409 - t488 * t487, r_i_i_C(1) * t406 - t407 * r_i_i_C(2); -t475 * t505 + t409 * pkin(5) + t407 * r_i_i_C(1) + t406 * r_i_i_C(2) - t431 * qJD(4) - t415 * t470 + t478 * t461 - t473 * t450 + t523 * t408 + (-t462 * t473 + t494 * t475) * qJD(1) (-t417 * t471 + t429 * t506) * r_i_i_C(1) + (-t417 * t474 - t429 * t507) * r_i_i_C(2) + t417 * t470 + t429 * qJD(4) + t524 * t418 + ((-t469 * t511 - t503) * qJD(2) + (-t469 * t503 - t511) * qJD(1)) * pkin(2) - t529 * t428, t501, -t418, t523 * t411 - (-t429 * t463 - t504) * t487 + t490 * t477, -t533 * r_i_i_C(1) + t532 * r_i_i_C(2) + (t530 * r_i_i_C(1) + t531 * r_i_i_C(2)) * qJD(6); 0 (-t440 * t471 + t446 * t506) * r_i_i_C(1) + (-t440 * t474 - t446 * t507) * r_i_i_C(2) + t440 * t470 + t446 * qJD(4) - t468 * t505 - t524 * t441 - t529 * t445, 0, t441, t523 * t421 - t491 * t487 + t490 * (-t434 * qJD(5) + t440 * t463) (-t421 * t471 + t441 * t474) * r_i_i_C(1) + (-t421 * t474 - t441 * t471) * r_i_i_C(2) + ((-t434 * t474 + t445 * t471) * r_i_i_C(1) + (t434 * t471 + t445 * t474) * r_i_i_C(2)) * qJD(6);];
JaD_transl  = t1;
