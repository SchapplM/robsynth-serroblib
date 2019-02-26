% Zeitableitung der translatorischen Teilmatrix der analytischen Jacobi-Matrix für Segment Nr. 6 (0=Basis) von
% S6RRRRRR6
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
%   pkin=[a2,a3,a4,a5,a6,alpha2,d1,d2,d3,d4,d5,d6]';
%
% Output:
% JaD_transl [3x6]
%   Zeitableitung der translatorischen Teilmatrix der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-02-26 22:50
% Revision: d75aae1ac561373cd3be920984c3df29a1c2ecc8 (2019-02-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JaD_transl = S6RRRRRR6_jacobiaD_transl_6_sym_varpar(qJ, qJD, r_i_i_C, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(3,1),zeros(12,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRRRR6_jacobiaD_transl_6_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRRRRR6_jacobiaD_transl_6_sym_varpar: qJD has to be [6x1] (double)');
assert(isa(r_i_i_C,'double') && isreal(r_i_i_C) && all(size(r_i_i_C) == [3 1]), ...
	'S6RRRRRR6_jacobiaD_transl_6_sym_varpar: Position vector r_i_i_C has to be [3x1] double');
assert(isreal(pkin) && all(size(pkin) == [12 1]), ...
  'S6RRRRRR6_jacobiaD_transl_6_sym_varpar: pkin has to be [12x1] (double)');

%% Symbolic Calculation
% From jacobiaD_transl_6_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 22:50:06
% EndTime: 2019-02-26 22:50:07
% DurationCPUTime: 0.95s
% Computational Cost: add. (1416->140), mult. (2311->224), div. (0->0), fcn. (2307->14), ass. (0->89)
t476 = sin(qJ(1));
t537 = cos(pkin(6));
t493 = qJD(2) * t537 + qJD(1);
t475 = sin(qJ(2));
t508 = t476 * t537;
t499 = t475 * t508;
t520 = qJD(2) * t475;
t479 = cos(qJ(2));
t480 = cos(qJ(1));
t525 = t480 * t479;
t443 = -qJD(1) * t499 - t476 * t520 + t493 * t525;
t471 = qJ(3) + qJ(4);
t465 = sin(t471);
t467 = cos(t471);
t469 = qJD(3) + qJD(4);
t507 = t480 * t537;
t452 = t475 * t507 + t476 * t479;
t472 = sin(pkin(6));
t521 = qJD(1) * t472;
t512 = t476 * t521;
t491 = -t452 * t469 + t512;
t527 = t472 * t480;
t513 = t467 * t527;
t431 = t443 * t467 + t491 * t465 - t469 * t513;
t498 = t479 * t507;
t526 = t476 * t475;
t451 = -t498 + t526;
t468 = qJD(5) + qJD(6);
t504 = t451 * t468 + t431;
t474 = sin(qJ(3));
t470 = qJ(5) + qJ(6);
t464 = sin(t470);
t466 = cos(t470);
t496 = r_i_i_C(1) * t464 + r_i_i_C(2) * t466;
t473 = sin(qJ(5));
t539 = t473 * pkin(5);
t488 = -qJD(5) * t539 - t496 * t468;
t538 = r_i_i_C(3) + pkin(12) + pkin(11);
t477 = cos(qJ(5));
t462 = t477 * pkin(5) + pkin(4);
t497 = t466 * r_i_i_C(1) - t464 * r_i_i_C(2);
t549 = t462 + t497;
t483 = qJD(3) * t474 * pkin(3) + t465 * t469 * t549 - (t538 * t469 + t488) * t467;
t453 = t480 * t475 + t479 * t508;
t442 = t453 * qJD(1) + t452 * qJD(2);
t446 = -t452 * t467 + t465 * t527;
t546 = -t446 * t468 - t442;
t454 = -t499 + t525;
t511 = t480 * t521;
t544 = -t454 * t469 + t511;
t478 = cos(qJ(3));
t463 = t478 * pkin(3) + pkin(2);
t542 = t538 * t465 + t467 * t549 + t463;
t534 = t464 * t468;
t532 = t466 * t468;
t531 = t472 * t475;
t530 = t472 * t476;
t529 = t472 * t478;
t528 = t472 * t479;
t519 = qJD(2) * t479;
t440 = -qJD(1) * t498 - t480 * t519 + t493 * t526;
t448 = t454 * t467 + t465 * t530;
t503 = -t448 * t468 - t440;
t441 = t452 * qJD(1) + t453 * qJD(2);
t494 = t469 * t530 - t441;
t429 = t465 * t544 + t494 * t467;
t506 = t453 * t468 + t429;
t416 = -t506 * t464 + t503 * t466;
t417 = t503 * t464 + t506 * t466;
t524 = t416 * r_i_i_C(1) - t417 * r_i_i_C(2);
t523 = (-t504 * t464 - t466 * t546) * r_i_i_C(1) + (t464 * t546 - t504 * t466) * r_i_i_C(2);
t514 = t467 * t531;
t450 = t537 * t465 + t514;
t510 = t472 * t520;
t490 = -t450 * t468 + t510;
t509 = t472 * t519;
t489 = t537 * t469 + t509;
t515 = t465 * t531;
t439 = t489 * t467 - t469 * t515;
t495 = t468 * t528 - t439;
t522 = (t495 * t464 + t490 * t466) * r_i_i_C(1) + (-t490 * t464 + t495 * t466) * r_i_i_C(2);
t518 = qJD(5) * t477;
t430 = t491 * t467 + (t469 * t527 - t443) * t465;
t428 = t494 * t465 - t467 * t544;
t486 = t488 * (-t454 * t465 + t467 * t530) + t538 * t429 - t549 * t428;
t485 = t488 * (-t452 * t465 - t513) + t538 * t431 + t549 * t430;
t484 = t488 * (t537 * t467 - t515) + t538 * t439 + t549 * (-t489 * t465 - t469 * t514);
t482 = -pkin(10) - pkin(9);
t1 = [-t431 * t462 + t442 * t482 - t443 * t463 + (-t504 * r_i_i_C(1) + r_i_i_C(2) * t546) * t466 + (r_i_i_C(1) * t546 + t504 * r_i_i_C(2)) * t464 + t538 * t430 + (-t480 * pkin(1) - pkin(8) * t530) * qJD(1) + (-t442 * t473 + (-t446 * t473 - t451 * t477) * qJD(5)) * pkin(5) + (-t474 * t512 + (t452 * t474 + t478 * t527) * qJD(3)) * pkin(3) (-t441 * t464 + t454 * t532) * r_i_i_C(1) + (-t441 * t466 - t454 * t534) * r_i_i_C(2) + t441 * t482 + (-t441 * t473 + t454 * t518) * pkin(5) + t542 * t440 + t483 * t453 (t478 * t511 + t441 * t474 + (-t454 * t478 - t474 * t530) * qJD(3)) * pkin(3) + t486, t486 (-t429 * t473 - t440 * t477 + (-t448 * t477 - t453 * t473) * qJD(5)) * pkin(5) + t524, t524; t417 * r_i_i_C(1) + t416 * r_i_i_C(2) + t429 * t462 + t440 * t482 - t441 * t463 + t538 * t428 + (-pkin(1) * t476 + pkin(8) * t527) * qJD(1) + (-t440 * t473 + (-t448 * t473 + t453 * t477) * qJD(5)) * pkin(5) + (t474 * t511 + (-t454 * t474 + t476 * t529) * qJD(3)) * pkin(3) (t443 * t464 + t452 * t532) * r_i_i_C(1) + (t443 * t466 - t452 * t534) * r_i_i_C(2) - t443 * t482 + (t443 * t473 + t452 * t518) * pkin(5) - t542 * t442 + t483 * t451 (t478 * t512 - t443 * t474 + (-t452 * t478 + t474 * t527) * qJD(3)) * pkin(3) + t485, t485 (-t431 * t473 + t442 * t477 + (t446 * t477 - t451 * t473) * qJD(5)) * pkin(5) + t523, t523; 0 ((pkin(5) * t518 - qJD(2) * t542 + t497 * t468) * t475 + ((-t482 + t496 + t539) * qJD(2) - t483) * t479) * t472 (-t474 * t509 + (-t537 * t474 - t475 * t529) * qJD(3)) * pkin(3) + t484, t484 (t477 * t510 - t439 * t473 + (-t450 * t477 + t473 * t528) * qJD(5)) * pkin(5) + t522, t522;];
JaD_transl  = t1;
