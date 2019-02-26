% Zeitableitung der translatorischen Teilmatrix der analytischen Jacobi-Matrix für Segment Nr. 5 (0=Basis) von
% S7RRRRRRR1
% Use Code from Maple symbolic Code Generation
%
% analytische Jacobi-Matrix: Differentieller Zusammenhang zwischen
% Endeffektorposition und verallgemeinerten Koordinaten.
%
% Input:
% qJ [7x1]
%   Generalized joint coordinates (joint angles)
% qJD [7x1]
%   Generalized joint velocities
% r_i_i_C [3x1]
%   Ortsvektor vom KörperKS-Ursprung zum gesuchten Punkt
% pkin [4x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[d1,d3,d5,d7]';
%
% Output:
% JaD_transl [3x7]
%   Zeitableitung der translatorischen Teilmatrix der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-02-26 22:54
% Revision: d75aae1ac561373cd3be920984c3df29a1c2ecc8 (2019-02-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JaD_transl = S7RRRRRRR1_jacobiaD_transl_5_sym_varpar(qJ, qJD, r_i_i_C, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(7,1),zeros(7,1),zeros(3,1),zeros(4,1)}
assert(isreal(qJ) && all(size(qJ) == [7 1]), ...
  'S7RRRRRRR1_jacobiaD_transl_5_sym_varpar: qJ has to be [7x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [7 1]), ...
  'S7RRRRRRR1_jacobiaD_transl_5_sym_varpar: qJD has to be [7x1] (double)');
assert(isa(r_i_i_C,'double') && isreal(r_i_i_C) && all(size(r_i_i_C) == [3 1]), ...
	'S7RRRRRRR1_jacobiaD_transl_5_sym_varpar: Position vector r_i_i_C has to be [3x1] double');
assert(isreal(pkin) && all(size(pkin) == [4 1]), ...
  'S7RRRRRRR1_jacobiaD_transl_5_sym_varpar: pkin has to be [4x1] (double)');

%% Symbolic Calculation
% From jacobiaD_transl_5_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 22:54:23
% EndTime: 2019-02-26 22:54:25
% DurationCPUTime: 1.16s
% Computational Cost: add. (553->140), mult. (1736->256), div. (0->0), fcn. (1741->10), ass. (0->96)
t460 = sin(qJ(4));
t465 = cos(qJ(4));
t466 = cos(qJ(3));
t461 = sin(qJ(3));
t468 = cos(qJ(1));
t522 = t468 * t461;
t463 = sin(qJ(1));
t467 = cos(qJ(2));
t524 = t463 * t467;
t448 = t466 * t524 + t522;
t462 = sin(qJ(2));
t516 = qJD(2) * t467;
t518 = qJD(1) * t468;
t479 = t462 * t518 + t463 * t516;
t475 = -qJD(4) * t448 + t479;
t519 = qJD(1) * t467;
t495 = qJD(3) + t519;
t512 = qJD(3) * t467;
t497 = t461 * t512;
t517 = qJD(2) * t462;
t501 = t463 * t517;
t520 = qJD(1) * t463;
t503 = t461 * t520;
t521 = t468 * t466;
t437 = -t463 * t497 - t466 * t501 + t495 * t521 - t503;
t511 = qJD(4) * t462;
t487 = t463 * t511 + t437;
t423 = t475 * t460 + t487 * t465;
t480 = t463 * t466 + t467 * t522;
t436 = t480 * qJD(1) + t448 * qJD(3) - t461 * t501;
t459 = sin(qJ(5));
t464 = cos(qJ(5));
t551 = t423 * t459 + t436 * t464;
t550 = -t423 * t464 + t436 * t459;
t494 = qJD(4) * t466 - qJD(2);
t514 = qJD(3) * t461;
t527 = t462 * t465;
t493 = qJD(2) * t466 - qJD(4);
t541 = t493 * t467;
t549 = t460 * (-t462 * t514 + t541) + t494 * t527;
t489 = t464 * r_i_i_C(1) - t459 * r_i_i_C(2);
t536 = r_i_i_C(3) + pkin(3);
t471 = (t489 * t460 - t536 * t465) * qJD(4);
t488 = t459 * r_i_i_C(1) + t464 * r_i_i_C(2);
t507 = qJD(5) * t465;
t548 = t488 * t507 + t471;
t505 = t536 * t460;
t547 = t489 * t465 + t505;
t528 = t462 * t460;
t439 = t448 * t465 + t463 * t528;
t447 = t461 * t524 - t521;
t545 = t439 * t459 + t447 * t464;
t544 = t439 * t464 - t447 * t459;
t470 = (t494 * t460 + t465 * t514) * t462 - t465 * t541;
t515 = qJD(2) * t468;
t478 = t462 * t520 - t467 * t515;
t538 = t493 * t462 + t497;
t529 = t461 * t465;
t526 = t462 * t466;
t525 = t462 * t468;
t523 = t467 * t465;
t513 = qJD(3) * t466;
t510 = qJD(5) * t459;
t509 = qJD(5) * t461;
t508 = qJD(5) * t464;
t506 = pkin(2) * t519;
t502 = t461 * t516;
t496 = qJD(1) + t512;
t491 = t459 * t502;
t490 = t464 * t502;
t484 = t494 * t467;
t486 = t460 * t484 + t538 * t465 + t467 * t509;
t485 = t462 * t509 + t470;
t451 = -t463 * t461 + t467 * t521;
t442 = t451 * t465 + t460 * t525;
t446 = -t467 * t460 + t465 * t526;
t481 = t460 * t526 + t523;
t444 = t446 * t468;
t477 = qJD(1) * t481;
t476 = t462 * t515 + t495 * t463;
t474 = qJD(5) * t446 + t462 * t513 + t502;
t473 = -qJD(5) * (t466 * t523 + t528) + t461 * t517 - t466 * t512;
t422 = -t487 * t460 + t475 * t465;
t469 = (t459 * t513 + t461 * t508) * r_i_i_C(1) + (-t459 * t509 + t464 * t513) * r_i_i_C(2) + qJD(2) * pkin(2);
t443 = t446 * t463;
t441 = -t451 * t460 + t465 * t525;
t438 = -t448 * t460 + t463 * t527;
t435 = t476 * t466 + t496 * t522;
t434 = t476 * t461 - t496 * t521;
t429 = -qJD(1) * t444 + t470 * t463;
t427 = t446 * t520 + t470 * t468;
t421 = (t468 * t511 - t435) * t465 + (-qJD(4) * t451 - t478) * t460;
t420 = t442 * qJD(4) - t435 * t460 + t478 * t465;
t419 = t421 * t464 + t434 * t459 + (-t442 * t459 - t464 * t480) * qJD(5);
t418 = -t421 * t459 + t434 * t464 + (-t442 * t464 + t459 * t480) * qJD(5);
t1 = [t550 * r_i_i_C(1) + t551 * r_i_i_C(2) + t536 * t422 + (r_i_i_C(1) * t545 + r_i_i_C(2) * t544) * qJD(5) + t479 * pkin(2) (t427 * t464 + t444 * t510 + t468 * t491) * r_i_i_C(1) + (-t427 * t459 + t444 * t508 + t468 * t490) * r_i_i_C(2) + t463 * t506 + t536 * (t463 * t477 - t468 * t549) + (t469 * t468 - t488 * t503) * t462 (t435 * t459 - t451 * t508) * r_i_i_C(1) + (t435 * t464 + t451 * t510) * r_i_i_C(2) + t547 * t434 + t548 * t480, t536 * t421 + (t420 * t459 - t441 * t508) * r_i_i_C(2) + (-t420 * t464 - t441 * t510) * r_i_i_C(1), t418 * r_i_i_C(1) - t419 * r_i_i_C(2), 0, 0; t478 * pkin(2) + t419 * r_i_i_C(1) + t418 * r_i_i_C(2) + t536 * t420 (t429 * t464 + t443 * t510 + t463 * t491) * r_i_i_C(1) + (-t429 * t459 + t443 * t508 + t463 * t490) * r_i_i_C(2) - t468 * t506 - t536 * (t549 * t463 + t468 * t477) + (t488 * t461 * t518 + t469 * t463) * t462 (-t437 * t459 - t448 * t508) * r_i_i_C(1) + (-t437 * t464 + t448 * t510) * r_i_i_C(2) - t547 * t436 + t548 * t447, t536 * t423 + (-t422 * t459 - t438 * t508) * r_i_i_C(2) + (t422 * t464 - t438 * t510) * r_i_i_C(1), -t551 * r_i_i_C(1) + t550 * r_i_i_C(2) + (-r_i_i_C(1) * t544 + r_i_i_C(2) * t545) * qJD(5), 0, 0; 0, -pkin(2) * t516 - t536 * (t538 * t460 - t465 * t484) + (-t486 * r_i_i_C(1) + t473 * r_i_i_C(2)) * t464 + (t473 * r_i_i_C(1) + t486 * r_i_i_C(2)) * t459 ((-t459 * t466 - t464 * t529) * r_i_i_C(1) + (t459 * t529 - t464 * t466) * r_i_i_C(2) - t461 * t505) * t516 + ((-qJD(3) * t547 - t489 * qJD(5)) * t466 + (t471 + t488 * (qJD(3) + t507)) * t461) * t462, -t536 * t470 + (t459 * t549 + t481 * t508) * r_i_i_C(2) + (-t464 * t549 + t481 * t510) * r_i_i_C(1) (-t474 * r_i_i_C(1) + t485 * r_i_i_C(2)) * t464 + (t485 * r_i_i_C(1) + t474 * r_i_i_C(2)) * t459, 0, 0;];
JaD_transl  = t1;
