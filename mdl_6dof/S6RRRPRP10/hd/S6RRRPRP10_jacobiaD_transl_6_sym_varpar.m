% Zeitableitung der translatorischen Teilmatrix der analytischen Jacobi-Matrix für Segment Nr. 6 (0=Basis) von
% S6RRRPRP10
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
%   pkin=[a2,a3,a4,a5,a6,alpha2,d1,d2,d3,d5,theta4]';
%
% Output:
% JaD_transl [3x6]
%   Zeitableitung der translatorischen Teilmatrix der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-02-26 22:14
% Revision: d75aae1ac561373cd3be920984c3df29a1c2ecc8 (2019-02-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JaD_transl = S6RRRPRP10_jacobiaD_transl_6_sym_varpar(qJ, qJD, r_i_i_C, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(3,1),zeros(11,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRPRP10_jacobiaD_transl_6_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRRPRP10_jacobiaD_transl_6_sym_varpar: qJD has to be [6x1] (double)');
assert(isa(r_i_i_C,'double') && isreal(r_i_i_C) && all(size(r_i_i_C) == [3 1]), ...
	'S6RRRPRP10_jacobiaD_transl_6_sym_varpar: Position vector r_i_i_C has to be [3x1] double');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6RRRPRP10_jacobiaD_transl_6_sym_varpar: pkin has to be [11x1] (double)');

%% Symbolic Calculation
% From jacobiaD_transl_6_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 22:14:27
% EndTime: 2019-02-26 22:14:28
% DurationCPUTime: 0.64s
% Computational Cost: add. (1062->122), mult. (2499->192), div. (0->0), fcn. (2581->12), ass. (0->77)
t472 = sin(qJ(3));
t475 = cos(qJ(3));
t470 = cos(pkin(6));
t473 = sin(qJ(2));
t477 = cos(qJ(1));
t516 = t477 * t473;
t474 = sin(qJ(1));
t476 = cos(qJ(2));
t517 = t474 * t476;
t451 = t470 * t516 + t517;
t469 = sin(pkin(6));
t514 = qJD(1) * t469;
t539 = -qJD(3) * t451 + t474 * t514;
t518 = t474 * t473;
t506 = t470 * t518;
t513 = qJD(2) * t473;
t515 = t477 * t476;
t440 = -qJD(1) * t506 - t474 * t513 + (qJD(2) * t470 + qJD(1)) * t515;
t519 = t469 * t477;
t540 = -qJD(3) * t519 + t440;
t430 = t472 * t539 + t540 * t475;
t452 = t470 * t517 + t516;
t439 = qJD(1) * t452 + qJD(2) * t451;
t467 = pkin(11) + qJ(5);
t465 = sin(t467);
t466 = cos(t467);
t444 = t451 * t475 - t472 * t519;
t450 = -t470 * t515 + t518;
t493 = t444 * t466 + t450 * t465;
t417 = qJD(5) * t493 + t430 * t465 - t439 * t466;
t492 = t444 * t465 - t450 * t466;
t541 = -qJD(5) * t492 + t430 * t466 + t439 * t465;
t464 = cos(pkin(11)) * pkin(4) + pkin(3);
t509 = qJD(6) * t465;
t530 = r_i_i_C(2) + pkin(10) + qJ(4);
t536 = (-t464 * t472 + t475 * t530) * qJD(3) + t472 * qJD(4) + t475 * t509;
t521 = t469 * t475;
t449 = t470 * t472 + t473 * t521;
t520 = t469 * t476;
t535 = -t449 * t466 + t465 * t520;
t534 = t475 * t464 + t472 * t530 + pkin(2);
t511 = qJD(3) * t472;
t533 = (qJD(2) * t475 - qJD(5)) * t473 + t476 * t511;
t531 = r_i_i_C(1) + pkin(5);
t529 = r_i_i_C(3) + qJ(6);
t522 = t469 * t474;
t510 = qJD(5) * t475;
t508 = t466 * qJD(6);
t505 = sin(pkin(11)) * pkin(4) + pkin(9);
t503 = t477 * t514;
t502 = t469 * t513;
t501 = qJD(2) * t520;
t438 = qJD(1) * t451 + qJD(2) * t452;
t495 = t452 * t510 - t438;
t494 = t450 * t510 + t440;
t486 = t506 - t515;
t447 = t472 * t522 - t475 * t486;
t491 = t447 * t466 + t452 * t465;
t490 = -t447 * t465 + t452 * t466;
t489 = (qJD(2) - t510) * t476;
t488 = t451 * t472 + t475 * t519;
t446 = t472 * t486 + t474 * t521;
t487 = -t469 * t473 * t472 + t470 * t475;
t429 = t540 * t472 - t539 * t475;
t482 = -t529 * t465 - t466 * t531 - t464;
t437 = qJD(1) * t450 + qJD(2) * t486;
t481 = -qJD(5) * t486 + t437 * t475 + t452 * t511;
t480 = qJD(5) * t451 - t439 * t475 + t450 * t511;
t479 = t509 + (-t465 * t531 + t529 * t466) * qJD(5);
t442 = qJD(3) * t487 + t475 * t501;
t441 = qJD(3) * t449 + t472 * t501;
t428 = qJD(3) * t446 - t438 * t475 + t472 * t503;
t427 = qJD(3) * t447 - t438 * t472 - t475 * t503;
t425 = -t535 * qJD(5) + t442 * t465 - t466 * t502;
t416 = qJD(5) * t490 + t428 * t466 - t437 * t465;
t415 = qJD(5) * t491 + t428 * t465 + t437 * t466;
t1 = [-t492 * qJD(6) - t430 * t464 - t488 * qJD(4) - t440 * pkin(2) - t505 * t439 - t530 * t429 - t531 * t541 - t529 * t417 + (-t477 * pkin(1) - pkin(8) * t522) * qJD(1), t486 * t508 - t505 * t438 + t531 * (t465 * t495 + t466 * t481) + t529 * (t465 * t481 - t466 * t495) - t536 * t452 + t534 * t437, t447 * qJD(4) + t427 * t482 + t428 * t530 + t446 * t479, t427, t491 * qJD(6) - t415 * t531 + t529 * t416, t415; -t490 * qJD(6) + t428 * t464 - t446 * qJD(4) - t438 * pkin(2) - t505 * t437 + t530 * t427 + t531 * t416 + t529 * t415 + (-t474 * pkin(1) + pkin(8) * t519) * qJD(1), -t451 * t508 + t505 * t440 + t531 * (t465 * t494 + t466 * t480) + t529 * (t465 * t480 - t466 * t494) - t536 * t450 - t534 * t439, t444 * qJD(4) + t429 * t482 + t430 * t530 - t479 * t488, t429, t493 * qJD(6) - t531 * t417 + t529 * t541, t417; 0 (t531 * (t465 * t489 - t533 * t466) - t529 * (t533 * t465 + t466 * t489) - t473 * t508 + t536 * t476 + (-t473 * t534 + t476 * t505) * qJD(2)) * t469, t449 * qJD(4) + t441 * t482 + t442 * t530 + t479 * t487, t441, -t535 * qJD(6) + t529 * (t465 * t502 + t442 * t466 + (-t449 * t465 - t466 * t520) * qJD(5)) - t531 * t425, t425;];
JaD_transl  = t1;
