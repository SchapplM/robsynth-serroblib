% Zeitableitung der translatorischen Teilmatrix der analytischen Jacobi-Matrix für Segment Nr. 6 (0=Basis) von
% S6RRRRPP7
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
%   pkin=[a2,a3,a4,a5,a6,alpha2,d1,d2,d3,d4,theta5]';
%
% Output:
% JaD_transl [3x6]
%   Zeitableitung der translatorischen Teilmatrix der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-02-26 22:29
% Revision: d75aae1ac561373cd3be920984c3df29a1c2ecc8 (2019-02-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JaD_transl = S6RRRRPP7_jacobiaD_transl_6_sym_varpar(qJ, qJD, r_i_i_C, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(3,1),zeros(11,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRRPP7_jacobiaD_transl_6_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRRRPP7_jacobiaD_transl_6_sym_varpar: qJD has to be [6x1] (double)');
assert(isa(r_i_i_C,'double') && isreal(r_i_i_C) && all(size(r_i_i_C) == [3 1]), ...
	'S6RRRRPP7_jacobiaD_transl_6_sym_varpar: Position vector r_i_i_C has to be [3x1] double');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6RRRRPP7_jacobiaD_transl_6_sym_varpar: pkin has to be [11x1] (double)');

%% Symbolic Calculation
% From jacobiaD_transl_6_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 22:28:58
% EndTime: 2019-02-26 22:28:59
% DurationCPUTime: 0.81s
% Computational Cost: add. (1131->149), mult. (2714->231), div. (0->0), fcn. (2795->12), ass. (0->81)
t474 = sin(qJ(3));
t478 = cos(qJ(3));
t475 = sin(qJ(2));
t476 = sin(qJ(1));
t479 = cos(qJ(2));
t531 = cos(pkin(6));
t535 = cos(qJ(1));
t498 = t531 * t535;
t452 = t475 * t498 + t476 * t479;
t471 = sin(pkin(6));
t525 = t471 * t476;
t544 = qJD(1) * t525 - qJD(3) * t452;
t505 = t476 * t531;
t499 = t475 * t505;
t507 = t535 * qJD(1);
t521 = qJD(2) * t475;
t441 = -qJD(1) * t499 - t476 * t521 + (qJD(2) * t498 + t507) * t479;
t513 = t471 * t535;
t545 = -qJD(3) * t513 + t441;
t431 = t544 * t474 + t545 * t478;
t453 = t535 * t475 + t479 * t505;
t440 = t453 * qJD(1) + t452 * qJD(2);
t470 = qJ(4) + pkin(11);
t468 = sin(t470);
t469 = cos(t470);
t445 = t452 * t478 - t474 * t513;
t490 = t479 * t498;
t522 = t476 * t475;
t451 = -t490 + t522;
t495 = t445 * t469 + t451 * t468;
t418 = t495 * qJD(4) + t431 * t468 - t440 * t469;
t494 = t445 * t468 - t451 * t469;
t546 = -t494 * qJD(4) + t431 * t469 + t440 * t468;
t477 = cos(qJ(4));
t467 = t477 * pkin(4) + pkin(3);
t516 = t468 * qJD(6);
t533 = r_i_i_C(2) + qJ(5) + pkin(10);
t473 = sin(qJ(4));
t534 = t473 * pkin(4);
t480 = -(t533 * qJD(3) - qJD(4) * t534 + t516) * t478 + (qJD(3) * t467 - qJD(5)) * t474;
t524 = t471 * t478;
t450 = t531 * t474 + t475 * t524;
t523 = t471 * t479;
t540 = -t450 * t469 + t468 * t523;
t539 = t478 * t467 + t533 * t474 + pkin(2);
t519 = qJD(3) * t474;
t538 = (qJD(2) * t478 - qJD(4)) * t475 + t479 * t519;
t536 = r_i_i_C(1) + pkin(5);
t532 = r_i_i_C(3) + qJ(6);
t518 = qJD(4) * t477;
t517 = qJD(4) * t478;
t515 = t469 * qJD(6);
t512 = t535 * t479;
t510 = t471 * t521;
t509 = qJD(2) * t523;
t500 = t471 * t507;
t439 = t452 * qJD(1) + t453 * qJD(2);
t497 = t453 * t517 - t439;
t496 = t451 * t517 + t441;
t454 = t512 - t499;
t448 = t454 * t478 + t474 * t525;
t493 = t448 * t469 + t453 * t468;
t492 = -t448 * t468 + t453 * t469;
t491 = (qJD(2) - t517) * t479;
t447 = -t454 * t474 + t476 * t524;
t486 = t452 * t474 + t478 * t513;
t485 = -t471 * t475 * t474 + t531 * t478;
t430 = t545 * t474 - t544 * t478;
t484 = -t532 * t468 - t536 * t469 - t467;
t438 = -qJD(1) * t490 - qJD(2) * t512 + (qJD(2) * t531 + qJD(1)) * t522;
t483 = qJD(4) * t454 + t438 * t478 + t453 * t519;
t482 = qJD(4) * t452 - t440 * t478 + t451 * t519;
t481 = t516 + (-t536 * t468 + t532 * t469 - t534) * qJD(4);
t443 = t485 * qJD(3) + t478 * t509;
t442 = t450 * qJD(3) + t474 * t509;
t429 = t447 * qJD(3) - t439 * t478 + t474 * t500;
t428 = t448 * qJD(3) - t439 * t474 - t478 * t500;
t426 = -t540 * qJD(4) + t443 * t468 - t469 * t510;
t417 = t492 * qJD(4) + t429 * t469 - t438 * t468;
t416 = t493 * qJD(4) + t429 * t468 + t438 * t469;
t1 = [-t494 * qJD(6) - t431 * t467 - t486 * qJD(5) - t441 * pkin(2) - t440 * pkin(9) - t533 * t430 - t536 * t546 - t532 * t418 + (-t535 * pkin(1) - pkin(8) * t525) * qJD(1) + (-t440 * t473 + (t445 * t473 - t451 * t477) * qJD(4)) * pkin(4), -t454 * t515 - t439 * pkin(9) + t536 * (t497 * t468 + t483 * t469) + t532 * (t483 * t468 - t497 * t469) + (-t439 * t473 + t454 * t518) * pkin(4) + t539 * t438 + t480 * t453, t448 * qJD(5) + t484 * t428 + t533 * t429 + t481 * t447, t493 * qJD(6) + t532 * t417 - t536 * t416 + (-t429 * t473 - t438 * t477 + (-t448 * t477 - t453 * t473) * qJD(4)) * pkin(4), t428, t416; -t492 * qJD(6) + t429 * t467 - t447 * qJD(5) - t439 * pkin(2) - t438 * pkin(9) + t533 * t428 + t536 * t417 + t532 * t416 + (-t476 * pkin(1) + pkin(8) * t513) * qJD(1) + (-t438 * t473 + (-t448 * t473 + t453 * t477) * qJD(4)) * pkin(4), -t452 * t515 + t441 * pkin(9) + t536 * (t496 * t468 + t482 * t469) + t532 * (t482 * t468 - t496 * t469) + (t441 * t473 + t452 * t518) * pkin(4) - t539 * t440 + t480 * t451, t445 * qJD(5) + t484 * t430 + t533 * t431 - t481 * t486, t495 * qJD(6) + t532 * t546 - t536 * t418 + (-t431 * t473 + t440 * t477 + (-t445 * t477 - t451 * t473) * qJD(4)) * pkin(4), t430, t418; 0 (t536 * (t468 * t491 - t538 * t469) - t532 * (t538 * t468 + t469 * t491) + (pkin(4) * t518 - qJD(2) * t539 - t515) * t475 + ((pkin(9) + t534) * qJD(2) - t480) * t479) * t471, t450 * qJD(5) + t484 * t442 + t533 * t443 + t481 * t485, -t540 * qJD(6) + t532 * (t468 * t510 + t443 * t469 + (-t450 * t468 - t469 * t523) * qJD(4)) - t536 * t426 + (t477 * t510 - t443 * t473 + (-t450 * t477 + t473 * t523) * qJD(4)) * pkin(4), t442, t426;];
JaD_transl  = t1;
