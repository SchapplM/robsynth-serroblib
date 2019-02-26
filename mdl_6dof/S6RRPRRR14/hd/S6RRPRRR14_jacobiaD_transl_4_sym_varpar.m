% Zeitableitung der translatorischen Teilmatrix der analytischen Jacobi-Matrix für Segment Nr. 4 (0=Basis) von
% S6RRPRRR14
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
% pkin [14x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,alpha3,alpha4,d1,d2,d4,d5,d6,theta3]';
%
% Output:
% JaD_transl [3x6]
%   Zeitableitung der translatorischen Teilmatrix der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox (ehem. IRT-Maple-Toolbox)
% Datum: 2019-01-03 10:25
% Revision: 5fdbc45bcf2cc60deefd7ac2d71d743ed41bf7e4 (2018-12-21)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für mechatronische Systeme, Universität Hannover

function JaD_transl = S6RRPRRR14_jacobiaD_transl_4_floatb_twist_sym_varpar(qJ, qJD, r_i_i_C, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(3,1),zeros(14,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPRRR14_jacobiaD_transl_4_floatb_twist_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRPRRR14_jacobiaD_transl_4_floatb_twist_sym_varpar: qJD has to be [6x1] (double)');
assert(isa(r_i_i_C,'double') && isreal(r_i_i_C) && all(size(r_i_i_C) == [3 1]), ...
	'S6RRPRRR14_jacobiaD_transl_4_floatb_twist_sym_varpar: Position vector r_i_i_C has to be [3x1] double');
assert(isreal(pkin) && all(size(pkin) == [14 1]), ...
  'S6RRPRRR14_jacobiaD_transl_4_floatb_twist_sym_varpar: pkin has to be [14x1] (double)');

%% Symbolic Calculation
% From jacobiaD_transl_4_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-01-03 10:25:36
% EndTime: 2019-01-03 10:25:37
% DurationCPUTime: 1.04s
% Computational Cost: add. (651->136), mult. (2094->255), div. (0->0), fcn. (2168->14), ass. (0->91)
t496 = sin(qJ(1));
t495 = sin(qJ(2));
t548 = cos(pkin(6));
t519 = t496 * t548;
t515 = t495 * t519;
t526 = qJD(2) * t495;
t498 = cos(qJ(2));
t499 = cos(qJ(1));
t529 = t499 * t498;
t466 = -qJD(1) * t515 - t496 * t526 + (qJD(2) * t548 + qJD(1)) * t529;
t487 = sin(pkin(14));
t491 = cos(pkin(14));
t518 = t499 * t548;
t480 = t495 * t518 + t496 * t498;
t501 = t499 * t495 + t498 * t519;
t465 = t501 * qJD(1) + t480 * qJD(2);
t489 = sin(pkin(7));
t493 = cos(pkin(7));
t490 = sin(pkin(6));
t528 = qJD(1) * t490;
t522 = t496 * t528;
t503 = -t465 * t493 + t489 * t522;
t442 = -t466 * t487 + t503 * t491;
t443 = t466 * t491 + t503 * t487;
t546 = t465 * t489;
t458 = t493 * t522 + t546;
t494 = sin(qJ(4));
t492 = cos(pkin(8));
t497 = cos(qJ(4));
t532 = t492 * t497;
t488 = sin(pkin(8));
t539 = t488 * t497;
t557 = t442 * t532 - t443 * t494 + t458 * t539;
t533 = t492 * t494;
t540 = t488 * t494;
t556 = -t442 * t533 - t443 * t497 - t458 * t540;
t479 = t496 * t495 - t498 * t518;
t536 = t490 * t499;
t507 = t479 * t493 + t489 * t536;
t452 = t480 * t487 + t507 * t491;
t453 = -t480 * t491 + t507 * t487;
t469 = -t479 * t489 + t493 * t536;
t555 = t452 * t532 - t453 * t494 + t469 * t539;
t554 = t452 * t533 + t453 * t497 + t469 * t540;
t549 = pkin(11) + r_i_i_C(3);
t514 = t549 * t492 + qJ(3);
t502 = t515 - t529;
t463 = t479 * qJD(1) + t502 * qJD(2);
t547 = t463 * t489;
t541 = t487 * t493;
t538 = t489 * t490;
t537 = t490 * t496;
t535 = t491 * t493;
t534 = t491 * t495;
t531 = t493 * t495;
t530 = t493 * t498;
t525 = qJD(2) * t498;
t524 = qJD(4) * t494;
t523 = qJD(4) * t497;
t521 = t499 * t528;
t520 = qJ(3) * t493 + pkin(10);
t517 = t548 * t489;
t516 = t526 * t538;
t513 = t488 * t516;
t464 = t480 * qJD(1) + t501 * qJD(2);
t504 = t463 * t493 + t489 * t521;
t440 = t464 * t487 + t504 * t491;
t456 = t493 * t521 - t547;
t512 = t440 * t492 + t456 * t488;
t471 = t489 * t501 + t493 * t537;
t506 = t489 * t537 - t493 * t501;
t509 = (t487 * t502 + t506 * t491) * t492 + t471 * t488;
t508 = t497 * r_i_i_C(1) - t494 * r_i_i_C(2) + pkin(3);
t505 = -t487 * t495 + t491 * t530;
t476 = (-t487 * t498 - t491 * t531) * t490;
t477 = (-t487 * t531 + t491 * t498) * t490;
t500 = (t494 * r_i_i_C(1) + t497 * r_i_i_C(2)) * t492 - t549 * t488;
t478 = t548 * t493 - t498 * t538;
t475 = qJD(2) * t477;
t474 = qJD(2) * t476;
t468 = t490 * t534 + (t490 * t530 + t517) * t487;
t467 = t505 * t490 + t491 * t517;
t462 = -t491 * t501 + t502 * t541;
t461 = t487 * t501 + t502 * t535;
t460 = -t479 * t491 - t480 * t541;
t459 = t479 * t487 - t480 * t535;
t455 = t506 * t487 - t491 * t502;
t441 = -t464 * t491 + t504 * t487;
t439 = t441 * t497 + t512 * t494 + (-t455 * t494 + t509 * t497) * qJD(4);
t438 = -t441 * t494 + t512 * t497 + (-t455 * t497 - t509 * t494) * qJD(4);
t1 = [t556 * r_i_i_C(1) - t557 * r_i_i_C(2) - t443 * pkin(3) - t466 * pkin(2) - qJ(3) * t546 + t469 * qJD(3) + (-t499 * pkin(1) - t520 * t537) * qJD(1) + (t555 * r_i_i_C(1) - t554 * r_i_i_C(2)) * qJD(4) + t549 * (t442 * t488 - t458 * t492) t463 * pkin(2) + t508 * (t463 * t491 + t464 * t541) + t500 * (-t463 * t487 + t464 * t535) + ((t461 * t532 - t462 * t494) * r_i_i_C(1) + (-t461 * t533 - t462 * t497) * r_i_i_C(2)) * qJD(4) + (-t502 * qJD(3) - t514 * t464 + ((-t464 * t494 - t502 * t523) * r_i_i_C(1) + (-t464 * t497 + t502 * t524) * r_i_i_C(2)) * t488) * t489, t456, r_i_i_C(1) * t438 - t439 * r_i_i_C(2), 0, 0; t439 * r_i_i_C(1) + t438 * r_i_i_C(2) + t441 * pkin(3) - t464 * pkin(2) - qJ(3) * t547 + t471 * qJD(3) + (-t496 * pkin(1) + t520 * t536) * qJD(1) + t549 * (-t440 * t488 + t456 * t492) -t465 * pkin(2) + t508 * (-t465 * t491 - t466 * t541) + t500 * (t465 * t487 - t466 * t535) + ((t459 * t532 - t460 * t494) * r_i_i_C(1) + (-t459 * t533 - t460 * t497) * r_i_i_C(2)) * qJD(4) + (t480 * qJD(3) + t514 * t466 + ((t466 * t494 + t480 * t523) * r_i_i_C(1) + (t466 * t497 - t480 * t524) * r_i_i_C(2)) * t488) * t489, t458, t557 * r_i_i_C(1) + t556 * r_i_i_C(2) + (t554 * r_i_i_C(1) + t555 * r_i_i_C(2)) * qJD(4), 0, 0; 0 ((t476 * t532 - t477 * t494) * r_i_i_C(1) + (-t476 * t533 - t477 * t497) * r_i_i_C(2)) * qJD(4) + (-pkin(2) * t526 + (t495 * qJD(3) + t514 * t525 + ((t494 * t525 + t495 * t523) * r_i_i_C(1) + (-t495 * t524 + t497 * t525) * r_i_i_C(2)) * t488) * t489 + (t508 * (-t487 * t530 - t534) - t500 * t505) * qJD(2)) * t490, t516 (t474 * t532 - t475 * t494 + t497 * t513) * r_i_i_C(1) + (-t474 * t533 - t475 * t497 - t494 * t513) * r_i_i_C(2) + ((-t467 * t533 - t468 * t497 - t478 * t540) * r_i_i_C(1) + (-t467 * t532 + t468 * t494 - t478 * t539) * r_i_i_C(2)) * qJD(4), 0, 0;];
JaD_transl  = t1;
