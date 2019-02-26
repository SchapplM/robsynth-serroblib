% Zeitableitung der translatorischen Teilmatrix der analytischen Jacobi-Matrix für Segment Nr. 4 (0=Basis) von
% S6RPRRRR12
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
%   pkin=[a2,a3,a4,a5,a6,alpha2,alpha3,alpha4,d1,d3,d4,d5,d6,theta2]';
%
% Output:
% JaD_transl [3x6]
%   Zeitableitung der translatorischen Teilmatrix der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-02-26 21:21
% Revision: d75aae1ac561373cd3be920984c3df29a1c2ecc8 (2019-02-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JaD_transl = S6RPRRRR12_jacobiaD_transl_4_sym_varpar(qJ, qJD, r_i_i_C, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(3,1),zeros(14,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRRRR12_jacobiaD_transl_4_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RPRRRR12_jacobiaD_transl_4_sym_varpar: qJD has to be [6x1] (double)');
assert(isa(r_i_i_C,'double') && isreal(r_i_i_C) && all(size(r_i_i_C) == [3 1]), ...
	'S6RPRRRR12_jacobiaD_transl_4_sym_varpar: Position vector r_i_i_C has to be [3x1] double');
assert(isreal(pkin) && all(size(pkin) == [14 1]), ...
  'S6RPRRRR12_jacobiaD_transl_4_sym_varpar: pkin has to be [14x1] (double)');

%% Symbolic Calculation
% From jacobiaD_transl_4_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 21:21:12
% EndTime: 2019-02-26 21:21:13
% DurationCPUTime: 0.86s
% Computational Cost: add. (609->114), mult. (2041->203), div. (0->0), fcn. (2200->14), ass. (0->86)
t489 = cos(pkin(14));
t492 = cos(pkin(6));
t485 = sin(pkin(14));
t495 = sin(qJ(1));
t531 = t495 * t485;
t520 = t492 * t531;
t498 = cos(qJ(1));
t527 = qJD(1) * t498;
t474 = -qJD(1) * t520 + t489 * t527;
t494 = sin(qJ(3));
t497 = cos(qJ(3));
t529 = t498 * t485;
t530 = t495 * t489;
t505 = t492 * t530 + t529;
t473 = t505 * qJD(1);
t491 = cos(pkin(7));
t487 = sin(pkin(7));
t488 = sin(pkin(6));
t537 = t488 * t495;
t519 = qJD(1) * t537;
t515 = t487 * t519;
t502 = -t473 * t491 + t515;
t536 = t488 * t498;
t523 = t487 * t536;
t528 = t498 * t489;
t476 = -t492 * t528 + t531;
t543 = t476 * t491;
t508 = t523 + t543;
t477 = t492 * t529 + t530;
t541 = t477 * t497;
t449 = (t508 * t494 - t541) * qJD(3) - t474 * t494 + t502 * t497;
t545 = t473 * t487;
t463 = t491 * t519 + t545;
t490 = cos(pkin(8));
t496 = cos(qJ(4));
t533 = t490 * t496;
t486 = sin(pkin(8));
t539 = t486 * t496;
t557 = t449 * t533 + t463 * t539;
t493 = sin(qJ(4));
t534 = t490 * t493;
t540 = t486 * t493;
t556 = -t449 * t534 - t463 * t540;
t455 = t477 * t494 + t508 * t497;
t532 = t491 * t494;
t456 = t476 * t532 + t494 * t523 - t541;
t466 = -t476 * t487 + t491 * t536;
t555 = t455 * t533 - t456 * t493 + t466 * t539;
t554 = t455 * t534 + t456 * t496 + t466 * t540;
t547 = pkin(11) + r_i_i_C(3);
t538 = t487 * t492;
t465 = (t485 * t497 + t489 * t532) * t488 + t494 * t538;
t504 = t520 - t528;
t507 = t487 * t537 - t491 * t505;
t458 = t507 * t494 - t504 * t497;
t499 = t547 * t486 - (t493 * r_i_i_C(1) + t496 * r_i_i_C(2)) * t490;
t471 = t476 * qJD(1);
t546 = t471 * t487;
t535 = t489 * t491;
t526 = qJD(3) * t494;
t525 = qJD(3) * t497;
t524 = t488 * qJD(2);
t521 = t497 * t538;
t518 = t488 * t527;
t517 = t488 * t525;
t516 = pkin(10) * t491 + qJ(2);
t472 = t477 * qJD(1);
t503 = t471 * t491 + t487 * t518;
t447 = -t458 * qJD(3) + t472 * t494 + t503 * t497;
t461 = t491 * t518 - t546;
t513 = t447 * t490 + t461 * t486;
t501 = t504 * t494;
t457 = t507 * t497 + t501;
t510 = t457 * t490 + (t487 * t505 + t491 * t537) * t486;
t509 = t496 * r_i_i_C(1) - t493 * r_i_i_C(2) + pkin(3);
t480 = t498 * t487 * t517;
t475 = -t488 * t489 * t487 + t492 * t491;
t464 = t521 + (-t485 * t494 + t497 * t535) * t488;
t460 = t465 * qJD(3);
t459 = t488 * t485 * t526 - qJD(3) * t521 - t517 * t535;
t452 = t480 + (qJD(3) * t543 - t474) * t497 + (qJD(3) * t477 - t502) * t494;
t450 = t494 * t515 + t474 * t497 - t477 * t526 - t480 + (-t473 * t494 - t476 * t525) * t491;
t448 = qJD(3) * t501 + t503 * t494 + (t507 * qJD(3) - t472) * t497;
t446 = t448 * t496 + t513 * t493 + (-t458 * t493 + t510 * t496) * qJD(4);
t445 = -t448 * t493 + t513 * t496 + (-t458 * t496 - t510 * t493) * qJD(4);
t1 = [(t452 * t496 + t556) * r_i_i_C(1) + (-t452 * t493 - t557) * r_i_i_C(2) + t452 * pkin(3) - t474 * pkin(2) - pkin(10) * t545 + t498 * t524 + (-t498 * pkin(1) - t516 * t537) * qJD(1) + (t555 * r_i_i_C(1) - t554 * r_i_i_C(2)) * qJD(4) + t547 * (t449 * t486 - t463 * t490) t518, t509 * t447 + t499 * t448 + ((-t457 * t493 - t458 * t533) * r_i_i_C(1) + (-t457 * t496 + t458 * t534) * r_i_i_C(2)) * qJD(4), t445 * r_i_i_C(1) - t446 * r_i_i_C(2), 0, 0; t446 * r_i_i_C(1) + t445 * r_i_i_C(2) + t448 * pkin(3) - t472 * pkin(2) - pkin(10) * t546 + t495 * t524 + (-t495 * pkin(1) + t516 * t536) * qJD(1) + t547 * (-t447 * t486 + t461 * t490) t519, t509 * t449 + t499 * t450 + ((t455 * t493 + t456 * t533) * r_i_i_C(1) + (t455 * t496 - t456 * t534) * r_i_i_C(2)) * qJD(4) (-t450 * t493 + t557) * r_i_i_C(1) + (-t450 * t496 + t556) * r_i_i_C(2) + (t554 * r_i_i_C(1) + t555 * r_i_i_C(2)) * qJD(4), 0, 0; 0, 0, -t509 * t460 - t499 * t459 + ((-t464 * t493 - t465 * t533) * r_i_i_C(1) + (-t464 * t496 + t465 * t534) * r_i_i_C(2)) * qJD(4) (t459 * t493 - t460 * t533) * r_i_i_C(1) + (t459 * t496 + t460 * t534) * r_i_i_C(2) + ((-t464 * t534 - t465 * t496 - t475 * t540) * r_i_i_C(1) + (-t464 * t533 + t465 * t493 - t475 * t539) * r_i_i_C(2)) * qJD(4), 0, 0;];
JaD_transl  = t1;
