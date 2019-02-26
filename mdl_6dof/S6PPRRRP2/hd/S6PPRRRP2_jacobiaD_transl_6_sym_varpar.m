% Zeitableitung der translatorischen Teilmatrix der analytischen Jacobi-Matrix für Segment Nr. 6 (0=Basis) von
% S6PPRRRP2
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
%   pkin=[a2,a3,a4,a5,a6,alpha2,alpha3,d3,d4,d5,theta1,theta2]';
%
% Output:
% JaD_transl [3x6]
%   Zeitableitung der translatorischen Teilmatrix der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-02-26 19:42
% Revision: d75aae1ac561373cd3be920984c3df29a1c2ecc8 (2019-02-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JaD_transl = S6PPRRRP2_jacobiaD_transl_6_sym_varpar(qJ, qJD, r_i_i_C, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(3,1),zeros(12,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PPRRRP2_jacobiaD_transl_6_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6PPRRRP2_jacobiaD_transl_6_sym_varpar: qJD has to be [6x1] (double)');
assert(isa(r_i_i_C,'double') && isreal(r_i_i_C) && all(size(r_i_i_C) == [3 1]), ...
	'S6PPRRRP2_jacobiaD_transl_6_sym_varpar: Position vector r_i_i_C has to be [3x1] double');
assert(isreal(pkin) && all(size(pkin) == [12 1]), ...
  'S6PPRRRP2_jacobiaD_transl_6_sym_varpar: pkin has to be [12x1] (double)');

%% Symbolic Calculation
% From jacobiaD_transl_6_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 19:42:09
% EndTime: 2019-02-26 19:42:09
% DurationCPUTime: 0.57s
% Computational Cost: add. (1160->104), mult. (3610->180), div. (0->0), fcn. (4263->14), ass. (0->81)
t538 = sin(pkin(12));
t539 = sin(pkin(11));
t521 = t539 * t538;
t542 = cos(pkin(12));
t543 = cos(pkin(11));
t528 = t543 * t542;
t545 = cos(pkin(6));
t498 = -t545 * t528 + t521;
t544 = cos(pkin(7));
t496 = t498 * t544;
t540 = sin(pkin(7));
t541 = sin(pkin(6));
t526 = t541 * t540;
t511 = t543 * t526;
t553 = t496 + t511;
t523 = t539 * t542;
t527 = t543 * t538;
t499 = t545 * t523 + t527;
t522 = t539 * t541;
t552 = t499 * t544 - t540 * t522;
t529 = t544 * t541;
t551 = t542 * t529 + t540 * t545;
t480 = t545 * t527 + t523;
t491 = sin(qJ(3));
t547 = cos(qJ(3));
t462 = t480 * t491 + t553 * t547;
t550 = t552 * t547;
t525 = t541 * t538;
t471 = t491 * t525 - t551 * t547;
t549 = pkin(10) + r_i_i_C(2);
t548 = -r_i_i_C(1) - pkin(5);
t546 = r_i_i_C(3) + qJ(6);
t489 = sin(qJ(5));
t493 = cos(qJ(4));
t537 = t489 * t493;
t536 = qJD(3) * t491;
t490 = sin(qJ(4));
t535 = qJD(4) * t490;
t534 = qJD(5) * t493;
t533 = t480 * t547;
t458 = t462 * qJD(3);
t532 = t462 * t534 - t458;
t481 = -t545 * t521 + t528;
t460 = qJD(3) * t550 + t481 * t536;
t464 = t481 * t491 + t550;
t531 = t464 * t534 - t460;
t469 = t471 * qJD(3);
t530 = t471 * t534 - t469;
t463 = -t553 * t491 + t533;
t473 = t498 * t540 - t543 * t529;
t453 = t463 * t493 + t473 * t490;
t492 = cos(qJ(5));
t520 = t453 * t492 + t462 * t489;
t465 = t481 * t547 - t552 * t491;
t474 = t499 * t540 + t544 * t522;
t455 = t465 * t493 + t474 * t490;
t519 = t455 * t492 + t464 * t489;
t518 = -t463 * t490 + t473 * t493;
t517 = -t465 * t490 + t474 * t493;
t472 = t551 * t491 + t547 * t525;
t479 = -t542 * t526 + t545 * t544;
t467 = t472 * t493 + t479 * t490;
t516 = t467 * t492 + t471 * t489;
t515 = -t472 * t490 + t479 * t493;
t513 = -t493 * pkin(4) - t549 * t490 - pkin(3);
t508 = qJD(4) * (pkin(4) * t490 - t549 * t493);
t507 = t546 * t489 - t548 * t492 + pkin(4);
t459 = -t511 * t536 + (-t491 * t496 + t533) * qJD(3);
t506 = qJD(5) * t463 - t459 * t493 + t462 * t535;
t461 = t465 * qJD(3);
t505 = qJD(5) * t465 - t461 * t493 + t464 * t535;
t470 = t472 * qJD(3);
t504 = qJD(5) * t472 - t470 * t493 + t471 * t535;
t501 = qJD(6) * t489 + (t548 * t489 + t546 * t492) * qJD(5);
t451 = t515 * qJD(4) - t469 * t493;
t449 = t517 * qJD(4) - t460 * t493;
t447 = t518 * qJD(4) - t458 * t493;
t442 = t516 * qJD(5) + t451 * t489 - t470 * t492;
t436 = t519 * qJD(5) + t449 * t489 - t461 * t492;
t434 = t520 * qJD(5) + t447 * t489 - t459 * t492;
t1 = [0, 0 -(t464 * t537 + t465 * t492) * qJD(6) - t460 * pkin(9) - t548 * (t531 * t489 + t505 * t492) + t546 * (t505 * t489 - t531 * t492) + t464 * t508 + t513 * t461, t549 * t449 + t501 * t517 + t507 * (-t455 * qJD(4) + t460 * t490) t519 * qJD(6) + t546 * (t449 * t492 + t461 * t489 + (-t455 * t489 + t464 * t492) * qJD(5)) + t548 * t436, t436; 0, 0 -(t462 * t537 + t463 * t492) * qJD(6) - t458 * pkin(9) - t548 * (t532 * t489 + t506 * t492) + t546 * (t506 * t489 - t532 * t492) + t462 * t508 + t513 * t459, t549 * t447 + t501 * t518 + t507 * (-t453 * qJD(4) + t458 * t490) t520 * qJD(6) + t546 * (t447 * t492 + t459 * t489 + (-t453 * t489 + t462 * t492) * qJD(5)) + t548 * t434, t434; 0, 0 -(t471 * t537 + t472 * t492) * qJD(6) - t469 * pkin(9) - t548 * (t530 * t489 + t504 * t492) + t546 * (t504 * t489 - t530 * t492) + t471 * t508 + t513 * t470, t549 * t451 + t501 * t515 + t507 * (-t467 * qJD(4) + t469 * t490) t516 * qJD(6) + t546 * (t451 * t492 + t470 * t489 + (-t467 * t489 + t471 * t492) * qJD(5)) + t548 * t442, t442;];
JaD_transl  = t1;
