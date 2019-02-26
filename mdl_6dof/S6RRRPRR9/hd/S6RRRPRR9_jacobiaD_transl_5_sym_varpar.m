% Zeitableitung der translatorischen Teilmatrix der analytischen Jacobi-Matrix für Segment Nr. 5 (0=Basis) von
% S6RRRPRR9
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
%   pkin=[a2,a3,a4,a5,a6,alpha2,alpha3,d1,d2,d3,d5,d6,theta4]';
%
% Output:
% JaD_transl [3x6]
%   Zeitableitung der translatorischen Teilmatrix der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-02-26 22:20
% Revision: d75aae1ac561373cd3be920984c3df29a1c2ecc8 (2019-02-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JaD_transl = S6RRRPRR9_jacobiaD_transl_5_sym_varpar(qJ, qJD, r_i_i_C, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(3,1),zeros(13,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRPRR9_jacobiaD_transl_5_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRRPRR9_jacobiaD_transl_5_sym_varpar: qJD has to be [6x1] (double)');
assert(isa(r_i_i_C,'double') && isreal(r_i_i_C) && all(size(r_i_i_C) == [3 1]), ...
	'S6RRRPRR9_jacobiaD_transl_5_sym_varpar: Position vector r_i_i_C has to be [3x1] double');
assert(isreal(pkin) && all(size(pkin) == [13 1]), ...
  'S6RRRPRR9_jacobiaD_transl_5_sym_varpar: pkin has to be [13x1] (double)');

%% Symbolic Calculation
% From jacobiaD_transl_5_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 22:20:31
% EndTime: 2019-02-26 22:20:32
% DurationCPUTime: 1.30s
% Computational Cost: add. (1116->193), mult. (3439->332), div. (0->0), fcn. (3690->14), ass. (0->103)
t530 = sin(qJ(2));
t531 = sin(qJ(1));
t533 = cos(qJ(2));
t534 = cos(qJ(1));
t583 = cos(pkin(6));
t558 = t534 * t583;
t507 = t530 * t558 + t531 * t533;
t559 = t531 * t583;
t541 = t534 * t530 + t533 * t559;
t487 = t541 * qJD(1) + t507 * qJD(2);
t554 = t530 * t559;
t572 = qJD(2) * t530;
t575 = t534 * t533;
t488 = -qJD(1) * t554 - t531 * t572 + (qJD(2) * t583 + qJD(1)) * t575;
t525 = sin(pkin(7));
t524 = sin(pkin(13));
t582 = cos(pkin(13));
t587 = cos(qJ(3));
t553 = t587 * t582;
t529 = sin(qJ(3));
t570 = qJD(3) * t529;
t590 = -qJD(3) * t553 + t524 * t570;
t493 = t590 * t525;
t527 = cos(pkin(7));
t495 = t590 * t527;
t560 = t529 * t582;
t540 = t587 * t524 + t560;
t498 = t540 * t525;
t500 = t540 * t527;
t561 = t587 * qJD(3);
t505 = -qJD(3) * t560 - t524 * t561;
t506 = t531 * t530 - t533 * t558;
t526 = sin(pkin(6));
t539 = -t529 * t524 + t553;
t574 = qJD(1) * t531;
t460 = -t487 * t500 + t488 * t539 + t506 * t495 + t507 * t505 + (t493 * t534 + t498 * t574) * t526;
t563 = t526 * t574;
t483 = t487 * t525 + t527 * t563;
t528 = sin(qJ(5));
t532 = cos(qJ(5));
t597 = t460 * t528 - t483 * t532;
t596 = -t460 * t532 - t483 * t528;
t549 = t500 * t530 - t533 * t539;
t593 = t583 * t493 + (t549 * qJD(2) + t495 * t533 - t505 * t530) * t526;
t577 = t526 * t534;
t474 = t498 * t577 + t506 * t500 - t507 * t539;
t566 = t527 * t577;
t489 = -t506 * t525 + t566;
t592 = t474 * t532 + t489 * t528;
t591 = -t474 * t528 + t489 * t532;
t542 = t554 - t575;
t485 = t506 * qJD(1) + t542 * qJD(2);
t486 = t507 * qJD(1) + t541 * qJD(2);
t573 = qJD(1) * t534;
t457 = t485 * t500 - t486 * t539 + t541 * t495 - t542 * t505 + (-t493 * t531 + t498 * t573) * t526;
t588 = r_i_i_C(3) + pkin(11);
t586 = pkin(3) * t529;
t584 = pkin(10) + qJ(4);
t585 = t525 * t586 + t584 * t527 + pkin(9);
t581 = t525 * t526;
t580 = t525 * t528;
t579 = t525 * t532;
t578 = t526 * t531;
t576 = t529 * t533;
t571 = qJD(2) * t533;
t569 = qJD(5) * t528;
t568 = qJD(5) * t532;
t567 = pkin(3) * t570;
t565 = t527 * t587;
t564 = t587 * t530;
t556 = pkin(3) * t561;
t555 = t572 * t581;
t499 = t539 * t527;
t551 = t499 * t533 - t530 * t540;
t550 = t500 * t533 + t530 * t539;
t548 = qJD(1) * t587 * t581;
t547 = t532 * r_i_i_C(1) - t528 * r_i_i_C(2) + pkin(4);
t545 = qJD(5) * (-t528 * r_i_i_C(1) - t532 * r_i_i_C(2));
t538 = qJD(3) * t540;
t494 = t525 * t538;
t496 = t527 * t538;
t497 = t539 * t525;
t504 = t539 * qJD(3);
t535 = t487 * t499 + t488 * t540 - t494 * t577 - t506 * t496 - t497 * t563 + t507 * t504;
t523 = t587 * pkin(3) + pkin(2);
t511 = -t525 * qJD(4) + t527 * t556;
t510 = t527 * qJD(4) + t525 * t556;
t503 = t583 * t527 - t533 * t581;
t502 = -t584 * t525 + t527 * t586;
t491 = t525 * t541 + t527 * t578;
t484 = t549 * t526;
t481 = qJD(1) * t566 - t485 * t525;
t480 = t500 * t542 - t539 * t541;
t479 = -t507 * t500 - t506 * t539;
t478 = t583 * t498 + t550 * t526;
t476 = t498 * t578 - t500 * t541 - t539 * t542;
t471 = (-t550 * qJD(2) + t495 * t530 + t505 * t533) * t526;
t466 = -t487 * t539 - t488 * t500 + t507 * t495 - t506 * t505;
t464 = t485 * t539 + t486 * t500 - t495 * t542 - t505 * t541;
t456 = t485 * t499 + t486 * t540 + t541 * t496 + t542 * t504 + (-t494 * t531 + t497 * t573) * t526;
t454 = t457 * t532 + t481 * t528 + (-t476 * t528 + t491 * t532) * qJD(5);
t453 = -t457 * t528 + t481 * t532 + (-t476 * t532 - t491 * t528) * qJD(5);
t1 = [t596 * r_i_i_C(1) + t597 * r_i_i_C(2) - t460 * pkin(4) - t488 * t523 + t507 * t567 + t487 * t502 + t506 * t511 + t510 * t577 - t588 * t535 + (t591 * r_i_i_C(1) - t592 * r_i_i_C(2)) * qJD(5) + (-t534 * pkin(1) - t585 * t578) * qJD(1) (t464 * t532 - t486 * t580) * r_i_i_C(1) + (-t464 * t528 - t486 * t579) * r_i_i_C(2) + t464 * pkin(4) + t485 * t523 + t541 * t567 + t486 * t502 + t542 * t511 - t588 * (-t485 * t540 + t486 * t499 - t496 * t542 + t504 * t541) + ((-t480 * t528 - t542 * t579) * r_i_i_C(1) + (-t480 * t532 + t542 * t580) * r_i_i_C(2)) * qJD(5), t588 * t457 + (t497 * t578 - t499 * t541 + t540 * t542) * t545 + t547 * t456 + (t534 * t548 + t485 * t565 + t486 * t529 + (t587 * t542 + (-t525 * t578 + t527 * t541) * t529) * qJD(3)) * pkin(3), t481, t453 * r_i_i_C(1) - t454 * r_i_i_C(2), 0; t542 * t567 + t510 * t578 + t457 * pkin(4) + t454 * r_i_i_C(1) + t453 * r_i_i_C(2) + t485 * t502 - t486 * t523 - t541 * t511 - t588 * t456 + (-pkin(1) * t531 + t585 * t577) * qJD(1) (t466 * t532 + t488 * t580) * r_i_i_C(1) + (-t466 * t528 + t488 * t579) * r_i_i_C(2) + t466 * pkin(4) - t487 * t523 + t506 * t567 - t488 * t502 - t507 * t511 - t588 * (t487 * t540 - t488 * t499 + t507 * t496 + t506 * t504) + ((-t479 * t528 + t507 * t579) * r_i_i_C(1) + (-t479 * t532 - t507 * t580) * r_i_i_C(2)) * qJD(5), t588 * t460 + (-t497 * t577 - t506 * t499 - t507 * t540) * t545 - t547 * t535 + (t531 * t548 - t487 * t565 - t488 * t529 + (-t587 * t507 + (t506 * t527 + t525 * t577) * t529) * qJD(3)) * pkin(3), t483, -t597 * r_i_i_C(1) + t596 * r_i_i_C(2) + (t592 * r_i_i_C(1) + t591 * r_i_i_C(2)) * qJD(5), 0; 0 (t471 * t532 + t484 * t569) * r_i_i_C(1) + (-t471 * t528 + t484 * t568) * r_i_i_C(2) + t471 * pkin(4) + (t588 * (t551 * qJD(2) - t496 * t530 + t504 * t533) - t533 * t567 - t530 * t511 + (-t533 * t502 - t530 * t523) * qJD(2) + ((t528 * t571 + t530 * t568) * r_i_i_C(1) + (-t530 * t569 + t532 * t571) * r_i_i_C(2)) * t525) * t526, -t588 * t593 + (t583 * t497 + t551 * t526) * t545 + t547 * (-t583 * t494 + (-t496 * t533 - t504 * t530 + (-t499 * t530 - t533 * t540) * qJD(2)) * t526) + (-t583 * t525 * t570 + ((-t527 * t576 - t564) * qJD(3) + (-t527 * t564 - t576) * qJD(2)) * t526) * pkin(3), t555 (t528 * t593 + t532 * t555) * r_i_i_C(1) + (-t528 * t555 + t532 * t593) * r_i_i_C(2) + ((-t478 * t532 - t503 * t528) * r_i_i_C(1) + (t478 * t528 - t503 * t532) * r_i_i_C(2)) * qJD(5), 0;];
JaD_transl  = t1;
