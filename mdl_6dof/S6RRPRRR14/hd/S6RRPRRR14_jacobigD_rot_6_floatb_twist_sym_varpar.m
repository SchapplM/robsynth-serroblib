% Zeitableitung der rotatorischen Teilmatrix der geometrischen Jacobi-Matrix für Segment Nr. 6 (0=Basis) von
% S6RRPRRR14
% Use Code from Maple symbolic Code Generation
%
% geometrische Jacobi-Matrix: Differentieller Zusammenhang zwischen
% Endeffektorposition und verallgemeinerten Koordinaten.
%
%
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [14x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,alpha3,alpha4,d1,d2,d4,d5,d6,theta3]';
%
% Output:
% JgD_rot [3x6]
%   Zeitableitung der rotatorischen Teilmatrix der geometrischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox (ehem. IRT-Maple-Toolbox)
% Datum: 2018-12-10 18:38
% Revision: bb42a8b95257d9bc83910d26e849f5825122f662 (2018-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für mechatronische Systeme, Universität Hannover

function JgD_rot = S6RRPRRR14_jacobigD_rot_6_floatb_twist_sym_varpar(qJ, qJD, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(14,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPRRR14_jacobigD_rot_6_floatb_twist_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRPRRR14_jacobigD_rot_6_floatb_twist_sym_varpar: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [14 1]), ...
  'S6RRPRRR14_jacobigD_rot_6_floatb_twist_sym_varpar: pkin has to be [14x1] (double)');

%% Symbolic Calculation
% From jacobigD_rot_6_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2018-12-10 18:38:25
% EndTime: 2018-12-10 18:38:25
% DurationCPUTime: 0.34s
% Computational Cost: add. (986->116), mult. (1196->182), div. (0->0), fcn. (999->28), ass. (0->107)
t598 = qJD(2) / 0.2e1;
t597 = qJD(4) / 0.2e1;
t565 = pkin(6) + qJ(2);
t559 = cos(t565);
t544 = t559 * t598;
t566 = pkin(6) - qJ(2);
t560 = cos(t566);
t591 = qJD(2) * t560;
t528 = t544 - t591 / 0.2e1;
t569 = sin(pkin(7));
t596 = t528 * t569;
t570 = sin(pkin(6));
t578 = sin(qJ(1));
t595 = t570 * t578;
t582 = cos(qJ(1));
t594 = t570 * t582;
t593 = qJD(1) * t570;
t555 = sin(t565);
t592 = qJD(2) * t555;
t590 = qJD(2) * t578;
t589 = qJD(2) * t582;
t563 = pkin(8) + qJ(4);
t553 = sin(t563);
t588 = qJD(4) * t553;
t564 = pkin(8) - qJ(4);
t558 = cos(t564);
t587 = qJD(4) * t558;
t576 = sin(qJ(4));
t586 = qJD(4) * t576;
t580 = cos(qJ(4));
t585 = qJD(4) * t580;
t584 = t578 * t593;
t583 = t582 * t593;
t548 = t555 / 0.2e1;
t556 = sin(t566);
t536 = t548 - t556 / 0.2e1;
t581 = cos(qJ(2));
t518 = t582 * t536 + t578 * t581;
t520 = -t578 * t536 + t582 * t581;
t550 = t560 / 0.2e1;
t539 = t550 + t559 / 0.2e1;
t577 = sin(qJ(2));
t517 = t582 * t539 - t578 * t577;
t519 = -t578 * t539 - t582 * t577;
t579 = cos(qJ(5));
t575 = sin(qJ(5));
t574 = cos(pkin(6));
t573 = cos(pkin(7));
t572 = cos(pkin(8));
t571 = cos(pkin(14));
t568 = sin(pkin(8));
t567 = sin(pkin(14));
t562 = pkin(7) - pkin(14);
t561 = pkin(7) + pkin(14);
t557 = cos(t563);
t554 = sin(t564);
t552 = cos(t561);
t551 = sin(t562);
t549 = t558 / 0.2e1;
t547 = t553 / 0.2e1;
t546 = cos(t562) / 0.2e1;
t545 = sin(t561) / 0.2e1;
t543 = t556 * t598;
t542 = t557 * t597;
t541 = t554 * t597;
t540 = t550 - t559 / 0.2e1;
t538 = t549 - t557 / 0.2e1;
t537 = t549 + t557 / 0.2e1;
t535 = t548 + t556 / 0.2e1;
t534 = t547 - t554 / 0.2e1;
t533 = t547 + t554 / 0.2e1;
t532 = t546 - t552 / 0.2e1;
t531 = t546 + t552 / 0.2e1;
t530 = t545 - t551 / 0.2e1;
t529 = t545 + t551 / 0.2e1;
t527 = t544 + t591 / 0.2e1;
t526 = t543 - t592 / 0.2e1;
t525 = t543 + t592 / 0.2e1;
t524 = t542 - t587 / 0.2e1;
t523 = t542 + t587 / 0.2e1;
t522 = t541 - t588 / 0.2e1;
t521 = t541 + t588 / 0.2e1;
t516 = -t535 * t569 + t574 * t573;
t515 = -t519 * t569 + t573 * t595;
t514 = -t517 * t569 - t573 * t594;
t513 = t525 * t571 + t528 * t530;
t512 = -t525 * t567 + t528 * t531;
t511 = qJD(1) * t520 + t582 * t527 - t577 * t590;
t510 = qJD(1) * t519 + t582 * t526 - t581 * t590;
t509 = -qJD(1) * t518 - t578 * t527 - t577 * t589;
t508 = -qJD(1) * t517 - t578 * t526 - t581 * t589;
t507 = t535 * t530 + t574 * t532 + t540 * t571;
t506 = t574 * t529 + t535 * t531 - t540 * t567;
t505 = -t510 * t569 + t573 * t584;
t504 = -t508 * t569 + t573 * t583;
t503 = -t512 * t568 - t572 * t596;
t502 = t519 * t530 + t520 * t571 + t532 * t595;
t501 = t519 * t531 - t520 * t567 + t529 * t595;
t500 = t517 * t530 + t518 * t571 - t532 * t594;
t499 = t517 * t531 - t518 * t567 - t529 * t594;
t498 = t510 * t530 + t511 * t571 + t532 * t584;
t497 = t510 * t531 - t511 * t567 + t529 * t584;
t496 = t508 * t530 + t509 * t571 + t532 * t583;
t495 = t508 * t531 - t509 * t567 + t529 * t583;
t494 = -t497 * t568 + t505 * t572;
t493 = -t495 * t568 + t504 * t572;
t1 = [0, t583, 0, t493, -t495 * t537 + t496 * t576 - t501 * t522 + t502 * t585 - t504 * t533 - t515 * t524 (t495 * t534 + t496 * t580 + t501 * t523 - t502 * t586 + t504 * t538 + t515 * t521) * t575 - t493 * t579 + ((t501 * t534 + t502 * t580 + t515 * t538) * t579 + (-t501 * t568 + t515 * t572) * t575) * qJD(5); 0, t584, 0, t494, -t497 * t537 + t498 * t576 - t499 * t522 + t500 * t585 - t505 * t533 - t514 * t524 (t497 * t534 + t498 * t580 + t499 * t523 - t500 * t586 + t505 * t538 + t514 * t521) * t575 - t494 * t579 + ((t499 * t534 + t500 * t580 + t514 * t538) * t579 + (-t499 * t568 + t514 * t572) * t575) * qJD(5); 0, 0, 0, t503, -t506 * t522 + t507 * t585 - t512 * t537 + t513 * t576 - t516 * t524 + t533 * t596 (t506 * t523 - t507 * t586 + t512 * t534 + t513 * t580 + t516 * t521 - t538 * t596) * t575 - t503 * t579 + ((t506 * t534 + t507 * t580 + t516 * t538) * t579 + (-t506 * t568 + t516 * t572) * t575) * qJD(5);];
JgD_rot  = t1;
