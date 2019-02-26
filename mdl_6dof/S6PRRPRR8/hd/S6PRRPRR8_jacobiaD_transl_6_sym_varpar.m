% Zeitableitung der translatorischen Teilmatrix der analytischen Jacobi-Matrix für Segment Nr. 6 (0=Basis) von
% S6PRRPRR8
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
%   pkin=[a2,a3,a4,a5,a6,alpha2,alpha3,d2,d3,d5,d6,theta1]';
%
% Output:
% JaD_transl [3x6]
%   Zeitableitung der translatorischen Teilmatrix der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-02-26 20:08
% Revision: d75aae1ac561373cd3be920984c3df29a1c2ecc8 (2019-02-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JaD_transl = S6PRRPRR8_jacobiaD_transl_6_sym_varpar(qJ, qJD, r_i_i_C, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(3,1),zeros(12,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRRPRR8_jacobiaD_transl_6_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6PRRPRR8_jacobiaD_transl_6_sym_varpar: qJD has to be [6x1] (double)');
assert(isa(r_i_i_C,'double') && isreal(r_i_i_C) && all(size(r_i_i_C) == [3 1]), ...
	'S6PRRPRR8_jacobiaD_transl_6_sym_varpar: Position vector r_i_i_C has to be [3x1] double');
assert(isreal(pkin) && all(size(pkin) == [12 1]), ...
  'S6PRRPRR8_jacobiaD_transl_6_sym_varpar: pkin has to be [12x1] (double)');

%% Symbolic Calculation
% From jacobiaD_transl_6_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 20:08:11
% EndTime: 2019-02-26 20:08:11
% DurationCPUTime: 0.89s
% Computational Cost: add. (1242->168), mult. (4014->298), div. (0->0), fcn. (4428->14), ass. (0->105)
t569 = sin(pkin(12));
t572 = cos(pkin(12));
t578 = sin(qJ(2));
t574 = cos(pkin(6));
t582 = cos(qJ(2));
t618 = t574 * t582;
t557 = -t569 * t578 + t572 * t618;
t635 = pkin(4) + pkin(9);
t634 = r_i_i_C(3) + pkin(11);
t619 = t574 * t578;
t558 = t569 * t582 + t572 * t619;
t577 = sin(qJ(3));
t633 = t558 * t577;
t581 = cos(qJ(3));
t632 = t558 * t581;
t590 = t569 * t619 - t572 * t582;
t631 = t590 * t577;
t570 = sin(pkin(7));
t571 = sin(pkin(6));
t629 = t570 * t571;
t576 = sin(qJ(5));
t628 = t570 * t576;
t627 = t570 * t577;
t580 = cos(qJ(5));
t626 = t570 * t580;
t625 = t570 * t581;
t624 = t570 * t582;
t623 = t571 * t572;
t573 = cos(pkin(7));
t622 = t571 * t573;
t621 = t573 * t577;
t620 = t573 * t581;
t617 = t577 * t578;
t616 = t577 * t582;
t615 = t578 * t581;
t614 = t581 * t582;
t613 = qJD(2) * t571;
t612 = qJD(5) * t576;
t611 = qJD(5) * t580;
t610 = t635 * t570;
t609 = t578 * t629;
t608 = t571 * t625;
t607 = t571 * t617;
t605 = t573 * t614;
t604 = t574 * t625;
t603 = t570 * t613;
t602 = qJD(3) * t627;
t601 = t571 * t605;
t600 = t578 * t603;
t599 = t582 * t603;
t575 = sin(qJ(6));
t579 = cos(qJ(6));
t598 = -t579 * r_i_i_C(1) + t575 * r_i_i_C(2);
t597 = -t575 * r_i_i_C(1) - t579 * r_i_i_C(2);
t529 = -t557 * t620 + t572 * t608 + t633;
t544 = -t557 * t570 - t572 * t622;
t520 = t529 * t576 + t544 * t580;
t591 = t569 * t618 + t572 * t578;
t531 = -t569 * t608 + t591 * t620 - t631;
t545 = t569 * t622 + t570 * t591;
t522 = t531 * t576 + t545 * t580;
t542 = -t601 - t604 + t607;
t556 = -t571 * t624 + t574 * t573;
t596 = t542 * t580 - t556 * t576;
t534 = t542 * t576 + t556 * t580;
t595 = pkin(5) - t598;
t594 = pkin(3) + pkin(10) - t597;
t537 = t557 * t577 + t558 * t620;
t523 = t537 * t576 + t558 * t626;
t539 = -t577 * t591 - t590 * t620;
t524 = t539 * t576 - t590 * t626;
t538 = t557 * t581 - t558 * t621;
t593 = t557 * t573 - t570 * t623;
t540 = -t581 * t591 + t590 * t621;
t592 = t569 * t629 - t573 * t591;
t589 = t573 * t615 + t616;
t588 = t573 * t616 + t615;
t587 = -t573 * t617 + t614;
t586 = qJD(6) * t598;
t585 = qJD(6) * t597;
t550 = t589 * t571;
t541 = t550 * t576 + t580 * t609;
t532 = t592 * t577 - t581 * t590;
t584 = t595 * t576 - t634 * t580 + qJ(4);
t583 = qJD(4) + t576 * t585 + (t634 * t576 + t595 * t580) * qJD(5);
t555 = t590 * qJD(2);
t554 = t591 * qJD(2);
t553 = t558 * qJD(2);
t552 = t557 * qJD(2);
t551 = t587 * t571;
t543 = t588 * t571 + t574 * t627;
t535 = -qJD(2) * t601 - t571 * qJD(3) * t614 + (qJD(3) * t573 + qJD(2)) * t607;
t530 = t593 * t577 + t632;
t528 = qJD(3) * t604 + ((t605 - t617) * qJD(3) + t587 * qJD(2)) * t571;
t527 = t574 * t602 + (t589 * qJD(2) + t588 * qJD(3)) * t571;
t517 = t540 * qJD(3) - t554 * t620 + t555 * t577;
t515 = t538 * qJD(3) + t552 * t620 - t553 * t577;
t512 = t555 * t621 - t554 * t581 + (t592 * t581 + t631) * qJD(3);
t511 = t532 * qJD(3) - t554 * t577 - t555 * t620;
t510 = -t553 * t621 + t552 * t581 + (t593 * t581 - t633) * qJD(3);
t509 = t552 * t577 + t553 * t620 - t602 * t623 + (t557 * t621 + t632) * qJD(3);
t506 = t596 * qJD(5) + t527 * t576 + t580 * t600;
t500 = -t511 * t576 - t531 * t611 + t545 * t612 + t555 * t626;
t498 = -t509 * t576 - t529 * t611 + t544 * t612 - t553 * t626;
t1 = [0, t555 * pkin(2) + t517 * qJ(4) + t539 * qJD(4) - t554 * t610 + t595 * (-t554 * t626 + t517 * t576 + (t539 * t580 + t590 * t628) * qJD(5)) + t594 * (-t539 * qJD(3) + t554 * t621 + t555 * t581) - t634 * (-t524 * qJD(5) + t517 * t580 + t554 * t628) + ((-t524 * t575 + t540 * t579) * r_i_i_C(1) + (-t524 * t579 - t540 * t575) * r_i_i_C(2)) * qJD(6), -t594 * t511 + t584 * t512 + t531 * t586 + t583 * t532, t511, -t634 * t500 + (t531 * t580 - t545 * t576) * t585 + t595 * (-t522 * qJD(5) + t511 * t580 + t555 * t628) (t500 * t575 + t512 * t579) * r_i_i_C(1) + (t500 * t579 - t512 * t575) * r_i_i_C(2) + ((-t522 * t579 - t532 * t575) * r_i_i_C(1) + (t522 * t575 - t532 * t579) * r_i_i_C(2)) * qJD(6); 0, -t553 * pkin(2) + t515 * qJ(4) + t537 * qJD(4) + t552 * t610 + t595 * (t552 * t626 + t515 * t576 + (t537 * t580 - t558 * t628) * qJD(5)) + t594 * (-t537 * qJD(3) - t552 * t621 - t553 * t581) - t634 * (-t523 * qJD(5) + t515 * t580 - t552 * t628) + ((-t523 * t575 + t538 * t579) * r_i_i_C(1) + (-t523 * t579 - t538 * t575) * r_i_i_C(2)) * qJD(6), -t594 * t509 + t584 * t510 + t529 * t586 + t583 * t530, t509, -t634 * t498 + (t529 * t580 - t544 * t576) * t585 + t595 * (-t520 * qJD(5) + t509 * t580 - t553 * t628) (t498 * t575 + t510 * t579) * r_i_i_C(1) + (t498 * t579 - t510 * t575) * r_i_i_C(2) + ((-t520 * t579 - t530 * t575) * r_i_i_C(1) + (t520 * t575 - t530 * t579) * r_i_i_C(2)) * qJD(6); 0, -t535 * qJ(4) + t550 * qJD(4) + t595 * (t580 * t599 - t535 * t576 + (t550 * t580 - t576 * t609) * qJD(5)) - t594 * (t588 * qJD(2) + t589 * qJD(3)) * t571 + t634 * (t541 * qJD(5) + t535 * t580 + t576 * t599) + ((-t541 * t575 + t551 * t579) * r_i_i_C(1) + (-t541 * t579 - t551 * t575) * r_i_i_C(2)) * qJD(6) + (-pkin(2) * t578 + t635 * t624) * t613, -t594 * t527 + t584 * t528 + t542 * t586 + t583 * t543, t527, t634 * t506 + t596 * t585 + t595 * (-t534 * qJD(5) + t527 * t580 - t576 * t600) (-t506 * t575 + t528 * t579) * r_i_i_C(1) + (-t506 * t579 - t528 * t575) * r_i_i_C(2) + ((-t534 * t579 - t543 * t575) * r_i_i_C(1) + (t534 * t575 - t543 * t579) * r_i_i_C(2)) * qJD(6);];
JaD_transl  = t1;
