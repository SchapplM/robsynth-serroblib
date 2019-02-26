% Zeitableitung der translatorischen Teilmatrix der analytischen Jacobi-Matrix für Segment Nr. 5 (0=Basis) von
% S6PRRRRR5
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
%   pkin=[a2,a3,a4,a5,a6,alpha2,alpha3,d2,d3,d4,d5,d6,theta1]';
%
% Output:
% JaD_transl [3x6]
%   Zeitableitung der translatorischen Teilmatrix der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-02-26 20:21
% Revision: d75aae1ac561373cd3be920984c3df29a1c2ecc8 (2019-02-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JaD_transl = S6PRRRRR5_jacobiaD_transl_5_sym_varpar(qJ, qJD, r_i_i_C, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(3,1),zeros(13,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRRRRR5_jacobiaD_transl_5_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6PRRRRR5_jacobiaD_transl_5_sym_varpar: qJD has to be [6x1] (double)');
assert(isa(r_i_i_C,'double') && isreal(r_i_i_C) && all(size(r_i_i_C) == [3 1]), ...
	'S6PRRRRR5_jacobiaD_transl_5_sym_varpar: Position vector r_i_i_C has to be [3x1] double');
assert(isreal(pkin) && all(size(pkin) == [13 1]), ...
  'S6PRRRRR5_jacobiaD_transl_5_sym_varpar: pkin has to be [13x1] (double)');

%% Symbolic Calculation
% From jacobiaD_transl_5_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 20:21:14
% EndTime: 2019-02-26 20:21:15
% DurationCPUTime: 0.93s
% Computational Cost: add. (1151->179), mult. (3722->319), div. (0->0), fcn. (4116->14), ass. (0->103)
t551 = sin(qJ(2));
t554 = cos(qJ(2));
t546 = cos(pkin(13));
t605 = cos(pkin(6));
t583 = t546 * t605;
t604 = sin(pkin(13));
t532 = -t604 * t551 + t554 * t583;
t548 = sin(qJ(5));
t552 = cos(qJ(5));
t567 = (r_i_i_C(1) * t548 + r_i_i_C(2) * t552) * qJD(5);
t608 = r_i_i_C(3) + pkin(11);
t607 = cos(qJ(3));
t544 = sin(pkin(7));
t606 = pkin(9) * t544;
t549 = sin(qJ(4));
t603 = t544 * t549;
t553 = cos(qJ(4));
t602 = t544 * t553;
t601 = t544 * t554;
t545 = sin(pkin(6));
t600 = t545 * t546;
t547 = cos(pkin(7));
t550 = sin(qJ(3));
t599 = t547 * t550;
t598 = t550 * t551;
t597 = t550 * t554;
t596 = qJD(2) * t545;
t595 = qJD(5) * t548;
t594 = qJD(5) * t552;
t593 = t544 * t600;
t592 = t544 * t545 * t551;
t591 = t545 * t598;
t561 = -t551 * t583 - t604 * t554;
t590 = t561 * t607;
t589 = t547 * t607;
t588 = t607 * t551;
t587 = t607 * t554;
t586 = t544 * t596;
t585 = t544 * t605;
t584 = t545 * t604;
t581 = t547 * t587;
t580 = t551 * t586;
t579 = t554 * t586;
t577 = t550 * t585;
t576 = t544 * t584;
t574 = t605 * t604;
t573 = t545 * t581;
t505 = -t590 + (t532 * t547 - t593) * t550;
t519 = -t532 * t544 - t547 * t600;
t497 = t505 * t553 + t519 * t549;
t572 = -t505 * t549 + t519 * t553;
t559 = t546 * t551 + t554 * t574;
t560 = -t546 * t554 + t551 * t574;
t507 = -t560 * t607 + (-t547 * t559 + t576) * t550;
t520 = t544 * t559 + t547 * t584;
t499 = t507 * t553 + t520 * t549;
t571 = -t507 * t549 + t520 * t553;
t563 = t547 * t597 + t588;
t518 = t563 * t545 + t577;
t531 = -t545 * t601 + t605 * t547;
t509 = t518 * t553 + t531 * t549;
t570 = -t518 * t549 + t531 * t553;
t569 = r_i_i_C(1) * t552 - r_i_i_C(2) * t548 + pkin(4);
t568 = t607 * t585;
t513 = t532 * t607 + t561 * t599;
t500 = t513 * t553 - t561 * t603;
t515 = -t559 * t607 + t560 * t599;
t501 = t515 * t553 - t560 * t603;
t564 = -t547 * t598 + t587;
t526 = t564 * t545;
t516 = t526 * t553 + t549 * t592;
t566 = -t532 * t550 + t561 * t589;
t565 = t550 * t559 + t560 * t589;
t562 = t547 * t588 + t597;
t558 = -t608 * t549 - t569 * t553 - pkin(3);
t557 = t532 * t589 + t550 * t561 - t607 * t593;
t556 = t550 * t560 - t559 * t589 + t607 * t576;
t555 = t553 * t567 + (t569 * t549 - t608 * t553) * qJD(4);
t530 = t560 * qJD(2);
t529 = t559 * qJD(2);
t528 = t561 * qJD(2);
t527 = t532 * qJD(2);
t525 = t562 * t545;
t517 = -t568 - t573 + t591;
t511 = (-t563 * qJD(2) - t562 * qJD(3)) * t545;
t510 = -qJD(2) * t573 - t545 * qJD(3) * t587 + (qJD(3) * t547 + qJD(2)) * t591;
t503 = qJD(3) * t568 + ((t581 - t598) * qJD(3) + t564 * qJD(2)) * t545;
t502 = qJD(3) * t577 + (t562 * qJD(2) + t563 * qJD(3)) * t545;
t495 = t565 * qJD(3) + t529 * t599 + t530 * t607;
t494 = t515 * qJD(3) - t529 * t589 + t530 * t550;
t493 = t566 * qJD(3) - t527 * t599 + t528 * t607;
t492 = t513 * qJD(3) + t527 * t589 + t528 * t550;
t491 = t549 * t579 + t511 * t553 + (-t526 * t549 + t553 * t592) * qJD(4);
t489 = t556 * qJD(3) - t529 * t607 + t530 * t599;
t488 = t507 * qJD(3) - t529 * t550 - t530 * t589;
t487 = t557 * qJD(3) + t527 * t607 + t528 * t599;
t486 = t527 * t550 - t528 * t589 + (t532 * t599 - t550 * t593 - t590) * qJD(3);
t485 = t570 * qJD(4) + t503 * t553 + t549 * t580;
t483 = -t529 * t603 + t495 * t553 + (-t515 * t549 - t560 * t602) * qJD(4);
t481 = t527 * t603 + t493 * t553 + (-t513 * t549 - t561 * t602) * qJD(4);
t479 = t571 * qJD(4) + t489 * t553 - t530 * t603;
t477 = t572 * qJD(4) + t487 * t553 - t528 * t603;
t1 = [0 (t483 * t552 + t494 * t548) * r_i_i_C(1) + (-t483 * t548 + t494 * t552) * r_i_i_C(2) + t483 * pkin(4) + t495 * pkin(3) + t494 * pkin(10) + t530 * pkin(2) - t529 * t606 + t608 * (t501 * qJD(4) + t495 * t549 + t529 * t602) + ((-t501 * t548 - t552 * t565) * r_i_i_C(1) + (-t501 * t552 + t548 * t565) * r_i_i_C(2)) * qJD(5) (t489 * t548 + t507 * t594) * r_i_i_C(1) + (t489 * t552 - t507 * t595) * r_i_i_C(2) + t489 * pkin(10) + t558 * t488 - t555 * t556, t608 * t479 - t571 * t567 + t569 * (-t499 * qJD(4) - t489 * t549 - t530 * t602) (-t479 * t548 + t488 * t552) * r_i_i_C(1) + (-t479 * t552 - t488 * t548) * r_i_i_C(2) + ((-t499 * t552 + t548 * t556) * r_i_i_C(1) + (t499 * t548 + t552 * t556) * r_i_i_C(2)) * qJD(5), 0; 0 (t481 * t552 + t492 * t548) * r_i_i_C(1) + (-t481 * t548 + t492 * t552) * r_i_i_C(2) + t481 * pkin(4) + t493 * pkin(3) + t492 * pkin(10) + t528 * pkin(2) + t527 * t606 + t608 * (t500 * qJD(4) + t493 * t549 - t527 * t602) + ((-t500 * t548 - t552 * t566) * r_i_i_C(1) + (-t500 * t552 + t548 * t566) * r_i_i_C(2)) * qJD(5) (t487 * t548 + t505 * t594) * r_i_i_C(1) + (t487 * t552 - t505 * t595) * r_i_i_C(2) + t487 * pkin(10) + t558 * t486 - t555 * t557, t608 * t477 - t572 * t567 + t569 * (-t497 * qJD(4) - t487 * t549 - t528 * t602) (-t477 * t548 + t486 * t552) * r_i_i_C(1) + (-t477 * t552 - t486 * t548) * r_i_i_C(2) + ((-t497 * t552 + t548 * t557) * r_i_i_C(1) + (t497 * t548 + t552 * t557) * r_i_i_C(2)) * qJD(5), 0; 0 (t491 * t552 - t510 * t548) * r_i_i_C(1) + (-t491 * t548 - t510 * t552) * r_i_i_C(2) + t491 * pkin(4) + t511 * pkin(3) - t510 * pkin(10) + t608 * (t516 * qJD(4) + t511 * t549 - t553 * t579) + ((-t516 * t548 + t525 * t552) * r_i_i_C(1) + (-t516 * t552 - t525 * t548) * r_i_i_C(2)) * qJD(5) + (-pkin(2) * t551 + pkin(9) * t601) * t596 (t503 * t548 + t518 * t594) * r_i_i_C(1) + (t503 * t552 - t518 * t595) * r_i_i_C(2) + t503 * pkin(10) + t558 * t502 + t555 * t517, t608 * t485 - t570 * t567 + t569 * (-t509 * qJD(4) - t503 * t549 + t553 * t580) (-t485 * t548 + t502 * t552) * r_i_i_C(1) + (-t485 * t552 - t502 * t548) * r_i_i_C(2) + ((-t509 * t552 - t517 * t548) * r_i_i_C(1) + (t509 * t548 - t517 * t552) * r_i_i_C(2)) * qJD(5), 0;];
JaD_transl  = t1;
