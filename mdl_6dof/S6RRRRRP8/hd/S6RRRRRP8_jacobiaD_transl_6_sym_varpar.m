% Zeitableitung der translatorischen Teilmatrix der analytischen Jacobi-Matrix für Segment Nr. 6 (0=Basis) von
% S6RRRRRP8
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
%   pkin=[a2,a3,a4,a5,a6,alpha2,d1,d2,d3,d4,d5]';
%
% Output:
% JaD_transl [3x6]
%   Zeitableitung der translatorischen Teilmatrix der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-02-26 22:43
% Revision: d75aae1ac561373cd3be920984c3df29a1c2ecc8 (2019-02-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JaD_transl = S6RRRRRP8_jacobiaD_transl_6_sym_varpar(qJ, qJD, r_i_i_C, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(3,1),zeros(11,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRRRP8_jacobiaD_transl_6_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRRRRP8_jacobiaD_transl_6_sym_varpar: qJD has to be [6x1] (double)');
assert(isa(r_i_i_C,'double') && isreal(r_i_i_C) && all(size(r_i_i_C) == [3 1]), ...
	'S6RRRRRP8_jacobiaD_transl_6_sym_varpar: Position vector r_i_i_C has to be [3x1] double');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6RRRRRP8_jacobiaD_transl_6_sym_varpar: pkin has to be [11x1] (double)');

%% Symbolic Calculation
% From jacobiaD_transl_6_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 22:43:52
% EndTime: 2019-02-26 22:43:53
% DurationCPUTime: 0.98s
% Computational Cost: add. (1606->147), mult. (2965->229), div. (0->0), fcn. (3031->12), ass. (0->94)
t546 = cos(pkin(6));
t549 = sin(qJ(2));
t554 = cos(qJ(1));
t599 = t554 * t549;
t550 = sin(qJ(1));
t553 = cos(qJ(2));
t600 = t550 * t553;
t529 = t546 * t599 + t600;
t530 = t546 * t600 + t599;
t514 = t530 * qJD(1) + t529 * qJD(2);
t551 = cos(qJ(5));
t544 = qJ(3) + qJ(4);
t541 = sin(t544);
t542 = cos(t544);
t545 = sin(pkin(6));
t602 = t545 * t554;
t518 = -t529 * t542 + t541 * t602;
t598 = t554 * t553;
t601 = t550 * t549;
t528 = -t546 * t598 + t601;
t547 = sin(qJ(5));
t626 = t518 * t551 - t528 * t547;
t635 = -t626 * qJD(5) - t514 * t551;
t627 = t518 * t547 + t528 * t551;
t634 = t627 * qJD(5) + t514 * t547;
t620 = r_i_i_C(1) + pkin(5);
t633 = t620 * t551 + pkin(4);
t619 = r_i_i_C(3) + qJ(6);
t621 = pkin(11) + r_i_i_C(2);
t593 = qJD(6) * t547;
t632 = -qJD(5) * t547 * t620 + t593;
t543 = qJD(3) + qJD(4);
t548 = sin(qJ(3));
t607 = t541 * t543;
t625 = -qJD(3) * t548 * pkin(3) - pkin(4) * t607 + (t621 * t543 + t593) * t542;
t552 = cos(qJ(3));
t540 = t552 * pkin(3) + pkin(2);
t624 = t542 * pkin(4) + t621 * t541 + t540;
t622 = (qJD(2) * t542 - qJD(5)) * t549 + t553 * t607;
t588 = t546 * t601;
t566 = t588 - t598;
t608 = t566 * t542;
t606 = t545 * t549;
t605 = t545 * t550;
t604 = t545 * t552;
t603 = t545 * t553;
t597 = qJD(1) * t545;
t596 = qJD(2) * t549;
t595 = qJD(5) * t542;
t594 = qJD(5) * t551;
t592 = t551 * qJD(6);
t591 = t541 * t606;
t590 = t542 * t606;
t589 = t547 * t603;
t586 = t542 * t602;
t585 = t550 * t597;
t584 = t554 * t597;
t583 = t545 * t596;
t582 = qJD(2) * t603;
t515 = -qJD(1) * t588 - t550 * t596 + (qJD(2) * t546 + qJD(1)) * t598;
t577 = -t515 * t542 + t543 * t586;
t513 = t529 * qJD(1) + t530 * qJD(2);
t574 = t543 * t605 - t513;
t573 = t530 * t595 - t513;
t572 = t528 * t595 + t515;
t520 = t541 * t605 - t608;
t569 = t520 * t551 + t530 * t547;
t568 = -t520 * t547 + t530 * t551;
t567 = (qJD(2) - t595) * t553;
t565 = -t529 * t543 + t585;
t564 = t543 * t546 + t582;
t512 = t528 * qJD(1) + t566 * qJD(2);
t561 = -qJD(5) * t566 + t512 * t542 + t530 * t607;
t560 = qJD(5) * t529 - t514 * t542 + t528 * t607;
t488 = t565 * t542 + (t543 * t602 - t515) * t541;
t486 = t574 * t541 - t542 * t584 - t543 * t608;
t487 = t541 * t584 + t574 * t542 + t566 * t607;
t519 = t541 * t566 + t542 * t605;
t558 = t621 * t487 + t632 * t519 + t619 * (-t486 * t547 + t519 * t594) - t633 * t486;
t489 = -t529 * t607 + t541 * t585 - t577;
t516 = -t529 * t541 - t586;
t557 = t621 * t489 + t632 * t516 + t619 * (t488 * t547 + t516 * t594) + t633 * t488;
t504 = -t564 * t541 - t543 * t590;
t505 = t564 * t542 - t543 * t591;
t526 = t546 * t542 - t591;
t556 = t621 * t505 + t632 * t526 + t619 * (t504 * t547 + t526 * t594) + t633 * t504;
t555 = -pkin(10) - pkin(9);
t527 = t546 * t541 + t590;
t492 = -qJD(5) * t589 + t505 * t547 + t527 * t594 - t551 * t583;
t491 = -t565 * t541 + t577;
t466 = t489 * t547 + t635;
t465 = t568 * qJD(5) + t487 * t551 - t512 * t547;
t464 = t569 * qJD(5) + t487 * t547 + t512 * t551;
t1 = [t627 * qJD(6) + t491 * pkin(4) - t515 * t540 + t514 * t555 + t621 * t488 + t620 * (t491 * t551 - t634) + t619 * (t491 * t547 - t635) + (-t554 * pkin(1) - pkin(8) * t605) * qJD(1) + (-t548 * t585 + (t529 * t548 + t552 * t602) * qJD(3)) * pkin(3), t566 * t592 + t513 * t555 + t620 * (t573 * t547 + t561 * t551) + t619 * (t561 * t547 - t573 * t551) - t625 * t530 + t624 * t512 (t552 * t584 + t513 * t548 + (-t548 * t605 + t552 * t566) * qJD(3)) * pkin(3) + t558, t558, t569 * qJD(6) - t620 * t464 + t619 * t465, t464; -t568 * qJD(6) + t487 * pkin(4) - t513 * t540 + t512 * t555 + t621 * t486 + t620 * t465 + t619 * t464 + (-t550 * pkin(1) + pkin(8) * t602) * qJD(1) + (t548 * t584 + (t548 * t566 + t550 * t604) * qJD(3)) * pkin(3), -t529 * t592 - t515 * t555 + t620 * (t572 * t547 + t560 * t551) + t619 * (t560 * t547 - t572 * t551) - t625 * t528 - t624 * t514 (t552 * t585 - t515 * t548 + (-t529 * t552 + t548 * t602) * qJD(3)) * pkin(3) + t557, t557, -t626 * qJD(6) + t619 * (t489 * t551 + t634) - t620 * t466, t466; 0 (t620 * (t547 * t567 - t622 * t551) - t619 * (t622 * t547 + t551 * t567) - t549 * t592 + t625 * t553 + (-t549 * t624 - t553 * t555) * qJD(2)) * t545 (-t548 * t582 + (-t546 * t548 - t549 * t604) * qJD(3)) * pkin(3) + t556, t556 -(-t527 * t551 + t589) * qJD(6) + t619 * (t547 * t583 + t505 * t551 + (-t527 * t547 - t551 * t603) * qJD(5)) - t620 * t492, t492;];
JaD_transl  = t1;
