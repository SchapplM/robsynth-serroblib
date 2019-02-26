% Zeitableitung der translatorischen Teilmatrix der analytischen Jacobi-Matrix für Segment Nr. 6 (0=Basis) von
% S6RPRPRR13
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
%   pkin=[a2,a3,a4,a5,a6,alpha2,alpha3,d1,d3,d5,d6,theta2]';
%
% Output:
% JaD_transl [3x6]
%   Zeitableitung der translatorischen Teilmatrix der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-02-26 20:55
% Revision: d75aae1ac561373cd3be920984c3df29a1c2ecc8 (2019-02-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JaD_transl = S6RPRPRR13_jacobiaD_transl_6_sym_varpar(qJ, qJD, r_i_i_C, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(3,1),zeros(12,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRPRR13_jacobiaD_transl_6_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RPRPRR13_jacobiaD_transl_6_sym_varpar: qJD has to be [6x1] (double)');
assert(isa(r_i_i_C,'double') && isreal(r_i_i_C) && all(size(r_i_i_C) == [3 1]), ...
	'S6RPRPRR13_jacobiaD_transl_6_sym_varpar: Position vector r_i_i_C has to be [3x1] double');
assert(isreal(pkin) && all(size(pkin) == [12 1]), ...
  'S6RPRPRR13_jacobiaD_transl_6_sym_varpar: pkin has to be [12x1] (double)');

%% Symbolic Calculation
% From jacobiaD_transl_6_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 20:55:47
% EndTime: 2019-02-26 20:55:48
% DurationCPUTime: 0.91s
% Computational Cost: add. (1269->129), mult. (4073->213), div. (0->0), fcn. (4514->14), ass. (0->91)
t580 = sin(qJ(3));
t584 = cos(qJ(3));
t573 = sin(pkin(12));
t576 = cos(pkin(12));
t581 = sin(qJ(1));
t585 = cos(qJ(1));
t635 = cos(pkin(6));
t612 = t585 * t635;
t591 = -t581 * t573 + t576 * t612;
t592 = -t573 * t612 - t581 * t576;
t574 = sin(pkin(7));
t575 = sin(pkin(6));
t626 = t575 * t585;
t617 = t574 * t626;
t577 = cos(pkin(7));
t625 = t577 * t580;
t534 = t580 * t617 + t584 * t592 - t591 * t625;
t613 = t581 * t635;
t593 = t585 * t573 + t576 * t613;
t555 = t593 * qJD(1);
t606 = t573 * t613;
t622 = qJD(1) * t585;
t556 = -qJD(1) * t606 + t576 * t622;
t623 = qJD(1) * t581;
t615 = t575 * t623;
t608 = t574 * t615;
t624 = t577 * t584;
t520 = t534 * qJD(3) - t555 * t624 - t556 * t580 + t584 * t608;
t633 = t555 * t574;
t543 = t577 * t615 + t633;
t579 = sin(qJ(5));
t583 = cos(qJ(5));
t653 = t520 * t579 - t543 * t583;
t652 = -t520 * t583 - t543 * t579;
t630 = t592 * t580;
t638 = (-t577 * t591 + t617) * t584 - t630;
t651 = t638 * qJD(3) - (-t555 * t577 + t608) * t580 - t556 * t584;
t578 = sin(qJ(6));
t650 = t534 * t578;
t582 = cos(qJ(6));
t649 = t534 * t582;
t546 = t591 * t574 + t577 * t626;
t646 = t546 * t579;
t645 = t546 * t583;
t644 = pkin(10) + pkin(3);
t637 = pkin(11) + r_i_i_C(3);
t643 = pkin(9) * t577 + qJ(2);
t589 = t591 * qJD(1);
t614 = t575 * t622;
t642 = t574 * t614 - t577 * t589;
t611 = t635 * t574;
t628 = t575 * t576;
t544 = t575 * t573 * t580 - t584 * t611 - t624 * t628;
t557 = -t574 * t628 + t635 * t577;
t641 = -t544 * t583 + t557 * t579;
t561 = t585 * t576 - t606;
t627 = t575 * t581;
t619 = t574 * t627;
t640 = -t561 * t580 + t584 * t619 - t593 * t624;
t604 = -t582 * r_i_i_C(1) + t578 * r_i_i_C(2);
t600 = pkin(5) - t604;
t587 = t600 * t579 - t637 * t583 + qJ(4);
t621 = t575 * qJD(2);
t603 = -t578 * r_i_i_C(1) - t582 * r_i_i_C(2);
t531 = t584 * t617 - t591 * t624 - t630;
t602 = t531 * t583 + t646;
t526 = t531 * t579 - t645;
t524 = -t579 * t638 + t645;
t548 = t574 * t593 + t577 * t627;
t601 = -t548 * t579 - t583 * t640;
t528 = t548 * t583 - t579 * t640;
t530 = t544 * t579 + t557 * t583;
t599 = -t603 + t644;
t596 = qJD(6) * t604;
t595 = qJD(6) * t603;
t536 = t561 * t584 + (-t577 * t593 + t619) * t580;
t545 = t580 * t611 + (t573 * t584 + t576 * t625) * t575;
t586 = qJD(4) + t579 * t595 + (t637 * t579 + t600 * t583) * qJD(5);
t554 = t592 * qJD(1);
t541 = t546 * qJD(1);
t540 = t545 * qJD(3);
t539 = t544 * qJD(3);
t522 = t641 * qJD(5) - t540 * t579;
t517 = t640 * qJD(3) + t554 * t584 + t642 * t580;
t516 = t536 * qJD(3) + t554 * t580 - t642 * t584;
t514 = t602 * qJD(5) - t653;
t510 = t601 * qJD(5) + t516 * t579 + t541 * t583;
t509 = t528 * qJD(5) - t516 * t583 + t541 * t579;
t508 = t510 * t582 + t517 * t578 + (-t528 * t578 + t536 * t582) * qJD(6);
t507 = -t510 * t578 + t517 * t582 + (-t528 * t582 - t536 * t578) * qJD(6);
t1 = [-pkin(9) * t633 + t585 * t621 - t556 * pkin(2) - t543 * pkin(4) + t520 * qJ(4) - t638 * qJD(4) + t600 * ((-t583 * t638 - t646) * qJD(5) + t653) + t599 * t651 + t637 * (t524 * qJD(5) + t652) + ((-t524 * t578 + t649) * r_i_i_C(1) + (-t524 * t582 - t650) * r_i_i_C(2)) * qJD(6) + (-t585 * pkin(1) - t643 * t627) * qJD(1), t614, -t599 * t516 + t587 * t517 + t586 * t536 - t596 * t640, t516, -t600 * t509 + t637 * t510 + t601 * t595, t507 * r_i_i_C(1) - t508 * r_i_i_C(2); t574 * pkin(9) * t589 - pkin(1) * t623 + t554 * pkin(2) + t541 * pkin(4) + t510 * pkin(5) + t508 * r_i_i_C(1) + t507 * r_i_i_C(2) + t516 * qJ(4) - qJD(4) * t640 + t637 * t509 + t644 * t517 + t581 * t621 + t643 * t614, t615, t520 * t599 + t531 * t596 - t534 * t586 - t587 * t651, -t520, t637 * t514 + t602 * t595 + t600 * (-t526 * qJD(5) + t652) (-t514 * t578 - t582 * t651) * r_i_i_C(1) + (-t514 * t582 + t578 * t651) * r_i_i_C(2) + ((-t526 * t582 + t650) * r_i_i_C(1) + (t526 * t578 + t649) * r_i_i_C(2)) * qJD(6); 0, 0, -t587 * t539 - t599 * t540 + t544 * t596 + t586 * t545, t540, -t637 * t522 - t641 * t595 + t600 * (-t530 * qJD(5) + t540 * t583) (t522 * t578 - t539 * t582) * r_i_i_C(1) + (t522 * t582 + t539 * t578) * r_i_i_C(2) + ((-t530 * t582 - t545 * t578) * r_i_i_C(1) + (t530 * t578 - t545 * t582) * r_i_i_C(2)) * qJD(6);];
JaD_transl  = t1;
