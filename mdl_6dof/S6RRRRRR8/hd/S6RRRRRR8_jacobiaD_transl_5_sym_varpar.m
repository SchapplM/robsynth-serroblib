% Zeitableitung der translatorischen Teilmatrix der analytischen Jacobi-Matrix für Segment Nr. 5 (0=Basis) von
% S6RRRRRR8
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
%   pkin=[a2,a3,a4,a5,a6,alpha2,alpha3,d1,d2,d3,d4,d5,d6]';
%
% Output:
% JaD_transl [3x6]
%   Zeitableitung der translatorischen Teilmatrix der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-02-26 22:51
% Revision: d75aae1ac561373cd3be920984c3df29a1c2ecc8 (2019-02-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JaD_transl = S6RRRRRR8_jacobiaD_transl_5_sym_varpar(qJ, qJD, r_i_i_C, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(3,1),zeros(13,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRRRR8_jacobiaD_transl_5_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRRRRR8_jacobiaD_transl_5_sym_varpar: qJD has to be [6x1] (double)');
assert(isa(r_i_i_C,'double') && isreal(r_i_i_C) && all(size(r_i_i_C) == [3 1]), ...
	'S6RRRRRR8_jacobiaD_transl_5_sym_varpar: Position vector r_i_i_C has to be [3x1] double');
assert(isreal(pkin) && all(size(pkin) == [13 1]), ...
  'S6RRRRRR8_jacobiaD_transl_5_sym_varpar: pkin has to be [13x1] (double)');

%% Symbolic Calculation
% From jacobiaD_transl_5_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 22:51:35
% EndTime: 2019-02-26 22:51:36
% DurationCPUTime: 1.02s
% Computational Cost: add. (1151->163), mult. (3096->278), div. (0->0), fcn. (3240->14), ass. (0->97)
t543 = sin(qJ(1));
t542 = sin(qJ(2));
t608 = cos(pkin(6));
t578 = t543 * t608;
t566 = t542 * t578;
t586 = qJD(2) * t542;
t546 = cos(qJ(2));
t547 = cos(qJ(1));
t591 = t547 * t546;
t514 = -qJD(1) * t566 - t543 * t586 + (qJD(2) * t608 + qJD(1)) * t591;
t577 = t547 * t608;
t523 = t542 * t577 + t543 * t546;
t541 = sin(qJ(3));
t545 = cos(qJ(3));
t553 = t547 * t542 + t546 * t578;
t513 = t553 * qJD(1) + t523 * qJD(2);
t537 = sin(pkin(7));
t539 = cos(pkin(7));
t538 = sin(pkin(6));
t587 = qJD(1) * t538;
t580 = t543 * t587;
t556 = -t513 * t539 + t537 * t580;
t522 = t543 * t542 - t546 * t577;
t598 = t538 * t547;
t562 = t522 * t539 + t537 * t598;
t487 = (-qJD(3) * t523 + t556) * t541 + (-t562 * qJD(3) + t514) * t545;
t517 = -t522 * t537 + t539 * t598;
t535 = qJD(4) + qJD(5);
t612 = t517 * t535 - t487;
t500 = -t523 * t545 + t562 * t541;
t606 = t513 * t537;
t505 = t539 * t580 + t606;
t611 = -t500 * t535 - t505;
t609 = r_i_i_C(3) + pkin(12) + pkin(11);
t554 = t566 - t591;
t511 = t522 * qJD(1) + t554 * qJD(2);
t607 = t511 * t537;
t536 = qJ(4) + qJ(5);
t533 = sin(t536);
t603 = t533 * t535;
t534 = cos(t536);
t602 = t534 * t535;
t601 = t535 * t542;
t600 = t537 * t538;
t599 = t538 * t543;
t597 = t539 * t541;
t596 = t539 * t545;
t595 = t541 * t542;
t594 = t541 * t546;
t593 = t542 * t545;
t592 = t545 * t546;
t561 = t537 * t599 - t539 * t553;
t502 = t561 * t541 - t545 * t554;
t579 = t547 * t587;
t503 = t539 * t579 - t607;
t569 = -t502 * t535 + t503;
t512 = t523 * qJD(1) + t553 * qJD(2);
t550 = t541 * t554 + t561 * t545;
t568 = t537 * t579;
t485 = -t512 * t545 + (t511 * t539 + t568) * t541 + t550 * qJD(3);
t519 = t537 * t553 + t539 * t599;
t575 = t519 * t535 + t485;
t482 = -t575 * t533 + t569 * t534;
t483 = t569 * t533 + t575 * t534;
t590 = t482 * r_i_i_C(1) - t483 * r_i_i_C(2);
t589 = (t533 * t612 - t534 * t611) * r_i_i_C(1) + (t533 * t611 + t534 * t612) * r_i_i_C(2);
t558 = t539 * t594 + t593;
t576 = t608 * t537;
t516 = t558 * t538 + t541 * t576;
t567 = t586 * t600;
t555 = -t516 * t535 + t567;
t557 = t539 * t595 - t592;
t560 = t539 * t592 - t595;
t565 = qJD(3) * t576;
t497 = t545 * t565 + (-t557 * qJD(2) + t560 * qJD(3)) * t538;
t521 = t608 * t539 - t546 * t600;
t572 = -t521 * t535 - t497;
t588 = (t572 * t533 + t555 * t534) * r_i_i_C(1) + (-t555 * t533 + t572 * t534) * r_i_i_C(2);
t585 = qJD(2) * t546;
t544 = cos(qJ(4));
t584 = qJD(4) * t544;
t540 = sin(qJ(4));
t583 = qJD(4) * t540 * pkin(4);
t581 = pkin(10) * t539 + pkin(9);
t532 = t544 * pkin(4) + pkin(3);
t563 = t534 * r_i_i_C(1) - t533 * r_i_i_C(2) + t532;
t509 = -t522 * t545 - t523 * t597;
t510 = -t545 * t553 + t554 * t597;
t559 = -t539 * t593 - t594;
t551 = -t583 + (-t533 * r_i_i_C(1) - t534 * r_i_i_C(2)) * t535;
t549 = t500 * qJD(3) - t514 * t541 + t556 * t545;
t520 = t557 * t538;
t507 = (-t558 * qJD(2) + t559 * qJD(3)) * t538;
t495 = -t514 * t597 - t513 * t545 + (t522 * t541 - t523 * t596) * qJD(3);
t493 = t512 * t597 + t511 * t545 + (t541 * t553 + t554 * t596) * qJD(3);
t484 = t502 * qJD(3) - t511 * t596 - t512 * t541 - t545 * t568;
t1 = [-pkin(10) * t606 - t514 * pkin(2) - t487 * t532 + (r_i_i_C(1) * t612 + r_i_i_C(2) * t611) * t534 + (r_i_i_C(1) * t611 - r_i_i_C(2) * t612) * t533 + t609 * t549 + (-t547 * pkin(1) - t581 * t599) * qJD(1) + (-t505 * t540 + (-t500 * t540 + t517 * t544) * qJD(4)) * pkin(4) (t493 * t534 - t510 * t603) * r_i_i_C(1) + (-t493 * t533 - t510 * t602) * r_i_i_C(2) + t493 * t532 - t510 * t583 + t511 * pkin(2) + t609 * (t510 * qJD(3) + t511 * t541 - t512 * t596) + ((-t512 * t533 - t554 * t602) * r_i_i_C(1) + (-t512 * t534 + t554 * t603) * r_i_i_C(2) - t512 * pkin(10) + (-t512 * t540 - t554 * t584) * pkin(4)) * t537, -t563 * t484 + t609 * t485 + t551 * t550 (-t485 * t540 + t503 * t544 + (-t502 * t544 - t519 * t540) * qJD(4)) * pkin(4) + t590, t590, 0; -pkin(10) * t607 - t512 * pkin(2) + t483 * r_i_i_C(1) + t482 * r_i_i_C(2) + t485 * t532 + t609 * t484 + (-pkin(1) * t543 + t581 * t598) * qJD(1) + (t503 * t540 + (-t502 * t540 + t519 * t544) * qJD(4)) * pkin(4) (t495 * t534 - t509 * t603) * r_i_i_C(1) + (-t495 * t533 - t509 * t602) * r_i_i_C(2) + t495 * t532 - t509 * t583 - t513 * pkin(2) + t609 * (t509 * qJD(3) - t513 * t541 + t514 * t596) + ((t514 * t533 + t523 * t602) * r_i_i_C(1) + (t514 * t534 - t523 * t603) * r_i_i_C(2) + t514 * pkin(10) + (t514 * t540 + t523 * t584) * pkin(4)) * t537, t609 * t487 + t551 * (-t523 * t541 - t562 * t545) + t563 * t549 (-t487 * t540 + t505 * t544 + (t500 * t544 + t517 * t540) * qJD(4)) * pkin(4) + t589, t589, 0; 0 (t507 * t534 + t520 * t603) * r_i_i_C(1) + (-t507 * t533 + t520 * t602) * r_i_i_C(2) + t507 * t532 + t520 * t583 + (-t609 * (-t560 * qJD(2) + t557 * qJD(3)) - pkin(2) * t586 + ((t533 * t585 + t534 * t601) * r_i_i_C(1) + (-t533 * t601 + t534 * t585) * r_i_i_C(2) + pkin(10) * t585 + (t540 * t585 + t542 * t584) * pkin(4)) * t537) * t538, t609 * t497 + t551 * (t560 * t538 + t545 * t576) + t563 * (-t541 * t565 + (t559 * qJD(2) - t558 * qJD(3)) * t538) (t544 * t567 - t497 * t540 + (-t516 * t544 - t521 * t540) * qJD(4)) * pkin(4) + t588, t588, 0;];
JaD_transl  = t1;
