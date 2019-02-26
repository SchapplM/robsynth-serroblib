% Zeitableitung der translatorischen Teilmatrix der analytischen Jacobi-Matrix für Segment Nr. 6 (0=Basis) von
% S6PPPRRR1
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
%   pkin=[a2,a3,a4,a5,a6,alpha2,alpha3,alpha4,d4,d5,d6,theta1,theta2,theta3]';
%
% Output:
% JaD_transl [3x6]
%   Zeitableitung der translatorischen Teilmatrix der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-02-26 19:38
% Revision: d75aae1ac561373cd3be920984c3df29a1c2ecc8 (2019-02-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JaD_transl = S6PPPRRR1_jacobiaD_transl_6_sym_varpar(qJ, qJD, r_i_i_C, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(3,1),zeros(14,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PPPRRR1_jacobiaD_transl_6_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6PPPRRR1_jacobiaD_transl_6_sym_varpar: qJD has to be [6x1] (double)');
assert(isa(r_i_i_C,'double') && isreal(r_i_i_C) && all(size(r_i_i_C) == [3 1]), ...
	'S6PPPRRR1_jacobiaD_transl_6_sym_varpar: Position vector r_i_i_C has to be [3x1] double');
assert(isreal(pkin) && all(size(pkin) == [14 1]), ...
  'S6PPPRRR1_jacobiaD_transl_6_sym_varpar: pkin has to be [14x1] (double)');

%% Symbolic Calculation
% From jacobiaD_transl_6_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 19:38:48
% EndTime: 2019-02-26 19:38:49
% DurationCPUTime: 0.68s
% Computational Cost: add. (1582->99), mult. (4747->188), div. (0->0), fcn. (6032->18), ass. (0->84)
t585 = sin(pkin(13));
t586 = sin(pkin(12));
t571 = t586 * t585;
t591 = cos(pkin(13));
t592 = cos(pkin(12));
t578 = t592 * t591;
t595 = cos(pkin(6));
t558 = t595 * t578 - t571;
t594 = cos(pkin(7));
t555 = t558 * t594;
t572 = t586 * t591;
t577 = t592 * t585;
t559 = t595 * t577 + t572;
t588 = sin(pkin(7));
t590 = cos(pkin(14));
t574 = t588 * t590;
t589 = sin(pkin(6));
t564 = t589 * t574;
t584 = sin(pkin(14));
t541 = -t590 * t555 + t559 * t584 + t592 * t564;
t579 = t594 * t589;
t550 = -t558 * t588 - t592 * t579;
t587 = sin(pkin(8));
t593 = cos(pkin(8));
t603 = t541 * t593 - t550 * t587;
t560 = -t595 * t572 - t577;
t556 = t560 * t594;
t561 = -t595 * t571 + t578;
t542 = -t590 * t556 + t561 * t584 - t586 * t564;
t551 = -t560 * t588 + t586 * t579;
t602 = t542 * t593 - t551 * t587;
t576 = t589 * t591;
t565 = t594 * t576;
t575 = t589 * t585;
t549 = -t590 * t565 - t595 * t574 + t584 * t575;
t557 = -t588 * t576 + t595 * t594;
t601 = t549 * t593 - t557 * t587;
t596 = cos(qJ(4));
t600 = t603 * t596;
t599 = t602 * t596;
t598 = t601 * t596;
t532 = sin(qJ(6));
t535 = cos(qJ(6));
t566 = qJD(6) * (t532 * r_i_i_C(1) + t535 * r_i_i_C(2));
t597 = pkin(11) + r_i_i_C(3);
t534 = sin(qJ(4));
t583 = qJD(4) * t534;
t582 = qJD(6) * t532;
t581 = qJD(6) * t535;
t573 = t588 * t584;
t563 = t589 * t573;
t520 = t584 * t555 + t559 * t590 - t592 * t563;
t507 = t520 * t596 - t603 * t534;
t514 = t541 * t587 + t550 * t593;
t533 = sin(qJ(5));
t536 = cos(qJ(5));
t497 = t507 * t536 + t514 * t533;
t570 = -t507 * t533 + t514 * t536;
t521 = t584 * t556 + t561 * t590 + t586 * t563;
t509 = t521 * t596 - t602 * t534;
t515 = t542 * t587 + t551 * t593;
t499 = t509 * t536 + t515 * t533;
t569 = -t509 * t533 + t515 * t536;
t529 = t584 * t565 + t595 * t573 + t590 * t575;
t513 = t529 * t596 - t601 * t534;
t522 = t549 * t587 + t557 * t593;
t505 = t513 * t536 + t522 * t533;
t568 = -t513 * t533 + t522 * t536;
t567 = r_i_i_C(1) * t535 - r_i_i_C(2) * t532 + pkin(5);
t562 = -t597 * t533 - t567 * t536 - pkin(4);
t552 = t536 * t566 + (t567 * t533 - t597 * t536) * qJD(5);
t512 = t529 * t534 + t598;
t511 = t513 * qJD(4);
t510 = t598 * qJD(4) + t529 * t583;
t508 = t521 * t534 + t599;
t506 = t520 * t534 + t600;
t503 = t509 * qJD(4);
t502 = t599 * qJD(4) + t521 * t583;
t501 = t507 * qJD(4);
t500 = t600 * qJD(4) + t520 * t583;
t495 = t568 * qJD(5) - t510 * t536;
t493 = t569 * qJD(5) - t502 * t536;
t491 = t570 * qJD(5) - t500 * t536;
t1 = [0, 0, 0 (-t502 * t532 + t509 * t581) * r_i_i_C(1) + (-t502 * t535 - t509 * t582) * r_i_i_C(2) - t502 * pkin(10) + t562 * t503 + t552 * t508, t597 * t493 - t569 * t566 + t567 * (-t499 * qJD(5) + t502 * t533) (-t493 * t532 + t503 * t535) * r_i_i_C(1) + (-t493 * t535 - t503 * t532) * r_i_i_C(2) + ((-t499 * t535 - t508 * t532) * r_i_i_C(1) + (t499 * t532 - t508 * t535) * r_i_i_C(2)) * qJD(6); 0, 0, 0 (-t500 * t532 + t507 * t581) * r_i_i_C(1) + (-t500 * t535 - t507 * t582) * r_i_i_C(2) - t500 * pkin(10) + t562 * t501 + t552 * t506, t597 * t491 - t570 * t566 + t567 * (-t497 * qJD(5) + t500 * t533) (-t491 * t532 + t501 * t535) * r_i_i_C(1) + (-t491 * t535 - t501 * t532) * r_i_i_C(2) + ((-t497 * t535 - t506 * t532) * r_i_i_C(1) + (t497 * t532 - t506 * t535) * r_i_i_C(2)) * qJD(6); 0, 0, 0 (-t510 * t532 + t513 * t581) * r_i_i_C(1) + (-t510 * t535 - t513 * t582) * r_i_i_C(2) - t510 * pkin(10) + t562 * t511 + t552 * t512, t597 * t495 - t568 * t566 + t567 * (-t505 * qJD(5) + t510 * t533) (-t495 * t532 + t511 * t535) * r_i_i_C(1) + (-t495 * t535 - t511 * t532) * r_i_i_C(2) + ((-t505 * t535 - t512 * t532) * r_i_i_C(1) + (t505 * t532 - t512 * t535) * r_i_i_C(2)) * qJD(6);];
JaD_transl  = t1;
