% Zeitableitung der translatorischen Teilmatrix der analytischen Jacobi-Matrix für Segment Nr. 6 (0=Basis) von
% S6RRRRRP10
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
% Datum: 2019-02-26 22:45
% Revision: d75aae1ac561373cd3be920984c3df29a1c2ecc8 (2019-02-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JaD_transl = S6RRRRRP10_jacobiaD_transl_6_sym_varpar(qJ, qJD, r_i_i_C, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(3,1),zeros(11,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRRRP10_jacobiaD_transl_6_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRRRRP10_jacobiaD_transl_6_sym_varpar: qJD has to be [6x1] (double)');
assert(isa(r_i_i_C,'double') && isreal(r_i_i_C) && all(size(r_i_i_C) == [3 1]), ...
	'S6RRRRRP10_jacobiaD_transl_6_sym_varpar: Position vector r_i_i_C has to be [3x1] double');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6RRRRRP10_jacobiaD_transl_6_sym_varpar: pkin has to be [11x1] (double)');

%% Symbolic Calculation
% From jacobiaD_transl_6_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 22:45:04
% EndTime: 2019-02-26 22:45:05
% DurationCPUTime: 1.06s
% Computational Cost: add. (1430->149), mult. (3116->234), div. (0->0), fcn. (3218->12), ass. (0->87)
t518 = qJ(4) + qJ(5);
t515 = sin(t518);
t523 = sin(qJ(1));
t582 = cos(pkin(6));
t543 = qJD(2) * t582 + qJD(1);
t522 = sin(qJ(2));
t556 = t523 * t582;
t548 = t522 * t556;
t570 = qJD(2) * t522;
t526 = cos(qJ(2));
t527 = cos(qJ(1));
t572 = t527 * t526;
t488 = -qJD(1) * t548 - t523 * t570 + t543 * t572;
t555 = t527 * t582;
t501 = t522 * t555 + t523 * t526;
t521 = sin(qJ(3));
t525 = cos(qJ(3));
t519 = sin(pkin(6));
t571 = qJD(1) * t519;
t563 = t523 * t571;
t574 = t519 * t527;
t564 = t525 * t574;
t474 = (-qJD(3) * t501 + t563) * t521 - qJD(3) * t564 + t488 * t525;
t547 = t526 * t555;
t573 = t523 * t522;
t500 = -t547 + t573;
t517 = qJD(4) + qJD(5);
t553 = t500 * t517 + t474;
t592 = t553 * t515;
t516 = cos(t518);
t591 = t553 * t516;
t586 = r_i_i_C(1) + pkin(5);
t583 = r_i_i_C(3) + qJ(6);
t524 = cos(qJ(4));
t514 = t524 * pkin(4) + pkin(3);
t520 = sin(qJ(4));
t585 = t520 * pkin(4);
t542 = -qJD(4) * t585 + t515 * qJD(6);
t568 = qJD(3) * t521;
t584 = r_i_i_C(2) + pkin(11) + pkin(10);
t590 = -t514 * t568 + (t584 * qJD(3) + t542) * t525;
t589 = t525 * t514 + t584 * t521 + pkin(2);
t588 = (qJD(2) * t525 - t517) * t522 + t526 * t568;
t534 = t583 * t515 + t586 * t516 + t514;
t580 = t515 * t517;
t579 = t516 * t517;
t578 = t517 * t525;
t577 = t519 * t523;
t576 = t519 * t525;
t575 = t519 * t526;
t569 = qJD(2) * t526;
t567 = qJD(4) * t524;
t566 = t516 * qJD(6);
t565 = t515 * t575;
t562 = t527 * t571;
t561 = t519 * t570;
t560 = t519 * t569;
t502 = t527 * t522 + t526 * t556;
t486 = t501 * qJD(1) + t502 * qJD(2);
t503 = -t548 + t572;
t541 = -t503 * t521 + t523 * t576;
t472 = t541 * qJD(3) - t486 * t525 + t521 * t562;
t554 = t502 * t517 + t472;
t487 = t502 * qJD(1) + t501 * qJD(2);
t494 = -t501 * t525 + t521 * t574;
t551 = t494 * t517 + t487;
t546 = t502 * t578 - t486;
t545 = t500 * t578 + t488;
t544 = (qJD(2) - t578) * t526;
t496 = t503 * t525 + t521 * t577;
t538 = -t519 * t522 * t521 + t582 * t525;
t499 = t582 * t521 + t522 * t576;
t485 = -qJD(1) * t547 - t527 * t569 + t543 * t573;
t456 = t485 * t516 + t496 * t579 + t554 * t515;
t457 = -t485 * t515 - t496 * t580 + t554 * t516;
t537 = -(-t496 * t516 - t502 * t515) * qJD(6) + t583 * t457 - t586 * t456;
t458 = -t487 * t516 - t494 * t579 + t592;
t536 = -(t494 * t516 - t500 * t515) * qJD(6) + t583 * (t487 * t515 + t494 * t580 + t591) - t586 * t458;
t491 = t538 * qJD(3) + t525 * t560;
t469 = t491 * t515 + t499 * t579 - t516 * t561 - t517 * t565;
t535 = -(-t499 * t516 + t565) * qJD(6) + t583 * (t491 * t516 - t499 * t580 + t515 * t561 - t575 * t579) - t586 * t469;
t533 = t485 * t525 + t502 * t568 + t503 * t517;
t532 = -t487 * t525 + t500 * t568 + t501 * t517;
t531 = (-t586 * t515 + t583 * t516) * t517 + t542;
t530 = t494 * qJD(3) - t488 * t521 + t525 * t563;
t471 = t496 * qJD(3) - t486 * t521 - t525 * t562;
t1 = [-(-t494 * t515 - t500 * t516) * qJD(6) - t474 * t514 - t488 * pkin(2) - t487 * pkin(9) + t584 * t530 + t586 * (-t551 * t515 - t591) + t583 * (t551 * t516 - t592) + (-t527 * pkin(1) - pkin(8) * t577) * qJD(1) + (-t487 * t520 + (-t494 * t520 - t500 * t524) * qJD(4)) * pkin(4), -t503 * t566 - t486 * pkin(9) + t586 * (t546 * t515 + t533 * t516) + t583 * (t533 * t515 - t546 * t516) + (-t486 * t520 + t503 * t567) * pkin(4) - t590 * t502 + t589 * t485, -t534 * t471 + t584 * t472 + t531 * t541 (-t472 * t520 - t485 * t524 + (-t496 * t524 - t502 * t520) * qJD(4)) * pkin(4) + t537, t537, t456; -(-t496 * t515 + t502 * t516) * qJD(6) + t472 * t514 - t486 * pkin(2) - t485 * pkin(9) + t584 * t471 + t586 * t457 + t583 * t456 + (-t523 * pkin(1) + pkin(8) * t574) * qJD(1) + (-t485 * t520 + (-t496 * t520 + t502 * t524) * qJD(4)) * pkin(4), -t501 * t566 + t488 * pkin(9) + t586 * (t545 * t515 + t532 * t516) + t583 * (t532 * t515 - t545 * t516) + (t488 * t520 + t501 * t567) * pkin(4) - t590 * t500 - t589 * t487, t584 * t474 + t534 * t530 + t531 * (-t501 * t521 - t564) (-t474 * t520 + t487 * t524 + (t494 * t524 - t500 * t520) * qJD(4)) * pkin(4) + t536, t536, t458; 0 (t586 * (t515 * t544 - t588 * t516) - t583 * (t588 * t515 + t516 * t544) + (pkin(4) * t567 - t566) * t522 + t590 * t526 + ((pkin(9) + t585) * t526 - t589 * t522) * qJD(2)) * t519, t584 * t491 + t534 * (-t499 * qJD(3) - t521 * t560) + t531 * t538 (t524 * t561 - t491 * t520 + (-t499 * t524 + t520 * t575) * qJD(4)) * pkin(4) + t535, t535, t469;];
JaD_transl  = t1;
