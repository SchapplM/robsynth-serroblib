% Zeitableitung der translatorischen Teilmatrix der analytischen Jacobi-Matrix für Segment Nr. 6 (0=Basis) von
% S6RRPRRR4
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
%   pkin=[a2,a3,a4,a5,a6,alpha2,d1,d2,d4,d5,d6,theta3]';
%
% Output:
% JaD_transl [3x6]
%   Zeitableitung der translatorischen Teilmatrix der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-02-26 21:55
% Revision: d75aae1ac561373cd3be920984c3df29a1c2ecc8 (2019-02-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JaD_transl = S6RRPRRR4_jacobiaD_transl_6_sym_varpar(qJ, qJD, r_i_i_C, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(3,1),zeros(12,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPRRR4_jacobiaD_transl_6_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRPRRR4_jacobiaD_transl_6_sym_varpar: qJD has to be [6x1] (double)');
assert(isa(r_i_i_C,'double') && isreal(r_i_i_C) && all(size(r_i_i_C) == [3 1]), ...
	'S6RRPRRR4_jacobiaD_transl_6_sym_varpar: Position vector r_i_i_C has to be [3x1] double');
assert(isreal(pkin) && all(size(pkin) == [12 1]), ...
  'S6RRPRRR4_jacobiaD_transl_6_sym_varpar: pkin has to be [12x1] (double)');

%% Symbolic Calculation
% From jacobiaD_transl_6_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 21:55:49
% EndTime: 2019-02-26 21:55:50
% DurationCPUTime: 1.04s
% Computational Cost: add. (1457->150), mult. (3222->250), div. (0->0), fcn. (3524->14), ass. (0->91)
t579 = r_i_i_C(3) + pkin(11);
t520 = cos(qJ(6));
t558 = qJD(6) * t520;
t516 = sin(qJ(6));
t559 = qJD(6) * t516;
t589 = -r_i_i_C(1) * t559 - t558 * r_i_i_C(2);
t512 = qJ(4) + qJ(5);
t509 = sin(t512);
t510 = cos(t512);
t515 = cos(pkin(6));
t513 = sin(pkin(12));
t518 = sin(qJ(2));
t575 = cos(pkin(12));
t547 = t518 * t575;
t578 = cos(qJ(2));
t534 = t513 * t578 + t547;
t493 = t534 * t515;
t519 = sin(qJ(1));
t522 = cos(qJ(1));
t543 = t578 * t575;
t533 = -t518 * t513 + t543;
t540 = t522 * t493 + t519 * t533;
t514 = sin(pkin(6));
t567 = t514 * t522;
t468 = t509 * t567 - t510 * t540;
t532 = t533 * t515;
t474 = -t519 * t534 + t522 * t532;
t588 = -t468 * t516 + t474 * t520;
t587 = t468 * t520 + t474 * t516;
t511 = qJD(4) + qJD(5);
t517 = sin(qJ(4));
t585 = -t520 * r_i_i_C(1) + t516 * r_i_i_C(2) - pkin(5);
t586 = (-t509 * t585 - t510 * t579) * t511 + (t516 * r_i_i_C(1) + t520 * r_i_i_C(2)) * t510 * qJD(6) + qJD(4) * t517 * pkin(4);
t584 = qJD(1) * t532 + t533 * qJD(2);
t561 = qJD(2) * t518;
t583 = -qJD(2) * t543 + t513 * t561;
t488 = t583 * t515;
t548 = t578 * qJD(2);
t495 = -qJD(2) * t547 - t513 * t548;
t563 = qJD(1) * t519;
t458 = t493 * t563 - t519 * t495 + t522 * (-qJD(1) * t533 + t488);
t521 = cos(qJ(4));
t507 = t521 * pkin(4) + pkin(3);
t580 = t509 * t579 - t510 * t585 + t507;
t577 = pkin(2) * t515;
t489 = qJD(2) * t493;
t562 = qJD(1) * t522;
t459 = -t522 * t489 - t584 * t519 - t534 * t562;
t574 = t459 * t516;
t573 = t459 * t520;
t570 = t509 * t511;
t569 = t510 * t511;
t568 = t514 * t519;
t565 = t518 * t519;
t564 = t518 * t522;
t557 = pkin(2) * t561;
t554 = t510 * t567;
t553 = t578 * t519;
t552 = t578 * t522;
t551 = t514 * t563;
t550 = t514 * t562;
t546 = t458 * t510 + t511 * t554;
t486 = t583 * t514;
t545 = t511 * t515 - t486;
t528 = -qJD(1) * t540 + t519 * t488 + t522 * t495;
t541 = t511 * t568 + t528;
t539 = -t519 * t493 + t522 * t533;
t535 = -t511 * t540 + t551;
t492 = t534 * t514;
t451 = t535 * t510 + (t511 * t567 + t458) * t509;
t452 = t509 * t551 - t540 * t570 - t546;
t527 = t589 * (-t509 * t540 - t554) + t579 * t452 - t585 * t451;
t449 = t509 * t541 - t510 * t550 + t539 * t569;
t450 = t509 * t550 + t510 * t541 - t539 * t570;
t526 = t589 * (-t509 * t539 + t510 * t568) + t579 * t450 + t585 * t449;
t465 = -t492 * t570 + t510 * t545;
t525 = t589 * (-t492 * t509 + t515 * t510) + t579 * t465 - t585 * (-t492 * t569 - t509 * t545);
t523 = -pkin(10) - pkin(9);
t508 = pkin(2) * t578 + pkin(1);
t496 = -t514 * qJD(3) + t548 * t577;
t494 = t518 * t577 + (-pkin(8) - qJ(3)) * t514;
t491 = t533 * t514;
t487 = qJD(2) * t492;
t480 = t492 * t510 + t515 * t509;
t477 = -t519 * t532 - t522 * t534;
t470 = t509 * t568 + t510 * t539;
t456 = -t519 * t489 + t584 * t522 - t534 * t563;
t454 = -t509 * t535 + t546;
t442 = t450 * t520 + t456 * t516 + (-t470 * t516 - t477 * t520) * qJD(6);
t441 = -t450 * t516 + t456 * t520 + (-t470 * t520 + t477 * t516) * qJD(6);
t1 = [(t454 * t520 + t574) * r_i_i_C(1) + (-t454 * t516 + t573) * r_i_i_C(2) + t454 * pkin(5) + t458 * t507 - t459 * t523 + t519 * t557 - t522 * t496 + t579 * t451 + (t588 * r_i_i_C(1) - t587 * r_i_i_C(2)) * qJD(6) + (t519 * t494 - t522 * t508) * qJD(1) + (-t517 * t551 + (t517 * t540 + t521 * t567) * qJD(4)) * pkin(4) (t516 * t528 + t539 * t558) * r_i_i_C(1) + (t520 * t528 - t539 * t559) * r_i_i_C(2) - t528 * t523 - t580 * t456 + ((t515 * t565 - t552) * qJD(2) + (-t515 * t552 + t565) * qJD(1)) * pkin(2) - t586 * t477, t550 (t521 * t550 - t528 * t517 + (-t517 * t568 - t521 * t539) * qJD(4)) * pkin(4) + t526, t526, t441 * r_i_i_C(1) - t442 * r_i_i_C(2); -t522 * t557 + t450 * pkin(5) + t442 * r_i_i_C(1) + t441 * r_i_i_C(2) - t456 * t523 + t528 * t507 - t519 * t496 + t579 * t449 + (-t494 * t522 - t508 * t519) * qJD(1) + (t517 * t550 + (-t517 * t539 + t521 * t568) * qJD(4)) * pkin(4) (-t458 * t516 + t540 * t558) * r_i_i_C(1) + (-t458 * t520 - t540 * t559) * r_i_i_C(2) + t458 * t523 + t580 * t459 + ((-t515 * t564 - t553) * qJD(2) + (-t515 * t553 - t564) * qJD(1)) * pkin(2) - t586 * t474, t551 (t521 * t551 + t458 * t517 + (t517 * t567 - t521 * t540) * qJD(4)) * pkin(4) + t527, t527 (-t452 * t516 - t573) * r_i_i_C(1) + (-t452 * t520 + t574) * r_i_i_C(2) + (t587 * r_i_i_C(1) + t588 * r_i_i_C(2)) * qJD(6); 0 (-t486 * t516 + t492 * t558) * r_i_i_C(1) + (-t486 * t520 - t492 * t559) * r_i_i_C(2) + t486 * t523 - t514 * t557 - t580 * t487 - t586 * t491, 0 (t486 * t517 + (-t492 * t521 - t515 * t517) * qJD(4)) * pkin(4) + t525, t525 (-t465 * t516 + t487 * t520) * r_i_i_C(1) + (-t465 * t520 - t487 * t516) * r_i_i_C(2) + ((-t480 * t520 + t491 * t516) * r_i_i_C(1) + (t480 * t516 + t491 * t520) * r_i_i_C(2)) * qJD(6);];
JaD_transl  = t1;
