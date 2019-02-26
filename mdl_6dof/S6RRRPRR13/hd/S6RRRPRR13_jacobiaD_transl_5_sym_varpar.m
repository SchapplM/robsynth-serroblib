% Zeitableitung der translatorischen Teilmatrix der analytischen Jacobi-Matrix für Segment Nr. 5 (0=Basis) von
% S6RRRPRR13
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
% Datum: 2019-02-26 22:23
% Revision: d75aae1ac561373cd3be920984c3df29a1c2ecc8 (2019-02-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JaD_transl = S6RRRPRR13_jacobiaD_transl_5_sym_varpar(qJ, qJD, r_i_i_C, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(3,1),zeros(13,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRPRR13_jacobiaD_transl_5_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRRPRR13_jacobiaD_transl_5_sym_varpar: qJD has to be [6x1] (double)');
assert(isa(r_i_i_C,'double') && isreal(r_i_i_C) && all(size(r_i_i_C) == [3 1]), ...
	'S6RRRPRR13_jacobiaD_transl_5_sym_varpar: Position vector r_i_i_C has to be [3x1] double');
assert(isreal(pkin) && all(size(pkin) == [13 1]), ...
  'S6RRRPRR13_jacobiaD_transl_5_sym_varpar: pkin has to be [13x1] (double)');

%% Symbolic Calculation
% From jacobiaD_transl_5_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 22:23:07
% EndTime: 2019-02-26 22:23:08
% DurationCPUTime: 0.96s
% Computational Cost: add. (904->146), mult. (2538->248), div. (0->0), fcn. (2654->14), ass. (0->97)
t520 = sin(qJ(1));
t519 = sin(qJ(2));
t575 = cos(pkin(6));
t551 = t520 * t575;
t545 = t519 * t551;
t559 = qJD(2) * t519;
t522 = cos(qJ(2));
t523 = cos(qJ(1));
t561 = t523 * t522;
t488 = -qJD(1) * t545 - t520 * t559 + (qJD(2) * t575 + qJD(1)) * t561;
t514 = sin(pkin(7));
t515 = sin(pkin(6));
t568 = t515 * t523;
t556 = t514 * t568;
t586 = -qJD(3) * t556 + t488;
t550 = t523 * t575;
t497 = t519 * t550 + t520 * t522;
t527 = t523 * t519 + t522 * t551;
t487 = t527 * qJD(1) + t497 * qJD(2);
t516 = cos(pkin(7));
t518 = sin(qJ(3));
t521 = cos(qJ(3));
t560 = qJD(1) * t515;
t554 = t520 * t560;
t549 = t514 * t554;
t496 = t520 * t519 - t522 * t550;
t572 = t496 * t516;
t462 = (-qJD(3) * t497 - t487 * t516 + t549) * t518 + (-qJD(3) * t572 + t586) * t521;
t573 = t487 * t514;
t478 = t516 * t554 + t573;
t512 = pkin(13) + qJ(5);
t510 = sin(t512);
t511 = cos(t512);
t585 = t462 * t510 - t478 * t511;
t584 = -t462 * t511 - t478 * t510;
t567 = t516 * t518;
t538 = t496 * t567 - t497 * t521;
t473 = t518 * t556 + t538;
t491 = -t496 * t514 + t516 * t568;
t581 = t473 * t511 + t491 * t510;
t580 = -t473 * t510 + t491 * t511;
t579 = (t556 + t572) * t521 + t497 * t518;
t566 = t516 * t521;
t578 = t538 * qJD(3) - t487 * t566 - t518 * t586 + t521 * t549;
t577 = pkin(4) * sin(pkin(13));
t576 = r_i_i_C(3) + pkin(11) + qJ(4);
t528 = t545 - t561;
t485 = t496 * qJD(1) + t528 * qJD(2);
t574 = t485 * t514;
t570 = t514 * t515;
t569 = t515 * t520;
t565 = t518 * t519;
t564 = t518 * t522;
t563 = t519 * t521;
t562 = t521 * t522;
t558 = qJD(5) * t510;
t557 = qJD(5) * t511;
t555 = pkin(10) * t516 + pkin(9);
t553 = t523 * t560;
t552 = t514 * t575;
t548 = t514 * t553;
t547 = t559 * t570;
t544 = qJD(3) * t552;
t543 = r_i_i_C(1) * t511 - r_i_i_C(2) * t510;
t542 = -r_i_i_C(1) * t510 - r_i_i_C(2) * t511;
t509 = cos(pkin(13)) * pkin(4) + pkin(3);
t540 = -t509 - t543;
t539 = t496 * t518 - t497 * t566;
t483 = -t496 * t521 - t497 * t567;
t536 = t518 * t527 + t528 * t566;
t484 = -t521 * t527 + t528 * t567;
t535 = t514 * t569 - t516 * t527;
t534 = t516 * t562 - t565;
t533 = t516 * t563 + t564;
t532 = t516 * t564 + t563;
t531 = t516 * t565 - t562;
t530 = qJD(5) * t543;
t529 = qJD(5) * t542;
t526 = pkin(10) - t542 + t577;
t474 = t518 * t528 + t535 * t521;
t475 = t535 * t518 - t521 * t528;
t495 = t575 * t516 - t522 * t570;
t494 = t531 * t515;
t493 = t514 * t527 + t516 * t569;
t490 = t532 * t515 + t518 * t552;
t486 = t497 * qJD(1) + t527 * qJD(2);
t480 = (-t532 * qJD(2) - t533 * qJD(3)) * t515;
t476 = t516 * t553 - t574;
t470 = t521 * t544 + (-t531 * qJD(2) + t534 * qJD(3)) * t515;
t469 = t518 * t544 + (t533 * qJD(2) + t532 * qJD(3)) * t515;
t468 = t539 * qJD(3) - t487 * t521 - t488 * t567;
t466 = t536 * qJD(3) + t485 * t521 + t486 * t567;
t460 = -t486 * t521 + (t485 * t516 + t548) * t518 + t474 * qJD(3);
t459 = t475 * qJD(3) - t485 * t566 - t486 * t518 - t521 * t548;
t458 = t460 * t511 + t476 * t510 + (-t475 * t510 + t493 * t511) * qJD(5);
t457 = -t460 * t510 + t476 * t511 + (-t475 * t511 - t493 * t510) * qJD(5);
t1 = [t584 * r_i_i_C(1) + t585 * r_i_i_C(2) - t462 * t509 - t579 * qJD(4) - t478 * t577 - t488 * pkin(2) - pkin(10) * t573 + t576 * t578 + (t580 * r_i_i_C(1) - t581 * r_i_i_C(2)) * qJD(5) + (-t523 * pkin(1) - t555 * t569) * qJD(1) (t466 * t511 - t484 * t558) * r_i_i_C(1) + (-t466 * t510 - t484 * t557) * r_i_i_C(2) + t466 * t509 - t536 * qJD(4) + t485 * pkin(2) + t576 * (t484 * qJD(3) + t485 * t518 - t486 * t566) + (-t526 * t486 - t528 * t530) * t514, qJD(4) * t475 + t540 * t459 + t576 * t460 + t474 * t529, t459, r_i_i_C(1) * t457 - t458 * r_i_i_C(2), 0; t476 * t577 - pkin(10) * t574 - t486 * pkin(2) + t458 * r_i_i_C(1) + t457 * r_i_i_C(2) - t474 * qJD(4) + t460 * t509 + t576 * t459 + (-pkin(1) * t520 + t555 * t568) * qJD(1) (t468 * t511 - t483 * t558) * r_i_i_C(1) + (-t468 * t510 - t483 * t557) * r_i_i_C(2) + t468 * t509 - t539 * qJD(4) - t487 * pkin(2) + t576 * (t483 * qJD(3) - t487 * t518 + t488 * t566) + (t526 * t488 + t497 * t530) * t514, -qJD(4) * t473 + t576 * t462 - t579 * t529 - t540 * t578, -t578, -t585 * r_i_i_C(1) + t584 * r_i_i_C(2) + (t581 * r_i_i_C(1) + t580 * r_i_i_C(2)) * qJD(5), 0; 0 (t480 * t511 + t494 * t558) * r_i_i_C(1) + (-t480 * t510 + t494 * t557) * r_i_i_C(2) + t480 * t509 + (-t576 * (-t534 * qJD(2) + t531 * qJD(3)) + t533 * qJD(4) - pkin(2) * t559 + (t526 * t522 * qJD(2) + t519 * t530) * t514) * t515, qJD(4) * t490 + t576 * t470 + (t534 * t515 + t521 * t552) * t529 + t540 * t469, t469 (-t470 * t510 + t511 * t547) * r_i_i_C(1) + (-t470 * t511 - t510 * t547) * r_i_i_C(2) + ((-t490 * t511 - t495 * t510) * r_i_i_C(1) + (t490 * t510 - t495 * t511) * r_i_i_C(2)) * qJD(5), 0;];
JaD_transl  = t1;
