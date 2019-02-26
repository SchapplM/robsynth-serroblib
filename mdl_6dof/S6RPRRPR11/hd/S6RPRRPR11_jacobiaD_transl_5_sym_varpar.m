% Zeitableitung der translatorischen Teilmatrix der analytischen Jacobi-Matrix für Segment Nr. 5 (0=Basis) von
% S6RPRRPR11
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
%   pkin=[a2,a3,a4,a5,a6,alpha2,alpha3,d1,d3,d4,d6,theta2,theta5]';
%
% Output:
% JaD_transl [3x6]
%   Zeitableitung der translatorischen Teilmatrix der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-02-26 21:06
% Revision: d75aae1ac561373cd3be920984c3df29a1c2ecc8 (2019-02-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JaD_transl = S6RPRRPR11_jacobiaD_transl_5_sym_varpar(qJ, qJD, r_i_i_C, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(3,1),zeros(13,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRRPR11_jacobiaD_transl_5_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RPRRPR11_jacobiaD_transl_5_sym_varpar: qJD has to be [6x1] (double)');
assert(isa(r_i_i_C,'double') && isreal(r_i_i_C) && all(size(r_i_i_C) == [3 1]), ...
	'S6RPRRPR11_jacobiaD_transl_5_sym_varpar: Position vector r_i_i_C has to be [3x1] double');
assert(isreal(pkin) && all(size(pkin) == [13 1]), ...
  'S6RPRRPR11_jacobiaD_transl_5_sym_varpar: pkin has to be [13x1] (double)');

%% Symbolic Calculation
% From jacobiaD_transl_5_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 21:06:43
% EndTime: 2019-02-26 21:06:44
% DurationCPUTime: 0.78s
% Computational Cost: add. (917->99), mult. (2979->171), div. (0->0), fcn. (3248->14), ass. (0->75)
t518 = sin(qJ(1));
t521 = cos(qJ(1));
t511 = sin(pkin(12));
t571 = cos(pkin(6));
t553 = t511 * t571;
t570 = cos(pkin(12));
t503 = -t518 * t553 + t521 * t570;
t498 = t503 * qJD(1);
t517 = sin(qJ(3));
t520 = cos(qJ(3));
t529 = -t518 * t570 - t521 * t553;
t542 = t571 * t570;
t530 = t521 * t511 + t518 * t542;
t497 = t530 * qJD(1);
t512 = sin(pkin(7));
t515 = cos(pkin(7));
t513 = sin(pkin(6));
t560 = qJD(1) * t518;
t555 = t513 * t560;
t532 = -t497 * t515 + t512 * t555;
t558 = qJD(3) * t520;
t563 = t513 * t512;
t548 = t558 * t563;
t528 = -t518 * t511 + t521 * t542;
t566 = t528 * t515;
t475 = (qJD(3) * t529 + t532) * t517 - (-qJD(3) * t566 - t498) * t520 - t521 * t548;
t516 = sin(qJ(4));
t519 = cos(qJ(4));
t567 = t497 * t512;
t533 = t515 * t555 + t567;
t561 = t513 * t521;
t535 = t512 * t561 - t566;
t482 = t535 * t517 + t520 * t529;
t492 = t528 * t512 + t515 * t561;
t581 = t482 * t519 + t492 * t516;
t586 = t581 * qJD(4) - t475 * t516 + t533 * t519;
t582 = t482 * t516 - t492 * t519;
t585 = t582 * qJD(4) + t475 * t519 + t533 * t516;
t572 = r_i_i_C(3) + qJ(5);
t526 = t528 * qJD(1);
t554 = qJD(1) * t561;
t576 = t512 * t554 - t515 * t526;
t562 = t513 * t518;
t564 = t530 * t515;
t534 = t512 * t562 - t564;
t484 = t503 * t520 + t534 * t517;
t494 = t512 * t530 + t515 * t562;
t575 = -t484 * t516 + t494 * t519;
t574 = pkin(9) * t515 + qJ(2);
t550 = t571 * t512;
t552 = t515 * t570;
t491 = (t511 * t520 + t517 * t552) * t513 + t517 * t550;
t510 = sin(pkin(13));
t514 = cos(pkin(13));
t537 = t514 * r_i_i_C(1) - t510 * r_i_i_C(2) + pkin(4);
t524 = t572 * t516 + t537 * t519 + pkin(3);
t559 = qJD(3) * t517;
t557 = t513 * qJD(2);
t545 = t520 * t550;
t544 = t520 * t552;
t539 = t484 * t519 + t494 * t516;
t499 = t571 * t515 - t570 * t563;
t538 = t491 * t519 + t499 * t516;
t536 = t510 * r_i_i_C(1) + t514 * r_i_i_C(2) + pkin(10);
t523 = t516 * qJD(5) + (-t537 * t516 + t572 * t519) * qJD(4);
t522 = t492 * qJD(1);
t476 = t482 * qJD(3) - t498 * t517 + t532 * t520;
t496 = t529 * qJD(1);
t487 = t513 * t511 * t559 + (-t513 * t544 - t545) * qJD(3);
t478 = t538 * qJD(4) - t487 * t516;
t473 = t496 * t520 - t503 * t559 + t576 * t517 + t518 * t548 - t558 * t564;
t472 = t484 * qJD(3) + t496 * t517 - t576 * t520;
t467 = t575 * qJD(4) + t473 * t519 + t516 * t522;
t466 = t539 * qJD(4) + t473 * t516 - t519 * t522;
t1 = [(t476 * t510 - t514 * t585) * r_i_i_C(1) + (t476 * t514 + t510 * t585) * r_i_i_C(2) - t585 * pkin(4) + t582 * qJD(5) - t475 * pkin(3) + t476 * pkin(10) - t498 * pkin(2) - pkin(9) * t567 + t521 * t557 + t572 * t586 + (-t521 * pkin(1) - t574 * t562) * qJD(1), t554, t536 * t473 + t523 * (-t503 * t517 + t534 * t520) - t524 * t472, t539 * qJD(5) - t537 * t466 + t572 * t467, t466, 0; (t467 * t514 + t472 * t510) * r_i_i_C(1) + (-t467 * t510 + t472 * t514) * r_i_i_C(2) + t467 * pkin(4) - t575 * qJD(5) + t473 * pkin(3) + t472 * pkin(10) + t496 * pkin(2) + t512 * pkin(9) * t526 - pkin(1) * t560 + t518 * t557 + t574 * t554 + t572 * t466, t555, t536 * t475 + t523 * (t517 * t529 - t535 * t520) + t524 * t476, -t581 * qJD(5) + t537 * t586 + t572 * t585, -t586, 0; 0, 0, -t536 * t487 + t523 * (t545 + (-t511 * t517 + t544) * t513) - t524 * t491 * qJD(3), t538 * qJD(5) + t572 * (-t487 * t519 + (-t491 * t516 + t499 * t519) * qJD(4)) - t537 * t478, t478, 0;];
JaD_transl  = t1;
