% Zeitableitung der translatorischen Teilmatrix der analytischen Jacobi-Matrix für Segment Nr. 5 (0=Basis) von
% S6PRRRPR7
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
%   pkin=[a2,a3,a4,a5,a6,alpha2,alpha3,d2,d3,d4,d6,theta1,theta5]';
%
% Output:
% JaD_transl [3x6]
%   Zeitableitung der translatorischen Teilmatrix der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-02-26 20:14
% Revision: d75aae1ac561373cd3be920984c3df29a1c2ecc8 (2019-02-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JaD_transl = S6PRRRPR7_jacobiaD_transl_5_sym_varpar(qJ, qJD, r_i_i_C, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(3,1),zeros(13,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRRRPR7_jacobiaD_transl_5_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6PRRRPR7_jacobiaD_transl_5_sym_varpar: qJD has to be [6x1] (double)');
assert(isa(r_i_i_C,'double') && isreal(r_i_i_C) && all(size(r_i_i_C) == [3 1]), ...
	'S6PRRRPR7_jacobiaD_transl_5_sym_varpar: Position vector r_i_i_C has to be [3x1] double');
assert(isreal(pkin) && all(size(pkin) == [13 1]), ...
  'S6PRRRPR7_jacobiaD_transl_5_sym_varpar: pkin has to be [13x1] (double)');

%% Symbolic Calculation
% From jacobiaD_transl_5_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 20:14:07
% EndTime: 2019-02-26 20:14:08
% DurationCPUTime: 0.62s
% Computational Cost: add. (893->134), mult. (2926->236), div. (0->0), fcn. (3175->14), ass. (0->87)
t485 = sin(pkin(12));
t489 = cos(pkin(12));
t494 = sin(qJ(2));
t491 = cos(pkin(6));
t497 = cos(qJ(2));
t529 = t491 * t497;
t475 = -t485 * t494 + t489 * t529;
t493 = sin(qJ(3));
t496 = cos(qJ(3));
t530 = t491 * t494;
t507 = t485 * t530 - t489 * t497;
t490 = cos(pkin(7));
t508 = t485 * t529 + t489 * t494;
t486 = sin(pkin(7));
t487 = sin(pkin(6));
t538 = t486 * t487;
t509 = t485 * t538 - t490 * t508;
t460 = t509 * t493 - t496 * t507;
t476 = t485 * t497 + t489 * t530;
t510 = -t475 * t490 + t489 * t538;
t544 = -t476 * t496 + t510 * t493;
t543 = pkin(9) * t486;
t542 = r_i_i_C(3) + qJ(5);
t537 = t486 * t491;
t492 = sin(qJ(4));
t536 = t486 * t492;
t495 = cos(qJ(4));
t535 = t486 * t495;
t534 = t486 * t497;
t533 = t487 * t490;
t532 = t490 * t493;
t531 = t490 * t496;
t528 = t493 * t494;
t527 = t493 * t497;
t526 = t494 * t496;
t525 = t496 * t497;
t524 = qJD(2) * t487;
t523 = t494 * t538;
t521 = t497 * t524;
t520 = qJD(3) * t537;
t519 = qJD(2) * t523;
t518 = t486 * t521;
t467 = -t475 * t486 - t489 * t533;
t517 = t467 * t492 - t495 * t544;
t468 = t485 * t533 + t486 * t508;
t516 = t460 * t495 + t468 * t492;
t504 = t490 * t527 + t526;
t466 = t504 * t487 + t493 * t537;
t474 = -t487 * t534 + t491 * t490;
t515 = t466 * t495 + t474 * t492;
t484 = sin(pkin(13));
t488 = cos(pkin(13));
t514 = t488 * r_i_i_C(1) - t484 * r_i_i_C(2) + pkin(4);
t513 = t484 * r_i_i_C(1) + t488 * r_i_i_C(2) + pkin(10);
t463 = t475 * t496 - t476 * t532;
t512 = -t463 * t492 + t476 * t535;
t464 = -t496 * t508 + t507 * t532;
t511 = -t464 * t492 - t507 * t535;
t506 = t490 * t525 - t528;
t505 = -t490 * t526 - t527;
t503 = -t490 * t528 + t525;
t469 = t503 * t487;
t502 = -t469 * t492 + t495 * t523;
t501 = -t476 * t493 - t510 * t496;
t500 = t493 * t507 + t509 * t496;
t499 = t542 * t492 + t514 * t495 + pkin(3);
t498 = t492 * qJD(5) + (-t514 * t492 + t542 * t495) * qJD(4);
t473 = t507 * qJD(2);
t472 = t508 * qJD(2);
t471 = t476 * qJD(2);
t470 = t475 * qJD(2);
t462 = (-t504 * qJD(2) + t505 * qJD(3)) * t487;
t461 = -t521 * t531 + (-qJD(3) * t525 + (qJD(3) * t490 + qJD(2)) * t528) * t487;
t456 = t496 * t520 + (t503 * qJD(2) + t506 * qJD(3)) * t487;
t454 = t472 * t532 + t473 * t496 + (t493 * t508 + t507 * t531) * qJD(3);
t453 = t464 * qJD(3) - t472 * t531 + t473 * t493;
t452 = -t470 * t532 - t471 * t496 + (-t475 * t493 - t476 * t531) * qJD(3);
t451 = t463 * qJD(3) + t470 * t531 - t471 * t493;
t450 = t502 * qJD(4) + t462 * t495 + t492 * t518;
t448 = t500 * qJD(3) - t472 * t496 + t473 * t532;
t446 = t501 * qJD(3) + t470 * t496 - t471 * t532;
t443 = t515 * qJD(4) + t456 * t492 - t495 * t519;
t442 = t511 * qJD(4) + t454 * t495 - t472 * t536;
t440 = t512 * qJD(4) + t452 * t495 + t470 * t536;
t437 = t516 * qJD(4) + t448 * t492 + t473 * t535;
t435 = t517 * qJD(4) + t446 * t492 - t471 * t535;
t1 = [0 (t442 * t488 + t453 * t484) * r_i_i_C(1) + (-t442 * t484 + t453 * t488) * r_i_i_C(2) + t442 * pkin(4) - t511 * qJD(5) + t454 * pkin(3) + t453 * pkin(10) + t473 * pkin(2) - t472 * t543 + t542 * (t472 * t535 + t454 * t492 + (t464 * t495 - t507 * t536) * qJD(4)) t513 * t448 + t498 * t500 + t499 * (-t460 * qJD(3) + t472 * t493 + t473 * t531) t516 * qJD(5) + t542 * (-t473 * t536 + t448 * t495 + (-t460 * t492 + t468 * t495) * qJD(4)) - t514 * t437, t437, 0; 0 (t440 * t488 + t451 * t484) * r_i_i_C(1) + (-t440 * t484 + t451 * t488) * r_i_i_C(2) + t440 * pkin(4) - t512 * qJD(5) + t452 * pkin(3) + t451 * pkin(10) - t471 * pkin(2) + t470 * t543 + t542 * (-t470 * t535 + t452 * t492 + (t463 * t495 + t476 * t536) * qJD(4)) t513 * t446 + t498 * t501 + t499 * (t544 * qJD(3) - t470 * t493 - t471 * t531) t517 * qJD(5) + t542 * (t471 * t536 + t446 * t495 + (t467 * t495 + t492 * t544) * qJD(4)) - t514 * t435, t435, 0; 0 (t450 * t488 - t461 * t484) * r_i_i_C(1) + (-t450 * t484 - t461 * t488) * r_i_i_C(2) + t450 * pkin(4) - t502 * qJD(5) + t462 * pkin(3) - t461 * pkin(10) + t542 * (-t495 * t518 + t462 * t492 + (t469 * t495 + t492 * t523) * qJD(4)) + (-pkin(2) * t494 + pkin(9) * t534) * t524, t513 * t456 + t498 * (t506 * t487 + t496 * t537) + t499 * (-t493 * t520 + (t505 * qJD(2) - t504 * qJD(3)) * t487) t515 * qJD(5) + t542 * (t492 * t519 + t456 * t495 + (-t466 * t492 + t474 * t495) * qJD(4)) - t514 * t443, t443, 0;];
JaD_transl  = t1;
