% Zeitableitung der translatorischen Teilmatrix der analytischen Jacobi-Matrix für Segment Nr. 5 (0=Basis) von
% S6RPRRPR12
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
%   pkin=[a2,a3,a4,a5,a6,alpha2,alpha3,d1,d3,d4,d6,theta2]';
%
% Output:
% JaD_transl [3x6]
%   Zeitableitung der translatorischen Teilmatrix der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-02-26 21:07
% Revision: d75aae1ac561373cd3be920984c3df29a1c2ecc8 (2019-02-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JaD_transl = S6RPRRPR12_jacobiaD_transl_5_sym_varpar(qJ, qJD, r_i_i_C, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(3,1),zeros(12,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRRPR12_jacobiaD_transl_5_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RPRRPR12_jacobiaD_transl_5_sym_varpar: qJD has to be [6x1] (double)');
assert(isa(r_i_i_C,'double') && isreal(r_i_i_C) && all(size(r_i_i_C) == [3 1]), ...
	'S6RPRRPR12_jacobiaD_transl_5_sym_varpar: Position vector r_i_i_C has to be [3x1] double');
assert(isreal(pkin) && all(size(pkin) == [12 1]), ...
  'S6RPRRPR12_jacobiaD_transl_5_sym_varpar: pkin has to be [12x1] (double)');

%% Symbolic Calculation
% From jacobiaD_transl_5_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 21:07:19
% EndTime: 2019-02-26 21:07:20
% DurationCPUTime: 0.60s
% Computational Cost: add. (746->86), mult. (2410->145), div. (0->0), fcn. (2622->12), ass. (0->70)
t468 = cos(pkin(12));
t470 = cos(pkin(6));
t465 = sin(pkin(12));
t473 = sin(qJ(1));
t510 = t473 * t465;
t500 = t470 * t510;
t476 = cos(qJ(1));
t506 = qJD(1) * t476;
t453 = -qJD(1) * t500 + t468 * t506;
t508 = t476 * t465;
t509 = t473 * t468;
t456 = t470 * t508 + t509;
t472 = sin(qJ(3));
t475 = cos(qJ(3));
t486 = t470 * t509 + t508;
t452 = t486 * qJD(1);
t466 = sin(pkin(7));
t469 = cos(pkin(7));
t467 = sin(pkin(6));
t513 = t467 * t473;
t497 = qJD(1) * t513;
t485 = -t452 * t469 + t466 * t497;
t512 = t467 * t476;
t499 = t466 * t512;
t496 = t475 * t499;
t507 = t476 * t468;
t487 = t470 * t507 - t510;
t528 = t487 * t469;
t429 = -qJD(3) * t496 + (-qJD(3) * t456 + t485) * t472 - (-qJD(3) * t528 - t453) * t475;
t518 = t452 * t466;
t444 = t469 * t497 + t518;
t471 = sin(qJ(4));
t474 = cos(qJ(4));
t490 = -t528 + t499;
t537 = t490 * t472;
t437 = -t456 * t475 + t537;
t527 = t487 * t466 + t469 * t512;
t533 = t437 * t474 + t527 * t471;
t539 = t533 * qJD(4) - t429 * t471 + t444 * t474;
t534 = t437 * t471 - t527 * t474;
t538 = t534 * qJD(4) + t429 * t474 + t444 * t471;
t511 = t468 * t469;
t514 = t466 * t470;
t526 = (-t465 * t472 + t475 * t511) * t467 + t475 * t514;
t458 = -t500 + t507;
t489 = t466 * t513 - t469 * t486;
t439 = t458 * t475 + t489 * t472;
t449 = t466 * t486 + t469 * t513;
t524 = -t439 * t471 + t449 * t474;
t446 = (t465 * t475 + t472 * t511) * t467 + t472 * t514;
t521 = r_i_i_C(3) + qJ(5);
t522 = r_i_i_C(2) - pkin(4);
t483 = t521 * t471 - t522 * t474 + pkin(3);
t523 = r_i_i_C(1) + pkin(10);
t504 = t467 * qJD(2);
t492 = t439 * t474 + t449 * t471;
t454 = -t467 * t468 * t466 + t470 * t469;
t491 = t446 * t474 + t454 * t471;
t481 = -t458 * t472 + t489 * t475;
t479 = qJD(5) * t471 + (t522 * t471 + t521 * t474) * qJD(4);
t478 = t527 * qJD(1);
t477 = t437 * qJD(3) - t453 * t472 + t485 * t475;
t451 = t456 * qJD(1);
t441 = t526 * qJD(3);
t432 = t491 * qJD(4) + t441 * t471;
t427 = qJD(1) * t537 + t481 * qJD(3) - t451 * t475;
t426 = t439 * qJD(3) - t451 * t472 + (t475 * t528 - t496) * qJD(1);
t421 = t524 * qJD(4) + t427 * t474 + t471 * t478;
t420 = t492 * qJD(4) + t427 * t471 - t474 * t478;
t1 = [t534 * qJD(5) - t429 * pkin(3) - t453 * pkin(2) - pkin(9) * t518 + t476 * t504 + t523 * t477 + t522 * t538 + t521 * t539 + (-t476 * pkin(1) + (-pkin(9) * t469 - qJ(2)) * t513) * qJD(1), t467 * t506, -t483 * t426 + t523 * t427 + t479 * t481, t492 * qJD(5) + t522 * t420 + t521 * t421, t420, 0; -t524 * qJD(5) + t427 * pkin(3) - t451 * pkin(2) + t473 * t504 + t523 * t426 - t522 * t421 + t521 * t420 + (-t473 * pkin(1) + pkin(9) * t527 + qJ(2) * t512) * qJD(1), t497, t523 * t429 + t479 * (-t456 * t472 - t490 * t475) + t483 * t477, -t533 * qJD(5) + t521 * t538 - t522 * t539, -t539, 0; 0, 0, -t483 * t446 * qJD(3) + t523 * t441 + t479 * t526, t491 * qJD(5) + t521 * (t441 * t474 + (-t446 * t471 + t454 * t474) * qJD(4)) + t522 * t432, t432, 0;];
JaD_transl  = t1;
