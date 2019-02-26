% Zeitableitung der translatorischen Teilmatrix der analytischen Jacobi-Matrix für Segment Nr. 6 (0=Basis) von
% S6PPRRRR1
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
%   pkin=[a2,a3,a4,a5,a6,alpha2,alpha3,d3,d4,d5,d6,theta1,theta2]';
%
% Output:
% JaD_transl [3x6]
%   Zeitableitung der translatorischen Teilmatrix der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-02-26 19:42
% Revision: d75aae1ac561373cd3be920984c3df29a1c2ecc8 (2019-02-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JaD_transl = S6PPRRRR1_jacobiaD_transl_6_sym_varpar(qJ, qJD, r_i_i_C, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(3,1),zeros(13,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PPRRRR1_jacobiaD_transl_6_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6PPRRRR1_jacobiaD_transl_6_sym_varpar: qJD has to be [6x1] (double)');
assert(isa(r_i_i_C,'double') && isreal(r_i_i_C) && all(size(r_i_i_C) == [3 1]), ...
	'S6PPRRRR1_jacobiaD_transl_6_sym_varpar: Position vector r_i_i_C has to be [3x1] double');
assert(isreal(pkin) && all(size(pkin) == [13 1]), ...
  'S6PPRRRR1_jacobiaD_transl_6_sym_varpar: pkin has to be [13x1] (double)');

%% Symbolic Calculation
% From jacobiaD_transl_6_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 19:42:46
% EndTime: 2019-02-26 19:42:46
% DurationCPUTime: 0.62s
% Computational Cost: add. (1128->102), mult. (2701->183), div. (0->0), fcn. (3130->16), ass. (0->76)
t530 = sin(pkin(13));
t531 = sin(pkin(12));
t511 = t531 * t530;
t534 = cos(pkin(13));
t535 = cos(pkin(12));
t517 = t535 * t534;
t537 = cos(pkin(6));
t499 = -t537 * t517 + t511;
t536 = cos(pkin(7));
t497 = t499 * t536;
t532 = sin(pkin(7));
t533 = sin(pkin(6));
t515 = t533 * t532;
t507 = t535 * t515;
t546 = t497 + t507;
t513 = t531 * t534;
t516 = t535 * t530;
t500 = t537 * t513 + t516;
t512 = t531 * t533;
t545 = t500 * t536 - t532 * t512;
t518 = t536 * t533;
t544 = t534 * t518 + t537 * t532;
t543 = pkin(11) + r_i_i_C(3);
t488 = cos(qJ(6));
t525 = qJD(6) * t488;
t485 = sin(qJ(6));
t526 = qJD(6) * t485;
t542 = -r_i_i_C(1) * t526 - t525 * r_i_i_C(2);
t510 = t488 * r_i_i_C(1) - t485 * r_i_i_C(2) + pkin(5);
t471 = t537 * t516 + t513;
t487 = sin(qJ(3));
t539 = cos(qJ(3));
t455 = t471 * t487 + t546 * t539;
t541 = t545 * t539;
t514 = t533 * t530;
t462 = t487 * t514 - t544 * t539;
t484 = qJ(4) + qJ(5);
t481 = sin(t484);
t483 = qJD(4) + qJD(5);
t529 = t481 * t483;
t482 = cos(t484);
t528 = t482 * t483;
t527 = qJD(3) * t487;
t523 = t471 * t539;
t449 = t455 * qJD(3);
t464 = t499 * t532 - t535 * t518;
t522 = t464 * t483 - t449;
t472 = -t537 * t511 + t517;
t451 = t541 * qJD(3) + t472 * t527;
t465 = t500 * t532 + t536 * t512;
t521 = t465 * t483 - t451;
t460 = t462 * qJD(3);
t470 = -t534 * t515 + t537 * t536;
t520 = t470 * t483 - t460;
t489 = cos(qJ(4));
t502 = -t489 * pkin(4) - t481 * t543 - t510 * t482 - pkin(3);
t456 = -t546 * t487 + t523;
t435 = -t456 * t529 + t522 * t482;
t496 = t542 * (-t456 * t481 + t464 * t482) + t543 * t435 + t510 * (-t456 * t528 - t522 * t481);
t458 = t472 * t539 - t545 * t487;
t437 = -t458 * t529 + t521 * t482;
t495 = t542 * (-t458 * t481 + t465 * t482) + t543 * t437 + t510 * (-t458 * t528 - t521 * t481);
t463 = t544 * t487 + t539 * t514;
t442 = -t463 * t529 + t520 * t482;
t494 = t542 * (-t463 * t481 + t470 * t482) + t543 * t442 + t510 * (-t463 * t528 - t520 * t481);
t486 = sin(qJ(4));
t491 = qJD(4) * t486 * pkin(4) + (t485 * r_i_i_C(1) + t488 * r_i_i_C(2)) * t482 * qJD(6) + (t510 * t481 - t482 * t543) * t483;
t490 = -pkin(10) - pkin(9);
t461 = t463 * qJD(3);
t457 = t472 * t487 + t541;
t454 = t463 * t482 + t470 * t481;
t452 = t458 * qJD(3);
t450 = -t507 * t527 + (-t487 * t497 + t523) * qJD(3);
t446 = t458 * t482 + t465 * t481;
t444 = t456 * t482 + t464 * t481;
t1 = [0, 0 (-t451 * t485 + t458 * t525) * r_i_i_C(1) + (-t451 * t488 - t458 * t526) * r_i_i_C(2) + t451 * t490 + t502 * t452 + t491 * t457 (t451 * t486 + (-t458 * t489 - t465 * t486) * qJD(4)) * pkin(4) + t495, t495 (-t437 * t485 + t452 * t488) * r_i_i_C(1) + (-t437 * t488 - t452 * t485) * r_i_i_C(2) + ((-t446 * t488 - t457 * t485) * r_i_i_C(1) + (t446 * t485 - t457 * t488) * r_i_i_C(2)) * qJD(6); 0, 0 (-t449 * t485 + t456 * t525) * r_i_i_C(1) + (-t449 * t488 - t456 * t526) * r_i_i_C(2) + t449 * t490 + t502 * t450 + t491 * t455 (t449 * t486 + (-t456 * t489 - t464 * t486) * qJD(4)) * pkin(4) + t496, t496 (-t435 * t485 + t450 * t488) * r_i_i_C(1) + (-t435 * t488 - t450 * t485) * r_i_i_C(2) + ((-t444 * t488 - t455 * t485) * r_i_i_C(1) + (t444 * t485 - t455 * t488) * r_i_i_C(2)) * qJD(6); 0, 0 (-t460 * t485 + t463 * t525) * r_i_i_C(1) + (-t460 * t488 - t463 * t526) * r_i_i_C(2) + t460 * t490 + t502 * t461 + t491 * t462 (t460 * t486 + (-t463 * t489 - t470 * t486) * qJD(4)) * pkin(4) + t494, t494 (-t442 * t485 + t461 * t488) * r_i_i_C(1) + (-t442 * t488 - t461 * t485) * r_i_i_C(2) + ((-t454 * t488 - t462 * t485) * r_i_i_C(1) + (t454 * t485 - t462 * t488) * r_i_i_C(2)) * qJD(6);];
JaD_transl  = t1;
