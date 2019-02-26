% Zeitableitung der translatorischen Teilmatrix der analytischen Jacobi-Matrix für Segment Nr. 6 (0=Basis) von
% S6RRRPRP7
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
%   pkin=[a2,a3,a4,a5,a6,alpha2,d1,d2,d3,d5,theta4]';
%
% Output:
% JaD_transl [3x6]
%   Zeitableitung der translatorischen Teilmatrix der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-02-26 22:12
% Revision: d75aae1ac561373cd3be920984c3df29a1c2ecc8 (2019-02-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JaD_transl = S6RRRPRP7_jacobiaD_transl_6_sym_varpar(qJ, qJD, r_i_i_C, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(3,1),zeros(11,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRPRP7_jacobiaD_transl_6_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRRPRP7_jacobiaD_transl_6_sym_varpar: qJD has to be [6x1] (double)');
assert(isa(r_i_i_C,'double') && isreal(r_i_i_C) && all(size(r_i_i_C) == [3 1]), ...
	'S6RRRPRP7_jacobiaD_transl_6_sym_varpar: Position vector r_i_i_C has to be [3x1] double');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6RRRPRP7_jacobiaD_transl_6_sym_varpar: pkin has to be [11x1] (double)');

%% Symbolic Calculation
% From jacobiaD_transl_6_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 22:12:43
% EndTime: 2019-02-26 22:12:44
% DurationCPUTime: 0.95s
% Computational Cost: add. (1181->144), mult. (2491->226), div. (0->0), fcn. (2555->12), ass. (0->80)
t471 = cos(pkin(6));
t476 = sin(qJ(1));
t475 = sin(qJ(2));
t519 = t476 * t475;
t509 = t471 * t519;
t514 = qJD(2) * t475;
t479 = cos(qJ(2));
t480 = cos(qJ(1));
t516 = t480 * t479;
t445 = -qJD(1) * t509 - t476 * t514 + (qJD(2) * t471 + qJD(1)) * t516;
t517 = t480 * t475;
t518 = t476 * t479;
t456 = t471 * t517 + t518;
t469 = qJ(3) + pkin(11);
t467 = sin(t469);
t468 = cos(t469);
t470 = sin(pkin(6));
t515 = qJD(1) * t470;
t507 = t476 * t515;
t521 = t470 * t480;
t508 = t468 * t521;
t433 = t467 * (-qJD(3) * t456 + t507) - qJD(3) * t508 + t445 * t468;
t457 = t471 * t518 + t517;
t444 = qJD(1) * t457 + qJD(2) * t456;
t473 = sin(qJ(5));
t477 = cos(qJ(5));
t450 = -t456 * t468 + t467 * t521;
t455 = -t471 * t516 + t519;
t539 = t450 * t477 - t455 * t473;
t544 = t539 * qJD(5) - t433 * t473 + t444 * t477;
t540 = t450 * t473 + t455 * t477;
t543 = t540 * qJD(5) + t433 * t477 + t444 * t473;
t525 = t470 * t475;
t454 = t471 * t467 + t468 * t525;
t520 = t473 * t479;
t538 = -t454 * t477 + t470 * t520;
t478 = cos(qJ(3));
t466 = t478 * pkin(3) + pkin(2);
t534 = pkin(10) + r_i_i_C(2);
t537 = t468 * pkin(4) + t467 * t534 + t466;
t474 = sin(qJ(3));
t536 = -pkin(3) * t474 - pkin(4) * t467 + t468 * t534;
t512 = qJD(3) * t479;
t535 = (qJD(2) * t468 - qJD(5)) * t475 + t467 * t512;
t532 = r_i_i_C(3) + qJ(6);
t533 = r_i_i_C(1) + pkin(5);
t486 = t473 * t532 + t477 * t533 + pkin(4);
t526 = t468 * t473;
t524 = t470 * t476;
t523 = t470 * t478;
t522 = t470 * t479;
t513 = qJD(3) * t467;
t511 = qJD(5) * t468;
t506 = t480 * t515;
t505 = qJD(2) * t522;
t504 = t470 * t514;
t443 = qJD(1) * t456 + qJD(2) * t457;
t498 = t457 * t511 - t443;
t497 = t455 * t511 + t445;
t489 = t509 - t516;
t452 = t467 * t524 - t468 * t489;
t494 = t452 * t477 + t457 * t473;
t493 = -t452 * t473 + t457 * t477;
t492 = (qJD(2) - t511) * t479;
t491 = t467 * t489 + t468 * t524;
t490 = -t467 * t525 + t471 * t468;
t442 = qJD(1) * t455 + qJD(2) * t489;
t485 = -qJD(5) * t489 + t442 * t468 + t457 * t513;
t484 = qJD(5) * t456 - t444 * t468 + t455 * t513;
t483 = qJD(3) * t536;
t482 = qJD(6) * t473 + (-t473 * t533 + t477 * t532) * qJD(5);
t481 = qJD(3) * t450 - t445 * t467 + t468 * t507;
t472 = -qJ(4) - pkin(9);
t447 = qJD(3) * t490 + t468 * t505;
t436 = -t538 * qJD(5) + t447 * t473 - t477 * t504;
t431 = qJD(3) * t491 - t443 * t468 + t467 * t506;
t430 = qJD(3) * t452 - t443 * t467 - t468 * t506;
t421 = qJD(5) * t493 + t431 * t477 - t442 * t473;
t420 = qJD(5) * t494 + t431 * t473 + t442 * t477;
t1 = [t540 * qJD(6) - t433 * pkin(4) - t445 * t466 + t444 * t472 - t455 * qJD(4) + t534 * t481 - t533 * t543 + t532 * t544 + (-t480 * pkin(1) - pkin(8) * t524) * qJD(1) + (-t474 * t507 + (t456 * t474 + t478 * t521) * qJD(3)) * pkin(3) -(t457 * t526 - t477 * t489) * qJD(6) + t443 * t472 - t489 * qJD(4) + t533 * (t473 * t498 + t477 * t485) + t532 * (t473 * t485 - t477 * t498) + t537 * t442 - t457 * t483, t534 * t431 + (t478 * t506 + t443 * t474 + (-t474 * t524 + t478 * t489) * qJD(3)) * pkin(3) + t482 * t491 - t486 * t430, -t442, qJD(6) * t494 - t420 * t533 + t421 * t532, t420; -t493 * qJD(6) + t431 * pkin(4) - t443 * t466 + t442 * t472 + t457 * qJD(4) + t534 * t430 + t533 * t421 + t532 * t420 + (-t476 * pkin(1) + pkin(8) * t521) * qJD(1) + (t474 * t506 + (t474 * t489 + t476 * t523) * qJD(3)) * pkin(3) -(t455 * t526 + t456 * t477) * qJD(6) - t445 * t472 + t456 * qJD(4) + t533 * (t473 * t497 + t477 * t484) + t532 * (t473 * t484 - t477 * t497) - t537 * t444 - t455 * t483, t534 * t433 + (t478 * t507 - t445 * t474 + (-t456 * t478 + t474 * t521) * qJD(3)) * pkin(3) + t482 * (-t456 * t467 - t508) + t486 * t481, t444, -t539 * qJD(6) + t532 * t543 + t533 * t544, -t544; 0 (t533 * (t473 * t492 - t535 * t477) - t532 * (t535 * t473 + t477 * t492) - (-t468 * t520 + t475 * t477) * qJD(6) + t475 * qJD(4) + t536 * t512 + (-t479 * t472 - t475 * t537) * qJD(2)) * t470, t534 * t447 + (-t474 * t505 + (-t471 * t474 - t475 * t523) * qJD(3)) * pkin(3) + t482 * t490 + t486 * (-qJD(3) * t454 - t467 * t505) t504, -t538 * qJD(6) + t532 * (t473 * t504 + t447 * t477 + (-t454 * t473 - t477 * t522) * qJD(5)) - t533 * t436, t436;];
JaD_transl  = t1;
