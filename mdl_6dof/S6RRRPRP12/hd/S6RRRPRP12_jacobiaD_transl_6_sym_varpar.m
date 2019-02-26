% Zeitableitung der translatorischen Teilmatrix der analytischen Jacobi-Matrix für Segment Nr. 6 (0=Basis) von
% S6RRRPRP12
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
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d1,d2,d3,d5]';
%
% Output:
% JaD_transl [3x6]
%   Zeitableitung der translatorischen Teilmatrix der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-02-26 22:15
% Revision: d75aae1ac561373cd3be920984c3df29a1c2ecc8 (2019-02-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JaD_transl = S6RRRPRP12_jacobiaD_transl_6_sym_varpar(qJ, qJD, r_i_i_C, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(3,1),zeros(10,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRPRP12_jacobiaD_transl_6_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRRPRP12_jacobiaD_transl_6_sym_varpar: qJD has to be [6x1] (double)');
assert(isa(r_i_i_C,'double') && isreal(r_i_i_C) && all(size(r_i_i_C) == [3 1]), ...
	'S6RRRPRP12_jacobiaD_transl_6_sym_varpar: Position vector r_i_i_C has to be [3x1] double');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RRRPRP12_jacobiaD_transl_6_sym_varpar: pkin has to be [10x1] (double)');

%% Symbolic Calculation
% From jacobiaD_transl_6_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 22:15:34
% EndTime: 2019-02-26 22:15:35
% DurationCPUTime: 0.72s
% Computational Cost: add. (894->122), mult. (2628->192), div. (0->0), fcn. (2705->10), ass. (0->77)
t466 = sin(qJ(1));
t523 = cos(pkin(6));
t482 = qJD(2) * t523 + qJD(1);
t465 = sin(qJ(2));
t494 = t466 * t523;
t490 = t465 * t494;
t509 = qJD(2) * t465;
t469 = cos(qJ(2));
t470 = cos(qJ(1));
t510 = t470 * t469;
t436 = -qJD(1) * t490 - t466 * t509 + t482 * t510;
t493 = t470 * t523;
t449 = t465 * t493 + t466 * t469;
t464 = sin(qJ(3));
t462 = sin(pkin(6));
t468 = cos(qJ(3));
t515 = t462 * t468;
t500 = qJD(1) * t515;
t512 = t464 * t470;
t501 = t462 * t512;
t507 = qJD(3) * t468;
t425 = -qJD(3) * t501 + t436 * t464 + t449 * t507 - t466 * t500;
t450 = t470 * t465 + t469 * t494;
t435 = qJD(1) * t450 + qJD(2) * t449;
t463 = sin(qJ(5));
t467 = cos(qJ(5));
t489 = t469 * t493;
t511 = t466 * t465;
t448 = -t489 + t511;
t513 = t462 * t470;
t531 = t449 * t464 + t468 * t513;
t532 = -t448 * t463 + t467 * t531;
t537 = t532 * qJD(5) + t425 * t463 + t435 * t467;
t486 = t448 * t467 + t463 * t531;
t536 = -qJD(5) * t486 + t425 * t467 - t435 * t463;
t491 = -t467 * qJD(6) + qJD(4);
t504 = r_i_i_C(2) + pkin(10) + pkin(3);
t471 = -qJ(4) * t507 + (qJD(3) * t504 - t491) * t464;
t530 = (qJD(2) * t464 + qJD(5)) * t465 - t469 * t507;
t529 = t464 * qJ(4) + t468 * t504 + pkin(2);
t516 = t462 * t466;
t502 = t464 * t516;
t528 = -qJD(1) * t502 + qJD(3) * t531 - t436 * t468;
t526 = pkin(4) + pkin(9);
t525 = -r_i_i_C(1) - pkin(5);
t524 = r_i_i_C(3) + qJ(6);
t451 = -t490 + t510;
t517 = t451 * t464;
t514 = t462 * t469;
t508 = qJD(2) * t469;
t506 = qJD(5) * t464;
t505 = t463 * qJD(6);
t503 = t467 * t514;
t499 = t462 * t509;
t498 = t462 * t508;
t434 = qJD(1) * t449 + qJD(2) * t450;
t488 = -t450 * t506 - t434;
t487 = -t448 * t506 + t436;
t443 = -t466 * t515 + t517;
t484 = t443 * t467 - t450 * t463;
t483 = t443 * t463 + t450 * t467;
t481 = (qJD(2) + t506) * t469;
t479 = t451 * t468 + t502;
t446 = t462 * t465 * t464 - t468 * t523;
t477 = t464 * t523 + t465 * t515;
t475 = -t463 * t525 - t467 * t524 + qJ(4);
t433 = -qJD(1) * t489 - t470 * t508 + t482 * t511;
t474 = qJD(5) * t451 - t433 * t464 + t450 * t507;
t473 = qJD(5) * t449 + t435 * t464 + t448 * t507;
t472 = (t463 * t524 - t467 * t525) * qJD(5) + t491;
t438 = qJD(3) * t477 + t464 * t498;
t429 = -t438 * t467 - qJD(5) * t503 + (qJD(5) * t446 + t499) * t463;
t424 = -t434 * t468 - qJD(3) * t517 + (qJD(1) * t512 + t466 * t507) * t462;
t423 = qJD(3) * t479 - t434 * t464 - t470 * t500;
t412 = qJD(5) * t484 + t423 * t463 - t433 * t467;
t411 = qJD(5) * t483 - t423 * t467 - t433 * t463;
t1 = [t532 * qJD(6) - t425 * qJ(4) - t531 * qJD(4) - t436 * pkin(2) - t526 * t435 + t525 * t537 + t524 * t536 + (-t470 * pkin(1) - pkin(8) * t516) * qJD(1) + t504 * t528, t451 * t505 - t526 * t434 - t525 * (-t463 * t474 + t467 * t488) + t524 * (t463 * t488 + t467 * t474) + t471 * t450 + t529 * t433, -t423 * t504 + t424 * t475 + t472 * t479, t423, qJD(6) * t483 + t411 * t525 + t412 * t524, t411; -t484 * qJD(6) + t423 * qJ(4) + t443 * qJD(4) - t434 * pkin(2) - t526 * t433 - t525 * t412 + t524 * t411 + (-t466 * pkin(1) + pkin(8) * t513) * qJD(1) + t504 * t424, t449 * t505 + t526 * t436 - t525 * (-t463 * t473 + t467 * t487) + t524 * (t463 * t487 + t467 * t473) + t471 * t448 - t529 * t435, -t504 * t425 - t475 * t528 + t472 * (t449 * t468 - t501) t425, t486 * qJD(6) + t524 * t537 - t525 * t536, -t536; 0 (-t525 * (-t530 * t463 + t467 * t481) + t524 * (t463 * t481 + t530 * t467) + t465 * t505 - t471 * t469 + (-t465 * t529 + t469 * t526) * qJD(2)) * t462, -t504 * t438 + t475 * (-qJD(3) * t446 + t468 * t498) + t472 * t477, t438 -(-t446 * t463 + t503) * qJD(6) + t524 * (t467 * t499 + t438 * t463 + (t446 * t467 + t463 * t514) * qJD(5)) + t525 * t429, t429;];
JaD_transl  = t1;
