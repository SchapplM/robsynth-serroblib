% Zeitableitung der translatorischen Teilmatrix der analytischen Jacobi-Matrix für Segment Nr. 6 (0=Basis) von
% S6RRPRRP14
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
%   pkin=[a2,a3,a4,a5,a6,alpha2,d1,d2,d4,d5]';
%
% Output:
% JaD_transl [3x6]
%   Zeitableitung der translatorischen Teilmatrix der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-02-26 21:53
% Revision: d75aae1ac561373cd3be920984c3df29a1c2ecc8 (2019-02-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JaD_transl = S6RRPRRP14_jacobiaD_transl_6_sym_varpar(qJ, qJD, r_i_i_C, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(3,1),zeros(10,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPRRP14_jacobiaD_transl_6_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRPRRP14_jacobiaD_transl_6_sym_varpar: qJD has to be [6x1] (double)');
assert(isa(r_i_i_C,'double') && isreal(r_i_i_C) && all(size(r_i_i_C) == [3 1]), ...
	'S6RRPRRP14_jacobiaD_transl_6_sym_varpar: Position vector r_i_i_C has to be [3x1] double');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RRPRRP14_jacobiaD_transl_6_sym_varpar: pkin has to be [10x1] (double)');

%% Symbolic Calculation
% From jacobiaD_transl_6_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 21:53:28
% EndTime: 2019-02-26 21:53:28
% DurationCPUTime: 0.76s
% Computational Cost: add. (827->124), mult. (2432->190), div. (0->0), fcn. (2503->10), ass. (0->76)
t468 = cos(pkin(6));
t471 = sin(qJ(2));
t476 = cos(qJ(1));
t514 = t476 * t471;
t472 = sin(qJ(1));
t475 = cos(qJ(2));
t515 = t472 * t475;
t456 = t468 * t514 + t515;
t457 = t468 * t515 + t514;
t444 = t457 * qJD(1) + t456 * qJD(2);
t470 = sin(qJ(4));
t474 = cos(qJ(4));
t513 = t476 * t475;
t516 = t472 * t471;
t455 = -t468 * t513 + t516;
t467 = sin(pkin(6));
t518 = t467 * t476;
t488 = t455 * t474 + t470 * t518;
t512 = qJD(1) * t472;
t505 = t467 * t512;
t431 = t488 * qJD(4) + t444 * t470 + t474 * t505;
t511 = qJD(2) * t471;
t500 = t472 * t511;
t503 = t471 * t512;
t445 = -t468 * t503 - t500 + (qJD(2) * t468 + qJD(1)) * t513;
t469 = sin(qJ(5));
t473 = cos(qJ(5));
t506 = t474 * t518;
t451 = -t455 * t470 + t506;
t532 = t451 * t473 - t456 * t469;
t537 = t532 * qJD(5) - t431 * t469 + t445 * t473;
t533 = t451 * t469 + t456 * t473;
t536 = t533 * qJD(5) + t431 * t473 + t445 * t469;
t531 = (qJD(4) * t455 + t505) * t470 - qJD(4) * t506 - t444 * t474;
t526 = r_i_i_C(2) + pkin(10);
t530 = -t470 * pkin(4) + t526 * t474 - qJ(3);
t525 = r_i_i_C(3) + qJ(6);
t527 = r_i_i_C(1) + pkin(5);
t481 = t525 * t469 + t527 * t473 + pkin(4);
t529 = pkin(3) + pkin(8);
t528 = pkin(9) + pkin(2);
t520 = t467 * t472;
t519 = t467 * t475;
t517 = t471 * t473;
t510 = qJD(4) * t474;
t509 = qJD(5) * t470;
t508 = qJD(6) * t469;
t507 = t473 * qJD(6);
t504 = qJD(1) * t518;
t502 = qJD(2) * t519;
t501 = t467 * t511;
t497 = qJD(2) + t509;
t442 = t468 * t500 + t503 + (-qJD(1) * t468 - qJD(2)) * t513;
t458 = -t468 * t516 + t513;
t496 = t458 * t509 - t442;
t495 = t456 * t509 + t444;
t449 = t457 * t470 + t474 * t520;
t494 = t449 * t473 + t458 * t469;
t493 = -t449 * t469 + t458 * t473;
t490 = (qJD(2) * t470 + qJD(5)) * t475;
t485 = -t468 * t474 + t470 * t519;
t489 = t467 * t471 * t469 - t473 * t485;
t487 = t457 * t474 - t470 * t520;
t486 = t468 * t470 + t474 * t519;
t443 = t456 * qJD(1) + t457 * qJD(2);
t480 = -qJD(5) * t457 - t443 * t470 + t458 * t510;
t479 = -qJD(5) * t455 + t445 * t470 + t456 * t510;
t478 = t508 + (-t527 * t469 + t525 * t473) * qJD(5);
t477 = t470 * t508 + qJD(3) + (pkin(4) * t474 + t526 * t470) * qJD(4);
t446 = t486 * qJD(4) - t470 * t501;
t435 = t489 * qJD(5) - t446 * t469 - t473 * t502;
t434 = t487 * qJD(4) - t442 * t470 + t474 * t504;
t433 = t449 * qJD(4) + t442 * t474 + t470 * t504;
t424 = t493 * qJD(5) + t434 * t473 - t443 * t469;
t423 = t494 * qJD(5) + t434 * t469 + t443 * t473;
t1 = [t533 * qJD(6) - t431 * pkin(4) - t444 * qJ(3) - t455 * qJD(3) - t528 * t445 - t526 * t531 - t527 * t536 + t525 * t537 + (-t476 * pkin(1) - t529 * t520) * qJD(1), t457 * t507 + t528 * t442 + t527 * (-t496 * t469 + t480 * t473) + t525 * (t480 * t469 + t496 * t473) + t477 * t458 + t530 * t443, -t442, -t481 * t433 + t526 * t434 + t478 * t487, t494 * qJD(6) - t527 * t423 + t525 * t424, t423; -t493 * qJD(6) + t434 * pkin(4) - t442 * qJ(3) + t457 * qJD(3) - t528 * t443 + t526 * t433 + t527 * t424 + t525 * t423 + (-t472 * pkin(1) + t529 * t518) * qJD(1), t455 * t507 - t528 * t444 + t527 * (-t495 * t469 + t479 * t473) + t525 * (t479 * t469 + t495 * t473) + t477 * t456 - t530 * t445, t444, t526 * t431 + t478 * t488 - t481 * t531, -t532 * qJD(6) + t525 * t536 + t527 * t537, -t537; 0 (t527 * (t473 * t490 + (-t497 * t469 + t473 * t510) * t471) + t525 * (t497 * t517 + (t471 * t510 + t490) * t469) - t475 * t507 + t477 * t471 + (-t528 * t471 - t475 * t530) * qJD(2)) * t467, t501, -t526 * t446 - t478 * t486 + t481 * (t485 * qJD(4) + t474 * t501) t489 * qJD(6) + t525 * (t469 * t502 - t446 * t473 + (t467 * t517 + t469 * t485) * qJD(5)) - t527 * t435, t435;];
JaD_transl  = t1;
