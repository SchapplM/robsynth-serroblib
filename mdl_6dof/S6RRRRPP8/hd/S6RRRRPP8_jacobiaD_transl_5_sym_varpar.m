% Zeitableitung der translatorischen Teilmatrix der analytischen Jacobi-Matrix für Segment Nr. 5 (0=Basis) von
% S6RRRRPP8
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
%   pkin=[a2,a3,a4,a5,a6,alpha2,d1,d2,d3,d4]';
%
% Output:
% JaD_transl [3x6]
%   Zeitableitung der translatorischen Teilmatrix der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-02-26 22:29
% Revision: d75aae1ac561373cd3be920984c3df29a1c2ecc8 (2019-02-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JaD_transl = S6RRRRPP8_jacobiaD_transl_5_sym_varpar(qJ, qJD, r_i_i_C, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(3,1),zeros(10,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRRPP8_jacobiaD_transl_5_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRRRPP8_jacobiaD_transl_5_sym_varpar: qJD has to be [6x1] (double)');
assert(isa(r_i_i_C,'double') && isreal(r_i_i_C) && all(size(r_i_i_C) == [3 1]), ...
	'S6RRRRPP8_jacobiaD_transl_5_sym_varpar: Position vector r_i_i_C has to be [3x1] double');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RRRRPP8_jacobiaD_transl_5_sym_varpar: pkin has to be [10x1] (double)');

%% Symbolic Calculation
% From jacobiaD_transl_5_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 22:29:38
% EndTime: 2019-02-26 22:29:39
% DurationCPUTime: 0.78s
% Computational Cost: add. (793->117), mult. (2339->192), div. (0->0), fcn. (2413->10), ass. (0->74)
t455 = cos(pkin(6));
t459 = sin(qJ(1));
t458 = sin(qJ(2));
t501 = t459 * t458;
t491 = t455 * t501;
t496 = qJD(2) * t458;
t462 = cos(qJ(2));
t463 = cos(qJ(1));
t498 = t463 * t462;
t433 = -qJD(1) * t491 - t459 * t496 + (qJD(2) * t455 + qJD(1)) * t498;
t499 = t463 * t458;
t500 = t459 * t462;
t444 = t455 * t499 + t500;
t457 = sin(qJ(3));
t461 = cos(qJ(3));
t454 = sin(pkin(6));
t497 = qJD(1) * t454;
t488 = t459 * t497;
t504 = t454 * t463;
t490 = t461 * t504;
t421 = (-qJD(3) * t444 + t488) * t457 - qJD(3) * t490 + t433 * t461;
t445 = t455 * t500 + t499;
t432 = t445 * qJD(1) + t444 * qJD(2);
t456 = sin(qJ(4));
t460 = cos(qJ(4));
t438 = -t444 * t461 + t457 * t504;
t443 = -t455 * t498 + t501;
t521 = t438 * t460 - t443 * t456;
t526 = t521 * qJD(4) - t421 * t456 + t432 * t460;
t522 = t438 * t456 + t443 * t460;
t525 = t522 * qJD(4) + t421 * t460 + t432 * t456;
t506 = t454 * t461;
t442 = t455 * t457 + t458 * t506;
t502 = t456 * t462;
t520 = -t442 * t460 + t454 * t502;
t515 = r_i_i_C(2) + pkin(10);
t519 = t461 * pkin(3) + t515 * t457 + pkin(2);
t494 = qJD(3) * t462;
t518 = (qJD(2) * t461 - qJD(4)) * t458 + t457 * t494;
t517 = -pkin(3) * t457 + t515 * t461;
t513 = r_i_i_C(3) + qJ(5);
t516 = r_i_i_C(1) + pkin(4);
t468 = t513 * t456 + t516 * t460 + pkin(3);
t507 = t454 * t459;
t505 = t454 * t462;
t503 = t456 * t461;
t495 = qJD(3) * t457;
t493 = qJD(4) * t461;
t487 = t463 * t497;
t486 = t454 * t496;
t485 = qJD(2) * t505;
t431 = t444 * qJD(1) + t445 * qJD(2);
t481 = t445 * t493 - t431;
t480 = t443 * t493 + t433;
t472 = t491 - t498;
t440 = t457 * t507 - t461 * t472;
t477 = t440 * t460 + t445 * t456;
t476 = -t440 * t456 + t445 * t460;
t475 = (qJD(2) - t493) * t462;
t474 = t457 * t472 + t459 * t506;
t473 = -t454 * t458 * t457 + t455 * t461;
t469 = qJD(3) * t517;
t430 = t443 * qJD(1) + t472 * qJD(2);
t467 = -qJD(4) * t472 + t430 * t461 + t445 * t495;
t466 = qJD(4) * t444 - t432 * t461 + t443 * t495;
t465 = qJD(5) * t456 + (-t516 * t456 + t513 * t460) * qJD(4);
t464 = t438 * qJD(3) - t433 * t457 + t461 * t488;
t435 = t473 * qJD(3) + t461 * t485;
t424 = -t520 * qJD(4) + t435 * t456 - t460 * t486;
t419 = t474 * qJD(3) - t431 * t461 + t457 * t487;
t418 = t440 * qJD(3) - t431 * t457 - t461 * t487;
t409 = t476 * qJD(4) + t419 * t460 - t430 * t456;
t408 = t477 * qJD(4) + t419 * t456 + t430 * t460;
t1 = [t522 * qJD(5) - t421 * pkin(3) - t433 * pkin(2) - t432 * pkin(9) + t515 * t464 - t516 * t525 + t513 * t526 + (-t463 * pkin(1) - pkin(8) * t507) * qJD(1) -(t445 * t503 - t460 * t472) * qJD(5) - t431 * pkin(9) + t516 * (t481 * t456 + t467 * t460) + t513 * (t467 * t456 - t481 * t460) - t445 * t469 + t519 * t430, -t468 * t418 + t515 * t419 + t465 * t474, t477 * qJD(5) - t516 * t408 + t513 * t409, t408, 0; -t476 * qJD(5) + t419 * pkin(3) - t431 * pkin(2) - t430 * pkin(9) + t515 * t418 + t516 * t409 + t513 * t408 + (-t459 * pkin(1) + pkin(8) * t504) * qJD(1) -(t443 * t503 + t444 * t460) * qJD(5) + t433 * pkin(9) + t516 * (t480 * t456 + t466 * t460) + t513 * (t466 * t456 - t480 * t460) - t443 * t469 - t519 * t432, t515 * t421 + t465 * (-t444 * t457 - t490) + t468 * t464, -t521 * qJD(5) + t513 * t525 + t516 * t526, -t526, 0; 0 (t516 * (t456 * t475 - t518 * t460) - t513 * (t518 * t456 + t460 * t475) - (t458 * t460 - t461 * t502) * qJD(5) + t517 * t494 + (t462 * pkin(9) - t458 * t519) * qJD(2)) * t454, t515 * t435 + t465 * t473 + t468 * (-t442 * qJD(3) - t457 * t485) -t520 * qJD(5) + t513 * (t456 * t486 + t435 * t460 + (-t442 * t456 - t460 * t505) * qJD(4)) - t516 * t424, t424, 0;];
JaD_transl  = t1;
