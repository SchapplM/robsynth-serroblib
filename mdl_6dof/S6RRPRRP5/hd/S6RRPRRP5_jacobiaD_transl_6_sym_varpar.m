% Zeitableitung der translatorischen Teilmatrix der analytischen Jacobi-Matrix für Segment Nr. 6 (0=Basis) von
% S6RRPRRP5
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
%   pkin=[a2,a3,a4,a5,a6,alpha2,d1,d2,d4,d5,theta3]';
%
% Output:
% JaD_transl [3x6]
%   Zeitableitung der translatorischen Teilmatrix der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-02-26 21:48
% Revision: d75aae1ac561373cd3be920984c3df29a1c2ecc8 (2019-02-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JaD_transl = S6RRPRRP5_jacobiaD_transl_6_sym_varpar(qJ, qJD, r_i_i_C, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(3,1),zeros(11,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPRRP5_jacobiaD_transl_6_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRPRRP5_jacobiaD_transl_6_sym_varpar: qJD has to be [6x1] (double)');
assert(isa(r_i_i_C,'double') && isreal(r_i_i_C) && all(size(r_i_i_C) == [3 1]), ...
	'S6RRPRRP5_jacobiaD_transl_6_sym_varpar: Position vector r_i_i_C has to be [3x1] double');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6RRPRRP5_jacobiaD_transl_6_sym_varpar: pkin has to be [11x1] (double)');

%% Symbolic Calculation
% From jacobiaD_transl_6_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 21:48:23
% EndTime: 2019-02-26 21:48:24
% DurationCPUTime: 0.85s
% Computational Cost: add. (1089->132), mult. (3144->219), div. (0->0), fcn. (3472->12), ass. (0->85)
t462 = sin(qJ(4));
t466 = cos(qJ(4));
t461 = sin(qJ(5));
t465 = cos(qJ(5));
t526 = pkin(5) + r_i_i_C(1);
t484 = t465 * r_i_i_C(2) + t526 * t461;
t477 = t484 * qJD(5);
t455 = pkin(5) * t465 + pkin(4);
t486 = r_i_i_C(1) * t465 - r_i_i_C(2) * t461 + t455;
t522 = r_i_i_C(3) + qJ(6) + pkin(10);
t468 = -(t486 * t462 - t522 * t466) * qJD(4) + t462 * qJD(6) - t466 * t477;
t459 = cos(pkin(6));
t457 = sin(pkin(11));
t521 = cos(pkin(11));
t525 = cos(qJ(2));
t495 = t525 * t521;
t463 = sin(qJ(2));
t511 = qJD(2) * t463;
t529 = -qJD(2) * t495 + t457 * t511;
t433 = t529 * t459;
t498 = t463 * t521;
t476 = t525 * t457 + t498;
t438 = t476 * t459;
t499 = qJD(2) * t525;
t440 = -qJD(2) * t498 - t457 * t499;
t464 = sin(qJ(1));
t467 = cos(qJ(1));
t475 = -t463 * t457 + t495;
t513 = qJD(1) * t464;
t409 = t438 * t513 - t464 * t440 + (-qJD(1) * t475 + t433) * t467;
t458 = sin(pkin(6));
t517 = t458 * t467;
t535 = qJD(4) * t517 + t409;
t488 = t438 * t467 + t464 * t475;
t503 = t458 * t513;
t534 = -qJD(4) * t488 + t503;
t415 = -t462 * t517 + t466 * t488;
t474 = t459 * t475;
t420 = -t464 * t476 + t467 * t474;
t533 = t415 * t461 + t420 * t465;
t532 = t415 * t465 - t420 * t461;
t530 = qJD(1) * t474 + qJD(2) * t475;
t403 = t534 * t462 - t535 * t466;
t527 = t522 * t462 + t486 * t466 + pkin(3);
t524 = pkin(2) * t459;
t518 = t458 * t464;
t515 = t463 * t464;
t514 = t463 * t467;
t512 = qJD(1) * t467;
t509 = qJD(5) * t461;
t508 = qJD(5) * t465;
t506 = pkin(2) * t511;
t505 = t525 * t464;
t504 = t525 * t467;
t502 = t458 * t512;
t434 = qJD(2) * t438;
t410 = -t467 * t434 - t530 * t464 - t476 * t512;
t494 = -t403 * t461 - t410 * t465;
t431 = t529 * t458;
t437 = t476 * t458;
t489 = -t437 * t462 + t459 * t466;
t413 = t489 * qJD(4) - t431 * t466;
t432 = qJD(2) * t437;
t493 = -t413 * t461 + t432 * t465;
t426 = t437 * t466 + t459 * t462;
t436 = t475 * t458;
t490 = -t426 * t465 + t436 * t461;
t487 = -t464 * t438 + t467 * t475;
t483 = t462 * t488 + t466 * t517;
t417 = -t462 * t487 + t466 * t518;
t418 = t462 * t518 + t466 * t487;
t402 = -t535 * t462 - t534 * t466;
t407 = -t464 * t434 + t530 * t467 - t476 * t513;
t423 = -t464 * t474 - t467 * t476;
t472 = t407 * t461 + (-t418 * t461 - t423 * t465) * qJD(5);
t469 = -t488 * qJD(1) + t464 * t433 + t467 * t440;
t401 = t417 * qJD(4) + t462 * t502 + t466 * t469;
t398 = -t401 * t461 + t407 * t465 + (-t418 * t465 + t423 * t461) * qJD(5);
t456 = t525 * pkin(2) + pkin(1);
t441 = -t458 * qJD(3) + t499 * t524;
t439 = t463 * t524 + (-pkin(8) - qJ(3)) * t458;
t412 = t426 * qJD(4) - t431 * t462;
t400 = t418 * qJD(4) + t462 * t469 - t466 * t502;
t399 = t401 * t465 + t472;
t1 = [-t483 * qJD(6) + t409 * pkin(3) + t464 * t506 - t467 * t441 - t486 * t403 + (pkin(9) + t484) * t410 - t522 * t402 + (t439 * t464 - t456 * t467) * qJD(1) + (t532 * r_i_i_C(2) + t526 * t533) * qJD(5) (t465 * t469 - t487 * t509) * r_i_i_C(2) + t469 * pkin(9) - t527 * t407 + ((t459 * t515 - t504) * qJD(2) + (-t459 * t504 + t515) * qJD(1)) * pkin(2) + t468 * t423 + t526 * (t461 * t469 + t487 * t508) t502, t418 * qJD(6) - t486 * t400 + t522 * t401 - t417 * t477, -t399 * r_i_i_C(2) + t526 * t398, t400; -t467 * t506 + t469 * pkin(3) + t407 * pkin(9) + t399 * r_i_i_C(1) + t398 * r_i_i_C(2) - t417 * qJD(6) + t401 * t455 - t464 * t441 + t522 * t400 + (-t439 * t467 - t456 * t464) * qJD(1) + t472 * pkin(5) (-t409 * t465 - t488 * t509) * r_i_i_C(2) - t409 * pkin(9) + t527 * t410 + ((-t459 * t514 - t505) * qJD(2) + (-t459 * t505 - t514) * qJD(1)) * pkin(2) + t468 * t420 + t526 * (-t409 * t461 + t488 * t508) t503, t415 * qJD(6) - t486 * t402 + t522 * t403 + t483 * t477, t494 * r_i_i_C(1) + (-t403 * t465 + t410 * t461) * r_i_i_C(2) + (-r_i_i_C(1) * t532 + t533 * r_i_i_C(2)) * qJD(5) + (-qJD(5) * t532 + t494) * pkin(5), t402; 0 (-t431 * t465 - t437 * t509) * r_i_i_C(2) - t431 * pkin(9) - t458 * t506 - t527 * t432 + t468 * t436 + t526 * (-t431 * t461 + t437 * t508) 0, t426 * qJD(6) - t486 * t412 + t522 * t413 - t489 * t477, t493 * r_i_i_C(1) + (-t413 * t465 - t432 * t461) * r_i_i_C(2) + (t490 * r_i_i_C(1) + (t426 * t461 + t436 * t465) * r_i_i_C(2)) * qJD(5) + (t490 * qJD(5) + t493) * pkin(5), t412;];
JaD_transl  = t1;
