% Zeitableitung der translatorischen Teilmatrix der analytischen Jacobi-Matrix für Segment Nr. 6 (0=Basis) von
% S6RRRPPR9
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
%   pkin=[a2,a3,a4,a5,a6,alpha2,d1,d2,d3,d6,theta4]';
%
% Output:
% JaD_transl [3x6]
%   Zeitableitung der translatorischen Teilmatrix der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-02-26 22:08
% Revision: d75aae1ac561373cd3be920984c3df29a1c2ecc8 (2019-02-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JaD_transl = S6RRRPPR9_jacobiaD_transl_6_sym_varpar(qJ, qJD, r_i_i_C, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(3,1),zeros(11,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRPPR9_jacobiaD_transl_6_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRRPPR9_jacobiaD_transl_6_sym_varpar: qJD has to be [6x1] (double)');
assert(isa(r_i_i_C,'double') && isreal(r_i_i_C) && all(size(r_i_i_C) == [3 1]), ...
	'S6RRRPPR9_jacobiaD_transl_6_sym_varpar: Position vector r_i_i_C has to be [3x1] double');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6RRRPPR9_jacobiaD_transl_6_sym_varpar: pkin has to be [11x1] (double)');

%% Symbolic Calculation
% From jacobiaD_transl_6_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 22:08:12
% EndTime: 2019-02-26 22:08:13
% DurationCPUTime: 1.08s
% Computational Cost: add. (1049->158), mult. (3138->268), div. (0->0), fcn. (3263->12), ass. (0->94)
t471 = cos(pkin(6));
t474 = sin(qJ(2));
t523 = cos(qJ(1));
t499 = t523 * t474;
t475 = sin(qJ(1));
t478 = cos(qJ(2));
t510 = t475 * t478;
t454 = t471 * t499 + t510;
t473 = sin(qJ(3));
t477 = cos(qJ(3));
t469 = sin(pkin(6));
t500 = t469 * t523;
t443 = t454 * t477 - t473 * t500;
t498 = t523 * t478;
t489 = t471 * t498;
t511 = t475 * t474;
t453 = -t489 + t511;
t468 = sin(pkin(11));
t470 = cos(pkin(11));
t420 = t443 * t468 - t453 * t470;
t421 = t443 * t470 + t453 * t468;
t472 = sin(qJ(6));
t476 = cos(qJ(6));
t530 = -t420 * t476 + t421 * t472;
t529 = t420 * t472 + t421 * t476;
t492 = t523 * qJD(2);
t493 = t523 * qJD(1);
t502 = t471 * t511;
t508 = qJD(2) * t474;
t437 = -qJD(1) * t502 - t475 * t508 + (t471 * t492 + t493) * t478;
t528 = -qJD(3) * t500 + t437;
t504 = -r_i_i_C(3) - pkin(10) + qJ(4);
t527 = qJD(3) * (-pkin(3) * t473 + t504 * t477) + t473 * qJD(4);
t526 = t477 * pkin(3) + t504 * t473 + pkin(2);
t524 = pkin(4) + pkin(5);
t455 = t471 * t510 + t499;
t436 = t455 * qJD(1) + t454 * qJD(2);
t521 = t436 * t468;
t520 = t436 * t470;
t517 = t468 * t477;
t516 = t469 * t475;
t515 = t469 * t478;
t514 = t470 * t477;
t513 = t470 * t478;
t512 = t474 * t477;
t509 = t477 * t478;
t507 = qJD(3) * t473;
t506 = qJD(3) * t477;
t503 = t469 * t474 * t473;
t497 = qJD(1) * t516;
t496 = t469 * t508;
t495 = qJD(2) * t515;
t494 = t478 * t507;
t491 = t528 * t477;
t490 = t477 * t500;
t487 = t472 * r_i_i_C(1) + t476 * r_i_i_C(2) + qJ(5);
t486 = t476 * r_i_i_C(1) - t472 * r_i_i_C(2) + t524;
t456 = t498 - t502;
t446 = t456 * t477 + t473 * t516;
t452 = t469 * t512 + t471 * t473;
t434 = -qJD(1) * t489 - t478 * t492 + (qJD(2) * t471 + qJD(1)) * t511;
t485 = t434 * t477 + t455 * t507;
t484 = -t436 * t477 + t453 * t507;
t483 = t454 * t473 + t490;
t416 = t454 * t506 + t528 * t473 - t477 * t497;
t480 = -t487 * t468 - t486 * t470 - pkin(3);
t479 = t468 * qJD(5) + ((t468 * t476 - t470 * t472) * r_i_i_C(1) + (-t468 * t472 - t470 * t476) * r_i_i_C(2)) * qJD(6);
t448 = (t468 * t474 + t470 * t509) * t469;
t447 = (t468 * t509 - t470 * t474) * t469;
t445 = -t456 * t473 + t477 * t516;
t441 = -qJD(3) * t503 + (qJD(3) * t471 + t495) * t477;
t440 = t452 * qJD(3) + t473 * t495;
t439 = t452 * t470 - t468 * t515;
t438 = t452 * t468 + t469 * t513;
t435 = t454 * qJD(1) + t455 * qJD(2);
t431 = -t455 * t514 + t456 * t468;
t430 = -t455 * t517 - t456 * t470;
t429 = -t453 * t514 + t454 * t468;
t428 = -t453 * t517 - t454 * t470;
t427 = t441 * t470 + t468 * t496;
t426 = t441 * t468 - t470 * t496;
t425 = t446 * t470 + t455 * t468;
t424 = t446 * t468 - t455 * t470;
t419 = (qJD(3) * t454 - t497) * t473 - t491;
t417 = -t454 * t507 + t473 * t497 + t491;
t415 = -t435 * t477 - t456 * t507 + (t473 * t493 + t475 * t506) * t469;
t414 = -qJD(1) * t490 + t446 * qJD(3) - t435 * t473;
t407 = t417 * t470 + t521;
t406 = t417 * t468 - t520;
t405 = t415 * t470 - t434 * t468;
t404 = t415 * t468 + t434 * t470;
t403 = t404 * t472 + t405 * t476 + (t424 * t476 - t425 * t472) * qJD(6);
t402 = t404 * t476 - t405 * t472 + (-t424 * t472 - t425 * t476) * qJD(6);
t1 = [-t420 * qJD(5) + t419 * pkin(3) - t483 * qJD(4) - t437 * pkin(2) - t436 * pkin(9) + t487 * (t419 * t468 + t520) + t486 * (t419 * t470 - t521) + (t530 * r_i_i_C(1) + t529 * r_i_i_C(2)) * qJD(6) + (-t523 * pkin(1) - pkin(8) * t516) * qJD(1) - t504 * t416, -t435 * pkin(9) + t430 * qJD(5) + t487 * (t435 * t470 + t485 * t468) + t486 * (-t435 * t468 + t485 * t470) + ((t430 * t476 - t431 * t472) * r_i_i_C(1) + (-t430 * t472 - t431 * t476) * r_i_i_C(2)) * qJD(6) - t527 * t455 + t526 * t434, t446 * qJD(4) + t480 * t414 + t504 * t415 + t479 * t445, t414, t404, r_i_i_C(1) * t402 - r_i_i_C(2) * t403; -t435 * pkin(2) + t415 * pkin(3) - t434 * pkin(9) + t403 * r_i_i_C(1) + t402 * r_i_i_C(2) + t404 * qJ(5) - t445 * qJD(4) + t424 * qJD(5) + t524 * t405 + (-pkin(1) * t475 + pkin(8) * t500) * qJD(1) + t504 * t414, t437 * pkin(9) + t428 * qJD(5) + t487 * (-t437 * t470 + t484 * t468) + t486 * (t437 * t468 + t484 * t470) + ((t428 * t476 - t429 * t472) * r_i_i_C(1) + (-t428 * t472 - t429 * t476) * r_i_i_C(2)) * qJD(6) - t527 * t453 - t526 * t436, t443 * qJD(4) + t480 * t416 + t504 * t417 - t479 * t483, t416, t406 (t406 * t476 - t407 * t472) * r_i_i_C(1) + (-t406 * t472 - t407 * t476) * r_i_i_C(2) + (-t529 * r_i_i_C(1) + t530 * r_i_i_C(2)) * qJD(6); 0, t447 * qJD(5) + ((t447 * t476 - t448 * t472) * r_i_i_C(1) + (-t447 * t472 - t448 * t476) * r_i_i_C(2)) * qJD(6) + (-t487 * (t468 * t494 + (t468 * t512 + t513) * qJD(2)) + t486 * (-t470 * t494 + (t468 * t478 - t470 * t512) * qJD(2)) + t527 * t478 + (pkin(9) * t478 - t474 * t526) * qJD(2)) * t469, t452 * qJD(4) + t504 * t441 + t479 * (t471 * t477 - t503) + t480 * t440, t440, t426 (t426 * t476 - t427 * t472) * r_i_i_C(1) + (-t426 * t472 - t427 * t476) * r_i_i_C(2) + ((-t438 * t472 - t439 * t476) * r_i_i_C(1) + (-t438 * t476 + t439 * t472) * r_i_i_C(2)) * qJD(6);];
JaD_transl  = t1;
