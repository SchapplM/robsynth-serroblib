% Zeitableitung der translatorischen Teilmatrix der analytischen Jacobi-Matrix für Segment Nr. 6 (0=Basis) von
% S6RRRPPP1
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
%   pkin=[a2,a3,a4,a5,a6,alpha4,d1,d2,d3,theta4]';
%
% Output:
% JaD_transl [3x6]
%   Zeitableitung der translatorischen Teilmatrix der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-02-26 22:02
% Revision: d75aae1ac561373cd3be920984c3df29a1c2ecc8 (2019-02-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JaD_transl = S6RRRPPP1_jacobiaD_transl_6_sym_varpar(qJ, qJD, r_i_i_C, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(3,1),zeros(10,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRPPP1_jacobiaD_transl_6_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRRPPP1_jacobiaD_transl_6_sym_varpar: qJD has to be [6x1] (double)');
assert(isa(r_i_i_C,'double') && isreal(r_i_i_C) && all(size(r_i_i_C) == [3 1]), ...
	'S6RRRPPP1_jacobiaD_transl_6_sym_varpar: Position vector r_i_i_C has to be [3x1] double');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RRRPPP1_jacobiaD_transl_6_sym_varpar: pkin has to be [10x1] (double)');

%% Symbolic Calculation
% From jacobiaD_transl_6_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 22:02:50
% EndTime: 2019-02-26 22:02:51
% DurationCPUTime: 1.05s
% Computational Cost: add. (647->159), mult. (2194->269), div. (0->0), fcn. (2123->10), ass. (0->104)
t429 = sin(qJ(2));
t425 = sin(pkin(6));
t428 = sin(qJ(3));
t508 = t425 * t428;
t431 = cos(qJ(3));
t519 = t431 * pkin(3);
t449 = qJ(4) * t508 + pkin(2) + t519;
t432 = cos(qJ(2));
t427 = cos(pkin(6));
t464 = qJ(4) * t427 + pkin(9);
t529 = t464 * t432;
t537 = -t449 * t429 + t529;
t433 = cos(qJ(1));
t488 = qJD(2) * t433;
t472 = t432 * t488;
t430 = sin(qJ(1));
t494 = qJD(1) * t430;
t445 = t429 * t494 - t472;
t485 = qJD(3) * t433;
t465 = t431 * t485;
t436 = t445 * t428 - t429 * t465;
t467 = t429 * t488;
t493 = qJD(1) * t432;
t447 = t430 * t493 + t467;
t536 = t447 * t425 + t436 * t427;
t486 = qJD(3) * t431;
t468 = t430 * t486;
t492 = qJD(1) * t433;
t443 = t428 * t492 + t468;
t466 = t428 * t485;
t444 = t431 * t494 + t466;
t490 = qJD(2) * t430;
t475 = t429 * t490;
t399 = -t428 * t475 + t443 * t432 - t444;
t495 = t432 * t433;
t499 = t430 * t428;
t408 = t431 * t495 + t499;
t469 = qJD(3) * t499;
t500 = t429 * t431;
t474 = qJD(2) * t500;
t400 = t408 * qJD(1) - t430 * t474 - t432 * t469 - t465;
t424 = sin(pkin(10));
t426 = cos(pkin(10));
t489 = qJD(2) * t432;
t473 = t430 * t489;
t478 = t429 * t492;
t446 = t473 + t478;
t535 = -t400 * t424 + t426 * (-t399 * t427 + t446 * t425);
t455 = qJD(5) * t426 - qJD(6) * t424;
t471 = t425 * t486;
t484 = t425 * qJD(4);
t491 = qJD(2) * t427;
t534 = (-qJD(3) * pkin(3) + t455 * t427 + t484) * t428 + (t471 + t491) * qJ(4) + (qJD(5) * t424 + qJD(6) * t426) * t431 + qJD(2) * pkin(9);
t438 = t427 * qJD(4) - t455 * t425;
t533 = qJD(2) * t529 + t429 * (-qJD(2) * pkin(2) + t438);
t527 = (t431 * t492 - t469) * t429;
t496 = t431 * t433;
t397 = (-qJD(3) * t432 + qJD(1)) * t496 + (t467 + (-qJD(3) + t493) * t430) * t428;
t526 = -t397 * t427 + t445 * t425;
t524 = t399 * t425 + t446 * t427;
t460 = t429 * t468;
t522 = t446 * t428 + t460;
t520 = r_i_i_C(1) + pkin(5);
t518 = r_i_i_C(2) + qJ(5);
t511 = t424 * t425;
t510 = t424 * t427;
t509 = t424 * t428;
t507 = t425 * t429;
t506 = t426 * t427;
t504 = t427 * t428;
t503 = t427 * t431;
t502 = t428 * t429;
t501 = t428 * t432;
t498 = t430 * t432;
t497 = t431 * t432;
t487 = qJD(3) * t429;
t482 = qJ(4) + t520;
t481 = r_i_i_C(3) + qJ(6) + pkin(4);
t480 = t427 * t502;
t479 = t426 * t497;
t477 = t432 * t492;
t470 = t429 * t486;
t463 = t424 * t478;
t462 = t491 * t509;
t461 = t489 * t511;
t453 = t426 * t503 - t509;
t452 = t424 * t503 + t426 * t428;
t451 = t427 * t501 - t507;
t448 = -t475 + t477;
t442 = t453 * t432;
t441 = t452 * t432;
t440 = t399 * t510 - t400 * t426 - t425 * t463 - t430 * t461;
t439 = -pkin(2) * t432 - t464 * t429 - pkin(1);
t437 = t444 * t429 - t431 * t472;
t435 = -t449 * qJD(2) + t438;
t434 = -t534 * t429 + t435 * t432;
t407 = -t428 * t495 + t430 * t431;
t406 = t428 * t433 - t430 * t497;
t405 = t428 * t498 + t496;
t398 = t447 * t431 + t432 * t466 - t443;
t389 = -t397 * t425 - t445 * t427;
t378 = -t398 * t426 - t424 * t526;
t377 = -t398 * t424 + t426 * t526;
t1 = [-(-t405 * t510 - t406 * t426) * qJD(6) - (t405 * t506 - t406 * t424) * qJD(5) - t400 * pkin(3) + (-qJ(4) * t399 - qJD(4) * t405) * t425 - t520 * t524 + t518 * t535 + t481 * t440 - t533 * t430 + (-t430 * pkin(8) + t439 * t433) * qJD(1), t520 * (t436 * t425 - t447 * t427) + t518 * (t437 * t424 + t536 * t426) + t481 * (-t536 * t424 + t437 * t426) - t537 * t494 + t434 * t433 -(-t407 * t426 + t408 * t510) * qJD(6) - (-t407 * t424 - t408 * t506) * qJD(5) + t397 * pkin(3) + t518 * (t397 * t424 - t398 * t506) + t481 * (t397 * t426 + t398 * t510) + (t408 * qJD(4) - t482 * t398) * t425, t389, t377, t378; -(-t407 * t510 - t408 * t426) * qJD(6) - (t407 * t506 - t408 * t424) * qJD(5) - t398 * pkin(3) + (-qJ(4) * t397 - qJD(4) * t407) * t425 + t520 * t389 + t518 * t377 + t481 * t378 + t533 * t433 + (t433 * pkin(8) + t439 * t430) * qJD(1), t520 * (-t425 * t522 + t448 * t427) - t518 * ((t431 * t473 + t527) * t424 + (t448 * t425 + t427 * t522) * t426) - t481 * (-t460 * t510 - t462 * t498 - t463 * t504 - t477 * t511 + (t424 * t507 + t479) * t490 + t426 * t527) + t537 * t492 + t434 * t430 -(t405 * t426 - t406 * t510) * qJD(6) - (t405 * t424 + t406 * t506) * qJD(5) - t399 * pkin(3) + t518 * (-t399 * t424 + t400 * t506) + t481 * (-t399 * t426 - t400 * t510) + (-t406 * qJD(4) + t482 * t400) * t425, t524, -t535, -t440; 0, t520 * (t432 * t471 + (-t425 * t502 + t427 * t432) * qJD(2)) - t518 * (-qJD(3) * t442 + (t424 * t500 + (t425 * t432 + t480) * t426) * qJD(2)) - t481 * (qJD(3) * t441 + t426 * t474 - t429 * t462 - t461) + t435 * t429 + t534 * t432, -t518 * ((t424 * t431 + t426 * t504) * t487 - qJD(2) * t442) - t481 * (-qJD(3) * t424 * t480 + qJD(2) * t441 + t426 * t470) + (t482 * t431 * t425 - pkin(3) * t428) * t489 + (-t452 * qJD(6) + t453 * qJD(5) + t431 * t484 + (-t482 * t508 - t519) * qJD(3)) * t429, t425 * t470 + (t425 * t501 + t427 * t429) * qJD(2), t453 * t487 + (t424 * t497 + t451 * t426) * qJD(2), -t452 * t487 + (-t451 * t424 + t479) * qJD(2);];
JaD_transl  = t1;
