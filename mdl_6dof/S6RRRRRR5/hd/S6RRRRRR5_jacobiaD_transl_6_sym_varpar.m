% Zeitableitung der translatorischen Teilmatrix der analytischen Jacobi-Matrix für Segment Nr. 6 (0=Basis) von
% S6RRRRRR5
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
%   pkin=[a2,a3,a4,a5,a6,alpha2,d1,d2,d3,d4,d5,d6]';
%
% Output:
% JaD_transl [3x6]
%   Zeitableitung der translatorischen Teilmatrix der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-02-26 22:49
% Revision: d75aae1ac561373cd3be920984c3df29a1c2ecc8 (2019-02-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JaD_transl = S6RRRRRR5_jacobiaD_transl_6_sym_varpar(qJ, qJD, r_i_i_C, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(3,1),zeros(12,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRRRR5_jacobiaD_transl_6_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRRRRR5_jacobiaD_transl_6_sym_varpar: qJD has to be [6x1] (double)');
assert(isa(r_i_i_C,'double') && isreal(r_i_i_C) && all(size(r_i_i_C) == [3 1]), ...
	'S6RRRRRR5_jacobiaD_transl_6_sym_varpar: Position vector r_i_i_C has to be [3x1] double');
assert(isreal(pkin) && all(size(pkin) == [12 1]), ...
  'S6RRRRRR5_jacobiaD_transl_6_sym_varpar: pkin has to be [12x1] (double)');

%% Symbolic Calculation
% From jacobiaD_transl_6_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 22:49:28
% EndTime: 2019-02-26 22:49:29
% DurationCPUTime: 0.76s
% Computational Cost: add. (1718->141), mult. (2202->217), div. (0->0), fcn. (2167->14), ass. (0->94)
t510 = pkin(12) + r_i_i_C(3);
t451 = sin(qJ(6));
t455 = cos(qJ(6));
t468 = t455 * r_i_i_C(1) - t451 * r_i_i_C(2);
t466 = pkin(5) + t468;
t482 = qJD(6) * t455;
t483 = qJD(6) * t451;
t518 = -r_i_i_C(1) * t483 - t482 * r_i_i_C(2);
t450 = cos(pkin(6));
t453 = sin(qJ(2));
t458 = cos(qJ(1));
t490 = t458 * t453;
t454 = sin(qJ(1));
t457 = cos(qJ(2));
t491 = t454 * t457;
t426 = t450 * t490 + t491;
t448 = qJ(3) + qJ(4);
t444 = qJ(5) + t448;
t439 = sin(t444);
t440 = cos(t444);
t449 = sin(pkin(6));
t493 = t449 * t458;
t414 = -t426 * t440 + t439 * t493;
t489 = t458 * t457;
t476 = t450 * t489;
t492 = t454 * t453;
t425 = -t476 + t492;
t517 = -t414 * t451 - t425 * t455;
t516 = t414 * t455 - t425 * t451;
t446 = qJD(3) + qJD(4);
t452 = sin(qJ(3));
t504 = pkin(3) * qJD(3);
t442 = sin(t448);
t509 = pkin(4) * t442;
t423 = -t446 * t509 - t452 * t504;
t441 = qJD(5) + t446;
t515 = -(t466 * t439 - t510 * t440) * t441 + t423;
t513 = t451 * r_i_i_C(1) + t455 * r_i_i_C(2);
t443 = cos(t448);
t456 = cos(qJ(3));
t432 = t456 * pkin(3) + pkin(4) * t443;
t430 = pkin(2) + t432;
t511 = t510 * t439 + t466 * t440 + t430;
t431 = t452 * pkin(3) + t509;
t505 = pkin(8) + t431;
t427 = t450 * t491 + t490;
t410 = t427 * qJD(1) + t426 * qJD(2);
t503 = t410 * t451;
t502 = t410 * t455;
t478 = t450 * t492;
t428 = -t478 + t489;
t499 = t428 * t440;
t498 = t439 * t441;
t497 = t443 * t446;
t496 = t449 * t453;
t495 = t449 * t454;
t494 = t449 * t457;
t488 = qJD(1) * t454;
t487 = qJD(1) * t458;
t486 = qJD(2) * t453;
t485 = qJD(2) * t457;
t484 = qJD(6) * t440;
t480 = t439 * t496;
t479 = t440 * t496;
t477 = t440 * t493;
t475 = t449 * t488;
t474 = t449 * t487;
t473 = t449 * t486;
t472 = t449 * t485;
t470 = qJD(2) * t450 + qJD(1);
t411 = -qJD(1) * t478 - t454 * t486 + t470 * t489;
t471 = -t411 * t440 + t441 * t477;
t409 = t426 * qJD(1) + t427 * qJD(2);
t467 = t441 * t495 - t409;
t465 = -t426 * t441 + t475;
t464 = t441 * t450 + t472;
t399 = t465 * t440 + (t441 * t493 - t411) * t439;
t400 = -t426 * t498 + t439 * t475 - t471;
t462 = t518 * (-t426 * t439 - t477) + t510 * t400 + t466 * t399;
t397 = t467 * t439 - t440 * t474 + t441 * t499;
t398 = -t428 * t498 + t439 * t474 + t467 * t440;
t461 = t518 * (-t428 * t439 + t440 * t495) + t510 * t398 - t466 * t397;
t407 = t464 * t440 - t441 * t480;
t460 = t518 * (t450 * t440 - t480) + t510 * t407 + t466 * (-t464 * t439 - t441 * t479);
t459 = t513 * t484 - t515;
t447 = -pkin(11) - pkin(10) - pkin(9);
t424 = pkin(4) * t497 + t456 * t504;
t420 = t450 * t439 + t479;
t416 = t439 * t495 + t499;
t408 = -qJD(1) * t476 - t458 * t485 + t470 * t492;
t402 = -t465 * t439 + t471;
t390 = t398 * t455 - t408 * t451 + (-t416 * t451 + t427 * t455) * qJD(6);
t389 = -t398 * t451 - t408 * t455 + (-t416 * t455 - t427 * t451) * qJD(6);
t1 = [(t402 * t455 - t503) * r_i_i_C(1) + (-t402 * t451 - t502) * r_i_i_C(2) + t402 * pkin(5) - t411 * t430 - t426 * t423 + t410 * t447 + t424 * t493 + t510 * t399 + (t517 * r_i_i_C(1) - t516 * r_i_i_C(2)) * qJD(6) + (-t458 * pkin(1) - t505 * t495) * qJD(1) (-t409 * t451 + t428 * t482) * r_i_i_C(1) + (-t409 * t455 - t428 * t483) * r_i_i_C(2) + t409 * t447 + t511 * t408 + t459 * t427, t409 * t431 - t428 * t424 + (t423 * t454 + t432 * t487) * t449 + t461 ((-t428 * t446 + t474) * t443 + (-t446 * t495 + t409) * t442) * pkin(4) + t461, t461, t389 * r_i_i_C(1) - t390 * r_i_i_C(2); t424 * t495 + t398 * pkin(5) + t390 * r_i_i_C(1) + t389 * r_i_i_C(2) + t408 * t447 - t409 * t430 + t428 * t423 + t510 * t397 + (-pkin(1) * t454 + t505 * t493) * qJD(1) (t411 * t451 + t426 * t482) * r_i_i_C(1) + (t411 * t455 - t426 * t483) * r_i_i_C(2) - t411 * t447 - t511 * t410 + t459 * t425, -t411 * t431 - t426 * t424 + (-t423 * t458 + t432 * t488) * t449 + t462 ((-t426 * t446 + t475) * t443 + (t446 * t493 - t411) * t442) * pkin(4) + t462, t462 (-t400 * t451 + t502) * r_i_i_C(1) + (-t400 * t455 - t503) * r_i_i_C(2) + (t516 * r_i_i_C(1) + t517 * r_i_i_C(2)) * qJD(6); 0 ((-qJD(2) * t511 + t468 * qJD(6)) * t453 + (-qJD(2) * t447 + t513 * (qJD(2) - t484) + t515) * t457) * t449, t450 * t423 + (-t424 * t453 - t431 * t485) * t449 + t460 (-t496 * t497 + (-t446 * t450 - t472) * t442) * pkin(4) + t460, t460 (-t407 * t451 + t455 * t473) * r_i_i_C(1) + (-t407 * t455 - t451 * t473) * r_i_i_C(2) + ((-t420 * t455 + t451 * t494) * r_i_i_C(1) + (t420 * t451 + t455 * t494) * r_i_i_C(2)) * qJD(6);];
JaD_transl  = t1;
