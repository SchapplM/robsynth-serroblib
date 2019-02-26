% Zeitableitung der translatorischen Teilmatrix der analytischen Jacobi-Matrix für Segment Nr. 5 (0=Basis) von
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

function JaD_transl = S6RRRPPP1_jacobiaD_transl_5_sym_varpar(qJ, qJD, r_i_i_C, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(3,1),zeros(10,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRPPP1_jacobiaD_transl_5_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRRPPP1_jacobiaD_transl_5_sym_varpar: qJD has to be [6x1] (double)');
assert(isa(r_i_i_C,'double') && isreal(r_i_i_C) && all(size(r_i_i_C) == [3 1]), ...
	'S6RRRPPP1_jacobiaD_transl_5_sym_varpar: Position vector r_i_i_C has to be [3x1] double');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RRRPPP1_jacobiaD_transl_5_sym_varpar: pkin has to be [10x1] (double)');

%% Symbolic Calculation
% From jacobiaD_transl_5_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 22:02:50
% EndTime: 2019-02-26 22:02:51
% DurationCPUTime: 0.93s
% Computational Cost: add. (484->131), mult. (1653->223), div. (0->0), fcn. (1578->10), ass. (0->85)
t407 = sin(qJ(2));
t406 = sin(qJ(3));
t409 = cos(qJ(3));
t403 = sin(pkin(6));
t475 = r_i_i_C(1) + qJ(4);
t438 = t403 * t475;
t427 = -t409 * pkin(3) - t406 * t438;
t418 = -pkin(2) + t427;
t410 = cos(qJ(2));
t405 = cos(pkin(6));
t484 = t475 * t405 + pkin(9);
t433 = t484 * t410;
t493 = t418 * t407 + t433;
t411 = cos(qJ(1));
t459 = t411 * t409;
t408 = sin(qJ(1));
t462 = t408 * t406;
t393 = t410 * t459 + t462;
t451 = qJD(3) * t411;
t439 = t409 * t451;
t444 = qJD(3) * t462;
t446 = t408 * qJD(2) * t407;
t389 = t393 * qJD(1) - t409 * t446 - t410 * t444 - t439;
t402 = sin(pkin(10));
t404 = cos(pkin(10));
t452 = qJD(3) * t409;
t443 = t408 * t452;
t456 = qJD(1) * t411;
t421 = t406 * t456 + t443;
t440 = t406 * t451;
t458 = qJD(1) * t408;
t422 = t409 * t458 + t440;
t388 = -t406 * t446 + t410 * t421 - t422;
t455 = qJD(2) * t410;
t445 = t408 * t455;
t424 = t407 * t456 + t445;
t482 = -t388 * t405 + t403 * t424;
t491 = -t389 * t402 + t482 * t404;
t490 = (t406 * t424 + t407 * t443) * t405 + (t410 * t456 - t446) * t403;
t454 = qJD(2) * t411;
t441 = t410 * t454;
t423 = t407 * t458 - t441;
t442 = t407 * t454;
t457 = qJD(1) * t410;
t425 = t408 * t457 + t442;
t489 = (t406 * t423 - t407 * t439) * t405 + t425 * t403;
t447 = t403 * t452;
t449 = t403 * qJD(4);
t450 = qJD(5) * t404;
t469 = t402 * t409;
t488 = t475 * (qJD(2) * t405 + t447) + (-qJD(3) * pkin(3) + t405 * t450 + t449) * t406 + qJD(2) * pkin(9) + qJD(5) * t469;
t428 = -t405 * qJD(4) + t403 * t450;
t485 = t428 * t407;
t386 = (-qJD(3) * t410 + qJD(1)) * t459 + (t442 + (-qJD(3) + t457) * t408) * t406;
t483 = -t386 * t405 + t403 * t423;
t478 = r_i_i_C(2) - pkin(4);
t477 = t407 * pkin(2);
t474 = r_i_i_C(3) + qJ(5);
t470 = t402 * t405;
t468 = t404 * t405;
t467 = t405 * t406;
t466 = t405 * t409;
t465 = t406 * t410;
t464 = t407 * t405;
t463 = t407 * t409;
t461 = t409 * t410;
t460 = t411 * t406;
t453 = qJD(3) * t407;
t448 = -pkin(2) * t410 - pkin(1);
t437 = t405 * qJ(4) + pkin(9);
t431 = -t402 * t406 + t404 * t466;
t430 = t403 * t410 + t406 * t464;
t420 = t431 * t410;
t419 = (t402 * t466 + t404 * t406) * t410;
t415 = t409 * t445 + (t409 * t456 - t444) * t407;
t414 = t407 * t422 - t409 * t441;
t413 = qJD(2) * t418 - t428;
t412 = -t488 * t407 + t413 * t410;
t392 = t408 * t409 - t410 * t460;
t391 = -t408 * t461 + t460;
t390 = t410 * t462 + t459;
t387 = t409 * t425 + t410 * t440 - t421;
t381 = -t386 * t403 - t405 * t423;
t369 = -t387 * t402 + t483 * t404;
t1 = [-(t390 * t468 - t391 * t402) * qJD(5) - t389 * pkin(3) - t478 * (-t389 * t404 - t402 * t482) + t474 * t491 + (-t390 * qJD(4) - t475 * t388) * t403 + (t485 + (-t433 + t477) * qJD(2)) * t408 + (-t408 * pkin(8) + (-t407 * t484 + t448) * t411) * qJD(1), -t478 * (-t489 * t402 + t414 * t404) + t474 * (t414 * t402 + t489 * t404) - t493 * t458 + t412 * t411 -(-t392 * t402 - t393 * t468) * qJD(5) + t386 * pkin(3) - t478 * (t386 * t404 + t387 * t470) + t474 * (t386 * t402 - t387 * t468) + (t393 * qJD(4) - t475 * t387) * t403, t381, t369, 0; t381 * r_i_i_C(1) - (t392 * t468 - t393 * t402) * qJD(5) - t387 * pkin(3) + (-t386 * qJ(4) - t392 * qJD(4)) * t403 - t478 * (-t387 * t404 - t483 * t402) + t474 * t369 + (-t485 + (t437 * t410 - t477) * qJD(2)) * t411 + (t411 * pkin(8) + (-t437 * t407 + t448) * t408) * qJD(1), t478 * (-t490 * t402 + t415 * t404) - t474 * (t415 * t402 + t490 * t404) + t493 * t456 + t412 * t408 -(t390 * t402 + t391 * t468) * qJD(5) - t388 * pkin(3) - t478 * (-t388 * t404 - t389 * t470) + t474 * (-t388 * t402 + t389 * t468) + (-t391 * qJD(4) + t475 * t389) * t403, t388 * t403 + t405 * t424, -t491, 0; 0, t478 * (qJD(3) * t419 + (-t430 * t402 + t404 * t463) * qJD(2)) - t474 * (-qJD(3) * t420 + (t402 * t463 + t430 * t404) * qJD(2)) + t413 * t407 + t488 * t410, t478 * ((-t402 * t467 + t404 * t409) * t453 + qJD(2) * t419) - t474 * ((t404 * t467 + t469) * t453 - qJD(2) * t420) + (-pkin(3) * t406 + t409 * t438) * t455 + (t427 * qJD(3) + t431 * qJD(5) + t409 * t449) * t407, t407 * t447 + (t403 * t465 + t464) * qJD(2), t431 * t453 + (t402 * t461 + (-t403 * t407 + t405 * t465) * t404) * qJD(2), 0;];
JaD_transl  = t1;
