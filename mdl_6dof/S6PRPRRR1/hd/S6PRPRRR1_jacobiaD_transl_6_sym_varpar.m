% Zeitableitung der translatorischen Teilmatrix der analytischen Jacobi-Matrix für Segment Nr. 6 (0=Basis) von
% S6PRPRRR1
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
%   pkin=[a2,a3,a4,a5,a6,alpha2,d2,d4,d5,d6,theta1,theta3]';
%
% Output:
% JaD_transl [3x6]
%   Zeitableitung der translatorischen Teilmatrix der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-02-26 19:53
% Revision: d75aae1ac561373cd3be920984c3df29a1c2ecc8 (2019-02-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JaD_transl = S6PRPRRR1_jacobiaD_transl_6_sym_varpar(qJ, qJD, r_i_i_C, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(3,1),zeros(12,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRPRRR1_jacobiaD_transl_6_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6PRPRRR1_jacobiaD_transl_6_sym_varpar: qJD has to be [6x1] (double)');
assert(isa(r_i_i_C,'double') && isreal(r_i_i_C) && all(size(r_i_i_C) == [3 1]), ...
	'S6PRPRRR1_jacobiaD_transl_6_sym_varpar: Position vector r_i_i_C has to be [3x1] double');
assert(isreal(pkin) && all(size(pkin) == [12 1]), ...
  'S6PRPRRR1_jacobiaD_transl_6_sym_varpar: pkin has to be [12x1] (double)');

%% Symbolic Calculation
% From jacobiaD_transl_6_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 19:53:45
% EndTime: 2019-02-26 19:53:46
% DurationCPUTime: 0.62s
% Computational Cost: add. (873->99), mult. (1890->183), div. (0->0), fcn. (2030->14), ass. (0->67)
t490 = pkin(10) + r_i_i_C(3);
t449 = cos(qJ(6));
t475 = qJD(6) * t449;
t446 = sin(qJ(6));
t476 = qJD(6) * t446;
t496 = -r_i_i_C(1) * t476 - t475 * r_i_i_C(2);
t440 = qJ(4) + qJ(5);
t437 = sin(t440);
t438 = cos(t440);
t439 = qJD(4) + qJD(5);
t447 = sin(qJ(4));
t462 = r_i_i_C(1) * t449 - r_i_i_C(2) * t446 + pkin(5);
t495 = (t462 * t437 - t490 * t438) * t439 + (t446 * r_i_i_C(1) + t449 * r_i_i_C(2)) * t438 * qJD(6) + qJD(4) * t447 * pkin(4);
t445 = cos(pkin(6));
t441 = sin(pkin(12));
t448 = sin(qJ(2));
t486 = cos(pkin(12));
t489 = cos(qJ(2));
t460 = t489 * t441 + t448 * t486;
t424 = t460 * t445;
t467 = t489 * t486;
t478 = qJD(2) * t448;
t493 = -qJD(2) * t467 + t441 * t478;
t421 = t493 * t445;
t426 = t460 * qJD(2);
t442 = sin(pkin(11));
t444 = cos(pkin(11));
t403 = t421 * t444 + t426 * t442;
t443 = sin(pkin(6));
t482 = t443 * t444;
t494 = t439 * t482 + t403;
t459 = -t448 * t441 + t467;
t450 = cos(qJ(4));
t456 = pkin(4) * t450 + t490 * t437 + t462 * t438 + pkin(3);
t487 = pkin(2) * qJD(2);
t485 = t437 * t439;
t484 = t438 * t439;
t483 = t442 * t443;
t481 = t443 * t447;
t480 = t445 * t448;
t419 = t493 * t443;
t468 = t439 * t445 - t419;
t405 = t442 * t421 - t426 * t444;
t465 = t439 * t483 + t405;
t464 = t424 * t444 + t442 * t459;
t463 = -t424 * t442 + t444 * t459;
t458 = t445 * t459;
t423 = t460 * t443;
t457 = qJD(2) * t424;
t387 = -t494 * t438 - t464 * t485;
t455 = t496 * (-t437 * t464 - t438 * t482) + t490 * t387 + t462 * (t494 * t437 - t464 * t484);
t389 = t465 * t438 - t463 * t485;
t454 = t496 * (-t437 * t463 + t438 * t483) + t490 * t389 + t462 * (-t465 * t437 - t463 * t484);
t394 = -t423 * t485 + t468 * t438;
t453 = t496 * (-t423 * t437 + t438 * t445) + t490 * t394 + t462 * (-t423 * t484 - t468 * t437);
t451 = -pkin(9) - pkin(8);
t425 = t459 * qJD(2);
t422 = t459 * t443;
t420 = qJD(2) * t423;
t414 = t423 * t438 + t437 * t445;
t411 = -t442 * t458 - t444 * t460;
t408 = -t442 * t460 + t444 * t458;
t404 = -t425 * t444 + t442 * t457;
t401 = -t442 * t425 - t444 * t457;
t398 = t437 * t483 + t438 * t463;
t396 = -t437 * t482 + t438 * t464;
t1 = [0 (t405 * t446 + t463 * t475) * r_i_i_C(1) + (t405 * t449 - t463 * t476) * r_i_i_C(2) - t405 * t451 + (t442 * t480 - t489 * t444) * t487 + t456 * t404 - t495 * t411, 0 (-t405 * t447 + (-t442 * t481 - t450 * t463) * qJD(4)) * pkin(4) + t454, t454 (-t389 * t446 - t404 * t449) * r_i_i_C(1) + (-t389 * t449 + t404 * t446) * r_i_i_C(2) + ((-t398 * t449 + t411 * t446) * r_i_i_C(1) + (t398 * t446 + t411 * t449) * r_i_i_C(2)) * qJD(6); 0 (-t403 * t446 + t464 * t475) * r_i_i_C(1) + (-t403 * t449 - t464 * t476) * r_i_i_C(2) + t403 * t451 + (-t489 * t442 - t444 * t480) * t487 + t456 * t401 - t495 * t408, 0 (t403 * t447 + (t444 * t481 - t450 * t464) * qJD(4)) * pkin(4) + t455, t455 (-t387 * t446 - t401 * t449) * r_i_i_C(1) + (-t387 * t449 + t401 * t446) * r_i_i_C(2) + ((-t396 * t449 + t408 * t446) * r_i_i_C(1) + (t396 * t446 + t408 * t449) * r_i_i_C(2)) * qJD(6); 0 (-t419 * t446 + t423 * t475) * r_i_i_C(1) + (-t419 * t449 - t423 * t476) * r_i_i_C(2) + t419 * t451 - t443 * pkin(2) * t478 - t456 * t420 - t495 * t422, 0 (t419 * t447 + (-t423 * t450 - t445 * t447) * qJD(4)) * pkin(4) + t453, t453 (-t394 * t446 + t420 * t449) * r_i_i_C(1) + (-t394 * t449 - t420 * t446) * r_i_i_C(2) + ((-t414 * t449 + t422 * t446) * r_i_i_C(1) + (t414 * t446 + t422 * t449) * r_i_i_C(2)) * qJD(6);];
JaD_transl  = t1;
