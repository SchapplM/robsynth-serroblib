% Zeitableitung der translatorischen Teilmatrix der analytischen Jacobi-Matrix für Segment Nr. 6 (0=Basis) von
% S6RRPRRR12
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
%   pkin=[a2,a3,a4,a5,a6,alpha2,d1,d2,d4,d5,d6]';
%
% Output:
% JaD_transl [3x6]
%   Zeitableitung der translatorischen Teilmatrix der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-02-26 22:00
% Revision: d75aae1ac561373cd3be920984c3df29a1c2ecc8 (2019-02-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JaD_transl = S6RRPRRR12_jacobiaD_transl_6_sym_varpar(qJ, qJD, r_i_i_C, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(3,1),zeros(11,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPRRR12_jacobiaD_transl_6_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRPRRR12_jacobiaD_transl_6_sym_varpar: qJD has to be [6x1] (double)');
assert(isa(r_i_i_C,'double') && isreal(r_i_i_C) && all(size(r_i_i_C) == [3 1]), ...
	'S6RRPRRR12_jacobiaD_transl_6_sym_varpar: Position vector r_i_i_C has to be [3x1] double');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6RRPRRR12_jacobiaD_transl_6_sym_varpar: pkin has to be [11x1] (double)');

%% Symbolic Calculation
% From jacobiaD_transl_6_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 22:00:27
% EndTime: 2019-02-26 22:00:28
% DurationCPUTime: 0.87s
% Computational Cost: add. (999->114), mult. (1906->182), div. (0->0), fcn. (1888->12), ass. (0->74)
t445 = sin(qJ(6));
t449 = cos(qJ(6));
t509 = r_i_i_C(1) * t445 + r_i_i_C(2) * t449;
t506 = t509 * qJD(6);
t500 = pkin(11) + r_i_i_C(3);
t467 = t449 * r_i_i_C(1) - t445 * r_i_i_C(2);
t464 = pkin(5) + t467;
t444 = cos(pkin(6));
t448 = sin(qJ(1));
t451 = cos(qJ(2));
t486 = t448 * t451;
t447 = sin(qJ(2));
t452 = cos(qJ(1));
t487 = t447 * t452;
t427 = t444 * t487 + t486;
t428 = t444 * t486 + t487;
t417 = t428 * qJD(1) + t427 * qJD(2);
t441 = qJD(4) + qJD(5);
t443 = sin(pkin(6));
t489 = t443 * t452;
t507 = t441 * t489 + t417;
t485 = t451 * t452;
t475 = t444 * t485;
t488 = t447 * t448;
t426 = -t475 + t488;
t442 = qJ(4) + qJ(5);
t439 = sin(t442);
t440 = cos(t442);
t422 = -t426 * t439 + t440 * t489;
t505 = -t422 * t445 - t427 * t449;
t504 = t422 * t449 - t427 * t445;
t468 = qJD(2) * t444 + qJD(1);
t482 = qJD(2) * t451;
t415 = -qJD(1) * t475 - t452 * t482 + t468 * t488;
t491 = t443 * t448;
t502 = t441 * t491 + t415;
t462 = t467 * qJD(6);
t484 = qJD(1) * t443;
t474 = t448 * t484;
t461 = t426 * t441 + t474;
t501 = t461 * t439 - t507 * t440;
t405 = t507 * t439 + t461 * t440;
t446 = sin(qJ(4));
t458 = t446 * pkin(4) + t464 * t439 - t500 * t440 + qJ(3);
t450 = cos(qJ(4));
t499 = pkin(4) * t450;
t496 = -pkin(2) - pkin(10) - pkin(9);
t495 = pkin(8) + pkin(3) + t499;
t492 = t443 * t447;
t490 = t443 * t451;
t483 = qJD(2) * t447;
t479 = t439 * t490;
t476 = t444 * t488;
t473 = t452 * t484;
t472 = t443 * t482;
t471 = t443 * t483;
t463 = t439 * t444 + t440 * t490;
t460 = t428 * t441 + t473;
t459 = t509 - t496;
t407 = t460 * t439 + t502 * t440;
t408 = -t502 * t439 + t460 * t440;
t457 = -t506 * (t428 * t440 - t439 * t491) + t500 * t408 - t464 * t407;
t456 = -t506 * (t426 * t440 + t439 * t489) + t500 * t405 - t464 * t501;
t413 = -t439 * t471 + t463 * t441;
t455 = t506 * t463 + t464 * (t441 * t479 + (-t441 * t444 + t471) * t440) - t500 * t413;
t454 = qJD(4) * t499 + qJD(3) - t439 * t506 + (t500 * t439 + t464 * t440) * t441;
t429 = -t476 + t485;
t425 = t440 * t444 - t479;
t420 = t428 * t439 + t440 * t491;
t418 = -qJD(1) * t476 - t448 * t483 + t468 * t485;
t416 = t427 * qJD(1) + t428 * qJD(2);
t396 = t408 * t449 - t416 * t445 + (-t420 * t445 + t429 * t449) * qJD(6);
t395 = -t408 * t445 - t416 * t449 + (-t420 * t449 - t429 * t445) * qJD(6);
t1 = [-t417 * qJ(3) - t426 * qJD(3) - t464 * t405 - t459 * t418 - t500 * t501 + (t505 * r_i_i_C(1) - t504 * r_i_i_C(2)) * qJD(6) + (-t452 * pkin(1) - t495 * t491) * qJD(1) + (-t417 * t446 + (-t426 * t450 - t446 * t489) * qJD(4)) * pkin(4), t459 * t415 - t458 * t416 - t428 * t462 + t454 * t429, -t415 (-t446 * t473 - t415 * t450 + (-t428 * t446 - t450 * t491) * qJD(4)) * pkin(4) + t457, t457, r_i_i_C(1) * t395 - t396 * r_i_i_C(2); t408 * pkin(5) + t396 * r_i_i_C(1) + t395 * r_i_i_C(2) - t415 * qJ(3) + t428 * qJD(3) + t496 * t416 + t500 * t407 + (-pkin(1) * t448 + t495 * t489) * qJD(1) + (-t415 * t446 + (t428 * t450 - t446 * t491) * qJD(4)) * pkin(4), -t459 * t417 + t458 * t418 - t426 * t462 + t454 * t427, t417 (-t446 * t474 + t417 * t450 + (-t426 * t446 + t450 * t489) * qJD(4)) * pkin(4) + t456, t456 (-t405 * t445 + t418 * t449) * r_i_i_C(1) + (-t405 * t449 - t418 * t445) * r_i_i_C(2) + (t504 * r_i_i_C(1) + t505 * r_i_i_C(2)) * qJD(6); 0 ((t458 * qJD(2) + t462) * t451 + (-t459 * qJD(2) + t454) * t447) * t443, t471 (t450 * t471 + (-t444 * t450 + t446 * t490) * qJD(4)) * pkin(4) + t455, t455 (t413 * t445 + t449 * t472) * r_i_i_C(1) + (t413 * t449 - t445 * t472) * r_i_i_C(2) + ((-t425 * t449 - t445 * t492) * r_i_i_C(1) + (t425 * t445 - t449 * t492) * r_i_i_C(2)) * qJD(6);];
JaD_transl  = t1;
