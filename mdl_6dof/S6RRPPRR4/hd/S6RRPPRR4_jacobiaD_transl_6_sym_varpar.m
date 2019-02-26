% Zeitableitung der translatorischen Teilmatrix der analytischen Jacobi-Matrix für Segment Nr. 6 (0=Basis) von
% S6RRPPRR4
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
%   pkin=[a2,a3,a4,a5,a6,alpha2,d1,d2,d5,d6,theta3]';
%
% Output:
% JaD_transl [3x6]
%   Zeitableitung der translatorischen Teilmatrix der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-02-26 21:30
% Revision: d75aae1ac561373cd3be920984c3df29a1c2ecc8 (2019-02-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JaD_transl = S6RRPPRR4_jacobiaD_transl_6_sym_varpar(qJ, qJD, r_i_i_C, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(3,1),zeros(11,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPPRR4_jacobiaD_transl_6_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRPPRR4_jacobiaD_transl_6_sym_varpar: qJD has to be [6x1] (double)');
assert(isa(r_i_i_C,'double') && isreal(r_i_i_C) && all(size(r_i_i_C) == [3 1]), ...
	'S6RRPPRR4_jacobiaD_transl_6_sym_varpar: Position vector r_i_i_C has to be [3x1] double');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6RRPPRR4_jacobiaD_transl_6_sym_varpar: pkin has to be [11x1] (double)');

%% Symbolic Calculation
% From jacobiaD_transl_6_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 21:30:13
% EndTime: 2019-02-26 21:30:14
% DurationCPUTime: 0.82s
% Computational Cost: add. (914->117), mult. (2682->194), div. (0->0), fcn. (2944->12), ass. (0->75)
t474 = sin(qJ(2));
t527 = sin(pkin(11));
t529 = cos(pkin(6));
t497 = t529 * t527;
t528 = cos(pkin(11));
t498 = t529 * t528;
t531 = cos(qJ(2));
t481 = -t474 * t497 + t531 * t498;
t451 = t481 * qJD(2);
t484 = t474 * t528 + t531 * t527;
t458 = t484 * qJD(2);
t475 = sin(qJ(1));
t478 = cos(qJ(1));
t460 = t474 * t527 - t531 * t528;
t518 = -t474 * t498 - t531 * t497;
t487 = t475 * t460 + t478 * t518;
t543 = t487 * qJD(1) - t475 * t451 - t478 * t458;
t538 = t475 * t518;
t439 = -t478 * t460 + t538;
t435 = -t475 * t484 + t478 * t481;
t473 = sin(qJ(5));
t477 = cos(qJ(5));
t471 = sin(pkin(6));
t522 = t471 * t478;
t515 = t477 * t522;
t432 = t435 * t473 + t515;
t472 = sin(qJ(6));
t476 = cos(qJ(6));
t542 = -t432 * t472 + t476 * t487;
t541 = t432 * t476 + t472 * t487;
t426 = qJD(1) * t538 - (qJD(1) * t460 - t451) * t478 - t475 * t458;
t499 = t472 * r_i_i_C(1) + t476 * r_i_i_C(2);
t489 = qJD(6) * t499;
t500 = t476 * r_i_i_C(1) - t472 * r_i_i_C(2);
t495 = pkin(5) + t500;
t532 = pkin(10) + r_i_i_C(3);
t540 = qJD(4) - t473 * t489 + (t532 * t473 + t495 * t477) * qJD(5);
t455 = t460 * t471;
t537 = -t455 * t477 + t529 * t473;
t490 = qJD(6) * t500;
t483 = t460 * qJD(2);
t452 = t518 * qJD(2);
t480 = t475 * t481;
t517 = qJD(1) * t478;
t425 = -qJD(1) * t480 + t478 * t452 + t475 * t483 - t484 * t517;
t523 = t471 * t475;
t514 = qJD(1) * t523;
t536 = (-qJD(5) * t435 + t514) * t473 - qJD(5) * t515 + t425 * t477;
t492 = -t435 * t477 + t473 * t522;
t416 = t492 * qJD(5) - t425 * t473 + t477 * t514;
t482 = t495 * t473 - t532 * t477 + qJ(4);
t533 = -pkin(3) - pkin(9);
t530 = pkin(2) * qJD(2);
t516 = t474 * t530;
t513 = t471 * t517;
t511 = t474 * t529;
t512 = -pkin(2) * t511 + (pkin(4) + pkin(8) + qJ(3)) * t471;
t503 = t529 * t531;
t493 = t499 - t533;
t438 = -t478 * t484 - t480;
t430 = -t438 * t473 + t477 * t523;
t491 = -t438 * t477 - t473 * t523;
t441 = t455 * t473 + t529 * t477;
t456 = t484 * t471;
t470 = t531 * pkin(2) + pkin(1);
t459 = -t471 * qJD(3) + t503 * t530;
t450 = t471 * t483;
t449 = qJD(2) * t456;
t427 = t537 * qJD(5) - t449 * t473;
t422 = t435 * qJD(1) + t475 * t452 - t478 * t483;
t419 = t491 * qJD(5) + t422 * t473 + t477 * t513;
t418 = t430 * qJD(5) - t422 * t477 + t473 * t513;
t413 = t419 * t476 + t543 * t472 + (-t430 * t472 + t439 * t476) * qJD(6);
t412 = -t419 * t472 + t543 * t476 + (-t430 * t476 - t439 * t472) * qJD(6);
t1 = [t475 * t516 + t425 * qJ(4) + t435 * qJD(4) - t478 * t459 - t495 * t416 - t493 * t426 - t532 * t536 + (t542 * r_i_i_C(1) - t541 * r_i_i_C(2)) * qJD(6) + (-t478 * t470 - t512 * t475) * qJD(1), t438 * t490 - t493 * t422 + t482 * t543 + ((t475 * t511 - t531 * t478) * qJD(2) + (t474 * t475 - t478 * t503) * qJD(1)) * pkin(2) + t540 * t439, t513, t422, -t495 * t418 + t532 * t419 - t491 * t489, t412 * r_i_i_C(1) - t413 * r_i_i_C(2); -t478 * t516 + t419 * pkin(5) + t413 * r_i_i_C(1) + t412 * r_i_i_C(2) + t422 * qJ(4) - t438 * qJD(4) - t475 * t459 - t533 * t543 + t532 * t418 + (-t470 * t475 + t512 * t478) * qJD(1), t435 * t490 + t493 * t425 + t482 * t426 + ((-t531 * t475 - t478 * t511) * qJD(2) + (-t474 * t478 - t475 * t503) * qJD(1)) * pkin(2) - t540 * t487, t514, -t425, t532 * t416 - t492 * t489 - t495 * t536 (-t416 * t472 + t426 * t476) * r_i_i_C(1) + (-t416 * t476 - t426 * t472) * r_i_i_C(2) + (t541 * r_i_i_C(1) + t542 * r_i_i_C(2)) * qJD(6); 0, -t493 * t449 - t482 * t450 - t455 * t490 + t540 * t456 - t471 * t516, 0, t449, -t532 * t427 + t537 * t489 + t495 * (-t441 * qJD(5) + t449 * t477) (t427 * t472 - t450 * t476) * r_i_i_C(1) + (t427 * t476 + t450 * t472) * r_i_i_C(2) + ((-t441 * t476 - t456 * t472) * r_i_i_C(1) + (t441 * t472 - t456 * t476) * r_i_i_C(2)) * qJD(6);];
JaD_transl  = t1;
