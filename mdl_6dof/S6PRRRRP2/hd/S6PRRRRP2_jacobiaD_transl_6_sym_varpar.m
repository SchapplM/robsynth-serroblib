% Zeitableitung der translatorischen Teilmatrix der analytischen Jacobi-Matrix für Segment Nr. 6 (0=Basis) von
% S6PRRRRP2
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
%   pkin=[a2,a3,a4,a5,a6,alpha2,d2,d3,d4,d5,theta1]';
%
% Output:
% JaD_transl [3x6]
%   Zeitableitung der translatorischen Teilmatrix der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-02-26 20:15
% Revision: d75aae1ac561373cd3be920984c3df29a1c2ecc8 (2019-02-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JaD_transl = S6PRRRRP2_jacobiaD_transl_6_sym_varpar(qJ, qJD, r_i_i_C, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(3,1),zeros(11,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRRRRP2_jacobiaD_transl_6_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6PRRRRP2_jacobiaD_transl_6_sym_varpar: qJD has to be [6x1] (double)');
assert(isa(r_i_i_C,'double') && isreal(r_i_i_C) && all(size(r_i_i_C) == [3 1]), ...
	'S6PRRRRP2_jacobiaD_transl_6_sym_varpar: Position vector r_i_i_C has to be [3x1] double');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6PRRRRP2_jacobiaD_transl_6_sym_varpar: pkin has to be [11x1] (double)');

%% Symbolic Calculation
% From jacobiaD_transl_6_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 20:15:45
% EndTime: 2019-02-26 20:15:46
% DurationCPUTime: 0.62s
% Computational Cost: add. (1048->111), mult. (1833->184), div. (0->0), fcn. (1881->12), ass. (0->79)
t462 = cos(qJ(5));
t519 = pkin(5) + r_i_i_C(1);
t528 = t519 * t462 + pkin(4);
t517 = r_i_i_C(3) + qJ(6);
t518 = pkin(10) + r_i_i_C(2);
t459 = sin(qJ(5));
t500 = qJD(6) * t459;
t527 = -qJD(5) * t459 * t519 + t500;
t456 = qJ(3) + qJ(4);
t454 = cos(t456);
t455 = qJD(3) + qJD(4);
t460 = sin(qJ(3));
t453 = sin(t456);
t508 = t453 * t455;
t522 = -qJD(3) * t460 * pkin(3) - pkin(4) * t508 + (t518 * t455 + t500) * t454;
t461 = sin(qJ(2));
t464 = cos(qJ(2));
t457 = sin(pkin(11));
t516 = cos(pkin(6));
t489 = t457 * t516;
t515 = cos(pkin(11));
t444 = -t461 * t489 + t515 * t464;
t520 = (qJD(2) * t454 - qJD(5)) * t461 + t464 * t508;
t507 = t454 * t455;
t458 = sin(pkin(6));
t506 = t457 * t458;
t505 = t458 * t461;
t504 = t458 * t464;
t503 = qJD(2) * t461;
t502 = qJD(5) * t454;
t501 = qJD(5) * t462;
t499 = qJD(6) * t462;
t498 = t453 * t505;
t497 = t454 * t505;
t496 = t459 * t504;
t494 = t458 * t503;
t493 = qJD(2) * t504;
t488 = t458 * t515;
t483 = t454 * t488;
t443 = t515 * t461 + t464 * t489;
t439 = t443 * qJD(2);
t482 = t455 * t506 - t439;
t479 = t516 * t515;
t475 = t464 * t479;
t437 = -qJD(2) * t475 + t457 * t503;
t441 = t457 * t461 - t475;
t481 = t441 * t502 - t437;
t480 = t443 * t502 - t439;
t442 = t457 * t464 + t461 * t479;
t425 = t442 * t454 - t453 * t488;
t478 = t425 * t462 + t441 * t459;
t427 = t444 * t454 + t453 * t506;
t477 = t427 * t462 + t443 * t459;
t476 = (qJD(2) - t502) * t464;
t463 = cos(qJ(3));
t473 = -t463 * pkin(3) - t454 * pkin(4) - t518 * t453 - pkin(2);
t472 = t516 * t455 + t493;
t438 = t442 * qJD(2);
t471 = qJD(5) * t442 - t438 * t454 + t441 * t508;
t440 = t444 * qJD(2);
t470 = qJD(5) * t444 - t440 * t454 + t443 * t508;
t406 = -t442 * t507 + (t455 * t488 + t437) * t453;
t407 = -t437 * t454 - t442 * t508 - t455 * t483;
t424 = -t442 * t453 - t483;
t468 = t407 * t518 + t424 * t527 + t517 * (t406 * t459 + t424 * t501) + t528 * t406;
t408 = -t444 * t507 - t482 * t453;
t409 = -t444 * t508 + t482 * t454;
t426 = -t444 * t453 + t454 * t506;
t467 = t409 * t518 + t426 * t527 + t517 * (t408 * t459 + t426 * t501) + t528 * t408;
t418 = -t472 * t453 - t455 * t497;
t419 = t472 * t454 - t455 * t498;
t435 = t516 * t454 - t498;
t466 = t419 * t518 + t435 * t527 + t517 * (t418 * t459 + t435 * t501) + t528 * t418;
t465 = -pkin(9) - pkin(8);
t436 = t516 * t453 + t497;
t390 = -qJD(5) * t496 + t419 * t459 + t436 * t501 - t462 * t494;
t384 = t477 * qJD(5) + t409 * t459 - t440 * t462;
t382 = t478 * qJD(5) + t407 * t459 - t438 * t462;
t1 = [0, -t444 * t499 + t439 * t465 + t519 * (t480 * t459 + t470 * t462) + t517 * (t470 * t459 - t480 * t462) - t522 * t443 + t473 * t440 (t439 * t460 + (-t444 * t463 - t460 * t506) * qJD(3)) * pkin(3) + t467, t467, t477 * qJD(6) + t517 * (t409 * t462 + t440 * t459 + (-t427 * t459 + t443 * t462) * qJD(5)) - t519 * t384, t384; 0, -t442 * t499 + t437 * t465 + t519 * (t481 * t459 + t471 * t462) + t517 * (t471 * t459 - t481 * t462) - t522 * t441 + t473 * t438 (t437 * t460 + (-t442 * t463 + t460 * t488) * qJD(3)) * pkin(3) + t468, t468, t478 * qJD(6) + t517 * (t407 * t462 + t438 * t459 + (-t425 * t459 + t441 * t462) * qJD(5)) - t519 * t382, t382; 0 (t519 * (t459 * t476 - t520 * t462) - t517 * (t520 * t459 + t462 * t476) - t461 * t499 + t522 * t464 + (t473 * t461 - t464 * t465) * qJD(2)) * t458 (-t460 * t493 + (-t516 * t460 - t463 * t505) * qJD(3)) * pkin(3) + t466, t466 -(-t436 * t462 + t496) * qJD(6) + t517 * (t459 * t494 + t419 * t462 + (-t436 * t459 - t462 * t504) * qJD(5)) - t519 * t390, t390;];
JaD_transl  = t1;
