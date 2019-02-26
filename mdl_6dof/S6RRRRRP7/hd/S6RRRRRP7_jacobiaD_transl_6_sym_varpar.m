% Zeitableitung der translatorischen Teilmatrix der analytischen Jacobi-Matrix für Segment Nr. 6 (0=Basis) von
% S6RRRRRP7
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
%   pkin=[a2,a3,a4,a5,a6,alpha2,d1,d2,d3,d4,d5]';
%
% Output:
% JaD_transl [3x6]
%   Zeitableitung der translatorischen Teilmatrix der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-02-26 22:43
% Revision: d75aae1ac561373cd3be920984c3df29a1c2ecc8 (2019-02-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JaD_transl = S6RRRRRP7_jacobiaD_transl_6_sym_varpar(qJ, qJD, r_i_i_C, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(3,1),zeros(11,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRRRP7_jacobiaD_transl_6_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRRRRP7_jacobiaD_transl_6_sym_varpar: qJD has to be [6x1] (double)');
assert(isa(r_i_i_C,'double') && isreal(r_i_i_C) && all(size(r_i_i_C) == [3 1]), ...
	'S6RRRRRP7_jacobiaD_transl_6_sym_varpar: Position vector r_i_i_C has to be [3x1] double');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6RRRRRP7_jacobiaD_transl_6_sym_varpar: pkin has to be [11x1] (double)');

%% Symbolic Calculation
% From jacobiaD_transl_6_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 22:43:15
% EndTime: 2019-02-26 22:43:16
% DurationCPUTime: 0.86s
% Computational Cost: add. (1196->133), mult. (2192->211), div. (0->0), fcn. (2194->12), ass. (0->83)
t440 = qJ(3) + qJ(4);
t437 = sin(t440);
t438 = cos(t440);
t439 = qJD(3) + qJD(4);
t444 = sin(qJ(3));
t443 = sin(qJ(5));
t447 = cos(qJ(5));
t505 = pkin(5) + r_i_i_C(1);
t465 = t447 * r_i_i_C(2) + t505 * t443;
t458 = t465 * qJD(5);
t435 = t447 * pkin(5) + pkin(4);
t503 = t443 * r_i_i_C(2);
t467 = t447 * r_i_i_C(1) + t435 - t503;
t500 = r_i_i_C(3) + qJ(6) + pkin(11);
t451 = -(t500 * t439 - t458) * t438 + qJD(3) * t444 * pkin(3) + (t467 * t439 - qJD(6)) * t437;
t450 = -pkin(10) - pkin(9);
t513 = t450 - t465;
t445 = sin(qJ(2));
t446 = sin(qJ(1));
t449 = cos(qJ(2));
t499 = cos(pkin(6));
t504 = cos(qJ(1));
t473 = t499 * t504;
t422 = t445 * t473 + t446 * t449;
t441 = sin(pkin(6));
t484 = t441 * t504;
t429 = t437 * t484;
t415 = t422 * t438 - t429;
t468 = t449 * t473;
t491 = t446 * t445;
t421 = -t468 + t491;
t511 = t415 * t443 - t421 * t447;
t510 = t415 * t447 + t421 * t443;
t478 = t446 * t499;
t474 = t445 * t478;
t479 = t504 * qJD(1);
t490 = qJD(2) * t445;
t413 = -qJD(1) * t474 - t446 * t490 + (qJD(2) * t473 + t479) * t449;
t476 = t438 * t484;
t494 = t441 * t446;
t482 = qJD(1) * t494;
t401 = (-t422 * t439 + t482) * t437 + t413 * t438 - t439 * t476;
t448 = cos(qJ(3));
t436 = t448 * pkin(3) + pkin(2);
t506 = t500 * t437 + t467 * t438 + t436;
t496 = t438 * t439;
t495 = t441 * t445;
t493 = t441 * t448;
t492 = t441 * t449;
t489 = qJD(5) * t443;
t488 = qJD(5) * t447;
t486 = t437 * t495;
t485 = t438 * t495;
t483 = t504 * t449;
t481 = t441 * t490;
t480 = qJD(2) * t492;
t475 = t441 * t479;
t423 = t504 * t445 + t449 * t478;
t411 = t422 * qJD(1) + t423 * qJD(2);
t472 = t439 * t494 - t411;
t412 = t423 * qJD(1) + t422 * qJD(2);
t471 = -t401 * t443 + t412 * t447;
t420 = t499 * t437 + t485;
t464 = -t420 * t447 + t443 * t492;
t460 = t422 * t437 + t476;
t457 = t499 * t439 + t480;
t409 = t457 * t438 - t439 * t486;
t459 = -t409 * t443 + t447 * t481;
t400 = t413 * t437 + t422 * t496 - t439 * t429 - t438 * t482;
t410 = -qJD(1) * t468 - qJD(2) * t483 + (qJD(2) * t499 + qJD(1)) * t491;
t424 = t483 - t474;
t418 = t424 * t438 + t437 * t494;
t456 = -t410 * t443 + (-t418 * t443 + t423 * t447) * qJD(5);
t399 = t472 * t438 + (-t424 * t439 + t475) * t437;
t388 = -t399 * t443 - t410 * t447 + (-t418 * t447 - t423 * t443) * qJD(5);
t398 = t424 * t496 + t472 * t437 - t438 * t475;
t417 = -t424 * t437 + t438 * t494;
t454 = t418 * qJD(6) - t467 * t398 + t500 * t399 - t417 * t458;
t453 = t415 * qJD(6) - t467 * t400 + t500 * t401 + t460 * t458;
t408 = t457 * t437 + t439 * t485;
t452 = t420 * qJD(6) - (t499 * t438 - t486) * t458 + t500 * t409 - t467 * t408;
t389 = t399 * t447 + t456;
t1 = [-t460 * qJD(6) - t413 * t436 - t467 * t401 + t513 * t412 - t500 * t400 + (-t504 * pkin(1) - pkin(8) * t494) * qJD(1) + (-t444 * t482 + (t422 * t444 + t448 * t484) * qJD(3)) * pkin(3) + (t510 * r_i_i_C(2) + t505 * t511) * qJD(5) (-t411 * t447 - t424 * t489) * r_i_i_C(2) + t411 * t450 + t506 * t410 + t451 * t423 + t505 * (-t411 * t443 + t424 * t488) (t448 * t475 + t411 * t444 + (-t424 * t448 - t444 * t494) * qJD(3)) * pkin(3) + t454, t454, -t389 * r_i_i_C(2) + t505 * t388, t398; t389 * r_i_i_C(1) + t388 * r_i_i_C(2) - t417 * qJD(6) + t399 * t435 + t410 * t450 - t411 * t436 + t500 * t398 + (-t446 * pkin(1) + pkin(8) * t484) * qJD(1) + t456 * pkin(5) + (t444 * t475 + (-t424 * t444 + t446 * t493) * qJD(3)) * pkin(3) (t413 * t447 - t422 * t489) * r_i_i_C(2) - t413 * t450 - t506 * t412 + t451 * t421 + t505 * (t413 * t443 + t422 * t488) (t448 * t482 - t413 * t444 + (-t422 * t448 + t444 * t484) * qJD(3)) * pkin(3) + t453, t453, t471 * r_i_i_C(1) + (-t401 * t447 - t412 * t443) * r_i_i_C(2) + (-r_i_i_C(1) * t510 + t511 * r_i_i_C(2)) * qJD(5) + (-qJD(5) * t510 + t471) * pkin(5), t400; 0 (((t505 * t447 - t503) * qJD(5) - t506 * qJD(2)) * t445 + (-t513 * qJD(2) - t451) * t449) * t441 (-t444 * t480 + (-t499 * t444 - t445 * t493) * qJD(3)) * pkin(3) + t452, t452, t459 * r_i_i_C(1) + (-t409 * t447 - t443 * t481) * r_i_i_C(2) + (t464 * r_i_i_C(1) + (t420 * t443 + t447 * t492) * r_i_i_C(2)) * qJD(5) + (t464 * qJD(5) + t459) * pkin(5), t408;];
JaD_transl  = t1;
