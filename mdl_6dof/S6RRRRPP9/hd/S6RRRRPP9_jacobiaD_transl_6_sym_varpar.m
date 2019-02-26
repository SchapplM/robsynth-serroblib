% Zeitableitung der translatorischen Teilmatrix der analytischen Jacobi-Matrix für Segment Nr. 6 (0=Basis) von
% S6RRRRPP9
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
%   pkin=[a2,a3,a4,a5,a6,alpha2,d1,d2,d3,d4]';
%
% Output:
% JaD_transl [3x6]
%   Zeitableitung der translatorischen Teilmatrix der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-02-26 22:30
% Revision: d75aae1ac561373cd3be920984c3df29a1c2ecc8 (2019-02-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JaD_transl = S6RRRRPP9_jacobiaD_transl_6_sym_varpar(qJ, qJD, r_i_i_C, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(3,1),zeros(10,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRRPP9_jacobiaD_transl_6_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRRRPP9_jacobiaD_transl_6_sym_varpar: qJD has to be [6x1] (double)');
assert(isa(r_i_i_C,'double') && isreal(r_i_i_C) && all(size(r_i_i_C) == [3 1]), ...
	'S6RRRRPP9_jacobiaD_transl_6_sym_varpar: Position vector r_i_i_C has to be [3x1] double');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RRRRPP9_jacobiaD_transl_6_sym_varpar: pkin has to be [10x1] (double)');

%% Symbolic Calculation
% From jacobiaD_transl_6_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 22:30:18
% EndTime: 2019-02-26 22:30:19
% DurationCPUTime: 0.93s
% Computational Cost: add. (1056->135), mult. (3096->219), div. (0->0), fcn. (3220->10), ass. (0->80)
t442 = sin(qJ(1));
t441 = sin(qJ(2));
t502 = cos(pkin(6));
t469 = t442 * t502;
t466 = t441 * t469;
t486 = qJD(2) * t441;
t445 = cos(qJ(2));
t446 = cos(qJ(1));
t488 = t446 * t445;
t414 = -qJD(1) * t466 - t442 * t486 + (qJD(2) * t502 + qJD(1)) * t488;
t468 = t446 * t502;
t426 = t441 * t468 + t442 * t445;
t440 = sin(qJ(3));
t444 = cos(qJ(3));
t438 = sin(pkin(6));
t487 = qJD(1) * t438;
t475 = t442 * t487;
t493 = t438 * t446;
t477 = t444 * t493;
t400 = (-qJD(3) * t426 + t475) * t440 - qJD(3) * t477 + t414 * t444;
t427 = t446 * t441 + t445 * t469;
t413 = t427 * qJD(1) + t426 * qJD(2);
t439 = sin(qJ(4));
t443 = cos(qJ(4));
t420 = -t426 * t444 + t440 * t493;
t425 = t442 * t441 - t445 * t468;
t461 = t420 * t443 - t425 * t439;
t510 = t461 * qJD(4) - t400 * t439 + t413 * t443;
t507 = -t420 * t439 - t425 * t443;
t508 = qJD(4) * t507 - t400 * t443 - t413 * t439;
t479 = r_i_i_C(1) + pkin(5) + pkin(10);
t506 = t444 * pkin(3) + t479 * t440 + pkin(2);
t505 = -pkin(3) * t440 + t479 * t444;
t478 = r_i_i_C(3) + qJ(6) + pkin(4);
t503 = r_i_i_C(2) + qJ(5);
t449 = t503 * t439 + t478 * t443 + pkin(3);
t497 = t438 * t441;
t496 = t438 * t442;
t495 = t438 * t444;
t494 = t438 * t445;
t492 = t439 * t444;
t491 = t443 * t444;
t490 = t443 * t445;
t489 = t444 * t445;
t485 = qJD(2) * t444;
t484 = qJD(3) * t440;
t483 = qJD(3) * t445;
t482 = qJD(4) * t439;
t481 = qJD(4) * t443;
t480 = qJD(4) * t444;
t474 = t446 * t487;
t473 = t438 * t486;
t472 = qJD(2) * t494;
t471 = t445 * t482;
t470 = t440 * t483;
t412 = t426 * qJD(1) + t427 * qJD(2);
t465 = t427 * t480 - t412;
t464 = t425 * t480 + t414;
t453 = t466 - t488;
t422 = t440 * t496 - t444 * t453;
t459 = t422 * t443 + t427 * t439;
t424 = t502 * t440 + t441 * t495;
t458 = t424 * t439 + t438 * t490;
t457 = t440 * t453 + t442 * t495;
t454 = -t440 * t497 + t502 * t444;
t452 = qJD(3) * t505;
t411 = t425 * qJD(1) + t453 * qJD(2);
t451 = -qJD(4) * t453 + t411 * t444 + t427 * t484;
t450 = qJD(4) * t426 - t413 * t444 + t425 * t484;
t448 = t420 * qJD(3) - t414 * t440 + t444 * t475;
t447 = t439 * qJD(5) + t443 * qJD(6) + (-t478 * t439 + t503 * t443) * qJD(4);
t417 = t454 * qJD(3) + t444 * t472;
t407 = t422 * t439 - t427 * t443;
t404 = -t458 * qJD(4) + t417 * t443 + t439 * t473;
t403 = t417 * t439 + t424 * t481 - t438 * t471 - t443 * t473;
t398 = t457 * qJD(3) - t412 * t444 + t440 * t474;
t397 = t422 * qJD(3) - t412 * t440 - t444 * t474;
t388 = -t411 * t439 - t422 * t482 + (qJD(4) * t427 + t398) * t443;
t387 = t459 * qJD(4) + t398 * t439 + t411 * t443;
t1 = [t461 * qJD(6) - t507 * qJD(5) - t400 * pkin(3) - t414 * pkin(2) - t413 * pkin(9) + t503 * t510 + (-pkin(1) * t446 - pkin(8) * t496) * qJD(1) + t479 * t448 + t478 * t508 -(t427 * t491 + t439 * t453) * qJD(6) - (t427 * t492 - t443 * t453) * qJD(5) - t412 * pkin(9) + t503 * (t451 * t439 - t465 * t443) + t478 * (t465 * t439 + t451 * t443) - t427 * t452 + t506 * t411, -t449 * t397 + t479 * t398 + t447 * t457, qJD(5) * t459 - t407 * qJD(6) - t478 * t387 + t503 * t388, t387, t388; -t412 * pkin(2) + t398 * pkin(3) - t411 * pkin(9) + t407 * qJD(5) + t459 * qJD(6) + t503 * t387 + (-pkin(1) * t442 + pkin(8) * t493) * qJD(1) + t479 * t397 + t478 * t388 -(t425 * t491 - t426 * t439) * qJD(6) - (t425 * t492 + t426 * t443) * qJD(5) + t414 * pkin(9) + t503 * (t450 * t439 - t464 * t443) + t478 * (t464 * t439 + t450 * t443) - t425 * t452 - t506 * t413, t479 * t400 + t449 * t448 + t447 * (-t426 * t440 - t477) -qJD(5) * t461 - qJD(6) * t507 + t478 * t510 - t503 * t508, -t510, -t508; 0, -t478 * (-t439 * t472 - t481 * t497) + (-t503 * ((qJD(2) - t480) * t490 + (t470 + (-qJD(4) + t485) * t441) * t439) - t478 * (t444 * t471 + (t441 * t485 + t470) * t443) - (-t439 * t441 - t443 * t489) * qJD(6) - (-t439 * t489 + t441 * t443) * qJD(5) + t505 * t483 + (t445 * pkin(9) - t441 * t506) * qJD(2)) * t438, t479 * t417 + t449 * (-t424 * qJD(3) - t440 * t472) + t447 * t454, -t458 * qJD(6) - (-t424 * t443 + t439 * t494) * qJD(5) + t503 * t404 - t478 * t403, t403, t404;];
JaD_transl  = t1;
