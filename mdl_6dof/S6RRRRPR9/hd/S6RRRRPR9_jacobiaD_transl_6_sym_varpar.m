% Zeitableitung der translatorischen Teilmatrix der analytischen Jacobi-Matrix für Segment Nr. 6 (0=Basis) von
% S6RRRRPR9
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
%   pkin=[a2,a3,a4,a5,a6,alpha2,d1,d2,d3,d4,d6,theta5]';
%
% Output:
% JaD_transl [3x6]
%   Zeitableitung der translatorischen Teilmatrix der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-02-26 22:35
% Revision: d75aae1ac561373cd3be920984c3df29a1c2ecc8 (2019-02-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JaD_transl = S6RRRRPR9_jacobiaD_transl_6_sym_varpar(qJ, qJD, r_i_i_C, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(3,1),zeros(12,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRRPR9_jacobiaD_transl_6_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRRRPR9_jacobiaD_transl_6_sym_varpar: qJD has to be [6x1] (double)');
assert(isa(r_i_i_C,'double') && isreal(r_i_i_C) && all(size(r_i_i_C) == [3 1]), ...
	'S6RRRRPR9_jacobiaD_transl_6_sym_varpar: Position vector r_i_i_C has to be [3x1] double');
assert(isreal(pkin) && all(size(pkin) == [12 1]), ...
  'S6RRRRPR9_jacobiaD_transl_6_sym_varpar: pkin has to be [12x1] (double)');

%% Symbolic Calculation
% From jacobiaD_transl_6_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 22:35:15
% EndTime: 2019-02-26 22:35:16
% DurationCPUTime: 0.82s
% Computational Cost: add. (1171->121), mult. (1953->193), div. (0->0), fcn. (1954->14), ass. (0->80)
t442 = pkin(12) + qJ(6);
t438 = sin(t442);
t439 = cos(t442);
t469 = t438 * r_i_i_C(1) + t439 * r_i_i_C(2);
t510 = t469 * qJD(6);
t501 = r_i_i_C(3) + pkin(11) + qJ(5);
t436 = cos(pkin(12)) * pkin(5) + pkin(4);
t470 = t439 * r_i_i_C(1) - t438 * r_i_i_C(2);
t467 = t436 + t470;
t447 = cos(pkin(6));
t450 = sin(qJ(2));
t454 = cos(qJ(1));
t490 = t454 * t450;
t451 = sin(qJ(1));
t453 = cos(qJ(2));
t491 = t451 * t453;
t423 = t447 * t490 + t491;
t444 = qJ(3) + qJ(4);
t440 = sin(t444);
t446 = sin(pkin(6));
t493 = t446 * t454;
t430 = t440 * t493;
t441 = cos(t444);
t416 = t423 * t441 - t430;
t489 = t454 * t453;
t478 = t447 * t489;
t492 = t451 * t450;
t422 = -t478 + t492;
t509 = t416 * t438 - t422 * t439;
t508 = t416 * t439 + t422 * t438;
t443 = qJD(3) + qJD(4);
t449 = sin(qJ(3));
t456 = (t467 * t440 - t501 * t441) * t443 + qJD(3) * t449 * pkin(3) + t441 * t510 - t440 * qJD(5);
t471 = qJD(2) * t447 + qJD(1);
t480 = t447 * t492;
t487 = qJD(2) * t450;
t414 = -qJD(1) * t480 - t451 * t487 + t471 * t489;
t488 = qJD(1) * t446;
t477 = t451 * t488;
t479 = t441 * t493;
t402 = (-t423 * t443 + t477) * t440 + t414 * t441 - t443 * t479;
t452 = cos(qJ(3));
t437 = t452 * pkin(3) + pkin(2);
t504 = t501 * t440 + t467 * t441 + t437;
t498 = t441 * t443;
t497 = t446 * t450;
t496 = t446 * t451;
t495 = t446 * t452;
t494 = t446 * t453;
t486 = qJD(2) * t453;
t482 = t440 * t497;
t481 = t441 * t497;
t476 = t454 * t488;
t475 = t446 * t487;
t474 = t446 * t486;
t473 = -sin(pkin(12)) * pkin(5) - pkin(10) - pkin(9);
t424 = t447 * t491 + t490;
t412 = t423 * qJD(1) + t424 * qJD(2);
t468 = t443 * t496 - t412;
t465 = t423 * t440 + t479;
t464 = t470 * qJD(6);
t462 = t443 * t447 + t474;
t461 = t469 - t473;
t401 = t414 * t440 + t423 * t498 - t443 * t430 - t441 * t477;
t459 = t416 * qJD(5) - t401 * t467 + t501 * t402 + t510 * t465;
t425 = -t480 + t489;
t399 = t425 * t498 + t468 * t440 - t441 * t476;
t400 = t468 * t441 + (-t425 * t443 + t476) * t440;
t418 = -t425 * t440 + t441 * t496;
t419 = t425 * t441 + t440 * t496;
t458 = t419 * qJD(5) - t399 * t467 + t501 * t400 - t510 * t418;
t409 = t462 * t440 + t443 * t481;
t410 = t462 * t441 - t443 * t482;
t421 = t447 * t440 + t481;
t457 = t421 * qJD(5) - t510 * (t447 * t441 - t482) + t501 * t410 - t467 * t409;
t413 = t424 * qJD(1) + t423 * qJD(2);
t411 = -qJD(1) * t478 - t454 * t486 + t471 * t492;
t390 = t400 * t439 - t411 * t438 + (-t419 * t438 + t424 * t439) * qJD(6);
t389 = -t400 * t438 - t411 * t439 + (-t419 * t439 - t424 * t438) * qJD(6);
t1 = [-t465 * qJD(5) - t414 * t437 - t467 * t402 - t461 * t413 - t501 * t401 + (t509 * r_i_i_C(1) + t508 * r_i_i_C(2)) * qJD(6) + (-t454 * pkin(1) - pkin(8) * t496) * qJD(1) + (-t449 * t477 + (t423 * t449 + t452 * t493) * qJD(3)) * pkin(3), t504 * t411 - t461 * t412 + t456 * t424 + t425 * t464 (t452 * t476 + t412 * t449 + (-t425 * t452 - t449 * t496) * qJD(3)) * pkin(3) + t458, t458, t399, t389 * r_i_i_C(1) - t390 * r_i_i_C(2); t390 * r_i_i_C(1) + t389 * r_i_i_C(2) - t418 * qJD(5) + t400 * t436 - t412 * t437 + t473 * t411 + t501 * t399 + (-pkin(1) * t451 + pkin(8) * t493) * qJD(1) + (t449 * t476 + (-t425 * t449 + t451 * t495) * qJD(3)) * pkin(3), -t413 * t504 + t461 * t414 + t456 * t422 + t423 * t464 (t452 * t477 - t414 * t449 + (-t423 * t452 + t449 * t493) * qJD(3)) * pkin(3) + t459, t459, t401 (-t402 * t438 + t413 * t439) * r_i_i_C(1) + (-t402 * t439 - t413 * t438) * r_i_i_C(2) + (-t508 * r_i_i_C(1) + t509 * r_i_i_C(2)) * qJD(6); 0 ((-qJD(2) * t504 + t464) * t450 + (t461 * qJD(2) - t456) * t453) * t446 (-t449 * t474 + (-t447 * t449 - t450 * t495) * qJD(3)) * pkin(3) + t457, t457, t409 (-t410 * t438 + t439 * t475) * r_i_i_C(1) + (-t410 * t439 - t438 * t475) * r_i_i_C(2) + ((-t421 * t439 + t438 * t494) * r_i_i_C(1) + (t421 * t438 + t439 * t494) * r_i_i_C(2)) * qJD(6);];
JaD_transl  = t1;
