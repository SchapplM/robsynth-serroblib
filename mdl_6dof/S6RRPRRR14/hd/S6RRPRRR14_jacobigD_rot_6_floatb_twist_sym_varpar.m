% Zeitableitung der rotatorischen Teilmatrix der geometrischen Jacobi-Matrix für Segment Nr. 6 (0=Basis) von
% S6RRPRRR14
% Use Code from Maple symbolic Code Generation
%
% geometrische Jacobi-Matrix: Differentieller Zusammenhang zwischen
% Endeffektorposition und verallgemeinerten Koordinaten.
%
%
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [14x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,alpha3,alpha4,d1,d2,d4,d5,d6,theta3]';
%
% Output:
% JgD_rot [3x6]
%   Zeitableitung der rotatorischen Teilmatrix der geometrischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox (ehem. IRT-Maple-Toolbox)
% Datum: 2019-01-03 10:26
% Revision: 5fdbc45bcf2cc60deefd7ac2d71d743ed41bf7e4 (2018-12-21)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für mechatronische Systeme, Universität Hannover

function JgD_rot = S6RRPRRR14_jacobigD_rot_6_floatb_twist_sym_varpar(qJ, qJD, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(14,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPRRR14_jacobigD_rot_6_floatb_twist_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRPRRR14_jacobigD_rot_6_floatb_twist_sym_varpar: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [14 1]), ...
  'S6RRPRRR14_jacobigD_rot_6_floatb_twist_sym_varpar: pkin has to be [14x1] (double)');

%% Symbolic Calculation
% From jacobigD_rot_6_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-01-03 10:25:36
% EndTime: 2019-01-03 10:25:37
% DurationCPUTime: 0.31s
% Computational Cost: add. (278->73), mult. (900->150), div. (0->0), fcn. (999->16), ass. (0->73)
t453 = sin(pkin(8));
t461 = sin(qJ(4));
t496 = t453 * t461;
t454 = sin(pkin(7));
t459 = cos(pkin(6));
t495 = t454 * t459;
t455 = sin(pkin(6));
t463 = sin(qJ(1));
t494 = t455 * t463;
t467 = cos(qJ(1));
t493 = t455 * t467;
t457 = cos(pkin(8));
t492 = t457 * t461;
t458 = cos(pkin(7));
t462 = sin(qJ(2));
t491 = t458 * t462;
t466 = cos(qJ(2));
t490 = t458 * t466;
t489 = t463 * t462;
t488 = t463 * t466;
t487 = t467 * t462;
t486 = t467 * t466;
t485 = qJD(1) * t455;
t484 = qJD(2) * t455;
t460 = sin(qJ(5));
t483 = qJD(4) * t460;
t482 = t463 * t485;
t481 = t467 * t485;
t480 = t454 * t462 * t484;
t479 = t453 * t480;
t449 = t459 * t487 + t488;
t452 = sin(pkin(14));
t456 = cos(pkin(14));
t448 = t459 * t486 - t489;
t475 = t448 * t458 - t454 * t493;
t430 = -t449 * t452 + t456 * t475;
t443 = -t448 * t454 - t458 * t493;
t478 = t430 * t457 + t443 * t453;
t473 = t459 * t489 - t486;
t450 = -t459 * t488 - t487;
t474 = t450 * t458 + t454 * t494;
t432 = t452 * t473 + t456 * t474;
t444 = -t450 * t454 + t458 * t494;
t477 = t432 * t457 + t444 * t453;
t441 = t456 * t495 + (-t452 * t462 + t456 * t490) * t455;
t447 = -t455 * t466 * t454 + t459 * t458;
t476 = t441 * t457 + t447 * t453;
t437 = -qJD(1) * t448 + qJD(2) * t473;
t472 = t437 * t458 + t454 * t481;
t439 = qJD(1) * t450 - qJD(2) * t449;
t471 = t439 * t458 + t454 * t482;
t431 = t449 * t456 + t452 * t475;
t465 = cos(qJ(4));
t470 = t431 * t465 + t461 * t478;
t433 = t452 * t474 - t456 * t473;
t469 = t433 * t465 + t461 * t477;
t442 = t455 * t462 * t456 + (t455 * t490 + t495) * t452;
t468 = t442 * t465 + t461 * t476;
t464 = cos(qJ(5));
t446 = (-t452 * t491 + t456 * t466) * t484;
t445 = (-t452 * t466 - t456 * t491) * t484;
t440 = -qJD(1) * t473 + qJD(2) * t448;
t438 = -qJD(1) * t449 + qJD(2) * t450;
t436 = -t445 * t453 + t457 * t480;
t435 = -t439 * t454 + t458 * t482;
t434 = -t437 * t454 + t458 * t481;
t429 = t440 * t456 + t452 * t471;
t428 = -t440 * t452 + t456 * t471;
t427 = t438 * t456 + t452 * t472;
t426 = -t438 * t452 + t456 * t472;
t425 = -t428 * t453 + t435 * t457;
t424 = -t426 * t453 + t434 * t457;
t1 = [0, t481, 0, t424, t427 * t461 + (-t426 * t457 - t434 * t453) * t465 + t469 * qJD(4) (t426 * t492 + t427 * t465 + t434 * t496) * t460 - t424 * t464 + (t469 * t464 + (-t432 * t453 + t444 * t457) * t460) * qJD(5) + (-t433 * t461 + t465 * t477) * t483; 0, t482, 0, t425, t429 * t461 + (-t428 * t457 - t435 * t453) * t465 + t470 * qJD(4) (t428 * t492 + t429 * t465 + t435 * t496) * t460 - t425 * t464 + (t470 * t464 + (-t430 * t453 + t443 * t457) * t460) * qJD(5) + (-t431 * t461 + t465 * t478) * t483; 0, 0, 0, t436, t446 * t461 + (-t445 * t457 - t479) * t465 + t468 * qJD(4) (t445 * t492 + t446 * t465 + t461 * t479) * t460 - t436 * t464 + (t468 * t464 + (-t441 * t453 + t447 * t457) * t460) * qJD(5) + (-t442 * t461 + t465 * t476) * t483;];
JgD_rot  = t1;
