% Zeitableitung der rotatorischen Teilmatrix der geometrischen Jacobi-Matrix für Segment Nr. 5 (0=Basis) von
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
% Datum: 2018-12-10 18:38
% Revision: bb42a8b95257d9bc83910d26e849f5825122f662 (2018-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für mechatronische Systeme, Universität Hannover

function JgD_rot = S6RRPRRR14_jacobigD_rot_5_floatb_twist_sym_varpar(qJ, qJD, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(14,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPRRR14_jacobigD_rot_5_floatb_twist_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRPRRR14_jacobigD_rot_5_floatb_twist_sym_varpar: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [14 1]), ...
  'S6RRPRRR14_jacobigD_rot_5_floatb_twist_sym_varpar: pkin has to be [14x1] (double)');

%% Symbolic Calculation
% From jacobigD_rot_5_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2018-12-10 18:38:23
% EndTime: 2018-12-10 18:38:23
% DurationCPUTime: 0.24s
% Computational Cost: add. (392->79), mult. (495->126), div. (0->0), fcn. (396->26), ass. (0->77)
t449 = pkin(8) - qJ(4);
t479 = sin(t449) / 0.2e1;
t448 = pkin(8) + qJ(4);
t478 = cos(t448) / 0.2e1;
t477 = qJD(2) / 0.2e1;
t450 = pkin(6) + qJ(2);
t444 = cos(t450);
t431 = t444 * t477;
t451 = pkin(6) - qJ(2);
t445 = cos(t451);
t471 = qJD(2) * t445;
t419 = t431 - t471 / 0.2e1;
t454 = sin(pkin(7));
t476 = t419 * t454;
t455 = sin(pkin(6));
t462 = sin(qJ(1));
t475 = t455 * t462;
t465 = cos(qJ(1));
t474 = t455 * t465;
t473 = qJD(1) * t455;
t440 = sin(t450);
t472 = qJD(2) * t440;
t470 = qJD(2) * t462;
t469 = qJD(2) * t465;
t468 = qJD(4) * cos(qJ(4));
t467 = t462 * t473;
t466 = t465 * t473;
t434 = t440 / 0.2e1;
t441 = sin(t451);
t426 = t434 - t441 / 0.2e1;
t464 = cos(qJ(2));
t411 = t465 * t426 + t462 * t464;
t413 = -t462 * t426 + t465 * t464;
t435 = t445 / 0.2e1;
t428 = t435 + t444 / 0.2e1;
t461 = sin(qJ(2));
t410 = t465 * t428 - t462 * t461;
t412 = -t462 * t428 - t465 * t461;
t460 = sin(qJ(4));
t459 = cos(pkin(6));
t458 = cos(pkin(7));
t457 = cos(pkin(8));
t456 = cos(pkin(14));
t453 = sin(pkin(8));
t452 = sin(pkin(14));
t447 = pkin(7) - pkin(14);
t446 = pkin(7) + pkin(14);
t443 = cos(t449);
t438 = sin(t448);
t437 = cos(t446);
t436 = sin(t447);
t433 = cos(t447) / 0.2e1;
t432 = sin(t446) / 0.2e1;
t430 = t441 * t477;
t429 = t435 - t444 / 0.2e1;
t427 = t443 / 0.2e1 + t478;
t425 = t434 + t441 / 0.2e1;
t424 = t438 / 0.2e1 + t479;
t423 = t433 - t437 / 0.2e1;
t422 = t433 + t437 / 0.2e1;
t421 = t432 - t436 / 0.2e1;
t420 = t432 + t436 / 0.2e1;
t418 = t431 + t471 / 0.2e1;
t417 = t430 - t472 / 0.2e1;
t416 = t430 + t472 / 0.2e1;
t415 = (t478 - t443 / 0.2e1) * qJD(4);
t414 = (t479 - t438 / 0.2e1) * qJD(4);
t409 = -t416 * t452 + t419 * t422;
t408 = t413 * qJD(1) + t465 * t418 - t461 * t470;
t407 = t412 * qJD(1) + t465 * t417 - t464 * t470;
t406 = -t411 * qJD(1) - t462 * t418 - t461 * t469;
t405 = -t410 * qJD(1) - t462 * t417 - t464 * t469;
t404 = -t407 * t454 + t458 * t467;
t403 = -t405 * t454 + t458 * t466;
t402 = t407 * t422 - t408 * t452 + t420 * t467;
t401 = t405 * t422 - t406 * t452 + t420 * t466;
t1 = [0, t466, 0, -t401 * t453 + t403 * t457 (t405 * t421 + t406 * t456 + t423 * t466) * t460 + (t412 * t421 + t413 * t456 + t423 * t475) * t468 - t401 * t427 - (t412 * t422 - t413 * t452 + t420 * t475) * t414 - t403 * t424 - (-t412 * t454 + t458 * t475) * t415, 0; 0, t467, 0, -t402 * t453 + t404 * t457 (t407 * t421 + t408 * t456 + t423 * t467) * t460 + (t410 * t421 + t411 * t456 - t423 * t474) * t468 - t402 * t427 - (t410 * t422 - t411 * t452 - t420 * t474) * t414 - t404 * t424 - (-t410 * t454 - t458 * t474) * t415, 0; 0, 0, 0, -t409 * t453 - t457 * t476 (t416 * t456 + t419 * t421) * t460 + (t425 * t421 + t459 * t423 + t429 * t456) * t468 - t409 * t427 - (t459 * t420 + t425 * t422 - t429 * t452) * t414 + t424 * t476 - (-t425 * t454 + t459 * t458) * t415, 0;];
JgD_rot  = t1;
