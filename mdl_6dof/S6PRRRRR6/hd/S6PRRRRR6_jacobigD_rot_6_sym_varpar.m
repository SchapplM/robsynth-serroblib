% Zeitableitung der rotatorischen Teilmatrix der geometrischen Jacobi-Matrix für Segment Nr. 6 (0=Basis) von
% S6PRRRRR6
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
%   pkin=[a2,a3,a4,a5,a6,alpha2,alpha3,alpha4,d2,d3,d4,d5,d6,theta1]';
%
% Output:
% JgD_rot [3x6]
%   Zeitableitung der rotatorischen Teilmatrix der geometrischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-02-26 20:22
% Revision: d75aae1ac561373cd3be920984c3df29a1c2ecc8 (2019-02-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JgD_rot = S6PRRRRR6_jacobigD_rot_6_sym_varpar(qJ, qJD, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(14,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRRRRR6_jacobigD_rot_6_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6PRRRRR6_jacobigD_rot_6_sym_varpar: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [14 1]), ...
  'S6PRRRRR6_jacobigD_rot_6_sym_varpar: pkin has to be [14x1] (double)');

%% Symbolic Calculation
% From jacobigD_rot_6_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 20:21:53
% EndTime: 2019-02-26 20:21:53
% DurationCPUTime: 0.36s
% Computational Cost: add. (276->74), mult. (919->150), div. (0->0), fcn. (1042->16), ass. (0->72)
t428 = sin(qJ(3));
t432 = cos(qJ(3));
t418 = sin(pkin(14));
t422 = cos(pkin(14));
t433 = cos(qJ(2));
t425 = cos(pkin(6));
t429 = sin(qJ(2));
t455 = t425 * t429;
t439 = t418 * t455 - t422 * t433;
t454 = t425 * t433;
t416 = -t418 * t454 - t422 * t429;
t424 = cos(pkin(7));
t420 = sin(pkin(7));
t421 = sin(pkin(6));
t462 = t420 * t421;
t440 = t416 * t424 + t418 * t462;
t404 = t440 * t428 - t432 * t439;
t415 = t418 * t433 + t422 * t455;
t414 = -t418 * t429 + t422 * t454;
t441 = -t414 * t424 + t422 * t462;
t466 = -t415 * t432 + t441 * t428;
t419 = sin(pkin(8));
t463 = t419 * t420;
t423 = cos(pkin(8));
t461 = t420 * t423;
t460 = t420 * t425;
t459 = t421 * t424;
t427 = sin(qJ(4));
t458 = t423 * t427;
t457 = t424 * t428;
t456 = t424 * t432;
t453 = t428 * t429;
t452 = t428 * t433;
t451 = t429 * t432;
t450 = t432 * t433;
t426 = sin(qJ(5));
t449 = qJD(4) * t426;
t448 = t427 * t463;
t447 = qJD(3) * t460;
t446 = qJD(2) * t429 * t462;
t445 = t419 * t446;
t401 = -t415 * t428 - t441 * t432;
t407 = -t414 * t420 - t422 * t459;
t444 = t401 * t423 + t407 * t419;
t403 = t428 * t439 + t440 * t432;
t408 = -t416 * t420 + t418 * t459;
t443 = t403 * t423 + t408 * t419;
t438 = t424 * t450 - t453;
t405 = t438 * t421 + t432 * t460;
t413 = t425 * t424 - t433 * t462;
t442 = t405 * t423 + t413 * t419;
t437 = t424 * t452 + t451;
t431 = cos(qJ(4));
t436 = t444 * t427 - t431 * t466;
t435 = t404 * t431 + t443 * t427;
t406 = t437 * t421 + t428 * t460;
t434 = t406 * t431 + t442 * t427;
t430 = cos(qJ(5));
t412 = t439 * qJD(2);
t411 = t416 * qJD(2);
t410 = t415 * qJD(2);
t409 = t414 * qJD(2);
t400 = t432 * t447 + (t438 * qJD(3) + (-t424 * t453 + t450) * qJD(2)) * t421;
t399 = -t428 * t447 + (-t437 * qJD(3) + (-t424 * t451 - t452) * qJD(2)) * t421;
t398 = -t399 * t419 + t423 * t446;
t397 = t403 * qJD(3) + t411 * t432 + t412 * t457;
t396 = -t404 * qJD(3) - t411 * t428 + t412 * t456;
t395 = t401 * qJD(3) + t409 * t432 - t410 * t457;
t394 = t466 * qJD(3) - t409 * t428 - t410 * t456;
t393 = -t396 * t419 - t412 * t461;
t392 = -t394 * t419 + t410 * t461;
t1 = [0, 0, -t412 * t420, t393, t397 * t427 + (-t396 * t423 + t412 * t463) * t431 + t435 * qJD(4) (t396 * t458 + t397 * t431 - t412 * t448) * t426 - t393 * t430 + (t435 * t430 + (-t403 * t419 + t408 * t423) * t426) * qJD(5) + (-t404 * t427 + t443 * t431) * t449; 0, 0, t410 * t420, t392, t395 * t427 + (-t394 * t423 - t410 * t463) * t431 + t436 * qJD(4) (t394 * t458 + t395 * t431 + t410 * t448) * t426 - t392 * t430 + (t436 * t430 + (-t401 * t419 + t407 * t423) * t426) * qJD(5) + (t427 * t466 + t444 * t431) * t449; 0, 0, t446, t398, t400 * t427 + (-t399 * t423 - t445) * t431 + t434 * qJD(4) (t399 * t458 + t400 * t431 + t427 * t445) * t426 - t398 * t430 + (t434 * t430 + (-t405 * t419 + t413 * t423) * t426) * qJD(5) + (-t406 * t427 + t442 * t431) * t449;];
JgD_rot  = t1;
