% Zeitableitung der rotatorischen Teilmatrix der geometrischen Jacobi-Matrix für Segment Nr. 6 (0=Basis) von
% S6RPRRRR12
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
%   pkin=[a2,a3,a4,a5,a6,alpha2,alpha3,alpha4,d1,d3,d4,d5,d6,theta2]';
%
% Output:
% JgD_rot [3x6]
%   Zeitableitung der rotatorischen Teilmatrix der geometrischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-02-26 21:21
% Revision: d75aae1ac561373cd3be920984c3df29a1c2ecc8 (2019-02-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JgD_rot = S6RPRRRR12_jacobigD_rot_6_sym_varpar(qJ, qJD, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(14,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRRRR12_jacobigD_rot_6_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RPRRRR12_jacobigD_rot_6_sym_varpar: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [14 1]), ...
  'S6RPRRRR12_jacobigD_rot_6_sym_varpar: pkin has to be [14x1] (double)');

%% Symbolic Calculation
% From jacobigD_rot_6_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 21:21:11
% EndTime: 2019-02-26 21:21:11
% DurationCPUTime: 0.31s
% Computational Cost: add. (281->68), mult. (934->140), div. (0->0), fcn. (1057->16), ass. (0->69)
t416 = sin(pkin(14));
t419 = sin(pkin(6));
t426 = sin(qJ(3));
t430 = cos(qJ(3));
t420 = cos(pkin(14));
t422 = cos(pkin(7));
t453 = t420 * t422;
t418 = sin(pkin(7));
t423 = cos(pkin(6));
t456 = t418 * t423;
t404 = (t416 * t430 + t426 * t453) * t419 + t426 * t456;
t431 = cos(qJ(1));
t448 = t431 * t420;
t427 = sin(qJ(1));
t451 = t427 * t416;
t415 = -t423 * t451 + t448;
t449 = t431 * t416;
t450 = t427 * t420;
t414 = -t423 * t450 - t449;
t455 = t419 * t427;
t438 = t414 * t422 + t418 * t455;
t398 = t415 * t430 + t438 * t426;
t413 = t423 * t449 + t450;
t412 = t423 * t448 - t451;
t454 = t419 * t431;
t439 = -t412 * t422 + t418 * t454;
t461 = -t413 * t430 + t439 * t426;
t400 = t404 * qJD(3);
t417 = sin(pkin(8));
t460 = t400 * t417;
t425 = sin(qJ(4));
t457 = t417 * t425;
t421 = cos(pkin(8));
t452 = t421 * t425;
t447 = qJD(1) * t419;
t424 = sin(qJ(5));
t446 = qJD(4) * t424;
t444 = t427 * t447;
t443 = t431 * t447;
t395 = -t413 * t426 - t439 * t430;
t405 = -t412 * t418 - t422 * t454;
t442 = t395 * t421 + t405 * t417;
t397 = -t415 * t426 + t438 * t430;
t406 = -t414 * t418 + t422 * t455;
t441 = t397 * t421 + t406 * t417;
t403 = t430 * t456 + (-t416 * t426 + t430 * t453) * t419;
t411 = -t419 * t420 * t418 + t423 * t422;
t440 = t403 * t421 + t411 * t417;
t407 = t412 * qJD(1);
t436 = -t407 * t422 + t418 * t443;
t409 = t414 * qJD(1);
t435 = t409 * t422 + t418 * t444;
t429 = cos(qJ(4));
t434 = t442 * t425 - t429 * t461;
t433 = t398 * t429 + t441 * t425;
t432 = t404 * t429 + t440 * t425;
t428 = cos(qJ(5));
t410 = t415 * qJD(1);
t408 = t413 * qJD(1);
t402 = -t409 * t418 + t422 * t444;
t401 = t407 * t418 + t422 * t443;
t399 = t403 * qJD(3);
t394 = t395 * qJD(3) + t410 * t430 + t435 * t426;
t393 = t461 * qJD(3) - t410 * t426 + t435 * t430;
t392 = t397 * qJD(3) - t408 * t430 + t436 * t426;
t391 = -t398 * qJD(3) + t408 * t426 + t436 * t430;
t390 = -t393 * t417 + t402 * t421;
t389 = -t391 * t417 + t401 * t421;
t1 = [0, 0, t401, t389, t392 * t425 + (-t391 * t421 - t401 * t417) * t429 + t433 * qJD(4) (t391 * t452 + t392 * t429 + t401 * t457) * t424 - t389 * t428 + (t433 * t428 + (-t397 * t417 + t406 * t421) * t424) * qJD(5) + (-t398 * t425 + t441 * t429) * t446; 0, 0, t402, t390, t394 * t425 + (-t393 * t421 - t402 * t417) * t429 + t434 * qJD(4) (t393 * t452 + t394 * t429 + t402 * t457) * t424 - t390 * t428 + (t434 * t428 + (-t395 * t417 + t405 * t421) * t424) * qJD(5) + (t425 * t461 + t442 * t429) * t446; 0, 0, 0, t460, t400 * t421 * t429 + t432 * qJD(4) + t399 * t425 (t399 * t429 - t400 * t452) * t424 - t428 * t460 + (t432 * t428 + (-t403 * t417 + t411 * t421) * t424) * qJD(5) + (-t404 * t425 + t440 * t429) * t446;];
JgD_rot  = t1;
