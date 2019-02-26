% Zeitableitung der rotatorischen Teilmatrix der geometrischen Jacobi-Matrix für Segment Nr. 7 (0=Basis) von
% S7RRRRRRR1
% Use Code from Maple symbolic Code Generation
%
% geometrische Jacobi-Matrix: Differentieller Zusammenhang zwischen
% Endeffektorposition und verallgemeinerten Koordinaten.
%
%
% Input:
% qJ [7x1]
%   Generalized joint coordinates (joint angles)
% qJD [7x1]
%   Generalized joint velocities
% pkin [4x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[d1,d3,d5,d7]';
%
% Output:
% JgD_rot [3x7]
%   Zeitableitung der rotatorischen Teilmatrix der geometrischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-02-26 22:54
% Revision: d75aae1ac561373cd3be920984c3df29a1c2ecc8 (2019-02-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JgD_rot = S7RRRRRRR1_jacobigD_rot_7_sym_varpar(qJ, qJD, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(7,1),zeros(7,1),zeros(4,1)}
assert(isreal(qJ) && all(size(qJ) == [7 1]), ...
  'S7RRRRRRR1_jacobigD_rot_7_sym_varpar: qJ has to be [7x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [7 1]), ...
  'S7RRRRRRR1_jacobigD_rot_7_sym_varpar: qJD has to be [7x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [4 1]), ...
  'S7RRRRRRR1_jacobigD_rot_7_sym_varpar: pkin has to be [4x1] (double)');

%% Symbolic Calculation
% From jacobigD_rot_7_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 22:54:32
% EndTime: 2019-02-26 22:54:33
% DurationCPUTime: 0.26s
% Computational Cost: add. (174->63), mult. (520->120), div. (0->0), fcn. (557->12), ass. (0->60)
t400 = sin(qJ(2));
t401 = sin(qJ(1));
t406 = cos(qJ(2));
t422 = qJD(1) * t406 + qJD(3);
t407 = cos(qJ(1));
t427 = qJD(2) * t407;
t441 = t400 * t427 + t422 * t401;
t399 = sin(qJ(3));
t440 = t399 * t400;
t404 = cos(qJ(4));
t439 = t400 * t404;
t405 = cos(qJ(3));
t438 = t400 * t405;
t437 = t400 * t407;
t436 = t401 * t405;
t435 = t401 * t406;
t434 = t407 * t399;
t433 = t407 * t405;
t432 = qJD(1) * t401;
t431 = qJD(1) * t407;
t430 = qJD(2) * t400;
t429 = qJD(2) * t405;
t428 = qJD(2) * t406;
t426 = qJD(3) * t400;
t425 = qJD(4) * t400;
t423 = -qJD(3) * t406 - qJD(1);
t421 = qJD(4) * t405 - qJD(2);
t415 = t423 * t407;
t420 = t399 * t415 - t441 * t405 + t407 * t425;
t419 = t422 * t433 + (t423 * t399 - t400 * t429 + t425) * t401;
t398 = sin(qJ(4));
t414 = (-qJD(4) + t429) * t406;
t418 = -qJD(5) * t440 + t404 * t414 + (-qJD(3) * t399 * t404 - t421 * t398) * t400;
t393 = t405 * t435 + t434;
t389 = t401 * t400 * t398 + t393 * t404;
t392 = -t399 * t435 + t433;
t397 = sin(qJ(5));
t403 = cos(qJ(5));
t417 = t389 * t403 + t392 * t397;
t395 = -t401 * t399 + t406 * t433;
t390 = t395 * t404 + t398 * t437;
t394 = -t406 * t434 - t436;
t416 = t390 * t403 + t394 * t397;
t413 = -t400 * t431 - t401 * t428;
t412 = t400 * t432 - t406 * t427;
t411 = -t399 * t428 - t405 * t426;
t410 = -qJD(4) * t393 - t413;
t409 = qJD(4) * t395 + t412;
t391 = -t406 * t398 + t404 * t438;
t408 = qJD(5) * t391 - t411;
t402 = cos(qJ(6));
t396 = sin(qJ(6));
t387 = t423 * t436 + (t401 * t430 - t422 * t407) * t399;
t385 = t441 * t399 + t405 * t415;
t383 = t421 * t439 + (-t399 * t426 + t414) * t398;
t382 = t410 * t398 + t419 * t404;
t381 = t419 * t398 - t410 * t404;
t380 = -t409 * t398 + t420 * t404;
t379 = t420 * t398 + t409 * t404;
t1 = [0, t431, t412, t385, t379, t416 * qJD(5) + t380 * t397 - t385 * t403 -(t380 * t403 + t385 * t397 + (-t390 * t397 + t394 * t403) * qJD(5)) * t396 + t379 * t402 + (-t416 * t402 - (t395 * t398 - t404 * t437) * t396) * qJD(6); 0, t432, t413, t387, t381, t417 * qJD(5) + t382 * t397 - t387 * t403 -(t382 * t403 + t387 * t397 + (-t389 * t397 + t392 * t403) * qJD(5)) * t396 + t381 * t402 + (-t417 * t402 - (t393 * t398 - t401 * t439) * t396) * qJD(6); 0, 0, -t430, t411, t383, t418 * t397 + t408 * t403 (t383 - (t391 * t403 - t397 * t440) * qJD(6)) * t402 + (-t418 * t403 + t408 * t397 + (-t398 * t438 - t406 * t404) * qJD(6)) * t396;];
JgD_rot  = t1;
