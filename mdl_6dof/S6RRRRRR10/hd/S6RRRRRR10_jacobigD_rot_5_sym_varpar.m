% Zeitableitung der rotatorischen Teilmatrix der geometrischen Jacobi-Matrix für Segment Nr. 5 (0=Basis) von
% S6RRRRRR10
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
%   pkin=[a2,a3,a4,a5,a6,alpha2,alpha3,alpha4,d1,d2,d3,d4,d5,d6]';
%
% Output:
% JgD_rot [3x6]
%   Zeitableitung der rotatorischen Teilmatrix der geometrischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-02-26 22:53
% Revision: d75aae1ac561373cd3be920984c3df29a1c2ecc8 (2019-02-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JgD_rot = S6RRRRRR10_jacobigD_rot_5_sym_varpar(qJ, qJD, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(14,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRRRR10_jacobigD_rot_5_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRRRRR10_jacobigD_rot_5_sym_varpar: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [14 1]), ...
  'S6RRRRRR10_jacobigD_rot_5_sym_varpar: pkin has to be [14x1] (double)');

%% Symbolic Calculation
% From jacobigD_rot_5_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 22:52:48
% EndTime: 2019-02-26 22:52:49
% DurationCPUTime: 0.31s
% Computational Cost: add. (157->57), mult. (522->120), div. (0->0), fcn. (559->14), ass. (0->56)
t387 = sin(qJ(3));
t391 = cos(qJ(3));
t385 = cos(pkin(6));
t392 = cos(qJ(2));
t393 = cos(qJ(1));
t409 = t393 * t392;
t388 = sin(qJ(2));
t389 = sin(qJ(1));
t413 = t389 * t388;
t398 = t385 * t413 - t409;
t410 = t393 * t388;
t412 = t389 * t392;
t378 = -t385 * t412 - t410;
t381 = sin(pkin(7));
t384 = cos(pkin(7));
t382 = sin(pkin(6));
t418 = t382 * t389;
t401 = t378 * t384 + t381 * t418;
t424 = t401 * t387 - t391 * t398;
t377 = t385 * t410 + t412;
t376 = t385 * t409 - t413;
t417 = t382 * t393;
t402 = -t376 * t384 + t381 * t417;
t423 = -t377 * t391 + t402 * t387;
t420 = t381 * t385;
t419 = t381 * t388;
t416 = t387 * t388;
t415 = t387 * t392;
t414 = t388 * t391;
t411 = t391 * t392;
t408 = qJD(1) * t382;
t386 = sin(qJ(4));
t407 = qJD(3) * t386;
t406 = t389 * t408;
t405 = t393 * t408;
t404 = qJD(3) * t420;
t403 = t382 * qJD(2) * t419;
t400 = t384 * t411 - t416;
t399 = t384 * t415 + t414;
t372 = -t376 * qJD(1) + t398 * qJD(2);
t397 = t372 * t384 + t381 * t405;
t374 = t378 * qJD(1) - t377 * qJD(2);
t396 = t374 * t384 + t381 * t406;
t395 = -t377 * t387 - t402 * t391;
t394 = t387 * t398 + t401 * t391;
t390 = cos(qJ(4));
t383 = cos(pkin(8));
t380 = sin(pkin(8));
t375 = -t398 * qJD(1) + t376 * qJD(2);
t373 = -t377 * qJD(1) + t378 * qJD(2);
t371 = -t374 * t381 + t384 * t406;
t370 = -t372 * t381 + t384 * t405;
t369 = -t387 * t404 + (-t399 * qJD(3) + (-t384 * t414 - t415) * qJD(2)) * t382;
t368 = t423 * qJD(3) - t375 * t387 + t396 * t391;
t367 = -t424 * qJD(3) - t373 * t387 + t397 * t391;
t1 = [0, t405, t370, -t367 * t380 + t370 * t383 (t373 * t391 + t397 * t387) * t386 + (-t367 * t383 - t370 * t380) * t390 + t394 * t407 + (t424 * t390 + (t394 * t383 + (-t378 * t381 + t384 * t418) * t380) * t386) * qJD(4), 0; 0, t406, t371, -t368 * t380 + t371 * t383 (t375 * t391 + t396 * t387) * t386 + (-t368 * t383 - t371 * t380) * t390 + t395 * t407 + (-t423 * t390 + (t395 * t383 + (-t376 * t381 - t384 * t417) * t380) * t386) * qJD(4), 0; 0, 0, t403, -t369 * t380 + t383 * t403, t391 * t386 * t404 - t369 * t383 * t390 + (t400 * t407 + ((-t384 * t416 + t411) * t386 - t380 * t390 * t419) * qJD(2)) * t382 + ((t399 * t382 + t387 * t420) * t390 + ((t400 * t382 + t391 * t420) * t383 + (-t382 * t392 * t381 + t385 * t384) * t380) * t386) * qJD(4), 0;];
JgD_rot  = t1;
