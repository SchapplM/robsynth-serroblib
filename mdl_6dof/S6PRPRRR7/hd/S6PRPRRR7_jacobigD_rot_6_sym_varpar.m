% Zeitableitung der rotatorischen Teilmatrix der geometrischen Jacobi-Matrix für Segment Nr. 6 (0=Basis) von
% S6PRPRRR7
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
%   pkin=[a2,a3,a4,a5,a6,alpha2,alpha3,alpha4,d2,d4,d5,d6,theta1,theta3]';
%
% Output:
% JgD_rot [3x6]
%   Zeitableitung der rotatorischen Teilmatrix der geometrischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-02-26 19:57
% Revision: d75aae1ac561373cd3be920984c3df29a1c2ecc8 (2019-02-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JgD_rot = S6PRPRRR7_jacobigD_rot_6_sym_varpar(qJ, qJD, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(14,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRPRRR7_jacobigD_rot_6_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6PRPRRR7_jacobigD_rot_6_sym_varpar: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [14 1]), ...
  'S6PRPRRR7_jacobigD_rot_6_sym_varpar: pkin has to be [14x1] (double)');

%% Symbolic Calculation
% From jacobigD_rot_6_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 19:57:32
% EndTime: 2019-02-26 19:57:32
% DurationCPUTime: 0.26s
% Computational Cost: add. (194->65), mult. (664->140), div. (0->0), fcn. (763->16), ass. (0->68)
t384 = sin(pkin(14));
t392 = cos(pkin(7));
t425 = t384 * t392;
t386 = sin(pkin(8));
t387 = sin(pkin(7));
t424 = t386 * t387;
t388 = sin(pkin(6));
t423 = t387 * t388;
t391 = cos(pkin(8));
t422 = t387 * t391;
t393 = cos(pkin(6));
t421 = t387 * t393;
t420 = t388 * t392;
t389 = cos(pkin(14));
t419 = t389 * t392;
t395 = sin(qJ(4));
t418 = t391 * t395;
t396 = sin(qJ(2));
t417 = t392 * t396;
t399 = cos(qJ(2));
t416 = t392 * t399;
t415 = t393 * t396;
t414 = t393 * t399;
t413 = qJD(2) * t388;
t394 = sin(qJ(5));
t412 = qJD(4) * t394;
t411 = t395 * t424;
t410 = t387 * t396 * t413;
t409 = t386 * t410;
t385 = sin(pkin(13));
t390 = cos(pkin(13));
t381 = t385 * t399 + t390 * t415;
t380 = -t385 * t396 + t390 * t414;
t405 = t380 * t392 - t390 * t423;
t360 = -t381 * t384 + t405 * t389;
t371 = -t380 * t387 - t390 * t420;
t408 = t360 * t391 + t371 * t386;
t403 = t385 * t415 - t390 * t399;
t382 = -t385 * t414 - t390 * t396;
t404 = t382 * t392 + t385 * t423;
t362 = t384 * t403 + t404 * t389;
t372 = -t382 * t387 + t385 * t420;
t407 = t362 * t391 + t372 * t386;
t369 = t389 * t421 + (-t384 * t396 + t389 * t416) * t388;
t379 = t393 * t392 - t399 * t423;
t406 = t369 * t391 + t379 * t386;
t361 = t381 * t389 + t405 * t384;
t398 = cos(qJ(4));
t402 = t361 * t398 + t408 * t395;
t363 = t404 * t384 - t389 * t403;
t401 = t363 * t398 + t407 * t395;
t370 = t388 * t396 * t389 + (t388 * t416 + t421) * t384;
t400 = t370 * t398 + t406 * t395;
t397 = cos(qJ(5));
t378 = t403 * qJD(2);
t377 = t382 * qJD(2);
t376 = t381 * qJD(2);
t375 = t380 * qJD(2);
t374 = (-t384 * t417 + t389 * t399) * t413;
t373 = (-t384 * t399 - t389 * t417) * t413;
t368 = -t373 * t386 + t391 * t410;
t367 = t377 * t389 + t378 * t425;
t366 = -t377 * t384 + t378 * t419;
t365 = t375 * t389 - t376 * t425;
t364 = -t375 * t384 - t376 * t419;
t359 = -t366 * t386 - t378 * t422;
t358 = -t364 * t386 + t376 * t422;
t1 = [0, 0, 0, t359, t367 * t395 + (-t366 * t391 + t378 * t424) * t398 + t401 * qJD(4) (t366 * t418 + t367 * t398 - t378 * t411) * t394 - t359 * t397 + (t401 * t397 + (-t362 * t386 + t372 * t391) * t394) * qJD(5) + (-t363 * t395 + t407 * t398) * t412; 0, 0, 0, t358, t365 * t395 + (-t364 * t391 - t376 * t424) * t398 + t402 * qJD(4) (t364 * t418 + t365 * t398 + t376 * t411) * t394 - t358 * t397 + (t402 * t397 + (-t360 * t386 + t371 * t391) * t394) * qJD(5) + (-t361 * t395 + t408 * t398) * t412; 0, 0, 0, t368, t374 * t395 + (-t373 * t391 - t409) * t398 + t400 * qJD(4) (t373 * t418 + t374 * t398 + t395 * t409) * t394 - t368 * t397 + (t400 * t397 + (-t369 * t386 + t379 * t391) * t394) * qJD(5) + (-t370 * t395 + t406 * t398) * t412;];
JgD_rot  = t1;
