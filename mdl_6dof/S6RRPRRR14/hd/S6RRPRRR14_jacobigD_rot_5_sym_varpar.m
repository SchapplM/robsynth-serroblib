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

% Quelle: HybrDyn-Toolbox
% Datum: 2019-02-26 22:56
% Revision: d75aae1ac561373cd3be920984c3df29a1c2ecc8 (2019-02-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JgD_rot = S6RRPRRR14_jacobigD_rot_5_sym_varpar(qJ, qJD, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(14,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPRRR14_jacobigD_rot_5_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRPRRR14_jacobigD_rot_5_sym_varpar: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [14 1]), ...
  'S6RRPRRR14_jacobigD_rot_5_sym_varpar: pkin has to be [14x1] (double)');

%% Symbolic Calculation
% From jacobigD_rot_5_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 22:55:47
% EndTime: 2019-02-26 22:55:47
% DurationCPUTime: 0.27s
% Computational Cost: add. (110->49), mult. (371->110), div. (0->0), fcn. (396->14), ass. (0->46)
t354 = sin(pkin(7));
t359 = cos(pkin(6));
t384 = t354 * t359;
t361 = sin(qJ(2));
t383 = t354 * t361;
t355 = sin(pkin(6));
t362 = sin(qJ(1));
t382 = t355 * t362;
t365 = cos(qJ(1));
t381 = t355 * t365;
t358 = cos(pkin(7));
t380 = t358 * t361;
t364 = cos(qJ(2));
t379 = t358 * t364;
t378 = t361 * t365;
t377 = t362 * t361;
t376 = t362 * t364;
t375 = t364 * t365;
t374 = qJD(1) * t355;
t373 = qJD(2) * t355;
t372 = t362 * t374;
t371 = t365 * t374;
t348 = t359 * t375 - t377;
t370 = t348 * t358 - t354 * t381;
t350 = -t359 * t376 - t378;
t369 = t350 * t358 + t354 * t382;
t349 = t359 * t378 + t376;
t368 = t359 * t377 - t375;
t343 = -t348 * qJD(1) + t368 * qJD(2);
t367 = t343 * t358 + t354 * t371;
t345 = t350 * qJD(1) - t349 * qJD(2);
t366 = t345 * t358 + t354 * t372;
t363 = cos(qJ(4));
t360 = sin(qJ(4));
t357 = cos(pkin(8));
t356 = cos(pkin(14));
t353 = sin(pkin(8));
t352 = sin(pkin(14));
t347 = (-t352 * t364 - t356 * t380) * t373;
t346 = -t368 * qJD(1) + t348 * qJD(2);
t344 = -t349 * qJD(1) + t350 * qJD(2);
t342 = -t345 * t354 + t358 * t372;
t341 = -t343 * t354 + t358 * t371;
t340 = -t346 * t352 + t366 * t356;
t339 = -t344 * t352 + t367 * t356;
t1 = [0, t371, 0, -t339 * t353 + t341 * t357 (t344 * t356 + t367 * t352) * t360 + (-t339 * t357 - t341 * t353) * t363 + ((t369 * t352 - t356 * t368) * t363 + ((t352 * t368 + t369 * t356) * t357 + (-t350 * t354 + t358 * t382) * t353) * t360) * qJD(4), 0; 0, t372, 0, -t340 * t353 + t342 * t357 (t346 * t356 + t366 * t352) * t360 + (-t340 * t357 - t342 * t353) * t363 + ((t349 * t356 + t370 * t352) * t363 + ((-t349 * t352 + t370 * t356) * t357 + (-t348 * t354 - t358 * t381) * t353) * t360) * qJD(4), 0; 0, 0, 0, t357 * t373 * t383 - t347 * t353, -t347 * t357 * t363 + ((t355 * t361 * t356 + (t355 * t379 + t384) * t352) * t363 + ((t356 * t384 + (-t352 * t361 + t356 * t379) * t355) * t357 + (-t354 * t355 * t364 + t358 * t359) * t353) * t360) * qJD(4) + ((-t352 * t380 + t356 * t364) * t360 - t353 * t363 * t383) * t373, 0;];
JgD_rot  = t1;
