% Zeitableitung der rotatorischen Teilmatrix der geometrischen Jacobi-Matrix für Segment Nr. 4 (0=Basis) von
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

% Quelle: HybrDyn-Toolbox (ehem. IRT-Maple-Toolbox)
% Datum: 2018-11-23 11:28
% Revision: 76f9d5e39f14dc242b53c0d9d3d9db48bd8f37c0
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für mechatronische Systeme, Universität Hannover

function JgD_rot = S6RRRRRR10_jacobigD_rot_4_floatb_twist_sym_varpar(qJ, qJD, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(14,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRRRR10_jacobigD_rot_4_floatb_twist_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRRRRR10_jacobigD_rot_4_floatb_twist_sym_varpar: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [14 1]), ...
  'S6RRRRRR10_jacobigD_rot_4_floatb_twist_sym_varpar: pkin has to be [14x1] (double)');

%% Symbolic Calculation
% From jacobigD_rot_4_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2018-11-23 11:27:18
% EndTime: 2018-11-23 11:27:19
% DurationCPUTime: 0.12s
% Computational Cost: add. (161->48), mult. (228->81), div. (0->0), fcn. (177->20), ass. (0->56)
t368 = pkin(7) - qJ(3);
t399 = sin(t368) / 0.2e1;
t367 = pkin(7) + qJ(3);
t398 = cos(t367) / 0.2e1;
t397 = qJD(2) / 0.2e1;
t369 = pkin(6) + qJ(2);
t365 = cos(t369);
t356 = t365 * t397;
t370 = pkin(6) - qJ(2);
t366 = cos(t370);
t391 = qJD(2) * t366;
t350 = t356 - t391 / 0.2e1;
t372 = sin(pkin(7));
t396 = t350 * t372;
t373 = sin(pkin(6));
t378 = sin(qJ(1));
t395 = t373 * t378;
t381 = cos(qJ(1));
t394 = t373 * t381;
t393 = qJD(1) * t373;
t361 = sin(t369);
t392 = qJD(2) * t361;
t390 = qJD(2) * t378;
t389 = qJD(2) * t381;
t388 = qJD(3) * cos(qJ(3));
t387 = t378 * t393;
t386 = t381 * t393;
t357 = t361 / 0.2e1;
t362 = sin(t370);
t352 = t357 - t362 / 0.2e1;
t380 = cos(qJ(2));
t385 = t381 * t352 + t378 * t380;
t384 = t378 * t352 - t381 * t380;
t358 = t366 / 0.2e1;
t354 = t358 + t365 / 0.2e1;
t377 = sin(qJ(2));
t383 = t381 * t354 - t378 * t377;
t382 = -t378 * t354 - t381 * t377;
t376 = sin(qJ(3));
t375 = cos(pkin(7));
t374 = cos(pkin(8));
t371 = sin(pkin(8));
t364 = cos(t368);
t359 = sin(t367);
t355 = t362 * t397;
t353 = t364 / 0.2e1 + t398;
t351 = t359 / 0.2e1 + t399;
t349 = t356 + t391 / 0.2e1;
t348 = t355 - t392 / 0.2e1;
t347 = (t398 - t364 / 0.2e1) * qJD(3);
t346 = (t399 - t359 / 0.2e1) * qJD(3);
t345 = t382 * qJD(1) + t381 * t348 - t380 * t390;
t344 = -t383 * qJD(1) - t378 * t348 - t380 * t389;
t343 = -t345 * t372 + t375 * t387;
t342 = -t344 * t372 + t375 * t386;
t1 = [0, t386, t342 -(-(-t378 * t349 - t377 * t389) * t376 + t384 * t388 + t344 * t353 + t382 * t346 + t347 * t395 + (t351 * t394 + t385 * t376) * qJD(1)) * t371 + t342 * t374, 0, 0; 0, t387, t343 -(-(t381 * t349 - t377 * t390) * t376 - t385 * t388 + t345 * t353 + t383 * t346 - t347 * t394 + (t351 * t395 + t384 * t376) * qJD(1)) * t371 + t343 * t374, 0, 0; 0, 0, -t396 -(-(t355 + t392 / 0.2e1) * t376 - (t358 - t365 / 0.2e1) * t388 + t350 * t353 + (t357 + t362 / 0.2e1) * t346 + cos(pkin(6)) * t347) * t371 - t374 * t396, 0, 0;];
JgD_rot  = t1;
