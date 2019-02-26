% Rotatorische Teilmatrix der geometrischen Jacobi-Matrix für Segment Nr. 6 (0=Basis) von
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
% pkin [14x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,alpha3,alpha4,d1,d2,d3,d4,d5,d6]';
%
% Output:
% Jg_rot [3x6]
%   Rotatorische Teilmatrix der geometrischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox (ehem. IRT-Maple-Toolbox)
% Datum: 2018-11-23 11:27
% Revision: 76f9d5e39f14dc242b53c0d9d3d9db48bd8f37c0
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für mechatronische Systeme, Universität Hannover

function Jg_rot = S6RRRRRR10_jacobig_rot_6_floatb_twist_sym_varpar(qJ, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(14,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRRRR10_jacobig_rot_6_floatb_twist_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [14 1]), ...
  'S6RRRRRR10_jacobig_rot_6_floatb_twist_sym_varpar: pkin has to be [14x1] (double)');

%% Symbolic Calculation
% From jacobig_rot_6_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2018-11-23 11:27:20
% EndTime: 2018-11-23 11:27:21
% DurationCPUTime: 0.16s
% Computational Cost: add. (404->56), mult. (410->82), div. (0->0), fcn. (432->28), ass. (0->65)
t373 = sin(pkin(6));
t381 = sin(qJ(1));
t388 = t381 * t373;
t386 = cos(qJ(1));
t387 = t386 * t373;
t385 = cos(qJ(2));
t384 = cos(qJ(3));
t383 = cos(qJ(4));
t382 = cos(qJ(5));
t380 = sin(qJ(2));
t379 = sin(qJ(3));
t378 = sin(qJ(4));
t377 = sin(qJ(5));
t376 = cos(pkin(6));
t375 = cos(pkin(7));
t374 = cos(pkin(8));
t372 = sin(pkin(7));
t371 = sin(pkin(8));
t370 = pkin(6) - qJ(2);
t369 = pkin(6) + qJ(2);
t368 = pkin(7) - qJ(3);
t367 = pkin(7) + qJ(3);
t366 = pkin(8) - qJ(4);
t365 = pkin(8) + qJ(4);
t364 = cos(t369);
t363 = cos(t367);
t362 = cos(t365);
t361 = sin(t370);
t360 = sin(t368);
t359 = sin(t366);
t358 = cos(t370) / 0.2e1;
t357 = cos(t368) / 0.2e1;
t356 = cos(t366) / 0.2e1;
t355 = sin(t369) / 0.2e1;
t354 = sin(t367) / 0.2e1;
t353 = sin(t365) / 0.2e1;
t352 = t358 - t364 / 0.2e1;
t351 = t358 + t364 / 0.2e1;
t350 = t357 - t363 / 0.2e1;
t349 = t357 + t363 / 0.2e1;
t348 = t356 - t362 / 0.2e1;
t347 = t356 + t362 / 0.2e1;
t346 = t355 - t361 / 0.2e1;
t345 = t355 + t361 / 0.2e1;
t344 = t354 - t360 / 0.2e1;
t343 = t354 + t360 / 0.2e1;
t342 = t353 - t359 / 0.2e1;
t341 = t353 + t359 / 0.2e1;
t340 = -t346 * t381 + t385 * t386;
t339 = -t351 * t381 - t380 * t386;
t338 = t346 * t386 + t381 * t385;
t337 = t351 * t386 - t380 * t381;
t336 = -t345 * t372 + t375 * t376;
t335 = -t339 * t372 + t375 * t388;
t334 = -t337 * t372 - t375 * t387;
t333 = t344 * t345 + t350 * t376 + t352 * t384;
t332 = t343 * t376 + t345 * t349 - t352 * t379;
t331 = t339 * t344 + t340 * t384 + t350 * t388;
t330 = t339 * t349 - t340 * t379 + t343 * t388;
t329 = t337 * t344 + t338 * t384 - t350 * t387;
t328 = t337 * t349 - t338 * t379 - t343 * t387;
t327 = -t332 * t371 + t336 * t374;
t326 = -t330 * t371 + t335 * t374;
t325 = -t328 * t371 + t334 * t374;
t1 = [0, t388, t335, t326, -t330 * t347 + t331 * t378 - t335 * t341 (t330 * t342 + t331 * t383 + t335 * t348) * t377 - t326 * t382; 0, -t387, t334, t325, -t328 * t347 + t329 * t378 - t334 * t341 (t328 * t342 + t329 * t383 + t334 * t348) * t377 - t325 * t382; 1, t376, t336, t327, -t332 * t347 + t333 * t378 - t336 * t341 (t332 * t342 + t333 * t383 + t336 * t348) * t377 - t327 * t382;];
Jg_rot  = t1;
