% Rotatorische Teilmatrix der geometrischen Jacobi-Matrix für Segment Nr. 6 (0=Basis) von
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
% pkin [14x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,alpha3,alpha4,d1,d2,d4,d5,d6,theta3]';
%
% Output:
% Jg_rot [3x6]
%   Rotatorische Teilmatrix der geometrischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox (ehem. IRT-Maple-Toolbox)
% Datum: 2018-12-10 18:38
% Revision: bb42a8b95257d9bc83910d26e849f5825122f662 (2018-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für mechatronische Systeme, Universität Hannover

function Jg_rot = S6RRPRRR14_jacobig_rot_6_floatb_twist_sym_varpar(qJ, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(14,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPRRR14_jacobig_rot_6_floatb_twist_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [14 1]), ...
  'S6RRPRRR14_jacobig_rot_6_floatb_twist_sym_varpar: pkin has to be [14x1] (double)');

%% Symbolic Calculation
% From jacobig_rot_6_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2018-12-10 18:38:23
% EndTime: 2018-12-10 18:38:24
% DurationCPUTime: 0.11s
% Computational Cost: add. (390->56), mult. (392->82), div. (0->0), fcn. (409->28), ass. (0->65)
t338 = sin(pkin(6));
t346 = sin(qJ(1));
t352 = t346 * t338;
t350 = cos(qJ(1));
t351 = t350 * t338;
t349 = cos(qJ(2));
t348 = cos(qJ(4));
t347 = cos(qJ(5));
t345 = sin(qJ(2));
t344 = sin(qJ(4));
t343 = sin(qJ(5));
t342 = cos(pkin(6));
t341 = cos(pkin(7));
t340 = cos(pkin(8));
t339 = cos(pkin(14));
t337 = sin(pkin(7));
t336 = sin(pkin(8));
t335 = sin(pkin(14));
t334 = pkin(6) - qJ(2);
t333 = pkin(6) + qJ(2);
t332 = pkin(8) - qJ(4);
t331 = pkin(8) + qJ(4);
t330 = pkin(7) - pkin(14);
t329 = pkin(7) + pkin(14);
t328 = cos(t333);
t327 = cos(t331);
t326 = sin(t334);
t325 = sin(t332);
t324 = cos(t329);
t323 = sin(t330);
t322 = cos(t334) / 0.2e1;
t321 = cos(t332) / 0.2e1;
t320 = sin(t333) / 0.2e1;
t319 = sin(t331) / 0.2e1;
t318 = cos(t330) / 0.2e1;
t317 = sin(t329) / 0.2e1;
t316 = t322 - t328 / 0.2e1;
t315 = t322 + t328 / 0.2e1;
t314 = t321 - t327 / 0.2e1;
t313 = t321 + t327 / 0.2e1;
t312 = t320 - t326 / 0.2e1;
t311 = t320 + t326 / 0.2e1;
t310 = t319 - t325 / 0.2e1;
t309 = t319 + t325 / 0.2e1;
t308 = t318 - t324 / 0.2e1;
t307 = t318 + t324 / 0.2e1;
t306 = t317 - t323 / 0.2e1;
t305 = t317 + t323 / 0.2e1;
t304 = -t346 * t312 + t350 * t349;
t303 = -t346 * t315 - t350 * t345;
t302 = t350 * t312 + t346 * t349;
t301 = t350 * t315 - t346 * t345;
t300 = -t311 * t337 + t342 * t341;
t299 = -t303 * t337 + t341 * t352;
t298 = -t301 * t337 - t341 * t351;
t297 = t311 * t306 + t342 * t308 + t316 * t339;
t296 = t342 * t305 + t311 * t307 - t316 * t335;
t295 = t303 * t306 + t304 * t339 + t308 * t352;
t294 = t303 * t307 - t304 * t335 + t305 * t352;
t293 = t301 * t306 + t302 * t339 - t308 * t351;
t292 = t301 * t307 - t302 * t335 - t305 * t351;
t291 = -t296 * t336 + t300 * t340;
t290 = -t294 * t336 + t299 * t340;
t289 = -t292 * t336 + t298 * t340;
t1 = [0, t352, 0, t290, -t294 * t313 + t295 * t344 - t299 * t309 (t294 * t310 + t295 * t348 + t299 * t314) * t343 - t290 * t347; 0, -t351, 0, t289, -t292 * t313 + t293 * t344 - t298 * t309 (t292 * t310 + t293 * t348 + t298 * t314) * t343 - t289 * t347; 1, t342, 0, t291, -t296 * t313 + t297 * t344 - t300 * t309 (t296 * t310 + t297 * t348 + t300 * t314) * t343 - t291 * t347;];
Jg_rot  = t1;
