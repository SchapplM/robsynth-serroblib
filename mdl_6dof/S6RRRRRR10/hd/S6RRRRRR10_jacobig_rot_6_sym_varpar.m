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

% Quelle: HybrDyn-Toolbox
% Datum: 2019-02-26 22:53
% Revision: d75aae1ac561373cd3be920984c3df29a1c2ecc8 (2019-02-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Jg_rot = S6RRRRRR10_jacobig_rot_6_sym_varpar(qJ, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(14,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRRRR10_jacobig_rot_6_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [14 1]), ...
  'S6RRRRRR10_jacobig_rot_6_sym_varpar: pkin has to be [14x1] (double)');

%% Symbolic Calculation
% From jacobig_rot_6_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 22:52:59
% EndTime: 2019-02-26 22:52:59
% DurationCPUTime: 0.12s
% Computational Cost: add. (107->33), mult. (311->69), div. (0->0), fcn. (432->16), ass. (0->46)
t308 = sin(pkin(7));
t312 = cos(pkin(6));
t335 = t308 * t312;
t311 = cos(pkin(7));
t321 = cos(qJ(2));
t334 = t311 * t321;
t309 = sin(pkin(6));
t317 = sin(qJ(1));
t333 = t317 * t309;
t316 = sin(qJ(2));
t332 = t317 * t316;
t331 = t317 * t321;
t322 = cos(qJ(1));
t330 = t322 * t309;
t329 = t322 * t316;
t328 = t322 * t321;
t304 = t312 * t329 + t331;
t315 = sin(qJ(3));
t320 = cos(qJ(3));
t303 = t312 * t328 - t332;
t324 = t303 * t311 - t308 * t330;
t294 = -t304 * t315 + t324 * t320;
t300 = -t303 * t308 - t311 * t330;
t307 = sin(pkin(8));
t310 = cos(pkin(8));
t327 = t294 * t310 + t300 * t307;
t306 = -t312 * t332 + t328;
t305 = -t312 * t331 - t329;
t323 = t305 * t311 + t308 * t333;
t296 = -t306 * t315 + t323 * t320;
t301 = -t305 * t308 + t311 * t333;
t326 = t296 * t310 + t301 * t307;
t298 = t320 * t335 + (-t315 * t316 + t320 * t334) * t309;
t302 = -t309 * t321 * t308 + t312 * t311;
t325 = t298 * t310 + t302 * t307;
t319 = cos(qJ(4));
t318 = cos(qJ(5));
t314 = sin(qJ(4));
t313 = sin(qJ(5));
t299 = t315 * t335 + (t315 * t334 + t316 * t320) * t309;
t297 = t306 * t320 + t323 * t315;
t295 = t304 * t320 + t324 * t315;
t293 = -t298 * t307 + t302 * t310;
t292 = -t296 * t307 + t301 * t310;
t291 = -t294 * t307 + t300 * t310;
t1 = [0, t333, t301, t292, t297 * t314 - t326 * t319 (t297 * t319 + t326 * t314) * t313 - t292 * t318; 0, -t330, t300, t291, t295 * t314 - t327 * t319 (t295 * t319 + t327 * t314) * t313 - t291 * t318; 1, t312, t302, t293, t299 * t314 - t325 * t319 (t299 * t319 + t325 * t314) * t313 - t293 * t318;];
Jg_rot  = t1;
