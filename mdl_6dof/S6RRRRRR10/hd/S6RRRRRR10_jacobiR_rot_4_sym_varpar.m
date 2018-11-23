% Rotatorische Teilmatrix der Rotationsmatrix-Jacobi-Matrix für Segment Nr. 4 (0=Basis) von
% S6RRRRRR10
% Use Code from Maple symbolic Code Generation
%
% Rotationsmatrix-Jacobi-Matrix: Differentieller Zusammenhang zwischen
% gestapelter Endeffektor-Rotationsmatrix und verallgemeinerten Koordinaten.
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
% JR_rot [9x6]
%   Jacobi-Matrix der Endeffektor-Rotationsmatrix

% Quelle: HybrDyn-Toolbox (ehem. IRT-Maple-Toolbox)
% Datum: 2018-11-23 11:28
% Revision: 76f9d5e39f14dc242b53c0d9d3d9db48bd8f37c0
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für mechatronische Systeme, Universität Hannover

function JR_rot = S6RRRRRR10_jacobiR_rot_4_sym_varpar(qJ, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(14,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRRRR10_jacobiR_rot_4_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [14 1]), ...
  'S6RRRRRR10_jacobiR_rot_4_sym_varpar: pkin has to be [14x1] (double)');

%% Symbolic Calculation
% From jacobiR_rot_4_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2018-11-23 11:27:18
% EndTime: 2018-11-23 11:27:18
% DurationCPUTime: 0.31s
% Computational Cost: add. (995->81), mult. (973->129), div. (0->0), fcn. (987->26), ass. (0->74)
t308 = sin(qJ(2));
t309 = sin(qJ(1));
t313 = cos(qJ(1));
t332 = pkin(6) + qJ(2);
t321 = cos(t332) / 0.2e1;
t333 = pkin(6) - qJ(2);
t327 = cos(t333);
t314 = t327 / 0.2e1 + t321;
t271 = t309 * t308 - t313 * t314;
t318 = sin(t332) / 0.2e1;
t324 = sin(t333);
t283 = t318 - t324 / 0.2e1;
t312 = cos(qJ(2));
t272 = t313 * t283 + t309 * t312;
t330 = pkin(7) + qJ(3);
t317 = sin(t330) / 0.2e1;
t331 = pkin(7) - qJ(3);
t323 = sin(t331);
t281 = t317 - t323 / 0.2e1;
t320 = cos(t330) / 0.2e1;
t326 = cos(t331);
t286 = t320 - t326 / 0.2e1;
t311 = cos(qJ(3));
t302 = sin(pkin(6));
t334 = t302 * t313;
t251 = -t271 * t281 + t272 * t311 + t286 * t334;
t280 = t317 + t323 / 0.2e1;
t287 = t326 / 0.2e1 + t320;
t307 = sin(qJ(3));
t253 = t271 * t287 + t272 * t307 + t280 * t334;
t301 = sin(pkin(7));
t304 = cos(pkin(7));
t267 = -t271 * t301 + t304 * t334;
t328 = pkin(8) + qJ(4);
t316 = sin(t328) / 0.2e1;
t329 = pkin(8) - qJ(4);
t322 = sin(t329);
t278 = t316 + t322 / 0.2e1;
t319 = cos(t328) / 0.2e1;
t325 = cos(t329);
t285 = t325 / 0.2e1 + t319;
t306 = sin(qJ(4));
t341 = t251 * t306 + t253 * t285 + t267 * t278;
t279 = t316 - t322 / 0.2e1;
t284 = t319 - t325 / 0.2e1;
t310 = cos(qJ(4));
t340 = -t251 * t310 + t253 * t279 - t267 * t284;
t339 = t278 * t301;
t338 = t284 * t301;
t288 = t321 - t327 / 0.2e1;
t337 = t288 * t301;
t303 = cos(pkin(8));
t336 = t301 * t303;
t335 = t302 * t309;
t276 = t309 * t283 - t313 * t312;
t274 = -t313 * t308 - t309 * t314;
t255 = t274 * t287 + t276 * t307 + t280 * t335;
t257 = -t274 * t281 + t276 * t311 + t286 * t335;
t269 = -t274 * t301 + t304 * t335;
t315 = t255 * t279 - t257 * t310 - t269 * t284;
t282 = t318 + t324 / 0.2e1;
t305 = cos(pkin(6));
t263 = t282 * t281 - t305 * t286 - t288 * t311;
t300 = sin(pkin(8));
t270 = -t282 * t301 + t305 * t304;
t266 = t288 * t281 + t282 * t311;
t265 = -t282 * t307 + t288 * t287;
t262 = t305 * t280 + t282 * t287 + t288 * t307;
t261 = t274 * t311 + t276 * t281;
t260 = -t274 * t307 + t276 * t287;
t259 = -t271 * t311 - t272 * t281;
t258 = t271 * t307 - t272 * t287;
t249 = t255 * t285 + t257 * t306 + t269 * t278;
t1 = [t340, t260 * t279 + t261 * t310 + t276 * t338, t255 * t310 + t257 * t279, t249, 0, 0; t315, t258 * t279 + t259 * t310 - t272 * t338, -t251 * t279 - t253 * t310, -t341, 0, 0; 0, t265 * t279 + t266 * t310 + t284 * t337, t262 * t310 - t263 * t279, t262 * t285 - t263 * t306 + t270 * t278, 0, 0; t341, t260 * t285 - t261 * t306 - t276 * t339, -t255 * t306 + t257 * t285, -t315, 0, 0; t249, t258 * t285 - t259 * t306 + t272 * t339, -t251 * t285 + t253 * t306, t340, 0, 0; 0, t265 * t285 - t266 * t306 - t278 * t337, -t262 * t306 - t263 * t285, -t262 * t279 - t263 * t310 + t270 * t284, 0, 0; -t253 * t300 + t267 * t303, -t260 * t300 - t276 * t336, -t257 * t300, 0, 0, 0; -t255 * t300 + t269 * t303, -t258 * t300 + t272 * t336, t251 * t300, 0, 0, 0; 0, -t265 * t300 - t288 * t336, t263 * t300, 0, 0, 0;];
JR_rot  = t1;
