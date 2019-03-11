% Calculate minimal parameter regressor of potential energy for
% S6RRRPPP1
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha4,d1,d2,d3,theta4]';
% 
% Output:
% U_reg [1x29]
%   minimal parameter regressor of Potential energy

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 15:20
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U_reg = S6RRRPPP1_energypot_fixb_regmin_slag_vp(qJ, g, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(10,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRPPP1_energypot_fixb_regmin_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RRRPPP1_energypot_fixb_regmin_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RRRPPP1_energypot_fixb_regmin_slag_vp: pkin has to be [10x1] (double)');

%% Symbolic Calculation
% From energy_potential_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 15:18:28
% EndTime: 2019-03-09 15:18:28
% DurationCPUTime: 0.18s
% Computational Cost: add. (249->71), mult. (624->105), div. (0->0), fcn. (761->10), ass. (0->46)
t308 = sin(qJ(1));
t309 = cos(qJ(3));
t310 = cos(qJ(2));
t306 = sin(qJ(3));
t311 = cos(qJ(1));
t328 = t311 * t306;
t285 = t308 * t309 - t310 * t328;
t304 = sin(pkin(6));
t305 = cos(pkin(6));
t307 = sin(qJ(2));
t330 = t307 * t311;
t277 = -t285 * t304 + t305 * t330;
t327 = t311 * t309;
t329 = t308 * t310;
t283 = -t306 * t329 - t327;
t332 = t307 * t308;
t276 = -t283 * t304 + t305 * t332;
t337 = cos(pkin(10));
t334 = t304 * t310;
t333 = t306 * t307;
t331 = t307 * t309;
t326 = t304 * t333;
t323 = t305 * t337;
t322 = t307 * t337;
t321 = t304 * t322;
t320 = g(1) * t311 + g(2) * t308;
t284 = t309 * t329 - t328;
t303 = sin(pkin(10));
t269 = -t283 * t323 + t284 * t303 - t308 * t321;
t286 = t308 * t306 + t310 * t327;
t271 = -t285 * t323 + t286 * t303 - t311 * t321;
t274 = t337 * t334 + (t303 * t309 + t306 * t323) * t307;
t319 = g(1) * t271 + g(2) * t269 + g(3) * t274;
t270 = t284 * t337 + (t283 * t305 + t304 * t332) * t303;
t272 = t286 * t337 + (t285 * t305 + t304 * t330) * t303;
t275 = t309 * t322 + (-t305 * t333 - t334) * t303;
t318 = g(1) * t272 + g(2) * t270 + g(3) * t275;
t317 = pkin(7) + qJ(4) * t326 + pkin(3) * t331 + t307 * pkin(2) + (-qJ(4) * t305 - pkin(9)) * t310;
t316 = t286 * pkin(3) + t308 * pkin(8) + pkin(9) * t330 + (pkin(2) * t310 + pkin(1)) * t311 + t277 * qJ(4);
t315 = t308 * pkin(1) + pkin(2) * t329 + t284 * pkin(3) - t311 * pkin(8) + pkin(9) * t332 + t276 * qJ(4);
t314 = t275 * pkin(4) + t274 * qJ(5) + t317;
t313 = t272 * pkin(4) + t271 * qJ(5) + t316;
t312 = t270 * pkin(4) + t269 * qJ(5) + t315;
t282 = -t310 * t305 + t326;
t266 = -g(1) * t277 - g(2) * t276 - g(3) * t282;
t1 = [0, -t320, g(1) * t308 - g(2) * t311, 0, 0, 0, 0, 0, -g(3) * t307 - t320 * t310, -g(3) * t310 + t320 * t307, 0, 0, 0, 0, 0, -g(1) * t286 - g(2) * t284 - g(3) * t331, -g(1) * t285 - g(2) * t283 + g(3) * t333, -t318, t319, t266, -g(1) * t316 - g(2) * t315 - g(3) * t317, t266, t318, -t319, -g(1) * t313 - g(2) * t312 - g(3) * t314, t266, -t319, -t318, -g(1) * (t277 * pkin(5) + t272 * qJ(6) + t313) - g(2) * (t276 * pkin(5) + t270 * qJ(6) + t312) - g(3) * (t282 * pkin(5) + t275 * qJ(6) + t314);];
U_reg  = t1;
