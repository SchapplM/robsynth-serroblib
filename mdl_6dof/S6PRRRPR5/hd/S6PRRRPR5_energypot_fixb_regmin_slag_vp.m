% Calculate minimal parameter regressor of potential energy for
% S6PRRRPR5
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [13x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,alpha3,d2,d3,d4,d6,theta1,theta5]';
% 
% Output:
% U_reg [1x27]
%   minimal parameter regressor of Potential energy

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-08 23:30
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U_reg = S6PRRRPR5_energypot_fixb_regmin_slag_vp(qJ, g, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(13,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRRRPR5_energypot_fixb_regmin_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6PRRRPR5_energypot_fixb_regmin_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [13 1]), ...
  'S6PRRRPR5_energypot_fixb_regmin_slag_vp: pkin has to be [13x1] (double)');

%% Symbolic Calculation
% From energy_potential_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-08 23:28:18
% EndTime: 2019-03-08 23:28:18
% DurationCPUTime: 0.18s
% Computational Cost: add. (221->75), mult. (550->127), div. (0->0), fcn. (707->16), ass. (0->52)
t294 = sin(pkin(12));
t297 = cos(pkin(12));
t304 = sin(qJ(2));
t299 = cos(pkin(6));
t308 = cos(qJ(2));
t311 = t299 * t308;
t283 = -t294 * t304 + t297 * t311;
t295 = sin(pkin(7));
t296 = sin(pkin(6));
t298 = cos(pkin(7));
t318 = t296 * t298;
t278 = -t283 * t295 - t297 * t318;
t302 = sin(qJ(4));
t324 = t278 * t302;
t285 = -t294 * t311 - t297 * t304;
t279 = -t285 * t295 + t294 * t318;
t323 = t279 * t302;
t315 = t296 * t308;
t282 = -t295 * t315 + t299 * t298;
t322 = t282 * t302;
t321 = t294 * t296;
t320 = t295 * t299;
t319 = t296 * t297;
t317 = t296 * t304;
t307 = cos(qJ(3));
t316 = t296 * t307;
t314 = t298 * t307;
t313 = t298 * t308;
t312 = t299 * t304;
t310 = t295 * t316;
t284 = t294 * t308 + t297 * t312;
t303 = sin(qJ(3));
t272 = -t283 * t314 + t284 * t303 + t297 * t310;
t286 = -t294 * t312 + t297 * t308;
t274 = -t285 * t314 + t286 * t303 - t294 * t310;
t276 = t303 * t317 - t307 * t320 - t313 * t316;
t309 = g(1) * t274 + g(2) * t272 + g(3) * t276;
t306 = cos(qJ(4));
t305 = cos(qJ(6));
t301 = sin(qJ(6));
t300 = -qJ(5) - pkin(10);
t293 = qJ(4) + pkin(13);
t292 = cos(t293);
t291 = sin(t293);
t290 = t306 * pkin(4) + pkin(3);
t277 = t303 * t320 + (t303 * t313 + t304 * t307) * t296;
t275 = t286 * t307 + (t285 * t298 + t295 * t321) * t303;
t273 = t284 * t307 + (t283 * t298 - t295 * t319) * t303;
t271 = t277 * t292 + t282 * t291;
t270 = t275 * t292 + t279 * t291;
t269 = t273 * t292 + t278 * t291;
t1 = [-g(3) * qJ(1), 0, -g(1) * t286 - g(2) * t284 - g(3) * t317, -g(1) * t285 - g(2) * t283 - g(3) * t315, 0, 0, 0, 0, 0, -g(1) * t275 - g(2) * t273 - g(3) * t277, t309, 0, 0, 0, 0, 0, -g(1) * (t275 * t306 + t323) - g(2) * (t273 * t306 + t324) - g(3) * (t277 * t306 + t322) -g(1) * (-t275 * t302 + t279 * t306) - g(2) * (-t273 * t302 + t278 * t306) - g(3) * (-t277 * t302 + t282 * t306) -t309, -g(1) * (t297 * pkin(1) + t286 * pkin(2) + pkin(4) * t323 + pkin(8) * t321 - t274 * t300 + t275 * t290) - g(2) * (t294 * pkin(1) + t284 * pkin(2) + pkin(4) * t324 - pkin(8) * t319 - t272 * t300 + t273 * t290) - g(3) * (pkin(2) * t317 + pkin(4) * t322 + t299 * pkin(8) - t276 * t300 + t277 * t290 + qJ(1)) + (-g(1) * t279 - g(2) * t278 - g(3) * t282) * pkin(9), 0, 0, 0, 0, 0, -g(1) * (t270 * t305 + t274 * t301) - g(2) * (t269 * t305 + t272 * t301) - g(3) * (t271 * t305 + t276 * t301) -g(1) * (-t270 * t301 + t274 * t305) - g(2) * (-t269 * t301 + t272 * t305) - g(3) * (-t271 * t301 + t276 * t305);];
U_reg  = t1;
