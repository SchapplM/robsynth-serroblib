% Calculate minimal parameter regressor of potential energy for
% S6PRRPRR3
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [13x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,alpha3,d2,d3,d5,d6,theta1,theta4]';
% 
% Output:
% U_reg [1x27]
%   minimal parameter regressor of Potential energy

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-08 22:09
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U_reg = S6PRRPRR3_energypot_fixb_regmin_slag_vp(qJ, g, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(13,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRRPRR3_energypot_fixb_regmin_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6PRRPRR3_energypot_fixb_regmin_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [13 1]), ...
  'S6PRRPRR3_energypot_fixb_regmin_slag_vp: pkin has to be [13x1] (double)');

%% Symbolic Calculation
% From energy_potential_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-08 22:07:42
% EndTime: 2019-03-08 22:07:42
% DurationCPUTime: 0.22s
% Computational Cost: add. (233->72), mult. (607->131), div. (0->0), fcn. (788->16), ass. (0->52)
t297 = sin(pkin(7));
t301 = cos(pkin(7));
t320 = pkin(9) + qJ(4);
t305 = sin(qJ(3));
t322 = pkin(3) * t305;
t323 = t297 * t322 + t320 * t301 + pkin(8);
t296 = sin(pkin(12));
t298 = sin(pkin(6));
t319 = t296 * t298;
t300 = cos(pkin(12));
t318 = t298 * t300;
t317 = t298 * t301;
t310 = cos(qJ(2));
t316 = t298 * t310;
t302 = cos(pkin(6));
t306 = sin(qJ(2));
t315 = t302 * t306;
t314 = t302 * t310;
t295 = sin(pkin(13));
t299 = cos(pkin(13));
t309 = cos(qJ(3));
t313 = t309 * t295 + t305 * t299;
t293 = -t305 * t295 + t309 * t299;
t289 = t296 * t310 + t300 * t315;
t291 = -t296 * t315 + t300 * t310;
t312 = g(3) * t298 * t306 + g(1) * t291 + g(2) * t289;
t288 = -t296 * t306 + t300 * t314;
t290 = -t296 * t314 - t300 * t306;
t311 = -g(1) * (t290 * t301 + t297 * t319) - g(2) * (t288 * t301 - t297 * t318) - g(3) * (t297 * t302 + t301 * t316);
t308 = cos(qJ(5));
t307 = cos(qJ(6));
t304 = sin(qJ(5));
t303 = sin(qJ(6));
t294 = t309 * pkin(3) + pkin(2);
t287 = -t297 * t316 + t302 * t301;
t286 = -t320 * t297 + t301 * t322;
t284 = t313 * t301;
t283 = t293 * t301;
t282 = t313 * t297;
t281 = t293 * t297;
t280 = -t290 * t297 + t296 * t317;
t279 = -t288 * t297 - t300 * t317;
t278 = t302 * t282 + (t284 * t310 + t293 * t306) * t298;
t277 = -t302 * t281 + (-t283 * t310 + t306 * t313) * t298;
t276 = t278 * t308 + t287 * t304;
t275 = t282 * t319 + t290 * t284 + t291 * t293;
t274 = -t281 * t319 - t290 * t283 + t291 * t313;
t273 = -t282 * t318 + t288 * t284 + t289 * t293;
t272 = t281 * t318 - t288 * t283 + t289 * t313;
t271 = t275 * t308 + t280 * t304;
t270 = t273 * t308 + t279 * t304;
t1 = [-g(3) * qJ(1), 0, -t312, -g(1) * t290 - g(2) * t288 - g(3) * t316, 0, 0, 0, 0, 0, t311 * t305 - t312 * t309, t312 * t305 + t311 * t309, -g(1) * t280 - g(2) * t279 - g(3) * t287, -g(1) * (t300 * pkin(1) + t290 * t286 + t291 * t294) - g(2) * (t296 * pkin(1) + t288 * t286 + t289 * t294) - g(3) * (t323 * t302 + qJ(1)) + (-g(3) * (t286 * t310 + t294 * t306) + (-g(1) * t296 + g(2) * t300) * t323) * t298, 0, 0, 0, 0, 0, -g(1) * t271 - g(2) * t270 - g(3) * t276, -g(1) * (-t275 * t304 + t280 * t308) - g(2) * (-t273 * t304 + t279 * t308) - g(3) * (-t278 * t304 + t287 * t308) 0, 0, 0, 0, 0, -g(1) * (t271 * t307 + t274 * t303) - g(2) * (t270 * t307 + t272 * t303) - g(3) * (t276 * t307 + t277 * t303) -g(1) * (-t271 * t303 + t274 * t307) - g(2) * (-t270 * t303 + t272 * t307) - g(3) * (-t276 * t303 + t277 * t307);];
U_reg  = t1;
