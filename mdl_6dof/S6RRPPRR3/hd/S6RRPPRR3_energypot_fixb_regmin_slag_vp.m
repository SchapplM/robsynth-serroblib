% Calculate minimal parameter regressor of potential energy for
% S6RRPPRR3
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [12x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d1,d2,d5,d6,theta3,theta4]';
% 
% Output:
% U_reg [1x30]
%   minimal parameter regressor of Potential energy

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 09:00
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U_reg = S6RRPPRR3_energypot_fixb_regmin_slag_vp(qJ, g, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(12,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPPRR3_energypot_fixb_regmin_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RRPPRR3_energypot_fixb_regmin_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [12 1]), ...
  'S6RRPPRR3_energypot_fixb_regmin_slag_vp: pkin has to be [12x1] (double)');

%% Symbolic Calculation
% From energy_potential_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 08:58:37
% EndTime: 2019-03-09 08:58:37
% DurationCPUTime: 0.20s
% Computational Cost: add. (183->67), mult. (378->119), div. (0->0), fcn. (477->14), ass. (0->43)
t309 = pkin(8) + qJ(3);
t285 = sin(pkin(6));
t290 = sin(qJ(2));
t308 = t285 * t290;
t291 = sin(qJ(1));
t307 = t285 * t291;
t294 = cos(qJ(1));
t306 = t285 * t294;
t305 = t291 * t290;
t293 = cos(qJ(2));
t304 = t291 * t293;
t303 = t294 * t290;
t302 = t294 * t293;
t288 = cos(pkin(6));
t271 = t288 * t290 * pkin(2) - t309 * t285;
t277 = t293 * pkin(2) + pkin(1);
t301 = t294 * t271 + t291 * t277;
t300 = -t291 * t271 + t294 * t277;
t299 = pkin(2) * t308 + t309 * t288 + pkin(7);
t298 = g(1) * t291 - g(2) * t294;
t284 = sin(pkin(11));
t287 = cos(pkin(11));
t297 = t293 * t284 + t290 * t287;
t296 = t290 * t284 - t293 * t287;
t295 = t296 * t288;
t292 = cos(qJ(6));
t289 = sin(qJ(6));
t286 = cos(pkin(12));
t283 = sin(pkin(12));
t282 = pkin(12) + qJ(5);
t279 = cos(t282);
t278 = sin(t282);
t270 = t297 * t288;
t269 = t297 * t285;
t268 = t296 * t285;
t265 = t269 * t279 + t288 * t278;
t264 = -t291 * t270 - t294 * t296;
t263 = t291 * t295 - t294 * t297;
t262 = t294 * t270 - t291 * t296;
t261 = -t291 * t297 - t294 * t295;
t260 = t264 * t279 + t278 * t307;
t259 = t262 * t279 - t278 * t306;
t1 = [0, -g(1) * t294 - g(2) * t291, t298, 0, 0, 0, 0, 0, -g(1) * (-t288 * t305 + t302) - g(2) * (t288 * t303 + t304) - g(3) * t308, -g(1) * (-t288 * t304 - t303) - g(2) * (t288 * t302 - t305) - g(3) * t285 * t293, -g(3) * t288 - t298 * t285, -g(1) * t300 - g(2) * t301 - g(3) * t299, -g(1) * (t264 * t286 + t283 * t307) - g(2) * (t262 * t286 - t283 * t306) - g(3) * (t269 * t286 + t288 * t283) -g(1) * (-t264 * t283 + t286 * t307) - g(2) * (-t262 * t283 - t286 * t306) - g(3) * (-t269 * t283 + t288 * t286) g(1) * t263 + g(2) * t261 - g(3) * t268, -g(1) * (t264 * pkin(3) - t263 * qJ(4) + t300) - g(2) * (t262 * pkin(3) - t261 * qJ(4) + t301) - g(3) * (t269 * pkin(3) + t268 * qJ(4) + t299) 0, 0, 0, 0, 0, -g(1) * t260 - g(2) * t259 - g(3) * t265, -g(1) * (-t264 * t278 + t279 * t307) - g(2) * (-t262 * t278 - t279 * t306) - g(3) * (-t269 * t278 + t288 * t279) 0, 0, 0, 0, 0, -g(1) * (t260 * t292 - t263 * t289) - g(2) * (t259 * t292 - t261 * t289) - g(3) * (t265 * t292 + t268 * t289) -g(1) * (-t260 * t289 - t263 * t292) - g(2) * (-t259 * t289 - t261 * t292) - g(3) * (-t265 * t289 + t268 * t292);];
U_reg  = t1;
