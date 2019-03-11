% Calculate minimal parameter regressor of potential energy for
% S6RRRRPP7
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [11x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d1,d2,d3,d4,theta5]';
% 
% Output:
% U_reg [1x30]
%   minimal parameter regressor of Potential energy

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 21:29
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U_reg = S6RRRRPP7_energypot_fixb_regmin_slag_vp(qJ, g, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(11,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRRPP7_energypot_fixb_regmin_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RRRRPP7_energypot_fixb_regmin_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6RRRRPP7_energypot_fixb_regmin_slag_vp: pkin has to be [11x1] (double)');

%% Symbolic Calculation
% From energy_potential_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 21:26:30
% EndTime: 2019-03-09 21:26:31
% DurationCPUTime: 0.15s
% Computational Cost: add. (220->73), mult. (446->110), div. (0->0), fcn. (556->12), ass. (0->47)
t288 = cos(pkin(6));
t296 = cos(qJ(2));
t297 = cos(qJ(1));
t302 = t297 * t296;
t292 = sin(qJ(2));
t293 = sin(qJ(1));
t305 = t293 * t292;
t271 = -t288 * t302 + t305;
t290 = sin(qJ(4));
t312 = t271 * t290;
t303 = t297 * t292;
t304 = t293 * t296;
t273 = t288 * t304 + t303;
t311 = t273 * t290;
t287 = sin(pkin(6));
t310 = t287 * t292;
t309 = t287 * t293;
t295 = cos(qJ(3));
t308 = t287 * t295;
t307 = t287 * t296;
t306 = t287 * t297;
t272 = t288 * t303 + t304;
t291 = sin(qJ(3));
t259 = t272 * t291 + t295 * t306;
t274 = -t288 * t305 + t302;
t261 = t274 * t291 - t293 * t308;
t269 = -t288 * t295 + t291 * t310;
t301 = g(1) * t261 + g(2) * t259 + g(3) * t269;
t262 = t274 * t295 + t291 * t309;
t294 = cos(qJ(4));
t280 = t294 * pkin(4) + pkin(3);
t289 = -qJ(5) - pkin(10);
t300 = t297 * pkin(1) + t274 * pkin(2) + pkin(4) * t311 + pkin(8) * t309 + t273 * pkin(9) - t261 * t289 + t262 * t280;
t260 = t272 * t295 - t291 * t306;
t299 = t293 * pkin(1) + t272 * pkin(2) + pkin(4) * t312 - pkin(8) * t306 + t271 * pkin(9) - t259 * t289 + t260 * t280;
t270 = t288 * t291 + t292 * t308;
t298 = pkin(7) + t270 * t280 - t269 * t289 + pkin(2) * t310 + t288 * pkin(8) + (-pkin(4) * t290 - pkin(9)) * t307;
t286 = qJ(4) + pkin(11);
t282 = cos(t286);
t281 = sin(t286);
t256 = t270 * t282 - t281 * t307;
t255 = t270 * t281 + t282 * t307;
t252 = t262 * t282 + t273 * t281;
t251 = t262 * t281 - t273 * t282;
t250 = t260 * t282 + t271 * t281;
t249 = t260 * t281 - t271 * t282;
t1 = [0, -g(1) * t297 - g(2) * t293, g(1) * t293 - g(2) * t297, 0, 0, 0, 0, 0, -g(1) * t274 - g(2) * t272 - g(3) * t310, g(1) * t273 + g(2) * t271 - g(3) * t307, 0, 0, 0, 0, 0, -g(1) * t262 - g(2) * t260 - g(3) * t270, t301, 0, 0, 0, 0, 0, -g(1) * (t262 * t294 + t311) - g(2) * (t260 * t294 + t312) - g(3) * (t270 * t294 - t290 * t307) -g(1) * (-t262 * t290 + t273 * t294) - g(2) * (-t260 * t290 + t271 * t294) - g(3) * (-t270 * t290 - t294 * t307) -t301, -g(1) * t300 - g(2) * t299 - g(3) * t298, -g(1) * t252 - g(2) * t250 - g(3) * t256, -t301, -g(1) * t251 - g(2) * t249 - g(3) * t255, -g(1) * (t252 * pkin(5) + t251 * qJ(6) + t300) - g(2) * (t250 * pkin(5) + t249 * qJ(6) + t299) - g(3) * (t256 * pkin(5) + t255 * qJ(6) + t298);];
U_reg  = t1;
