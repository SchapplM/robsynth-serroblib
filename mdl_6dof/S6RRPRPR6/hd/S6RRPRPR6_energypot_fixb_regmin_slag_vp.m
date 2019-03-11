% Calculate minimal parameter regressor of potential energy for
% S6RRPRPR6
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [11x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d1,d2,d4,d6,theta3]';
% 
% Output:
% U_reg [1x30]
%   minimal parameter regressor of Potential energy

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 10:44
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U_reg = S6RRPRPR6_energypot_fixb_regmin_slag_vp(qJ, g, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(11,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPRPR6_energypot_fixb_regmin_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RRPRPR6_energypot_fixb_regmin_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6RRPRPR6_energypot_fixb_regmin_slag_vp: pkin has to be [11x1] (double)');

%% Symbolic Calculation
% From energy_potential_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 10:42:54
% EndTime: 2019-03-09 10:42:54
% DurationCPUTime: 0.19s
% Computational Cost: add. (185->64), mult. (440->107), div. (0->0), fcn. (559->12), ass. (0->45)
t296 = pkin(8) + qJ(3);
t269 = sin(pkin(6));
t274 = sin(qJ(2));
t295 = t269 * t274;
t275 = sin(qJ(1));
t294 = t269 * t275;
t279 = cos(qJ(1));
t293 = t269 * t279;
t292 = t275 * t274;
t278 = cos(qJ(2));
t291 = t275 * t278;
t290 = t279 * t274;
t289 = t279 * t278;
t271 = cos(pkin(6));
t257 = t271 * t274 * pkin(2) - t296 * t269;
t265 = t278 * pkin(2) + pkin(1);
t288 = t279 * t257 + t275 * t265;
t287 = -t275 * t257 + t279 * t265;
t286 = pkin(2) * t295 + t296 * t271 + pkin(7);
t285 = g(1) * t275 - g(2) * t279;
t268 = sin(pkin(11));
t270 = cos(pkin(11));
t284 = t278 * t268 + t274 * t270;
t283 = t274 * t268 - t278 * t270;
t282 = t283 * t271;
t256 = t284 * t271;
t247 = t279 * t256 - t275 * t283;
t273 = sin(qJ(4));
t277 = cos(qJ(4));
t242 = t247 * t273 + t277 * t293;
t249 = -t275 * t256 - t279 * t283;
t244 = t249 * t273 - t277 * t294;
t255 = t284 * t269;
t250 = t255 * t273 - t271 * t277;
t281 = g(1) * t244 + g(2) * t242 + g(3) * t250;
t243 = t247 * t277 - t273 * t293;
t245 = t249 * t277 + t273 * t294;
t251 = t255 * t277 + t271 * t273;
t280 = g(1) * t245 + g(2) * t243 + g(3) * t251;
t276 = cos(qJ(6));
t272 = sin(qJ(6));
t254 = t283 * t269;
t248 = t275 * t282 - t279 * t284;
t246 = -t275 * t284 - t279 * t282;
t1 = [0, -g(1) * t279 - g(2) * t275, t285, 0, 0, 0, 0, 0, -g(1) * (-t271 * t292 + t289) - g(2) * (t271 * t290 + t291) - g(3) * t295, -g(1) * (-t271 * t291 - t290) - g(2) * (t271 * t289 - t292) - g(3) * t269 * t278, -g(3) * t271 - t285 * t269, -g(1) * t287 - g(2) * t288 - g(3) * t286, 0, 0, 0, 0, 0, -t280, t281, g(1) * t248 + g(2) * t246 - g(3) * t254, t280, -t281, -g(1) * (t249 * pkin(3) + t245 * pkin(4) - t248 * pkin(9) + t244 * qJ(5) + t287) - g(2) * (t247 * pkin(3) + t243 * pkin(4) - t246 * pkin(9) + t242 * qJ(5) + t288) - g(3) * (t255 * pkin(3) + t251 * pkin(4) + t254 * pkin(9) + t250 * qJ(5) + t286) 0, 0, 0, 0, 0, -g(1) * (t244 * t272 - t248 * t276) - g(2) * (t242 * t272 - t246 * t276) - g(3) * (t250 * t272 + t254 * t276) -g(1) * (t244 * t276 + t248 * t272) - g(2) * (t242 * t276 + t246 * t272) - g(3) * (t250 * t276 - t254 * t272);];
U_reg  = t1;
