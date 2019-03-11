% Calculate minimal parameter regressor of potential energy for
% S6RRRPPR10
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [11x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d1,d2,d3,d6,theta5]';
% 
% Output:
% U_reg [1x32]
%   minimal parameter regressor of Potential energy

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 16:29
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U_reg = S6RRRPPR10_energypot_fixb_regmin_slag_vp(qJ, g, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(11,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRPPR10_energypot_fixb_regmin_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RRRPPR10_energypot_fixb_regmin_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6RRRPPR10_energypot_fixb_regmin_slag_vp: pkin has to be [11x1] (double)');

%% Symbolic Calculation
% From energy_potential_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 16:27:53
% EndTime: 2019-03-09 16:27:53
% DurationCPUTime: 0.17s
% Computational Cost: add. (184->71), mult. (408->106), div. (0->0), fcn. (508->12), ass. (0->40)
t290 = pkin(4) + pkin(9);
t266 = sin(pkin(6));
t270 = sin(qJ(2));
t289 = t266 * t270;
t271 = sin(qJ(1));
t288 = t266 * t271;
t272 = cos(qJ(3));
t287 = t266 * t272;
t273 = cos(qJ(2));
t286 = t266 * t273;
t274 = cos(qJ(1));
t285 = t266 * t274;
t284 = t271 * t270;
t283 = t271 * t273;
t282 = t274 * t270;
t281 = t274 * t273;
t268 = cos(pkin(6));
t269 = sin(qJ(3));
t247 = -t268 * t272 + t269 * t289;
t248 = t268 * t269 + t270 * t287;
t280 = pkin(2) * t289 + t248 * pkin(3) + t268 * pkin(8) + t247 * qJ(4) + pkin(7);
t252 = -t268 * t284 + t281;
t242 = t252 * t269 - t271 * t287;
t243 = t252 * t272 + t269 * t288;
t279 = t274 * pkin(1) + t252 * pkin(2) + t243 * pkin(3) + pkin(8) * t288 + t242 * qJ(4);
t250 = t268 * t282 + t283;
t240 = t250 * t269 + t272 * t285;
t278 = g(1) * t242 + g(2) * t240 + g(3) * t247;
t241 = t250 * t272 - t269 * t285;
t277 = g(1) * t243 + g(2) * t241 + g(3) * t248;
t249 = -t268 * t281 + t284;
t251 = t268 * t283 + t282;
t276 = -g(1) * t251 - g(2) * t249 + g(3) * t286;
t275 = t271 * pkin(1) + t250 * pkin(2) + t241 * pkin(3) - pkin(8) * t285 + t240 * qJ(4);
t267 = cos(pkin(11));
t265 = sin(pkin(11));
t264 = pkin(11) + qJ(6);
t260 = cos(t264);
t259 = sin(t264);
t1 = [0, -g(1) * t274 - g(2) * t271, g(1) * t271 - g(2) * t274, 0, 0, 0, 0, 0, -g(1) * t252 - g(2) * t250 - g(3) * t289, -t276, 0, 0, 0, 0, 0, -t277, t278, t276, t277, -t278, -g(1) * (t251 * pkin(9) + t279) - g(2) * (t249 * pkin(9) + t275) - g(3) * (-pkin(9) * t286 + t280) -g(1) * (t242 * t265 + t251 * t267) - g(2) * (t240 * t265 + t249 * t267) - g(3) * (t247 * t265 - t267 * t286) -g(1) * (t242 * t267 - t251 * t265) - g(2) * (t240 * t267 - t249 * t265) - g(3) * (t247 * t267 + t265 * t286) -t277, -g(1) * (t243 * qJ(5) + t251 * t290 + t279) - g(2) * (t241 * qJ(5) + t249 * t290 + t275) - g(3) * (t248 * qJ(5) - t286 * t290 + t280) 0, 0, 0, 0, 0, -g(1) * (t242 * t259 + t251 * t260) - g(2) * (t240 * t259 + t249 * t260) - g(3) * (t247 * t259 - t260 * t286) -g(1) * (t242 * t260 - t251 * t259) - g(2) * (t240 * t260 - t249 * t259) - g(3) * (t247 * t260 + t259 * t286);];
U_reg  = t1;
