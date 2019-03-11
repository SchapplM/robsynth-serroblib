% Calculate minimal parameter regressor of potential energy for
% S6RRRRPR11
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [12x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d1,d2,d3,d4,d6,theta5]';
% 
% Output:
% U_reg [1x33]
%   minimal parameter regressor of Potential energy

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 23:28
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U_reg = S6RRRRPR11_energypot_fixb_regmin_slag_vp(qJ, g, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(12,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRRPR11_energypot_fixb_regmin_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RRRRPR11_energypot_fixb_regmin_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [12 1]), ...
  'S6RRRRPR11_energypot_fixb_regmin_slag_vp: pkin has to be [12x1] (double)');

%% Symbolic Calculation
% From energy_potential_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 23:23:33
% EndTime: 2019-03-09 23:23:33
% DurationCPUTime: 0.12s
% Computational Cost: add. (145->63), mult. (284->102), div. (0->0), fcn. (356->12), ass. (0->38)
t268 = cos(pkin(6));
t276 = cos(qJ(2));
t277 = cos(qJ(1));
t280 = t277 * t276;
t272 = sin(qJ(2));
t273 = sin(qJ(1));
t283 = t273 * t272;
t256 = -t268 * t280 + t283;
t270 = sin(qJ(4));
t289 = t256 * t270;
t281 = t277 * t272;
t282 = t273 * t276;
t258 = t268 * t282 + t281;
t288 = t258 * t270;
t267 = sin(pkin(6));
t287 = t267 * t272;
t275 = cos(qJ(3));
t286 = t267 * t275;
t285 = t267 * t276;
t284 = t267 * t277;
t279 = g(1) * t273 - g(2) * t277;
t257 = t268 * t281 + t282;
t271 = sin(qJ(3));
t250 = t257 * t271 + t275 * t284;
t259 = -t268 * t283 + t280;
t252 = t259 * t271 - t273 * t286;
t254 = -t268 * t275 + t271 * t287;
t278 = g(1) * t252 + g(2) * t250 + g(3) * t254;
t274 = cos(qJ(4));
t269 = -qJ(5) - pkin(10);
t266 = qJ(4) + pkin(12) + qJ(6);
t265 = t274 * pkin(4) + pkin(3);
t264 = cos(t266);
t263 = sin(t266);
t255 = t268 * t271 + t272 * t286;
t253 = t273 * t267 * t271 + t259 * t275;
t251 = t257 * t275 - t271 * t284;
t1 = [0, -g(1) * t277 - g(2) * t273, t279, 0, 0, 0, 0, 0, -g(1) * t259 - g(2) * t257 - g(3) * t287, g(1) * t258 + g(2) * t256 - g(3) * t285, 0, 0, 0, 0, 0, -g(1) * t253 - g(2) * t251 - g(3) * t255, t278, 0, 0, 0, 0, 0, -g(1) * (t253 * t274 + t288) - g(2) * (t251 * t274 + t289) - g(3) * (t255 * t274 - t270 * t285) -g(1) * (-t253 * t270 + t258 * t274) - g(2) * (-t251 * t270 + t256 * t274) - g(3) * (-t255 * t270 - t274 * t285) -t278, -g(1) * (t277 * pkin(1) + t259 * pkin(2) + pkin(4) * t288 + t258 * pkin(9) - t252 * t269 + t253 * t265) - g(2) * (t273 * pkin(1) + t257 * pkin(2) + pkin(4) * t289 + t256 * pkin(9) - t250 * t269 + t251 * t265) - g(3) * (t268 * pkin(8) - t254 * t269 + t255 * t265 + pkin(7)) + (-g(3) * (pkin(2) * t272 + (-pkin(4) * t270 - pkin(9)) * t276) - t279 * pkin(8)) * t267, 0, 0, 0, 0, 0, -g(1) * (t253 * t264 + t258 * t263) - g(2) * (t251 * t264 + t256 * t263) - g(3) * (t255 * t264 - t263 * t285) -g(1) * (-t253 * t263 + t258 * t264) - g(2) * (-t251 * t263 + t256 * t264) - g(3) * (-t255 * t263 - t264 * t285);];
U_reg  = t1;
