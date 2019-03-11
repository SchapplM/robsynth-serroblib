% Calculate minimal parameter regressor of potential energy for
% S6RRPPRR11
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [11x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d1,d2,d5,d6,theta4]';
% 
% Output:
% U_reg [1x32]
%   minimal parameter regressor of Potential energy

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 09:44
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U_reg = S6RRPPRR11_energypot_fixb_regmin_slag_vp(qJ, g, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(11,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPPRR11_energypot_fixb_regmin_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RRPPRR11_energypot_fixb_regmin_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6RRPPRR11_energypot_fixb_regmin_slag_vp: pkin has to be [11x1] (double)');

%% Symbolic Calculation
% From energy_potential_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 09:42:43
% EndTime: 2019-03-09 09:42:43
% DurationCPUTime: 0.15s
% Computational Cost: add. (138->66), mult. (268->104), div. (0->0), fcn. (322->12), ass. (0->37)
t255 = sin(qJ(1));
t274 = g(1) * t255;
t258 = cos(qJ(1));
t273 = g(2) * t258;
t250 = sin(pkin(6));
t254 = sin(qJ(2));
t272 = t250 * t254;
t271 = t250 * t255;
t257 = cos(qJ(2));
t270 = t250 * t257;
t269 = t250 * t258;
t268 = t255 * t254;
t267 = t255 * t257;
t266 = t258 * t254;
t265 = t258 * t257;
t252 = cos(pkin(6));
t264 = pkin(2) * t272 + t252 * pkin(8) + pkin(7);
t263 = -t273 + t274;
t235 = -t252 * t265 + t268;
t236 = t252 * t266 + t267;
t262 = t255 * pkin(1) + t236 * pkin(2) + t235 * qJ(3);
t237 = t252 * t267 + t266;
t238 = -t252 * t268 + t265;
t261 = t258 * pkin(1) + t238 * pkin(2) + pkin(8) * t271 + t237 * qJ(3);
t260 = -g(1) * t237 - g(2) * t235 + g(3) * t270;
t259 = g(1) * t238 + g(2) * t236 + g(3) * t272;
t256 = cos(qJ(6));
t253 = sin(qJ(6));
t251 = cos(pkin(11));
t249 = sin(pkin(11));
t248 = pkin(11) + qJ(5);
t244 = cos(t248);
t243 = sin(t248);
t232 = -t243 * t270 + t252 * t244;
t231 = t235 * t243 - t244 * t269;
t230 = t237 * t243 + t244 * t271;
t1 = [0, -g(1) * t258 - g(2) * t255, t263, 0, 0, 0, 0, 0, -t259, -t260, -g(3) * t252 - t263 * t250, t259, t260, -g(1) * t261 - g(2) * (-pkin(8) * t269 + t262) - g(3) * (-qJ(3) * t270 + t264) -g(1) * (t237 * t249 + t251 * t271) - g(2) * (t235 * t249 - t251 * t269) - g(3) * (-t249 * t270 + t252 * t251) -g(1) * (t237 * t251 - t249 * t271) - g(2) * (t235 * t251 + t249 * t269) - g(3) * (-t252 * t249 - t251 * t270) -t259, -g(1) * (t238 * qJ(4) + t261) - g(2) * (t236 * qJ(4) + t262) - g(3) * (t252 * pkin(3) + t264) + (-pkin(3) * t274 - g(3) * (-qJ(3) * t257 + qJ(4) * t254) - (-pkin(3) - pkin(8)) * t273) * t250, 0, 0, 0, 0, 0, -g(1) * t230 - g(2) * t231 - g(3) * t232, -g(1) * (t237 * t244 - t243 * t271) - g(2) * (t235 * t244 + t243 * t269) - g(3) * (-t252 * t243 - t244 * t270) 0, 0, 0, 0, 0, -g(1) * (t230 * t256 + t238 * t253) - g(2) * (t231 * t256 + t236 * t253) - g(3) * (t232 * t256 + t253 * t272) -g(1) * (-t230 * t253 + t238 * t256) - g(2) * (-t231 * t253 + t236 * t256) - g(3) * (-t232 * t253 + t256 * t272);];
U_reg  = t1;
