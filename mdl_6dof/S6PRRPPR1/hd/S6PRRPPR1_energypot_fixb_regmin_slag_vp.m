% Calculate minimal parameter regressor of potential energy for
% S6PRRPPR1
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [12x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d2,d3,d6,theta1,theta4,theta5]';
% 
% Output:
% U_reg [1x24]
%   minimal parameter regressor of Potential energy

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-08 21:02
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U_reg = S6PRRPPR1_energypot_fixb_regmin_slag_vp(qJ, g, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(12,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRRPPR1_energypot_fixb_regmin_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6PRRPPR1_energypot_fixb_regmin_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [12 1]), ...
  'S6PRRPPR1_energypot_fixb_regmin_slag_vp: pkin has to be [12x1] (double)');

%% Symbolic Calculation
% From energy_potential_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-08 21:02:11
% EndTime: 2019-03-08 21:02:12
% DurationCPUTime: 0.16s
% Computational Cost: add. (200->74), mult. (327->119), div. (0->0), fcn. (399->14), ass. (0->43)
t251 = sin(pkin(10));
t252 = sin(pkin(6));
t274 = t251 * t252;
t254 = cos(pkin(10));
t273 = t252 * t254;
t257 = sin(qJ(3));
t272 = t252 * t257;
t258 = sin(qJ(2));
t271 = t252 * t258;
t259 = cos(qJ(3));
t270 = t252 * t259;
t260 = cos(qJ(2));
t269 = t252 * t260;
t255 = cos(pkin(6));
t268 = t255 * t257;
t267 = t255 * t258;
t266 = t255 * t260;
t265 = t251 * t272;
t232 = t251 * t266 + t254 * t258;
t233 = -t251 * t267 + t254 * t260;
t240 = t259 * pkin(3) + pkin(2);
t256 = -qJ(4) - pkin(8);
t264 = t254 * pkin(1) + pkin(3) * t265 + pkin(7) * t274 - t232 * t256 + t233 * t240;
t263 = pkin(3) * t268 + t255 * pkin(7) + t240 * t271 + t256 * t269 + qJ(1);
t230 = t251 * t258 - t254 * t266;
t262 = -g(1) * t232 - g(2) * t230 + g(3) * t269;
t231 = t251 * t260 + t254 * t267;
t261 = t231 * t240 - t230 * t256 + t251 * pkin(1) + (-pkin(3) * t257 - pkin(7)) * t273;
t253 = cos(pkin(12));
t250 = sin(pkin(12));
t249 = qJ(3) + pkin(11);
t248 = pkin(12) + qJ(6);
t244 = cos(t249);
t243 = cos(t248);
t242 = sin(t249);
t241 = sin(t248);
t227 = t255 * t242 + t244 * t271;
t226 = t242 * t271 - t255 * t244;
t223 = t233 * t244 + t242 * t274;
t222 = t233 * t242 - t244 * t274;
t221 = t231 * t244 - t242 * t273;
t220 = t231 * t242 + t244 * t273;
t1 = [-g(3) * qJ(1), 0, -g(1) * t233 - g(2) * t231 - g(3) * t271, -t262, 0, 0, 0, 0, 0, -g(1) * (t233 * t259 + t265) - g(2) * (t231 * t259 - t254 * t272) - g(3) * (t258 * t270 + t268) -g(1) * (-t233 * t257 + t251 * t270) - g(2) * (-t231 * t257 - t254 * t270) - g(3) * (t255 * t259 - t257 * t271) t262, -g(1) * t264 - g(2) * t261 - g(3) * t263, -g(1) * (t223 * t253 + t232 * t250) - g(2) * (t221 * t253 + t230 * t250) - g(3) * (t227 * t253 - t250 * t269) -g(1) * (-t223 * t250 + t232 * t253) - g(2) * (-t221 * t250 + t230 * t253) - g(3) * (-t227 * t250 - t253 * t269) -g(1) * t222 - g(2) * t220 - g(3) * t226, -g(1) * (t223 * pkin(4) + t222 * qJ(5) + t264) - g(2) * (t221 * pkin(4) + t220 * qJ(5) + t261) - g(3) * (t227 * pkin(4) + t226 * qJ(5) + t263) 0, 0, 0, 0, 0, -g(1) * (t223 * t243 + t232 * t241) - g(2) * (t221 * t243 + t230 * t241) - g(3) * (t227 * t243 - t241 * t269) -g(1) * (-t223 * t241 + t232 * t243) - g(2) * (-t221 * t241 + t230 * t243) - g(3) * (-t227 * t241 - t243 * t269);];
U_reg  = t1;
