% Calculate minimal parameter regressor of potential energy for
% S6RRRPPR8
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d1,d2,d3,d6]';
% 
% Output:
% U_reg [1x32]
%   minimal parameter regressor of Potential energy

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 16:10
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U_reg = S6RRRPPR8_energypot_fixb_regmin_slag_vp(qJ, g, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(10,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRPPR8_energypot_fixb_regmin_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RRRPPR8_energypot_fixb_regmin_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RRRPPR8_energypot_fixb_regmin_slag_vp: pkin has to be [10x1] (double)');

%% Symbolic Calculation
% From energy_potential_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 16:09:29
% EndTime: 2019-03-09 16:09:29
% DurationCPUTime: 0.15s
% Computational Cost: add. (159->62), mult. (374->90), div. (0->0), fcn. (461->10), ass. (0->35)
t263 = pkin(9) - qJ(5);
t262 = cos(pkin(6));
t240 = sin(pkin(6));
t243 = sin(qJ(2));
t261 = t240 * t243;
t244 = sin(qJ(1));
t260 = t240 * t244;
t246 = cos(qJ(3));
t259 = t240 * t246;
t247 = cos(qJ(2));
t258 = t240 * t247;
t248 = cos(qJ(1));
t257 = t240 * t248;
t256 = t244 * t262;
t255 = t248 * t262;
t242 = sin(qJ(3));
t226 = t242 * t261 - t262 * t246;
t227 = t262 * t242 + t243 * t259;
t254 = pkin(2) * t261 + t227 * pkin(3) + t262 * pkin(8) + t226 * qJ(4) + pkin(7);
t231 = -t243 * t256 + t248 * t247;
t221 = t231 * t242 - t244 * t259;
t222 = t231 * t246 + t242 * t260;
t253 = t248 * pkin(1) + t231 * pkin(2) + t222 * pkin(3) + pkin(8) * t260 + t221 * qJ(4);
t229 = t243 * t255 + t244 * t247;
t219 = t229 * t242 + t246 * t257;
t252 = g(1) * t221 + g(2) * t219 + g(3) * t226;
t220 = t229 * t246 - t242 * t257;
t251 = g(1) * t222 + g(2) * t220 + g(3) * t227;
t228 = t244 * t243 - t247 * t255;
t230 = t248 * t243 + t247 * t256;
t250 = -g(1) * t230 - g(2) * t228 + g(3) * t258;
t249 = t244 * pkin(1) + t229 * pkin(2) + t220 * pkin(3) - pkin(8) * t257 + t219 * qJ(4);
t245 = cos(qJ(6));
t241 = sin(qJ(6));
t1 = [0, -g(1) * t248 - g(2) * t244, g(1) * t244 - g(2) * t248, 0, 0, 0, 0, 0, -g(1) * t231 - g(2) * t229 - g(3) * t261, -t250, 0, 0, 0, 0, 0, -t251, t252, -t251, t250, -t252, -g(1) * (t230 * pkin(9) + t253) - g(2) * (t228 * pkin(9) + t249) - g(3) * (-pkin(9) * t258 + t254) -t252, t251, -t250, -g(1) * (t222 * pkin(4) + t263 * t230 + t253) - g(2) * (t220 * pkin(4) + t263 * t228 + t249) - g(3) * (t227 * pkin(4) - t263 * t258 + t254) 0, 0, 0, 0, 0, -g(1) * (t221 * t245 - t230 * t241) - g(2) * (t219 * t245 - t228 * t241) - g(3) * (t226 * t245 + t241 * t258) -g(1) * (-t221 * t241 - t230 * t245) - g(2) * (-t219 * t241 - t228 * t245) - g(3) * (-t226 * t241 + t245 * t258);];
U_reg  = t1;
