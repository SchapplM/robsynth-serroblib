% Calculate minimal parameter regressor of potential energy for
% S6RRRPPR1
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [11x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d3,d6,theta4,theta5]';
% 
% Output:
% U_reg [1x30]
%   minimal parameter regressor of Potential energy

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 15:23
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U_reg = S6RRRPPR1_energypot_fixb_regmin_slag_vp(qJ, g, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(11,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRPPR1_energypot_fixb_regmin_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RRRPPR1_energypot_fixb_regmin_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6RRRPPR1_energypot_fixb_regmin_slag_vp: pkin has to be [11x1] (double)');

%% Symbolic Calculation
% From energy_potential_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 15:23:05
% EndTime: 2019-03-09 15:23:05
% DurationCPUTime: 0.08s
% Computational Cost: add. (127->44), mult. (107->63), div. (0->0), fcn. (109->12), ass. (0->33)
t232 = qJ(2) + qJ(3);
t226 = pkin(10) + t232;
t221 = sin(t226);
t252 = g(3) * t221;
t231 = pkin(11) + qJ(6);
t224 = sin(t231);
t236 = sin(qJ(1));
t251 = t236 * t224;
t225 = cos(t231);
t250 = t236 * t225;
t233 = sin(pkin(11));
t249 = t236 * t233;
t234 = cos(pkin(11));
t248 = t236 * t234;
t238 = cos(qJ(1));
t247 = t238 * t224;
t246 = t238 * t225;
t245 = t238 * t233;
t244 = t238 * t234;
t228 = cos(t232);
t237 = cos(qJ(2));
t218 = t237 * pkin(2) + pkin(3) * t228 + pkin(1);
t230 = -qJ(4) - pkin(8) - pkin(7);
t243 = t236 * t218 + t238 * t230;
t227 = sin(t232);
t235 = sin(qJ(2));
t242 = t235 * pkin(2) + pkin(3) * t227 + pkin(6);
t241 = t238 * t218 - t236 * t230;
t240 = g(1) * t238 + g(2) * t236;
t222 = cos(t226);
t239 = pkin(4) * t222 + qJ(5) * t221;
t219 = g(1) * t236 - g(2) * t238;
t1 = [0, -t240, t219, 0, 0, 0, 0, 0, -g(3) * t235 - t240 * t237, -g(3) * t237 + t240 * t235, 0, 0, 0, 0, 0, -g(3) * t227 - t240 * t228, -g(3) * t228 + t240 * t227, -t219, -g(1) * t241 - g(2) * t243 - g(3) * t242, -g(1) * (t222 * t244 + t249) - g(2) * (t222 * t248 - t245) - t234 * t252, -g(1) * (-t222 * t245 + t248) - g(2) * (-t222 * t249 - t244) + t233 * t252, g(3) * t222 - t240 * t221, -g(1) * (t239 * t238 + t241) - g(2) * (t239 * t236 + t243) - g(3) * (t221 * pkin(4) - t222 * qJ(5) + t242) 0, 0, 0, 0, 0, -g(1) * (t222 * t246 + t251) - g(2) * (t222 * t250 - t247) - t225 * t252, -g(1) * (-t222 * t247 + t250) - g(2) * (-t222 * t251 - t246) + t224 * t252;];
U_reg  = t1;
