% Calculate minimal parameter regressor of potential energy for
% S6PRRRPR1
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [12x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d2,d3,d4,d6,theta1,theta5]';
% 
% Output:
% U_reg [1x27]
%   minimal parameter regressor of Potential energy

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-08 23:04
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U_reg = S6PRRRPR1_energypot_fixb_regmin_slag_vp(qJ, g, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(12,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRRRPR1_energypot_fixb_regmin_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6PRRRPR1_energypot_fixb_regmin_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [12 1]), ...
  'S6PRRRPR1_energypot_fixb_regmin_slag_vp: pkin has to be [12x1] (double)');

%% Symbolic Calculation
% From energy_potential_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-08 23:03:21
% EndTime: 2019-03-08 23:03:22
% DurationCPUTime: 0.18s
% Computational Cost: add. (140->61), mult. (212->104), div. (0->0), fcn. (258->14), ass. (0->36)
t231 = qJ(3) + qJ(4);
t228 = sin(t231);
t237 = sin(qJ(3));
t252 = t237 * pkin(3) + pkin(4) * t228 + pkin(7);
t232 = sin(pkin(11));
t233 = sin(pkin(6));
t250 = t232 * t233;
t234 = cos(pkin(11));
t249 = t233 * t234;
t248 = t233 * t237;
t238 = sin(qJ(2));
t247 = t233 * t238;
t240 = cos(qJ(3));
t246 = t233 * t240;
t241 = cos(qJ(2));
t245 = t233 * t241;
t235 = cos(pkin(6));
t244 = t235 * t238;
t243 = t235 * t241;
t218 = t232 * t238 - t234 * t243;
t220 = t232 * t243 + t234 * t238;
t242 = -g(1) * t220 - g(2) * t218 + g(3) * t245;
t239 = cos(qJ(6));
t236 = sin(qJ(6));
t230 = -qJ(5) - pkin(9) - pkin(8);
t229 = cos(t231);
t227 = pkin(12) + t231;
t226 = cos(t227);
t225 = sin(t227);
t222 = t240 * pkin(3) + pkin(4) * t229 + pkin(2);
t221 = -t232 * t244 + t234 * t241;
t219 = t232 * t241 + t234 * t244;
t217 = t235 * t225 + t226 * t247;
t216 = t221 * t226 + t225 * t250;
t215 = t219 * t226 - t225 * t249;
t1 = [-g(3) * qJ(1), 0, -g(1) * t221 - g(2) * t219 - g(3) * t247, -t242, 0, 0, 0, 0, 0, -g(1) * (t221 * t240 + t232 * t248) - g(2) * (t219 * t240 - t234 * t248) - g(3) * (t235 * t237 + t238 * t246) -g(1) * (-t221 * t237 + t232 * t246) - g(2) * (-t219 * t237 - t234 * t246) - g(3) * (t235 * t240 - t237 * t247) 0, 0, 0, 0, 0, -g(1) * (t221 * t229 + t228 * t250) - g(2) * (t219 * t229 - t228 * t249) - g(3) * (t235 * t228 + t229 * t247) -g(1) * (-t221 * t228 + t229 * t250) - g(2) * (-t219 * t228 - t229 * t249) - g(3) * (-t228 * t247 + t235 * t229) t242, -g(1) * (t234 * pkin(1) - t220 * t230 + t221 * t222) - g(2) * (t232 * pkin(1) - t218 * t230 + t219 * t222) - g(3) * (t252 * t235 + qJ(1)) + (-g(3) * (t222 * t238 + t230 * t241) + (-g(1) * t232 + g(2) * t234) * t252) * t233, 0, 0, 0, 0, 0, -g(1) * (t216 * t239 + t220 * t236) - g(2) * (t215 * t239 + t218 * t236) - g(3) * (t217 * t239 - t236 * t245) -g(1) * (-t216 * t236 + t220 * t239) - g(2) * (-t215 * t236 + t218 * t239) - g(3) * (-t217 * t236 - t239 * t245);];
U_reg  = t1;
