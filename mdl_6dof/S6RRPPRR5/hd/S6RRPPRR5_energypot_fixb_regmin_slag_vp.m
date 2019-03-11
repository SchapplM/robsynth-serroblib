% Calculate minimal parameter regressor of potential energy for
% S6RRPPRR5
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d1,d2,d5,d6]';
% 
% Output:
% U_reg [1x32]
%   minimal parameter regressor of Potential energy

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 09:12
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U_reg = S6RRPPRR5_energypot_fixb_regmin_slag_vp(qJ, g, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(10,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPPRR5_energypot_fixb_regmin_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RRPPRR5_energypot_fixb_regmin_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RRPPRR5_energypot_fixb_regmin_slag_vp: pkin has to be [10x1] (double)');

%% Symbolic Calculation
% From energy_potential_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 09:10:58
% EndTime: 2019-03-09 09:10:58
% DurationCPUTime: 0.11s
% Computational Cost: add. (106->54), mult. (247->87), div. (0->0), fcn. (293->10), ass. (0->36)
t233 = sin(qJ(1));
t252 = g(1) * t233;
t237 = cos(qJ(1));
t251 = g(2) * t237;
t228 = sin(pkin(6));
t232 = sin(qJ(2));
t250 = t228 * t232;
t249 = t228 * t233;
t235 = cos(qJ(5));
t248 = t228 * t235;
t236 = cos(qJ(2));
t247 = t228 * t236;
t246 = t228 * t237;
t245 = t233 * t232;
t244 = t233 * t236;
t243 = t237 * t232;
t242 = t237 * t236;
t229 = cos(pkin(6));
t241 = pkin(2) * t250 + t229 * pkin(8) + pkin(7);
t240 = -t251 + t252;
t218 = -t229 * t242 + t245;
t219 = t229 * t243 + t244;
t239 = t233 * pkin(1) + t219 * pkin(2) + t218 * qJ(3);
t220 = t229 * t244 + t243;
t221 = -t229 * t245 + t242;
t238 = t237 * pkin(1) + t221 * pkin(2) + pkin(8) * t249 + t220 * qJ(3);
t210 = -g(1) * t220 - g(2) * t218 + g(3) * t247;
t234 = cos(qJ(6));
t231 = sin(qJ(5));
t230 = sin(qJ(6));
t217 = -t229 * t231 + t232 * t248;
t214 = g(3) * t229 + t240 * t228;
t213 = t221 * t235 - t231 * t249;
t212 = t219 * t235 + t231 * t246;
t211 = -g(1) * t221 - g(2) * t219 - g(3) * t250;
t1 = [0, -g(1) * t237 - g(2) * t233, t240, 0, 0, 0, 0, 0, t211, -t210, t211, -t214, t210, -g(1) * t238 - g(2) * (-pkin(8) * t246 + t239) - g(3) * (-qJ(3) * t247 + t241) t211, t210, t214, -g(1) * (t221 * pkin(3) + t238) - g(2) * (t219 * pkin(3) + t239) - g(3) * (-t229 * qJ(4) + t241) + (qJ(4) * t252 - g(3) * (pkin(3) * t232 - qJ(3) * t236) - (-pkin(8) + qJ(4)) * t251) * t228, 0, 0, 0, 0, 0, -g(1) * t213 - g(2) * t212 - g(3) * t217, -g(1) * (-t221 * t231 - t233 * t248) - g(2) * (-t219 * t231 + t235 * t246) - g(3) * (-t229 * t235 - t231 * t250) 0, 0, 0, 0, 0, -g(1) * (t213 * t234 - t220 * t230) - g(2) * (t212 * t234 - t218 * t230) - g(3) * (t217 * t234 + t230 * t247) -g(1) * (-t213 * t230 - t220 * t234) - g(2) * (-t212 * t230 - t218 * t234) - g(3) * (-t217 * t230 + t234 * t247);];
U_reg  = t1;
