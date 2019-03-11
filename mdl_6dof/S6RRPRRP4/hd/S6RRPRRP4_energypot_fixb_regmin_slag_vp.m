% Calculate minimal parameter regressor of potential energy for
% S6RRPRRP4
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d4,d5,theta3]';
% 
% Output:
% U_reg [1x30]
%   minimal parameter regressor of Potential energy

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 11:56
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U_reg = S6RRPRRP4_energypot_fixb_regmin_slag_vp(qJ, g, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(10,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPRRP4_energypot_fixb_regmin_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RRPRRP4_energypot_fixb_regmin_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RRPRRP4_energypot_fixb_regmin_slag_vp: pkin has to be [10x1] (double)');

%% Symbolic Calculation
% From energy_potential_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 11:55:39
% EndTime: 2019-03-09 11:55:39
% DurationCPUTime: 0.12s
% Computational Cost: add. (138->50), mult. (138->67), div. (0->0), fcn. (148->10), ass. (0->38)
t225 = qJ(2) + pkin(10);
t220 = sin(t225);
t249 = g(3) * t220;
t229 = sin(qJ(2));
t248 = t229 * pkin(2) + pkin(6);
t231 = cos(qJ(4));
t218 = t231 * pkin(4) + pkin(3);
t221 = cos(t225);
t247 = t218 * t221;
t226 = qJ(4) + qJ(5);
t222 = sin(t226);
t230 = sin(qJ(1));
t246 = t230 * t222;
t223 = cos(t226);
t245 = t230 * t223;
t228 = sin(qJ(4));
t244 = t230 * t228;
t243 = t230 * t231;
t233 = cos(qJ(1));
t242 = t233 * t222;
t241 = t233 * t223;
t240 = t233 * t228;
t239 = t233 * t231;
t232 = cos(qJ(2));
t219 = t232 * pkin(2) + pkin(1);
t227 = -pkin(7) - qJ(3);
t238 = t230 * t219 + t233 * t227;
t237 = t233 * t219 - t230 * t227;
t236 = g(1) * t233 + g(2) * t230;
t209 = t221 * t246 + t241;
t211 = t221 * t242 - t245;
t235 = g(1) * t211 + g(2) * t209 + t222 * t249;
t234 = -pkin(9) - pkin(8);
t213 = g(1) * t230 - g(2) * t233;
t212 = t221 * t241 + t246;
t210 = t221 * t245 - t242;
t208 = -g(1) * t212 - g(2) * t210 - t223 * t249;
t1 = [0, -t236, t213, 0, 0, 0, 0, 0, -g(3) * t229 - t236 * t232, -g(3) * t232 + t236 * t229, -t213, -g(1) * t237 - g(2) * t238 - g(3) * t248, 0, 0, 0, 0, 0, -g(1) * (t221 * t239 + t244) - g(2) * (t221 * t243 - t240) - t231 * t249, -g(1) * (-t221 * t240 + t243) - g(2) * (-t221 * t244 - t239) + t228 * t249, 0, 0, 0, 0, 0, t208, t235, t208, g(3) * t221 - t236 * t220, -t235, -g(1) * (pkin(4) * t244 + t212 * pkin(5) + t211 * qJ(6) + t233 * t247 + t237) - g(2) * (-pkin(4) * t240 + t210 * pkin(5) + t209 * qJ(6) + t230 * t247 + t238) - g(3) * (t221 * t234 + t248) + (-g(3) * (pkin(5) * t223 + qJ(6) * t222 + t218) + t236 * t234) * t220;];
U_reg  = t1;
