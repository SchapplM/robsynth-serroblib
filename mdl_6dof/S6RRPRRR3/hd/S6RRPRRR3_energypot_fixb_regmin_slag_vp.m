% Calculate minimal parameter regressor of potential energy for
% S6RRPRRR3
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [11x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d4,d5,d6,theta3]';
% 
% Output:
% U_reg [1x33]
%   minimal parameter regressor of Potential energy

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 13:26
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U_reg = S6RRPRRR3_energypot_fixb_regmin_slag_vp(qJ, g, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(11,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPRRR3_energypot_fixb_regmin_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RRPRRR3_energypot_fixb_regmin_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6RRPRRR3_energypot_fixb_regmin_slag_vp: pkin has to be [11x1] (double)');

%% Symbolic Calculation
% From energy_potential_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 13:25:18
% EndTime: 2019-03-09 13:25:18
% DurationCPUTime: 0.08s
% Computational Cost: add. (88->40), mult. (86->60), div. (0->0), fcn. (95->12), ass. (0->32)
t233 = qJ(2) + pkin(11);
t255 = g(3) * sin(t233);
t234 = qJ(4) + qJ(5);
t232 = qJ(6) + t234;
t225 = sin(t232);
t238 = sin(qJ(1));
t254 = t238 * t225;
t226 = cos(t232);
t253 = t238 * t226;
t230 = sin(t234);
t252 = t238 * t230;
t231 = cos(t234);
t251 = t238 * t231;
t236 = sin(qJ(4));
t250 = t238 * t236;
t239 = cos(qJ(4));
t249 = t238 * t239;
t241 = cos(qJ(1));
t248 = t241 * t225;
t247 = t241 * t226;
t246 = t241 * t230;
t245 = t241 * t231;
t244 = t241 * t236;
t243 = t241 * t239;
t242 = g(1) * t241 + g(2) * t238;
t240 = cos(qJ(2));
t237 = sin(qJ(2));
t235 = -pkin(7) - qJ(3);
t229 = cos(t233);
t227 = t240 * pkin(2) + pkin(1);
t224 = g(1) * t238 - g(2) * t241;
t1 = [0, -t242, t224, 0, 0, 0, 0, 0, -g(3) * t237 - t242 * t240, -g(3) * t240 + t242 * t237, -t224, -g(1) * (t241 * t227 - t238 * t235) - g(2) * (t238 * t227 + t241 * t235) - g(3) * (t237 * pkin(2) + pkin(6)) 0, 0, 0, 0, 0, -g(1) * (t229 * t243 + t250) - g(2) * (t229 * t249 - t244) - t239 * t255, -g(1) * (-t229 * t244 + t249) - g(2) * (-t229 * t250 - t243) + t236 * t255, 0, 0, 0, 0, 0, -g(1) * (t229 * t245 + t252) - g(2) * (t229 * t251 - t246) - t231 * t255, -g(1) * (-t229 * t246 + t251) - g(2) * (-t229 * t252 - t245) + t230 * t255, 0, 0, 0, 0, 0, -g(1) * (t229 * t247 + t254) - g(2) * (t229 * t253 - t248) - t226 * t255, -g(1) * (-t229 * t248 + t253) - g(2) * (-t229 * t254 - t247) + t225 * t255;];
U_reg  = t1;
