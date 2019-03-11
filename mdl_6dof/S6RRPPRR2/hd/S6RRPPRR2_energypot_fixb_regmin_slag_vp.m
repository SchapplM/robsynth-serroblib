% Calculate minimal parameter regressor of potential energy for
% S6RRPPRR2
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [11x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d5,d6,theta3,theta4]';
% 
% Output:
% U_reg [1x30]
%   minimal parameter regressor of Potential energy

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 08:53
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U_reg = S6RRPPRR2_energypot_fixb_regmin_slag_vp(qJ, g, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(11,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPPRR2_energypot_fixb_regmin_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RRPPRR2_energypot_fixb_regmin_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6RRPPRR2_energypot_fixb_regmin_slag_vp: pkin has to be [11x1] (double)');

%% Symbolic Calculation
% From energy_potential_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 08:52:35
% EndTime: 2019-03-09 08:52:35
% DurationCPUTime: 0.13s
% Computational Cost: add. (114->48), mult. (111->71), div. (0->0), fcn. (117->12), ass. (0->37)
t231 = qJ(2) + pkin(10);
t225 = sin(t231);
t256 = g(3) * t225;
t235 = sin(qJ(2));
t255 = t235 * pkin(2) + pkin(6);
t230 = pkin(11) + qJ(5);
t228 = qJ(6) + t230;
t220 = sin(t228);
t236 = sin(qJ(1));
t254 = t236 * t220;
t221 = cos(t228);
t253 = t236 * t221;
t224 = sin(t230);
t252 = t236 * t224;
t226 = cos(t230);
t251 = t236 * t226;
t232 = sin(pkin(11));
t250 = t236 * t232;
t233 = cos(pkin(11));
t249 = t236 * t233;
t238 = cos(qJ(1));
t248 = t238 * t220;
t247 = t238 * t221;
t246 = t238 * t224;
t245 = t238 * t226;
t244 = t238 * t232;
t243 = t238 * t233;
t237 = cos(qJ(2));
t223 = t237 * pkin(2) + pkin(1);
t234 = -pkin(7) - qJ(3);
t242 = t236 * t223 + t238 * t234;
t241 = t238 * t223 - t236 * t234;
t240 = g(1) * t238 + g(2) * t236;
t227 = cos(t231);
t239 = pkin(3) * t227 + qJ(4) * t225;
t217 = g(1) * t236 - g(2) * t238;
t1 = [0, -t240, t217, 0, 0, 0, 0, 0, -g(3) * t235 - t240 * t237, -g(3) * t237 + t240 * t235, -t217, -g(1) * t241 - g(2) * t242 - g(3) * t255, -g(1) * (t227 * t243 + t250) - g(2) * (t227 * t249 - t244) - t233 * t256, -g(1) * (-t227 * t244 + t249) - g(2) * (-t227 * t250 - t243) + t232 * t256, g(3) * t227 - t240 * t225, -g(1) * (t239 * t238 + t241) - g(2) * (t239 * t236 + t242) - g(3) * (t225 * pkin(3) - t227 * qJ(4) + t255) 0, 0, 0, 0, 0, -g(1) * (t227 * t245 + t252) - g(2) * (t227 * t251 - t246) - t226 * t256, -g(1) * (-t227 * t246 + t251) - g(2) * (-t227 * t252 - t245) + t224 * t256, 0, 0, 0, 0, 0, -g(1) * (t227 * t247 + t254) - g(2) * (t227 * t253 - t248) - t221 * t256, -g(1) * (-t227 * t248 + t253) - g(2) * (-t227 * t254 - t247) + t220 * t256;];
U_reg  = t1;
