% Calculate minimal parameter regressor of potential energy for
% S6RRPRRR12
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [11x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d1,d2,d4,d5,d6]';
% 
% Output:
% U_reg [1x35]
%   minimal parameter regressor of Potential energy

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 14:43
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U_reg = S6RRPRRR12_energypot_fixb_regmin_slag_vp(qJ, g, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(11,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPRRR12_energypot_fixb_regmin_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RRPRRR12_energypot_fixb_regmin_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6RRPRRR12_energypot_fixb_regmin_slag_vp: pkin has to be [11x1] (double)');

%% Symbolic Calculation
% From energy_potential_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 14:41:42
% EndTime: 2019-03-09 14:41:42
% DurationCPUTime: 0.10s
% Computational Cost: add. (111->55), mult. (213->93), div. (0->0), fcn. (262->12), ass. (0->32)
t246 = sin(pkin(6));
t250 = sin(qJ(2));
t266 = t246 * t250;
t251 = sin(qJ(1));
t265 = t246 * t251;
t254 = cos(qJ(2));
t264 = t246 * t254;
t255 = cos(qJ(1));
t263 = t246 * t255;
t262 = t251 * t250;
t261 = t251 * t254;
t260 = t255 * t250;
t259 = t255 * t254;
t258 = g(1) * t251 - g(2) * t255;
t247 = cos(pkin(6));
t238 = -t247 * t259 + t262;
t240 = t247 * t261 + t260;
t257 = -g(1) * t240 - g(2) * t238 + g(3) * t264;
t239 = t247 * t260 + t261;
t241 = -t247 * t262 + t259;
t256 = g(1) * t241 + g(2) * t239 + g(3) * t266;
t253 = cos(qJ(4));
t252 = cos(qJ(6));
t249 = sin(qJ(4));
t248 = sin(qJ(6));
t245 = qJ(4) + qJ(5);
t244 = cos(t245);
t243 = sin(t245);
t237 = -t243 * t264 + t247 * t244;
t236 = t238 * t243 - t244 * t263;
t235 = t240 * t243 + t244 * t265;
t1 = [0, -g(1) * t255 - g(2) * t251, t258, 0, 0, 0, 0, 0, -t256, -t257, -g(3) * t247 - t258 * t246, t256, t257, -g(1) * (t255 * pkin(1) + t241 * pkin(2) + pkin(8) * t265 + t240 * qJ(3)) - g(2) * (t251 * pkin(1) + t239 * pkin(2) - pkin(8) * t263 + t238 * qJ(3)) - g(3) * (t247 * pkin(8) + pkin(7) + (pkin(2) * t250 - qJ(3) * t254) * t246) 0, 0, 0, 0, 0, -g(1) * (t240 * t249 + t253 * t265) - g(2) * (t238 * t249 - t253 * t263) - g(3) * (t247 * t253 - t249 * t264) -g(1) * (t240 * t253 - t249 * t265) - g(2) * (t238 * t253 + t249 * t263) - g(3) * (-t247 * t249 - t253 * t264) 0, 0, 0, 0, 0, -g(1) * t235 - g(2) * t236 - g(3) * t237, -g(1) * (t240 * t244 - t243 * t265) - g(2) * (t238 * t244 + t243 * t263) - g(3) * (-t247 * t243 - t244 * t264) 0, 0, 0, 0, 0, -g(1) * (t235 * t252 + t241 * t248) - g(2) * (t236 * t252 + t239 * t248) - g(3) * (t237 * t252 + t248 * t266) -g(1) * (-t235 * t248 + t241 * t252) - g(2) * (-t236 * t248 + t239 * t252) - g(3) * (-t237 * t248 + t252 * t266);];
U_reg  = t1;
