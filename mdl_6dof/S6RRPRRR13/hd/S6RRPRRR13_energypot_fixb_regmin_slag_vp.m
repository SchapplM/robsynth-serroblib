% Calculate minimal parameter regressor of potential energy for
% S6RRPRRR13
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
% Datum: 2019-03-09 14:56
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U_reg = S6RRPRRR13_energypot_fixb_regmin_slag_vp(qJ, g, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(11,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPRRR13_energypot_fixb_regmin_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RRPRRR13_energypot_fixb_regmin_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6RRPRRR13_energypot_fixb_regmin_slag_vp: pkin has to be [11x1] (double)');

%% Symbolic Calculation
% From energy_potential_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 14:53:18
% EndTime: 2019-03-09 14:53:19
% DurationCPUTime: 0.10s
% Computational Cost: add. (109->55), mult. (239->93), div. (0->0), fcn. (298->12), ass. (0->32)
t244 = sin(pkin(6));
t248 = sin(qJ(2));
t264 = t244 * t248;
t249 = sin(qJ(1));
t263 = t244 * t249;
t252 = cos(qJ(2));
t262 = t244 * t252;
t253 = cos(qJ(1));
t261 = t244 * t253;
t260 = t249 * t248;
t259 = t249 * t252;
t258 = t253 * t248;
t257 = t253 * t252;
t256 = g(1) * t249 - g(2) * t253;
t245 = cos(pkin(6));
t236 = -t245 * t257 + t260;
t238 = t245 * t259 + t258;
t255 = -g(1) * t238 - g(2) * t236 + g(3) * t262;
t237 = t245 * t258 + t259;
t239 = -t245 * t260 + t257;
t254 = g(1) * t239 + g(2) * t237 + g(3) * t264;
t251 = cos(qJ(4));
t250 = cos(qJ(5));
t247 = sin(qJ(4));
t246 = sin(qJ(5));
t243 = qJ(5) + qJ(6);
t242 = cos(t243);
t241 = sin(t243);
t235 = t245 * t251 - t247 * t262;
t234 = t236 * t247 - t251 * t261;
t233 = t238 * t247 + t251 * t263;
t1 = [0, -g(1) * t253 - g(2) * t249, t256, 0, 0, 0, 0, 0, -t254, -t255, -g(3) * t245 - t256 * t244, t254, t255, -g(1) * (t253 * pkin(1) + t239 * pkin(2) + pkin(8) * t263 + t238 * qJ(3)) - g(2) * (t249 * pkin(1) + t237 * pkin(2) - pkin(8) * t261 + t236 * qJ(3)) - g(3) * (t245 * pkin(8) + pkin(7) + (pkin(2) * t248 - qJ(3) * t252) * t244) 0, 0, 0, 0, 0, -g(1) * t233 - g(2) * t234 - g(3) * t235, -g(1) * (t238 * t251 - t247 * t263) - g(2) * (t236 * t251 + t247 * t261) - g(3) * (-t245 * t247 - t251 * t262) 0, 0, 0, 0, 0, -g(1) * (t233 * t250 + t239 * t246) - g(2) * (t234 * t250 + t237 * t246) - g(3) * (t235 * t250 + t246 * t264) -g(1) * (-t233 * t246 + t239 * t250) - g(2) * (-t234 * t246 + t237 * t250) - g(3) * (-t235 * t246 + t250 * t264) 0, 0, 0, 0, 0, -g(1) * (t233 * t242 + t239 * t241) - g(2) * (t234 * t242 + t237 * t241) - g(3) * (t235 * t242 + t241 * t264) -g(1) * (-t233 * t241 + t239 * t242) - g(2) * (-t234 * t241 + t237 * t242) - g(3) * (-t235 * t241 + t242 * t264);];
U_reg  = t1;
