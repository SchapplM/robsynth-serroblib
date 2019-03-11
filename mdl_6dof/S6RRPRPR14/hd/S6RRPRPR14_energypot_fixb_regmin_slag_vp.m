% Calculate minimal parameter regressor of potential energy for
% S6RRPRPR14
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d1,d2,d4,d6]';
% 
% Output:
% U_reg [1x32]
%   minimal parameter regressor of Potential energy

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 11:38
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U_reg = S6RRPRPR14_energypot_fixb_regmin_slag_vp(qJ, g, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(10,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPRPR14_energypot_fixb_regmin_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RRPRPR14_energypot_fixb_regmin_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RRPRPR14_energypot_fixb_regmin_slag_vp: pkin has to be [10x1] (double)');

%% Symbolic Calculation
% From energy_potential_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 11:36:56
% EndTime: 2019-03-09 11:36:56
% DurationCPUTime: 0.16s
% Computational Cost: add. (130->63), mult. (304->92), div. (0->0), fcn. (368->10), ass. (0->39)
t250 = sin(qJ(1));
t272 = g(1) * t250;
t254 = cos(qJ(1));
t271 = g(2) * t254;
t245 = sin(pkin(6));
t249 = sin(qJ(2));
t270 = t245 * t249;
t269 = t245 * t250;
t253 = cos(qJ(2));
t268 = t245 * t253;
t267 = t245 * t254;
t266 = t250 * t249;
t265 = t250 * t253;
t264 = t254 * t249;
t263 = t254 * t253;
t246 = cos(pkin(6));
t262 = pkin(2) * t270 + t246 * pkin(8) + pkin(7);
t261 = -t271 + t272;
t234 = -t246 * t263 + t266;
t235 = t246 * t264 + t265;
t260 = t250 * pkin(1) + t235 * pkin(2) + t234 * qJ(3);
t236 = t246 * t265 + t264;
t237 = -t246 * t266 + t263;
t259 = t254 * pkin(1) + t237 * pkin(2) + pkin(8) * t269 + t236 * qJ(3);
t248 = sin(qJ(4));
t252 = cos(qJ(4));
t224 = -t236 * t252 + t248 * t269;
t226 = t234 * t252 + t248 * t267;
t232 = t246 * t248 + t252 * t268;
t258 = g(1) * t224 - g(2) * t226 + g(3) * t232;
t225 = t236 * t248 + t252 * t269;
t227 = -t234 * t248 + t252 * t267;
t233 = t246 * t252 - t248 * t268;
t257 = g(1) * t225 - g(2) * t227 + g(3) * t233;
t256 = -g(1) * t236 - g(2) * t234 + g(3) * t268;
t255 = g(1) * t237 + g(2) * t235 + g(3) * t270;
t251 = cos(qJ(6));
t247 = sin(qJ(6));
t1 = [0, -g(1) * t254 - g(2) * t250, t261, 0, 0, 0, 0, 0, -t255, -t256, -g(3) * t246 - t261 * t245, t255, t256, -g(1) * t259 - g(2) * (-pkin(8) * t267 + t260) - g(3) * (-qJ(3) * t268 + t262) 0, 0, 0, 0, 0, -t257, t258, -t255, t257, -t258, -g(1) * (t225 * pkin(4) + t237 * pkin(9) + t224 * qJ(5) + t259) - g(2) * (-t227 * pkin(4) + t235 * pkin(9) - t226 * qJ(5) + t260) - g(3) * (t246 * pkin(3) + t233 * pkin(4) + t232 * qJ(5) + t262) + (-pkin(3) * t272 - g(3) * (pkin(9) * t249 - qJ(3) * t253) - (-pkin(3) - pkin(8)) * t271) * t245, 0, 0, 0, 0, 0, -g(1) * (t224 * t247 + t237 * t251) - g(2) * (-t226 * t247 + t235 * t251) - g(3) * (t232 * t247 + t251 * t270) -g(1) * (t224 * t251 - t237 * t247) - g(2) * (-t226 * t251 - t235 * t247) - g(3) * (t232 * t251 - t247 * t270);];
U_reg  = t1;
