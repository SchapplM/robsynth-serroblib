% Calculate minimal parameter regressor of potential energy for
% S6RRPRPP1
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d4,theta3,theta5]';
% 
% Output:
% U_reg [1x25]
%   minimal parameter regressor of Potential energy

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 09:48
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U_reg = S6RRPRPP1_energypot_fixb_regmin_slag_vp(qJ, g, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(10,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPRPP1_energypot_fixb_regmin_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RRPRPP1_energypot_fixb_regmin_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RRPRPP1_energypot_fixb_regmin_slag_vp: pkin has to be [10x1] (double)');

%% Symbolic Calculation
% From energy_potential_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 09:47:51
% EndTime: 2019-03-09 09:47:51
% DurationCPUTime: 0.14s
% Computational Cost: add. (148->51), mult. (150->70), div. (0->0), fcn. (153->10), ass. (0->40)
t237 = cos(qJ(4));
t223 = t237 * pkin(4) + pkin(3);
t231 = qJ(2) + pkin(9);
t226 = sin(t231);
t228 = cos(t231);
t232 = -qJ(5) - pkin(8);
t258 = t223 * t228 - t226 * t232;
t257 = g(3) * t226;
t235 = sin(qJ(2));
t256 = t235 * pkin(2) + pkin(6);
t230 = qJ(4) + pkin(10);
t225 = sin(t230);
t236 = sin(qJ(1));
t253 = t236 * t225;
t227 = cos(t230);
t252 = t236 * t227;
t234 = sin(qJ(4));
t251 = t236 * t234;
t250 = t236 * t237;
t239 = cos(qJ(1));
t249 = t239 * t225;
t248 = t239 * t227;
t247 = t239 * t234;
t246 = t239 * t237;
t238 = cos(qJ(2));
t224 = t238 * pkin(2) + pkin(1);
t233 = -pkin(7) - qJ(3);
t245 = t236 * t224 + t239 * t233;
t244 = t239 * t224 - t236 * t233;
t243 = t226 * t223 + t228 * t232 + t256;
t242 = g(1) * t239 + g(2) * t236;
t241 = pkin(4) * t251 + t258 * t239 + t244;
t240 = -pkin(4) * t247 + t258 * t236 + t245;
t217 = g(1) * t236 - g(2) * t239;
t213 = t228 * t248 + t253;
t212 = t228 * t249 - t252;
t211 = t228 * t252 - t249;
t210 = t228 * t253 + t248;
t209 = g(3) * t228 - t242 * t226;
t1 = [0, -t242, t217, 0, 0, 0, 0, 0, -g(3) * t235 - t242 * t238, -g(3) * t238 + t242 * t235, -t217, -g(1) * t244 - g(2) * t245 - g(3) * t256, 0, 0, 0, 0, 0, -g(1) * (t228 * t246 + t251) - g(2) * (t228 * t250 - t247) - t237 * t257, -g(1) * (-t228 * t247 + t250) - g(2) * (-t228 * t251 - t246) + t234 * t257, t209, -g(1) * t241 - g(2) * t240 - g(3) * t243, -g(1) * t213 - g(2) * t211 - t227 * t257, t209, -g(1) * t212 - g(2) * t210 - t225 * t257, -g(1) * (t213 * pkin(5) + t212 * qJ(6) + t241) - g(2) * (t211 * pkin(5) + t210 * qJ(6) + t240) - g(3) * ((pkin(5) * t227 + qJ(6) * t225) * t226 + t243);];
U_reg  = t1;
