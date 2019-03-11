% Calculate minimal parameter regressor of potential energy for
% S6RRRPPR4
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d3,d6,theta4]';
% 
% Output:
% U_reg [1x30]
%   minimal parameter regressor of Potential energy

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 15:37
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U_reg = S6RRRPPR4_energypot_fixb_regmin_slag_vp(qJ, g, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(10,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRPPR4_energypot_fixb_regmin_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RRRPPR4_energypot_fixb_regmin_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RRRPPR4_energypot_fixb_regmin_slag_vp: pkin has to be [10x1] (double)');

%% Symbolic Calculation
% From energy_potential_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 15:35:59
% EndTime: 2019-03-09 15:36:00
% DurationCPUTime: 0.17s
% Computational Cost: add. (130->55), mult. (178->80), div. (0->0), fcn. (198->10), ass. (0->31)
t230 = -qJ(4) - pkin(8);
t233 = sin(qJ(2));
t251 = -t230 * t233 + pkin(1);
t250 = g(3) * t233;
t232 = sin(qJ(3));
t234 = sin(qJ(1));
t248 = t234 * t232;
t237 = cos(qJ(2));
t247 = t234 * t237;
t229 = qJ(3) + pkin(10);
t223 = sin(t229);
t238 = cos(qJ(1));
t246 = t238 * t223;
t224 = cos(t229);
t245 = t238 * t224;
t244 = t238 * t232;
t236 = cos(qJ(3));
t243 = t238 * t236;
t222 = t236 * pkin(3) + pkin(2);
t242 = t233 * t222 + t237 * t230 + pkin(6);
t241 = g(1) * t238 + g(2) * t234;
t240 = pkin(3) * t248 + t234 * pkin(7) + (t222 * t237 + t251) * t238;
t239 = t222 * t247 + (-pkin(3) * t232 - pkin(7)) * t238 + t251 * t234;
t235 = cos(qJ(6));
t231 = sin(qJ(6));
t216 = -g(3) * t237 + t241 * t233;
t215 = t234 * t223 + t237 * t245;
t214 = -t234 * t224 + t237 * t246;
t213 = t224 * t247 - t246;
t212 = t223 * t247 + t245;
t1 = [0, -t241, g(1) * t234 - g(2) * t238, 0, 0, 0, 0, 0, -t241 * t237 - t250, t216, 0, 0, 0, 0, 0, -g(1) * (t237 * t243 + t248) - g(2) * (t236 * t247 - t244) - t236 * t250, -g(1) * (t234 * t236 - t237 * t244) - g(2) * (-t232 * t247 - t243) + t232 * t250, -t216, -g(1) * t240 - g(2) * t239 - g(3) * t242, -g(1) * t215 - g(2) * t213 - t224 * t250, -t216, -g(1) * t214 - g(2) * t212 - t223 * t250, -g(1) * (t215 * pkin(4) + t214 * qJ(5) + t240) - g(2) * (t213 * pkin(4) + t212 * qJ(5) + t239) - g(3) * ((pkin(4) * t224 + qJ(5) * t223) * t233 + t242) 0, 0, 0, 0, 0, -g(1) * (t214 * t231 + t215 * t235) - g(2) * (t212 * t231 + t213 * t235) - (t223 * t231 + t224 * t235) * t250, -g(1) * (t214 * t235 - t215 * t231) - g(2) * (t212 * t235 - t213 * t231) - (t223 * t235 - t224 * t231) * t250;];
U_reg  = t1;
