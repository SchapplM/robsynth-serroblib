% Calculate minimal parameter regressor of potential energy for
% S6RRRRPP4
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d3,d4,theta5]';
% 
% Output:
% U_reg [1x30]
%   minimal parameter regressor of Potential energy

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 21:04
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U_reg = S6RRRRPP4_energypot_fixb_regmin_slag_vp(qJ, g, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(10,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRRPP4_energypot_fixb_regmin_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RRRRPP4_energypot_fixb_regmin_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RRRRPP4_energypot_fixb_regmin_slag_vp: pkin has to be [10x1] (double)');

%% Symbolic Calculation
% From energy_potential_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 21:02:34
% EndTime: 2019-03-09 21:02:34
% DurationCPUTime: 0.16s
% Computational Cost: add. (156->55), mult. (162->77), div. (0->0), fcn. (172->10), ass. (0->34)
t231 = qJ(3) + qJ(4);
t224 = sin(t231);
t232 = sin(qJ(3));
t252 = t232 * pkin(3) + pkin(4) * t224 + pkin(7);
t230 = -qJ(5) - pkin(9) - pkin(8);
t233 = sin(qJ(2));
t251 = -t230 * t233 + pkin(1);
t250 = g(3) * t233;
t234 = sin(qJ(1));
t236 = cos(qJ(2));
t248 = t234 * t236;
t223 = pkin(10) + t231;
t221 = sin(t223);
t237 = cos(qJ(1));
t247 = t237 * t221;
t222 = cos(t223);
t246 = t237 * t222;
t245 = t237 * t224;
t225 = cos(t231);
t244 = t237 * t225;
t243 = t237 * t232;
t235 = cos(qJ(3));
t242 = t237 * t235;
t218 = t235 * pkin(3) + pkin(4) * t225 + pkin(2);
t241 = t233 * t218 + t236 * t230 + pkin(6);
t240 = g(1) * t237 + g(2) * t234;
t239 = t252 * t234 + (t218 * t236 + t251) * t237;
t238 = t218 * t248 + t251 * t234 - t237 * t252;
t213 = -g(3) * t236 + t240 * t233;
t212 = t234 * t221 + t236 * t246;
t211 = -t234 * t222 + t236 * t247;
t210 = t222 * t248 - t247;
t209 = t221 * t248 + t246;
t1 = [0, -t240, g(1) * t234 - g(2) * t237, 0, 0, 0, 0, 0, -t240 * t236 - t250, t213, 0, 0, 0, 0, 0, -g(1) * (t234 * t232 + t236 * t242) - g(2) * (t235 * t248 - t243) - t235 * t250, -g(1) * (t234 * t235 - t236 * t243) - g(2) * (-t232 * t248 - t242) + t232 * t250, 0, 0, 0, 0, 0, -g(1) * (t234 * t224 + t236 * t244) - g(2) * (t225 * t248 - t245) - t225 * t250, -g(1) * (t234 * t225 - t236 * t245) - g(2) * (-t224 * t248 - t244) + t224 * t250, -t213, -g(1) * t239 - g(2) * t238 - g(3) * t241, -g(1) * t212 - g(2) * t210 - t222 * t250, -t213, -g(1) * t211 - g(2) * t209 - t221 * t250, -g(1) * (t212 * pkin(5) + t211 * qJ(6) + t239) - g(2) * (t210 * pkin(5) + t209 * qJ(6) + t238) - g(3) * ((pkin(5) * t222 + qJ(6) * t221) * t233 + t241);];
U_reg  = t1;
