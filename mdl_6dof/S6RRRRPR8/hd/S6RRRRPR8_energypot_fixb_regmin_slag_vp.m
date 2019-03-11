% Calculate minimal parameter regressor of potential energy for
% S6RRRRPR8
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d3,d4,d6]';
% 
% Output:
% U_reg [1x35]
%   minimal parameter regressor of Potential energy

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 22:46
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U_reg = S6RRRRPR8_energypot_fixb_regmin_slag_vp(qJ, g, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(10,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRRPR8_energypot_fixb_regmin_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RRRRPR8_energypot_fixb_regmin_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RRRRPR8_energypot_fixb_regmin_slag_vp: pkin has to be [10x1] (double)');

%% Symbolic Calculation
% From energy_potential_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 22:44:47
% EndTime: 2019-03-09 22:44:48
% DurationCPUTime: 0.10s
% Computational Cost: add. (128->51), mult. (169->76), div. (0->0), fcn. (196->10), ass. (0->30)
t227 = sin(qJ(2));
t243 = g(3) * t227;
t228 = sin(qJ(1));
t231 = cos(qJ(2));
t242 = t228 * t231;
t224 = qJ(3) + qJ(4);
t222 = sin(t224);
t232 = cos(qJ(1));
t241 = t232 * t222;
t223 = cos(t224);
t240 = t232 * t223;
t226 = sin(qJ(3));
t239 = t232 * t226;
t230 = cos(qJ(3));
t238 = t232 * t230;
t237 = pkin(3) * t226 + pkin(7);
t236 = g(1) * t232 + g(2) * t228;
t221 = t230 * pkin(3) + pkin(2);
t233 = -pkin(9) - pkin(8);
t235 = t221 * t231 - t227 * t233 + pkin(1);
t215 = t222 * t242 + t240;
t217 = -t228 * t223 + t231 * t241;
t234 = g(1) * t217 + g(2) * t215 + t222 * t243;
t229 = cos(qJ(6));
t225 = sin(qJ(6));
t219 = -g(3) * t231 + t236 * t227;
t218 = t228 * t222 + t231 * t240;
t216 = t223 * t242 - t241;
t214 = -g(1) * t218 - g(2) * t216 - t223 * t243;
t1 = [0, -t236, g(1) * t228 - g(2) * t232, 0, 0, 0, 0, 0, -t236 * t231 - t243, t219, 0, 0, 0, 0, 0, -g(1) * (t228 * t226 + t231 * t238) - g(2) * (t230 * t242 - t239) - t230 * t243, -g(1) * (t228 * t230 - t231 * t239) - g(2) * (-t226 * t242 - t238) + t226 * t243, 0, 0, 0, 0, 0, t214, t234, t214, -t219, -t234, -g(1) * (t218 * pkin(4) + t217 * qJ(5)) - g(2) * (t216 * pkin(4) + t215 * qJ(5)) - g(3) * (t231 * t233 + pkin(6)) - (pkin(4) * t223 + qJ(5) * t222 + t221) * t243 + (-g(1) * t235 + g(2) * t237) * t232 + (-g(1) * t237 - g(2) * t235) * t228, 0, 0, 0, 0, 0, -g(1) * (t217 * t225 + t218 * t229) - g(2) * (t215 * t225 + t216 * t229) - (t222 * t225 + t223 * t229) * t243, -g(1) * (t217 * t229 - t218 * t225) - g(2) * (t215 * t229 - t216 * t225) - (t222 * t229 - t223 * t225) * t243;];
U_reg  = t1;
