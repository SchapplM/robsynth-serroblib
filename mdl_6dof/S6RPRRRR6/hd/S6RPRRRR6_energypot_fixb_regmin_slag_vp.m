% Calculate minimal parameter regressor of potential energy for
% S6RPRRRR6
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [11x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d3,d4,d5,d6,theta2]';
% 
% Output:
% U_reg [1x35]
%   minimal parameter regressor of Potential energy

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 07:16
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U_reg = S6RPRRRR6_energypot_fixb_regmin_slag_vp(qJ, g, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(11,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRRRR6_energypot_fixb_regmin_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RPRRRR6_energypot_fixb_regmin_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6RPRRRR6_energypot_fixb_regmin_slag_vp: pkin has to be [11x1] (double)');

%% Symbolic Calculation
% From energy_potential_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 07:14:32
% EndTime: 2019-03-09 07:14:33
% DurationCPUTime: 0.10s
% Computational Cost: add. (93->39), mult. (93->61), div. (0->0), fcn. (102->12), ass. (0->31)
t223 = pkin(11) + qJ(3);
t218 = sin(t223);
t244 = g(3) * t218;
t224 = qJ(4) + qJ(5);
t222 = qJ(6) + t224;
t216 = sin(t222);
t228 = sin(qJ(1));
t243 = t228 * t216;
t217 = cos(t222);
t242 = t228 * t217;
t220 = sin(t224);
t241 = t228 * t220;
t221 = cos(t224);
t240 = t228 * t221;
t227 = sin(qJ(4));
t239 = t228 * t227;
t229 = cos(qJ(4));
t238 = t228 * t229;
t230 = cos(qJ(1));
t237 = t230 * t216;
t236 = t230 * t217;
t235 = t230 * t220;
t234 = t230 * t221;
t233 = t230 * t227;
t232 = t230 * t229;
t231 = g(1) * t230 + g(2) * t228;
t226 = cos(pkin(11));
t225 = sin(pkin(11));
t219 = cos(t223);
t215 = g(1) * t228 - g(2) * t230;
t1 = [0, -t231, t215, -g(3) * t225 - t231 * t226, -g(3) * t226 + t231 * t225, -t215, -g(1) * (t230 * pkin(1) + t228 * qJ(2)) - g(2) * (t228 * pkin(1) - t230 * qJ(2)) - g(3) * pkin(6), 0, 0, 0, 0, 0, -t231 * t219 - t244, -g(3) * t219 + t231 * t218, 0, 0, 0, 0, 0, -g(1) * (t219 * t232 + t239) - g(2) * (t219 * t238 - t233) - t229 * t244, -g(1) * (-t219 * t233 + t238) - g(2) * (-t219 * t239 - t232) + t227 * t244, 0, 0, 0, 0, 0, -g(1) * (t219 * t234 + t241) - g(2) * (t219 * t240 - t235) - t221 * t244, -g(1) * (-t219 * t235 + t240) - g(2) * (-t219 * t241 - t234) + t220 * t244, 0, 0, 0, 0, 0, -g(1) * (t219 * t236 + t243) - g(2) * (t219 * t242 - t237) - t217 * t244, -g(1) * (-t219 * t237 + t242) - g(2) * (-t219 * t243 - t236) + t216 * t244;];
U_reg  = t1;
