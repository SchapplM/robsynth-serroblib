% Calculate minimal parameter regressor of potential energy for
% S6RPRRRR4
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
% Datum: 2019-03-09 07:07
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U_reg = S6RPRRRR4_energypot_fixb_regmin_slag_vp(qJ, g, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(11,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRRRR4_energypot_fixb_regmin_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RPRRRR4_energypot_fixb_regmin_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6RPRRRR4_energypot_fixb_regmin_slag_vp: pkin has to be [11x1] (double)');

%% Symbolic Calculation
% From energy_potential_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 07:05:53
% EndTime: 2019-03-09 07:05:53
% DurationCPUTime: 0.06s
% Computational Cost: add. (85->27), mult. (73->41), div. (0->0), fcn. (74->12), ass. (0->23)
t213 = pkin(11) + qJ(3);
t212 = qJ(4) + t213;
t209 = qJ(5) + t212;
t205 = sin(t209);
t225 = g(3) * t205;
t216 = sin(qJ(6));
t217 = sin(qJ(1));
t224 = t217 * t216;
t218 = cos(qJ(6));
t223 = t217 * t218;
t219 = cos(qJ(1));
t222 = t219 * t216;
t221 = t219 * t218;
t220 = g(1) * t219 + g(2) * t217;
t215 = cos(pkin(11));
t214 = sin(pkin(11));
t211 = cos(t213);
t210 = sin(t213);
t208 = cos(t212);
t207 = sin(t212);
t206 = cos(t209);
t204 = g(1) * t217 - g(2) * t219;
t1 = [0, -t220, t204, -g(3) * t214 - t220 * t215, -g(3) * t215 + t220 * t214, -t204, -g(1) * (t219 * pkin(1) + t217 * qJ(2)) - g(2) * (t217 * pkin(1) - t219 * qJ(2)) - g(3) * pkin(6), 0, 0, 0, 0, 0, -g(3) * t210 - t220 * t211, -g(3) * t211 + t220 * t210, 0, 0, 0, 0, 0, -g(3) * t207 - t220 * t208, -g(3) * t208 + t220 * t207, 0, 0, 0, 0, 0, -t220 * t206 - t225, -g(3) * t206 + t220 * t205, 0, 0, 0, 0, 0, -g(1) * (t206 * t221 + t224) - g(2) * (t206 * t223 - t222) - t218 * t225, -g(1) * (-t206 * t222 + t223) - g(2) * (-t206 * t224 - t221) + t216 * t225;];
U_reg  = t1;
