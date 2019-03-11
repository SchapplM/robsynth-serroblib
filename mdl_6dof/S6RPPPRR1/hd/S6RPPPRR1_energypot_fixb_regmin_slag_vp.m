% Calculate minimal parameter regressor of potential energy for
% S6RPPPRR1
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d5,d6,theta2]';
% 
% Output:
% U_reg [1x24]
%   minimal parameter regressor of Potential energy

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 01:30
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U_reg = S6RPPPRR1_energypot_fixb_regmin_slag_vp(qJ, g, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(9,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPPPRR1_energypot_fixb_regmin_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RPPPRR1_energypot_fixb_regmin_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S6RPPPRR1_energypot_fixb_regmin_slag_vp: pkin has to be [9x1] (double)');

%% Symbolic Calculation
% From energy_potential_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 01:30:13
% EndTime: 2019-03-09 01:30:13
% DurationCPUTime: 0.08s
% Computational Cost: add. (68->29), mult. (67->43), div. (0->0), fcn. (62->8), ass. (0->20)
t131 = qJ(2) + pkin(6);
t144 = g(3) * t131;
t136 = cos(qJ(5));
t143 = g(3) * t136;
t132 = sin(qJ(6));
t133 = sin(qJ(5));
t142 = t132 * t133;
t135 = cos(qJ(6));
t141 = t133 * t135;
t130 = qJ(1) + pkin(9);
t126 = sin(t130);
t127 = cos(t130);
t137 = cos(qJ(1));
t140 = pkin(1) * t137 + pkin(2) * t127 + qJ(3) * t126;
t122 = g(1) * t127 + g(2) * t126;
t134 = sin(qJ(1));
t139 = -g(1) * t137 - g(2) * t134;
t138 = pkin(1) * t134 + pkin(2) * t126 - qJ(3) * t127;
t121 = -g(1) * t126 + g(2) * t127;
t1 = [0, t139, g(1) * t134 - g(2) * t137, pkin(1) * t139 - t144, t122, t121, -g(1) * t140 - g(2) * t138 - t144, t121, -t122, -g(1) * (qJ(4) * t127 + t140) - g(2) * (qJ(4) * t126 + t138) - g(3) * (pkin(3) + t131) 0, 0, 0, 0, 0, -t122 * t133 - t143, g(3) * t133 - t122 * t136, 0, 0, 0, 0, 0, -g(1) * (-t126 * t132 + t127 * t141) - g(2) * (t126 * t141 + t127 * t132) - t135 * t143, -g(1) * (-t126 * t135 - t127 * t142) - g(2) * (-t126 * t142 + t127 * t135) + t132 * t143;];
U_reg  = t1;
