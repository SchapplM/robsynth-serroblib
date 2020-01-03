% Calculate minimal parameter regressor of potential energy for
% S5RPRPR9
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d3,d5,theta2]';
% 
% Output:
% U_reg [1x22]
%   minimal parameter regressor of Potential energy

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 18:25
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U_reg = S5RPRPR9_energypot_fixb_regmin_slag_vp(qJ, g, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(3,1),zeros(8,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPRPR9_energypot_fixb_regmin_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RPRPR9_energypot_fixb_regmin_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RPRPR9_energypot_fixb_regmin_slag_vp: pkin has to be [8x1] (double)');

%% Symbolic Calculation
% From energy_potential_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 18:24:37
% EndTime: 2019-12-31 18:24:37
% DurationCPUTime: 0.05s
% Computational Cost: add. (61->28), mult. (70->41), div. (0->0), fcn. (68->8), ass. (0->20)
t135 = cos(qJ(3));
t144 = g(3) * t135;
t143 = qJ(2) + pkin(5);
t131 = sin(qJ(5));
t132 = sin(qJ(3));
t142 = t131 * t132;
t134 = cos(qJ(5));
t141 = t132 * t134;
t130 = qJ(1) + pkin(8);
t128 = sin(t130);
t129 = cos(t130);
t140 = g(1) * t129 + g(2) * t128;
t133 = sin(qJ(1));
t136 = cos(qJ(1));
t139 = -g(1) * t136 - g(2) * t133;
t138 = pkin(3) * t135 + qJ(4) * t132 + pkin(2);
t137 = t139 * pkin(1);
t127 = g(3) * t132 + t140 * t135;
t126 = t140 * t132 - t144;
t1 = [0, t139, g(1) * t133 - g(2) * t136, -g(3) * t143 + t137, 0, 0, 0, 0, 0, -t127, t126, -g(1) * t128 + g(2) * t129, t127, -t126, -g(3) * (t132 * pkin(3) - t135 * qJ(4) + t143) + t137 + (g(2) * pkin(6) - g(1) * t138) * t129 + (-g(1) * pkin(6) - g(2) * t138) * t128, 0, 0, 0, 0, 0, -g(1) * (t128 * t134 + t129 * t142) - g(2) * (t128 * t142 - t129 * t134) + t131 * t144, -g(1) * (-t128 * t131 + t129 * t141) - g(2) * (t128 * t141 + t129 * t131) + t134 * t144;];
U_reg = t1;
