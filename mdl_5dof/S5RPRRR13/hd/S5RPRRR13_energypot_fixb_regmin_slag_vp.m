% Calculate minimal parameter regressor of potential energy for
% S5RPRRR13
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d3,d4,d5]';
% 
% Output:
% U_reg [1x27]
%   minimal parameter regressor of Potential energy

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 19:16
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U_reg = S5RPRRR13_energypot_fixb_regmin_slag_vp(qJ, g, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(3,1),zeros(8,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPRRR13_energypot_fixb_regmin_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RPRRR13_energypot_fixb_regmin_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RPRRR13_energypot_fixb_regmin_slag_vp: pkin has to be [8x1] (double)');

%% Symbolic Calculation
% From energy_potential_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 19:15:27
% EndTime: 2019-12-31 19:15:27
% DurationCPUTime: 0.05s
% Computational Cost: add. (38->27), mult. (65->43), div. (0->0), fcn. (70->8), ass. (0->21)
t139 = cos(qJ(3));
t149 = g(3) * t139;
t134 = qJ(4) + qJ(5);
t132 = sin(t134);
t137 = sin(qJ(1));
t148 = t137 * t132;
t133 = cos(t134);
t147 = t137 * t133;
t135 = sin(qJ(4));
t146 = t137 * t135;
t138 = cos(qJ(4));
t145 = t137 * t138;
t140 = cos(qJ(1));
t144 = t140 * t132;
t143 = t140 * t133;
t142 = t140 * t135;
t141 = t140 * t138;
t130 = g(1) * t137 - g(2) * t140;
t136 = sin(qJ(3));
t131 = g(1) * t140 + g(2) * t137;
t1 = [0, -t131, t130, t131, -t130, -g(1) * (t140 * pkin(1) + t137 * qJ(2)) - g(2) * (t137 * pkin(1) - t140 * qJ(2)) - g(3) * pkin(5), 0, 0, 0, 0, 0, -t130 * t136 - t149, g(3) * t136 - t130 * t139, 0, 0, 0, 0, 0, -g(1) * (t136 * t145 + t142) - g(2) * (-t136 * t141 + t146) - t138 * t149, -g(1) * (-t136 * t146 + t141) - g(2) * (t136 * t142 + t145) + t135 * t149, 0, 0, 0, 0, 0, -g(1) * (t136 * t147 + t144) - g(2) * (-t136 * t143 + t148) - t133 * t149, -g(1) * (-t136 * t148 + t143) - g(2) * (t136 * t144 + t147) + t132 * t149;];
U_reg = t1;
