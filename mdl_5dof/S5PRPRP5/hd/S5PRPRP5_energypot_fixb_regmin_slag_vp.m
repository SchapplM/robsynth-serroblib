% Calculate minimal parameter regressor of potential energy for
% S5PRPRP5
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d2,d4,theta1,theta3]';
% 
% Output:
% U_reg [1x19]
%   minimal parameter regressor of Potential energy

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-05 15:39
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U_reg = S5PRPRP5_energypot_fixb_regmin_slag_vp(qJ, g, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(3,1),zeros(8,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PRPRP5_energypot_fixb_regmin_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5PRPRP5_energypot_fixb_regmin_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5PRPRP5_energypot_fixb_regmin_slag_vp: pkin has to be [8x1] (double)');

%% Symbolic Calculation
% From energy_potential_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 15:38:40
% EndTime: 2019-12-05 15:38:40
% DurationCPUTime: 0.10s
% Computational Cost: add. (101->52), mult. (142->72), div. (0->0), fcn. (151->8), ass. (0->28)
t133 = sin(pkin(7));
t147 = g(2) * t133;
t137 = sin(qJ(2));
t146 = g(3) * t137;
t132 = sin(pkin(8));
t145 = t133 * t132;
t138 = cos(qJ(2));
t144 = t133 * t138;
t135 = cos(pkin(7));
t143 = t135 * t138;
t142 = t135 * pkin(1) + t133 * pkin(5);
t141 = g(1) * t135 + t147;
t140 = pkin(2) * t138 + qJ(3) * t137;
t131 = pkin(8) + qJ(4);
t125 = sin(t131);
t126 = cos(t131);
t118 = t125 * t144 + t135 * t126;
t120 = t125 * t143 - t133 * t126;
t139 = g(1) * t120 + g(2) * t118 + t125 * t146;
t136 = -pkin(6) - qJ(3);
t134 = cos(pkin(8));
t128 = t133 * pkin(1);
t124 = t134 * pkin(3) + pkin(2);
t122 = -g(3) * t138 + t141 * t137;
t121 = t133 * t125 + t126 * t143;
t119 = -t135 * t125 + t126 * t144;
t117 = -g(1) * t121 - g(2) * t119 - t126 * t146;
t1 = [-g(3) * qJ(1), 0, -t141 * t138 - t146, t122, -g(1) * (t134 * t143 + t145) - g(2) * (-t135 * t132 + t134 * t144) - t134 * t146, -g(1) * (-t132 * t143 + t133 * t134) - g(2) * (-t132 * t144 - t135 * t134) + t132 * t146, -t122, -g(1) * (t140 * t135 + t142) - g(2) * (-t135 * pkin(5) + t140 * t133 + t128) - g(3) * (t137 * pkin(2) - t138 * qJ(3) + qJ(1)), 0, 0, 0, 0, 0, t117, t139, t117, -t122, -t139, -g(1) * (pkin(3) * t145 + t121 * pkin(4) + t120 * qJ(5) + t142) - g(2) * (t119 * pkin(4) + t118 * qJ(5) + t124 * t144 + t128) - g(3) * (t138 * t136 + qJ(1)) + (t136 * t147 - g(3) * (pkin(4) * t126 + qJ(5) * t125 + t124)) * t137 + (-g(1) * (t124 * t138 - t136 * t137) - g(2) * (-pkin(3) * t132 - pkin(5))) * t135;];
U_reg = t1;
