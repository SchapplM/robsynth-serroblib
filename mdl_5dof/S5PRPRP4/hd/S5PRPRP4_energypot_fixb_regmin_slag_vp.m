% Calculate minimal parameter regressor of potential energy for
% S5PRPRP4
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
% U_reg [1x16]
%   minimal parameter regressor of Potential energy

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-05 15:36
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U_reg = S5PRPRP4_energypot_fixb_regmin_slag_vp(qJ, g, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(3,1),zeros(8,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PRPRP4_energypot_fixb_regmin_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5PRPRP4_energypot_fixb_regmin_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5PRPRP4_energypot_fixb_regmin_slag_vp: pkin has to be [8x1] (double)');

%% Symbolic Calculation
% From energy_potential_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 15:36:05
% EndTime: 2019-12-05 15:36:05
% DurationCPUTime: 0.07s
% Computational Cost: add. (84->35), mult. (106->49), div. (0->0), fcn. (111->8), ass. (0->28)
t130 = qJ(2) + pkin(8);
t128 = cos(t130);
t148 = pkin(3) * t128;
t127 = sin(t130);
t147 = g(3) * t127;
t131 = sin(pkin(7));
t134 = sin(qJ(4));
t146 = t131 * t134;
t136 = cos(qJ(4));
t145 = t131 * t136;
t132 = cos(pkin(7));
t144 = t132 * t134;
t143 = t132 * t136;
t137 = cos(qJ(2));
t126 = t137 * pkin(2) + pkin(1);
t133 = -qJ(3) - pkin(5);
t142 = t131 * t126 + t132 * t133;
t135 = sin(qJ(2));
t141 = t135 * pkin(2) + qJ(1);
t140 = t132 * t126 - t131 * t133;
t139 = g(1) * t132 + g(2) * t131;
t118 = t128 * t146 + t143;
t120 = t128 * t144 - t145;
t138 = g(1) * t120 + g(2) * t118 + t134 * t147;
t121 = t128 * t143 + t146;
t119 = t128 * t145 - t144;
t117 = -g(1) * t121 - g(2) * t119 - t136 * t147;
t1 = [-g(3) * qJ(1), 0, -g(3) * t135 - t139 * t137, -g(3) * t137 + t139 * t135, -g(1) * t140 - g(2) * t142 - g(3) * t141, 0, 0, 0, 0, 0, t117, t138, t117, g(3) * t128 - t139 * t127, -t138, -g(1) * (t121 * pkin(4) + t120 * qJ(5) + t132 * t148 + t140) - g(2) * (t119 * pkin(4) + t118 * qJ(5) + t131 * t148 + t142) - g(3) * (-t128 * pkin(6) + t141) + (-g(3) * (pkin(4) * t136 + qJ(5) * t134 + pkin(3)) - t139 * pkin(6)) * t127;];
U_reg = t1;
