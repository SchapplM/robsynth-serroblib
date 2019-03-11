% Calculate inertial parameters regressor of potential energy for
% S6RRPRPR1
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [11x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d4,d6,theta3,theta5]';
% 
% Output:
% U_reg [1x(6*10)]
%   inertial parameter regressor of Potential energy

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 10:11
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U_reg = S6RRPRPR1_energypot_fixb_reg2_slag_vp(qJ, g, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(11,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPRPR1_energypot_fixb_reg2_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RRPRPR1_energypot_fixb_reg2_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6RRPRPR1_energypot_fixb_reg2_slag_vp: pkin has to be [11x1] (double)');

%% Symbolic Calculation
% From energy_potential_fixb_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 10:09:56
% EndTime: 2019-03-09 10:09:56
% DurationCPUTime: 0.15s
% Computational Cost: add. (211->72), mult. (174->92), div. (0->0), fcn. (166->12), ass. (0->42)
t94 = qJ(2) + pkin(10);
t88 = qJ(4) + t94;
t80 = sin(t88);
t81 = cos(t88);
t96 = cos(pkin(11));
t82 = t96 * pkin(5) + pkin(4);
t97 = -pkin(9) - qJ(5);
t121 = -t80 * t97 + t81 * t82;
t120 = g(3) * pkin(6);
t119 = g(3) * t80;
t99 = sin(qJ(2));
t118 = t99 * pkin(2) + pkin(6);
t101 = cos(qJ(2));
t83 = t101 * pkin(2) + pkin(1);
t98 = -pkin(7) - qJ(3);
t100 = sin(qJ(1));
t102 = cos(qJ(1));
t87 = cos(t94);
t75 = pkin(3) * t87 + t83;
t92 = -pkin(8) + t98;
t115 = t100 * t75 + t102 * t92;
t93 = pkin(11) + qJ(6);
t84 = sin(t93);
t114 = t100 * t84;
t86 = cos(t93);
t113 = t100 * t86;
t95 = sin(pkin(11));
t112 = t100 * t95;
t111 = t100 * t96;
t110 = t102 * t84;
t109 = t102 * t86;
t108 = t102 * t95;
t107 = t102 * t96;
t85 = sin(t94);
t106 = pkin(3) * t85 + t118;
t74 = t102 * t75;
t105 = -t100 * t92 + t74;
t104 = pkin(4) * t81 + qJ(5) * t80;
t103 = g(1) * t102 + g(2) * t100;
t76 = g(1) * t100 - g(2) * t102;
t72 = -g(3) * t81 + t103 * t80;
t1 = [0, 0, 0, 0, 0, 0, -t103, t76, -g(3), -t120, 0, 0, 0, 0, 0, 0, -g(3) * t99 - t103 * t101, -g(3) * t101 + t103 * t99, -t76, -g(1) * (t102 * pkin(1) + t100 * pkin(7)) - g(2) * (t100 * pkin(1) - t102 * pkin(7)) - t120, 0, 0, 0, 0, 0, 0, -g(3) * t85 - t103 * t87, -g(3) * t87 + t103 * t85, -t76, -g(1) * (-t100 * t98 + t102 * t83) - g(2) * (t100 * t83 + t102 * t98) - g(3) * t118, 0, 0, 0, 0, 0, 0, -t103 * t81 - t119, t72, -t76, -g(1) * t105 - g(2) * t115 - g(3) * t106, 0, 0, 0, 0, 0, 0, -g(1) * (t81 * t107 + t112) - g(2) * (t81 * t111 - t108) - t96 * t119, -g(1) * (-t81 * t108 + t111) - g(2) * (-t81 * t112 - t107) + t95 * t119, -t72, -g(1) * (t104 * t102 + t105) - g(2) * (t104 * t100 + t115) - g(3) * (t80 * pkin(4) - t81 * qJ(5) + t106) 0, 0, 0, 0, 0, 0, -g(1) * (t81 * t109 + t114) - g(2) * (t81 * t113 - t110) - t86 * t119, -g(1) * (-t81 * t110 + t113) - g(2) * (-t81 * t114 - t109) + t84 * t119, -t72, -g(1) * (t121 * t102 + t74) - g(2) * (-pkin(5) * t108 + t115) - g(3) * (t80 * t82 + t81 * t97 + t106) + (-g(1) * (pkin(5) * t95 - t92) - g(2) * t121) * t100;];
U_reg  = t1;
