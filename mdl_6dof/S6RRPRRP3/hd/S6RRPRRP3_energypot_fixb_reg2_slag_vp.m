% Calculate inertial parameters regressor of potential energy for
% S6RRPRRP3
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d4,d5,theta3]';
% 
% Output:
% U_reg [1x(6*10)]
%   inertial parameter regressor of Potential energy

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 11:51
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U_reg = S6RRPRRP3_energypot_fixb_reg2_slag_vp(qJ, g, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(10,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPRRP3_energypot_fixb_reg2_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RRPRRP3_energypot_fixb_reg2_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RRPRRP3_energypot_fixb_reg2_slag_vp: pkin has to be [10x1] (double)');

%% Symbolic Calculation
% From energy_potential_fixb_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 11:50:33
% EndTime: 2019-03-09 11:50:33
% DurationCPUTime: 0.15s
% Computational Cost: add. (202->77), mult. (200->96), div. (0->0), fcn. (196->10), ass. (0->44)
t117 = g(3) * pkin(6);
t98 = -pkin(9) - pkin(8);
t89 = qJ(2) + pkin(10);
t81 = sin(t89);
t116 = g(3) * t81;
t92 = sin(qJ(4));
t115 = t92 * pkin(4);
t93 = sin(qJ(2));
t114 = t93 * pkin(2) + pkin(6);
t95 = cos(qJ(4));
t79 = t95 * pkin(4) + pkin(3);
t88 = -qJ(6) + t98;
t113 = t81 * t88;
t112 = t81 * t98;
t82 = cos(t89);
t97 = cos(qJ(1));
t111 = t82 * t97;
t90 = qJ(4) + qJ(5);
t83 = sin(t90);
t94 = sin(qJ(1));
t110 = t94 * t83;
t84 = cos(t90);
t109 = t94 * t84;
t108 = t94 * t92;
t107 = t94 * t95;
t106 = t97 * t83;
t105 = t97 * t84;
t104 = t97 * t92;
t103 = t97 * t95;
t96 = cos(qJ(2));
t80 = t96 * pkin(2) + pkin(1);
t91 = -pkin(7) - qJ(3);
t102 = t94 * t80 + t97 * t91;
t76 = t97 * t80;
t101 = -t94 * t91 + t76;
t100 = pkin(3) * t82 + pkin(8) * t81;
t99 = g(1) * t97 + g(2) * t94;
t74 = g(1) * t94 - g(2) * t97;
t73 = pkin(5) * t83 + t115;
t72 = pkin(5) * t84 + t79;
t71 = -g(3) * t82 + t99 * t81;
t70 = -g(1) * (t82 * t105 + t110) - g(2) * (t82 * t109 - t106) - t84 * t116;
t69 = -g(1) * (-t82 * t106 + t109) - g(2) * (-t82 * t110 - t105) + t83 * t116;
t1 = [0, 0, 0, 0, 0, 0, -t99, t74, -g(3), -t117, 0, 0, 0, 0, 0, 0, -g(3) * t93 - t99 * t96, -g(3) * t96 + t99 * t93, -t74, -g(1) * (t97 * pkin(1) + t94 * pkin(7)) - g(2) * (t94 * pkin(1) - t97 * pkin(7)) - t117, 0, 0, 0, 0, 0, 0, -t99 * t82 - t116, t71, -t74, -g(1) * t101 - g(2) * t102 - g(3) * t114, 0, 0, 0, 0, 0, 0, -g(1) * (t82 * t103 + t108) - g(2) * (t82 * t107 - t104) - t95 * t116, -g(1) * (-t82 * t104 + t107) - g(2) * (-t82 * t108 - t103) + t92 * t116, -t71, -g(1) * (t100 * t97 + t101) - g(2) * (t100 * t94 + t102) - g(3) * (t81 * pkin(3) - t82 * pkin(8) + t114) 0, 0, 0, 0, 0, 0, t70, t69, -t71, -g(1) * (t79 * t111 - t97 * t112 + t76) - g(2) * (-pkin(4) * t104 + t102) - g(3) * (t81 * t79 + t82 * t98 + t114) + (-g(1) * (-t91 + t115) - g(2) * (t79 * t82 - t112)) * t94, 0, 0, 0, 0, 0, 0, t70, t69, -t71, -g(1) * (t72 * t111 - t97 * t113 + t76) - g(2) * (-t97 * t73 + t102) - g(3) * (t81 * t72 + t82 * t88 + t114) + (-g(1) * (t73 - t91) - g(2) * (t72 * t82 - t113)) * t94;];
U_reg  = t1;
