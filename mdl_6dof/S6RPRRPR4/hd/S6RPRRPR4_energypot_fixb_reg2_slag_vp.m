% Calculate inertial parameters regressor of potential energy for
% S6RPRRPR4
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [11x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d3,d4,d6,theta2,theta5]';
% 
% Output:
% U_reg [1x(6*10)]
%   inertial parameter regressor of Potential energy

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 05:11
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U_reg = S6RPRRPR4_energypot_fixb_reg2_slag_vp(qJ, g, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(11,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRRPR4_energypot_fixb_reg2_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RPRRPR4_energypot_fixb_reg2_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6RPRRPR4_energypot_fixb_reg2_slag_vp: pkin has to be [11x1] (double)');

%% Symbolic Calculation
% From energy_potential_fixb_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 05:10:15
% EndTime: 2019-03-09 05:10:15
% DurationCPUTime: 0.15s
% Computational Cost: add. (211->72), mult. (174->92), div. (0->0), fcn. (166->12), ass. (0->42)
t90 = pkin(10) + qJ(3);
t84 = qJ(4) + t90;
t76 = sin(t84);
t77 = cos(t84);
t93 = cos(pkin(11));
t78 = t93 * pkin(5) + pkin(4);
t95 = -pkin(9) - qJ(5);
t117 = -t76 * t95 + t77 * t78;
t116 = g(3) * pkin(6);
t115 = g(3) * t76;
t92 = sin(pkin(10));
t114 = t92 * pkin(2) + pkin(6);
t94 = cos(pkin(10));
t79 = t94 * pkin(2) + pkin(1);
t89 = pkin(11) + qJ(6);
t80 = sin(t89);
t97 = sin(qJ(1));
t111 = t97 * t80;
t82 = cos(t89);
t110 = t97 * t82;
t91 = sin(pkin(11));
t109 = t97 * t91;
t108 = t97 * t93;
t98 = cos(qJ(1));
t107 = t98 * t80;
t106 = t98 * t82;
t105 = t98 * t91;
t104 = t98 * t93;
t96 = -pkin(7) - qJ(2);
t83 = cos(t90);
t71 = pkin(3) * t83 + t79;
t88 = -pkin(8) + t96;
t103 = t97 * t71 + t98 * t88;
t81 = sin(t90);
t102 = pkin(3) * t81 + t114;
t70 = t98 * t71;
t101 = -t97 * t88 + t70;
t100 = g(1) * t98 + g(2) * t97;
t99 = pkin(4) * t77 + qJ(5) * t76;
t72 = g(1) * t97 - g(2) * t98;
t68 = -g(3) * t77 + t100 * t76;
t1 = [0, 0, 0, 0, 0, 0, -t100, t72, -g(3), -t116, 0, 0, 0, 0, 0, 0, -g(3) * t92 - t100 * t94, -g(3) * t94 + t100 * t92, -t72, -g(1) * (t98 * pkin(1) + t97 * qJ(2)) - g(2) * (t97 * pkin(1) - t98 * qJ(2)) - t116, 0, 0, 0, 0, 0, 0, -g(3) * t81 - t100 * t83, -g(3) * t83 + t100 * t81, -t72, -g(1) * (t98 * t79 - t97 * t96) - g(2) * (t97 * t79 + t98 * t96) - g(3) * t114, 0, 0, 0, 0, 0, 0, -t100 * t77 - t115, t68, -t72, -g(1) * t101 - g(2) * t103 - g(3) * t102, 0, 0, 0, 0, 0, 0, -g(1) * (t77 * t104 + t109) - g(2) * (t77 * t108 - t105) - t93 * t115, -g(1) * (-t77 * t105 + t108) - g(2) * (-t77 * t109 - t104) + t91 * t115, -t68, -g(1) * (t99 * t98 + t101) - g(2) * (t99 * t97 + t103) - g(3) * (t76 * pkin(4) - t77 * qJ(5) + t102) 0, 0, 0, 0, 0, 0, -g(1) * (t77 * t106 + t111) - g(2) * (t77 * t110 - t107) - t82 * t115, -g(1) * (-t77 * t107 + t110) - g(2) * (-t77 * t111 - t106) + t80 * t115, -t68, -g(1) * (t117 * t98 + t70) - g(2) * (-pkin(5) * t105 + t103) - g(3) * (t76 * t78 + t77 * t95 + t102) + (-g(1) * (pkin(5) * t91 - t88) - g(2) * t117) * t97;];
U_reg  = t1;
