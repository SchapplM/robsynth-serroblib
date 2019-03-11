% Calculate inertial parameters regressor of potential energy for
% S6RPRRPR8
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d3,d4,d6,theta5]';
% 
% Output:
% U_reg [1x(6*10)]
%   inertial parameter regressor of Potential energy

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 05:25
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U_reg = S6RPRRPR8_energypot_fixb_reg2_slag_vp(qJ, g, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(10,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRRPR8_energypot_fixb_reg2_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RPRRPR8_energypot_fixb_reg2_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RPRRPR8_energypot_fixb_reg2_slag_vp: pkin has to be [10x1] (double)');

%% Symbolic Calculation
% From energy_potential_fixb_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 05:24:50
% EndTime: 2019-03-09 05:24:50
% DurationCPUTime: 0.16s
% Computational Cost: add. (151->82), mult. (188->103), div. (0->0), fcn. (184->10), ass. (0->47)
t106 = g(3) * pkin(6);
t105 = pkin(2) + pkin(6);
t82 = cos(qJ(1));
t104 = g(2) * t82;
t81 = cos(qJ(3));
t103 = g(3) * t81;
t77 = sin(qJ(4));
t102 = t77 * pkin(4);
t80 = cos(qJ(4));
t64 = t80 * pkin(4) + pkin(3);
t78 = sin(qJ(3));
t79 = sin(qJ(1));
t101 = t78 * t79;
t75 = qJ(4) + pkin(10);
t67 = qJ(6) + t75;
t62 = sin(t67);
t100 = t79 * t62;
t63 = cos(t67);
t99 = t79 * t63;
t65 = sin(t75);
t98 = t79 * t65;
t66 = cos(t75);
t97 = t79 * t66;
t96 = t79 * t77;
t95 = t79 * t80;
t94 = t79 * t81;
t93 = t82 * t62;
t92 = t82 * t63;
t91 = t82 * t65;
t90 = t82 * t66;
t89 = t82 * t77;
t88 = t82 * t80;
t76 = -qJ(5) - pkin(8);
t69 = t79 * pkin(7);
t70 = t79 * pkin(1);
t87 = t69 + t70;
t86 = t82 * pkin(1) + t79 * qJ(2);
t85 = t82 * pkin(7) + t86;
t84 = -t82 * qJ(2) + t70;
t83 = pkin(3) * t78 - pkin(8) * t81;
t60 = g(1) * t79 - t104;
t74 = -pkin(9) + t76;
t61 = g(1) * t82 + g(2) * t79;
t59 = pkin(5) * t65 + t102;
t58 = pkin(5) * t66 + t64;
t57 = -g(3) * t78 + t60 * t81;
t1 = [0, 0, 0, 0, 0, 0, -t61, t60, -g(3), -t106, 0, 0, 0, 0, 0, 0, -g(3), t61, -t60, -g(1) * t86 - g(2) * t84 - t106, 0, 0, 0, 0, 0, 0, -t60 * t78 - t103, -t57, -t61, -g(1) * t85 - g(2) * (t69 + t84) - g(3) * t105, 0, 0, 0, 0, 0, 0, -g(1) * (t78 * t95 + t89) - g(2) * (-t78 * t88 + t96) - t80 * t103, -g(1) * (-t78 * t96 + t88) - g(2) * (t78 * t89 + t95) + t77 * t103, t57, -g(1) * (t83 * t79 + t85) - g(2) * t87 - g(3) * (t81 * pkin(3) + t78 * pkin(8) + t105) - (-qJ(2) - t83) * t104, 0, 0, 0, 0, 0, 0, -g(1) * (t78 * t97 + t91) - g(2) * (-t78 * t90 + t98) - t66 * t103, -g(1) * (-t78 * t98 + t90) - g(2) * (t78 * t91 + t97) + t65 * t103, t57, -g(1) * (t64 * t101 + t76 * t94 + t85) - g(2) * (pkin(4) * t96 + t87) - g(3) * (t81 * t64 - t78 * t76 + t105) + (-g(1) * t102 - g(2) * (-t64 * t78 - t76 * t81 - qJ(2))) * t82, 0, 0, 0, 0, 0, 0, -g(1) * (t78 * t99 + t93) - g(2) * (-t78 * t92 + t100) - t63 * t103, -g(1) * (-t78 * t100 + t92) - g(2) * (t78 * t93 + t99) + t62 * t103, t57, -g(1) * (t58 * t101 + t74 * t94 + t85) - g(2) * (t79 * t59 + t87) - g(3) * (t81 * t58 - t78 * t74 + t105) + (-g(1) * t59 - g(2) * (-t58 * t78 - t74 * t81 - qJ(2))) * t82;];
U_reg  = t1;
