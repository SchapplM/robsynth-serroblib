% Calculate inertial parameters regressor of potential energy for
% S5RPRPR5
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d3,d5,theta2,theta4]';
% 
% Output:
% U_reg [1x(5*10)]
%   inertial parameter regressor of Potential energy

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-05 17:58
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U_reg = S5RPRPR5_energypot_fixb_reg2_slag_vp(qJ, g, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(3,1),zeros(9,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPRPR5_energypot_fixb_reg2_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RPRPR5_energypot_fixb_reg2_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5RPRPR5_energypot_fixb_reg2_slag_vp: pkin has to be [9x1] (double)');

%% Symbolic Calculation
% From energy_potential_fixb_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 17:57:11
% EndTime: 2019-12-05 17:57:11
% DurationCPUTime: 0.17s
% Computational Cost: add. (132->73), mult. (169->98), div. (0->0), fcn. (168->10), ass. (0->42)
t102 = g(1) * pkin(5);
t75 = sin(pkin(8));
t101 = g(1) * t75;
t79 = sin(qJ(1));
t100 = g(2) * t79;
t78 = sin(qJ(3));
t99 = t78 * pkin(3);
t80 = cos(qJ(3));
t65 = t80 * pkin(3) + pkin(2);
t81 = cos(qJ(1));
t98 = t75 * t81;
t76 = cos(pkin(8));
t97 = t76 * t81;
t74 = qJ(3) + pkin(9);
t68 = qJ(5) + t74;
t63 = sin(t68);
t96 = t79 * t63;
t64 = cos(t68);
t95 = t79 * t64;
t66 = sin(t74);
t94 = t79 * t66;
t67 = cos(t74);
t93 = t79 * t67;
t92 = t79 * t78;
t91 = t79 * t80;
t90 = t81 * t63;
t89 = t81 * t64;
t88 = t81 * t66;
t87 = t81 * t67;
t86 = t81 * t78;
t85 = t81 * t80;
t77 = -qJ(4) - pkin(6);
t84 = t81 * pkin(1) + t79 * qJ(2);
t83 = pkin(2) * t76 + pkin(6) * t75;
t82 = -g(3) * t81 + t100;
t73 = -pkin(7) + t77;
t70 = t81 * qJ(2);
t62 = g(2) * t81 + g(3) * t79;
t61 = pkin(4) * t66 + t99;
t60 = pkin(4) * t67 + t65;
t59 = g(1) * t76 + t82 * t75;
t1 = [0, 0, 0, 0, 0, 0, t82, t62, -g(1), -t102, 0, 0, 0, 0, 0, 0, t82 * t76 - t101, -t59, -t62, -t102 - g(2) * (-t79 * pkin(1) + t70) - g(3) * t84, 0, 0, 0, 0, 0, 0, -t80 * t101 - g(2) * (-t76 * t91 + t86) - g(3) * (t76 * t85 + t92), t78 * t101 - g(2) * (t76 * t92 + t85) - g(3) * (-t76 * t86 + t91), t59, -g(1) * (t75 * pkin(2) - t76 * pkin(6) + pkin(5)) - g(2) * t70 - g(3) * (t83 * t81 + t84) - (-pkin(1) - t83) * t100, 0, 0, 0, 0, 0, 0, -t67 * t101 - g(2) * (-t76 * t93 + t88) - g(3) * (t76 * t87 + t94), t66 * t101 - g(2) * (t76 * t94 + t87) - g(3) * (-t76 * t88 + t93), t59, -g(1) * (t75 * t65 + t76 * t77 + pkin(5)) - g(2) * (pkin(3) * t86 + t70) - g(3) * (t65 * t97 - t77 * t98 + t84) + (-g(2) * (-t65 * t76 + t75 * t77 - pkin(1)) - g(3) * t99) * t79, 0, 0, 0, 0, 0, 0, -t64 * t101 - g(2) * (-t76 * t95 + t90) - g(3) * (t76 * t89 + t96), t63 * t101 - g(2) * (t76 * t96 + t89) - g(3) * (-t76 * t90 + t95), t59, -g(1) * (t75 * t60 + t76 * t73 + pkin(5)) - g(2) * (t81 * t61 + t70) - g(3) * (t60 * t97 - t73 * t98 + t84) + (-g(2) * (-t60 * t76 + t73 * t75 - pkin(1)) - g(3) * t61) * t79;];
U_reg = t1;
