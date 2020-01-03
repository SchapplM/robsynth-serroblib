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
% Datum: 2020-01-03 11:44
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
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
% StartTime: 2020-01-03 11:43:00
% EndTime: 2020-01-03 11:43:01
% DurationCPUTime: 0.24s
% Computational Cost: add. (132->71), mult. (169->96), div. (0->0), fcn. (168->10), ass. (0->38)
t76 = sin(qJ(3));
t84 = pkin(3) * t76 + qJ(2);
t98 = g(1) * pkin(5);
t73 = sin(pkin(8));
t96 = g(1) * t73;
t77 = sin(qJ(1));
t69 = t77 * pkin(1);
t95 = g(2) * t69;
t78 = cos(qJ(3));
t65 = t78 * pkin(3) + pkin(2);
t74 = cos(pkin(8));
t79 = cos(qJ(1));
t94 = t74 * t79;
t93 = t76 * t79;
t72 = qJ(3) + pkin(9);
t68 = qJ(5) + t72;
t63 = sin(t68);
t92 = t77 * t63;
t64 = cos(t68);
t91 = t77 * t64;
t66 = sin(t72);
t90 = t77 * t66;
t67 = cos(t72);
t89 = t77 * t67;
t88 = t77 * t76;
t87 = t77 * t78;
t86 = t78 * t79;
t75 = -qJ(4) - pkin(6);
t85 = pkin(4) * t66 + t84;
t83 = pkin(2) * t74 + pkin(6) * t73;
t82 = -g(2) * t77 + g(3) * t79;
t60 = pkin(4) * t67 + t65;
t71 = -pkin(7) + t75;
t81 = t60 * t74 - t71 * t73;
t80 = t65 * t74 - t73 * t75;
t62 = g(2) * t79 + g(3) * t77;
t59 = g(1) * t74 + t82 * t73;
t1 = [0, 0, 0, 0, 0, 0, t82, -t62, -g(1), -t98, 0, 0, 0, 0, 0, 0, t82 * t74 - t96, -t59, t62, -t98 - g(2) * (-t79 * qJ(2) + t69) - g(3) * (-t79 * pkin(1) - t77 * qJ(2)), 0, 0, 0, 0, 0, 0, -t78 * t96 - g(2) * (t74 * t87 - t93) - g(3) * (-t74 * t86 - t88), t76 * t96 - g(2) * (-t74 * t88 - t86) - g(3) * (t74 * t93 - t87), t59, -g(1) * (pkin(2) * t73 - pkin(6) * t74 + pkin(5)) - t95 + (-g(2) * t83 + g(3) * qJ(2)) * t77 + (g(2) * qJ(2) - g(3) * (-pkin(1) - t83)) * t79, 0, 0, 0, 0, 0, 0, -t67 * t96 - g(2) * (-t66 * t79 + t74 * t89) - g(3) * (-t67 * t94 - t90), t66 * t96 - g(2) * (-t67 * t79 - t74 * t90) - g(3) * (t66 * t94 - t89), t59, -g(1) * (t65 * t73 + t74 * t75 + pkin(5)) - t95 + (-g(2) * t80 + g(3) * t84) * t77 + (g(2) * t84 - g(3) * (-pkin(1) - t80)) * t79, 0, 0, 0, 0, 0, 0, -t64 * t96 - g(2) * (-t63 * t79 + t74 * t91) - g(3) * (-t64 * t94 - t92), t63 * t96 - g(2) * (-t64 * t79 - t74 * t92) - g(3) * (t63 * t94 - t91), t59, -g(1) * (t60 * t73 + t71 * t74 + pkin(5)) - t95 + (-g(2) * t81 + g(3) * t85) * t77 + (g(2) * t85 - g(3) * (-pkin(1) - t81)) * t79;];
U_reg = t1;
