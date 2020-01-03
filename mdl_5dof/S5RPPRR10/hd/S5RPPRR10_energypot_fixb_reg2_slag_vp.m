% Calculate inertial parameters regressor of potential energy for
% S5RPPRR10
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d4,d5,theta2]';
% 
% Output:
% U_reg [1x(5*10)]
%   inertial parameter regressor of Potential energy

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 18:04
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U_reg = S5RPPRR10_energypot_fixb_reg2_slag_vp(qJ, g, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(3,1),zeros(8,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPPRR10_energypot_fixb_reg2_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RPPRR10_energypot_fixb_reg2_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RPPRR10_energypot_fixb_reg2_slag_vp: pkin has to be [8x1] (double)');

%% Symbolic Calculation
% From energy_potential_fixb_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 18:04:16
% EndTime: 2019-12-31 18:04:17
% DurationCPUTime: 0.15s
% Computational Cost: add. (100->54), mult. (168->66), div. (0->0), fcn. (167->8), ass. (0->34)
t69 = sin(qJ(1));
t66 = sin(pkin(8));
t83 = qJ(3) * t66;
t67 = cos(pkin(8));
t86 = t67 * t69;
t92 = pkin(2) * t86 + t69 * t83;
t71 = cos(qJ(1));
t77 = g(1) * t71 + g(2) * t69;
t91 = g(3) * pkin(5);
t88 = t66 * pkin(2) + pkin(5);
t68 = sin(qJ(4));
t87 = t66 * t68;
t85 = t67 * t71;
t84 = t71 * pkin(1) + t69 * qJ(2);
t82 = pkin(4) * t87;
t62 = t69 * pkin(1);
t81 = t62 + t92;
t80 = -t71 * qJ(2) + t62;
t79 = pkin(2) * t85 + t71 * t83 + t84;
t78 = -t67 * qJ(3) + t88;
t65 = qJ(4) + qJ(5);
t58 = sin(t65);
t59 = cos(t65);
t76 = t67 * t58 - t66 * t59;
t75 = t66 * t58 + t67 * t59;
t70 = cos(qJ(4));
t74 = t66 * t70 - t67 * t68;
t73 = t67 * t70 + t87;
t72 = -pkin(7) - pkin(6);
t57 = t70 * pkin(4) + pkin(3);
t52 = g(1) * t69 - g(2) * t71;
t51 = -g(3) * t66 - t67 * t77;
t50 = -g(3) * t67 + t66 * t77;
t1 = [0, 0, 0, 0, 0, 0, -t77, t52, -g(3), -t91, 0, 0, 0, 0, 0, 0, t51, t50, -t52, -g(1) * t84 - g(2) * t80 - t91, 0, 0, 0, 0, 0, 0, t51, -t52, -t50, -g(1) * t79 - g(2) * (t80 + t92) - g(3) * t78, 0, 0, 0, 0, 0, 0, -g(3) * t74 - t77 * t73, g(3) * t73 - t77 * t74, t52, -g(1) * (pkin(3) * t85 - t69 * pkin(6) + t79) - g(2) * (pkin(3) * t86 + (pkin(6) - qJ(2)) * t71 + t81) - g(3) * (t66 * pkin(3) + t78), 0, 0, 0, 0, 0, 0, g(3) * t76 - t77 * t75, g(3) * t75 + t77 * t76, t52, -g(1) * (t69 * t72 + t79) - g(2) * (t57 * t86 + t69 * t82 + t81) - g(3) * (t66 * t57 + (-pkin(4) * t68 - qJ(3)) * t67 + t88) + (-g(1) * (t57 * t67 + t82) - g(2) * (-qJ(2) - t72)) * t71;];
U_reg = t1;
