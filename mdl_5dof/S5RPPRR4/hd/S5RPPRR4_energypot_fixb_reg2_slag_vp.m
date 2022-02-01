% Calculate inertial parameters regressor of potential energy for
% S5RPPRR4
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d4,d5,theta2,theta3]';
% 
% Output:
% U_reg [1x(5*10)]
%   inertial parameter regressor of Potential energy

% Quelle: HybrDyn-Toolbox
% Datum: 2022-01-23 09:17
% Revision: fd3771346c4aea32fdeb66112c511235427c26a7 (2022-01-20)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U_reg = S5RPPRR4_energypot_fixb_reg2_slag_vp(qJ, g, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(3,1),zeros(9,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPPRR4_energypot_fixb_reg2_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RPPRR4_energypot_fixb_reg2_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5RPPRR4_energypot_fixb_reg2_slag_vp: pkin has to be [9x1] (double)');

%% Symbolic Calculation
% From energy_potential_fixb_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2022-01-23 09:16:50
% EndTime: 2022-01-23 09:16:51
% DurationCPUTime: 0.18s
% Computational Cost: add. (132->73), mult. (159->94), div. (0->0), fcn. (158->10), ass. (0->40)
t80 = cos(pkin(9));
t67 = t80 * pkin(3) + pkin(2);
t77 = pkin(9) + qJ(4);
t69 = cos(t77);
t60 = pkin(4) * t69 + t67;
t82 = qJ(3) + pkin(6);
t76 = -pkin(7) - t82;
t79 = sin(pkin(8));
t81 = cos(pkin(8));
t101 = t60 * t81 - t76 * t79;
t100 = g(3) * pkin(5);
t78 = sin(pkin(9));
t99 = pkin(3) * t78;
t98 = g(3) * t79;
t84 = cos(qJ(1));
t95 = t81 * t84;
t70 = qJ(5) + t77;
t65 = sin(t70);
t83 = sin(qJ(1));
t94 = t83 * t65;
t66 = cos(t70);
t93 = t83 * t66;
t68 = sin(t77);
t92 = t83 * t68;
t91 = t83 * t69;
t90 = t83 * t78;
t89 = t83 * t80;
t88 = t84 * t80;
t71 = t83 * qJ(2);
t87 = t84 * pkin(1) + t71;
t86 = qJ(2) * t84;
t85 = g(1) * t84 + g(2) * t83;
t74 = t83 * pkin(1);
t64 = qJ(2) + t99;
t63 = g(1) * t83 - g(2) * t84;
t62 = pkin(2) * t81 + qJ(3) * t79 + pkin(1);
t61 = pkin(4) * t68 + t99;
t59 = t67 * t81 + t79 * t82 + pkin(1);
t58 = -g(3) * t81 + t85 * t79;
t1 = [0, 0, 0, 0, 0, 0, -t85, t63, -g(3), -t100, 0, 0, 0, 0, 0, 0, -t85 * t81 - t98, t58, -t63, -g(1) * t87 - g(2) * (t74 - t86) - t100, 0, 0, 0, 0, 0, 0, -g(1) * (t81 * t88 + t90) - g(2) * (-t78 * t84 + t81 * t89) - t80 * t98, -g(1) * (-t78 * t95 + t89) - g(2) * (-t81 * t90 - t88) + t78 * t98, -t58, -g(1) * (t62 * t84 + t71) - g(2) * (t62 * t83 - t86) - g(3) * (pkin(2) * t79 - qJ(3) * t81 + pkin(5)), 0, 0, 0, 0, 0, 0, -g(1) * (t69 * t95 + t92) - g(2) * (-t68 * t84 + t81 * t91) - t69 * t98, -g(1) * (-t68 * t95 + t91) - g(2) * (-t69 * t84 - t81 * t92) + t68 * t98, -t58, -g(1) * (t59 * t84 + t64 * t83) - g(2) * (t59 * t83 - t64 * t84) - g(3) * (t67 * t79 - t81 * t82 + pkin(5)), 0, 0, 0, 0, 0, 0, -g(1) * (t66 * t95 + t94) - g(2) * (-t65 * t84 + t81 * t93) - t66 * t98, -g(1) * (-t65 * t95 + t93) - g(2) * (-t66 * t84 - t81 * t94) + t65 * t98, -t58, -g(1) * (t83 * t61 + t87) - g(2) * (t101 * t83 + t74) - g(3) * (t60 * t79 + t76 * t81 + pkin(5)) + (-g(1) * t101 - g(2) * (-qJ(2) - t61)) * t84;];
U_reg = t1;
