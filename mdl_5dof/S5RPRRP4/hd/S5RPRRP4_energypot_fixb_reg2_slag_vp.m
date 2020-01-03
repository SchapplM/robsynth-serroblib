% Calculate inertial parameters regressor of potential energy for
% S5RPRRP4
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d3,d4,theta2]';
% 
% Output:
% U_reg [1x(5*10)]
%   inertial parameter regressor of Potential energy

% Quelle: HybrDyn-Toolbox
% Datum: 2020-01-03 11:51
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U_reg = S5RPRRP4_energypot_fixb_reg2_slag_vp(qJ, g, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(3,1),zeros(8,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPRRP4_energypot_fixb_reg2_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RPRRP4_energypot_fixb_reg2_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RPRRP4_energypot_fixb_reg2_slag_vp: pkin has to be [8x1] (double)');

%% Symbolic Calculation
% From energy_potential_fixb_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2020-01-03 11:50:23
% EndTime: 2020-01-03 11:50:23
% DurationCPUTime: 0.21s
% Computational Cost: add. (122->62), mult. (169->81), div. (0->0), fcn. (168->8), ass. (0->36)
t73 = sin(qJ(3));
t82 = t73 * pkin(3) + qJ(2);
t95 = g(1) * pkin(5);
t77 = -pkin(7) - pkin(6);
t71 = sin(pkin(8));
t94 = g(1) * t71;
t74 = sin(qJ(1));
t67 = t74 * pkin(1);
t93 = g(2) * t67;
t75 = cos(qJ(3));
t64 = t75 * pkin(3) + pkin(2);
t70 = qJ(3) + qJ(4);
t65 = sin(t70);
t91 = t74 * t65;
t66 = cos(t70);
t90 = t74 * t66;
t89 = t74 * t73;
t88 = t74 * t75;
t76 = cos(qJ(1));
t87 = t76 * t65;
t86 = t76 * t66;
t85 = t76 * t73;
t84 = t76 * t75;
t83 = pkin(4) * t65 + t82;
t72 = cos(pkin(8));
t81 = pkin(2) * t72 + pkin(6) * t71;
t80 = -g(2) * t74 + g(3) * t76;
t61 = pkin(4) * t66 + t64;
t69 = -qJ(5) + t77;
t79 = t61 * t72 - t69 * t71;
t78 = t64 * t72 - t71 * t77;
t63 = g(2) * t76 + g(3) * t74;
t60 = g(1) * t72 + t80 * t71;
t59 = -t66 * t94 - g(2) * (t72 * t90 - t87) - g(3) * (-t72 * t86 - t91);
t58 = t65 * t94 - g(2) * (-t72 * t91 - t86) - g(3) * (t72 * t87 - t90);
t1 = [0, 0, 0, 0, 0, 0, t80, -t63, -g(1), -t95, 0, 0, 0, 0, 0, 0, t80 * t72 - t94, -t60, t63, -t95 - g(2) * (-t76 * qJ(2) + t67) - g(3) * (-t76 * pkin(1) - t74 * qJ(2)), 0, 0, 0, 0, 0, 0, -t75 * t94 - g(2) * (t72 * t88 - t85) - g(3) * (-t72 * t84 - t89), t73 * t94 - g(2) * (-t72 * t89 - t84) - g(3) * (t72 * t85 - t88), t60, -g(1) * (t71 * pkin(2) - t72 * pkin(6) + pkin(5)) - t93 + (-g(2) * t81 + g(3) * qJ(2)) * t74 + (g(2) * qJ(2) - g(3) * (-pkin(1) - t81)) * t76, 0, 0, 0, 0, 0, 0, t59, t58, t60, -g(1) * (t71 * t64 + t72 * t77 + pkin(5)) - t93 + (-g(2) * t78 + g(3) * t82) * t74 + (g(2) * t82 - g(3) * (-pkin(1) - t78)) * t76, 0, 0, 0, 0, 0, 0, t59, t58, t60, -g(1) * (t71 * t61 + t72 * t69 + pkin(5)) - t93 + (-g(2) * t79 + g(3) * t83) * t74 + (g(2) * t83 - g(3) * (-pkin(1) - t79)) * t76;];
U_reg = t1;
