% Calculate inertial parameters regressor of potential energy for
% S5RPRRR10
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d3,d4,d5,theta2]';
% 
% Output:
% U_reg [1x(5*10)]
%   inertial parameter regressor of Potential energy

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 19:11
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U_reg = S5RPRRR10_energypot_fixb_reg2_slag_vp(qJ, g, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(3,1),zeros(9,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPRRR10_energypot_fixb_reg2_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RPRRR10_energypot_fixb_reg2_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5RPRRR10_energypot_fixb_reg2_slag_vp: pkin has to be [9x1] (double)');

%% Symbolic Calculation
% From energy_potential_fixb_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 19:10:43
% EndTime: 2019-12-31 19:10:43
% DurationCPUTime: 0.18s
% Computational Cost: add. (131->61), mult. (143->79), div. (0->0), fcn. (138->10), ass. (0->36)
t74 = cos(qJ(4));
t60 = t74 * pkin(4) + pkin(3);
t67 = pkin(9) + qJ(3);
t61 = sin(t67);
t62 = cos(t67);
t76 = -pkin(8) - pkin(7);
t94 = t60 * t62 - t61 * t76;
t93 = g(3) * pkin(5);
t92 = g(3) * t61;
t69 = sin(pkin(9));
t91 = t69 * pkin(2) + pkin(5);
t68 = qJ(4) + qJ(5);
t63 = sin(t68);
t73 = sin(qJ(1));
t88 = t73 * t63;
t64 = cos(t68);
t87 = t73 * t64;
t72 = sin(qJ(4));
t86 = t73 * t72;
t85 = t73 * t74;
t75 = cos(qJ(1));
t84 = t75 * t63;
t83 = t75 * t64;
t82 = t75 * t72;
t81 = t75 * t74;
t70 = cos(pkin(9));
t58 = t70 * pkin(2) + pkin(1);
t71 = -pkin(6) - qJ(2);
t80 = t73 * t58 + t75 * t71;
t55 = t75 * t58;
t79 = -t73 * t71 + t55;
t78 = pkin(3) * t62 + pkin(7) * t61;
t77 = g(1) * t75 + g(2) * t73;
t56 = g(1) * t73 - g(2) * t75;
t53 = -g(3) * t62 + t77 * t61;
t1 = [0, 0, 0, 0, 0, 0, -t77, t56, -g(3), -t93, 0, 0, 0, 0, 0, 0, -g(3) * t69 - t77 * t70, -g(3) * t70 + t77 * t69, -t56, -g(1) * (t75 * pkin(1) + t73 * qJ(2)) - g(2) * (t73 * pkin(1) - t75 * qJ(2)) - t93, 0, 0, 0, 0, 0, 0, -t77 * t62 - t92, t53, -t56, -g(1) * t79 - g(2) * t80 - g(3) * t91, 0, 0, 0, 0, 0, 0, -g(1) * (t62 * t81 + t86) - g(2) * (t62 * t85 - t82) - t74 * t92, -g(1) * (-t62 * t82 + t85) - g(2) * (-t62 * t86 - t81) + t72 * t92, -t53, -g(1) * (t78 * t75 + t79) - g(2) * (t78 * t73 + t80) - g(3) * (t61 * pkin(3) - t62 * pkin(7) + t91), 0, 0, 0, 0, 0, 0, -g(1) * (t62 * t83 + t88) - g(2) * (t62 * t87 - t84) - t64 * t92, -g(1) * (-t62 * t84 + t87) - g(2) * (-t62 * t88 - t83) + t63 * t92, -t53, -g(1) * (t75 * t94 + t55) - g(2) * (-pkin(4) * t82 + t80) - g(3) * (t61 * t60 + t62 * t76 + t91) + (-g(1) * (pkin(4) * t72 - t71) - g(2) * t94) * t73;];
U_reg = t1;
