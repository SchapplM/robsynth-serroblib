% Calculate inertial parameters regressor of potential energy for
% S5RPPPR6
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d5,theta2,theta4]';
% 
% Output:
% U_reg [1x(5*10)]
%   inertial parameter regressor of Potential energy

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 17:48
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U_reg = S5RPPPR6_energypot_fixb_reg2_slag_vp(qJ, g, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(3,1),zeros(8,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPPPR6_energypot_fixb_reg2_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RPPPR6_energypot_fixb_reg2_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RPPPR6_energypot_fixb_reg2_slag_vp: pkin has to be [8x1] (double)');

%% Symbolic Calculation
% From energy_potential_fixb_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 17:47:50
% EndTime: 2019-12-31 17:47:50
% DurationCPUTime: 0.12s
% Computational Cost: add. (96->62), mult. (192->80), div. (0->0), fcn. (199->8), ass. (0->40)
t74 = sin(qJ(1));
t70 = sin(pkin(7));
t85 = qJ(3) * t70;
t72 = cos(pkin(7));
t93 = t72 * t74;
t98 = pkin(2) * t93 + t74 * t85;
t97 = g(3) * pkin(5);
t96 = g(3) * t72;
t95 = t70 * pkin(2) + pkin(5);
t69 = sin(pkin(8));
t94 = t69 * t72;
t75 = cos(qJ(5));
t92 = t72 * t75;
t76 = cos(qJ(1));
t91 = t72 * t76;
t90 = t74 * t69;
t71 = cos(pkin(8));
t89 = t74 * t71;
t88 = t76 * t69;
t87 = t76 * t71;
t86 = t76 * pkin(1) + t74 * qJ(2);
t84 = qJ(4) * t72;
t66 = t74 * pkin(1);
t83 = -t76 * qJ(2) + t66;
t82 = pkin(2) * t91 + t76 * t85 + t86;
t81 = -t72 * qJ(3) + t95;
t80 = g(1) * t76 + g(2) * t74;
t79 = t74 * pkin(3) + t76 * t84 + t82;
t78 = t74 * t84 + t66 + (-pkin(3) - qJ(2)) * t76 + t98;
t50 = -t70 * t87 + t90;
t52 = t70 * t89 + t88;
t77 = g(1) * t50 - g(2) * t52 + t71 * t96;
t73 = sin(qJ(5));
t62 = t70 * qJ(4);
t55 = g(1) * t74 - g(2) * t76;
t53 = t70 * t90 - t87;
t51 = t70 * t88 + t89;
t49 = g(3) * t70 + t80 * t72;
t48 = t80 * t70 - t96;
t1 = [0, 0, 0, 0, 0, 0, -t80, t55, -g(3), -t97, 0, 0, 0, 0, 0, 0, -t49, t48, -t55, -g(1) * t86 - g(2) * t83 - t97, 0, 0, 0, 0, 0, 0, -t55, t49, -t48, -g(1) * t82 - g(2) * (t83 + t98) - g(3) * t81, 0, 0, 0, 0, 0, 0, -g(1) * t51 - g(2) * t53 + g(3) * t94, t77, -t49, -g(1) * t79 - g(2) * t78 - g(3) * (t62 + t81), 0, 0, 0, 0, 0, 0, -g(1) * (t51 * t75 + t73 * t91) - g(2) * (t53 * t75 + t73 * t93) - g(3) * (-t69 * t92 + t70 * t73), -g(1) * (-t51 * t73 + t75 * t91) - g(2) * (-t53 * t73 + t74 * t92) - g(3) * (t70 * t75 + t73 * t94), -t77, -g(1) * (t51 * pkin(4) + t50 * pkin(6) + t79) - g(2) * (t53 * pkin(4) - t52 * pkin(6) + t78) - g(3) * (t62 + t95) - (-pkin(4) * t69 + pkin(6) * t71 - qJ(3)) * t96;];
U_reg = t1;
