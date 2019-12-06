% Calculate inertial parameters regressor of potential energy for
% S5PRRPR7
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,alpha2,d2,d3,d5,theta1,theta4]';
% 
% Output:
% U_reg [1x(5*10)]
%   inertial parameter regressor of Potential energy

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-05 16:38
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U_reg = S5PRRPR7_energypot_fixb_reg2_slag_vp(qJ, g, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(3,1),zeros(10,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PRRPR7_energypot_fixb_reg2_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5PRRPR7_energypot_fixb_reg2_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S5PRRPR7_energypot_fixb_reg2_slag_vp: pkin has to be [10x1] (double)');

%% Symbolic Calculation
% From energy_potential_fixb_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 16:37:32
% EndTime: 2019-12-05 16:37:33
% DurationCPUTime: 0.16s
% Computational Cost: add. (223->75), mult. (525->117), div. (0->0), fcn. (651->12), ass. (0->49)
t83 = sin(pkin(5));
t112 = pkin(6) * t83;
t88 = sin(qJ(3));
t111 = t83 * t88;
t89 = sin(qJ(2));
t110 = t83 * t89;
t91 = cos(qJ(3));
t109 = t83 * t91;
t92 = cos(qJ(2));
t108 = t83 * t92;
t86 = cos(pkin(5));
t107 = t86 * t89;
t106 = t86 * t92;
t82 = sin(pkin(9));
t85 = cos(pkin(9));
t105 = t85 * pkin(1) + t82 * t112;
t104 = t86 * pkin(6) + qJ(1);
t103 = t82 * pkin(1) - t85 * t112;
t102 = g(1) * t82 - g(2) * t85;
t69 = t82 * t106 + t85 * t89;
t70 = -t82 * t107 + t85 * t92;
t101 = t70 * pkin(2) + t69 * pkin(7) + t105;
t100 = pkin(2) * t110 - pkin(7) * t108 + t104;
t68 = t85 * t107 + t82 * t92;
t59 = -t85 * t111 + t68 * t91;
t67 = -t85 * t106 + t82 * t89;
t81 = sin(pkin(10));
t84 = cos(pkin(10));
t50 = t59 * t81 - t67 * t84;
t61 = t82 * t111 + t70 * t91;
t52 = t61 * t81 - t69 * t84;
t72 = t89 * t109 + t86 * t88;
t56 = t84 * t108 + t72 * t81;
t99 = g(1) * t52 + g(2) * t50 + g(3) * t56;
t58 = t85 * t109 + t68 * t88;
t60 = -t82 * t109 + t70 * t88;
t71 = t88 * t110 - t86 * t91;
t98 = g(1) * t60 + g(2) * t58 + g(3) * t71;
t97 = t68 * pkin(2) + t67 * pkin(7) + t103;
t96 = -g(1) * t69 - g(2) * t67 + g(3) * t108;
t95 = t61 * pkin(3) + t60 * qJ(4) + t101;
t94 = t72 * pkin(3) + t71 * qJ(4) + t100;
t93 = t59 * pkin(3) + t58 * qJ(4) + t97;
t90 = cos(qJ(5));
t87 = sin(qJ(5));
t57 = -t81 * t108 + t72 * t84;
t53 = t61 * t84 + t69 * t81;
t51 = t59 * t84 + t67 * t81;
t1 = [0, 0, 0, 0, 0, 0, -g(1) * t85 - g(2) * t82, t102, -g(3), -g(3) * qJ(1), 0, 0, 0, 0, 0, 0, -g(1) * t70 - g(2) * t68 - g(3) * t110, -t96, -g(3) * t86 - t102 * t83, -g(1) * t105 - g(2) * t103 - g(3) * t104, 0, 0, 0, 0, 0, 0, -g(1) * t61 - g(2) * t59 - g(3) * t72, t98, t96, -g(1) * t101 - g(2) * t97 - g(3) * t100, 0, 0, 0, 0, 0, 0, -g(1) * t53 - g(2) * t51 - g(3) * t57, t99, -t98, -g(1) * t95 - g(2) * t93 - g(3) * t94, 0, 0, 0, 0, 0, 0, -g(1) * (t53 * t90 + t60 * t87) - g(2) * (t51 * t90 + t58 * t87) - g(3) * (t57 * t90 + t71 * t87), -g(1) * (-t53 * t87 + t60 * t90) - g(2) * (-t51 * t87 + t58 * t90) - g(3) * (-t57 * t87 + t71 * t90), -t99, -g(1) * (t53 * pkin(4) + t52 * pkin(8) + t95) - g(2) * (t51 * pkin(4) + t50 * pkin(8) + t93) - g(3) * (t57 * pkin(4) + t56 * pkin(8) + t94);];
U_reg = t1;
