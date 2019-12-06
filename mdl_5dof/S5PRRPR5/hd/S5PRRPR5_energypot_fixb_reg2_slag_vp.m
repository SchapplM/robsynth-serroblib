% Calculate inertial parameters regressor of potential energy for
% S5PRRPR5
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
% Datum: 2019-12-05 16:28
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U_reg = S5PRRPR5_energypot_fixb_reg2_slag_vp(qJ, g, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(3,1),zeros(10,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PRRPR5_energypot_fixb_reg2_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5PRRPR5_energypot_fixb_reg2_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S5PRRPR5_energypot_fixb_reg2_slag_vp: pkin has to be [10x1] (double)');

%% Symbolic Calculation
% From energy_potential_fixb_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 16:27:53
% EndTime: 2019-12-05 16:27:53
% DurationCPUTime: 0.24s
% Computational Cost: add. (207->81), mult. (360->124), div. (0->0), fcn. (422->12), ass. (0->46)
t82 = sin(pkin(9));
t83 = sin(pkin(5));
t110 = t82 * t83;
t84 = cos(pkin(9));
t109 = t83 * t84;
t88 = sin(qJ(3));
t108 = t83 * t88;
t89 = sin(qJ(2));
t107 = t83 * t89;
t91 = cos(qJ(3));
t106 = t83 * t91;
t92 = cos(qJ(2));
t105 = t83 * t92;
t85 = cos(pkin(5));
t104 = t85 * t88;
t103 = t85 * t89;
t102 = t85 * t92;
t101 = t84 * pkin(1) + pkin(6) * t110;
t100 = t85 * pkin(6) + qJ(1);
t99 = t82 * t108;
t78 = t82 * pkin(1);
t98 = -pkin(6) * t109 + t78;
t65 = t82 * t102 + t84 * t89;
t66 = -t82 * t103 + t84 * t92;
t75 = t91 * pkin(3) + pkin(2);
t86 = -qJ(4) - pkin(7);
t97 = pkin(3) * t99 - t65 * t86 + t66 * t75 + t101;
t96 = g(1) * t82 - g(2) * t84;
t95 = pkin(3) * t104 + t86 * t105 + t75 * t107 + t100;
t64 = t84 * t103 + t82 * t92;
t81 = qJ(3) + pkin(10);
t76 = sin(t81);
t77 = cos(t81);
t53 = t77 * t109 + t64 * t76;
t55 = -t77 * t110 + t66 * t76;
t59 = t76 * t107 - t85 * t77;
t94 = g(1) * t55 + g(2) * t53 + g(3) * t59;
t63 = -t84 * t102 + t82 * t89;
t52 = -g(1) * t65 - g(2) * t63 + g(3) * t105;
t93 = t64 * t75 - t63 * t86 + t78 + (-pkin(3) * t88 - pkin(6)) * t109;
t90 = cos(qJ(5));
t87 = sin(qJ(5));
t60 = t77 * t107 + t85 * t76;
t56 = t76 * t110 + t66 * t77;
t54 = -t76 * t109 + t64 * t77;
t1 = [0, 0, 0, 0, 0, 0, -g(1) * t84 - g(2) * t82, t96, -g(3), -g(3) * qJ(1), 0, 0, 0, 0, 0, 0, -g(1) * t66 - g(2) * t64 - g(3) * t107, -t52, -g(3) * t85 - t96 * t83, -g(1) * t101 - g(2) * t98 - g(3) * t100, 0, 0, 0, 0, 0, 0, -g(1) * (t66 * t91 + t99) - g(2) * (-t84 * t108 + t64 * t91) - g(3) * (t89 * t106 + t104), -g(1) * (t82 * t106 - t66 * t88) - g(2) * (-t84 * t106 - t64 * t88) - g(3) * (-t88 * t107 + t85 * t91), t52, -g(1) * (t66 * pkin(2) + t65 * pkin(7) + t101) - g(2) * (t64 * pkin(2) + t63 * pkin(7) + t98) - g(3) * ((pkin(2) * t89 - pkin(7) * t92) * t83 + t100), 0, 0, 0, 0, 0, 0, -g(1) * t56 - g(2) * t54 - g(3) * t60, t94, t52, -g(1) * t97 - g(2) * t93 - g(3) * t95, 0, 0, 0, 0, 0, 0, -g(1) * (t56 * t90 + t65 * t87) - g(2) * (t54 * t90 + t63 * t87) - g(3) * (-t87 * t105 + t60 * t90), -g(1) * (-t56 * t87 + t65 * t90) - g(2) * (-t54 * t87 + t63 * t90) - g(3) * (-t90 * t105 - t60 * t87), -t94, -g(1) * (t56 * pkin(4) + t55 * pkin(8) + t97) - g(2) * (t54 * pkin(4) + t53 * pkin(8) + t93) - g(3) * (t60 * pkin(4) + t59 * pkin(8) + t95);];
U_reg = t1;
