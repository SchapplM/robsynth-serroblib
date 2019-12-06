% Calculate inertial parameters regressor of potential energy for
% S5PRRRR8
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,alpha2,d2,d3,d4,d5,theta1]';
% 
% Output:
% U_reg [1x(5*10)]
%   inertial parameter regressor of Potential energy

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-05 17:17
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U_reg = S5PRRRR8_energypot_fixb_reg2_slag_vp(qJ, g, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(3,1),zeros(10,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PRRRR8_energypot_fixb_reg2_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5PRRRR8_energypot_fixb_reg2_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S5PRRRR8_energypot_fixb_reg2_slag_vp: pkin has to be [10x1] (double)');

%% Symbolic Calculation
% From energy_potential_fixb_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 17:16:39
% EndTime: 2019-12-05 17:16:39
% DurationCPUTime: 0.24s
% Computational Cost: add. (207->81), mult. (360->124), div. (0->0), fcn. (422->12), ass. (0->46)
t81 = sin(pkin(10));
t82 = sin(pkin(5));
t109 = t81 * t82;
t83 = cos(pkin(10));
t108 = t82 * t83;
t86 = sin(qJ(3));
t107 = t82 * t86;
t87 = sin(qJ(2));
t106 = t82 * t87;
t89 = cos(qJ(3));
t105 = t82 * t89;
t90 = cos(qJ(2));
t104 = t82 * t90;
t84 = cos(pkin(5));
t103 = t84 * t86;
t102 = t84 * t87;
t101 = t84 * t90;
t100 = t83 * pkin(1) + pkin(6) * t109;
t99 = t84 * pkin(6) + qJ(1);
t98 = t81 * t107;
t77 = t81 * pkin(1);
t97 = -pkin(6) * t108 + t77;
t64 = t81 * t101 + t83 * t87;
t65 = -t81 * t102 + t83 * t90;
t74 = t89 * pkin(3) + pkin(2);
t91 = -pkin(8) - pkin(7);
t96 = pkin(3) * t98 - t64 * t91 + t65 * t74 + t100;
t95 = g(1) * t81 - g(2) * t83;
t94 = pkin(3) * t103 + t91 * t104 + t74 * t106 + t99;
t63 = t83 * t102 + t81 * t90;
t80 = qJ(3) + qJ(4);
t75 = sin(t80);
t76 = cos(t80);
t52 = t76 * t108 + t63 * t75;
t54 = -t76 * t109 + t65 * t75;
t58 = t75 * t106 - t84 * t76;
t93 = g(1) * t54 + g(2) * t52 + g(3) * t58;
t62 = -t83 * t101 + t81 * t87;
t51 = -g(1) * t64 - g(2) * t62 + g(3) * t104;
t92 = t63 * t74 - t62 * t91 + t77 + (-pkin(3) * t86 - pkin(6)) * t108;
t88 = cos(qJ(5));
t85 = sin(qJ(5));
t59 = t76 * t106 + t84 * t75;
t55 = t75 * t109 + t65 * t76;
t53 = -t75 * t108 + t63 * t76;
t1 = [0, 0, 0, 0, 0, 0, -g(1) * t83 - g(2) * t81, t95, -g(3), -g(3) * qJ(1), 0, 0, 0, 0, 0, 0, -g(1) * t65 - g(2) * t63 - g(3) * t106, -t51, -g(3) * t84 - t95 * t82, -g(1) * t100 - g(2) * t97 - g(3) * t99, 0, 0, 0, 0, 0, 0, -g(1) * (t65 * t89 + t98) - g(2) * (-t83 * t107 + t63 * t89) - g(3) * (t87 * t105 + t103), -g(1) * (t81 * t105 - t65 * t86) - g(2) * (-t83 * t105 - t63 * t86) - g(3) * (-t86 * t106 + t84 * t89), t51, -g(1) * (t65 * pkin(2) + t64 * pkin(7) + t100) - g(2) * (t63 * pkin(2) + t62 * pkin(7) + t97) - g(3) * ((pkin(2) * t87 - pkin(7) * t90) * t82 + t99), 0, 0, 0, 0, 0, 0, -g(1) * t55 - g(2) * t53 - g(3) * t59, t93, t51, -g(1) * t96 - g(2) * t92 - g(3) * t94, 0, 0, 0, 0, 0, 0, -g(1) * (t55 * t88 + t64 * t85) - g(2) * (t53 * t88 + t62 * t85) - g(3) * (-t85 * t104 + t59 * t88), -g(1) * (-t55 * t85 + t64 * t88) - g(2) * (-t53 * t85 + t62 * t88) - g(3) * (-t88 * t104 - t59 * t85), -t93, -g(1) * (t55 * pkin(4) + t54 * pkin(9) + t96) - g(2) * (t53 * pkin(4) + t52 * pkin(9) + t92) - g(3) * (t59 * pkin(4) + t58 * pkin(9) + t94);];
U_reg = t1;
