% Calculate inertial parameters regressor of potential energy for
% S5PRRRP7
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,alpha2,d2,d3,d4,theta1]';
% 
% Output:
% U_reg [1x(5*10)]
%   inertial parameter regressor of Potential energy

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-05 16:57
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U_reg = S5PRRRP7_energypot_fixb_reg2_slag_vp(qJ, g, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(3,1),zeros(9,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PRRRP7_energypot_fixb_reg2_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5PRRRP7_energypot_fixb_reg2_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5PRRRP7_energypot_fixb_reg2_slag_vp: pkin has to be [9x1] (double)');

%% Symbolic Calculation
% From energy_potential_fixb_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 16:56:25
% EndTime: 2019-12-05 16:56:25
% DurationCPUTime: 0.21s
% Computational Cost: add. (187->71), mult. (424->104), div. (0->0), fcn. (511->10), ass. (0->45)
t75 = sin(pkin(5));
t104 = pkin(6) * t75;
t80 = sin(qJ(3));
t103 = t75 * t80;
t81 = sin(qJ(2));
t102 = t75 * t81;
t83 = cos(qJ(3));
t101 = t75 * t83;
t84 = cos(qJ(2));
t100 = t75 * t84;
t77 = cos(pkin(5));
t99 = t77 * t81;
t98 = t77 * t84;
t74 = sin(pkin(9));
t76 = cos(pkin(9));
t97 = t76 * pkin(1) + t74 * t104;
t96 = t77 * pkin(6) + qJ(1);
t62 = -t74 * t99 + t76 * t84;
t95 = t62 * pkin(2) + t97;
t94 = pkin(2) * t102 + t96;
t79 = sin(qJ(4));
t93 = pkin(4) * t79 + pkin(7);
t92 = t74 * pkin(1) - t76 * t104;
t91 = g(1) * t74 - g(2) * t76;
t60 = t74 * t84 + t76 * t99;
t90 = t60 * pkin(2) + t92;
t61 = t74 * t98 + t76 * t81;
t89 = t61 * pkin(7) + t95;
t88 = -pkin(7) * t100 + t94;
t53 = t76 * t101 + t60 * t80;
t55 = -t74 * t101 + t62 * t80;
t63 = t80 * t102 - t77 * t83;
t87 = g(1) * t55 + g(2) * t53 + g(3) * t63;
t59 = t74 * t81 - t76 * t98;
t86 = t59 * pkin(7) + t90;
t85 = -g(1) * t61 - g(2) * t59 + g(3) * t100;
t82 = cos(qJ(4));
t78 = -qJ(5) - pkin(8);
t70 = t82 * pkin(4) + pkin(3);
t64 = t81 * t101 + t77 * t80;
t56 = t74 * t103 + t62 * t83;
t54 = -t76 * t103 + t60 * t83;
t51 = -g(1) * (t56 * t82 + t61 * t79) - g(2) * (t54 * t82 + t59 * t79) - g(3) * (-t79 * t100 + t64 * t82);
t50 = -g(1) * (-t56 * t79 + t61 * t82) - g(2) * (-t54 * t79 + t59 * t82) - g(3) * (-t82 * t100 - t64 * t79);
t1 = [0, 0, 0, 0, 0, 0, -g(1) * t76 - g(2) * t74, t91, -g(3), -g(3) * qJ(1), 0, 0, 0, 0, 0, 0, -g(1) * t62 - g(2) * t60 - g(3) * t102, -t85, -g(3) * t77 - t91 * t75, -g(1) * t97 - g(2) * t92 - g(3) * t96, 0, 0, 0, 0, 0, 0, -g(1) * t56 - g(2) * t54 - g(3) * t64, t87, t85, -g(1) * t89 - g(2) * t86 - g(3) * t88, 0, 0, 0, 0, 0, 0, t51, t50, -t87, -g(1) * (t56 * pkin(3) + t55 * pkin(8) + t89) - g(2) * (t54 * pkin(3) + t53 * pkin(8) + t86) - g(3) * (t64 * pkin(3) + t63 * pkin(8) + t88), 0, 0, 0, 0, 0, 0, t51, t50, -t87, -g(1) * (-t55 * t78 + t56 * t70 + t93 * t61 + t95) - g(2) * (-t53 * t78 + t54 * t70 + t93 * t59 + t90) - g(3) * (-t93 * t100 - t63 * t78 + t64 * t70 + t94);];
U_reg = t1;
