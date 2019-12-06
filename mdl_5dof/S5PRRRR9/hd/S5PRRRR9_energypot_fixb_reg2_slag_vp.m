% Calculate inertial parameters regressor of potential energy for
% S5PRRRR9
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
% Datum: 2019-12-05 17:22
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U_reg = S5PRRRR9_energypot_fixb_reg2_slag_vp(qJ, g, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(3,1),zeros(10,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PRRRR9_energypot_fixb_reg2_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5PRRRR9_energypot_fixb_reg2_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S5PRRRR9_energypot_fixb_reg2_slag_vp: pkin has to be [10x1] (double)');

%% Symbolic Calculation
% From energy_potential_fixb_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 17:21:09
% EndTime: 2019-12-05 17:21:10
% DurationCPUTime: 0.24s
% Computational Cost: add. (199->82), mult. (424->122), div. (0->0), fcn. (511->12), ass. (0->46)
t80 = sin(pkin(5));
t109 = pkin(6) * t80;
t84 = sin(qJ(3));
t108 = t80 * t84;
t85 = sin(qJ(2));
t107 = t80 * t85;
t87 = cos(qJ(3));
t106 = t80 * t87;
t88 = cos(qJ(2));
t105 = t80 * t88;
t82 = cos(pkin(5));
t104 = t82 * t85;
t103 = t82 * t88;
t79 = sin(pkin(10));
t81 = cos(pkin(10));
t102 = t81 * pkin(1) + t79 * t109;
t101 = t82 * pkin(6) + qJ(1);
t64 = -t104 * t79 + t81 * t88;
t100 = t64 * pkin(2) + t102;
t99 = pkin(2) * t107 + t101;
t83 = sin(qJ(4));
t98 = pkin(4) * t83 + pkin(7);
t97 = t79 * pkin(1) - t109 * t81;
t96 = g(1) * t79 - g(2) * t81;
t62 = t104 * t81 + t79 * t88;
t95 = t62 * pkin(2) + t97;
t63 = t103 * t79 + t81 * t85;
t94 = pkin(7) * t63 + t100;
t93 = -pkin(7) * t105 + t99;
t55 = t106 * t81 + t62 * t84;
t57 = -t106 * t79 + t64 * t84;
t65 = t107 * t84 - t82 * t87;
t92 = g(1) * t57 + g(2) * t55 + g(3) * t65;
t61 = -t103 * t81 + t79 * t85;
t91 = pkin(7) * t61 + t95;
t90 = -g(1) * t63 - g(2) * t61 + g(3) * t105;
t89 = -pkin(9) - pkin(8);
t86 = cos(qJ(4));
t78 = qJ(4) + qJ(5);
t74 = cos(t78);
t73 = sin(t78);
t72 = pkin(4) * t86 + pkin(3);
t66 = t106 * t85 + t82 * t84;
t58 = t108 * t79 + t64 * t87;
t56 = -t108 * t81 + t62 * t87;
t1 = [0, 0, 0, 0, 0, 0, -g(1) * t81 - g(2) * t79, t96, -g(3), -g(3) * qJ(1), 0, 0, 0, 0, 0, 0, -g(1) * t64 - g(2) * t62 - g(3) * t107, -t90, -g(3) * t82 - t80 * t96, -g(1) * t102 - g(2) * t97 - g(3) * t101, 0, 0, 0, 0, 0, 0, -g(1) * t58 - g(2) * t56 - g(3) * t66, t92, t90, -g(1) * t94 - g(2) * t91 - g(3) * t93, 0, 0, 0, 0, 0, 0, -g(1) * (t58 * t86 + t63 * t83) - g(2) * (t56 * t86 + t61 * t83) - g(3) * (-t105 * t83 + t66 * t86), -g(1) * (-t58 * t83 + t63 * t86) - g(2) * (-t56 * t83 + t61 * t86) - g(3) * (-t105 * t86 - t66 * t83), -t92, -g(1) * (pkin(3) * t58 + pkin(8) * t57 + t94) - g(2) * (pkin(3) * t56 + pkin(8) * t55 + t91) - g(3) * (pkin(3) * t66 + pkin(8) * t65 + t93), 0, 0, 0, 0, 0, 0, -g(1) * (t58 * t74 + t63 * t73) - g(2) * (t56 * t74 + t61 * t73) - g(3) * (-t105 * t73 + t66 * t74), -g(1) * (-t58 * t73 + t63 * t74) - g(2) * (-t56 * t73 + t61 * t74) - g(3) * (-t105 * t74 - t66 * t73), -t92, -g(1) * (-t57 * t89 + t58 * t72 + t63 * t98 + t100) - g(2) * (-t55 * t89 + t56 * t72 + t61 * t98 + t95) - g(3) * (-t105 * t98 - t65 * t89 + t66 * t72 + t99);];
U_reg = t1;
