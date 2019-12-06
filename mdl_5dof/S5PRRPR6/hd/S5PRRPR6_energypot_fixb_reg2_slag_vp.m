% Calculate inertial parameters regressor of potential energy for
% S5PRRPR6
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
% Datum: 2019-12-05 16:33
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U_reg = S5PRRPR6_energypot_fixb_reg2_slag_vp(qJ, g, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(3,1),zeros(10,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PRRPR6_energypot_fixb_reg2_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5PRRPR6_energypot_fixb_reg2_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S5PRRPR6_energypot_fixb_reg2_slag_vp: pkin has to be [10x1] (double)');

%% Symbolic Calculation
% From energy_potential_fixb_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 16:32:45
% EndTime: 2019-12-05 16:32:46
% DurationCPUTime: 0.27s
% Computational Cost: add. (199->82), mult. (424->122), div. (0->0), fcn. (511->12), ass. (0->46)
t80 = sin(pkin(5));
t108 = pkin(6) * t80;
t85 = sin(qJ(3));
t107 = t80 * t85;
t86 = sin(qJ(2));
t106 = t80 * t86;
t87 = cos(qJ(3));
t105 = t80 * t87;
t88 = cos(qJ(2));
t104 = t80 * t88;
t83 = cos(pkin(5));
t103 = t83 * t86;
t102 = t83 * t88;
t79 = sin(pkin(9));
t82 = cos(pkin(9));
t101 = t82 * pkin(1) + t79 * t108;
t100 = t83 * pkin(6) + qJ(1);
t63 = -t79 * t103 + t82 * t88;
t99 = t63 * pkin(2) + t101;
t98 = pkin(2) * t106 + t100;
t78 = sin(pkin(10));
t97 = pkin(4) * t78 + pkin(7);
t96 = t79 * pkin(1) - t82 * t108;
t95 = g(1) * t79 - g(2) * t82;
t61 = t82 * t103 + t79 * t88;
t94 = t61 * pkin(2) + t96;
t62 = t79 * t102 + t82 * t86;
t93 = t62 * pkin(7) + t99;
t92 = -pkin(7) * t104 + t98;
t54 = t82 * t105 + t61 * t85;
t56 = -t79 * t105 + t63 * t85;
t64 = t85 * t106 - t83 * t87;
t91 = g(1) * t56 + g(2) * t54 + g(3) * t64;
t60 = -t82 * t102 + t79 * t86;
t90 = t60 * pkin(7) + t94;
t89 = -g(1) * t62 - g(2) * t60 + g(3) * t104;
t84 = -pkin(8) - qJ(4);
t81 = cos(pkin(10));
t77 = pkin(10) + qJ(5);
t73 = cos(t77);
t72 = sin(t77);
t71 = t81 * pkin(4) + pkin(3);
t65 = t86 * t105 + t83 * t85;
t57 = t79 * t107 + t63 * t87;
t55 = -t82 * t107 + t61 * t87;
t1 = [0, 0, 0, 0, 0, 0, -g(1) * t82 - g(2) * t79, t95, -g(3), -g(3) * qJ(1), 0, 0, 0, 0, 0, 0, -g(1) * t63 - g(2) * t61 - g(3) * t106, -t89, -g(3) * t83 - t95 * t80, -g(1) * t101 - g(2) * t96 - g(3) * t100, 0, 0, 0, 0, 0, 0, -g(1) * t57 - g(2) * t55 - g(3) * t65, t91, t89, -g(1) * t93 - g(2) * t90 - g(3) * t92, 0, 0, 0, 0, 0, 0, -g(1) * (t57 * t81 + t62 * t78) - g(2) * (t55 * t81 + t60 * t78) - g(3) * (-t78 * t104 + t65 * t81), -g(1) * (-t57 * t78 + t62 * t81) - g(2) * (-t55 * t78 + t60 * t81) - g(3) * (-t81 * t104 - t65 * t78), -t91, -g(1) * (t57 * pkin(3) + t56 * qJ(4) + t93) - g(2) * (t55 * pkin(3) + t54 * qJ(4) + t90) - g(3) * (t65 * pkin(3) + t64 * qJ(4) + t92), 0, 0, 0, 0, 0, 0, -g(1) * (t57 * t73 + t62 * t72) - g(2) * (t55 * t73 + t60 * t72) - g(3) * (-t72 * t104 + t65 * t73), -g(1) * (-t57 * t72 + t62 * t73) - g(2) * (-t55 * t72 + t60 * t73) - g(3) * (-t73 * t104 - t65 * t72), -t91, -g(1) * (-t56 * t84 + t57 * t71 + t97 * t62 + t99) - g(2) * (-t54 * t84 + t55 * t71 + t97 * t60 + t94) - g(3) * (-t97 * t104 - t64 * t84 + t65 * t71 + t98);];
U_reg = t1;
