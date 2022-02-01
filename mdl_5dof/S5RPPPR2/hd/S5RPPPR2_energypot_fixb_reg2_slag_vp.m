% Calculate inertial parameters regressor of potential energy for
% S5RPPPR2
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d5,theta2,theta3,theta4]';
% 
% Output:
% U_reg [1x(5*10)]
%   inertial parameter regressor of Potential energy

% Quelle: HybrDyn-Toolbox
% Datum: 2022-01-23 09:00
% Revision: fd3771346c4aea32fdeb66112c511235427c26a7 (2022-01-20)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U_reg = S5RPPPR2_energypot_fixb_reg2_slag_vp(qJ, g, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(3,1),zeros(9,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPPPR2_energypot_fixb_reg2_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RPPPR2_energypot_fixb_reg2_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5RPPPR2_energypot_fixb_reg2_slag_vp: pkin has to be [9x1] (double)');

%% Symbolic Calculation
% From energy_potential_fixb_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2022-01-23 08:59:31
% EndTime: 2022-01-23 08:59:32
% DurationCPUTime: 0.21s
% Computational Cost: add. (138->81), mult. (252->121), div. (0->0), fcn. (280->10), ass. (0->50)
t110 = g(3) * pkin(5);
t81 = sin(pkin(7));
t109 = g(3) * t81;
t108 = t81 * qJ(3) + pkin(1);
t80 = sin(pkin(8));
t107 = t80 * qJ(4) + pkin(2);
t85 = sin(qJ(5));
t106 = t80 * t85;
t87 = cos(qJ(5));
t105 = t80 * t87;
t83 = cos(pkin(8));
t104 = t81 * t83;
t86 = sin(qJ(1));
t103 = t81 * t86;
t88 = cos(qJ(1));
t102 = t81 * t88;
t82 = cos(pkin(9));
t84 = cos(pkin(7));
t101 = t84 * t82;
t100 = t86 * t80;
t99 = t86 * t83;
t98 = t88 * t80;
t97 = t88 * t83;
t96 = t88 * qJ(2);
t95 = -t84 * qJ(3) + pkin(5);
t94 = qJ(4) * t83 - qJ(2);
t93 = g(1) * t88 + g(2) * t86;
t79 = sin(pkin(9));
t92 = t79 * pkin(4) - t82 * pkin(6) + qJ(3);
t74 = t82 * pkin(4) + t79 * pkin(6) + pkin(3);
t91 = -t74 * t80 + t94;
t68 = t84 * t99 - t98;
t70 = t84 * t97 + t100;
t90 = g(1) * (-t82 * t102 + t70 * t79) + g(2) * (-t82 * t103 + t68 * t79) + g(3) * (t79 * t104 + t101);
t67 = t84 * t100 + t97;
t89 = g(1) * (t84 * t98 - t99) + g(2) * t67 + t80 * t109;
t78 = t86 * qJ(2);
t75 = g(1) * t86 - g(2) * t88;
t73 = t83 * pkin(3) + t107;
t72 = pkin(2) * t84 + t108;
t71 = -t80 * pkin(3) + t94;
t66 = t82 * t106 + t87 * t83;
t65 = t83 * t101 + t81 * t79;
t64 = t82 * t104 - t84 * t79;
t62 = -g(3) * t84 + t93 * t81;
t61 = t74 * t83 + t107;
t60 = t73 * t84 + t108;
t57 = t84 * t105 - t65 * t85;
t56 = t61 * t84 + t92 * t81 + pkin(1);
t1 = [0, 0, 0, 0, 0, 0, -t93, t75, -g(3), -t110, 0, 0, 0, 0, 0, 0, -t93 * t84 - t109, t62, -t75, -g(1) * (t88 * pkin(1) + t78) - g(2) * (t86 * pkin(1) - t96) - t110, 0, 0, 0, 0, 0, 0, -g(1) * t70 - g(2) * t68 - g(3) * t104, t89, -t62, -g(1) * (t72 * t88 + t78) - g(2) * (t72 * t86 - t96) - g(3) * (t81 * pkin(2) + t95), 0, 0, 0, 0, 0, 0, -g(1) * (t79 * t102 + t70 * t82) - g(2) * (t79 * t103 + t68 * t82) - g(3) * t64, t90, -t89, -g(1) * (t60 * t88 - t71 * t86) - g(2) * (t60 * t86 + t71 * t88) - g(3) * (t73 * t81 + t95), 0, 0, 0, 0, 0, 0, -g(1) * ((t84 * t106 + t65 * t87) * t88 + t86 * (t82 * t105 - t85 * t83)) - g(2) * ((t65 * t86 - t82 * t98) * t87 + t67 * t85) - g(3) * (t81 * t106 + t64 * t87), -g(1) * (t57 * t88 - t86 * t66) - g(2) * (t57 * t86 + t88 * t66) - g(3) * (t81 * t105 - t64 * t85), -t90, -g(1) * (t56 * t88 - t91 * t86) - g(2) * (t56 * t86 + t91 * t88) - g(3) * (t61 * t81 - t92 * t84 + pkin(5));];
U_reg = t1;
