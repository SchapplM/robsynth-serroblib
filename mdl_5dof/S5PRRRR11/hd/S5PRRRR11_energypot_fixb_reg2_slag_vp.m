% Calculate inertial parameters regressor of potential energy for
% S5PRRRR11
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,alpha3,d2,d3,d4,d5,theta1]';
% 
% Output:
% U_reg [1x(5*10)]
%   inertial parameter regressor of Potential energy

% Quelle: HybrDyn-Toolbox
% Datum: 2024-09-27 21:46
% Revision: 4128c0dfad36ceb08eaa7a9fbc8f797c888f59be (2024-07-29)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U_reg = S5PRRRR11_energypot_fixb_reg2_slag_vp(qJ, g, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(3,1),zeros(10,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PRRRR11_energypot_fixb_reg2_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5PRRRR11_energypot_fixb_reg2_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S5PRRRR11_energypot_fixb_reg2_slag_vp: pkin has to be [10x1] (double)');

%% Symbolic Calculation
% From energy_potential_fixb_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2024-09-27 21:45:48
% EndTime: 2024-09-27 21:45:48
% DurationCPUTime: 0.05s
% Computational Cost: add. (223->80), mult. (167->110), div. (0->0), fcn. (154->22), ass. (0->49)
t97 = pkin(7) + pkin(8);
t91 = sin(pkin(5));
t105 = pkin(7) * t91;
t96 = cos(qJ(3));
t76 = t96 * pkin(3) + pkin(2);
t95 = sin(qJ(3));
t104 = t91 * t95;
t93 = cos(pkin(5));
t103 = t93 * t95;
t102 = t93 * t96;
t101 = pkin(6) + qJ(1);
t89 = qJ(3) + qJ(4);
t80 = pkin(5) - t89;
t79 = pkin(5) + t89;
t87 = pkin(10) + qJ(2);
t77 = sin(t87);
t78 = cos(t87);
t100 = g(1) * t77 - g(2) * t78;
t90 = sin(pkin(10));
t92 = cos(pkin(10));
t99 = -g(1) * t92 - g(2) * t90;
t98 = pkin(4) * sin(qJ(4)) * t96 + (cos(qJ(4)) * pkin(4) + pkin(3)) * t95;
t88 = pkin(9) + t97;
t85 = qJ(5) + t89;
t84 = t92 * pkin(1);
t83 = t90 * pkin(1);
t82 = cos(t89);
t81 = sin(t89);
t74 = cos(t85);
t73 = sin(t85);
t72 = -qJ(5) + t80;
t71 = qJ(5) + t79;
t70 = cos(t79);
t69 = sin(t80);
t68 = cos(t71);
t67 = sin(t72);
t66 = cos(t80) / 0.2e1;
t65 = sin(t79) / 0.2e1;
t64 = cos(t72) / 0.2e1;
t63 = sin(t71) / 0.2e1;
t62 = pkin(4) * t82 + t76;
t61 = pkin(3) * t103 - t91 * t97;
t60 = t66 + t70 / 0.2e1;
t59 = t65 - t69 / 0.2e1;
t58 = t64 + t68 / 0.2e1;
t57 = t63 - t67 / 0.2e1;
t56 = -g(3) * t93 - t100 * t91;
t55 = -t91 * t88 + t98 * t93;
t1 = [0, 0, 0, 0, 0, 0, t99, g(1) * t90 - g(2) * t92, -g(3), -g(3) * qJ(1), 0, 0, 0, 0, 0, 0, -g(1) * t78 - g(2) * t77, t100, -g(3), t99 * pkin(1) - g(3) * t101, 0, 0, 0, 0, 0, 0, -g(1) * (-t77 * t103 + t78 * t96) - g(2) * (t78 * t103 + t77 * t96) - g(3) * t104, -g(1) * (-t102 * t77 - t78 * t95) - g(2) * (t102 * t78 - t77 * t95) - g(3) * t91 * t96, t56, -g(1) * (t78 * pkin(2) + t105 * t77 + t84) - g(2) * (t77 * pkin(2) - t105 * t78 + t83) - g(3) * (t93 * pkin(7) + t101), 0, 0, 0, 0, 0, 0, -g(1) * (-t77 * t59 + t78 * t82) - g(2) * (t78 * t59 + t77 * t82) - g(3) * (t66 - t70 / 0.2e1), -g(1) * (-t77 * t60 - t78 * t81) - g(2) * (t78 * t60 - t77 * t81) - g(3) * (t65 + t69 / 0.2e1), t56, -g(1) * (-t77 * t61 + t78 * t76 + t84) - g(2) * (t78 * t61 + t77 * t76 + t83) - g(3) * (pkin(3) * t104 + t93 * t97 + t101), 0, 0, 0, 0, 0, 0, -g(1) * (-t77 * t57 + t78 * t74) - g(2) * (t78 * t57 + t77 * t74) - g(3) * (t64 - t68 / 0.2e1), -g(1) * (-t77 * t58 - t78 * t73) - g(2) * (t78 * t58 - t77 * t73) - g(3) * (t63 + t67 / 0.2e1), t56, -g(1) * (-t77 * t55 + t78 * t62 + t84) - g(2) * (t78 * t55 + t77 * t62 + t83) - g(3) * (t93 * t88 + t91 * t98 + t101);];
U_reg = t1;
