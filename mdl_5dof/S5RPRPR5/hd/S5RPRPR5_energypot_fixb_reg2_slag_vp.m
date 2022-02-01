% Calculate inertial parameters regressor of potential energy for
% S5RPRPR5
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d3,d5,theta2,theta4]';
% 
% Output:
% U_reg [1x(5*10)]
%   inertial parameter regressor of Potential energy

% Quelle: HybrDyn-Toolbox
% Datum: 2022-01-23 09:26
% Revision: fd3771346c4aea32fdeb66112c511235427c26a7 (2022-01-20)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U_reg = S5RPRPR5_energypot_fixb_reg2_slag_vp(qJ, g, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(3,1),zeros(9,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPRPR5_energypot_fixb_reg2_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RPRPR5_energypot_fixb_reg2_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5RPRPR5_energypot_fixb_reg2_slag_vp: pkin has to be [9x1] (double)');

%% Symbolic Calculation
% From energy_potential_fixb_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2022-01-23 09:25:52
% EndTime: 2022-01-23 09:25:52
% DurationCPUTime: 0.18s
% Computational Cost: add. (132->73), mult. (162->94), div. (0->0), fcn. (161->10), ass. (0->47)
t84 = cos(qJ(3));
t106 = t84 * pkin(3);
t78 = qJ(3) + pkin(9);
t69 = cos(t78);
t61 = pkin(4) * t69 + pkin(2) + t106;
t81 = qJ(4) + pkin(6);
t77 = -pkin(7) - t81;
t79 = sin(pkin(8));
t80 = cos(pkin(8));
t110 = t61 * t80 - t77 * t79;
t109 = g(3) * pkin(5);
t108 = g(3) * t79;
t82 = sin(qJ(3));
t107 = t82 * pkin(3);
t105 = pkin(2) * t80 + pkin(1);
t104 = t79 * pkin(2) + pkin(5);
t101 = t79 * t84;
t70 = qJ(5) + t78;
t65 = sin(t70);
t83 = sin(qJ(1));
t100 = t83 * t65;
t66 = cos(t70);
t99 = t83 * t66;
t68 = sin(t78);
t98 = t83 * t68;
t97 = t83 * t69;
t96 = t83 * t82;
t95 = t83 * t84;
t85 = cos(qJ(1));
t94 = t85 * t65;
t93 = t85 * t66;
t92 = t85 * t68;
t91 = t85 * t69;
t90 = t85 * t82;
t89 = t85 * t84;
t72 = t83 * qJ(2);
t88 = t85 * pkin(1) + t72;
t87 = t85 * qJ(2);
t86 = g(1) * t85 + g(2) * t83;
t75 = t83 * pkin(1);
t67 = qJ(2) + t107;
t64 = g(1) * t83 - g(2) * t85;
t63 = t79 * pkin(6) + t105;
t62 = pkin(4) * t68 + t107;
t60 = -g(3) * t80 + t86 * t79;
t59 = t80 * t106 + t81 * t79 + t105;
t1 = [0, 0, 0, 0, 0, 0, -t86, t64, -g(3), -t109, 0, 0, 0, 0, 0, 0, -t86 * t80 - t108, t60, -t64, -g(1) * t88 - g(2) * (t75 - t87) - t109, 0, 0, 0, 0, 0, 0, -g(1) * (t80 * t89 + t96) - g(2) * (t80 * t95 - t90) - g(3) * t101, -g(1) * (-t80 * t90 + t95) - g(2) * (-t80 * t96 - t89) + t82 * t108, -t60, -g(1) * (t63 * t85 + t72) - g(2) * (t63 * t83 - t87) - g(3) * (-t80 * pkin(6) + t104), 0, 0, 0, 0, 0, 0, -g(1) * (t80 * t91 + t98) - g(2) * (t80 * t97 - t92) - t69 * t108, -g(1) * (-t80 * t92 + t97) - g(2) * (-t80 * t98 - t91) + t68 * t108, -t60, -g(1) * (t59 * t85 + t67 * t83) - g(2) * (t59 * t83 - t67 * t85) - g(3) * (pkin(3) * t101 - t80 * t81 + t104), 0, 0, 0, 0, 0, 0, -g(1) * (t80 * t93 + t100) - g(2) * (t80 * t99 - t94) - t66 * t108, -g(1) * (-t80 * t94 + t99) - g(2) * (-t80 * t100 - t93) + t65 * t108, -t60, -g(1) * (t83 * t62 + t88) - g(2) * (t110 * t83 + t75) - g(3) * (t79 * t61 + t80 * t77 + pkin(5)) + (-g(1) * t110 - g(2) * (-qJ(2) - t62)) * t85;];
U_reg = t1;
