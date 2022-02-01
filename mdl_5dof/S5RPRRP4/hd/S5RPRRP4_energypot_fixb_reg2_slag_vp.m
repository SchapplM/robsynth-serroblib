% Calculate inertial parameters regressor of potential energy for
% S5RPRRP4
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d3,d4,theta2]';
% 
% Output:
% U_reg [1x(5*10)]
%   inertial parameter regressor of Potential energy

% Quelle: HybrDyn-Toolbox
% Datum: 2022-01-23 09:33
% Revision: fd3771346c4aea32fdeb66112c511235427c26a7 (2022-01-20)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U_reg = S5RPRRP4_energypot_fixb_reg2_slag_vp(qJ, g, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(3,1),zeros(8,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPRRP4_energypot_fixb_reg2_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RPRRP4_energypot_fixb_reg2_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RPRRP4_energypot_fixb_reg2_slag_vp: pkin has to be [8x1] (double)');

%% Symbolic Calculation
% From energy_potential_fixb_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2022-01-23 09:32:46
% EndTime: 2022-01-23 09:32:46
% DurationCPUTime: 0.18s
% Computational Cost: add. (122->64), mult. (162->80), div. (0->0), fcn. (161->8), ass. (0->42)
t81 = cos(qJ(3));
t100 = t81 * pkin(3);
t76 = qJ(3) + qJ(4);
t68 = cos(t76);
t62 = pkin(4) * t68 + pkin(2) + t100;
t83 = pkin(7) + pkin(6);
t75 = -qJ(5) - t83;
t77 = sin(pkin(8));
t78 = cos(pkin(8));
t104 = t62 * t78 - t75 * t77;
t103 = g(3) * pkin(5);
t102 = g(3) * t77;
t79 = sin(qJ(3));
t101 = t79 * pkin(3);
t99 = pkin(2) * t78 + pkin(1);
t98 = t77 * pkin(2) + pkin(5);
t95 = t77 * t81;
t67 = sin(t76);
t80 = sin(qJ(1));
t94 = t80 * t67;
t93 = t80 * t68;
t92 = t80 * t79;
t91 = t80 * t81;
t82 = cos(qJ(1));
t90 = t82 * t67;
t89 = t82 * t68;
t88 = t82 * t79;
t87 = t82 * t81;
t70 = t80 * qJ(2);
t86 = t82 * pkin(1) + t70;
t85 = t82 * qJ(2);
t84 = g(1) * t82 + g(2) * t80;
t73 = t80 * pkin(1);
t66 = qJ(2) + t101;
t65 = g(1) * t80 - g(2) * t82;
t64 = pkin(6) * t77 + t99;
t63 = pkin(4) * t67 + t101;
t61 = t100 * t78 + t77 * t83 + t99;
t60 = -g(3) * t78 + t77 * t84;
t59 = -g(1) * (t78 * t89 + t94) - g(2) * (t78 * t93 - t90) - t68 * t102;
t58 = -g(1) * (-t78 * t90 + t93) - g(2) * (-t78 * t94 - t89) + t67 * t102;
t1 = [0, 0, 0, 0, 0, 0, -t84, t65, -g(3), -t103, 0, 0, 0, 0, 0, 0, -t78 * t84 - t102, t60, -t65, -g(1) * t86 - g(2) * (t73 - t85) - t103, 0, 0, 0, 0, 0, 0, -g(1) * (t78 * t87 + t92) - g(2) * (t78 * t91 - t88) - g(3) * t95, -g(1) * (-t78 * t88 + t91) - g(2) * (-t78 * t92 - t87) + t79 * t102, -t60, -g(1) * (t64 * t82 + t70) - g(2) * (t64 * t80 - t85) - g(3) * (-pkin(6) * t78 + t98), 0, 0, 0, 0, 0, 0, t59, t58, -t60, -g(1) * (t61 * t82 + t66 * t80) - g(2) * (t61 * t80 - t66 * t82) - g(3) * (pkin(3) * t95 - t78 * t83 + t98), 0, 0, 0, 0, 0, 0, t59, t58, -t60, -g(1) * (t80 * t63 + t86) - g(2) * (t104 * t80 + t73) - g(3) * (t62 * t77 + t75 * t78 + pkin(5)) + (-g(1) * t104 - g(2) * (-qJ(2) - t63)) * t82;];
U_reg = t1;
