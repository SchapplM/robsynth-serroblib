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
% Datum: 2020-01-03 11:23
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
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
% StartTime: 2020-01-03 11:22:53
% EndTime: 2020-01-03 11:22:54
% DurationCPUTime: 0.24s
% Computational Cost: add. (138->69), mult. (302->98), div. (0->0), fcn. (348->10), ass. (0->44)
t80 = sin(pkin(7));
t83 = cos(pkin(7));
t109 = pkin(2) * t83 + qJ(3) * t80;
t108 = g(1) * pkin(5);
t79 = sin(pkin(8));
t106 = t79 * t80;
t82 = cos(pkin(8));
t105 = t80 * t82;
t85 = sin(qJ(1));
t104 = t80 * t85;
t87 = cos(qJ(1));
t103 = t80 * t87;
t102 = t85 * t79;
t101 = t85 * t82;
t100 = t87 * t79;
t99 = t87 * t82;
t97 = t85 * qJ(2);
t77 = t85 * pkin(1);
t96 = t109 * t85 + t77;
t95 = t80 * pkin(2) - t83 * qJ(3) + pkin(5);
t94 = -g(2) * t85 + g(3) * t87;
t67 = t83 * t100 - t101;
t68 = -t83 * t99 - t102;
t93 = t68 * pkin(3) - t67 * qJ(4) - t97;
t92 = pkin(3) * t105 + qJ(4) * t106 + t95;
t65 = t83 * t102 + t99;
t66 = t83 * t101 - t100;
t91 = t66 * pkin(3) + t65 * qJ(4) + t96;
t78 = sin(pkin(9));
t81 = cos(pkin(9));
t56 = -t81 * t104 + t66 * t78;
t58 = t81 * t103 + t68 * t78;
t63 = t78 * t105 + t83 * t81;
t90 = g(1) * t63 + g(2) * t56 + g(3) * t58;
t89 = g(1) * t106 + g(2) * t65 - g(3) * t67;
t88 = (g(2) * qJ(2) - g(3) * (-pkin(1) - t109)) * t87;
t86 = cos(qJ(5));
t84 = sin(qJ(5));
t70 = g(2) * t87 + g(3) * t85;
t64 = t81 * t105 - t83 * t78;
t60 = g(1) * t83 + t94 * t80;
t59 = -t78 * t103 + t68 * t81;
t57 = t78 * t104 + t66 * t81;
t1 = [0, 0, 0, 0, 0, 0, t94, -t70, -g(1), -t108, 0, 0, 0, 0, 0, 0, -g(1) * t80 + t94 * t83, -t60, t70, -t108 - g(2) * (-t87 * qJ(2) + t77) - g(3) * (-t87 * pkin(1) - t97), 0, 0, 0, 0, 0, 0, -g(1) * t105 - g(2) * t66 - g(3) * t68, t89, t60, -g(1) * t95 - g(2) * t96 + g(3) * t97 + t88, 0, 0, 0, 0, 0, 0, -g(1) * t64 - g(2) * t57 - g(3) * t59, t90, -t89, -g(1) * t92 - g(2) * t91 - g(3) * t93 + t88, 0, 0, 0, 0, 0, 0, -g(1) * (t84 * t106 + t64 * t86) - g(2) * (t57 * t86 + t65 * t84) - g(3) * (t59 * t86 - t67 * t84), -g(1) * (t86 * t106 - t64 * t84) - g(2) * (-t57 * t84 + t65 * t86) - g(3) * (-t59 * t84 - t67 * t86), -t90, -g(1) * (t64 * pkin(4) + t63 * pkin(6) + t92) - g(2) * (t57 * pkin(4) + t56 * pkin(6) + t91) - g(3) * (t59 * pkin(4) + t58 * pkin(6) + t93) + t88;];
U_reg = t1;
