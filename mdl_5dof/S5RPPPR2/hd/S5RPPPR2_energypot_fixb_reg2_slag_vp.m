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
% Datum: 2019-12-05 17:32
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
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
% StartTime: 2019-12-05 17:31:30
% EndTime: 2019-12-05 17:31:30
% DurationCPUTime: 0.22s
% Computational Cost: add. (138->68), mult. (302->96), div. (0->0), fcn. (348->10), ass. (0->45)
t82 = sin(pkin(7));
t85 = cos(pkin(7));
t112 = pkin(2) * t85 + qJ(3) * t82;
t111 = g(1) * pkin(5);
t87 = sin(qJ(1));
t109 = g(2) * t87;
t81 = sin(pkin(8));
t108 = t81 * t82;
t84 = cos(pkin(8));
t107 = t82 * t84;
t106 = t82 * t87;
t89 = cos(qJ(1));
t105 = t82 * t89;
t104 = t87 * t81;
t103 = t87 * t84;
t102 = t89 * t81;
t101 = t89 * t84;
t100 = t89 * pkin(1) + t87 * qJ(2);
t98 = t112 * t89 + t100;
t97 = t82 * pkin(2) - qJ(3) * t85 + pkin(5);
t65 = t104 * t85 + t101;
t66 = -t103 * t85 + t102;
t78 = t89 * qJ(2);
t96 = t66 * pkin(3) - t65 * qJ(4) + t78;
t95 = -g(3) * t89 + t109;
t94 = pkin(3) * t107 + qJ(4) * t108 + t97;
t80 = sin(pkin(9));
t83 = cos(pkin(9));
t56 = t106 * t83 + t66 * t80;
t68 = t101 * t85 + t104;
t58 = -t105 * t83 + t68 * t80;
t63 = t107 * t80 + t83 * t85;
t93 = g(1) * t63 + g(2) * t56 + g(3) * t58;
t67 = t102 * t85 - t103;
t92 = t68 * pkin(3) + qJ(4) * t67 + t98;
t91 = g(1) * t108 - g(2) * t65 + g(3) * t67;
t90 = (-pkin(1) - t112) * t109;
t88 = cos(qJ(5));
t86 = sin(qJ(5));
t70 = g(2) * t89 + g(3) * t87;
t64 = t107 * t83 - t80 * t85;
t60 = g(1) * t85 + t82 * t95;
t59 = t105 * t80 + t68 * t83;
t57 = -t106 * t80 + t66 * t83;
t1 = [0, 0, 0, 0, 0, 0, t95, t70, -g(1), -t111, 0, 0, 0, 0, 0, 0, -g(1) * t82 + t85 * t95, -t60, -t70, -t111 - g(2) * (-pkin(1) * t87 + t78) - g(3) * t100, 0, 0, 0, 0, 0, 0, -g(1) * t107 - g(2) * t66 - g(3) * t68, t91, t60, -g(1) * t97 - g(2) * t78 - g(3) * t98 - t90, 0, 0, 0, 0, 0, 0, -g(1) * t64 - g(2) * t57 - g(3) * t59, t93, -t91, -g(1) * t94 - g(2) * t96 - g(3) * t92 - t90, 0, 0, 0, 0, 0, 0, -g(1) * (t108 * t86 + t64 * t88) - g(2) * (t57 * t88 - t65 * t86) - g(3) * (t59 * t88 + t67 * t86), -g(1) * (t108 * t88 - t64 * t86) - g(2) * (-t57 * t86 - t65 * t88) - g(3) * (-t59 * t86 + t67 * t88), -t93, -g(1) * (pkin(4) * t64 + pkin(6) * t63 + t94) - g(2) * (t57 * pkin(4) + t56 * pkin(6) + t96) - g(3) * (pkin(4) * t59 + pkin(6) * t58 + t92) - t90;];
U_reg = t1;
