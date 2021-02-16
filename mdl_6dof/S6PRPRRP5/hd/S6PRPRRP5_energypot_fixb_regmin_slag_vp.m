% Calculate minimal parameter regressor of potential energy for
% S6PRPRRP5
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d2,d4,d5,theta1]';
% 
% Output:
% U_reg [1x25]
%   minimal parameter regressor of Potential energy

% Quelle: HybrDyn-Toolbox
% Datum: 2021-01-16 01:52
% Revision: 24b2e7d74a0c1a3b64fa2f8f5ad758691ad61af3 (2021-01-15)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U_reg = S6PRPRRP5_energypot_fixb_regmin_slag_vp(qJ, g, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(10,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRPRRP5_energypot_fixb_regmin_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6PRPRRP5_energypot_fixb_regmin_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6PRPRRP5_energypot_fixb_regmin_slag_vp: pkin has to be [10x1] (double)');

%% Symbolic Calculation
% From energy_potential_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2021-01-16 01:51:49
% EndTime: 2021-01-16 01:51:49
% DurationCPUTime: 0.18s
% Computational Cost: add. (151->86), mult. (335->125), div. (0->0), fcn. (403->10), ass. (0->56)
t79 = sin(pkin(6));
t85 = sin(qJ(2));
t105 = t79 * t85;
t81 = cos(pkin(6));
t87 = cos(qJ(4));
t84 = sin(qJ(4));
t88 = cos(qJ(2));
t97 = t84 * t88;
t70 = -t79 * t97 + t81 * t87;
t83 = sin(qJ(5));
t86 = cos(qJ(5));
t112 = g(3) * (t86 * t105 - t70 * t83);
t111 = g(3) * (t83 * t105 + t70 * t86);
t78 = sin(pkin(10));
t110 = t78 * pkin(2);
t80 = cos(pkin(10));
t109 = t80 * pkin(2);
t101 = t81 * t85;
t66 = t80 * t101 + t78 * t88;
t62 = t66 * t83;
t68 = -t78 * t101 + t80 * t88;
t108 = t68 * t83;
t107 = t68 * t86;
t106 = t78 * t79;
t104 = t79 * t87;
t103 = t79 * t88;
t102 = t80 * t79;
t100 = t81 * t88;
t90 = pkin(2) + pkin(8);
t99 = t81 * t90;
t98 = t84 * t85;
t96 = qJ(3) * t88;
t74 = t78 * qJ(3);
t75 = t80 * qJ(3);
t95 = t81 * t75;
t94 = t81 * t97 + t104;
t67 = t78 * t100 + t80 * t85;
t56 = t84 * t106 - t67 * t87;
t65 = -t80 * t100 + t78 * t85;
t58 = t84 * t102 + t65 * t87;
t69 = t87 * t103 + t81 * t84;
t93 = g(1) * t56 - g(2) * t58 + g(3) * t69;
t92 = -g(1) * t67 - g(2) * t65 + g(3) * t103;
t91 = g(1) * t68 + g(2) * t66 + g(3) * t105;
t89 = pkin(3) + pkin(7);
t82 = -qJ(6) - pkin(9);
t77 = t80 * pkin(1);
t76 = t78 * pkin(1);
t73 = t86 * pkin(5) + pkin(4);
t72 = t81 * t74;
t63 = t66 * t86;
t59 = -t87 * t102 + t65 * t84;
t57 = t78 * t104 + t67 * t84;
t55 = -t78 * t98 + t94 * t80;
t54 = t94 * t78 + t80 * t98;
t1 = [-g(3) * qJ(1), 0, -t91, -t92, t91, t92, -g(1) * ((t72 + t109) * t88 + (-t81 * t110 + t75) * t85 + pkin(7) * t106 + t77) - g(2) * ((-t95 + t110) * t88 + (t81 * t109 + t74) * t85 - pkin(7) * t102 + t76) - g(3) * (t81 * pkin(7) + qJ(1) + (pkin(2) * t85 - t96) * t79), 0, 0, 0, 0, 0, -g(1) * t57 - g(2) * t59 - g(3) * t70, t93, 0, 0, 0, 0, 0, -g(1) * (t54 * t86 + t108) - g(2) * (-t55 * t86 + t62) - t111, -g(1) * (-t54 * t83 + t107) - g(2) * (t55 * t83 + t63) - t112, -g(1) * (t57 * t86 + t108) - g(2) * (t59 * t86 + t62) - t111, -g(1) * (-t57 * t83 + t107) - g(2) * (-t59 * t83 + t63) - t112, -t93, -g(1) * (t57 * t73 - t56 * t82 + pkin(5) * t108 + (t80 * t90 + t72) * t88 + (-t78 * t99 + t75) * t85 + t77) - g(2) * (t59 * t73 + t58 * t82 + pkin(5) * t62 + (t78 * t90 - t95) * t88 + (t80 * t99 + t74) * t85 + t76) - g(3) * (-t69 * t82 + t70 * t73 + t89 * t81 + qJ(1)) + (-g(3) * (-t96 + (pkin(5) * t83 + t90) * t85) + (-g(1) * t78 + g(2) * t80) * t89) * t79;];
U_reg = t1;
