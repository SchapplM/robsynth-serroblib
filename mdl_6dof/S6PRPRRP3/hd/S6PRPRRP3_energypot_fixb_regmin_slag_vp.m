% Calculate minimal parameter regressor of potential energy for
% S6PRPRRP3
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [11x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d2,d4,d5,theta1,theta3]';
% 
% Output:
% U_reg [1x25]
%   minimal parameter regressor of Potential energy

% Quelle: HybrDyn-Toolbox
% Datum: 2021-01-16 01:40
% Revision: 24b2e7d74a0c1a3b64fa2f8f5ad758691ad61af3 (2021-01-15)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U_reg = S6PRPRRP3_energypot_fixb_regmin_slag_vp(qJ, g, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(11,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRPRRP3_energypot_fixb_regmin_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6PRPRRP3_energypot_fixb_regmin_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6PRPRRP3_energypot_fixb_regmin_slag_vp: pkin has to be [11x1] (double)');

%% Symbolic Calculation
% From energy_potential_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2021-01-16 01:38:30
% EndTime: 2021-01-16 01:38:30
% DurationCPUTime: 0.18s
% Computational Cost: add. (213->81), mult. (347->127), div. (0->0), fcn. (418->12), ass. (0->48)
t83 = sin(pkin(10));
t86 = cos(pkin(10));
t91 = sin(qJ(2));
t87 = cos(pkin(6));
t93 = cos(qJ(2));
t98 = t87 * t93;
t69 = t83 * t91 - t86 * t98;
t90 = sin(qJ(5));
t108 = t69 * t90;
t71 = t83 * t98 + t86 * t91;
t107 = t71 * t90;
t84 = sin(pkin(6));
t106 = t83 * t84;
t105 = t83 * t87;
t104 = t84 * t91;
t103 = t84 * t93;
t102 = t86 * t84;
t101 = t86 * t87;
t89 = qJ(3) + pkin(8);
t100 = t87 * t89;
t99 = t87 * t91;
t97 = t83 * qJ(3);
t96 = t86 * qJ(3);
t68 = t83 * t99 - t86 * t93;
t81 = pkin(11) + qJ(4);
t77 = sin(t81);
t78 = cos(t81);
t62 = t78 * t106 + t68 * t77;
t70 = t83 * t93 + t86 * t99;
t63 = t78 * t102 + t70 * t77;
t66 = t77 * t104 - t87 * t78;
t95 = g(1) * t62 - g(2) * t63 - g(3) * t66;
t94 = -g(1) * t71 - g(2) * t69 + g(3) * t103;
t92 = cos(qJ(5));
t88 = -qJ(6) - pkin(9);
t85 = cos(pkin(11));
t82 = sin(pkin(11));
t80 = t86 * pkin(1);
t79 = t83 * pkin(1);
t76 = t92 * pkin(5) + pkin(4);
t75 = t85 * pkin(3) + pkin(2);
t74 = t82 * pkin(3) + pkin(7);
t67 = t78 * t104 + t87 * t77;
t65 = t77 * t106 - t68 * t78;
t64 = -t77 * t102 + t70 * t78;
t61 = -g(1) * (t65 * t92 + t107) - g(2) * (t64 * t92 + t108) - g(3) * (-t90 * t103 + t67 * t92);
t60 = -g(1) * (-t65 * t90 + t71 * t92) - g(2) * (-t64 * t90 + t69 * t92) - g(3) * (-t92 * t103 - t67 * t90);
t1 = [-g(3) * qJ(1), 0, g(1) * t68 - g(2) * t70 - g(3) * t104, -t94, -g(1) * (t82 * t106 - t68 * t85) - g(2) * (-t82 * t102 + t70 * t85) - g(3) * (t85 * t104 + t87 * t82), t94, -g(1) * ((t86 * pkin(2) + t87 * t97) * t93 + (-pkin(2) * t105 + t96) * t91 + pkin(7) * t106 + t80) - g(2) * ((t83 * pkin(2) - t87 * t96) * t93 + (pkin(2) * t101 + t97) * t91 - pkin(7) * t102 + t79) - g(3) * (t87 * pkin(7) + qJ(1) + (pkin(2) * t91 - qJ(3) * t93) * t84), 0, 0, 0, 0, 0, -g(1) * t65 - g(2) * t64 - g(3) * t67, -t95, 0, 0, 0, 0, 0, t61, t60, t61, t60, t95, -g(1) * (t65 * t76 + t62 * t88 + pkin(5) * t107 + (t83 * t100 + t86 * t75) * t93 + (-t75 * t105 + t86 * t89) * t91 + t80) - g(2) * (t64 * t76 - t63 * t88 + pkin(5) * t108 + (-t86 * t100 + t83 * t75) * t93 + (t75 * t101 + t83 * t89) * t91 + t79) - g(3) * (-t66 * t88 + t67 * t76 + t74 * t87 + qJ(1)) + (-g(3) * (t75 * t91 + (-pkin(5) * t90 - t89) * t93) + (-g(1) * t83 + g(2) * t86) * t74) * t84;];
U_reg = t1;
