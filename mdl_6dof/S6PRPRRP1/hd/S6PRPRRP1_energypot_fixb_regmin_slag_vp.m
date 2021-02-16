% Calculate minimal parameter regressor of potential energy for
% S6PRPRRP1
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
% U_reg [1x23]
%   minimal parameter regressor of Potential energy

% Quelle: HybrDyn-Toolbox
% Datum: 2021-01-16 01:25
% Revision: 24b2e7d74a0c1a3b64fa2f8f5ad758691ad61af3 (2021-01-15)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U_reg = S6PRPRRP1_energypot_fixb_regmin_slag_vp(qJ, g, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(11,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRPRRP1_energypot_fixb_regmin_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6PRPRRP1_energypot_fixb_regmin_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6PRPRRP1_energypot_fixb_regmin_slag_vp: pkin has to be [11x1] (double)');

%% Symbolic Calculation
% From energy_potential_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2021-01-16 01:24:16
% EndTime: 2021-01-16 01:24:17
% DurationCPUTime: 0.16s
% Computational Cost: add. (281->76), mult. (324->118), div. (0->0), fcn. (383->18), ass. (0->53)
t88 = qJ(2) + pkin(11);
t84 = pkin(6) + t88;
t85 = pkin(6) - t88;
t77 = cos(t84) + cos(t85);
t86 = sin(t88);
t90 = sin(pkin(10));
t93 = cos(pkin(10));
t69 = t90 * t86 - t93 * t77 / 0.2e1;
t97 = sin(qJ(5));
t116 = t69 * t97;
t70 = t93 * t86 + t90 * t77 / 0.2e1;
t115 = t70 * t97;
t75 = -sin(t85) / 0.2e1 - sin(t84) / 0.2e1;
t114 = t75 * t97;
t94 = cos(pkin(6));
t113 = t90 * t94;
t91 = sin(pkin(6));
t96 = pkin(7) + qJ(3);
t112 = t91 * t96;
t98 = sin(qJ(4));
t111 = t91 * t98;
t99 = sin(qJ(2));
t110 = t91 * t99;
t109 = t93 * t94;
t108 = t94 * t99;
t101 = cos(qJ(4));
t107 = t101 * t91;
t102 = cos(qJ(2));
t106 = t90 * t102;
t105 = t93 * t102;
t104 = t94 * t96 + qJ(1);
t87 = cos(t88);
t71 = t86 * t113 - t93 * t87;
t65 = t90 * t107 + t71 * t98;
t72 = t86 * t109 + t90 * t87;
t66 = t93 * t107 + t72 * t98;
t73 = -t94 * t101 + t86 * t111;
t103 = g(1) * t65 - g(2) * t66 - g(3) * t73;
t100 = cos(qJ(5));
t95 = -qJ(6) - pkin(9);
t92 = cos(pkin(11));
t89 = sin(pkin(11));
t83 = t102 * pkin(2) + pkin(1);
t82 = t100 * pkin(5) + pkin(4);
t79 = -t89 * pkin(3) + t92 * pkin(8);
t78 = t92 * pkin(3) + t89 * pkin(8) + pkin(2);
t76 = pkin(2) * t108 - t112;
t74 = t86 * t107 + t94 * t98;
t68 = -t71 * t101 + t90 * t111;
t67 = t72 * t101 - t93 * t111;
t64 = -g(1) * (t68 * t100 + t115) - g(2) * (t67 * t100 + t116) - g(3) * (t74 * t100 + t114);
t63 = -g(1) * (t70 * t100 - t68 * t97) - g(2) * (t69 * t100 - t67 * t97) - g(3) * (t75 * t100 - t74 * t97);
t1 = [-g(3) * qJ(1), 0, -g(1) * (-t90 * t108 + t105) - g(2) * (t93 * t108 + t106) - g(3) * t110, -g(1) * (-t94 * t106 - t93 * t99) - g(2) * (t94 * t105 - t90 * t99) - g(3) * t91 * t102, -g(1) * (-t90 * t76 + t93 * t83) - g(2) * (t93 * t76 + t90 * t83) - g(3) * (pkin(2) * t110 + t104), 0, 0, 0, 0, 0, -g(1) * t68 - g(2) * t67 - g(3) * t74, -t103, 0, 0, 0, 0, 0, t64, t63, t64, t63, t103, -g(1) * (t68 * t82 + t65 * t95 + pkin(5) * t115 + (t79 * t113 + t93 * t78) * t102 + (-t78 * t113 + t93 * t79) * t99 + t90 * t112 + t93 * pkin(1)) - g(2) * (t67 * t82 - t66 * t95 + pkin(5) * t116 + (-t79 * t109 + t90 * t78) * t102 + (t78 * t109 + t90 * t79) * t99 - t93 * t112 + t90 * pkin(1)) - g(3) * (pkin(5) * t114 - t73 * t95 + t74 * t82 + (-t102 * t79 + t78 * t99) * t91 + t104);];
U_reg = t1;
