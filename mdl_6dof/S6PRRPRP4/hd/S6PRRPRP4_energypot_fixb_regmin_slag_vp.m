% Calculate minimal parameter regressor of potential energy for
% S6PRRPRP4
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d2,d3,d5,theta1]';
% 
% Output:
% U_reg [1x26]
%   minimal parameter regressor of Potential energy

% Quelle: HybrDyn-Toolbox
% Datum: 2021-01-16 03:14
% Revision: 24b2e7d74a0c1a3b64fa2f8f5ad758691ad61af3 (2021-01-15)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U_reg = S6PRRPRP4_energypot_fixb_regmin_slag_vp(qJ, g, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(10,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRRPRP4_energypot_fixb_regmin_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6PRRPRP4_energypot_fixb_regmin_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6PRRPRP4_energypot_fixb_regmin_slag_vp: pkin has to be [10x1] (double)');

%% Symbolic Calculation
% From energy_potential_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2021-01-16 03:12:46
% EndTime: 2021-01-16 03:12:46
% DurationCPUTime: 0.33s
% Computational Cost: add. (192->82), mult. (400->128), div. (0->0), fcn. (464->10), ass. (0->52)
t87 = cos(qJ(5));
t102 = t87 * pkin(5) + pkin(4) + pkin(8);
t89 = cos(qJ(2));
t121 = t102 * t89;
t86 = sin(qJ(2));
t84 = sin(qJ(5));
t73 = t84 * pkin(5) + qJ(4);
t79 = qJ(6) + pkin(3) + pkin(9);
t85 = sin(qJ(3));
t88 = cos(qJ(3));
t97 = t73 * t85 + t79 * t88;
t123 = t121 - (pkin(2) + t97) * t86;
t81 = sin(pkin(6));
t106 = t88 * t81;
t113 = t81 * t85;
t83 = cos(pkin(6));
t110 = t83 * t86;
t80 = sin(pkin(10));
t82 = cos(pkin(10));
t62 = t80 * t110 - t82 * t89;
t64 = t82 * t110 + t80 * t89;
t93 = -g(2) * (-t82 * t113 + t64 * t88) - g(3) * (t86 * t106 + t83 * t85) + g(1) * (-t80 * t113 + t62 * t88);
t100 = t88 * pkin(3) + qJ(4) * t85;
t114 = t89 * pkin(8);
t120 = (pkin(2) + t100) * t86 - t114;
t118 = g(1) * t80;
t117 = g(2) * t82;
t76 = t82 * pkin(2);
t112 = t81 * t89;
t111 = t82 * t86;
t109 = t83 * t89;
t108 = t85 * t86;
t107 = t85 * t89;
t105 = t80 * t81 * pkin(7) + t82 * pkin(1);
t103 = t83 * pkin(7) + qJ(1);
t99 = t85 * pkin(3) - qJ(4) * t88;
t98 = -t73 * t88 + t79 * t85;
t67 = t81 * t108 - t83 * t88;
t94 = g(1) * (t80 * t106 + t62 * t85) - g(2) * (t82 * t106 + t64 * t85) - g(3) * t67;
t63 = -t82 * t109 + t80 * t86;
t65 = t80 * t109 + t111;
t92 = -g(1) * t65 - g(2) * t63 + g(3) * t112;
t91 = t102 * t86 + t97 * t89;
t75 = t80 * pkin(1);
t74 = t80 * pkin(2);
t71 = t83 * t76;
t66 = t83 * t108 + t106;
t61 = t82 * t107 - t80 * t66;
t60 = t80 * t107 + t82 * t66;
t55 = -g(1) * (t61 * t87 - t65 * t84) - g(2) * (t60 * t87 - t63 * t84) - g(3) * (t84 * t112 + t67 * t87);
t54 = -g(1) * (t61 * t84 + t65 * t87) - g(2) * (t60 * t84 + t63 * t87) - g(3) * (-t87 * t112 + t67 * t84);
t1 = [-g(3) * qJ(1), 0, -g(3) * t81 * t86 + g(1) * t62 - g(2) * t64, -t92, 0, 0, 0, 0, 0, t93, -t94, t92, -t93, t94, -g(1) * ((t100 * t82 + t76) * t89 + pkin(8) * t111 + t105) - g(2) * ((t100 * t80 + t74) * t89 + (t80 * pkin(8) + t71) * t86 + t75) - g(3) * t103 + (-g(3) * t99 - (t100 * t86 - t114) * t117 + t120 * t118) * t83 + (-g(3) * t120 - t99 * t118 - (-pkin(7) - t99) * t117) * t81, 0, 0, 0, 0, 0, t54, t55, t54, t55, t93, -g(1) * (t76 * t89 + t105) - g(2) * (t71 * t86 + t74 * t89 + t75) - g(3) * (-t123 * t81 + t98 * t83 + t103) + (-g(1) * t91 + (-(-pkin(7) - t98) * t81 - (t97 * t86 - t121) * t83) * g(2)) * t82 + (-g(2) * t91 + (-t123 * t83 - t98 * t81) * g(1)) * t80;];
U_reg = t1;
