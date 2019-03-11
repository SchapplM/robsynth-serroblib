% Calculate inertial parameters regressor of potential energy for
% S6RRPPRP2
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d5,theta3]';
% 
% Output:
% U_reg [1x(6*10)]
%   inertial parameter regressor of Potential energy

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 08:32
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U_reg = S6RRPPRP2_energypot_fixb_reg2_slag_vp(qJ, g, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(9,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPPRP2_energypot_fixb_reg2_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RRPPRP2_energypot_fixb_reg2_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S6RRPPRP2_energypot_fixb_reg2_slag_vp: pkin has to be [9x1] (double)');

%% Symbolic Calculation
% From energy_potential_fixb_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 08:31:56
% EndTime: 2019-03-09 08:31:56
% DurationCPUTime: 0.13s
% Computational Cost: add. (172->67), mult. (192->77), div. (0->0), fcn. (184->8), ass. (0->40)
t87 = qJ(2) + pkin(9);
t83 = sin(t87);
t103 = qJ(4) * t83;
t84 = cos(t87);
t95 = cos(qJ(1));
t109 = t84 * t95;
t115 = pkin(3) * t109 + t95 * t103;
t114 = g(3) * pkin(6);
t90 = sin(qJ(5));
t113 = pkin(5) * t90;
t112 = g(3) * t84;
t91 = sin(qJ(2));
t111 = t91 * pkin(2) + pkin(6);
t92 = sin(qJ(1));
t110 = t84 * t92;
t108 = t92 * t90;
t93 = cos(qJ(5));
t107 = t92 * t93;
t106 = t95 * t90;
t105 = t95 * t93;
t94 = cos(qJ(2));
t82 = t94 * pkin(2) + pkin(1);
t89 = -pkin(7) - qJ(3);
t104 = t92 * t82 + t95 * t89;
t102 = t83 * pkin(3) + t111;
t101 = t83 * t106;
t78 = t95 * t82;
t100 = t78 + t115;
t99 = -t92 * t89 + t78;
t98 = pkin(3) * t110 + t92 * t103 + t104;
t97 = g(1) * t95 + g(2) * t92;
t96 = -t84 * qJ(4) + t102;
t88 = -qJ(6) - pkin(8);
t81 = t93 * pkin(5) + pkin(4);
t74 = g(1) * t92 - g(2) * t95;
t71 = g(3) * t83 + t97 * t84;
t70 = t97 * t83 - t112;
t69 = -g(1) * (t101 + t107) - g(2) * (t83 * t108 - t105) + t90 * t112;
t68 = -g(1) * (t83 * t105 - t108) - g(2) * (t83 * t107 + t106) + t93 * t112;
t1 = [0, 0, 0, 0, 0, 0, -t97, t74, -g(3), -t114, 0, 0, 0, 0, 0, 0, -g(3) * t91 - t97 * t94, -g(3) * t94 + t97 * t91, -t74, -g(1) * (t95 * pkin(1) + t92 * pkin(7)) - g(2) * (t92 * pkin(1) - t95 * pkin(7)) - t114, 0, 0, 0, 0, 0, 0, -t71, t70, -t74, -g(1) * t99 - g(2) * t104 - g(3) * t111, 0, 0, 0, 0, 0, 0, -t74, t71, -t70, -g(1) * (t99 + t115) - g(2) * t98 - g(3) * t96, 0, 0, 0, 0, 0, 0, t69, t68, -t71, -g(1) * (pkin(8) * t109 + (pkin(4) - t89) * t92 + t100) - g(2) * (-t95 * pkin(4) + pkin(8) * t110 + t98) - g(3) * (t83 * pkin(8) + t96) 0, 0, 0, 0, 0, 0, t69, t68, -t71, -g(1) * (pkin(5) * t101 - t88 * t109 + t100) - g(2) * (-t95 * t81 + t98) - g(3) * (-t83 * t88 + (-qJ(4) - t113) * t84 + t102) + (-g(1) * (t81 - t89) - g(2) * (t83 * t113 - t84 * t88)) * t92;];
U_reg  = t1;
