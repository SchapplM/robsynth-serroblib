% Calculate inertial parameters regressor of potential energy for
% S6RRPPRP1
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d5,theta3,theta4]';
% 
% Output:
% U_reg [1x(6*10)]
%   inertial parameter regressor of Potential energy

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 08:28
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U_reg = S6RRPPRP1_energypot_fixb_reg2_slag_vp(qJ, g, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(10,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPPRP1_energypot_fixb_reg2_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RRPPRP1_energypot_fixb_reg2_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RRPPRP1_energypot_fixb_reg2_slag_vp: pkin has to be [10x1] (double)');

%% Symbolic Calculation
% From energy_potential_fixb_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 08:27:46
% EndTime: 2019-03-09 08:27:46
% DurationCPUTime: 0.11s
% Computational Cost: add. (215->71), mult. (215->90), div. (0->0), fcn. (215->10), ass. (0->42)
t120 = g(3) * pkin(6);
t94 = qJ(2) + pkin(9);
t88 = sin(t94);
t119 = g(3) * t88;
t99 = sin(qJ(2));
t118 = t99 * pkin(2) + pkin(6);
t97 = -pkin(8) - qJ(4);
t117 = t88 * t97;
t100 = sin(qJ(1));
t102 = cos(qJ(1));
t101 = cos(qJ(2));
t86 = t101 * pkin(2) + pkin(1);
t98 = -pkin(7) - qJ(3);
t116 = t100 * t86 + t102 * t98;
t90 = cos(t94);
t115 = t100 * t90;
t95 = sin(pkin(10));
t114 = t100 * t95;
t96 = cos(pkin(10));
t113 = t100 * t96;
t112 = t102 * t90;
t111 = t102 * t95;
t110 = t102 * t96;
t84 = t96 * pkin(4) + pkin(3);
t109 = t88 * t84 + t90 * t97 + t118;
t108 = -t100 * t98 + t102 * t86;
t107 = pkin(3) * t90 + qJ(4) * t88;
t106 = g(1) * t102 + g(2) * t100;
t93 = pkin(10) + qJ(5);
t87 = sin(t93);
t89 = cos(t93);
t70 = t102 * t89 + t87 * t115;
t72 = -t100 * t89 + t87 * t112;
t105 = g(1) * t72 + g(2) * t70 + t87 * t119;
t104 = pkin(4) * t114 - t102 * t117 + t84 * t112 + t108;
t103 = -pkin(4) * t111 - t100 * t117 + t84 * t115 + t116;
t77 = g(1) * t100 - g(2) * t102;
t73 = t100 * t87 + t89 * t112;
t71 = -t102 * t87 + t89 * t115;
t69 = -g(3) * t90 + t106 * t88;
t68 = -g(1) * t73 - g(2) * t71 - t89 * t119;
t1 = [0, 0, 0, 0, 0, 0, -t106, t77, -g(3), -t120, 0, 0, 0, 0, 0, 0, -g(3) * t99 - t106 * t101, -g(3) * t101 + t106 * t99, -t77, -g(1) * (t102 * pkin(1) + t100 * pkin(7)) - g(2) * (t100 * pkin(1) - t102 * pkin(7)) - t120, 0, 0, 0, 0, 0, 0, -t106 * t90 - t119, t69, -t77, -g(1) * t108 - g(2) * t116 - g(3) * t118, 0, 0, 0, 0, 0, 0, -g(1) * (t90 * t110 + t114) - g(2) * (t90 * t113 - t111) - t96 * t119, -g(1) * (-t90 * t111 + t113) - g(2) * (-t90 * t114 - t110) + t95 * t119, -t69, -g(1) * (t107 * t102 + t108) - g(2) * (t107 * t100 + t116) - g(3) * (t88 * pkin(3) - t90 * qJ(4) + t118) 0, 0, 0, 0, 0, 0, t68, t105, -t69, -g(1) * t104 - g(2) * t103 - g(3) * t109, 0, 0, 0, 0, 0, 0, t68, -t69, -t105, -g(1) * (t73 * pkin(5) + t72 * qJ(6) + t104) - g(2) * (t71 * pkin(5) + t70 * qJ(6) + t103) - g(3) * ((pkin(5) * t89 + qJ(6) * t87) * t88 + t109);];
U_reg  = t1;
