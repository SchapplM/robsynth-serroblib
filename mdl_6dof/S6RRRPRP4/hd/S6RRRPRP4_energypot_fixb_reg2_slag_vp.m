% Calculate inertial parameters regressor of potential energy for
% S6RRRPRP4
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d3,d5]';
% 
% Output:
% U_reg [1x(6*10)]
%   inertial parameter regressor of Potential energy

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 16:46
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U_reg = S6RRRPRP4_energypot_fixb_reg2_slag_vp(qJ, g, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(9,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRPRP4_energypot_fixb_reg2_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RRRPRP4_energypot_fixb_reg2_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S6RRRPRP4_energypot_fixb_reg2_slag_vp: pkin has to be [9x1] (double)');

%% Symbolic Calculation
% From energy_potential_fixb_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 16:45:55
% EndTime: 2019-03-09 16:45:55
% DurationCPUTime: 0.12s
% Computational Cost: add. (177->62), mult. (202->72), div. (0->0), fcn. (198->8), ass. (0->41)
t121 = g(3) * pkin(6);
t94 = qJ(2) + qJ(3);
t90 = cos(t94);
t120 = g(3) * t90;
t96 = sin(qJ(2));
t119 = t96 * pkin(2) + pkin(6);
t97 = sin(qJ(1));
t118 = t90 * t97;
t95 = sin(qJ(5));
t117 = t97 * t95;
t98 = cos(qJ(5));
t116 = t97 * t98;
t100 = cos(qJ(1));
t101 = -pkin(8) - pkin(7);
t99 = cos(qJ(2));
t87 = t99 * pkin(2) + pkin(1);
t115 = t100 * t101 + t97 * t87;
t89 = sin(t94);
t114 = qJ(4) * t89;
t113 = t100 * t90;
t112 = t100 * t95;
t111 = t100 * t98;
t110 = t89 * pkin(3) + t119;
t109 = t100 * t87 - t97 * t101;
t108 = pkin(3) * t118 + t97 * t114 + t115;
t107 = g(1) * t100 + g(2) * t97;
t106 = -t90 * qJ(4) + t110;
t105 = pkin(3) * t113 + t100 * t114 + t109;
t104 = -t100 * pkin(4) + pkin(9) * t118 + t108;
t103 = t97 * pkin(4) + pkin(9) * t113 + t105;
t71 = -t89 * t111 + t117;
t73 = t89 * t116 + t112;
t102 = g(1) * t71 - g(2) * t73 + t98 * t120;
t85 = t89 * pkin(9);
t76 = g(1) * t97 - g(2) * t100;
t74 = t89 * t117 - t111;
t72 = t89 * t112 + t116;
t70 = g(3) * t89 + t107 * t90;
t69 = t107 * t89 - t120;
t68 = -g(1) * t72 - g(2) * t74 + t95 * t120;
t1 = [0, 0, 0, 0, 0, 0, -t107, t76, -g(3), -t121, 0, 0, 0, 0, 0, 0, -g(3) * t96 - t107 * t99, -g(3) * t99 + t107 * t96, -t76, -g(1) * (t100 * pkin(1) + t97 * pkin(7)) - g(2) * (t97 * pkin(1) - t100 * pkin(7)) - t121, 0, 0, 0, 0, 0, 0, -t70, t69, -t76, -g(1) * t109 - g(2) * t115 - g(3) * t119, 0, 0, 0, 0, 0, 0, -t76, t70, -t69, -g(1) * t105 - g(2) * t108 - g(3) * t106, 0, 0, 0, 0, 0, 0, t68, t102, -t70, -g(1) * t103 - g(2) * t104 - g(3) * (t106 + t85) 0, 0, 0, 0, 0, 0, t68, -t70, -t102, -g(1) * (t72 * pkin(5) + t71 * qJ(6) + t103) - g(2) * (t74 * pkin(5) - t73 * qJ(6) + t104) - g(3) * (t85 + t110) - (-pkin(5) * t95 + qJ(6) * t98 - qJ(4)) * t120;];
U_reg  = t1;
