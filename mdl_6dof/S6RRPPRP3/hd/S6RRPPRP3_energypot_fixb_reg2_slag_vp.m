% Calculate inertial parameters regressor of potential energy for
% S6RRPPRP3
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d5]';
% 
% Output:
% U_reg [1x(6*10)]
%   inertial parameter regressor of Potential energy

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 08:36
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U_reg = S6RRPPRP3_energypot_fixb_reg2_slag_vp(qJ, g, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(8,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPPRP3_energypot_fixb_reg2_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RRPPRP3_energypot_fixb_reg2_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S6RRPPRP3_energypot_fixb_reg2_slag_vp: pkin has to be [8x1] (double)');

%% Symbolic Calculation
% From energy_potential_fixb_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 08:35:32
% EndTime: 2019-03-09 08:35:32
% DurationCPUTime: 0.16s
% Computational Cost: add. (118->67), mult. (214->73), div. (0->0), fcn. (206->6), ass. (0->41)
t83 = sin(qJ(1));
t85 = cos(qJ(2));
t103 = t83 * t85;
t82 = sin(qJ(2));
t98 = qJ(3) * t82;
t111 = pkin(2) * t103 + t83 * t98;
t86 = cos(qJ(1));
t110 = pkin(3) * t103 + t86 * qJ(4);
t109 = g(3) * pkin(6);
t81 = sin(qJ(5));
t108 = pkin(5) * t81;
t107 = g(3) * t85;
t106 = t82 * pkin(2) + pkin(6);
t105 = t83 * t81;
t84 = cos(qJ(5));
t104 = t83 * t84;
t102 = t85 * t86;
t101 = t86 * t81;
t100 = t86 * t84;
t99 = t86 * pkin(1) + t83 * pkin(7);
t73 = t82 * pkin(3);
t97 = t73 + t106;
t76 = t83 * pkin(1);
t96 = -t86 * pkin(7) + t76;
t95 = pkin(2) * t102 + t86 * t98 + t99;
t94 = -t85 * qJ(3) + t106;
t93 = pkin(3) * t102 + t95;
t92 = pkin(4) * t82 + pkin(8) * t85;
t91 = g(1) * t86 + g(2) * t83;
t71 = t84 * pkin(5) + pkin(4);
t80 = -qJ(6) - pkin(8);
t90 = t71 * t82 - t80 * t85;
t89 = t96 + t111;
t88 = t89 + t110;
t87 = -t83 * qJ(4) + t93;
t64 = g(1) * t83 - g(2) * t86;
t63 = g(3) * t82 + t91 * t85;
t62 = t91 * t82 - t107;
t61 = -g(1) * (t82 * t100 - t105) - g(2) * (t82 * t104 + t101) + t84 * t107;
t60 = -g(1) * (-t82 * t101 - t104) - g(2) * (-t82 * t105 + t100) - t81 * t107;
t1 = [0, 0, 0, 0, 0, 0, -t91, t64, -g(3), -t109, 0, 0, 0, 0, 0, 0, -t63, t62, -t64, -g(1) * t99 - g(2) * t96 - t109, 0, 0, 0, 0, 0, 0, -t63, -t64, -t62, -g(1) * t95 - g(2) * t89 - g(3) * t94, 0, 0, 0, 0, 0, 0, -t62, t63, t64, -g(1) * t87 - g(2) * t88 - g(3) * (t73 + t94) 0, 0, 0, 0, 0, 0, t61, t60, -t63, -g(1) * (t92 * t86 + t87) - g(2) * (t92 * t83 + t88) - g(3) * (t82 * pkin(8) + (-pkin(4) - qJ(3)) * t85 + t97) 0, 0, 0, 0, 0, 0, t61, t60, -t63, -g(1) * t93 - g(2) * (t76 + t110 + t111) - g(3) * (-t82 * t80 + (-qJ(3) - t71) * t85 + t97) + (-g(1) * t90 - g(2) * (-pkin(7) + t108)) * t86 + (-g(1) * (-qJ(4) - t108) - g(2) * t90) * t83;];
U_reg  = t1;
