% Calculate inertial parameters regressor of potential energy for
% S6RPRPRP2
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d3,d5,theta2,theta4]';
% 
% Output:
% U_reg [1x(6*10)]
%   inertial parameter regressor of Potential energy

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 03:06
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U_reg = S6RPRPRP2_energypot_fixb_reg2_slag_vp(qJ, g, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(10,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRPRP2_energypot_fixb_reg2_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RPRPRP2_energypot_fixb_reg2_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RPRPRP2_energypot_fixb_reg2_slag_vp: pkin has to be [10x1] (double)');

%% Symbolic Calculation
% From energy_potential_fixb_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 03:06:12
% EndTime: 2019-03-09 03:06:12
% DurationCPUTime: 0.13s
% Computational Cost: add. (215->59), mult. (173->69), div. (0->0), fcn. (169->10), ass. (0->41)
t85 = qJ(3) + pkin(10);
t78 = sin(t85);
t80 = cos(t85);
t112 = pkin(4) * t80 + pkin(8) * t78;
t109 = g(3) * t78;
t88 = qJ(2) + pkin(6);
t108 = g(3) * t88;
t86 = qJ(1) + pkin(9);
t79 = sin(t86);
t89 = sin(qJ(5));
t107 = t79 * t89;
t92 = cos(qJ(5));
t106 = t79 * t92;
t81 = cos(t86);
t105 = t81 * t89;
t104 = t81 * t92;
t93 = cos(qJ(3));
t77 = t93 * pkin(3) + pkin(2);
t91 = sin(qJ(1));
t83 = t91 * pkin(1);
t87 = -qJ(4) - pkin(7);
t103 = t79 * t77 + t81 * t87 + t83;
t90 = sin(qJ(3));
t102 = t90 * pkin(3) + t88;
t94 = cos(qJ(1));
t84 = t94 * pkin(1);
t101 = t81 * t77 - t79 * t87 + t84;
t100 = t112 * t79 + t103;
t99 = g(1) * t81 + g(2) * t79;
t98 = -g(1) * t94 - g(2) * t91;
t97 = t78 * pkin(4) - t80 * pkin(8) + t102;
t96 = t112 * t81 + t101;
t61 = t80 * t107 + t104;
t63 = t80 * t105 - t106;
t95 = g(1) * t63 + g(2) * t61 + t89 * t109;
t65 = g(1) * t79 - g(2) * t81;
t64 = t80 * t104 + t107;
t62 = t80 * t106 - t105;
t60 = -g(3) * t80 + t99 * t78;
t59 = -g(1) * t64 - g(2) * t62 - t92 * t109;
t1 = [0, 0, 0, 0, 0, 0, t98, g(1) * t91 - g(2) * t94, -g(3), -g(3) * pkin(6), 0, 0, 0, 0, 0, 0, -t99, t65, -g(3), t98 * pkin(1) - t108, 0, 0, 0, 0, 0, 0, -g(3) * t90 - t99 * t93, -g(3) * t93 + t99 * t90, -t65, -g(1) * (t81 * pkin(2) + t79 * pkin(7) + t84) - g(2) * (t79 * pkin(2) - t81 * pkin(7) + t83) - t108, 0, 0, 0, 0, 0, 0, -t99 * t80 - t109, t60, -t65, -g(1) * t101 - g(2) * t103 - g(3) * t102, 0, 0, 0, 0, 0, 0, t59, t95, -t60, -g(1) * t96 - g(2) * t100 - g(3) * t97, 0, 0, 0, 0, 0, 0, t59, -t60, -t95, -g(1) * (t64 * pkin(5) + t63 * qJ(6) + t96) - g(2) * (t62 * pkin(5) + t61 * qJ(6) + t100) - g(3) * ((pkin(5) * t92 + qJ(6) * t89) * t78 + t97);];
U_reg  = t1;
