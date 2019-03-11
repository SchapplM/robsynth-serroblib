% Calculate inertial parameters regressor of potential energy for
% S6RPRRRP1
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d3,d4,d5,theta2]';
% 
% Output:
% U_reg [1x(6*10)]
%   inertial parameter regressor of Potential energy

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 05:57
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U_reg = S6RPRRRP1_energypot_fixb_reg2_slag_vp(qJ, g, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(10,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRRRP1_energypot_fixb_reg2_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RPRRRP1_energypot_fixb_reg2_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RPRRRP1_energypot_fixb_reg2_slag_vp: pkin has to be [10x1] (double)');

%% Symbolic Calculation
% From energy_potential_fixb_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 05:57:05
% EndTime: 2019-03-09 05:57:05
% DurationCPUTime: 0.13s
% Computational Cost: add. (215->59), mult. (173->71), div. (0->0), fcn. (169->10), ass. (0->39)
t85 = qJ(3) + qJ(4);
t79 = sin(t85);
t80 = cos(t85);
t109 = pkin(4) * t80 + pkin(9) * t79;
t106 = g(3) * t79;
t86 = qJ(2) + pkin(6);
t105 = g(3) * t86;
t87 = sin(qJ(5));
t104 = t80 * t87;
t90 = cos(qJ(5));
t103 = t80 * t90;
t91 = cos(qJ(3));
t76 = pkin(3) * t91 + pkin(2);
t84 = qJ(1) + pkin(10);
t77 = sin(t84);
t78 = cos(t84);
t89 = sin(qJ(1));
t82 = t89 * pkin(1);
t93 = -pkin(8) - pkin(7);
t102 = t77 * t76 + t78 * t93 + t82;
t88 = sin(qJ(3));
t101 = t88 * pkin(3) + t86;
t92 = cos(qJ(1));
t83 = t92 * pkin(1);
t100 = t78 * t76 - t77 * t93 + t83;
t99 = t109 * t77 + t102;
t98 = g(1) * t78 + g(2) * t77;
t97 = -g(1) * t92 - g(2) * t89;
t96 = t79 * pkin(4) - t80 * pkin(9) + t101;
t95 = t109 * t78 + t100;
t60 = t77 * t104 + t78 * t90;
t62 = t78 * t104 - t77 * t90;
t94 = g(1) * t62 + g(2) * t60 + t87 * t106;
t64 = g(1) * t77 - g(2) * t78;
t63 = t78 * t103 + t77 * t87;
t61 = t77 * t103 - t78 * t87;
t59 = -g(3) * t80 + t98 * t79;
t58 = -g(1) * t63 - g(2) * t61 - t90 * t106;
t1 = [0, 0, 0, 0, 0, 0, t97, g(1) * t89 - g(2) * t92, -g(3), -g(3) * pkin(6), 0, 0, 0, 0, 0, 0, -t98, t64, -g(3), t97 * pkin(1) - t105, 0, 0, 0, 0, 0, 0, -g(3) * t88 - t98 * t91, -g(3) * t91 + t98 * t88, -t64, -g(1) * (pkin(2) * t78 + pkin(7) * t77 + t83) - g(2) * (pkin(2) * t77 - pkin(7) * t78 + t82) - t105, 0, 0, 0, 0, 0, 0, -t98 * t80 - t106, t59, -t64, -g(1) * t100 - g(2) * t102 - g(3) * t101, 0, 0, 0, 0, 0, 0, t58, t94, -t59, -g(1) * t95 - g(2) * t99 - g(3) * t96, 0, 0, 0, 0, 0, 0, t58, -t59, -t94, -g(1) * (t63 * pkin(5) + t62 * qJ(6) + t95) - g(2) * (pkin(5) * t61 + qJ(6) * t60 + t99) - g(3) * ((pkin(5) * t90 + qJ(6) * t87) * t79 + t96);];
U_reg  = t1;
