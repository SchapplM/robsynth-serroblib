% Calculate inertial parameters regressor of potential energy for
% S6RPPRRP2
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d4,d5,theta2,theta3]';
% 
% Output:
% U_reg [1x(6*10)]
%   inertial parameter regressor of Potential energy

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 02:01
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U_reg = S6RPPRRP2_energypot_fixb_reg2_slag_vp(qJ, g, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(10,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPPRRP2_energypot_fixb_reg2_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RPPRRP2_energypot_fixb_reg2_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RPPRRP2_energypot_fixb_reg2_slag_vp: pkin has to be [10x1] (double)');

%% Symbolic Calculation
% From energy_potential_fixb_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 02:01:16
% EndTime: 2019-03-09 02:01:16
% DurationCPUTime: 0.16s
% Computational Cost: add. (215->59), mult. (173->69), div. (0->0), fcn. (169->10), ass. (0->41)
t82 = pkin(10) + qJ(4);
t75 = sin(t82);
t77 = cos(t82);
t109 = pkin(4) * t77 + pkin(8) * t75;
t106 = g(3) * t75;
t86 = qJ(2) + pkin(6);
t105 = g(3) * t86;
t83 = qJ(1) + pkin(9);
t76 = sin(t83);
t88 = sin(qJ(5));
t104 = t76 * t88;
t90 = cos(qJ(5));
t103 = t76 * t90;
t78 = cos(t83);
t102 = t78 * t88;
t101 = t78 * t90;
t85 = cos(pkin(10));
t74 = pkin(3) * t85 + pkin(2);
t89 = sin(qJ(1));
t80 = t89 * pkin(1);
t87 = -pkin(7) - qJ(3);
t100 = t76 * t74 + t78 * t87 + t80;
t84 = sin(pkin(10));
t99 = t84 * pkin(3) + t86;
t91 = cos(qJ(1));
t81 = t91 * pkin(1);
t98 = t78 * t74 - t76 * t87 + t81;
t97 = t109 * t76 + t100;
t96 = g(1) * t78 + g(2) * t76;
t95 = -g(1) * t91 - g(2) * t89;
t94 = t75 * pkin(4) - t77 * pkin(8) + t99;
t93 = t109 * t78 + t98;
t58 = t77 * t104 + t101;
t60 = t77 * t102 - t103;
t92 = g(1) * t60 + g(2) * t58 + t88 * t106;
t62 = g(1) * t76 - g(2) * t78;
t61 = t77 * t101 + t104;
t59 = t77 * t103 - t102;
t57 = -g(3) * t77 + t96 * t75;
t56 = -g(1) * t61 - g(2) * t59 - t90 * t106;
t1 = [0, 0, 0, 0, 0, 0, t95, g(1) * t89 - g(2) * t91, -g(3), -g(3) * pkin(6), 0, 0, 0, 0, 0, 0, -t96, t62, -g(3), t95 * pkin(1) - t105, 0, 0, 0, 0, 0, 0, -g(3) * t84 - t96 * t85, -g(3) * t85 + t96 * t84, -t62, -g(1) * (pkin(2) * t78 + qJ(3) * t76 + t81) - g(2) * (pkin(2) * t76 - qJ(3) * t78 + t80) - t105, 0, 0, 0, 0, 0, 0, -t96 * t77 - t106, t57, -t62, -g(1) * t98 - g(2) * t100 - g(3) * t99, 0, 0, 0, 0, 0, 0, t56, t92, -t57, -g(1) * t93 - g(2) * t97 - g(3) * t94, 0, 0, 0, 0, 0, 0, t56, -t57, -t92, -g(1) * (pkin(5) * t61 + qJ(6) * t60 + t93) - g(2) * (pkin(5) * t59 + qJ(6) * t58 + t97) - g(3) * ((pkin(5) * t90 + qJ(6) * t88) * t75 + t94);];
U_reg  = t1;
