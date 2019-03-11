% Calculate inertial parameters regressor of potential energy for
% S6RPRRRR8
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d3,d4,d5,d6]';
% 
% Output:
% U_reg [1x(6*10)]
%   inertial parameter regressor of Potential energy

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 07:22
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U_reg = S6RPRRRR8_energypot_fixb_reg2_slag_vp(qJ, g, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(10,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRRRR8_energypot_fixb_reg2_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RPRRRR8_energypot_fixb_reg2_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RPRRRR8_energypot_fixb_reg2_slag_vp: pkin has to be [10x1] (double)');

%% Symbolic Calculation
% From energy_potential_fixb_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 07:21:12
% EndTime: 2019-03-09 07:21:12
% DurationCPUTime: 0.14s
% Computational Cost: add. (150->73), mult. (168->89), div. (0->0), fcn. (160->10), ass. (0->41)
t102 = g(3) * pkin(6);
t101 = pkin(2) + pkin(6);
t75 = sin(qJ(3));
t100 = pkin(3) * t75;
t73 = qJ(3) + qJ(4);
t67 = cos(t73);
t99 = pkin(9) * t67;
t98 = g(3) * t67;
t72 = qJ(5) + qJ(6);
t64 = sin(t72);
t79 = cos(qJ(1));
t97 = t64 * t79;
t65 = sin(t73);
t76 = sin(qJ(1));
t96 = t65 * t76;
t66 = cos(t72);
t95 = t66 * t76;
t94 = t66 * t79;
t74 = sin(qJ(5));
t93 = t74 * t76;
t92 = t74 * t79;
t77 = cos(qJ(5));
t91 = t76 * t77;
t90 = t77 * t79;
t89 = t79 * pkin(1) + t76 * qJ(2);
t78 = cos(qJ(3));
t88 = t78 * pkin(3) + t101;
t87 = t76 * t100 + t89;
t81 = -pkin(8) - pkin(7);
t86 = pkin(5) * t74 - t81;
t85 = -qJ(2) - t100;
t69 = t76 * pkin(1);
t84 = -t76 * t81 + t69;
t83 = -t79 * qJ(2) + t69;
t60 = g(1) * t76 - g(2) * t79;
t63 = pkin(5) * t77 + pkin(4);
t80 = -pkin(10) - pkin(9);
t82 = t63 * t65 + t67 * t80;
t61 = g(1) * t79 + g(2) * t76;
t59 = -g(3) * t65 + t60 * t67;
t1 = [0, 0, 0, 0, 0, 0, -t61, t60, -g(3), -t102, 0, 0, 0, 0, 0, 0, -g(3), t61, -t60, -g(1) * t89 - g(2) * t83 - t102, 0, 0, 0, 0, 0, 0, -g(3) * t78 - t60 * t75, g(3) * t75 - t60 * t78, -t61, -g(1) * (pkin(7) * t79 + t89) - g(2) * (pkin(7) * t76 + t83) - g(3) * t101, 0, 0, 0, 0, 0, 0, -t60 * t65 - t98, -t59, -t61, -g(1) * (-t79 * t81 + t87) - g(2) * (t85 * t79 + t84) - g(3) * t88, 0, 0, 0, 0, 0, 0, -g(1) * (t65 * t91 + t92) - g(2) * (-t65 * t90 + t93) - t77 * t98, -g(1) * (-t65 * t93 + t90) - g(2) * (t65 * t92 + t91) + t74 * t98, t59, -g(1) * (pkin(4) * t96 - t76 * t99 + t87) - g(2) * t84 - g(3) * (pkin(4) * t67 + pkin(9) * t65 + t88) + (g(1) * t81 - g(2) * (-pkin(4) * t65 + t85 + t99)) * t79, 0, 0, 0, 0, 0, 0, -g(1) * (t65 * t95 + t97) - g(2) * (t64 * t76 - t65 * t94) - t66 * t98, -g(1) * (-t64 * t96 + t94) - g(2) * (t65 * t97 + t95) + t64 * t98, t59, -g(1) * t87 - g(2) * t69 - g(3) * (t63 * t67 - t65 * t80 + t88) + (-g(1) * t82 - g(2) * t86) * t76 + (-g(1) * t86 - g(2) * (-t82 + t85)) * t79;];
U_reg  = t1;
