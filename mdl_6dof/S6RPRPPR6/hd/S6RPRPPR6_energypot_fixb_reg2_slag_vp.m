% Calculate inertial parameters regressor of potential energy for
% S6RPRPPR6
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d3,d6,theta4,theta5]';
% 
% Output:
% U_reg [1x(6*10)]
%   inertial parameter regressor of Potential energy

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 02:55
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U_reg = S6RPRPPR6_energypot_fixb_reg2_slag_vp(qJ, g, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(10,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRPPR6_energypot_fixb_reg2_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RPRPPR6_energypot_fixb_reg2_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RPRPPR6_energypot_fixb_reg2_slag_vp: pkin has to be [10x1] (double)');

%% Symbolic Calculation
% From energy_potential_fixb_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 02:54:37
% EndTime: 2019-03-09 02:54:37
% DurationCPUTime: 0.16s
% Computational Cost: add. (150->72), mult. (168->87), div. (0->0), fcn. (160->10), ass. (0->41)
t73 = qJ(3) + pkin(9);
t65 = sin(t73);
t67 = cos(t73);
t104 = pkin(4) * t65 - qJ(5) * t67;
t103 = g(3) * pkin(6);
t102 = pkin(2) + pkin(6);
t78 = sin(qJ(3));
t101 = pkin(3) * t78;
t99 = g(3) * t67;
t72 = pkin(10) + qJ(6);
t64 = sin(t72);
t81 = cos(qJ(1));
t98 = t64 * t81;
t66 = cos(t72);
t97 = t66 * t81;
t74 = sin(pkin(10));
t96 = t74 * t81;
t75 = cos(pkin(10));
t95 = t75 * t81;
t79 = sin(qJ(1));
t94 = t79 * t64;
t93 = t79 * t66;
t92 = t79 * t74;
t91 = t79 * t75;
t90 = t81 * pkin(1) + t79 * qJ(2);
t80 = cos(qJ(3));
t88 = t80 * pkin(3) + t102;
t87 = t79 * t101 + t90;
t76 = -qJ(4) - pkin(7);
t86 = pkin(5) * t74 - t76;
t85 = -qJ(2) - t101;
t69 = t79 * pkin(1);
t84 = -t79 * t76 + t69;
t83 = -qJ(2) * t81 + t69;
t60 = g(1) * t79 - g(2) * t81;
t63 = pkin(5) * t75 + pkin(4);
t77 = -pkin(8) - qJ(5);
t82 = t63 * t65 + t67 * t77;
t61 = g(1) * t81 + g(2) * t79;
t59 = -g(3) * t65 + t60 * t67;
t1 = [0, 0, 0, 0, 0, 0, -t61, t60, -g(3), -t103, 0, 0, 0, 0, 0, 0, -g(3), t61, -t60, -g(1) * t90 - g(2) * t83 - t103, 0, 0, 0, 0, 0, 0, -g(3) * t80 - t60 * t78, g(3) * t78 - t60 * t80, -t61, -g(1) * (pkin(7) * t81 + t90) - g(2) * (t79 * pkin(7) + t83) - g(3) * t102, 0, 0, 0, 0, 0, 0, -t60 * t65 - t99, -t59, -t61, -g(1) * (-t76 * t81 + t87) - g(2) * (t85 * t81 + t84) - g(3) * t88, 0, 0, 0, 0, 0, 0, -g(1) * (t65 * t91 + t96) - g(2) * (-t65 * t95 + t92) - t75 * t99, -g(1) * (-t65 * t92 + t95) - g(2) * (t65 * t96 + t91) + t74 * t99, t59, -g(1) * (t104 * t79 + t87) - g(2) * t84 - g(3) * (pkin(4) * t67 + qJ(5) * t65 + t88) + (g(1) * t76 - g(2) * (t85 - t104)) * t81, 0, 0, 0, 0, 0, 0, -g(1) * (t65 * t93 + t98) - g(2) * (-t65 * t97 + t94) - t66 * t99, -g(1) * (-t65 * t94 + t97) - g(2) * (t65 * t98 + t93) + t64 * t99, t59, -g(1) * t87 - g(2) * t69 - g(3) * (t63 * t67 - t65 * t77 + t88) + (-g(1) * t82 - g(2) * t86) * t79 + (-g(1) * t86 - g(2) * (-t82 + t85)) * t81;];
U_reg  = t1;
