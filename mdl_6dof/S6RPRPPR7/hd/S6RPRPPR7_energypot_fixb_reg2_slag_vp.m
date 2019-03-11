% Calculate inertial parameters regressor of potential energy for
% S6RPRPPR7
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d3,d6,theta4]';
% 
% Output:
% U_reg [1x(6*10)]
%   inertial parameter regressor of Potential energy

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 02:57
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U_reg = S6RPRPPR7_energypot_fixb_reg2_slag_vp(qJ, g, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(9,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRPPR7_energypot_fixb_reg2_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RPRPPR7_energypot_fixb_reg2_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S6RPRPPR7_energypot_fixb_reg2_slag_vp: pkin has to be [9x1] (double)');

%% Symbolic Calculation
% From energy_potential_fixb_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 02:57:25
% EndTime: 2019-03-09 02:57:25
% DurationCPUTime: 0.14s
% Computational Cost: add. (132->66), mult. (155->73), div. (0->0), fcn. (143->8), ass. (0->39)
t102 = g(3) * pkin(6);
t101 = pkin(2) + pkin(6);
t76 = sin(qJ(3));
t100 = pkin(3) * t76;
t73 = qJ(3) + pkin(9);
t67 = sin(t73);
t99 = pkin(8) * t67;
t80 = cos(qJ(1));
t98 = g(2) * t80;
t97 = g(3) * t67;
t74 = -qJ(4) - pkin(7);
t96 = pkin(5) - t74;
t77 = sin(qJ(1));
t95 = t67 * t77;
t75 = sin(qJ(6));
t94 = t77 * t75;
t78 = cos(qJ(6));
t93 = t77 * t78;
t92 = t80 * t75;
t91 = t80 * t78;
t90 = t80 * pkin(1) + t77 * qJ(2);
t68 = cos(t73);
t89 = qJ(5) * t68;
t79 = cos(qJ(3));
t88 = t79 * pkin(3) + t101;
t87 = t77 * t100 + t90;
t86 = -qJ(2) - t100;
t70 = t77 * pkin(1);
t85 = -t77 * t74 + t70;
t84 = -t80 * qJ(2) + t70;
t83 = pkin(4) * t95 + t87;
t82 = t68 * pkin(4) + t67 * qJ(5) + t88;
t60 = g(1) * t77 - t98;
t81 = -pkin(4) * t67 + t86;
t61 = g(1) * t80 + g(2) * t77;
t59 = t80 * t89;
t58 = t60 * t68 - t97;
t57 = g(1) * t95 + g(3) * t68 - t67 * t98;
t1 = [0, 0, 0, 0, 0, 0, -t61, t60, -g(3), -t102, 0, 0, 0, 0, 0, 0, -g(3), t61, -t60, -g(1) * t90 - g(2) * t84 - t102, 0, 0, 0, 0, 0, 0, -g(3) * t79 - t60 * t76, g(3) * t76 - t60 * t79, -t61, -g(1) * (t80 * pkin(7) + t90) - g(2) * (t77 * pkin(7) + t84) - g(3) * t101, 0, 0, 0, 0, 0, 0, -t57, -t58, -t61, -g(1) * (-t80 * t74 + t87) - g(2) * (t86 * t80 + t85) - g(3) * t88, 0, 0, 0, 0, 0, 0, -t61, t57, t58, -g(1) * (-t77 * t89 + t83) - g(2) * (t59 + t85) - g(3) * t82 + (g(1) * t74 - g(2) * t81) * t80, 0, 0, 0, 0, 0, 0, -g(1) * (-t68 * t94 + t91) - g(2) * (t68 * t92 + t93) - t75 * t97, -g(1) * (-t68 * t93 - t92) - g(2) * (t68 * t91 - t94) - t78 * t97, -t57, -g(1) * t83 - g(2) * (t59 + t70) - g(3) * (t68 * pkin(8) + t82) + (-g(1) * (-t89 + t99) - g(2) * t96) * t77 + (-g(1) * t96 - g(2) * (t81 - t99)) * t80;];
U_reg  = t1;
