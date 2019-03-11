% Calculate inertial parameters regressor of potential energy for
% S6RPRPRR8
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d3,d5,d6,theta4]';
% 
% Output:
% U_reg [1x(6*10)]
%   inertial parameter regressor of Potential energy

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 04:00
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U_reg = S6RPRPRR8_energypot_fixb_reg2_slag_vp(qJ, g, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(10,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRPRR8_energypot_fixb_reg2_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RPRPRR8_energypot_fixb_reg2_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RPRPRR8_energypot_fixb_reg2_slag_vp: pkin has to be [10x1] (double)');

%% Symbolic Calculation
% From energy_potential_fixb_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 03:59:30
% EndTime: 2019-03-09 03:59:30
% DurationCPUTime: 0.16s
% Computational Cost: add. (150->72), mult. (168->87), div. (0->0), fcn. (160->10), ass. (0->41)
t72 = qJ(3) + pkin(10);
t64 = sin(t72);
t65 = cos(t72);
t104 = pkin(4) * t64 - pkin(8) * t65;
t103 = g(3) * pkin(6);
t102 = pkin(2) + pkin(6);
t76 = sin(qJ(3));
t101 = pkin(3) * t76;
t98 = g(3) * t65;
t73 = qJ(5) + qJ(6);
t66 = sin(t73);
t77 = sin(qJ(1));
t97 = t66 * t77;
t80 = cos(qJ(1));
t96 = t66 * t80;
t67 = cos(t73);
t95 = t67 * t77;
t94 = t67 * t80;
t75 = sin(qJ(5));
t93 = t75 * t77;
t92 = t75 * t80;
t78 = cos(qJ(5));
t91 = t77 * t78;
t90 = t78 * t80;
t89 = t80 * pkin(1) + t77 * qJ(2);
t79 = cos(qJ(3));
t88 = t79 * pkin(3) + t102;
t87 = t77 * t101 + t89;
t74 = -qJ(4) - pkin(7);
t86 = pkin(5) * t75 - t74;
t85 = -qJ(2) - t101;
t69 = t77 * pkin(1);
t84 = -t77 * t74 + t69;
t83 = -t80 * qJ(2) + t69;
t60 = g(1) * t77 - g(2) * t80;
t63 = pkin(5) * t78 + pkin(4);
t81 = -pkin(9) - pkin(8);
t82 = t63 * t64 + t65 * t81;
t61 = g(1) * t80 + g(2) * t77;
t59 = -g(3) * t64 + t60 * t65;
t1 = [0, 0, 0, 0, 0, 0, -t61, t60, -g(3), -t103, 0, 0, 0, 0, 0, 0, -g(3), t61, -t60, -g(1) * t89 - g(2) * t83 - t103, 0, 0, 0, 0, 0, 0, -g(3) * t79 - t60 * t76, g(3) * t76 - t60 * t79, -t61, -g(1) * (pkin(7) * t80 + t89) - g(2) * (pkin(7) * t77 + t83) - g(3) * t102, 0, 0, 0, 0, 0, 0, -t60 * t64 - t98, -t59, -t61, -g(1) * (-t80 * t74 + t87) - g(2) * (t85 * t80 + t84) - g(3) * t88, 0, 0, 0, 0, 0, 0, -g(1) * (t64 * t91 + t92) - g(2) * (-t64 * t90 + t93) - t78 * t98, -g(1) * (-t64 * t93 + t90) - g(2) * (t64 * t92 + t91) + t75 * t98, t59, -g(1) * (t104 * t77 + t87) - g(2) * t84 - g(3) * (pkin(4) * t65 + pkin(8) * t64 + t88) + (g(1) * t74 - g(2) * (t85 - t104)) * t80, 0, 0, 0, 0, 0, 0, -g(1) * (t64 * t95 + t96) - g(2) * (-t64 * t94 + t97) - t67 * t98, -g(1) * (-t64 * t97 + t94) - g(2) * (t64 * t96 + t95) + t66 * t98, t59, -g(1) * t87 - g(2) * t69 - g(3) * (t65 * t63 - t64 * t81 + t88) + (-g(1) * t82 - g(2) * t86) * t77 + (-g(1) * t86 - g(2) * (-t82 + t85)) * t80;];
U_reg  = t1;
