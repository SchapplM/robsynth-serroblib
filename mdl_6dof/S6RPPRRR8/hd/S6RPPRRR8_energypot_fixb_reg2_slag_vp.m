% Calculate inertial parameters regressor of potential energy for
% S6RPPRRR8
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d4,d5,d6,theta3]';
% 
% Output:
% U_reg [1x(6*10)]
%   inertial parameter regressor of Potential energy

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 02:36
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U_reg = S6RPPRRR8_energypot_fixb_reg2_slag_vp(qJ, g, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(10,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPPRRR8_energypot_fixb_reg2_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RPPRRR8_energypot_fixb_reg2_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RPPRRR8_energypot_fixb_reg2_slag_vp: pkin has to be [10x1] (double)');

%% Symbolic Calculation
% From energy_potential_fixb_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 02:36:19
% EndTime: 2019-03-09 02:36:19
% DurationCPUTime: 0.16s
% Computational Cost: add. (150->72), mult. (168->87), div. (0->0), fcn. (160->10), ass. (0->41)
t69 = pkin(10) + qJ(4);
t61 = sin(t69);
t62 = cos(t69);
t101 = pkin(4) * t61 - pkin(8) * t62;
t100 = g(3) * pkin(6);
t99 = pkin(2) + pkin(6);
t71 = sin(pkin(10));
t98 = pkin(3) * t71;
t95 = g(3) * t62;
t70 = qJ(5) + qJ(6);
t63 = sin(t70);
t75 = sin(qJ(1));
t94 = t75 * t63;
t64 = cos(t70);
t93 = t75 * t64;
t74 = sin(qJ(5));
t92 = t75 * t74;
t76 = cos(qJ(5));
t91 = t75 * t76;
t77 = cos(qJ(1));
t90 = t77 * t63;
t89 = t77 * t64;
t88 = t77 * t74;
t87 = t77 * t76;
t86 = t77 * pkin(1) + t75 * qJ(2);
t72 = cos(pkin(10));
t85 = t72 * pkin(3) + t99;
t84 = t75 * t98 + t86;
t73 = -pkin(7) - qJ(3);
t83 = pkin(5) * t74 - t73;
t82 = -qJ(2) - t98;
t67 = t75 * pkin(1);
t81 = -t75 * t73 + t67;
t80 = -t77 * qJ(2) + t67;
t57 = g(1) * t75 - g(2) * t77;
t60 = t76 * pkin(5) + pkin(4);
t78 = -pkin(9) - pkin(8);
t79 = t60 * t61 + t62 * t78;
t58 = g(1) * t77 + g(2) * t75;
t56 = -g(3) * t61 + t57 * t62;
t1 = [0, 0, 0, 0, 0, 0, -t58, t57, -g(3), -t100, 0, 0, 0, 0, 0, 0, -g(3), t58, -t57, -g(1) * t86 - g(2) * t80 - t100, 0, 0, 0, 0, 0, 0, -g(3) * t72 - t57 * t71, g(3) * t71 - t57 * t72, -t58, -g(1) * (t77 * qJ(3) + t86) - g(2) * (t75 * qJ(3) + t80) - g(3) * t99, 0, 0, 0, 0, 0, 0, -t57 * t61 - t95, -t56, -t58, -g(1) * (-t77 * t73 + t84) - g(2) * (t82 * t77 + t81) - g(3) * t85, 0, 0, 0, 0, 0, 0, -g(1) * (t61 * t91 + t88) - g(2) * (-t61 * t87 + t92) - t76 * t95, -g(1) * (-t61 * t92 + t87) - g(2) * (t61 * t88 + t91) + t74 * t95, t56, -g(1) * (t101 * t75 + t84) - g(2) * t81 - g(3) * (t62 * pkin(4) + t61 * pkin(8) + t85) + (g(1) * t73 - g(2) * (t82 - t101)) * t77, 0, 0, 0, 0, 0, 0, -g(1) * (t61 * t93 + t90) - g(2) * (-t61 * t89 + t94) - t64 * t95, -g(1) * (-t61 * t94 + t89) - g(2) * (t61 * t90 + t93) + t63 * t95, t56, -g(1) * t84 - g(2) * t67 - g(3) * (t62 * t60 - t61 * t78 + t85) + (-g(1) * t79 - g(2) * t83) * t75 + (-g(1) * t83 - g(2) * (-t79 + t82)) * t77;];
U_reg  = t1;
