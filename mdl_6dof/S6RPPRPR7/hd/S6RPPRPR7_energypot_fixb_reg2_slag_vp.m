% Calculate inertial parameters regressor of potential energy for
% S6RPPRPR7
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d4,d6,theta3,theta5]';
% 
% Output:
% U_reg [1x(6*10)]
%   inertial parameter regressor of Potential energy

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 01:54
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U_reg = S6RPPRPR7_energypot_fixb_reg2_slag_vp(qJ, g, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(10,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPPRPR7_energypot_fixb_reg2_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RPPRPR7_energypot_fixb_reg2_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RPPRPR7_energypot_fixb_reg2_slag_vp: pkin has to be [10x1] (double)');

%% Symbolic Calculation
% From energy_potential_fixb_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 01:53:43
% EndTime: 2019-03-09 01:53:43
% DurationCPUTime: 0.16s
% Computational Cost: add. (150->72), mult. (168->87), div. (0->0), fcn. (160->10), ass. (0->41)
t70 = pkin(9) + qJ(4);
t62 = sin(t70);
t64 = cos(t70);
t101 = pkin(4) * t62 - qJ(5) * t64;
t100 = g(3) * pkin(6);
t99 = pkin(2) + pkin(6);
t72 = sin(pkin(9));
t98 = pkin(3) * t72;
t96 = g(3) * t64;
t69 = pkin(10) + qJ(6);
t61 = sin(t69);
t78 = cos(qJ(1));
t95 = t61 * t78;
t63 = cos(t69);
t94 = t63 * t78;
t71 = sin(pkin(10));
t93 = t71 * t78;
t73 = cos(pkin(10));
t92 = t73 * t78;
t77 = sin(qJ(1));
t91 = t77 * t61;
t90 = t77 * t63;
t89 = t77 * t71;
t88 = t77 * t73;
t87 = t78 * pkin(1) + t77 * qJ(2);
t74 = cos(pkin(9));
t85 = t74 * pkin(3) + t99;
t84 = t77 * t98 + t87;
t76 = -pkin(7) - qJ(3);
t83 = pkin(5) * t71 - t76;
t82 = -qJ(2) - t98;
t67 = t77 * pkin(1);
t81 = -t77 * t76 + t67;
t80 = -t78 * qJ(2) + t67;
t57 = g(1) * t77 - g(2) * t78;
t60 = pkin(5) * t73 + pkin(4);
t75 = -pkin(8) - qJ(5);
t79 = t60 * t62 + t64 * t75;
t58 = g(1) * t78 + g(2) * t77;
t56 = -g(3) * t62 + t57 * t64;
t1 = [0, 0, 0, 0, 0, 0, -t58, t57, -g(3), -t100, 0, 0, 0, 0, 0, 0, -g(3), t58, -t57, -g(1) * t87 - g(2) * t80 - t100, 0, 0, 0, 0, 0, 0, -g(3) * t74 - t57 * t72, g(3) * t72 - t57 * t74, -t58, -g(1) * (qJ(3) * t78 + t87) - g(2) * (t77 * qJ(3) + t80) - g(3) * t99, 0, 0, 0, 0, 0, 0, -t57 * t62 - t96, -t56, -t58, -g(1) * (-t78 * t76 + t84) - g(2) * (t82 * t78 + t81) - g(3) * t85, 0, 0, 0, 0, 0, 0, -g(1) * (t62 * t88 + t93) - g(2) * (-t62 * t92 + t89) - t73 * t96, -g(1) * (-t62 * t89 + t92) - g(2) * (t62 * t93 + t88) + t71 * t96, t56, -g(1) * (t101 * t77 + t84) - g(2) * t81 - g(3) * (pkin(4) * t64 + qJ(5) * t62 + t85) + (g(1) * t76 - g(2) * (t82 - t101)) * t78, 0, 0, 0, 0, 0, 0, -g(1) * (t62 * t90 + t95) - g(2) * (-t62 * t94 + t91) - t63 * t96, -g(1) * (-t62 * t91 + t94) - g(2) * (t62 * t95 + t90) + t61 * t96, t56, -g(1) * t84 - g(2) * t67 - g(3) * (t60 * t64 - t62 * t75 + t85) + (-g(1) * t79 - g(2) * t83) * t77 + (-g(1) * t83 - g(2) * (-t79 + t82)) * t78;];
U_reg  = t1;
