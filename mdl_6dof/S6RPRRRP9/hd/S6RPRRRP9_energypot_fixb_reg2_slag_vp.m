% Calculate inertial parameters regressor of potential energy for
% S6RPRRRP9
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d3,d4,d5]';
% 
% Output:
% U_reg [1x(6*10)]
%   inertial parameter regressor of Potential energy

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 06:29
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U_reg = S6RPRRRP9_energypot_fixb_reg2_slag_vp(qJ, g, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(9,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRRRP9_energypot_fixb_reg2_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RPRRRP9_energypot_fixb_reg2_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S6RPRRRP9_energypot_fixb_reg2_slag_vp: pkin has to be [9x1] (double)');

%% Symbolic Calculation
% From energy_potential_fixb_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 06:28:30
% EndTime: 2019-03-09 06:28:30
% DurationCPUTime: 0.14s
% Computational Cost: add. (141->73), mult. (188->88), div. (0->0), fcn. (184->8), ass. (0->43)
t99 = g(3) * pkin(6);
t98 = pkin(2) + pkin(6);
t78 = -pkin(9) - pkin(8);
t77 = cos(qJ(1));
t97 = g(2) * t77;
t76 = cos(qJ(3));
t96 = g(3) * t76;
t72 = sin(qJ(4));
t95 = t72 * pkin(4);
t75 = cos(qJ(4));
t61 = t75 * pkin(4) + pkin(3);
t70 = -qJ(6) + t78;
t94 = t70 * t76;
t73 = sin(qJ(3));
t74 = sin(qJ(1));
t93 = t73 * t74;
t71 = qJ(4) + qJ(5);
t62 = sin(t71);
t92 = t74 * t62;
t63 = cos(t71);
t91 = t74 * t63;
t90 = t74 * t72;
t89 = t74 * t75;
t88 = t76 * t78;
t87 = t77 * t62;
t86 = t77 * t63;
t85 = t77 * t72;
t84 = t77 * t75;
t65 = t74 * pkin(7);
t66 = t74 * pkin(1);
t83 = t65 + t66;
t82 = t77 * pkin(1) + t74 * qJ(2);
t81 = t77 * pkin(7) + t82;
t80 = -t77 * qJ(2) + t66;
t79 = pkin(3) * t73 - pkin(8) * t76;
t59 = g(1) * t74 - t97;
t60 = g(1) * t77 + g(2) * t74;
t58 = pkin(5) * t62 + t95;
t57 = pkin(5) * t63 + t61;
t56 = -g(3) * t73 + t59 * t76;
t55 = -g(1) * (t73 * t91 + t87) - g(2) * (-t73 * t86 + t92) - t63 * t96;
t54 = -g(1) * (-t73 * t92 + t86) - g(2) * (t73 * t87 + t91) + t62 * t96;
t1 = [0, 0, 0, 0, 0, 0, -t60, t59, -g(3), -t99, 0, 0, 0, 0, 0, 0, -g(3), t60, -t59, -g(1) * t82 - g(2) * t80 - t99, 0, 0, 0, 0, 0, 0, -t59 * t73 - t96, -t56, -t60, -g(1) * t81 - g(2) * (t65 + t80) - g(3) * t98, 0, 0, 0, 0, 0, 0, -g(1) * (t73 * t89 + t85) - g(2) * (-t73 * t84 + t90) - t75 * t96, -g(1) * (-t73 * t90 + t84) - g(2) * (t73 * t85 + t89) + t72 * t96, t56, -g(1) * (t79 * t74 + t81) - g(2) * t83 - g(3) * (t76 * pkin(3) + t73 * pkin(8) + t98) - (-qJ(2) - t79) * t97, 0, 0, 0, 0, 0, 0, t55, t54, t56, -g(1) * (t61 * t93 + t74 * t88 + t81) - g(2) * (pkin(4) * t90 + t83) - g(3) * (t76 * t61 - t73 * t78 + t98) + (-g(1) * t95 - g(2) * (-t61 * t73 - qJ(2) - t88)) * t77, 0, 0, 0, 0, 0, 0, t55, t54, t56, -g(1) * (t57 * t93 + t74 * t94 + t81) - g(2) * (t74 * t58 + t83) - g(3) * (t76 * t57 - t73 * t70 + t98) + (-g(1) * t58 - g(2) * (-t57 * t73 - qJ(2) - t94)) * t77;];
U_reg  = t1;
