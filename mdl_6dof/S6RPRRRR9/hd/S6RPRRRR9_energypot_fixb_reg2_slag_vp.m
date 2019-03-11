% Calculate inertial parameters regressor of potential energy for
% S6RPRRRR9
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
% Datum: 2019-03-09 07:26
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U_reg = S6RPRRRR9_energypot_fixb_reg2_slag_vp(qJ, g, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(10,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRRRR9_energypot_fixb_reg2_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RPRRRR9_energypot_fixb_reg2_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RPRRRR9_energypot_fixb_reg2_slag_vp: pkin has to be [10x1] (double)');

%% Symbolic Calculation
% From energy_potential_fixb_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 07:25:17
% EndTime: 2019-03-09 07:25:18
% DurationCPUTime: 0.15s
% Computational Cost: add. (151->82), mult. (188->103), div. (0->0), fcn. (184->10), ass. (0->40)
t99 = g(3) * pkin(6);
t98 = pkin(2) + pkin(6);
t82 = -pkin(9) - pkin(8);
t76 = sin(qJ(4));
t97 = pkin(4) * t76;
t81 = cos(qJ(1));
t96 = g(2) * t81;
t80 = cos(qJ(3));
t95 = g(3) * t80;
t79 = cos(qJ(4));
t64 = t79 * pkin(4) + pkin(3);
t74 = -pkin(10) + t82;
t94 = t74 * t80;
t78 = sin(qJ(1));
t93 = t76 * t78;
t77 = sin(qJ(3));
t92 = t77 * t78;
t91 = t77 * t81;
t90 = t78 * t79;
t89 = t79 * t81;
t88 = t80 * t82;
t69 = t78 * pkin(7);
t70 = t78 * pkin(1);
t87 = t69 + t70;
t86 = t81 * pkin(1) + t78 * qJ(2);
t75 = qJ(4) + qJ(5);
t85 = t81 * pkin(7) + t86;
t84 = -t81 * qJ(2) + t70;
t83 = pkin(3) * t77 - pkin(8) * t80;
t60 = g(1) * t78 - t96;
t68 = qJ(6) + t75;
t66 = cos(t75);
t65 = sin(t75);
t63 = cos(t68);
t62 = sin(t68);
t61 = g(1) * t81 + g(2) * t78;
t59 = pkin(5) * t65 + t97;
t58 = pkin(5) * t66 + t64;
t57 = -g(3) * t77 + t60 * t80;
t1 = [0, 0, 0, 0, 0, 0, -t61, t60, -g(3), -t99, 0, 0, 0, 0, 0, 0, -g(3), t61, -t60, -g(1) * t86 - g(2) * t84 - t99, 0, 0, 0, 0, 0, 0, -t60 * t77 - t95, -t57, -t61, -g(1) * t85 - g(2) * (t69 + t84) - g(3) * t98, 0, 0, 0, 0, 0, 0, -g(1) * (t76 * t81 + t77 * t90) - g(2) * (-t77 * t89 + t93) - t79 * t95, -g(1) * (-t76 * t92 + t89) - g(2) * (t76 * t91 + t90) + t76 * t95, t57, -g(1) * (t83 * t78 + t85) - g(2) * t87 - g(3) * (pkin(3) * t80 + pkin(8) * t77 + t98) - (-qJ(2) - t83) * t96, 0, 0, 0, 0, 0, 0, -g(1) * (t65 * t81 + t66 * t92) - g(2) * (t65 * t78 - t66 * t91) - t66 * t95, -g(1) * (-t65 * t92 + t66 * t81) - g(2) * (t65 * t91 + t66 * t78) + t65 * t95, t57, -g(1) * (t64 * t92 + t78 * t88 + t85) - g(2) * (pkin(4) * t93 + t87) - g(3) * (t80 * t64 - t77 * t82 + t98) + (-g(1) * t97 - g(2) * (-t64 * t77 - qJ(2) - t88)) * t81, 0, 0, 0, 0, 0, 0, -g(1) * (t62 * t81 + t63 * t92) - g(2) * (t62 * t78 - t63 * t91) - t63 * t95, -g(1) * (-t62 * t92 + t63 * t81) - g(2) * (t62 * t91 + t63 * t78) + t62 * t95, t57, -g(1) * (t58 * t92 + t78 * t94 + t85) - g(2) * (t78 * t59 + t87) - g(3) * (t58 * t80 - t74 * t77 + t98) + (-g(1) * t59 - g(2) * (-t58 * t77 - qJ(2) - t94)) * t81;];
U_reg  = t1;
